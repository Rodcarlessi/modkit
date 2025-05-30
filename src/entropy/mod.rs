use std::collections::{BTreeSet, HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::ops::Range;
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{anyhow, bail, Context};
use derive_new::new;
use itertools::{Itertools, MinMaxResult};
use log::{debug, info};
use nom::character::complete::multispace1;
use nom::IResult;
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{self, FetchDefinition, Read};
use rustc_hash::FxHashMap;

use crate::entropy::methylation_entropy::calc_me_entropy;
use crate::errs::{MkError, MkResult};
use crate::mod_bam::{BaseModCall, ModBaseInfo};
use crate::mod_base_code::{DnaBase, ModCodeRepr};
use crate::motifs::motif_bed::RegexMotif;
use crate::read_ids_to_base_mod_probs::{PositionModCalls, ReadBaseModProfile};
use crate::reads_sampler::sampling_schedule::ReferenceSequencesLookup;
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::thresholds::percentile_linear_interp;
use crate::util::{record_is_not_primary, ReferenceRecord, Strand};

mod methylation_entropy;
pub mod subcommand;
mod writers;

type BaseAndPosition = (DnaBase, u64);

#[derive(Debug)]
pub(super) enum GenomeWindow {
    CombineStrands {
        interval: Range<u64>,
        neg_to_pos_positions: FxHashMap<BaseAndPosition, BaseAndPosition>,
        read_patterns: Vec<Vec<BaseModCall>>,
        position_valid_coverages: Vec<u32>,
    },
    Stranded {
        // todo instead of having pos/neg for everything, make one struct and
        // have  an optional for all of it
        pos_interval: Option<Range<u64>>,
        neg_interval: Option<Range<u64>>,
        pos_positions: Option<Vec<BaseAndPosition>>,
        neg_positions: Option<Vec<BaseAndPosition>>,
        pos_read_patterns: Vec<Vec<BaseModCall>>,
        neg_read_patterns: Vec<Vec<BaseModCall>>,
        pos_position_valid_coverages: Vec<u32>,
        neg_position_valid_coverages: Vec<u32>,
    },
}

impl GenomeWindow {
    fn new_combine_strands(
        interval: Range<u64>,
        num_positions: usize,
        neg_to_pos_positions: FxHashMap<BaseAndPosition, BaseAndPosition>,
    ) -> Self {
        let position_valid_coverages = vec![0u32; num_positions];
        Self::CombineStrands {
            interval,
            neg_to_pos_positions,
            read_patterns: Vec::new(),
            position_valid_coverages,
        }
    }

    fn new_stranded(
        pos_positions: Option<Vec<BaseAndPosition>>,
        neg_positions: Option<Vec<BaseAndPosition>>,
        num_positions: usize,
    ) -> Self {
        let pos_interval = pos_positions.as_ref().map(|positions| {
            match positions.iter().map(|(_, p)| p).minmax() {
                MinMaxResult::MinMax(s, t) => *s..*t,
                MinMaxResult::OneElement(x) => *x..(*x + 1u64),
                MinMaxResult::NoElements => {
                    unreachable!("should have >0 elements")
                }
            }
        });
        let neg_interval = neg_positions.as_ref().map(|positions| {
            match positions.iter().map(|(_, p)| p).minmax() {
                MinMaxResult::MinMax(s, t) => *s..*t,
                MinMaxResult::OneElement(x) => *x..(*x + 1u64),
                MinMaxResult::NoElements => {
                    unreachable!("should have >0 elements")
                }
            }
        });

        #[cfg(debug_assertions)]
        let check = |positions: Option<&Vec<BaseAndPosition>>| {
            if let Some(ps) = positions {
                ps.iter().skip(1).fold(ps[0].1, |last, (_, next)| {
                    assert!(last < *next, "needs to be sorted");
                    *next
                });
            }
        };

        #[cfg(debug_assertions)]
        check(pos_positions.as_ref());
        #[cfg(debug_assertions)]
        check(neg_positions.as_ref());

        let pos_position_valid_coverages = vec![0u32; num_positions];
        let neg_position_valid_coverages = vec![0u32; num_positions];
        // debug!(
        //     "interval {pos_interval:?}, {neg_interval:?} \n\t> pos: \
        //      {pos_positions:?} neg {neg_positions:?}"
        // );
        Self::Stranded {
            pos_interval,
            neg_interval,
            pos_positions,
            neg_positions,
            pos_read_patterns: Vec::new(),
            neg_read_patterns: Vec::new(),
            pos_position_valid_coverages,
            neg_position_valid_coverages,
        }
    }

    #[inline]
    fn inc_coverage(&mut self, pos: usize, strand: &Strand) {
        match self {
            Self::CombineStrands { position_valid_coverages, .. } => {
                assert!(
                    pos < position_valid_coverages.len(),
                    "pos is larger than the window size?"
                );
                position_valid_coverages[pos] += 1;
            }
            Self::Stranded {
                pos_position_valid_coverages,
                neg_position_valid_coverages,
                ..
            } => match strand {
                Strand::Positive => {
                    assert!(
                        pos < pos_position_valid_coverages.len(),
                        "pos is larger than the window size?"
                    );
                    pos_position_valid_coverages[pos] += 1;
                }
                Strand::Negative => {
                    assert!(
                        pos < neg_position_valid_coverages.len(),
                        "pos is larger than the window size?"
                    );
                    neg_position_valid_coverages[pos] += 1;
                }
            },
        };
    }

    fn add_pattern(&mut self, strand: &Strand, pattern: Vec<BaseModCall>) {
        match self {
            Self::Stranded { pos_read_patterns, neg_read_patterns, .. } => {
                match strand {
                    Strand::Positive => pos_read_patterns.push(pattern),
                    Strand::Negative => neg_read_patterns.push(pattern),
                }
            }
            Self::CombineStrands { read_patterns, .. } => {
                read_patterns.push(pattern);
            }
        }
    }

    fn leftmost(&self) -> u64 {
        match (self.start(&Strand::Positive), self.start(&Strand::Negative)) {
            (Some(x), Some(y)) => std::cmp::min(x, y),
            (Some(x), None) => x,
            (None, Some(x)) => x,
            _ => unreachable!(
                "should always have either a positive or negative interval!"
            ),
        }
    }

    fn rightmost(&self) -> u64 {
        match (self.end(&Strand::Positive), self.end(&Strand::Negative)) {
            (Some(x), Some(y)) => std::cmp::max(x, y),
            (Some(x), None) => x,
            (None, Some(x)) => x,
            _ => unreachable!(
                "should always have either a positive or negative interval!"
            ),
        }
    }

    fn start(&self, strand: &Strand) -> Option<u64> {
        match self {
            Self::CombineStrands { interval, .. } => Some(interval.start),
            Self::Stranded { pos_interval, neg_interval, .. } => match strand {
                Strand::Positive => pos_interval.as_ref().map(|x| x.start),
                Strand::Negative => neg_interval.as_ref().map(|x| x.start),
            },
        }
    }

    fn end(&self, strand: &Strand) -> Option<u64> {
        match self {
            Self::CombineStrands { interval, .. } => Some(interval.end),
            Self::Stranded { pos_interval, neg_interval, .. } => match strand {
                Strand::Positive => pos_interval.as_ref().map(|x| x.end),
                Strand::Negative => neg_interval.as_ref().map(|x| x.end),
            },
        }
    }

    fn add_read_to_patterns(
        &mut self,
        ref_pos_to_basemod_call: &FxHashMap<BaseAndPosition, BaseModCall>,
        reference_start: i64,
        reference_end: i64,
        strand: Strand,
        max_filtered_positions: usize,
    ) {
        // check that the read fully covers the interval
        let reference_start = if reference_start >= 0 {
            Some(reference_start as u64)
        } else {
            None
        };
        let reference_end = if reference_start
            .map(|x| reference_end > x as i64)
            .unwrap_or(false)
        {
            Some(reference_end as u64)
        } else {
            None
        };

        let overlaps = reference_start
            .and_then(|s| reference_end.map(|t| (s, t)))
            .map(|(s, t)| match (self.start(&strand), self.end(&strand)) {
                (Some(wind_start), Some(wind_end)) => {
                    s <= wind_start && t >= wind_end
                }
                _ => false,
            })
            // .map(|(s, t)| s <= self.start() && t >= self.end())
            .unwrap_or(false);
        if !overlaps {
            return;
        }

        let pattern: Vec<BaseModCall> = match strand {
            Strand::Positive => match &self {
                Self::Stranded { pos_positions: Some(positions), .. } => {
                    positions
                        .iter()
                        .map(|p| {
                            ref_pos_to_basemod_call
                                .get(p)
                                .copied()
                                .unwrap_or(BaseModCall::Filtered)
                        })
                        .collect()
                }
                Self::CombineStrands { neg_to_pos_positions, .. } => {
                    neg_to_pos_positions
                        .values()
                        .map(|p| {
                            let call = ref_pos_to_basemod_call
                                .get(p)
                                .copied()
                                .unwrap_or(BaseModCall::Filtered);
                            (p, call)
                        })
                        .sorted_by(|((_, a), _), ((_, b), _)| a.cmp(b))
                        .map(|(_, call)| call)
                        .collect()
                }
                _ => return,
            },
            Strand::Negative => match &self {
                Self::Stranded { neg_positions: Some(positions), .. } => {
                    positions
                        .iter()
                        .map(|p| {
                            ref_pos_to_basemod_call
                                .get(p)
                                .copied()
                                .unwrap_or(BaseModCall::Filtered)
                        })
                        .collect()
                }
                Self::CombineStrands { neg_to_pos_positions, .. } => {
                    neg_to_pos_positions
                        .iter()
                        .map(|(neg_position, positive_position)| {
                            let call = ref_pos_to_basemod_call
                                .get(neg_position)
                                .copied()
                                .unwrap_or(BaseModCall::Filtered);
                            (positive_position, call)
                        })
                        .sorted_by(|((_, a), _), ((_, b), _)| a.cmp(b))
                        .map(|(_, call)| call)
                        .collect()
                }
                _ => return,
            },
        };

        if pattern.iter().filter(|&bmc| bmc == &BaseModCall::Filtered).count()
            > max_filtered_positions
        {
            // skip when too many filtered positions
            return;
        }

        for (i, call) in pattern.iter().enumerate() {
            match call {
                BaseModCall::Filtered => {}
                _ => self.inc_coverage(i, &strand),
            }
        }
        self.add_pattern(&strand, pattern);
    }

    fn get_mod_code_lookup(&self) -> FxHashMap<ModCodeRepr, char> {
        // looks complicated, but it just iterates over either the positive and
        // negative read patterns or the positive-combined read patterns
        let read_patterns: Box<dyn Iterator<Item = &Vec<BaseModCall>>> =
            match self {
                Self::Stranded {
                    pos_read_patterns, neg_read_patterns, ..
                } => {
                    Box::new(pos_read_patterns.iter().chain(neg_read_patterns))
                }
                Self::CombineStrands { read_patterns, .. } => {
                    Box::new(read_patterns.iter())
                }
            };

        // todo this could be done more simply with a set, but the idea is to
        // make  a single char code (e.g. '1', '2', '3', etc. for each
        // modification code
        read_patterns
            .flat_map(|pattern| {
                pattern.iter().filter_map(|call| match call {
                    BaseModCall::Modified(_, code) => Some(*code),
                    _ => None,
                })
            })
            .collect::<BTreeSet<ModCodeRepr>>()
            .into_iter()
            .enumerate()
            .map(|(id, code)| {
                // save 0 for canonical
                let id = id.saturating_add(1);
                let encoded = format!("{id}").parse::<char>().unwrap();
                (code, encoded)
            })
            .collect::<FxHashMap<ModCodeRepr, char>>()
    }

    fn encode_patterns(
        &self,
        chrom_id: u32,
        strand: Strand,
        patterns: &Vec<Vec<BaseModCall>>,
        mod_code_lookup: &FxHashMap<ModCodeRepr, char>,
        position_valid_coverages: &[u32],
        min_coverage: u32,
    ) -> MkResult<Vec<String>> {
        // todo remove these checks after testing
        assert!(
            self.start(&strand).is_some(),
            "start should be Some when encoding pattern for strand \
             {strand:?}, {patterns:?}"
        );
        assert!(
            self.end(&strand).is_some(),
            "end should be Some when encoding pattern for strand {strand:?}, \
             {patterns:?}"
        );

        if position_valid_coverages.iter().all(|x| *x >= min_coverage) {
            let encoded = patterns
                .iter()
                .map(|pat| {
                    let pattern = pat
                        .iter()
                        .map(|call| match call {
                            BaseModCall::Canonical(_) => '0',
                            BaseModCall::Modified(_, code) => {
                                *mod_code_lookup.get(code).unwrap()
                            }
                            BaseModCall::Filtered => '*',
                        })
                        .collect::<String>();
                    // todo remove after testing
                    assert_eq!(
                        pattern.len(),
                        position_valid_coverages.len(),
                        "pattern {pattern} is the wrong size? \
                         {position_valid_coverages:?}"
                    );
                    pattern
                })
                .collect();
            Ok(encoded)
        } else {
            let zero_coverage =
                position_valid_coverages.iter().all(|&cov| cov == 0);
            if zero_coverage {
                return Err(MkError::EntropyZeroCoverage {
                    chrom_id,
                    start: self.start(&strand).unwrap(),
                    end: self.end(&strand).unwrap(),
                });
            } else {
                let err = MkError::EntropyInsufficientCoverage {
                    chrom_id,
                    start: self.start(&strand).unwrap(),
                    end: self.end(&strand).unwrap(),
                };
                return Err(err);
            }
        }
    }

    fn into_entropy(
        &self,
        chrom_id: u32,
        min_valid_coverage: u32,
    ) -> WindowEntropy {
        let window_size = self.size();
        let constant = 1f32 / window_size as f32; // todo make this configurable

        let mod_code_lookup = self.get_mod_code_lookup();
        let positive_encoded_patterns = match &self {
            Self::CombineStrands {
                read_patterns,
                position_valid_coverages,
                ..
            } => Some(self.encode_patterns(
                chrom_id,
                Strand::Positive,
                read_patterns,
                &mod_code_lookup,
                position_valid_coverages,
                min_valid_coverage,
            )),
            Self::Stranded {
                pos_interval: Some(_),
                pos_read_patterns,
                pos_position_valid_coverages,
                ..
            } => Some(self.encode_patterns(
                chrom_id,
                Strand::Positive,
                pos_read_patterns,
                &mod_code_lookup,
                &pos_position_valid_coverages,
                min_valid_coverage,
            )),
            _ => None,
        };
        let negative_patterns = match &self {
            Self::Stranded {
                neg_interval: Some(_),
                neg_read_patterns,
                neg_position_valid_coverages,
                ..
            } => Some(self.encode_patterns(
                chrom_id,
                Strand::Negative,
                neg_read_patterns,
                &mod_code_lookup,
                neg_position_valid_coverages,
                min_valid_coverage,
            )),
            _ => None,
        };
        // left for debugging
        // debug!(
        //     "{}:{}-{} (+), {:?}",
        //     chrom,
        //     self.leftmost(),
        //     self.rightmost(),
        //     &positive_encoded_patterns
        // );
        // if let Some(nps) = negative_patterns.as_ref() {
        //     debug!(
        //         "{}:{}-{} (-), {:?}",
        //         chrom,
        //         self.leftmost(),
        //         self.rightmost(),
        //         &nps
        //     );
        // }

        // TODO: make sure there is a proper entropy test
        #[cfg(debug_assertions)]
        {
            if let Some(Ok(patterns)) = positive_encoded_patterns.as_ref() {
                debug_assert!(
                    patterns.iter().all(|x| x.len() == window_size),
                    "patterns are the wrong size {positive_encoded_patterns:?}"
                );
            }
            if let Some(Ok(neg_patterns)) = negative_patterns.as_ref() {
                debug_assert!(neg_patterns
                    .iter()
                    .all(|x| x.len() == window_size));
            }
        }

        let pos_me_entropy = positive_encoded_patterns.map(|maybe_patterns| {
            maybe_patterns.map(|patterns| {
                let me_entropy =
                    calc_me_entropy(&patterns, window_size, constant);
                let num_reads = patterns.len();
                let interval = self.start(&Strand::Positive).unwrap()
                    ..self.end(&Strand::Positive).unwrap().saturating_add(1);
                MethylationEntropy::new(me_entropy, num_reads, interval)
            })
        });

        let neg_me_entropy = negative_patterns.map(|maybe_patterns| {
            maybe_patterns.map(|patterns| {
                let me_entropy =
                    calc_me_entropy(&patterns, window_size, constant);
                let num_reads = patterns.len();
                let interval = self.start(&Strand::Negative).unwrap()
                    ..self.end(&Strand::Negative).unwrap().saturating_add(1);
                MethylationEntropy::new(me_entropy, num_reads, interval)
            })
        });

        WindowEntropy::new(chrom_id, pos_me_entropy, neg_me_entropy)
    }

    #[inline]
    fn size(&self) -> usize {
        match self {
            Self::Stranded { pos_position_valid_coverages, .. } => {
                pos_position_valid_coverages.len()
            }
            Self::CombineStrands { position_valid_coverages, .. } => {
                position_valid_coverages.len()
            }
        }
    }
}

pub(super) struct GenomeWindows {
    chrom_id: u32,
    entropy_windows: Vec<GenomeWindow>,
    region_name: Option<String>,
}

pub(super) enum EntropyCalculation {
    Windows(Vec<WindowEntropy>),
    Region(RegionEntropy),
}

impl GenomeWindows {
    fn new(
        chrom_id: u32,
        entropy_windows: Vec<GenomeWindow>,
        region_name: Option<String>,
    ) -> Self {
        assert!(!entropy_windows.is_empty());
        Self { chrom_id, entropy_windows, region_name }
    }

    fn get_range(&self) -> Range<u64> {
        // these expects are checked in a few places, make them .unwrap()s
        let start = self
            .entropy_windows
            .first()
            .expect("self.entropy_windows should not be empty")
            .leftmost();
        let end = self
            .entropy_windows
            .last()
            .expect("self.entropy_windows should not be empty")
            .rightmost();
        start..end
    }

    fn get_fetch_definition(&self) -> FetchDefinition {
        let range = self.get_range();
        let start = range.start as i64;
        let end = range.end as i64;
        let chrom_id = self.chrom_id;
        FetchDefinition::Region(chrom_id as i32, start, end)
    }

    fn into_entropy_calculation(
        self,
        chrom_id: u32,
        min_coverage: u32,
    ) -> EntropyCalculation {
        // to appease the bC we have to get the interval
        // here, but it's only used if we're summarizing a region
        let interval = self.get_range();
        let window_entropies = self
            .entropy_windows
            .par_iter()
            .map(|ew| ew.into_entropy(chrom_id, min_coverage))
            .collect::<Vec<_>>();
        let chrom_id = self.chrom_id;
        if let Some(region_name) = self.region_name {
            let mut pos_entropies = Vec::with_capacity(window_entropies.len());
            let mut pos_num_reads = Vec::with_capacity(window_entropies.len());
            let mut pos_num_fails = 0usize;
            let mut neg_entropies = Vec::with_capacity(window_entropies.len());
            let mut neg_num_reads = Vec::with_capacity(window_entropies.len());
            let mut neg_num_fails = 0usize;

            for window_entropy in window_entropies.iter() {
                match window_entropy.pos_me_entropy.as_ref() {
                    Some(Ok(me_entropy)) => {
                        pos_entropies.push(me_entropy.me_entropy);
                        pos_num_reads.push(me_entropy.num_reads);
                    }
                    Some(Err(_e)) => {
                        pos_num_fails += 1;
                    }
                    None => {}
                }
                match window_entropy.neg_me_entropy.as_ref() {
                    Some(Ok(me_entropy)) => {
                        neg_entropies.push(me_entropy.me_entropy);
                        neg_num_reads.push(me_entropy.num_reads);
                    }
                    Some(Err(_e)) => {
                        neg_num_fails += 1;
                    }
                    // this means it was combine strands
                    None => {}
                }
            }

            // todo make sure the semantics here are what I want,
            //  should pos_entropy_stats be an Option?
            let pos_entropy_stats = DescriptiveStats::new(
                &pos_entropies,
                &pos_num_reads,
                pos_num_fails,
                chrom_id,
                &interval,
            );
            // if neg_entropies is empty and there are no fails, we never saw
            // any negative strand me entropies
            let neg_entropy_stats = if neg_entropies.is_empty()
                && neg_num_fails == 0
            {
                assert!(
                    neg_num_reads.is_empty(),
                    "neg num reads and window entropies should both be empty"
                );
                None
            } else {
                // this will fail correctly if there are neg_entropies is empty
                // but there are fails
                Some(DescriptiveStats::new(
                    &neg_entropies,
                    &neg_num_reads,
                    neg_num_fails,
                    chrom_id,
                    &interval,
                ))
            };

            let region_entropy = RegionEntropy::new(
                chrom_id,
                interval,
                pos_entropy_stats,
                neg_entropy_stats,
                region_name,
                window_entropies,
            );
            EntropyCalculation::Region(region_entropy)
        } else {
            EntropyCalculation::Windows(window_entropies)
        }
    }
}

#[derive(new)]
struct MotifHit {
    pos: u64,
    neg_position: Option<u64>,
    strand: Strand,
    base: DnaBase,
}

struct SlidingWindows {
    motifs: Vec<RegexMotif>,
    work_queue: VecDeque<(ReferenceRecord, Vec<char>)>,
    region_names: VecDeque<String>,
    window_size: usize,
    num_positions: usize,
    batch_size: usize,
    curr_position: usize,
    curr_contig: ReferenceRecord,
    curr_seq: Vec<char>,
    curr_region_name: Option<String>,
    combine_strands: bool,
    /// the longest motif length, so we find motifs that are in the window, but
    /// reach outside the window
    motif_search_adj: usize,
    done: bool,
}

impl SlidingWindows {
    fn new_with_regions(
        reference_sequences_lookup: ReferenceSequencesLookup,
        regions_bed_fp: &PathBuf,
        motifs: Vec<RegexMotif>,
        combine_strands: bool,
        num_positions: usize,
        window_size: usize,
        batch_size: usize,
    ) -> anyhow::Result<Self> {
        let regions_iter =
            BufReader::new(File::open(regions_bed_fp).with_context(|| {
                format!("failed to load regions at {regions_bed_fp:?}")
            })?)
            .lines()
            // change the lines into Errors
            .map(|r| r.map_err(|e| anyhow!("failed to read line, {e}")))
            // Parse the lines
            .map(|r| r.and_then(|l| BedRegion::parse_str(&l)))
            // grab the subsequences, also collect up the errors for invalid BED
            // lines
            .map(|r| {
                r.and_then(|bed_region| {
                    let start = bed_region.interval.start;
                    let end = bed_region.interval.end;
                    let interval = start..end;
                    reference_sequences_lookup
                        .get_subsequence_by_name(
                            bed_region.chrom.as_str(),
                            interval,
                        )
                        .map(|seq| (bed_region, seq))
                })
            })
            .map_ok(|(bed_region, seq)| {
                let tid = reference_sequences_lookup
                    .name_to_chrom_id(bed_region.chrom.as_str())
                    .unwrap();
                let start = bed_region.interval.start as u32;
                let length = bed_region.length() as u32;
                let chrom_name = bed_region.chrom;
                let region_name = bed_region.name;
                let reference_record =
                    ReferenceRecord::new(tid, start, length, chrom_name);
                (reference_record, region_name, seq)
            });

        // accumulators for the above iterator, could have done this all in a
        // fold, but with 3 accumulators this is easier to look at and
        // ends up being the same thing
        let mut work_queue = VecDeque::new();
        let mut region_queue = VecDeque::new();
        let mut failures = HashMap::new();

        let mut add_failure = |cause: String| {
            *failures.entry(cause).or_insert(0) += 1;
        };

        for res in regions_iter {
            match res {
                Ok((reference_record, region_name, subseq)) => {
                    work_queue.push_back((reference_record, subseq));
                    region_queue.push_back(region_name);
                }
                Err(e) => {
                    add_failure(e.to_string());
                }
            }
        }

        if !failures.is_empty() {
            debug!("failure reasons while parsing regions BED file");
            for (cause, count) in
                failures.iter().sorted_by(|(_, a), (_, b)| a.cmp(b))
            {
                debug!("\t {cause}: {count}")
            }
        }

        if work_queue.is_empty() {
            bail!("no valid regions parsed");
        }

        assert_eq!(region_queue.len(), work_queue.len());
        let (curr_contig, curr_seq, curr_position, curr_region_name) = loop {
            let (ref_record, subseq, region_name) =
                match (work_queue.pop_front(), region_queue.pop_front()) {
                    (Some((rr, subseq)), Some(region_name)) => {
                        anyhow::Ok((rr, subseq, region_name))
                    }
                    _ => bail!(
                        "didn't find at least 1 sequence with valid start \
                         position"
                    ),
                }?;
            if let Some(start_position) =
                Self::find_start_position(&subseq, &motifs)
            {
                info!(
                    "starting with region {region_name} at 0-based position \
                     {} on contig {}",
                    start_position + ref_record.start as usize,
                    &ref_record.name
                );
                break (ref_record, subseq, start_position, region_name);
            } else {
                info!("region {region_name} has no valid positions, skipping");
                continue;
            }
        };
        debug!(
            "parsed {} regions, starting with {} on contig {}",
            region_queue.len() + 1usize,
            &curr_region_name,
            curr_contig.name
        );
        let motif_search_adj = motifs
            .iter()
            .map(|motif| motif.length())
            .filter(|l| *l > 1)
            .max()
            .unwrap_or(0);

        Ok(Self {
            motifs,
            work_queue,
            region_names: region_queue,
            window_size,
            num_positions,
            batch_size,
            curr_position,
            curr_contig,
            curr_seq,
            curr_region_name: Some(curr_region_name),
            combine_strands,
            motif_search_adj,
            done: false,
        })
    }

    fn new(
        reference_sequence_lookup: ReferenceSequencesLookup,
        motifs: Vec<RegexMotif>,
        combine_strands: bool,
        num_positions: usize,
        window_size: usize,
        batch_size: usize,
    ) -> anyhow::Result<Self> {
        let mut work_queue =
            reference_sequence_lookup.into_reference_sequences();

        let (curr_contig, curr_seq, curr_position) = loop {
            let (curr_record, curr_seq) =
                work_queue.pop_front().ok_or_else(|| {
                    anyhow!(
                        "didn't find at least 1 sequence with a valid start \
                         position"
                    )
                })?;
            if let Some(pos) = Self::find_start_position(&curr_seq, &motifs) {
                info!(
                    "starting with contig {} at 0-based position {pos}",
                    &curr_record.name
                );
                break (curr_record, curr_seq, pos);
            } else {
                info!(
                    "contig {} had no valid motif positions, skipping..",
                    curr_record.name
                );
            }
        };
        let motif_search_adj = motifs
            .iter()
            .map(|motif| motif.length())
            .filter(|l| *l > 1)
            .max()
            .unwrap_or(0);

        Ok(Self {
            motifs,
            work_queue,
            region_names: VecDeque::new(),
            window_size,
            num_positions,
            batch_size,
            curr_position,
            curr_contig,
            curr_seq,
            curr_region_name: None,
            combine_strands,
            motif_search_adj,
            done: false,
        })
    }

    #[inline]
    fn take_hits_if_enough(
        &self,
        motif_hits: &[MotifHit],
    ) -> Option<Vec<BaseAndPosition>> {
        let positions = motif_hits
            .into_iter()
            .take(self.num_positions)
            .map(|mh| (mh.base, mh.pos))
            .sorted_by(|(_, a), (_, b)| a.cmp(b))
            .collect::<Vec<BaseAndPosition>>();
        if positions.len() == self.num_positions {
            Some(positions)
        } else {
            None
        }
    }

    #[inline]
    fn enough_hits_for_window(
        &self,
        pos_hits: &[MotifHit],
        neg_hits: &[MotifHit],
    ) -> Option<GenomeWindow> {
        if self.combine_strands {
            let neg_to_pos = pos_hits
                .into_iter()
                .filter(|x| x.strand == Strand::Positive)
                .take(self.num_positions)
                .filter_map(|motif_hit| {
                    assert_eq!(
                        motif_hit.strand,
                        Strand::Positive,
                        "logic error!"
                    );
                    motif_hit.neg_position.map(|np| {
                        ((motif_hit.base, np), (motif_hit.base, motif_hit.pos))
                    })
                })
                .collect::<FxHashMap<BaseAndPosition, BaseAndPosition>>();
            if neg_to_pos.len() < self.num_positions {
                None
            } else {
                let (start, end) = match neg_to_pos
                    .keys()
                    .chain(neg_to_pos.values())
                    .map(|(_, x)| x)
                    .minmax()
                {
                    MinMaxResult::MinMax(s, t) => (*s, *t),
                    MinMaxResult::OneElement(x) => (*x, *x + 1u64), /* should probably fail here too? */
                    _ => unreachable!("there must be more than 1 element"),
                };
                let interval = start..end;
                Some(GenomeWindow::new_combine_strands(
                    interval,
                    self.num_positions,
                    neg_to_pos,
                ))
            }
        } else {
            if pos_hits.len() >= self.num_positions
                || neg_hits.len() >= self.num_positions
            {
                let pos_positions = self.take_hits_if_enough(pos_hits);
                let neg_positions = self.take_hits_if_enough(neg_hits);
                match (pos_positions, neg_positions) {
                    (Some(p), Some(n)) => {
                        assert_eq!(p.len(), self.num_positions);
                        assert!(!p.is_empty());
                        assert_eq!(n.len(), self.num_positions);
                        assert!(!n.is_empty());
                        let leftmost_positive_ref_pos = p
                            .iter()
                            .min_by(|(_, a), (_, b)| a.cmp(b))
                            .map(|(_, p)| *p)
                            .unwrap();
                        let leftmost_negative_ref_pos = n
                            .iter()
                            .min_by(|(_, a), (_, b)| a.cmp(b))
                            .map(|(_, p)| *p)
                            .unwrap();
                        if leftmost_positive_ref_pos < leftmost_negative_ref_pos
                        {
                            // debug!("(+) is lefter, using {p:?}");
                            Some(GenomeWindow::new_stranded(
                                Some(p),
                                None,
                                self.num_positions,
                            ))
                        } else if leftmost_negative_ref_pos
                            < leftmost_positive_ref_pos
                        {
                            // debug!("(-) is lefter, using {n:?}");
                            Some(GenomeWindow::new_stranded(
                                None,
                                Some(n),
                                self.num_positions,
                            ))
                        } else {
                            assert_eq!(
                                leftmost_positive_ref_pos,
                                leftmost_negative_ref_pos
                            );
                            // debug!("they are the same, using {p:?} and
                            // {n:?}");
                            Some(GenomeWindow::new_stranded(
                                Some(p),
                                Some(n),
                                self.num_positions,
                            ))
                        }
                    }
                    (Some(p), None) => {
                        // debug!("(+) only, using {p:?}");
                        Some(GenomeWindow::new_stranded(
                            Some(p),
                            None,
                            self.num_positions,
                        ))
                    }
                    (None, Some(n)) => {
                        // debug!("(-) only, using {n:?}");
                        Some(GenomeWindow::new_stranded(
                            None,
                            Some(n),
                            self.num_positions,
                        ))
                    }
                    _ => None,
                }
            } else {
                None
            }
        }
    }

    fn next_window(&mut self) -> Option<GenomeWindow> {
        while !self.at_end_of_contig() {
            // search forward for hits
            let end = std::cmp::min(
                self.curr_position.saturating_add(self.window_size),
                self.curr_seq.len(),
            );
            // todo optimize?
            // debug!(
            //     "genome space position at top {}, {}, {}",
            //     self.curr_position + self.curr_contig.start as usize,
            //     self.curr_position,
            //     self.motif_search_adj
            // );
            let subseq_start =
                self.curr_position.saturating_sub(self.motif_search_adj);
            let offset = self.curr_position.checked_sub(subseq_start).expect(
                "curr_position should always be greater than subset_start",
            );
            let subseq = self.curr_seq[subseq_start..end]
                .iter()
                .map(|x| *x)
                .collect::<String>();
            // debug!("subseq at the top {subseq}");
            // N.B. the 'position' in these tuples are  _genome coordinates_!
            // this is because when we fetch reads we need to do it with the
            // proper genome coordinates. when we're using normal
            // sliding windows, the relative coordinates and the
            // genome coordinates _should_ be the same however when
            // using regions, we slice the reference genome, so the
            // relative (to the sequence) and genome coordinates will _not_ be
            // the same
            let (pos_hits, neg_hits): (Vec<MotifHit>, Vec<MotifHit>) = self
                .motifs
                .iter()
                .flat_map(|motif| {
                    motif
                        .find_hits(&subseq)
                        .into_iter()
                        // this filter removes positions found before
                        // self.curr-position
                        .filter_map(|(pos, strand)| {
                            pos.checked_sub(offset).map(|p| (p, strand))
                        })
                        .map(|(pos, strand)| {
                            let adjusted_position = pos
                                .saturating_add(self.curr_position)
                                .saturating_add(
                                    self.curr_contig.start as usize,
                                );
                            let dna_base = DnaBase::parse(
                                self.curr_seq[pos + self.curr_position],
                            )
                            .unwrap();
                            let base = if strand == Strand::Negative {
                                dna_base.complement()
                            } else {
                                dna_base
                            };
                            let neg_position = motif
                                .motif_info
                                .negative_strand_position(
                                    adjusted_position as u32,
                                )
                                .map(|x| x as u64);
                            MotifHit::new(
                                adjusted_position as u64,
                                neg_position,
                                strand,
                                base,
                            )
                        })
                        .collect::<Vec<MotifHit>>()
                })
                .sorted_by(|a, b| a.pos.cmp(&b.pos))
                .partition(|x| x.strand == Strand::Positive);
            if let Some(entropy_window) =
                self.enough_hits_for_window(&pos_hits, &neg_hits)
            {
                let new_genome_space_position =
                    (entropy_window.leftmost() as usize).saturating_add(1usize);
                // info!("new genome position {new_genome_space_position}");
                // need to re-adjust to relative coordinates instead of genome
                // coordinates
                self.curr_position = new_genome_space_position
                    .checked_sub(self.curr_contig.start as usize)
                    .expect(
                        "should be able to subtract contig start from position",
                    );

                return Some(entropy_window);
            } else {
                // not enough on (+) or (-)
                let hits = pos_hits
                    .into_iter()
                    .chain(neg_hits)
                    .map(|mh| mh.pos as usize)
                    .map(|p| {
                        // need to re-adjust to relative coordinates instead of
                        // genome coordinates
                        p.checked_sub(self.curr_contig.start as usize)
                            .expect("should be able to re-adjust position")
                    })
                    .collect::<BTreeSet<usize>>();
                if let Some(&first) = hits.first() {
                    // at least 1
                    if self.curr_position == first {
                        match hits.iter().nth(1) {
                            Some(&second_hit) => {
                                self.curr_position = second_hit
                            }
                            None => {
                                // there was only 1
                                self.curr_position = end;
                            }
                        }
                    } else {
                        self.curr_position = first;
                    }
                } else {
                    // hits was empty, set to end
                    self.curr_position = end;
                }
                continue;
            }
        }
        None
    }

    fn find_start_position(
        seq: &[char],
        motifs: &[RegexMotif],
    ) -> Option<usize> {
        seq.par_chunks(10_000).find_map_first(|c| {
            let s = c.iter().collect::<String>();
            let min_pos = motifs
                .iter()
                .flat_map(|motif| {
                    motif.find_hits(&s).into_iter().nth(0).map(|(pos, _)| pos)
                })
                .min();
            min_pos
        })
    }

    #[inline]
    fn at_end_of_contig(&self) -> bool {
        self.curr_position >= self.curr_contig.length as usize
    }

    fn update_current_contig(&mut self) {
        'search: loop {
            if let Some((record, seq)) = self.work_queue.pop_front() {
                match Self::find_start_position(&seq, &self.motifs) {
                    Some(start_pos) => {
                        self.curr_contig = record;
                        self.curr_position = start_pos;
                        self.curr_seq = seq;
                        let region_name = self.region_names.pop_front();
                        self.curr_region_name = region_name;
                        break 'search;
                    }
                    None => {
                        if let Some(region_name) = self.region_names.pop_front()
                        {
                            debug!(
                                "skipping region {region_name}, no valid \
                                 positions for motifs {:?}",
                                &self.motifs
                            )
                        } else {
                            debug!(
                                "skipping {}, no valid positions for motifs \
                                 {:?}",
                                &record.name, &self.motifs
                            )
                        }
                        continue;
                    }
                }
            } else {
                assert!(self.region_names.is_empty());
                self.done = true;
                break 'search;
            }
        }
    }

    pub(super) fn total_length(&self) -> usize {
        self.work_queue.iter().map(|(_, s)| s.len()).sum::<usize>()
            + self.curr_seq.len()
    }
}

impl Iterator for SlidingWindows {
    type Item = Vec<GenomeWindows>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut batch = Vec::with_capacity(self.batch_size);
        let mut windows = Vec::new();
        loop {
            // stopping conditions
            if self.done || batch.len() >= self.batch_size {
                break;
            }

            // grab the next window
            if let Some(entropy_window) = self.next_window() {
                windows.push(entropy_window);
            }

            // update conditions
            if self.at_end_of_contig() {
                // need to rotate the windows since we're moving on to another
                // contig
                let finished_windows =
                    std::mem::replace(&mut windows, Vec::new());
                let finished_region =
                    std::mem::replace(&mut self.curr_region_name, None);
                if !finished_windows.is_empty() {
                    let entropy_windows = GenomeWindows::new(
                        self.curr_contig.tid,
                        finished_windows,
                        finished_region,
                    );
                    batch.push(entropy_windows);
                }
                self.update_current_contig();
                continue;
            }

            // N.B. semantics, if the current region name is None, we're just
            // doing sliding windows over the genome and we can cut
            // this batch once the window size is the batch size.
            // otoh, if current_region is Some, we never cut the batch until
            // we've finished the contig so that an entire region ends up
            // in a single batch
            if self.curr_region_name.is_none()
                && windows.len() > self.batch_size
            {
                assert!(
                    self.region_names.is_empty(),
                    "region names should be empty here!"
                );
                let finished_windows =
                    std::mem::replace(&mut windows, Vec::new());
                if !finished_windows.is_empty() {
                    let entropy_windows = GenomeWindows::new(
                        self.curr_contig.tid,
                        finished_windows,
                        None,
                    );
                    batch.push(entropy_windows);
                }
            }
        }

        if !windows.is_empty() {
            assert!(
                self.region_names.is_empty(),
                "region names should be empty here also!"
            );
            let entropy_windows =
                GenomeWindows::new(self.curr_contig.tid, windows, None);
            batch.push(entropy_windows)
        }

        if batch.is_empty() {
            None
        } else {
            Some(batch)
        }
    }
}

#[derive(new, Debug)]
pub(super) struct MethylationEntropy {
    me_entropy: f32,
    num_reads: usize,
    interval: Range<u64>,
}

// todo make this an enum, one for regions
#[derive(new, Debug)]
pub(super) struct WindowEntropy {
    chrom_id: u32,
    pos_me_entropy: Option<MkResult<MethylationEntropy>>,
    neg_me_entropy: Option<MkResult<MethylationEntropy>>,
}

struct DescriptiveStats {
    mean_entropy: f32,
    median_entropy: f32,
    max_entropy: f32,
    min_entropy: f32,
    mean_num_reads: f32,
    max_num_reads: usize,
    min_num_reads: usize,
    failed_count: usize,
    successful_count: usize,
}

impl DescriptiveStats {
    fn mean(xs: &[f32]) -> f32 {
        xs.iter().sum::<f32>() / (xs.len() as f32)
    }

    fn new(
        measurements: &[f32],
        n_reads: &[usize],
        n_fails: usize,
        chrom_id: u32,
        interval: &Range<u64>,
    ) -> MkResult<Self> {
        if measurements.is_empty() {
            debug_assert!(
                n_reads.is_empty(),
                "measurements and reads should be empty together"
            );
            Err(MkError::EntropyZeroCoverage {
                chrom_id,
                start: interval.start,
                end: interval.end,
            })
        } else {
            debug_assert_eq!(
                measurements.len(),
                n_reads.len(),
                "measurements and n_reads should be the same length"
            );
            let mean_entropy = Self::mean(measurements);
            let median_entropy =
                percentile_linear_interp(measurements, 0.5f32)?;
            // safe because of above check
            let (min_entropy, max_entropy) = match measurements.iter().minmax()
            {
                MinMaxResult::OneElement(x) => (*x, *x),
                MinMaxResult::MinMax(m, x) => (*m, *x),
                MinMaxResult::NoElements => {
                    unreachable!("checked for empty above")
                }
            };

            let mean_num_reads = Self::mean(
                &n_reads.iter().map(|&x| x as f32).collect::<Vec<_>>(),
            );
            let (min_num_reads, max_num_reads) = match n_reads.iter().minmax() {
                MinMaxResult::OneElement(x) => (*x, *x),
                MinMaxResult::MinMax(m, x) => (*m, *x),
                MinMaxResult::NoElements => {
                    unreachable!("checked for empty above")
                }
            };

            let success_count = measurements.len();

            Ok(Self {
                mean_entropy,
                median_entropy,
                max_entropy,
                min_entropy,
                mean_num_reads,
                max_num_reads,
                min_num_reads,
                successful_count: success_count,
                failed_count: n_fails,
            })
        }
    }

    pub(super) fn to_row(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        strand: Strand,
        region_name: &str,
    ) -> String {
        use crate::util::TAB;

        format!(
            "\
            {chrom}{TAB}\
            {start}{TAB}\
            {end}{TAB}\
            {region_name}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}\n",
            self.mean_entropy,
            strand.to_char(),
            self.median_entropy,
            self.min_entropy,
            self.max_entropy,
            self.mean_num_reads,
            self.min_num_reads,
            self.max_num_reads,
            self.successful_count,
            self.failed_count
        )
    }
}

#[derive(new)]
pub(super) struct RegionEntropy {
    chrom_id: u32,
    interval: Range<u64>,
    pos_entropy_stats: MkResult<DescriptiveStats>,
    neg_entropy_stats: Option<MkResult<DescriptiveStats>>,
    region_name: String,
    window_entropies: Vec<WindowEntropy>,
}

#[derive(new)]
struct Message {
    mod_calls: FxHashMap<BaseAndPosition, BaseModCall>,
    reference_start: i64,
    reference_end: i64,
    strand: Strand,
    // _name: String,
}

fn process_bam_fp(
    bam_fp: &PathBuf,
    fetch_definition: FetchDefinition,
    caller: Arc<MultipleThresholdModCaller>,
    io_threads: usize,
) -> anyhow::Result<Vec<Message>> {
    let mut reader = bam::IndexedReader::from_path(bam_fp)?;
    reader.set_threads(io_threads)?;
    reader.fetch(fetch_definition)?;

    let record_iter = reader
        .records()
        .filter_map(|r| r.ok())
        .filter(|record| {
            !record.is_unmapped()
                && !(record_is_not_primary(&record) || record.seq_len() == 0)
        })
        .filter_map(|record| {
            String::from_utf8(record.qname().to_vec())
                .ok()
                .map(|name| (record, name))
        })
        .filter_map(|(record, name)| {
            match ModBaseInfo::new_from_record(&record) {
                Ok(modbase_info) => Some((modbase_info, record, name)),
                Err(run_error) => {
                    debug!(
                        "read {name}, failed to parse modbase info, \
                         {run_error}"
                    );
                    None
                }
            }
        });

    let mut messages = Vec::new();
    for (modbase_info, record, name) in record_iter {
        match ReadBaseModProfile::process_record(
            &record,
            &name,
            modbase_info,
            None,
            None,
            1,
        ) {
            Ok(profile) => {
                let position_calls = PositionModCalls::from_profile(&profile);
                let strands = position_calls
                    .iter()
                    .map(|p| p.mod_strand)
                    .collect::<HashSet<Strand>>();
                if strands.len() > 1 {
                    debug!("duplex not yet supported");
                } else {
                    let strand = if record.is_reverse() {
                        Strand::Negative
                    } else {
                        Strand::Positive
                    };
                    let mod_calls = position_calls
                        .into_iter()
                        .filter_map(|p| {
                            match (p.ref_position, p.alignment_strand) {
                                (Some(ref_pos), Some(aln_strand)) => {
                                    Some((p, ref_pos, aln_strand))
                                }
                                _ => None,
                            }
                        })
                        .map(|(p, ref_pos, _alignment_strand)| {
                            let mod_base_call = caller
                                .call(&p.canonical_base, &p.base_mod_probs);
                            ((p.canonical_base, ref_pos as u64), mod_base_call)
                        })
                        .collect::<FxHashMap<BaseAndPosition, BaseModCall>>();
                    let msg = Message::new(
                        mod_calls,
                        record.reference_start(),
                        record.reference_end(),
                        strand,
                    );
                    messages.push(msg);
                }
            }
            Err(e) => {
                debug!("read {name} failed to extract modbase info, {e}");
            }
        };
    }
    Ok(messages)
}

pub(super) fn process_entropy_window(
    mut entropy_windows: GenomeWindows,
    min_coverage: u32,
    max_filtered_positions: usize,
    io_threads: usize,
    caller: Arc<MultipleThresholdModCaller>,
    bam_fps: &[PathBuf],
) -> anyhow::Result<EntropyCalculation> {
    let bam_fp = &bam_fps[0];
    let reader = bam::IndexedReader::from_path(bam_fp)?;
    let chrom_id = entropy_windows.chrom_id;
    drop(reader);

    let results = bam_fps
        .into_par_iter()
        .map(|fp| {
            process_bam_fp(
                fp,
                entropy_windows.get_fetch_definition(),
                caller.clone(),
                io_threads,
            )
        })
        .collect::<Vec<anyhow::Result<Vec<Message>>>>();

    for message_result in results {
        match message_result {
            Ok(messages) => {
                for message in messages {
                    entropy_windows.entropy_windows.par_iter_mut().for_each(
                        |window| {
                            window.add_read_to_patterns(
                                &message.mod_calls,
                                message.reference_start,
                                message.reference_end,
                                message.strand,
                                max_filtered_positions,
                            )
                        },
                    );
                }
            }
            Err(e) => {
                debug!("failed to run bam {e}");
            }
        }
    }

    Ok(entropy_windows.into_entropy_calculation(chrom_id, min_coverage))
}

#[derive(new, Debug)]
struct BedRegion {
    chrom: String,
    interval: Range<usize>,
    name: String,
}

impl BedRegion {
    fn length(&self) -> usize {
        self.interval.end - self.interval.start
    }

    fn parser(raw: &str) -> IResult<&str, Self> {
        let n_parts = raw.split('\t').count();
        let (rest, chrom) = crate::parsing_utils::consume_string(raw)?;
        let (rest, start) = crate::parsing_utils::consume_digit(rest)?;
        let (rest, stop) = crate::parsing_utils::consume_digit(rest)?;
        let (rest, name) = if n_parts == 3 {
            (rest, format!("{chrom}:{start}-{stop}"))
        } else {
            let (rest, _leading_tab) = multispace1(rest)?;
            crate::parsing_utils::consume_string_spaces(rest)?
        };

        let interval = (start as usize)..(stop as usize);
        let this = Self { chrom, interval, name };
        Ok((rest, this))
    }

    fn parse_str(raw: &str) -> anyhow::Result<Self> {
        Self::parser(raw)
            .map_err(|e| anyhow!("failed to parse {raw} into BED3 line, {e}"))
            .and_then(|(_, this)| {
                if this.interval.end > this.interval.start {
                    Ok(this)
                } else {
                    bail!("end must be after start")
                }
            })
    }
}

#[cfg(test)]
mod entropy_mod_tests {
    use crate::entropy::BedRegion;

    #[test]
    fn test_bed_region_parsing() {
        let raw = "chr1\t100\t101\tfoo\n";
        let bed_region = BedRegion::parse_str(raw).expect("should parse");
        assert_eq!(&bed_region.chrom, "chr1");
        assert_eq!(bed_region.interval, 100usize..101);
        assert_eq!(&bed_region.name, "foo");
        let raw = "chr1\t100\t101\tfoo\t400\t.\tmorestuff\n";
        let bed_region = BedRegion::parse_str(raw).expect("should parse");
        assert_eq!(&bed_region.chrom, "chr1");
        assert_eq!(bed_region.interval, 100usize..101);
        assert_eq!(&bed_region.name, "foo");

        let raw = "chr20\t279148\t279507\tCpG: 39";
        let bed_region = BedRegion::parse_str(raw).expect("should parse");
        assert_eq!(&bed_region.chrom, "chr20");
        assert_eq!(bed_region.interval, 279148usize..279507);
        assert_eq!(&bed_region.name, "CpG: 39");
    }
}
