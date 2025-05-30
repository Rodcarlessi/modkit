use std::collections::{BTreeMap, HashMap, HashSet};
use std::ops::ControlFlow;

use bio::alphabets::dna::revcomp;
use derive_new::new;
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use log::debug;
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{self, Read, Records};
use rustc_hash::{FxHashMap, FxHashSet};

use crate::errs::{MkError, MkResult};
use crate::mod_bam::{
    prob_to_qual, BaseModCall, BaseModProbs, CollapseMethod, EdgeFilter,
    ModBaseInfo, SeqPosBaseModProbs, SkipMode, TrackingModRecordIter,
    WithModBaseInfos,
};
use crate::mod_base_code::{
    BaseAndState, BaseState, DnaBase, ModCodeRepr, ProbHistogram,
};
use crate::monoid::Moniod;
use crate::motifs::motif_bed::MotifPositionLookup;
use crate::position_filter::StrandedPositionFilter;
use crate::reads_sampler::record_sampler::{Indicator, RecordSampler};
use crate::record_processor::{RecordProcessor, WithRecords};
use crate::util::{
    self, get_aligned_pairs_forward, get_master_progress_bar,
    get_query_name_string, get_reference_mod_strand, get_ticker,
    record_is_primary, Kmer, Strand, MISSING_SYMBOL, TAB,
};

/// Read IDs mapped to their base modification probabilities, organized
/// by the canonical base. This data structure contains essentially all
/// of the same data as in the records themselves, but with the query
/// position and the alternative probabilities removed (i.e. it only has
/// the probability of the called modification).
pub(crate) struct ReadIdsToBaseModProbs {
    // mapping of read id to canonical base mapped to a vec
    // of base mod calls on that canonical base
    pub(crate) inner: HashMap<String, HashMap<DnaBase, Vec<BaseModProbs>>>,
}

impl ReadIdsToBaseModProbs {
    fn add_read_without_probs(&mut self, read_id: &str) {
        self.inner.entry(read_id.to_owned()).or_insert(HashMap::new());
    }

    fn add_mod_probs_for_read(
        &mut self,
        read_id: &str,
        canonical_base: DnaBase,
        mod_probs: Vec<BaseModProbs>,
    ) {
        self.inner
            .entry(read_id.to_owned())
            .or_insert(HashMap::new())
            .entry(canonical_base)
            .or_insert(Vec::new())
            .extend(mod_probs)
    }

    #[inline]
    /// Returns most likely probabilities for base modifications predicted for
    /// each canonical base.
    pub(crate) fn mle_probs_per_base(
        &self,
        suppress_progress: bool,
    ) -> HashMap<DnaBase, Vec<f32>> {
        let pb = get_master_progress_bar(self.inner.len());
        if suppress_progress {
            pb.set_draw_target(indicatif::ProgressDrawTarget::hidden())
        }
        pb.set_message("aggregating per-base modification probabilities");
        self.inner
            .par_iter()
            .progress_with(pb)
            .map(|(_, can_base_to_base_mod_probs)| {
                can_base_to_base_mod_probs
                    .iter()
                    .map(|(canonical_base, base_mod_probs)| {
                        let probs = base_mod_probs
                            .iter()
                            .map(|bmc| match bmc.argmax_base_mod_call() {
                                BaseModCall::Modified(f, _) => f,
                                BaseModCall::Canonical(f) => f,
                                BaseModCall::Filtered => {
                                    unreachable!(
                                        "argmax base mod call should not \
                                         return Filtered"
                                    )
                                }
                            })
                            .collect::<Vec<f32>>();
                        (*canonical_base, probs)
                    })
                    .collect::<HashMap<DnaBase, Vec<f32>>>()
            })
            .reduce(|| HashMap::zero(), |a, b| a.op(b))
    }

    /// return argmax probs for each mod-code
    pub(crate) fn mle_probs_per_base_mod(
        &self,
        suppress_progress: bool,
    ) -> HashMap<BaseAndState, Vec<f64>> {
        // todo(arand) should really aggregate per mod-code
        let pb = get_master_progress_bar(self.inner.len());
        if suppress_progress {
            pb.set_draw_target(indicatif::ProgressDrawTarget::hidden())
        }
        pb.set_message("aggregating per-mod probabilities");
        self.inner
            .par_iter()
            .progress_with(pb)
            .filter_map(|(_, base_mod_probs)| {
                let grouped = base_mod_probs
                    .iter()
                    .map(|(dna_base, base_mod_probs)| {
                        base_mod_probs
                            .iter()
                            // can make this .base_mod_call
                            .map(|bmc| match bmc.argmax_base_mod_call() {
                                BaseModCall::Modified(p, code) => (
                                    (*dna_base, BaseState::Modified(code)),
                                    p as f64,
                                ),
                                BaseModCall::Canonical(p) => (
                                    (
                                        *dna_base,
                                        BaseState::Canonical(*dna_base),
                                    ),
                                    p as f64,
                                ),
                                BaseModCall::Filtered => {
                                    unreachable!(
                                        "argmax base mod call should not \
                                         return Filtered"
                                    )
                                }
                            })
                            .fold(
                                HashMap::<BaseAndState, Vec<f64>>::new(),
                                |mut acc, (base, p)| {
                                    acc.entry(base)
                                        .or_insert(Vec::new())
                                        .push(p);
                                    acc
                                },
                            )
                    })
                    .reduce(|a, b| a.op(b));
                grouped
            })
            .reduce(|| HashMap::zero(), |a, b| a.op(b))
    }

    pub(crate) fn get_per_mod_histograms(
        &self,
        suppress_progress: bool,
    ) -> ProbHistogram {
        let base_state_probs = self.mle_probs_per_base_mod(suppress_progress);
        let pb = get_master_progress_bar(base_state_probs.len());
        pb.set_message("preparing histograms");
        let prob_counts = base_state_probs
            .into_par_iter()
            .progress_with(pb)
            .map(|(base_state, probs)| {
                let max_p = probs
                    .iter()
                    .copied()
                    .max_by(|a, b| a.partial_cmp(b).unwrap())
                    .unwrap();
                let counts = probs
                    .into_iter()
                    .map(|x| prob_to_qual(x as f32))
                    .counts()
                    .into_iter()
                    .collect::<BTreeMap<u8, usize>>();
                let max_q = counts.keys().max().unwrap();
                debug!("{base_state:?} {max_p} {max_q}");
                (base_state, counts)
            })
            .collect::<HashMap<BaseAndState, BTreeMap<u8, usize>>>();
        ProbHistogram::new(prob_counts)
    }

    pub(crate) fn seen(&self, record_name: &str) -> bool {
        self.inner.contains_key(record_name)
    }
}

impl Moniod for ReadIdsToBaseModProbs {
    fn zero() -> Self {
        Self { inner: HashMap::new() }
    }

    fn op(self, other: Self) -> Self {
        let mut this = self;
        this.op_mut(other);
        this
    }

    fn op_mut(&mut self, other: Self) {
        for (read_id, base_mod_calls) in other.inner {
            if self.inner.contains_key(&read_id) {
                continue;
            } else {
                self.inner.insert(read_id, base_mod_calls);
            }
        }
    }

    fn len(&self) -> usize {
        self.inner.len()
    }
}

impl RecordProcessor for ReadIdsToBaseModProbs {
    type Output = Self;

    fn process_records<T: Read>(
        records: Records<T>,
        with_progress: bool,
        mut record_sampler: RecordSampler,
        collapse_method: Option<&CollapseMethod>,
        edge_filter: Option<&EdgeFilter>,
        position_filter: Option<&StrandedPositionFilter<()>>,
        only_mapped: bool,
        allow_non_primary: bool,
        _cut: Option<u32>,
        _kmer_size: Option<usize>,
    ) -> anyhow::Result<Self::Output> {
        let spinner = if with_progress {
            Some(record_sampler.get_progress_bar())
        } else {
            None
        };
        let mod_base_info_iter = records
            .with_mod_base_info()
            .filter(|(record, _)| {
                if only_mapped || edge_filter.is_some() {
                    !record.is_unmapped()
                } else {
                    true
                }
            })
            .filter(|(record, _)| {
                if !allow_non_primary {
                    record_is_primary(record)
                } else {
                    true
                }
            });
        let mut read_ids_to_mod_base_probs = Self::zero();
        for (record, mod_base_info) in mod_base_info_iter {
            match record_sampler.ask() {
                Indicator::Use(token) => {
                    let record_name = get_query_name_string(&record);
                    let aligned_pairs = if only_mapped {
                        get_aligned_pairs_forward(&record)
                            .filter_map(|pair| pair.ok())
                            .collect::<FxHashMap<usize, u64>>()
                    } else {
                        FxHashMap::default()
                    };
                    if record_name.is_err() {
                        debug!("record name failed UTF-8 decode");
                        continue;
                    }
                    let record_name = record_name.unwrap();
                    if read_ids_to_mod_base_probs.seen(&record_name) {
                        debug!(
                            "record: {record_name}, already processed, \
                             consider de-duplicating alignments."
                        );
                        continue;
                    }
                    if mod_base_info.is_empty() {
                        // the current iterator should filter these out, but
                        // leaving this check
                        // here in case that changes..
                        // add count of unused/no calls
                        // debug!("record {record_name} contains no mod-base
                        // information");
                        read_ids_to_mod_base_probs
                            .add_read_without_probs(&record_name);
                        continue;
                    }

                    let (_, base_mod_probs_iter) =
                        mod_base_info.into_iter_base_mod_probs();
                    let mut added_probs_for_record = false;
                    for (dna_base, strand, seq_pos_base_mod_probs) in
                        base_mod_probs_iter
                    {
                        let canonical_base = match strand {
                            Strand::Positive => dna_base,
                            Strand::Negative => dna_base.complement(),
                        };

                        let seq_pos_base_mod_probs = seq_pos_base_mod_probs
                            .filter_positions(
                                edge_filter,
                                position_filter,
                                only_mapped,
                                &aligned_pairs,
                                strand,
                                &record,
                            );

                        // must stay such that mod_probs will not be empty if
                        // seq_pos_base_mod_probs
                        // is Some otherwise added_mod_probs_for_record should
                        // not be flipped to true
                        if let Some(seq_pos_base_mod_probs) =
                            seq_pos_base_mod_probs
                        {
                            let mod_probs = seq_pos_base_mod_probs
                                .pos_to_base_mod_probs
                                .into_iter()
                                .map(|(_q_pos, base_mod_probs)| {
                                    if let Some(method) = collapse_method {
                                        base_mod_probs.into_collapsed(method)
                                    } else {
                                        base_mod_probs
                                    }
                                })
                                .collect::<Vec<BaseModProbs>>();
                            read_ids_to_mod_base_probs.add_mod_probs_for_read(
                                &record_name,
                                canonical_base,
                                mod_probs,
                            );
                            added_probs_for_record = true;
                        } else {
                            // trace!("all base mod positions were removed by
                            // filtering \
                            //     for {record_name} and base
                            // {raw_canonical_base}");
                            continue;
                        }
                    }
                    if let Some(pb) = &spinner {
                        pb.inc(1);
                    }
                    if added_probs_for_record {
                        record_sampler.used(token);
                    }
                }
                Indicator::Skip => continue,
                Indicator::Done => break,
            }
        }

        if let Some(pb) = &spinner {
            pb.finish_and_clear();
        }

        Ok(read_ids_to_mod_base_probs)
    }
}

impl WithRecords for ReadIdsToBaseModProbs {
    fn size(&self) -> u64 {
        let s = self
            .inner
            .iter()
            .map(|(_, base_mod_calls)| {
                base_mod_calls.values().map(|vs| vs.len()).sum::<usize>()
            })
            .sum::<usize>();
        s as u64
    }

    fn num_reads(&self) -> usize {
        self.inner.len()
    }
}

#[derive(new, Debug)]
pub(crate) struct ModProfile {
    pub(crate) query_position: usize,
    pub(crate) ref_position: Option<i64>,
    pub(crate) num_soft_clipped_start: usize,
    pub(crate) num_soft_clipped_end: usize,
    pub(crate) read_length: usize,
    pub(crate) q_mod: f32,
    pub(crate) raw_mod_code: ModCodeRepr,
    pub(crate) q_base: u8,
    pub(crate) query_kmer: Kmer,
    pub(crate) mod_strand: Strand,
    pub(crate) alignment_strand: Option<Strand>,
    pub(crate) canonical_base: DnaBase,
    pub(crate) inferred: bool,
}

impl ModProfile {
    pub(crate) fn header(with_motifs: bool) -> String {
        let mut fields = vec![
            "read_id",
            "forward_read_position",
            "ref_position",
            "chrom",
            "mod_strand",
            "ref_strand",
            "ref_mod_strand",
            "fw_soft_clipped_start",
            "fw_soft_clipped_end",
            "alignment_start",
            "alignment_end",
            "read_length",
            "mod_qual",
            "mod_code",
            "base_qual",
            "ref_kmer",
            "query_kmer",
            "canonical_base",
            "modified_primary_base",
            "inferred",
            "flag",
        ];
        if with_motifs {
            fields.push("motifs")
        }
        fields.join(&TAB.to_string())
    }

    pub(crate) fn within_alignment(&self) -> bool {
        util::within_alignment(
            self.query_position,
            self.num_soft_clipped_start,
            self.num_soft_clipped_end,
            self.read_length,
        )
    }

    pub(crate) fn to_row(
        &self,
        read_id: &str,
        chrom_name: &str,
        tid: Option<u32>,
        alignment_start: Option<u64>,
        alignment_end: Option<u64>,
        reference_seqs: &HashMap<String, Vec<u8>>,
        flag: u16,
        motif_positions_lookup: Option<&MotifPositionLookup>,
        with_motifs: bool,
    ) -> String {
        let query_kmer = format!("{}", self.query_kmer);
        let motif_hits = motif_positions_lookup.and_then(|lu| {
            match (self.ref_position, tid, self.alignment_strand) {
                (Some(i), Some(tid), Some(strand)) if i > 0i64 => {
                    let pos = i as usize;
                    let motif_hits = lu.get_motif_hits(
                        tid,
                        pos,
                        get_reference_mod_strand(self.mod_strand, strand),
                    );
                    motif_hits
                }
                _ => None,
            }
        });

        let kmer_size = self.query_kmer.size;
        let ref_kmer = if let Some(ref_pos) = self.ref_position {
            if ref_pos < 0 {
                ".".to_string()
            } else {
                reference_seqs
                    .get(chrom_name)
                    .map(|s| {
                        Kmer::from_seq(s, ref_pos as usize, kmer_size)
                            .to_string()
                    })
                    .unwrap_or(".".to_string())
            }
        } else {
            ".".to_string()
        };
        let modified_primary_base = if self.mod_strand == Strand::Negative {
            self.canonical_base.complement().char()
        } else {
            self.canonical_base.char()
        };

        let _within_alignment = self.within_alignment();
        let mut s = format!(
            "\
            {read_id}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {chrom_name}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}",
            self.query_position,
            self.ref_position.unwrap_or(-1),
            self.mod_strand.to_char(),
            self.alignment_strand.map(|s| s.to_char()).unwrap_or('.'),
            self.alignment_strand
                .map(|s| get_reference_mod_strand(self.mod_strand, s).to_char())
                .unwrap_or('.'),
            self.num_soft_clipped_start,
            self.num_soft_clipped_end,
            alignment_start.map(|x| x as i64).unwrap_or(-1i64),
            alignment_end.map(|x| x as i64).unwrap_or(-1i64),
            self.read_length,
            self.q_mod,
            self.raw_mod_code,
            self.q_base,
            ref_kmer,
            query_kmer,
            self.canonical_base.char(),
            modified_primary_base,
            self.inferred,
            flag,
            // motif_hits.unwrap_or_else(|| MISSING_SYMBOL.to_string())
        );

        if with_motifs {
            s.push(TAB);
            if let Some(x) = motif_hits.as_ref() {
                s.push_str(x.as_str());
            } else {
                s.push_str(MISSING_SYMBOL);
            }
        }

        s.push_str("\n");
        s
    }
}

#[derive(new, Debug)]
pub(crate) struct ReadBaseModProfile {
    pub(crate) record_name: String,
    pub(crate) chrom_id: Option<u32>,
    pub(crate) flag: u16,
    pub(crate) alignment_start: Option<u64>,
    pub(crate) alignment_end: Option<u64>,
    pub(crate) profile: Vec<ModProfile>,
}

impl ReadBaseModProfile {
    #[cfg(test)]
    pub(crate) fn from_record(
        record: &bam::Record,
        collapse_method: Option<&CollapseMethod>,
        edge_filter: Option<&EdgeFilter>,
        kmer_size: usize,
    ) -> MkResult<Self> {
        let mod_base_info = ModBaseInfo::new_from_record(record)?;
        let record_name = get_query_name_string(record)?;
        Self::process_record(
            record,
            &record_name,
            mod_base_info,
            collapse_method,
            edge_filter,
            kmer_size,
        )
    }

    fn get_kmer_from_sequence(
        forward_sequence: &[u8],
        forward_position: usize,
        mod_strand: Strand,
        kmer_size: usize,
    ) -> Kmer {
        let kmer = Kmer::new(&forward_sequence, forward_position, kmer_size);
        if mod_strand == Strand::Negative {
            kmer.reverse_complement()
        } else {
            kmer
        }
    }

    #[inline]
    fn base_mod_probs_to_mod_profile(
        query_pos_forward: usize,
        primary_base: DnaBase,
        mod_strand: Strand,
        base_mod_probs: BaseModProbs,
        base_qual: u8,
        kmer: Kmer,
        read_length: usize,
        ref_pos: Option<i64>,
        alignment_strand: Option<Strand>,
        num_clip_start: usize,
        num_clip_end: usize,
    ) -> Vec<ModProfile> {
        let inferred = base_mod_probs.inferred_unmodified;
        base_mod_probs
            .iter_probs()
            .map(|(raw_mod_code, prob)| {
                ModProfile::new(
                    query_pos_forward,
                    ref_pos,
                    num_clip_start,
                    num_clip_end,
                    read_length,
                    *prob,
                    *raw_mod_code,
                    base_qual,
                    kmer,
                    mod_strand,
                    alignment_strand,
                    primary_base,
                    inferred,
                )
            })
            .collect::<Vec<ModProfile>>()
    }

    pub(crate) fn process_record(
        record: &bam::Record,
        record_name: &str,
        mod_base_info: ModBaseInfo,
        collapse_method: Option<&CollapseMethod>,
        edge_filter: Option<&EdgeFilter>,
        kmer_size: usize,
    ) -> MkResult<Self> {
        let read_length = record.seq_len();
        let (num_clip_start, num_clip_end) =
            match ReadsBaseModProfile::get_soft_clipped(&record) {
                Ok((sc_start, sc_end)) => {
                    if record.is_reverse() {
                        (sc_end, sc_start)
                    } else {
                        (sc_start, sc_end)
                    }
                }
                Err(e) => {
                    debug!(
                        "record: {record_name}, has improper CIGAR, {}",
                        e.to_string()
                    );
                    return Err(MkError::InvalidCigar);
                }
            };

        let (alignment_strand, chrom_tid) = if record.is_unmapped() {
            (None, None)
        } else {
            if record.is_reverse() {
                (Some(Strand::Negative), Some(record.tid() as u32))
            } else {
                (Some(Strand::Positive), Some(record.tid() as u32))
            }
        };
        // mapping of _forward_ query position to (alignment_qpos, ref_pos)
        let forward_query_pos_to_ref_pos = if record.is_unmapped() {
            HashMap::new()
        } else {
            record
                .aligned_pairs_full()
                .filter_map(|pair| {
                    match (pair[0], pair[1]) {
                        // aligned or insert (r_pos is None)
                        (Some(qpos), r_pos) => {
                            if qpos < 0 {
                                None
                            } else {
                                let qpos = qpos as usize;
                                if record.is_reverse() {
                                    // shouldn't _really_ need to perform this
                                    // checked_sub
                                    // but better to do it this way than to
                                    // panic when there
                                    // is some bug/invalid CIGAR in a dependency
                                    read_length
                                        .checked_sub(qpos + 1usize)
                                        // todo make sure you dont need to check
                                        // that r_pos is < 0
                                        .map(|qpos_adj| {
                                            (qpos_adj, (qpos, r_pos))
                                        })
                                } else {
                                    Some((qpos, (qpos, r_pos)))
                                }
                            }
                        }
                        // delete
                        (None, _) => None,
                    }
                })
                .collect::<HashMap<usize, (usize, Option<i64>)>>()
        };

        let quals = if record.is_reverse() {
            record.qual().to_vec().into_iter().rev().collect()
        } else {
            record.qual().to_vec()
        };
        let seq_len = record.seq_len();
        let forward_sequence = if record.is_reverse() {
            revcomp(record.seq().as_bytes())
        } else {
            record.seq().as_bytes()
        };
        let (_, iter) = mod_base_info.into_iter_base_mod_probs();
        let base_mod_probs_iter = iter
            .into_iter()
            .filter_map(|(base, strand, probs)| {
                let filtered = if let Some(edge_filter) = edge_filter {
                    let x = probs
                        .edge_filter_positions(edge_filter, record.seq_len());
                    if x.is_none() {
                        debug!(
                            "\
                        {record_name}: all positions for primary base {base} \
                             were removed by edge filter."
                        )
                    }
                    x
                } else {
                    Some(probs)
                };
                filtered.map(|probs| (base, strand, probs))
            })
            .map(|(base, strand, mut probs)| {
                if let Some(collapse_method) = collapse_method {
                    probs = probs.into_collapsed(collapse_method);
                }
                (base, strand, probs)
            });

        let mut mod_profiles = base_mod_probs_iter
            .flat_map(|(primary_base, mod_strand, seq_pos_base_mod_probs)| {
                seq_pos_base_mod_probs
                    .pos_to_base_mod_probs
                    .into_iter()
                    .flat_map(|(forward_pos, base_mod_probs)| {
                        let ref_pos = forward_query_pos_to_ref_pos
                            .get(&forward_pos)
                            .and_then(|(_query_aligned_pos, ref_pos)| *ref_pos);
                        let seq_kmer = Self::get_kmer_from_sequence(
                            &forward_sequence,
                            forward_pos,
                            mod_strand,
                            kmer_size,
                        );
                        let base_qual = quals
                            .get(forward_pos)
                            .map(|q| *q)
                            .unwrap_or_else(|| {
                                debug!(
                                    "record: {record_name}, didn't find base \
                                     quality for position {forward_pos}"
                                );
                                0u8
                            });
                        Self::base_mod_probs_to_mod_profile(
                            forward_pos,
                            primary_base,
                            mod_strand,
                            base_mod_probs,
                            base_qual,
                            seq_kmer,
                            seq_len,
                            ref_pos,
                            alignment_strand,
                            num_clip_start,
                            num_clip_end,
                        )
                    })
                    .collect::<Vec<ModProfile>>()
            })
            .collect::<Vec<ModProfile>>();
        mod_profiles.par_sort_by(|a, b| {
            if record.is_reverse() {
                b.query_position.cmp(&a.query_position)
            } else {
                a.query_position.cmp(&b.query_position)
            }
        });
        let flag = record.flags();
        let alignment_start = record.reference_start();
        let alignment_end = record.reference_end();
        let alignment_start = if alignment_start >= 0 {
            Some(alignment_start as u64)
        } else {
            None
        };
        let alignment_end =
            if alignment_end >= 0 { Some(alignment_end as u64) } else { None };

        Ok(Self {
            record_name: record_name.to_owned(),
            chrom_id: chrom_tid,
            flag,
            alignment_start,
            alignment_end,
            profile: mod_profiles,
        })
    }

    pub(crate) fn remove_inferred(self) -> Self {
        let profile =
            self.profile.into_iter().filter(|p| !p.inferred).collect();
        Self::new(
            self.record_name,
            self.chrom_id,
            self.flag,
            self.alignment_start,
            self.alignment_end,
            profile,
        )
    }

    fn primary_alignment(&self) -> bool {
        self.flag == 0 || self.flag == 16
    }

    fn unmapped_alignment(&self) -> bool {
        self.flag == 4
    }

    pub(crate) fn iter_profiles(
        &self,
    ) -> Box<dyn Iterator<Item = &ModProfile> + '_> {
        if self.unmapped_alignment() || self.primary_alignment() {
            Box::new(self.profile.iter())
        } else {
            Box::new(self.profile.iter().filter(|p| p.within_alignment()))
        }
    }
}

#[derive(new, Debug)]
pub(crate) struct ReadsBaseModProfile {
    pub(crate) profiles: Vec<ReadBaseModProfile>,
    pub(crate) num_skips: usize,
    pub(crate) num_fails: usize,
}

impl ReadsBaseModProfile {
    fn get_soft_clipped(record: &bam::Record) -> MkResult<(usize, usize)> {
        if record.is_unmapped() {
            return Ok((0, 0));
        }
        let cigar = &record.cigar().0;
        let sc_start = cigar.iter().try_fold(0, |acc, op| match op {
            Cigar::SoftClip(l) => ControlFlow::Continue(acc + *l),
            _ => ControlFlow::Break(acc),
        });
        let sc_end = cigar.iter().rev().try_fold(0, |acc, op| match op {
            Cigar::SoftClip(l) => ControlFlow::Continue(acc + *l),
            _ => ControlFlow::Break(acc),
        });

        match (sc_start, sc_end) {
            (ControlFlow::Break(s), ControlFlow::Break(e)) => {
                Ok((s as usize, e as usize))
            }
            _ => return Err(MkError::InvalidCigar),
        }
    }

    pub(crate) fn remove_inferred(self) -> Self {
        let profiles =
            self.profiles.into_iter().map(|p| p.remove_inferred()).collect();
        Self::new(profiles, self.num_skips, self.num_fails)
    }
}

impl Moniod for ReadsBaseModProfile {
    fn zero() -> Self {
        Self { profiles: Vec::new(), num_skips: 0, num_fails: 0 }
    }

    fn op(self, other: Self) -> Self {
        let seen = self
            .profiles
            .iter()
            .map(|p| p.record_name.as_str())
            .collect::<HashSet<&str>>();
        let to_add = other
            .profiles
            .into_iter()
            .filter(|p| !seen.contains(p.record_name.as_str()))
            .collect::<Vec<ReadBaseModProfile>>();
        drop(seen);
        let mut profiles = self.profiles;
        profiles.extend(to_add.into_iter());

        let num_skips = self.num_skips + other.num_skips;
        let num_fails = self.num_fails + other.num_fails;
        Self { profiles, num_skips, num_fails }
    }

    fn op_mut(&mut self, other: Self) {
        let seen = self
            .profiles
            .iter()
            .map(|p| p.record_name.as_str())
            .collect::<HashSet<&str>>();
        let to_add = other
            .profiles
            .into_iter()
            .filter(|p| !seen.contains(p.record_name.as_str()))
            .collect::<Vec<ReadBaseModProfile>>();
        drop(seen);
        self.profiles.extend(to_add.into_iter());

        self.num_skips += other.num_skips;
        self.num_fails += other.num_fails;
    }

    fn len(&self) -> usize {
        self.profiles.len()
    }
}

impl RecordProcessor for ReadsBaseModProfile {
    type Output = Self;

    fn process_records<T: Read>(
        records: Records<T>,
        with_progress: bool,
        mut record_sampler: RecordSampler,
        collapse_method: Option<&CollapseMethod>,
        edge_filter: Option<&EdgeFilter>,
        _position_filter: Option<&StrandedPositionFilter<()>>,
        _only_mapped: bool,
        allow_non_primary: bool,
        cut: Option<u32>,
        kmer_size: Option<usize>,
    ) -> anyhow::Result<Self::Output> {
        let mut mod_iter =
            TrackingModRecordIter::new(records, false, allow_non_primary);
        let mut agg = Vec::new();
        let mut seen = HashSet::new();
        let pb = if with_progress { Some(get_ticker()) } else { None };

        let mut n_fails = 0usize;
        for (record, record_name, modbase_info) in &mut mod_iter {
            if let Some(cut) = cut {
                if record.reference_start() < cut as i64 {
                    continue;
                }
            }

            match record_sampler.ask() {
                Indicator::Use(token) => {
                    match ReadBaseModProfile::process_record(
                        &record,
                        &record_name,
                        modbase_info,
                        collapse_method,
                        edge_filter,
                        kmer_size.unwrap_or(5),
                    ) {
                        Ok(read_base_mod_profile) => {
                            if seen.contains(&record_name) {
                                debug!(
                                    "record: {record_name}, added more than \
                                     once"
                                );
                            } else {
                                seen.insert(record_name);
                            }
                            agg.push(read_base_mod_profile);

                            if let Some(pb) = &pb {
                                pb.inc(1);
                            }
                            record_sampler.used(token);
                        }
                        Err(_) => n_fails += 1,
                    }
                }
                Indicator::Skip => continue,
                Indicator::Done => break,
            }
        }

        let num_failed = mod_iter.num_failed + n_fails;
        let num_skipped = mod_iter.num_skipped;

        Ok(ReadsBaseModProfile {
            profiles: agg,
            num_skips: num_skipped,
            num_fails: num_failed,
        })
    }
}

impl WithRecords for ReadsBaseModProfile {
    fn size(&self) -> u64 {
        self.profiles.iter().map(|p| p.profile.len() as u64).sum::<u64>()
    }

    fn num_reads(&self) -> usize {
        self.profiles.len()
    }
}

impl SeqPosBaseModProbs {
    fn filter_positions(
        self,
        edge_filter: Option<&EdgeFilter>,
        position_filter: Option<&StrandedPositionFilter<()>>,
        only_mapped: bool,
        aligned_pairs: &FxHashMap<usize, u64>,
        mod_strand: Strand,
        record: &bam::Record,
    ) -> Option<Self> {
        let read_length = record.seq_len();
        let read_can_be_trimmed = edge_filter
            .map(|ef| ef.read_can_be_trimmed(read_length))
            .unwrap_or(true);
        if !read_can_be_trimmed {
            return None;
        }

        let starting_positions = self.pos_to_base_mod_probs.len();
        let starting_skip_mode = self.get_skip_mode();
        let probs = self
            .pos_to_base_mod_probs
            .into_iter()
            .filter(|(q_pos, _)| {
                // use edge filter, if provided
                let edge_keep = edge_filter
                    .map(|ef| match ef.keep_position(*q_pos, read_length) {
                        Ok(b) => b,
                        Err(e) => {
                            let read_name = get_query_name_string(record)
                                .unwrap_or_else(|e| {
                                    format!(
                                        "UTF-8 DECODE ERROR, {}",
                                        e.to_string()
                                    )
                                });
                            debug!(
                                "{read_name}, error when trying to filter \
                                 edge positions, {}",
                                e.to_string()
                            );
                            false
                        }
                    })
                    .unwrap_or(true);

                // only mapped, if asked for
                let only_mapped_keep = if only_mapped {
                    aligned_pairs.contains_key(q_pos)
                } else {
                    true
                };

                // "bedtools intersect" keep only positions in interval
                let position_keep = match position_filter {
                    Some(position_filter) => aligned_pairs
                        .get(q_pos)
                        .map(|ref_pos| {
                            let reference_strand =
                                match (mod_strand, record.is_reverse()) {
                                    (Strand::Positive, false) => {
                                        Strand::Positive
                                    }
                                    (Strand::Positive, true) => {
                                        Strand::Negative
                                    }
                                    (Strand::Negative, false) => {
                                        Strand::Negative
                                    }
                                    (Strand::Negative, true) => {
                                        Strand::Positive
                                    }
                                };

                            position_filter.contains(
                                record.tid(),
                                *ref_pos,
                                reference_strand,
                            )
                        })
                        .unwrap_or(false),
                    None => true,
                };

                edge_keep && only_mapped_keep && position_keep
            })
            .collect::<FxHashMap<usize, BaseModProbs>>();
        if probs.is_empty() {
            None
        } else {
            let skip_mode = if probs.len() == starting_positions {
                starting_skip_mode
            } else {
                // change to Explicit if we filtered any positions. This
                // is a little unnecessary since this method is private to
                // this module and it does not write MM/ML tags, but just
                // in case that isn't always true.
                // N.B. If you filter out calls, you _must_ change to Explicit
                // mode otherwise the trimmed calls could be considered
                // canonical.
                SkipMode::Explicit
            };
            Some(Self::new(skip_mode, probs))
        }
    }
}

#[derive(new, Debug)]
pub(crate) struct PositionModCalls {
    pub(crate) query_position: usize,
    pub(crate) ref_position: Option<i64>,
    pub(crate) num_soft_clipped_start: usize,
    pub(crate) num_soft_clipped_end: usize,
    pub(crate) read_length: usize,
    pub(crate) base_mod_probs: BaseModProbs,
    pub(crate) q_base: u8,
    pub(crate) query_kmer: Kmer,
    pub(crate) mod_strand: Strand,
    pub(crate) alignment_strand: Option<Strand>,
    pub(crate) canonical_base: DnaBase,
}

impl PositionModCalls {
    pub(crate) fn from_profile(
        read_base_mod_profile: &ReadBaseModProfile,
    ) -> Vec<Self> {
        type Key = (usize, Strand, DnaBase);
        let (grouped, mod_codes): (
            FxHashMap<Key, Vec<&ModProfile>>,
            FxHashMap<DnaBase, FxHashSet<ModCodeRepr>>,
        ) = read_base_mod_profile.iter_profiles().fold(
            (FxHashMap::default(), FxHashMap::default()),
            |(mut acc, mut codes), x| {
                let k = (x.query_position, x.mod_strand, x.canonical_base);
                acc.entry(k).or_insert(Vec::new()).push(x);
                codes
                    .entry(x.canonical_base)
                    .or_insert_with(FxHashSet::default)
                    .insert(x.raw_mod_code);
                (acc, codes)
            },
        );
        let mod_codes = mod_codes
            .into_iter()
            .map(|(base, codes)| {
                (base, codes.into_iter().sorted().collect::<Vec<ModCodeRepr>>())
            })
            .collect::<FxHashMap<DnaBase, Vec<ModCodeRepr>>>();

        grouped
            .into_iter()
            .fold(
                Vec::<Self>::new(),
                |mut acc, ((query_pos, strand, base), mod_profile)| {
                    let codes_for_base = mod_codes.get(&base).unwrap();
                    let base_mod_probs =
                        if mod_profile.iter().any(|x| x.inferred) {
                            BaseModProbs::new_inferred_canonical(
                                codes_for_base.iter(),
                            )
                        } else {
                            let mut probs = mod_profile
                                .iter()
                                .map(|x| (x.raw_mod_code, x.q_mod))
                                .collect::<FxHashMap<ModCodeRepr, f32>>();
                            for code in codes_for_base.iter() {
                                if !probs.contains_key(&code) {
                                    probs.insert(*code, 0f32);
                                }
                            }

                            BaseModProbs::new(probs, false)
                        };
                    let template = &mod_profile[0];
                    let ref_position = template.ref_position;
                    let num_clip_start = template.num_soft_clipped_start;
                    let num_clip_end = template.num_soft_clipped_end;
                    let q_base = template.q_base;
                    let kmer = template.query_kmer;
                    let alignment_strand = template.alignment_strand;

                    let pos_mod_calls = PositionModCalls::new(
                        query_pos,
                        ref_position,
                        num_clip_start,
                        num_clip_end,
                        template.read_length,
                        base_mod_probs,
                        q_base,
                        kmer,
                        strand,
                        alignment_strand,
                        base,
                    );
                    acc.push(pos_mod_calls);

                    acc
                },
            )
            .into_iter()
            .sorted_by(|a, b| {
                if a.alignment_strand
                    .map(|s| s == Strand::Negative)
                    .unwrap_or(false)
                {
                    b.query_position.cmp(&a.query_position)
                } else {
                    a.query_position.cmp(&b.query_position)
                }
            })
            .collect()
    }

    pub(crate) fn within_alignment(&self) -> bool {
        util::within_alignment(
            self.query_position,
            self.num_soft_clipped_start,
            self.num_soft_clipped_end,
            self.read_length,
        )
    }
}

#[cfg(test)]
mod read_ids_to_base_mod_probs_tests {
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::Path;

    use rust_htslib::bam::{self, Read};
    use rustc_hash::{FxHashMap, FxHashSet};

    use crate::mod_bam::filter_records_iter;
    use crate::position_filter::StrandedPositionFilter;
    use crate::util::get_aligned_pairs_forward;

    #[test]
    fn test_seq_pos_base_mod_probs_filter_positions() {
        let mut reader = bam::Reader::from_path(
            "tests/resources/bc_anchored_10_reads.sorted.bam",
        )
        .unwrap();
        let header = reader.header().to_owned();
        let records = reader.records();
        let chrom_to_tid = (0..header.target_count())
            .map(|tid| {
                (String::from_utf8(header.tid2name(tid).to_vec()).unwrap(), tid)
            })
            .collect::<HashMap<String, u32>>();

        let position_bed_fp = "tests/resources/CGI_ladder_3.6kb_ref_CG.bed";
        let position_filter = StrandedPositionFilter::from_bed_file(
            &Path::new(position_bed_fp).to_path_buf(),
            &chrom_to_tid.iter().map(|(k, v)| (k.as_str(), *v)).collect(),
            true,
        )
        .unwrap();

        let mut pos_positions = FxHashSet::default();
        let mut neg_positions = FxHashSet::default();
        for line in BufReader::new(File::open(position_bed_fp).unwrap())
            .lines()
            .map(|r| r.unwrap())
        {
            let parts = line.split_whitespace().collect::<Vec<&str>>();
            assert_eq!(parts.len(), 6);
            if parts[0] != "oligo_1512_adapters" {
                continue;
            }
            // assert_eq!(parts[0], );
            let pos = parts[1].parse::<u64>().unwrap();
            match parts[5] {
                "+" => assert!(pos_positions.insert(pos)),
                "-" => assert!(neg_positions.insert(pos)),
                _ => panic!("illegal strand in BED"),
            }
        }

        let mod_base_info_iter = filter_records_iter(records);
        for (record, mod_base_info) in mod_base_info_iter {
            let aligned_pairs = get_aligned_pairs_forward(&record)
                .filter_map(|pair| pair.ok())
                .collect::<FxHashMap<usize, u64>>();

            let (_converters, base_mod_probs_iter) =
                mod_base_info.into_iter_base_mod_probs();
            for (_primary_base, mod_strand, seq_pos_mod_base_probs) in
                base_mod_probs_iter
            {
                let seq_pos_mod_base_probs = seq_pos_mod_base_probs
                    .filter_positions(
                        None,
                        Some(&position_filter),
                        true,
                        &aligned_pairs,
                        mod_strand,
                        &record,
                    )
                    .unwrap();
                let positions_to_check = if record.is_reverse() {
                    &neg_positions
                } else {
                    &pos_positions
                };
                for (q_pos, _) in seq_pos_mod_base_probs.pos_to_base_mod_probs {
                    let r_pos = aligned_pairs.get(&q_pos).unwrap();
                    assert!(positions_to_check.contains(r_pos));
                }
            }
        }
    }
}
