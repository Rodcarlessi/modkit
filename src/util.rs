use std::collections::{HashMap, HashSet};
use std::fmt::{Debug, Display, Formatter};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::str;

use anyhow::{anyhow, bail};
use anyhow::{Context, Result as AnyhowResult};
use bio::alphabets::dna::complement;
use clap::ValueEnum;
use derive_new::new;
use indexmap::IndexMap;
use indicatif::{ProgressBar, ProgressStyle};
use itertools::Itertools;
use lazy_static::lazy_static;
use linear_map::LinearMap;
use log::{debug, error, info};
use nom::bytes::complete::tag;
use nom::character::complete::one_of;
use nom::combinator::map_res;
use nom::multi::many0;
use nom::IResult;
use prettytable::row;
use regex::Regex;
use rust_htslib::bam::{
    self, ext::BamRecordExtensions, header::HeaderRecord, record::Aux,
    HeaderView, Read,
};
use rustc_hash::FxHashMap;
use substring::Substring;

use crate::errs::{MkError, MkResult};
use crate::mod_base_code::{DnaBase, ParseChar};
use crate::monoid::Moniod;
use crate::parsing_utils::{
    consume_char, consume_digit, consume_dot, consume_float, consume_string,
    consume_string_spaces,
};

pub(crate) const TAB: char = '\t';
pub(crate) const MISSING_SYMBOL: &'static str = ".";

pub(crate) fn create_out_directory<T: AsRef<std::ffi::OsStr>>(
    raw_path: T,
) -> anyhow::Result<()> {
    if let Some(p) = Path::new(&raw_path).parent() {
        if !p.exists() && p != Path::new("") {
            info!("creating directory at {p:?}");
            std::fs::create_dir_all(p)?;
        }
    }
    Ok(())
}

pub(crate) fn get_ticker() -> ProgressBar {
    let ticker = ProgressBar::new_spinner();
    ticker.set_style(ProgressStyle::with_template("> {pos} {msg}").unwrap());
    ticker
}
pub(crate) fn get_spinner() -> ProgressBar {
    let spinner = ProgressBar::new_spinner();
    spinner.set_style(
        ProgressStyle::with_template(
            "{spinner:.blue} [{elapsed_precise}] {pos} {msg}",
        )
        .unwrap()
        .tick_strings(&[
            "▹▹▹▹▹",
            "▸▹▹▹▹",
            "▹▸▹▹▹",
            "▹▹▸▹▹",
            "▹▹▹▸▹",
            "▹▹▹▹▸",
            "▪▪▪▪▪",
        ]),
    );
    spinner
}

fn get_master_progress_bar_style() -> ProgressStyle {
    ProgressStyle::with_template(
        "[{elapsed_precise}] {bar:40.green/yellow} {pos:>7}/{len:7} {msg}",
    )
    .unwrap()
    .progress_chars("##-")
}

fn get_subroutine_progress_bar_style() -> ProgressStyle {
    ProgressStyle::with_template(
        "[{elapsed_precise}] {bar:40.blue/cyan} {pos:>7}/{len:7} {msg}",
    )
    .unwrap()
    .progress_chars("##-")
}

fn get_gauge_style() -> ProgressStyle {
    ProgressStyle::with_template("{bar:40.red/blue} {pos:>7}/{len:7} {msg}")
        .unwrap()
        .progress_chars("||-")
}

pub(crate) fn get_master_progress_bar<
    T: num_traits::Num + num_traits::cast::AsPrimitive<u64>,
>(
    n: T,
) -> ProgressBar {
    ProgressBar::new(n.as_()).with_style(get_master_progress_bar_style())
}

pub(crate) fn get_subroutine_progress_bar<
    T: num_traits::Num + num_traits::cast::AsPrimitive<u64>,
>(
    n: T,
) -> ProgressBar {
    ProgressBar::new(n.as_()).with_style(get_subroutine_progress_bar_style())
}

pub(crate) fn get_guage<
    T: num_traits::Num + num_traits::cast::AsPrimitive<u64>,
>(
    n_max: T,
) -> ProgressBar {
    ProgressBar::new(n_max.as_()).with_style(get_gauge_style())
}

pub(crate) fn get_aligned_pairs_forward(
    record: &bam::Record,
) -> impl Iterator<Item = AnyhowResult<(usize, u64)>> + '_ {
    let read_length = record.seq_len();
    record.aligned_pairs().map(move |pair| {
        let q_pos = pair[0] as usize;
        let q_pos = if record.is_reverse() {
            read_length.checked_sub(q_pos).and_then(|x| x.checked_sub(1))
        } else {
            Some(q_pos)
        };
        if q_pos.is_none() || pair[1] < 0 {
            let read_id = get_query_name_string(&record)
                .unwrap_or("failed-to-parse-utf8".to_owned());
            debug!("record {read_id} has invalid aligned pair {:?}", pair);
            return Err(anyhow!("pair {:?} is invalid", pair));
        }

        let r_pos = pair[1];
        assert!(r_pos >= 0);
        let r_pos = r_pos as u64;
        Ok((q_pos.unwrap(), r_pos))
    })
}

pub(crate) fn get_query_name_string(record: &bam::Record) -> MkResult<String> {
    String::from_utf8(record.qname().to_vec())
        .map_err(|_e| MkError::InvalidRecordName)
}

#[inline]
pub(crate) fn get_forward_sequence(record: &bam::Record) -> Vec<u8> {
    if record.is_reverse() {
        bio::alphabets::dna::revcomp(record.seq().as_bytes())
    } else {
        record.seq().as_bytes()
    }
}

#[inline]
pub(crate) fn get_forward_sequence_str(
    record: &bam::Record,
) -> MkResult<String> {
    let seq_bs = get_forward_sequence(record);
    let seq = String::from_utf8(seq_bs)
        .map_err(|e| MkError::InvalidReadSequence(e))?;
    if seq.len() == 0 {
        return Err(MkError::EmptyReadSequence);
    }
    Ok(seq)
}

pub(crate) fn get_tag<T>(
    record: &bam::Record,
    tag_keys: &[&'static str; 2],
    parser: &dyn Fn(&Aux) -> MkResult<T>,
) -> MkResult<(T, &'static str)> {
    let tag_new = record.aux(tag_keys[0].as_bytes());
    let tag_old = record.aux(tag_keys[1].as_bytes());

    let (tag, t) = match (tag_new, tag_old) {
        (Ok(aux), _) => Ok((aux, tag_keys[0])),
        (Err(_), Ok(aux)) => Ok((aux, tag_keys[1])),
        _ => Err(MkError::AuxMissing),
    }?;
    parser(&tag).map(|v| (v, t))
}

pub(crate) fn parse_nm(record: &bam::Record) -> anyhow::Result<u32> {
    let nm_tag = record.aux("NM".as_bytes())?;
    match nm_tag {
        Aux::U8(x) => Ok(x as u32),
        Aux::U16(x) => Ok(x as u32),
        Aux::U32(x) => Ok(x),
        Aux::I8(x) => Ok(x as u32),
        Aux::I16(x) => Ok(x as u32),
        Aux::I32(x) => Ok(x as u32),
        _ => bail!("invalid NM tag {nm_tag:?}"),
    }
}

// Regex split into three possible elements
// (\d+) - matches
// (\^[A-Z]+) - deletions
// ([A-Z]) - mismatch
lazy_static! {
    pub static ref MDTAG_REGEX: Regex =
        Regex::new(r"(\d+)|(\^[A-Za-z]+)|([A-Za-z])").unwrap();
}

#[allow(dead_code)]
pub(crate) enum MdTag {
    // Number of matches
    Match(usize),
    // Mismatch base
    Mismatch(DnaBase),
    // Base deletions with length
    Deletion(Vec<DnaBase>),
}

// Parse BAM tags
// returns a vector of Option<MdTag> in the event the BAM tag has invalid
// elements
#[allow(dead_code)]
pub(crate) fn parse_md(record: &bam::Record) -> anyhow::Result<Vec<MdTag>> {
    let md_tag = record.aux("MD".as_bytes()).context("missing MD tag")?;
    let Aux::String(md_tag) = md_tag else { bail!("MD tag isn't a String") };

    MDTAG_REGEX
        .captures_iter(&md_tag)
        .map(|op| {
            if let Some(md_match) = op.get(1) {
                md_match
                    .as_str()
                    .parse::<usize>()
                    .map_err(|e| anyhow!("invalid match number, {e}"))
                    .map(|n| MdTag::Match(n))
            } else if let Some(md_deletion) = op.get(2) {
                md_deletion
                    .as_str()
                    .to_uppercase()
                    .chars()
                    .map(|b| DnaBase::parse_char(b).map_err(|e| e.into()))
                    .collect::<anyhow::Result<Vec<DnaBase>>>()
                    .map(|bases| MdTag::Deletion(bases))
            } else if let Some(md_mismatch) = op.get(3) {
                md_mismatch
                    .as_str()
                    .parse::<char>()
                    .map_err(|e| anyhow!("invalid mismatch char, {e}"))
                    .and_then(|b| DnaBase::parse_char(b).map_err(|e| e.into()))
                    .map(|b| MdTag::Mismatch(b))
            } else {
                bail!("invalid MD, should match one of the groups")
            }
        })
        .collect::<anyhow::Result<Vec<MdTag>>>()
}

#[derive(Debug, PartialEq, Eq, Copy, Clone, Hash, Default, PartialOrd, Ord)]
pub enum Strand {
    #[default]
    Positive,
    Negative,
}

impl Strand {
    pub fn parse_char(x: char) -> MkResult<Self> {
        match x {
            '+' => Ok(Self::Positive),
            '-' => Ok(Self::Negative),
            _ => Err(MkError::InvalidStrand),
        }
    }
    pub fn to_char(&self) -> char {
        match self {
            Strand::Positive => '+',
            Strand::Negative => '-',
        }
    }

    pub fn opposite(&self) -> Self {
        match self {
            Strand::Positive => Strand::Negative,
            Strand::Negative => Strand::Positive,
        }
    }
}

impl Display for Strand {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_char())
    }
}

#[derive(
    Debug, Copy, Clone, Eq, PartialEq, Hash, ValueEnum, PartialOrd, Ord,
)]
pub enum StrandRule {
    #[clap(name = "positive")]
    Positive,
    #[clap(name = "negative")]
    Negative,
    #[clap(name = "both")]
    Both,
}

impl StrandRule {
    pub fn overlaps(&self, other: &Self) -> bool {
        (self == &StrandRule::Both || other == &StrandRule::Both) || {
            match self {
                StrandRule::Positive => other == &StrandRule::Positive,
                StrandRule::Negative => other == &StrandRule::Negative,
                _ => unreachable!(),
            }
        }
    }

    pub fn covers(&self, strand: Strand) -> bool {
        match &self {
            StrandRule::Positive => strand == Strand::Positive,
            StrandRule::Negative => strand == Strand::Negative,
            StrandRule::Both => true,
        }
    }

    pub fn same_as(&self, strand: Strand) -> bool {
        match &self {
            StrandRule::Positive => strand == Strand::Positive,
            StrandRule::Negative => strand == Strand::Negative,
            StrandRule::Both => false,
        }
    }

    pub fn absorb(self, strand: Strand) -> Self {
        if self.same_as(strand) {
            self
        } else {
            // self is either both or they are opposite strands
            // so that means to "absorb" the rule is now both
            StrandRule::Both
        }
    }

    pub fn combine(self, other: Self) -> Self {
        if self == other {
            self
        } else {
            Self::Both
        }
    }
}

impl From<Strand> for StrandRule {
    fn from(value: Strand) -> Self {
        match value {
            Strand::Positive => Self::Positive,
            Strand::Negative => Self::Negative,
        }
    }
}

impl TryFrom<char> for StrandRule {
    type Error = anyhow::Error;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value {
            '+' => Ok(Self::Positive),
            '-' => Ok(Self::Negative),
            '.' => Ok(Self::Both),
            _ => bail!("illegal strand rule {value}"),
        }
    }
}

impl From<StrandRule> for char {
    fn from(value: StrandRule) -> Self {
        match value {
            StrandRule::Positive => '+',
            StrandRule::Negative => '-',
            StrandRule::Both => '.',
        }
    }
}

impl Display for StrandRule {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let lab = match self {
            Self::Positive => '+',
            Self::Negative => '-',
            Self::Both => '.',
        };

        write!(f, "{}", lab)
    }
}

#[inline]
pub fn record_is_primary(record: &bam::Record) -> bool {
    !record_is_not_primary(record)
}

#[inline]
pub fn record_is_not_primary(record: &bam::Record) -> bool {
    record.is_supplementary() || record.is_secondary() || record.is_duplicate()
}

pub(crate) fn get_targets(
    header: &HeaderView,
    region: Option<&Region>,
) -> Vec<ReferenceRecord> {
    (0..header.target_count())
        .filter_map(|tid| {
            let chrom_name = String::from_utf8(header.tid2name(tid).to_vec())
                .unwrap_or("???".to_owned());
            if let Some(region) = &region {
                if chrom_name == region.name {
                    Some(ReferenceRecord::new(
                        tid,
                        region.start,
                        region.length(),
                        chrom_name,
                    ))
                } else {
                    None
                }
            } else {
                match header.target_len(tid) {
                    Some(size) => Some(ReferenceRecord::new(
                        tid,
                        0,
                        size as u32,
                        chrom_name,
                    )),
                    None => {
                        debug!(
                            "no size information for {chrom_name} (tid: {tid})"
                        );
                        None
                    }
                }
            }
        })
        .collect::<Vec<ReferenceRecord>>()
}

#[derive(Debug, new)]
pub struct ReferenceRecord {
    // todo make this usize and unify all of the "Genome types"
    pub tid: u32,
    pub start: u32,
    pub(crate) length: u32,
    pub name: String,
}

impl ReferenceRecord {
    pub fn end(&self) -> u32 {
        self.start + self.length
    }
}

#[derive(new, Debug, Eq, PartialEq)]
pub struct Region {
    pub name: String,
    pub start: u32,
    pub end: u32,
}

impl Region {
    pub fn length(&self) -> u32 {
        self.end - self.start
    }

    fn parse_start_stop(raw: &str) -> Option<(u32, u32)> {
        fn parse_coordinates(input: &str) -> IResult<&str, (u32, u32)> {
            let (rest, start) = nom::character::complete::u32(input)?;
            let (rest, _) = tag("-")(rest)?;
            let (rest, stop) = nom::character::complete::u32(rest)?;
            Ok((rest, (start as u32, stop as u32)))
        }

        let is_coordinates = raw
            .chars()
            .filter(|c| (*c != '-') && (*c != ','))
            .all(|c| c.is_numeric());
        let has_sep = raw.contains('-');
        if is_coordinates && has_sep {
            let no_commas = raw.replace(",", "");
            parse_coordinates(&no_commas).map(|(_, (s, t))| (s, t)).ok()
        } else {
            None
        }
    }

    fn get_region_subsection(
        contig: &str,
        start: u32,
        stop: u32,
        header: &HeaderView,
    ) -> MkResult<Self> {
        let target_id = (0..header.target_count()).find_map(|tid| {
            String::from_utf8(header.tid2name(tid).to_vec()).ok().and_then(
                |sq_record| if &sq_record == contig { Some(tid) } else { None },
            )
        });

        let target_length = target_id.and_then(|tid| header.target_len(tid));
        if let Some(len) = target_length {
            let end = std::cmp::min(stop as u64, len) as u32;
            Ok(Self { name: contig.to_owned(), start, end })
        } else {
            Err(MkError::ContigMissing(contig.to_string()))
        }
    }

    pub fn parse_str(raw: &str, header: &HeaderView) -> MkResult<Self> {
        let final_colon_pos = raw
            .rfind(":")
            // add one to remove the ":"
            .map(|x| std::cmp::min(x.saturating_add(1), raw.len()));
        if let Some(final_col_pos) = final_colon_pos {
            let start_stop = raw.substring(final_col_pos, raw.len());
            let contig = raw.substring(0, final_col_pos.saturating_sub(1));
            if let Some((start, stop)) = Self::parse_start_stop(start_stop) {
                Self::get_region_subsection(contig, start, stop, header)
            } else {
                Self::get_region_subsection(raw, 0, u32::MAX, header)
            }
        } else {
            Self::get_region_subsection(raw, 0, u32::MAX, header)
        }
    }

    pub fn get_fetch_definition(
        &self,
        header: &HeaderView,
    ) -> AnyhowResult<bam::FetchDefinition> {
        let tid = (0..header.target_count())
            .find_map(|tid| {
                String::from_utf8(header.tid2name(tid).to_vec()).ok().and_then(
                    |chrom| {
                        if &chrom == &self.name {
                            Some(tid)
                        } else {
                            None
                        }
                    },
                )
            })
            .ok_or(anyhow!(
                "failed to find target ID for chrom {}",
                self.name.as_str()
            ))?;
        let tid = tid as i32;
        Ok(bam::FetchDefinition::Region(
            tid,
            self.start as i64,
            self.end as i64,
        ))
    }

    pub(crate) fn to_string(&self) -> String {
        format!("{}:{}-{}", self.name, self.start, self.end)
    }
}

// shouldn't need this once it's fixed in rust-htslib or the repo moves to
// noodles..
fn header_to_hashmap(
    header: &bam::Header,
) -> anyhow::Result<HashMap<String, Vec<LinearMap<String, String>>>> {
    let mut header_map = HashMap::default();
    let record_type_regex = Regex::new(r"@([A-Z][A-Z])").unwrap();
    let tag_regex = Regex::new(r"([A-Za-z][A-Za-z0-9]):([ -~]*)").unwrap();

    if let Ok(header_string) = String::from_utf8(header.to_bytes()) {
        for line in header_string.split('\n').filter(|x| !x.is_empty()) {
            let parts: Vec<_> =
                line.split('\t').filter(|x| !x.is_empty()).collect();
            if parts.is_empty() {
                continue;
            }
            let record_type = record_type_regex
                .captures(parts[0])
                .and_then(|captures| captures.get(1))
                .map(|m| m.as_str().to_owned());

            if let Some(record_type) = record_type {
                if record_type == "CO" {
                    continue;
                }
                let mut field = LinearMap::default();
                for part in parts.iter().skip(1) {
                    if let Some(cap) = tag_regex.captures(part) {
                        let tag = cap.get(1).unwrap().as_str().to_owned();
                        let value = cap.get(2).unwrap().as_str().to_owned();
                        field.insert(tag, value);
                    } else {
                        debug!("encounted illegal record line {line}");
                    }
                }
                header_map
                    .entry(record_type)
                    .or_insert_with(Vec::new)
                    .push(field);
            } else {
                debug!("encountered illegal record type in line {line}");
            }
        }
        Ok(header_map)
    } else {
        bail!("failed to parse header string")
    }
}

pub fn add_modkit_pg_records(header: &mut bam::Header) {
    let header_map = match header_to_hashmap(&header) {
        Ok(hm) => hm,
        Err(_) => {
            error!(
                "failed to parse input BAM header, not adding PG header \
                 record for modkit"
            );
            return;
        }
    };
    let (id, pp) = if let Some(pg_tags) = header_map.get("PG") {
        let modkit_invocations = pg_tags.iter().filter_map(|tags| {
            tags.get("ID").and_then(|v| {
                if v.contains("modkit") {
                    let last_run = v.split('.').nth(1).unwrap_or("0");
                    last_run.parse::<usize>().ok()
                } else {
                    None
                }
            })
        });
        if let Some(latest_run_number) = modkit_invocations.max() {
            let pp = if latest_run_number > 0 {
                Some(format!("modkit.{}", latest_run_number))
            } else {
                Some(format!("modkit"))
            };
            (format!("modkit.{}", latest_run_number + 1), pp)
        } else {
            (format!("modkit"), None)
        }
    } else {
        (format!("modkit"), None)
    };

    let command_line = std::env::args().collect::<Vec<String>>();
    let command_line = command_line.join(" ");
    let version = env!("CARGO_PKG_VERSION");
    let mut modkit_header_record = HeaderRecord::new("PG".as_bytes());
    modkit_header_record.push_tag("ID".as_bytes(), &id);
    modkit_header_record.push_tag("PN".as_bytes(), &"modkit".to_owned());
    modkit_header_record.push_tag("VN".as_bytes(), &version.to_owned());
    if let Some(pp) = pp {
        modkit_header_record.push_tag("PP".as_bytes(), &pp);
    }
    modkit_header_record.push_tag("CL".as_bytes(), &command_line);

    header.push_record(&modkit_header_record);
}

#[derive(new, Debug, Eq, PartialEq, Hash, Ord, PartialOrd, Copy, Clone)]
pub struct SamTag {
    inner: [u8; 2],
}

#[cfg(test)]
impl SamTag {
    pub(crate) fn parse(chars: [char; 2]) -> Self {
        Self { inner: [chars[0] as u8, chars[1] as u8] }
    }
}

pub(crate) fn get_stringable_aux(
    record: &bam::Record,
    sam_tag: &SamTag,
) -> Option<String> {
    record.aux(&sam_tag.inner).ok().and_then(|aux| match aux {
        Aux::String(s) => Some(s.to_string()),
        Aux::Char(c) => Some(format!("{}", c)),
        Aux::Double(f) => Some(format!("{}", f)),
        Aux::Float(f) => Some(format!("{}", f)),
        Aux::HexByteArray(a) => Some(a.to_string()),
        Aux::I8(i) => Some(format!("{}", i)),
        Aux::I16(i) => Some(format!("{}", i)),
        Aux::I32(i) => Some(format!("{}", i)),
        Aux::U8(u) => Some(format!("{}", u)),
        Aux::U16(u) => Some(format!("{}", u)),
        Aux::U32(u) => Some(format!("{}", u)),
        _ => None,
    })
}

pub(crate) fn parse_partition_tags(
    raw_tags: &[String],
) -> anyhow::Result<Vec<SamTag>> {
    let mut tags_seen = HashSet::with_capacity(raw_tags.len());
    let mut tags = Vec::with_capacity(raw_tags.len());
    for raw_tag in raw_tags {
        if raw_tag.len() != 2 {
            bail!("illegal tag {raw_tag} should be length 2")
        }
        let raw_tag_parts = raw_tag.chars().collect::<Vec<char>>();
        assert_eq!(raw_tag_parts.len(), 2);
        let inner = [raw_tag_parts[0] as u8, raw_tag_parts[1] as u8];
        let tag = SamTag::new(inner);

        let inserted = tags_seen.insert(tag);
        if inserted {
            tags.push(tag);
        } else {
            bail!("cannot repeat partition-tags, got {raw_tag} twice")
        }
    }

    Ok(tags)
}

#[inline]
pub fn get_reference_mod_strand(
    read_mod_strand: Strand,
    alignment_strand: Strand,
) -> Strand {
    match (read_mod_strand, alignment_strand) {
        (Strand::Positive, Strand::Positive) => Strand::Positive,
        (Strand::Positive, Strand::Negative) => Strand::Negative,
        (Strand::Negative, Strand::Positive) => Strand::Negative,
        (Strand::Negative, Strand::Negative) => Strand::Positive,
    }
}

#[inline]
pub(crate) fn reader_is_bam(reader: &bam::IndexedReader) -> bool {
    unsafe {
        (*reader.htsfile()).format.format
            == rust_htslib::htslib::htsExactFormat_bam
    }
}

pub(crate) const KMER_SIZE: usize = 50;

#[derive(Copy, Clone)]
pub(crate) struct Kmer {
    inner: [u8; KMER_SIZE],
    pub(crate) size: usize,
}

impl Kmer {
    pub(crate) fn from_seq(seq: &[u8], pos: usize, kmer_size: usize) -> Kmer {
        Kmer::new(seq, pos, kmer_size)
    }

    // kinda risky, size needs to be < 12
    pub(crate) fn new(seq: &[u8], position: usize, size: usize) -> Self {
        if size > KMER_SIZE {
            debug!("kmers greater that size {KMER_SIZE} will be corrupted");
        }
        let get_back_base_safe = |i| -> Option<u8> {
            position.checked_sub(i).and_then(|idx| seq.get(idx).map(|b| *b))
        };
        let before = if size % 2 == 0 { size / 2 - 1 } else { size / 2 };

        let after = size / 2;
        let mut buffer = [Some(45u8); KMER_SIZE];
        let mut i = 0;
        let mut assign = |b: Option<u8>| {
            buffer[i] = b;
            i = std::cmp::min(i + 1, KMER_SIZE - 1);
        };

        for offset in (1..=before).rev() {
            let b = get_back_base_safe(offset);
            assign(b);
        }
        assign(seq.get(position).map(|b| *b));
        for offset in 1..=after {
            assign(seq.get(position + offset).map(|b| *b))
        }
        let inner = buffer.map(|b| b.unwrap_or(45));
        Self { inner, size }
    }

    pub(crate) fn reverse_complement(self) -> Self {
        let mut inner = [45u8; KMER_SIZE];
        for (i, p) in (0..self.size).rev().enumerate() {
            let mut b = self.inner[p];
            if b != 45 {
                b = complement(b)
            }
            inner[i] = b
        }
        Self { inner, size: self.size }
    }

    #[cfg(test)]
    pub(crate) fn get_nt(&self, pos: usize) -> Option<u8> {
        if pos > self.size || pos > KMER_SIZE {
            None
        } else {
            Some(self.inner[pos])
        }
    }
}

impl Debug for Kmer {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let s = self.inner.iter().take(self.size).map(|b| *b as char).join("");
        write!(f, "{s}")
    }
}

impl Display for Kmer {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[inline]
pub fn within_alignment(
    query_position: usize,
    num_soft_clipped_start: usize,
    num_soft_clipped_end: usize,
    read_length: usize,
) -> bool {
    read_length
        .checked_sub(num_soft_clipped_end)
        .map(|x| query_position >= num_soft_clipped_start && query_position < x)
        .unwrap_or_else(|| {
            debug!(
                "read_length ({read_length}) is less than \
                 num_soft_clipped_end ({num_soft_clipped_end})"
            );
            false
        })
}

pub fn format_int_with_commas(val: isize) -> String {
    let mut num = val
        .abs()
        .to_string()
        .as_bytes()
        .rchunks(3)
        .rev()
        .map(str::from_utf8)
        .collect::<Result<Vec<&str>, _>>()
        .unwrap()
        .join(",");
    if val < 0 {
        num = format!("-{num}")
    }
    num
}

#[derive(new, Debug, Clone, PartialEq, Eq, Hash)]
pub struct GenomeRegion {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: StrandRule,
    pub name: Option<String>,
}

impl GenomeRegion {
    pub fn midpoint(&self) -> u64 {
        (self.start + self.end) / 2
    }

    #[inline]
    fn parse_bed_line(l: &str) -> IResult<&str, Self> {
        let (rest, chrom) = consume_string(l)?;
        let (rest, start) = consume_digit(rest)?;
        let (rest, stop) = consume_digit(rest)?;

        let (rest, name) = many0(one_of(" \t\r\n"))(rest)
            .and_then(|(rest, _)| consume_string_spaces(rest))
            .map(|(rest, name)| (rest, Some(name)))
            .unwrap_or_else(|_| (rest, None));
        Ok((
            rest,
            Self { chrom, start, end: stop, name, strand: StrandRule::Both },
        ))
    }

    pub(super) fn parse_unstranded_bed_line(
        line: &str,
    ) -> anyhow::Result<Self> {
        Self::parse_bed_line(line)
            .map_err(|e| {
                anyhow!(
                    "failed to parse un-stranded (bed3/4/5) line: {line}, {e}"
                )
            })
            .map(|(_, this)| this)
    }

    pub(super) fn parse_stranded_bed_line(line: &str) -> anyhow::Result<Self> {
        fn inner(l: &str) -> IResult<&str, GenomeRegion> {
            let mut parse_strand =
                map_res(consume_char, |x| StrandRule::try_from(x));
            let (rest, mut this) = GenomeRegion::parse_bed_line(l)?;
            let (rest, _score) = consume_float(rest)
                .or_else(|_| consume_dot(rest).map(|(r, _)| (r, 0f32)))?;

            let (rest, strand) = parse_strand(rest)?;
            this.strand = strand;
            Ok((rest, this))
        }
        inner(line)
            .map_err(|e| {
                anyhow!("failed to parse stranded (bed5+) line: {line}, {e}")
            })
            .map(|(_, this)| this)
    }
}

// todo could make this a trait and have some of the other structs implement it,
// like BedMethylLine
#[derive(new, Debug)]
pub(crate) struct ModPositionInfo<T> {
    pub n_valid: T,
    pub n_mod: T,
}

impl<T: num_traits::Num + num_traits::cast::AsPrimitive<f32> + Debug>
    ModPositionInfo<T>
{
    pub(crate) fn frac_modified(&self) -> f32 {
        if self.n_valid == T::zero() {
            0f32
        } else {
            let n_mod: f32 = self.n_mod.as_();
            let n_valid: f32 = self.n_valid.as_();
            n_mod / n_valid
        }
    }

    pub(crate) fn percent_modified(&self) -> f32 {
        self.frac_modified() * 100f32
    }
}

impl<T> Moniod for ModPositionInfo<T>
where
    T: num_traits::Num
        + num_traits::cast::AsPrimitive<f32>
        + num_traits::cast::AsPrimitive<usize>,
{
    fn zero() -> Self {
        Self { n_valid: T::zero(), n_mod: T::zero() }
    }

    fn op(self, other: Self) -> Self {
        Self {
            n_mod: self.n_mod + other.n_mod,
            n_valid: self.n_valid + other.n_valid,
        }
    }

    fn op_mut(&mut self, other: Self) {
        let n_mod = self.n_mod + other.n_mod;
        let n_valid = self.n_valid + other.n_valid;
        self.n_mod = n_mod;
        self.n_valid = n_valid;
    }

    fn len(&self) -> usize {
        let l: usize = self.n_valid.as_();
        l
    }
}

/// Read a "contig sizes" tab-separated file.
pub(crate) fn read_sequence_lengths_file(
    p: &PathBuf,
) -> anyhow::Result<IndexMap<String, u64>> {
    fn parse_line(line: &str) -> IResult<&str, (String, u64)> {
        let (rest, chrom) = consume_string(line)?;
        let (rest, length) = consume_digit(rest)?;
        Ok((rest, (chrom, length)))
    }

    BufReader::new(std::fs::File::open(p)?)
        .lines()
        .map(|l| {
            l.map_err(|e| anyhow!("failed to read from sizes, {e}")).and_then(
                |l| {
                    parse_line(&l)
                        .map(|(_, this)| this)
                        .map_err(|e| anyhow!("failed to parse sizes {l}, {e}"))
                },
            )
        })
        .collect::<anyhow::Result<IndexMap<_, _>>>()
}

pub(crate) fn format_errors_table(
    error_counts: &FxHashMap<String, usize>,
) -> prettytable::Table {
    let mut tab = get_human_readable_table();
    tab.set_titles(row!["error", "count"]);
    error_counts.iter().sorted_by(|(_, a), (_, b)| a.cmp(b)).for_each(
        |(er, c)| {
            tab.add_row(row![er, c]);
        },
    );
    tab
}

pub(crate) fn get_human_readable_table() -> prettytable::Table {
    let mut tab = prettytable::Table::new();
    tab.set_format(*prettytable::format::consts::FORMAT_NO_LINESEP_WITH_TITLE);
    tab
}

#[cfg(test)]
mod utils_tests {
    use anyhow::Context;
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use similar_asserts::assert_eq;

    use crate::errs::MkError;
    use crate::util::{
        get_query_name_string, get_stringable_aux, parse_partition_tags,
        GenomeRegion, Region, SamTag, StrandRule,
    };

    use super::Kmer;

    #[test]
    fn test_util_get_stringable_tag() {
        let bam_fp = "tests/resources/bc_anchored_10_reads.sorted.bam";
        let mut reader = bam::Reader::from_path(bam_fp).unwrap();
        let mut checked = false;
        for record in reader.records().filter_map(|r| r.ok()) {
            let read_id = get_query_name_string(&record).unwrap();
            if read_id == "068ce426-129e-4870-bd34-16cd78edaa43".to_string() {
                let tag = SamTag {
                    inner: [102, 98], // 'fb' not a tag that's there
                };
                assert!(get_stringable_aux(&record, &tag).is_none());
                let expected_rg =
                    "5598049b1b3264566b162bf035344e7ec610d608_dna_r10.4.1_e8.\
                     2_400bps_hac@v3.5.2"
                        .to_string();
                let tag = SamTag { inner: [82, 71] };
                assert_eq!(
                    get_stringable_aux(&record, &tag),
                    Some(expected_rg)
                );
                let tag = SamTag { inner: [114, 110] };
                assert_eq!(
                    get_stringable_aux(&record, &tag),
                    Some("6335".to_string())
                );
                checked = true
            }
        }
        assert!(checked)
    }

    #[test]
    fn test_util_parse_partition_tags() {
        let raw_tags = ["HP".to_string(), "RG".to_string()];
        let parsed = parse_partition_tags(&raw_tags)
            .context("should have parsed raw tags")
            .unwrap();
        let expected =
            vec![SamTag::parse(['H', 'P']), SamTag::parse(['R', 'G'])];
        assert_eq!(parsed, expected);
    }

    #[test]
    fn test_strand_rule_semantics() {
        let pos = StrandRule::Positive;
        let neg = StrandRule::Negative;
        let both = StrandRule::Both;
        assert!(both.overlaps(&pos));
        assert!(both.overlaps(&neg));
        assert!(pos.overlaps(&both));
        assert!(neg.overlaps(&both));
        assert!(pos.overlaps(&StrandRule::Positive));
        assert!(!pos.overlaps(&StrandRule::Negative));
        assert!(!neg.overlaps(&StrandRule::Positive));
        assert!(neg.overlaps(&StrandRule::Negative));
    }

    #[test]
    fn test_genome_region_parse_bedlines() {
        let line = "chr1\t938169\t938373\tmerged_peak1\t.\t.\n";
        let gr = GenomeRegion::parse_stranded_bed_line(line).unwrap();
        let expected = GenomeRegion {
            chrom: "chr1".to_string(),
            start: 938169,
            end: 938373,
            name: Some("merged_peak1".to_string()),
            strand: StrandRule::Both,
        };
        assert_eq!(gr, expected);
        let line = "chr1\t938169\t938373\tmerged_peak1\t.\t+\n";
        let gr = GenomeRegion::parse_stranded_bed_line(line).unwrap();
        let expected = GenomeRegion {
            chrom: "chr1".to_string(),
            start: 938169,
            end: 938373,
            name: Some("merged_peak1".to_string()),
            strand: StrandRule::Positive,
        };
        assert_eq!(gr, expected);
        let line = "chr1\t938169\t938373\tmerged_peak1\n";
        let gr = GenomeRegion::parse_unstranded_bed_line(line).unwrap();
        let expected = GenomeRegion {
            chrom: "chr1".to_string(),
            start: 938169,
            end: 938373,
            name: Some("merged_peak1".to_string()),
            strand: StrandRule::Both,
        };
        assert_eq!(gr, expected);
        let line = "chr1\t938169\t938373\tmerged_peak1\t1000\t+\n";
        let gr = GenomeRegion::parse_stranded_bed_line(line).unwrap();
        let expected = GenomeRegion {
            chrom: "chr1".to_string(),
            start: 938169,
            end: 938373,
            name: Some("merged_peak1".to_string()),
            strand: StrandRule::Positive,
        };
        assert_eq!(gr, expected);
        let line = "chr20\t9838623\t9839213\tCpG: 47\n";
        let gr = GenomeRegion::parse_unstranded_bed_line(line).unwrap();
        let expected = GenomeRegion {
            chrom: "chr20".to_string(),
            start: 9838623,
            end: 9839213,
            name: Some("CpG: 47".to_string()),
            strand: StrandRule::Both,
        };
        assert_eq!(gr, expected);
    }

    #[test]
    fn test_kmer_get_nt() {
        let seq = "GATTACA".as_bytes();
        let kmer = Kmer::from_seq(seq, 2, 5);
        assert_eq!(format!("{kmer}"), "GATTA".to_string());
        let nt = kmer.get_nt(0).unwrap();
        assert_eq!(nt, 'G' as u8);
        assert!(kmer.get_nt(6).is_none());
    }

    #[test]
    fn test_parse_coordinates() {
        let raw = "1-2,000";
        let (start, stop) = Region::parse_start_stop(raw).unwrap();
        assert_eq!(start, 1);
        assert_eq!(stop, 2000);
        let raw = "1,200-2,000";
        let (start, stop) = Region::parse_start_stop(raw).unwrap();
        assert_eq!(start, 1200);
        assert_eq!(stop, 2000);
        let raw = "000,1-2,000";
        let (start, stop) = Region::parse_start_stop(raw).unwrap();
        assert_eq!(start, 1);
        assert_eq!(stop, 2000);
        let raw = "1200";
        let x = Region::parse_start_stop(raw);
        assert!(x.is_none());
    }

    #[test]
    #[rustfmt::skip]
    fn test_parse_transcript_region() {
        let reader =
            bam::Reader::from_path("tests/resources/transcriptome_header.sam")
                .unwrap();
        let raw = "ENST00000616016.5|ENSG00000187634.13|OTTHUMG00000040719.11|OTTHUMT00000316521.3|SAMD11-209|SAMD11|3465|UTR5:1-509|CDS:510-3044|UTR3:3045-3465|:1-3,200";
        let region = Region::parse_str(raw, &reader.header()).unwrap();
        assert_eq!(
            &region.name,
            "ENST00000616016.5|ENSG00000187634.13|OTTHUMG00000040719.11|OTTHUMT00000316521.3|SAMD11-209|SAMD11|3465|UTR5:1-509|CDS:510-3044|UTR3:3045-3465|"
        );
        assert_eq!(region.start, 1);
        assert_eq!(region.end, 3200);
        let raw = "ENST00000616016.5|ENSG00000187634.13|OTTHUMG00000040719.11|OTTHUMT00000316521.3|SAMD11-209|SAMD11|3465|UTR5:1-509|CDS:510-3044|UTR3:3045-3465|:1-99,999";
        let region = Region::parse_str(raw, &reader.header()).unwrap();
        assert_eq!(
            &region.name,
            "ENST00000616016.5|ENSG00000187634.13|OTTHUMG00000040719.11|OTTHUMT00000316521.3|SAMD11-209|SAMD11|3465|UTR5:1-509|CDS:510-3044|UTR3:3045-3465|"
        );
        assert_eq!(region.start, 1);
        assert_eq!(region.end, 3465);
        let raw = "ENST00000616016.5|ENSG00000187634.13|OTTHUMG00000040719.11|OTTHUMT00000316521.3|SAMD11-209|SAMD11|3465|UTR5:1-509|CDS:510-3044|UTR3:3045-3465|";
        let region = Region::parse_str(raw, &reader.header()).unwrap();
        assert_eq!(
            &region.name,
            "ENST00000616016.5|ENSG00000187634.13|OTTHUMG00000040719.11|OTTHUMT00000316521.3|SAMD11-209|SAMD11|3465|UTR5:1-509|CDS:510-3044|UTR3:3045-3465|"
        );
        assert_eq!(region.start, 0);
        assert_eq!(region.end, 3465);
        let raw = "tag|anothertag|decoy:100-200|:2-3,200";
        let region = Region::parse_str(raw, &reader.header());
        match region.unwrap_err() {
            MkError::ContigMissing(s) => {
                assert_eq!(s, "tag|anothertag|decoy:100-200|".to_string())
            }
            e @ _ => assert!(false, "incorrect error {e}"),
        }
        let reader = bam::Reader::from_path("tests/resources/bc_anchored_10_reads.sorted.bam").unwrap();
        let raw = "oligo_1512_adapters:1-10";
        let region = Region::parse_str(raw, reader.header()).unwrap();
        assert_eq!(&region.name, "oligo_1512_adapters");
        assert_eq!(region.start, 1);
        assert_eq!(region.end, 10);
        let raw = "oligo_1512_adapters:-1-10";
        let region = Region::parse_str(raw, reader.header());
        match region.unwrap_err() {
            MkError::ContigMissing(s) => {
                assert_eq!(s, "oligo_1512_adapters:-1-10".to_string())
            }
            e @ _ => assert!(false, "incorrect error {e}"),
        }
    }
}
