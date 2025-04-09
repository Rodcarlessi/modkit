use anyhow::{anyhow, bail, Result as AnyhowResult};
use derive_new::new;
use mod_kit::mod_bam::{CollapseMethod, EdgeFilter};
use mod_kit::position_filter::StrandedPositionFilter;
use mod_kit::summarize::{summarize_modbam, ModSummary};
use mod_kit::threshold_mod_caller::MultipleThresholdModCaller;
use serde::Deserialize;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::process::Output;

pub fn run_modkit(args: &[&str]) -> AnyhowResult<Output> {
    let exe = Path::new(env!("CARGO_BIN_EXE_modkit"));
    assert!(exe.exists());

    let output = std::process::Command::new(exe)
        .args(args)
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .spawn()?
        .wait_with_output()?;
    if output.status.success() {
        Ok(output)
    } else {
        Err(anyhow!("failed to run {:?}", args.join(" ")))
    }
}

fn run_summary<'a>(
    bam_fp: &str,
    interval_size: u32,
    collapse_method: Option<&CollapseMethod>,
    edge_filter: Option<&EdgeFilter>,
) -> AnyhowResult<ModSummary<'a>> {
    let threads = 1usize;
    let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build()?;
    pool.install(|| {
        summarize_modbam(
            &Path::new(bam_fp).to_path_buf(),
            threads,
            interval_size,
            None,
            None,
            None,
            None,
            0.1, // doesn't matter
            Some(MultipleThresholdModCaller::new_passthrough()),
            None,
            collapse_method,
            edge_filter,
            None,
            false,
            true,
        )
    })
}

pub fn run_simple_summary(
    bam_fp: &str,
    interval_size: u32,
) -> AnyhowResult<ModSummary> {
    run_summary(bam_fp, interval_size, None, None)
}

pub fn run_summary_with_include_positions<'a>(
    bam_fp: &PathBuf,
    include_bed_fp: &PathBuf,
) -> anyhow::Result<ModSummary<'a>> {
    let threads = 1usize;
    let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build()?;
    let position_filter =
        StrandedPositionFilter::from_bam_and_bed(bam_fp, include_bed_fp, true)?;
    let caller = MultipleThresholdModCaller::new_passthrough();
    pool.install(|| {
        summarize_modbam(
            &Path::new(bam_fp).to_path_buf(),
            threads,
            32,
            None,
            None,
            None,
            None,
            0.1, // doesn't matter
            Some(caller),
            None,
            None,
            None,
            Some(&position_filter),
            true,
            true,
        )
    })
}

pub fn run_simple_summary_with_collapse_method<'a>(
    bam_fp: &str,
    interval_size: u32,
    collapse_method: &CollapseMethod,
) -> AnyhowResult<ModSummary<'a>> {
    run_summary(bam_fp, interval_size, Some(collapse_method), None)
}

pub fn run_simple_summary_with_edge_filter<'a>(
    bam_fp: &str,
    interval_size: u32,
    edge_filter: &EdgeFilter,
) -> AnyhowResult<ModSummary<'a>> {
    run_summary(bam_fp, interval_size, None, Some(edge_filter))
}

pub fn check_against_expected_text_file(output_fp: &str, expected_fp: &str) {
    assert_ne!(output_fp, expected_fp, "cannot check a file against itself");
    let test = {
        let mut fh = File::open(output_fp).unwrap();
        let mut buff = String::new();
        fh.read_to_string(&mut buff).unwrap();
        buff
    };
    let expected = {
        // this file was hand-checked for correctness or should be equivalent
        // due to orthogonal process
        let mut fh = File::open(expected_fp).unwrap();
        let mut buff = String::new();
        fh.read_to_string(&mut buff).unwrap();
        buff
    };

    similar_asserts::assert_eq!(
        test,
        expected,
        "{output_fp} is not the same as {expected_fp}"
    );
}

#[derive(Deserialize)]
pub struct ExtractFullRecord {
    read_id: String,
    forward_read_position: usize,
    ref_position: i64,
    mod_code: char,
    #[serde(rename(deserialize = "ref_strand"))]
    strand: char,
    read_length: usize,
    chrom: String,
}

#[derive(new, Eq, PartialEq, Debug)]
pub struct ModData {
    pub q_pos: usize,
    pub ref_pos: i64,
    pub mod_code: char,
    pub strand: char,
    pub read_length: usize,
    pub contig: String,
    // pub data: String,
}

impl PartialOrd for ModData {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ModData {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.q_pos.cmp(&other.q_pos) {
            Ordering::Equal => match self.mod_code.cmp(&other.mod_code) {
                Ordering::Equal => self.strand.cmp(&other.strand),
                ord => ord,
            },
            ord => ord,
        }
    }
}

pub fn parse_mod_profile(
    fp: &PathBuf,
) -> anyhow::Result<HashMap<String, Vec<ModData>>> {
    let mut agg = HashMap::new();
    let mut reader = csv::ReaderBuilder::new()
        .delimiter('\t' as u8)
        .has_headers(true)
        .from_path(fp)
        .unwrap();

    for record in reader.deserialize() {
        let record: ExtractFullRecord = record.unwrap();
        let read_id = record.read_id;
        agg.entry(read_id).or_insert_with(Vec::new).push(ModData::new(
            record.forward_read_position,
            record.ref_position,
            record.mod_code,
            record.strand,
            record.read_length,
            record.chrom,
        ))
    }
    // let mut reader =
    //     BufReader::new(File::open(fp)?).lines().map(|l| l.unwrap());
    // while let Some(line) = reader.next() {
    //     let parts = line.split_ascii_whitespace().collect::<Vec<&str>>();
    //     let read_id = parts[0].to_owned();
    //     let q_pos = parts[1].parse::<usize>().unwrap();
    //     let ref_pos = parts[2].parse::<i64>().unwrap();
    //     let mod_code = parts[13].parse::<char>().unwrap();
    //     let strand = parts[5].parse::<char>().unwrap();
    //     let read_length = parts[11].parse::<usize>().unwrap();
    //     let contig = parts[3].to_owned();
    //     agg.entry(read_id).or_insert(Vec::new()).push(ModData::new(
    //         q_pos,
    //         ref_pos,
    //         mod_code,
    //         strand,
    //         read_length,
    //         contig,
    //         // line,
    //     ));
    // }
    for (_, dat) in agg.iter_mut() {
        dat.sort()
    }

    Ok(agg)
}

pub fn check_legal_csv<const SEP: u8>(fp: &PathBuf) -> anyhow::Result<()> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(SEP)
        .from_reader(File::open(fp).expect("should open file"));
    let mut i = 1;
    for record in reader.records() {
        match record {
            Ok(_) => {}
            Err(e) => {
                bail!("failed to parse line at {i}, {e}")
            }
        }
        i += 1;
    }

    Ok(())
}
