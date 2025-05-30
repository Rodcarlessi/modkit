use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{anyhow, bail, Context};
use clap::{Args, Subcommand};
use indicatif::{MultiProgress, ProgressIterator};
use itertools::Itertools;
use log::{debug, error, info};
use prettytable::row;
use rustc_hash::FxHashMap;

use crate::dmr::bedmethyl::BedMethylLine;
use crate::dmr::pairwise::run_pairwise_dmr;
use crate::dmr::single_site::SingleSiteDmrAnalysis;
use crate::dmr::tabix::MultiSampleIndex;
use crate::dmr::util::{parse_roi_bed, HandleMissing, RoiIter};
use crate::errs::MkResult;
use crate::genome_positions::GenomePositions;
use crate::logging::init_logging;
use crate::mod_base_code::{DnaBase, ModCodeRepr, MOD_CODE_TO_DNA_BASE};
use crate::monoid::Moniod;
use crate::tabix::{BedMethylTbxIndex, HtsTabixHandler};
use crate::util::{
    create_out_directory, format_errors_table, get_master_progress_bar,
    get_subroutine_progress_bar, get_ticker,
};

#[derive(Subcommand)]
pub enum BedMethylDmr {
    /// Compare regions in a pair of samples (for example, tumor and normal or
    /// control and experiment). A sample is input as a bgzip pileup bedMethyl
    /// (produced by pileup, for example) that has an associated tabix index.
    /// Output is a BED file with the score column indicating the magnitude of
    /// the difference in methylation between the two samples. See the online
    /// documentation for additional details.
    Pair(PairwiseDmr),
    /// Compare regions between all pairs of samples (for example a trio sample
    /// set or haplotyped trio sample set). As with `pair` all inputs must be
    /// bgzip compressed bedMethyl files with associated tabix indices.
    /// Each sample must be assigned a name. Output is a directory of BED
    /// files with the score column indicating the magnitude of the
    /// difference in methylation between the two samples indicated in the
    /// file name. See the online documentation for additional details.
    Multi(MultiSampleDmr),
}

impl BedMethylDmr {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            Self::Pair(x) => x.run(),
            Self::Multi(x) => x.run(),
        }
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct PairwiseDmr {
    /// Bgzipped bedMethyl file for the first (usually control) sample. There
    /// should be a tabix index with the same name and .tbi next to this
    /// file or the --index-a option must be provided.
    #[clap(help_heading = "Sample Options")]
    #[arg(short = 'a')]
    control_bed_methyl: Vec<PathBuf>,
    /// Bgzipped bedMethyl file for the second (usually experimental) sample.
    /// There should be a tabix index with the same name and .tbi next to
    /// this file or the --index-b option must be provided.
    #[clap(help_heading = "Sample Options")]
    #[arg(short = 'b')]
    exp_bed_methyl: Vec<PathBuf>,
    /// Path to file to direct output, optional, no argument will direct output
    /// to stdout.
    #[clap(help_heading = "Output Options")]
    #[arg(short = 'o', long)]
    out_path: Option<String>,
    /// Include header in output
    #[clap(help_heading = "Output Options")]
    #[arg(long, alias = "with-header", default_value_t = false)]
    header: bool,
    /// BED file of regions over which to compare methylation levels. Should be
    /// tab-separated (spaces allowed in the "name" column). Requires
    /// chrom, chromStart and chromEnd. The Name column is optional. Strand
    /// is currently ignored. When omitted, methylation levels are compared at
    /// each site.
    #[arg(long, short = 'r', alias = "regions")]
    regions_bed: Option<PathBuf>,
    /// Path to reference fasta for used in the pileup/alignment.
    #[arg(long = "ref")]
    reference_fasta: PathBuf,
    /// Run segmentation, output segmented differentially methylated regions to
    /// this file.
    #[clap(help_heading = "Segmentation Options")]
    #[arg(long = "segment", conflicts_with = "regions_bed")]
    segmentation_fp: Option<PathBuf>,

    /// Maximum number of base pairs between modified bases for them to be
    /// segmented together.
    #[clap(help_heading = "Segmentation Options")]
    #[arg(long, requires = "segmentation_fp", default_value_t = 5000)]
    max_gap_size: u64,
    /// Prior probability of a differentially methylated position
    #[clap(help_heading = "Segmentation Options")]
    #[arg(
        long,
        requires = "segmentation_fp",
        default_value_t = 0.1,
        hide_short_help = true
    )]
    dmr_prior: f64,
    /// Maximum probability of continuing a differentially methylated block,
    /// decay will be dynamic based on proximity to the next position.
    #[clap(help_heading = "Segmentation Options")]
    #[arg(
        long,
        requires = "segmentation_fp",
        default_value_t = 0.9,
        hide_short_help = true
    )]
    diff_stay: f64,
    /// Significance factor, effective p-value necessary to favor the
    /// "Different" state.
    #[clap(help_heading = "Segmentation Options")]
    #[arg(
        long,
        requires = "segmentation_fp",
        default_value_t = 0.01,
        hide_short_help = true
    )]
    significance_factor: f64,
    /// Use logarithmic decay for "Different" stay probability
    #[clap(help_heading = "Segmentation Options")]
    #[arg(
        long,
        requires = "segmentation_fp",
        default_value_t = false,
        hide_short_help = true
    )]
    log_transition_decay: bool,
    /// After this many base pairs, the transition probability will become the
    /// prior probability of encountering a differentially modified
    /// position.
    #[clap(help_heading = "Segmentation Options")]
    #[arg(
        long,
        requires = "segmentation_fp",
        default_value_t = 500,
        hide_short_help = true
    )]
    decay_distance: u32,
    /// Preset HMM segmentation parameters for higher propensity to switch from
    /// "Same" to "Different" state. Results will be shorter segments, but
    /// potentially higher sensitivity.
    #[clap(help_heading = "Segmentation Options")]
    #[arg(
        long,
        requires = "segmentation_fp",
        conflicts_with_all=["log_transition_decay", "significance_factor", "diff_stay", "dmr_prior"],
        default_value_t=false
    )]
    fine_grained: bool,
    /// Bases to use to calculate DMR, may be multiple. For example, to
    /// calculate differentially methylated regions using only cytosine
    /// modifications use --base C.
    #[clap(help_heading = "Sample Options")]
    #[arg(short, long="base", alias = "modified-bases", action=clap::ArgAction::Append)]
    modified_bases: Vec<char>,
    /// Extra assignments of modification codes to their respective primary
    /// bases. In general, modkit dmr will use the SAM specification to
    /// know which modification codes are appropriate to use for a given
    /// primary base. For example "h" is the code for 5hmC, so is appropriate
    /// for cytosine bases, but not adenine bases. However, if your
    /// bedMethyl file contains custom codes or codes that are not part of
    /// the specification, you can specify which primary base they
    /// belong to here with --assign-code x:C meaning associate modification
    /// code "x" with cytosine (C) primary sequence bases. If a code is
    /// encountered that is not part of the specification, the bedMethyl
    /// record will not be used, this will be logged.
    #[clap(help_heading = "Sample Options")]
    #[arg(long="assign-code", action=clap::ArgAction::Append)]
    mod_code_assignments: Option<Vec<String>>,

    /// Log out which sequences are in common between the samples and the
    /// reference FASTA, useful for debugging
    #[clap(help_heading = "Logging Options")]
    #[arg(
        long = "careful",
        requires = "log_filepath",
        hide_short_help = true,
        default_value_t = false
    )]
    careful: bool,
    /// File to write logs to, it's recommended to use this option.
    #[clap(help_heading = "Logging Options")]
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Number of threads to use.
    #[clap(help_heading = "Compute Options")]
    #[arg(short = 't', long, default_value_t = 4)]
    threads: usize,
    /// Number of threads to use when for decompression.
    #[clap(help_heading = "Compute Options")]
    #[arg(long, default_value_t = 4)]
    io_threads: usize,
    /// Control the  batch size. The batch size is the number of regions to
    /// load at a time. Each region will be processed concurrently. Loading
    /// more regions at a time will decrease IO to load data, but will use
    /// more memory. Default will be 50% more than the number of
    /// threads assigned.
    #[clap(help_heading = "Compute Options")]
    #[arg(long, alias = "batch")]
    batch_size: Option<usize>,
    /// Respect soft masking in the reference FASTA.
    #[clap(help_heading = "Sample Options")]
    #[arg(long, short = 'k', default_value_t = false)]
    mask: bool,
    /// Don't show progress bars
    #[clap(help_heading = "Logging Options")]
    #[arg(long, default_value_t = false)]
    suppress_progress: bool,
    /// Force overwrite of output file, if it already exists.
    #[clap(help_heading = "Compute Options")]
    #[arg(short = 'f', long, default_value_t = false)]
    force: bool,
    /// How to handle regions found in the `--regions` BED file.
    /// quiet => ignore regions that are not found in the tabix header
    /// warn => log (debug) regions that are missing
    /// fatal => log (error) and exit the program when a region is missing.
    #[clap(help_heading = "Logging Options")]
    #[arg(long="missing", requires = "regions_bed", default_value_t=HandleMissing::quiet)]
    handle_missing: HandleMissing,
    /// Minimum valid coverage required to use an entry from a bedMethyl. See
    /// the help for pileup for the specification and description of valid
    /// coverage.
    #[clap(help_heading = "Sample Options")]
    #[arg(long, alias = "min-coverage", default_value_t = 0)]
    min_valid_coverage: u64,
    /// Prior distribution for estimating MAP-based p-value. Should be two
    /// arguments for alpha and beta (e.g. 1.0 1.0). See
    /// `dmr_scoring_details.md` for additional details on how the metric
    /// is calculated.

    #[clap(help_heading = "Single-site Options")]
    #[arg(
        long,
        num_args = 2,
        conflicts_with = "regions_bed",
        hide_short_help = true
    )]
    prior: Option<Vec<f64>>,
    /// Consider only effect sizes greater than this when calculating the
    /// MAP-based p-value.
    #[clap(help_heading = "Single-site Options")]
    #[arg(
        long,
        default_value_t = 0.05,
        conflicts_with = "regions_bed",
        hide_short_help = true
    )]
    delta: f64,
    /// Sample this many reads when estimating the max coverage thresholds.
    #[clap(help_heading = "Single-site Options")]
    #[arg(
        long,
        short='N',
        default_value_t = 10_042,
        conflicts_with_all = ["max_coverages", "regions_bed"],
    )]
    n_sample_records: usize,
    /// Max coverages to enforce when calculating estimated MAP-based p-value.
    #[clap(help_heading = "Single-site Options")]
    #[arg(long, num_args = 2, conflicts_with = "regions_bed")]
    max_coverages: Option<Vec<usize>>,
    /// When using replicates, cap coverage to be equal to the maximum coverage
    /// for a single sample. For example, if there are 3 replicates with
    /// max_coverage of 30, the total coverage would normally be 90. Using
    /// --cap-coverages will down sample the data to 30X.
    #[clap(help_heading = "Single-site Options")]
    #[arg(
        long,
        conflicts_with = "regions_bed",
        default_value_t = false,
        hide_short_help = true
    )]
    cap_coverages: bool,
    /// Interval chunk size in base pairs to process concurrently. Smaller
    /// interval chunk sizes will use less memory but incur more overhead.
    #[clap(help_heading = "Compute Options")]
    #[arg(
        short = 'i',
        long,
        default_value_t = 100_000,
        hide_short_help = true
    )]
    interval_size: u64,
}

impl PairwiseDmr {
    fn check_modified_bases(
        &self,
    ) -> anyhow::Result<FxHashMap<ModCodeRepr, DnaBase>> {
        Self::validate_modified_bases(
            &self.modified_bases,
            self.mod_code_assignments.as_ref(),
        )
    }

    fn is_single_site(&self) -> bool {
        self.regions_bed.is_none()
    }

    fn parse_raw_assignments(
        raw_mod_code_assignments: Option<&Vec<String>>,
    ) -> anyhow::Result<FxHashMap<ModCodeRepr, DnaBase>> {
        if let Some(raw_assignments) = raw_mod_code_assignments {
            let user_assignments = raw_assignments.iter().try_fold(
                FxHashMap::default(),
                |mut acc, next| {
                    if next.contains(':') {
                        let parts = next.split(':').collect::<Vec<&str>>();
                        if parts.len() != 2 {
                            bail!(
                                "invalid assignment {next}, should be \
                                 <code>:<DNA>, such as m:C"
                            )
                        } else {
                            let dna_base = parts[1]
                                .parse::<char>()
                                .map_err(|e| {
                                    anyhow!(
                                        "invalid DNA base, should be single \
                                         letter, {e}"
                                    )
                                })
                                .and_then(|raw| {
                                    DnaBase::parse(raw).map_err(|e| e.into())
                                })?;
                            let mod_code = ModCodeRepr::parse(parts[0])?;
                            debug!(
                                "assigning modification code {mod_code:?} to \
                                 {dna_base:?}"
                            );
                            acc.insert(mod_code, dna_base);
                            Ok(acc)
                        }
                    } else {
                        bail!(
                            "invalid assignment {next}, should be \
                             <code>:<DNA>, such as m:C"
                        )
                    }
                },
            )?;
            Ok(MOD_CODE_TO_DNA_BASE
                .clone()
                .into_iter()
                .chain(user_assignments.into_iter())
                .collect())
        } else {
            Ok(MOD_CODE_TO_DNA_BASE.clone())
        }
    }

    fn validate_modified_bases(
        bases: &[char],
        raw_mod_code_assignments: Option<&Vec<String>>,
    ) -> anyhow::Result<FxHashMap<ModCodeRepr, DnaBase>> {
        if bases.is_empty() {
            bail!("need to specify at least 1 modified base")
        }
        for b in bases.iter() {
            match *b {
                'A' | 'C' | 'G' | 'T' => {
                    debug!("using primary sequence base {b}");
                }
                _ => bail!("modified base needs to be A, C, G, or T."),
            }
        }

        Self::parse_raw_assignments(raw_mod_code_assignments)
    }

    pub fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()?;
        if self.control_bed_methyl.is_empty() || self.exp_bed_methyl.is_empty()
        {
            bail!("need to provide at least 1 'a' sample and 'b' sample")
        }
        let code_lookup = self.check_modified_bases()?;

        let mpb = MultiProgress::new();
        if self.suppress_progress {
            mpb.set_draw_target(indicatif::ProgressDrawTarget::hidden());
        }

        let modified_bases = self
            .modified_bases
            .iter()
            .map(|c| DnaBase::parse(*c))
            .collect::<MkResult<Vec<DnaBase>>>()?;

        if self.regions_bed.is_some()
            & (self.control_bed_methyl.len() > 1
                || self.exp_bed_methyl.len() > 1)
        {
            info!(
                "multiple samples will be combined and DMR will be performed \
                 over regions"
            );
        }

        let a_handlers = self
            .control_bed_methyl
            .iter()
            .map(|fp| BedMethylTbxIndex::from_path(fp))
            .collect::<anyhow::Result<Vec<BedMethylTbxIndex>>>()?;
        let b_handlers = self
            .exp_bed_methyl
            .iter()
            .map(|fp| HtsTabixHandler::<BedMethylLine>::from_path(fp))
            .collect::<anyhow::Result<Vec<BedMethylTbxIndex>>>()?;
        let handlers = a_handlers
            .into_iter()
            .chain(b_handlers)
            .collect::<Vec<BedMethylTbxIndex>>();

        let sample_index = MultiSampleIndex::new(
            handlers,
            code_lookup,
            self.min_valid_coverage,
            self.io_threads,
        );
        let total = self.control_bed_methyl.len() + self.exp_bed_methyl.len();
        let control_idxs =
            (0..self.control_bed_methyl.len()).collect::<Vec<usize>>();
        let exp_idxs =
            (self.control_bed_methyl.len()..total).collect::<Vec<usize>>();

        let writer: Box<dyn Write> = {
            match self.out_path.as_ref() {
                None => Box::new(BufWriter::new(std::io::stdout())),
                Some(fp) => {
                    let p = Path::new(fp);
                    create_out_directory(p)?;
                    if p.exists() && !self.force {
                        bail!("refusing to overwrite existing file {}", fp)
                    } else {
                        let fh = File::create(p)?;
                        Box::new(BufWriter::new(fh))
                    }
                }
            }
        };

        info!("reading reference FASTA at {:?}", self.reference_fasta);
        let genome_positions = GenomePositions::new_from_sequences(
            &modified_bases,
            &self.reference_fasta,
            self.mask,
            &sample_index.all_contigs(),
            &mpb,
        )?;
        let mut tab = prettytable::Table::new();
        tab.set_format(
            *prettytable::format::consts::FORMAT_NO_LINESEP_WITH_TITLE,
        );
        tab.set_titles(row!["contig", "a_contains", "b_contains", "both"]);
        let mut common_contigs = 0usize;
        for (name, _) in genome_positions.contig_sizes() {
            let a_contains =
                control_idxs.iter().any(|i| sample_index.has_contig(*i, name));
            let b_contains =
                exp_idxs.iter().any(|i| sample_index.has_contig(*i, name));
            tab.add_row(row![
                name,
                a_contains,
                b_contains,
                a_contains && b_contains
            ]);
            if a_contains && b_contains {
                common_contigs += 1;
            }
        }
        if self.careful || common_contigs == 0 {
            debug!("contig breakdown:\n{tab}");
        }
        mpb.suspend(|| {
            info!(
                "{common_contigs} common sequence(s) between FASTA and both \
                 samples"
            );
        });

        let batch_size =
            self.batch_size.as_ref().map(|x| *x).unwrap_or_else(|| {
                (self.threads as f32 * 1.5f32).floor() as usize
            });

        if self.is_single_site() {
            info!("running single-site analysis");
            let linear_transitions = if self.fine_grained {
                false
            } else {
                !self.log_transition_decay
            };
            return SingleSiteDmrAnalysis::new(
                sample_index,
                genome_positions,
                self.cap_coverages,
                self.control_bed_methyl.len(),
                self.exp_bed_methyl.len(),
                batch_size,
                self.interval_size,
                self.prior.as_ref(),
                self.max_coverages.as_ref(),
                self.delta,
                self.n_sample_records,
                self.header,
                self.segmentation_fp.as_ref(),
                mpb.clone(),
                &pool,
            )?
            .run(
                pool,
                self.max_gap_size,
                self.dmr_prior,
                self.diff_stay,
                self.significance_factor,
                self.decay_distance,
                linear_transitions,
                writer,
            );
        }

        let sample_index = Arc::new(sample_index);
        let genome_positions = Arc::new(genome_positions);

        let regions_of_interest =
            if let Some(roi_bed) = self.regions_bed.as_ref() {
                let rois = parse_roi_bed(roi_bed).with_context(|| {
                    format!("failed to parse supplied regions at {roi_bed:?}")
                })?;
                info!("loaded {} regions", rois.len());
                rois
            } else {
                unreachable!(
                    "regions should always be available unless we're doing \
                     single-site analysis"
                )
            };

        info!("loading {batch_size} regions at a time");

        let pb = mpb.add(get_master_progress_bar(regions_of_interest.len()));
        pb.set_message("regions processed");
        let failures = mpb.add(get_ticker());
        failures.set_message("regions failed to process");
        let batch_failures = mpb.add(get_ticker());
        batch_failures.set_message("failed batches");

        let dmr_interval_iter = RoiIter::new(
            control_idxs.as_slice(),
            exp_idxs.as_slice(),
            "a",
            "b",
            sample_index.clone(),
            regions_of_interest,
            batch_size,
            self.handle_missing,
            genome_positions.clone(),
            &mpb,
        )?;

        let (success_count, region_errors) = run_pairwise_dmr(
            dmr_interval_iter,
            sample_index.clone(),
            pool,
            writer,
            pb,
            self.header,
            "a",
            "b",
            failures.clone(),
            batch_failures.clone(),
            mpb.clone(),
        )?;

        mpb.suspend(|| {
            info!(
                "{} regions processed successfully and {} regions failed",
                success_count,
                failures.position()
            );
            if !region_errors.is_empty() {
                let tab = format_errors_table(&region_errors);
                error!("region errors:\n{tab}");
            }
        });

        Ok(())
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct MultiSampleDmr {
    /// Two or more named samples to compare. Two arguments are required <path>
    /// <name>. This option should be repeated at least two times. When two
    /// samples have the same name, they will be combined.
    #[clap(help_heading = "Sample Options")]
    #[arg(short = 's', long = "sample", num_args = 2)]
    samples: Vec<String>,
    /// BED file of regions over which to compare methylation levels. Should be
    /// tab-separated (spaces allowed in the "name" column). Requires
    /// chrom, chromStart and chromEnd. The Name column is optional. Strand
    /// is currently ignored.
    #[clap(help_heading = "Sample Options")]
    #[arg(long, short = 'r', alias = "regions")]
    regions_bed: PathBuf,
    /// Include header in output
    #[clap(help_heading = "Output Options")]
    #[arg(long, alias = "with-header", default_value_t = false)]
    header: bool,
    /// Directory to place output DMR results in BED format.
    #[clap(help_heading = "Output Options")]
    #[arg(short = 'o', long)]
    out_dir: PathBuf,
    /// Prefix files in directory with this label
    #[clap(help_heading = "Output Options")]
    #[arg(short = 'p', long)]
    prefix: Option<String>,
    /// Path to reference fasta for the pileup.
    #[clap(help_heading = "Sample Options")]
    #[arg(long = "ref")]
    reference_fasta: PathBuf,
    /// Bases to use to calculate DMR, may be multiple. For example, to
    /// calculate differentially methylated regions using only cytosine
    /// modifications use --base C.
    #[clap(help_heading = "Sample Options")]
    #[arg(short, long="base", alias = "modified-bases", action=clap::ArgAction::Append)]
    modified_bases: Vec<char>,
    /// Extra assignments of modification codes to their respective primary
    /// bases. In general, modkit dmr will use the SAM specification to
    /// know which modification codes are appropriate to use for a given
    /// primary base. For example "h" is the code for 5hmC, so is appropriate
    /// for cytosine bases, but not adenine bases. However, if your
    /// bedMethyl file contains custom codes or codes that are not part of
    /// the specification, you can specify which primary base they
    /// belong to here with --assign-code x:C meaning associate modification
    /// code "x" with cytosine (C) primary sequence bases. If a code is
    /// encountered that is not part of the specification, the bedMethyl
    /// record will not be used, this will be logged.
    #[clap(help_heading = "Sample Options")]
    #[arg(long="assign-code", action=clap::ArgAction::Append)]
    mod_code_assignments: Option<Vec<String>>,
    /// File to write logs to, it's recommended to use this option.
    #[clap(help_heading = "Logging Options")]
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Number of threads to use.
    #[clap(help_heading = "Compute Options")]
    #[arg(short = 't', long, default_value_t = 4)]
    threads: usize,
    /// Number of threads to use when for decompression.
    #[clap(help_heading = "Compute Options")]
    #[arg(long, default_value_t = 4)]
    io_threads: usize,
    /// Respect soft masking in the reference FASTA.
    #[clap(help_heading = "Sample Options")]
    #[arg(long, short = 'k', default_value_t = false)]
    mask: bool,
    /// Don't show progress bars
    #[clap(help_heading = "Logging Options")]
    #[arg(long, default_value_t = false)]
    suppress_progress: bool,
    /// Force overwrite of output file, if it already exists.
    #[clap(help_heading = "Output Options")]
    #[arg(short = 'f', long, default_value_t = false)]
    force: bool,
    /// How to handle regions found in the `--regions` BED file.
    /// quiet => ignore regions that are not found in the tabix header
    /// warn => log (debug) regions that are missing
    /// fatal => log (error) and exit the program when a region is missing.
    #[clap(help_heading = "Logging Options")]
    #[arg(long="missing", requires = "regions_bed", default_value_t=HandleMissing::quiet)]
    handle_missing: HandleMissing,
    /// Minimum valid coverage required to use an entry from a bedMethyl. See
    /// the help for pileup for the specification and description of valid
    /// coverage.
    #[clap(help_heading = "Sample Options")]
    #[arg(long, alias = "min-coverage", default_value_t = 0)]
    min_valid_coverage: u64,
}

impl MultiSampleDmr {
    fn get_writer(
        &self,
        a_name: &str,
        b_name: &str,
    ) -> anyhow::Result<Box<BufWriter<File>>> {
        let fp = if let Some(p) = self.prefix.as_ref() {
            self.out_dir.join(format!("{}_{}_{}.bed", p, a_name, b_name))
        } else {
            self.out_dir.join(format!("{}_{}.bed", a_name, b_name))
        };
        if fp.exists() && !self.force {
            bail!(
                "refusing to overwrite {:?}",
                fp.to_str().unwrap_or("failed decode")
            )
        } else {
            let fh = File::create(fp.clone()).with_context(|| {
                format!("failed to make output file at {fp:?}")
            })?;
            Ok(Box::new(BufWriter::new(fh)))
        }
    }

    pub fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        if !self.out_dir.exists() {
            info!("creating directory at {:?}", &self.out_dir);
            std::fs::create_dir_all(&self.out_dir)?;
        }
        let code_lookup = PairwiseDmr::validate_modified_bases(
            &self.modified_bases,
            self.mod_code_assignments.as_ref(),
        )?;

        let handlers = self
            .samples
            .chunks(2)
            .enumerate()
            .filter_map(|(i, raw)| {
                if raw.len() != 2 {
                    error!(
                        "illegal sample pair {:?}, should be length 2 of the \
                         form <path> <name>",
                        raw
                    );
                    None
                } else {
                    let fp = Path::new(raw[0].as_str()).to_path_buf();
                    let name = raw[1].to_string();
                    if fp.exists() {
                        match BedMethylTbxIndex::from_path(&fp) {
                            Ok(handler) => Some((i, name, handler)),
                            Err(e) => {
                                error!("failed to load {name}, {e}");
                                None
                            }
                        }
                    } else {
                        error!("bedMethyl for {name} at {} not found", &raw[0]);
                        None
                    }
                }
            })
            .collect::<Vec<(usize, String, BedMethylTbxIndex)>>();

        let mpb = MultiProgress::new();

        let motifs = self
            .modified_bases
            .iter()
            .map(|c| DnaBase::parse(*c))
            .collect::<MkResult<Vec<DnaBase>>>()
            .context("failed to parse modified base")?;

        let (names, handlers) = handlers.into_iter().fold(
            (HashMap::new(), Vec::new()),
            |(mut names, mut handlers), (sample_id, name, handler)| {
                names.entry(name).or_insert_with(Vec::new).push(sample_id);
                handlers.push(handler);
                (names, handlers)
            },
        );
        for (name, ids) in &names {
            if ids.len() > 1 {
                info!(
                    "sample {name} has {} replicates, they will be combined",
                    ids.len()
                );
            }
        }

        let sample_index = MultiSampleIndex::new(
            handlers,
            code_lookup,
            self.min_valid_coverage,
            self.io_threads,
        );

        let genome_positions = GenomePositions::new_from_sequences(
            &motifs,
            &self.reference_fasta,
            self.mask,
            &sample_index.all_contigs(),
            &mpb,
        )?;

        let regions_of_interest = parse_roi_bed(&self.regions_bed)?;

        let sample_index = Arc::new(sample_index);
        let genome_positions = Arc::new(genome_positions);

        info!("loaded {} regions", regions_of_interest.len());

        let chunk_size = (self.threads as f32 * 1.5f32).floor() as usize;
        info!("processing {chunk_size} regions concurrently");

        let sample_pb =
            mpb.add(get_master_progress_bar(sample_index.num_combinations()?));

        let samples = names.keys().sorted().collect::<Vec<&String>>();
        for pair in
            samples.into_iter().combinations(2).progress_with(sample_pb.clone())
        {
            let a_name = pair[0];
            let b_name = pair[1];
            let a_idxs = names.get(a_name).unwrap();
            let b_idxs = names.get(b_name).unwrap();

            sample_pb
                .set_message(format!("comparing {} and {}", a_name, b_name));
            let pb =
                mpb.add(get_subroutine_progress_bar(regions_of_interest.len()));
            pb.set_message("regions processed");
            let failures = mpb.add(get_ticker());
            failures.set_message("regions failed to process");
            let batch_failures = mpb.add(get_ticker());
            batch_failures.set_message("failed batches");

            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(self.threads)
                .build()?;

            debug!("running {a_name} as control and {b_name} as experiment");
            let mut all_region_errors = FxHashMap::default();
            match RoiIter::new(
                a_idxs,
                b_idxs,
                a_name,
                b_name,
                sample_index.clone(),
                regions_of_interest.clone(),
                chunk_size,
                self.handle_missing,
                genome_positions.clone(),
                &mpb,
            ) {
                Ok(dmr_interval_iter) => {
                    let writer = self.get_writer(a_name, b_name)?;
                    let (success_count, region_errors) = run_pairwise_dmr(
                        dmr_interval_iter,
                        sample_index.clone(),
                        pool,
                        writer,
                        pb,
                        self.header,
                        a_name,
                        b_name,
                        failures.clone(),
                        batch_failures.clone(),
                        mpb.clone(),
                    )?;
                    mpb.suspend(|| {
                        info!(
                            "{} regions processed successfully and {} regions \
                             failed for pair {} {}",
                            success_count,
                            failures.position(),
                            &a_name,
                            &b_name,
                        );
                        if !region_errors.is_empty() {
                            let tab = format_errors_table(&region_errors);
                            error!("region errors:\n{tab}");
                            all_region_errors.op_mut(region_errors);
                        }
                    });
                }
                Err(e) => {
                    mpb.suspend(|| {
                        error!(
                            "pair {} {} failed to process, {}",
                            &a_name,
                            &b_name,
                            e.to_string()
                        );
                    });
                    if self.handle_missing == HandleMissing::fail {
                        return Err(e);
                    }
                }
            }
        }

        Ok(())
    }
}
