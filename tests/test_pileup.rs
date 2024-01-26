use anyhow::Context;
use itertools::Itertools;
use rust_htslib::bam;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use common::{check_against_expected_text_file, run_modkit};
use mod_kit::dmr::bedmethyl::BedMethylLine;
use mod_kit::mod_base_code::{ModCodeRepr, METHYL_CYTOSINE};

mod common;

#[test]
fn test_pileup_help() {
    let pileup_help_args = ["pileup", "--help"];
    let _out = run_modkit(&pileup_help_args).unwrap();
}

#[test]
fn test_pileup_no_filt() {
    let temp_file = std::env::temp_dir().join("test_pileup_nofilt.bed");
    let args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--no-filtering",
        "--only-tabs",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ];

    run_modkit(&args).unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "tests/resources/modbam.modpileup_nofilt.methyl.bed",
    );
}

#[test]
fn test_pileup_with_filt() {
    let temp_file = std::env::temp_dir().join("test_pileup_withfilt.bed");
    let args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "-f",
        "1.0",
        "-p",
        "0.25",
        "--only-tabs",
        "--seed",
        "42",
        "--include-unmapped",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ];

    run_modkit(&args).unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "tests/resources/modbam.modpileup_filt025.methyl.bed",
    );
}

#[test]
fn test_pileup_combine() {
    let test_adjusted_bam = std::env::temp_dir().join("test_combined.bed");
    let pileup_args = [
        "pileup",
        "--combine-mods",
        "--no-filtering",
        "--only-tabs",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        test_adjusted_bam.to_str().unwrap(),
    ];
    run_modkit(&pileup_args).unwrap();
    assert!(test_adjusted_bam.exists());

    check_against_expected_text_file(
        test_adjusted_bam.to_str().unwrap(),
        "tests/resources/modbam.modpileup_combined.methyl.bed",
    );
}

#[test]
fn test_pileup_collapse() {
    let test_collapsed_bam = std::env::temp_dir().join("test_collapsed.bam");
    let test_collapsed_bed = std::env::temp_dir().join("test_collapsed.bed");
    let test_restricted_bed = std::env::temp_dir().join("test_restricted.bed");

    let collapse_args = [
        "adjust-mods",
        "--ignore",
        "h",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        test_collapsed_bam.to_str().unwrap(),
    ];
    run_modkit(&collapse_args).unwrap();
    assert!(test_collapsed_bam.exists());
    bam::index::build(
        test_collapsed_bam.clone(),
        None,
        bam::index::Type::Bai,
        1,
    )
    .unwrap();

    let pileup_args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--no-filtering",
        test_collapsed_bam.to_str().unwrap(),
        test_collapsed_bed.to_str().unwrap(),
    ];
    run_modkit(&pileup_args).unwrap();
    assert!(test_collapsed_bed.exists());

    let pileup_args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--ignore",
        "h",
        "--no-filtering",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        test_restricted_bed.to_str().unwrap(),
    ];
    run_modkit(&pileup_args).unwrap();
    assert!(test_restricted_bed.exists());
    check_against_expected_text_file(
        test_restricted_bed.to_str().unwrap(),
        test_collapsed_bed.to_str().unwrap(),
    );
}

#[test]
fn test_pileup_no_mod_calls() {
    let empty_bedfile =
        std::env::temp_dir().join("test_pileup_no_mod_calls_outbed.bed");
    let args = [
        "pileup",
        "--no-filtering",
        "tests/resources/empty-tags.sorted.bam",
        empty_bedfile.to_str().unwrap(),
    ];

    run_modkit(&args).unwrap();

    let reader = BufReader::new(File::open(empty_bedfile).unwrap());
    let lines = reader.lines().collect::<Vec<Result<String, _>>>();
    assert_eq!(lines.len(), 0);
}

#[test]
fn test_pileup_old_tags() {
    let updated_file =
        std::env::temp_dir().join("test_pileup_old_tags_updated.bam");
    run_modkit(&[
        "update-tags",
        "tests/resources/HG002_small.ch20._other.sorted.bam",
        "--mode",
        "ambiguous",
        updated_file.to_str().unwrap(),
    ])
    .unwrap();
    assert!(updated_file.exists());
    bam::index::build(updated_file.clone(), None, bam::index::Type::Bai, 1)
        .unwrap();

    let out_file = std::env::temp_dir().join("test_pileup_old_tags.bed");
    run_modkit(&[
        "pileup",
        "--no-filtering",
        "--only-tabs",
        updated_file.to_str().unwrap(),
        out_file.to_str().unwrap(),
    ])
    .unwrap();
    assert!(out_file.exists());
    check_against_expected_text_file(
        out_file.to_str().unwrap(),
        "tests/resources/pileup-old-tags-regressiontest.methyl.bed",
    );
}

#[test]
fn test_pileup_with_region() {
    let temp_file = std::env::temp_dir().join("test_pileup_with_region.bed");
    let args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--no-filtering",
        "--region",
        "oligo_1512_adapters:0-50",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ];

    run_modkit(&args).unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "tests/resources/modbam.modpileup_nofilt_oligo_1512_adapters_10_50.bed",
    );
}

#[test]
fn test_pileup_duplex_reads() {
    let temp_file = std::env::temp_dir().join("test_pileup_duplex_reads.bed");
    run_modkit(&[
        "pileup",
        "tests/resources/duplex_modbam.sorted.bam",
        temp_file.to_str().unwrap(),
        "--region",
        "chr17",
        "--no-filtering",
    ])
    .unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "tests/resources/duplex_modbam_pileup_nofilt.bed",
    );
}

#[test]
fn test_pileup_cpg_motif_filtering() {
    let temp_file = std::env::temp_dir().join("test_cpg_motif_filtering.bed");
    run_modkit(&[
        "pileup",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
        "--no-filtering",
        "--cpg",
        "--ref",
        "tests/resources/CGI_ladder_3.6kb_ref.fa",
    ])
    .unwrap();
    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "tests/resources/bc_anchored_10_reads_nofilt_cg_motif.bed",
    );
}

#[test]
fn test_pileup_cpg_motif_filtering_strand_combine() {
    let temp_file = std::env::temp_dir()
        .join("test_cpg_motif_filtering_strand_combine.bed");
    let interval_sizes =
        ["10", "88", "89", "90", "91", "92", "93", "94", "10000"];
    for interval_size in interval_sizes {
        run_modkit(&[
            "pileup",
            "tests/resources/bc_anchored_10_reads.sorted.bam",
            temp_file.to_str().unwrap(),
            "--no-filtering",
            "-i",
            interval_size,
            "--cpg",
            "--combine-strands",
            "--ref",
            "tests/resources/CGI_ladder_3.6kb_ref.fa",
        ])
        .unwrap();
        check_against_expected_text_file(
            temp_file.to_str().unwrap(),
            "tests/resources/\
             bc_anchored_10_reads_nofilt_cg_motif_strand_combine.bed",
        );
    }
}

#[test]
fn test_pileup_presets_traditional_same_as_options() {
    let preset_temp_file = std::env::temp_dir()
        .join("test_presets_traditional_same_as_options.bed");
    let options_temp_file = std::env::temp_dir()
        .join("test_presets_traditional_same_as_options2.bed");

    run_modkit(&[
        "pileup",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        preset_temp_file.to_str().unwrap(),
        "--no-filtering",
        "--preset",
        "traditional",
        "--ref",
        "tests/resources/CGI_ladder_3.6kb_ref.fa",
    ])
    .unwrap();

    run_modkit(&[
        "pileup",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        options_temp_file.to_str().unwrap(),
        "--cpg",
        "--no-filtering",
        "--ignore",
        "h",
        "--combine-strands",
        "--ref",
        "tests/resources/CGI_ladder_3.6kb_ref.fa",
    ])
    .unwrap();
    check_against_expected_text_file(
        preset_temp_file.to_str().unwrap(),
        options_temp_file.to_str().unwrap(),
    );
}

#[test]
fn test_pileup_duplicated_reads_ignored() {
    let control_fp =
        std::env::temp_dir().join("test_duplicated_reads_ignored_control.bed");
    let test_fp =
        std::env::temp_dir().join("test_duplicated_reads_ignored_marked.bed");
    run_modkit(&[
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--no-filtering",
        "--only-tabs",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        control_fp.to_str().unwrap(),
    ])
    .unwrap();

    run_modkit(&[
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--no-filtering",
        "--only-tabs",
        "tests/resources/duplicated.marked.fixed.bam",
        test_fp.to_str().unwrap(),
    ])
    .unwrap();

    check_against_expected_text_file(
        control_fp.to_str().unwrap(),
        test_fp.to_str().unwrap(),
    );
}

#[test]
fn test_pileup_edge_filter_regression() {
    let adjusted_bam =
        std::env::temp_dir().join("test_pileup_edge_filter_adjusted.bam");
    let edge_filter_bed = std::env::temp_dir()
        .join("test_pileup_edge_filter_edge_filter_50.pileup.bed");
    let edge_filter_bed_2 = std::env::temp_dir()
        .join("test_pileup_edge_filter_edge_filter_50_2.pileup.bed");

    run_modkit(&[
        "pileup",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        edge_filter_bed.to_str().unwrap(),
        "--no-filtering",
        "--edge-filter",
        "50",
    ])
    .context("test_pileup_edge_filter_regression failed to make bedMethyl")
    .unwrap();
    check_against_expected_text_file(
        edge_filter_bed.to_str().unwrap(),
        "tests/resources/bc_anchored_10_reads_edge_filter50.bed",
    );

    run_modkit(&[
        "adjust-mods",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        adjusted_bam.to_str().unwrap(),
        "--edge-filter",
        "50",
    ])
    .context("test_pileup_edge_filter_regression failed to run adjust-mods")
    .unwrap();

    bam::index::build(adjusted_bam.clone(), None, bam::index::Type::Bai, 1)
        .unwrap();

    run_modkit(&[
        "pileup",
        adjusted_bam.to_str().unwrap(),
        edge_filter_bed_2.to_str().unwrap(),
        "--no-filtering",
        "--edge-filter",
        "50",
    ])
    .context(
        "test_pileup_edge_filter_regression failed to make bedMethyl on \
         adjusted bam",
    )
    .unwrap();
    check_against_expected_text_file(
        edge_filter_bed.to_str().unwrap(),
        edge_filter_bed_2.to_str().unwrap(),
    );
}

#[test]
fn test_pileup_edge_filter_asymmetric_regression() {
    let adjusted_bam = std::env::temp_dir()
        .join("test_pileup_edge_filter_asymmetric_regression.bam");
    let edge_filter_bed = std::env::temp_dir().join(
        "test_pileup_edge_filter_asymmetric_regression_filter_50.pileup.bed",
    );
    let edge_filter_bed_2 = std::env::temp_dir()
        .join("test_pileup_edge_filter_asymmetric_regression_50_2.pileup.bed");

    run_modkit(&[
        "pileup",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        edge_filter_bed.to_str().unwrap(),
        "--no-filtering",
        "--edge-filter",
        "50,50",
    ])
    .context(
        "test_pileup_edge_filter_asymmetric_regression failed to make \
         bedMethyl with 50,50",
    )
    .unwrap();
    check_against_expected_text_file(
        edge_filter_bed.to_str().unwrap(),
        "tests/resources/bc_anchored_10_reads_edge_filter50.bed",
    );

    run_modkit(&[
        "pileup",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        edge_filter_bed.to_str().unwrap(),
        "--no-filtering",
        "--edge-filter",
        "50,0",
    ])
    .context(
        "test_pileup_edge_filter_asymmetric_regression failed to make \
         bedMethyl with 50,0",
    )
    .unwrap();
    check_against_expected_text_file(
        edge_filter_bed.to_str().unwrap(),
        "tests/resources/bc_anchored_10_reads_edge_filter50-0.bed",
    );

    run_modkit(&[
        "adjust-mods",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        adjusted_bam.to_str().unwrap(),
        "--edge-filter",
        "50,0",
    ])
    .context(
        "test_pileup_edge_filter_asymmetric_regression failed to make adjust \
         mods BAM",
    )
    .unwrap();

    bam::index::build(adjusted_bam.clone(), None, bam::index::Type::Bai, 1)
        .unwrap();

    run_modkit(&[
        "pileup",
        adjusted_bam.to_str().unwrap(),
        edge_filter_bed_2.to_str().unwrap(),
        "--no-filtering",
    ])
    .context(
        "test_pileup_edge_filter_asymmetric_regression failed to make pileup \
         on adjusted BAM",
    )
    .unwrap();

    check_against_expected_text_file(
        edge_filter_bed_2.to_str().unwrap(),
        "tests/resources/bc_anchored_10_reads_edge_filter50-0.bed",
    );
}

#[test]
fn test_pileup_partition_tags_partitioned() {
    let tmp_dir =
        std::env::temp_dir().join("test_pileup_partition_tags_partitioned");
    let control_file =
        std::env::temp_dir().join("test_pileup_partition_tags_control.bed");

    // control BED, all of the partitioned BED files should be the same as this
    // one
    run_modkit(&[
        "pileup",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        control_file.to_str().unwrap(),
        "--no-filtering",
    ])
    .context("failed to run modkit on control")
    .unwrap();

    // run partitioned on HP and RG tags. This test file has 2 HP tags {1, 2}
    // and 3 read groups {A, B, C}. So we expect 6 files, all the same as the
    // control
    run_modkit(&[
        "pileup",
        "tests/resources/bc_anchored_10_reads.haplotyped.sorted.bam",
        tmp_dir.to_str().unwrap(),
        "--partition-tag",
        "RG",
        "--partition-tag",
        "HP",
        "--no-filtering",
    ])
    .context("failed to run modkit with partition tags")
    .unwrap();

    let mut count = 0;
    for result in tmp_dir.read_dir().unwrap() {
        let dir_entry = result.unwrap().path();
        check_against_expected_text_file(
            dir_entry.to_str().unwrap(),
            control_file.to_str().unwrap(),
        );
        count += 1;
    }
    assert_eq!(count, 6);
}

#[test]
fn test_pileup_partition_tags_bedgraph() {
    let tmp_dir = std::env::temp_dir()
        .join("test_pileup_partition_tags_bedgraph_partitioned");
    let control_dir = std::env::temp_dir()
        .join("test_pileup_partition_tags_bedgraph_control");

    let collect_bedgraph_files =
        |dir_path: &PathBuf| -> std::io::Result<Vec<PathBuf>> {
            dir_path.read_dir().map(|read_dir| {
                read_dir
                    .filter_map(|dir| match dir {
                        Ok(dir) => {
                            if dir.path().extension().and_then(|fp| fp.to_str())
                                == Some("bedgraph")
                            {
                                Some(dir.path())
                            } else {
                                None
                            }
                        }
                        Err(_) => None,
                    })
                    .collect::<Vec<PathBuf>>()
            })
        };

    // control BED, all of the partitioned BED files should be the same as this
    // one
    run_modkit(&[
        "pileup",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        control_dir.to_str().unwrap(),
        "--no-filtering",
        "--bedgraph",
    ])
    .context("failed to run modkit on control bedgraph")
    .unwrap();

    let control_bedgraph_files = collect_bedgraph_files(&control_dir)
        .unwrap()
        .into_iter()
        .map(|fp| {
            let file_name = fp.file_name().unwrap().to_str().unwrap();
            match (file_name.starts_with("h"), file_name.contains("positive")) {
                (true, true) => (('h', "positive"), fp),
                (true, false) => (('h', "negative"), fp),
                (false, true) => (('m', "positive"), fp),
                (false, false) => (('m', "negative"), fp),
            }
        })
        .collect::<HashMap<(char, &str), PathBuf>>();

    // run partitioned on HP and RG tags. This test file has 2 HP tags {1, 2}
    // and 3 read groups {A, B, C}. So we expect 6 files, all the same as the
    // control
    run_modkit(&[
        "pileup",
        "tests/resources/bc_anchored_10_reads.haplotyped.sorted.bam",
        tmp_dir.to_str().unwrap(),
        "--partition-tag",
        "RG",
        "--partition-tag",
        "HP",
        "--no-filtering",
        "--bedgraph",
    ])
    .context("failed to run modkit with partition tags")
    .unwrap();

    let mut count = 0;
    for result in tmp_dir.read_dir().unwrap() {
        let dir_entry = result.unwrap().path();
        if dir_entry.extension().and_then(|s| s.to_str()) != Some("bedgraph") {
            continue;
        }
        let file_name = dir_entry.file_name().unwrap().to_str().unwrap();
        let stripped = file_name.replace(".bedgraph", "");
        let parts = stripped.split('_').collect::<Vec<&str>>();
        let mod_code = parts[2].parse::<char>().unwrap();
        let strand = parts[3];
        let key = (mod_code, strand);
        let file_to_compare_to = control_bedgraph_files.get(&key).unwrap();
        check_against_expected_text_file(
            dir_entry.to_str().unwrap(),
            file_to_compare_to.to_str().unwrap(),
        );
        count += 1;
    }
    assert_eq!(count, 24);
}

#[test]
fn test_pileup_with_filt_position_filter() {
    let temp_file =
        std::env::temp_dir().join("test_pileup_with_filt_position_filter.bed");
    run_modkit(&[
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "-p",
        "0.25",
        "--include-positions",
        "tests/resources/CGI_ladder_3.6kb_ref_include_positions.bed",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ])
    .unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "tests/resources/modbam.modpileup_filt_positions_025.methyl.bed",
    );
}

#[test]
fn test_pileup_with_filter_positions_and_traditional() {
    let temp_file = std::env::temp_dir()
        .join("test_pileup_with_filter_positions_and_traditional.bed");
    run_modkit(&[
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "-p",
        "0.25",
        "--preset",
        "traditional",
        "--ref",
        "tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--include-positions",
        "tests/resources/CGI_ladder_3.6kb_ref_include_positions.bed",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ])
    .unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "tests/resources/modbam.modpileup_filt_positions_025_traditional.\
         methyl.bed",
    );
}

#[test]
fn test_pileup_partition_tags_combine_strands() {
    let exp_dir = std::env::temp_dir()
        .join("test_pileup_partition_tags_combine_strands_partitioned");
    let control_file = std::env::temp_dir()
        .join("test_pileup_partition_tags_combine_strands_control.bed");
    run_modkit(&[
        "pileup",
        "tests/resources/bc_anchored_10_reads.sorted.bam",
        control_file.to_str().unwrap(),
        "--combine-strands",
        "--ref",
        "tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--cpg",
        "--no-filtering",
    ])
    .context("failed to run modkit on control")
    .unwrap();
    run_modkit(&[
        "pileup",
        "tests/resources/bc_anchored_10_reads.haplotyped.sorted.bam",
        exp_dir.to_str().unwrap(),
        "--partition-tag",
        "RG",
        "--partition-tag",
        "HP",
        "--combine-strands",
        "--ref",
        "tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--cpg",
        "--no-filtering",
    ])
    .context("failed to run modkit with partition tags")
    .unwrap();
    let mut count = 0;
    for result in exp_dir.read_dir().unwrap() {
        let dir_entry = result.unwrap().path();
        check_against_expected_text_file(
            dir_entry.to_str().unwrap(),
            control_file.to_str().unwrap(),
        );
        count += 1;
    }
    assert_eq!(count, 6);
}

#[test]
fn test_pileup_motifs_cg0_cgcg2() {
    let temp_file =
        std::env::temp_dir().join("test_pileup_motifs_cg0_cgcg2.bed");
    run_modkit(&[
        "pileup",
        "tests/resources/CG_5mC_20230207_1700_6A_PAG66026_3c0abf27_oligo_741_adapters_modcalls_0th_sort_10_reads.bam",
        temp_file.to_str().unwrap(),
        "--motif", "CG", "0",
        "--motif", "CGCG", "2",
        "--no-filtering",
        "--ref", "tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--region", "oligo_741_adapters:22-62",
    ])
        .unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "tests/resources/cgcg2_cg0_test1.bed",
    );

    run_modkit(&[
        "pileup",
        "tests/resources/CG_5mC_20230207_1700_6A_PAG66026_3c0abf27_oligo_741_adapters_modcalls_0th_sort_10_reads-2.bam",
        temp_file.to_str().unwrap(),
        "--motif", "CG", "0",
        "--motif", "CGCG", "2",
        "--no-filtering",
        "--ref", "tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--region", "oligo_741_adapters:22-62",
    ])
        .unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "tests/resources/cgcg2_cg0_test2.bed",
    );
}

#[test]
fn test_pileup_motifs_cg0_cgcg2_combined() {
    let temp_file =
        std::env::temp_dir().join("test_pileup_motifs_cg0_cgcg2_combined.bed");
    run_modkit(&[
        "pileup",
        "tests/resources/CG_5mC_20230207_1700_6A_PAG66026_3c0abf27_oligo_741_adapters_modcalls_0th_sort_10_reads.bam",
        temp_file.to_str().unwrap(),
        "--motif", "CG", "0",
        "--motif", "CGCG", "2",
        "--no-filtering",
        "--combine-strands",
        "--ref", "tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--region", "oligo_741_adapters:22-62",
    ])
        .unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "tests/resources/cgcg2_cg0_test1_combine_strands.bed",
    );

    run_modkit(&[
        "pileup",
        "tests/resources/CG_5mC_20230207_1700_6A_PAG66026_3c0abf27_oligo_741_adapters_modcalls_0th_sort_10_reads-2.bam",
        temp_file.to_str().unwrap(),
        "--motif", "CG", "0",
        "--motif", "CGCG", "2",
        "--no-filtering",
        "--combine-strands",
        "--ref", "tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--region", "oligo_741_adapters:22-62",
    ])
        .unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "tests/resources/cgcg2_cg0_test2_combine_strands.bed",
    );
}

#[test]
fn test_pileup_chebi_code_same_output() {
    let ord_bm_line = |a: &BedMethylLine, b: &BedMethylLine| -> Ordering {
        match a.chrom.cmp(&b.chrom) {
            Ordering::Equal => match a.start().cmp(&b.start()) {
                Ordering::Equal => a.raw_mod_code.cmp(&b.raw_mod_code),
                o @ _ => o,
            },
            o @ _ => o,
        }
    };
    let adjusted_bam =
        std::env::temp_dir().join("test_chebi_code_same_output_hmc2chEBI.bam");
    let pileup =
        std::env::temp_dir().join("test_chebi_code_same_output_pileup.bed");
    let expected_fp = "tests/resources/modbam.modpileup_nofilt.methyl.bed";
    let expected = BufReader::new(File::open(expected_fp).unwrap())
        .lines()
        .map(|l| BedMethylLine::parse(&l.unwrap()).unwrap())
        .sorted_by(|a, b| ord_bm_line(a, b))
        .collect::<Vec<BedMethylLine>>();
    for to_code in [ModCodeRepr::ChEbi(76792), ModCodeRepr::Code('c')] {
        run_modkit(&[
            "adjust-mods",
            "tests/resources/bc_anchored_10_reads.sorted.bam",
            adjusted_bam.to_str().unwrap(),
            "--convert",
            "h",
            &to_code.to_string(),
        ])
        .with_context(|| format!("failed to change 5hmC to {to_code}"))
        .unwrap();

        bam::index::build(adjusted_bam.clone(), None, bam::index::Type::Bai, 1)
            .unwrap();

        run_modkit(&[
            "pileup",
            adjusted_bam.to_str().unwrap(),
            pileup.to_str().unwrap(),
            "-i",
            "25", // use small interval to make sure chunking works
            "--no-filtering",
            "--only-tabs",
        ])
        .context("failed to generate pileup")
        .unwrap();

        let observed = BufReader::new(File::open(pileup.clone()).unwrap())
            .lines()
            .map(|l| {
                let bm = BedMethylLine::parse(&l.unwrap()).unwrap();
                if bm.raw_mod_code != METHYL_CYTOSINE {
                    assert_eq!(bm.raw_mod_code, to_code);
                    BedMethylLine::new(
                        bm.chrom,
                        bm.interval,
                        ModCodeRepr::Code('h'), /* change back so we can
                                                 * compare */
                        bm.strand,
                        bm.count_methylated,
                        bm.valid_coverage,
                    )
                } else {
                    bm
                }
            })
            .sorted_by(|a, b| ord_bm_line(a, b))
            .collect::<Vec<BedMethylLine>>();
        assert_eq!(expected, observed);
    }
}
