use crate::common::{
    check_against_expected_text_file, check_legal_csv, run_modkit,
};

mod common;

#[test]
fn test_dmr_helps() {
    let _ = run_modkit(&["dmr", "pair", "--help"])
        .expect("failed to run modkit dmr pair help");
    let _ = run_modkit(&["dmr", "multi", "--help"])
        .expect("failed to run modkit dmr multi help");
}

#[test]
fn test_dmr_regression() {
    let out_bed = std::env::temp_dir().join("test_dmr_regression.bed");
    let _ = run_modkit(&[
        "dmr",
        "pair",
        "-a",
        "tests/resources/\
         lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "-b",
        "tests/resources/\
         lung_00733-m_primary-tumour_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "-o",
        out_bed.to_str().unwrap(),
        "-r",
        "tests/resources/cpg_chr20_with_orig_names_selection.bed",
        "--ref",
        "tests/resources/GRCh38_chr20.fa",
        "--header",
        "-f",
        "--base",
        "C",
    ])
    .expect("failed to run modkit dmr");

    check_legal_csv::<{ '\t' as u8 }>(&out_bed);
    check_against_expected_text_file(
        out_bed.to_str().unwrap(),
        "tests/resources/test_output_chr20-2.bed",
    );

    let out_bed =
        std::env::temp_dir().join("foo").join("test_dmr_regression_2.bed");

    let _ = run_modkit(&[
        "dmr",
        "pair",
        "-a",
        "tests/resources/\
         lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "-b",
        "tests/resources/\
         lung_00733-m_primary-tumour_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        "-o",
        out_bed.to_str().unwrap(),
        "-r",
        "tests/resources/cpg_chr20_with_orig_names_selection.bed",
        "--ref",
        "tests/resources/GRCh38_chr20.fa",
        "-f",
        "--header",
        "--base",
        "C",
    ])
    .expect("failed to run modkit dmr");

    check_legal_csv::<{ '\t' as u8 }>(&out_bed);
    check_against_expected_text_file(
        out_bed.to_str().unwrap(),
        "tests/resources/test_output_chr20-2.bed",
    );
}

// todo
//  test pair with explicit index
//  test multi
