use std::collections::HashMap;
use std::io::Write;

use crate::mod_bam::BaseModCall;
use crate::motifs::motif_bed::MotifPositionLookup;
use crate::read_ids_to_base_mod_probs::{
    PositionModCalls, ReadBaseModProfile, ReadsBaseModProfile,
};
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::util::{
    get_reference_mod_strand, Kmer, Strand, MISSING_SYMBOL, TAB,
};
use crate::writers::TsvWriter;

impl PositionModCalls {
    pub(super) fn header(with_motifs: bool) -> String {
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
            "call_prob",
            "call_code",
            "base_qual",
            "ref_kmer",
            "query_kmer",
            "canonical_base",
            "modified_primary_base",
            "fail",
            "inferred",
            "within_alignment",
            "flag",
        ];
        if with_motifs {
            fields.push("motifs")
        }
        fields.join("\t")
    }

    pub(crate) fn to_row(
        &self,
        profile: &ReadBaseModProfile,
        chrom_name: Option<&String>,
        caller: &MultipleThresholdModCaller,
        reference_seqs: &HashMap<String, Vec<u8>>,
        pass_only: bool,
        skip_inferred: bool,
        motif_position_lookup: Option<&MotifPositionLookup>,
        with_motifs: bool,
    ) -> Option<String> {
        let filtered = caller.call(&self.canonical_base, &self.base_mod_probs)
            == BaseModCall::Filtered;
        let inferred = self.base_mod_probs.inferred_unmodified;
        let motif_hits = motif_position_lookup.and_then(|lu| {
            match (self.ref_position, profile.chrom_id, self.alignment_strand) {
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
        if filtered && pass_only {
            return None;
        }
        if inferred && skip_inferred {
            return None;
        }

        let missing = ".".to_string();
        let chrom_name_label = chrom_name.unwrap_or(&missing).to_owned();
        let forward_read_position = self.query_position;
        let ref_position = self.ref_position.unwrap_or(-1);
        let mod_strand = self.mod_strand.to_char();
        let ref_strand =
            self.alignment_strand.map(|x| x.to_char()).unwrap_or('.');
        let ref_mod_strand = self
            .alignment_strand
            .map(|x| get_reference_mod_strand(self.mod_strand, x).to_char())
            .unwrap_or('.');
        let fw_soft_clipped_start = self.num_soft_clipped_start;
        let fw_soft_clipped_end = self.num_soft_clipped_end;
        let (mod_call_prob, mod_call_code) =
            match self.base_mod_probs.argmax_base_mod_call() {
                BaseModCall::Canonical(p) => (p, "-".to_string()),
                BaseModCall::Modified(p, code) => (p, code.to_string()),
                BaseModCall::Filtered => {
                    unreachable!("argmax should not output filtered calls")
                }
            };
        let read_length = self.read_length;
        let base_qual = self.q_base;
        let query_kmer = format!("{}", self.query_kmer);
        let ref_kmer = if let Some(ref_pos) = self.ref_position {
            if ref_pos < 0 {
                None
            } else {
                reference_seqs.get(&chrom_name_label).map(|s| {
                    Kmer::from_seq(s, ref_pos as usize, self.query_kmer.size)
                        .to_string()
                })
            }
        } else {
            None
        };
        let ref_kmer_rep = ref_kmer.as_ref().unwrap_or(&missing);
        let canonical_base = self.canonical_base.char();
        let modified_primary_base = if self.mod_strand == Strand::Negative {
            self.canonical_base.complement().char()
        } else {
            self.canonical_base.char()
        };
        let within_alignment = chrom_name.is_some() && self.within_alignment();

        let mut s = format!(
            "\
            {}{TAB}\
            {forward_read_position}{TAB}\
            {ref_position}{TAB}\
            {chrom_name_label}{TAB}\
            {mod_strand}{TAB}\
            {ref_strand}{TAB}\
            {ref_mod_strand}{TAB}\
            {fw_soft_clipped_start}{TAB}\
            {fw_soft_clipped_end}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {read_length}{TAB}\
            {mod_call_prob}{TAB}\
            {mod_call_code}{TAB}\
            {base_qual}{TAB}\
            {ref_kmer_rep}{TAB}\
            {query_kmer}{TAB}\
            {canonical_base}{TAB}\
            {modified_primary_base}{TAB}\
            {filtered}{TAB}\
            {inferred}{TAB}\
            {within_alignment}{TAB}\
            {}",
            &profile.record_name,
            profile.alignment_start.map(|x| x as i64).unwrap_or(-1i64),
            profile.alignment_end.map(|x| x as i64).unwrap_or(-1i64),
            &profile.flag,
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
        Some(s)
    }
}

pub(crate) trait OutwriterWithMemory<T> {
    fn write(
        &mut self,
        item: T,
        motif_position_lookup: Option<&MotifPositionLookup>,
    ) -> anyhow::Result<u64>;
    fn num_reads(&self) -> usize;
}

pub struct TsvWriterWithContigNames<W: Write, C> {
    tsv_writer: TsvWriter<W>,
    tid_to_name: HashMap<u32, String>,
    name_to_seq: HashMap<String, Vec<u8>>,
    number_of_written_reads: usize,
    caller: C,
    pass_only: bool,
    with_motifs: bool,
}

impl<W: Write> TsvWriterWithContigNames<W, ()> {
    pub(crate) fn new(
        output_writer: TsvWriter<W>,
        tid_to_name: HashMap<u32, String>,
        name_to_seq: HashMap<String, Vec<u8>>,
        with_motifs: bool,
    ) -> anyhow::Result<Self> {
        Ok(Self {
            tsv_writer: output_writer,
            tid_to_name,
            name_to_seq,
            number_of_written_reads: 0,
            caller: (),
            pass_only: false,
            with_motifs,
        })
    }
}

impl<W: Write> OutwriterWithMemory<ReadsBaseModProfile>
    for TsvWriterWithContigNames<W, ()>
{
    fn write(
        &mut self,
        item: ReadsBaseModProfile,
        motif_position_lookup: Option<&MotifPositionLookup>,
    ) -> anyhow::Result<u64> {
        let mut rows_written = 0u64;
        for profile in item.profiles.iter() {
            let chrom_name = profile
                .chrom_id
                .and_then(|chrom_id| self.tid_to_name.get(&chrom_id))
                .map(|x| x.as_str())
                .unwrap_or(MISSING_SYMBOL);
            for mod_profile in profile.iter_profiles() {
                let row = mod_profile.to_row(
                    &profile.record_name,
                    chrom_name,
                    profile.chrom_id,
                    profile.alignment_start,
                    profile.alignment_end,
                    &self.name_to_seq,
                    profile.flag,
                    motif_position_lookup,
                    self.with_motifs,
                );
                self.tsv_writer.write(row.as_bytes())?;
                rows_written += 1;
            }
            self.number_of_written_reads += 1;
        }
        Ok(rows_written)
    }

    fn num_reads(&self) -> usize {
        self.number_of_written_reads
    }
}

impl<W: Write> TsvWriterWithContigNames<W, MultipleThresholdModCaller> {
    pub(crate) fn new_with_caller(
        output_writer: TsvWriter<W>,
        tid_to_name: HashMap<u32, String>,
        name_to_seq: HashMap<String, Vec<u8>>,
        caller: MultipleThresholdModCaller,
        pass_only: bool,
        with_motifs: bool,
    ) -> anyhow::Result<Self> {
        Ok(Self {
            tsv_writer: output_writer,
            tid_to_name,
            name_to_seq,
            number_of_written_reads: 0,
            caller,
            pass_only,
            with_motifs,
        })
    }
}

impl<W: Write> OutwriterWithMemory<ReadsBaseModProfile>
    for TsvWriterWithContigNames<W, MultipleThresholdModCaller>
{
    fn write(
        &mut self,
        item: ReadsBaseModProfile,
        motif_position_lookup: Option<&MotifPositionLookup>,
    ) -> anyhow::Result<u64> {
        let mut rows_written = 0u64;
        for profile in item.profiles.iter() {
            let chrom_name = profile
                .chrom_id
                .and_then(|chrom_id| self.tid_to_name.get(&chrom_id));
            let position_calls = PositionModCalls::from_profile(&profile);
            for call in position_calls {
                call.to_row(
                    profile,
                    chrom_name,
                    &self.caller,
                    &self.name_to_seq,
                    self.pass_only,
                    false,
                    motif_position_lookup,
                    self.with_motifs,
                )
                .map(|s| self.tsv_writer.write(s.as_bytes()))
                .transpose()?;
                rows_written += 1;
            }
            self.number_of_written_reads += 1;
        }
        Ok(rows_written)
    }

    fn num_reads(&self) -> usize {
        self.number_of_written_reads
    }
}
