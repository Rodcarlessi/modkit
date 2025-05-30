# Extracting base modification information

The `modkit extract full` sub-commands will produce a table containing the base modification probabilities, the read sequence context, and optionally aligned reference information.
For `extract full` and `extract calls`, if a correct `MN` tag is found, secondary and supplementary alignments may be output with the `--allow-non-primary` flag. 
See [troubleshooting](./troubleshooting.md) for details.

The table will by default contain unmapped sections of the read (soft-clipped sections, for example). 
To only include mapped bases use the `--mapped` flag. To only include sites of interest, pass a 
BED-formatted file to the `--include-bed` option. Similarly, to exclude sites, pass a BED-formatted
file to the `--exclude` option. One caution, the files generated by `modkit extract` can be large (2-2.5x
the size of the BAM). You may want to either use the `--num-reads` option, the `--region` option, or
pre-filter the modBAM ahead of time. You can also stream the output to stdout by setting the output to `-`
or `stdout` and filter the columns before writing to disk.

## Description of output table for `extract full`

| column | name                  | description                                                                                                             | type |
|--------|-----------------------|-------------------------------------------------------------------------------------------------------------------------|------|
| 1      | read_id               | name of the read                                                                                                        | str  |
| 2      | forward_read_position | 0-based position on the forward-oriented read sequence                                                                  | int  |
| 3      | ref_position          | aligned 0-based reference sequence position, -1 means unmapped                                                          | int  |
| 4      | chrom                 | name of aligned contig, or '.' if the read is Gunmapped                                                                 | str  |
| 5      | mod_strand            | strand of the molecule the base modification is on                                                                      | str  |
| 6      | ref_strand            | strand of the reference the read is aligned to, or '.' if unmapped                                                      | str  |
| 7      | ref_mod_strand        | strand of the reference with the base modification, or '.' if unmapped                                                  | str  |
| 8      | fw_soft_clipped_start | number of bases soft clipped from the start of the forward-oriented read                                                | int  |
| 9      | fw_soft_clipped_end   | number of bases soft clipped from the end of the forward-oriented read                                                  | int  |
| 10     | alignment_start       | leftmost (i.e. smallest) aligned reference position                                                                     | int  |
| 11     | alignment_end         | rightmost (i.e. largest) aligned reference position                                                                     | int  |
| 12     | read_length           | total length of the read                                                                                                | int  |
| 13     | mod_qual              | probability of the base modification in the next column                                                                 | int  |
| 14     | mod_code              | base modification code from the MM tag                                                                                  | str  |
| 15     | base_qual             | basecall quality score (phred)                                                                                          | int  |
| 16     | ref_kmer              | reference 5-mer sequence context (center base is aligned base), '.' if unmapped                                         | str  |
| 17     | query_kmer            | read 5-mer sequence context (center base is aligned base)                                                               | str  |
| 18     | canonical_base        | canonical base from the query sequence, from the MM tag                                                                 | str  |
| 19     | modified_primary_base | primary sequence base with the modification                                                                             | str  |
| 20     | inferred              | whether the base modification call is implicit canonical                                                                | str  |
| 21     | flag                  | FLAG from alignment record                                                                                              | str  |
| 22     | motifs                | comma-separated list of reference motifs matching at this position, **only present when `--motifs` or `--cpg` is used** | str  |


# Tabulating base modification _calls_ for each read position with `extract calls`
The `modkit extract calls` command will generate a table of read-level base modification calls using the same [thresholding](./filtering.md) algorithm employed by `modkit pileup`.
The resultant table has, for each read, one row for each base modification call in that read.
If a base is called as modified then `call_code` will be the code in the `MM` tag. If the base is called as canonical the `call_code` will be `-` (`A`, `C`, `G`, and `T` are
reserved for "any modification"). The full schema of the table is below:

| column | name                  | description                                                                                                             | type |
|--------|-----------------------|-------------------------------------------------------------------------------------------------------------------------|------|
| 1      | read_id               | name of the read                                                                                                        | str  |
| 2      | forward_read_position | 0-based position on the forward-oriented read sequence                                                                  | int  |
| 3      | ref_position          | aligned 0-based reference sequence position, -1 means unmapped                                                          | int  |
| 4      | chrom                 | name of aligned contig, or '.' if unmapped                                                                              | str  |
| 5      | mod_strand            | strand of the molecule the base modification is on                                                                      | str  |
| 6      | ref_strand            | strand of the reference the read is aligned to, or '.' if unmapped                                                      | str  |
| 7      | ref_mod_strand        | strand of the reference with the base modification, or '.' if unmapped                                                  | str  |
| 8      | fw_soft_clipped_start | number of bases soft clipped from the start of the forward-oriented read                                                | int  |
| 9      | fw_soft_clipped_end   | number of bases soft clipped from the end of the forward-oriented read                                                  | int  |
| 10     | alignment_start       | leftmost (i.e. smallest) aligned reference position                                                                     | int  |
| 11     | alignment_end         | rightmost (i.e. largest) aligned reference position                                                                     | int  |
| 12     | read_length           | total length of the read                                                                                                | int  |
| 13     | call_prob             | probability of the base modification call in the next column                                                            | int  |
| 14     | call_code             | base modification call, `-` indicates a canonical call                                                                  | str  |
| 15     | base_qual             | basecall quality score (phred)                                                                                          | int  |
| 16     | ref_kmer              | reference 5-mer sequence context (center base is aligned base), '.' if unmapped                                         | str  |
| 17     | query_kmer            | read 5-mer sequence context (center base is aligned base)                                                               | str  |
| 18     | canonical_base        | canonical base from the query sequence, from the MM tag                                                                 | str  |
| 19     | modified_primary_base | primary sequence base with the modification                                                                             | str  |
| 20     | fail                  | true if the base modification call fell below the pass threshold                                                        | str  |
| 21     | inferred              | whether the base modification call is implicit canonical                                                                | str  |
| 22     | within_alignment      | when alignment information is present, is this base aligned to the reference                                            | str  |
| 23     | flag                  | FLAG from alignment record                                                                                              | str  |
| 24     | motifs                | comma-separated list of reference motifs matching at this position, **only present when `--motifs` or `--cpg` is used** | str  |


## Note on implicit base modification calls.
The `.` MM flag indicates that primary sequence bases without an associated base modification probability 
should be inferred to be canonical. By default, when this flag is encountered in a modBAM, `modkit extract` will 
output rows with the `inferred` column set to `true` and a `mod_qual` value of `0.0` for the base modifications
called on that read. For example, if you have a `A+a.` MM tag, and there are `A` bases in the read for which 
there aren't base modification calls (identifiable as non-0s in the MM tag) will be rows where the `mod_code` 
is `a` and the `mod_qual` is 0.0.

## Note on non-primary alignments
If a valid `MN` tag is found, secondary and supplementary alignments can be output in the `modkit extract` tables above.
See [troubleshooting](./troubleshooting.md) for details on how to get valid `MN` tags.
To have non-primary alignments appear in the output, the `--allow-non-primary` flag must be passed. 
By default, the primary alignment will have all base modification information contained on the read, including soft-clipped and unaligned read positions. 
If the `--mapped-only` flag is used, soft clipped sections of the read will not be included. 
For secondary and supplementary alignments, soft-clipped positions are not repeated. See [advanced usage](./advanced_usage.md) for more details.

## Example usages:

### Extract a table of base modification probabilities from an aligned and indexed BAM 
```
modkit extract full <input.bam> <output.tsv> [--bgzf]
```
If the index `input.bam.bai` can be found, intervals along the aligned genome can be performed
in parallel. The optional `--bgzf` flag will emit compressed output.

### Extract a table from a region of a large modBAM
The below example will extract reads from only chr20, and include reference sequence context
```
modkit extract full <intput.bam> <output.tsv> --region chr20 --ref <ref.fasta>
```

### Extract only sites aligned to a CG motif
```
modkit motif bed <reference.fasta> CG 0 > CG_motifs.bed
modkit extract full <in.bam> <out.tsv> --ref <ref.fasta> --include-bed CG_motifs.bed
```

### Extract only sites that are at least 50 bases from the ends of the reads
```
modkit extract full <in.bam> <out.tsv> --edge-filter 50
```

### Extract read-level base modification calls

```
modkit extract calls <input.bam> <calls.tsv>
```

Use `--allow-non-primary` to get secondary and supplementary mappings in the output.

```
modkit extract calls <input.bam> <output.tsv> --allow-non-primary
```

See the help string and/or [advanced_usage](./advanced_usage.md) for more details and [performace considerations](./perf_considerations.md) if you encounter issues with memory usage.
