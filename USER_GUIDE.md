# genoflowSTATS

`genoflowSTATS` is the standalone genome-assembly statistics and QC engine developed for the genoflow framework. It replaces the uses of QUAST and Qualimap in genoflowHiFi with deterministic FASTA metrics, direct parsing of `samtools stats`, robust coverage summaries, optional BAM interrogation through pysam, and comprehensive multi-assembly comparison.

A central goal of `genoflowSTATS` is to evaluate **multiple genome assemblies from different assembly, polishing, correction, scaffolding, and gap-closing stages in a single analysis**. This makes it possible to quantify how each processing step changes contiguity, assembly size, gap content, read support, coverage, and other QC metrics.

## Design

The main evidence layers are:

- **Assembly sequence:** scaffold/contig counts, N50/N90, NG50/NG90, auN/auNG, gap statistics, GC, N content, assembly/genome-size agreement.
- **Alignment summary:** direct parsing of `samtools stats` `SN` records.
- **Coverage:** reconstruction from `COV` plus `RFS`/target length, including zero-depth sites, median/mode/MAD/quantiles and coverage-threshold fractions.
- **Mapping confidence:** exact `MAPQ` histogram summaries.
- **Indels:** `ID` event and base-length distributions.
- **GC-depth behavior:** `GCD` quantiles and robust GC-depth heterogeneity metrics.
- **Reference integrity cross-check:** `RFS` summary and per-reference rows.
- **Direct BAM QC:** optional one-pass pysam mode for metrics not present in saved `samtools stats`.
- **Comparison:** one command evaluates multiple assemblies and creates wide/long TSV reports and comparative SVG plots.

Plots are exported as SVG. Use `--palette colorblind` (default), `--palette ocean`, or `--palette heat` with `fasta`, `stats`, and `compare`. Comparison curves use distinct categorical colors; the comparison heatmap uses a diverging scale.

## Recommended multi-assembly workflow

For genome-assembly development, the recommended workflow is to retain the representative FASTA produced after each major stage.

For example:

```text
01_hifiasm_primary.fasta
02_fcs_cleaned.fasta
03_prepolished.fasta
04_inspector_corrected.fasta
05_scaffolded.fasta
06_final_refined.fasta
```

The same HiFi reads can then be mapped independently against each assembly to generate corresponding `samtools stats` reports:

```text
01_hifiasm_primary.stats
02_fcs_cleaned.stats
03_prepolished.stats
04_inspector_corrected.stats
05_scaffolded.stats
06_final_refined.stats
```

These FASTA and statistics files can be provided together through a tab-delimited manifest and analyzed simultaneously with:

```bash
genoflowSTATS compare \
  --manifest assemblies.tsv \
  --genome-size 612904761 \
  --threads 8 \
  --outdir comparison
```

This longitudinal approach makes it possible to measure improvements or regressions introduced by each assembly-processing stage.

## Why `samtools stats COV` is preferred for global coverage

When `samtools stats` is run with unit coverage bins, for example:

```bash
samtools stats \
  --coverage 1,1000,1 \
  --ref-seq assembly.fasta \
  assembly.bam > assembly.stats
```

`COV` is already a compact depth histogram. `genoflowSTATS` combines it with the target/reference length from `SN` or `RFS` to reconstruct zero-depth sequence.

This avoids generating and parsing genome-scale `samtools depth -aa` text when only global coverage statistics are needed.

The overflow bin (`[1000<]` in this example) is retained. If it is non-zero, the reported histogram mean is explicitly labeled a lower bound because actual depths in that bin exceed the bin boundary.

## Commands

### 1. FASTA metrics for an individual assembly

The `fasta` command evaluates one FASTA at a time and is useful for inspecting an individual assembly:

```bash
genoflowSTATS fasta \
  --fasta 01_hifiasm_primary.fasta \
  --genome-size 612904761 \
  --outdir assess_01_hifiasm
```

For another assembly:

```bash
genoflowSTATS fasta \
  --fasta 04_inspector_corrected.fasta \
  --genome-size 612904761 \
  --outdir assess_04_inspector
```

And for the final assembly:

```bash
genoflowSTATS fasta \
  --fasta 06_final_refined.fasta \
  --genome-size 612904761 \
  --outdir assess_06_final
```

Principal outputs for each assembly are:

```text
<assembly>.assembly_metrics.tsv
<assembly>.assembly_metrics.json
<assembly>.scaffolds.tsv
<assembly>.scaffold_cumulative.svg
<assembly>.nx_ngx.svg
```

For simultaneous comparison of multiple FASTA files, use the `compare` command described below.

### 2. Generate `samtools stats` for multiple assemblies

Each assembly should be evaluated against the same HiFi read dataset whenever read-support comparisons are desired.

For example:

```bash
minimap2 -ax map-hifi \
  --secondary=no \
  --cs=long \
  01_hifiasm_primary.fasta \
  hifi_reads.fastq.gz |
samtools sort \
  -@ 16 \
  -o 01_hifiasm_primary.bam
```

```bash
samtools index \
  -@ 16 \
  01_hifiasm_primary.bam
```

```bash
samtools stats \
  --coverage 1,1000,1 \
  --ref-seq 01_hifiasm_primary.fasta \
  01_hifiasm_primary.bam \
  > 01_hifiasm_primary.stats
```

Repeat the same procedure for the other assembly stages:

```text
02_fcs_cleaned.fasta
03_prepolished.fasta
04_inspector_corrected.fasta
05_scaffolded.fasta
06_final_refined.fasta
```

producing:

```text
02_fcs_cleaned.stats
03_prepolished.stats
04_inspector_corrected.stats
05_scaffolded.stats
06_final_refined.stats
```

Using the same read dataset and the same minimap2/SAMtools parameters across all assemblies is important for direct comparison.

### 3. Parse individual `samtools stats` files

Each statistics file can also be examined independently:

```bash
genoflowSTATS stats \
  --stats 01_hifiasm_primary.stats \
  --outdir assess_stats/01_hifiasm
```

```bash
genoflowSTATS stats \
  --stats 04_inspector_corrected.stats \
  --outdir assess_stats/04_inspector
```

```bash
genoflowSTATS stats \
  --stats 06_final_refined.stats \
  --reference-fasta 06_final_refined.fasta \
  --outdir assess_stats/06_final
```

If total input HiFi bases are known, provide the same value for every assembly generated from that dataset:

```bash
genoflowSTATS stats \
  --stats 06_final_refined.stats \
  --hifi-bases 20000000000 \
  --outdir assess_stats/06_final
```

This allows nominal input coverage to be compared with mapped-CIGAR coverage and empirical median coverage.

Detailed `RFS` rows are not retained by default because target-restricted SAMtools reports can contain hundreds of thousands of interval rows.

Write them explicitly when needed:

```bash
genoflowSTATS stats \
  --stats 06_final_refined.stats \
  --outdir assess_stats/06_final \
  --write-reference-table
```

Principal outputs include:

```text
<assembly>.samtools_metrics.tsv
<assembly>.samtools_metrics.json
<assembly>.coverage_histogram.tsv
<assembly>.mapq_histogram.tsv
<assembly>.indel_distribution.tsv
<assembly>.gc_depth.tsv
<assembly>.reference_sequences.tsv
<assembly>.coverage_histogram.svg
<assembly>.coverage_cumulative.svg
<assembly>.mapq.svg
<assembly>.indels.svg
<assembly>.gc_depth.svg
```

`<assembly>.reference_sequences.tsv` is generated only when `--write-reference-table` is enabled.

### 4. Direct BAM QC

Direct BAM analysis remains available for individual assemblies:

```bash
genoflowSTATS bam \
  --bam 01_hifiasm_primary.bam \
  --outdir assess_bam/01_hifiasm
```

```bash
genoflowSTATS bam \
  --bam 06_final_refined.bam \
  --outdir assess_bam/06_final
```

This is optional for global statistics when detailed `samtools stats` files already exist. It remains useful for metrics that are unavailable from the saved summary and for future spatial/junction analyses.

## 5. Compare multiple assemblies

Multi-assembly comparison is the recommended mode for genoflowSTATS.

Create a **tab-delimited** manifest such as:

```text
name	fasta	samtools_stats	hifi_bases
01_hifiasm_primary	/path/01_hifiasm_primary.fasta	/path/01_hifiasm_primary.stats	20000000000
02_fcs_cleaned	/path/02_fcs_cleaned.fasta	/path/02_fcs_cleaned.stats	20000000000
03_prepolished	/path/03_prepolished.fasta	/path/03_prepolished.stats	20000000000
04_inspector_corrected	/path/04_inspector_corrected.fasta	/path/04_inspector_corrected.stats	20000000000
05_scaffolded	/path/05_scaffolded.fasta	/path/05_scaffolded.stats	20000000000
06_final_refined	/path/06_final_refined.fasta	/path/06_final_refined.stats	20000000000
```

The columns are:

```text
name
fasta
samtools_stats
hifi_bases
```

Only `name` and `fasta` are mandatory.

Therefore, FASTA-only comparison is also possible:

```text
name	fasta
01_hifiasm_primary	/path/01_hifiasm_primary.fasta
02_fcs_cleaned	/path/02_fcs_cleaned.fasta
03_prepolished	/path/03_prepolished.fasta
04_inspector_corrected	/path/04_inspector_corrected.fasta
05_scaffolded	/path/05_scaffolded.fasta
06_final_refined	/path/06_final_refined.fasta
```

Run the comparison:

```bash
genoflowSTATS compare \
  --manifest assemblies.tsv \
  --genome-size 612904761 \
  --threads 8 \
  --outdir comparison
```

Principal outputs are:

```text
comparison.tsv
comparison.long.tsv
comparison.json
comparison.nx_ngx.svg
comparison.heatmap.svg
```

### Wide comparison table

`comparison.tsv` is designed for direct inspection and follows a QUAST-like structure:

```text
metric	unit	direction	01_hifiasm_primary	02_fcs_cleaned	03_prepolished	04_inspector_corrected	05_scaffolded	06_final_refined
```

Metrics are rows and assemblies are columns.

The additional `direction` field describes the expected interpretation of each metric, for example:

```text
higher
lower
target_1
descriptive
```

This permits downstream ranking or stepwise-improvement calculations without embedding metric interpretation in the Bash pipeline.

### Long comparison table

`comparison.long.tsv` uses a tidy format:

```text
assembly	metric	value	unit	direction
```

This representation is suitable for:

- Polars;
- pandas;
- R;
- custom plotting;
- longitudinal stage analysis;
- automated candidate ranking.

## Longitudinal assembly assessment

When assemblies correspond to sequential processing stages, the comparison should be interpreted as a trajectory rather than as independent samples.

For example:

```text
Hifiasm
   ↓
FCS cleaning
   ↓
DeepVariant polishing
   ↓
Inspector correction
   ↓
scaffolding
   ↓
gap filling / final refinement
```

Useful changes to monitor include:

### Assembly size

```text
total_length
assembly_genome_size_ratio
```

Large unexpected changes can indicate sequence removal, duplication, contamination removal, or changes in haplotig representation.

### Contiguity

```text
N50
N90
NG50
NG90
auN
auNG
```

Scaffolding is expected primarily to improve this metric group.

### Fragmentation

```text
scaffold_count
contig_count
```

### Gap content

```text
gap_count
gap_bp
gap_fraction
median_gap
p95_gap
largest_gap
```

This is particularly important when evaluating scaffolding and gap-closing stages.

### Read support

When corresponding `samtools stats` files are available:

```text
mapped_fraction
MAPQ distribution
mapped_CIGAR_coverage
median empirical depth
zero_coverage_fraction
SAMtools edit/error rate
indel burden
GC-depth behavior
```

A structurally more contiguous assembly should ideally improve contiguity without causing deterioration in read support.

## Coverage quantities

`genoflowSTATS` intentionally separates different coverage definitions.

### Nominal input coverage

\[
C_{\mathrm{nominal}}
=
\frac{\text{input HiFi bases}}
{\text{target or assembly length}}
\]

### Mapped-CIGAR coverage

\[
C_{\mathrm{mapped}}
=
\frac{\text{samtools bases mapped (cigar)}}
{\text{target or assembly length}}
\]

### Empirical coverage distribution

Reconstructed from `COV`, including zero-depth positions whenever a reference/target denominator is available.

### Empirical median positive depth

Median depth among positions with coverage greater than zero.

### Depth-implied sequence space

\[
G_{\mathrm{depth}}
=
\frac{\text{mapped CIGAR bases}}
{\text{median positive depth}}
\]

This is a consistency check and is not intended to replace GenomeScope2.

## Interpretation of the SAMtools error metric

The `SN error rate` is retained exactly as reported by SAMtools.

`genoflowSTATS` converts it to:

```text
samtools_nm_edit_qv
```

for convenient interpretation:

\[
QV=-10\log_{10}(e)
\]

where \(e\) is the SAMtools error rate.

Because SAM `NM` represents edit distance and can include indels, genoflowSTATS reports insertion and deletion burdens separately and does **not** add those indels again to `SN mismatches`. This avoids double counting.

## Example complete multi-assembly workflow

Given:

```text
01_hifiasm_primary.fasta
02_fcs_cleaned.fasta
03_prepolished.fasta
04_inspector_corrected.fasta
05_scaffolded.fasta
06_final_refined.fasta
```

generate an alignment and statistics file for each FASTA using the same HiFi reads.

For example:

```bash
for fasta in \
  01_hifiasm_primary.fasta \
  02_fcs_cleaned.fasta \
  03_prepolished.fasta \
  04_inspector_corrected.fasta \
  05_scaffolded.fasta \
  06_final_refined.fasta
do
    sample="${fasta%.fasta}"

    minimap2 -ax map-hifi \
      --secondary=no \
      --cs=long \
      "${fasta}" \
      hifi_reads.fastq.gz |
    samtools sort \
      -@ 16 \
      -o "${sample}.bam"

    samtools index \
      -@ 16 \
      "${sample}.bam"

    samtools stats \
      --coverage 1,1000,1 \
      --ref-seq "${fasta}" \
      "${sample}.bam" \
      > "${sample}.stats"
done
```

Create:

```text
assemblies.tsv
```

containing:

```text
name	fasta	samtools_stats	hifi_bases
01_hifiasm_primary	01_hifiasm_primary.fasta	01_hifiasm_primary.stats	20000000000
02_fcs_cleaned	02_fcs_cleaned.fasta	02_fcs_cleaned.stats	20000000000
03_prepolished	03_prepolished.fasta	03_prepolished.stats	20000000000
04_inspector_corrected	04_inspector_corrected.fasta	04_inspector_corrected.stats	20000000000
05_scaffolded	05_scaffolded.fasta	05_scaffolded.stats	20000000000
06_final_refined	06_final_refined.fasta	06_final_refined.stats	20000000000
```

Then run:

```bash
genoflowSTATS compare \
  --manifest assemblies.tsv \
  --genome-size 612904761 \
  --threads 8 \
  --outdir genoflowSTATS_comparison
```

This produces one integrated comparison of the complete assembly-development trajectory.

## Current scope

The package implements:

- deterministic FASTA assembly metrics;
- multiple-assembly comparison;
- direct `samtools stats` parsing;
- global coverage reconstruction;
- MAPQ distributions;
- indel statistics;
- GC-depth summaries;
- optional BAM QC;
- comparative TSV/JSON outputs;
- comparative SVG visualizations.

Planned extensions can include:

- spatial genomic-window analysis;
- scaffold-junction HiFi support;
- multivariate anomaly detection;
- GenomeScope2 ingestion;
- Merqury ingestion;
- automated Pareto filtering;
- direct genoflowHiFi Step-5 candidate ranking;
- explicit per-stage improvement tables.

These extensions can use the existing TSV/JSON schemas without changing the definitions of the current core metrics.

### Coverage domain and GC-depth interpretation (0.2.0)

For a recognized whole-reference `samtools stats` command, genoflowSTATS automatically locates the recorded `--ref-seq` FASTA and sums its sequence-line lengths, including N bases. The stats file and reference must be accessible; GVFS SFTP paths are resolved on the same remote mount. Supply `stats --reference-fasta assembly.fasta` when the recorded reference path is unavailable or the report has no command header. `compare` reuses the length already calculated for the manifest FASTA when the saved command identifies a whole-reference run. Each FASTA must match its stats report; length alone does not establish sequence identity.

SAMtools target or RFS denominators take precedence. A restricted or unrecognized recorded command does not establish that a full FASTA is the histogram domain, so no FASTA-derived denominator is used in that case. If a restricted report lacks its target denominator, provide the matching target size with `--reference-length`; without a valid denominator, breadth remains unknown and fractions are labeled `_covered_sites`. A conflicting or too-small denominator is rejected. COV bins must have unit width (for example, `--coverage 1,1000,1`); the overflow bin makes the mean a lower bound. Thresholds beyond its lower bound are unknown.

The GCD unique-sequence column is cumulative. The weighted mean of per-bin medians uses differences between consecutive cumulative percentages; it is a summary of bin medians, not a global median. The plot retains the SAMtools GCD depth scale.

Comparison heatmaps show relative, direction-aware scores, not absolute quality or statistical significance. Assembly length is descriptive; only the assembly-to-genome ratio is scored against 1. Missing measurements are blank. An NGx curve ends when an assembly does not reach the requested percentage of estimated genome size.


### FASTA length and QUAST comparisons

`total_length` sums the lengths of every FASTA record, including N and other ambiguity symbols. It is an assembly length, not an independent genome-size estimate. Empty records, sequence before a header, and internal sequence whitespace are rejected. Compressed `.gz` FASTA files are supported.

For a plain, consistently formatted FASTA with Unix line endings and no whitespace in sequence lines, the same count is available with:

```bash
LC_ALL=C awk '!/^>/ {sum += length($0)} END {printf "%.0f\n", sum}' input.fasta
```

The `printf` avoids scientific notation on typical awk implementations. Windows CRLF line endings and whitespace can inflate the simple AWK count. AWK counts every non-header character and does not validate FASTA records.

QUAST defaults to a 500 bp minimum sequence length for its main report. To compare all records with genoflowSTATS, run QUAST with `--min-contig 0` and compare its `Total length (>= 0 bp)` or its main `Total length`. QUAST calls FASTA records contigs even when the records are scaffolds; genoflowSTATS separately reports scaffold metrics and gap-split contig metrics. A supplied genome-size estimate controls NGx and auNG, not total assembly length. When the assembly is too short to reach an NGx threshold, genoflowSTATS reports `NA`.
