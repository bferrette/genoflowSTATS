<div align="center">
  <img src="genoflowSTATS_badge.png" alt="genoflowSTATS hexagonal badge" width="340">

# genoflowSTATS

**Genome assembly statistics, read-support quality control, and comparative visualisation**

Part of the [genoflowVerse](https://github.com/bferrette) ecosystem · Python ≥ 3.14 · Version 0.2.0

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24483249.svg)](https://doi.org/10.5281/zenodo.24483249)

</div>

---

genoflowSTATS is a command-line toolkit for assessing genome assemblies and comparing processing stages. It combines FASTA-based contiguity and gap metrics with summaries from `samtools stats`, optional BAM-level alignment statistics, and SVG visualisations. It can be used independently or alongside the [genoflowHiFi](https://github.com/bferrette/genoflowHiFi) assembly workflow.

| Input | What it describes | Command |
|---|---|---|
| Assembly FASTA | Length, sequence composition, contiguity and gaps | `fasta` |
| `samtools stats` report | Coverage, mapping quality, indels, GC–depth and read lengths | `stats` |
| BAM alignment | One-pass alignment and CIGAR summaries | `bam` |
| Assembly manifest | Comparative metrics and plots across assemblies | `compare` |

```text
FASTA assemblies ───────────────┐
                                ├──> genoflowSTATS ──> TSV + JSON + SVG
SAMtools stats / BAM ───────────┘                       individual or comparative
```

## Highlights

- **Assembly quality:** scaffold and gap-split contig N50/N90, L50/L90, auN, gap statistics, GC content, and N content. An optional estimated genome size enables NGx, auNG, and assembly-to-genome comparisons.
- **Read support:** depth distributions and threshold fractions, MAPQ distributions, read-length summaries, and indel-size summaries from existing SAMtools reports.
- **Comparisons:** a tab-delimited manifest brings assemblies and optional read-support reports together in wide and long metric tables, Nx/NGx curves, and a relative metric heatmap.
- **Plots:** SVG output with `colorblind`, `ocean`, and `heat` palettes. `colorblind` is the default.
- **Flexible execution:** evaluate one assembly, inspect one alignment report, or compare assembly stages without remapping reads inside genoflowSTATS.

## Installation

The package requires Python 3.14 or later. Create the supplied environment and install the local source:

```bash
git clone https://github.com/bferrette/genoflowSTATS.git
cd genoflowSTATS
mamba env create -f environment.yml
mamba activate GENOFLOWSTATSenv
python -m pip install --no-deps --no-build-isolation .
genoflowSTATS --version
```

The environment contains the Python dependencies. Tools used to generate inputs, such as SAMtools and an aligner, are installed separately. See [INSTALL.md](INSTALL.md) for environment setup and upgrade instructions.

## Quick start

### 1. Assess an assembly

```bash
genoflowSTATS fasta \
  --fasta assembly.fasta \
  --genome-size 700000000 \
  --outdir results/assembly \
  --prefix assembly \
  --palette ocean
```

`--genome-size` is optional and must be an independent estimate in base pairs. Total FASTA length is measured from the assembly; it is not an independent genome-size estimate. By default, runs of at least 10 Ns split scaffold records into contigs for contig statistics; change this with `--gap-min`.

### 2. Summarise read support

Generate a coordinate-sorted BAM and a SAMtools report against the same assembly. The following example uses unit-width coverage bins:

```bash
samtools stats --coverage 1,1000,1 \
  --ref-seq assembly.fasta aligned.sorted.bam > assembly.stats

genoflowSTATS stats \
  --stats assembly.stats \
  --reference-fasta assembly.fasta \
  --outdir results/read_support \
  --prefix assembly
```

The matching reference FASTA allows the package to calculate its length from sequence lines when the report lacks a usable reference denominator. It can also resolve the reference recorded in a recognised whole-reference SAMtools command. Target-restricted reports require the matching target denominator, not the full FASTA length. Use `--reference-length` for a known target length when the report does not provide it. If the applicable denominator cannot be established, zero-depth sites and reference-wide breadth remain unknown, and fractions are labelled as conditional on covered sites.

SAMtools COV input must use unit-width bins. An overflow bin limits the precision of the mean depth; genoflowSTATS reports it as a lower bound when overflow observations are present. Read the [coverage interpretation](USER_GUIDE.md#coverage-domain-and-gc-depth-interpretation-020) before comparing reports generated with different settings or target domains.

### 3. Compare assemblies

Create a tab-separated file with unique assembly names and paths. `samtools_stats` and `hifi_bases` are optional; the supplied FASTA and stats for each row must describe the same assembly.

```text
name	fasta	samtools_stats
primary	/path/to/primary.fasta	/path/to/primary.stats
scaffolded	/path/to/scaffolded.fasta	/path/to/scaffolded.stats
final	/path/to/final.fasta	/path/to/final.stats
```

```bash
genoflowSTATS compare \
  --manifest assemblies.tsv \
  --genome-size 700000000 \
  --threads 4 \
  --outdir results/comparison \
  --palette colorblind
```

Omit `--genome-size` if no independent estimate is available. The comparison heatmap uses relative, direction-aware scores; it is not an absolute quality score or a statistical test. A whole-reference FASTA length is reused only when the report identifies a recognised whole-reference SAMtools command.

### 4. Inspect a BAM directly

```bash
genoflowSTATS bam \
  --bam aligned.sorted.bam \
  --outdir results/bam \
  --prefix assembly
```

BAM mode requires `pysam`. See the [user guide](USER_GUIDE.md) for command options, output fields, and additional examples.

## Outputs

| Analysis | Main outputs |
|---|---|
| FASTA | `*.assembly_metrics.tsv`, `*.assembly_metrics.json`, `*.scaffolds.tsv`, cumulative scaffold and Nx/NGx SVGs |
| SAMtools stats | `*.samtools_metrics.tsv`, `*.samtools_metrics.json`, coverage, MAPQ, indel and GC–depth tables and SVGs when available |
| BAM | `*.bam_metrics.tsv` and `*.bam_metrics.json` |
| Comparison | `comparison.tsv`, `comparison.long.tsv`, `comparison.json`, Nx/NGx and heatmap SVGs |

Use `--no-plots` with `fasta`, `stats`, or `compare` to write tables without SVGs. Use `genoflowSTATS <command> --help` for the complete options.

## Interpretation

- Coverage is specific to the reads, alignment filters, target domain, and histogram settings used to produce the SAMtools report. Zero-depth reconstruction is meaningful only with the matching reference or target length.
- The `samtools_nm_edit_qv` is derived from the SAMtools alignment error rate. It describes read-to-reference edit burden and is not an independent assembly consensus QV.
- GC–depth summaries describe the SAMtools GCD bins. The weighted mean of bin medians is not a global coverage median.
- `flagstat.txt` is not parsed by this version. It can be inspected separately alongside the generated summaries.
- QUAST defaults to filtering short sequences in its main report. Set its minimum contig length to zero when comparing all FASTA records with genoflowSTATS.

## Validation and documentation

The project includes focused tests:

```bash
python -m unittest discover -s tests -v
```

[VALIDATION.md](VALIDATION.md) describes synthetic and QUAST comparisons, coverage reconstruction checks, and their limits. Additional workflows, metric definitions, and complete output details are in [USER_GUIDE.md](USER_GUIDE.md).

## Citation

If you use genoflowSTATS v0.2.0, cite:

bferrette. (2026). *bferrette/genoflowSTATS: genoflowSTATS v0.2.0* (Version v0.2.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.24483249

Cite [genoflowHiFi](https://github.com/bferrette/genoflowHiFi) separately when using its workflow.

## License

genoflowSTATS is distributed under the [GNU General Public License v3.0](LICENSE).
