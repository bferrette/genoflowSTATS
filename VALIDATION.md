# genoflowSTATS validation

The initial implementation was validated against two local SAMtools 1.24 reports and synthetic assembly FASTAs. The reports are not distributed in this repository.

## SAMtools coverage reconstruction

### `SFA01.stats`

From `COV` plus the `RFS` reference length:

- reference positions: **612,904,761**
- reconstructed zero-depth positions: **3,422,444**
- fraction with depth ≥1×: **0.9944160264** (99.4416%)
- median positive depth: **11×**
- coverage mean from the finite histogram: **11.5937× lower bound**
- mapped-CIGAR coverage: **11.9294×**
- SAMtools NM/edit QV: **21.3130**

The histogram mean is labeled a lower bound because the source contains a non-empty `[1000<]` overflow bin.

### `SFA01.nonrep.stats`

From `COV` plus `SN bases inside the target`:

- target positions: **476,311,082**
- reconstructed zero-depth positions: **1,515,895**
- fraction with depth ≥1×: **0.9968174266** (99.6817%)
- median positive depth: **11×**

SAMtools reports `percentage of target genome with coverage > 0 (%) = 99.68`, matching the reconstructed value after rounding.

## Multi-assembly comparison

A two-assembly synthetic manifest was run through `genoflowSTATS compare` with process-level parallel FASTA evaluation. It successfully produced:

- `comparison.tsv`
- `comparison.long.tsv`
- `comparison.json`
- `comparison.nx_ngx.svg`
- `comparison.heatmap.svg`

## Syntax/import validation

```bash
python -m py_compile src/genoflowstats/__init__.py src/genoflowstats/cli.py
PYTHONPATH=src python -m genoflowstats.cli --version
```

Expected:

```text
genoflowSTATS 0.2.0
```

## Version 0.2.0 checks

Focused unit tests cover cumulative GCD weights, coverage denominator reconstruction and rejection, missing breadth, missing comparison values, and unreachable NGx values. The local `SFA01.stats` was reprocessed and retains 3,422,444 reconstructed zero-depth sites and 0.9944160264 breadth. The Snook `samtools.stats.txt` was read remotely without modification and processed locally using the corresponding reported assembly length, 693,292,604 bp: 460,818 zero-depth positions and 0.99933532 breadth. Its corrected GCD weighted mean of bin medians is 11.79419. A synthetic two-assembly comparison and FASTA run generated SVGs; the comparison was run with one worker because process creation was denied in the sandbox. These checks do not constitute a remote production rerun.

## QUAST and AWK comparison

QUAST 5.3.0 at `/home/bferrette/quast` was run without a reference, using `--min-contig 0`. A synthetic assembly of records 1,000, 600, and 300 bp had total length 1,900 bp, N50 1,000 bp, N90 300 bp, L50 1, L90 3, and auN 763.1579 in genoflowSTATS; QUAST agreed, reporting auN rounded to 763.2. The AWK count was 1,900. QUAST's bundled `contigs_1.fasta` had length 6,710, N50 3,980, N90 1,120, and auN 2,933.9642; QUAST agreed to its displayed precision and AWK counted 6,710. With QUAST's default 500 bp minimum, the synthetic main total was 1,600 while its >=0 bp total remained 1,900. With `--est-ref-size 5000`, QUAST omitted unreachable NG50 and genoflowSTATS now reports NA rather than zero; both calculate auNG 290.0. These are controlled local comparisons, not a production assembly rerun.

## Automatic FASTA reference length

Focused tests cover sequence-line counting with plain and gzip FASTA, automatic resolution from a whole-reference SAMtools command, reuse of the manifest assembly length in comparison mode, explicit FASTA input for reports lacking a command header, and preservation of target-specific denominators. The Snook report command was parsed and its recorded reference resolved to the correct GVFS path. The complete remote FASTA scan was interrupted before completion, so end-to-end automatic counting of the 693 Mb remote file remains unverified; the previously recorded 693,292,604 bp result used an explicitly supplied denominator.
