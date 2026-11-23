# genoflowSTATS — Mamba environment and installation

## Requirements

- Linux x86_64
- Mamba/Miniforge
- Python 3.14
- For producing input alignment statistics: external `samtools` (genoflowHiFi currently targets SAMtools 1.24) and `minimap2`.

The Python environment intentionally contains the statistical/QC layer only. `samtools`, `minimap2`, GenomeScope2, Merqury, BUSCO, and BlobToolKit can remain independently versioned external tools/environments.

## Option 1 — create the environment directly

```bash
mamba create -n GENOFLOWSTATSenv \
  -c conda-forge -c bioconda \
  --override-channels \
  --strict-channel-priority \
  -y \
  python=3.14 pip "setuptools>=77" wheel \
  numpy scipy pysam polars pyarrow \
  statsmodels scikit-learn matplotlib \
  rich psutil jinja2 pydantic

mamba activate GENOFLOWSTATSenv
```

## Option 2 — create from `environment.yml`

From the project directory:

```bash
mamba env create -f environment.yml -y
mamba activate GENOFLOWSTATSenv
```

The channels are already declared in `environment.yml`. To verify strict channel priority for the active environment:

```bash
mamba config --env --set channel_priority strict
```

## Install genoflowSTATS without resolving dependencies again

All Python dependencies are already installed by Mamba. Install only the local project:

```bash
PIP_REQUIRE_VIRTUALENV=false \
python -m pip install --no-deps --no-build-isolation .
```

Verify:

```bash
genoflowSTATS --version
```

Expected for this release:

```text
genoflowSTATS 0.2.0
```

## Environment integrity test

```bash
python - <<'PY'
import sys
import numpy, scipy, pysam, polars, pyarrow
import statsmodels, sklearn, matplotlib, rich, psutil, jinja2, pydantic

print("Python", sys.version.split()[0])
for mod in (
    numpy, scipy, pysam, polars, pyarrow,
    statsmodels, sklearn, matplotlib, rich,
    psutil, jinja2, pydantic,
):
    print(f"{mod.__name__}: {getattr(mod, '__version__', 'OK')}")
PY

genoflowSTATS --version
```

## Upgrade after editing the source

```bash
mamba activate GENOFLOWSTATSenv
PIP_REQUIRE_VIRTUALENV=false \
python -m pip install --no-deps --no-build-isolation --upgrade .
```

## Remove the environment

```bash
mamba deactivate
mamba env remove -n GENOFLOWSTATSenv -y
```
