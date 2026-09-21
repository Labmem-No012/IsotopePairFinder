# HDPairFinder
[![Generic badge](https://img.shields.io/badge/HDPairFinder-ver_1.0-<COLOR>.svg)](https://github.com/HuanLab/HDPairFinder)
![Maintainer](https://img.shields.io/badge/maintainer-Tingting_Zhao,_Tao_Huan-blue)

`HDPairFinder` is a bioinformatics software to effectively recognize hydrogen and deuterium labeled compounds in untargeted LC-MS  analysis.

## Part 1: Quick start
### 1.1 Installation
The source code of `HDPairFinder` can be freely downloaded on this [GitHub release page](https://github.com/HuanLab/HDPairFinder/releases/tag/v1.0).

`Demo data` can be freely downloaded on this [Demo data page](https://pan.baidu.com/s/1wuEh4VMJyHmsTQxwvbsRNQ?pwd=nxyy)

To run `HDPairFinder`, R version 4.2.0 or above is required and we recommend using RStudio.
### 1.2 Task

- [Extraction of H/D-labeled chemical features](https://github.com/HuanLab/HDPairFinder)
- [Alignment](https://github.com/HuanLab/HDPairFinder)
- [Evidence-based missing value imputation](https://github.com/HuanLab/HDPairFinder)
- [Putative compound annotation](https://github.com/HuanLab/HDPairFinder)

### 1.3 File import
- Single import (a single mzML or mzXML file)
- Batch import (mzML or mzXML files)
### 1.4 Result export
- Single export (a single CSV pair table) 
- Batch export (aligned CSV pair table)

## Part 2: User manual
Detailed instructions on `HDPairFinder` can be found in [HDPairFinder user manual](https://github.com/HuanLab/HDPairFinder).

## Part 3: Citation

## Parallel launcher

`run_parallel.R` processes independent mzML/mass-accuracy CSV pairs in parallel
without adding parallel code to the HDPairFinder workflow itself. The launcher
only examines files directly inside the input directory (not its subfolders)
and matches files by name:

```text
sample_A.mzML
sample_A_mass_accuracy.csv
```

Install the launcher dependency once:

```r
install.packages("future")
```

Run with the default input directory (`mzML/`). The launcher uses up to 90% of
the available CPU cores when enough input pairs exist:

```sh
Rscript run_parallel.R
```

Or provide the input directory, worker count, and output directory:

```sh
Rscript run_parallel.R /path/to/input 2 /path/to/output
```

Each pair runs in an isolated directory under `.hdpairfinder_runs`, with its own
`HDPairFinder.log`. Files without a matching mass-accuracy CSV are skipped. The
sixth column of the CSV's first data row (spreadsheet cell F2) supplies that
sample's pair-picking and alignment m/z tolerances in ppm.

ISFrag runs only when the mzML file contains MS2 spectra. For MS1-only files,
HDPairFinder records that ISFrag was skipped and continues pair picking from the
raw feature table.

The launcher parallelizes independent per-sample processing. Cross-sample
alignment should be run once after all sample jobs have completed.

### Resource limits

The launcher defaults to a maximum of four concurrent jobs as well as a 90%
CPU-capacity ceiling. The four-job cap is based on an observed peak of about
4.1 GiB per job and is intended to keep combined job RAM below 17 GiB. Every
sample is limited to one BLAS/OpenMP thread, and Unix jobs run at niceness 10.
A command-line worker count above either the four-job cap or CPU budget is
reduced automatically. RAM is monitored; because memory use varies with input,
the job cap is a conservative scheduling limit rather than an OS-enforced RAM
boundary.

The limits can be changed explicitly:

```sh
HDPAIRFINDER_CPU_FRACTION=0.90 \
HDPAIRFINDER_MAX_WORKERS=4 \
HDPAIRFINDER_RESERVED_CORES=0 \
HDPAIRFINDER_THREADS_PER_JOB=1 \
HDPAIRFINDER_NICE=10 \
HDPAIRFINDER_TELEMETRY_SECONDS=5 \
Rscript run_parallel.R /path/to/input 4 /path/to/output
```

### Telemetry

Each run directory contains:

- `telemetry.csv`: timestamped whole-system CPU utilization, load averages,
  active jobs, RAM, and swap usage.
- `telemetry_summary.csv`: peak CPU, load, and RAM values for the run.
- `settings.csv`: detected capacity and the limits applied to the run.
- `summary.csv`: per-job elapsed time, CPU time, average CPU utilization, and
  peak resident memory.
- `<sample>/resources.txt`: the complete GNU Time report for that sample.

Whole-system telemetry includes activity from programs other than
HDPairFinder. Per-job measurements describe the individual HDPairFinder
subprocesses.
