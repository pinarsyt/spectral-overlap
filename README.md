# spectral-overlap

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

`spectral-overlap` compares Gaussian TD-DFT output files with the CIE
LED-B4 reference spectrum.

It reads `.out` files and writes the results as CSV and JSON.

- Manual: [`docs/MANUAL.md`](docs/MANUAL.md)

## Requirements

- Python 3.10+
- `numpy`
- `pandas`

Install:

```bash
pip install numpy pandas
```

## Usage

From the repository root:

```bash
python run.py <input_dir> <output_dir>
```

Example:

```bash
python run.py input output
```

The script scans `<input_dir>` for `.out` files, processes valid
TD-DFT outputs, and writes the results into `<output_dir>`.

## Output

The output directory contains:

- `results.csv`
- `results.json`
- `skipped.json`

## Main Parameters

The main settings are defined at the top of [`run.py`](run.py):

```python
SIGMA_EV         = 0.3
WL_MIN           = 200.0
WL_MAX           = 800.0
NUM_POINTS       = 4000
CONCENTRATION_M  = 1.0e-5
PATH_CM          = 1.0
LED_COLUMN       = 4
```

Main settings:

- `SIGMA_EV` controls broadening
- `CONCENTRATION_M` controls reference concentration
- `PATH_CM` controls path length
- `LED_COLUMN` selects the LED spectrum

## Examples

Included example workflows:

- Single dye: run one TD-DFT output and inspect one result row
- Library screen: run the full 13-dye set and rank by `absorbed_fraction`

Single dye example:

```bash
mkdir -p case_studies/01_single_dye/input
cp input/slurm-5473089.out case_studies/01_single_dye/input/
python run.py case_studies/01_single_dye/input case_studies/01_single_dye/output
```

Library screen example:

```bash
python run.py input case_studies/02_library_screen/output
```

Useful columns in `results.csv`:

- `sample_id`
- `lambda_max_nm`
- `absorbed_fraction`
- `shape_overlap`
- `absorptance_max`

## Limitations

- Works with Gaussian-style TD-DFT `.out` files
- Uses Gaussian broadening
- Uses bundled CIE LED data

## Third-Party Data

The bundled [`CIE_illum_LEDs_1nm.csv`](CIE_illum_LEDs_1nm.csv) file
comes from CIE reference data.

- Website: <https://cie.co.at/>
