import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
LED_CSV = SCRIPT_DIR / "CIE_illum_LEDs_1nm.csv"
LED_COLUMN = 4

SIGMA_EV = 0.3
WL_MIN = 200.0
WL_MAX = 800.0
NUM_POINTS = 4000
CONCENTRATION_M = 1.0e-5
PATH_CM = 1.0
EPSILON_PREFACTOR = 2.315e8
CM1_PER_EV = 8065.544

_ROUTE_TD = re.compile(r"(?i)#\s*.*\btd\b")
_STATE_RE = re.compile(
    r"Excited State\s+(\d+)\s*:[^\n]*?([\d.]+)\s*eV\s+([\d.]+)\s*nm\s+f=([\d.]+)",
    re.IGNORECASE,
)
_CHK_RE = re.compile(r"(?i)^\s*%chk\s*=\s*(\S+)")


def extract_chk_name(path):
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                m = _CHK_RE.match(line)
                if m:
                    name = Path(m.group(1)).name
                    if name.lower().endswith(".chk"):
                        name = name[:-4]
                    return name or None
    except OSError:
        return None
    return None


def is_tddft_output(path):
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if _ROUTE_TD.search(line) or "Excited State" in line:
                    return True
    except OSError:
        return False
    return False


def parse_td_output(path):
    states = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            m = _STATE_RE.search(line)
            if m:
                states.append(
                    {
                        "state": int(m.group(1)),
                        "energy_ev": float(m.group(2)),
                        "wavelength_nm": float(m.group(3)),
                        "f": float(m.group(4)),
                    }
                )
    return states


def build_epsilon_gaussian(states):
    wavelengths = np.linspace(WL_MIN, WL_MAX, NUM_POINTS)
    nu_grid = 1.0e7 / wavelengths
    sigma_cm1 = SIGMA_EV * CM1_PER_EV
    epsilon = np.zeros_like(wavelengths)
    for s in states:
        if s["wavelength_nm"] <= 0:
            continue
        nu_center = 1.0e7 / s["wavelength_nm"]
        profile = np.exp(-0.5 * ((nu_grid - nu_center) / sigma_cm1) ** 2) / (
            sigma_cm1 * np.sqrt(2.0 * np.pi)
        )
        epsilon += EPSILON_PREFACTOR * s["f"] * profile
    return wavelengths, epsilon


def load_ledb4_on_grid(wavelengths):
    df = pd.read_csv(LED_CSV, header=None)
    wl = pd.to_numeric(df.iloc[:, 0], errors="coerce").to_numpy(dtype=float)
    intensity = pd.to_numeric(df.iloc[:, LED_COLUMN], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(wl) & np.isfinite(intensity)
    wl = wl[valid]
    intensity = intensity[valid]
    order = np.argsort(wl)
    return np.interp(wavelengths, wl[order], intensity[order], left=0.0, right=0.0)


def analyze_states(sample_id, states):
    wavelengths, epsilon = build_epsilon_gaussian(states)
    absorbance = epsilon * CONCENTRATION_M * PATH_CM
    absorptance = 1.0 - np.power(10.0, -absorbance)
    intensity = load_ledb4_on_grid(wavelengths)

    absorbed_flux = float(np.trapezoid(absorptance * intensity, wavelengths))
    light_flux_total = float(np.trapezoid(intensity, wavelengths))
    absorbed_fraction = absorbed_flux / light_flux_total if light_flux_total > 0.0 else 0.0

    absorptance_max = float(np.max(absorptance))
    epsilon_max = float(np.max(epsilon))
    lambda_max_nm = float(wavelengths[int(np.argmax(epsilon))])

    a_norm = absorptance / max(absorptance_max, 1e-12)
    i_peak = max(float(np.max(intensity)), 1e-12)
    i_norm = intensity / i_peak
    i_norm_area = float(np.trapezoid(i_norm, wavelengths))
    shape_overlap = (
        float(np.trapezoid(a_norm * i_norm, wavelengths)) / i_norm_area
        if i_norm_area > 0.0
        else 0.0
    )

    return {
        "sample_id": sample_id,
        "light_source": "LEDB4",
        "light_source_unit": "relative",
        "broadening": "gaussian",
        "sigma_ev": SIGMA_EV,
        "reference_concentration_molar": CONCENTRATION_M,
        "reference_path_cm": PATH_CM,
        "lambda_max_nm": round(lambda_max_nm, 4),
        "molar_absorptivity_max_M_1_cm_1": round(epsilon_max, 4),
        "absorptance_max": round(absorptance_max, 6),
        "absorbed_flux": round(absorbed_flux, 6),
        "light_flux_total": round(light_flux_total, 6),
        "absorbed_fraction": round(absorbed_fraction, 6),
        "shape_overlap": round(shape_overlap, 6),
        "num_excited_states": len(states),
    }


def main(in_dir, out_dir):
    in_dir = Path(in_dir)
    out_dir = Path(out_dir)
    if not in_dir.exists():
        print(f"ERROR: input directory not found: {in_dir}")
        return 1
    if not LED_CSV.exists():
        print(f"ERROR: LED data file not found: {LED_CSV}")
        return 1
    out_dir.mkdir(parents=True, exist_ok=True)

    out_files = sorted(p for p in in_dir.rglob("*.out") if p.is_file())
    if not out_files:
        print(f"No .out files found under: {in_dir}")
        return 1

    print(f"Scanning {len(out_files)} .out file(s) under {in_dir} ...")
    results = []
    skipped = []
    for path in out_files:
        sample_id = extract_chk_name(path) or path.stem
        if not is_tddft_output(path):
            skipped.append({"sample_id": sample_id, "reason": "not a TD-DFT output (likely OPT)"})
            continue
        states = parse_td_output(path)
        if not states:
            skipped.append({"sample_id": sample_id, "reason": "no Excited State lines parsed"})
            continue
        row = analyze_states(sample_id, states)
        results.append(row)
        print(
            f"  [ok] {sample_id}: states={row['num_excited_states']}, "
            f"lambda_max={row['lambda_max_nm']} nm, "
            f"eps_max={row['molar_absorptivity_max_M_1_cm_1']:.0f} M^-1 cm^-1, "
            f"absorbed_fraction={row['absorbed_fraction']:.4f}, "
            f"absorbed_flux={row['absorbed_flux']:.4f} (relative)"
        )

    if not results:
        print("No analyzable samples produced; see skipped.json if any.")

    if results:
        pd.DataFrame(results).to_csv(out_dir / "results.csv", index=False)
        (out_dir / "results.json").write_text(
            json.dumps(results, indent=2), encoding="utf-8"
        )
    if skipped:
        (out_dir / "skipped.json").write_text(
            json.dumps(skipped, indent=2), encoding="utf-8"
        )

    print()
    print(f"Produced: {len(results)} result(s)   Skipped: {len(skipped)}")
    print(f"Output directory: {out_dir.resolve()}")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python standalone/run.py <input_dir> <output_dir>")
        sys.exit(2)
    sys.exit(main(sys.argv[1], sys.argv[2]))
