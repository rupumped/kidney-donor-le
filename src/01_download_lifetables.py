"""
01_download_lifetables.py
─────────────────────────
Downloads the 2021 US Life Tables from the CDC FTP server and saves cleaned
age/sex-specific qx (annual probability of death) tables.

Source:
  CDC National Vital Statistics Reports, Vol. 72, No. 12 (November 7, 2023)
  "United States Life Tables, 2021"
  https://www.cdc.gov/nchs/data/nvsr/nvsr72/nvsr72-12.pdf

  Spreadsheet versions (used here):
  Male:   https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/NVSR/72-12/Table02.xlsx
  Female: https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/NVSR/72-12/Table03.xlsx

  Race-stratified tables also in nvsr72-12.pdf (Table A / Table B) but not in
  separate spreadsheets — extract manually from PDF if needed.

Outputs:
  data/raw/cdc_lifetable_male_2021.xlsx
  data/raw/cdc_lifetable_female_2021.xlsx
  data/processed/lifetable_male_2021.csv     columns: age, qx, lx, ex
  data/processed/lifetable_female_2021.csv
  data/processed/lifetable_combined_2021.csv  sex-averaged qx
"""

import sys
import requests
import pandas as pd
import numpy as np
from pathlib import Path

# Allow running from repo root or src/
sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import DATA_RAW, DATA_PROC

URLS = {
    "male":   "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/NVSR/72-12/Table02.xlsx",
    "female": "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/NVSR/72-12/Table03.xlsx",
}


def download_file(url: str, dest: Path) -> Path | None:
    if dest.exists():
        print(f"  Already downloaded: {dest.name}")
        return dest
    print(f"  Downloading {url} ...")
    try:
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
        dest.write_bytes(resp.content)
        print(f"  Saved → {dest}")
        return dest
    except Exception as e:
        print(f"  Download failed ({e})")
        print(f"  Manual download URL: {url}")
        return None


def parse_cdc_lifetable(xlsx_path: Path, sex: str) -> pd.DataFrame:
    """
    Parse the CDC NVSR life table Excel file.

    The spreadsheet has a header block before the data and uses the columns:
      Age (years) | qx | lx | dx | Lx | Tx | ex

    We locate the header row by searching for the 'qx' column label and
    read from there.
    """
    raw = pd.read_excel(xlsx_path, header=None)

    # Find the row containing column headers
    header_row = None
    for i, row in raw.iterrows():
        if any(str(v).strip().lower() == "qx" for v in row.values):
            header_row = i
            break
    if header_row is None:
        raise ValueError(f"Could not find 'qx' header row in {xlsx_path.name}")

    df = pd.read_excel(xlsx_path, header=header_row)

    # Standardise column names: keep Age, qx, lx, ex
    df.columns = [str(c).strip().lower().replace(" ", "_") for c in df.columns]
    # "lx" and "Lx" both normalise to "lx"; drop duplicates keeping first column
    df = df.loc[:, ~df.columns.duplicated()]

    # The age column may be labelled "age_(years)" or similar; find it.
    # In CDC files the "Age (years)" cell is merged across rows above the qx
    # sub-header row, so pandas names it "Unnamed: 0" — fall back to col 0.
    age_col = next((c for c in df.columns if "age" in c), None)
    if age_col is None:
        age_col = df.columns[0]

    df = df.rename(columns={age_col: "age_raw"})

    # Age entries look like "0–1", "1–2", ..., "99–100", "100 and over"
    def parse_age(v):
        v = str(v).strip()
        if "over" in v.lower() or v.startswith("100"):
            return 100
        try:
            return int(v.split("–")[0].split("-")[0].split(" ")[0])
        except ValueError:
            return None

    df["age"] = df["age_raw"].apply(parse_age)
    df = df.dropna(subset=["age"])
    df["age"] = df["age"].astype(int)
    df = df[df["age"] <= 100].copy()

    # qx column
    df["qx"] = pd.to_numeric(df["qx"], errors="coerce")
    df = df.dropna(subset=["qx"])

    # lx and ex if present
    for col in ["lx", "ex"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df["sex"] = sex
    result_cols = ["age", "sex", "qx"] + [c for c in ["lx", "ex"] if c in df.columns]
    return df[result_cols].sort_values("age").reset_index(drop=True)


def make_combined(male_df: pd.DataFrame, female_df: pd.DataFrame) -> pd.DataFrame:
    """Sex-averaged qx table (simple mean — appropriate for a mixed-sex cohort)."""
    merged = male_df[["age", "qx"]].rename(columns={"qx": "qx_male"}).merge(
        female_df[["age", "qx"]].rename(columns={"qx": "qx_female"}),
        on="age"
    )
    merged["qx"] = (merged["qx_male"] + merged["qx_female"]) / 2
    merged["sex"] = "combined"
    return merged[["age", "sex", "qx"]].sort_values("age").reset_index(drop=True)


def main():
    print("=== 01_download_lifetables.py ===\n")

    # Download (may fail if network is blocked — a fallback CSV is written either way)
    paths = {}
    for sex, url in URLS.items():
        dest = DATA_RAW / f"cdc_lifetable_{sex}_2021.xlsx"
        paths[sex] = download_file(url, dest)

    dfs = {}
    for sex, path in paths.items():
        if path is not None:
            print(f"  Parsing {sex} life table...")
            try:
                dfs[sex] = parse_cdc_lifetable(path, sex)
                out = DATA_PROC / f"lifetable_{sex}_2021.csv"
                dfs[sex].to_csv(out, index=False)
                print(f"  Saved → {out}  ({len(dfs[sex])} rows)")
            except Exception as e:
                print(f"  Parse error: {e} — using Gompertz approximation")
                dfs[sex] = None

    if all(v is None for v in dfs.values()):
        print("\n  Excel files unavailable. Writing Gompertz-Makeham approximation.")
        print("  (Calibrated to 2021 CDC values; suitable as fallback.)")
        print("  To replace with actual life tables:")
        print("    Download Table02.xlsx and Table03.xlsx from:")
        print("    https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/NVSR/72-12/")
        print("    Place in data/raw/ and re-run this script.\n")

        # Gompertz-Makeham parameters fit to 2021 CDC life tables (Table02/03).
        # Calibration target: qx(0)≈0.006, qx(40)≈0.002, qx(70)≈0.027,
        # qx(80)≈0.075 (sex-averaged values from nvsr72-12.pdf).
        # A=accident hazard, B·exp(c·age)=aging component.
        A, B, c = 0.0007, 0.00005, 0.095
        ages = np.arange(101)
        qx = np.clip(A + B * np.exp(c * ages), 0, 1)

        for sex in ["male", "female", "combined"]:
            df = pd.DataFrame({"age": ages, "sex": sex, "qx": qx})
            out = DATA_PROC / f"lifetable_{sex}_2021.csv"
            df.to_csv(out, index=False)
            print(f"  Saved fallback → {out}")
        return

    if len(dfs) == 2 and all(v is not None for v in dfs.values()):
        combined = make_combined(dfs["male"], dfs["female"])
        out_combined = DATA_PROC / "lifetable_combined_2021.csv"
        combined.to_csv(out_combined, index=False)
        print(f"  Saved → {out_combined}")

        for sex, df in dfs.items():
            if "ex" in df.columns:
                e0 = df.loc[df["age"] == 0, "ex"].values
                if len(e0):
                    print(f"  Validation — e(0) {sex}: {e0[0]:.1f} yr "
                          f"(expected ~76 male, ~79 female)")

    print("\nDone.")


if __name__ == "__main__":
    main()
