from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"
AURORA_DIR = DATA_DIR / "aurora"
TCGA_RAW_DIR = DATA_DIR / "tcga_raw"
TCGA_PROCESSED_DIR = DATA_DIR / "tcga_processed"
RESULTS_DIR = REPO_ROOT / "results"
TABLES_DIR = RESULTS_DIR / "tables"
FIGURES_DIR = RESULTS_DIR / "figures"


def ensure_output_dirs():
    for path in (AURORA_DIR, TCGA_RAW_DIR, TCGA_PROCESSED_DIR, TABLES_DIR, FIGURES_DIR):
        path.mkdir(parents=True, exist_ok=True)


def require_columns(frame, columns, source):
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"{source} is missing required columns: {', '.join(missing)}")


def tcga_sample_type(sample_ids):
    return sample_ids.astype(str).str.split("-").str[3].str[:2]
