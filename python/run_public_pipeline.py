"""Packages: see requirements.txt. Prerequisites: repository data files; R is run separately."""

import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent


for script in ("03_panel_construction.py", "04_tcga_cleaning.py", "05_manuscript_descriptives.py"):
    print(f"Running {script}", flush=True)
    subprocess.run([sys.executable, str(SCRIPT_DIR / script)], check=True)

subprocess.run(
    [sys.executable, str(SCRIPT_DIR / "06_validate_manuscript_values.py"), "--require-figures"],
    check=True,
)
