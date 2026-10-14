"""
Run the whole analysis from the anonymised data and check the results against expected/.

    python run_all.py                     # everything, about 5 minutes plus the invariance models (10 to 20 minutes)
    python run_all.py --skip-invariance   # reuse the invariance results already in output/ (copied from expected/
                                          # when output/ has none)
    python run_all.py --rscript /path/to/Rscript

The last step prints "REPRODUCTION OK" when every result table matches the published results.
"""
import argparse
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("--skip-invariance", action="store_true", help="do not refit the invariance models (script 03)")
ap.add_argument("--rscript", default=os.environ.get("RSCRIPT", "Rscript"), help="path to Rscript")
a = ap.parse_args()

if not os.path.exists(os.path.join(HERE, "data", "trauma_survey_anonymised.csv")):
    sys.exit("Missing data/trauma_survey_anonymised.csv (see data/README.md).")
for folder in ["output", "figures"]:
    os.makedirs(os.path.join(HERE, folder), exist_ok=True)


def run(*cmd):
    t = time.time()
    print(f"\n>>> {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, cwd=HERE, check=True)
    print(f"    done in {time.time() - t:.0f} s", flush=True)


py = sys.executable
run(py, "scripts/01_build_dataset.py")                     # score every questionnaire, define the groups
run(py, "scripts/02_factor_structure.py")                  # factor analysis, group congruence, item bias, pooling check
if a.skip_invariance:
    for ref in ["control", "an_with_imagery"]:
        for kind in ["fit", "freed", "latent_means"]:
            f = f"03_{'invariance_' if kind != 'latent_means' else ''}{kind}_{ref}.csv"
            if not os.path.exists(os.path.join(HERE, "output", f)):
                shutil.copy2(os.path.join(HERE, "expected", f), os.path.join(HERE, "output", f))
else:                                                      # ordinal multi-group CFA, both comparisons in parallel
    print("\n>>> measurement invariance, both comparisons in parallel (10 to 20 minutes)", flush=True)
    procs = [subprocess.Popen([a.rscript, "scripts/03_invariance.R", ref], cwd=HERE) for ref in ["control", "an_with_imagery"]]
    if any(p.wait() for p in procs):
        sys.exit("script 03 failed")
run(py, "scripts/04_trauma_analyses.py")                   # every trauma analysis in the paper
run(py, "scripts/05_tables.py")                            # the paper's tables and the mediation estimates
run(py, "scripts/06_figure_factor_structure.py")           # Supplementary Figure S1
run(py, "scripts/07_figure_invariance.py")                 # Supplementary Figure S2
run(py, "scripts/08_figures_trauma.py")                    # Figures 1 to 4
run(py, "scripts/09_check_reproduction.py")                # every result table against expected/
