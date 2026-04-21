"""Nesting script: HurryWave (detail) nested within HurryWave (overall).

Workflow:
  1. Load overall and detail HurryWave models.
  2. nest1 — add spectral observation points to the overall model at the
     detail model's boundary locations, then write the updated overall model.
  3. Run the overall HurryWave model (manual step).
  4. nest2 — read spectral output from the overall model and apply it as
     boundary conditions to the detail model, then write the detail model.
"""
#%%
from pathlib import Path

from hydromt_hurrywave import HurrywaveModel
from cht_nesting import nest1, nest2

# ---------------------------------------------------------------------------
# Configuration — edit these paths
# ---------------------------------------------------------------------------

OVERALL_ROOT = Path(r"p:\11212696-fhics-bes\1_models\HurryWave\02_modelruns\nesting_test_amu_amv\overall")
DETAIL_ROOT   = Path(r"p:\11212696-fhics-bes\1_models\HurryWave\02_modelruns\nesting_test_amu_amv\detail")

# Prefix used to identify observation points in the overall model output.
# Must be unique and match what was used in nest1.
OBS_PREFIX = "detail"

# Directory that contains the overall model's output after it has been run.
# Defaults to OVERALL_ROOT if left as None.
OVERALL_OUTPUT_PATH = OVERALL_ROOT / "input"

# Directory to write detail model boundary condition files.
# Set to None to skip writing to disk (BCs are still set on the model object).
DETAIL_BC_PATH = DETAIL_ROOT / "input_0_05"

# ---------------------------------------------------------------------------
# Step 1 — add observation points to the overall model
# ---------------------------------------------------------------------------
#%%
print("Loading models ...")
overall = HurrywaveModel(root=str(OVERALL_ROOT/ "input"), mode="r+")
detail  = HurrywaveModel(root=str(DETAIL_ROOT/ "input_0_05"),  mode="r+")

overall.read()
detail.read()

#%%
# print("nest1: adding spectral observation points to overall model ...")
# success = nest1(overall, detail, obs_point_prefix=OBS_PREFIX)
# if not success:
#     raise RuntimeError("nest1 failed — HurryWave-in-HurryWave combination not supported.")

# overall.write()
# print(f"Overall model written to: {OVERALL_ROOT}")
# print()
# print(">>> Run the overall HurryWave model now, then re-run this script from Step 2.")
# print("    (Comment out everything above the Step 2 block once done.)")

#%% ---------------------------------------------------------------------------
# Step 2 — apply overall spectral output as BCs on the detail model
# ---------------------------------------------------------------------------

# Uncomment the block below after the overall model has been run.

print("nest2: reading spectral output and setting detail model BCs ...")

print(f"Overall output path: {OVERALL_OUTPUT_PATH}")
print(f"Detail BC output path: {DETAIL_BC_PATH}")

nest2(
    overall=overall,
    detail=detail,
    obs_point_prefix=OBS_PREFIX,
    output_path=str(OVERALL_OUTPUT_PATH) if OVERALL_OUTPUT_PATH else None,
    bc_path=str(DETAIL_BC_PATH) if DETAIL_BC_PATH else None,
)
detail.write()
print(f"Detail model BCs updated and written to: {DETAIL_BC_PATH}")

# %%
