from pathlib import Path

import numpy as np
from pylt import adjust_ylt

# Change to MIDAPRIL for forecasts initiated in April, or keep as MIDJUNE for forecasts
# initiated in June. By default, June forecasts are used.
FCT_TIMING = "MIDJUNE"

# change this to the correct path to your data
proj_dir = Path("/Users/ianbolliger/git-repos/ilw-paper")

data_dir = proj_dir / "DATA"
restricted_dir = data_dir / "RESTRICTED"
INPUT_YLT_PATH = restricted_dir / "ylt_Verisk_STD_cat.csv"

COUNTS_PATH = restricted_dir / "Reask_landfall-data.parquet"
METRICS_PATH = data_dir / "Reask_ClimateIndices.parquet"
GATES_PATH = data_dir / "Reask_gates.parquet"
SMOOTHING_MAP_PATH = data_dir / "Reask_smoothing-mapper.parquet"

# Adjust to seasonal forecasts:
# change fct_timing to either "MIDJUNE" or "MIDAPRIL" for forecasts initiated in June
# or April, respectively. By default, June forecasts are used.

for y in range(1985, 2025):
    try:
        print(y)
        sdir = data_dir / "CBRA_YLT" / f"STD_1951-2020_{FCT_TIMING[3:]}" / str(y)
        sdir.mkdir(parents=True, exist_ok=True)
        adjust_ylt(
            INPUT_YLT_PATH,
            "SEAS5",
            baseline_yrs=np.arange(1951, 2021),
            target_yrs=y,
            fct_timing=FCT_TIMING,
            gates_path=GATES_PATH,
            metrics_path=METRICS_PATH,
            ldf_data_path=COUNTS_PATH,
            smoothing_mapper=SMOOTHING_MAP_PATH,
            ylt_bypassing_col=None,
            max_point_to_gate_dist_deg=2,
            intensity_units="category",
            save_dir=sdir,
        )
        break
    except:
        print("Error in YLT resampling")
        raise
