
# This notebook calls the pylt method to generate ylimate-conditioned year loss tables for the period 1985-2024.


%load_ext autoreload
%autoreload 2

from pathlib import Path

import numpy as np
import polars_st as st

from pylt import adjust_ylt, get_optimal_baseline
from pylt.io import load_gates
from pylt.settings import CLIENT_DATA_DELIVERY_DIR

import os
    
import pylt
pylt.__version__

import cartopy
cartopy.config["data_dir"]


# change this to the correct path to your data

proj_dir = Path("./")
data_dir =  proj_dir / "../DATA/RESTRICTED/"
INPUT_YLT_PATH = data_dir / "Verisk_STD_YLT.csv"

COUNTS_PATH = data_dir / "Reask_landfall-data.parquet"
METRICS_PATH = data_dir / "Reask_climate-metrics.parquet"
GATES_PATH = data_dir / "Reask_gates.parquet"
SMOOTHING_MAP_PATH = data_dir / "Reask_smoothing-mapper.parquet"

# ## Adjust to seasonal forecasts:
# ## change fct_timing to either "july" or "may" for forecasts initiated in June or April, respectively. By default, June forecasts are used.

for y in range(1985,2025):
    try:
        print(y)
        sdir = proj_dir / "../DATA/CBRA_YLT" / "STD_1951-2020_JUNE" / str(y)
        if not os.path.exists(sdir):
            os.makedirs(sdir)
        adjust_ylt(
            INPUT_YLT_PATH,
            "SEAS5",
            baseline_yrs=np.arange(1951, 2021),
            target_yrs=y,
            fct_timing="july",
            gates_path=GATES_PATH,
            metrics_path=METRICS_PATH,
            ldf_data_path=COUNTS_PATH,
            smoothing_mapper=SMOOTHING_MAP_PATH,
            ylt_bypassing_col=None,
            max_point_to_gate_dist_deg=2,
            intensity_units='category',
            save_dir=sdir
        )
    except:
        print("Error in YLT resampling")
