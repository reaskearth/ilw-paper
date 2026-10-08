# Climate-Conditioned Catastrophe Modeling for Dynamic Risk Assessment

This repository contains code and data necessary to reproduce the results presented in the article "Climate-Conditioned Catastrophe Modeling for Dynamic Risk Assessment (preprint available [here](https://doi.org/10.21203/rs.3.rs-9124834/v2)). The following files are needed for this:

## Data

### Publicly Available

The `DATA/` folder contains the publicly available input data:

* **DATA/ILW_prices.csv**: Prices quoted for different Industry Loss Warranties (ILWs).
* **DATA/indexation-factors.csv**: Loss indexation factors, used to adjust historical losses to present-day values based on inflation and estimated changes in exposure. This follows the method used in [Comola et al., 2024](https://www.nature.com/articles/s43247-024-01824-7#Sec5)
* **DATA/landfall.csv**: Landfall region of historical storms. This file can be used to benchmark modeled landfall frequencies against historical data.
* **DATA/Reask_climate-indices.parquet**: Season-by-ensemble values of August-October means for various indices calculated from the gridded weather model output used to force the Reask Unified Tropical Cyclone (UTC) model. These are used to generate Figure 2b-d.
* **DATA/Reask_gates.parquet**: Landfall gates used in CBRA.
* **DATA/Reask_smoothing-mapper.parquet**: By default, nearby gates are smoothed to improve convergence in CBRA. This crosswalk file defines which gates should be smoothed to which other gates.
* **DATA/CBRA_YLT/STD_1951-2020_[APRIL,JUNE]**: These folders contain the seasonally adjusted loss tables for each seasonal forecast from 1985-2024, for both the April and June forecast initialization dates. They are produced by the Climate-Based Risk Adjuster (CBRA), using data from `Reask_landfall-data.parquet` and `Verisk_STD_YLT.csv`.

### Restricted Data

This data in the `DATA/RESTRICTED` folder is subject to commercial confidentiality agreements and cannot be publicly disclosed. Contact the corresponding author ([francesco.comola@lgtcp.com](mailto:francesco.comola@lgtcp.com)) for additional information.

* **DATA/RESTRICTED/PCS_regional-split.csv**: The regional distribution of historical losses from the industry-standard source for claims data: Verisk's Property Claim Services (PCS).
* **DATA/RESTRICTED/Reask_landfall-data.parquet**: Long-term and seasonal landfall simulations: Reask's generative Unified Tropical Cyclone (UTC) model generates synthetic event sets forced by weather model output. The file holds an event set forced by the ERA5 reanalysis (seasons 1940-2024, 2,500 simulations per season) and event sets forced by the European Center for Medium-range Weather Forecasting (ECMWF)'s SEAS5 ensemble forecast initialized in April and June of each year from 1981-2024 (1,000 simulations per ensemble member; 25 members through 2016 and 51 from 2017). `adjust-ylt.py` uses the 1951-2020 ERA5 seasons to define our "baseline" climate and adjusts the loss tables for the 1985-2024 forecast seasons used in the investment simulation; `likelihood-test-figure.py` scores all 1981-2024 seasons of the SEAS5 and ERA5 event sets. This represents a backend database for the Climate-Based Risk Adjuster (CBRA) tool used to adjust loss tables to reflect forecasted risk in each season (see `adjust-ylt.py`). While this database is not publicly available, the anonymized outputs of the CBRA process, e.g. loss tables for each season, are included in the public data associated with this repository.
* **DATA/RESTRICTED/Reask_ibtracs-landfall-data.parquet**: Observed tropical cyclone landfalls at the CONUS landfall gates in `Reask_gates.parquet`, 1981-2024, from IBTrACS best tracks (USA agency), linearly interpolated to hourly track points and processed by Reask through the same gate-crossing code as the synthetic event sets. One row per landfall, with the maximum sustained wind at the last hourly point before landfall. Used by `likelihood-test-figure.py` to score the seasonal and static event sets against the observed record.
* **DATA/RESTRICTED/Verisk_STD_regional-split.csv**: The regional distribution of modeled losses in the static loss model used in this analysis.
* **DATA/RESTRICTED/Verisk_STD_YLT.csv**: The static Year-Event Loss Table (YELT) used in this analysis (AIR North Atlantic Hurricane Model, as implemented in TouchstoneRe v13).
* **DATA/RESTRICTED/event-mapping.csv**: Mapping file to match the event IDs in the adjusted YLTs with those in the original YLT. This mapping makes it possible to split the industry losses in the adjusted YLTs geographically in the same proportions as in the original YLT.

## Code

* **environment.yml**: A conda environment that will load the packages you need to execute `adjust-ylt.py`, `likelihood-test-figure.py` and `investment-simulation.R`.
* **likelihood-test-figure.py**: A python script that produces the landfall likelihood-ratio figure. It scores the observed record of CONUS landfalls (`Reask_ibtracs-landfall-data.parquet`) under the seasonally forced SEAS5 event sets and under static event sets pooled over 1981-2024 (both from `Reask_landfall-data.parquet`), and tests the difference in skill. The figure and its data tables are written to `OUTPUT/`.
* **adjust-ylt.py**: A python script that adjusts the long-term View of Risk (VoR), represented by the industry standard catastrophe model output in the form of a Year-Event Loss Table (YELT) (`Verisk_STD_YLT.csv`). This output is adjusted to season-specific YELTs for two different forecast initialization months. June (used for the main results of the paper) and April (shown in the supplementary information).
* **investment-simulation.R**: the R script that simulates the investment strategies based on the seasonally-adjusted risk model output and the long-term risk model output. The script reads in the output of `adjust-ylt.py` and other input data contained in the `DATA/` folder, structures the data using data.tables, performes the investment calculations, and generates the plots shown in the manuscript (with the exception of Fig. 2a, c, e, which are generated from the CBRA package).

## Replication

To replicate the findings in the paper, you will need access both to CBRA (if you wish to run the optional steps 3 and 4 below) and to the restricted data (see above). CBRA is made available to Reask clients and for additional noncommercial research purposes on a case-by-case basis. To request access to CBRA, email [contact@reask.earth](mailto:contact@reask.earth) with a description of your use case. The restricted data is subject to third-party commercial confidentiality agreements and may not be available. The scripts in `CODE/` look for `DATA/` and write `OUTPUT/` one level up, so run them from the `CODE/` directory (or edit the project path at the top of each script).

1. Obtain all data listed above.

2. Install and activate the conda environment with the necessary packages

   a. First, put your Reask account credentials in a file at `~/.netrc`. This will allow pip to install CBRA from the private pyPI server on which the package is hosted. If you do not wish to run the optional steps 3 and 4, CBRA is not needed and you may remove the corresponding line from `environment.yml` and skip to step 2b. Once obtained, your `~/.netrc` file should look like this:

      ```bash
      machine pypi.reask.earth
      login your_username
      password your_password
      ```

   b. Next, create and activate the conda environment:

      ```bash
      conda env create -f environment.yml
      conda activate ilwpaper
      ```

3. (optional) Run `python adjust-ylt.py`. To do this, you will need `Reask_landfall-data.parquet`. This step is not necessary, as the outputs for this step are already provided in the publicly available data under the `DATA/CBRA_YLT/` folder.

4. (optional) Run `python likelihood-test-figure.py`. To do this, you will need `Reask_landfall-data.parquet` and `Reask_ibtracs-landfall-data.parquet`. This generates the landfall likelihood-ratio figure, written to `OUTPUT/`. Step 5 does not depend on it.

5. Run `Rscript investment-simulation.R`. This will generate all of the other figures and tables presented in the manuscript.

## License

See the `LICENSE` file for the license associated with all source code provided within this repository. It does not apply to CBRA, which has its own license.
