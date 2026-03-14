# Climate-Conditioned Catastrophe Modeling for Dynamic Risk Assessment

This repository contains code and data necessary to reproduce the results presented in the article "Climate-Conditioned Catastrophe Modeling for Dynamic Risk Assessment (preprint available [here](add url)). The following files are needed for this:

## Data

### Publicly Available

The `Data/` folder contains the publicly available input data:

* **Data/ILW_prices.csv**: Prices quoted for different Industry Loss Warranties (ILWs).
* **Data/IndexationFactors.csv**: Loss indexation factors, used to adjust historical losses to present-day values based on inflation and estimated changes in exposure. This follows the method used in [Comola et al., 2024](https://www.nature.com/articles/s43247-024-01824-7#Sec5)
* **Data/Landfall.csv**: Landfall region of historical storms
* **Data/Reask_ClimateIndices.parquet**: Season-by-ensemble valuse of August-October mean values for various indices calculated from the gridded weather model output used to force the Reask Unified Tropical Cyclone (UTC) model. These are used to generate Figure 2b-d.
* **CBRA/Output/STD_1951-2020_[APRIL,JUNE]**: These folders contain the seasonally adjusted loss tables for each seasonal forecast from 1985-2024, for both the April and June forecast initialization dates. They are produced by the Climate-Based Risk Adjuster (CBRA), using data from `landfall-data.parquet` and `ylt_Verisk_STD_cat.csv`.

### Restricted Data

This data in the `Restricted_Data` folder is subject to commercial confidentiality agreements and cannot be publicly disclosed. Contact the corresponding author ([francesco.comola@lgtcp.com](mailto:francesco.comola@lgtcp.com)) for additional information.

* **Data/Restricted_Data/PCS_RegionalSplit.csv**: The regional distribution of historical losses from the industry-standard source for claims data: Verisk's Property Claim Services (PCS).
* **Data/Restricted_Data/Verisk_STD_RegionalSplit.csv**: The regional distribution of modeled losses in the static loss model used in this analysis.
* **Data/Restricted_Data/ylt_Verisk_STD_cat.csv**: The static Year-Event Loss Table (YELT) used in this analysis.
* **Data/Restricted_Data/landfall-data.parquet**: Long-term and seasonal landfall simulations: Reask's generative Unified Tropical Cyclone (UTC) model generates synthetic event sets forced by weather model output. For this analysis, we use an event set forced by the ERA5 reanalysis over 1950-2024 (used to define our "baseline" climate) and event sets forced by the European Center for Medium-range Weather Forecasting (ECMWF)'s 5-member SEAS5 ensemble forecast initialized in April and June of each year from 1985-2024. This represents a backend database for the Climate-Based Risk Adjuster (CBRA) tool used to adjust loss tables reflect forecasted risk in each season (see `adjust-ylt.py`). While this database is not publicly available, the outputs of the CBRA process, e.g. loss tables for each season, are included in the public data associated with this repository.

## Code

* **environment.yml**: A conda environment that will load the packages you need to execute `Investment Simulation.R`.
* **adjust-ylt.py**: A python script that adjusts the long-term View of Risk (VoR), represented by the industry standard catastrophe model output in the form of a Year-Event Loss Table (YELT) (`ylt_Verisk_STD_cat.csv`). This output is adjusted to season-specific YELTs for two different forecast initialization months. June (used for the main results of the paper) and April (shown in the supplementary information). Note that this code is provided for transparency but is not runnable with the
* **Investment_Simulation.R**: the R script that simulates the investment strategies based on the seasonally-adjusted risk model output and the long-term risk model output. The script reads in the output of `adjust-ylt.py` and other input data contained in the `Data/` folder, structures the data using data.tables, performes the investment calculations, and generates the plots shown in the manuscript (with the exception of Fig. 2a, c, e, which are generated from the CBRA package).

## Replication

To replicate the findings in the paper, execute the code with the following steps:

1. Obtain all data listed above. You will need both the publicly available and the restricted data to fully replicate the results.

2. Install and activate the conda environment with the necessary packages

   a. First, put your Reask account credentials in a file at `~/.netrc`. This will allow pip to install CBRA from the private pyPI server on which the package is hosted. Your file should look like this:

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

3. (optional) Run `python adjust-ylt.py`. To do this, you will need `landfall-data.parquet`. This step is not necessary, as the outputs for this step are already provided in the publicly available data under the `CBRA/` folder.

4. Run `Rscript Investment_Simulation.R`. This will generate all of the figures and tables presented in the manuscript.

## License

See the `LICENSE` file for the license associated with all source code provided within this repository. It does not apply to CBRA, which has its own license.
