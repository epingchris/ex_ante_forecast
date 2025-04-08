# Generating _ex ante_ forecasts for planned REDD+ projects

This project generates ex ante forecasts of carbon outcomes of REDD+ projects and evaluates methods of counterfactual estimation. As input, it uses the output of the implementation code of the Canopy PACT 2.0 methodology for tropical forest carbon accreditation (https://github.com/quantifyearth/tmf-implementation) code, which are parquet files for a given set of projects and parquet files for their matched pixels.

This project adopts the method judged to be the most suitable based on the analysis in (epingchris/placebo_evaluation) to generate _ex ante_ forecasts of counterfactual carbon loss rates (MgC ha-1 yr-1) for ongoing REDD+ projects, and compares them against ex post estimated additionality.

## Requirements

This project is developed under R 4.2, and requires the packages `tidyverse`, `magrittr`, `units`, `sf`, `arrow`, `MatchIt`, `boot`, `scales`, `Metrics`, and `patchwork`. The `sf` package runs on GDAL 3.10.

## Structure

The project contains a core script (forecast.R), which is divided into the following five parts:
0. Setup
A. Obtain observed additionality
B. Predict deforestation probability of baseline pixels using logistic regression
C. Calculate boostrapped baseline C loss
D. Generate additionality forecast and estimate overclaiming risk

It requires the following input variables to read TMF implementation output and other data.
All variables are vectors containing one value for each project to be analysed:

1. analysis_type: "ongoing" for real, ongoing REDD+ projects; "control" for randomly selected placebo "projects"

2. projects: an index of all projects to be analysed
This should correspond to the filenames of the shapefiles and to the -p argument in the implementation code
It is usually the ongoing projects' VCS ID or customised (e.g. prefixed series of integers)

3. pair_dirs: absolute paths of the directories containing all matched pair pixel sets (typically "/pairs/xxx.parquet" and  "/pairs/xxx_matchless.parquet")
The directory should containing pairs of parquet files with the same file name, with and without the "_matchless" suffix.
This is used to calculate estimated observed additionality.

4. k_paths: absolute paths of the set K (typically "k.parquet")
5. m_paths: absolute paths of the set M (typically "matches.parquet")
Both should be in parquet format, containing the following columns:
"lat", "lng" (degrees), "slope" (degrees), "elevation" (metres), "access" (remoteness, in minutes), "cpc10_u", "cpc10_d", "cpc0_u", "cpc0_d" (from 0 to 1), "luc_[t-10]" to "luc_2021" (categorical, 1-6, based on the JRC-TMF dataset, Vancutsem et al. 2021), "ecoregion" (categorical, based on the RESOLVE dataset, Dinerstein et al. 2017)

6. acd_paths: absolute paths of the carbon density per LUC (typically "carbon-density.csv")
This should be an csv file (although the script could be modiified in the future to support txt format) containing columns "land.use.class" (categorical, 1-6) and "carbon.density" (MgC/ha) for all six LUCs, although the script checks and fill missing LUC with NAs

7. polygon_paths: absolute paths of the shapefile of the project extent
This should be a geojson file containing valid geometries in WGS84 (EPSG: 4326), although the script checks for both conditions.
This is currently only used to calculate project area (ha), but could be useful for other purposes in the future.

8. country: country of the project
9. t0: year of start of the project (real or hypothetical)

10. out_path: absolute path of the directory where outputs are to be saved; include file prefix if desired
11. OPTIONAL: fig_path: absolute path of the directory where output figures are to be saved


It generates the following output:

A.
1. project_var: project-level variables
This is a csv file containings N rows (one for each project), and the columns "ID" (project code), "COUNTRY" (country where the project is located), "t0" (project start year), "area_ha" (project area in hectares), and "acd_undisturbed" (aboveground carbon density estimated for the land class of undisturbed forests).

2. acd_df: aboveground carbon density per land class
This is a csv file containing N rows (one for each project), and six columns (class_1 to class_6) for the aboveground carbon density estimated for the each of of six land classes (undisturbed forest, degraded forest, deforested land, forest regrowth, water surface, and others).

3. additionality: estimated additionality time series

One csv file is generated for each project ("_additionality_[project ID].csv"). Each row is a year (starting from 10 years before t0) for one of the 100 sample pairs. The csv file contains the following columns:
- year: year of the estimation
- c_loss: counterfactual carbon loss (MgC ha-1 yr-1)
- t_loss: project carbon loss (MgC ha-1 yr-1)
- Additionality (c_loss - t_loss): additionality (amount of carbon credits) (MgC ha-1 yr-1)
- pair: sample pair ID
- started: boolean, TRUE if the year > t0, FALSE otherwise
- project: project ID

4. Estimated baseline carbon loss (MgC ha-1 yr-1)
- baseline_best: "close matching" baseline
- baseline_loose: "loose matching" baseline
- baseline_lagged: "time-lagged matching" baseline

One csv file is generated for each project and each baseline type ("_baseline_[type]_[project ID].csv"). Close matching and time-lagged matching baselines contain 1000 rows (100 sample pairs x 10 years before the project start), whereas the loose matching baseline has the same amount of rows as the set M pixels in the implementation code, downsampled to 250000 if set M size is larger than that. The csv file contains the following columns:
- c_loss: baseline carbon loss (MgC ha-1)
- project: project ID


It is currently a list of N element, saved as an RDS file that can be read as an R object.
Each element is a data frame containing the columns "year", "c_loss" (counterfactual carbon loss), "t_loss" (project carbon loss), "additionality", "pair" (index of the sampled pairs), "started" (FALSE for years before and including t0, TRUE for years after t0). Each row is a combination of a year and a sampled pair.

B.
1. baseline: predicted baseline deforestation probability
It is currently a list of N elements, saved as an RDS file that can be read as an R object.
Each element is a data frame containing the columns "lat", "lng", "slope", "elevation", "access", "luc10", "luc0", "cpc10_u", "cpc10_d", "cpc0_u", "cpc0_d", "defor" (TRUE if land use class has changed from 1 to 2, 3, or 4 between t-10 and t0, FALSE otherwise), "defor_prob" (predicted deforestation probability based on logistic regression, between 1 an 0), "risk" ("low" if defor_prob < 0.01, "high" otherwise), "acd10", "acd0", "c_loss" (annual carbon loss rate). Each row is a baseline pixel.

2. #project_defor_prob: predicted project deforestation probability
It is currently a list of N elements, saved as an RDS file that can be read as an R object.
Each element is a data frame containing the columns "lat", "lng", "slope", "elevation", "access", "luc10", "luc0", "cpc10_u", "cpc10_d", "cpc0_u", "cpc0_d", "defor" (TRUE if land use class has changed from 1 to 2, 3, or 4 between t-10 and t0, FALSE otherwise), "defor_prob" (predicted deforestation probability based on logistic regression, between 1 an 0), "risk" ("low" if defor_prob < 0.01, "high" otherwise), "acd10", "acd0", "c_loss" (annual carbon loss rate). Each row is a project pixel.

3. range_defor_prob: total range of predicted baseline deforestation probability in baseline pixels
It is a data frame of N columns, each for a project. Each column contains the minimum and maximum of predicted baseline deforestation probability.

5. baseline_summary: basic information about the baseline
This is a csv file containing N rows (one for each project), and the columns "project", "baseline_area" (in number of pixels), "low_risk_ratio" (ratio of baseline pixels classified as "low-risk" based on a threshold of 1% predicted deforestation probability), and "slope", "elevation", "access" (indicating what is the logistic regression result for each of the environmental variables: "Neg." for negative effect, "Pos." for positive effect, and "N.S." for non-significant effect)

C.
1. df_ses: standardised effect size of change in bootstrapped carbon loss by using only high-risk pixels instead of all pixels in 
This is a csv file containing N rows (one for each project), and the columns "project", "ses" (standardised effect size).
The standardised effect size is calculated as ("mean carbon loss rate in high-risk pixels" - "mean carbon loss rate in all pixels") / "standard deviation of carbon loss rate in all pixels"

2. c_loss_boot: bootstrapped baseline annual carbon loss rates for all projects
This is a csv file containing a long-form data frame, containing the columns "project", "type" ("all" for all pixels, "high_risk" for only high-risk pixels), and "val" (each bootstrapped carbon loss rate). The number of values generated (number of rows in each project-type combination) is defined by the variable "boot_n"

D.
1. forecast_summ: additionality forecast under different scenarios
This is a csv file containing N rows (one for each project), and the columns "mean_100", "ci_100", "mean_75", "ci_75", "mean_50", "ci_50", "mean_25", "ci_25", and "project". The columns prefixed "mean_" indicate the forecasted mean annual additionality (MgC/ha/yr), and the columns prefixed "ci_" indicate its confidence interval, under the scenarios of 100%, 75%, 50%, and 25% project effectiveness, respectively.

Apart from these outputs, the script also generates a number of plots for visualisation. Currently they are all used for figures in E-Ping's manuscript, and can be deactivated by setting the variable "visualise" to FALSE.