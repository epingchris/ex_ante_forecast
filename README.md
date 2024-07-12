This R repository receives the output of the TMF implementation (https://github.com/quantifyearth/tmf-implementation) for N REDD+ projects to be analysed,
and calculates the baseline deforestation and carbon loss, generates additionality forecast and evaluates overclaiming risk. The baseline is currently based on the set M (matches.parquet), and is subject to change in the future.


The core script, forecast.R, contains the following five parts:
0. Setup
A. Obtain observed additionality
B. Predict deforestation probability of baseline pixels using logistic regression
C. Calculate boostrapped baseline C loss
D. Generate additionality forecast and estimate overclaiming risk


It requires the following input variables to read TMF implementation output and other data.
All variables are vectors of length N (number of projects to be analysed):

1. projects: character, an index of all projects to be analysed; it could be the projects' VCS ID or customised (e.g. simply a series of integers)

2. pair_dirs: character, absolute paths of the directories containing all matched pair pixel sets (typically "/pairs/xxx.parquet" and  "/pairs/xxx_matchless.parquet")
The directory should containing pairs of parquet files with the same file name, with and without the "_matchless" suffix.
This is used to calculate estimated observed additionality.

3. k_paths: character, absolute paths of the set K (typically "k.parquet")
4. m_paths: character, absolute paths of the set M (typically "matches.parquet")
Both should be in parquet format, containing the following columns:
"lat", "lng" (degrees), "slope" (degrees), "elevation" (metres), "access" (remoteness, in minutes), "cpc10_u", "cpc10_d", "cpc0_u", "cpc0_d" (from 0 to 1), "luc_[t-10]" to "luc_2021" (categorical, 1-6, based on the JRC-TMF dataset, Vancutsem et al. 2021), "ecoregion" (categorical, based on the RESOLVE dataset, Dinerstein et al. 2017)

5. acd_paths: character, absolute paths of the carbon density per LUC (typically "carbon-density.csv")
This should be an csv file (although the script could be modiified in the future to support txt format) containing columns "land.use.class" (categorical, 1-6) and "carbon.density" (MgC/ha) for all six LUCs, although the script checks and fill missing LUC with NAs

6. polygon_paths: character, absolute paths of the shapefile of the project extent
This should be a geojson file containing valid geometries in WGS84 (EPSG: 4326), although the script checks for both conditions.
This is currently only used to calculate project area (ha), but could be useful for other purposes in the future.

7. country: character, country of the project
8. t0: numerical, year of start of the project (real or hypothetical)
9. OPTIONAL: proj_name: character, full name of the project for readability (if unspecified, the projects variable will be used)
10. out_path: character, absolute paths of the directory where outputs are to be saved; include file prefix if desired


It generates the following output:

A.
1. project_var: project-level variables
This is a csv file containings N rows (one for each project), and the columns "project", "t0", "country", "area_ha", "acd_undisturbed", and the min/max/median over 100 samples of the following variables: "slope", "elevation", "accessibility", "cpc{0, 5, 10}_ {u, d}" (coarse proportional cover at t0/t-5/t-10 of undisturbed/deforested pixels), "defor_ {5_0, 10_0, 10_5}" (deforestation defined as CPC change over the period of t-5 to t0, t-10 to t0, or t-10 to t-5).

2. OPTIONAL output: only basic project-level variables "project", "t0", "country", "area_ha"

3. additonality_estimates: observed additionality time series
It is currently a list of N element, saved as an RDS file that can be read as an R object.
Each element is a data frame containing the columns "year", "c_loss" (counterfactual carbon loss), "t_loss" (project carbon loss), "additionality", "pair" (index of the sampled pairs), "started" (FALSE for years before and including t0, TRUE for years after t0). Each row is a combination of a year and a sampled pair.

B.
1. baseline: predicted baseline deforestation probability
It is currently a list of N elements, saved as an RDS file that can be read as an R object.
Each element is a data frame containing the columns "lat", "lng", "slope", "elevation", "access", "luc10", "luc0", "cpc10_u", "cpc10_d", "cpc0_u", "cpc0_d", "defor" (TRUE if land use class has changed from 1 to 2, 3, or 4 between t-10 and t0, FALSE otherwise), "defor_prob" (predicted deforestation probability based on logistic regression, between 1 an 0), "risk" ("low" if defor_prob < 0.01, "high" otherwise), "acd10", "acd0", "c_loss" (annual carbon loss rate). Each row is a baseline pixel.

2. #project_defor_prob: predicted project deforestation probability
It is currently a list of N elements, saved as an RDS file that can be read as an R object.
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
This is a csv file containing N rows (one for each project)
 [1] "X"        "mean_100" "ci_100"   "mean_75"  "ci_75"    "mean_50" 
 [7] "ci_50"    "mean_25"  "ci_25"    "project"
