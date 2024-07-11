This R script receives the output of the TMF implementation (https://github.com/quantifyearth/tmf-implementation),
and calculates the baseline deforestation and carbon loss. The baseline is currently based on the set M (matches.parquet), and is subject to change in the future.

It requires the following input variables to read TMF implementation output and other data.
All variables are vectors containing one value for each project to be analysed:

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


Its outputs are:
1. project_var: a data frame of basic information of every project
2. additionality_estimates: a list of the datas frames containing observed additionality estimates for each project
3. baseline_list: a list of the data frames containing each baseline pixel, its predicted deforestation probability, and its C loss
4. forecast_summ: a data frame of additionality forecast based on bootstrapped baseline C loss estimates and its confidence interval under different assumptions of project effectiveness
