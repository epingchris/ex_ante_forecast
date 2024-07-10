This R code receives the output of the TMF implementation (https://github.com/quantifyearth/tmf-implementation),
and calculates the baseline deforestation and carbon loss. The baseline is currently based on the set M (matches.parquet), and is subject to change in the future.

Its inputs are:
1. Sampled project pixel sets (set K, "directory_path/file_prefix_k.parquet")
2. Potentially matched pixel sets (set M, "directory_path/file_prefix_matches.parquet")
3. Matched pair pixel sets ("directory_path/pairs/file_prefix_xxx.parquet") - for observed additionality calculation only
4. Carbon density per LUC ("directory_path/file_prefix_carbon-density.csv")
5. Project polygon paths: a vector of the full paths of the project polygons (used to calculate project area)
6. Year of start of project: real or hypothetical, depending whether ongoing projects are being analysed or not
7. Country of project

Its outputs are:
1. project_var: a data frame of basic information of every project
2. additionality_estimates: a list of the datas frames containing observed additionality estimates for each project
3. baseline_list: a list of the data frames containing each baseline pixel, its predicted deforestation probability, and its C loss
4. forecast_summ: a data frame of additionality forecast based on bootstrapped baseline C loss estimates and its confidence interval under different assumptions of project effectiveness
