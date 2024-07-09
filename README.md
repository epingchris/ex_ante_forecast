This R code receives the output of the TMF implementation (https://github.com/quantifyearth/tmf-implementation),
and calculates the baseline deforestation and carbon loss. The baseline is currently based on the set M (matches.parquet), and is subject to change in the future.

Its inputs are:
1. Project directory paths: a vector of the full paths of directories containing the TMF implementation output of each project
2. Project polygon paths: a vector of the full paths of the project polygons (used to calculate project area)
3. Year of start of project: real or hypothetical, depending whether ongoing projects are being analysed or not

Its outputs are:

