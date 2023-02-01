# Antimicrobial resistance in patients with COVID-19: a systematic review and meta-analysis

## Purpose of this repository

This repository contains all data and code necessary to reproduce the analyses presented in the manuscript "[Antimicrobial resistance in patients with COVID-19: a systematic review and meta-analysis](https://doi.org/10.1016/S2666-5247(22)00355-X)" by Langford et al. (doi: [10.1016/j.cmi.2022.12.006](https://doi.org/10.1016/S2666-5247(22)00355-X)). **This includes every table and figure except for figure 1, supplementary table 1 and supplementary table 3. Supplementary table 1 is a subset of the dataset found in the file `data.csv`. Supplementary table 3 is also calculated from this dataset.** Note that the included scripts produces additional figures beyond those directly presented in the manuscript.

For more details on the methodological approach used in this meta-analysis, please see the manuscript [Meta-analysis of Proportions Using Generalized Linear Mixed Models](https://doi.org/10.1097/EDE.0000000000001232) by Lin & Chu (2020). Sample code is also provided.

## Requirements

All code is written in the programming language [R](https://www.r-project.org/). The easiest way to run it is to use the [RStudio](https://rstudio.com/) IDE. An `.Rproj` file is included with this repository for ease of use with RStudio. The scripts should run with any modern version of R.

The R packages required to reproduce the tables and figures at listed at the top of their respective scripts. They must be installed using `install.packages` or similar functionality within RStudio prior to running the script.

## Reproducing tables and figures

Run `analysis.R` to create the output tables and figures.

## Tables

- [Table 1](tables/summary_table_6.csv)
- [Table 2](tables/summary_table_4.csv)
- [Supplementary table 2](tables/summary_table_2.csv)
- [Supplementary table 4 (row 1)](figures/1_plot_1_coinfection.png), [Supplementary table 4 (row 2)](figures/1_plot_1_secondary_infection.png), [Supplementary table 4 (row 3)](figures/1_plot_1_bacterial_infection_unspecified.png), [Supplementary table 4 (rows 4–6)](figures/3_plot_8_Resistant_organisms.png)

## Figures

- [Figure 2A](figures/1_plot_1_coinfection.png), [Figure 2B](figures/1_plot_1_secondary_infection.png), [Figure 2C](figures/1_plot_1_bacterial_infection_unspecified.png)
- [Supplementary figure 1A](figures/3_plot_1_Resistant_pts_total.png), [Supplementary figure 1B](figures/3_plot_1_Resistant_organisms.png)
- [Supplementary figure 2A](figures/3_plot_6_Resistant_pts_total.png), [Supplementary figure 2B](figures/3_plot_6_Resistant_organisms.png)
- [Supplementary figure 3A](figures/3_scatterplots_unadjusted_Resistant_pts_total.png), [Supplementary figure 3B](figures/3_scatterplots_unadjusted_Resistant_organisms.png)