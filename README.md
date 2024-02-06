## General

This repository contains all the scripts and data used for the sorghum project.

## Content of Repository

-   Raw data used in the study from WEST sorghum diversity pannel

-   Relationship matrices created for different model

-   All output in a folder

-   R-script for the analysis for stage-1 and stage-2

-   figures are created through pre-processing script

-   Figures are stored in figures folder

# Methodology for the analysis based on script:

1.  Pre-processing script

-   Read the field design and wavelength data and filter year 16

-   exploratory data analysis on the design data.

2.  stage-1 script (phenotypic data)

-    calculate heritability for wavelength for location EF, MW,and joint for 350-2500

-   calculate BLUEs for Loc EF,MW and joint for 350-2500 wave number

-   calculate heritability for trait of interest(Narea and SLA)\

    |       | EF    | MW    | joint |
    |-------|-------|-------|-------|
    | SLA   | 0.345 | 0.33  | 0.58  |
    | Narea | 0.355 | 0.587 | 0.59  |

    : Table:1- Heritability of trait of inetrest

-   calculate BLUEs for Narea and SLA for loc EF, MW, and Joint and included weight and this is the first stage model

-   filter name2 based on corrected name from name_west file

3.  NA Imputation script:
