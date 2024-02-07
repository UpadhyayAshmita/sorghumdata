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

-   calculate heritability for wavelength for location EF, MW,and joint for 350-2500

-   calculate BLUEs for Loc EF,MW and joint for 350-2500 wave number

-   calculate heritability for trait of interest(Narea and SLA)\

    |       | EF    | MW    | joint |
    |-------|-------|-------|-------|
    | SLA   | 0.345 | 0.33  | 0.58  |
    | Narea | 0.355 | 0.587 | 0.59  |

    : Table:1- Heritability of trait of inetrest

-   calculate BLUEs for Narea and SLA for loc EF, MW, and Joint and included weight and this is the first stage model

-   filter name2 based on corrected name from name_west file

    3.  NA Imputation script

```{=html}
<!-- -->
```
-   We have NA for several wavelength so we pass simpler model for those wavelength and impute the NA

-   Likewise, we impute NA with the simpler heritability model for wavelength and have complete output from first stage of analysis

    4.  Relationship matrix script

-   We calculate kinship matrix from the genotype information which was in vcf file format

-   We load the wavelength blues for three different location "wavebluesEF", "wavebluesMW", and "wavejoint" file after imputing NA in a first stage output and then correct it by left joining corrected name with name2 from wavelength blues to filter individual that are also present in genotypic data.

-   With the corrected waveblues we calculate relationship matrix which are: NIRS, Wholewave and highly heritable matrix.

-   NIRS relationship matrix was created by filtering wavenumber that lies in nirs region (wave number 860 to 1660), which was further created a square matrix calling it NIRS matrix using function calculate_relationship

-   Wholewave relationship matrix was created using whole 350-2500 wavenumber into a square matrix using function similar as above

-   We calculated heritability in first stage, so we plotted the ggplot for all wavenumber and then we set threshold of 0.4 to filter the wavelength and then further created relationship matrix with the wavelength that has heritability greater than 0.4 calling the square matrix highly heritable(highh2)matrix.

-   We calculated NIRS, Wholewave (WW) and highh2 relationship matrix for three different location "EF", "MW" and "Joint".

    5.  Second stage script:

-   From the blues data we obtained from first stage model, we filtered them using selected line because selected line were the line that had genotype information so we filter Narea blues and SLA blues for three different location

-   We fit second stage model where we fit GBLUP model using kinship matrix for three different location for Narea and SLA and here we performed 5- fold cross-validation to assess the model performance where 20 repetition was performed and each repetition has 5 fold.

-    We obtained Genomic estimated breeding value(gebv) from the GBLUP model and the accuracy i.e., correlation between the predicted value and observed value in a validation set.

-   After, GBLUP we fitted NIRS, WW and highh2 model where we used blues from first stage for Narea and SLA and here we used relationship matrix we generated from waveblues. Then we obtained accuracy and gebv from the second stage NIRS, WW and highh2 model.

-   At the end of this analysis ,we have total 4 different types of model i.e., GBLUP, NIRS, Wholewave, and highlyheritable model.

    6.  Combination scheme

        Relationship matrix script:

-   We then create different relationship matrix i.e., different combined matrix using kinship + NIRS, highh2, wholewave

-   We first remove 10%, 25% and 50% individual from kinship matrix and then calculated kinship matrix called G10 which has 10% less individual than kinship and G25 and G50 respectively.

-   We generated a combined relationship matrix called Gnirs where 10% removed individual from kinship were imputed using NIRS matrix and then 25% and 50% removed individual info was imputed using NIRS matrix

-   Likewise, we created GWW where 10% removed individual from kinship were imputed using wholewave relationship matrix and then 25% and 50% removed individual info was imputed using wholewave matrix.

-   Finally we generated Gh2 relationship matrix where 10% , 25% and 50% removed individual from kinship were imputed using highly heritable relationship matrix information.

-   By this time we have four new kind of relationship matrix. The first one is called GBLUP_reduced relationship matrix, and three are combined relationship matrix called Gnirs, GWW, and Gh2 for 10%, 25% and 50% scheme for three different location

    7.  Second stage script

-    Using kinship matrix that had 10%, 25% and 50% individual removed we fit second stage model called GBLUP_reduced model (10%, 25% and 50% ) for three different location using blues from first stage for Narea and SLA where 5 fold cross validation was performed with 20 repetition.

-   We fit second stage model called Gh2, Gnirs and GWW model using combined relationship matrix for 10%, 25% and 50% scheme for three different location

-   By this time we have new combined model called GBLUP reduced, Gnirs, GWW, and Gh2 for 10%, 25% and 50% scheme for three different location, their accuracy and gebv for all individual.
