# SKAT R-package

SKAT is an R-package for gene/region-based multiple variant tests. It has functions for Burden test, SKAT, and SKAT-O with adjusting for covariates and kinship. 

## Main Functions and key references
Following are main functions and key references. For details, please refer the package [manual](https://cran.r-project.org/web/packages/SKAT/SKAT.pdf) and [vignettes](https://cran.r-project.org/web/packages/SKAT/vignettes/SKAT.pdf).

1. SKAT function: Burden test, SKAT, and SKAT-O 
    * Lee, S., Emond, M.J., ..., and Lin, X. (2012). Optimal unified approach for rare variant association testing with application to small sample case-control whole-exome sequencing studies. *AJHG*, 91, 224-237.
    * Lee, S., Wu, M. and Lin, X. (2012). Optimal tests for rare variant effects in sequencing association studies. *Biostatistics*, 13, 762-775 
    * Wu, M., Lee, S., Cai, T., Li, Y., Boehnke, M. and Lin, X. (2011). Rare Variant Association Testing for Sequencing Data Using the Sequence Kernel Association Test (SKAT). *AJHG*, 89, 82-93. 

2. SKATBinary function: Burden test, SKAT, and SKAT-O with efficient resampling for binary traits

    * Lee, S., Fuchsberger, C., Kim, S., Scott, L. (2016) An efficient resampling method for calibrating single and gene-based rare variant association analysis in case-control studies, *Biostatistics*, 17, 1-15.

3. SKAT_CommonRare function: joint test for common and rare variants

    * Ionita-Laza, I., Lee, S., Makarov, V., Buxbaum, J. Lin, X. (2013). Sequence kernel association tests for the combined effect of rare and common variants. *AJHG*, 92, 841-853. 

4. SKAT_ChrX function: X-chromosome test

    * Ma, C., Boehnke, M., Lee, S. and the GoT2D investigators (2015) Evaluating the calibration and power of three gene-based association tests for the X chromosome, *Genetic Epidemiology*, 39, 499-508.
 
5. SKAT_NULL_emmaX function: Kinship adjustment 

6. SSD functions: plink binary file related functions

7. Power_Continuous, ...: power calculation functions 

## Link
* SKAT CRAN: [Link](https://cran.r-project.org/web/packages/SKAT/index.html)
* SKAT google group: [Link](https://groups.google.com/forum/#!forum/skat_slee)
* Example dataset: [Download](https://github.com/leeshawn/SKAT/blob/master/vignettes/Example.zip)  

