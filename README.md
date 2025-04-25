# EZbakRsuite_paper_code
Scripts for reproducing figures from EZbakR-suite paper.

Data is available on Zenodo (DOI: [10.5281/zenodo.13929898](https://zenodo.org/records/13946128)).

All scripts are written in R.

## Package versions used for preprint figures:

* R 4.3.1
* EZbakR 0.0.0.9000 ([Preprint release on Github](https://github.com/isaacvock/EZbakR/tree/v0.0.0.9000))
* dplyr 1.1.4
* ggplot 3.4.4
* rtracklayer 1.62.0
* GenomicFeatures 1.54.1
* MASS 7.3.60
* data.table 1.14.10
* tidyr 1.3.0
* readr 2.1.5
* bakR 1.0.1
* readxl 1.4.3
* arrow 14.0.0.2


## Full session info when making preprint figures:

```
R version 4.3.1 (2023-06-16 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 22631)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.utf8 
[2] LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8
[4] LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices
[5] utils     datasets  methods   base     

other attached packages:
 [1] bakR_1.0.1            
 [2] readr_2.1.5           
 [3] tidyr_1.3.0           
 [4] data.table_1.14.10    
 [5] arrow_14.0.0.2        
 [6] readxl_1.4.3          
 [7] GenomicFeatures_1.54.1
 [8] AnnotationDbi_1.64.1  
 [9] Biobase_2.62.0        
[10] rtracklayer_1.62.0    
[11] GenomicRanges_1.54.1  
[12] GenomeInfoDb_1.38.5   
[13] IRanges_2.36.0        
[14] S4Vectors_0.40.1      
[15] BiocGenerics_0.48.1   
[16] MASS_7.3-60           
[17] EZbakR_0.0.0.9000     
[18] dplyr_1.1.4           
[19] ggplot2_3.4.4         

loaded via a namespace (and not attached):
 [1] DBI_1.2.1                  
 [2] bitops_1.0-7               
 [3] gridExtra_2.3              
 [4] inline_0.3.19              
 [5] biomaRt_2.58.0             
 [6] rlang_1.1.2                
 [7] magrittr_2.0.3             
 [8] matrixStats_1.1.0          
 [9] compiler_4.3.1             
[10] RSQLite_2.3.5              
[11] loo_2.6.0                  
[12] png_0.1-8                  
[13] vctrs_0.6.4                
[14] stringr_1.5.1              
[15] pkgconfig_2.0.3            
[16] crayon_1.5.2               
[17] fastmap_1.1.1              
[18] dbplyr_2.4.0               
[19] XVector_0.42.0             
[20] utf8_1.2.4                 
[21] Rsamtools_2.18.0           
[22] tzdb_0.4.0                 
[23] purrr_1.0.2                
[24] bit_4.0.5                  
[25] zlibbioc_1.48.0            
[26] cachem_1.0.8               
[27] jsonlite_1.8.8             
[28] progress_1.2.3             
[29] blob_1.2.4                 
[30] DelayedArray_0.28.0        
[31] BiocParallel_1.36.0        
[32] parallel_4.3.1             
[33] prettyunits_1.2.0          
[34] R6_2.5.1                   
[35] StanHeaders_2.32.5         
[36] stringi_1.8.1              
[37] cellranger_1.1.0           
[38] Rcpp_1.0.12                
[39] assertthat_0.2.1           
[40] rstan_2.32.5               
[41] SummarizedExperiment_1.32.0
[42] Matrix_1.6-3               
[43] tidyselect_1.2.0           
[44] rstudioapi_0.15.0          
[45] abind_1.4-5                
[46] yaml_2.3.7                 
[47] codetools_0.2-19           
[48] curl_5.2.0                 
[49] pkgbuild_1.4.3             
[50] lattice_0.22-5             
[51] tibble_3.2.1               
[52] withr_3.0.0                
[53] KEGGREST_1.42.0            
[54] RcppParallel_5.1.7         
[55] BiocFileCache_2.10.1       
[56] xml2_1.3.6                 
[57] Biostrings_2.70.1          
[58] pillar_1.9.0               
[59] filelock_1.0.3             
[60] MatrixGenerics_1.14.0      
[61] generics_0.1.3             
[62] RCurl_1.98-1.13            
[63] hms_1.1.3                  
[64] munsell_0.5.0              
[65] scales_1.3.0               
[66] glue_1.6.2                 
[67] tools_4.3.1                
[68] BiocIO_1.12.0              
[69] GenomicAlignments_1.38.0   
[70] XML_3.99-0.15              
[71] grid_4.3.1                 
[72] QuickJSR_1.1.0             
[73] colorspace_2.1-0           
[74] GenomeInfoDbData_1.2.11    
[75] restfulr_0.0.15            
[76] cli_3.6.1                  
[77] rappdirs_0.3.3             
[78] fansi_1.0.5                
[79] S4Arrays_1.2.0             
[80] V8_4.4.1                   
[81] gtable_0.3.4               
[82] digest_0.6.34              
[83] SparseArray_1.2.2          
[84] rjson_0.2.21               
[85] memoise_2.0.1              
[86] lifecycle_1.0.4            
[87] httr_1.4.7                 
[88] bit64_4.0.5    
```


