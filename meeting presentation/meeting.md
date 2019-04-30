Comparison of batch effect correction methods
========================================================
author: Guillaume Heger
date: May 3rd, 2019
autosize: true

Datasets
========================================================
- E-MTAB-513
- E-MTAB-2836
- E-MTAB-3716
- E-MTAB-4344
- E-MTAB-3871

RNASeq counts matrix from the Human transcriptome

Used methods
========================================================
- Mean-centering (pamr package)
- ComBat (with and without variables of interest) (sva package)
- Ratio-based methods (A & G) (bapred package)
- RUVg (RUVSeq package)

Slide With Code
========================================================


```r
summary(cars)
```

```
     speed           dist       
 Min.   : 4.0   Min.   :  2.00  
 1st Qu.:12.0   1st Qu.: 26.00  
 Median :15.0   Median : 36.00  
 Mean   :15.4   Mean   : 42.98  
 3rd Qu.:19.0   3rd Qu.: 56.00  
 Max.   :25.0   Max.   :120.00  
```

Slide With Plot
========================================================

![plot of chunk unnamed-chunk-2](meeting-figure/unnamed-chunk-2-1.png)
