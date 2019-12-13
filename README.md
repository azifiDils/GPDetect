# Description of the package

GPDetect is meant to showcase the methodology used in the paper by Ramzan et al.

The function `GPDetect()` calculates cubic smoothing splines on the test statistic values to obtain a smoothed curve. The inflection points of the curve are determined and the region between two consecutive inflection points having a downward concave form is regarded as a peak. The maximum test statistic value within a peak after smoothing is recorded as the height of the peak. To determine the optimal level of smoothing, cross validation is implemented. The GWAS results must be obtained beforehand and supplied to the function. Any test statictic that represents the strength of association between a SNP and the phenotype can be used.

# Installation

To get the package:

```R
devtools::install_github("azifidils/GPDetect")
```

# Usage

To check out the results of GPDetect on a simulated data set:

```R
library("GPDetect")
#run the example given in the package
example("GPDetect", ask = FALSE)
```

The results will be saved in the data frame 'result' and scatter plots with the applied smoothed spline for each chromosome will be saved in the current directory in the folder 'results'.

Each row in the data frame defines a peak with its start (Pos.right) and end position (Pos.right), their SNP ids (InitialSNP/lastSNP), the total number of SNPS within the peak (NSNP) as well as the height of this peak:

| Peak  | NSNP  | Pos.left  | Pos.right | InitialSNP    | lastSNP   |  Height   | Chr   |
|:-----:|:-----:|:---------:|:---------:|:-------------:|:---------:|:---------:|:-----:|
| 1     | 171   | 2.73804   | 4.55606   | M274          | M444      | 1.0982116 | 1     |
| 2     | 169   | 6.44007   | 8.15209   | M623          | M791      | 1.3910353 | 1     |
| 3     | 122   | 9.21610   | 10.37611  | M910          | M1031     | 2.2179283 | 1     |

For further information refer to the documentation of the `GPDetect()`

# Prerequisites

GPDetect requires single-SNP based GWAS results as input. The format may look like this:

| SNP       | Chromosome    | Position  | Pvalue    | Qvalue | N    | NullLogLike   | AltLogLike    | SNPWeight | SNPWeightSE   | OddsRatio | WaldStat  | NullLogDelta  | NullGeneticVar    | NullResidualVar   | NullBias  |
|:---------:|:-------------:|:---------:|:---------:|:------:|:----:|:-------------:|:-------------:|:---------:|:-------------:|:---------:|:---------:| :------------:|:-----------------:|:-----------------:|:---------:| 
| M41855    | 3             | 15.04016  | 4.61784507639991E-43 | 9.23569015279983E-38 | 5000 | -18084.0632491666 | -17991.0043438445 | -2.0212977536 | 0.1455316896 | 0.0446536889 | 192.9058702374 | 0.1113069532 | 42.0518232423 | 47.0029185959 | 27.9081213446 |
|M41863 | 3 | 15.10616 | 1.64702166601606E-25 | 1.29677373810485E-20 | 5000 | -18084.0632491666 | -18031.0279927869 | -1.5211844846 | 0.1449217496 | 0.0844246696 | 110.1784446639 | 0.1113069532 | 42.0518232423 | 47.0029185959 | 27.9081213446 |
M3943 | 1 | 39.81441 | 1.94516060715727E-25 | 1.29677373810485E-20 | 5000 | -18084.0632491666 | -18031.0901319964 | -1.5621555522 | 0.1490531951 | 0.0998462992 | 109.841378134 | 0.1113069532 | 42.0518232423 | 47.0029185959 | 27.9081213446 |

The columns SNP, Chromosome, Position and a test statistic value (e.g. WaldStat) are mandatory for the function to work.
