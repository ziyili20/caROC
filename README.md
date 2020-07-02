# Evaluate continuous biomarkers with caROC

`caROC` is an R package devoted to the assessment of continuous biomarkers. The metrics considered include **specificity at contolled sensitivity level**, **sensitivity at controlled specificity level**, **receiver operating characteristic** (ROC) curve. We allow both categorical and continuous covariates to be adjusted for in computing these metrics.

## Installation and quick start

### Install caROC

Before installing caROC, be sure you have installed the `cequre` package, which you can download from [here](http://web1.sph.emory.edu/users/yhuang5/software/index.html). You can download the `cequre` package to your computer. If you use R studio, install it through Tools->Install Packages->Install from local directory. If you use terminal, navigate to the folder containing the cequre gz file, you can install cequre through command line

```
R CMD INSTALL cequre_1.3.tar.gz
```

After `cequre` package has been installed, install `caROC` through

```{r install, message=FALSE, warning=FALSE}
library(devtools)
install_github("ziyili20/caROC")
```

### How to get help for caROC

Any caROC questions should be posted
to the GitHub Issue section of caROC 
homepage at https://github.com/ziyili20/caROC/issues.

### Quick start on evaluating continuous biomarker with covariates adjusted

```{r quick_start, eval = FALSE}
## get specificity at controlled sensitivity levels 0.2, 0.8, 0.9
caROC(diseaseData,controlData,formula,
      control_sensitivity = c(0.2,0.8, 0.9),
      control_specificity = NULL)
      
## get covariate-adjusted ROC curve with curve-based monotonizing method
curveROC <- caROC(diseaseData,controlData,formula,
            mono_resp_method = "curve", 
            verbose = FALSE)
```

### Illustrating the usage of caROC in details

The tutorial is based on a simulation dataset:

```{r getdata, eval = FALSE}
### n1: number of cases
### n0: number of controls
n1 = n0 = 1000

## Z_D and Z_C are the covariates in the disease and control groups
Z_D <- rbinom(n1, size = 1, prob = 0.3)
Z_C <- rbinom(n0, size = 1, prob = 0.7)

Y_C_Z0 <- rnorm(n0, 0.1, 1)
Y_D_Z0 <- rnorm(n1, 1.1, 1)
Y_C_Z1 <- rnorm(n0, 0.2, 1)
Y_D_Z1 <- rnorm(n1, 0.9, 1)

## M0 and M1 are the outcome of interest (biomarker to be evaluated) in the control and disease groups
M0 <- Y_C_Z0 * (Z_C == 0) + Y_C_Z1 * (Z_C == 1)
M1 <- Y_D_Z0 * (Z_D == 0) + Y_D_Z1 * (Z_D == 1)

diseaseData <- data.frame(M = M1, Z = Z_D)
controlData <- data.frame(M = M0, Z = Z_C)

## we are interested in evaluating biomarker M while adjusting for covariate Z
formula = "M~Z"
```

## 1. Covariate-adjusted sensitivity at controlled specificity level (or the reverse)


One can easily compute covariate-adjusted specificity at controlled sensitivity levels by specifying `control_sensitivity` and leaving `control_specificity` NULL. 

`mono_resp_method` is to choose which monotonicity restoration method to use, "none", "mono" or "curve". `whichSE` is to choose how to compute standard error. It could be "boostrap" or "numerical", i.e. boostrap-based or sample-based SE. Try ?caROC to see more details of these arguments.

```{r controlspec}
caROC(diseaseData,controlData,formula,
      control_sensitivity = c(0.2,0.8, 0.9),
      control_specificity = NULL,
      mono_resp_method = "curve",
      whichSE = "bootstrap",nbootstrap = 100,
      CI_alpha = 0.95, logit_CI = TRUE)
```

To compute covariate-adjusted sensitivity at controlled specificity levels by specifying `control_specificity` and leaving `control_sensitivity` NULL. 

```{r controlsens}
caROC(diseaseData,controlData,formula,
      control_sensitivity = NULL,
      control_specificity = c(0.7,0.8, 0.9),
      mono_resp_method = "none",
      whichSE = "numerical",nbootstrap = 100,
      CI_alpha = 0.95, logit_CI = TRUE)
```

## 2. Covariate-adjusted ROC curve

Obtaining the covariate-adjusted ROC curve with sensitivity controlled through the whole spectrum is very easy. You can choose different ways of restoring monotonicity, or no restoration when constructing ROC through argument `mono_resp_method`. It could be "none" (no monotonicity restoration), "mono" (method based on Huang, 2017), "curve" (curve-based monotonicity restoration). 

```{r ROC}
### ROC with curve-based monotonicity restoration
curveROC <- caROC(diseaseData,controlData,formula,
                 mono_resp_method = "curve", 
                 verbose = FALSE)
```

Plot the ROC curves:

```{r plotROC}
plot_caROC(curveROC)
```

Construct confidence-band for the ROC curve:

```{r ROC}
curveROC_CB <- caROC_CB(diseaseData,controlData,
						formula, 
						mono_resp_method = "curve",
						CB_alpha = 0.95,
						nbin = 100,verbose = FALSE)
 ```   
 
Plot the confidence band:

```{r plotROCband}
plot_caROC_CB(curveROC_CB, add = FALSE, lty = 2, col = "blue")                   
```

or plot the ROC and confidence band on the same plot:

```{r plotROCband}
plot_caROC(curveROC)
plot_caROC_CB(curveROC_CB, add = TRUE, lty = 2, col = "blue")
```

## 3. Threshold at controlled sensitivity/specificity for given covariate values

In clinical setting, it is useful to know the specific thresholds of biomarkers at controlled sensitivity or specificity level for given covariate values.

```{r treshold}
### this is the given covariates of interest
new_covariates <- data.frame(Z = controlData[2:11, 2])

### controlling sensitivity levels
caThreshold(formula, new_covariates,
            diseaseData = diseaseData,
            controlData = NULL,
            control_sensitivity = c(0.7,0.8,0.9),
            control_specificity = NULL)
            
### controlling specificity levels
caThreshold(formula,new_covariates,
            diseaseData = NULL,
            controlData = controlData,
            control_sensitivity = NULL,
            control_specificity = c(0.7,0.8,0.9))
```




