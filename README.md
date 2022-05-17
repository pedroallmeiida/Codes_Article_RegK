# Codes: K-Bessel Regression Model for Speckled Data

## Project description

The objective of this project is to build functions and algorithms for the K-Bessel regression model proposed in the Article ```K-Bessel Regression Model for Speckled Data```. Several functions have been created for using the model. We applied the model in different parts of the Image of São Francisco (Ocean, Forest and Urban). In the file ```Code_regK.R``` we provide the code needed to run the application. 

## Datasets

```data_ocean.csv```: dataset of the ocean region of the San Francisco image

```data_forest.csv```: dataset of the forest region of the San Francisco image

```data_urban.csv```: dataset of the urban region of the San Francisco image


## Functions

```meijer.py -- formato .py```: function created to calculate Meijer function in python

``` func_MEIJER.R -- formato .R ```: function used to call the created function ``meijer.py`` to .R format using the reticulate library

``` Functions_RegK.R -- formato .R ```: function with the functions created from the model: Generator, likelihood functions and optimization algorithms, etc.

``` Functions_regressaoK.R -- formato .R ```: function with the functions created from the model: Generator, likelihood functions and optimization algorithms, etc.

## How to run
#### Step (1): Installation of dependencies

```
install.packages(dplyr)
install.packages(parallel)
install.packages(pbmcapply)
install.packages(magrittr)
install.packages(ggplot2)
install.packages(reshape2)
install.packages(robustbase)
install.packages(maxLik)
install.packages(reticulate)
install.packages(nortest)
```  
 
#### Step (2): load the functions
 
``` 
source(func_MEIJER.R) 
source(Functions_RegK.R) 
source(Functions_regressaoK.R)  
```

#### Step (3): File for execution
  
 ``Code_regK.R``: The ``glm.K()`` function defined in the file is responsible for executing the model. It depends on the response variable, covariate, optimization method and number of looks;

Below is an example of the output of the function ``glm.K(Y=variable_response, X=covariate, method = "NM", L=4 )``:

```
$estimate
    alpha       b11       b12
12.916549 -4.283202 44.836116 

$residual
 [1] -0.469932802  0.388750276  1.583027762  1.113346938  1.784391669
 [6]  0.372897335  0.857168118 -0.964570009  0.233094002  0.706256891
[11] -0.075567417  1.212633450  0.006968036  2.611626940  0.271688085
[16] -0.530191225  0.421818410 -0.630984397 -0.674811873 -0.340070297
[21] -0.142433453 -0.803176809  1.226112771 -0.970526529 -0.530511305
[26]  0.202838266  0.899918175  1.083397967 -1.430884767 -0.391027436
[31]  0.637085550 -0.118946257 -1.740669808 -1.288912794 -0.778620700
[36]  0.893494276  0.068170879 -0.573987789 -1.264469829 -0.622365645
[41]  0.300039319 -0.746774460  1.200334674  1.188845158  0.635669083
[46] -0.305209917 -0.733972635 -2.595207229 -1.103807635

$AIC
-263.1548

$BIC
-261.3712

$fitted.values
 [1] 0.03495813 0.03437094 0.02446610 0.03203000 0.02862450 0.03003447
 [7] 0.02523675 0.04563472 0.02940072 0.03162495 0.03346641 0.03216943
[13] 0.04006289 0.04692924 0.01680074 0.02091545 0.02822633 0.05006250
[19] 0.05210805 0.07178651 0.10507865 0.05158560 0.01692342 0.04544088
[25] 0.03475375 0.03431724 0.03413735 0.03846679 0.01636678 0.01647054
[31] 0.01998212 0.04226513 0.10412002 0.02984331 0.03877040 0.02268912
[37] 0.01570074 0.01532045 0.01571473 0.01682082 0.02493277 0.06441928
[43] 0.02726532 0.03484596 0.05505844 0.03401007 0.01917499 0.01507662
[49] 0.01690230

$r.squared
0.3715021

$mae
[1] 0.01625003
```

## Article link

.....

## MIT License
*Copyright (c) 2022 Pedro Almeida*
