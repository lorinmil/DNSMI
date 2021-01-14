# DNSMI
Identifying supervised multi-omic networks by DNSMI


# PRELIMINARIES: 
To perform this method, you need to have two datasets with the same number of rows and a univariate outcome. You may also have covariates to additionally adjust for. The method will return a set of variables from each dataset that both relate to each other AND relate to the outcome. This method is based off of Generalized Linear Models, so the distributions may be set if a different link should be used (defaults to guassian, but you may also specify "binomial" for logistic regression). Note that there is a directionality to the DNSMI method, so Dataset G -> Dataset X -> outcome y. For more details on the algorithm, please visit 

    F. Zhang, J. C. Miecznikowski, and D. L. Tritchler. Identification of Supervised and Sparse Functional Genomic Pathways. Statistical Applications in Genetics and Molecular       Biology, 19(1), 2020.


# STEP 1: 
Load in supporting R packages

    library(PMA)
    library(stringr)
    library(MASS)
    library(ROCR)
    
    
# STEP 2: 
Load in all of the supporting R functions from the files within this repository


# STEP 3: 
Now you may run the DNSMI method using the R function InstabilityMax_naive. An example is shown below:

    beta <- 0.03	# This is the instability threshold (specifies how sparse the solution is)
    repeatsN <- 100 #This is how many times to subsample when callibrating stability of the penalty thresholds
    level.precision <- 3 #This is how many decimal points for convergence when callibrating stability (to get it as close to beta as possible)
    max.iter <- 1500 #This is the maximum number of iterations to go out when callibrating the stability threshold
    k <- 1				# Do not change this number, keep it as 1
    set.seed(2255)
    sparseout <- InstabilityMax_naive(G=Dataset1, X=Dataset2, Y=outcomeVar,
                                      covariatesY = DatasetWithCovariates,
                                      distG="gaussian", distX="gaussian", distY="binomial",
                                      beta = beta, repeatsN = repeatsN, level.precision = level.precision, max.iter = max.iter)
                
                
# STEP 4: 
Extract the ideitified network features and other information

    sparseout$w       #The NSM matrix using all of the data. High values correspond to likely to be apart of the network.    
    sparseout$U_spas  #Indicators for each Dataset1 variable for whether it is part of the network (0=no, otherwise=yes)    
    sparseout$V_spas  #Indicators for each Dataset2 variable for whether it is part of the network (0=no, otherwise=yes)    
    sparseout$opt_c1  #Tuned sparsity threshold for Dataset1    
    sparseout$opt_c2  #Tuned sparsity threshold for Dataset2
