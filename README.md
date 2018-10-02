# Joint Indirect Standardization when Only Marginal Distributions are Observed in the Index Population

# Author Contributions Checklist Form

## Data

### Abstract 
Dataset describing one binary outcome, whether an abdominal computed tomography scan has
gone over a fixed radiation dose pre-determined to be high. This outcome is influenced by four
covariates, each of which is also categorical. These covariates are: patient’s age, patient’s
gender, patient’s abdominal size, and whether the scan was single phase or multi-phase.

### Availability 
To protect patient confidentiality, the hospitals providing the example data used in this paper
have not given permission for the data to be made publicly available. For public use, an
alternative set of data summaries has been generated to run with provided code, but this
dataset should be not expected to produce identical results to the example data used in this
paper.

### Description 
Data contained in file “RSCdataObject.Rdata”. This is an Rdata file containing the R object
“full.out”, which is a list containing five objects. These objects are as follows:

```
1) An array with a number of rows equal to the number of covariate combinations, and
columns equal to the number of reference hospitals. Each column is the observed joint
covariate probability vector in a reference hospital.
2) A vector containing the observed prevalence of the outcome in each reference hospital.
3) An array, with each row being a covariate combination.
4) A vector containing the expected prevalence of the outcome for each covariate
combination.
5) A vector containing the sample size of each hospital.
```
## Code

### Abstract 
A set of files with functions that can be used to conduct the methods described in the
manuscript, as well as two additional files which apply the methods to assess performance, as
also described in the paper. One of such files was used to apply the example data, another to
run simulations. Code written and run in R version 3.4.3.

### Description 
All files are delivered within a single zip file (RSCfiles.zip). Within this zip file is contained the
following files:


1. RSCfiles/example_data.R: Contains code used to apply our algorithm to example data.
2. RSCfiles/simulated.R: Contains code used to apply our algorithm to simulations.
3. RSCfiles/Package: Contains R files of all functions used to perform the methods
    describes in this paper.
4. RSCfiles/Package/RSCcond.R: Back-end function which produces a matrix of all unique
    combinations of levels of multiple categorical variables.
5. RSCfiles/Package/RSCform.R: Back-end function which produces the matrices used in
    equations used in raking.
6. RSCfiles/Package/RSCpoint.R: Rakes one index hospital against one reference hospital.
7. RSCfiles/Package/RSCconcentrationBase.R: Given a reference dataset, predicts the
    coverage rate given Dirichlet concentration parameter.
8. RSCfiles/Package/RSCconcentration.R: Given a reference dataset, optimizes Dirichlet
    concentration parameter with respect to coverage rate.
9. RSCfiles/Package/RSC.R: Profiles a single index hospital against a reference dataset.
10. RSCfiles/forPublic/: To protect patient confidentiality, the hospitals providing the example
    data used in this paper have not given permission for the data to be made publicly
    available. We have thus instead provided a simulated “fake” dataset. The results from
    applying our proposed methods to this dataset produces very different SIR estimates,
    due to very different index and reference hospitals in terms of outcome prevalence,
    covariate distribution, and sample size. However, the performance of our methods, in
    terms of coverage rate, prediction rate, and bias remain just as good as when using real
    data.
11. RSCfiles/forPublic/RSCdataObject.Rdata: Dataset in the correct format to apply to
    “RSCfiles/Package/RSC.R”.
12. RSCfiles/forPublic/RSCrakedObject.Rdata: All hospitals raked against all other hospitals.
13. RSCfiles/forPublic/RSCshapeObject.Rdata: Correct Dirichlet concentration parameter
    for the example data application, given our hospitals randomly-selected as referential.
14. RSCfiles/forPublic/RSCexampleObject.Rdata: Primary output of example data
    application.
15. RSCfiles/forPublic/RSCsimulationObject.Rdata: Primary output of simulation study.
16. RSCfiles/make_fake_data.R: This file was used to generate the “fake” dataset contained
    in the “forPublic" folder. As the generation of the “forPublic” fake dataset uses some
    information from the real, unpublished, dataset, this file will not run for the public
    audience, but has been included for completeness.


These files can also be found on my personal GitHub account at
[http://www.github.com/dryifeiwang.](http://www.github.com/dryifeiwang.)

## Instructions for Use

### Reproducibility
Extract the entirety of RSCfiles.zip into the working directory. Outcomes akin to those shown in
the paper can be produced by running the entireties of “example_data.R” and “simulated.R”.

Running “example_data.R” will produce
1) An R object called “test.outputs” in the global environment, which is an array used to
compute the numerical summaries of the Example Data section. Said numerical
summaries will be presented in the R console immediately after creation of
“test.outputs”.
2) Image file “plotSMR.pdf” in the working directory, which is Figure 1.

Running “simulated.R” will produce
1) An R object called “bunch.of.outputs”, which is a list of arrays, each array being
analogous to “test.outputs” object described in “example_data.R”, but under different
simulated hospitals. This object is used to compute the numerical summaries of the
Simulation section, presented in the R console immediately after creation of
“bunch.of.outputs”.
2) Image files “sumplotSIRsim1.png” and “sumplotSIRsim2.png” in the working directory,
which are Figures 2 and 3.

Since the data available for public use is not the same as the one used in the manuscript, note
that these outcomes will not show identical numerical results to the paper.


