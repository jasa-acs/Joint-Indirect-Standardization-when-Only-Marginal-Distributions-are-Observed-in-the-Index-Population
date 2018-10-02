Code
Abstract
A set of files with functions that can be used to conduct the methods described in the manuscript, as well as two additional files which apply the methods to assess performance, as also described in the paper. One of such files was used to apply the example data, another to run simulations. Code written and run in R version 3.4.3.

Description
All files are delivered within a single zip file (RSCfiles.zip). Within this zip file is contained the following files:

1.	RSCfiles/example_data.R: Contains code used to apply our algorithm to example data.

2.	RSCfiles/simulated.R: Contains code used to apply our algorithm to simulations.

3.	RSCfiles/Package: Contains R files of all functions used to perform the methods describes in this paper.

4.	RSCfiles/Package/RSCcond.R: Back-end function which produces a matrix of all unique combinations of levels of multiple categorical variables.

5.	RSCfiles/Package/RSCform.R: Back-end function which produces the matrices used in equations used in raking.

6.	RSCfiles/Package/RSCpoint.R: Rakes one index hospital against one reference hospital.

7.	RSCfiles/Package/RSCconcentrationBase.R: Given a reference dataset, predicts the coverage rate given Dirichlet concentration parameter.

8.	RSCfiles/Package/RSCconcentration.R: Given a reference dataset, optimizes Dirichlet concentration parameter with respect to coverage rate.

9.	RSCfiles/Package/RSC.R: Profiles a single index hospital against a reference dataset.

10.	RSCfiles/forPublic/: To protect patient confidentiality, the hospitals providing the example data used in this paper have not given permission for the data to be made publicly available. We have thus instead provided a simulated ÅgfakeÅh dataset. The results from applying our proposed methods to this dataset produces very different SIR estimates, due to very different index and reference hospitals in terms of outcome prevalence, covariate distribution, and sample size. However, the performance of our methods, in terms of coverage rate, prediction rate, and bias remain just as good as when using real data.

11.	RSCfiles/forPublic/RSCdataObject.Rdata: Dataset in the correct format to apply to ÅgRSCfiles/Package/RSC.RÅh.

12.	RSCfiles/forPublic/RSCrakedObject.Rdata: All hospitals raked against all other hospitals.

13.	RSCfiles/forPublic/RSCshapeObject.Rdata: Correct Dirichlet concentration parameter for the example data application, given our hospitals randomly-selected as referential.

14.	RSCfiles/forPublic/RSCexampleObject.Rdata: Primary output of example data application.

15.	RSCfiles/forPublic/RSCsimulationObject.Rdata: Primary output of simulation study.

16.	RSCfiles/make_fake_data.R: This file was used to generate the ÅgfakeÅh dataset contained in the ÅgforPublic" folder. As the generation of the ÅgforPublicÅh fake dataset uses some information from the real, unpublished, dataset, this file does not run for the public audience, but has been included for completeness.

These files can also be found on my personal GitHub account at http://www.github.com/dryifeiwang.
