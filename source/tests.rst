.. _tests:

Association Tests
==================================

Common Single Variant Association Test
----------------------------------

For both quantitative and binary trait analyses, a common variant test refers to a score test which has a Chi-squared distribution with 1 degree of freedom under no association hypothesis. The general form of the score test appears as follows:

    :math:`T={{S^2} \over var(S)}`

Where \\(T\\) is the test statistic following a Chi-squared distribution and \\(S\\) is the score. This test is used to perform a genetic association analysis between the phenotype \\(Y\\) and a single variant \\(G_j\\). For testing variant \\(j\\) given \\(n\\) indivduals and phenotype vector \\(Y\\) and genotype matrix \\(G\\),

    :math:`S_j=\sum_{i=1}^n(Y_i-E(Y_i))G_{ij}`

    :math:`E(Y_i)=Y_i - \hat{Y}`

\\(E(Y_i)\\) is estimated from a vector of fitted values \\(\\hat{Y}\\) which is dependent on the underlying distribution of \\(Y\\) (ex. case-control vs quantatitive). With no covariates, \\(\\hat{Y}=\\bar{Y}\\) which is the simple average of the observed phenotypes. 

Under a case-control setting with no covariates, the score is an indication of how often the tested genotype appears in one group over the other. When coded as \\(Y_i=1\\) for cases and \\(Y_i=0\\) for controls, and genotypes coded as {0,1,2} corresponding to the number of alleles a particular individual possesses. Given this framwork, cases with the allele of interest contribute positively to the overall score and controls contribute negatively. Therefore, the more a particular allele is associated with one group, the larger the magnitude of the score.

For genotypes hardcoded as {0,1,2}, the conventional variance formula is used to calculate \\(var(S_j)\\). To produce the test statistic \\(T_j\\), the square of the score \\(S_j\\) is normalized by the variance \\(var(S_j)\\) and a p-value is produced by evaluating \\(T_j\\) with respect to a Chi-squared distribution with 1 degree of freedom. In general, a large score and a small variance will result in a small (more significant) p-value.

In the vRVS methodology available in VikNGS, the genotype value \\(G_{ij}\\) is replaced with the expected genotype value calulated from the sequence read data \\(E(G_{ikj}\\mid D_{ikj})\\). When integrating data from an arbitrary number of cohorts, the variance is calculated for each group separately and summed together to produce \\(var(S_j)\\). The details of the derivation of \\(var(S_j)\\) are given in the Supplementary document of the VikNGS paper *VIKNGS: A C++ Variant Integration Kit for next generation sequencing association analysis*.

Rare Variant Association Test
----------------------------------

For joint variant analysis, the score statistics for \\(J\\) variants, \\(\\boldsymbol{S}=[S_1,...,S_J]\\). Please review the common variant section above to review the general structure of a score test. In VikNGS, multiple differen genetic association tests are available which are described below.

For the CAST- and SKAT-like tests, we recommend the use of permutation to calculate p-values. This involves shuffling the phenotype vector \\(Y\\) and recalculating the p-value many times for every variant. After iteratively calculating a set of p-values, the final p-value is calculated based on the number of values that are less than or equal to the value that was calculated for the unshuffled data set divided by the number of iterations plus 1. 

.. note::
    Using permutation, the smallest p-value obtainable is \\(1/{(# iterations + 1}\\). Since this method can be very computationally expensive, an an early stopping procedure is avaliable to terminate the calculation early if the p-value appears to be > 0.05. This uses the method designed by  Jiang and Salzman (2012 `🔗 <https://www.ncbi.nlm.nih.gov/pubmed/23843675>`).
    
    In VikNGS, these tests can be run by assuming the asymptoic distribution by setting the number of iterations to 1. Based on our testing, the results appear to be behaved but we offer no statistical guarantees.

Linear Test (CAST-like)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This test related to the CAST method described by Morgenthaler and Thilly (2007 `🔗 <https://www.ncbi.nlm.nih.gov/pubmed/17101154>`). In this  "Linear" refers to the fact that the score in this test is a linear combination of 


\item \textbf{When using true genotypes or genotype calls in conventional score test:}\\
``CAST-like" and ``"SKAT-like" refer to the CAST \citep{morgenthaler:2007} and SKAT with weights $w^{1/2}=1/[MAF(1-MAF)]^{1/2}$ \citep{wu:2011} respectively, when the true genotypes (in simulation part) and genotype calls are used.  The P-values can be obtain through permutation. 
\textit{\textbf{Without covariates}}, the $Y_i$s are permuted and a test statistics is calculated for each permuted data set where the user defines the number of iterations, e.g. 1,000. At the end of the iterations, P-value is calculated based on the number of values that are less than or equal to the value that was calculated for the initial data set divided by the number of iterations plus 1. \textit{\textbf{When covariates are added}}, we implement the bootstrap methodology defined in \citep{lin:2011} for for binary and quantitative trait analyses.
The user can also use asymptotic distribution to compute the P-values. In that case, the user should define the number of iterations as 1 in VikNGS (e.g. see Figure 2 in the paper ``VIKNGS: A C++ Variant Integration Kit for next generation sequencing association analysis")


Quadratic Test (SKAT-like)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

\item \textbf{When using expected genotypes:}\\
\begin{itemize}
\item \textbf{vRVS}: ``CAST-like" and ``"SKAT-like" refer to CAST \citep{morgenthaler:2007} and SKAT with weights $w^{1/2}=1/[MAF(1-MAF)]^{1/2}$ \citep{wu:2011}, respectively, where $E(G_{ij} \mid {D_{ij}})$ are used instead of genotype calls (MAF=minor allele frequency). Since the distribution of the expected genotypes given sequence data, $E(G_{ij} \mid {D_{ij}})$, depends on read depth, permutation is not valid. We adopted the bootstrap approach defined in \citet{derkach:2014} for binary trait analysis. We basically use centered expected genotypes for $J$ variants,  $\left[ E(G_{i1} \mid {D_{i1}})- \overline{E(G_{i1} \mid {D_{i1}})} ... E(G_{iJ} \mid {D_{iJ}})- \overline{ E(G_{iJ} \mid {D_{iJ}})} \right.]$ and sample these with replacement, separately for each read depth group. \textit{\textbf{Without covariates}}, we bootstrap on the expected genotypes only. \textit{\textbf{When covariates are added}}, we also bootstrap the added covariates, independently from expected genotypes. For quantitative trait analysis, we implement the permutation methodology defined in \citet{lin:2011} within each read depth group.



Likelihood Method (Coming soon)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This method refers to the test described in *Association testing for next-generation sequencing data using score statistics* _ from Skotte and Albrechtsen (2012 `🔗 <https://www.ncbi.nlm.nih.gov/pubmed/22570057>`) Their method provides a score test where genotype calls are substituted by their expected values, \\(E(G_{ikj}\\mid D_{ikj})\\). The variance of the score test is obtained from the second derivative of the joint likelihood of the observed \\(Y_i\\) and the observed sequencing data, \\(D_{ij}\\) individual \\(i\\) at locus \\(j\\). The p-values are calculated using the asymptotic distribution of the score test. For a joint rare analysis of \\(J\\) variants, the score test is distributed as a chi-square distribution with \\(J\\) degrees of freedom.  This can also be used for common single variant association test which is distributed as chi-squared with one degree of freedom.
