.. _tests:

Association Tests
==================================

Common Variant Association Test
----------------------------------

For both quantitative and binary trait analyses, a common variant test refers to a score test which has a chi-square distribution with 1 degree of freedom under no association hypothesis. For testing the variant     
:math: $j$, let  $T_j=S^2_j/var(S_j)$ where
    :math: $S_j=\sum_{i=1}^n(Y_i-E(Y_i))G_{ij}$. $E(Y_i)$ is estimated from the fitted values $\hat{Y}$ of $Y$ based on the model assumed. e.g. $\hat{Y}=\bar{Y}$ when there are no covariates. This test is used to perform a genetic association analysis between the phenotype $Y$ and a single variant $G$. For ``True common" and ``Call common", the conventional variance formula of the score statistics is used.  

Note that the score test for vRVS is  $T_j=S^2_j/var(S_j)$ where $S_j=\sum_{i=1}^n(Y_i-E(Y_i))E(G_{ij}\mid D_{ij})$. When we consider integrating sequence data from a arbitrary number of cohorts$(k)$,   $S_j=\sum_{i=1}^{n_k} (Y_{ik}-E(Y_{ik}))E(G_{ikj}\mid D_{ikj})$. The details of the derivation of $Var(S_j)$ are given in the Supplementary document of the paper ``VIKNGS: A C++ Variant Integration Kit for next generation sequencing association analysis".



Rare Variant Association Test
----------------------------------

For joint variant analysis, the score statistics for $J$ joint variants, $\boldsymbol{S}=(S_1,...,S_J)$ is considered and  one can choose a linear (CAST-like) test or a quadratic test (SKAT-like) or Likelihood method (in progress) \citep{skotte:2012} to perform an genetic association test. 

\item \textbf{When using true genotypes or genotype calls in conventional score test:}\\
``CAST-like" and ``"SKAT-like" refer to the CAST \citep{morgenthaler:2007} and SKAT with weights $w^{1/2}=1/[MAF(1-MAF)]^{1/2}$ \citep{wu:2011} respectively, when the true genotypes (in simulation part) and genotype calls are used.  The P-values can be obtain through permutation. 
\textit{\textbf{Without covariates}}, the $Y_i$s are permuted and a test statistics is calculated for each permuted data set where the user defines the number of iterations, e.g. 1,000. At the end of the iterations, P-value is calculated based on the number of values that are less than or equal to the value that was calculated for the initial data set divided by the number of iterations plus 1. \textit{\textbf{When covariates are added}}, we implement the bootstrap methodology defined in \citep{lin:2011} for for binary and quantitative trait analyses.
The user can also use asymptotic distribution to compute the P-values. In that case, the user should define the number of iterations as 1 in VikNGS (e.g. see Figure 2 in the paper ``VIKNGS: A C++ Variant Integration Kit for next generation sequencing association analysis")

\item \textbf{When using expected genotypes:}\\
\begin{itemize}
\item \textbf{vRVS}: ``CAST-like" and ``"SKAT-like" refer to CAST \citep{morgenthaler:2007} and SKAT with weights $w^{1/2}=1/[MAF(1-MAF)]^{1/2}$ \citep{wu:2011}, respectively, where $E(G_{ij} \mid {D_{ij}})$ are used instead of genotype calls (MAF=minor allele frequency). Since the distribution of the expected genotypes given sequence data, $E(G_{ij} \mid {D_{ij}})$, depends on read depth, permutation is not valid. We adopted the bootstrap approach defined in \citet{derkach:2014} for binary trait analysis. We basically use centered expected genotypes for $J$ variants,  $\left[ E(G_{i1} \mid {D_{i1}})- \overline{E(G_{i1} \mid {D_{i1}})} ... E(G_{iJ} \mid {D_{iJ}})- \overline{ E(G_{iJ} \mid {D_{iJ}})} \right.]$ and sample these with replacement, separately for each read depth group. \textit{\textbf{Without covariates}}, we bootstrap on the expected genotypes only. \textit{\textbf{When covariates are added}}, we also bootstrap the added covariates, independently from expected genotypes. For quantitative trait analysis, we implement the permutation methodology defined in \citet{lin:2011} within each read depth group.

To reduce running time, an early stopping procedure can be chosen. It terminates the iterations early if calculations suggest a big P-value \citep{jiang:2012}.

\item \textbf{Likelihood}: The ``Likelihood" method refers to the method described in  \citep{skotte:2012}. Their method provides a score test where genotype calls are substituted by their expected values, $E(G_{ij} \mid {D_{ij}})$. The variance of the score test is obtained from the second derivative of the joint likelihood of the observed $Y_i$ and the observed sequencing data, $D_{ij}$ for individual $i$ at locus $j$. The p-values are calculated using the asymptotic distribution of the score test. For a joint rare aanalysis of $J$ variants, the score test is distributed as a chi-square distribution with $J$ degrees of freedom.  This can also be used for common variant association test which is distributed as chi-squared with one degree of freedom. 


















