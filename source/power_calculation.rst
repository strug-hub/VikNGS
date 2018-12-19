Power Calculator
==================================


\section{Power simulation}
Power simulation package can be used for two purposes:
\begin{enumerate}
\item \textbf{Type I error rare:} The user can test the performance of the association tests with respect to the control of type I error under different sequencing settings, e.g. different combinations of read depths with varying sample sizes, different base calling errors. In this setting odds ratio ($OR$) for a binary trait analysis is set to 1, and the proportion of variation explained by the genetic effect ($R^2$-coefficient of determination) for quantitative trait analysis is set to 0.
\item \textbf{Power analysis:} The user can calculate the minimum sample size required to detect a prespecified effect size, e.g.$OR=1.2$ for a binary trait analysis and $R^2=0.1$  for quantitative trait analysis. 
\end{enumerate}
\textbf{Other Input parameters:}
\begin{itemize}
\item On the left upper corner of the simulation window, the user chooses the type of the analysis, Case-control or Quantitative.
\item \textit{Variant Parameters}: This part includes information on the number of variants, effect size ($OR$ or $R^2$), the minimum and the maximum MAF for the variants.
\item \textit{Run Parameters}:  \textit{High/low cutoff} button defines the cut off value which defines high or low read depth. e.g. high/low cutoff=30 indicates the the cohorts with mean read depth smaller than $30x$ will be considered low read depth cohort. This part is important for the application of the vRVS test as the robust variance estimate in vRVS depends on whether the cohort is high read depth or low read depth. \textit{Steps} button is related to the power analysis. The user can define increment sample sizes on the cohort table, e.g. 500:1500 in controls and Steps=3, this indicates the simulations will be run for number of controls at three different sample size, 500, 1000 and 1500. The results will be saved and plotted on the sample size versus power graph. \textit{\# of threads} defines the ...
\item \textit{Test Parameters}: The tests here are describe in the \textbf{Association Test} part. \textit{\# iteration} defines the number of resampling for permutation/bootstrap step. \textit{Collapse} defines the number of rare variant collapsed. \textit{Stop Early} is the option for stopping the iterations early if calculations suggest a big P-value.
\item \textit{\# Individuals}: This is the sample size corresponding to each group, e.g. case-control, or normal (quantitative) with a particular read depth.
\item \textit{Cohort}: This part indicates cohorts that are combined for the study. The user can add/delete cohorts by clicking on \textit{Add Group} and \textit{Remove Selected} on the menu above the cohort table.
\item \textit{Mean Depth} and  \textit{Depth SD}: The user can define the average read depth and its standard deviation that belong to the corresponding cohort. e.g. mean depth=100 and Depth SD=10 mean the read depth is normally distributed with a mean of 100 and a standard deviation of 10 for that cohort.
\item \textit{Error Rate}: This part defines the base calling error rate. e.g. Error rate=0.01 mean one percent of bases are incorrectly called as one of the three other possible bases.
\item \textit{Quantitative Phenotype}: At the bottom of the cohort table, the user can define the mean and the standard deviation of the quantitative phenotype.



\end{itemize}
\subsection{Output}
The program generates a Q-Q plot and histograms of the p-values when Type I error rate is of interest.  For power analysis, the user can define ... 



Parameters
------------------------------

Please check the :ref:`relevant section <explain_param>` for information regarding parameters not specific to the power calculator. 

Output
------------------------------


