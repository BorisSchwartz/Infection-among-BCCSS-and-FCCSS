SAS program created by Boris Schwartz, CESP-INSERM-U1018 Radiation Epidemiology Team (France) - 2024

Late Infection-Related Risk among Childhood Solid Cancer Survivors: A binational study from the French and British Childhood Cancer Survivor Studies	

Free space required on your computer : around 30 Go.
Note: Some warnings may appear in the log, but this has no effect on the program.
Note: results appear as sas tables and/or rtf files (Word format). These are crude internal working documents. Result tables in the article have been formatted then from these documents.

Global program with:

1.Set up (change folders and library according to your computer to execute the following program)

2.Code list for infections (ICD-9-10)

3.frequently used macros

4.Sankey diagram (/*Save sankey scripts (rawtosankey.sas, sankey.sas and saneybarchart.sas in the folder before). SAS Sankey macro created by Shane Rosanbalm of Rho, Inc. 2015 https://github.com/RhoInc/sas-sankeybarchart*/)

5.Calculation of PY : from a base with one line per patient to a base with one line per patient per year

6.SHR/AER Calculation: example for inf_grave endpoint (To adapt and reproduce for each endpoint)

7.Macro for relative risks from Poisson regression and example for inf_grave endpoint

8.Cumulative incidence function (CIF): example for inf_grave endpoint

9.Additionnal code for sensitivity: COX ph, Andersen Gill method for recurrent events and Marginal means and rate model for recurrent events analysis (/*_See tutorial in appendix of:
Amorim, L. D., & Cai, J. (2015). Modelling recurrent events: a tutorial for analysis in epidemiology. 
International journal of epidemiology, 44(1), 324–333. https://doi-org.proxy.insermbiblio.inist.fr/10.1093/ije/dyu222 _*/)

