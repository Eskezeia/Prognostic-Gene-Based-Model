Paper "Construction and Validation of a Prognostic Gene-Based Model for Overall Survival Prediction in Hepatocellular Carcinoma Using an Integrated Statistical and Bioinformatic Approach"

##**Introdduction** 


##**Data preparation**

Identifying differentially expressed genes 

Input files
Frist we need to organize the input file for normalized gene expression and overall survival time and status usingTCGA dataset. For example, the input file of survival related gene selection have this kind matrix format.

ID	os	event	GHR	ADH4	LCAT	…	FAM83D
TCGA.WX.AA46.01A.11R.A39D.07	24.84	0	11.68619	16.60054	12.11408	…	5.088862
TCGA.BC.A10X.01A.11R.A131.07	25.3	1	13.36989	16.53799	11.884	…	6.215566
TCGA.G3.AAV0.01A.11R.A37K.07	15.64	0	12.4672	15.21704	12.72918	…	6.026993
TCGA.G3.A5SK.01A.11R.A27V.07	24.44	0	13.11411	15.94487	11.9492	…	6.160402
TCGA.NI.A4U2.01A.11R.A28V.07	58.84	1	12.90194	17.29704	12.0351	…	6.523212
TCGA.MR.A520.01A.11R.A266.07	7.52	0	12.04749	15.69578	11.81886	…	5.510978
TCGA.K7.A5RF.01A.11R.A28V.07	20.73	0	12.50312	16.17981	12.63151	…	6.550564
TCGA.DD.A4NL.01A.11R.A28V.07	56.21	0	12.31974	16.77668	12.15429	…	6.276694
TCGA.BC.A110.01A.11R.A131.07	69.51	1	12.05064	15.32384	12.51753	…	6.291554
TCGA.UB.AA0V.01A.11R.A38B.07	10.32	0	12.29882	16.16648	11.78306	…	6.36227
TCGA.2Y.A9GV.01A.11R.A38B.07	83.18	1	13.06243	14.13687	11.85227	…	6.529362
TCGA.DD.A3A2.01A.11R.A213.07	70.01	1	12.28312	16.82152	10.98399	…	6.098186
TCGA.G3.A3CI.01A.11R.A213.07	5.91	0	11.59614	15.01456	11.79008	…	6.117502
TCGA.G3.A7M8.01A.11R.A33R.07	14.13	0	10.21954	15.01972	12.28718	…	5.522605
TCGA.LG.A9QD.01A.11R.A38B.07	12.02	0	10.87814	15.32892	11.07987	…	5.302374
TCGA.DD.A4NK.01A.11R.A28V.07	39.75	1	11.32406	14.00305	12.48706	…	6.3055



Datasets used in this study
1. TCGA-LIHC dataset(http://gdac.broadinstitute.org/)
2. GEne expression profiles from the gene expression omnibus (GEO)(https://www.ncbi.nlm.nih.gov/geo/) inccluding GSE112790, GSE84402, and GSE45267.
3. International Cancer Genome Consortium (ICGC LIRI-JP dataset) dataset(https://icgc.org/)
4. China Medical University Hospital HCC dataset




