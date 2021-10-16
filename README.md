# LSM3241 - Microarray project to find differentially expressed genes after addition of miR-203 to cells. 
Done in collaboration with CHUA YUEN SIONG

## Introduction
Metastasis is the leading cause of death from cancer. miR-203 is a known micro-RNA that
attenuates the migration capabilities of cells by inhibiting the epithelial-mesenchymal transition
(EMT). Yet, the transcriptomic response to miR-203 relevant to EMT was not well elucidated.
Here, we investigated two Geoquery series (GSE50697 and GSE45121) microarray datasets
consisting of treatment samples overexpressing miR-203 and control samples. The objective of
microarray analysis was to identify the differentially expressed genes due to miR-203
expression. The tools used to manipulate and statistically evaluate the datasets were RStudio
and relevant R packages (GEOquery, affy, limma, hgu133plus2.db and mouse4302.db). Genes
identified in RStudio were further evaluated for biological relevance in cancer, EMT and their
possible use as cancer treatments using DAVID Bioinformatics and literature review.
Subsequently the genes were ranked in accordance to our recommendation for their qPCR
validation. In descending order, we recommend the following eight genes: (1) interleukin 6 (IL-6),
(2) nidogen 1, (3) LAPTM5, (4) SH3 and SYLF domain containing 1 (SH3YL1), (5) SERPINB3,
(6) cyclin D2, (7) chondroitin sulfate N-acetylgalactosaminyltransferase 1 (CSGALNACT1) and
(8) parathyroid hormone like hormone (PTHLH). 

## Contents 
1. R Script (GSE45121, mouse).R
2. R Script (GSE50697, human).R
3. Written report
