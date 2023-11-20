Digital Spatial Profiling of Psychosis in Alzheimers Disease : A Pilot
Study
================
Jacob Thomas-Hegarty

Code used in differential gene expression and pathway analyses of
digital spatial profiling pilot data.

## Study

This study compared transcription in 6 post-mortem brain slices from
participants who had Alzheimer’s disease with psychosis (AD+P) and
Alzheimer’s disease without psychosis (AD-P). Digital spatial profiling
allows manual selection of regions of interest (ROIs) based on cell
morphology and areas of pathology. Our pilot study selected regions
based off the presence/absence of amyloid plaque pathology and
inner/outer cortex regions.

## Analysis

Dimensional reduction techniques were used to visualise gene count data.
Differential gene expression analysis was performed using DESeq2 package
to find differentially expressed genes (DEGs). All analyses compared
transcription in AD+P to transcription in AD-P. Analyses were performed
using all ROIs, only β amyloid plaque positive ROIs, only β amyloid
plaque negative ROIs, only inner cortex ROIs and only outer cortex ROIs.
Volcano plots and DEG lists (padj\<0.05, \|LFC\|\>1) are created for
each analysis. Venn diagrams were created to look at similarities in
DEGs found betrween analysis Pathway enrichment analysis was performed
on DEG lists using the pathfindR package. This script has been altered
to ensure no phenotype or transcriptional information is included to
ensure data protection is adhered to.

This is not a perfect analysis, and the field of spatial transcriptions
is rapidly evolving. This code is published as an aid to individuals
learning R for transcriptional analysis in the hope that it can be of
some use. A report on the findings of this pilot study is available on
request.
