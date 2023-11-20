library(tidyverse)
library(DESeq2)
library(pathfindR)
library(zeallot)
library(ggvenn)

##################################################
#                                                #
#    function for creating DESeq2 object         #
#                                                #
#                                                #
##################################################
##function to produce DESeq2 object from count data and phenotype file, 
##filtering for specific value of specific variable

#countData = normalized gene count matrix
#pheno = phenotype/region of interest type data for each region of interest
#var=variable to filter pheno by
#val=value to filter for

obj<-function(countData,pheno,var='',val=''){
  
#make countData rownames genes
targetNames<-countData$TargetName
countData<-countData[,-1]

#round all count data to be compatible with DESeq
countData<-round(countData)
rownames(countData)<-targetNames

rownames(pheno)<-pheno$sample
#model variables to factors
pheno$psychosis<-factor(pheno$psychosis)
pheno$cortex<-factor(pheno$cortex)
pheno$plaque<-factor(pheno$plaque)


#filter pheno and count
if (var==''){
  DESeqDataSetFromMatrix(countData,pheno , design = ~ age+ PMD+psychosis  )->des
  return(des)
  }else {phenoFilt<-pheno[pheno[[var]]==val,]
countFilt<-countData[,colnames(countData) %in% phenoFilt$sample]
rownames(countFilt)<-targetNames

stopifnot(identical(phenoFilt$sample,colnames(countFilt)))

DESeqDataSetFromMatrix(countFilt,phenoFilt , design = ~ age+PMD+psychosis  )->des

#Filter low gene counts (genes with counts lower than 10 in more than 80% of ROIs)
keep <- rowSums(counts(des) >= 10) >= ncol(counts(des))*0.8
des <- des[keep,]

return(des)
}

}


##################################################
#                                                #
#          dimensionality reduction              #
#                                                #
#                                                #
##################################################
#des=DESeq2 object
#var = variable to assign shapes in UMPA and PCA plots
dim_red<-function(des,var){
  
##UMAP analysis
library(umap)
# Set seed 
set.seed(12345)
#normalisation
vst(des)->de_norm
# retrieve the normalized data and transpose
normalized_counts <- assay(de_norm) %>%
  t() 

# perform UMAP
umap_results <- umap(normalized_counts)

#pull values into dataframe for ploting
umap_plot_df <- data.frame(umap_results$layout) %>%
  tibble::rownames_to_column("sample")

# Add the pheno data into data frame
pheno%>%
  dplyr::select(sample,psychosis,plaque,cortex,BBNID,PMD,sex)->pheno2
pheno2$psychosis<-ifelse(pheno2$psychosis == 1, "AD+P", ifelse(pheno2$psychosis == 0, "AD-P", pheno2$psychosis))
umap_plot_df<-inner_join(umap_plot_df,pheno2,by='sample')

#plot
UMAP<-ggplot(
  umap_plot_df,
  aes(
    x = X1,
    y = X2
  )
) +
  geom_point(aes(color=psychosis,shape=var,group=factor(BBNID)),size=2.5)



#principal component analysis
plotPCA(de_norm,intgroup=c('psychosis',var,'age','PMD'),returnData=T)->PCAdata

PCA<-ggplot(PCAdata,aes(PC1,PC2))+
  geom_point(aes(color=psychosis,shape=var),size=2)

return(list(UMAP,PCA))
}

##################################################
#                                                #
#     Differential gene expression analysis      #
#                                                #
#                                                #
##################################################
##des = DESeq2 object
##filters low counts, performs analysis (AD+P vs AD-P)
##outputs filtered DEGs: abs(LFC)>log2(1.1), adjusted p < 0.05 and volcano plots
DEGA<-function(des){

  #perfoem differential gene expression analysis
  DESeq(des)->des
  
  #extract results and filter for adjusted p-value <0.05 and |LFC| >1, printing after each filter
  results(des,contrast=c('psychosis',1, 0))->res
  as.data.frame(subset(res, padj < 0.05))->DEGs
  print(paste(nrow(DEGs),' DEGs bellow adjusted p-value threshold of 0.05'))
  as.data.frame(subset(DEGs, abs(log2FoldChange) > 1))->DEGs
  print(paste(nrow(DEGs),'DEGs above LFC threshold of 1'))

  
#volcano plot
  
  
 if (nrow(DEGs)>0){
   #label top 50
   data.frame(res) %>% 
    mutate(threshold = padj < 0.05 & abs(log2FoldChange)>1)%>% 
    arrange(padj,desc=F) %>%
    mutate(labthresh = row_number() <= 50 & threshold) ->res_volc
  for (n in 1:nrow(res_volc)){
    if (res_volc$labthresh[n] == T){
      res_volc$labthresh[n]<-rownames(res_volc)[n]
    } else{
      res_volc$labthresh[n]<-''
    }
  }
  
  #plot
  library(ggrepel)
  volcplot <-ggplot(res_volc,aes(x = log2FoldChange, y = -log10(padj), color = threshold,label=labthresh)) + 
    geom_point() + 
    xlab("log2 fold change") + 
    ylab("-log10 adjusted p-value") + 
    theme(legend.position = "none", 
          plot.title = element_text(size = rel(1.5), hjust = 0.5), 
          axis.title = element_text(size = rel(1.25)))+
    geom_vline(xintercept=c(-1, 1), col="red") +
    geom_hline(yintercept=-log10(0.1), col="red")+
    theme_minimal()+
    geom_text_repel(color='black')
  
  return(list(res,DEGs,volcplot))
 }
  #only plots if >1 DEG to avoid error
  else{return(list(res,DEGs))}
}


c(resCortexInner,DEGsCortexInner,VolcCortexInner)%<-%DEGA(obj('cortex','inner'))

c(resCortexOuter,DEGsCortexOuter,VolcCortexOuter)%<-%DEGA(obj('cortex','outer'))

c(resPlaquePresent,DEGsPlaquePresent,VolcPlaquePresent)%<-%DEGA(obj('plaque',T))

c(resPlaqueAbsent,DEGsPlaqueAbsent,VolcPlaqueAbsent)%<-%DEGA(obj('plaque',F))


c(resAll,DEGsAll,VolcAll)%<-%DEGA(obj())

##################################################
#                                                #
#                Venn Diagram                    #
#                                                #
#                                                #
##################################################
#Venn for all groups
D <- list('Inner Cortex'=rownames(DEGsCortexInner),'Outer Cortex'=rownames(DEGsCortexOuter),
          'Plaque Present'=rownames(DEGsPlaquePresent),'Plaque Absent'=rownames(DEGsPlaqueAbsent),
          'All'=rownames(DEGsAll))
ggvenn(D,show_percentage = F,fill_alpha=0.1)+
  ggtitle("AD+P vs AD-P")

#Venn for plaque/plaque free DEGs
D <- list('Plaque Present'=rownames(DEGsPlaquePresent),'Plaque Absent'=rownames(DEGsPlaqueAbsent))
ggvenn(D,show_percentage = F,fill_alpha=0.1)+
  ggtitle("AD+P vs AD-P")

#Venn for inner/outer cortex DEGs
D <- list('Inner Cortex'=rownames(DEGsCortexInner),'Outer Cortex'=rownames(DEGsCortexOuter))
ggvenn(D,show_percentage = F,fill_alpha=0.1)+
  ggtitle("AD+P vs AD-P")

#all ROI, inner cortex, outer cortex, plaque positive venn
D <- list('All'=rownames(DEGsAll),'Inner Cortex'=rownames(DEGsCortexInner),
          'BA+'=rownames(DEGsPlaquePresent),'Outer Cortex'=rownames(DEGsCortexOuter))
ggvenn(D,show_percentage = F,fill_alpha=0.1)+
  ggtitle("AD+P vs AD-P")

##################################################
#                                                #
#               pathway analysis                 #
#                                                #
#                                                #
##################################################

#DEGs = dataframe of DEGs with adjusted p-value and log fold change
#set = gene set list to use
path<-function(DEGs,set='KEGG'){
  #start with ADDH vs C----
  dsDEG<-data.frame(Gene=rownames(DEGs),log2FoldChange=DEGs$log2FoldChange,
                    padj=DEGs$padj)
  
  #pathway analysis----
  run_pathfindR(dsDEG,'KEGG')->pathwayKegg
  
  ggplot(pathwayKegg[1:20,], aes(y = Term_Description, x = Fold_Enrichment,fill=support)) +
    geom_bar(stat = "identity")->bar
  
  enrichment_chart(
    pathwayKegg,
    top_terms =20,
    num_bubbles = 4
  )->enrichment_chart
  return(list(pathwayKegg,bar,enrichment_chart))
}

path(Plaque_signature)->plaqueSignatureKegg
path(DEGsAll)->AllDEGKegg
path(DEGsPlaquePresent)->PlaquePresentKegg
path(DEGsPlaqueAbsent)->PlaqueAbsentKegg
path(DEGsCortexInner)->CortexInnerKegg
path(DEGsCortexOuter)->CortexOuterKegg

##function to search for like terms from pathway analyses output
like_terms<- function(list1,list2){list1[[1]]$Term_Description[list1[[1]]$Term_Description %in% 
                                                                 list2[[1]]$Term_Description]
  
}
#eg of like_terms
like_terms(CortexInnerKegg,CortexOuterKegg)
