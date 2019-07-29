setwd("/Users/rudnevav/Documents/Projects/Uveal_melanoma/Results")
library("dplyr")
library("circlize")
library("RColorBrewer")
library("wesanderson")
library("ComplexHeatmap")
library("CNTools")
library("ComplexHeatmap")
library("dendextend")
library("dendsort")
library("gplots")
library("ggplot2")
options(stringsAsFactors = FALSE)
source("/Users/rudnevav/scripts/Funtions/ProcessChromData.R")

#load("Oncoprint_data.RData")
load("tmp.RData")
load("IMPACT-Annots-info.md2.RData")

## MSKCC-IMPACT
impact.annots.sel=md2; impact.annots.sel$Merged_Sample_ID=unlist(lapply(impact.annots.sel$Merged_Sample_ID, function(x) paste0(unlist(strsplit(x, split = "-"))[1:3], collapse = "-")[[1]]))
impact.annots.sel = t(impact.annots.sel); colnames(impact.annots.sel)=impact.annots.sel[1,]; impact.annots.sel=impact.annots.sel[-1,]; rownames(impact.annots.sel)[dim(impact.annots.sel)[1]]="Ploidy"
TCGA_Cluster=rep("MSKCC", dim(impact.annots.sel)[2]); impact.annots.sel=rbind(impact.annots.sel, TCGA_Cluster)
## TCGA
tcga.md2=read.table("../TCGA-UVM_MolClinData_processed.txt", header = T, sep = "\t")
ReplaceVariantsWithMUTs3 = function(tmp){
  tmp=as.matrix(tmp)
  tmp[is.na(tmp)]="na"
  tmp[!tmp %in% c("0", "B;", "gain", "loss", "Q209P", "Q209L", "R183H", "R183Q", "R183C", "na", "AMP;", "HOMDEL;", "LOH;", "LOH")]="MUT;"
  tmp[tmp %in% c("0")]=" "
  tmp[tmp %in% c("LOH")]="LOH;"
  tmp[tmp %in% c("B;")]=" "
  return(tmp)
}
tmp=ReplaceVariantsWithMUTs3(tcga.md2[,4:12])
tcga.md2=cbind(tcga.md2[,1:3], tmp, tcga.md2[,13:15])
tcga.annots.sel=t(tcga.md2); colnames(tcga.annots.sel)=tcga.annots.sel[1,]; tcga.annots.sel=tcga.annots.sel[-1,]; rownames(tcga.annots.sel)[1]="TCGA_Cluster"
Recurrent=rep(0, dim(tcga.annots.sel)[2]); tcga.annots.sel=rbind(tcga.annots.sel, Recurrent)
## Combine
annots.sel=cbind(impact.annots.sel[rownames(tcga.annots.sel),], tcga.annots.sel)

mydata.nr2=t(mydata.nr); colnames(mydata.nr2)=unlist(lapply(colnames(mydata.nr2), function(x) paste(unlist(strsplit(x, split = "-"))[1:3], collapse = "-")))
chr.label=MakeSeparatorForChromHeatmap2(t(mydata.nr2))

## Define colors for annotations
col = c("LOH;"="lightblue", "HOMDEL;" = "blue", "AMP;" = "red", "MUT;" = "#008000", "Q209L"="darkolivegreen1", "Q209P"="green", "R183Q"="coral", "R183H"="coral", "R183C"="coral", "na"="gray90")
col_fun_ploidy = colorRamp2(breaks = c(0, 4), colors = c("white", "deeppink3"))
col_yes_no = c("0"="darkgreen", "1"="coral")
col.tmp=c("#F1BB7B", "red", "lightblue", "blue", "white"); names(col.tmp)=c("1","2", "3", "4", "MSKCC")

### Clustering
col_dend1 = hclust(dist(t(mydata.nr2), method = "manhattan"), method = "ward.D"); pt.order1=colnames(mydata.nr2[,col_dend1$order])
#annots.sel1=PrepareAnnots.Sel(annots.sel, pt.order1, mydata.nr2); clust1=PrepareClust(mydata.nr2, tcga.molclin, pt.order1)

#ha = HeatmapAnnotation("TCGA Clusters"=clust1, col=list("TCGA Clusters"=col.tmp), gp = gpar(fontsize = 18, col="black"))

annots.sel=t(annots.sel)
annots.sel=annots.sel[rownames(annots.sel)[rownames(annots.sel) %in% pt.order1],]
annots.sel=annots.sel[pt.order1,]

ha = HeatmapAnnotation("TCGA Clusters"=annots.sel[pt.order1,"TCGA_Cluster"], col=list("TCGA Clusters"=col.tmp), gp = gpar(fontsize = 18, col="black"))
hb = HeatmapAnnotation(
  #gp = gpar(fontsize = 12, col="black"),
  show_annotation_name = T, 
  show_legend=T,
  #column_order=pt.order,
  annotation_name_gp = gpar(fontsize = 12),
  Recurrence = annots.sel[,"Recurrent"],  
  GNAQ = annots.sel[,"GNAQ"],
  GNA11 = annots.sel[,"GNA11"],
  EIF1AX = annots.sel[,"EIF1AX"],
  SF3B1 = annots.sel[,"SF3B1"],
  BAP1 = annots.sel[,"BAP1"],
  chr3p = annots.sel[,"chr3p"],
  chr6p = annots.sel[,"chr6p"],
  chr6q = annots.sel[,"chr6q"],
  chr8q = annots.sel[,"chr8q"],
  Ploidy = as.numeric(annots.sel[,"Ploidy"]),
  col = list(
    Recurrence = col_yes_no,
    GNAQ = col,
    GNA11 = col,
    EIF1AX = col,
    SF3B1 = col,
    BAP1 = col,
    chr3p = col,
    chr6p = col,
    chr6q = col,
    chr8q = col,
    Ploidy = col_fun_ploidy),
  gp = gpar(fontsize = 12, col="black"),
  na_col = "gray90")

h3=Heatmap(mydata.nr2[,pt.order1], 
           width = unit(100, "cm"), height = unit(40, "cm"),   
           name = "Copy Number", #title of legend
           row_split = factor(chr.label[-length(chr.label)], levels = seq(1:22)),
           cluster_rows = FALSE,
           cluster_row_slices = FALSE,
           show_row_names = FALSE,
           
           column_dend_height = unit(3, "cm"),
           column_dend_reorder =T,
           clustering_distance_columns = "manhattan",
           clustering_method_columns = "ward.D",
        
           #column_split = 4,
           column_km =4,
           column_km_repeats = 100,
           #cluster_column_slices = F,
           column_title_gp = gpar(fill = c("lightblue", "red", "blue", "#F1BB7B"), font = 2, fontsize = 18),
           column_names_gp = gpar(col = c("lightblue", "red","blue","#F1BB7B"), fontsize = 18),

           gap = unit(2, "mm"),
           bottom_annotation = hb,
           top_annotation = ha
)
pdf("19.TCGA_IMPACT_Clustering.pdf",useDingbats=F, height=30, width = 45)
draw(h3,
     ht_gap = unit(3, "mm"),
     #row_title = "Alterations",
     #column_title = column_title.sel,
     #column_title_gp = gpar(fontsize = 12, fontface = "bold"),
     merge_legend = TRUE,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     adjust_annotation_extension = FALSE,
     row_title = "Chromosome",
     column_title = "Manhattan/Ward.D Clustering"
)
dev.off()
