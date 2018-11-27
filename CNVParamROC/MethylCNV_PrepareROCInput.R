setwd(dir="/Users/rudneva/Documents")

print(paste("Amp: ", amp_formula, " Del: ", del_formula, sep=""))
options(warn=-1)

suppressMessages(library("minfi")) 
suppressMessages(library("IlluminaHumanMethylation450kanno.ilmn12.hg19")) 
suppressMessages(library("IlluminaHumanMethylation450kmanifest")) 
suppressMessages(library("minfiData"))
suppressMessages(library("stringr"))
suppressMessages(library("CopyNumber450kData"))
suppressMessages(library("conumee"))

data(RGcontrolSetEx)
data(exclude_regions)
data(detail_regions)

require(data.table)
masterTable<-t(as.data.frame(fread("heatmap.tsv")))
masterTable=as.data.frame(masterTable, stringsAsFactors = F)
rownames(masterTable)=gsub("0;;", "", rownames(masterTable))
colnames(masterTable)=masterTable[1,]
masterTable=masterTable[-1,]
masterTable$PID=rownames(masterTable)
masterTable=masterTable[,-1]

load("query.data.RData")
load("cytoBand.RData")

samples=names(query.data)

ROC_input=as.data.frame(cbind(rep(NA, length(samples)), rep(NA, length(samples)), rep(NA, length(samples)), rep(NA, length(samples)), rep(NA, length(samples)), rep(NA, length(samples))), ncol=6)
colnames(ROC_input)=c("TP", "TN", "FP", "FN", "P", "N"); rownames(ROC_input)=masterTable[masterTable$IDAT %in% samples,]$PID

logFile = paste(amp_formula,"_", del_formula,".log_file.txt", sep="")
cat("Excluded samples:", file=logFile, append=FALSE, sep = "\n")


for (a2.filename in samples)
{
  # The chromosomes that were studied
  chrs=eval(parse(text=paste("chr_",masterTable[masterTable$IDAT==a2.filename,]$Subgroup, sep="")))
  x <- CNV.fit(query.data[a2.filename], controls.data, anno)
  
  x <- CNV.bin(x)
  x <- CNV.detail(x)
  x <- CNV.segment(x)
  
  #CNV.detailplot_wrap(x)
  out=""
  res=CNV.write(x, what="bins")
  
    ext_low=res[which(res[,5]< -0.4),]
    ext_high=res[which(res[,5]> 0.4),]
    
    qc_out=vector()
    
    e=ext_low
    print("Extreme low...")
    for (c in paste("chr", seq(1:22), sep="")){
      qc_out=c(qc_out, mean(c(sort(e[e$Chromosome==c,]$Start),0)-c(0,sort(e[e$Chromosome==c,]$Start))[1:length(e[e$Chromosome==c,]$Start)]))
      
    }
    e=ext_high
    print("Extreme high...")
    for (c in paste("chr", seq(1:22), sep="")){
      qc_out=c(qc_out, mean(c(sort(e[e$Chromosome==c,]$Start),0)-c(0,sort(e[e$Chromosome==c,]$Start))[1:length(e[e$Chromosome==c,]$Start)]))
    }
    
    if (length(qc_out[qc_out> 0]) <= 15){
      
      amp_cutoff=eval(parse(text=amp_formula))
      del_cutoff=eval(parse(text=del_formula))
      
      for (i in 1:22)
      {
        chr=paste("chr", i, sep="")
        p.start=cyto.mappings[i,1];p.end=cyto.mappings[i,2]
        q.start=cyto.mappings[i,3];q.end=cyto.mappings[i,4]
        
        this_res=res[res$Chromosome==chr,]
        med.p=median(as.numeric(this_res[(this_res$Start>=p.start)&(this_res$End<=p.end),][,5]))
        med.q=median(as.numeric(this_res[(this_res$Start>=q.start)&(this_res$End<q.end),][,5]))
        
        if(is.na(med.p)){out=paste(out, chr, "p=NA; ", sep="")} else {
          if (med.p > amp_cutoff) {out=paste(out, chr, "p=AMP; ", sep="")} else {
            if (med.p < del_cutoff) {out=paste(out, chr, "p=DEL; ", sep="")} else {out=paste(out, chr, "p=0; ", sep="")}
          }
        }
        if(is.na(med.q)){out=paste(out, chr, "q=NA; ", sep="")} else {
          if (med.q > amp_cutoff) {out=paste(out, chr, "q=AMP; ", sep="")} else {
            if (med.q < del_cutoff) {out=paste(out, chr, "q=DEL; ", sep="")} else {out=paste(out, chr, "q=0; ", sep="")}
          }
        }
      }
            
      WGS=masterTable[masterTable$IDAT==a2.filename,51:94]
      WGS=WGS[which(!is.na(WGS))]; 
      if (dim(WGS)[2]!=0){WGS=WGS[which(WGS!=0)]} 
      if (dim(WGS)[2]==0){out1_a="No"; out1_d="No"} else {
        WGS_A=WGS[which(WGS=="AMP")]; WGS_D=WGS[which(WGS=="DEL")]
        out1_a=""; for (c in sort(names(WGS_A))){out1_a=paste(out1_a, c, ";", sep="")}
        out1_d=""; for (c in sort(names(WGS_D))){out1_d=paste(out1_d, c, ";", sep="")}
      }
      
      Meth_CNV=unlist(lapply(unlist(strsplit(out, split = "; ")), function(x) strsplit(x,"=")[[1]][2])); Meth_CNV[Meth_CNV=="NA"]=NA
      names(Meth_CNV)=names(masterTable[masterTable$IDAT==a2.filename,51:99])[1:44]; Meth_CNV=Meth_CNV[which(!is.na(Meth_CNV))]; Meth_CNV=Meth_CNV[which(Meth_CNV!=0)]
      if (length(Meth_CNV)==0){out2_a="No"; out2_d="No"} else {
        Meth_CNV_A=Meth_CNV[which(Meth_CNV=="AMP")]; Meth_CNV_D=Meth_CNV[which(Meth_CNV=="DEL")]
        out2_a=""; for (c in sort(names(Meth_CNV_A))){out2_a=paste(out2_a, c, ";", sep="")}
        out2_d=""; for (c in sort(names(Meth_CNV_D))){out2_d=paste(out2_d, c, ";", sep="")}
      }
      
      pdf(paste(masterTable[masterTable$IDAT %in% a2.filename,]$PID, "._", amp_formula,"_", del_formula, "_.CNV.WGS_Meth450K_Comparison.pdf", sep=""), height=20, width=20, onefile=TRUE)
      CNV.genomeplot(x)
      legend("topleft",legend = paste("WGS:\nAMP: ", out1_a, "\nDEL: ",out1_d,"\n\nMeth_450K:\nAMP: ", out2_a, "\nDEL: ", out2_d, "\nChrs studied: ", paste(chrs, collapse = ";"), sep=""), cex=0.75, x.intersp=1, y.intersp=7)
      dev.off()
      
      # Calculate sensitivity and specificity
      Meth_CNV=t(as.data.frame(Meth_CNV))
      P=dim(WGS)[2]; N=length(chrs)-P; TP=0; TN=0; FN=0; FP=0 
      
      for (chr in chrs)
      {
        if (chr %in% names(WGS)){
          #print(paste("true alteration in:", chr, sep=""))
          if (chr %in% colnames(Meth_CNV)){
            if (WGS[,chr]==Meth_CNV[,chr]) {TP=TP+1} else {FN=FN+1}
          } else {FN=FN+1}
        } else {
          if (!chr %in% colnames(Meth_CNV)){TN=TN+1} else {FP=FP+1}
        }
      }
      ROC_input[masterTable[masterTable$IDAT %in% a2.filename,]$PID,]$TP=TP
      ROC_input[masterTable[masterTable$IDAT %in% a2.filename,]$PID,]$TN=TN
      ROC_input[masterTable[masterTable$IDAT %in% a2.filename,]$PID,]$FP=FP
      ROC_input[masterTable[masterTable$IDAT %in% a2.filename,]$PID,]$FN=FN
      ROC_input[masterTable[masterTable$IDAT %in% a2.filename,]$PID,]$P=P
      ROC_input[masterTable[masterTable$IDAT %in% a2.filename,]$PID,]$N=N
    } else {
    cat(masterTable[masterTable$IDAT %in% a2.filename,]$PID, file=logFile, append=TRUE, sep = "\n")
    
    pdf(paste("Excluded.", masterTable[masterTable$IDAT %in% a2.filename,]$PID, "._", amp_formula,"_", del_formula, "_.CNV.WGS_Meth450K_Comparison.pdf", sep=""), height=20, width=20, onefile=TRUE)
    CNV.genomeplot(x)
    dev.off()
  }
}

write.table(x = ROC_input, file = paste("ROC_input.", amp_formula,"_", del_formula, ".txt", sep=""), sep="\t", quote=F)
