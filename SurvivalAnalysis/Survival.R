library(survival)
library(survminer)

dataset=read.table("/Users/vrudneva/Dropbox/Vasilisa-shared/GeneExpression/U133p2-R2-MB-series-April13-17.txt", header=T, stringsAsFactors = F)
annot=t(dataset[1:8,]); colnames(annot)=annot[2,]; annot=annot[-1,]; annot=annot[-1,]; annot=as.data.frame(annot, stringsAsFactors=F)
dataset=dataset[-1,]; dataset=dataset[-1,]; dataset=dataset[-1,]; dataset=dataset[-1,]; dataset=dataset[-1,]; dataset=dataset[-1,]; dataset=dataset[-1,]; dataset=dataset[-1,]

gene="CCR2"

tmp=t(dataset[dataset$H.hugo==gene,]); md=as.data.frame(cbind(rownames(tmp), tmp[,1]), stringsAsFactors=F); colnames(md)=md[1,];md=md[-1,]; md=md[-1,]
md=merge(md, annot, by=0, all=T)[,-1]; md$death=as.numeric(md$death); md$CCR2=as.numeric(md$CCR2); md$follow.up_.months.=as.numeric(md$overall_survival)
md=md[md$follow.up_.months.<=120,]; md=na.omit(md)


MakeKaplanMeier=function(data, DataSet, LowHighMetod){
  fit <- survfit(Surv(follow.up_.months., death) ~ Expr, data = data, conf.type = "log-log")
  res=survdiff(Surv(follow.up_.months., death) ~ Expr, data = data)
  
  plot(fit, lty = 1, col = c("darkred", "darkblue"), xlab = "Follow-up time in month", ylab="Overall survival", mark.time=T)
  legend("topright", legend = c(paste("High (n=", dim(data[data$Expr=="High",])[1], ")", sep=""), paste("Low (n=", dim(data[data$Expr=="Low",])[1], ")", sep="")), 
         lty = 1, col = c("darkred", "darkblue"))
  title(paste("Kaplan-Meier Curves for SHH-MB, ", DataSet, " data\n(", LowHighMetod, ")", sep="")) 
  txt=1 - pchisq(res$chisq, length(res$n) - 1)
  text(round(range(data$follow.up_.months.)[2]/7),0, paste("P-value = ",format(txt, digits = 2), sep=""))
}

MakeKaplanMeier2=function(data, DataSet, LowHighMetod){
  fit <- survfit(Surv(follow.up_.months., death) ~ Expr, data = data)

  g=ggsurvplot(fit, data=data, risk.table = T, pval=TRUE, 
             main=paste("Kaplan-Meier Curves for SHH-MB, ", DataSet, " data\n(", LowHighMetod, ")", sep=""), 
             legend.labs = c(paste("High (n=", dim(data[data$Expr=="High",])[1], ")", sep=""), paste("Low (n=", dim(data[data$Expr=="Low",])[1], ")", sep="")),
             legend.title=NULL)
  return(g)
  
}

md1=read.table("/Users/vrudneva/Dropbox/Vasilisa-shared/AIF1_Survival/GeneArray.txt", header=T, stringsAsFactors = F)
md1$death=as.numeric(md1$death)
md1$Expr <- ifelse(md1$AIF1 >=mean(md1$AIF1), "High", "Low")

md2=read.table("/Users/vrudneva/Dropbox/Vasilisa-shared/AIF1_Survival/Plus2.txt", header=T, stringsAsFactors = F)
md2$Expr <- ifelse(md2$AIF1 >=mean(md2$AIF1), "High", "Low")


md1=md1[md1$follow.up_.months.<=120,]; md1=na.omit(md1)
md2=md2[md2$follow.up_.months.<=120,]; md2=na.omit(md2)

dim(md1)[1]; dim(md2)[1]


pdf("/Users/vrudneva/Dropbox/Vasilisa-shared/AIF1_Survival/K-M_Curves.UpTo10years.pdf")  

MakeKaplanMeier(md1, "GeneArray", "compared to mean")
summary(survfit(Surv(follow.up_.months., death) ~ Expr, data = md1), times=c(60,120))
MakeKaplanMeier(md2, "Plus2", "compared to mean")
summary(survfit(Surv(follow.up_.months., death) ~ Expr, data = md2), times=c(60,120))


md1$Expr=rep(NA, dim(md1)[1])
md1[md1$AIF1>(mean(md1$AIF1)+sd(md1$AIF1)),]$Expr="High"
md1[md1$AIF1<(mean(md1$AIF1)-sd(md1$AIF1)),]$Expr="Low"
md1=na.omit(md1)

MakeKaplanMeier(md1, "GeneArray", "compared to 1st/3rd quantiles")
summary(survfit(Surv(follow.up_.months., death) ~ Expr, data = md1), times=c(60,120))


md2$Expr=rep(NA, dim(md2)[1])
md2[md2$AIF1>(mean(md2$AIF1)+sd(md2$AIF1)),]$Expr="High"
md2[md2$AIF1<(mean(md2$AIF1)-sd(md2$AIF1)),]$Expr="Low"
md2=na.omit(md2)

MakeKaplanMeier(md2, "Plus2", "compared to 1st/3rd quantiles")
summary(survfit(Surv(follow.up_.months., death) ~ Expr, data = md2), times=c(60,120))


md=rbind(md1,md2)

MakeKaplanMeier(md, "Combined GeneArray/Plus2", "compared to 1st/3rd quantiles")
summary(survfit(Surv(follow.up_.months., death) ~ Expr, data = md), times=c(60,120))


dev.off()


# Using the new gene
pdf(paste("/Users/vrudneva/Dropbox/Vasilisa-shared/AIF1_Survival/", gene,".K-M_Curves.UpTo10years.pdf", sep=""))  

md$Expr <- ifelse(md$CCR2 >=mean(md$CCR2), "High", "Low")
MakeKaplanMeier(md, "Plus2", "compared to mean")
summary(survfit(Surv(follow.up_.months., death) ~ Expr, data = md), times=c(60,120))


md$Expr <- ifelse(md$CCR2 >=median(md$CCR2), "High", "Low")
MakeKaplanMeier(md, "Plus2", "compared to median")
summary(survfit(Surv(follow.up_.months., death) ~ Expr, data = md), times=c(60,120))


md$Expr=rep(NA, dim(md)[1])
md[md$CCR2>quantile(md$CCR2)[4],]$Expr="High"
md[md$CCR2<quantile(md$CCR2)[2],]$Expr="Low"
md=na.omit(md)



MakeKaplanMeier(md, "Plus2", "compared to 1st/3rd quantiles")
summary(survfit(Surv(follow.up_.months., death) ~ Expr, data = md), times=c(60,120))


md$Expr=rep(NA, dim(md)[1])
md[md$CCR2>(mean(md$CCR2)+sd(md$CCR2)),]$Expr="High"
md[md$CCR2<(mean(md$CCR2)-sd(md$CCR2)),]$Expr="Low"
md=na.omit(md)

MakeKaplanMeier(md, "Plus2", "mean +- SD")
summary(survfit(Surv(follow.up_.months., death) ~ Expr, data = md), times=c(60,120))

dev.off()



