dir="/Users/rudneva/Documents/"
files=list.files(dir, pattern = "ROC_input.*.txt")

i=1
for (f in files){
  if (i!=2){
    assign(paste("a",i,sep=""),read.table(paste(dir, f, sep=""),  header=T))
  }
  i=i+1
}

pdf(paste(dir, "ROC_Curves.pdf", sep=""), height=20, width=20, onefile=TRUE)
par(mfrow=c(4,4))
for (j in 1:(i-1)){
  if (j!=2){
    tmp=eval(parse(text=paste("a",j,sep="")))
    plot(tmp$FP/(tmp$FP+tmp$TN), tmp$TP/(tmp$FN+tmp$TP), col="blue", 
         ylab = "Sensiitivity", 
         xlab="1 - Specififcity", 
         pch=16, 
         main=paste(paste("AMP: ",strsplit(gsub("ROC_input.","", files[j]), "_")[[1]][1], sep=""), "\n", paste("DEL: ", gsub(".txt", "", strsplit(gsub("ROC_input.","", files[j]), "_")[[1]][2]), sep=""), sep=""))
    #legend("topright", 
    #       legend = c(paste("AMP: ",strsplit(gsub("ROC_input.","", files[j]), "_")[[1]][1], sep=""), paste("DEL: ", gsub(".txt", "", strsplit(gsub("ROC_input.","", files[j]), "_")[[1]][2]), sep="")),
    #       col="red", border = F )
  }
  abline(h = 0.6, col="red")
  abline(v = 0.2, col="red")
  j=j+1
}
dev.off()
