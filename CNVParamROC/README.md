##### Making ROC Curves

module load R/3.3.1
cd /home/vrudneva/data

a_f="quantile(res[which(res[,5]>0),][,5])[4] median(res[,5])+sd(res[,5]) mean(res[,5])+sd(res[,5]) median(res[which(res[,5]>0),][,5])+sd(res[which(res[,5]>0),][,5])"
d_f="quantile(res[which(res[,5]<0),][,5])[4] median(res[,5])-sd(res[,5]) mean(res[,5])-sd(res[,5]) median(res[which(res[,5]<0),][,5])+sd(res[which(res[,5]<0),][,5])"

for a in $a_f
do 
    for d in $d_f
    do
        bsub -M 46000 -R "rusage[mem=46000]" -P methCNV -e AllCNVs.lsferr -o AllCNVs.lsfout "bash ../InfantMB_MethylCNV_PrepareROCInput_ClusterWrapper.sh 'mean(res[,5])+sd(res[,5])' 'mean(res[,5])-sd(res[,5])'"
    done
done


# Moving files
for a in $a_f; do      for d in $d_f;     do mkdir "${a}"_"${d}"; done; done
for a in $a_f; do      for d in $d_f;     do mv *_"${a}"*"${d}"*.pdf "${a}"_"${d}"; done; done





cat CNV_Call_450K.All.mean\(res\[\,5\]\)+sd\(res\[\,5\]\)_mean\(res\[\,5\]\)-sd\(res\[\,5\]\).txt | cut -f100 | perl -ne '
chomp($_); 
my $out=""; 
if ($_ eq "NA") {$out="NA"} 
else 
{
   my @all=split(" ", $_); 
   my %hash=(); 
   foreach my $e (@all) 
   {
       my @info=split("=", $e); 
       my $g=$info[0]; 
       my $r=$info[1]; 
       $hash{$g}=$r;
   } 
   foreach my $k (keys %hash)
   {
       $out=$out.$k."=".$hash{$k};
   } 
} 
print "$out\n";'
