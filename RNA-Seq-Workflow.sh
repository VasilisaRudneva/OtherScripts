### Quality control

cd /QC/mapped

# To get only the mapped reads:
samtools view -b -F 4 /data/others/vasilisa/runs/150728_SN999_0239_BC778EACXX/3NINACp5WGA/lane215s012249_sequence.bam > lane215s012249_sequence.mapped.bam &
samtools view -b -F 4 /data/others/vasilisa/runs/150728_SN999_0239_BC778EACXX/3NINACs1WGA/lane215s012247_sequence.bam > lane215s012247_sequence.mapped.bam &
samtools view -b -F 4 /data/others/vasilisa/runs/150728_SN999_0239_BC778EACXX/3NINACs2WGA/lane215s012248_sequence.bam > lane215s012248_sequence.mapped.bam &
samtools view -b -F 4 /data/others/vasilisa/runs/150728_SN999_0239_BC778EACXX/3NIp5WGA/lane215s012246_sequence.bam > lane215s012246_sequence.mapped.bam &
samtools view -b -F 4 /data/others/vasilisa/runs/150728_SN999_0239_BC778EACXX/3NIs1WGA/lane215s012244_sequence.bam > 3NIs1WGA.mapped.bam &
samtools view -b -F 4 /data/others/vasilisa/runs/150728_SN999_0239_BC778EACXX/3NIs2WGA/lane215s012245_sequence.bam > 3NIs2WGA.mapped.bam
samtools view -b -F 4 /data/others/vasilisa/runs/150728_SN999_0239_BC778EACXX/3NIs2WGA/lane215s012245_sequence.bam > 3NIs2WGA.mapped.bam &
samtools view -b -F 4 /data/others/vasilisa/runs/150728_SN999_0239_BC778EACXX/4OnOffNACp5WGA/lane215s012255_sequence.bam > 4OnOffNACp5WGA.mapped.bam &
samtools view -b -F 4 /data/others/vasilisa/runs/150728_SN999_0239_BC778EACXX/4OnOffNACs2WGA/lane215s012254_sequence.bam > 4OnOffNACs2WGA.mapped.bam &
samtools view -b -F 4 /data/others/vasilisa/runs/150728_SN999_0239_BC778EACXX/4OnOffp5WGA/lane215s012252_sequence.bam > 4OnOffp5WGA.mapped.bam &
samtools view -b -F 4 /data/others/vasilisa/runs/150728_SN999_0239_BC778EACXX/4OnOffs1WGA/lane215s012250_sequence.bam > 4OnOffs1WGA.mapped.bam &

for file_name in *.mapped.bam
do
    echo "infile='${file_name}'" | cat - /scripts/QC_plot_coverage_per_chr.R | R --vanilla --slave &
done 

### Controls
c1="/data/others/chris/runs/150617_SN999_0232_AC5HD3ACXX/HC1/lane3HC1_sequence.bam"
c2="/data/others/chris/runs/150617_SN999_0232_AC5HD3ACXX/HC3/lane3HC3_sequence.bam"
c3="/data/others/chris/runs/150617_SN999_0232_AC5HD3ACXX/HC4/lane3HC4_sequence.bam"

samtools view -b -F 4 /data/others/chris/runs/150617_SN999_0232_AC5HD3ACXX/HC1/lane3HC1_sequence.bam > control_lane3HC1_sequence.mapped.bam
samtools view -b -F 4 /data/others/chris/runs/150617_SN999_0232_AC5HD3ACXX/HC3/lane3HC3_sequence.bam > control_lane3HC3_sequence.mapped.bam
samtools view -b -F 4 /data/others/chris/runs/150617_SN999_0232_AC5HD3ACXX/HC4/lane3HC4_sequence.bam > control_lane3HC4_sequence.mapped.bam
samtools view -b -F 4 /data/others/chris/runs/150617_SN999_0232_AC5HD3ACXX/HC2/lane3HC2_sequence.bam > control_lane3HC2.mapped.bam



for file_name in control_*.mapped.bam
do
    echo "infile='${file_name}'" | cat - /scripts/QC_plot_coverage_per_chr.R | R --vanilla --slave &
done 



#### Read Depth Plots

REFGENOME="/datasets/refgenomes/mouse/Mus_musculus_mm10.fa"
COV_TOOL="/software/tools/seqtools/cov"
GCNORM="/readDepth/chris_gcNormCov.R"
REFGENOME_GC="/datasets/refgenomes/mouse/Mus_musculus_mm10.10kb.GC"

for file_name in *.mapped.bam
do
    bsub -M 46000 -R "rusage[mem=46000]" -e ReadDepth_Cov.${file_name}.${this_window_size}.lsferr -o ReadDepth_Cov.${file_name}.${this_window_size}.lsfout "/software/tools/seqtools/cov -d -q 5 -w 10000 -o 10000 -g ${REFGENOME} ${file_name} > ${file_name}_10kb.cov; /software/bin/Rscript-2.15.0 $GCNORM ${file_name}_10kb.cov $REFGENOME_GC ${file_name}_10kb.gcnorm.cov"
done


SAMPLES="3NINACp5WGA.mapped.bam_10kb.gcnorm.cov
3NINACs1WGA.mapped.bam_10kb.gcnorm.cov
3NINACs2WGA.mapped.bam_10kb.gcnorm.cov
3NIp5WGA.mapped.bam_10kb.gcnorm.cov
3NIs1WGA.mapped.bam_10kb.gcnorm.cov
3NIs2WGA.mapped.bam_10kb.gcnorm.cov
4OnOffNACp5WGA.mapped.bam_10kb.gcnorm.cov
4OnOffNACs2WGA.mapped.bam_10kb.gcnorm.cov
4OnOffp5WGA.mapped.bam_10kb.gcnorm.cov
4OnOffs1WGA.mapped.bam_10kb.gcnorm.cov"
CONTROLS="control_lane3HC1_sequence.mapped.bam_10kb.gcnorm.cov
control_lane3HC3_sequence.mapped.bam_10kb.gcnorm.cov
control_lane3HC4_sequence.mapped.bam_10kb.gcnorm.cov"


for control_file in $CONTROLS
do
    for sample_file in $SAMPLES
    do
	echo "Sample:" ${sample_file} "; Control:" ${control_file}
        echo "sample_file='${sample_file}'; control_file='${control_file}'" | cat - /scripts/QC_plot_ReadDepth_per_chr.R | R --vanilla --slave &
    done
done

# Seqstat
bsub < /scripts/seqstat_mouse.sh


# Preparing VCF file

cd /mouse_vcf
strains=`zcat SC_MOUSE_GENOMES.genotype.vcf.gz | grep -v "##" | head -1 | cut -f10- | sed 's/\//_/g'`

count=10
for this_strain in $strains
do
    echo $count $this_strain
    bash /scripts/make_vcf_for_each_mouse_strain.sh $this_strain $count &
    let count++
done


#### FVB mouse strain
cat 2012-0612-snps_only_FVBNJ_annotated.vcf | grep "#" > header
cat 2012-0612-snps_only_FVBNJ_annotated.vcf | grep -v "#" | grep -v "" > input
perl /scripts/fix_contigs_in_vcf_file_for_GATK.pl input /datasets/refgenomes/mouse/Mus_musculus_mm10.fa.fai  > output
java -jar /soft/GenomeAnalysisTK.jar -T LiftoverVariants -R /datasets/refgenomes/mouse/Mus_musculus_mm10.fa -V header -chain /soft/mm9ToMm10.over.chain -dict /datasets/refgenomes/mouse/Mus_musculus_mm10.dict -o 2012-0612-snps_only_FVBNJ_annotated.liftover.mm10.vcf

# LiftOver coordinates for VCF file from mm9 to mm10
java -jar /soft/GenomeAnalysisTK.jar -T LiftoverVariants -R /datasets/refgenomes/mouse/Mus_musculus_mm10.fa -V 2012-0612-snps_only_FVBNJ_annotated.vcf -chain /g/korbel/rudneva/soft/mm9ToMm10.over.chain -dict /g/korbel/shared/datasets/refgenomes/mouse/Mus_musculus_mm10.dict -o 2012-0612-snps_only_FVBNJ_annotated.liftover.mm10.vcf

cat 2012-0612-snps_only_FVBNJ_annotated.liftover.mm10.vcf | grep "#" > 2012-0612-snps_only_FVBNJ_annotated.liftover.mm10.sorted.vcf; grep -v "#" 2012-0612-snps_only_FVBNJ_annotated.liftover.mm10.vcf | sort -k1,1 -k2,2n >> 2012-0612-snps_only_FVBNJ_annotated.liftover.mm10.sorted.vcf
bgzip 2012-0612-snps_only_FVBNJ_annotated.liftover.mm10.sorted.vcf
tabix -p vcf 2012-0612-snps_only_FVBNJ_annotated.liftover.mm10.sorted.vcf.gz


# Freebayes regenotyping for my samples against FVB_NJ strain SNPs
cd /analysis/FVB_NJ

SAMPLES="
3NINACp5WGA
3NINACs1WGA
3NINACs2WGA
3NIp5WGA
3NIs1WGA
3NIs2WGA
4OnOffNACp5WGA
4OnOffNACs2WGA
4OnOffp5WGA
4OnOffs1WGA
control_lane3HC1_sequence
control_lane3HC3_sequence
control_lane3HC4_sequence
"

FREEBAYES="/soft/freebayes/bin/freebayes"
VCF_FILE="/mouse_vcf/2012-0612-snps_only_FVBNJ_annotated.liftover.mm10.sorted.vcf.gz"
REFGENOME="/datasets/refgenomes/mouse/Mus_musculus_mm10.fa"

for sample in $SAMPLES
do
    echo $sample
    $FREEBAYES -f $REFGENOME -@ ${VCF_FILE} /QC/mapped/${sample}.mapped.bam > ${sample}.mapped.var.vcf &
done

#### 129S1_SvImJ mouse strain for "control samples" from Chris

cd /analysis/129S1_SvImJ

SAMPLES="
control_lane3HC1_sequence
control_lane3HC3_sequence
control_lane3HC4_sequence
"

FREEBAYES="/soft/freebayes/bin/freebayes"
VCF_FILE="/mouse_vcf/SC_MOUSE_GENOMES.genotype.vcf.129S1_SvImJ.gz"
REFGENOME="/datasets/refgenomes/mouse/Mus_musculus_mm10.fa"

for sample in $SAMPLES
do
    echo $sample
    $FREEBAYES -f $REFGENOME -@ ${VCF_FILE} -l /QC/mapped/${sample}.mapped.bam | vt normalize -r $REFGENOME - > ${sample}.mapped.var.vcf &
done

### Quality Control

for f in *.mapped.var.vcf
do
    echo $f
    cat $f | grep -v "#" |  perl -ne 'my @all=split("\t", $_); print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]"; 
my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^AO/){my @res=split("=",$el); print "\t$res[1]\t";}if($el=~/^RO/){my @res2=split("=",$el); print "\t$res2[1]\n";}}' > ${f}.processed &
done


# 2. How many SNPs were identified out of total number of SNPs
# select the lines in which AO>0
cat 3NINACp5WGA.mapped.var.vcf.processed | awk '($7>0){print}' | grep -v ","  | wc -l

# 3. Estimate genetic background of the sample: heterozygous vs homozygous


##### Merge BAM files
cd /QC/mapped
FILES=`ls -1 *.mapped.bam |grep -v "control" `
samtools merge merged.bam -b $FILES &

# Genotyping the merged 1x file (11 * 0.1x)
cd /analysis/merged_bams
FREEBAYES="/soft/freebayes/bin/freebayes"
VCF_FILE="/mouse_vcf/2012-0612-snps_only_FVBNJ_annotated.liftover.mm10.sorted.vcf.gz"
REFGENOME="/datasets/refgenomes/mouse/Mus_musculus_mm10.fa"

$FREEBAYES -f $REFGENOME -@ ${VCF_FILE} -l /QC/mapped/merged.bam > merged.var.vcf &
bgzip merged.var.vcf; tabix -p vcf merged.var.vcf.gz
/soft/vt/vt normalize -r $REFGENOME merged.var.vcf.gz  > merged.var.vt.vcf

# Select only positions with more than 4 reads 
cat merged.var.vt.vcf | grep "#" | sed 's/\tunknown//g' > merged.var.vt.vcf.DP_filtered.processed
cat merged.var.vt.vcf | grep -v "#" |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^DP=/){my @res=split("=",$el); if ($res[1]>3) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t.\t.\n"}}}' >>  merged.var.vt.vcf.DP_filtered.processed &

cat merged.var.vt.vcf | grep "#" | sed 's/\tunknown//g' > merged.var.vt.vcf.DP_filtered.processed_AO
cat merged.var.vt.vcf | grep -v "#" |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^DP=/){my @res=split("=",$el); if ($res[1]>3) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t$all[7]\t.\n"}}}' >>  merged.var.vt.vcf.DP_filtered.processed_AO &

cat merged.var.vt.vcf.DP_filtered.processed_AO |  grep -v "#" |  perl -ne 'my @all=split("\t", $_); print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]"; 
my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^AO/){my @res=split("=",$el); print "\t$res[1]\t";}if($el=~/^RO/){my @res2=split("=",$el); print "\t$res2[1]\n";}}' > merged.var.vt.vcf.DP_filtered.processed_AO_filtered

cat merged.var.vt.vcf | grep "#" | sed 's/\tunknown//g' > merged.var.vt.vcf.DP_filtered.processed.final
cat merged.var.vt.vcf | grep -v "#" |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^DP=/){my @res=split("=",$el); if ($res[1]>3) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t$all[7]\t.\n"}}}' |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^AO=/){my @res=split("=",$el); if ($res[1]>0) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t.\t.\n"}}}' >> merged.var.vt.vcf.DP_filtered.processed.final



cd /analysis/merged_bams/individual_bams/

SAMPLES="
3NINACp5WGA
3NINACs1WGA
3NINACs2WGA
3NIp5WGA
3NIs1WGA
3NIs2WGA
4OnOffNACp5WGA
4OnOffNACs2WGA
4OnOffp5WGA
4OnOffs1WGA
"

FREEBAYES="/soft/freebayes/bin/freebayes"
VCF_FILE="/analysis/merged_bams/merged.var.vt.vcf.DP_filtered.processed.final.gz"
REFGENOME="/datasets/refgenomes/mouse/Mus_musculus_mm10.fa"

for sample in $SAMPLES
do
    echo $sample
    $FREEBAYES -f $REFGENOME -@ ${VCF_FILE} -l /QC/mapped/${sample}.mapped.bam > ${sample}.mapped.var.vcf &
done

# RECALL

touch stats.txt
for f in *.mapped.var.vcf
do
   echo $f
   cat $f | grep -v "#" |  perl -ne 'my @all=split("\t", $_); print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]"; 
my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^AO/){my @res=split("=",$el); print "\t$res[1]\t";}if($el=~/^RO/){my @res2=split("=",$el); print "\t$res2[1]\n";}}' > ${f}.processed &
   res=`cat ${f}.processed | awk '($7>0){print}' | wc -l`
   echo $f $res >> stats.txt
done


### Individual BAMs

cd /analysis/merged_bams/individual_bams/
for f in *.mapped.vsr.vcf
do

echo $f
cat merged.var.vt.vcf | grep -v "#" |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^DP=/){my @res=split("=",$el); if ($res[1]>3) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t.\t.\n"}}}' >>  merged.var.vt.vcf.DP_filtered.processed &

cat merged.var.vt.vcf | grep "#" | sed 's/\tunknown//g' > merged.var.vt.vcf.DP_filtered.processed_AO
cat merged.var.vt.vcf | grep -v "#" |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^DP=/){my @res=split("=",$el); if ($res[1]>3) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t$all[7]\t.\n"}}}' >>  merged.var.vt.vcf.DP_filtered.processed_AO &

cat merged.var.vt.vcf.DP_filtered.processed_AO |  grep -v "#" |  perl -ne 'my @all=split("\t", $_); print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]"; 
my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^AO/){my @res=split("=",$el); print "\t$res[1]\t";}if($el=~/^RO/){my @res2=split("=",$el); print "\t$res2[1]\n";}}' > merged.var.vt.vcf.DP_filtered.processed_AO_filtered

cat merged.var.vt.vcf | grep "#" | sed 's/\tunknown//g' > merged.var.vt.vcf.DP_filtered.processed.final
cat merged.var.vt.vcf | grep -v "#" |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^DP=/){my @res=split("=",$el); if ($res[1]>3) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t$all[7]\t.\n"}}}' |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^AO=/){my @res=split("=",$el); if ($res[1]>0) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t.\t.\n"}}}' >> merged.var.vt.vcf.DP_filtered.processed.final

# 2
cd /analysis/merged_bams

C_FILES=`ls -1 /QC/mapped/control*.mapped.bam `
samtools merge control_merged.bam $C_FILES &

# Genotyping the merged 1x file (11 * 0.1x)
cd /learning/analysis/merged_bams
FREEBAYES="/soft/freebayes/bin/freebayes"
VCF_FILE="/mouse_vcf/SC_MOUSE_GENOMES.genotype.vcf.129S1_SvImJ.gz"
REFGENOME="/datasets/refgenomes/mouse/Mus_musculus_mm10.fa"

samtools index control_merged.bam

$FREEBAYES -f $REFGENOME -@ ${VCF_FILE} -l /analysis/merged_bams/control_merged.bam > control_merged.var.vcf &
bgzip control_merged.var.vcf; tabix -p vcfcontrol_merged.var.vcf.gz
/g/korbel/rudneva/soft/vt/vt normalize -r $REFGENOME control_merged.var.vcf.gz  > control_merged.var.vt.vcf

# Select only positions with more than 4 reads 
cat merged.var.vt.vcf | grep "#" | sed 's/\tunknown//g' > merged.var.vt.vcf.DP_filtered.processed
cat merged.var.vt.vcf | grep -v "#" |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^DP=/){my @res=split("=",$el); if ($res[1]>3) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t.\t.\n"}}}' >>  merged.var.vt.vcf.DP_filtered.processed &

### FINAL
cd /analysis/merged_bams

# Select only positions with more than 0 reads 
cat merged.var.vt.vcf | grep "#" | sed 's/\tunknown//g' > final/merged.var.vt.vcf.DP_more_than_zero.processed
cat merged.var.vt.vcf | grep -v "#" |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^DP=/){my @res=split("=",$el); if ($res[1]>0) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t$all[7]\t$all[8]\n"}}}' >>  final/merged.var.vt.vcf.DP_more_than_zero.processed &

# Subsample from the merged control bam file to get an appropriate coverage
/software/bin/samtools-0.1.19 view -s 0.96 -b control_merged.bam > control_merged.subsampled.bam &
samtools index control_merged.subsampled.bam

bgzip final/merged.var.vt.vcf.DP_more_than_zero.processed
tabix -p vcf final/merged.var.vt.vcf.DP_more_than_zero.processed.gz

FREEBAYES="/soft/freebayes/bin/freebayes"
VCF_FILE="/mouse_vcf/2012-0612-snps_only_FVBNJ_annotated.liftover.mm10.sorted.vcf.gz"
REFGENOME="/datasets/refgenomes/mouse/Mus_musculus_mm10.fa"


$FREEBAYES -f $REFGENOME -@ ${VCF_FILE} -l /learning/analysis/merged_bams/control_merged.subsampled.bam > control_merged.subsample.var.vcf &

bgzip control_merged.subsample.var.vcf; tabix -p vcf control_merged.subsample.var.vcf.gz
/soft/vt/vt normalize -r $REFGENOME control_merged.subsample.var.vcf.gz  > control_merged.subsample.var.vt.vcf

# Select only positions with more than 0 reads 
cat control_merged.subsample.var.vt.vcf | grep "#" | sed 's/\tunknown//g' > final/control_merged.subsample.var.vt.vcf.DP_more_than_zero.processed
cat control_merged.subsample.var.vt.vcf | grep -v "#" |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^DP=/){my @res=split("=",$el); if ($res[1]>0) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t$all[7]\t$all[8]\n"}}}' >>  final/control_merged.subsample.var.vt.vcf.DP_more_than_zero.processed &

# Preparing for making genetic background plots
zcat /analysis/merged_bams/final/merged.var.vt.vcf.DP_more_than_zero.processed.gz |  grep -v "#" |  perl -ne 'my @all=split("\t", $_); print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]"; my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^AO/){my @res=split("=",$el); print "\t$res[1]\t";}if($el=~/^RO/){my @res2=split("=",$el); print "\t$res2[1]\n";}}' > /g/korbel/rudneva/martin/learning/analysis/merged_bams/final/merged.var.vt.vcf.DP_more_than_zero.processed_AO_filtered &

cat /analysis/merged_bams/final/control_merged.subsample.var.vt.vcf.DP_more_than_zero.processed | grep -v "#" |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^DP=/){my @res=split("=",$el); if ($res[1]>3) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t$all[7]\t.\n"}}}' |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^AO=/){my @res=split("=",$el); if ($res[1]>0) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t.\t.\n"}}}' > /g/korbel/rudneva/martin/learning/analysis/merged_bams/final/control_merged.subsample.var.vt.vcf.DP_more_than_zero.processed_AO_filtered &


# Coverage vs the recall
cd /analysis/merged_bams/individual_bams/

SAMPLES="
3NINACp5WGA
3NINACs1WGA
3NINACs2WGA
3NIp5WGA
3NIs1WGA
3NIs2WGA
4OnOffNACp5WGA
4OnOffNACs2WGA
4OnOffp5WGA
4OnOffs1WGA
"

FREEBAYES="/soft/freebayes/bin/freebayes"
VCF_FILE="/analysis/merged_bams/final/merged.var.vt.vcf.DP_more_than_zero.processed.gz"
REFGENOME="/datasets/refgenomes/mouse/Mus_musculus_mm10.fa"

for sample in $SAMPLES
do
    echo $sample
    $FREEBAYES -f $REFGENOME -@ ${VCF_FILE} -l /QC/mapped/${sample}.mapped.bam > ${sample}.mapped.var.vcf &
done

# RECALL

SAMPLES="
3NINACp5WGA
3NINACs1WGA
3NINACs2WGA
3NIp5WGA
3NIs1WGA
3NIs2WGA
4OnOffNACp5WGA
4OnOffNACs2WGA
4OnOffp5WGA
4OnOffs1WGA
"

for s in $SAMPLES
do
   f=${s}.mapped.var.vcf
   echo $f
   cat $f | grep -v "#" |  perl -ne 'my @all=split("\t", $_); print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]"; 
my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^AO/){my @res=split("=",$el); print "\t$res[1]\t";}if($el=~/^RO/){my @res2=split("=",$el); print "\t$res2[1]\n";}}' > ${f}.processed &
done


touch stats.txt
for s in $SAMPLES
do
   f=${s}.mapped.var.vcf
   echo $f
   res=`cat ${f}.processed | awk '($7>0){print}' | wc -l`
   cov=`cat /QC/mapped/seqstat/${s}.mapped.bam.r2.seqstats.pdf.stats | grep "Haploid sequencing coverage" | awk '{print $4}'`
   echo $s $res $cov >> stats.txt
done

merged_recall=`zcat /analysis/merged_bams/final/merged.var.vt.vcf.DP_more_than_zero.processed.gz | grep -v "#" | wc -l`
total=`zcat /mouse_vcf/2012-0612-snps_only_FVBNJ_annotated.liftover.mm10.sorted.vcf.gz | grep -v "#" | wc -l`
touch total_stats.txt
echo merged_recall total_SNVs >> total_stats.txt
echo $merged_recall $total >> total_stats.txt 

total_cov=`cat /QC/mapped/seqstat/merged.bam.r2.seqstats.pdf.stats | grep "Haploid sequencing" | awk '{print $4}'`
echo $total_cov > /analysis/merged_bams/individual_bams/total_cov.txt

# Recall for more than 4 reads
cd /analysis/merged_bams
cat control_merged.subsample.var.vt.vcf | grep "#" | sed 's/\tunknown//g' > final/control_merged.subsample.var.vt.vcf.DP_more_than_4.processed
cat control_merged.subsample.var.vt.vcf | grep -v "#" |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^DP=/){my @res=split("=",$el); if ($res[1]>4) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t$all[7]\t$all[8]\n"}}}' >>  final/control_merged.subsample.var.vt.vcf.DP_more_than_4.processed &
cat merged.var.vt.vcf | grep "#" | sed 's/\tunknown//g' > final/merged.var.vt.vcf.DP_more_than_4.processed
cat merged.var.vt.vcf | grep -v "#" |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^DP=/){my @res=split("=",$el); if ($res[1]>4) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t$all[7]\t$all[8]\n"}}}' >>  final/merged.var.vt.vcf.DP_more_than_4.processed &
  
# Recall for more than 2 reads
cd /analysis/merged_bams
cat control_merged.subsample.var.vt.vcf | grep "#" | sed 's/\tunknown//g' > final/control_merged.subsample.var.vt.vcf.DP_more_than_2.processed
cat control_merged.subsample.var.vt.vcf | grep -v "#" |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^DP=/){my @res=split("=",$el); if ($res[1]>2) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t$all[7]\t$all[8]\n"}}}' >>  final/control_merged.subsample.var.vt.vcf.DP_more_than_2.processed &
cat merged.var.vt.vcf | grep "#" | sed 's/\tunknown//g' > final/merged.var.vt.vcf.DP_more_than_2.processed
cat merged.var.vt.vcf | grep -v "#" |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^DP=/){my @res=split("=",$el); if ($res[1]>2) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t$all[7]\t$all[8]\n"}}}' >>  final/merged.var.vt.vcf.DP_more_than_2.processed &

# Repeat the analysis for "s" samples only
##### Merge BAM files
cd /QC/mapped
FILES=`ls -1 *.mapped.bam |grep -v "control"| grep "s"`
samtools merge merged.s.bam $FILES &
samtools index merged.s.bam 

# Genotyping the merged 1x file (11 * 0.1x)
cd /analysis/merged_bams
FREEBAYES="/soft/freebayes/bin/freebayes"
VCF_FILE="/mouse_vcf/2012-0612-snps_only_FVBNJ_annotated.liftover.mm10.sorted.vcf.gz"
REFGENOME="/datasets/refgenomes/mouse/Mus_musculus_mm10.fa"

$FREEBAYES -f $REFGENOME -@ ${VCF_FILE} -l /QC/mapped/merged.s.bam > merged.s.var.vcf &

bgzip merged.s.var.vcf; tabix -p vcf merged.s.var.vcf.gz
/soft/vt/vt normalize -r $REFGENOME merged.s.var.vcf.gz  > merged.s.var.vt.vcf

# Select only positions with more than 4 reads 
cat merged.s.var.vt.vcf | grep "#" | sed 's/\tunknown//g' > merged.s.var.vt.vcf.DP_more_than_4.processed
cat merged.s.var.vt.vcf | grep -v "#" |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^DP=/){my @res=split("=",$el); if ($res[1]>4) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t.\t.\n"}}}' >>  merged.s.var.vt.vcf.DP_more_than_4.processed &
# 646 positions

# Select only positions with more than 0 reads 
cat merged.s.var.vt.vcf | grep "#" | sed 's/\tunknown//g' > merged.s.var.vt.vcf.DP_more_than_zero.processed
cat merged.s.var.vt.vcf | grep -v "#" |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^DP=/){my @res=split("=",$el); if ($res[1]>0) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t.\t.\n"}}}' >>  merged.s.var.vt.vcf.DP_more_than_zero.processed; bgzip merged.s.var.vt.vcf.DP_more_than_zero.processed; tabix -p vcf merged.s.var.vt.vcf.DP_more_than_zero.processed.gz &

# Coverage vs the recall
cd /analysis/merged_bams/individual_bams/
bgzip merged.s.var.vt.vcf.DP_more_than_4.processed; tabix -p vcf merged.s.var.vt.vcf.DP_more_than_4.processed.gz &

SAMPLES=`ls -1 /QC/mapped/*.mapped.bam |grep -v "control"| sed 's/\//\t/g' | awk '{print $NF}' | sed 's/\.mapped\.bam//g' | grep "s"`

FREEBAYES="/soft/freebayes/bin/freebayes"
VCF_FILE="/analysis/merged_bams/merged.s.var.vt.vcf.DP_more_than_4.processed.gz"
REFGENOME="/datasets/refgenomes/mouse/Mus_musculus_mm10.fa"

for sample in $SAMPLES
do
    echo $sample
    $FREEBAYES -f $REFGENOME -@ ${VCF_FILE} -l /QC/mapped/${sample}.mapped.bam > new_more_than_DP_zero/${sample}.s.mapped.var.vcf &
done

# The same for more than zero reads
cd /analysis/merged_bams/individual_bams/

FREEBAYES="/soft/freebayes/bin/freebayes"
VCF_FILE="/analysis/merged_bams/merged.s.var.vt.vcf.DP_more_than_zero.processed.gz"
REFGENOME="/datasets/refgenomes/mouse/Mus_musculus_mm10.fa"

for sample in $SAMPLES
do
    echo $sample
    $FREEBAYES -f $REFGENOME -@ ${VCF_FILE} -l /QC/mapped/${sample}.mapped.bam > ${sample}.s.mapped.var.vcf &
done


##### To Do

# RECALL

SAMPLES=`ls -1 /QC/mapped/*.mapped.bam |grep -v "control"| sed 's/\//\t/g' | awk '{print $NF}' | sed 's/\.mapped\.bam//g' | grep "s"`
for s in $SAMPLES
do
   f=${s}.s.mapped.var.vcf
   echo $f
   cat $f | grep -v "#" |  perl -ne 'my @all=split("\t", $_); print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]"; 
my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^AO/){my @res=split("=",$el); print "\t$res[1]\t";}if($el=~/^RO/){my @res2=split("=",$el); print "\t$res2[1]\n";}}' > ${f}.processed &
done


touch stats.txt
for s in $SAMPLES
do
   f=${s}.s.mapped.var.vcf
   echo $f
   res=`cat ${f}.processed | awk '($7>0){print}' | wc -l`
   cov=`cat /QC/mapped/seqstat/${s}.mapped.bam.r2.seqstats.pdf.stats | grep "Haploid sequencing coverage" | awk '{print $4}'`
   echo $s $res $cov >> stats.txt
done

merged_recall=`zcat /analysis/merged_bams/final/merged.s.var.vt.vcf.DP_more_than_4.processed.gz | grep -v "#" | wc -l`
total=`zcat /mouse_vcf/2012-0612-snps_only_FVBNJ_annotated.liftover.mm10.sorted.vcf.gz | grep -v "#" | wc -l`
touch total_stats.txt
echo merged_recall total_SNVs >> total_stats.txt
echo $merged_recall $total >> total_stats.txt 

total_cov=`cat /QC/mapped/seqstat/merged.bam.r2.seqstats.pdf.stats | grep "Haploid sequencing" | awk '{print $4}'`
echo $total_cov > /analysis/merged_bams/individual_bams/total_cov.txt


### Analyze using Gingko

#bamToBed -i reads.bam > reads.bed

samples=`ls -1 *.bam | grep -v "control"`
for s in $samples; do echo $s; /software/bin/bamToBed -i $s > BED/$s.bed; done


### Distribution of base coverage

cd /analysis/merged_bams/final
cat merged.var.vt.vcf.DP_more_than_zero.processed | grep -v "#" | perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^DP=/){my @res=split("=",$el); if ($res[1]>0) {print "$all[0]\t$all[1]\t$res[1]\n"}}}' > merged.var.vt.vcf.DP_more_than_zero.processed.short
cat control_merged.subsample.var.vt.vcf.DP_more_than_zero.processed | grep -v "#" | perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^DP=/){my @res=split("=",$el); if ($res[1]>0) {print "$all[0]\t$all[1]\t$res[1]\n"}}}' > control_merged.subsample.var.vt.vcf.DP_more_than_zero.processed.short



########################################
##### Comparing different mouse strains
########################################

VCF_DIR="/g/korbel/rudneva/martin/learning/mouse_vcf/REL-1505-SNPs_Indels"
FREEBAYES="/g/korbel/rudneva/soft/freebayes/bin/freebayes"
REFGENOME="/g/korbel/shared/datasets/refgenomes/mouse/Mus_musculus_mm10.fa"


/g/korbel/rudneva/soft/vt/vt normalize -r $REFGENOME ${VCF_DIR}/BALB_cJ.mgp.v5.snps.dbSNP142.vcf.gz  > ${VCF_DIR}/BALB_cJ.mgp.v5.snps.dbSNP142.var.vcf
cat ${VCF_DIR}/BALB_cJ.mgp.v5.snps.dbSNP142.var.vcf | grep "#" | sed 's/\tBALB_cJ//g' > ${VCF_DIR}/BALB_cJ.mgp.v5.snps.dbSNP142.var.processed.vcf
cat ${VCF_DIR}/BALB_cJ.mgp.v5.snps.dbSNP142.var.vcf | grep -v "#" | awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t.\t.\t."}' >>${VCF_DIR}/BALB_cJ.mgp.v5.snps.dbSNP142.var.processed.vcf &
bgzip ${VCF_DIR}/BALB_cJ.mgp.v5.snps.dbSNP142.var.processed.vcf; tabix -p vcf ${VCF_DIR}/BALB_cJ.mgp.v5.snps.dbSNP142.var.processed.vcf.gz


/g/korbel/rudneva/soft/vt/vt normalize -r $REFGENOME ${VCF_DIR}/FVB_NJ.mgp.v5.snps.dbSNP142.vcf.gz  > ${VCF_DIR}/FVB_NJ.mgp.v5.snps.dbSNP142.var.vcf
cat ${VCF_DIR}/FVB_NJ.mgp.v5.snps.dbSNP142.var.vcf | grep "#" | sed 's/\tFVB_NJ//g' > ${VCF_DIR}/FVB_NJ.mgp.v5.snps.dbSNP142.var.processed.vcf
cat ${VCF_DIR}/FVB_NJ.mgp.v5.snps.dbSNP142.var.vcf | grep -v "#" | awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t.\t.\t."}' >> ${VCF_DIR}/FVB_NJ.mgp.v5.snps.dbSNP142.var.processed.vcf &
bgzip ${VCF_DIR}/FVB_NJ.mgp.v5.snps.dbSNP142.var.processed.vcf; tabix -p vcf ${VCF_DIR}/FVB_NJ.mgp.v5.snps.dbSNP142.var.processed.vcf.gz



FVB="FVB_NJ.mgp.v5.snps.dbSNP142.var.processed.vcf.gz"
BALB="BALB_cJ.mgp.v5.snps.dbSNP142.var.processed.vcf.gz"

cd /g/korbel/rudneva/martin/learning/analysis/merged_bams

VCF_FILE1=${VCF_DIR}/${FVB}
VCF_FILE2=${VCF_DIR}/${BALB}

mkdir compare

# samples
$FREEBAYES -f $REFGENOME -@ ${VCF_FILE1} -l /g/korbel/rudneva/martin/learning/QC/mapped/merged.s.bam > compare/merged.FVB.var.vcf &
$FREEBAYES -f $REFGENOME -@ ${VCF_FILE2} -l /g/korbel/rudneva/martin/learning/QC/mapped/merged.s.bam > compare/merged.BALB.var.vcf &

# controls
$FREEBAYES -f $REFGENOME -@ ${VCF_FILE1} -l /g/korbel/rudneva/martin/learning/analysis/merged_bams/control_merged.subsampled.bam > compare/control_merged.subsample.FVB.var.vcf &
$FREEBAYES -f $REFGENOME -@ ${VCF_FILE2} -l /g/korbel/rudneva/martin/learning/analysis/merged_bams/control_merged.subsampled.bam > compare/control_merged.subsample.BALB.var.vcf &


# To Do:

cd compare/

bgzip merged.FVB.var.vcf; tabix -p vcf merged.FVB.var.vcf.gz
/g/korbel/rudneva/soft/vt/vt normalize -r $REFGENOME merged.FVB.var.vcf.gz  > merged.FVB.var.vt.vcf

bgzip merged.BALB.var.vcf; tabix -p vcf merged.BALB.var.vcf.gz
/g/korbel/rudneva/soft/vt/vt normalize -r $REFGENOME merged.BALB.var.vcf.gz  > merged.BALB.var.vt.vcf

bgzip control_merged.subsample.FVB.var.vcf; tabix -p vcf control_merged.subsample.FVB.var.vcf.gz
/g/korbel/rudneva/soft/vt/vt normalize -r $REFGENOME control_merged.subsample.FVB.var.vcf.gz  > control_merged.subsample.FVB.var.vt.vcf

bgzip control_merged.subsample.BALB.var.vcf; tabix -p vcf control_merged.subsample.BALB.var.vcf.gz
/g/korbel/rudneva/soft/vt/vt normalize -r $REFGENOME control_merged.subsample.BALB.var.vcf.gz  > control_merged.subsample.BALB.var.vt.vcf


# Select only positions with more than 0 reads 
cat merged.FVB.var.vt.vcf | grep "#" | sed 's/\tunknown//g' > merged.FVB.var.vt.vcf.DP_more_than_zero.processed
cat merged.FVB.var.vt.vcf | grep -v "#" |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^DP=/){my @res=split("=",$el); if ($res[1]>0) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t.\t.\n"}}}' >>  merged.FVB.var.vt.vcf.DP_more_than_zero.processed; bgzip merged.FVB.var.vt.vcf.DP_more_than_zero.processed; tabix -p vcf merged.FVB.var.vt.vcf.DP_more_than_zero.processed.gz &

cat merged.BALB.var.vt.vcf | grep "#" | sed 's/\tunknown//g' > merged.BALB.var.vt.vcf.DP_more_than_zero.processed
cat merged.BALB.var.vt.vcf | grep -v "#" |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^DP=/){my @res=split("=",$el); if ($res[1]>0) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t.\t.\n"}}}' >>  merged.BALB.var.vt.vcf.DP_more_than_zero.processed; bgzip merged.BALB.var.vt.vcf.DP_more_than_zero.processed; tabix -p vcf merged.BALB.var.vt.vcf.DP_more_than_zero.processed.gz &


cat control_merged.subsample.FVB.var.vt.vcf | grep "#" | sed 's/\tunknown//g' > control_merged.FVB.var.vt.vcf.DP_more_than_zero.processed
cat control_merged.subsample.FVB.var.vt.vcf | grep -v "#" |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^DP=/){my @res=split("=",$el); if ($res[1]>0) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t.\t.\n"}}}' >>  control_merged.FVB.var.vt.vcf.DP_more_than_zero.processed; bgzip control_merged.FVB.var.vt.vcf.DP_more_than_zero.processed; tabix -p vcf control_merged.FVB.var.vt.vcf.DP_more_than_zero.processed.gz &

cat control_merged.BALB.var.vt.vcf | grep "#" | sed 's/\tunknown//g' > control_merged.BALB.var.vt.vcf.DP_more_than_zero.processed
cat control_merged.BALB.var.vt.vcf | grep -v "#" |  perl -ne 'my @all=split("\t", $_); my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^DP=/){my @res=split("=",$el); if ($res[1]>0) {print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t.\t.\n"}}}' >>  control_merged.BALB.var.vt.vcf.DP_more_than_zero.processed; bgzip control_merged.BALB.var.vt.vcf.DP_more_than_zero.processed; tabix -p vcf control_merged.BALB.var.vt.vcf.DP_more_than_zero.processed.gz &


# Recall 1.2 m shared SNPs on control samples


FREEBAYES="/g/korbel/rudneva/soft/freebayes/bin/freebayes"
VCF_FILE="/g/korbel/rudneva/martin/learning/analysis/merged_bams/final/merged.var.vt.vcf.DP_more_than_zero.processed.gz"
REFGENOME="/g/korbel/shared/datasets/refgenomes/mouse/Mus_musculus_mm10.fa"

for sample in $SAMPLES
do
    echo $sample
    $FREEBAYES -f $REFGENOME -@ ${VCF_FILE} -l /g/korbel/rudneva/martin/learning/QC/mapped/${sample}.mapped.bam > ${sample}.mapped.var.vcf &
done

# To Do:

for s in $SAMPLES
do
   f=${s}.mapped.var.vcf
   echo $f
   cat $f | grep -v "#" |  perl -ne 'my @all=split("\t", $_); print "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]"; 
my @info=split(";",$all[7]); foreach my $el (@info){if ($el=~/^AO/){my @res=split("=",$el); print "\t$res[1]\t";}if($el=~/^RO/){my @res2=split("=",$el); print "\t$res2[1]\n";}}' > ${f}.processed &
done
touch control_stats.txt
for s in $SAMPLES
do
   f=${s}.mapped.var.vcf
   echo $f
   res=`cat ${f}.processed | awk '($7>0){print}' | wc -l`
   cov=`cat /g/korbel/rudneva/martin/learning/QC/mapped/seqstat/${s}.mapped.bam.r2.seqstats.pdf.stats | grep "Haploid sequencing coverage" | awk '{print $4}'`
   echo $s $res $cov >> control_stats.txt
done
