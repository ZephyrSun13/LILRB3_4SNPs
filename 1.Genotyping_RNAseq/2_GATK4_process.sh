#!/bin/bash -l
## working directory

wk_dir=$(pwd)
db_dir=$(pwd)
bin_dir=$(pwd)

mkdir -p ${wk_dir}/gatk_lsf
mkdir -p ${wk_dir}/gatk_log
mkdir -p ${wk_dir}/gatk_err
mkdir -p ${wk_dir}/gatk_res

FilesToRun=${wk_dir}/Sample.lst

picard=${bin_dir}/picard/1.93/bin
fasta=${db_dir}/GATK/human_g1k_v37.sorted.fasta
vcfmills=${db_dir}/GATK/Mills_and_1000G_gold_standard.indels.b37.vcf
vcf1000g=${db_dir}/GATK/1000G_phase1.indels.b37.sorted.vcf
vcfdbsnp=${db_dir}/GATK/dbsnp_138.b37.excluding_sites_after_129.sorted.vcf

while read line ; do
    sample_id="$(echo ${line} | cut -d' ' -f1)"

    sample_dir=${wk_dir}/gatk_res/${sample_id}
    mkdir -p ${sample_dir}

    echo ${sample_id}
    echo ${sample_dir}

echo "

set -e

#!/bin/bash
#BSUB -J gatk_${sample_id}
#BSUB -eo ${wk_dir}/gatk_err/${sample_id}.e
#BSUB -oo ${wk_dir}/gatk_log/${sample_id}.o
#BSUB -q premium
#BSUB -P acc_***
#BSUB -W 24:00
#BSUB -n 4
#BSUB -R "rusage[mem=8000]" 
#BSUB -R "span[hosts=1]"
#BSUB -L /bin/bash

mkdir -p ${sample_dir}
cd ${sample_dir}

module purge
module load picard/1.93 samtools/1.9 R
ml gatk/4.2.0.0

## remove dups
java -jar ${picard}/MarkDuplicates.jar \
I=${wk_dir}/star_res/${sample_id}/${sample_id}_aln.bam \
O=${sample_dir}/${sample_id}_aln.dedup.bam \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT \
M=${sample_dir}/output.metrics 

## need to sort bam in karyotypic order for GATK input
java -jar ${picard}/ReorderSam.jar \
I=${sample_dir}/${sample_id}_aln.dedup.bam \
O=${sample_dir}/${sample_id}_aln.dedup.resort.bam \
REFERENCE=${fasta} \
CREATE_INDEX=true

## add group information
java -jar ${picard}/AddOrReplaceReadGroups.jar I=${sample_dir}/${sample_id}_aln.dedup.resort.bam  O=${sample_dir}/${sample_id}_aln.dedup.resort.group.bam  RGLB=Lane1 RGPL=illumina RGPU=unit1 RGSM=${sample_id}

## index bam file
java -jar ${picard}/BuildBamIndex.jar I=${sample_dir}/${sample_id}_aln.dedup.resort.group.bam

## split reads

gatk SplitNCigarReads \
-R ${fasta} \
-I ${sample_dir}/${sample_id}_aln.dedup.resort.group.bam \
-O ${sample_dir}/${sample_id}_aln.dedup.split.bam

## Base recalibration

gatk BaseRecalibrator \
   -I ${sample_dir}/${sample_id}_aln.dedup.split.bam \
   -R ${fasta} \
   --known-sites ${vcfmills} \
   --known-sites ${vcf1000g} \
   --known-sites ${vcfdbsnp} \
   -O ${sample_dir}/recal_data.table

gatk ApplyBQSR \
   -R ${fasta} \
   -I ${sample_dir}/${sample_id}_aln.dedup.split.bam \
   --bqsr-recal-file ${sample_dir}/recal_data.table \
   -O ${sample_dir}/${sample_id}_aln.dedup.split.realn.recal.bam
 
gatk AnalyzeCovariates \
     -bqsr ${sample_dir}/recal_data.table \
     -plots ${sample_dir}/AnalyzeCovariates.pdf

## final bam output: ${sample_id}_aln.dedup.split.realn.recal.bam
rm -f ${sample_dir}/${sample_id}_aln.dedup.bam
rm -f ${sample_dir}/${sample_id}_aln.dedup.resort.bam
rm -f ${sample_dir}/${sample_id}_aln.dedup.resort.group.bam
rm -f ${sample_dir}/${sample_id}_aln.dedup.split.bam

## call variants ( phred-scaled confidence at calling >= 20.0 and emitting at >= 20.0; HQ is >= 30.0);  , 
## filter variants ( Phred-scaled p-value > 30.0, Quality by Depth >20 ). 

gatk --java-options "-Xmx8g" HaplotypeCaller  \
   -R ${fasta} \
   -I ${sample_dir}/${sample_id}_aln.dedup.split.realn.recal.bam \
   -O ${sample_dir}/${sample_id}_variants.g.vcf.gz \
   -ERC GVCF

gatk --java-options "-Xmx8g" HaplotypeCaller  \
   -R ${fasta} \
   -I ${sample_dir}/${sample_id}_aln.dedup.split.realn.recal.bam \
   -O ${sample_dir}/${sample_id}_variants.vcf.gz

gatk --java-options "-Xmx8g" VariantFiltration  \
   -R ${fasta} \
   -V ${sample_dir}/${sample_id}_variants.vcf.gz \
   -O ${sample_dir}/${sample_id}_variants.filter.vcf.gz \
   --filter-expression \"FS > 30.0\" \
   --filter-name \"FS\" \
   --filter-expression \"QD < 2.0\" \
   --filter-name \"QD\"

echo "\"analysis completed.\""

" > ${wk_dir}/gatk_lsf/${sample_id}.lsf

done < ${FilesToRun}

echo "analysis completed."

ls ${wk_dir}/gatk_lsf/*.lsf | awk '{printf("bsub <%s\n", $1)}' > bsub.gatk.sh

