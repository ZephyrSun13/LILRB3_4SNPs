#!/bin/bash -l

wk_dir=$(pwd)

rm -rf ${wk_dir}/star_lsf

mkdir -p ${wk_dir}/star_lsf
mkdir -p ${wk_dir}/star_log
mkdir -p ${wk_dir}/star_err
mkdir -p ${wk_dir}/star_res

FilesToRun=${wk_dir}/Sample.lst
Genome=${GenomeDir}
stargtf=${StarGTFDir}


while read line ; do
    sample_id="$(echo ${line} | cut -d' ' -f1)"
    fq1_dir="$(echo ${line} | cut -d' ' -f2)"
    fq2_dir="$(echo ${line} | cut -d' ' -f3)"
    
    sample_dir=${wk_dir}/star_res/${sample_id}
    mkdir -p ${sample_dir}

    ## code need to be adjusted for single end reads
    echo ${sample_id}
    echo ${fq1_dir}
    echo ${fq2_dir}

echo "
#!/bin/bash
#BSUB -J star_${sample_id}
#BSUB -eo ${wk_dir}/star_err/${sample_id}.e
#BSUB -oo ${wk_dir}/star_log/${sample_id}.o
#BSUB -q long 
#BSUB -P acc_XXX
#BSUB -W 150:00
#BSUB -n 1
#BSUB -R "rusage[mem=40000]" 
#BSUB -R "span[ptile=1]"
#BSUB -L /bin/bash

##### STAR needs >32GB RAM, ~40 to be safe, and spanned across just one host.
##### this sets us up with 12*3.5=42GB

mkdir -p ${sample_dir}
cd ${sample_dir}

module purge
module load star/2.6.1d samtools/1.9

STAR --genomeDir ${Genome} \
--sjdbGTFfile ${stargtf} \
--readFilesIn ${fq1_dir} ${fq2_dir} \
--twopassMode Basic \
--readFilesCommand zcat \
--runThreadN 12 \
--outReadsUnmapped Fastx \
--chimSegmentMin 15 \
--chimJunctionOverhangMin 15 \
--outSAMstrandField intronMotif \
--outSAMmapqUnique 60 \
--outStd BAM_SortedByCoordinate \
--outSAMtype BAM SortedByCoordinate >${sample_id}_aln.bam

samtools index ${sample_id}_aln.bam

## counts

module load subread/1.6.3

featureCounts -T 12 -t exon -g gene_id -a ${stargtf} -o ${sample_id}_counts.txt ${sample_id}_aln.bam

echo "\"analysis completed.\""

" > ${wk_dir}/star_lsf/${sample_id}.lsf

done < ${FilesToRun}

echo "analysis completed."

ls star_lsf/*.lsf | awk '{printf("bsub <%s\n", $1)}' > bsub.star.sh

