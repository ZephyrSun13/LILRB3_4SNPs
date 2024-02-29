
Wk_dir=$(pwd)
Bin_dir=$(pwd)

mkdir -p $Wk_dir

mkdir -p $Wk_dir/Separate

awk -v Wk_dir="$Wk_dir" '{

        printf("echo Start ");

	printf("&& module load star/2.6.1d && mkdir -p %s/Separate/%s && STAR --genomeDir %s/1.MaskGenome/StarDefault.all/ --sjdbGTFfile Homo_sapiens.GRCh37.75.gtf --readFilesIn %s %s --readFilesCommand zcat --runThreadN 12 --outReadsUnmapped Fastx --outSAMstrandField intronMotif --outSAMmapqUnique 60 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix %s/Separate/%s/%s", Wk_dir, Wk_dir, $1, $2, $3, Wk_dir, $1, $1);

	printf(" && module load samtools/1.9 && samtools index %s/Separate/%s/%sAligned.sortedByCoord.out.bam", Wk_dir, $1, $1);

	printf(" && module load picard/1.93 && module load java/1.8.0_211 && java -jar /hpc/packages/minerva-common/picard/1.93/bin/MarkDuplicates.jar I=%s/Separate/%s/%sAligned.sortedByCoord.out.bam O=%s/Separate/%s/%s_aln.dedup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=%s/Separate/%s/output.metrics", Wk_dir, $1, $1, Wk_dir, $1, $1, Wk_dir, $1);
	
	printf("\n");

}' Sample.lst > run_Align_all.sh

nohup ${Bin_dir}/qsub-sge.pl --queue premium --convert no --pro_code acc_xxxx --reqsub --jobprefix StarPair --resource 40000 --time 3600 --verbose run_Align_all.sh &

