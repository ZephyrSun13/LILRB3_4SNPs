
Bin_dir=$(pwd)
Wk_dir=$(pwd)
Db_dir=$(pwd)

ls ${Wk_dir}/Separate/*/*dedup.bam | awk '{

	printf("echo Start")

	printf(" && ml picard/2.22.3 && ml java/1.8.0_211 && java -jar $PICARD AddOrReplaceReadGroups I=%s O=%s.RG.bam RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20", $1, $1)

	printf(" && ml gatk/4.2.0.0 && gatk ASEReadCounter -R Homo_sapiens.GRCh37.dna.primary_assembly.fa -I %s.RG.bam -V GOCAR.Exonic.GT.snp.bed.vcf.recode.corr.dedup.vcf.gz -O %s.RG.bam.ase.output.table\n", $1, $1)}' > ase_all.sh

nohup ${Bin_dir}/qsub-sge.pl --queue premium --convert no --pro_code acc_xxx --jobprefix ASE --resource 20000 --time 800 --maxjob 300 --verbose ase_all.sh &

