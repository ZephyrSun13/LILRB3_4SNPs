
Wk_dir=$(pwd)
Bin_dir=$(pwd)

Ref=${Wk_dir}/ref.fa

mkdir -p $Wk_dir

awk -v Wk_dir="$Wk_dir" -v Ref="$Ref" '{

	printf("echo Start && mkdir -p %s/%s", Wk_dir, $1);

	printf(" && ml bwa/0.7.15 && ml samtools/1.9 && bwa mem -t 5 %s %s %s | samtools view -bS - > %s/%s/%s.bam && samtools sort -o %s/%s/%s.sorted.bam %s/%s/%s.bam && samtools index %s/%s/%s.sorted.bam", Ref, $2, $3, Wk_dir, $1, $1, Wk_dir, $1, $1, Wk_dir, $1, $1, Wk_dir, $1, $1);

	printf(" && ml samtools/1.9 && samtools stats %s/%s/%s.sorted.bam > %s/%s/%s.sorted.bam.stats && samtools flagstat %s/%s/%s.sorted.bam > %s/%s/%s.sorted.bam.flagstat", Wk_dir, $1, $1, Wk_dir, $1, $1, Wk_dir, $1, $1, Wk_dir, $1, $1);

	printf(" && ml bcftools/1.9 && bcftools mpileup --max-depth 300000 -f %s %s/%s/%s.sorted.bam | bcftools call -m -Ov -p 1 -o %s/%s/%s.sorted.bam.vcf", Ref, Wk_dir, $1, $1, Wk_dir, $1, $1);

	printf("\n");

}' Sam > All.sh

nohup ${Bin_dir}/qsub-sge.pl --queue premium --convert no --pro_code acc_*** --jobprefix Target --resource 50000 --time 1440 --verbose All.sh &

