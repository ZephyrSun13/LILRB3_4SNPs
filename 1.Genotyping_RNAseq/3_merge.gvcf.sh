
set -e

Wk_dir=$(pwd)
Db_dir=$(pwd)

ml gatk/4.2.0.0

ls ${Wk_dir}/gatk_res/*/*.g.vcf.gz | awk '{split($1, tt, "/"); printf("%s\t%s\n", tt[10], $1)}' > cohort.sample_map

awk '{split($1, tt, "="); gsub(",length", "", tt[3]); gsub(">", "", tt[4]);printf("%s\t1\t%s\n", tt[3], tt[4])}' contig.lst > Interval.lst.bed

rm -rf GBase

gatk --java-options "-Xmx4g -Xms4g" \
       GenomicsDBImport \
       --genomicsdb-workspace-path ${Wk_dir}/GBase \
       --intervals Interval.lst.bed \
       --sample-name-map cohort.sample_map \
       --tmp-dir /sc/arion/scratch/sunz04/Sun \
       --reader-threads 1
  
gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R ${Db_dir}/GATK/human_g1k_v37.sorted.fasta \
   -V gendb://${Wk_dir}/GBase/ \
   -O Merged.vcf.gz


gatk --java-options "-Xmx8g" VariantFiltration -R ${Db_dir}/GATK/human_g1k_v37.sorted.fasta \
  -V Merged.vcf.gz \
  -O Merged.filter.vcf.gz \
  --filter-expression "FS > 30.0" \
  --filter-name "FS" \
  --filter-expression "QD < 2.0" \
  --filter-name "QD"

ml bcftools/1.9

        bcftools filter -i 'FILTER="PASS"' Merged.filter.vcf.gz > Merged.filter.vcf.gz.clean.vcf

	grep '#CHROM' Merged.filter.vcf.gz.clean.vcf | cut -f 2,10- > Merged.filter.vcf.gz.clean.vcf.tab

bcftools query -f '%CHROM %POS %REF %ALT[\t%GT]\n' Merged.filter.vcf.gz.clean.vcf >> Merged.filter.vcf.gz.clean.vcf.tab

