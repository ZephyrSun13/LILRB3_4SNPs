
Wk_dir=$(pwd)

ml bcftools/1.9

for vv in ${Wk_dir}/LILR/Target/*/*.vcf
do

	bgzip $vv

	tabix ${vv}.gz

done

ls ${Wk_dir}/LILR/Target/*/*.vcf.gz > VCF.lst

bcftools merge -l VCF.lst > ${Wk_dir}/Target.merge.vcf 


for vv in ${Wk_dir}/LILR/Target/*/*.vcf.gz
do

	zgrep '#CHROM' $vv | cut -f 2,10- > ${vv}.tab 

	bcftools query -f '%CHROM %POS %REF %ALT[\t%INFO/DP4][\t%GT]\n' $vv >> ${vv}.tab 

done


ls ${Wk_dir}/Target/*/*.tab > Tabs

ml R/4.0.3

Rscript ${Wk_dir}/summary.r

