
gsutil -u $GOOGLE_PROJECT ls gs://fc-aou-datasets-controlled

gsutil -u $GOOGLE_PROJECT ls gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/exome_v7.1/vcf/ | less -S

gsutil -u $GOOGLE_PROJECT du -h gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/exome_v7.1/vcf/exome.chr19.vcf.bgz

gsutil -u $GOOGLE_PROJECT du -h gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/exome_v7.1/vcf/exome.chr22.vcf.bgz


gsutil -u $GOOGLE_PROJECT -m cp gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/exome_v7.1/vcf/exome.chr19.vcf.bgz .

gsutil -u $GOOGLE_PROJECT -m cp gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/exome_v7.1/vcf/exome.chr19.vcf.bgz.tbi .

gsutil -u $GOOGLE_PROJECT -m cp gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/exome_v7.1/vcf/exome.chr22.vcf.bgz .

gsutil -u $GOOGLE_PROJECT -m cp gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/exome_v7.1/vcf/exome.chr22.vcf.bgz.tbi .


bcftools view -R LILR.hg38.region exome.chr19.vcf.bgz > LILR.WES.all.vcf

bcftools view -H -R LILR.hg38.region exome.chr22.vcf.bgz >> LILR.WES.all.vcf


grep '#CHROM' LILR.WES.all.vcf | cut -f 2,10- > LILR.WES.vcf.tab

bcftools query -f '%CHROM %POS %REF %ALT[\t%GT]\n' LILR.WES.all.vcf >> LILR.WES.vcf.tab

grep '#CHROM' LILR.WES.all.vcf | cut -f 2,10- > LILR.WES.vcf.GQ.tab

bcftools query -f '%CHROM %POS %REF %ALT[\t%GQ]\n' LILR.WES.all.vcf >> LILR.WES.vcf.GQ.tab

