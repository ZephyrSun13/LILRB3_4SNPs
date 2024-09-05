
ml bcftools/1.9


grep '#CHROM' LILR.WES.vcf | cut -f 2,10- > LILR.WES.vcf.tab 

bcftools query -f '%CHROM %POS %REF %ALT[\t%GT]\n' LILR.WES.vcf >> LILR.WES.vcf.tab 

grep '#CHROM' LILR.WES.vcf | cut -f 2,10- > LILR.WES.vcf.GQ.tab

bcftools query -f '%CHROM %POS %REF %ALT[\t%GQ]\n' LILR.WES.vcf >> LILR.WES.vcf.GQ.tab


ml R/4.0.3

Rscript stat.r


head -1 TabMerge.xls > TabMerge.AA.xls

awk -F "\t" '($9 == "African American" || $9 == "AFRICAN AMERICAN (BLACK)"){print $0}' TabMerge.xls  >> TabMerge.AA.xls

head -1 TabMerge.xls > TabMerge.AAHis.xls

awk -F "\t" '($9 == "African American" || $9 == "AFRICAN AMERICAN (BLACK)" || $9 == "Hispanic"){print $0}' TabMerge.xls  >> TabMerge.AAHis.xls


Rscript genoClassify.r TabMerge.AAHis.xls

Rscript genoClassify.r TabMerge.AA.xls


## Diagnosis

awk '(NR==FNR){tt[$1] = 1; next}{split($0, tt2, "|"); if(tt[tt2[2]] != ""){print $0}}' TabMerge.AA.xls Pheno2/Encounter_Diagnosis.txt > Encounter_Diagnosis.AA.txt

awk -F "|" '{print $3}' Encounter_Diagnosis.AA.txt | sort | uniq -c | awk '{printf("%s\t%s\n", $2, $1)}' | sort -nr -k 2,2 > Encounter_Diagnosis.AA.txt.code

awk -F "|" '{printf("%s\t%s\n", $3, $5)}' Encounter_Diagnosis.AA.txt | sort | uniq > Code.map

awk -F "|" '($3=="Z94.0" || $3 == "V42.0"){print $2}' Encounter_Diagnosis.AA.txt | sort | uniq -c | awk '{printf("%s\t%s\n", $2, $1)}' | sort -nr -k 2,2 > Encounter_Diagnosis.AA.txt.code.KTR 

awk -F "|" '$2=="SINAI_35135_AB42709878"' Encounter_Diagnosis.AA.txt | sed 's/|/\t/g' > SINAI_35135_AB42709878.xls


Rscript biomeDiaGN.r LILRAll/WES.rds AA CodeTab.all.rds.filt.rds Dom "19_54217126_C_T" No

Rscript biomeDiaGN.r LILRAll/WES.rds AA CodeTab.all.rds.filt.rds Add "19_54217126_C_T" No

## Blood work

awk '(NR==FNR){tt[$1] = 1; next}{split($0, tt2, "|"); if(tt[tt2[2]] != ""){print $0}}' TabMerge.AA.xls Pheno2/Order_results.txt > Order_results.AA.txt

zcat Order_results.AA.txt.gz | awk -F '|' '{print $4}' | sort | uniq -c | sort -k '1,1' -nr > BloodWorks

awk '($1 >= 1000){$1=""; print $0}' BloodWorks | sed 's/^ //g' | awk '{tt = $0; gsub(" ", "_", tt); gsub("-", "_", tt); gsub("%", "P", tt); gsub("#", "C", tt); gsub(",", "", tt); gsub("\\.", "", tt); gsub("\\(", "_", tt); gsub("\\)", "_", tt); gsub("\\/", "_", tt); printf("%s\t%s\n", $0, tt)}' > BloodWorks.sele


## AA Hispanic eGFR

Rscript eGFR.suv.Steve.r 15 "LILRRisk" "Transplant;ESRD;CKD;AKI" AA "WT;LILRB3-4SNPs" 90 365 200

Rscript eGFR.suv.Steve.r 15 "APOL1Geno" "Transplant;ESRD;CKD;AKI" AA "WT;APOL1-single;APOL1-double" 90 365 200

Rscript eGFR.suv.Steve.r 15 "Combine2" "Transplant;ESRD;CKD;AKI" AA "WT;APOL1-double;LILRB3-4SNPs;Combined" 90 365 200

