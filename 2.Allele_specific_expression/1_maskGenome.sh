
## Mask exonic SNPs as N to avoid reference biased alignment

Db_dir=${pwd}

ml BEDTools/2.29.0

## change 0-based coordinates to 1-based coordinates

cat ../0.Data/GOCAR.Exonic.snp.bed | awk '{printf("%s\t%s\t%s\n", $1, $2-1, $3)}' > GOCAR.Exonic.snp.bed 

cat ../0.Data/GOCAR.RNA.Exonic.snp.bed | awk '{printf("%s\t%s\t%s\n", $1, $2-1, $3)}' >> GOCAR.Exonic.snp.bed

sort GOCAR.Exonic.snp.bed | uniq > tt && mv tt GOCAR.Exonic.snp.bed

bedtools maskfasta -fi ${Db_dir}/UCSC/hg19/Homo_sapiens.GRCh37.dna.primary_assembly.fa -bed GOCAR.Exonic.snp.bed -fo Homo_sapiens.GRCh37.dna.primary_assembly.exonMasked.all.fa


## build STAR reference and picard dictionary for masked genome

module load star/2.6.1d

mkdir -p StarDefault.all

STAR --runThreadN 1 --runMode genomeGenerate --genomeDir StarDefault --genomeFastaFiles Homo_sapiens.GRCh37.dna.primary_assembly.exonMasked.all.fa --sjdbGTFfile ${Db_dir}/UCSC/hg19/Homo_sapiens.GRCh37.75.gtf --sjdbOverhang 100 


ml picard/2.22.3

java -jar $PICARD CreateSequenceDictionary \
      R=Homo_sapiens.GRCh37.dna.primary_assembly.exonMasked.fa \
      O=Homo_sapiens.GRCh37.dna.primary_assembly.exonMasked.dict

