## LILRB3_4SNPs

Codes to run LILRB3_4SNPs and generate figures for publication[1].

### 1.Genotyping_RNAseq: Genotyping with RNAseq sequencing data.
    Data can be downloaded from GSE252272.
    
    Dependent softwares: star/2.6.1d; samtools/1.9; subread/1.6.3; picard/1.93; gatk/4.2.0.0; R/4.0.3; bcftools/1.9; 

    Running steps:
        1. Samples need to be configured in the format showed in "Sample.lst", single-ended reads can eliminate the third column. 
        2. LSF enrioment lsf files for STAR alignment for each sample can be generated in "star_lsf" sub-folder by running "sh 1_star_mapping_2pass.pair.sh". The LSF parameters need to be modified as needed. The alignment results can be found in "star_res" sub-folder after running "sh bsub.star.sh".
        3. LSF enrioment lsf files for GATK best practice of RNAseq for each sample can be generated in "gatk_lsf" sub-folder by running "sh 2_GATK4_process.sh". The LSF parameters need to be modified as needed. The SNP calling results can be found in "gatk_res" sub-folder after running "sh bsub.gatk.sh".
        4. The gvcf files for each sample can be further merged into one vcf file and filtered by running "sh 3_merge.gvcf.sh".

### 2.Allele_specific_expression: quantifying allele-specific expression with GATK ASEReadCounter.
    Data can be downloaded from GSE252272.
    
    Dependent softwares:  BEDTools/2.29.0; star/2.6.1d; picard/2.22.3; samtools/1.9; R/4.0.3

    Running steps:
        1. Human genome can be masked by the bed file of SNPs by running "sh 1_maskGenome.sh".
        2. The RNAseq reads configured as the format of "Sample.lst" will be aligned to masked genome with STAR, by running "sh 2_align.sh".
        3. The allele specific expression for each SNP will be counted with GATK ASEReadCounter by running "sh 3_allele_sepcific_expression.sh".
        4. The allele counts will be summarized for each SNP in each sample by running "sh 4_summary.sh" with the R scripts in this folder.
  
### 3.Target_DNAseq: Aligning and genotyping with targeted DNAseq data around the LILRB3-4SNPs region.
     Data available upon request.
     
     Dependent softwares: bwa/0.7.15; samtools/1.9; bcftools/1.9 

     Running steps:
         1. The targeted sequencing reads will be aligned to LILRB3-4SNPs gene region (ref.fa in the folder) and SNPs will be called by bcftools mileup by running "sh 1_align.sh", with samples configured as "Sam" file.
         2. Each SNP in each sample will be summarized as a table by running "sh 2_summary.sh", with the R scripts in the folder. 

### 4.SingleCell: Single-cell analysis on the PBMC data.
    Data can be downloaded from GSE252273.
    
    Dependent softwares:  R/4.0.3; Seurat/4.1.1; 

    Running steps:
        1. Ambient RNA will be imputed and cleaned by SoupX for each sample configured in "List" by running "sh 1_ambient.sh" with the R scripts in the folder.
        2. QC, doublet removal, and clustering and annotation by Seurat package can be run by "sh 2_seurat.sh", with the R scripts in the folder.
        3. DEG and enrichment can be called by running "sh 3_deg.sh" with the R scripts in the folder.

### 5. Biobank: scripts for processing the biobank data.

    Data is available upon request for BioMe or All-of-Us commitee.

    Dependent software: R/4.0.3 (Metafor, Survival, )

### 5.Figures: Codes to generate figures for publication, the corresonding data is availabe upon request to the corresponding author.

[1] Sun, Z., Yi, Z., Wei, C., Wang, W., Cravedi, P., Tedla, F., Ward, S.C., Azeloglu, E., Schrider, D.R., Li, Y. and Ali, S., 2024. Genetic polymorphisms of Leukocyte Immunoglobulin-Like Receptor B3 (LILRB3) gene in African American kidney transplant recipients are associated with post-transplant graft failure. bioRxiv, pp.2024-02.
