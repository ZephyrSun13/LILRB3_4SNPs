
ml R/3.5.3

Rscript limmaVoom.r Expr.xls Factor.xls unpaired X007SNPs "Yes;No" 0 0.05 Ensembl Group007 No voom 1 /sc/arion/projects/zhangw09a/PANDA/db_ZS/Ensembl/Ann_ensembl_refseq_gene.xls

awk -F "\t" '($5 <= 0.1){printf("%s\t%s\n", $8, $2)}' DEG_diff.ann.xls | sort | uniq > AllGene

awk -F "\t" '($5 <= 0.05 && $2 > 0){printf("%s\t%s\n", $8, $2)}' DEG_diff.ann.xls | sort | uniq > PosGene

awk -F "\t" '($5 <= 0.05 && $2 < 0){printf("%s\t%s\n", $8, $2)}' DEG_diff.ann.xls | sort | uniq > NegGene

sed '1d' DEG_diff.ann.xls | sort -k 4,4 -gr | awk '($4!="NA"){printf("%s\t%s\n", $8, $4)}' > DEG.lst

Rscript fgsea.r DEG.lst GSEA/c2.all.v7.4.symbols.kegg.gmt C3KEGG 1 1

Rscript fgsea.plot.r C3KEGG_fgseaRes.rds C3KEGG.plot DEG.lst GSEA/c2.all.v7.4.symbols.kegg.gmt C3KEGG 1

