
Wk_dir=$(pwd)

ls ${Wk_dir}/Separate/*/*_aln.dedup.bam.RG.bam.ase.output.table > Tabs

ls ${Wk_dir}/Separate/*/*_aln.dedup.bam.RG.bam.ase.output.table >> Tabs

ml R/4.0.3

Rscript ${Wk_dir}/summary.r

Rscript ${Wk_dir}/exprRefAlt.r

Rscript ${Wk_dir}/rate.r

