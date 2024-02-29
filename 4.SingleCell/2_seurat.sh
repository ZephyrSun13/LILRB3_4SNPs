
set -e

Wk_dir=$(pwd)

module load R/4.0.3

Rscript ${Wk_dir}/seurat4_integrate_Har.r List 10x ${Wk_dir}/Markers.pbmc 0.8 30 30 3 200 7000 30 All 2000 &> seurat.log

awk 'BEGIN{FS="\t"; OFS="\t"} {print $0, $2}' ${Wk_dir}/cluster.lst > ${Wk_dir}/cluster.lst.ori

Rscript ${Wk_dir}/doubletFinder.r ${Wk_dir}/cluster.lst.ori_Obj_post.rds

awk '(NR==FNR){if($14=="Doublet")(tt[$1] = 1); next} {if(tt[$1]!=1){print $0}}' ${Wk_dir}/DBTab.xls ${Wk_dir}/cluster.lst.ori > tt && mv tt ${Wk_dir}/cluster.lst.ori

Rscript ${Wk_dir}/seurat4_integrate_post.r ${Wk_dir}/Obj.rds No ${Wk_dir}/cluster.lst.ori ${Wk_dir}/Sample.info Null

Rscript ${Wk_dir}/seurat4_checkMarkers.r ${Wk_dir}/cluster.lst.ori_Obj_post.rds ${Wk_dir}/Markers.pbmc cluster.lst.ori ${Wk_dir}/Markers.pbmc RNA Null


cut -f 2 ${Wk_dir}/cluster.lst | sort | uniq | sort -n > ${Wk_dir}/anno.lst

awk 'BEGIN{FS = "\t"; OFS = "\t";} (NR==FNR){if($2!="RM"){tt[$1] = $2} next;} {if(tt[$2] != ""){print $0, tt[$2]}}' ${Wk_dir}/anno.lst ${Wk_dir}/cluster.lst > ${Wk_dir}/cluster.lst.ann

Rscript ${Wk_dir}/seurat4_integrate_post.r ${Wk_dir}/Obj.rds No ${Wk_dir}/cluster.lst.ann ${Wk_dir}/Sample.info "NK_CD56Bright;NK_CD56Dim;CD8T;CD4CD8;Proliferating;CD4T;CD14Mono;CD16Mono;DC;B;Platelet;CD14CD16Mono"

