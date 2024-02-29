
Wk_dir=$(pwd)

ml R/4.0.3

Rscript ${Wk_dir}/seurat4_integrate_deg.r ${Wk_dir}/cluster.lst.ann_Obj_post.rds ${Wk_dir}/Sample.info "Risk;Ctrl" "NK_CD56Bright;NK_CD56Dim;CD8T;CD4CD8;Proliferating;CD4T;CD14CD16Mono;CD14Mono;CD16Mono;DC;B;Platelet"

for dd in ${Wk_dir}/*DEG.xls
do

        awk '($2<=0.001 && ($3>=0.25 || $3<=-0.25) && ($4>=0.1 && $5>=0.1) && NR>1){printf("%s\t%s\n", $1, $3)}' $dd > ${Wk_dir}/${dd}.all.lst

        awk '($2<=0.001 && $3>=0.25 && ($4>=0.1 && $5>=0.1) && NR>1){printf("%s\t%s\n", $1, $3)}' $dd > ${Wk_dir}/${dd}.pos.lst

        awk '($2<=0.001 && $3<=-0.25 && ($4>=0.1 && $5>=0.1) && NR>1){printf("%s\t%s\n", $1, $3)}' $dd > ${Wk_dir}/${dd}.neg.lst

done


ml R/4.0.3

for ll in ${Wk_dir}/*DEG*.lst
do

        Num=$(wc -l $ll | awk '{print $1}')

        if [ "$Num" -gt 20 ]
        then

                Rscript ${Wk_dir}/enrichR.r $ll

        fi

done

