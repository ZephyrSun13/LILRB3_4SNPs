
Wk_dir=$(pwd)

ls -d ${Wk_dir}/cellRanger.sh.138377.qsub/*/outs | awk '{split($0, tt, "/"); printf("%s\t%s\n", $1, tt[9])}' > Files

ml R/4.0.3

Rscript ${Wk_dir}/ambient_soupX.r Files

