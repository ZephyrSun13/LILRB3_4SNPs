
ml R/4.0.3

awk '($10>0){printf("%s\t%s\n", $1, $10)}' meta_effect_size_p_selected2.txt > DEG.up.lst

awk '($10<0){printf("%s\t%s\n", $1, $10)}' meta_effect_size_p_selected2.txt > DEG.down.lst

Rscript enrichR.r DEG.up.lst

Rscript enrichR.r DEG.down.lst


Rscript plotEnrichNP.r DEG.up.lst_EnrichR.sig.xls UpFunc "Blue;Red"

Rscript plotEnrichNP.r DEG.down.lst_EnrichR.sig.xls DownFunc "Blue;Red"


Rscript plotEnrichNP.r DEG.up.lst_EnrichR.sig.xls UpFunc "Blue;Red"

Rscript plotEnrichNP.r DEG.down.lst_EnrichR.sig.xls DownFunc "Blue;Red"

