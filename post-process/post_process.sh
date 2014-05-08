
## compute mpear label
Rscript compute_mpear_label.R $1

## compute v-measure for traces and mpear label
python compute_vmeasures.py $1 $2

## plot vmeasure
Rscript plot_vmeasure.R $1
