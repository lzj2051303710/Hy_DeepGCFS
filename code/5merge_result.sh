ls ../Quantification/*.count >genes.quant_files.txt

perl script/abundance_estimates_to_matrix.pl --est_method featureCounts --quant_files genes.quant_files.txt --out_prefix genes
