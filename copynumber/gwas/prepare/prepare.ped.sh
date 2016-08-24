outfile="/home/fs01/rz253/project/uk10k/analysis/copynumber/gwas/prepare/pheno.ped"
single="/home/fs01/rz253/project/uk10k/analysis/copynumber/calculate_copy_number/single.copy.txt"
multi="/home/fs01/rz253/project/uk10k/analysis/copynumber/calculate_copy_number/multi.copy.txt"

echo "fid iid fatid matid sex copynumber" > $outfile

awk -F ";" '{if (NR!=1) print $1,$1,0,0,2,$2}' $single >> $outfile
awk -F ";" '{if (NR!=1) print $1,$1,0,0,2,$2}' $multi >> $outfile






