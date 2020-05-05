## Use vcftools to generate 012 matrix, 0=homoRef 1=hetro 2=homoAlt. -1= missing data (Matrix012loop.sh)
vcf=~/vcftools_0.1.13/bin/vcftools

for f in /STORAGE/George/GTEx/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/split_by_chr/*.vcf
do
 out="${f%.*}"
 $vcf --vcf $f --012 --out $out'_matrix012.txt'
done
