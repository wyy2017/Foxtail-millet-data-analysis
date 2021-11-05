## user defined software path ## 
#-- you must install these softwares, and test by yourself to make sure they work well
perl=/usr/bin/perl
Rscript=/dellfsqd2/ST_OCEAN/USER/sunshuai/software/miniconda3/bin/Rscript
plink=/dellfsqd2/ST_OCEAN/USER/sunshuai/software/variation/plink_1.9/plink
gemma=/hwfssz1/ST_OCEAN/USER/sunshuai/software/gwas/gemma_v0.98/gemma


#-- data used in this analysis
snp=./millet.snp.gwas   # prefix of genotype data (PLINK binary format which has three files including: .bed, .bim and .fam)
phe=./growth_yiled.phe
outdir=./host_phe

for i in MSPD MSW PGW MSPL; do
    spt_dir=$outdir/$i
    mkdir $spt_dir

    # kinship to correct the relatedness of indiviuals
    $gemma --bfile $snp -p $phe -miss 0.1 -maf 0.01 -gk -outdir $spt_dir -o snp.kinship
    
    # first ten PC to represent the population structure/strafation.
    $plink --bfile $snp --pca 10 --out $spt_dir/PC10
    cut -d ' ' -f 3- $spt_dir/PC10.eigenvec | sed 's/^/1 /' > $spt_dir/snp_PC10

    # Genome-wide association study (GWAS) by GEMMA
    $gemma --bfile $snp -p $phe -miss 0.1 -n $i -k $spt_dir/snp.kinship.cXX.txt -c $spt_dir/snp_PC10 -lmm -outdir $spt_dir -o $i

    # format type of p-value
    echo "SNP CHROM BP P" > $spt_dir/$i.assoc.for_plot.txt
    tail -n +2 $spt_dir/$i.assoc.txt | perl -ne '@tmp=split; print "si$tmp[0]:$tmp[2] $tmp[0] $tmp[2] $tmp[11]\n"' >> $spt_dir/$i.assoc.for_plot.txt
    
    # plot GWAS result
    cat > $spt_dir/$i.CMplot.r <<EOF
library(CMplot)
df <- read.table("$spt_dir/$i.assoc.for_plot.txt", header=True)
CMplot(df, plot.type = "m", LOG10=TRUE, cex=0.4, threshold=2.2E-5,, threshold.lty=2, threshold.col="black",  amplify=TRUE, signal.col="red", signal.cex=1, width=14, bin.size=1e6, chr.den.col=c("darkgreen", "yellow", "red"), verbose=TRUE, height=5, file="pdf")
EOF
    $Rscript $spt_dir/$i.CMplot.r
done
