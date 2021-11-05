## user defined software path ##
#-- you must install these softwares, and test by yourself to make sure they work well
perl=/usr/bin/perl
Rscript=/dellfsqd2/ST_OCEAN/USER/sunshuai/software/miniconda3/bin/Rscript
plink=/dellfsqd2/ST_OCEAN/USER/sunshuai/software/variation/plink_1.9/plink
gemma=/hwfssz1/ST_OCEAN/USER/sunshuai/software/gwas/gemma_v0.98/gemma

#-- data used in this analysis
snp=./millet.snp.gwas   # prefix of genotype data (PLINK binary format which has three files including: .bed, .bim and .fam)
phe=./OTU.rank.phe
outdir=./host_OTU

# GWAS analysis for all available OTU 
##-- total OTU number
total=`wc -l ${phe}.id | cut -d ' ' -f 1`
##-- run one by one
for i in $(seq 1 $total); do
    name=`head -n $i ${phe}.id | awk '{print $1}' | tail -n 1`
    spt_dir=$outdir/$name
    mkdir $spt_dir
    
    # kinship to correct the relatedness of indiviuals
    $gemma --bfile $snp -p $phe -miss 0.1 -maf 0.01 -gk -outdir $spt_dir -o snp.kinship

    # first ten PC to represent the population structure/strafation.
    $plink --bfile $snp --pca 10 --out $spt_dir/PC10
    cut -d ' ' -f 3- $spt_dir/PC10.eigenvec | sed 's/^/1 /' > $spt_dir/snp_PC10

    # Genome-wide association study (GWAS) by GEMMA
    $gemma --bfile $snp -p $phe -miss 0.1 -n $i -k $spt_dir/snp.kinship.cXX.txt -c $spt_dir/snp_PC10 -lmm -outdir $spt_dir -o $name
    
    # plot GWAS result
    cat > $spt_dir/$name.CMplot.r <<EOF
library(CMplot)
df <- read.table("$spt_dir/$name.assoc.for_plot.txt", header=True)
CMplot(df, plot.type = "m", LOG10=TRUE, cex=0.4, threshold=2.2E-5,, threshold.lty=2, threshold.col="black",  amplify=TRUE, signal.col="red", signal.cex=1,width=14, bin.size=1e6, chr.den.col=c("darkgreen", "yellow", "red"), verbose=TRUE, height=5, file="pdf")
EOF
    $Rscript $spt_dir/$name.CMplot.r
done

# combine interested GWAS result of some OTUs
##-- select only one p-value of the most significance for each SNP
for name in `cat ./interested_OTU.id.list`; do
    echo "$outdir/$name.assoc.txt" >> $outdir/combined.assoc.list
done
perl ./min_p.from.multi.pl $outdir/combined.assoc.list $outdir/combined.p_fmt.lst

# create R script to plot Manhantun figure
cat > $outdir/combined.assoc.plot.r <<EOF
library(CMplot)
df <- read.table("$outdir/combined.p_fmt.lst", header=F)
CMplot(df, type="p", plot.type = "m", LOG10=TRUE, threshold.lty=2, threshold.col="black", threshold=2.2E-5, amplify=FALSE, width=14, height=5, cex=0.4, file="pdf")
EOF

# run this R script
Rscript $outdir/combined.assoc.plot.r
