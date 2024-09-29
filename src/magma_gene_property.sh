#!/bin/bash

# # Run matlab make_gene_covar.m to prepare data for magma gene-property analysis
matlab -nodisplay -nosplash -nodesktop -r "run('./make_gene_covar.m');exit;"


# SCZ PGC wav-3
rawfile='../data/SCZ_GWAS/magma.genes.raw'
genecovarfile='../data/genecovar.txt'
outfileprefix='../results/gsa_scz_log2fc'

./magma/magma \
--gene-results $rawfile \
--gene-covar $genecovarfile --model direction=pos  \
--out $outfileprefix


# SCZ PGC eas
rawfile='../data/SCZ_EAS_GWAS/magma.genes.raw'
outfileprefix='../results/gsa_scz_eas_log2fc'

./magma/magma \
--gene-results $rawfile \
--gene-covar $genecovarfile --model direction=pos  \
--out $outfileprefix


# PGC BD
rawfile='../data/BD_GWAS/magma.genes.raw'
outfileprefix='../results/gsa_bd_log2fc'

./magma/magma \
--gene-results $rawfile \
--gene-covar $genecovarfile --model direction=pos  \
--out $outfileprefix


# PGC ADHD
rawfile='../data/ADHD_GWAS/magma.genes.raw'
outfileprefix='../results/gsa_adhd_log2fc'

./magma/magma \
--gene-results $rawfile \
--gene-covar $genecovarfile --model direction=pos  \
--out $outfileprefix


# PGC ASD
rawfile='../data/ASD_GWAS/magma.genes.raw'
outfileprefix='../results/gsa_asd_log2fc'

./magma/magma \
--gene-results $rawfile \
--gene-covar $genecovarfile --model direction=pos  \
--out $outfileprefix


# PGC MDD
rawfile='../data/MDD_GWAS/magma.genes.raw'
outfileprefix='../results/gsa_mdd_log2fc'

./magma/magma \
--gene-results $rawfile \
--gene-covar $genecovarfile --model direction=pos  \
--out $outfileprefix


# PGC Insomnia
rawfile='../data/INS_GWAS/magma.genes.raw'
outfileprefix='../results/gsa_ins_log2fc'

./magma/magma \
--gene-results $rawfile \
--gene-covar $genecovarfile --model direction=pos  \
--out $outfileprefix


