# 10X Genomics singlecell genomic analysis
# run one sample at a time

cellranger count --id=pyk \
--fastqs=/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/fastq/pyk \
--transcriptome=/home/bmrc/Public/phil_ubuntu/sc/gex/ref_gex/refdata-gex-mm10-2020-A \
--sample=PYK # or PYK
