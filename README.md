# TaxDistill

CUDA_VISIBLE_DEVICES=0,1,2,3 torchrun \
     --nproc_per_node=4 \
     --master_addr="127.0.0.1" \
     --master_port=29500 \
     vamb/__main__.py taxometer \
     --outdir your_dir \
     --fasta your_data/contigs_oral.fna.gz \
     --abundance_tsv your_data/abundances.tsv \
     --taxonomy your_data/disk1/yerongye/TaxDistill/Taxnometer/vamb/data/Marine/mmseq2_results_taxometer_fixed.tsv \
     --genomeocean GenomeOcean/100M  \
     --kd_alpha 0.3   --kd_temp 4 -pt 64 --cuda
