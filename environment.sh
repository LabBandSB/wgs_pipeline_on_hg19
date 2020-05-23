#
conda install -c bioconda git python fastqc bwa samtools picard gatk perl-vcftools-vcf tabix
# conda install -c bioconda https://anaconda.org/bioconda/samtools/1.2/download/linux-64/samtools-1.2-2.tar.bz2

ref=""
bwa index ${ref}
samtools faidx ${ref}