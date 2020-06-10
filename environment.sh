#
# conda install -c anaconda perl # for perl-vcftools-vcf
# conda install -c dan_blanchard glibc # for bwa



conda install -c bioconda git python fastqc bwa samtools picard gatk perl-vcftools-vcf tabix
# conda install -c bioconda https://anaconda.org/bioconda/samtools/1.2/download/linux-64/samtools-1.2-2.tar.bz2
# conda install -c bioconda https://anaconda.org/bioconda/bwa/0.7.12/download/linux-64/bwa-0.7.12-1.tar.bz2
ref=""
bwa index ${ref}
samtools faidx ${ref}