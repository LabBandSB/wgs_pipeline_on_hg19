# create environment variables for the project
conda create -y -n wgs 

# activate the environment
conda activate wgs 

conda install -c bioconda git python=3.9 # supplementary
conda install -y -c anaconda perl # for perl-vcftools-vcf
# broken # conda install -y -c dan_blanchard glibc # for bwa # seems it causes segmentation fault error
yum install glibc glibc-devel # use instead 

conda install -c bioconda https://anaconda.org/bioconda/bwa/0.7.12/download/linux-64/bwa-0.7.12-1.tar.bz2
conda install -c bioconda https://anaconda.org/bioconda/samtools/1.2/download/linux-64/samtools-1.2-2.tar.bz2

conda install -c bioconda fastqc picard gatk perl-vcftools-vcf tabix

# index the reference genome
ref=""
bwa index ${ref}
samtools faidx ${ref}
