#!/bin/bash

# sample aliases
read1=""
read2=""
sample=""
sample_dir=""
sample_prefix=${sample_dir}/${sample}

# databases list
ref=""
dbsnp138=""
Mills_indels=""
thousand_genomes_indels=""
hapmap_snps=""
omni_snps=""
thousand_genomes_snps=""

# tools
picard=""
gatk=""

# pipeline
mkdir -p ${sample_dir}

bwa mem -M -t 4 ${ref} ${read1} ${read2} > ${sample_prefix}.mem.sam

samtools view -bT ${ref} ${sample_prefix}.mem.sam > ${sample_prefix}.view.bam

samtools sort -l 9 -O bam -T ${sample_prefix}.sorted.tmp ${sample_prefix}.view.bam > ${sample_prefix}.sorted.bam

${picard} AddOrReplaceReadGroups \
  INPUT=${sample_prefix}.sorted.bam \
  OUTPUT=${sample_prefix}.ARRG.bam \
  SORT_ORDER=coordinate \
  RGID=${sample} \
  RGLB=${sample} \
  RGPL=ILLUMINA \
  RGPU=SureSelectV4 \
  RGSM=${sample} \
  RGCN=NLA \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=LENIENT \
  MAX_RECORDS_IN_RAM=1000000

${picard} MarkDuplicates \
  INPUT=${sample_prefix}.ARRG.bam \
  OUTPUT=${sample_prefix}.MD.bam \
  METRICS_FILE=${sample_prefix}.MD_${picard}_metrics.txt \
  ASSUME_SORTED=true \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=LENIENT \
  MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

${gatk} -T BaseRecalibrator \
  -R ${ref} \
  -I ${sample_prefix}.MD.bam \
  -knownSites ${dbsnp138} \
  -knownSites ${Mills_indels} \
  -knownSites ${thousand_genomes_indels} \
  -o ${sample_prefix}.BR_table.txt

${gatk} -T BaseRecalibrator \
  -R ${ref} \
  -I ${sample_prefix}.MD.bam \
  -knownSites ${dbsnp138} \
  -knownSites ${Mills_indels} \
  -knownSites ${thousand_genomes_indels} \
  -BQSR ${sample_prefix}.BR_table.txt \
  -o ${sample_prefix}.BQSR_BR_table.txt

# warning this step may fail due to absence of R modules
# you can safely skip this step
${gatk} -T AnalyzeCovariates \
  -R ${ref} \
  -before ${sample_prefix}.BR_table.txt \
  -after ${sample_prefix}.BQSR_BR_table.txt \
  -plots ${sample_prefix}.AC_plot.pdf

${gatk} -T PrintReads \
  -R ${ref} \
  -I ${sample_prefix}.MD.bam \
  -BQSR ${sample_prefix}.BQSR_BR_table.txt \
  -o ${sample_prefix}.BQSR_BR.bam

${gatk} -T HaplotypeCaller \
  -R ${ref} \
  -I ${sample_prefix}.BQSR_BR.bam \
  -D ${dbsnp138} \
  --genotyping_mode DISCOVERY \
  -stand_call_conf 30 \
  -o ${sample_prefix}.HC.vcf \
  -nct 4

${gatk} -T VariantRecalibrator \
  -R ${ref} \
  -input ${sample_prefix}.HC.vcf \
  -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap_snps} \
  -resource:omni,known=false,training=true,truth=true,prior=12.0 ${omni_snps} \
  -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${thousand_genomes_snps} \
  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp138} \
  -an DP \
  -an QD \
  -an FS \
  -an MQRankSum \
  -an ReadPosRankSum \
  -mode SNP \
  -tranche 100.0 \
  -tranche 99.9 \
  -tranche 99.0 \
  -tranche 90.0 \
  --maxGaussians 4 \
  -recalFile ${sample_prefix}.HC.VR_SNP.recal \
  -tranchesFile ${sample_prefix}.HC.VR_SNP.tranches \
  -rscriptFile ${sample_prefix}.HC.VR_SNP_plots.R

${gatk} -T ApplyRecalibration \
  -R ${ref} \
  -input ${sample_prefix}.HC.vcf \
  -mode SNP \
  --ts_filter_level 99.0 \
  -recalFile ${sample_prefix}.HC.VR_SNP.recal \
  -tranchesFile ${sample_prefix}.HC.VR_SNP.tranches \
  -o ${sample_prefix}.HC.VQSR_AR_SNP.vcf

${gatk} -T VariantRecalibrator \
  -R ${ref} \
  -input ${sample_prefix}.HC.VQSR_AR_SNP.vcf \
  -resource:mills,known=true,training=true,truth=true,prior=12.0 ${Mills_indels} \
  -an DP \
  -an FS \
  -an MQRankSum \
  -an ReadPosRankSum \
  -mode INDEL \
  -tranche 100.0 \
  -tranche 99.9 \
  -tranche 99.0 \
  -tranche 90.0 \
  --maxGaussians 4 \
  -recalFile ${sample_prefix}.HC.VR_INDEL.recal \
  -tranchesFile ${sample_prefix}.HC.VR_INDEL.tranches \
  -rscriptFile ${sample_prefix}.HC.VR_INDEL_plots.R

${gatk} -T ApplyRecalibration \
  -R ${ref} \
  -input ${sample_prefix}.HC.VQSR_AR_SNP.vcf \
  -mode INDEL \
  --ts_filter_level 99.0 \
  -recalFile ${sample_prefix}.HC.VR_INDEL.recal \
  -tranchesFile ${sample_prefix}.HC.VR_INDEL.tranches \
  -o ${sample_prefix}.HC.VQSR_AR_SNP.VQSR_AR_INDEL.vcf

${gatk} -T SelectVariants \
  -R ${ref} \
  -V ${sample_prefix}.HC.VQSR_AR_SNP.VQSR_AR_INDEL.vcf \
  -selectType SNP \
  -o ${sample_prefix}.HC.VQSR.raw_snp.vcf

${gatk} -T VariantFiltration \
  -R ${ref} \
  -V ${sample_prefix}.HC.VQSR.raw_snp.vcf \
  --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
  --filterName "SNP_FAIL" \
  -o ${sample_prefix}.HC.VQSR.fil_snp.vcf

${gatk} -T SelectVariants \
  -R ${ref} \
  -V ${sample_prefix}.HC.VQSR_AR_SNP.VQSR_AR_INDEL.vcf \
  -selectType INDEL \
  -o ${sample_prefix}.HC.VQSR.raw_indel.vcf

${gatk} -T VariantFiltration \
  -R ${ref} \
  -V ${sample_prefix}.HC.VQSR.raw_indel.vcf \
  --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
  --filterName "INDEL_FAIL" \
  -o ${sample_prefix}.HC.VQSR.fil_indel.vcf

vcf-concat ${sample_prefix}.HC.VQSR.fil_snp.vcf ${sample_prefix}.HC.VQSR.fil_indel.vcf > ${sample_prefix}.FINAL.concat.vcf

cat ${sample_prefix}.FINAL.concat.vcf | vcf-sort > ${sample_prefix}.FINAL.sorted.vcf
