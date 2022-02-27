

mkdir -p ${alignment_dir}


# lock sample to prevent runnig from 2 servers 
[ ! -f ${alignment_dir}/token.${sample}.__first_run_lock__ ] && \
touch ${alignment_dir}/token.${sample}.__first_run_lock__ \
|| exit 1


# bwa alignment
token="${alignment_dir}/token.${sample}.fastq_2_sam_bwa_mem"
output_file="${alignment_dir}/${sample}.bwa_mem.sam"
[ ! -f ${token} ] && \
[ -f ${read1} ] && \
[ -f ${read2} ] && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${bwa} mem -M -t 4 ${ref} ${read1} ${read2} > ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"



# samtools convert sam to bam
token="${alignment_dir}/token.${sample}.sam_2_bam_samtools_view"
input_file="${alignment_dir}/${sample}.bwa_mem.sam"
output_file="${alignment_dir}/${sample}.samtools_view.bam"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${samtools} view -bT ${ref} ${input_file} > ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"



# samtools sort bam to bam
token="${alignment_dir}/token.${sample}.bam_2_bam_samtools_sort"
input_file="${alignment_dir}/${sample}.samtools_view.bam"
output_file="${alignment_dir}/${sample}.samtools_sort.bam"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${samtools} sort -l 9 -O bam -T ${alignment_dir}/${sample}.sorted.tmp ${input_file} > ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# picard_ARRG
token="${alignment_dir}/token.${sample}.bam_2_bam_picard_ARRG"
input_file="${alignment_dir}/${sample}.samtools_sort.bam"
output_file="${alignment_dir}/${sample}.picard_ARRG.bam"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${picard} AddOrReplaceReadGroups \
  INPUT=${input_file} \
  OUTPUT=${output_file} \
  SORT_ORDER=coordinate \
  RGID=${sample} \
  RGLB=${sample} \
  RGPL=ILLUMINA \
  RGPU=SureSelectV4 \
  RGSM=${sample} \
  RGCN=NLA \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=LENIENT \
  MAX_RECORDS_IN_RAM=1000000 && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# #######
token="${alignment_dir}/token.${sample}.bam_2_bam_picard_MD"
input_file="${alignment_dir}/${sample}.picard_ARRG.bam"
output_file="${alignment_dir}/${sample}.picard_MD.bam"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${picard} MarkDuplicates \
  INPUT=${input_file} \
  OUTPUT=${output_file} \
  METRICS_FILE=${sample_prefix}.poicard_MD_metrics.txt \
  ASSUME_SORTED=true \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=LENIENT \
  MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"

# #######
token="${alignment_dir}/token.${sample}.bam_2_txt_gatk_BR"
input_file="${alignment_dir}/${sample}.picard_MD.bam"
output_file="${alignment_dir}/${sample}.BR_table.txt"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${gatk} -T BaseRecalibrator \
  -R ${ref} \
  -I ${input_file} \
  -knownSites ${dbsnp138} \
  -knownSites ${Mills_indels} \
  -knownSites ${thousand_genomes_indels} \
  -o ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"

# #######
token="${alignment_dir}/token.${sample}.txt_2_txt_gatk_BR"
input_file="${alignment_dir}/${sample}.picard_MD.bam"
input_file_2="${alignment_dir}/${sample}.BR_table.txt"
output_file="${alignment_dir}/${sample}.BQSR_BR_table.txt"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${gatk} -T BaseRecalibrator \
  -R ${ref} \
  -I ${input_file} \
  -knownSites ${dbsnp138} \
  -knownSites ${Mills_indels} \
  -knownSites ${thousand_genomes_indels} \
  -BQSR ${input_file_2} \
  -o {output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# #######
token="${alignment_dir}/token.${sample}.txt_2_pdf_gatk_AC"
input_file="${alignment_dir}/${sample}.BR_table.txt"
input_file_2="${alignment_dir}/${sample}.BQSR_BR_table.txt"
output_file="${alignment_dir}/${sample}.gatk_AC_plot.pdf"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${gatk} -T AnalyzeCovariates \
  -R ${ref} \
  -before ${input_file} \
  -after ${input_file_2} \
  -plots ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# #######
token="${alignment_dir}/token.${sample}.bam_2_bam_gatk_PR"
input_file="${alignment_dir}/${sample}.picard_MD.bam"
input_file_2="${alignment_dir}/${sample}.BQSR_BR_table.txt"
output_file="${alignment_dir}/${sample}.gatk_PR_BQSR_BR_table.bam"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${gatk} -T PrintReads \
  -R ${ref} \
  -I ${input_file} \
  -BQSR ${input_file_2} \
  -o ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# #######
token="${alignment_dir}/token.${sample}.bam_2_vcf_gatk_HC"
input_file="${alignment_dir}/${sample}.gatk_PR_BQSR_BR_table.bam"
output_file="${alignment_dir}/${sample}.gatk_HC.vcf"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${gatk} -T HaplotypeCaller \
  -R ${ref} \
  -I ${input_file} \
  -D ${dbsnp138} \
  --genotyping_mode DISCOVERY \
  -stand_call_conf 30 \
  -o ${output_file} \
  -nct 4 && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# #######
token="${alignment_dir}/token.${sample}.vcf_2_recal_tranches_gatk_VR_SNP"
input_file="${alignment_dir}/${sample}.gatk_HC.vcf"
output_file="${alignment_dir}/${sample}.gatk_VR_SNP.recal"
output_file_2="${alignment_dir}/${sample}.gatk_VR_SNP.tranches"
output_file_3="${alignment_dir}/${sample}.VR_SNP.R"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${gatk} -T VariantRecalibrator \
  -R ${ref} \
  -input ${input_file} \
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
  -recalFile ${output_file} \
  -tranchesFile ${output_file_2} \
  -rscriptFile ${output_file_3} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# #######
token="${alignment_dir}/token.${sample}.vcf_2_vcf_gatk_AR_SNP_VQSR"
input_file="${alignment_dir}/${sample}.gatk_HC.vcf"
input_file_2="${alignment_dir}/${sample}.gatk_VR_SNP.recal"
input_file_3="${alignment_dir}/${sample}.gatk_VR_SNP.tranches"
output_file="${alignment_dir}/${sample}.gatk_AR_SNP_VQSR.vcf"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
[ -f ${input_file_2} ] && \
[ -f ${input_file_3} ] && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${gatk} -T ApplyRecalibration \
  -R ${ref} \
  -input ${sample_prefix}.HC.vcf \
  -mode SNP \
  --ts_filter_level 99.0 \
  -recalFile ${input_file_2} \
  -tranchesFile ${input_file_3} \
  -o ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# #######
token="${alignment_dir}/token.${sample}.vcf_2_recal_tranches_gatk_VR_INDEL"
input_file="${alignment_dir}/${sample}.gatk_AR_SNP_VQSR.vcf"
output_file="${alignment_dir}/${sample}.gatk_VR_INDEL.recal"
output_file_2="${alignment_dir}/${sample}.gatk_VR_INDEL.tranches"
output_file_3="${alignment_dir}/${sample}.VR_INDEL.R"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${gatk} -T VariantRecalibrator \
  -R ${ref} \
  -input ${input_file} \
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
  -recalFile ${output_file} \
  -tranchesFile ${output_file_2} \
  -rscriptFile ${output_file_3} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# #######
token="${alignment_dir}/token.${sample}.vcf_2_vcf_gatk_AR_INDEL_VQSR"
input_file="${alignment_dir}/${sample}.gatk_AR_SNP_VQSR.vcf"
input_file_2="${alignment_dir}/${sample}.gatk_VR_INDEL.recal"
input_file_3="${alignment_dir}/${sample}.gatk_VR_INDEL.tranches"
output_file="${alignment_dir}/${sample}.gatk_AR_SNP_VQSR.gatk_AR_INDEL_VQSR.vcf"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
[ -f ${input_file_2} ] && \
[ -f ${input_file_3} ] && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${gatk} -T ApplyRecalibration \
  -R ${ref} \
  -input ${input_file} \
  -mode INDEL \
  --ts_filter_level 99.0 \
  -recalFile ${input_file_2} \
  -tranchesFile ${input_file_3} \
  -o ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# #######
token="${alignment_dir}/token.${sample}.vcf_2_vcf_gatk_SV_SNP"
input_file="${alignment_dir}/${sample}.gatk_AR_SNP_VQSR.gatk_AR_INDEL_VQSR.vcf"
output_file="${alignment_dir}/${sample}.gatk_SV_SNP.vcf"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${gatk} -T SelectVariants \
  -R ${ref} \
  -V ${input_file} \
  -selectType SNP \
  -o ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# #######
token="${alignment_dir}/token.${sample}.vcf_2_vcf_gatk_VF_SNP"
input_file="${alignment_dir}/${sample}.gatk_SV_SNP.vcf"
output_file="${alignment_dir}/${sample}.gatk_VF_SNP.vcf"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${gatk} -T VariantFiltration \
  -R ${ref} \
  -V ${input_file} \
  --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
  --filterName "SNP_FAIL" \
  -o ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# #######
token="${alignment_dir}/token.${sample}.vcf_2_vcf_gatk_SV_INDEL"
input_file="${alignment_dir}/${sample}.gatk_AR_SNP_VQSR.gatk_AR_INDEL_VQSR.vcf"
output_file="${alignment_dir}/${sample}.gatk_SV_INDEL.vcf"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${gatk} -T SelectVariants \
  -R ${ref} \
  -V ${input_file} \
  -selectType INDEL \
  -o ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# #######
token="${alignment_dir}/token.${sample}.vcf_2_vcf_gatk_VF_INDEL"
input_file="${alignment_dir}/${sample}.gatk_SV_INDEL.vcf"
output_file="${alignment_dir}/${sample}.gatk_VF_INDEL.vcf"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
${gatk} -T VariantFiltration \
  -R ${ref} \
  -V ${input_file} \
  --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
  --filterName "INDEL_FAIL" \
  -o ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# #######
token="${alignment_dir}/token.${sample}.vcf_vcf_2_vcf_vcf_concat"
input_file="${alignment_dir}/${sample}.gatk_VF_SNP.vcf"
input_file_2="${alignment_dir}/${sample}.gatk_VF_INDEL.vcf"
output_file="${alignment_dir}/${sample}.vcf_concat.vcf"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
[ -f ${input_file_2} ] && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
vcf-concat ${input_file} ${input_file_2} > ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# #######
token="${alignment_dir}/token.${sample}.vcf_2_vcf_vcf_sort"
input_file="${alignment_dir}/${sample}.vcf_concat.vcf"
output_file="${alignment_dir}/${sample}.vcf_sort.vcf"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
cat ${input_file} | vcf-sort > ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"


# #######
token="${alignment_dir}/token.${sample}.vcf_2_FINAL"
input_file="${alignment_dir}/${sample}.vcf_sort.vcf"
output_file="${alignment_dir}/${sample}.FINAL.vcf"
[ ! -f ${token} ] && \
[ -f ${input_file} ] && \
dt1=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${token} && \
cp ${input_file} ${output_file} && \
du ${output_file} > ${output_file}.${dt1}.du && \
md5sum ${output_file} > ${output_file}.${dt1}.md5 && \
dt2=`date +%y%m%d_%H%M%S` && \
echo ${dt1} ${dt2} > ${token} \
|| echo "TOKEN SKIPPED ${token}"

# THE END 
echo "DONE!"