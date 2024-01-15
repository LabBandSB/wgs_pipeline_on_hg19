import argparse
import json
import os
from collections import defaultdict

d = {
    "RGID": "__sample__",
    "RGLB": "__sample__",
    "RGPL": "ILLUMINA",
    "RGPU": "SureSelectV4",
    "RGSM": "__sample__",
    "RGCN": "NLA",
    #
    "fastqc": "fastqc",
    "bwa": "bwa",
    "samtools": "samtools",
    "bcftools": "bcftools",
    "java": "java",
    "picard": "picard",
    "gatk": "gatk3",
    "vcf_concat": "vcf-concat",
    "vcf_sort": "vcf-sort",
    "vcf_merge": "vcf-merge",
    "bgzip": "bgzip",
    "tabix": "tabix",
    #
    "dbsnp": "/home/PublicData/broadinstitute/2.8/hg19/dbsnp_138.hg19.vcf",
    "gold_indel": "/home/PublicData/broadinstitute/2.8/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf",
    "hapmap_snp": "/home/PublicData/broadinstitute/2.8/hg19/hapmap_3.3.hg19.sites.vcf",
    "oneKG_indel": "/home/PublicData/broadinstitute/2.8/hg19/1000G_phase1.indels.hg19.sites.vcf",
    "oneKG_snp": "/home/PublicData/broadinstitute/2.8/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf",
    "onmi_snp": "/home/PublicData/broadinstitute/2.8/hg19/1000G_omni2.5.hg19.sites.vcf",
    "ref": "/home/PublicData/broadinstitute/2.8/hg19/ucsc.hg19.ref/ucsc.hg19.fasta",
    "ref_mtDNA": "/home/PublicData/h.sapiens_mtDNA/HS_mtDNA.fa",
    "ref_ucsc_hg19": "/home/PublicData/broadinstitute/2.8/hg19/ucsc.hg19.ref/ucsc.hg19.fasta",
    "target_region": "/home/PublicData/Agilent_v4_71m_reduced.bed",
}


def parse_arguments_to_settings():
    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--settings_json", default=None, required=True)
    parser.add_argument("-s", "--samples_txt", default=None, required=True)

    args = parser.parse_args()
    if args.settings_json:
        settings = json.load(open(args.settings_json))
    else:
        settings = dict()
    if args.samples_txt:
        sample = os.path.splitext(os.path.basename(args.samples_txt))[0]
        settings["sample"] = f"MSAP_{sample}"

        samples = [line.strip() for line in open(args.samples_txt)]
        settings["samples"] = samples

    else:
        raise Exception("No --samples_txt file provided")

    settings["read1"] = ""
    settings["read2"] = ""

    return settings


def prepare_msap_sh(settings):
    d.update(settings)
    sample = d["sample"]
    project_dir = d["project_dir"]
    script_dir = d["script_dir"]

    d["alignment_dir"] = os.path.join(project_dir, sample)
    bam_files = [
        os.path.join(project_dir, s, f"{s}.gatk_PR_BQSR_BR_table.bam")
        for s in d["samples"]
    ]
    d["input_bams"] = " ".join([f"-I {bam}" for bam in bam_files])

    script_name = f"{sample}.sh"
    script_file = os.path.join(script_dir, script_name)
    with open(script_file, "w") as f:
        for line in open("bash_sample_template.sh"):
            new_line = line.format(**d)
            f.write(new_line)
        for line in open("bash_input_bams_template.sh"):
            new_line = line.format(**d)
            f.write(new_line)
        for line in open("bash_msap_template.sh"):
            f.write(line)


def main():
    settings = parse_arguments_to_settings()

    script_dir = settings["script_dir"]
    print(f"script_dir: {script_dir}")
    os.makedirs(script_dir, exist_ok=True)

    prepare_msap_sh(settings)


if __name__ == "__main__":
    main()
