#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 1 -R "rusage[mem=8]"
#BSUB -W 48:00

module load java/1.8.0_31
ref_fasta="/data/ldiaz/luol2/gatk_resource_bundle/hg38/Homo_sapiens_assembly38.fasta"
ref_dbsnps="/data/ldiaz/luol2/gatk_resource_bundle/hg38/dbsnp_146.hg38.vcf.gz"

project=PROJECT
subject=SUBJECT
sample=TEST_SAMPLE

data_path=DATA_PATH

#unfiltered_vcf=${data_path}/${sample}/${sample}.aligned.duplicates_marked.recalibrated.unfiltered.vcf
unfiltered_vcf=${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicates_marked.recalibrated.targeted_sequencing.tumor_only.gatk4_mutect2.raw_somatic_mutation.vcf.gz

command_mem="8000"

/data/ldiaz/luol2/gatk-4.0.1.2/gatk  --java-options "-Xmx${command_mem}m" SelectVariants -R ${ref_fasta} -V ${unfiltered_vcf} --discordance ${ref_dbsnps} -O ${unfiltered_vcf/.vcf.gz/.rm_dbsnps.vcf.gz}

