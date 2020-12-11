#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 8 -R "rusage[mem=8]"
#BSUB -W 48:00

### Refer to  ~/lingqi_workspace/Projects/Mitesh/WGS/intel-gatk4-somatic-with-preprocessing/mutect2_nodocker.wdl

module load java/1.8.0_31

ref_fasta="/data/ldiaz/luol2/gatk_resource_bundle/hg38/Homo_sapiens_assembly38.fasta"
data_path=DATA_PATH
project=PROJECT
subject=SUBJECT
sample=TEST_SAMPLE


tumor_bam=${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicates_marked.recalibrated.bam

command_mem="8000"

# index the bam
if [ ! -e ${tumor_bam/.bam/.bai} ];then
module load samtools/1.7
samtools index ${tumor_bam}
fi

# 1. Generate OXOG metrics:

/data/ldiaz/luol2/gatk-4.1.2.0/gatk --java-options "-Xmx${command_mem}m" \
CollectSequencingArtifactMetrics \
-I ${tumor_bam} \
-O ${data_path}/${project}/${subject}/${sample}/${sample} \
--FILE_EXTENSION .txt \
-R ${ref_fasta}

#2. Generate pileup summaries on tumor sample:
/data/ldiaz/luol2/gatk-4.1.2.0/gatk --java-options "-Xmx${command_mem}m" \
GetPileupSummaries \
-I ${tumor_bam} \
-O ${data_path}/${project}/${subject}/${sample}/${sample}.targeted_sequencing.table \
-V /data/ldiaz/luol2/gatk_resource_bundle/mutect2/af-only-gnomad-common-biallelic.grch38.main.vcf.gz \
-L /home/luol2/lingqi_workspace/SureSelect_All_Exon_BED/hg38/S31285117_Regions.bed \
-R ${ref_fasta}

#-L intervals.bed \ ## Only chr1-22 + XYM

## 3. Calculate contamination on tumor sample

/data/ldiaz/luol2/gatk-4.1.2.0/gatk --java-options "-Xmx${command_mem}m" \
CalculateContamination \
-I ${data_path}/${project}/${subject}/${sample}/${sample}.targeted_sequencing.table \
-O ${data_path}/${project}/${subject}/${sample}/${sample}.targeted_sequencing.contamination.table





