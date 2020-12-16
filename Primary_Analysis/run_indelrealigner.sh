#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 4 -R "rusage[mem=8]"
#BSUB -W 48:00

### Refer to  ~/lingqi_workspace/Projects/Mitesh/WGS/intel-gatk4-somatic-with-preprocessing/mutect2_nodocker.wdl

module load java/1.8.0_31
ref_fasta="/data/ldiaz/luol2/gatk_resource_bundle/hg38/Homo_sapiens_assembly38.fasta"
ref="/data/ldiaz/luol2/gatk_resource_bundle/hg38"
tool_path="/data/ldiaz/luol2/GATK_workflows/up-to-date-broad-genomes-in-the-cloud-tools"
compression_level="1"
java_heap_memory_initial="2g"
mem_limit="2G"
samtools_threads="32"


data_path=DATA_PATH
project=PROJECT
subject=SUBJECT
sample=TEST_SAMPLE

#tumor_bam=${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicate_marked.sorted.bam
tumor_bam=${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicate_marked.bam

SORT_BAM="${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicate_marked.IndelRealigned.sorted.bam"
SORT_BAI="${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicate_marked.IndelRealigned.sorted.bai"


# index the bam
if [ ! -e ${tumor_bam/.bam/.bai} ];then
module load samtools/1.7
samtools index ${tumor_bam}
fi


# Apply Indelrealign
java -jar /data/ldiaz/luol2/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R ${ref_fasta} \
-I ${tumor_bam} \
-o ${tumor_bam/.bam/.realign_target.intervals} 



java -jar /data/ldiaz/luol2/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R ${ref_fasta} \
-targetIntervals ${tumor_bam/.bam/.realign_target.intervals} \
--noOriginalAlignmentTags \
-I ${tumor_bam} \
-o ${tumor_bam/.bam/.IndelRealigned.bam}

# sort deduplicated & realigned bam
${tool_path}/samtools sort -m ${mem_limit} --threads ${samtools_threads} -l ${compression_level} ${tumor_bam/.bam/.IndelRealigned.bam} -o ${SORT_BAM} 
${tool_path}/samtools index ${SORT_BAM} ${SORT_BAI}




