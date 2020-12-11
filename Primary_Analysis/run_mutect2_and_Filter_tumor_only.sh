#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 8 -R "rusage[mem=8]"
#BSUB -W 48:00

### Refer to  ~/lingqi_workspace/Projects/Mitesh/WGS/intel-gatk4-somatic-with-preprocessing/mutect2_nodocker.wdl

module load java/1.8.0_31
ref_fasta="/data/ldiaz/luol2/gatk_resource_bundle/hg38/Homo_sapiens_assembly38.fasta"
ref="/data/ldiaz/luol2/gatk_resource_bundle/hg38"
mutect_ref="/data/ldiaz/luol2/gatk_resource_bundle/mutect2"
tool_path="/data/ldiaz/luol2/GATK_workflows/up-to-date-broad-genomes-in-the-cloud-tools"

data_path=DATA_PATH
project=PROJECT
subject=SUBJECT
sample=TEST_SAMPLE

tumor_bam=${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicates_marked.recalibrated.bam
tumor_name=${data_path}/${project}/${subject}/${sample}/tumor_name.txt


unfiltered_output_vcf=${tumor_bam/.bam/.unfiltered.vcf}
output_vcf=${tumor_bam/.bam/.targeted_sequencing.mutect2.tumor_only.contFiltered.vcf.gz}
command_mem="8000"

# index the bam
if [ ! -e ${tumor_bam/.bam/.bai} ];then
module load samtools/1.7
samtools index ${tumor_bam}
fi


# Run Mutect2

        # We need to create these files regardless, even if they stay empty
        touch ${data_path}/${project}/${subject}/${sample}/bamout.bam


        /data/ldiaz/luol2/gatk-4.1.2.0/gatk --java-options "-Xmx${command_mem}m" GetSampleName \
	    -R ${ref_fasta} \
	    -I ${tumor_bam} \
	    -O ${tumor_name} \
	    -encode
        
	tumor_command_line="-I ${tumor_bam} -tumor `cat ${tumor_name}`"


        /data/ldiaz/luol2/gatk-4.1.2.0/gatk --java-options "-Xmx${command_mem}m" Mutect2 \
            -R ${ref_fasta} \
            $tumor_command_line \
            -L ${ref}/S31285117_Regions.interval_list \
            -O "${unfiltered_output_vcf}" \
	    -pairHMM AVX_LOGLESS_CACHING \
            --native-pair-hmm-threads 1 \
            --smith-waterman AVX_ENABLED \
	    --af-of-alleles-not-in-resource 2.5e-06 \
	    --germline-resource ${mutect_ref}/af-only-gnomad.hg38.vcf.gz \
	    -pon ${mutect_ref}/gatk4_mutect2_4136_pon.vcf.gz
                        
           
## Sort VCF with Picard

java -Xmx${command_mem}G -jar ${tool_path}/picard.jar \
SortVcf \
SEQUENCE_DICTIONARY=${ref}/Homo_sapiens_assembly38.dict \
OUTPUT=${unfiltered_output_vcf/.vcf/.targeted_sequencing.mutect2.tumor_only.sorted.vcf.gz} \
I=${unfiltered_output_vcf} \
CREATE_INDEX=true


cp ${unfiltered_output_vcf}.stats ${unfiltered_output_vcf/.vcf/.targeted_sequencing.mutect2.tumor_only.sorted.vcf.gz}.stats

# Run variant Filtration
        /data/ldiaz/luol2/gatk-4.1.2.0/gatk --java-options "-Xmx${command_mem}m" FilterMutectCalls \
	     -V ${unfiltered_output_vcf/.vcf/.targeted_sequencing.mutect2.tumor_only.sorted.vcf.gz} \
             -R ${ref_fasta} \
             -O ${output_vcf} \
	     --contamination-table ${data_path}/${project}/${subject}/${sample}/${sample}.targeted_sequencing.contamination.table \
	     -L ${ref}/S31285117_Regions.interval_list

# Filter variants by orientation bias
/data/ldiaz/luol2/gatk-4.1.2.0/gatk --java-options "-Xmx${command_mem}m" FilterByOrientationBias \
-O ${tumor_bam/.bam/.targeted_sequencing.tumor_only.gatk4_mutect2.raw_somatic_mutation.vcf.gz} \
-P ${data_path}/${project}/${subject}/${sample}/${sample}.pre_adapter_detail_metrics.txt \
-V ${output_vcf} \
-R ${ref_fasta} \
-L ${ref}/S31285117_Regions.interval_list \
-AM G/T \
-AM C/T

#-L ${ref}/S31285117_Regions.interval_list \ 

