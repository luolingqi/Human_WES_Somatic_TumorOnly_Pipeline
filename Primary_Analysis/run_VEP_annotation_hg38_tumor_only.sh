#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 1 -R "rusage[mem=8]"
#BSUB -W 48:00

module load singularity/3.1.1
module load java/1.8.0_31

project=PROJECT
subject=SUBJECT
sample=TEST_SAMPLE

#data_path="/data/ldiaz/luol2/Projects/Jones_Lee_Exercise_Project"
data_path=DATA_PATH
ref_fasta="/data/ldiaz/luol2/gatk_resource_bundle/hg38/Homo_sapiens_assembly38.fasta"


VEPInputVcf=${sample}.aligned.duplicates_marked.recalibrated.targeted_sequencing.tumor_only.gatk4_mutect2.raw_somatic_mutation.rm_dbsnps.singleSample.PASSOnly.vcf
VEPOutputVcf=${VEPInputVcf/.vcf/.VEP.ann.vcf}
VEPOutputStat=${VEPInputVcf/.vcf/VEP.ann.summary.html}

#Prepare single sample (tumor only) vcf for VEP; Also subset to get only the "PASS" variants
tumor_sample_name=$(cat ${data_path}/${project}/${subject}/${sample}/tumor_name.txt)

# add one more filter to:
# 1) remove variants with DP in both tumor and normal (e.g.: DP >= 10)

#TumorID=$(bcftools query -l ${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicates_marked.recalibrated.targeted_sequencing.tumor_only.gatk4_mutect2.raw_somatic_mutation.rm_dbsnps.vcf.gz | nl -v 0 | grep "${tumor_sample_name}" | cut -f1)
#echo -e "TumorID is \t${TumorID}"

cp ${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicates_marked.recalibrated.targeted_sequencing.tumor_only.gatk4_mutect2.raw_somatic_mutation.rm_dbsnps.vcf.gz ${data_path}/${project}/${subject}/${sample}/temp.vcf.gz

#bcftools view -i "FORMAT/DP >= 10" -Oz -o ${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicates_marked.recalibrated.targeted_sequencing.tumor_only.gatk4_mutect2.raw_somatic_mutation.rm_dbsnps.vcf.gz ${data_path}/${project}/${subject}/${sample}/temp.vcf.gz 


java -jar /data/ldiaz/luol2/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T SelectVariants -R ${ref_fasta} -V ${data_path}/${project}/${subject}/${sample}/${sample}.aligned.duplicates_marked.recalibrated.targeted_sequencing.tumor_only.gatk4_mutect2.raw_somatic_mutation.rm_dbsnps.vcf.gz -o ${data_path}/${project}/${subject}/${sample}/${VEPInputVcf} -sn ${tumor_sample_name} --excludeFiltered


#Perform VEP using local cache (Mouse)
singularity run --bind ${data_path}/${project}/${subject}/${sample}:/data --bind /home/luol2/lingqi_workspace/vep_data:/vep_data \
/data/ldiaz/luol2/ensembl-vep.sif \
vep -i /data/${VEPInputVcf} \
    -o /data/${VEPOutputVcf} \
    --vcf \
    --force_overwrite \
    --cache \
    --dir_cache /vep_data \
    --species homo_sapiens \
    --stats_file /data/${VEPOutputStat}

#Perform VEP summarization after manual filtration
bcftools view -i 'FORMAT/DP[0:0] > 10 & FORMAT/AD[0:1] > 4 & FORMAT/AF[0:0] > 0.1 & INFO/MBQ[1] > 20 & INFO/MMQ[1] > 50 & INFO/MPOS > 5' ${data_path}/${project}/${subject}/${sample}/${VEPInputVcf} > ${data_path}/${project}/${subject}/${sample}/${VEPInputVcf/.vcf/.altAD4_AF0.1_MBQ20_MMQ50_MPOS5.vcf} 
singularity exec --bind ${data_path}/${project}/${subject}/${sample}:/data --bind /home/luol2/lingqi_workspace/vep_data:/vep_data \
/data/ldiaz/luol2/ensembl-vep.sif \
vep -i /data/${VEPInputVcf/.vcf/.altAD4_AF0.1_MBQ20_MMQ50_MPOS5.vcf} \
    -o /data/${VEPInputVcf/.vcf/.altAD4_AF0.1_MBQ20_MMQ50_MPOS5.VEP.ann.vcf} \
    --vcf \
    --force_overwrite \
    --cache \
    --dir_cache /vep_data \
    --species mouse \
    --stats_file /data/${VEPInputVcf/.vcf/.altAD4_AF0.1_MBQ20_MMQ50_MPOS5.VEP.ann.summary.html}

# VEP output in tab format
singularity exec --bind ${data_path}/${project}/${subject}/${sample}:/data --bind /home/luol2/lingqi_workspace/vep_data:/vep_data \
/data/ldiaz/luol2/ensembl-vep.sif \
vep -i /data/${VEPInputVcf/.vcf/.altAD4_AF0.1_MBQ20_MMQ50_MPOS5.vcf} \
    -o /data/${VEPInputVcf/.vcf/.altAD4_AF0.1_MBQ20_MMQ50_MPOS5.VEP.ann.txt} \
    --tab \
    --fields "Uploaded_variation,Location,Allele,Gene,CANONICAL,SYMBOL,SYMBOL_SOURCE,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,DISTANCE,STRAND" \
    --force_overwrite \
    --cache \
    --dir_cache /vep_data \
    --species homo_sapiens \
    --symbol


# Perform VEP filtrations:
singularity exec --bind ${data_path}/${project}/${subject}/${sample}:/data \
/data/ldiaz/luol2/ensembl-vep.sif \
filter_vep -i /data/${VEPOutputVcf} --filter "Consequence in missense_variant,frameshift_variant,inframe_deletion,coding_sequence_variant,stop_gained,inframe_insertion,splice_region_variant" | \
bcftools view -i 'FORMAT/DP[0:0] > 10 & FORMAT/AD[0:1] > 4 & FORMAT/AF[0:0] > 0.1 & INFO/MBQ[1] > 20 & INFO/MMQ[1] > 50 & INFO/MPOS > 5' > ${data_path}/${project}/${subject}/${sample}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.1_MBQ20_MMQ50_MPOS5.non_synonymous_variant_only.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t[%GT\t%AD{0}\t%AD{1}\t%AF\t]\t%MBQ[1]\t%MMQ[1]\t%MPOS\t\n" ${data_path}/${project}/${subject}/${sample}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.1_MBQ20_MMQ50_MPOS5.non_synonymous_variant_only.vcf} > ${data_path}/${project}/${subject}/${sample}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.1_MBQ20_MMQ50_MPOS5.non_synonymous_variant_only.query}

singularity exec --bind ${data_path}/${project}/${subject}/${sample}:/data \
/data/ldiaz/luol2/ensembl-vep.sif \
filter_vep -i /data/${VEPOutputVcf} --filter "Consequence is missense_variant" | \
bcftools view -i 'FORMAT/DP[0:0] > 10 & FORMAT/AD[0:1] > 4 & FORMAT/AF[0:0] > 0.1 & INFO/MBQ[1] > 20 & INFO/MMQ[1] > 50 & INFO/MPOS > 5' > ${data_path}/${project}/${subject}/${sample}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.1_MBQ20_MMQ50_MPOS5.missense_variant_only.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t[%GT\t%AD{0}\t%AD{1}\t%AF\t]\t%MBQ[1]\t%MMQ[1]\t%MPOS\t\n" ${data_path}/${project}/${subject}/${sample}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.1_MBQ20_MMQ50_MPOS5.missense_variant_only.vcf} > ${data_path}/${project}/${subject}/${sample}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.1_MBQ20_MMQ50_MPOS5.missense_variant_only.query}

singularity exec --bind ${data_path}/${project}/${subject}/${sample}:/data \
/data/ldiaz/luol2/ensembl-vep.sif \
filter_vep -i /data/${VEPOutputVcf} --filter "Consequence is frameshift_variant" | \
bcftools view -i 'FORMAT/DP[0:0] > 10 & FORMAT/AD[0:1] > 4 & FORMAT/AF[0:0] > 0.1 & INFO/MBQ[1] > 20 & INFO/MMQ[1] > 50 & INFO/MPOS > 5' > ${data_path}/${project}/${subject}/${sample}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.1_MBQ20_MMQ50_MPOS5.frameshift_variant_only.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t[%GT\t%AD{0}\t%AD{1}\t%AF\t]\t%MBQ[1]\t%MMQ[1]\t%MPOS\t\n" ${data_path}/${project}/${subject}/${sample}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.1_MBQ20_MMQ50_MPOS5.frameshift_variant_only.vcf} > ${data_path}/${project}/${subject}/${sample}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.1_MBQ20_MMQ50_MPOS5.frameshift_variant_only.query}

singularity exec --bind ${data_path}/${project}/${subject}/${sample}:/data \
/data/ldiaz/luol2/ensembl-vep.sif \
filter_vep -i /data/${VEPOutputVcf} --filter "Consequence is synonymous_variant" | \
bcftools view -i 'FORMAT/DP[0:0] > 10 & FORMAT/AD[0:1] > 4 & FORMAT/AF[0:0] > 0.1 & INFO/MBQ[1] > 20 & INFO/MMQ[1] > 50 & INFO/MPOS > 5' > ${data_path}/${project}/${subject}/${sample}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.1_MBQ20_MMQ50_MPOS5.synonymous_variant_only.vcf}
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT{0}\t%DP\t[%GT\t%AD{0}\t%AD{1}\t%AF\t]\t%MBQ[1]\t%MMQ[1]\t%MPOS\t\n" ${data_path}/${project}/${subject}/${sample}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.1_MBQ20_MMQ50_MPOS5.synonymous_variant_only.vcf} > ${data_path}/${project}/${subject}/${sample}/${VEPOutputVcf/.ann.vcf/.ann.altAD4_AF0.1_MBQ20_MMQ50_MPOS5.synonymous_variant_only.query}
