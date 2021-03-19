# Diaz Lab Human Tumor-only WES Somatic Variant Calling Pipeline (SomVaRIUS)


This is a version-controlled code repository for **Human Tumor-only Somatic Variant Calling Pipeline** in development by the Diaz group at MSKCC. The pipeline implementation is specific for execution by the members in the **Diaz group only!!!**  

The entire pipeline consists of two parts: **`Primary Analysis`** and **`Secondary Analysis`**.  
1. **Primary analysis** takes raw FASTQ files as inputs, follow the procedures described below to generate well-filtered & annotated mutation callset.  
2. **Secondary analysis** picks up the callset and produce multiple summary analyses in terms of both table and plot to cover the routine tasks like betwee-condition variant gain/loss comparison, variant number, VAF distribution by type and mutation signaure deconvolution. In more detail:  




> The **`primary analysis`** was built on top of both **GATK Best Practice** for somatic calling and a high performance Tumor-only variant calling algorithm, **SomVaRIUS**. The sample FASTQ files are aligned and called against the GRCh38 reference genome and its relevant resource bundle. The result set of variants are further manually filtered and VEP annotated to reach a clean, well function-annotated set. 
>
> The **`secondary analysis`** performs downstream summary analysis with the annotated callset. Before its operation, ***RESEARCHERS*** need to *MANUALLY* prepare 3 ***MANIFEST FILES*** to define the analytical scope for 1) `variant comparison by type`; 2) `variant allele frequency (VAF) plot`; 3) `mutation signature`. The template for these manifest files can be found in the tables shown below. A R script summarizes all the secondary analysis outputs and generates a markdown notebook with tables and plots of the final results. 
	  
	  

![GitHub Logo](/images/Mouse_WES_Somatic_Mutation_Calling_Pipeline.png)

_________________________
### **Primary Lilac Locations for the Pipeline and associated Resource Bundles**

* Top-level Github directory of the pipeline: **_/luolingqi/Human_WES_Somatic_TumorOnly_SomVaRIUS_Pipeline_** w/sub-directories: 
  - Main pipeline implementation scripts: **_`/Primary`_**
  - Downstream R analysis for variant gain/loss, MAF distribution plots, Mutation signature, etc.: **_`/Seconcary`_**
  
* Primary Mouse Reference Genome (GRCh38): **_/home/luol2/lingqi_workspace/gatk_resource_bundle/hg38/_** w/contents:
  - Reference Genome files in diverse format: fasta, dict, bwamem index
  - The Agilent SureSelect exome enrichment target bed file
  - Germline SNP/INDEL recorded by various sources (dbsnp)
  
* ENSEMBL VEP resources for both Mouse and Human: **_/home/luol2/lingqi_workspace/vep_data_**

_________________________
### **Primary Scripts & manifest files for Automatic Pipeline Running**
  * `step1_preprocessing_simple.sh` -- it takes fastq file and readgroup info as inputs, and runs the following as listed in the table below

  * `step2_SomVarIUS_VEP_simple.sh` -- it takes the duplicate-removed and BQSR-recalibrated file outputs from `step1_preprocessing_simple.sh`, and runs the following as listed in the table below
  
  * `step3_gain_loss_VEP_simple.sh` -- it takes a manifest sample comparison file (as shown in the template table below) and all the well annotated/filtered variant files from `step2_SomVarIUS_VEP_simple.sh`, output the gain/loss for all the pairwise comparisons (e.g. Parental vs Treatment)
  
  * `step4_MutSig_deconstruction.sh` -- it takes a manifest file (as shown in the template table below) recording vcf files for the individual samples and gain/loss comparisons, deconstruct and plot the mutation signatures
    
  * `items_to_compare.txt` -- as shown in the templeate `Table 1`
  
**Table 1. Description of the workflow on the above 4 steps**
Step1: Data Quality Checking & Preprocessing  |  Step2: Variant Calling, filtering & Annotation | Step3: Variant Gain/Loss Comparison (e.g. Parental v.s. Treated) | Step4: Mutation Signature Deconvolution
-------------------------------------------   |  ---------------------------------------------- |  --------------------------------------------------------------- | ----------------------------------------
Quality checking of raw fastq files <br/> **(Fastqc - run_fastqc.sh)**  |  Somatic Variant Calling and filtration <br/> **(SomVaRIUS - run_somvarius_and_Filter_tumor_only.sh)** | Variant Gain/Loss (Parental vs Treated) <br/> **(run_variants_gain_loss_treatment.vs.Parental.AF.0.05.sh)** | Mutation Signature Deconvolution <br/> **(run_deconstructsigshg38.sh)**
Adapter & low quality reads trimming <br/> **(Trimgalore - run_trim_galore.sh)** |  Removing Germline SNP/INDEL variants <br/> **(For human - run_remove_hg38_germline_dbsnp.sh)**
Trimmed fastq to uBAM format conversion <br/> **(required by GATK pipeline - run_fastq_to_uBAM.sh)**  |  ENSEMBL VEP variant annotation & type filtration <br/> **(missense, frameshit, nonsynonymous, etc. - run_VEP_annotation_hg38_tumor_only_AF_0.05.sh)**
BWA MEM alignment to GRCh38 <br/> **(BWA MEM - run_bwa_mem.sh)**  |  Extra manual filtrations by quality <br/> **(AD, MBQ, MMQ, MPOS5, etc. - run_VEP_annotation_hg38_tumor_only_AF_0.05.sh)**
Alignment quality metrics collection <br/> **(Qualimap - run_qualimap.sh)**  |  
Hybrid selection quality metrics collection <br/> **(GATK HsMetrics - run_CollectHsMetrics.sh)**  |  
Estimate and Apply MarkDuplicate and <br/> Base Quality Score recalibration <br/> **(GATK MarkDuplicate, BQSR  - run_markduplicate.sh)**  |  


**Table 2. Sample Comparison Manifest file as a template**
**Parental Sample** | **Treated Sample**
------------------- | ------------------
Sample_SW480_P8W_IGO_10212_G_2 | Sample_SW480_P0W_IGO_10212_G_1
Sample_SW480_P8W_IGO_10212_G_2 | Sample_SW480_TMZ8W_IGO_10212_G_3
Sample_SW480_P8W_IGO_10212_G_2 | Sample_SW480_CDDP8W_IGO_10212_G_4
Sample_SW480_P8W_IGO_10212_G_2 | Sample_SW480_COMBO8W_IGO_10212_G_5

**Table 3. Mutation Signature Manifest file as a template**
**cats** | **comparison** | **files**
-------- | -------------- | ---------
Sample_SW480_COMBO8W_IGO_10212_G_5 | loss | PITT_0522/Sample_SW480_COMBO8W_IGO_10212_G_5_vs_Sample_SW480_P8W_IGO_10212_G_2.VEP.ann.AF0.05_BQ20_MQ50.vcf
Sample_SW480_TMZ8W_IGO_10212_G_3 | loss | PITT_0522/Sample_SW480_TMZ8W_IGO_10212_G_3_vs_Sample_SW480_P8W_IGO_10212_G_2.VEP.ann.AF0.05_BQ20_MQ50.vcf
Sample_SW480_CDDP8W_IGO_10212_G_4 | gain | PITT_0522/Sample_SW480_P8W_IGO_10212_G_2_vs_Sample_SW480_CDDP8W_IGO_10212_G_4.VEP.ann.AF0.05_BQ20_MQ50.vcf
Sample_SW480_COMBO8W_IGO_10212_G_5 | gain | PITT_0522/Sample_SW480_P8W_IGO_10212_G_2_vs_Sample_SW480_COMBO8W_IGO_10212_G_5.VEP.ann.AF0.05_BQ20_MQ50.vcf
Sample_SW480_CDDP8W_IGO_10212_G_4 | total | PITT_0522/Sample_SW480_CDDP8W_IGO_10212_G_4/Sample_SW480_CDDP8W_IGO_10212_G_4.aligned.duplicates_marked.recalibrated.targeted_sequencing.somvarius.tumor_only.sorted.rm_dbsnps.AF0.05_BQ20_MQ50.VEP.ann.vcf

_________________________
### **Prerequisites for Running the Pipeline**<br/>

* The entire pipeline was built on MSK High Performance Computing (HPC) platform with all the individual building blocks/tools developed in worry-free encapsulated enviroment (Singularity). So, there is little dependency to the system we log in on Lilac, which means, **_anyone with an active Lilac account and basic skill of linux_** can easily run it without any bothering of environment/parameters tuning.
* The input data structure needs to be organized as following, so that the pipeline can locate the pair-end fastq files in gz format in each sample folder.

```
DATA_PATH/
|-- PROJECT/
|   |-- SUBJECT/ # can be missed
|       |-- SAMPLE/
|           |-- *R1*fastq.gz
|           |-- *R2*fastq.gz

```

_________________________
### **Main Pipeline Usage (Primary and Secondary Analysis)**

&nbsp;&nbsp;&nbsp;&nbsp;
The entire pipeline is split into the following four steps for the sake of efficient debugging and implementation. The operator will need to supervise the successful execution of each step by lauching subsequent one. There is one batch script linked to each of the 4 steps respectively as shown below.  

* step1_preprocessing_simple.sh
* step2_SomVarIUS_VEP_simple.sh
* step3_gain_loss_VEP_simple.sh
* step4_MutSig_deconstruction.sh

&nbsp;&nbsp;&nbsp;&nbsp;
Please first copy the above shell scripts into a folder (e.g. data\_analysis) as your working directory to launch these 4 steps! Each step consists of multiple heavy-load tasks, which take time to accomplish. Please be patient to wait till one step ends successfully before kicking off next step. While I don't expect end users to debug any errors which interrupt the pipeline, the pipeline does log the job status and errors in a log file named like "nohup\_step\_\*.log". A note message "Mission Accomplished!" at the end of the log file indicates the success of the step. Make sure you see it before you go to next step.
  
```
  You can write your own loop to run the pipeline for multiple samples. The following USAGE example is for single-sample only.
  
  # preprocessing
  .USAGE.
  nohup sh step1_preprocessing_simple.sh DATA_PATH PROJECT SUBJECT SAMPLE 2>&1 >nohup_step1_SAMPLE.log &
  
  .OPTIONS.
  DATA_PATH  a root directory of the entire study, required.             e.g. /home/luol2/lingqi_workspace/Projects/Ben_Projects
  PROJECT    a project name, required.                                   e.g. WES_mouse_Project_10212_E
  SUBJECT    a subject name if any.                                      e.g. any name here, if no, just use '.'
  SAMPLE     a sample name, required.                                    e.g. Sample_CT26CDDP_M1_IGO_10212_E_13
  
  
  # SomVaRIUS calling & VEP annotation
  nohup sh step2_SomVarIUS_VEP_simple.sh DATA_PATH PROJECT SUBJECT SAMPLE 2>&1 >nohup_step2_SAMPLE.log &
  
  .OPTIONS.
  Same as the options in step 1

  # Variant Gain/Loss analysis (e.g. Parental vs Treated) & VEP annotation 
  nohup sh step3_gain_loss_VEP_simple.sh DATA_PATH PROJECT SUBJECT PARENTAL TREATED 2>&1 >nohup_step3_PARENTAL_vs_TREATED.log &
  
  .OPTIONS.
  Same as the options in step 1, except that you need to specify both PARENTAL and TREATED sample for comparison. Loop through a manifest comparison file as Table 2 above to run the annotated gain/loss analysis as many as you want
  
  # Mutation Signature analysis (deconstructSigs)
  nohup sh step4_MutSig_deconstruction.sh DATA_PATH PROJECT MANIFEST 2>&1 >nohup_step4_MutSig.log
  
  .OPTIONS.
  Same as the option in step 1, except
  MANIFEST      a file recording the samples, variant types whose signatures are to be analyzed, together with file locations to the VCF files           e.g. data_analysis/items_to_plots_MutSig.txt Table 3 above shows an example
  
```
_________________________
## Versioning
For the versions available, see the [tags on this repository](https://github.com/luolingqi/Mouse_WES_Somatic_Primary_Analysis/releases/tag/v0.2-alpha).

## Authors
* **Lingqi Luo, PhD** - initial drafting - [luolingqi](https://github.com/luolingqi) <br/>
See also the list of [contributors](https://github.com/luolingqi/Mouse_WES_Somatic_Primary_Analysis/contributors) who participated in this project.

## License
This project is licensed by the Diaz laboratory at MSKCC. Only Diaz group is authorized to use this pipeline.

## Acknowledgements
* Team DiazLab @ MSKCC

