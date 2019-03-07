# prod-wgs-germline-snps-indels
### Purpose : 
Workflows used in production at Broad for germline short variant discovery in WGS data.

This repo contains a version modified by the *Monash Bioinformatics Platform* for running on
local HPC resources.

### Installing the workflow and cromwell

```bash
git clone --recursive https://github.com/MonashBioinformaticsPlatform/broad-prod-wgs-germline-snps-indels.git
# Create a conda environment called 'cromwell', and install cromwell into it
conda create -n cromwell -c bioconda -c conda-forge cromwell
conda activate cromwell
pip install gsutil
pip install j2cli
```

### Getting the references and example data
```bash
export REF_BASE=/scratch/pl41/references
mkdir -p $REF_BASE/references/broad-references/
cd $REF_BASE/references/broad-references/
gsutil -m rsync -r -x "CrossSpeciesContamination/*" gs://broad-references/hg38 hg38
```

For PairedEndSingleSampleWf:
```bash
mkdir -p $REF_BASE/references/broad-public-datasets
cd $REF_BASE
gsutil -m rsync -r gs://broad-public-datasets/NA12878_downsampled_for_testing broad-public-datasets/NA12878_downsampled_for_testing

mkdir -p $REF_BASE/references/dsde-data-na12878-public
wget -O dsde-data-na12878-public/NA12878.hg38.reference.fingerprint.vcf \
     https://storage.googleapis.com/dsde-data-na12878-public/NA12878.hg38.reference.fingerprint.vcf
wget -O dsde-data-na12878-public/NA12878.hg38.reference.fingerprint.vcf.idx \
     https://storage.googleapis.com/dsde-data-na12878-public/NA12878.hg38.reference.fingerprint.vcf.idx
```

For JointGenotypingWf:
```bash
mkdir -p $REF_BASE/references/gatk-test-data
cd $REF_BASE
gsutil -m rsync -r gs://gatk-test-data/joint_discovery/NA12878.sample_map gatk-test-data/joint_discovery/
gsutil -m rsync -r gs://gatk-test-data/intervals/hg38.even.handcurated.20k.intervals gatk-test-data/intervals/hg38.even.handcurated.20k.intervals
```

### PairedSingleSampleWF :
This WDL pipeline implements data pre-processing and initial variant calling (GVCF
generation) according to the GATK Best Practices (June 2016) for germline SNP and
Indel discovery in human whole-genome sequencing (WGS) data.

*Note: PairedEndSingleSampleWf-fc is a modified version of the original workflow
that is easier to run on FireCloud/Terra*

#### Requirements/expectations
- Human whole-genome pair-end sequencing data in unmapped BAM (uBAM) format
- One or more read groups, one per uBAM file, all belonging to a single sample (SM)
- Input uBAM files must additionally comply with the following requirements:
- - filenames all have the same suffix (we use ".unmapped.bam")
- - files must pass validation by ValidateSamFile
- - reads are provided in query-sorted order
- - all reads must have an RG tag
- Reference genome must be Hg38 with ALT contigs

#### Outputs 
- Cram, cram index, and cram md5 
- GVCF and its gvcf index 
- BQSR Report
- Several Summary Metrics 

#### Running (example)

Generate an `inputs.json` file based on the template:
```bash
export REF_BASE=/scratch/pl41/references
export SLURM_ACCOUNT=pl41
export SINGULARITY_CACHE=/scratch/${SLURM_ACCOUNT}/singularity_cache
export SINGULARITY_BIND_PATH=/scratch
export TMPDIR=/scratch/${SLURM_ACCOUNT}/tmp

j2 broad-prod-wgs-germline-snps-indels/PairedEndSingleSampleWf.hg38.inputs.json.j2 \
   >broad-prod-wgs-germline-snps-indels/PairedEndSingleSampleWf.hg38.local.inputs.json

j2 broad-prod-wgs-germline-snps-indels/PairedEndSingleSampleWf.gatk4.0.options.json.j2 \
   >broad-prod-wgs-germline-snps-indels/PairedEndSingleSampleWf.local.options.json
```

This example uses the provided SLURM backend and specifies some paths for Singularity,
and the SLURM account to submit jobs under.

```bash
cromwell \
  -Dconfig.file=$(pwd)/broad-prod-wgs-germline-snps-indels/slurm-backend.conf \
  -Dsystem.input-read-limits.lines=1000000 \
  run broad-prod-wgs-germline-snps-indels/PairedEndSingleSampleWf.gatk4.0.wdl \
  --inputs broad-prod-wgs-germline-snps-indels/PairedEndSingleSampleWf.hg38.local.inputs.json \
  --options broad-prod-wgs-germline-snps-indels/PairedEndSingleSampleWf.local.options.json
```

### JointGenotypingWf :
The second WDL implements the joint discovery and VQSR 
filtering portion of the GATK Best Practices (June 2016) for germline SNP and Indel 
discovery in human whole-genome sequencing (WGS) and exome sequencing data.

#### Requirements/expectations
- One or more GVCFs produced by HaplotypeCaller in GVCF mode.
- Bare minimum 1 WGS sample or 30 Exome samples. Gene panels are not supported.

#### Outputs 
- VCF  and its vcf index
 *Note: The gvcf is filtered using variant quality score recalibration  
  (VQSR) with genotypes for all samples present in the input VCF. All sites that  
  are present in the input VCF are retained; filtered sites are annotated as such  
  in the FILTER field.*
- Summary Metrics

#### Running (example)

Generate an `inputs.json` and `options.json` file based on the template:
```bash
export REF_BASE=/scratch/pl41/references
export SLURM_ACCOUNT=pl41
export SINGULARITY_CACHE=/scratch/${SLURM_ACCOUNT}/singularity_cache
export SINGULARITY_BIND_PATH=/scratch
export TMPDIR=/scratch/${SLURM_ACCOUNT}/tmp

j2 broad-prod-wgs-germline-snps-indels/JointGenotypingWf.hg38.inputs.json.j2 \
   >broad-prod-wgs-germline-snps-indels/JointGenotypingWf.hg38.local.inputs.json

j2 broad-prod-wgs-germline-snps-indels/JointGenotypingWf.options.json.j2 \
   >broad-prod-wgs-germline-snps-indels/JointGenotypingWf.local.options.json
```

This example uses the provided SLURM backend and specifies some paths for Singularity,
and the SLURM account to submit jobs under.

```bash
cromwell -Dconfig.file=$(pwd)/broad-prod-wgs-germline-snps-indels/slurm-backend.conf \
         -Dsystem.input-read-limits.lines=1000000 \
         run broad-prod-wgs-germline-snps-indels/JointGenotypingWf.wdl \
         --inputs broad-prod-wgs-germline-snps-indels/JointGenotypingWf.hg38.local.inputs.json \
         --options broad-prod-wgs-germline-snps-indels/JointGenotypingWf.local.options.json
         
# --options DOES seem to work for JointGenotypingWf
```

### Software version requirements :
These workflows user Docker/Singularity containers pulled on demand, so specific local installations aren't required 
(you will need `cromwell` and `singularity` or `docker`).

- GATK 4.beta.5 for PairedSingleSampleWF and GATK 4.0 for PairedSingleSampleWF-fc. GATK4.0.11.0 for JointGenotypingWf.
- Picard 2.16.0-SNAPSHOT
- Samtools 1.3.1
- Python 2.7
- Cromwell version support 
  - Successfully tested on v37
  - Does not work on versions < v23 due to output syntax

### Important Note :
- VQSR wiring. The SNP and INDEL models are built in parallel, but then the corresponding 
  recalibrations are applied in series. Because the INDEL model is generally ready 
  first (because there are fewer indels than SNPs) we set INDEL recalibration to 
  be applied first to the input VCF, while the SNP model is still being built. By 
  the time the SNP model is available, the indel-recalibrated file is available to 
  serve as input to apply the SNP recalibration. If we did it the other way around, 
  we would have to wait until the SNP recal file was available despite the INDEL 
  recal file being there already, then apply SNP recalibration, then apply INDEL 
  recalibration. This would lead to a longer wall clock time for complete workflow 
  execution. Wiring the INDEL recalibration to be applied first solves the problem.
- The current version of the posted "Generic germline short variant joint genotyping" 
  is derived from the Broad production version of the workflow, which was adapted for 
  large WGS callsets of up to 20K samples.  We believe the results of this workflow run 
  on a single WGS sample are equally accurate, but there may be some shortcomings when 
  the workflow is modified and run on small cohorts.  Specifically, modifying the SNP 
  ApplyRecalibration step for higher specificity may not be effective.  The user can verify 
  if this is an issue by consulting the gathered SNP tranches file.  If the listed 
  truthSensitivity in the rightmost column is not well matched to the targetTruthSensitivity 
  in the leftmost column, then requesting that targetTruthSensitivity from ApplyVQSR will 
  not use an accurate filtering threshold.  This workflow has not been tested on exomes.  
  The dynamic scatter interval creating was optimized for genomes.  The scattered SNP 
  VariantRecalibration may fail because of too few "bad" variants to build the negative model. 
  Also, apologies that the logging for SNP recalibration is overly verbose.
- The provided JSON is meant to be a ready to use example JSON template of the workflow. It is the user’s responsibility to correctly set the reference and resource input variables using the [GATK Tool and Tutorial Documentations](https://software.broadinstitute.org/gatk/documentation/).
- Relevant reference and resources bundles can be accessed in [Resource Bundle](https://software.broadinstitute.org/gatk/download/bundle).
- Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
- For help running workflows on the Google Cloud Platform or locally please
view the following tutorial [(How to) Execute Workflows from the gatk-workflows Git Organization](https://software.broadinstitute.org/gatk/documentation/article?id=12521).
- The following material is provided by the GATK Team. Please post any questions or concerns to one of our forum sites : [GATK](https://gatkforums.broadinstitute.org/gatk/categories/ask-the-team/) , [FireCloud](https://gatkforums.broadinstitute.org/firecloud/categories/ask-the-firecloud-team) or [Terra](https://broadinstitute.zendesk.com/hc/en-us/community/topics/360000500432-General-Discussion) , [WDL/Cromwell](https://gatkforums.broadinstitute.org/wdl/categories/ask-the-wdl-team).
- Please visit the [User Guide](https://software.broadinstitute.org/gatk/documentation/) site for further documentation on our workflows and tools.

### LICENSING :
Copyright Broad Institute, 2019 | BSD-3
This script is released under the WDL open source code license (BSD-3) (full license text at https://github.com/openwdl/wdl/blob/master/LICENSE). Note however that the programs it calls may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running this script.
- [GATK](https://software.broadinstitute.org/gatk/download/licensing.php)
- [BWA](http://bio-bwa.sourceforge.net/bwa.shtml#13)
- [Picard](https://broadinstitute.github.io/picard/)
- [Samtools](http://www.htslib.org/terms/)
