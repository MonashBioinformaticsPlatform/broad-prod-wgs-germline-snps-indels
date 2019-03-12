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
   >broad-prod-wgs-germline-snps-indels/PairedEndSingleSampleWf.gatk4.0.local.options.json
```

This example uses the provided SLURM backend and specifies some paths for Singularity,
and the SLURM account to submit jobs under.

```bash
cromwell \
  -Dconfig.file=$(pwd)/broad-prod-wgs-germline-snps-indels/slurm-backend.conf \
  -Dsystem.input-read-limits.lines=1000000 \
  run broad-prod-wgs-germline-snps-indels/PairedEndSingleSampleWf.gatk4.0.wdl \
  --inputs broad-prod-wgs-germline-snps-indels/PairedEndSingleSampleWf.hg38.local.inputs.json \
  --options broad-prod-wgs-germline-snps-indels/PairedEndSingleSampleWf.gatk4.0.local.options.json
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

- GATK 4.beta.5 or later for PairedSingleSampleWF. GATK4.0.1.0 or later for JointGenotypingWf.
- Picard 2.x
- Samtools (see gotc docker)
- Python 2.7

Cromwell version support 
- Successfully tested on v31
- Does not work on versions < v23 due to output syntax
