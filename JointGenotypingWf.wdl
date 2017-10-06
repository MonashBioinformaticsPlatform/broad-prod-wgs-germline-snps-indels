## Copyright Broad Institute, 2017
##
## This WDL pipeline implements joint calling across multiple samples according to the 
## GATK Best Practices (June 2016) for germline SNP and Indel discovery in human whole-genome 
## sequencing (WGS) data.
##
## Requirements/expectations :
## - One or more human whole-genome per-sample GVCF files
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

# WORKFLOW DEFINITION
workflow JointGenotyping {
  File unpadded_intervals_file
  Array[String] unpadded_intervals = read_lines(unpadded_intervals_file)

  String callset_name
  File sample_name_map
  
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  File dbsnp_vcf
  File dbsnp_vcf_index

  Int small_disk
  Int medium_disk
  Int huge_disk

  Array[String] snp_recalibration_tranche_values
  Array[String] snp_recalibration_annotation_values
  Array[String] indel_recalibration_tranche_values
  Array[String] indel_recalibration_annotation_values

  File eval_interval_list
  File hapmap_resource_vcf
  File hapmap_resource_vcf_index
  File omni_resource_vcf
  File omni_resource_vcf_index
  File one_thousand_genomes_resource_vcf
  File one_thousand_genomes_resource_vcf_index
  File mills_resource_vcf
  File mills_resource_vcf_index
  File axiomPoly_resource_vcf
  File axiomPoly_resource_vcf_index
  File dbsnp_resource_vcf = dbsnp_vcf
  File dbsnp_resource_vcf_index = dbsnp_vcf_index

  # ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
  # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
  Float excess_het_threshold = 54.69
  Float snp_filter_level
  Float indel_filter_level
  Int SNP_VQSR_downsampleFactor

  scatter (idx in range(length(unpadded_intervals))) {

    # the batch_size value was carefully chosen here as it
    # is the optimal value for the amount of memory allocated
    # within the task; please do not change it without consulting
    # the Hellbender (GATK engine) team!
    call ImportGVCFs {
      input:
        sample_name_map = sample_name_map,
        interval = unpadded_intervals[idx],
        workspace_dir_name = "genomicsdb",
        disk_size = medium_disk,
        batch_size = 50
    }

    call GenotypeGVCFs {
      input:
        workspace_tar = ImportGVCFs.output_genomicsdb,
        interval = unpadded_intervals[idx],
        output_vcf_filename = "output.vcf.gz",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        dbsnp_vcf = dbsnp_vcf,
        disk_size = medium_disk
    }

    call HardFilterAndMakeSitesOnlyVcf {
      input:
        vcf = GenotypeGVCFs.output_vcf,
        vcf_index = GenotypeGVCFs.output_vcf_index,
        excess_het_threshold = excess_het_threshold,
        variant_filtered_vcf_filename = callset_name + "." + idx + ".variant_filtered.vcf.gz",
        sites_only_vcf_filename = callset_name + "." + idx + ".sites_only.variant_filtered.vcf.gz",
        disk_size = medium_disk
    }
  }

  call GatherVcfs as SitesOnlyGatherVcf {
    input:
      input_vcfs_fofn = write_lines(HardFilterAndMakeSitesOnlyVcf.sites_only_vcf),
      input_vcf_indexes_fofn = write_lines(HardFilterAndMakeSitesOnlyVcf.sites_only_vcf_index),
      output_vcf_name = callset_name + ".sites_only.vcf.gz",
      disk_size = medium_disk
  }

  call SNPsVariantRecalibratorCreateModel {
    input:
      sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
      sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
      recalibration_filename = callset_name + ".snps.recal",
      tranches_filename = callset_name + ".snps.tranches",
      recalibration_tranche_values = snp_recalibration_tranche_values,
      recalibration_annotation_values = snp_recalibration_annotation_values,
      downsampleFactor = SNP_VQSR_downsampleFactor,
      model_report_filename = callset_name + ".snps.model.report",
      hapmap_resource_vcf = hapmap_resource_vcf,
      hapmap_resource_vcf_index = hapmap_resource_vcf_index,
      omni_resource_vcf = omni_resource_vcf,
      omni_resource_vcf_index = omni_resource_vcf_index,
      one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
      one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
      dbsnp_resource_vcf = dbsnp_resource_vcf,
      dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
      disk_size = small_disk
  }

  call IndelsVariantRecalibrator {
    input:
      sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
      sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
      recalibration_filename = callset_name + ".indels.recal",
      tranches_filename = callset_name + ".indels.tranches",
      recalibration_tranche_values = indel_recalibration_tranche_values,
      recalibration_annotation_values = indel_recalibration_annotation_values,
      mills_resource_vcf = mills_resource_vcf,
      mills_resource_vcf_index = mills_resource_vcf_index,
      axiomPoly_resource_vcf = axiomPoly_resource_vcf,
      axiomPoly_resource_vcf_index = axiomPoly_resource_vcf_index,
      dbsnp_resource_vcf = dbsnp_resource_vcf,
      dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
      disk_size = small_disk
  }

  scatter (idx in range(length(HardFilterAndMakeSitesOnlyVcf.sites_only_vcf))) {

    call SNPsVariantRecalibratorScattered {
        input:
          sites_only_variant_filtered_vcf = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf[idx],
          sites_only_variant_filtered_vcf_index = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf_index[idx],
          recalibration_filename = callset_name + ".snps." + idx + ".recal",
          tranches_filename = callset_name + ".snps." + idx + ".tranches",
          recalibration_tranche_values = snp_recalibration_tranche_values,
          recalibration_annotation_values = snp_recalibration_annotation_values,
          model_report = SNPsVariantRecalibratorCreateModel.model_report,
          hapmap_resource_vcf = hapmap_resource_vcf,
          hapmap_resource_vcf_index = hapmap_resource_vcf_index,
          omni_resource_vcf = omni_resource_vcf,
          omni_resource_vcf_index = omni_resource_vcf_index,
          one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
          one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
          dbsnp_resource_vcf = dbsnp_resource_vcf,
          dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
          disk_size = small_disk
      }
  }

  call GatherTranches as SNPGatherTranches {
  input:
       input_fofn = write_lines(SNPsVariantRecalibratorScattered.tranches),
       output_filename = callset_name + ".snps.gathered.tranches",
       disk_size = small_disk
  }

  scatter (idx in range(length(HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf))) {
    call ApplyRecalibration {
      input:
        recalibrated_vcf_filename = callset_name + ".filtered." + idx + ".vcf.gz",
        input_vcf = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf[idx],
        input_vcf_index = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index[idx],
        indels_recalibration = IndelsVariantRecalibrator.recalibration,
        indels_recalibration_index = IndelsVariantRecalibrator.recalibration_index,
        indels_tranches = IndelsVariantRecalibrator.tranches,
        snps_recalibration = SNPsVariantRecalibratorScattered.recalibration[idx],
        snps_recalibration_index = SNPsVariantRecalibratorScattered.recalibration_index[idx],
        snps_tranches = SNPGatherTranches.tranches,
        indel_filter_level = indel_filter_level,
        snp_filter_level = snp_filter_level,
        disk_size = medium_disk
    }

    call CollectVariantCallingMetrics {
      input:
	input_vcf = ApplyRecalibration.recalibrated_vcf,
      	input_vcf_index = ApplyRecalibration.recalibrated_vcf_index,
      	metrics_filename_prefix = callset_name + "." + idx,
      	dbsnp_vcf = dbsnp_vcf,
      	dbsnp_vcf_index = dbsnp_vcf_index,
      	interval_list = eval_interval_list,
      	ref_dict = ref_dict,
      	disk_size = small_disk
    }
  }

  call GatherVcfs as FinalGatherVcf {
    input:
      input_vcfs_fofn = write_lines(ApplyRecalibration.recalibrated_vcf),
      input_vcf_indexes_fofn = write_lines(ApplyRecalibration.recalibrated_vcf_index),
      output_vcf_name = callset_name + ".vcf.gz",
      disk_size = huge_disk
  }

  call GatherMetrics {
    input:
      input_details_fofn = write_lines(CollectVariantCallingMetrics.detail_metrics_file),
      input_summaries_fofn = write_lines(CollectVariantCallingMetrics.summary_metrics_file),
      output_prefix = callset_name,
      disk_size = medium_disk
  }

  output {
    FinalGatherVcf.*
    GatherMetrics.*
  }
}

task ImportGVCFs {
  File sample_name_map
  String interval

  String workspace_dir_name

  Int disk_size
  Int batch_size

  command <<<

    set -e

    rm -rf ${workspace_dir_name}

    # The memory setting here is very important and must be several GB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    /usr/gitc/gatk-launch --javaOptions "-Xmx4g -Xms4g" \
    GenomicsDBImport \
    --genomicsDBWorkspace ${workspace_dir_name} \
    --batchSize ${batch_size} \
    -L ${interval} \
    --sampleNameMap ${sample_name_map} \
    --readerThreads 5 \
    -ip 500

    tar -cf ${workspace_dir_name}.tar ${workspace_dir_name}

  >>>
  runtime {
    memory: "7 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_genomicsdb = "${workspace_dir_name}.tar"
  }
}

task GenotypeGVCFs {
  File workspace_tar
  String interval

  String output_vcf_filename

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  String dbsnp_vcf

  Int disk_size

  command <<<

    set -e

    tar -xf ${workspace_tar}
    WORKSPACE=$( basename ${workspace_tar} .tar)

    /usr/gitc/gatk-launch --javaOptions "-Xmx5g -Xms5g" \
     GenotypeGVCFs \
     -R ${ref_fasta} \
     -O ${output_vcf_filename} \
     -D ${dbsnp_vcf} \
     -G StandardAnnotation \
     --onlyOutputCallsStartingInIntervals \
     -newQual \
     -V gendb://$WORKSPACE \
     -L ${interval}
  >>>
  runtime {
    memory: "7 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
  }
  output {
    File output_vcf = "${output_vcf_filename}"
    File output_vcf_index = "${output_vcf_filename}.tbi"
  }
}

task HardFilterAndMakeSitesOnlyVcf {
  File vcf
  File vcf_index
  Float excess_het_threshold

  String variant_filtered_vcf_filename
  String sites_only_vcf_filename

  Int disk_size

  command {
    set -e

    /usr/gitc/gatk-launch --javaOptions "-Xmx3g -Xms3g" \
      VariantFiltration \
      --filterExpression "ExcessHet > ${excess_het_threshold}" \
      --filterName ExcessHet \
      -O ${variant_filtered_vcf_filename} \
      -V ${vcf}

    java -Xmx3g -Xms3g -jar /usr/gitc/picard.jar \
      MakeSitesOnlyVcf \
      INPUT=${variant_filtered_vcf_filename} \
      OUTPUT=${sites_only_vcf_filename}

  }
  runtime {
    memory: "3.5 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
  }
  output {
    File variant_filtered_vcf = "${variant_filtered_vcf_filename}"
    File variant_filtered_vcf_index = "${variant_filtered_vcf_filename}.tbi"
    File sites_only_vcf = "${sites_only_vcf_filename}"
    File sites_only_vcf_index = "${sites_only_vcf_filename}.tbi"
  }
}

task IndelsVariantRecalibrator {
  String recalibration_filename
  String tranches_filename

  Array[String] recalibration_tranche_values
  Array[String] recalibration_annotation_values

  File sites_only_variant_filtered_vcf
  File sites_only_variant_filtered_vcf_index

  File mills_resource_vcf
  File axiomPoly_resource_vcf
  File dbsnp_resource_vcf
  File mills_resource_vcf_index
  File axiomPoly_resource_vcf_index
  File dbsnp_resource_vcf_index

  Int disk_size

  command {
    /usr/gitc/gatk-launch --javaOptions "-Xmx24g -Xms24g" \
      VariantRecalibrator \
      -V ${sites_only_variant_filtered_vcf} \
      -O ${recalibration_filename} \
      --tranchesFile ${tranches_filename} \
      -allPoly \
      -tranche ${sep=' -tranche ' recalibration_tranche_values} \
      -an ${sep=' -an ' recalibration_annotation_values} \
      -mode INDEL \
      --maxGaussians 4 \
      -resource mills,known=false,training=true,truth=true,prior=12:${mills_resource_vcf} \
      -resource axiomPoly,known=false,training=true,truth=false,prior=10:${axiomPoly_resource_vcf} \
      -resource dbsnp,known=true,training=false,truth=false,prior=2:${dbsnp_resource_vcf}
  }
  runtime {
    memory: "26 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
  }
  output {
    File recalibration = "${recalibration_filename}"
    File recalibration_index = "${recalibration_filename}.idx"
    File tranches = "${tranches_filename}"
  }
}

task SNPsVariantRecalibratorCreateModel {
  String recalibration_filename
  String tranches_filename
  Int downsampleFactor
  String model_report_filename

  Array[String] recalibration_tranche_values
  Array[String] recalibration_annotation_values

  File sites_only_variant_filtered_vcf
  File sites_only_variant_filtered_vcf_index

  File hapmap_resource_vcf
  File omni_resource_vcf
  File one_thousand_genomes_resource_vcf
  File dbsnp_resource_vcf
  File hapmap_resource_vcf_index
  File omni_resource_vcf_index
  File one_thousand_genomes_resource_vcf_index
  File dbsnp_resource_vcf_index

  Int disk_size

  command {

    /usr/gitc/gatk-launch --javaOptions "-Xmx100g -Xms100g" \
      VariantRecalibrator \
      -V ${sites_only_variant_filtered_vcf} \
      -O ${recalibration_filename} \
      --tranchesFile ${tranches_filename} \
      -allPoly \
      -tranche ${sep=' -tranche ' recalibration_tranche_values} \
      -an ${sep=' -an ' recalibration_annotation_values} \
      -mode SNP \
      -sampleEvery ${downsampleFactor} \
      --output_model ${model_report_filename} \
      --maxGaussians 6 \
      -resource hapmap,known=false,training=true,truth=true,prior=15:${hapmap_resource_vcf} \
      -resource omni,known=false,training=true,truth=true,prior=12:${omni_resource_vcf} \
      -resource 1000G,known=false,training=true,truth=false,prior=10:${one_thousand_genomes_resource_vcf} \
      -resource dbsnp,known=true,training=false,truth=false,prior=7:${dbsnp_resource_vcf}
  }
  runtime {
    memory: "104 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
  }
  output {
    File model_report = "${model_report_filename}"
  }
}

task SNPsVariantRecalibratorScattered {
  String recalibration_filename
  String tranches_filename
  File model_report

  Array[String] recalibration_tranche_values
  Array[String] recalibration_annotation_values

  File sites_only_variant_filtered_vcf
  File sites_only_variant_filtered_vcf_index

  File hapmap_resource_vcf
  File omni_resource_vcf
  File one_thousand_genomes_resource_vcf
  File dbsnp_resource_vcf
  File hapmap_resource_vcf_index
  File omni_resource_vcf_index
  File one_thousand_genomes_resource_vcf_index
  File dbsnp_resource_vcf_index

  Int disk_size

  command {

    /usr/gitc/gatk-launch --javaOptions "-Xmx3g -Xms3g" \
      VariantRecalibrator \
      -V ${sites_only_variant_filtered_vcf} \
      -O ${recalibration_filename} \
      --tranchesFile ${tranches_filename} \
      -allPoly \
      -tranche ${sep=' -tranche ' recalibration_tranche_values} \
      -an ${sep=' -an ' recalibration_annotation_values} \
      -mode SNP \
      --input_model ${model_report} \
      -scatterTranches \
      --maxGaussians 6 \
      -resource hapmap,known=false,training=true,truth=true,prior=15:${hapmap_resource_vcf} \
      -resource omni,known=false,training=true,truth=true,prior=12:${omni_resource_vcf} \
      -resource 1000G,known=false,training=true,truth=false,prior=10:${one_thousand_genomes_resource_vcf} \
      -resource dbsnp,known=true,training=false,truth=false,prior=7:${dbsnp_resource_vcf}
  }
  runtime {
    memory: "3.5 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
  }
  output {
    File recalibration = "${recalibration_filename}"
    File recalibration_index = "${recalibration_filename}.idx"
    File tranches = "${tranches_filename}"
  }
}

task GatherTranches {
  File input_fofn
  String output_filename

  Int disk_size

  command <<<
    set -e
    set -o pipefail

    # this is here to deal with the JES bug where commands may be run twice
    rm -rf tranches

    mkdir tranches
    RETRY_LIMIT=5

    count=0
    until cat ${input_fofn} | /root/google-cloud-sdk/bin/gsutil -m cp -L cp.log -c -I tranches/; do
        sleep 1
        ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
        echo 'Could not copy all the tranches from the cloud' && exit 1
    fi

    cat ${input_fofn} | rev | cut -d '/' -f 1 | rev | awk '{print "tranches/" $1}' > inputs.list

      /usr/gitc/gatk-launch --javaOptions "-Xmx6g -Xms6g" \
      GatherTranches \
      --input inputs.list \
      --output ${output_filename}
  >>>
  runtime {
    memory: "7 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
  }
  output {
    File tranches = "${output_filename}"
  }
}

task ApplyRecalibration {
  String recalibrated_vcf_filename
  File input_vcf
  File input_vcf_index
  File indels_recalibration
  File indels_recalibration_index
  File indels_tranches
  File snps_recalibration
  File snps_recalibration_index
  File snps_tranches

  Float indel_filter_level
  Float snp_filter_level

  Int disk_size

  command {
    set -e

    /usr/gitc/gatk-launch --javaOptions "-Xmx5g -Xms5g" \
      ApplyVQSR \
      -O tmp.indel.recalibrated.vcf \
      -V ${input_vcf} \
      --recalFile ${indels_recalibration} \
      -tranchesFile ${indels_tranches} \
      -ts_filter_level ${indel_filter_level} \
      --createOutputVariantIndex true \
      -mode INDEL
      
    /usr/gitc/gatk-launch --javaOptions "-Xmx5g -Xms5g" \
      ApplyVQSR \
      -O ${recalibrated_vcf_filename} \
      -V tmp.indel.recalibrated.vcf \
      --recalFile ${snps_recalibration} \
      -tranchesFile ${snps_tranches} \
      -ts_filter_level ${snp_filter_level} \
      --createOutputVariantIndex true \
      -mode SNP
  }
  runtime {
    memory: "7 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
  }
  output {
    File recalibrated_vcf = "${recalibrated_vcf_filename}"
    File recalibrated_vcf_index = "${recalibrated_vcf_filename}.tbi"
  }
}

task GatherVcfs {

  File input_vcfs_fofn
  File input_vcf_indexes_fofn

  String output_vcf_name

  Int disk_size

  command <<<
    set -e
    set -o pipefail

    # this is here to deal with the JES bug where commands may be run twice
    rm -rf vcfs

    mkdir vcfs
    RETRY_LIMIT=5

    count=0
    until cat ${input_vcfs_fofn} | /root/google-cloud-sdk/bin/gsutil -m cp -L cp.log -c -I vcfs/; do
        sleep 1
        ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
        echo 'Could not copy all the vcfs from the cloud' && exit 1
    fi

    count=0
    until cat ${input_vcf_indexes_fofn} | /root/google-cloud-sdk/bin/gsutil -m cp -L cp.log -c -I vcfs/; do
        sleep 1
        ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
        echo 'Could not copy all the indexes from the cloud' && exit 1
    fi

    cat ${input_vcfs_fofn} | rev | cut -d '/' -f 1 | rev | awk '{print "vcfs/" $1}' > inputs.list

    # the various flags to disable safety checks are not ideal but are
    # (temporarily) necessary until Intel fixes a bug in GenomicsDB
    # which outputs contigs out of order in the VCF header
    /usr/gitc/gatk-launch --javaOptions "-Xmx6g -Xms6g" \
    GatherVcfs \
    --ignoreSafetyChecks \
    --disableContigOrderingCheck \
    --gatherType BLOCK \
    --input inputs.list \
    --output ${output_vcf_name}

    tabix ${output_vcf_name}
  >>>
  runtime {
    memory: "7 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
  }
  output {
    File output_vcf = "${output_vcf_name}"
    File output_vcf_index = "${output_vcf_name}.tbi"
  }
}

task CollectVariantCallingMetrics {
  File input_vcf
  File input_vcf_index
  
  String metrics_filename_prefix
  File dbsnp_vcf
  File dbsnp_vcf_index
  File interval_list
  File ref_dict

  Int disk_size

  command {
    java -Xmx6g -Xms6g -jar /usr/gitc/picard.jar \
      CollectVariantCallingMetrics \
      INPUT=${input_vcf} \
      DBSNP=${dbsnp_vcf} \
      SEQUENCE_DICTIONARY=${ref_dict} \
      OUTPUT=${metrics_filename_prefix} \
      THREAD_COUNT=8 \
      TARGET_INTERVALS=${interval_list}
  }
  output {
    File detail_metrics_file = "${metrics_filename_prefix}.variant_calling_detail_metrics"
    File summary_metrics_file = "${metrics_filename_prefix}.variant_calling_summary_metrics"
  }
  runtime {
    memory: "7 GB"
    cpu: 2
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
  }
}

task GatherMetrics {

  File input_details_fofn
  File input_summaries_fofn

  String output_prefix

  Int disk_size

  command <<<
    set -e
    set -o pipefail

    # this is here to deal with the JES bug where commands may be run twice
    rm -rf metrics

    mkdir metrics
    RETRY_LIMIT=5

    count=0
    until cat ${input_details_fofn} | /root/google-cloud-sdk/bin/gsutil -m cp -L cp.log -c -I metrics/; do
        sleep 1
        ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
        echo 'Could not copy all the metrics from the cloud' && exit 1
    fi

    count=0
    until cat ${input_summaries_fofn} | /root/google-cloud-sdk/bin/gsutil -m cp -L cp.log -c -I metrics/; do
        sleep 1
        ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
        echo 'Could not copy all the metrics from the cloud' && exit 1
    fi

    INPUT=`cat ${input_details_fofn} | rev | cut -d '/' -f 1 | rev | sed s/.variant_calling_detail_metrics//g | awk '{printf("I=metrics/%s ", $1)}'`

    java -Xmx14g -Xms14g -jar /usr/gitc/picard.jar \
    AccumulateVariantCallingMetrics \
    $INPUT \
    O= ${output_prefix}
  >>>
  runtime {
    memory: "15 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
  }
  output {
    File detail_metrics_file = "${output_prefix}.variant_calling_detail_metrics"
    File summary_metrics_file = "${output_prefix}.variant_calling_summary_metrics"
  }
}
