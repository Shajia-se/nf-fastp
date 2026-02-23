#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def fastp_output = params.fastp_output ?: "fastp_output"
def raw_data_dir = params.fastp_raw_data ?: params.fastqc_raw_data

process fastp {
  stageInMode 'symlink'
  stageOutMode 'move'
  tag "${sample_id}"

  input:
    tuple val(sample_id), path(read1), path(read2)

  publishDir "${params.project_folder}/${fastp_output}", mode: 'copy'

  output:
    path "${sample_id}_R1.fastp.trimmed.fastq.gz"
    path "${sample_id}_R2.fastp.trimmed.fastq.gz"
    path "${sample_id}.fastp.json"
    path "${sample_id}.fastp.html"

  script:
  """
  fastp \\
      -i ${read1} \\
      -I ${read2} \\
      -o ${sample_id}_R1.fastp.trimmed.fastq.gz \\
      -O ${sample_id}_R2.fastp.trimmed.fastq.gz \\
      -j ${sample_id}.fastp.json \\
      -h ${sample_id}.fastp.html \\
      -w ${task.cpus}
  """
}

workflow {

  def data
  if (params.samples_master) {
    data = Channel
      .fromPath(params.samples_master, checkIfExists: true)
      .splitCsv(header: true)
      .filter { row ->
        def enabled = row.enabled?.toString()?.trim()?.toLowerCase()
        enabled == null || enabled == '' || enabled == 'true'
      }
      .map { row ->
        def sample = row.sample_id?.toString()?.trim()
        def r1 = row.fastq_r1?.toString()?.trim()
        def r2 = row.fastq_r2?.toString()?.trim()

        assert sample : "samples_master row is missing sample_id"
        assert r1 : "samples_master row '${sample}' is missing fastq_r1"
        assert r2 : "samples_master row '${sample}' is missing fastq_r2 (nf-fastp expects paired-end data)"

        def f1 = file(r1)
        def f2 = file(r2)
        assert f1.exists() : "FASTQ not found for sample '${sample}': ${f1}"
        assert f2.exists() : "FASTQ not found for sample '${sample}': ${f2}"

        tuple(sample, f1, f2)
      }
      .ifEmpty { exit 1, "ERROR: No enabled paired-end rows found in samples_master: ${params.samples_master}" }
  } else {
    def pattern = params.fastp_pattern ?: "*_R{1,2}_001.fastq.gz"
    data = Channel.fromFilePairs(
      "${raw_data_dir}/${pattern}",
      flat: true,
      checkIfExists: true
    ).ifEmpty { exit 1, "ERROR: No FASTQ pairs found for pattern: ${raw_data_dir}/${pattern}" }
  }

  data = data.filter { sample_id, read1, read2 ->
      def report = file("${params.project_folder}/${fastp_output}/${sample_id}.fastp.html")
      ! report.exists()
  }

  fastp(data)
}
