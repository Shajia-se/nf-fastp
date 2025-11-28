#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def fastp_output = params.fastp_output ?: "fastp_output"

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
      -i ${reads[0]} \\
      -I ${reads[1]} \\
      -o ${sample_id}_R1.fastp.trimmed.fastq.gz \\
      -O ${sample_id}_R2.fastp.trimmed.fastq.gz \\
      -j ${sample_id}.fastp.json \\
      -h ${sample_id}.fastp.html \\
      -w ${task.cpus}
  """
}

workflow {

  def data = Channel.fromFilePairs(
    "${params.fastqc_raw_data}/*_{R1,R2}_001.fastq.gz",
    flat: true
  )

  data = data.filter { sample_id, read1, read2 ->
      def report = new File("${params.project_folder}/${fastp_output}/${sample_id}.fastp.html")
      ! report.exists()
  }

  fastp(data)
}
