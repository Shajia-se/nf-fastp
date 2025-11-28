#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def fastp_output = params.fastp_output ?: "fastp_output"

process fastp {
  tag "${sample_id}"
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    tuple val(sample_id), val(reads) 

  publishDir "${params.project_folder}/${fastp_output}", mode: 'copy'

  output:
    path "${sample_id}_R1.fastp.trimmed.fastq.gz"
    path("${sample_id}_R2.fastp.trimmed.fastq.gz"), optional: true
    path "${sample_id}.fastp.json"
    path "${sample_id}.fastp.html"

  script:
    def is_se = (reads instanceof Path)
    def r1    = is_se ? reads : reads[0]
    def r2    = is_se ? null  : reads[1]

    if (is_se) {
      """
      fastp \\
        -i ${r1} \\
        -o ${sample_id}_R1.fastp.trimmed.fastq.gz \\
        -j ${sample_id}.fastp.json \\
        -h ${sample_id}.fastp.html \\
        -w ${task.cpus}
      """
    }
    else {
      """
      fastp \\
        -i ${r1} \\
        -I ${r2} \\
        -o ${sample_id}_R1.fastp.trimmed.fastq.gz \\
        -O ${sample_id}_R2.fastp.trimmed.fastq.gz \\
        -j ${sample_id}.fastp.json \\
        -h ${sample_id}.fastp.html \\
        -w ${task.cpus}
      """
    }
}

workflow {

  def outdir = "${params.project_folder}/${fastp_output}"

  def ch = Channel
    .fromFilePairs("${params.fastqc_raw_data}/*_{R1,R2}.fastq.gz", flat: true)
    .ifEmpty {
      Channel
        .fromPath("${params.fastqc_raw_data}/*fastq.gz")
        .map { f ->
          tuple(f.simpleName.replaceAll(/.fastq.gz$/, ''), f)
        }
    }

  ch = ch.filter { sample_id, reads ->
    ! file("${outdir}/${sample_id}.fastp.html").exists()
  }

  fastp(ch)
}
