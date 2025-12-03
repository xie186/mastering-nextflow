nextflow.enable.dsl=2

// If you use typed params, enable the strict parser when running:
// NXF_SYNTAX_PARSER=v2 nextflow run sra2fastq.nf ...
params {
  meta:   Path = "sample_meta.csv"
  outdir: Path = "fastq_test"   // we'll point -output-dir at this
}

process SRA2FASTQ {
  tag "${sample_id}"

  input:
  tuple val(sample_id), val(sra_id)

  output:
  // emit sample_id along with both FASTQs as a single pattern
  //[id: val(sample_id), read1: file("${sample_id}_1.fastq.gz"), read2: file("${sample_id}_2.fastq.gz")]
  tuple val(sample_id), path("${sample_id}_*.fastq.gz")

  script:
  """
  echo "Processing ${sample_id} with SRA ID ${sra_id} R1\n" |gzip - > ${sample_id}_1.fastq.gz
  echo "Processing ${sample_id} with SRA ID ${sra_id} R2\n" |gzip - > ${sample_id}_2.fastq.gz
  """
}

workflow {
  main:
  ch_fastqs = Channel
               .fromPath(params.meta)
               .splitCsv(header: false)
               .map { row -> tuple(row[0], row[1]) }
               | SRA2FASTQ

  ch_fastqs.view()
  println "Total samples processed: ${ch_fastqs.count().next()}"
  // println "Sample IDs: ${ch_fastqs.map{ it[0] }.join(', ')}"
  println "Meta data file used: ${params.meta.toString()}"
  println "Output directory: ${params.outdir.toString()}"
  publish:
  fastqs = ch_fastqs
  
}


// ---- workflow outputs block (v25.10+) ----
output {
  // Publish to <outputDir>/<sample_id>/<files>
  fastqs {
    mode 'copy'
    // Use a closure so we can include the sample_id and keep filenames
    path { fastq -> 
      fastq[1][0] >> "${params.outdir}/${fastq[0]}/" 
      fastq[1][1] >> "${params.outdir}/${fastq[0]}/"
    }
  }
}