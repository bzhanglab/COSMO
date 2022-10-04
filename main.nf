#!/usr/bin/env nextflow

/* Prints help when asked for and exits */
def helpMessage() {
    log.info"""
    =========================================
    COSMO => COrrection of Sample Mislabeling by Omics
    =========================================
    Usage:
    nextflow run cosmo.nf
    Arguments:
      --d1_file               Dataset with quantification data at gene level.
      --d2_file               Dataset with quantification data at gene level.
      --cli_file              Sample annotation data.
      --cli_attribute         Sample attribute(s) for prediction. Multiple attributes 
                              must be separated by ",".
      --outdir                Output folder, default is "results".
      --help                  Print help message.
    """.stripIndent()
}

// Show help emssage
if (params.help) {
    helpMessage()
    exit 0
}

checkPathParamList = [params.d1_file, params.d2_file, params.cli_file]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

if (params.d1_file) { d1_file     = file(params.d1_file)  } else { exit 1, 'No file specified with --d1_file'  }
if (params.d1_file) { d2_file     = file(params.d2_file)  } else { exit 1, 'No file specified with --d2_file'  }
if (params.d1_file) { sample_file = file(params.cli_file) } else { exit 1, 'No file specified with --cli_file' }

outdir       = file(params.outdir)

log.info "Sample attribute will be used: $params.cli_attribute \n"

if(!outdir.isDirectory()){
    outdir = outdir.mkdirs()
    println outdir ? "Create folder: $outdir!" : "Cannot create directory: $outdir!"
}

process PREPROCESS {
    label 'process_low'

    input:
    path(d1_file)
    path(d2_file)
    path(sample_file)

    output:
    tuple path("out/${d1_file.name}"), path("out/${d2_file.name}"), path("out/${sample_file.name}")

    script:
    "format_input_data --d1 $d1_file --d2 $d2_file --samples $sample_file --out out"
}


process METHOD1 {
    label 'process_medium'

    input:
    tuple path(d1_file), path(d2_file), path(samplefile)
    path(gene_tsv)

    output:
    path("out/*")

    script:
    """
    cosmo one \\
      --d1 $d1_file \\
      --d2 $d2_file \\
      --samples $samplefile \\
      --out out \\
      --genes $gene_tsv \\
      --attributes ${params.cli_attribute} \\
      --cpus ${task.cpus}
    """
}


process METHOD2 {
    label 'process_medium'

    input:
    tuple path(d1_file), path(d2_file), path(samplefile)

    output:
    path("out/*")

    script:
    """
    python /opt/cosmo/method2_function.py \\
        -d1 ${d1_file} \\
        -d2 ${d2_file} \\
        -s ${samplefile} \\
        -l ${params.cli_attribute} \\
        -o out
    """
}

process COMBINE {
    label 'process_medium'

    input:
    path(method1_out_folder)
    path(method2_out_folder)
    path(sample_file)

    output:
    path("cosmo*")

    script:
    """
    cosmo combine \\
      --method-one-out $method1_out_folder \\
      --method-two-out $method2_out_folder \\
      --samples $sample_file \\
      --attributes ${params.cli_attribute} \\
      --prefix cosmo \\
      --cpus ${task.cpus} \\
      --out .
    """
}

workflow {
    genes = Channel.fromPath(params.genes)

    PREPROCESS(d1_file, d2_file, sample_file)

    METHOD1(PREPROCESS.out, genes.first())
    METHOD2(PREPROCESS.out)

    COMBINE(METHOD1.out, METHOD2.out, sample_file)
}