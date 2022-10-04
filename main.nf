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
    input:
    tuple path(d1_file), path(d2_file), path(samplefile)
    path(gene_tsv)

    output:
    path("results_method1")

    script:
    """
    method1 \
      --d1 $d1_file \
      --d2 $d2_file \
      --samples $samplefile \
      --out results_method1 \
      --genes $gene_tsv \
      --attributes ${params.cli_attribute}
    """
}


process METHOD2 {
    input:
    tuple path(d1_file), path(d2_file), path(samplefile)

    output:
    path("method2_folder")

    script:
    """
    python /opt/cosmo/method2_function.py \
        -d1 ${d1_file} \
        -d2 ${d2_file} \
        -s ${samplefile} \
        -l ${params.cli_attribute} \
        -o method2_folder
    """
}

process COMBINE {
    input:
    path(method1_out_folder)
    path(method2_out_folder)
    path(sample_file)

    output:
    path("cosmo*")

    script:
    """
    #!/usr/bin/env Rscript
    source("/opt/cosmo/method1_function.R")
    source("/opt/cosmo/combine_methods.R")
    method1_folder <- "${method1_out_folder}"
    method2_folder <- "${method2_out_folder}"
    sample_annotation_file <- "${sample_file}"
    clinical_attributes <- unlist(strsplit(x="${params.cli_attribute}",split=","))
    combine_methods(method1_folder, method2_folder, 
                    sample_annotation_file,
                    clinical_attributes = clinical_attributes, 
                    out_dir = "./", prefix = "cosmo")
    """
}

workflow {
    genes = Channel.fromPath(params.genes)

    PREPROCESS(d1_file, d2_file, sample_file)
    METHOD1(PREPROCESS.out, genes)
    METHOD2(PREPROCESS.out)

    COMBINE( METHOD1.out, METHOD2.out, sample_file)
}