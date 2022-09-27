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

sample_label = params.cli_attribute
outdir       = file(params.outdir)

log.info "Sample attribute will be used: $sample_label \n"

if(!outdir.isDirectory()){
    outdir = outdir.mkdirs()
    println outdir ? "Create folder: $outdir!" : "Cannot create directory: $outdir!"
}

process pre_process {
    container "proteomics/cosmo:latest"
    publishDir "${outdir}/", mode: "copy", overwrite: true

    input:
    file d1_file
    file d2_file
    file sample_file

    output:
    file "data_use/${d1_file.name}" into d1_file_use_1,d1_file_use_2
    file "data_use/${d2_file.name}" into d2_file_use_1,d2_file_use_2
    file "data_use/${sample_file.name}" into sample_file_use_1,sample_file_use_2


    script:
    """
    #!/usr/bin/env Rscript
    source("/opt/cosmo/tools.R")
    d1_file <- "${d1_file}"
    d2_file <- "${d2_file}"
    sample_file <- "${sample_file}"
    outdir <- "data_use"
    dir.create(outdir)
    format_input_data(d1_file, d2_file, sample_file, out_dir = outdir)

    """

}


process run_method_1 {
    container "proteomics/cosmo:latest"
    publishDir "${outdir}/method1_folder/", mode: "copy", overwrite: true

    input:
    file d1_file_use_1
    file d2_file_use_1
    file sample_file_use_1

    output:
    file "method1_folder" into method1_out_folder

    script:
    """
    #!/usr/bin/env Rscript
    source("/opt/cosmo/method1_function.R")
    d1_file <- "${d1_file_use_1}"
    d2_file <- "${d2_file_use_1}"
    sample_file <- "${sample_file_use_1}"
    gene_file <- "/opt/cosmo/genes.tsv"
    outdir <- "method1_folder"
    clinical_attributes <- unlist(strsplit(x="${sample_label}",split=","))
    run_2b(d1_file, d2_file, sample_file, gene_file, out_dir=outdir, clinical_attributes=clinical_attributes)
    """
}


process run_method_2 {
    container "proteomics/cosmo:latest"
    publishDir "${outdir}/method2_folder/", mode: "copy", overwrite: true

    input:
    file d1_file_use_2
    file d2_file_use_2
    file sample_file_use_2

    output:
    file "method2_folder" into method2_out_folder

    script:
    """
    python /opt/cosmo/method2_function.py \
        -d1 ${d1_file_use_2} \
        -d2 ${d2_file_use_2} \
        -s ${sample_file_use_2} \
        -l ${sample_label} \
        -o method2_folder

    """
}

process combine_methods {
    container "proteomics/cosmo:latest"
    publishDir "${outdir}/final_res_folder/", mode: "copy", overwrite: true

    input:
    file method1_out_folder
    file method2_out_folder
    file sample_file

    output:
    file "cosmo*" into final_res_folder

    script:
    """
    #!/usr/bin/env Rscript
    source("/opt/cosmo/method1_function.R")
    source("/opt/cosmo/combine_methods.R")
    method1_folder <- "${method1_out_folder}"
    method2_folder <- "${method2_out_folder}"
    sample_annotation_file <- "${sample_file}"
    clinical_attributes <- unlist(strsplit(x="${sample_label}",split=","))
    combine_methods(method1_folder, method2_folder, 
                    sample_annotation_file,
                    clinical_attributes = clinical_attributes, 
                    out_dir = "./", prefix = "cosmo")
    """
}