#!/usr/bin/env nextflow

params.help          = false 
params.d1_file       = "-"
params.d2_file       = "-"
params.d1_type       = "-"
params.d2_type       = "-"
params.cli_file      = "-"
params.cli_attribute = "-"
params.threads       = 4
params.out_dir       = "./output"


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
      --d1_type               The type for dataset d1. This is used to label the dataset in the output files.
      --d2_type               The type for dataset d2. This is used to label the dataset in the output files.
      --cli_file              Sample annotation data.
      --cli_attribute         Sample attribute(s) for prediction. Multiple attributes 
                              must be separated by ",".
      --out_dir               Output folder, default is "./output".
      --threads               The number of threads.
      --help                  Print help message.
    """.stripIndent()
}


// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}



d1_file     = file(params.d1_file)
d2_file     = file(params.d2_file)
d1_type     = params.d1_type
d2_type     = file(params.d2_type)
sample_file = file(params.cli_file)
sample_label= params.cli_attribute
out_dir     = file(params.out_dir)
threads     = params.threads



if("${d1_file}" == "-" || !d1_file.exists()){
    if("${d1_file}" =~ /-$/){
        println "There is no file provided to --d1_file!"
        helpMessage()
        exit 0
    }else{    
        exit 1, "\n${d1_file} does not exist!\n"
    }
}

if("${d2_file}" == "-" || !d2_file.exists()){
    if("${d2_file}" =~ /-$/){
        println "There is no file provided to --d2_file!"
        helpMessage()
        exit 0
    }else{    
        exit 1, "\n${d2_file} does not exist!\n"
    }
}

if("${sample_file}" == "-" || !sample_file.exists()){
    if("${sample_file}" =~ /-$/){
        println "There is no file provided to --sample_file!"
        helpMessage()
        exit 0
    }else{    
        exit 1, "\n${sample_file} does not exist!\n"
    }
}

if(params.threads <= 0 || threads <= 0){
    threads = 4
}

println "sample attribute will be used: $sample_label \n"

if(!out_dir.isDirectory()){
    out_dir_result = out_dir.mkdirs()
    println out_dir_result ? "Create folder: $out_dir!" : "Cannot create directory: $myDir!"
}




process pre_process {
    tag "preprocessing"
    
    echo true
    
    container "proteomics/cosmo:latest"
    
    publishDir "${out_dir}/", mode: "copy", overwrite: true

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
    #!/usr/bin/env /usr/local/bin/Rscript
    source("/opt/cosmo/tools.R")
    d1_file <- "${d1_file}"
    d2_file <- "${d2_file}"
    sample_file <- "${sample_file}"
    out_dir <- "data_use"
    dir.create(out_dir)
    format_input_data(d1_file, d2_file, sample_file, out_dir = out_dir)

    """

}


process run_method_1 {

    tag "run_method_1"

    //cpus ${threads}

    echo true

    container "proteomics/cosmo:latest"

    publishDir "${out_dir}/method1_folder/", mode: "copy", overwrite: true

    input:
    file d1_file_use_1
    file d2_file_use_1
    file sample_file_use_1

    output:
    file "method1_folder" into method1_out_folder

    script:
    """
    #!/usr/bin/env /usr/local/bin/Rscript
    source("/opt/cosmo/method1_function.R")
    d1_file <- "${d1_file_use_1}"
    d2_file <- "${d2_file_use_1}"
    sample_file <- "${sample_file_use_1}"
    gene_file <- "${baseDir}/bin/genes.tsv"
    out_dir <- "method1_folder"
    clinical_attributes <- unlist(strsplit(x="${sample_label}",split=","))
    run_2b(d1_file, d2_file, sample_file, gene_file, out_dir=out_dir, clinical_attributes=clinical_attributes)

    """
}


process run_method_2 {

    tag "run_method_2"

    //cpus ${threads}

    echo true

    container "proteomics/cosmo:latest"

    publishDir "${out_dir}/method2_folder/", mode: "copy", overwrite: true

    input:
    file d1_file_use_2
    file d2_file_use_2
    file sample_file_use_2

    output:
    file "method2_folder" into method2_out_folder

    script:
    """
    python /opt/cosmo/method2_function.py \
        -pro ${d1_file_use_2} \
        -rna ${d2_file_use_2} \
        -s ${sample_file_use_2} \
        -l ${sample_label} \
        -o method2_folder

    """
    
}

process combine_methods {

    tag "combine_methods"

    //cpus ${threads}

    echo true

    container "proteomics/cosmo:latest"

    publishDir "${out_dir}/final_res_folder/", mode: "copy", overwrite: true

    input:
    file method1_out_folder
    file method2_out_folder
    file sample_file

    output:
    file "cosmo*" into final_res_folder

    script:
    """
    #!/usr/bin/env /usr/local/bin/Rscript
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



