#!/usr/bin/env nextflow

params.help        = false 
params.pro_file    = "-"
params.rna_file    = "-"
params.cli_file    = ""
params.cli_attribute = "-"
params.cpu         = 0
params.out_dir     = "./output"


/* Prints help when asked for and exits */

def helpMessage() {
    log.info"""
    =========================================
    COSMO => COrrection of Sample Mislabeling by Omics
    =========================================
    Usage:
    nextflow run cosmo.nf
    Arguments:
      --pro_file              Protein expression data at gene level.
      --rna_file              RNA expressio data at gene level.
      --cli_file              Sample annotation data.
	  --cli_attribute         Sample attribute(s) for prediction. Multiple attributes 
	                          must be separated by ",".
      --out_dir               Output folder, default is "./output".
      --cpu                   The number of CPUs.
      --help                  Print help message.
    """.stripIndent()
}


// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}



pro_file    = file(params.pro_file)
rna_file    = file(params.rna_file)
sample_file = file(params.cli_file)
sample_label= params.cli_attribute
out_dir     = file(params.out_dir)
cpus        = params.cpu

println "sample attribute will be used: $sample_label \n"

if(!out_dir.isDirectory()){
    out_dir_result = out_dir.mkdirs()
    println out_dir_result ? "Create folder: $out_dir!" : "Cannot create directory: $myDir!"
}


process run_method_1 {

    tag "run_method_1"

    echo true

    container "cosmo:latest"

    publishDir "${out_dir}/method1_folder/", mode: "copy", overwrite: true

    input:
    file pro_file
    file rna_file
    file sample_file

    output:
    file "method1_folder" into method1_out_folder

    script:
    """
    #!/usr/bin/env /usr/local/bin/Rscript
    source("${baseDir}/bin/SoonJye_function.R")
    pro_file <- "${pro_file}"
    rna_file <- "${rna_file}"
    sample_file <- "${sample_file}"
    out_dir <- "method1_folder"
    run_2b(pro_file, rna_file, sample_file, out_dir=out_dir)

    """
}


process run_method_2 {

    tag "run_method_2"

    echo true

    container "cosmo:latest"

    publishDir "${out_dir}/method2_folder/", mode: "copy", overwrite: true

    input:
    file pro_file
    file rna_file
    file sample_file

    output:
    file "method2_folder" into method2_out_folder

    script:
    """
    python ${baseDir}/bin/sentieon.py \
        -pro ${pro_file} \
        -rna ${rna_file} \
        -s ${sample_file} \
        -l ${sample_label} \
        -o method2_folder

    """
    
}

process combine_methods {

    tag "combine_methods"

    echo true

    container "cosmo:latest"

    publishDir "${out_dir}/${final_res_folder}/", mode: "copy", overwrite: true

    input:
    file method1_out_folder
    file method2_out_folder
    file sample_file

    output:
    file "*cosmo*" into final_res_folder

    script:
    """
    #!/usr/bin/env /usr/local/bin/Rscript
    source("${baseDir}/bin/SoonJye_function.R")
    source("${baseDir}/bin/combine_methods.R")
    method1_folder <- "${method1_out_folder}"
    method2_folder <- "${method2_out_folder}"
    sample_annotation_file <- "${sample_file}"
    clinical_attributes <- "${sample_label}"
    combine_methods(method1_folder, method2_folder, 
                    sample_annotation_file,
                    clinical_attributes = clinical_attributes, 
                    out_dir = "./", prefix = "cosmo")
    """
    
}



