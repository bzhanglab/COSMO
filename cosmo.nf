#!/usr/bin/env nextflow

params.help        = false 
params.pro_file    = "-"
params.rna_file    = "-"
params.sample_file = ""
params.method_id   = 1
params.task_id     = "2b"
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
      --sample_file           Sample annotation data.
      --method_id             1:SoonJye, 2:Sentieon. Default is 1.
      --task_id               The task ID, 2b or 2c, default is 2b.
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
sample_file = file(params.sample_file)
method_id   = params.method_id
task_id     = params.task_id
out_dir     = file(params.out_dir)
cpus        = params.cpu



if(!out_dir.isDirectory()){
    out_dir_result = out_dir.mkdirs()
    println out_dir_result ? "Create folder: $out_dir!" : "Cannot create directory: $myDir!"
}


process main_process {
    tag "test"

    echo true

    container "cosmo:latest"

    publishDir "${out_dir}/", mode: "copy", overwrite: true

    input:
    file pro_file
    file rna_file
    file sample_file

    output:
    file "out_dir" into final_res_folder

    script:
    if (method_id == 1 && task_id == "2b"){
        println "Use method ${method_id} for task ${task_id}"
        """
        Rscript ${baseDir}/bin/SoonJye_2b.R ${pro_file} ${rna_file} ${sample_file} out_dir
        """
    
    }else if(method_id == 1 && task_id == "2c") {
        println "Use method ${method_id} for task ${task_id}"
        """
        Rscript ${baseDir}/bin/SoonJye_2c.R ${pro_file} ${rna_file} ${sample_file} out_dir
        """
    } else if(method_id == 2 && task_id == "2b"){
        println "Use method ${method_id} for task ${task_id}"
        """
        python ${baseDir}/bin/sentieon_2b.py \
            -pro ${pro_file} \
            -rna ${rna_file} \
            -s ${sample_file} \
            -o out_dir

        """
    } else if(method_id == 2 && task_id == "2c"){
        println "Use method ${method_id} for task ${task_id}"
        """
        python ${baseDir}/bin/sentieon_2c.py \
            -pro ${pro_file} \
            -rna ${rna_file} \
            -s ${sample_file} \
            -o out_dir
        """
    } else {
        println "Invalid method ${method_id} or task ${task_id}"
    }
}


