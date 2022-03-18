/*
#########################################################################################
##                                                                                     ##
##     main.nf                                                                         ##
##     Analysis of homopolymers in a batch of sequences in a fasta file                ##
##                                                                                     ##
##     Gael A. Millot                                                                  ##
##     Bioinformatics and Biostatistics Hub                                            ##
##     Computational Biology Department                                                ##
##     Institut Pasteur Paris                                                          ##
##                                                                                     ##
#########################################################################################
*/


//////// Arguments of nextflow run

params.modules = ""

//////// end Arguments of nextflow run


//////// Variables

// from the nextflow.config file
config_file = file("${projectDir}/nextflow.config")
log_file = file("${launchDir}/.nextflow.log")
file_name = file("${in_path}/${fasta_file}").baseName



// end from the nextflow.config file

// from parameters
modules = params.modules // remove the dot -> can be used in bash scripts
// end from parameters

//////// end Variables


//////// Variables from config.file that need to be checked

fasta_ch_test = file("${in_path}/${fasta_file}") // to test if exist below

//////// end Variables from config.file that need to be checked


//////// Channels

fasta_ch = Channel.fromPath("${in_path}/${fasta_file}", checkIfExists: false) // I could use true, but I prefer to perform the check below, in order to have a more explicit error message

//////// end Channels

//////// Checks

if(system_exec == 'local' || system_exec == 'slurm' || system_exec == 'slurm_local'){
    def file_exists1 = fasta_ch_test.exists()
    if( ! file_exists1){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID in_path AND fasta_file PARAMETERS IN nextflow.config FILE: ${in_path}/${fasta_file}\n\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
}else{
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID system_exec PARAMETER IN nextflow.config FILE: ${system_exec}\nTHE ONLY POSSIBLE VALUES ARE local, slurm OR slurm_local\n\n========\n\n"
}

//////// end Checks



//////// Processes

process init {
    label 'bash' // see the withLabel: bash in the nextflow config file 
    cache 'false'

    output:
    file "report.rmd" into log_ch0

    script:
    """
    echo "---
    title: 'Insertion Sites Report'
    author: 'Gael Millot'
    date: '`r Sys.Date()`'
    output:
      html_document:
        toc: TRUE
        toc_float: TRUE
    ---

    \\n\\n<br /><br />\\n\\n
    " >> report.rmd
    """
}



process homopolymer {
    label 'python'
    cache 'true'

    input:
    set  id, seqString from fasta_ch.splitFasta(record:[header: true, seqString: true]) // warning: several data -> parall expected
    // splitFasta(record:[header: true, seqString: true]) split a fasta files in records (one per sequence), each record with two items: the first is the header and the second the seq string
    // see https://www.nextflow.io/docs/latest/operator.html#splitfasta
    // set takes the first and second items of the record

    output:
    file "homopol_report.tsv" into homopol_report_ch

    script:
    """
    homopolymer.py "${seqString}" "${id}" "homopol_report.tsv"
    """
}

homopol_report_ch.collectFile(name: "${file_name}.tsv").subscribe{it -> it.copyTo("${out_path}/files")} // concatenate all the homopol_report.tsv files in channel homopol_report_ch into a single file published into 





process backup {
    label 'bash' // see the withLabel: bash in the nextflow config file 
    publishDir "${out_path}/reports", mode: 'copy', pattern: "{*.config,*.log}", overwrite: false // since I am in mode copy, all the output files will be copied into the publishDir. See \\wsl$\Ubuntu-20.04\home\gael\work\aa\a0e9a739acae026fb205bc3fc21f9b
    cache 'false'

    input:
    file config_file
    file log_file

    output:
    file "${config_file}" // warning message if we use file config_file
    file "${log_file}" // warning message if we use file log_file
    file "report.rmd" into log_ch2

    script:
    """
    echo -e "\\n\\n<br /><br />\\n\\n###  Backup\\n\\n" > report.rmd
    echo -e "See the [reports](./reports) folder for all the details of the analysis, including the parameters used in the .config file" >> report.rmd
    echo -e "\\n\\nFull .nextflow.log is in: ${launchDir}<br />The one in the [reports](./reports) folder is not complete (miss the end)" >> report.rmd
    """
}


process workflowVersion { // create a file with the workflow version in out_path
    label 'bash' // see the withLabel: bash in the nextflow config file 
    cache 'false'

    output:
    file "report.rmd" into log_ch3

    script:
    """
    modules=$modules # this is just to deal with variable interpretation during the creation of the .command.sh file by nextflow. See also \$modules below
    echo -e "\\n\\n<br /><br />\\n\\n###  Workflow Version\\n\\n" > report.rmd
    echo -e "\\n\\n#### General\\n\\n
| Variable | Value |
| :-- | :-- |
| Project<br />(empty means no .git folder where the main.nf file is present) | \$(git -C ${projectDir} remote -v | head -n 1) | # works only if the main script run is located in a directory that has a .git folder, i.e., that is connected to a distant repo
| Git info<br />(empty means no .git folder where the main.nf file is present) | \$(git -C ${projectDir} describe --abbrev=10 --dirty --always --tags) | # idem. Provide the small commit number of the script and nextflow.config used in the execution
| Cmd line | ${workflow.commandLine} |
| execution mode | ${system_exec} |" >> report.rmd 

    if [[ ! -z \$modules ]] ; then
        echo "| loaded modules (according to specification by the user thanks to the --modules argument of main.nf) | ${modules} |" >> report.rmd
    fi
    
    echo "| Manifest's pipeline version | ${workflow.manifest.version} |
| result path | ${out_path} |
| nextflow version | ${nextflow.version} |
    " >> report.rmd

    echo -e "\\n\\n<br /><br />\\n\\n#### Implicit variables\\n\\n
| Name | Description | Value | 
| :-- | :-- | :-- |
| launchDir | Directory where the workflow is run | ${launchDir} |
| nprojectDir | Directory where the main.nf script is located | ${projectDir} |
| workDir | Directory where tasks temporary files are created | ${workDir} |
    " >> report.rmd

    echo -e "\\n\\n<br /><br />\\n\\n#### User variables\\n\\n
| Name | Description | Value | 
| :-- | :-- | :-- |
| out_path | output folder path | ${out_path} |
| in_path | input folder path | ${in_path} |
    " >> report.rmd

    echo -e "\\n\\n<br /><br />\\n\\n#### Workflow diagram\\n\\nSee the [nf_dag.png](./reports/nf_dag.png) file" >> report.rmd
    """
}
//${projectDir} nextflow variable
//${workflow.commandLine} nextflow variable
//${workflow.manifest.version} nextflow variable
//Note that variables like ${out_path} are interpreted in the script block


process print_report { // section 8.8 of the labbook 20200520
    label 'r_ext' // see the withLabel: bash in the nextflow config file 
    publishDir "${out_path}", mode: 'copy', pattern: "{*.html}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    publishDir "${out_path}/reports", mode: 'copy', pattern: "{*.txt,*.rmd}", overwrite: true // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    cache 'false'

    input:
    val cute_path
    file report from log_ch0.concat(log_ch2, log_ch3).collectFile(name: 'report.rmd', sort: false)
    //file png1 from fig_ch1

    output:
    file "report.html"
    file "print_report.txt"
    file "report.rmd"

    script:
    """
    cp ${report} report_file.rmd # this is to get hard files, not symlinks
    mkdir figures
    mkdir files
    mkdir reports
    echo "" > ./reports/nf_dag.png # trick to delude the knitting during the print report
    print_report.R "${cute_path}" "report_file.rmd" "print_report.txt"
    """
}


//////// end Processes
