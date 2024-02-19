nextflow.enable.dsl=1
/*
#########################################################################################
##                                                                                     ##
##     homopolymer.nf                                                                  ##
##     Analysis of homopolymers in a batch of sequences in a fasta file                ##
##                                                                                     ##
##     Gael A. Millot                                                                  ##
##     Bioinformatics and Biostatistics Hub                                            ##
##     Computational Biology Department                                                ##
##     Institut Pasteur Paris                                                          ##
##                                                                                     ##
#########################################################################################
*/


print("\n\nRESULT DIRECTORY: ${out_path}")
print("\n\nWARNING: PARAMETERS ALREADY INTERPRETED IN THE .config FILE:")
print("    system_exec: ${system_exec}")
print("    out_path: ${out_path_ini}")
if("${system_exec}" != "local"){
    print("    queue: ${queue}")
    print("    qos: ${qos}")
    print("    add_options: ${add_options}")
}
print("\n\n")


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
    title: 'Homopolymer Report'
    author: 'Gael Millot'
    date: '`r Sys.Date()`'
    output:
      html_document:
        toc: TRUE
        toc_float: TRUE
    ---

    \\n\\n<br /><br />\\n\\n
    " > report.rmd
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
    val min_length

    output:
    file "homopol_df.tsv" into homopol_df_ch

    script:
    """
    homopolymer.py "${seqString}" "${id}" "${min_length}" "homopol_df.tsv"
    """
}

// homopol_df_ch1.collectFile(name: "${file_name}_homopol_summary.tsv").subscribe{it -> it.copyTo("${out_path}/files/${file_name}_homopol_summary.tsv")} // concatenate all the homopol_df.tsv files in channel homopol_df_ch into a single file published into the indicated path. STRONG WARNING: copyTo(${out_path}/files/) works only if the files folder already exists. Ohterwise, the saved file becomes "files"
homopol_df_ch.collectFile(name: "${file_name}_homopol_summary.tsv")into{tsv_ch1 ; tsv_ch2}




process graph_stat {
    label 'r_ext'
    publishDir "${out_path}/figures", mode: 'copy', pattern: "{*.png}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    publishDir "${out_path}/reports", mode: 'copy', pattern: "{graph_stat_report.txt}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    publishDir "${out_path}/files", mode: 'copy', pattern: "{*.tsv,*.RData}", overwrite: false  // warning,: no space after the comma in the pattern
    cache 'deep'

    input:
    val file_name
    val cute_path
    val min_length
    file  tsv from tsv_ch1

    output:
    file "*.png" into fig_ch1
    file "*.tsv" into table_ch1
    file "graph_stat_report.txt"
    file "graph_stat.RData"
    file "report.rmd" into log_ch1

    script:
    """
    graph_stat.R "${tsv}" "${file_name}" "${cute_path}" "graph_stat_report.txt"

    echo -e "\\n\\n<br /><br />\\n\\n###  Results\\n\\n" > report.rmd
    echo -e "The minimal polymer length of ${min_length}, set in the min_lenght paremeter of the nextflow.config file, is for the homopol_summary.tsv file only (except the two last columns of this file).\\n\\nRandomisation of each sequence was performed 10,000 times without any constrain.\\nThen means were computed for each homopolymer length category." >> report.rmd
    echo -e "\\n\\n<br /><br />\\n\\n#### Dot plot\\n\\n<br />" >> report.rmd
    echo -e "Each dot is a value obtained for one sequence.<br />See the [plot_raw_values.tsv](./files/) file for the plot raw values.<br /><br />\\n\\n" >> report.rmd
    echo -e "
\\n\\n</center>\\n\\n
![Figure 1: Proportions of homopolymer lengths. See the [scatterplot_stat.tsv](./files/) file for values](./figures/plot_${file_name}.png){width=600}
\\n\\n</center>\\n\\n
\\n\\n</center>\\n\\n
![Figure 2: Proportions of homopolymer lengths (log10)](./figures/plot_${file_name}_log10.png){width=600}
\\n\\n</center>\\n\\n
    " >> report.rmd
    echo -e "\\n\\n<br /><br />\\n\\n#### Boxplot plot\\n\\n<br />" >> report.rmd
    echo -e "Each dot is a value obtained for one sequence.<br /><br />\\n\\n" >> report.rmd
    echo -e "
\\n\\n</center>\\n\\n
![Figure 3: Proportions of homopolymer lengths: diamond, mean; whiskers, 1.5 x Inter Quartile Range; horizontal bars, quartiles; number at the top, mean. See the [boxplot_stat.tsv](./files/) file for values](./figures/boxplot_${file_name}.png){width=600}
\\n\\n</center>\\n\\n
\\n\\n</center>\\n\\n
![Figure 4: Proportions of homopolymer lengths (log10)](./figures/boxplot_${file_name}_log10.png){width=600}
\\n\\n</center>\\n\\n
    " >> report.rmd
    echo -e "\n\\n<br /><br />\\n\\n#### T test \"Obs versus Theo\" for each homopolymer length (see also the [t_test.tsv](./files/) file)<br />\\n\\n" >> report.rmd
    echo -e "See also the [t_test.tsv](./files/) file.<br />\\n\\n" >> report.rmd
    echo "
\\`\\`\\`{r, echo = FALSE}
tempo <- read.table('./files/t_test.tsv', header = TRUE, colClasses = 'character', sep = '\\t', check.names = FALSE) ; 
kableExtra::kable_styling(knitr::kable(tempo, row.names = FALSE, digits = 2, caption = NULL, format='html'), c('striped', 'bordered', 'responsive', 'condensed'), font_size=10, full_width = FALSE, position = 'left')
\\`\\`\\`
    \n\n
    " >> report.rmd
    """
}



process backup {
    label 'bash' // see the withLabel: bash in the nextflow config file 
    publishDir "${out_path}/reports", mode: 'copy', pattern: "{*.config,*.log}", overwrite: false // since I am in mode copy, all the output files will be copied into the publishDir. See \\wsl$\Ubuntu-20.04\home\gael\work\aa\a0e9a739acae026fb205bc3fc21f9b  // warning,: no space after the comma in the pattern
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
| Project<br />(empty means no .git folder where the homopolymer.nf file is present) | \$(git -C ${projectDir} remote -v | head -n 1) | # works only if the main script run is located in a directory that has a .git folder, i.e., that is connected to a distant repo
| Git info<br />(empty means no .git folder where the homopolymer.nf file is present) | \$(git -C ${projectDir} describe --abbrev=10 --dirty --always --tags) | # idem. Provide the small commit number of the script and nextflow.config used in the execution
| Cmd line | ${workflow.commandLine} |
| execution mode | ${system_exec} |" >> report.rmd 

    if [[ ! -z \$modules ]] ; then
        echo "| loaded modules (according to specification by the user thanks to the --modules argument of homopolymer.nf) | ${modules} |" >> report.rmd
    fi
    
    echo "| Manifest's pipeline version | ${workflow.manifest.version} |
| result path | ${out_path} |
| nextflow version | ${nextflow.version} |
    " >> report.rmd

    echo -e "\\n\\n<br /><br />\\n\\n#### Implicit variables\\n\\n
| Name | Description | Value | 
| :-- | :-- | :-- |
| launchDir | Directory where the workflow is run | ${launchDir} |
| nprojectDir | Directory where the homopolymer.nf script is located | ${projectDir} |
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
    publishDir "${out_path}/reports", mode: 'copy', pattern: "{*.txt,*.rmd}", overwrite: true // warning,: no space after the comma in the pattern
    cache 'false'

    input:
    val cute_path
    file report from log_ch0.concat(log_ch1, log_ch2, log_ch3).collectFile(name: 'report.rmd', sort: false)
    file png1 from fig_ch1.collect() // warning: several files
    file table from tsv_ch2
    file table2 from table_ch1.collect() // warning: several files

    output:
    file "report.html"
    file "print_report.txt"
    file "report.rmd"

    script:
    """
    cp ${report} report_file.rmd # this is to get hard files, not symlinks, for knitting
    mkdir figures
    mkdir files
    mkdir reports
    cp ${png1} ./figures/ # this is to get hard files, not symlinks, for knitting
    cp ${table} ${table2} ./files/ # this is to get hard files, not symlinks, for knitting
    echo "" > ./reports/nf_dag.png # trick to delude the knitting during the print report
    print_report.R "${cute_path}" "report_file.rmd" "print_report.txt"
    """
}


//////// end Processes
