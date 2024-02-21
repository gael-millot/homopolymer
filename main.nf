nextflow.enable.dsl=2
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





//////// Processes



process workflowParam { // create a file with the workflow parameters in out_path
    label 'bash'
    publishDir "${out_path}/reports", mode: 'copy', overwrite: false
    cache 'false'

    input: // does not work if modules if not declared
    val modules

    output:
    path "Run_info.txt"

    script:
    """
    echo "Project (empty means no .git folder where the main.nf file is present): " \$(git -C ${projectDir} remote -v | head -n 1) > Run_info.txt # works only if the main script run is located in a directory that has a .git folder, i.e., that is connected to a distant repo
    echo "Git info (empty means no .git folder where the main.nf file is present): " \$(git -C ${projectDir} describe --abbrev=10 --dirty --always --tags) >> Run_info.txt # idem. Provide the small commit number of the script and nextflow.config used in the execution
    echo "Cmd line: ${workflow.commandLine}" >> Run_info.txt
    echo "execution mode": ${system_exec} >> Run_info.txt
    modules=$modules # this is just to deal with variable interpretation during the creation of the .command.sh file by nextflow. See also \$modules below
    if [[ ! -z \$modules ]] ; then
        echo "loaded modules (according to specification by the user thanks to the --modules argument of main.nf): ${modules}" >> Run_info.txt
    fi
    echo "Manifest's pipeline version: ${workflow.manifest.version}" >> Run_info.txt
    echo "result path: ${out_path}" >> Run_info.txt
    echo "nextflow version: ${nextflow.version}" >> Run_info.txt
    echo -e "\\n\\nIMPLICIT VARIABLES:\\n\\nlaunchDir (directory where the workflow is run): ${launchDir}\\nprojectDir (directory where the main.nf script is located): ${projectDir}\\nworkDir (directory where tasks temporary files are created): ${workDir}" >> Run_info.txt
    echo -e "\\n\\nUSER VARIABLES:\\n\\nout_path: ${out_path}\\nsample_path: ${sample_path}" >> Run_info.txt
    """
}
//${projectDir} nextflow variable
//${workflow.commandLine} nextflow variable
//${workflow.manifest.version} nextflow variable
//Note that variables like ${out_path} are interpreted in the script block


process initial { // init name does not work !
    label 'r_ext'
    cache 'false'

    output:
    path "report.rmd", emit: log_ch0

    script:
    """
    echo "---
    title: 'Homopolymer Report'
    author: 'Gael Millot'
    date: \"\$(Rscript -e 'cat(format(Sys.Date(), "%Y-%m-%d"))')\"
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
    tuple val(id), val(seqString) // warning: several data -> parall expected
    val min_length

    output:
    path "homopol_df.tsv", emit: homopol_df_ch

    script:
    """
    homopolymer.py "${seqString}" "${id}" "${min_length}" "homopol_df.tsv"
    """
}


process graph_stat {
    label 'r_ext'
    publishDir "${out_path}/figures/png", mode: 'copy', pattern: "{*.png}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    publishDir "${out_path}/figures/svg", mode: 'copy', pattern: "{*.svg}", overwrite: false
    publishDir "${out_path}/reports", mode: 'copy', pattern: "{graph_stat_report.txt}", overwrite: false // https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
    publishDir "${out_path}/files", mode: 'copy', pattern: "{*.tsv,*.RData}", overwrite: false  // warning,: no space after the comma in the pattern
    cache 'deep'

    input:
    val file_name
    path cute_file
    val min_length
    path  tsv

    output:
    path "*.png", emit: fig_ch1
    path "*.svg"
    path "*.tsv", emit: table_ch1
    path "graph_stat_report.txt"
    path "graph_stat.RData"
    path "report.rmd", emit: log_ch1

    script:
    """
    graph_stat.R "${tsv}" "${file_name}" "${cute_file}" "graph_stat_report.txt"

    echo -e "\\n\\n<br /><br />\\n\\n###  Results\\n\\n" > report.rmd
    echo -e "The minimal polymer length of ${min_length}, set in the min_lenght paremeter of the nextflow.config file, is for the homopol_summary.tsv file only (except the two last columns of this file).\\n\\nRandomisation of each sequence was performed 10,000 times without any constrain.\\nThen means were computed for each homopolymer length category." >> report.rmd
    echo -e "\\n\\n<br /><br />\\n\\n#### Dot plot\\n\\n<br />" >> report.rmd
    echo -e "Each dot is a value obtained for one sequence.<br />See the [plot_raw_values.tsv](./files/) file for the plot raw values.<br /><br />\\n\\n" >> report.rmd
    echo -e "
\\n\\n</center>\\n\\n
![Figure 1: Proportions of homopolymer lengths. See the [scatterplot_stat.tsv](./files/) file for values](./figures/png/plot_${file_name}.png){width=600}
\\n\\n</center>\\n\\n
\\n\\n</center>\\n\\n
![Figure 2: Proportions of homopolymer lengths (log10)](./figures/png/plot_${file_name}_log10.png){width=600}
\\n\\n</center>\\n\\n
    " >> report.rmd
    echo -e "\\n\\n<br /><br />\\n\\n#### Boxplot plot\\n\\n<br />" >> report.rmd
    echo -e "Each dot is a value obtained for one sequence.<br /><br />\\n\\n" >> report.rmd
    echo -e "
\\n\\n</center>\\n\\n
![Figure 3: Proportions of homopolymer lengths: diamond, mean; whiskers, 1.5 x Inter Quartile Range; horizontal bars, quartiles; number at the top, mean. See the [boxplot_stat.tsv](./files/) file for values](./figures/png/boxplot_${file_name}.png){width=600}
\\n\\n</center>\\n\\n
\\n\\n</center>\\n\\n
![Figure 4: Proportions of homopolymer lengths (log10)](./figures/png/boxplot_${file_name}_log10.png){width=600}
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
    path config_file
    path log_file

    output:
    path "${config_file}" // warning message if we use file config_file
    path "${log_file}" // warning message if we use file log_file
    path "report.rmd", emit: log_ch2

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

    input: // does not work if modules if not declared
    val modules

    output:
    path "report.rmd", emit: log_ch3

    script:
    """
    modules=${modules} # this is just to deal with variable interpretation during the creation of the .command.sh file by nextflow. See also \$modules below
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
| in_path | input folder path | ${sample_path} |
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
    path cute_file
    path report
    path png1 // warning: several files
    path table
    path table2 // warning: several files

    output:
    path "report.html"
    path "print_report.txt"
    path "report.rmd"

    script:
    """
    cp ${report} report_file.rmd # this is to get hard files, not symlinks, for knitting
    mkdir -p figures/png
    mkdir files
    mkdir reports
    cp ${png1} ./figures/png/ # this is to get hard files, not symlinks, for knitting
    cp ${table} ${table2} ./files/ # this is to get hard files, not symlinks, for knitting
    echo "" > ./reports/nf_dag.png # trick to delude the knitting during the print report
    print_report.R "${cute_file}" "report_file.rmd" "print_report.txt"
    """
}


//////// end Processes



//////// Workflow


workflow {

    //////// Options of nextflow run

    print("\n\nINITIATION TIME: ${workflow.start}")

    //////// end Options of nextflow run


    //////// Options of nextflow run

    // --modules (it is just for the process workflowParam)
    params.modules = "" // if --module is used, this default value will be overridden
    // end --modules (it is just for the process workflowParam)

    //////// end Options of nextflow run


    //////// Variables


    config_file = workflow.configFiles[0] // better to use this than config_file = file("${projectDir}/nextflow.config") because the latter is not good if -c option of nextflow run is used // file() create a path object necessary o then create the file
    log_file = file("${launchDir}/.nextflow.log")

    // from parameters (options of the nexflow command line)
    modules = params.modules // remove the dot -> can be used in bash scripts
    // end from parameters (options of the nexflow command line)

    //////// end Variables


    //////// Checks

    if( ! (sample_path.class == String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN repertoire_profiler.config FILE:\n${sample_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (file(sample_path).exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN repertoire_profiler.config FILE (DOES NOT EXIST): ${sample_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    if( ! min_length.class == String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID thread_nb PARAMETER IN nextflow.config FILE:\n${min_length}\nMUST BE A SINGLE CHARACTER STRING OF AN INTEGER VALUE\n\n========\n\n"
    }else if( ! (min_length =~ /^[0123456789]+$/ || min_length == "NULL")){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID thread_nb PARAMETER IN nextflow.config FILE:\n${min_length}\nMUST BE A SINGLE CHARACTER STRING OF AN INTEGER VALUE\n\n========\n\n"
    }

    // below : those variable are already used in the config file. Thus, to late to check them. And not possible to check inside the config file
    // system_exec
    // out_ini
    print("\n\nRESULT DIRECTORY: ${out_path}")
    if("${system_exec}" != "local"){
        print("    queue: ${queue}")
        print("    qos: ${qos}")
        print("    add_options: ${add_options}")
    }
    print("\n\n")


    //////// end Checks


    //////// Variables from config.file that need to be modified


    //////// end Variables from config.file that need to be modified



    //////// Channels

    fasta_ch = Channel.fromPath(sample_path) 

    //////// end Channels


    //////// files import

    // in variable because a single file. If "NULL", will create a empty file, present in work folders, but that cannot be correctly linked. Thus, if the file has to be redirected into a channel inside a process, it will not work. Thus, in the first process using meta_file, I hard copy the NULL file if required (see below)
    cute_file = file(cute_path) // in variable because a single file
    sample_file = file(sample_path) // in variable because a single file, to get the name of the file

    //////// end files import



    //////// Main


    workflowParam(
        modules
    )

    initial()

    // fasta_ch.splitFasta(record:[header: true, seqString: true]).view()
    homopolymer(
        fasta_ch.splitFasta(record:[header: true, seqString: true]), // parallel
        min_length
    )
    // splitFasta(record:[header: true, seqString: true]) split a fasta files in records (one per sequence), each record with two items: the first is the header and the second the seq string
    // see https://www.nextflow.io/docs/latest/operator.html#splitfasta
    // set takes the first and second items of the record
    // homopol_df_ch1.collectFile(name: "${file_name}_homopol_summary.tsv").subscribe{it -> it.copyTo("${out_path}/files/${file_name}_homopol_summary.tsv")} // concatenate all the homopol_df.tsv files in channel homopol_df_ch into a single file published into the indicated path. STRONG WARNING: copyTo(${out_path}/files/) works only if the files folder already exists. Ohterwise, the saved file becomes "files"
    tsv_ch = homopolymer.out.homopol_df_ch.collectFile(name: "${sample_file.baseName}_homopol_summary.tsv")

    graph_stat(
        sample_file.baseName,
        cute_file, 
        min_length, 
        tsv_ch
    )

    backup(
        config_file,
        log_file
    )

    workflowVersion(
        modules
    )

    print_report(
        cute_file,
        initial.out.log_ch0.concat(graph_stat.out.log_ch1, backup.out.log_ch2, workflowVersion.out.log_ch3).collectFile(name: 'report.rmd', sort: false),
        graph_stat.out.fig_ch1.collect(),
        tsv_ch,
        graph_stat.out.table_ch1.collect()
    )




}