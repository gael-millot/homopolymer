#!/usr/bin/env Rscript

#########################################################################
##                                                                     ##
##     graph_stat.R                                                    ##
##                                                                     ##
##     Gael A. Millot                                                  ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Computational Biology Department                                ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################




################################ Aim


################################ End Aim


################################ Introduction


################################ End Introduction


################################ Acknowlegments


################################ End Acknowlegments


################################ Initialization


# R version checking
if(version$version.string != "R version 4.1.2 (2021-11-01)"){
    stop(paste0("\n\n================\n\nERROR IN plot_read_length.R\n", version$version.string, " IS NOT THE 4.1.2 RECOMMANDED\n\n================\n\n"))
}
# other initializations
erase.objects = TRUE # write TRUE to erase all the existing objects in R before starting the algorithm and FALSE otherwise. Beginners should use TRUE
if(erase.objects == TRUE){
    rm(list = ls(all.names = TRUE))
    erase.objects = TRUE
}
erase.graphs = TRUE # write TRUE to erase all the graphic windows in R before starting the algorithm and FALSE otherwise
script <- "graph_stat"


################################ End Initialization


################################ Parameters that need to be set by the user


################################ End Parameters that need to be set by the user


################################ Config import


tempo.cat <- "KIND OF RUN (SCRIPT, COPY-PASTE OR SOURCE): "
if(interactive() == FALSE){ # if(grepl(x = commandArgs(trailingOnly = FALSE), pattern = "R\\.exe$|\\/R$|Rcmd\\.exe$|Rcmd$|Rgui\\.exe$|Rgui$|Rscript\\.exe$|Rscript$|Rterm\\.exe$|Rterm$")){ # detection of script usage
    run.way <- "SCRIPT"
    cat(paste0("\n\n", tempo.cat, run.way, "\n"))
    command <- paste0(commandArgs(trailingOnly = FALSE), collapse = ",") # recover the full command
    args <- commandArgs(trailingOnly = TRUE) # recover arguments written after the call of the R script
    if(any(is.na(args))){
        stop(paste0("\n\n================\n\nERROR IN graph_stat.R\nTHE args OBJECT HAS NA\n\n================\n\n"), call. = FALSE)
    }
    tempo.arg.names <- c(
        "tsv", 
        "file_name", 
        "cute", 
        "log"
    ) # objects names exactly in the same order as in the bash code and recovered in args. Here only one, because only the path of the config file to indicate after the graph_stat.R script execution
    if(length(args) != length(tempo.arg.names)){
        stop(paste0("\n\n================\n\nERROR IN graph_stat.R\nTHE NUMBER OF ELEMENTS IN args (", length(args),") IS DIFFERENT FROM THE NUMBER OF ELEMENTS IN tempo.arg.names (", length(tempo.arg.names),")\nargs:", paste0(args, collapse = ","), "\ntempo.arg.names:", paste0(tempo.arg.names, collapse = ","), "\n\n================\n\n"), call. = FALSE)
    }
    for(i1 in 1:length(tempo.arg.names)){
        assign(tempo.arg.names[i1], args[i1])
    }
    rm(tempo.arg.names, args, i1)
}else if(sys.nframe() == 0L){ # detection of copy-paste/direct execution (for debugging). With script it is also 0, with source, it is 4
    run.way <- "COPY-PASTE"
    cat(paste0("\n\n", tempo.cat, run.way, "\n"))
}else{
    run.way <- "SOURCE" # using source(), sys.nframe() is 4
    cat(paste0("\n\n", tempo.cat, run.way, "\n"))
}
rm(tempo.cat)


################################ End Config import

################################ Test

# tsv <- "C:\\Users\\gael\\Documents\\Git_projects\\homopolymer\\example of results\\integrases.tsv"
# file_name <- "caca"
# cute <- "https://gitlab.pasteur.fr/gmillot/cute_little_R_functions/-/raw/v11.2.0/cute_little_R_functions.R"
# log <- "graph_stat_report.txt"



################################ end Test

################################ Recording of the initial parameters


param.list <- c(
    "erase.objects", 
    "erase.graphs", 
    "script", 
    "run.way",
    if(run.way == "SCRIPT"){"command"}, 
    "tsv", 
    "file_name", 
    "cute", 
    "log"
)
if(any(duplicated(param.list))){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 1 IN graph_stat.R\nTHE param.list OBJECT CONTAINS DUPLICATED ELEMENTS:\n", paste(param.list[duplicated(param.list)], collapse = " "), "\n\n================\n\n"), call. = FALSE) # message for developers
}
if(erase.objects == TRUE){
    created.object.control <- ls()[ ! ls() %in% "param.list"]
    if( ! (all(created.object.control %in% param.list) & all(param.list %in% created.object.control))){
        stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 2 IN graph_stat.R\nINCONSISTENCIES BETWEEN THE ARGUMENTS USED AND THE PARAMETERS REQUIRED IN THE EXECUTABLE CODE FILE\nTHE ARGUMENTS NOT PRESENT IN THE EXECUTABLE FILE (graph_stat.R) ARE:\n", paste(created.object.control[ ! created.object.control %in% param.list], collapse = " "), "\nTHE PARAMETERS OF THE EXECUTABLE FILE (graph_stat.R) NOT PRESENT IN THE ARGUMENTS ARE:\n", paste(param.list[ ! param.list %in% created.object.control], collapse = " "), "\n\n================\n\n"), call. = FALSE) # message for developers
    }
}
char.length <- nchar(param.list)
space.add <- max(char.length) - char.length + 5
param.ini.settings <- character(length = length(param.list))
for(i in 1:length(param.list)){
    param.ini.settings[i] <- paste0("\n", param.list[i], paste0(rep(" ", space.add[i]), collapse = ""), paste0(get(param.list[i]), collapse = ",")) # no env = sys.nframe(), inherit = FALSE in get() because look for function in the classical scope
}


################################ End Recording of the initial parameters


################################ Functions


# Functions are built such that they should have no direct use of Global objects (going through the R scope), and only use function arguments
# 1) Cute little function is sourced for the moment into the .GlobalEnv environment, but may be interesting to put it into a new environement just above .GlobalEnv environment. See https://stackoverflow.com/questions/9002544/how-to-add-functions-in-an-existing-environment
# 2) Argument names of each function must not be a name of Global objects (error message otherwise)
# 3) Argument name of each function ends with "_fun" in the first function, "_2fun" in the second, etc. This prevent conflicts with the argument partial names when using these functions, notably when they are imbricated


################ import functions from cute little functions toolbox


if(length(cute) != 1){
    stop(paste0("\n\n============\n\nERROR IN graph_stat.R\ncute PARAMETER MUST BE LENGTH 1: ", paste(cute, collapse = " "), "\n\n============\n\n"), call. = FALSE)
}else if(grepl(x = cute, pattern = "^http")){
    tempo.try <- try(suppressWarnings(suppressMessages(source(cute, local = .GlobalEnv))), silent = TRUE)
    if(any(grepl(x = tempo.try, pattern = "^[Ee]rror"))){
        stop(paste0("\n\n============\n\nERROR IN graph_stat.R\nHTTP INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else if( ! grepl(x = cute, pattern = "^http")){
    if( ! file.exists(cute)){
        stop(paste0("\n\n============\n\nERROR IN graph_stat.R\nFILE INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else{
    tempo.cat <- paste0("\n\n================\n\nINTERNAL CODE ERROR 3 IN graph_stat.R: CODE HAS TO BE MODIFIED\n\n============\n\n")
    stop(tempo.cat, call. = FALSE)
}


# required cute function checking
req.function <- c(
    "fun_check",
    "fun_pack", 
    "fun_df_remod", 
    "fun_gg_scatter", 
    "fun_gg_palette", 
    "fun_gg_empty_graph", 
    "fun_report"
)
tempo <- NULL
for(i1 in req.function){
    if(length(find(i1, mode = "function")) == 0L){
        tempo <- c(tempo, i1)
    }
}
if( ! is.null(tempo)){
    tempo.cat <- paste0("ERROR IN graph_stat.R\nREQUIRED cute FUNCTION", ifelse(length(tempo) > 1, "S ARE", " IS"), " MISSING IN THE R ENVIRONMENT:\n", paste0(tempo, collapse = "()\n"))
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end required function checking


################ local function: package import


# R Packages required
req.package.list <- c(
    "lubridate", 
    "ggplot2", 
    "lemon"
)
for(i in 1:length(req.package.list)){suppressMessages(library(req.package.list[i], character.only = TRUE))}
# fun_pack(req.package = req.package.list, load = TRUE, lib.path = NULL) # packages are imported even if inside functions are written as package.name::function() in the present code


################################ End Functions


################################ Pre-ignition checking


# reserved words
# end reserved words
# argument primary checking
arg.check <- NULL #
text.check <- NULL #
checked.arg.names <- NULL # for function debbuging: used by r_debugging_tools
ee <- expression(arg.check <- c(arg.check, tempo$problem) , text.check <- c(text.check, tempo$text) , checked.arg.names <- c(checked.arg.names, tempo$object.name))
tempo <- fun_check(data = tsv, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = file_name, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = cute, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = log, class = "vector", typeof = "character", length = 1) ; eval(ee)
if(any(arg.check) == TRUE){ # normally no NA
    stop(paste0("\n\n================\n\n", paste(text.check[arg.check], collapse = "\n"), "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between == #
}
# end argument primary checking
# second round of checking and data preparation
# management of NA arguments
# end management of NA arguments
# management of NULL arguments
tempo.arg <-c(
    "tsv", 
    "file_name", 
    "cute", 
    "log"
)
tempo.log <- sapply(lapply(tempo.arg, FUN = get, env = sys.nframe(), inherit = FALSE), FUN = is.null)
if(any(tempo.log) == TRUE){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN graph_stat.R:\n", ifelse(sum(tempo.log, na.rm = TRUE) > 1, "THESE ARGUMENTS\n", "THIS ARGUMENT\n"), paste0(tempo.arg[tempo.log], collapse = "\n"),"\nCANNOT BE NULL")
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end management of NULL arguments
# code that protects set.seed() in the global environment
# end code that protects set.seed() in the global environment
# warning initiation
ini.warning.length <- options()$warning.length
options(warning.length = 8170)
warn <- NULL
# warn.count <- 0 # not required
# end warning initiation
# other checkings
# end other checkings
# reserved word checking
# end reserved word checking
# end second round of checking and data preparation
# package checking
# end package checking


################################ End pre-ignition checking


################################ Main code


################ Ignition


fun_report(data = paste0("\n\n################################################################ graph_stat PROCESS\n\n"), output = log, path = "./", overwrite = FALSE)
ini.date <- Sys.time()
ini.time <- as.numeric(ini.date) # time of process begin, converted into seconds
fun_report(data = paste0("\n\n################################ RUNNING DATE AND STARTING TIME\n\n"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0(ini.date, "\n\n"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0("\n\n################################ RUNNING\n\n"), output = log, path = "./", overwrite = FALSE)


################ End ignition


################ Graphical parameter initialization


pdf(file = NULL)
par.ini <- par(no.readonly = TRUE) # to recover the initial graphical parameters if required (reset)
invisible(dev.off()) # close the new window
zone.ini <- matrix(1, ncol=1)
if(erase.graphs == TRUE){
    graphics.off()
}else{
    tempo.warn <- paste0("GRAPHICS HAVE NOT BEEN ERASED. GRAPHICAL PARAMETERS MAY HAVE NOT BEEN REINITIALIZED")
    fun_report(data = paste0("WARNING\n", tempo.warn), output = log, path = "./", overwrite = FALSE)
    warn <- paste0(ifelse(is.null(warn), tempo.warn, paste0(warn, "\n\n", tempo.warn)))
}


################ End graphical parameter initialization


################ Data import


tsv <- read.table(tsv, stringsAsFactors = FALSE, header = FALSE, sep = "\t")


################ end Data import


############ modifications of imported tables

# export tsv table with header
tempo <- tsv
tempo_names <- c('name', 'seq_length', 'max_homopol_size', 'nucleotide', 'starting_position', 'relative_position', "nb", 'mean_size', 'homopol_obs_distrib', 'homopol_theo_distrib')
if(length(tempo) != length(tempo_names)){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN graph_stat.R:\nLENGTH OF THE IMPORTED TABLE IS NOT 10: ", length(tempo))
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}else{
    names(tempo) <- tempo_names
    write.table(tempo, file = paste0("./homopol_summary.tsv"), row.names = FALSE, col.names = TRUE, append = FALSE, quote = FALSE, sep = "\t")
}
tsv <- tempo

# recovering of the freq distributions
obs <- tsv[ , 'homopol_obs_distrib']
theo <- tsv[ , 'homopol_theo_distrib']
obs <- strsplit(obs, split = ";") # list
theo <- strsplit(theo, split = ";") # list
max_cells <- max(sapply(X = c(obs, theo), FUN = length), na.rm = TRUE)
obs <- lapply(X = obs, FUN = function(x = X, y = max_cells){c(x, rep("0", y - length(x)))})
theo <- lapply(X = theo, FUN = function(x = X, y = max_cells){c(x, rep("0", y - length(x)))})
obs <- as.matrix(as.data.frame(obs))
theo <- as.matrix(as.data.frame(theo))
colnames(obs) <- tsv[ , 'name']
colnames(theo) <- tsv[ , 'name']
mode(obs) <- "double"
mode(theo) <- "double"
# end recovering of the freq distributions
# recovering of the prop distributions
obs.prop <- apply(obs, MARGIN = 2, FUN = function(x){x / sum(x, na.rm = TRUE)})
theo.prop <- apply(theo, MARGIN = 2, FUN = function(x){x / sum(x, na.rm = TRUE)})
# end recovering of the prop distributions

# data for barplot and scatterplot
obs2 <- apply(X = obs, MARGIN = 1, FUN = sum)
theo2 <- apply(X = theo, MARGIN = 1, FUN = sum)
max1 <- which.max(diff(obs2 == 0)) # to determine where there are only zero after the last number
max2 <- which.max(diff(theo2 == 0))
obs2 <- obs2[1:max(max1, max2)]
theo2 <- theo2[1:max(max1, max2)]
obs3 <- obs[1:max(max1, max2), ]
theo3 <- theo[1:max(max1, max2), ]
obs3 <- data.frame(fun_df_remod(as.data.frame(t(obs3))), kind = "obs")
theo3 <- data.frame(fun_df_remod(as.data.frame(t(theo3))), kind = "theo")
final3 <- rbind(obs3, theo3)
names(final3) <- c("freq", "length", "name" , "kind")
final3$length <- as.numeric(gsub(x = as.character(final3$length), pattern = "V", replacement = ""))
stat <- aggregate(final3$freq, by = list(final3$length, final3$kind), FUN = mean)
names(stat) <- c("length", "kind", "mean")
stat2 <- aggregate(final3$freq, by = list(final3$length, final3$kind), FUN = sd)
stat <- data.frame(stat, sd = stat2$x, CI95.inf = stat$mean - 1.96 * stat2$x, CI95.sup = stat$mean + 1.96 * stat2$x)
# write.table(stat, file = paste0("./scatterplot_stat.tsv"), row.names = FALSE, col.names = TRUE, append = FALSE, quote = FALSE, sep = "\t")
final3$graph.length <- sprintf("%02.0f", final3$length)
chi2.table <- matrix(c(obs2, theo2), ncol = 2, byrow = FALSE)
chi2.table <- chi2.table[ ! apply(chi2.table, 1, FUN = function(x){all(x == 0, na.rm = TRUE)}), ]
tempo <- suppressWarnings(suppressMessages(prop.test(chi2.table, correct = TRUE)[]))
tempo2 <- data.frame(Parameter = c("Statistic", "Df", "P value", "Method", "Alternative"), Value = c(tempo$statistic, tempo$parameter, tempo$p.value, tempo$method, tempo$alternative))
# write.table(tempo2, file = paste0("./chi2.tsv"), row.names = FALSE, col.names = TRUE, append = FALSE, quote = FALSE, sep = "\t")
# end data for barplot and scatterplot

# data for prop barplot and scatterplot
obs2.prop <- apply(X = obs.prop, MARGIN = 1, FUN = sum)
theo2.prop <- apply(X = theo.prop, MARGIN = 1, FUN = sum)
max1.prop <- which.max(diff(obs2.prop == 0)) # to determine where there are only zero after the last number
max2.prop <- which.max(diff(theo2.prop == 0))
obs2.prop <- obs2.prop[1:max(max1.prop, max2.prop)]
theo2.prop <- theo2.prop[1:max(max1.prop, max2.prop)]
obs3.prop <- obs.prop[1:max(max1.prop, max2.prop), ]
theo3.prop <- theo.prop[1:max(max1.prop, max2.prop), ]
obs3.prop <- data.frame(fun_df_remod(as.data.frame(t(obs3.prop))), kind = "obs")
theo3.prop <- data.frame(fun_df_remod(as.data.frame(t(theo3.prop))), kind = "theo")
final3.prop <- rbind(obs3.prop, theo3.prop)
names(final3.prop) <- c("freq", "length", "name" , "kind")
final3.prop$length <- as.numeric(gsub(x = as.character(final3.prop$length), pattern = "V", replacement = ""))
stat.prop <- aggregate(final3.prop$freq, by = list(final3.prop$length, final3.prop$kind), FUN = mean)
names(stat.prop) <- c("length", "kind", "mean")
stat2.prop <- aggregate(final3.prop$freq, by = list(final3.prop$length, final3.prop$kind), FUN = sd)
stat.prop <- data.frame(stat.prop, sd = stat2.prop$x, CI95.inf = stat.prop$mean - 1.96 * stat2.prop$x, CI95.sup = stat.prop$mean + 1.96 * stat2.prop$x)
write.table(stat.prop, file = paste0("./scatterplot_stat.tsv"), row.names = FALSE, col.names = TRUE, append = FALSE, quote = FALSE, sep = "\t")
final3.prop$graph.length <- sprintf("%02.0f", final3.prop$length)
# end data for prop barplot and scatterplot



# data for t.test
categ <- NULL
obs.mean <- NULL
theo.mean <- NULL
obs.sd <- NULL
theo.sd <- NULL
df <- NULL
t <- NULL
p.value<-NULL
for(i0 in unique(sort(final3$graph.length))){
    categ <- c(categ, i0)
    obs.tempo <- final3$freq[final3$graph.length == i0 & final3$kind == "obs"]
    theo.tempo <- final3$freq[final3$graph.length == i0 & final3$kind == "theo"]
    obs.mean <- c(obs.mean, mean(obs.tempo, na.rm = TRUE))
    theo.mean <- c(theo.mean, mean(theo.tempo, na.rm = TRUE))
    obs.sd <- c(obs.sd, sd(obs.tempo, na.rm = TRUE))
    theo.sd <- c(theo.sd, sd(theo.tempo, na.rm = TRUE))
    t.test.tempo <- t.test(obs.tempo, theo.tempo, var.equal = FALSE)
    df <- c(df, t.test.tempo$parameter)
    t <- c(t, t.test.tempo$statistic)
    p.value <- c(p.value, t.test.tempo$p.value)
}
p.mult<-data.frame(categ, obs.mean, theo.mean, obs.sd, theo.sd, df, t, p.value, BH.adj.p.value = p.adjust(p.value, method = "BH"))
# write.table(p.mult, file = paste0("./t_test.tsv"), row.names = FALSE, col.names = TRUE, append = FALSE, quote = FALSE, sep = "\t")
# end data for t.test

# data for prop t.test
categ.prop <- NULL
obs.mean.prop <- NULL
theo.mean.prop <- NULL
obs.sd.prop <- NULL
theo.sd.prop <- NULL
df.prop <- NULL
t.prop <- NULL
p.value.prop <-NULL
for(i0 in unique(sort(final3.prop$graph.length))){
    categ.prop <- c(categ.prop, i0)
    obs.tempo.prop <- final3.prop$freq[final3.prop$graph.length == i0 & final3.prop$kind == "obs"]
    theo.tempo.prop <- final3.prop$freq[final3.prop$graph.length == i0 & final3.prop$kind == "theo"]
    obs.mean.prop <- c(obs.mean.prop, mean(obs.tempo.prop, na.rm = TRUE))
    theo.mean.prop <- c(theo.mean.prop, mean(theo.tempo.prop, na.rm = TRUE))
    obs.sd.prop <- c(obs.sd.prop, sd(obs.tempo.prop, na.rm = TRUE))
    theo.sd.prop <- c(theo.sd.prop, sd(theo.tempo.prop, na.rm = TRUE))
    t.test.tempo.prop <- t.test(obs.tempo.prop, theo.tempo.prop, var.equal = FALSE)
    df.prop <- c(df.prop, t.test.tempo.prop$parameter)
    t.prop <- c(t.prop, t.test.tempo.prop$statistic)
    p.value.prop <- c(p.value.prop, t.test.tempo.prop$p.value)
}
p.mult.prop <- data.frame(categ.prop, obs.mean.prop, theo.mean.prop, obs.sd.prop, theo.sd.prop, df.prop, t.prop, p.value.prop, BH.adj.p.value = p.adjust(p.value.prop, method = "BH"))
names(p.mult.prop) <- c("length", "obs_mean", "theo_mean", "obs_sd", "theo_sd", "df", "t", "p.value", "BH.adj.p.value")
write.table(p.mult.prop, file = paste0("./t_test.tsv"), row.names = FALSE, col.names = TRUE, append = FALSE, quote = FALSE, sep = "\t")
# end data for prop t.test




############ end modifications of imported tables


############ plotting


png(filename = paste0("plot_", file_name, ".png"), width = 5000, height = 1800, units = "px", res = 300)
if(nrow(final3.prop) > 0){
    fun_gg_scatter(
        data1 = list(final3.prop, stat.prop, stat.prop, stat.prop), # res # res[res$KIND == "obs.freq", ]
        x = list("length", "length", "length", "length"), 
        y = list("freq", "mean", "CI95.inf", "CI95.sup"), 
        categ = list("kind", "kind", "kind", "kind"), 
        geom = list("geom_point", "geom_line", "geom_line", "geom_line"), 
        line.size = list(NULL, 2, 1, 1), 
        color = list(fun_gg_palette(n = 2, kind = "dark")[c(1, 2)], fun_gg_palette(n = 2, kind = "std")[c(1, 2)], fun_gg_palette(n = 2, kind = "light")[c(1, 2)], fun_gg_palette(n = 2, kind = "light")[c(1, 2)]), # fun_gg_palette(n = 2) # fun_gg_palette(n = 2)[1]
        dot.size = 4, 
        dot.shape = 21, 
        dot.border.size = 0.5, 
        dot.border.color = NULL, 
        alpha = list(0.1, 0.5, 0.5, 0.5), 
        line.type = "solid",
        legend.width = 0.2, 
        legend.name = list("Kind", "Means", "Lower CI95", "Upper CI95"), 
        title = "", 
        x.lab = "Homopolymer length", 
        x.left.extra.margin = 0.05, 
        y.top.extra.margin = 0.05, 
        y.bottom.extra.margin = 0, 
        x.right.extra.margin = 0.05, 
        x.second.tick.nb = NULL, 
        y.lab = "Proportion", 
        y.log = "no", 
        y.second.tick.nb = 5, 
        text.size = 24, 
        title.text.size = 16
    )
}else{
    fun_gg_empty_graph(text = "EMPTY .tsv FILE: NO PLOT DRAWN")
}


png(filename = paste0("boxplot_", file_name, ".png"), width = 5000, height = 1800, units = "px", res = 300)
# tempo <- data.frame(length = c((1:length(obs2.prop)) - 0.1, (1:length(obs2)) + 0.1), freq = c(obs2.prop, theo2.prop), kind = rep(c("obs", "theo"), each = length(obs2.prop)))
# write.table(tempo, file = paste0("./boxplot.tsv"), row.names = FALSE, col.names = TRUE, append = FALSE, quote = FALSE, sep = "\t")
if(nrow(final3.prop) > 0){
    tempo <- fun_gg_boxplot(
        data1 = final3.prop, # res # res[res$KIND == "obs.freq", ]
        categ = c("graph.length", "kind"),  
        y = "freq", 
        legend.width = 0.2, 
        title = "", 
        x.lab = "Homopolymer length", 
        y.top.extra.margin = 0.05, 
        y.bottom.extra.margin = 0, 
        y.lab = "Proportion", 
        y.log = "no", 
        y.second.tick.nb = 5, 
        text.size = 24, 
        title.text.size = 16,
        return = TRUE
    )
    tempo <- tempo[, c("BOX", "MIN", "QUART1", "MEDIAN", "MEAN", "QUART3", "MAX", "WHISK_INF", "BOX_INF", "NOTCH_INF", "NOTCH_SUP", "BOX_SUP", "WHISK_SUP", "OUTLIERS")]
    write.table(as.matrix(tempo$stat), file = paste0("./boxplot_stat.tsv"), row.names = FALSE, col.names = TRUE, append = FALSE, quote = FALSE, sep = "\t")
}else{
    fun_gg_empty_graph(text = "EMPTY .tsv FILE: NO PLOT DRAWN")
    write.table("", file = paste0("./boxplot_stat.tsv"), row.names = FALSE, col.names = TRUE, append = FALSE, quote = FALSE, sep = "\t")
}



############ end plotting


################ Pdf window closing


graphics.off()


################ end Pdf window closing


################ Seeding inactivation


set.seed(NULL)


################ end Seeding inactivation


################ Environment saving


save(list = ls(), file = "graph_stat.RData")
fun_report(data = paste0("\n\n################################ RUNNING END"), output = log, path = "./", overwrite = FALSE)
end.date <- Sys.time()
end.time <- as.numeric(end.date)
total.lapse <- round(lubridate::seconds_to_period(end.time - ini.time))
fun_report(data = paste0("\n\nEND TIME: ", end.date), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0("\n\nTOTAL TIME LAPSE: ", total.lapse), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0("\n\nALL DATA SAVED IN graph_stat.RData"), output = log, path = "./", overwrite = FALSE)


################ end Environment saving


################ Warning messages


fun_report(data = paste0("\n\n################################ RECAPITULATION OF WARNING MESSAGES"), output = log, path = "./", overwrite = FALSE)
if( ! is.null(warn)){
    fun_report(data = paste0("\n\n", warn), output = log, path = "./", overwrite = FALSE)
}else{
    fun_report(data = paste0("\n\nNO WARNING MESSAGE TO REPORT"), output = log, path = "./", overwrite = FALSE)
}


################ end Warning messages


################ Parameter printing


fun_report(data = paste0("\n\n################################ INITIAL SETTINGS OF PARAMETERS"), output = log, path = "./", overwrite = FALSE)
fun_report(data = param.ini.settings, output = log, path = "./", overwrite = FALSE, , vector.cat = TRUE)
fun_report(data = paste0("\n\n################################ R SYSTEM AND PACKAGES"), output = log, path = "./", overwrite = FALSE)
tempo <- sessionInfo()
tempo$otherPkgs <- tempo$otherPkgs[order(names(tempo$otherPkgs))] # sort the packages
tempo$loadedOnly <- tempo$loadedOnly[order(names(tempo$loadedOnly))] # sort the packages
fun_report(data = tempo, output = log, path = "./", overwrite = FALSE, , vector.cat = TRUE)
fun_report(data = paste0("\n\n################################ JOB END\n\nTIME: ", end.date, "\n\nTOTAL TIME LAPSE: ", total.lapse, "\n"), output = log, path = "./", overwrite = FALSE)


################ end Parameter printing


################################ End Main code







