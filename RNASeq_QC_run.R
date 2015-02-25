#!/usr/bin/env Rscript

## Copyright (c) 2015 by Joseph Tran. This work is made available under the terms of the Creative Commons Attribution-ShareAlike 3.0 license, http://creativecommons.org/licenses/by-sa/3.0/.

################################################
##
## @script: RNASeq_QC_run.R
##
## @author: Joseph Tran (Joseph.Tran@versailles.inra.fr)
##
## @version: 0.0.1.0
##
## @date: 2015-01-20
##
## @description: This R script is a wrapper to compile RNASeq Quality Control (QC) Rmd documents to pdf documents containing QC figures. 
##
###############################################

version <- "0.0.1.0"

copyright <- "Copyright (c) 2015 by Joseph Tran. This work is made available under the terms of the Creative Commons Attribution-ShareAlike 3.0 license, http://creativecommons.org/licenses/by-sa/3.0/."

## execution time: start
T1<-Sys.time()


##
## Options and usage
##
### optparse
if (is.element('optparse', installed.packages()[,1]))
{
  suppressPackageStartupMessages(require('optparse'));
} else
{
  install.packages('optparse');
  suppressPackageStartupMessages(library('optparse'));
}





option_list <- list(
  make_option(c('-v', '--verbosity-level'), type = "integer", dest = "verbosity", default = 1, help = "Verbosity threshold (5=DEBUG, 4=INFO 3=WARN, 2=ERROR, 1=FATAL)"),
  make_option(c("-c", "--count"), type="character", help="The count dataset input file: 2 formats are accepted, BBRIC or generic (see notes) [mandatory]", metavar="count_input_file"),
   make_option(c("-f", "--format"), type="character", default = "BBRIC", help="The count dataset input file format: 2 values are accepted, BBRIC or generic (see notes) [mandatory]", metavar="count_input_file_format"),
  make_option(c("-s", "--stats"), type="character", help="The BBRIC mapping statistics dataset input file", metavar="stats_input_file"),
  make_option(c("-d", "--design"), type="character", help="The experimental design input file describing the different factors and modalities in the data (see notes)", metavar="design_input_file"),
  make_option(c("-o", "--outdir"), type="character", default="rnaseq_qc_out", help="The output directory [default %default]", metavar="OUTPUT_DIRECTORY")
) 
parser <- OptionParser(usage = "usage:  %prog [options]", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE, print_help_and_exit=FALSE)
opt <- arguments$options
if (opt$help) {
  cat(paste("author:  Joseph Tran <Joseph.Tran@versailles.inra.fr>\n\nversion: ", version, "\n\n"))
  print_help(parser)
  cat(paste("notes:\n\t1. ", copyright, "\n\n"), sep="")
  cat(paste("\n\t2. ", "BBRIC format (output from Expression Measure protocol):", "\n\t", 
            "example: with 2 libs Lib1 and Lib2 having each 2 columns (-count and -rpkm)",
            "\n\t", "id  0_seqid	1_start	2_end	3_strand	4_length	5_type	6_Note	Lib1-count	Lib1-rpkm	Lib2-count	Lib2-rpkm",
            "\n\t", "gene1	Contig_1.5	1294	2901	-	1608	gene		0	0.00	0	0.00",
            "\n\t", "gene2	Contig_1.5	13364	15576	-	2213	gene		314	14.79	497	16.75",
            "\n\t...",
            "with as many libs as needed.", "\n\n",
            sep=""))
  cat(paste("\n\t3. ", "Generic format (tabulated):", "\n\t",
            "example: with 2 column-libs Lib1 and Lib2",
            "\n\t", "id  Lib1	Lib2",
            "\n\t", "gene1	0	0",
            "\n\t", "gene2	314	497",
            "\n\t", "gene3	10991	11553",
            "\n\t", "gene4	135	192",
            "\n\t", "gene5	294	449",
            "\n\t...", 
            "with as many libs as needed.", "\n\n",
            sep=""))
  cat(paste("\n\t4. ", "Experimental design file format (tabulated):", "\n\t",
            "example: 2 column-factors (genotype and repet) with thier modalities genotype: Mutant1, Mutant2, etc. ; repet: M1, M2, etc.",
            "\n\t", "lib genotype  repet",
            "\n\t", "Mutant1_M2  Mutant1 M2",
            "\n\t", "Mutant1_M3  Mutant1 M3",
            "\n\t", "Mutant1_M4  Mutant1 M4",
            "\n\t", "Mutant2_M1  Mutant2 M1",
            "\n\t", "Mutant2_M2  Mutant2 M2",
            "\n\t...", "\n\n",
            sep=""))
  quit()
}

##
## Load and install libraries
##

### RNASeq_QC_lib.R
#initial.options <- commandArgs(trailingOnly = FALSE)
#file.arg.name <- "--file="
#script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
#script.dir <- dirname(script.name)
#source(normalizePath(file.path(script.dir, "lib/RNASeq_QC_lib.R")), chdir=T)
suppressPackageStartupMessages(library(qc4rnaseq))
suppressPackageStartupMessages(library(log4r))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(markdown))
suppressPackageStartupMessages(library(rmarkdown))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(FactoMineR))



##
## Logger
##
logger <- create.logger(logfile = "", level = verbosity(opt$verbosity))

##
## Input files validation
##

###
### Check for count data input file
###
if (!file.exists(opt$count) || (file.access(opt$count, mode=4) == -1)) {
  error(logger, paste("Count input file ", opt$count, " does not exist or is not readable.", sep=""))
  stop(sprintf("Specified count file ( %s ) does not exist or is not readable.", opt$count))      
} else
{
  debug(logger, paste("OK Count input file", opt$count, " does exist and is readable", sep=""))
}

###
### Check for count data input file format
###
formats <- c("BBRIC", "generic")
if (!(opt$format %in% formats)) {
  error(logger, paste("Count input file format", opt$format, " is not a valid value.", sep=""))
  stop(sprintf("Specified count file format ( %s ) is not valid.", opt$format))      
} else
{
  debug(logger, paste("OK Count input file format", opt$format, " is valid.", sep=""))
}


###
### Check for stats data input file
###
if (!is.null(opt$stats)) {
  if (!file.exists(opt$stats) || (file.access(opt$stats, mode=4) == -1)) {
    warn(logger, paste("Mapping statistics input file ", opt$stats, " does not exist or is not readable.", sep=""))
  } else
  {
    debug(logger, paste("OK Mapping statistics input file ", opt$stats, " does exist and is readable", sep=""))
  }
}

###
### Check for design data input file
###
if (!is.null(opt$design)) {
  if (!file.exists(opt$design) || (file.access(opt$design, mode=4) == -1)) {
    warn(logger, paste("Experimental design input file ", opt$design, " does not exist or is not readable.", sep=""))
  } else
  {
    debug(logger, paste("OK Experimental design input file ", opt$design, " does exist and is readable", sep=""))
  }
}

##
## Create output directory
##
if (!file.exists(opt$outdir)){
  info(logger, paste("Creating output directory, ", opt$outdir, ", ...", sep=""))
  dir.create(file.path(opt$outdir))
  info(logger, paste("OK Output directory ", opt$outdir, " created.", sep=""))
} else
{
  info(logger, paste("Output directory, ", opt$outdir, ", already exists.", sep=""))
  if (!file.access(opt$outdir, mode=4) == 0) {
    error(logger, paste("Output directory ", opt$outdir, " is not readable.", sep=""))
    stop(sprintf("Specified output directory ( %s ) is not readable.", opt$outdir))
  } else
  {
    debug(logger, paste("OK output directory ", opt$outdir, " is readable", sep="")) 
  }
}

##
## Check count data format
##

### load count data
info(logger, paste("Loading count data, ", opt$count, ", ...", sep=""))
count.df <- loadCountData(opt$count)
info(logger, paste("Loading count data, ", opt$count, " done", sep=""))
debug(logger, paste("input count data dimensions: ", dim(count.df)[1], " x ", dim(count.df)[2], sep=""))

### check count format
# defaults
is_count_format = FALSE
warn_err <- tryCatch.W.E(isCountDataFormat(count.df, format=opt$format))
if (is.null(warn_err$warning)) {
	is_count_format <- warn_err$value
} else
	{
		stop(paste(geterrmessage(), warn_err$warning$message))
	}


# set count data frames
is_bbric_format <- FALSE
is_generic_format <- FALSE
if (is_count_format && opt$format == "BBRIC") {
  # if bbric
  info(logger, "Set bbric count variable")
  is_bbric_format <- TRUE
  
} else if (is_count_format && opt$format == "generic") {
  # if generic
  info(logger, "Set generic count variable")
  is_generic_format <- TRUE
  
} else {
  # count format not supported
  error(logger, "Count input file format is not supported.")
  stop("Count input file format is not supported.")
}  


##
## Check stats data format(optional)
##
# default
is_stats_format <- TRUE
if (!is.null(opt$stats)) {
	### load stats data 
	info(logger, paste("Loading stats data, ", opt$stats, ", ...", sep=""))
	stats.df <- loadStatsData(opt$stats)
	info(logger, paste("Loading stats data, ", opt$stats, " done", sep=""))
	debug(logger, paste("input stats data dimensions: ", dim(stats.df)[1], " x ", dim(stats.df)[2], sep=""))
	### check stats format
  stats_warn_err <- tryCatch.W.E(isStatsDataFormat(stats.df))
  if (is.null(stats_warn_err$warning)) {
	  is_stats_format <- stats_warn_err$value
  } else
	{
		stop(paste(geterrmessage(), stats_warn_err$warning$message))
	}
}

##
## Check design data format (optional)
##
# default
is_design_format <- TRUE
if (!is.null(opt$design)) {
	### load design data
	info(logger, paste("Loading design data, ", opt$design, ", ...", sep=""))
	design.df <- loadDesignData(opt$design)
	info(logger, paste("Loading design data, ", opt$design, " done", sep=""))
	debug(logger, paste("input design data dimensions: ", dim(design.df)[1], " x ", dim(design.df)[2], sep=""))
	### check design format
  design_warn_err <- tryCatch.W.E(isDesignDataFormat(design.df))
  if (is.null(design_warn_err$warning)) {
	  is_design_format <- design_warn_err$value
  } else
	{
		stop(paste(geterrmessage(), design_warn_err$warning$message))
	}
}

##
## Check count data and design data consistency (optional)
##
# default
is_count_design <- TRUE 
if (!is.null(opt$design)) {
	if (is_bbric_format || is_generic_format) {
    count_design_warn_err <- tryCatch.W.E(isCountDesign(count.df, design.df, format=opt$format))
    if (is.null(count_design_warn_err$warning)) {
	    is_count_design <- count_design_warn_err$value
    } else
	  {
		  stop(paste(geterrmessage(), count_design_warn_err$warning$message))
	  }
	
	}
}

##
## Check stats data and design data consistency (optional)
##
# default
is_stats_design <- TRUE
if ((!is.null(opt$stats)) && (!is.null(opt$design))) {
  stats_design_warn_err <- tryCatch.W.E(isStatsDesign(stats.df, design.df))
  if (is.null(stats_design_warn_err$warning)) {
	  is_stats_design <- stats_design_warn_err$value
  } else
	{
		stop(paste(geterrmessage(), stats_design_warn_err$warning$message))
	}
}



##
## Execute the .Rmd (prepare the plots in a .pdf file)
##


# Gather all checks and launch the appropriate.Rmd script
if (all(is_stats_format, is_design_format, is_count_design, is_stats_design)) {
	# BBRIC or generic format
	if (is_bbric_format){
		bbric_report_path <- system.file("report", "QC_RNASeq_Count_BBRIC.Rmd", package="qc4rnaseq")
		bbric_report_file <- paste(unlist(strsplit(basename(bbric_report_path),".Rmd")), ".pdf", sep="")
		outdir_abs_path <- normalizePath(opt$outdir)
		render(input=bbric_report_path, output_format="pdf_document", output_file=bbric_report_file, output_dir=outdir_abs_path, intermediates_dir=outdir_abs_path)
}	 else if (is_generic_format) {
		generic_report_path <- system.file("report", "QC_RNASeq_Count_generic.Rmd", package="qc4rnaseq")
		render(input=generic_report_path, output_format="pdf_document", output_dir=opt$outdir)
	}
	if (!is.null(opt$stats)) {
		stats_report_path <- system.file("report", "QC_RNASeq_Stats_BBRIC.Rmd", package="qc4rnaseq")
		render(input=stats_report_path, output_format="pdf_document", output_dir=opt$outdir)
	}
} else {
	error(logger, "The stats and/or design format is inappropriate.")
	stop("The stats and/or design format is inappropriate.")
}



## execution time: end
T2<-Sys.time()
Tdiff= difftime(T2, T1)
cat("Execution time : ", Tdiff,"seconds\n")
