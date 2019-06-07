#' @title Utility functions for R
#' @author Jonatan Thompson, Perslab, rkm916 at ku dot dk

############################################################################################################################################################
####################################################### GET SCRIPT DIRECTORY ###############################################################################
############################################################################################################################################################

LocationOfThisScript = function() # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
  #' @usage returns the current location of the script
  #' @value directory of the script, character
  
  if (interactive()) {
    stop("LocationOfThisScript does not work in interactive sessions")
    }
    
  this.file = NULL
  # This file may be 'sourced'
  for (i in -(1:sys.nframe())) {
    if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
  }
  
  if (!is.null(this.file)) return(dirname(this.file))
  
  # But it may also be called from the command line
  cmd.args = commandArgs(trailingOnly = FALSE)
  cmd.args.trailing = commandArgs(trailingOnly = TRUE)
  cmd.args = cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
  res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)
  
  # If multiple --file arguments are given, R uses the last one
  res = tail(res[res != ""], 1)
  if (0 < length(res)) return(dirname(res))
  
  # Both are not the case. Maybe we are in an R GUI?
  return(NULL)
}

############################################################################################################################################################
########################################################### INSTALL PACKAGES ###############################################################################
############################################################################################################################################################

# install and require packages using a function (allows for automation)
ipak <- function(pkgs){
  #' @usage attempt to install and load packages 
  #' @example 
  #' packages <- c("ggplot2", "plyr", "reshape2", "RColorBrewer", "scales", "grid")
  #  ipak(packages)
  #' @param pkgs: vector of packages to install/load; character
  #' @value none 
  
  new.pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new.pkgs)) {
    sapply(new.pkgs, function(pkg) {
      tryCatch({
        install.packages(pkg, dependencies = TRUE)
      }, warning = function(war) {
        tryCatch({
          source("https://bioconductor.org/biocLite.R")
          biocLite(pkg, suppressUpdates = T)
        }, error = function(err1) {
          warning(paste0(pkg, " encountered the error: ", err1))
          dependency <- gsub("\\W|ERROR: dependency | is not available for package.*", "", err)
          ipak(dependency)
        })
      } ,
      error = function(err)
      {
        dependency <- gsub("\\W|ERROR: dependency | is not available for package.*", "", err)
        ipak(dependency)
      })
    })
  }
  suppressPackageStartupMessages(sapply(pkgs, require, character.only = TRUE))
  failed <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(failed)>0) warning(paste0(paste0(failed, collapse = " "), " failed to install"))
}

############################################################################################################################################################
####################################################### IF INSTALL_GITHUB FAILS .. ####################################################################
############################################################################################################################################################
# download the file, unzip, then
# install.packages("/projects/jonatan/tools/SoupX/SoupX-master/", repos=NULL,type="source")
#install.packages("/projects/jonatan/tools/LTMGSCA/LTMGSCA-master/", repos=NULL,type="source")
#install.packages("/projects/jonatan/tools/Rpackages3.5/Seurat/", repos=NULL,type="source", lib = "/projects/jonatan/tools/Rpackages3.4/")
#install.packages("/projects/jonatan/R_repository/seurat-release-3.0/", repos=NULL,type="source")

rmPkgs <- function() {
  #' @usage detacha and unload all packages
  #' @value none
  #' # https://stackoverflow.com/questions/7505547/detach-all-packages-while-working-in-r
  while (!is.null(names(sessionInfo()$otherPkgs))) lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
}
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
# Remove packages installed under 3.5.0
# Command line: R CMD REMOVE [options] [-l lib] pkgs
# pak3.5 <- installed.packages()[,"Package"][grep("3.5.0", installed.packages()[,"Built"])]

############################################################################################################################################################
############################################################### DOWNLOAD FILES AND CHECK INTEGRITY #########################################################
############################################################################################################################################################

dlfile <- function(url, destfile, correct_checksum) {
  #' @usage: check md5sum and download file if needed
  #' @param url: passed to download.file(). 
  #' @param destfile: passed to download.file().
  #' @param correct_checksum: to check file integrity
  #' @value: none 
  
  if(!file.exists(destfile)){
    print("Downloading file")
    download.file(url, destfile, quiet = FALSE)
  } else{
    print("Verifying file integrity...")
    checksum = md5sum(destfile)
    if(checksum != correct_checksum){
      print("Existing file looks corrupted or is out of date, downloading again.")
      try(download.file(url, destfile, quiet = FALSE))
    } else{
      print("Latest file already exists.")
    }
  }
}

############################################################################################################################################################
################################################################ READ FILES INTO R SESSION #################################################################
############################################################################################################################################################

load_obj <- function(f) {
  #' @usage Loads (gzip compressed) file from .RData, .RDS, .loom, .csv, .txt, .tab, .delim
  #' @param f: path to file
  #' @value RObject
  
  compressed = F
  
  if (grepl(pattern = "\\.gz|\\.gzip", x=f))  {
    compressed <- T
    f = paste0("gzfile('",f,"')")
  }
  
  if (grepl(pattern = "\\.RDS", x = f, ignore.case = T)) {
    out <- readRDS(file=if(compressed) eval(parse(text=f)) else f)
  } else if (grepl(pattern="\\.RData|\\.Rda", x=f, ignore.case = T)) { 
    env <- new.env()
    nm <- load(f, env)[1]
    out <- env[[nm]]
  } else if (grepl(pattern="\\.loom", x=f)) {
    out <- connect(filename=f, mode = "r+")
  } else if (grepl(pattern = "\\.csv", x=f)) {
    out <- read.csv(file=if(compressed) eval(parse(text=f)) else f, stringsAsFactors = F, quote="", header=T)
  } else if (grepl(pattern = "\\.tab|\\.tsv", x=f)) {
    out <- read.table(file=if(compressed) eval(parse(text=f)) else f, sep="\t", stringsAsFactors = F, quote="", header=T) 
  } else if (grepl(pattern = "\\.txt", x=f)) {
    out <- read.delim(file=if(compressed) eval(parse(text=f)) else f, stringsAsFactors = F, quote="", header=T)
  } 
  out
}

############################################################################################################################################################
####################################################### DISPLAY BIGGEST OBJECTS IN CURRENT ENV #############################################################
############################################################################################################################################################

if (FALSE) sapply(ls(), function(x) object.size(eval(parse(text=x)))) %>% sort(., decreasing=T) %>% head(., n=10)

############################################################################################################################################################
################################################## n_cores plus #############################################################
############################################################################################################################################################

detectCores_plus <- function(Gb_max=250, 
                             additional_Gb=1) {
  # args: 
  #  Gb_max: ceiling on session memory usage in Gb, assuming that each worker duplicates the session memory
  #  additional_Gb: max additional memory requirement for new (temporary) objects created within a parallel session 
  # value:
  #   n_cores (integer)
  obj_size_Gb <- as.numeric(sum(sapply(ls(envir = .GlobalEnv), function(x) object.size(x=eval(parse(text=x))))) / 1024^3)
  max(1, min(detectCores(), Gb_max %/% (obj_size_Gb + additional_Gb))-1) 
}

############################################################################################################################################################
############################################################## safepar #####################################################################################
############################################################################################################################################################

safeParallel = function(fun, list_iterable, simplify=F, MARGIN=NULL, n_cores=NULL, Gb_max=NULL, outfile=NULL,  ...) {
  #' @usage: calls the appropriate parallel computing function, with load balancing, 
  #'          and falls back on vectorised equivalent if makeCluster hangs or the parallelised computation fails.
  #' @param fun: function to iterate over each list or vector in list_iterable
  #' @param list_iterable: a named list of vectors or lists to iterate over with fun
  #' @param simplify: whether to attempt to simplify output to a vector or matrix
  #' @param MARGIN: when iterating over a multi-dimensional object, which ones to use
  #' @param n_cores: number of FORK cores to use for parallel computation. If NULL, a safe estimate is 
  #' made based on size of objects in the global environment and the length of the iterables in list_iterable
  #' @param outfile: path to output log file, passed to the parallelising function 
  #' @param Gb_max: Gb of RAM available; if not provided, estimated
  #' @param ... : additional arguments for fun to evaluate, but not iterate over
  #' @value: named list (or if args contains "SIMPLIFY" = T, a named vector or matrix); in case of failure, NULL
  #' names will be taken from first component of list_iterable
  
  if (length(list_iterable)==1) names(list_iterable) <- "X"
  
  list_out <- NULL
  
  if (is.null(n_cores)) {
    if (is.null(Gb_max)) Gb_max=200
    additional_Gb = max(as.numeric(sapply(list_iterable, FUN = function(x) object.size(x), simplify = T)))/1024^3
    obj_size_Gb <- as.numeric(sum(sapply(ls(envir = .GlobalEnv), function(x) object.size(x=eval(parse(text=x)))))) / 1024^3
    n_cores <- min(max(sapply(list_iterable, length)), min(detectCores()%/%3, Gb_max %/% (obj_size_Gb + additional_Gb))-1)
  }
  
  if (n_cores >=2) {
    cl <-  if (!is.null(outfile)) try(makeCluster(spec=max(1,n_cores), type="FORK", timeout=30, outfile = outfile)) else try(makeCluster(spec=max(1,n_cores), type="FORK", timeout=30))
  } else {
    cl <- "none"
    class(cl) <- "try-error"
  }
  
  list_args <- list_iterable
  
  if (!"try-error" %in% class(cl)) {
    if (length(list_iterable)>1) {
      fnc <- "clusterMap"
      list_args[[".scheduling"]] = c("dynamic")
      list_args[["SIMPLIFY"]] <- simplify
    } else {
      fnc <- "parLapplyLB"
      if (simplify) {
        fnc <- "parSapplyLB"
        list_args[["simplify"]] <- simplify
      }
      if (!is.null(MARGIN)) {
        fnc <- "parApplyLB"
        list_args[["MARGIN"]] <- MARGIN
        list_args[["simplify"]] <- NULL
      }
    }
    list_args[["fun"]] <- fun
    list_args[["cl"]] = cl
    
  } else if ("try-error" %in% class(cl)) {
    
    if (length(list_iterable)>1) {
      
      fnc = "mapply"
      list_args[["SIMPLIFY"]] <- simplify
      
    } else {
      
      fnc <- "lapply"
      if (simplify) {
        fnc <- "sapply"
        list_args[["simplify"]] = simplify
      }
      if (!is.null(MARGIN)) {
        fnc <- "apply"
        list_args[["MARGIN"]] <- MARGIN
        list_args[["simplify"]] <- NULL
      }
    }
    
    list_args[["FUN"]] <- fun
    
  }
  
  # pass ... through to the parallel or vectorising function
  if (length(list(...))) for (name in names(list(...))) list_args[[name]] <- list(...)[[name]]
  
  list_out <- tryCatch({ 
    do.call(what=fnc, args = list_args)
  }, error = function(err) {
    invisible(gc())
    
    if (fnc=="clusterMap") {
      fnc <- "mapply" 
    } else if (fnc=="parLapplyLB") {
      fnc <- "lapply"
    } else if (fnc=="parSapplyLB") {
      fnc <- "sapply"
      list_args[["SIMPLIFY"]] <- NULL
      list_args[["simplify"]] <- simplify
    } else if (fnc=="parApplyLB") {
      fnc <- "apply"
    }
    
    list_args[["cl"]] <- list_args[["fun"]] <- list_args[[".scheduling"]] <- NULL
    list_args[["FUN"]] <- fun
    do.call(what=fnc, args = list_args)
  })
  
  if (!"try-error" %in% class(cl)) try(stopCluster(cl))
  invisible(gc())
  
  names(list_out) <- names(list_iterable[[1]])
  
  list_out 
}

############################################################################################################################################################
################################################################ verboseFnc ################################################################################
############################################################################################################################################################

verboseFnc <- function(fnc,
                   args) {
  #' @usage Function wrapper to show head of inputs and outputs
  #' @param fnc: a function object
  #' @param args: a list of named arguments
  #' @value value returned by fnc, if any   
  
  # TODO: fix output to outfile. Currently NULL
  # TODO: deparse(substritute(fnc)) just prints fnc
  
  message(paste0("FUNCTION CALL: ", deparse(substitute(fnc))))
  
  #if (!is.null(path_outFile)) cat(text = paste0("\nFUNCTION CALL: ", deparse(substitute(fnc))), file=path_outFile, append=T, sep="\n")
  
  message("ARGUMENTS:")
  #if (!is.null(path_outFile)) cat(text = "\nARGUMENTS:", file=path_outFile, append=T, sep="\n")

  # first print head / str / first cols and rows of args
  for (i in 1:length(args)) {
    anArg = args[[i]]
    print(paste0(names(args)[i]))
    #if (!is.null(path_outFile)) cat(text = paste0("\n",names(args)[i]) , file =  path_outFile, append=T, sep="\n")
     if (any(class(anArg) %in% c("data.frame", "data.table", "matrix", "Matrix", "data.frame", "dgCMatrix"))) {
       #if(ncol(anArg)>8) print(anArg[1:5,1:(min(nrow(anArg), 8))]) else head(anArg, n=5)
       print(anArg[1:3,1:(min(nrow(anArg), 5))])
       #if (!is.null(path_outFile)) cat(text = anArg[1:3,1:(min(nrow(anArg), 5))] , file =  path_outFile, append=T, sep="\n")
       
     } else if (any(c("standardGeneric", "function") %in% class(anArg))) {
       print(deparse(anArg, nlines=3))
      # TODO: can we concatenate something printed?
      #if (!is.null(path_outFile)) {
    #  cat(text = format(str(anArg, vec.len=4, nchar.max = 40, list.len=3)), file =  path_outFile, append=T, sep="\n") 
      } else {
        str(anArg, vec.len=4, nchar.max = 40, list.len=3)
      }
    } 
  

  # call fnc 
  out <- do.call(what=fnc, args = args)

  # print value
  message("VALUE:") 

  #if (!is.null(path_outFile)) cat(text = "\nVALUE: ", file=path_outFile, append=T, sep="\n")   

  if (any(class(out) %in% c("data.frame", "data.table", "matrix", "Matrix", "data.frame", "dgCMatrix"))) {
    #if(ncol(out)>8) print(out[1:5,1:(min(nrow(out), 8))]) else head(out, n=5)
    #if (!is.null(path_outFile)) cat(text = out[1:3,1:(min(nrow(out), 5))], file =  path_outFile, append=T, sep="\n") 
    print(out[1:3,1:(min(nrow(out), 5))]) 
  } else if (any(c("standardGeneric", "function") %in% class(out))) {
    print(deparse(out, nlines=3))
    # TODO: can we concatenate something printed?
    #if (!is.null(path_outFile)) {
    #  cat(text = format(str(anArg, vec.len=4, nchar.max = 40, list.len=3)), file =  path_outFile, append=T, sep="\n") 
  } else {
    #if (!is.null(path_outFile)) { cat(text = format(str(out, vec.len=4, nchar.max = 40, list.len=4)), file =  path_outFile, append=T, sep="\n") 
    #  } else { 
    str(out, vec.len=4, nchar.max = 40, list.len=4)
   #   }
  } 
  
  # print any warnings
  if (!is.null(warnings())) {
    #if (doPrint) {
      message("\nWARNINGS:")
      print(format(warnings()))
    #}
    #if (!is.null(path_outFile)) {
    #  cat(text="WARNINGS: ", file=path_outFile, append=T, sep="\n")
    #  cat(text = format(warnings()), file=path_outFile, append=T, sep="\n")  
    #}
  }
  # return
  return(out)

}

############################################################################################################################################################
############################################################## saveMeta ####################################################################################
############################################################################################################################################################

# TODO: add pander pandoc markdown output option?
#require(pander) 

saveMeta <- function(savefnc=NULL, doPrint=F, path_log=NULL, ...) {
  #' @usage If savefnc provided, write a timestamped metadata file along with file; else write log anyway
  #'        If path_log is not given, i.e. NULL, makes a file in the working directory
  #'        Optionally print to screen
  #'        Includes
  #'         filename (if savefnc is not NULL) or else 
  #'         current date and time
  #'         sessionInfo() #devtools::session_info()
  #'         github log, if current or parent dir is a github dir
  #' @param savefnc a function to write some file to disk. If given as object is converted to character, default NULL
  #' @param doPrint print output to screen? useful if directin Rscript stdout to a log file (&>), default F
  #' @param path_log specify log file path; defaults to creating a file in current working dir, default NULL
  #' @param ... arguments to pass on to savefnc   
  #' @value NULL
  #' @examples 
  #' saveMeta(savefnc= ggsave, plot=p, filename= "/projects/myname/RNAplot.pdf", height=10,width=12)
  #' saveMeta(savefnc= save.image,  file= "/projects/myname/session.image.Rdata.gz", compress="gzip")
  #' saveMeta(doPrint=T)
  
  # packages
  require(magrittr)
  require(utils) 
  require(devtools) # for devtools::session_info()
  
  # check args
  if (is.null(savefnc) & is.null(path_log) & !doPrint) stop("saveMeta: savefnc and path_log are NULL and doPrint is FALSE, no output")
    
  args = list(...)
  
  gsub("\\:", ".", gsub("\\ ", "_",as.character(Sys.time())))-> flag_date
  
  if (!is.null(savefnc)) {
    # convert savefnc arg to character
    if (!"character" %in% class(savefnc)) fnc <- as.character(quote(savefnc))
    
    # call savefnc to save the file
    do.call(what=savefnc, args=args, quote = F)
    
    # infer the path or file argument given to savefnc and make a modified path for the log file
    name_path_log <- sapply(c("con", "file", "path"), function(str) grep(str, names(args), value=T))
    name_path_log <- Filter(x=name_path_log,f=length)[[1]]
    path_log <- if ("character"%in%class(args[[name_path_log]])) {
      args[[name_path_log]]
    } else { #if not given as a character capture it otherwise
      args[[name_path_log]] %>% capture.output %>% '['(2) %>% gsub('description|\\ |\\"',"",.)
    }
    # replace the file extention suffix with datestamp and .txt file extention
    path_log <- gsub("\\.[a-z]+$", paste0("_meta_", substr(flag_date,1,10),".txt"),path_log, ignore.case = T)
  }  

  # If savefnc is NULL and no path_log provided, save to a file in the working dir
  if (is.null(path_log)) path_log <- paste0(getwd(), "log_meta_", substr(flag_date,1,10),".txt")
  
  # Write to log file 
  # session_info()
  utils::capture.output(devtools::session_info()) %>% writeLines(text=.,con = path_log) #this has to come first since writeLines doesn't append
  #utils::capture.output(devtools::session_info()) %>% pander %>% cat(., file=path_log)
  cat("\n", file=path_log, append=T)
  cat(text="-Environment-----------------------------------", file=path_log, sep = "\n", append=T)
  capture.output(ls.str(envir = .GlobalEnv)) %>% cat(... = ., file=path_log, sep="\n", append=T)
  
  # git commit: check parent directory for .git file
  path_parent <- gsub("[^/]*$","", path_log) 
  path_parent <- substr(x=path_parent,1,nchar(path_parent)-1)
  path_parent <- gsub("[^/]*$","", path_parent) 
  
  if (length(dir(path = path_parent, pattern="\\.git", all.files = T, recursive = F))>1) {
    try({
    gitCommitEntry <- try(system2(command="git", args=c("log", "-n 1 --oneline"), stdout=F, stderr = F))
      if (!"try-error" %in% class(gitCommitEntry)) {
        cat("\n", file=path_log, append=T)
        cat(text = "-Git commit------------------------------------", file =  path_log,  sep = "\n", append=T)
        cat(text = gitCommitEntry, file =  path_log, append=T, sep = "\n")
      }
    })
  } else {
    gitCommitEntry <- NULL
  }
  
  #date 
  cat("\n", file=path_log, append=T)
  cat(text="-Date, time------------------------------------", file=path_log, sep = "\n", append=T)
  cat(text = flag_date, file = path_log, append=T, sep = "\n")
  cat(text="-----------------------------------------------", file=path_log, sep = "\n", append=T)
  
  
  if (doPrint) {
    # same as above, but print
    
    #date 
    message(x = flag_date, appendLF = T)
    
    devtools::session_info() %>% print
    
    message("")
    print(x="-Environment-----------------------------------")
    ls.str(envir = .GlobalEnv) %>% print 
  
    if (!"try-error" %in% class(gitCommitEntry) & !is.null(gitCommitEntry)) {
      message("")
      print(x = "-Git commit------------------------------------")
      print(x = gitCommitEntry)
    }
  }
}
