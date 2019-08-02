# Primarily used by tutorials, but also from dmc_hierarchical
load_model <- function (dir_name,model_name) {
  temp_wd <- getwd ()
  tmp<-try({setwd(file.path (temp_wd, "dmc/models",dir_name))},silent = TRUE)
  if (class(tmp)=="try-error")
  {print(paste("There is no model folder called",dir_name,"in dmc/models"))
    }else
    {tmp2<- try({source(model_name)}, silent = TRUE)
            if (class(tmp2)=="try-error")
            {print(paste("There is no model called ",model_name," in dmc/models/",dir_name,sep=""))}
    else
    { tmp3<-try({source("dists.R")})
    if (class(tmp3)=="try-error")
    {print(paste("There is no file dists.R in in dmc/models",dir_name,dir_name,sep="/"))}
      }
    }
 
  setwd(temp_wd)
}

download_file <- function (source, dest) {
  # TODO: Find an improvement to download.file 
  # (it creates an empty doc if it can't find the original)
  download.file (source, dest)
}

load_data <- function (file, local=TRUE, attempt_download=TRUE)
{
  local_data_dir <- "tutorial/data/"
  remote_data_dir <- "ftp://35.160.204.70/"
  path <- paste (local_data_dir,file,sep="")
  
#  cat ("Attempting to load file ", filename)
  if (local) {
    if (!file.exists(path)) {
#      cat ("File not found at ", path)
      if (attempt_download) {
#        cat ("Downloading from server")
        download_file (paste(remote_data_dir, file, sep=""),
                       paste(local_data_dir,  file, sep=""))
      } else { 
        warning (paste ("Could not load file ", path)) 
      }
    }
  } else {
    download_file (paste(remote_data_dir, file, sep=""),
                   paste(local_data_dir,  file, sep=""))
  }
  
  # Load the local copy
  tmp <- load (path, envir=.GlobalEnv)   
  invisible (tmp)
}


save_data <- function (..., file)
{
  local_data_dir <- "tutorial/data/"  # For local files  (Always) 
  
  save (..., file=paste(local_data_dir,file,sep=""))
}