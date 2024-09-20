if(packageVersion("Base") < "4.4.0"){stop("The minimum required R version for this protocol is 4.4.0! Please reinstall R with the latest version!")}
if(! "Rcode" %in% list.files()){stop("The Rcode folder is not found in the current working directory!")}
if(! "Input" %in% list.files()){stop("The Input directory is not found in the current working directory!")}
if(! "extdata" %in% list.files()){stop("The extdata directory is not found in the current working directory!")}
