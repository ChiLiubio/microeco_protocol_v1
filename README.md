
# A comprehensive protocol for statistical analysis and visualization of microbiome omics data with R microeco and file2meco packages



## Requirements

The minimum required R version for this protocol is 4.4.0, and it is recommended to install the latest version.
Windows, macOS, or Linux systems can all be used to run the process.


## File description

"Rcode": each script for the workflow

"Input": containing all the required input files for the protocol. These omics data come from a published study (https://doi.org/10.1007/s11104-024-06847-9).

"extdata": other files


## Working directory

First, create a folder to run this protocol as R working directory.
Download the compressed package from this GitHub page (click on the Code button) and extract the contents of the archive.
Then move all the downloaded folders (including Rcode, Input and extdata) to this working directory.
Run the following command to check whether the working directory setting are correct.

```r
source("Rcode/check_working_dir.R")
```


## Install packages

Please first open the script "Step1_Install_packages.R" using a text editor or Rstudio,
run the commands line by line to install the required R packages. 
The user should not run all commands at once but should run them one line at a time in an interactive R environment, 
because some packages have interactive prompts during installation.
The latter half of the code in this step is to load the packages, thereby testing whether each package has been correctly installed.  
If it is failed to download a package with the message related with time (e.g., 'Error in download.file; Timeout of 60 seconds was reached'),
please try to increase the allowable download time with the following code:

```r
options(timeout = 1000)
```

or change the resource download mirror. First, find the nearest mirror to you at "https://cran.r-project.org/mirrors.html", 
and then modify the CRAN and Bioconductor mirrors in R with the following example code：

```r
mirror = "https://mirrors.tuna.tsinghua.edu.cn"
options("repos" = c(CRAN = paste0(mirror, "/CRAN/")))
options(BioC_mirror = paste0(mirror, "/bioconductor/"))
```

Installing these dependency packages is greatly related to internet speed, hence we recommend users to select mirrors closer to their location.
It usually takes about half an hour to install these packages in a desktop computer.


## Demo

All explanations of the data and expected output results have been added as annotations in the corresponding sections of the code.
It will take approximately six hours to complete all the steps with the example data. 
Users can also choose to run the steps that interest them.
We strongly recommend that users take some time to understand the purpose of each step in the codes and 
appropriately review the help documents of functions in the corresponding class of microeco package. 









