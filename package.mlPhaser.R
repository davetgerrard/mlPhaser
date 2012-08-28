library(roxygen2)

#setwd("C:/Users/dave/HalfStarted/mlPhaser")
setwd("C:/Users/dave/HalfStarted/mlPhaser/package")


# FIRST TIME ONLY # package.skeleton(name = "mlPhaser", code_files="C:/Users/dave/HalfStarted/mlPhaser/mlPhaser.R", force=TRUE)
# Update DESCRIPTION FILE for VERSION
# delete contents of /man 
# Set the correct dates in mlPhaser.R and DESCRIPTION
# copy latest .R code to /R

roxygenize("mlPhaser", unlink.target=T)

#edit the namespace. Remove "import(seqinr)"

# In DOS, from same directory (/package)
# R CMD check mlPhaser
#R CMD build mlPhaser	# to get tar.gz
#R CMD INSTALL --build mlPhaser	# to get windows zip






## Install the new package version
###Generate the README 

setwd("C:/Users/dave/HalfStarted/mlPhaser/testREADME")
Sweave("../mlPhaser_README")	# from sub-directory up to main mlPhaser directory
#R CMD texify --pdf mlPhaser_README.tex
Stangle("../mlPhaser_README")
##Sweave("../mlPhaser_README_0_16")	# from sub-directory up to main mlPhaser directory
###R CMD texify --pdf mlPhaser_README_0_16.tex
##Stangle("../mlPhaser_README_0_16")