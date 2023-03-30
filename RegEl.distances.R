#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option ( c("--input"),
                help = "input gene distance matrix" ),
    
  make_option ( c("--output"), help = "Output file name") 
  
)


parser <- OptionParser(
  usage = "%prog [options]", 
  option_list=option_list,
  description = "\nReports mean and median values of distance matrix."
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options


#********
# BEGIN *
#********

# 1. read distance matrix
distances.matrix <- read.table( file = opt$input, sep="\t", header=F)

# 2. calculate mean and median
mean <- round(mean(distances.matrix[,3]),2)
median <- round(median(distances.matrix[,3]),2) 

# 3. print mean and media
print(paste("Mean: ",mean))
print(paste("Median: ",median))
