###### Axel Rosendahl Huber
# Aim: Generate Mutational Matrix for Blood samples 
# This function finds vcf_files and splits the filename to get the sample name
# Returns a list of filenames (paths), and sample names
get_vcfs <- function(directory, pattern, recursive = T) { 
  vcf_files <- dir(path = directory, pattern = ".vcf", full.names = T, recursive = recursive)
  sample_names <- list()
  sample_names <- c()
  count <- 0
  for (i in vcf_files){  #Create sample names
    count <- count + 1
    id <- strsplit(i, '/', fixed = T)
    id <- tail(unlist(id), 1)
    id <- unlist(strsplit(tail(id, 1), '_', fixed = T))[1]
    sample_names <- c(sample_names, id)
    }
  output <- data.frame(sample_names = sample_names, vcf_files = vcf_files)
  return(output)
}
