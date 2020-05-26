library(gkmSVM)
suppressMessages(suppressWarnings(library(tools)))

# Function to compute delta 
do_gkmsvm_delta <- function(population,trait="FMall_v2"){
  gkmsvm_delta(seqfile_allele1 = paste0('../../output/fastas_svm/',trait,"_ref.fasta"),
               seqfile_allele2 =  paste0('../../output/fastas_svm/',trait,"_alt.fasta"),
               svmfnprfx = paste0('../../output/kernel_svm/', population),
               outfile=paste0('../../output/variant_predictions_svm/',population,'-try1_',trait,'.out'))
  
  population
}

args <- commandArgs(trailingOnly = TRUE)
if(file_path_sans_ext(basename(args[1])) == "R"){
  i <- 2
} else { # Rscript
  i <- 0
}
idx <- as.numeric(args[i+1])

pops <- gsub(".bed", "", list.files("../../output/processed_peaks_svm/"))
do_gkmsvm_delta(pops[idx])

