library(gkmSVM)
suppressMessages(suppressWarnings(library(tools)))

do_gkmSVM <- function(pop){
  
  if(TRUE){
    genNullSeqs(paste0('..//',pop,'.bed'),
                nMaxTrials=10,xfold=1,
                genomeVersion='hg19',
                outputPosFastaFN=paste0('../../output/fastas_svm/',pop,'-positive.fa'),
                outputBedFN=paste0('../../output/fastas_svm/',pop,'-negative.bed'),
                outputNegFastaFN=paste0('../../output/fastas_svm/',pop,'-negative.fa'))
    
    gkmsvm_kernel(posfile = paste0('../../output/fastas_svm/',pop,'-positive.fa'),
                  negfile = paste0('../../output/fastas_svm/',pop,'-negative.fa'),
                  outfile = paste0('../../output/kernel_svm/', pop, ".kernel.out"))
  }
  gkmsvm_trainCV(kernelfn = paste0('../../output/kernel_svm/', pop, ".kernel.out"),
                 posfn = paste0('../../output/fastas_svm/',pop,'-positive.fa'),
                 negfn = paste0('../../output/fastas_svm/',pop,'-negative.fa'), 
                 svmfnprfx=paste0('../../output/kernel_svm/', pop),
                 outputCVpredfn=paste0('../../output/kernel_svm/', pop, ".cvPred.out"),
                 outputROCfn=paste0('../../output/kernel_svm/', pop, ".roc.out"))
  
  gkmsvm_classify('../../output/nr10mers.fa',
                  svmfnprfx=paste0('../../output/kernel_svm/', pop),
                  paste0('../../output/kernel_svm/', pop, ".weights.10mer.out"))
  
  pop
  
}

args <- commandArgs(trailingOnly = TRUE)
if(file_path_sans_ext(basename(args[1])) == "R"){
  i <- 2
} else { # Rscript
  i <- 0
}
idx <- as.numeric(args[i+1])

pops <- gsub(".bed", "", list.files("../../output/processed_peaks_svm/"))

do_gkmSVM(pops[idx])
