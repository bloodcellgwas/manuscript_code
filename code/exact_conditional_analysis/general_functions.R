#!/software/R-3.3.0/bin/Rscript
OUT_DIR <- Sys.getenv('OUT_DIR')
get_info_score <- function(variants) {
  library(data.table)
  all_min <- fread(paste0(OUT_DIR, "/all_min.tsv"), select=c("VARIANT", "INFO"));
  resp <- merge(data.frame(VARIANT=variants), all_min, by="VARIANT", all.x=TRUE)
  resp <- resp[match(variants, resp$VARIANT),]
  return(resp$INFO)
}

user_friendly_file_name <- function(traits) {
  resp <- c()
  for (trait in traits) {
    resp <- c(resp,
              switch(trait,
                     mcv="mcv", hgb="hgb", hct="hct", mch="mch", baso="baso", wbc="wbc", eo="eo", mrv="mrv",
                     lymph="lymph", mpv="mpv", pdw="pdw", neut="neut", "lymph_p"="lymph_p", rbc="rbc",
                     ret="ret", eo_p="eo_p", mono="mono", plt="plt", rdw_cv="rdw", pct="pct", hlr="hlsr",
                     hlr_p="hlsr_p", irf="irf", mscv="mscv", neut_p="neut_p", ret_p="ret_p", mono_p="mono_P",
                     mchc="mchc", baso_p="baso_p", stop("error")))
  }
  return(resp)
}

user_friendly_trait_name <- function(traits) {
  resp <- c()
  for (trait in traits) {
    resp <- c(resp,
              switch(trait,
                     mcv="MCV", hgb="HGB", hct="HCT", mch="MCH", baso="BASO#", wbc="WBC#", eo="EO#", mrv="MRV",
                     lymph="LYMPH#", mpv="MPV", pdw="PDW", neut="NEUT#", "lymph_p"="LYMPH%", rbc="RBC#",
                     ret="RET#", eo_p="EO%", mono="MONO#", plt="PLT#", rdw_cv="RDW", pct="PCT", hlr="HLSR#",
                     hlr_p="HLSR%", irf="IRF", mscv="MSCV", neut_p="NEUT%", ret_p="RET%", mono_p="MONO%",
                     mchc="MCHC", baso_p="BASO%", stop("error")))
  }
  return(resp)
}

get_trait_cell_type <- function(traits, user_friendly=FALSE) {
  resp <- c()
  for (trait in traits) {
    if (trait %in% c("plt", "mpv", "pdw", "pct")) { resp <- c(resp, "platelet") }
    else if (trait %in% c("rbc", "mcv", "hct", "mch", "mchc", "hgb", "rdw_cv", "mscv")) { resp <- c(resp, "red_cell") }
    else if (trait %in% c("ret", "ret_p", "irf", "hlr", "hlr_p", "mrv")) { resp <- c(resp, "red_cell") }
    else if (trait %in% c("mono", "mono_p")) { resp <- c(resp, "monocyte") }
    else if (trait %in% c("neut_p", "neut")) { resp <- c(resp, "neutrophil") }
    else if (trait %in% c("eo", "eo_p")) { resp <- c(resp, "eosinophil") }
    else if (trait %in% c("baso", "baso_p")) { resp <- c(resp, "basophil") }
    else if (trait %in% c("lymph", "lymph_p")) { resp <- c(resp, "lymphocyte") }
    else if (trait %in% c("wbc")) { resp <- c(resp, "white_cell") }
  }
  if (user_friendly==TRUE) {
    resp <- tools::toTitleCase(gsub("_", " ", resp))
  }
  return(resp)
}
get_trait_cell_type_order = function(traits, user_friendly=FALSE) {
    resp = c()
    for (trait in get_trait_cell_type(traits)) {
    if (trait %in% c("plt", "mpv", "pdw", "pct")) { resp <- c(resp, 1) }
    else if (trait %in% c("rbc", "mcv", "hct", "mch", "mchc", "hgb", "rdw_cv", "mscv")) { resp <- c(resp, 2) }
    else if (trait %in% c("ret", "ret_p", "irf", "hlr", "hlr_p", "mrv")) { resp <- c(resp, 2) }
    else if (trait %in% c("mono", "mono_p")) { resp <- c(resp, "monocyte") }
    else if (trait %in% c("neut_p", "neut")) { resp <- c(resp, "neutrophil") }
    else if (trait %in% c("eo", "eo_p")) { resp <- c(resp, "eosinophil") }
    else if (trait %in% c("baso", "baso_p")) { resp <- c(resp, "basophil") }
    else if (trait %in% c("lymph", "lymph_p")) { resp <- c(resp, "lymphocyte") }
    else if (trait %in% c("wbc")) { resp <- c(resp, "white_cell") }        
    }
    return(resp)

}
get_traits_in_eqtl = function(traits) {
    resp = c()
    for (trait in traits) {
        if (is.null(get_trait_eqtl_type(trait))) {
            resp = c(resp, FALSE)
        } else {
            resp = c(resp, TRUE)
        }
    }
    return(resp)
}
get_trait_eqtl_type = function(trait) {
  if (trait %in% c("plt", "mpv", "pdw", "pct")) { resp <- c("PLA") }
  else if (trait %in% c("mono", "mono_p")) { resp <- c("CD14") }
  else if (trait %in% c("neut_p", "neut")) { resp <- c("CD15") }
  else if (trait %in% c("lymph", "lymph_p")) { resp <- c("CD4", "CD8", "CD14", "CD19") }
  else {resp = c()}
  return(resp)
}

h <- function(w) if( any( grepl("running command \\'grep \\-r \\'.* had status 1", w) ) ) invokeRestart( "muffleWarning" )
h1 <- function(w) if( any( grepl("invalid factor level, NA generated", w) ) ) invokeRestart( "muffleWarning" )
# split output of grep search, required some postprocessing especially if
# multiple SNPs exist at that basepair position
process_grep_output <- function(header, output, variant) {
  output <- strsplit(output, "\t")
  snp_output <- ""
  if (length(output) == 1) {
    snp_output <- output[[1]]
    names(snp_output) <- header
    snp_output <- as.list(snp_output)
    return(snp_output)
  } else if (length(output) > 1) {
    # loop through all the outputs and find the one which is relevant
    for (item in 1:length(output)) {
      this_snp_output <- output[[item]]
      names(this_snp_output) <- header
      this_snp_output <- as.list(this_snp_output)
      bp <- as.integer(gsub(".+:(.+)_.+_.+","\\1", variant))
      ref <- gsub(".+:.+_(.+)_.+","\\1", variant)
      alt <- gsub(".+:.+_.+_(.+)","\\1", variant)
      chr <- gsub("(.+):.+_.+_.+","\\1", variant)
      if (this_snp_output$CHR == chr && this_snp_output$BP == bp
          && this_snp_output$ALLELE1 == ref && this_snp_output$ALLELE0 == alt) {
        snp_output <- this_snp_output
        return(snp_output)
      }
    }
  } else {
    stop(paste("for variant", variant, "grep output is empty", output))
  }
}

# retrieve the BETA, SE, PVAL and LOGP for this trait
get_variant_trait_boltlmm <- function(variant, trait) {
   chr <- gsub("(.+):.+_.+_.+","\\1", variant)
   bp <- as.integer(gsub(".+:(.+)_.+_.+","\\1", variant))
   filename <- paste(Sys.getenv('BOLT_OUT_DIR'), "/", trait, "_gwas_normalised_imputed.out", sep="")

   # get headers first (this is inefficient but I want my code to adapt if BOLT-LMM changes order of cols)
   cmd <- paste("head -n1", filename)
   output <- system(cmd, intern=TRUE)
   header <- strsplit(output, "\t")[[1]]
   cmd <- paste("grep -r '", bp, "' ", filename, sep="")
   withCallingHandlers(output <- system(cmd, intern=TRUE), warning=h)
   split_output <- process_grep_output(header, output, variant)

   # check output to ensure it's what we expect
   ref <- gsub(".+:.+_(.+)_.+","\\1", variant)
   alt <- gsub(".+:.+_.+_(.+)","\\1", variant)
   if (split_output$CHR != chr || split_output$BP != bp
       || split_output$ALLELE1 != ref || split_output$ALLELE0 != alt)
   {
     stop(paste("Error in get_variant_trait_boltlmm(), with", variant, "for trait,", trait, "ref is", split_output$ALLELE0, "alt is", split_output$ALLELE1, "output length", length(split_output)))
   }
   
   # flip the beta, because BOLT-LMM was given the ref and alt in reverse
   split_output$BETA <- as.numeric(split_output$BETA) * -1
   
   # calculate MLOG10P
   split_output$MLOG10P=-pchisq((split_output$BETA/as.numeric(split_output$SE))^2,df=1,lower.tail=FALSE,log.p=TRUE)/log(10)
   # calculate MAF
   split_output$ALT_FREQ=1-as.numeric(split_output$A1FREQ)
   split_output$ALT_MINOR <- as.numeric(split_output$ALT_FREQ)<0.5
   split_output$MA_FREQ <- as.numeric(split_output$ALT_FREQ)
   split_output$MA_FREQ[!split_output$ALT_MINOR] <- 1 - split_output$MA_FREQ[!split_output$ALT_MINOR]
   out_list <- list(variant, as.numeric(split_output$BETA), as.numeric(split_output$SE), as.numeric(split_output$P_BOLT_LMM_INF), as.numeric(split_output$MLOG10P), as.numeric(split_output$MA_FREQ))
   names(out_list) <- c("VARIANT", "BETA", "SE", "P", "MLOG10P", "MAF") 
   # get a df of the variant, beta, se and p value for this variant
   return(as.data.frame(out_list, stringsAsFactors=F))
}

get_rsid_from_variantid <- function(variantids, mapping_file=NULL) {
    variantids = gsub("XY:", "X:", variantids)
  if (is.null(mapping_file)) {
    mapping_file <- fread(Sys.getenv('MAPPING_FILE_RSID'))
    if (sum(is.element(variantids, mapping_file$COORDID)) != length(variantids)) {
      mapping_file <- fread(mapping_file)
    }
  }
  return(mapping_file$dbSNP[match(variantids, mapping_file$COORDID)])
}

get_variantid_from_rsid <- function(rsids, mapping_file=NULL) {
  library(data.table)
  if (is.null(mapping_file)) {
    mapping_file <- fread(Sys.getenv('MAPPING_FILE_RSID'))
    if (sum(is.element(rsids, mapping_file$dbSNP_49)) != length(rsids)) {
      mapping_file <- fread(paste(Sys.getenv('SUPPORT_FILES'), "/id_mapping_table_rsids.tsv", sep=""), drop=c("dbSNP_47"))
    }
  }
  return(mapping_file$COORDID[match(rsids, mapping_file$dbSNP_49)])
}

user_friendly_disease_name_short = function(diseases) {
  resp <- c()
  for (disease in diseases) {
      if (disease %in% c( "IBS_delange_2017", "IBS_IC_liu_2015", "IBS_LIU_2015")) {
          resp = c(resp, "IBD")
      } else if (disease %in% c("IBS_CD_delange_2017", "IBS_CD_IC_liu_2015", "IBS_CD_liu_2015")) {
          resp = c(resp, "Crohn's")
      } else if (disease %in% c("IBS_UC_delange_2017", "IBS_UC_IC_liu_2015", "IBS_UC_liu_2015")) {
          resp = c(resp, "Ulcerative Colitis")
      } else if (disease %in% c("hayfever_or_rhinitis")) {
          resp = c(resp, "Hayfever/Rhinitis")
      } else if (disease == "alzheimers_lambert_2013") {
          resp = c(resp, "ALZ")
      } else if (disease == "cad_nikpay_2015") {
          resp = c(resp, "CAD")
      } else if (disease %in% c("celiac_disease_dubois_2010", "celiac_disease_IC_trynka_2011")) {
          resp = c(resp, "Celiac")
      } else if (disease %in% c("multiple_sclerosis_patsopoulos_2017", "multiple_sclerosis_IC_beecham_2013", "multiple_sclerosis_sawcer_2011")) {
          resp = c(resp, "MS")
      } else if (disease %in% c("primary_biliary_cirrhosis_cordell_2015", "primary_biliary_cirrhosis_IC_liu_2012")) {
          resp = c(resp, "Primary Biliary Cirrhosis")
      } else if (disease %in% c("primary_sclerosing_cholangitis_ji_2016")) {
          resp = c(resp, "Primary Sclerosing Cholangitis")
      } else if (disease %in% c("rheumatoid_arthritis_okada_2014")) {
          resp = c(resp, "Rheumatoid Arthritis")
      } else if (disease %in% c("systemic_lupus_erythematosus_bentham_2015")) {
          resp = c(resp, "Lupus")
      } else if (disease %in% c("type_1_diabetes_IC_gumuscu_2015", "type_1_diabetes_meta_IC_gumuscu_2015")) {
          resp = c(resp, "T1D")
      } else if (disease %in% c("allergic_disease_EUR_ferreira_2017")) {
          resp = c(resp, "Allergic Disease")
      } else if (disease %in% c("asthma_EUR_tagc_2018")) {
          resp = c(resp, "Asthma")
      } else if (disease %in% c("eczema_eagle_2015")) {
          resp = c(resp, "Eczema")
      } else if (disease == "neut") {
          resp = c(resp, "NEUT#")
      }
  }
  return(resp)
}
convert_disease_name <- function(disease_names) {
  converted_names = c()
  disease_table=read.csv(Sys.getenv('DISEASE_NAMES'), header=T, stringsAsFactors=FALSE)
  for (disease in disease_names) {
    disease_category = disease_table[disease_table$full_name == disease,]
    if (nrow(disease_category) > 1) {stop("duplicates in disease table")}
    converted_names = c(converted_names,  disease_category$category)
  }
  return(converted_names)
}
theme_publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold",size = rel(1.5)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(size = rel(1)),
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "bottom",
           legend.direction = "horizontal",
           legend.key.size= unit(rel(1), "cm"),
           legend.margin = unit(0, "cm"),
           legend.title = element_text(face="italic"),
           legend.text=element_text(size=rel(1.2)),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#f0f0f0", fill="#f0f0f0"),
           strip.text = element_text(face="bold")
   ))
  
}
internal_get_trait_cell_type_order = function(traits) {
  cord = c()
  for (trait in traits) {
    cord = c(cord, switch(get_trait_cell_type(trait),
      "platelet"=1,
      "red_cell"=2,
      "white_cell"=3,
      "neutrophil"=4,
      "eosinophil"=5,
      "basophil"=6,
      "monocyte"=7,
      "lymphocyte"=8))
  }
  return(cord)
}
get_trait_cell_type_order = function(traits) {
    t = data.frame(traits=traits, ord = internal_get_trait_cell_type_order(traits), stringsAsFactors=F)
    return(t[order(-t$ord),]$traits)
    
}

alias_to_symbol = function(ids) {
  library(limma)
  resp = c()
  for (id in ids) {
    hgnc = alias2Symbol(id)
    if (length(hgnc) == 0) {
      resp = c(resp, NA)
    } else if (length(hgnc) == 1) {
      resp = c(resp, hgnc)
    } else {
      stop(paste(id, 'more than 1 ids alias_to_symbol():', hgnc))
    }
  }
  return(resp)
}
