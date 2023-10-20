
# Implimentation of parallel processing within core DGCA functions.
# Author: Max Segnitz, msegnitz@uw.edu
# Modified from DGCA  
# Zhang B, McKenzie A (2023). _DGCA: Differential Gene Correlation Analysis_. R package version 1.0.3, <https://CRAN.R-project.org/package=DGCA>.

# Started October 2023
#

# © Richard M Segnitz 2023 (modified from DGCA)
# License: This software is licensed under GNU General Public License, and may
# be modified and/or redistributed under GNU GPL version 3 or later. License details
# can be found in the accompanying this script, or at  (http://www.gnu.org/licenses/).
#


# DESCRIPTION:
# This script and the functions within attempt to accomplish several key aims. 
# 1) Improve the efficiency of calculations made within DGCA, especially when permutations are run.
# 2) Allow the user to return gene-level results used for calculating module-level statistics. 
#   This reduces redundancy and ensures that any permutation based statistics are calculated from 
#    same randomizations of the data.
#



moduleDC_par<-
  function (inputMat, design, compare, genes, labels, corr_cutoff = 0.99, 
            signType = "none", corrType = "pearson", nPerms = 50, oneSidedPVal = FALSE, 
            gene_avg_signif = 0.05, number_DC_genes = 3, dCorAvgMethod = "median", num_cores=4,
            save_gene_level=F, seed=NULL) {
    
    require(DGCA)
    
    if (!length(genes) == length(labels)) 
      stop("Genes and labels vectors must be the same length.")
    labels_names = unique(labels)
    mdc_vector = vector()
    mdc_signif = vector()
    module_size = vector()
    goc_genes = vector()
    loc_genes = vector()
    
    if(save_gene_level){
    dcPair_byMod<-list() # Add ability to save gene-level results
    dcGeneAvg_byMod<-list()}
    
    for (i in 1:length(labels_names)) {
      message(paste0("Calculating MDC for module #", i, ", which is called ", 
                     labels_names[i]))
      genes_tmp = genes[labels == labels_names[i]]
      genes_tmp = genes_tmp[genes_tmp %in% rownames(inputMat)]
      if (length(genes_tmp) == 0) 
        next
      module_size[i] = length(genes_tmp)
      inputMat_tmp = inputMat[genes_tmp, ]
      if(!is.null(seed)){set.seed(null)} # set random seed if provided.
      ddcor_res = suppressWarnings(
        ddcorAll_par(inputMat = inputMat_tmp, 
                               design = design, 
                               compare = compare, 
                               corrType = corrType, 
                               signType = signType, 
                               adjust = "none", 
                               nPerms = nPerms, 
                               corr_cutoff=corr_cutoff,
                               getDCorAvg = TRUE, 
                               dCorAvgType = "both", 
                               classify = FALSE,
                               dCorAvgMethod = dCorAvgMethod, 
                               num_cores=num_cores))
      mdc_vector[i] = ddcor_res[["total_avg_dcor"]][["total_zdiff"]]
      mdc_signif[i] = ddcor_res[["total_avg_dcor"]][["pVal"]]
      gene_avg = ddcor_res[["gene_avg_dcor"]]
      gene_avg_sig = gene_avg[gene_avg$pVal_adj < gene_avg_signif, ]
      gene_avg_goc = head(gene_avg_sig[gene_avg_sig$avgZDiff > 
                                         0, ]$Gene, number_DC_genes)
      gene_avg_loc = head(gene_avg_sig[gene_avg_sig$avgZDiff < 
                                         0, ]$Gene, number_DC_genes)
      goc_genes[i] = paste(gene_avg_goc, collapse = ", ")
      loc_genes[i] = paste(gene_avg_loc, collapse = ", ")
      
      if(save_gene_level){
      
        if(dCorAvgMethod == "median"){ddcor_res$gene_avg_dcor<-dplyr::rename(ddcor_res$gene_avg_dcor, medianZDiff=avgZDiff)}
        
      dcPair_byMod[[labels_names[i]]]<-  dplyr::mutate(ddcor_res$dcPair, module=labels_names[i], contrast_ref=compare[1], contrast_lvl=compare[2])
      dcGeneAvg_byMod[[labels_names[i]]]<-  dplyr::mutate(ddcor_res$gene_avg_dcor, module=labels_names[i], contrast_ref=compare[1], contrast_lvl=compare[2])
      
      
      }
      
    }
    res_df = data.frame(Module = labels_names, Size = module_size, 
                        MeDC = mdc_vector, pVal = mdc_signif, Top_GOC = goc_genes, 
                        Top_LOC = loc_genes)
    
    if(save_gene_level){
      return(list(moduleResults=res_df, dcPair_byMod=dcPair_byMod, dcGeneAvg_byMod=dcGeneAvg_byMod))
    }else{
    return(moduleResults=res_dfres_df)}
  }



#### 

ddcorAll_par<- function(inputMat, design, compare, inputMatB = NULL, splitSet = NULL, 
            impute = FALSE, corrType = "pearson", nPairs = "all", sortBy = "zScoreDiff", 
            adjust = "perm", nPerms = nPerms, classify = TRUE, sigThresh = 1, 
            corSigThresh = 0.05, heatmapPlot = FALSE, color_palette = NULL, 
            verbose = FALSE, plotFdr = FALSE, corr_cutoff = 0.99, signType = "none", 
            getDCorAvg = FALSE, dCorAvgType = "gene_average", dCorAvgMethod = "median", 
            oneSidedPVal = F, customize_heatmap = FALSE, heatmapClassic = FALSE, 
            corPower = 2, num_cores=4) {
    
    require(DGCA)
    
    if (!is.null(splitSet)) {
      if (!mode(splitSet) == "character") 
        stop("splitSet must be character type.\n")
    }
    if (!is.null(splitSet) & !is.null(inputMatB)) {
      stop("You cannot input both a character vector to split the first input matrix into two as well as a second matrix of identifiers.")
    }
    if (is.null(rownames(inputMat))) {
      stop("inputMat must specify identifiers as rownames.")
    }
    if (adjust == "perm" & nPerms == 0) {
      stop("If you choose permutation for p-value adjustment, then you need to generate at least one permutation sample by setting nPerms > 0.")
    }
    if (adjust != "perm" & nPerms > 0 & getDCorAvg == FALSE) {
      warning("If you are not choosing permutation for p-value adjustment or calculating the differential correlation average, then you may be wasting time by generating permutation samples. Consider setting nPerms to 0.")
    }
    if (!(is.numeric(nPairs) | nPairs == "all")) {
      stop("nPairs must either be numeric or be a character vector \"all\" to specify that all pairs should be returned.")
    }
    if (nPairs == "all" & is.null(inputMatB)) {
      nPairs = (nrow(inputMat)^2)/2 - nrow(inputMat)/2
    }
    if (nPairs == "all" & !is.null(inputMatB)) {
      nPairs = nrow(inputMat) * nrow(inputMatB)
    }
    if (!dCorAvgMethod %in% c("median", "mean")) {
      stop("The differential correlation average method chosen must be one of median or mean.")
    }
    SAF = getOption("stringsAsFactors", FALSE)
    on.exit(options(stringsAsFactors = SAF))
    options(stringsAsFactors = FALSE)
    if (!is.null(splitSet)) {
      splitSetFound = splitSet %in% rownames(inputMat)
      nPairs = (nrow(inputMat) - sum(splitSetFound)) * sum(splitSetFound)
      if (sum(splitSetFound) == 0) {
        stop("None of the splitSet identifiers were found in the rownames of the input matrix.")
      }
      else if (sum(splitSetFound) > 0) {
        splitSetRows = rownames(inputMat) %in% splitSet
        inputMatB = inputMat[splitSetRows, , drop = FALSE]
        inputMat = inputMat[!splitSetRows, , drop = FALSE]
        if (verbose) {
          message(sum(splitSetRows), " row(s) corresponding to the ", 
                  length(splitSet), " identifier(s) found in splitSet were found in the input matrix.")
        }
      }
    }
    secondMat = FALSE
    if (!is.null(splitSet) | !is.null(inputMatB)) {
      secondMat = TRUE
    }
    print("Now running initial DCors")
    ddcor_res = getDCors(inputMat = inputMat, design = design, 
                         compare = compare, inputMatB = inputMatB, impute = impute, 
                         corrType = corrType, corr_cutoff = corr_cutoff, signType = signType)
    if (nPerms > 0) {
      print("Now running  DCors permuation")
      ddcor_perm = getDCorPerm_par(inputMat = inputMat, design = design, 
                               compare = compare, inputMatB = inputMatB, impute = impute, 
                               corrType = corrType, nPerms = nPerms, corr_cutoff = corr_cutoff, 
                               signType = signType, num_cores=num_cores)
      close(pb)
    }
    if (adjust != "perm") {
      ddcor_table = dcTopPairs(dcObject = ddcor_res, nPairs = nPairs, 
                               adjust = adjust, plotFdr = plotFdr, classify = classify, 
                               compare = compare, sigThresh = sigThresh, corSigThresh = corSigThresh, 
                               verbose = verbose, secondMat = secondMat)
    }
    if (adjust == "perm") {
      ddcor_table = dcTopPairs(dcObject = ddcor_res, nPairs = nPairs, 
                               adjust = adjust, plotFdr = plotFdr, classify = classify, 
                               compare = compare, sigThresh = sigThresh, corSigThresh = corSigThresh, 
                               zScorePerm = ddcor_perm, verbose = verbose, secondMat = secondMat)
    }
    ddcor_table = ddcor_table[order(-abs(ddcor_table[, sortBy])), 
    ]
    if (heatmapPlot) {
      ddplot_plot = ddplot(dcObject = ddcor_res, color_palette = color_palette, 
                           customize_heatmap = customize_heatmap, heatmapClassic = heatmapClassic, 
                           corPower = corPower, ...)
    }
    if (getDCorAvg) {
      ZDiffs = slot(ddcor_res, "ZDiff")
      if (!dCorAvgType == "both") {
        avg_dcor = dCorAvg_mod(zDiff = ZDiffs, zDiffPerm = ddcor_perm, 
                           dCorAvgType = dCorAvgType, secondMat = secondMat, 
                           dCorAvgMethod = dCorAvgMethod)
        ddcor_res = list(dcPair = ddcor_table, avg_dcor = avg_dcor)
      }
      else {
        gene_avg_dcor = dCorAvg_mod(zDiff = ZDiffs, zDiffPerm = ddcor_perm, 
                                dCorAvgType = "gene_average", secondMat = secondMat, 
                                dCorAvgMethod = dCorAvgMethod)
        total_avg_dcor = dCorAvg_mod(zDiff = ZDiffs, zDiffPerm = ddcor_perm, 
                                 dCorAvgType = "total_average", secondMat = secondMat, 
                                 dCorAvgMethod = dCorAvgMethod)
        ddcor_res = list(dcPair = ddcor_table, gene_avg_dcor = gene_avg_dcor, 
                         total_avg_dcor = total_avg_dcor)
      }
    }
    if (!getDCorAvg) {
      ddcor_res = ddcor_table
    }
    
 
    return(ddcor_res)
  }

###

getDCorPerm_par<-
  function (inputMat, design, compare, inputMatB = NULL, impute = FALSE, 
          nPerms = 10, corrType = "pearson", corr_cutoff = 0.99, signType = "none", num_cores=num_cores) {
  
  # Load the required libraries
  require(foreach)
  require(doParallel)
  require(DGCA)
  require(doSNOW)
  
  secondMat = FALSE
  
  if (!is.null(inputMatB)) {
    zPermMat = array(dim = c(nrow(inputMat), nrow(inputMatB), 
                             nPerms))
    secondMat = TRUE
  } else {
    zPermMat = array(dim = c(nrow(inputMat), nrow(inputMat), 
                             nPerms))
  }


# Define a function that performs the permutation and returns zscores
permute_and_get_zscores <- function(i) {
  message("Calculating permutation number ", i, ".")
  inputMat_perm = inputMat[, sample(ncol(inputMat)), drop = FALSE]

  if (secondMat) {
    inputMatB_perm = inputMatB[, sample(ncol(inputMatB)), drop = FALSE]
    corMats_res = DGCA::getCors(inputMat_perm, design = design, inputMatB = inputMatB_perm, corrType = corrType, impute = impute)
  } else {
    corMats_res = DGCA::getCors(inputMat_perm, design = design, corrType = corrType, impute = impute)
  }

  dcPairs_res =
    DGCA::pairwiseDCor(corMats_res, compare, corr_cutoff = corr_cutoff,
                 secondMat = secondMat, signType = signType)

  zscores = slot(dcPairs_res, "ZDiff")

  return(zscores)
}

# Create a parallel foreach loop to perform the permutations
# Set the number of cores to use
# Initialize a parallel backend
cl <- makeCluster(num_cores)
#registerDoParallel(cl)
registerDoSNOW(cl)
pb <- txtProgressBar(max = nPerms, style = 3)

progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

print("Permuting in DCors in Parallel")
results <- foreach(i = 1:nPerms,
                   .packages=c('dplyr', 'DGCA', 'Matrix',"foreach","doParallel"),
                   .export =c('permute_and_get_zscores', 'getDCorPerm_par', 'ddcorAll_par', 'moduleDC_par', 'dCorAvg_mod'),
                   .options.snow = opts) %dopar% {
  print(paste0("Permutation: ", i))
  permute_and_get_zscores(i)
}

# Close the parallel backend
close(pb)
stopCluster(cl)

# Combine the results into a zPermMat array
zPermMat <- abind::abind(results, along=3) 

return(zPermMat)
}



###




dCorAvg_mod<-
function (zDiff, zDiffPerm, dCorAvgType, oneSidedPVal = FALSE, 
          secondMat = FALSE, dCorAvgMethod = "median") 
{
  message("Calculating the differential correlation average.")
  if (!dCorAvgMethod %in% c("median", "mean")) {
    stop("The differential correlation average method chosen must be one of median or mean.")
  }
  if (dCorAvgType == "gene_average") {
    if (!secondMat) {
      for (i in 1:dim(zDiffPerm)[3]) { 
        tmp = zDiffPerm[, , i]
        #tmp[lower.tri(tmp)] = t(tmp)[lower.tri(t(tmp))]
        tmp<-as.matrix(Matrix::forceSymmetric(tmp, uplo="U"))
        zDiffPerm[, , i] = tmp
      }
      # zDiff[lower.tri(zDiff)] = t(zDiff)[lower.tri(t(zDiff))]
      zDiff<-as.matrix(Matrix::forceSymmetric(zDiff, uplo="U"))
    }
    zdiff_medians = numeric(ncol(zDiff))
    empirical_pval = numeric(ncol(zDiff))
    nGenes = ncol(zDiff)
    
    for (i in 1:nGenes) {
      zdiffs_gene_nonself = zDiff[-i, i] # Pull all dCor Zscores for gene
      zdiff_perm_gene_nonself = zDiffPerm[-i, i, ] # Pull dCor Z scores for each permutation
      if (dCorAvgMethod == "median") {
        zdiff_medians[i] = median(zdiffs_gene_nonself)
        zdiff_perm_gene_medians = matrixStats::colMedians(zdiff_perm_gene_nonself)
      } else if (dCorAvgMethod == "mean") {
        zdiff_medians[i] = mean(zdiffs_gene_nonself)
        zdiff_perm_gene_medians = colMeans(zdiff_perm_gene_nonself)
      }
      
      empirical_pval[i] = 1 - sum(abs(zdiff_medians[i]) > 
                                    abs(zdiff_perm_gene_medians))/length(zdiff_perm_gene_medians)
    }
    
  
    avg_dcor_df = 
       data.frame(Gene = colnames(zDiff), avgZDiff = zdiff_medians, 
                             empirical_pVal = empirical_pval)
    
       avg_dcor_df$pVal_adj = p.adjust(avg_dcor_df$empirical_pVal, 
                                    method = "BH")
       avg_dcor_df = avg_dcor_df[order(abs(avg_dcor_df$avgZDiff), 
                                    decreasing = TRUE), ]
       
    
    return(avg_dcor_df)
  }
  if (dCorAvgType == "total_average") {
    if (dCorAvgMethod == "median") {
      zdiff_median = median(zDiff, na.rm = TRUE) # take median of all pairwise Z scores
      zdiff_median_permwise = apply(zDiffPerm, 3, median, # take median of all pairwise Z scores for permutations
                                    na.rm = TRUE)
    } else if (dCorAvgMethod == "mean") {
      zdiff_median = mean(zDiff, na.rm = TRUE)
      zdiff_median_permwise = apply(zDiffPerm, 3, mean, 
                                    na.rm = TRUE)
    }
    nPerms = dim(zDiffPerm)[3]
    if (oneSidedPVal) {
      if (zdiff_median > 0) {
        pVal = 1 - (sum(zdiff_median < zdiff_median_permwise)/nPerms)
      }
      if (zdiff_median < 0) {
        pVal = 1 - sum(zdiff_median > zdiff_median_permwise)/nPerms
      }
    }
    else if (!oneSidedPVal) {
      pVal = 1 - sum(abs(zdiff_median) > abs(zdiff_median_permwise))/nPerms
    }
    return(list(total_zdiff = zdiff_median, pVal = pVal))
  }
}
