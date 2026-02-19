# fit=
# Rsq_min = NULL 
# sft_value = 10
# minModuleSize = 50
# maxBlockSize = 500
# deepSplit = 4
# networkType = "signed" 
# TOMType = "signed"
# mods_mean = FALSE
# mods_eigen = FALSE
# 
# david = FALSE
# nThread = 8
# Rsq.min = NULL
# sft.value = NULL
# mods.mean = FALSE
# mods.eigen = FALSE

make_modules<- function (fit, Rsq_min = NULL, sft_value = NULL, minModuleSize = 20, 
          maxBlockSize = 500, deepSplit = 3, networkType = "signed", 
          TOMType = "signed", mods_mean = FALSE, mods_eigen = FALSE, 
          david = FALSE, nThread = 2, Rsq.min = NULL, sft.value = NULL, 
          mods.mean = FALSE, mods.eigen = FALSE) 
{
  SFT.R.sq <- Power <- module <- geneName <- module.char <- NULL
  if (!is.null(Rsq.min)) {
    Rsq_min <- Rsq.min
  }
  if (!is.null(sft.value)) {
    sft_value <- sft.value
  }
  if (mods.mean) {
    mods_mean <- mods.mean
  }
  if (mods.eigen) {
    mods_eigen <- mods.eigen
  }
  if (!is.null(Rsq_min)) {
    sft.select <- dplyr::filter(fit$sft, SFT.R.sq >= Rsq_min)
    if (nrow(sft.select) > 0) {
      power.t <- min(sft.select$Power)
    }
    else {
      stop("R-squared minimum not reached. Please input lower Rsq_min or set sft_value instead.")
    }
  } else if (!is.null(sft_value)) {
    sft.select <- dplyr::filter(fit$sft, Power == sft_value)
    power.t <- unique(sft.select$Power)
  } else {
    stop("Please set Rsq_min or sft_value.")
  }
  
  cor <- WGCNA::cor
  mod.net <- WGCNA::blockwiseModules(t(fit$dat$E), power = power.t, 
                                     networkType = networkType, TOMType = TOMType, maxBlockSize = maxBlockSize, 
                                     minModuleSize = minModuleSize, deepSplit = deepSplit, 
                                     numericLabels = TRUE, saveTOMFileBase = "TOM-blockwise", 
                                     nthreads = nThread)
  cor <- stats::cor
  mods <- as.data.frame(mod.net$colors) %>% tibble::rownames_to_column("geneName") %>% 
    dplyr::rename(module = "mod.net$colors") %>% dplyr::mutate(module.char = ifelse(module <= 
                                                                                      9, paste("0", module, sep = ""), module)) %>% dplyr::mutate(mod.color = WGCNA::labels2colors(mod.net$colors))
  if (!is.null(fit$dat$genes)) {
    mods <- mods %>% dplyr::left_join(fit$dat$genes, by = "geneName")
  }
  dat.mods <- list()
  dat.mods[["genes"]] <- fit$genes
  dat.mods[["mods"]] <- mods
  dat.mods[["sft"]] <- sft.select
  dat.mods[["top.plot"]] <- fit$top.plot + ggplot2::geom_hline(ggplot2::aes(yintercept = sft.select$SFT.R.sq[1]), 
                                                               color = "red") + ggplot2::geom_vline(xintercept = power.t, 
                                                                                                    color = "red")
  dat.mods[["connect.plot"]] <- fit$connect.plot + ggplot2::geom_hline(ggplot2::aes(yintercept = sft.select$mean.k.[1]), 
                                                                       color = "red") + ggplot2::geom_vline(xintercept = power.t, 
                                                                                                            color = "red")
  if (mods_mean) {
    mods.voom <- mods %>% dplyr::select(geneName, module.char) %>% 
      dplyr::left_join(tibble::rownames_to_column(as.data.frame(fit$dat$E), 
                                                  "geneName"), by = "geneName") %>% dplyr::group_by(module.char) %>% 
      dplyr::summarise_if(is.numeric, mean, na.rm = TRUE) %>% 
      tibble::rownames_to_column() %>% dplyr::mutate(rowname = paste("module", 
                                                                     module.char, sep = "_")) %>% tibble::column_to_rownames()
    dat.mods[["mods.mean"]] <- mods.voom
  }
  if (mods_eigen) {
    mods.E <- mod.net$MEs %>% t() %>% as.data.frame() %>% 
      tibble::rownames_to_column("module.char") %>% dplyr::mutate(module.char = as.numeric(gsub("ME", 
                                                                                                "", module.char))) %>% dplyr::mutate(module.char = ifelse(module.char <= 
                                                                                                                                                            9, paste("0", module.char, sep = ""), module.char)) %>% 
      dplyr::arrange(module.char)
    rownames(mods.E) <- paste0("module_", mods.E$module.char)
    dat.mods[["mods.eigen"]] <- mods.E
  }
  if (david) {
    mod.names <- sort(rownames(mods.voom))
    max.mod <- max(table(mod.net$colors))
    david.df <- data.frame(rowname = 1:max.mod)
    for (i in 1:length(mod.names)) {
      mod.name <- mod.names[i]
      gene.list <- mods %>% dplyr::filter(module == i - 
                                            1) %>% dplyr::select(geneName)
      add.genes <- max.mod - nrow(gene.list)
      gene.list <- c(gene.list$geneName, rep(NA, times = add.genes))
      gene.list <- as.data.frame(gene.list)
      colnames(gene.list) <- mod.name
      david.df <- david.df %>% dplyr::bind_cols(gene.list)
    }
    dat.mods[["david"]] <- david.df
  }
  return(dat.mods)
}


