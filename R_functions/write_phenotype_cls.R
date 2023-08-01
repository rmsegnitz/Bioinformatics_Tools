write_phenotype.cls<-
  function (file_prefix, counts, design, var_to_test, libID_col = "lib.id") 
  {
    counts <- countSubsetNorm::extract_counts(counts)
    if (is.matrix(design)) 
      design <- as.data.frame(design)
    if (libID_col == "row.names" & !("row.names" %in% colnames(design))) 
      design$row.names <- rownames(design)
    if (!(libID_col %in% colnames(design))) 
      stop(paste0("Design object is missing column ", libID_col, 
                  ", where I expected to find library identifiers."))
    if (!setequal(colnames(counts), design[, libID_col])) 
      stop("Library identifiers in the counts and design objects do not match. Check that the objects contain the same libraries.")
    if (any(duplicated(colnames(counts))) | any(duplicated(design[, 
                                                                  libID_col]))) 
      stop("Duplicated library identifiers in counts or design object.")
    sink(paste(file_prefix, var_to_test, "cls", sep = "."))
    if (is.numeric(design[, var_to_test])) {
      cat("#numeric\n")
      cat("#", var_to_test, "\n", sep = "")
      cat(design[match(colnames(counts), design[, libID_col]), 
                 var_to_test])
      cat("\n")
    }
    else {

      ref_level<-design%>%pull(var_to_test)%>%head(1)%>%as.character()
      
      design[, var_to_test]<-as.factor(as.character(design[, var_to_test]))
      design[, var_to_test]<-relevel(design[, var_to_test], ref=ref_level)
      
      cat(ncol(counts), length(levels(design[, var_to_test])), 
          "1\n", sep = " ")
      cat("#", levels(design[, var_to_test]), sep = " ")
      cat("\n")
      # cat(as.numeric(design[match(colnames(counts), design[, 
      #                                                      libID_col]), var_to_test]) - 1)
      cat(as.character(design[match(colnames(counts), design[, libID_col]), 
                 var_to_test]))
      cat("\n")
    }
    sink()
  }
