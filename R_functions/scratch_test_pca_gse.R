load("~/Documents/WORK/UW/ALTMAN/GITHUB/GRADUATE_batch2/data_clean/P662_2_voom_meta.RData")
temp_data<-dat_meta


expr=NULL        # input normalized expression matrix with genes/features as rows and samples/libraries as columns. rownames must be ensembl gene ids or HGNC symbols
gene_key =NULL   # input dataframe defining genes in expr. must include a column "symbol" denoting HGNC symbols
norm_data=dat_meta     # optionally an EList object which includes normalized expression matrix "E" and gene key "genes" which conform to specifications above.
min_geneset_size = 10   # min gene-set size to include
max_geneset_size = 2000 # max gene-set size to include
species = "Homo sapiens"
nPC=NULL
topN=250

rm(expr, gene_key, norm_data, min_geneset_size, max_geneset_size, species, nPC)


# test function

test_PC_enrichment<-pca_gse_light(norm_data = dat_meta)


# add libIDs as rownames (this is bug workaround)
dat_meta$targets<-
  dat_meta$targets%>%
  mutate(rowname=libID)%>%
  column_to_rownames("rowname")

dat_meta$genes <- dat_meta$genes %>% 
  mutate(geneName=ensembl_gene_id)

meta_NAC<- # filter to NAC samples
  filter(dat_meta$targets, !is.na(visit_nac_yr))

dat_NAC<-# filter to NAC samples
  RNAetc::subset_voom(dat_meta, lib_keep = meta_NAC$libID, libraryID = "libID")


NAC_PC_enrichment<-pca_gse_light(norm_data = dat_NAC)

NAC_PC_enrichment$plot_PC_enrichments

NAC_PC_enrichment$plot_PC_techscores_scatter

# refine dev

#Next Steps

#- reduce number of genesets rub by default (ie reduce redundancy)






# test on different data #

load("~/Documents/WORK/UW/ALTMAN/GITHUB/CATNIP/ANALYSIS/data_clean/CATNIP_346.1_dat_clean.RData")

dat_cat<-dat.voom
dat_cat$genes<-dat_cat$genes%>%mutate(symbol=hgnc_symbol, ensembl_gene_id=geneName)

cat_PC_enrichment<-pca_gse_light(norm_data = dat_cat)

cat_PC_enrichment$plot_PC_enrichments
cat_PC_enrichment$plot_PC_techscores

View(cat_PC_enrichment$PC_enrichments_top)
