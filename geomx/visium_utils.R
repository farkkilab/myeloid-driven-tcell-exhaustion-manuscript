# TODO prepare documentation

####################################################
# def
# inp
# args
# outp
run_unless_exists <- function(step_name, expected_output, script){
  if(!file.exists(expected_output)){
    source(script)
    print('$$$$$$$$$$')
    print(paste0(step_name, ' succeeded!'))
  } else{
    print(paste0(step_name, ' have been already run'))
  }
}
####################################################
# def
# inp
# args
# outp
visium_qc <- function(visium_obj){
  # count nr of mito and hemoglobine genes
  visium_obj <- PercentageFeatureSet(visium_obj, "^MT-", col.name = "percent_mito")
  visium_obj <- PercentageFeatureSet(visium_obj, "^HB.*", col.name = "percent_hb") # sth may be wrong with it
  
  # rmv spots with high nr of mit + hemoglobine genes + nr of detected genes < 500
  visium_obj <- visium_obj[, visium_obj$nFeature_Spatial > 500 & visium_obj$percent_mito < 25 & 
                             visium_obj$percent_hb < 20]
  
  # remove mit genes
  visium_obj <- visium_obj[!grepl("^MT-", rownames(visium_obj)), ]
  
  # remove non-expressed genes
  visium_obj <- visium_obj[which(rowSums(visium_obj@assays$Spatial$counts) != 0), ]
  
  return(visium_obj)
}
####################################################
# def
# inp
# args
# outp
draw_ncounts_plot <- function(visium_obj, out_dir){
  visium_obj <- SetIdent(visium_obj, value = 'histology')
  
  plot1 <- VlnPlot(visium_obj, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
  plot2 <- SpatialFeaturePlot(visium_obj, features = "nFeature_Spatial") + theme(legend.position = "right")
  plot3 <- wrap_plots(plot1, plot2)
  ggsave(file.path(out_dir, 'features_per_spot.png'))
}

####################################################
# def
# inp
# args
# outp
umap_and_plot <- function(visium_obj_norm, out_dir, save = F){
  # do PCA and UMAP reduction
  visium_obj_umap <- RunPCA(visium_obj_norm, assay = "SCT", verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30) %>%
    FindClusters(verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:30)
  
  # visualise UMAP
  p1 <- DimPlot(visium_obj_umap, reduction = "umap", label=T)
  p2 <- SpatialDimPlot(visium_obj_umap)
  
  hist_labs <- unique(visium_obj_umap@meta.data$histology)
  
  # visualise histology
  if (length(hist_labs) > 1){
    visium_obj_umap <- SetIdent(visium_obj_umap, value = 'histology')
    
    p3 <- DimPlot(visium_obj_umap, reduction = "umap", label=T)
    p4 <- SpatialDimPlot(visium_obj_umap)
    
    p1 + p2 + p3 + p4
    ggsave(file.path(out_dir, paste0('umap_', paste(hist_labs, collapse='_'), '.png')), width=10, height=10, units='in')
  } else if(length(hist_labs) == 1){
    p1 + p2
    ggsave(file.path(out_dir, paste0('umap_', paste(hist_labs, collapse='_'), '.png')))
  } else{
    p1 + p2
    ggsave(file.path(out_dir, paste0('umap_plain.png')))
  }
  
  # if needed for saving additionally subsetted data
  if(save){
    saveRDS(visium_obj_umap, 
            file.path(out_dir, 'umap', paste0('umap_', paste(hist_labs, collapse='_'), '.rds')))
  } else{
    return(visium_obj_umap)
  }
}
#################################################################
# for a given barcode and its position (rown, col in the array)
# finds all adjacent barcodes (closer spots vertically+horizontally+diagonally)
# def
# inp
# args
# outp
get_adj_barcodes <- function(barcode_pos_df, rown, coln){
  # TODO a bit inefficient. better to work on number pairs and filter at the end
  
  # 6-based
  barcode_adj <- barcode_pos_df %>%
    filter((array_row %in% seq(rown-1, rown+1)) & (array_col %in% seq(coln-2, coln+2))) %>%
    arrange(array_row, array_col) %>%
    filter(!((array_row == rown) & (array_col == coln)))
  
  # 8-based
  # barcode_adj <- barcode_pos %>%
  #   filter((array_row %in% seq(rown-2, rown+2)) & (array_col %in% seq(coln-2, coln+2))) %>%
  #   arrange(array_row, array_col) %>%
  #   filter(!((array_row %in% c(rown-2, rown+2)) & (array_col %in% c(coln-2, coln+2)))) %>%
  #   filter(!((array_row == rown) & (array_col == coln)))
  
  return(barcode_adj$barcode)
}

######################################################
# def
# inp
# args
# outp
# for a given condition in annotation df, find adjacent barcodes and
# save as new column wit 0/1 adjacency vector 
get_adj_barcodes_for_cond <- function(anno_df, barcode_pos_df, adj_anno_type, inp_colnames, outp_colnames){
  # get conditions according to type of anno
  
  
  # iterate trough conditions, return 0/1 vectors with adjacency information 
  adj_list <- lapply(inp_colnames, function(anno_cat){
    # get barcodes for given category
    if(adj_anno_type %in% c('histology', 'clone')){
      anno_cat_barc <- anno_df$barcode[anno_df[adj_anno_type] == anno_cat & !is.na(anno_df[adj_anno_type])] 
    } else{
      anno_cat_barc <- anno_df$barcode[anno_df[anno_cat] == 1] 
    }
    
    # get positions of barcodes from pos table
    anno_cat_barc_pos <- filter(barcode_pos_df, barcode %in% anno_cat_barc)
    
    # get all adjacent barcodes (unique)
    anno_cat_barc_adjacent <- apply(anno_cat_barc_pos, 1, function(x){
      adj_barc <- get_adj_barcodes(barcode_pos_df, as.numeric(x[['array_row']]), as.numeric(x[['array_col']]))
    })
    
    anno_cat_barc_adjacent <- unique(unlist(anno_cat_barc_adjacent))
    
    # change into 0/1 vector
    adj_vec <- ifelse(anno_df$barcode %in% anno_cat_barc_adjacent, 1, 0)
    return(adj_vec)
  })

  # combine with anno df table
  names(adj_list) <- outp_colnames
  anno_df <- cbind(anno_df, adj_list)
  return(anno_df)
}

################################################################
# make sorted factor of zonenames for plotting
# def
# inp
# args
# outp
sort_zone_names <- function(inp_df, zone_colname){
  zone_names_all <- unique(inp_df[[zone_colname]])
  zone_names <- mixedsort(zone_names_all[grepl('^tsi', zone_names_all)])

  if('tumortsi' %in% zone_names_all){
    inp_df[[zone_colname]] <- factor(x = inp_df[[zone_colname]],
                                     levels = c('tumor', 'tumortsi', zone_names, 'stroma'))
  } else{
    inp_df[[zone_colname]] <- factor(x = inp_df[[zone_colname]],
                                     levels = c('tumor', zone_names, 'stroma'))
  }
  return(inp_df)
}

###############################################################
###############################################################
# change ensembl/entrez into gene names for nested list of genes
# def
# inp
# args
# outp
gene_2names <- function(gene_inp_list, conv = c('ens', 'entrez'), type = 'list'){
  
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  
  #TODO make it for df if needed and also the other way around
  # gene_df_names <- getBM(attributes=c('external_gene_name', 'ensembl_gene_id'),
  #                           filters = 'ensembl_gene_id',
  #                           values = as.character(unlist(gene_list)),
  #                           mart = ensembl)
  # 
  # 
  # 
  # marker_ind <- left_join(marker_ind, marker_ind_names)
  # rm(marker_ind_names)
  # 
  if(conv == 'ens'){
    conv_name <- 'ensembl_gene_id'
  } else if(conv == 'entrez'){
    conv_name <- 'entrezgene_id'
  } else{stop()}
  
  if(type == 'list'){

    gene_list <- lapply(gene_inp_list, function(x){
      gene_names <- getBM(attributes=c('external_gene_name', conv_name),
                          filters = conv_name,
                          values = x,
                          mart = ensembl)
      
      gene_names <- gene_names$external_gene_name
    })
    return(gene_list)
  }
}

##############################################################
# no filtering, just parsing the df
# def
# inp
# args
# outp
process_signatures <- function(sign_path, use_names = T){
  
  # df into list
  marker_sign <- fread(sign_path, header = T)
  
  sign_list_orig <- lapply(colnames(marker_sign), function(x){
    sign <- marker_sign[[x]][marker_sign[[x]] != '']
    sign <- sign[!is.na(sign)]
  })
  names(sign_list_orig) <- colnames(marker_sign)
  
  # change to gene names if needed and list is with ensids or entrez
  if(use_names & all(grepl('^ENS', unlist(sign_list_orig)))){
    sign_list_orig <- gene_2names(sign_list_orig, 'ens', 'list')
  } else if(use_names & all(grepl('[1-9]', unlist(sign_list_orig)))){
    sign_list_orig <- gene_2names(sign_list_orig, 'entrez', 'list')
  }
  
  # remove null pathways from list
  sign_list_orig <- sign_list_orig[lengths(sign_list_orig) != 0]
  
  return(sign_list_orig)
}

##################################################3
# def
# inp
# args
# outp
plot_signatures_stats <- function(sign_list, visium_obj, output_dir){
  # calculate gene expression relative to sum of all transcripts expression
  rel_expression <- t( t(as.matrix(visium_obj@assays$SCT@counts)) /
                         colSums(as.matrix(visium_obj@assays$SCT@counts))) * 100

  dir.create(file.path(output_dir, 'signatures_stats'), showWarnings = TRUE)

  # plot relative expression for all genes
  for(sg_nr in seq(1:length(sign_list))){

    inter_genes <- intersect(rownames(rel_expression), sign_list[[sg_nr]])
    plot_data <- as.matrix(t(rel_expression[inter_genes,]))

    png(file= file.path(output_dir,'signatures_stats', paste0('sign_exp_boxp_', names(sign_list)[sg_nr], '.png')),
        width=2000, height=1200)
    boxplot(plot_data, cex=0.1, las=1, xlab="% total count per cell", horizontal=TRUE, main=names(sign_list)[sg_nr])
    dev.off()
  }
}

#############################################################
# process signatures df into list and filter based on relative expression
# removes pathways with too litlle genes expressed
# def
# inp
# args
# outp
filter_signatures <- function(sign_list_orig, visium_obj, output_dir, 
                               min_expr_thr = 0.005, min_nr_above = 10, 
                               filter_abovemin = F){
  
  # calculate gene expression relative to sum of all transcripts expression
  rel_expression <- t( t(as.matrix(visium_obj@assays$SCT@counts)) /
                         colSums(as.matrix(visium_obj@assays$SCT@counts))) * 100
  
  # rmv genes from signatures: min_expr_thr of transcripts per cell
  row_sub <- apply(rel_expression, 1, function(row) mean(row)>min_expr_thr)
  rel_expression_abovemin <- rel_expression[row_sub,]
  
  # update signatures list
  sign_list_abovemin <- lapply(sign_list_orig, function(x){
    sign_filt <- x[x %in% rownames(rel_expression_abovemin)]
  })
  
  #  remove pathways with < min nr of genes abovemin
  sign_list <- sapply(seq(1:length(sign_list_orig)), function(n){
    
    if(filter_abovemin){ # return only genes abovemin (better no for UCell)
      sign_list_return <- sign_list_abovemin[n]
    } else{
      sign_list_return <- sign_list_orig[n]
    }
    
    sign_filt <- ifelse((length(sign_list_abovemin[[n]]) > min_nr_above), sign_list_return, list())
  })
  
  names(sign_list) <- names(sign_list_orig)

  # remove null pathways from list
  sign_list <- sign_list[lengths(sign_list) != 0]
  
  return(sign_list)
}

######################################################
#!!! min_val and max_val not used for viridis
get_treatment_heatmap <- function(dt_treat_mt, treat_name, min_val, max_val, legend_name,
                                  keep_legend, color_scale='rb', clust_rows = F){
  if(color_scale == 'rb'){
    col_fun <- colorRamp2(c(min_val, 0, max_val), hcl_palette = "RdBu", reverse = TRUE)
  } else if(color_scale == 'viridis'){
    col_fun <- viridis(100)
  }
  
  

  row_size <- ifelse(nrow(dt_treat_mt) > 14, 4, ifelse(nrow(dt_treat_mt) > 12, 7, 10 ))
  
  treat_heat <- Heatmap(dt_treat_mt, height = unit(4.5, "cm") , width = unit(3, "cm"),border="white",
                        rect_gp = gpar(col = "white", lwd = 2), name=legend_name,
                        cluster_columns = F, cluster_rows= clust_rows, col = col_fun,
                        column_title = treat_name, row_names_gp = gpar(fontsize = row_size),
                        column_names_gp = gpar(fontsize = 6), column_title_gp = gpar(fontsize = 10),
                        show_heatmap_legend = keep_legend, column_names_rot = 45)
  
  
  return(treat_heat)
}

#########################################################
calculate_ora <- function(gene_vect, bcg_gene_vect, msigdb_df, padj = 0.1){
  ora <- enricher(
    gene = gene_vect,
    pvalueCutoff = padj, # Can choose a FDR cutoff
    pAdjustMethod = "BH", 
    universe = bcg_gene_vect, 
    TERM2GENE = dplyr::select(msigdb_df, gs_name, human_gene_symbol)
  )
  
  if(!is.null(ora)){
    ora_df <- data.frame(ora@result) %>%
      filter(p.adjust <= padj)
  } else{
    ora_df <- data.frame()
  }
  return(ora_df)
}
