#!/usr/bin/env Rscript

library(ggplot2) # for plotting
library(viridis) # for colours
library(PhyloProfile)
library(data.table)


library(yaml) # for reading config file
cfg <- yaml.load_file("config.yml")

plot_pdf_andor_png <- function(theplot, thepath, spec) {
  # spec: if dimensions and resolution should be specified
  if (spec) {
    if (cfg$output_png) {
      png(paste(c(thepath,".png"), collapse=""), units="in", width=6, height=6, res=500)
      print(theplot)
      dev.off()
    }
    if (cfg$output_pdf) {
      pdf(paste(c(thepath,".pdf"), collapse=""), width=6, height=6)
      print(theplot)
      # theplot
      dev.off()
    }
  }

  # if not (bc e.g. fviz plots may be impaired by this)
  else {
    if (cfg$output_png) {
      png(paste(c(thepath,".png"), collapse=""))
      print(theplot)
      dev.off()
    }
    if (cfg$output_pdf) {
      pdf(paste(c(thepath,".pdf"), collapse=""))
      print(theplot)
      # theplot
      dev.off()
    }
  }
}

id_to_name <- function(idList) {
    # get pre-processed taxonomy file
    ncbiFilein <- paste0(
        find.package("PhyloProfile"),
        "/PhyloProfile/data/preProcessedTaxonomy.txt"
    )
    if (file.exists(ncbiFilein)) {
        # print("Loading NCBI taxonomy file...")
        preProcessedTaxonomy <- data.table::fread(ncbiFilein)
    } else {
        # downloading & processing ncbi taxonomy file if needed
        # print("Download and process NCBI taxonomy files...")
        preProcessedTaxonomy <- processNcbiTaxonomy()
        # save to text (tab-delimited) file
        write.table(
            preProcessedTaxonomy,
            file = ncbiFilein,
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t"
        )
    }
    # print("done")
    # get taxon names
    nameList <- preProcessedTaxonomy[
        preProcessedTaxonomy$ncbiID %in% idList,
    ][,c("ncbiID","fullName")]
    return(nameList)
}

saveWidgetFix <- function (widget,file,...) {
  # workaround to saveWidget for relative paths
  wd<-getwd()
  on.exit(setwd(wd))
  outDir<-dirname(file)
  file<-basename(file)
  setwd(outDir);
  saveWidget(widget,file=file,...)
}


# _____________________DATA ASSEMBLY_______________________________

# read whole raw_gene_table_coords
gene_table_coords_w_nan <- data.table::fread(paste(c(cfg$output_path,'PCA_and_clustering/gene_table_coords.csv'), collapse=''), sep=",", header=TRUE)
# extract genes that where included in the PCA (coordinates available; excluded genes have coordinates = NaN)
gene_table_coords <- gene_table_coords_w_nan[complete.cases(gene_table_coords_w_nan[ , "Dim.1"]),]
# extract coordinates only



# _____________________READ LCA AND BEST HITS_______________________________
# read lca taxonomic assignment file and assign header
prots_lca_hits <- data.table::fread(cfg$taxon_hits_lca_path, sep="\t", header=FALSE, stringsAsFactors = FALSE)
setnames(prots_lca_hits, c("queryID", "lcaID", "score"))
# read blast tabular taxonomic assignment file with best hit for each query
prots_best_hit <- data.table::fread(cfg$best_taxon_hit_path, sep="\t", header=FALSE, stringsAsFactors = FALSE)
setnames(prots_best_hit, c("queryID", "sseqid", "evalue", "bitscore", "score", "pident", "best_hitID", "best_hit"))
prots_best_hit <- subset(prots_best_hit, select = c("queryID", "best_hitID", "best_hit", "evalue"))
prots_best_hit <- subset(prots_best_hit, best_hitID!="")
# take only the first ID and fullName (the others are often specific strains)
prots_best_hit$best_hitID <- gsub("?;.*", "", prots_best_hit$best_hitID)
prots_best_hit$best_hit <- gsub("?;.*", "", prots_best_hit$best_hit)

# _____________________IDENTIFY GENE ID FOR HITS_______________________________
# read protein ID to gene ID matching files
prot_gene_matching <- data.table::fread(paste0(cfg$output_path,'tmp.prot_gene_matching.txt'), sep="\t", header=FALSE, stringsAsFactors = FALSE)
setnames(prot_gene_matching, c("protID", "geneID"))
if ("TRUE" %in% (prots_lca_hits$queryID %in% prot_gene_matching$protID)){
  if ("FALSE" %in% (prots_lca_hits$queryID %in% prot_gene_matching$protID)){
    print("IDs in LCA hit file may be erroneous")
  }
  genes_lca_hits <- merge(prot_gene_matching, prots_lca_hits, by.x = "protID", by.y = "queryID",  all.y = TRUE)
}else{
  genes_lca_hits <- merge(prot_gene_matching, prots_lca_hits, by.x = "geneID", by.y = "queryID",  all.y = TRUE)

}
if ("TRUE" %in% (prots_best_hit$queryID %in% prot_gene_matching$protID)){
  if ("FALSE" %in% (prots_best_hit$queryID %in% prot_gene_matching$protID)){
    print("IDs in best hit file may be erroneous")
  }
  genes_best_hit <- merge(prot_gene_matching, prots_best_hit, by.x = "protID", by.y = "queryID",  all.y = TRUE)
}else{
  genes_best_hit <- merge(prot_gene_matching, prots_best_hit, by.x = "geneID", by.y = "queryID",  all.y = TRUE)
}

# remove tables from memory as they are no longer needed
rm(prot_gene_matching)
rm(prots_lca_hits)
rm(prots_best_hit)



# ______________________TRANSLATE NCBI ID FOR LCA_________________________________
## translate LCAs ncbiID to name and add the merging_rank information ##
# retrieve names for NCBI IDs and merge to taxon hits table

lca_names <- id_to_name(genes_lca_hits$lcaID)
colnames(lca_names) <- c("lcaID", "lca_hit")



# merge data
transl_taxon_hits_part <- merge(genes_lca_hits, lca_names, by = "lcaID", all.x = TRUE)
transl_taxon_hits <- merge(transl_taxon_hits_part, genes_best_hit, by = c("geneID", "protID"), all = TRUE)
# merge taxon info with gene coordinates (raw gene table + PCA coordinates)
genes_coords_taxon <- merge(gene_table_coords, transl_taxon_hits, by.x = "g_name", by.y="geneID", all.x = TRUE)
genes_coords_taxon <- transform(genes_coords_taxon, best_hitID = as.numeric(best_hitID))



# ____________PERFORM RANK REPLACEMENTS FOR LCA AND QUERY GROUP___________________

genes_coords_taxon$labelID <- as.numeric(genes_coords_taxon$lcaID)
genes_coords_taxon$labelID[is.na(genes_coords_taxon$labelID)] <- 0


ncbiFilein <- paste0(find.package("PhyloProfile"),"/PhyloProfile/data/preProcessedTaxonomy.txt")
if (file.exists(ncbiFilein)) {
    # print("Loading NCBI taxonomy file...")
    preProcessedTaxonomy <- data.table::fread(ncbiFilein)
} else {
    # downloading & processing ncbi taxonomy file if needed
    # print("Download and process NCBI taxonomy files...")
    preProcessedTaxonomy <- processNcbiTaxonomy()
    # save to text (tab-delimited) file
    write.table(
        preProcessedTaxonomy,
        file = ncbiFilein,
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE,
        sep = "\t"
    )
}

# get names of all descendants of query species
query_taxonomy <- PhyloProfile::getTaxonomyInfo(cfg$tax_id,preProcessedTaxonomy)[[1]]
# retrieve ranks of interest for query
query_plot_rank <- query_taxonomy[rank == cfg$plot_grouping_rank]$ncbiID
query_label <- query_plot_rank #paste0("Query ", cfg$plot_grouping_rank) # label for this group in the plots
query_replace_rank <- query_taxonomy[rank == cfg$lca_replacement_rank]$ncbiID

# retrieve <lca_replacement_rank> of best hits to determine if the best hit should replace LCA
best_hits_taxonomy <- PhyloProfile::getTaxonomyInfo(unique(genes_best_hit$best_hitID),preProcessedTaxonomy)
best_hit_replace_rank <- lapply(best_hits_taxonomy, function(x) x$ncbiID[x$rank == cfg$lca_replacement_rank])
best_hit_replace_rank <- lapply(best_hit_replace_rank, function(x) if(identical(x, numeric(0))) "higherrank" else x)
best_hit_replace_rank_name <- lapply(best_hits_taxonomy, function(x) x$fullName[x$rank == cfg$lca_replacement_rank])
best_hit_replace_rank_name <- lapply(best_hit_replace_rank_name, function(x) if(identical(x, character(0))) "higherrank" else x)
best_hit_id <- lapply(best_hits_taxonomy, function(x) x$ncbiID[1])
best_hit_replace_matching <- do.call(rbind, Map(data.frame, "best_hitID"=best_hit_id, "lca_replacement_rank"=best_hit_replace_rank, "lca_replacement_rank_name"=best_hit_replace_rank_name))
genes_coords_taxon <- merge(genes_coords_taxon, best_hit_replace_matching, by="best_hitID", all.x = TRUE)
genes_coords_taxon$lca_replacement_rank[is.na(genes_coords_taxon$lca_replacement_rank)] <- 0
genes_coords_taxon$labelID <- ifelse(genes_coords_taxon$lca_replacement_rank == query_replace_rank, genes_coords_taxon$best_hitID, genes_coords_taxon$labelID)


# retrieve <plot_grouping_rank> of current label to determine if it's in same <plot_grouping_rank> as query (and thus be merged into one group/label)
lca_hits_taxonomy <- PhyloProfile::getTaxonomyInfo(unique(genes_coords_taxon$labelID[genes_coords_taxon$labelID != "0"]),preProcessedTaxonomy)
lca_hit_plot_rank <- lapply(lca_hits_taxonomy, function(x) x$ncbiID[x$rank == cfg$plot_grouping_rank])
lca_hit_plot_rank <- lapply(lca_hit_plot_rank, function(x) if(identical(x, numeric(0))) "higherrank" else x)
lca_hit_plot_rank_name <- lapply(lca_hits_taxonomy, function(x) x$fullName[x$rank == cfg$plot_grouping_rank])
lca_hit_plot_rank_name <- lapply(lca_hit_plot_rank_name, function(x) if(identical(x, character(0))) "higherrank" else x)
lca_hit_id <- lapply(lca_hits_taxonomy, function(x) x$ncbiID[1])
lca_hit_group_matching <- do.call(rbind, Map(data.frame, "labelID"=lca_hit_id, "plot_grouping_rank"=lca_hit_plot_rank, "plot_grouping_rank_name"=lca_hit_plot_rank_name))
genes_coords_taxon <- merge(genes_coords_taxon, lca_hit_group_matching, by="labelID", all.x = TRUE)
genes_coords_taxon$labelID <- ifelse(genes_coords_taxon$plot_grouping_rank == query_plot_rank, query_plot_rank, genes_coords_taxon$labelID)
# get full name of ncbiID of labels
label_names <- id_to_name(genes_coords_taxon$labelID)
colnames(label_names) <- c("labelID", "label")
genes_coords_taxon <- merge(genes_coords_taxon, label_names, by = "labelID", all.x = TRUE)

# lca_hit equals NaN for genes excluded from assignment; replace NaN with "Unassigned"
genes_coords_taxon$label[is.na(genes_coords_taxon$label)] <- "Unassigned"
#genes_coords_taxon$label[genes_coords_taxon$lca_hit == "root"] <- "NCBI"


# __________________TEXT OUTPUT________________________

# csv output of genes with orphan + taxon info and PCA coordinates
if (cfg$include_coverage == "TRUE") {
  num_cov_cols <- length(cfg$pbc_paths)*10 # number of coverage related columns in genes_coords_taxon; information needed for indexing
}else {
  num_cov_cols <- 1  # by default there is one set of mock coverage data produced
}

# indexes of columns to be subsetted and reordered for output
print(colnames(genes_coords_taxon))
cols <- c(3:(30+num_cov_cols),32+num_cov_cols,34+num_cov_cols,33+num_cov_cols,2,(35+num_cov_cols),(36+num_cov_cols),(41+num_cov_cols)) #,(40+num_cov_cols),(41+num_cov_cols))
raw_gene_table_orphan_info <- data.frame(genes_coords_taxon)[,cols]
names(raw_gene_table_orphan_info)[names(raw_gene_table_orphan_info) == "rank"] <- cfg$plot_grouping_rank
write.csv(raw_gene_table_orphan_info, file=paste(c(cfg$output_path, "taxonomic_assignment/raw_gene_table_taxon.csv"), collapse=""), row.names=FALSE, quote=FALSE)

# _______DATA PREPARATION FOR PLOT GENERATION_________


#count how often each taxon was matched
taxon_count_table <- table(genes_coords_taxon$label) # look in label as these are the groups to be plotted
if (cfg$taxon_hit_threshold == "min") {
  #transform it into dataframe to sort it
  taxon_count_df <- data.frame(taxon_count_table)
  taxon_count_sort <- taxon_count_df[order(-taxon_count_df$Freq),]
  if (nrow(taxon_count_sort) > 18) {
    # can display 18 different groups: threshold equals the abundancy of the 18th least occuring taxon
    cfg$taxon_hit_threshold <- taxon_count_sort[18,2]+1
  }else{
    cfg$taxon_hit_threshold <- 1
  }
}
#print(cfg$taxon_hit_threshold)
# genes with taxonomic assignemnt of low abundance get label "Otherwise assigned"
genes_coords_taxon$label[genes_coords_taxon$label %in% names(taxon_count_table[taxon_count_table < cfg$taxon_hit_threshold])] <- "Otherwise assigned"


genes_coords_taxon_query <- genes_coords_taxon[genes_coords_taxon$label == query_label,]
genes_coords_taxon_rest <- genes_coords_taxon[genes_coords_taxon$label != query_label,]



# label for plots
label_count_table <- data.frame(table(genes_coords_taxon$label)) # how often does each label appear?
label_count_table$label <- as.character(paste0(label_count_table$Var1," (",label_count_table$Freq,")"))
genes_coords_taxon <- merge(genes_coords_taxon, label_count_table, by.x = "label", by.y = "Var1", all.x = TRUE)
label <- sort(unique(label_count_table$label))


# 1D density of 1 element groups is not possible
# for 1D density only use taxons with abundancy of at least 2
#genes_taxon_1dim <- genes_coords_taxon_plot[genes_coords_taxon_plot$taxon_hit %in% names(taxon_count_table[taxon_count_table >= 2]),]
genes_taxon_1dim <- genes_coords_taxon[genes_coords_taxon$label %in% names(taxon_count_table[taxon_count_table > 1]),]

# label for density plots
label_1d <- c(sort(unique(genes_taxon_1dim$label)))

# _______________PLOT GENERATION_____________________

# taxon_hit column as factor for marking in plots
taxon_rest_factor <- as.factor(genes_coords_taxon_rest$label)
taxon_query_factor <- as.factor(genes_coords_taxon_query$label)
taxon_factor <- as.factor(genes_coords_taxon$label)

taxons_1d <- as.factor(genes_taxon_1dim$label)



plot_o <- ggplot() +
                  scale_color_viridis(name="species", discrete=TRUE, begin=0.1, end=0.9, labels=label) +
                  scale_shape_manual(name="species", labels=label, values=c(1,2,15,16,17,3,7,8,9,23,18,4,10,11,12,13,14,0)) +
                  geom_point(data=genes_coords_taxon_query, aes(x=Dim.1, y=Dim.2, shape=taxon_query_factor, color=taxon_query_factor), alpha = 0.2) +
                  geom_point(data=genes_coords_taxon_rest, aes(x=Dim.1, y=Dim.2, shape=taxon_rest_factor, color=taxon_rest_factor)) +
                  #geom_point(data=genes_taxon_grouped, aes(x=Dim.1, y=Dim.2, shape=1, color=1)) +
                  theme_bw() +
                  #theme(legend.title = element_text(size = 16), legend.text = element_text(size = 13)) + #enlarge legend
                  #guides(shape = guide_legend(override.aes = list(size = 3))) +
                  ggtitle("all genes with orphan and taxon info")

plot_pdf_andor_png(plot_o, paste(c(cfg$output_path, "taxonomic_assignment/dot_plot_1_2"), collapse=""), TRUE)

plot_x <- ggplot() +
                  scale_color_viridis(name="species", discrete=TRUE, begin=0.1, end=0.9, labels=label_1d) +
                  #geom_density(data=genes_taxon_grouped, aes(x=Dim.1, color=taxons_group), linetype="dashed", alpha=0.25) +
                  geom_density(data=genes_taxon_1dim, aes(x=Dim.1, color=taxons_1d)) +
                  theme_bw() +
                  #expand_limits(x = c(0, max(coords$Dim.1))) +
                  ggtitle("density of Dim.1")

plot_pdf_andor_png(plot_x, paste(c(cfg$output_path, "taxonomic_assignment/density_x"), collapse=""), TRUE)

plot_y <- ggplot() +
                  scale_color_viridis(name="species", discrete=TRUE, begin=0.1, end=0.9, labels=label_1d) +
                  #geom_density(data=genes_taxon_grouped, aes(x=Dim.2, color=taxons_group), linetype="dashed", alpha=0.25) +
                  geom_density(data=genes_taxon_1dim, aes(x=Dim.2, color=taxons_1d)) +
                  theme_bw() +
                  #expand_limits(x = c(0, max(coords$Dim.2))) +
                  ggtitle("density of Dim.2")

plot_pdf_andor_png(plot_y, paste(c(cfg$output_path, "taxonomic_assignment/density_y"), collapse=""), TRUE)

library(plotly)
library(htmlwidgets)

prot_fasta <- Biostrings::readAAStringSet(cfg$proteins_path)
protID <- names(prot_fasta)
query_seq <- paste(prot_fasta)
prot_sequences <- data.frame(protID, query_seq)
prot_sequences$query_seq <- gsub("(.{70}?)", "\\1</br>", prot_sequences$query_seq)
genes_coords_taxon <- merge(genes_coords_taxon, prot_sequences, by="protID", all.x=TRUE)
genes_coords_taxon$g_terminal[genes_coords_taxon$g_terminal == "0"] <- "no"
genes_coords_taxon$g_terminal[genes_coords_taxon$g_terminal == "1"] <- "yes"

g_cov_vars <- as.factor(grep("g_cov_[0-9]",colnames(genes_coords_taxon), value=TRUE))
genes_coords_taxon$g_coverages <- apply(subset(genes_coords_taxon, select=g_cov_vars),1,paste,collapse="; ")
g_covdev_vars <- as.factor(grep("g_covdev_c_[0-9]",colnames(genes_coords_taxon), value=TRUE))
genes_coords_taxon$g_covdeviations <- apply(subset(genes_coords_taxon, select=g_covdev_vars),1,paste,collapse="; ")


fig <- plot_ly(genes_coords_taxon, type="scatter", mode="markers", x = ~Dim.1, y = ~Dim.2,
    hoverinfo = 'text',
    text = ~paste('</br>ID:',g_name,
                  '</br>Coverage:',g_coverages, '(SD from contig mean:',g_covdeviations,')',
                  '</br>Terminal:',g_terminal,
                  '</br>LCA:',lca_hit,
                  '</br>Best hit:',best_hit,'(e-value:',evalue,')',
                  '</br>Seq:',prot_sequences$query_seq),
    symbol=~as.factor(label.y),
    symbols=c('circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond'), # additional options: 'circle-open', 'square-open', 'diamond-open', 'x'
    color=~as.factor(label.y),
    colors=viridis_pal(option="D")(3),
    size=8
  )
fig <- fig %>% layout(title = 'taxonomic assignment')
saveWidgetFix(as_widget(fig),paste(c(cfg$output_path, "taxonomic_assignment/2D_plot.html"), collapse=""), selfcontained=FALSE)

if (as.numeric(cfg$num_pcs) >= 3){

  plot_1_3 <- ggplot() +
                    scale_color_viridis(name="species", discrete=TRUE, begin=0.1, end=0.9, labels=label) +
                    scale_shape_manual(name="species", labels=label, values=c(1,2,15,16,17,3,7,8,9,23,18,4,10,11,12,13,14,0)) +
                    geom_point(data=genes_coords_taxon_query, aes(x=Dim.1, y=Dim.3, shape=taxon_query_factor, color=taxon_query_factor), alpha = 0.1) +
                    geom_point(data=genes_coords_taxon_rest, aes(x=Dim.1, y=Dim.3, shape=taxon_rest_factor, color=taxon_rest_factor)) +
                    theme_bw() +
                    #theme(legend.title = element_text(size = 16), legend.text = element_text(size = 13)) + #enlarge, legend
                    #guides(shape = guide_legend(override.aes = list(size = 3))) +
                    ggtitle("all genes with orphan and taxon info")

  plot_pdf_andor_png(plot_1_3, paste(c(cfg$output_path, "taxonomic_assignment/dot_plot_1_3"), collapse=""), TRUE)

  plot_2_3 <- ggplot() +
                    scale_color_viridis(name="species", discrete=TRUE, begin=0.1, end=0.9, labels=label) +
                    scale_shape_manual(name="species", labels=label, values=c(1,2,15,16,17,3,7,8,9,23,18,4,10,11,12,13,14,0)) +
                    geom_point(data=genes_coords_taxon_query, aes(x=Dim.2, y=Dim.3, shape=taxon_query_factor, color=taxon_query_factor), alpha = 0.1) +
                    geom_point(data=genes_coords_taxon_rest, aes(x=Dim.2, y=Dim.3, shape=taxon_rest_factor, color=taxon_rest_factor)) +
                    theme_bw() +
                    #theme(legend.title = element_text(size = 16), legend.text = element_text(size = 13)) + #enlarge legend
                    #guides(shape = guide_legend(override.aes = list(size = 3))) +
                    ggtitle("all genes with orphan and taxon info")

  plot_pdf_andor_png(plot_2_3, paste(c(cfg$output_path, "taxonomic_assignment/dot_plot_2_3"), collapse=""), TRUE)


  plot_z <- ggplot() +
                    scale_color_viridis(discrete=TRUE, begin=0.1, end=0.9, labels=label_1d) +
                    #geom_density(data=genes_taxon_grouped, aes(x=Dim.3, color=taxons_group), linetype="dashed", alpha=0.25) +
                    geom_density(data=genes_taxon_1dim, aes(x=Dim.3, color=taxons_1d)) +
                    labs(color="color (species)") +
                    theme_bw() +
                    scale_x_continuous(limits = c(-5, 5)) +
                    ggtitle("density of Dim.3")

  plot_pdf_andor_png(plot_z, paste(c(cfg$output_path, "taxonomic_assignment/density_z"), collapse=""), TRUE)

  genes_coords_taxon_query <- genes_coords_taxon[genes_coords_taxon$labelID == query_label,]
  genes_coords_taxon_rest <- genes_coords_taxon[genes_coords_taxon$labelID != query_label | is.na(genes_coords_taxon$labelID),]


  fig <- plot_ly(type="scatter3d", mode="markers",
        symbols=c('circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond', 'circle', 'square', 'diamond'), # additional options: 'circle-open', 'square-open', 'diamond-open', 'x')
        colors=viridis_pal(option="D")(3),
        size=8)

  fig <- fig %>% add_markers(fig, data=genes_coords_taxon_query, x = ~Dim.1, y = ~Dim.2, z = ~Dim.3,
      hoverinfo = 'text',
      text = ~paste('</br>ID:',g_name,
                    '</br>Coverage:',g_coverages, '(SD from contig mean:',g_covdeviations,')',
                    '</br>Terminal:',g_terminal,
                    '</br>LCA:',lca_hit,
                    '</br>Best hit:',best_hit,'(e-value:',evalue,')',
                    '</br>Seq:',query_seq),
      #'</br>Coordinates: (',round(Dim.1,2),',',round(Dim.2,2),')',
      opacity=0.5,
      symbol=~as.factor(label.y),
      color=~as.factor(label.y),
      textposition="bottom right"
    )
    fig <- fig %>% add_markers(fig, data=genes_coords_taxon_rest, x = ~Dim.1, y = ~Dim.2, z = ~Dim.3,
        hoverinfo = 'text',
        text = ~paste('</br>ID:',g_name,
                      '</br>Coverage:',g_coverages, '(SD from contig mean:',g_covdeviations,')',
                      '</br>Terminal:',g_terminal,
                      '</br>LCA:',lca_hit,
                      '</br>Best hit:',best_hit,'(e-value:',evalue,')',
                      '</br>Seq:',query_seq),
        #'</br>Coordinates: (',round(Dim.1,2),',',round(Dim.2,2),')',
        symbol=~as.factor(label.y),
        color=~as.factor(label.y),
        textposition="bottom right"
      )

  fig <- fig %>% layout(title = 'taxonomic assignment')
  saveWidgetFix(as_widget(fig),paste(c(cfg$output_path, "taxonomic_assignment/3D_plot.html"), collapse=""), selfcontained=FALSE)



}
