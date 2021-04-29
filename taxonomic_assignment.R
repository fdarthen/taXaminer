#!/usr/bin/env Rscript

library(PhyloProfile)
library(data.table)
library(yaml) # for reading config file

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[1]

cfg <- yaml.load_file(config_path)


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


# _____________________DATA ASSEMBLY_______________________________

# read whole raw_gene_table_coords
gene_table_coords_w_nan <- data.table::fread(paste0(cfg$output_path,'PCA_and_clustering/gene_table_coords.csv'), sep=",", header=TRUE)
# extract genes that where included in the PCA (coordinates available; excluded genes have coordinates = NaN)
gene_table_coords <- gene_table_coords_w_nan[complete.cases(gene_table_coords_w_nan[ , "Dim.1"]),]



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
prot_gene_matching <- data.table::fread(paste0(cfg$output_path,'tmp.prot_gene_matching.txt'), sep="\t", header=FALSE, stringsAsFactors = FALSE)[,1:2]
setnames(prot_gene_matching, c("protID", "geneID"))

# if the GFF has the prefixes "transcript:" and "gene:" for the IDs and the protein FASTA doesn't
# delete the "transcript:" prefix from the matching table to allow it to work
if (startsWith(prot_gene_matching$protID[1],"transcript:")) {
    if (!startsWith(prots_best_hit$queryID[1],"transcript:")) {
        prot_gene_matching$protID <- gsub("transcript:", "", prot_gene_matching$protID)
    }
}

if ("TRUE" %in% (prots_lca_hits$queryID %in% prot_gene_matching$protID)){ # if queryID is protein ID
    if ("FALSE" %in% (prots_lca_hits$queryID %in% prot_gene_matching$protID)){# no clear indication whether it is protein or gene ID
        print("IDs in LCA hit file can not be mapped to gene IDs in GFF")
    }
    genes_lca_hits <- merge(prot_gene_matching, prots_lca_hits, by.x = "protID", by.y = "queryID",  all.y = TRUE)
}else{ # else queryID is gene ID
  genes_lca_hits <- merge(prot_gene_matching, prots_lca_hits, by.x = "geneID", by.y = "queryID",  all.y = TRUE)

}
if ("TRUE" %in% (prots_best_hit$queryID %in% prot_gene_matching$protID)){
    if ("FALSE" %in% (prots_best_hit$queryID %in% prot_gene_matching$protID)){
        print("IDs in best hit file can not be mapped to gene IDs in GFF")
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

#genes_lca_hits$lcaID[is.na(genes_lca_hits$lcaID)] <- 0
#lca_names_ranks <- PhyloProfile::getIDsRank(unique(genes_lca_hits$lcaID[genes_lca_hits$lcaID != "0"]),preProcessedTaxonomy)[[3]][,c("ncbiID","fullName","rank")]
#colnames(lca_names_ranks) <- c("lcaID", "lca_hit","lca_rank")

# merge data
transl_taxon_hits_part <- merge(genes_lca_hits, lca_names, by = "lcaID", all.x = TRUE)
transl_taxon_hits <- merge(transl_taxon_hits_part, genes_best_hit, by = c("geneID", "protID"), all = TRUE)
# merge taxon info with gene coordinates (raw gene table + PCA coordinates)
genes_coords_taxon <- merge(gene_table_coords, transl_taxon_hits, by.x = "g_name", by.y="geneID", all.x = TRUE)

# required for further data processing
genes_coords_taxon <- transform(genes_coords_taxon, best_hitID = as.numeric(best_hitID))
genes_coords_taxon$lcaID[genes_coords_taxon$lcaID == "0"] <- NA


# ____________PERFORM RANK REPLACEMENTS FOR LCA AND QUERY GROUP___________________

ncbiFilein <- paste0(find.package("PhyloProfile"),"/PhyloProfile/data/preProcessedTaxonomy.txt")
preProcessedTaxonomy <- data.table::fread(ncbiFilein)


# get names of all descendants of query species
query_taxonomy <- PhyloProfile::getTaxonomyInfo(cfg$tax_id,preProcessedTaxonomy)[[1]]

# rank at which labels are merged in the plots
# if cfg$plot_grouping_rank is numeric, it is already the ncbiID of the rank at which should be grouped
# else it is the rank (like phylum/class/etc.) at which should be grouped and the ID must be retained from the taxonomy first
query_plot_rank <- ifelse(grepl("\\d", cfg$plot_grouping_rank),as.numeric(cfg$plot_grouping_rank),query_taxonomy[rank == cfg$plot_grouping_rank]$ncbiID)
query_label <- query_plot_rank
query_label_name <- query_taxonomy$fullName[query_taxonomy$ncbiID == query_plot_rank]
# retrieve lineage information for query
query_lineage <- c(query_taxonomy$ncbiID,"1") # 1 equals the root node, which is per default not included in the taxonomy


# when the LCA is part of the lineage of the query, it is a candidate to be be replaced by the corrected LCA
genes_coords_taxon$LCA_replace <- (genes_coords_taxon$lcaID %in% query_lineage)

# compute (corrected) LCA of best hit and query species for genes where LCA_replace is TRUE,
# i.e. LCA is within query species lineage

best_hits_taxonomy <- PhyloProfile::getTaxonomyInfo(unique(genes_coords_taxon$best_hitID[!is.na(genes_coords_taxon$best_hitID) & genes_coords_taxon$LCA_replace == TRUE]),preProcessedTaxonomy)
best_hit_lineage <- lapply(best_hits_taxonomy, function(x) x$ncbiID)
corrected_lcaID <- lapply(best_hit_lineage, function(x) (intersect(x,query_lineage))[1]) # ID of LCA of best hit and query species, else NA
counter <- 0
corrected_lca <- lapply(best_hits_taxonomy, function(x) {counter <<- counter + 1; x$fullName[x$ncbiID == corrected_lcaID[counter][1]]})
best_hit_id <- lapply(best_hits_taxonomy, function(x) x$ncbiID[1]) # get for each position in the list the corresponding best hit ID
corrected_lca_for_best_hits <- do.call(rbind, Map(data.frame, "best_hitID"=best_hit_id, "corrected_lcaID"=as.numeric(corrected_lcaID), "corrected_lca"=corrected_lca))
genes_coords_taxon <- merge(genes_coords_taxon, corrected_lca_for_best_hits, by="best_hitID", all.x = TRUE)
# if the corrected LCA is NA, it means that there was no interesection between the lineages of the best hit and the query
# thus the best hit shares no lineage with the query
genes_coords_taxon$taxon_assignmentID <- ifelse(genes_coords_taxon$LCA_replace == TRUE, ifelse(is.na(genes_coords_taxon$corrected_lcaID), genes_coords_taxon$lcaID,genes_coords_taxon$corrected_lcaID), genes_coords_taxon$lcaID)


# __________MERGING OF TAXONOMIC ASSIGNMENTS_________ #
# taxonomic assignments thats are closely related to the query species are merged in the plots
# also, the least abundant taxonomic assignments are merged, so that there are in total 18 groups left

# get full name of ncbiID of taxonomic assignments
label_names <- id_to_name(genes_coords_taxon$taxon_assignmentID)
colnames(label_names) <- c("taxon_assignmentID", "taxon_assignment")
genes_coords_taxon <- merge(genes_coords_taxon, label_names, by = "taxon_assignmentID", all.x = TRUE)

# retrieve <query_plot_rank> of current taxon_assignment to determine if it's in same <query_plot_rank> as query (and thus be merged into one group/taxon_assignment)
labels_taxonomy <- PhyloProfile::getTaxonomyInfo(unique(genes_coords_taxon$taxon_assignmentID[!is.na(genes_coords_taxon$taxon_assignmentID)]),preProcessedTaxonomy)
label_grouping <- lapply(labels_taxonomy, function(x) ifelse(query_plot_rank %in% x$ncbiID, TRUE, FALSE))
label_id <- lapply(labels_taxonomy, function(x) x$ncbiID[1])
label_group_matching <- do.call(rbind, Map(data.frame, "taxon_assignmentID"=as.numeric(label_id), "label_grouping"=label_grouping))

genes_coords_taxon <- merge(genes_coords_taxon, label_group_matching, by="taxon_assignmentID", all.x = TRUE)
genes_coords_taxon$plot_label <- ifelse(genes_coords_taxon$label_grouping == TRUE, query_label_name, genes_coords_taxon$taxon_assignment)



# label equals NaN for genes excluded from assignment; replace NaN with "Unassigned"
genes_coords_taxon$plot_label[is.na(genes_coords_taxon$plot_label)] <- "Unassigned"


#count how often each taxon was matched
taxon_count_table <- table(genes_coords_taxon$plot_label) # look in label as these are the groups to be plotted
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
# genes with taxonomic assignemnt of low abundance get plot_label "Otherwise assigned"
genes_coords_taxon$plot_label[genes_coords_taxon$plot_label %in% names(taxon_count_table[taxon_count_table < cfg$taxon_hit_threshold])] <- "Otherwise assigned"


# __________________TEXT OUTPUT________________________ #


# csv output of genes with taxon info and PCA coordinates
num_cov_cols <- length(cfg$pbc_paths)*10 # number of coverage related columns in genes_coords_taxon; information needed for indexing
# indexes of columns to be subsetted and reordered for output
cols <- c(3:(32+num_cov_cols),(34+num_cov_cols),(33+num_cov_cols),2,(35+num_cov_cols),(36+num_cov_cols),(38+num_cov_cols),(39+num_cov_cols),1,(40+num_cov_cols),(42+num_cov_cols))
raw_gene_table_orphan_info <- data.frame(genes_coords_taxon)[,cols]
write.csv(raw_gene_table_orphan_info, file=paste(c(cfg$output_path, "taxonomic_assignment/gene_table_taxon_assignment.csv"), collapse=""), row.names=FALSE, quote=FALSE)
