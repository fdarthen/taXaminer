#!/usr/bin/env Rscript

#' reads output from gene_info and produces PCA_and_clustering
#' @usage Rscript perform_PCA_and_clustering.R config.yml
#' @param config.yml path to preprocessed config file
#' @author Simonida Zehr, Freya Arthen

library(factoextra) # for pca
library(ggplot2) # for plotting
library(viridis) # for colours
library(stats) # for hclust and kmeans clustering
library(yaml) # for reading config file


args <- commandArgs(trailingOnly = TRUE)
config_path <- args[1]

cfg <- yaml.load_file(config_path)

#' save ggplots to pdf or png
#' @param theplot ggplot plot object
#' @param thepath path for plot file
#' @param spec bool whether dimensions and resolution of should be specified
plot_pdf_andor_png <- function(theplot, thepath, spec) {
    # spec: if dimensions and resolution should be specified
    if (spec) {
        if (cfg$output_png) {
            png(paste0(thepath,".png"), units="in", width=3.5, height=3.5, res=500)
            print(theplot)
            dev.off()
        }
        if (cfg$output_pdf) {
            pdf(paste0(thepath,".pdf"), width=3.5, height=3.5)
            print(theplot)
            # theplot
            dev.off()
        }
    }

    # if not (bc e.g. fviz plots may be impaired by this)
    else {
        if (cfg$output_png) {
            png(paste0(thepath,".png"))
            print(theplot)
            dev.off()
        }
        if (cfg$output_pdf) {
            pdf(paste0(thepath,".pdf"))
            print(theplot)
            # theplot
            dev.off()
        }
    }
}

#' create output for clustering algorithm
#' @param coords pca coordinates
#' @param cluster cluster assignment for data points as factor
#' @param cluster_df cluster assignment for data points as data.frame
#' @param clustering_algo which clustering algorithm was performed
#' @param num_cluster number of clusters (eps for DBSCAN)
#' @param num_cluster2 (minPts for DBSCAN)
clustering_output <- function(coords, cluster, cluster_df, clustering_algo, num_cluster, num_cluster2=NULL) {
    if (clustering_algo == "kmeans"){
        plot_title <- paste0("k-means clustering, k=",num_cluster)
        dir_name <- "k-means_clustering"
        sub_name <- paste0("k-means_",num_cluster,"groups")
    }else if (clustering_algo == "hclust"){
        plot_title <- paste0("hierarchical clustering, k=",num_cluster)
        dir_name <- "hierarchical_clustering"
        sub_name <- paste0("hclust_",num_cluster,"groups")
    }else if (clustering_algo == "dbscan"){
        plot_title <- paste0("DB SCAN, eps=",num_cluster,", minPts=",num_cluster2)
        dir_name <- "DBSCAN_clustering"
        sub_name <- paste0("dbscan_eps",num_cluster,"_minPts",num_cluster2)
    }else if (clustering_algo == "mclust"){
        if (is.null(num_cluster)){
            num_cluster <- "opt_" #label when optimal number of groups is chosen
        }
        plot_title <- paste0("model-based clustering, k=",num_cluster)
        dir_name <- "model-based_clustering"
        sub_name <- paste0("mclust_",num_cluster,"groups")
    }else {
        print("Invalid clustering algorithm supplied")
    }

    plot <- ggplot(coords, aes(x=Dim.1, y=Dim.2)) + scale_color_viridis(discrete=TRUE, begin=0.1, end=0.9) + geom_point(aes(shape=cluster, color=cluster)) + labs(color="cluster", shape="cluster") + theme_bw() + ggtitle(plot_title)
    # visual output
    plot_pdf_andor_png(plot, paste0(cfg$output_path, "PCA_and_clustering/",dir_name,"/",sub_name), TRUE)

    # text output
    # merge PCA coordinates and cluster assignments into one data frame
    genecoords <- transform(merge(coords, cluster_df, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # merge original gene data and cluster assignments into one data frame
    raw_gene_table_cluster <- transform(merge(raw_gene_table, genecoords, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)

    # output 1:
    # store g_names as an individual column (for output)
    raw_gene_table_cluster_with_g_names <- cbind(g_name = rownames(raw_gene_table_cluster), raw_gene_table_cluster)
    # and write this info to a .csv file
    write.csv(raw_gene_table_cluster_with_g_names, file=paste0(cfg$output_path, "PCA_and_clustering/",dir_name,"/",sub_name,".taXaminer.csv"), row.names=FALSE, quote=FALSE)

    # output 2:
    # create separate directory to store the files holding all genes by cluster assignment
    if (!dir.exists(paste0(cfg$output_path, "PCA_and_clustering/",dir_name,"/genes_by_cluster/",sub_name))) {
        dir.create(paste0(cfg$output_path, "PCA_and_clustering/",dir_name,"/genes_by_cluster/",sub_name), recursive=TRUE)
    }
    # split the genes by cluster membership
    cluster_lists <- split(rownames(raw_gene_table_cluster), raw_gene_table_cluster$cluster) # --> this results in separate lists, stored in a superlist (cluster_lists)
    # for each list in the super list: print the gene names (row names) to a separate file
    lapply(names(cluster_lists), function(x){write.table(cluster_lists[[x]], file=paste0(cfg$output_path, "PCA_and_clustering/",dir_name,"/genes_by_cluster/",sub_name,"/cluster_", x, ".txt"),col.names=FALSE, row.names=FALSE, quote=FALSE)})
}

#' perform kmeans clustering
#' @param coords pca coordinates
#' @param num_cluster number of clusters
clustering_kmeans <- function(coords, num_cluster) {
    kmeans_clusters <- kmeans(coords, num_cluster)
    cluster <- as.factor(kmeans_clusters$cluster)
    # store cluster assignments in data frame for text output
    cluster_df <- as.data.frame(kmeans_clusters$cluster)
    colnames(cluster_df) <- "cluster"
    clustering_output(coords, cluster, cluster_df, "kmeans", num_cluster)
}

#' perform hierarchical clustering
#' @param coords pca coordinates
#' @param num_cluster number of clusters
clustering_hclust <- function(coords, num_cluster) {
    cluster <- as.factor(cutree(hclust_clustering, k=num_cluster))
    # store cluster assignments in data frame for text output
    cluster_df <- as.data.frame(cutree(hclust_clustering, k=num_cluster))
    colnames(cluster_df) <- "cluster"
    clustering_output(coords, cluster, cluster_df, "hclust", num_cluster)
}

#' perform kmeans clustering
#' @param coords pca coordinates
#' @param eps_val value for eps
#' @param minPts_val value for minPts
clustering_dbscan <- function(coords, eps_val, minPts_val) {
    run <- dbscan(coords, eps=eps_val, minPts=minPts_val)
    cluster <- as.factor(run$cluster+1L)
    # store cluster assignments in data frame for text output
    cluster_df <- as.data.frame(run$cluster+1L)
    colnames(cluster_df) <- "cluster"
    rownames(cluster_df) <- rownames(coords)
    clustering_output(coords, cluster, cluster_df, "dbscan", eps_val, minPts_val)
}

#' perform model based clustering
#' @param coords pca coordinates
#' @param num_cluster number of clusters
clustering_mclust <- function(coords, num_cluster) {
    model <- Mclust(coords, G=num_cluster)
    cluster <- as.factor(model$classification)
    # store cluster assignments in data frame for text output
    cluster_df <- as.data.frame(model$classification)
    colnames(cluster_df) <- "cluster"
    clustering_output(coords, cluster, cluster_df, "mclust", num_cluster)
}

# ::::::::::::::::: DATA PROCESSING ::::::::::::::::::::
rawdata <- read.table(paste0(cfg$output_path,'gene_info/imputed_gene_table.txt'), row.names=1, header=TRUE)
# WITHOUT GC AND ONLY WITH PEARSON R NUCLEOTIDE FREQS
myvars <- unlist(strsplit(cfg$input_variables,",")) # variable to be used for PCA specified by user in config file
# select only the subset of variables
subsetted_data <- subset(rawdata, select=myvars)

# check if all genes are on one contig
if (dim(table(subsetted_data["c_name"])) == 1) { # i.e. is there only one dimension (value) in the column specifying the contig name
    contig_columns <- c("c_num_of_genes", "c_len", "c_pct_assemby_len", "c_cov", "c_covsd", "c_covdev", "c_pearson_r", "c_pearson_p", "c_gc_cont",  "c_gcdev",  "c_genecovm",  "c_genecovsd",  "c_genelenm",  "c_genelensd")
data_columnfiltered1 <- subsetted_data[, !(names(subsetted_data) %in% contig_columns)]
# if so: drop the contig columns
# (since the contig columns have a variance of 0, i.e. contain no information)
} else {
    data_columnfiltered1 <- subsetted_data
}

# delete directory for PCA and clustering results and create new empty one (avoids warning messages)
unlink(paste0(cfg$output_path, "PCA_and_clustering/PCA_results"), recursive=TRUE)

if (!dir.exists(paste0(cfg$output_path, "PCA_and_clustering/PCA_results"))) {
    dir.create(paste0(cfg$output_path, "PCA_and_clustering/PCA_results"), recursive=TRUE)
}

# CLEAN DATA (FROM NaNS) -- COLUMNWISE
# if there are any columns with at least 30% NaN values, save them to a new data frame
columns_with_nans <- data_columnfiltered1[colSums(is.na(data_columnfiltered1))/nrow(data_columnfiltered1) >= .3]
# and print them to a file
write.table(colnames(columns_with_nans), file=paste0(cfg$output_path, "PCA_and_clustering/variables_excluded_from_PCA_and_clustering.txt"),col.names=FALSE, row.names=FALSE, quote=FALSE)
# keep only those columns with less than 30% NaN values
data_columnfiltered2 <- data_columnfiltered1[colSums(is.na(data_columnfiltered1))/nrow(data_columnfiltered1) < .3]

if (dim(table(subsetted_data["c_name"])) == 1) {
    cat("\nexcluded due to a variance of 0 (since only one contig is in the input set):\n", file=paste0(cfg$output_path, "PCA_and_clustering/variables_excluded_from_PCA_and_clustering.txt"), append=TRUE)
    write.table(colnames(subsetted_data[, (names(subsetted_data) %in% contig_columns)]), file=paste0(cfg$output_path, "PCA_and_clustering/variables_excluded_from_PCA_and_clustering.txt"), col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)
}

# CLEAN DATA FROM NaNS -- LISTWISE (i.e. by row, i.e. for gene)
# save the genes containing NaN values to a specific data frame
mydata_nans <- data_columnfiltered2[!complete.cases(data_columnfiltered2),]
# and print them to a file
write.csv(mydata_nans, file=paste0(cfg$output_path, "PCA_and_clustering/genes_excluded_from_PCA_and_clustering.csv"), row.names=TRUE, quote=FALSE)
# keep working only with the genes without NaNs (complete rows / cases)
complete_data <- data_columnfiltered2[complete.cases(data_columnfiltered2),]
# output which genes have been used as an input for the PCA
# write.csv(mydata, file=paste0(cfg$output_path, "PCA_and_clustering/genes_used_for_PCA_and_clustering.csv"), row.names=TRUE, quote=FALSE)
mydata <- complete_data

#TODO: adjust for multiple coverage sets
# FILTER GENES BASED ON COVERAGE
if (cfg$include_coverage){
  if (cfg$coverage_cutoff_mode=="transposons") {
    g_cov_median <- median(complete_data[,"g_cov_1"])
    mydata <- complete_data[complete_data$g_cov_1 >= (g_cov_median),]
    #c_cov_median <- median(complete_data[,"c_cov"])
    #mydata <- g_filtered_data[g_filtered_data$c_cov >= (c_cov_median),]
  } else if (cfg$coverage_cutoff_mode=="contamination") {
    g_cov_median <- median(complete_data[,"g_cov_1"])
    mydata <- complete_data[complete_data$g_cov_1 <= (g_cov_median),]
    #c_cov_median <- median(complete_data[,"c_cov"])
    #mydata <- g_filtered_data[g_filtered_data$c_cov <= (c_cov_median),]
  }
}

# :::::::::::::::::::: PCA ::::::::::::::::::::::::::::::
print("PCA: ")
ptm <- proc.time()
mypca <- prcomp(mydata[2:length(mydata)], scale = TRUE)
print(proc.time() - ptm)
# get new coordinates of genes on ALL principal components
genecoords <- get_pca_ind(mypca)$coord
genecoordsDf <- data.frame(genecoords)

# write summary (standard deviation, explained variance, ...) and loadings to the
write.csv(summary(mypca)$importance, quote=FALSE, file=paste0(cfg$output_path, "PCA_and_clustering/PCA_results/pca_summary.csv"))
write.csv(mypca$rotation, quote=FALSE, file=paste0(cfg$output_path, "PCA_and_clustering/PCA_results/pca_loadings.csv"))


# scree plot
myscree <- fviz_eig(mypca) + labs(title="Scree Plot")
plot_pdf_andor_png(myscree, paste0(cfg$output_path, "PCA_and_clustering/PCA_results/scree_plot"), FALSE)


# contributions of variables
mycontrib <- fviz_pca_var(mypca, col.var = "contrib", gradient.cols = viridis(length(mydata)), repel = TRUE) + labs(title="Contributions of Variables")
plot_pdf_andor_png(mycontrib, paste0(cfg$output_path, "PCA_and_clustering/PCA_results/contribution_of_variables"), FALSE)


# genes and variables
genes_and_vars <- fviz_pca_biplot(mypca, label="var", col.var="black",  alpha.var="contrib", col.ind="skyblue2", alpha.ind="contrib", addEllipses=FALSE)  + labs(title="Genes and Variables in PCA")
plot_pdf_andor_png(genes_and_vars, paste0(cfg$output_path, "PCA_and_clustering/PCA_results/genes_and_variables"), FALSE)


# HORNS PARALLEL ANALYSIS TO EXTRACT NUMBER OF COMPONENTS THAT (TOGETHER) EXPLAIN MORE VARIANCE THAN RANDOM
# perform only if requested
if (cfg$perform_parallel_analysis) {
  print("parallel analysis: ")
  ptm <- proc.time()
  library(paran) # for parallel analysis
  # open documents to store graph of parallel analysis in it
  png(paste0(cfg$output_path, "PCA_and_clustering/PCA_results/parallel_analysis.png"))
  a<-dev.cur()
  pdf(paste0(cfg$output_path, "PCA_and_clustering/PCA_results/parallel_analysis.pdf"))
  dev.control("enable")
  # conduct Horn's parallel analysis
  parallel_analysis <- paran(mydata[2:length(mydata)], centile=95, graph=TRUE)
  # close graph to store it
  dev.copy(which=a)
  dev.off()
  dev.off()
  # get number of components to perform clustering on
  ncomponents <- parallel_analysis$Retained
  # extract ncomponents columns from the original PCA output
  # these will define the new coordinates of the genes
  new_coords <- genecoordsDf[,c(1:ncomponents)]
  print(proc.time() - ptm)
} else { # IF NOT REQUESTED
  # retain 2 dimensions (default)
  new_coords <- genecoordsDf[,c(1:cfg$num_pcs)]
}


# :::::::::::::::::::: CLUSTERING ::::::::::::::::::::::::::::::

# read raw_gene_table to add text output to it
raw_gene_table <- read.table(paste0(cfg$output_path,'gene_info/raw_gene_table.txt'), row.names=1, header=TRUE)

# txt output of raw_gene_table with PCA coordinates
gene_table_coords <- transform(merge(raw_gene_table, new_coords, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
# add gene names as extra column for csv output (csv has per default no header for frist column)
gene_table_coords_w_g_names <- cbind(g_name = rownames(gene_table_coords), gene_table_coords)
write.csv(gene_table_coords_w_g_names, file=paste0(cfg$output_path, "PCA_and_clustering/gene_table_coords.csv"), row.names=FALSE, quote=FALSE)



# _______________K-MEANS CLUSTERING_____________________
# (if requested by the user)
if (cfg$perform_kmeans) {
  print("k means clustering: ")
  ptm <- proc.time()
  # create directory to store the results in
  if (!dir.exists(paste0(cfg$output_path, "PCA_and_clustering/k-means_clustering"))) {
      dir.create(paste0(cfg$output_path, "PCA_and_clustering/k-means_clustering"))
  }
  if (cfg$kmeans_k=="default") {
    # clusterin for default values of k
    clustering_kmeans(new_coords, 2)
    clustering_kmeans(new_coords, 3)
    clustering_kmeans(new_coords, 4)
  }

  # if the user did not choose the "default" output
  # but specified a desired number of clusters (in kmeans_k)
  else {
    kmeans_k <- as.numeric(cfg$kmeans_k)
    clustering_kmeans(new_coords, kmeans_k)
  }
  print(proc.time() - ptm)
}
# end of k-means clustering


# _______________HIERARCHICAL CLUSTERING_____________________
# (if requested by the user)
if (cfg$perform_hclust) {
  print("hierarchical clustering: ")
  ptm <- proc.time()
  # create directory to store the results in
  if (!dir.exists(paste0(cfg$output_path, "PCA_and_clustering/hierarchical_clustering"))) {
      dir.create(paste0(cfg$output_path, "PCA_and_clustering/hierarchical_clustering"))
  }

  d <- dist(new_coords, method="euclidean")
  hclust_clustering <- hclust(d, method="ward.D2")

  if (cfg$hclust_k=="default") {
    # clusterin for default values of k
    clustering_hclust(new_coords, 2)
    clustering_hclust(new_coords, 3)
    clustering_hclust(new_coords, 4)
  }

  # if the user did not choose the "default" output
  # but specified a desired number of clusters (in hclust_k)
  else {
    hclust_k <- as.numeric(cfg$hclust_k)
    clustering_hclust(new_coords, kmeans_k)
  }
  print(proc.time() - ptm)
}
# end of hierarchical clustering


# _______________DBSCAN CLUSTERING_____________________
# (if requested by the user)
if (cfg$perform_dbscan) {
  print("DBSCAN clustering: ")
  ptm <- proc.time()
  # create directory to store the results in
  if (!dir.exists(paste0(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering"))) {
      dir.create(paste0(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering"))
  }
  library(dbscan)

  # IF DEFAULT
  if(cfg$dbscan_groups == "default") {
    # clusterin for default values
    clustering_dbscan(new_coords, 0.15, 5)
    clustering_dbscan(new_coords, 0.20, 5)
    clustering_dbscan(new_coords, 0.20, 10)
    clustering_dbscan(new_coords, 0.25, 10)
    clustering_dbscan(new_coords, 0.25, 5)
    clustering_dbscan(new_coords, 0.30, 5)
  }
  else if ((is.numeric(cfg$custom_eps))&&(is.numeric(cfg$custom_minPts))) {
  clustering_dbscan(new_coords, cfg$custom_eps, cfg$custom_minPts)
  }
  else {
    print("DBSCAN clustering could not be performed. Either choose option ''default'' or specify custom_eps and custom_minPts. E.g. custom_eps <- 0.01, custom_minPts <- 20 (i.e. WITHOUT quote marks!)")
  }
  print(proc.time() - ptm)
}
# end of DBSCAN clustering


# _______________MODEL BASED CLUSTERING_____________________
# (if requested by the user)
if (cfg$perform_mclust) {
  print("model-based clustering: ")
  ptm <- proc.time()
  # create directory to store the results in
  if (!dir.exists(paste0(cfg$output_path, "PCA_and_clustering/model-based_clustering"))) {
      dir.create(paste0(cfg$output_path, "PCA_and_clustering/model-based_clustering"))
  }
  library(mclust)

  if (cfg$mclust_k=="default") {
    # clusterin for default values of k
    clustering_mclust(new_coords, 2)
    clustering_mclust(new_coords, 3)
    clustering_mclust(new_coords, 4)
    # ---- OPTIMAL (number clusters defined by BIC) ----
    clustering_mclust(new_coords, NULL)
  } else if (cfg$mclust_k=="BIC") {
    # ---- OPTIMAL (number clusters defined by BIC) ----
    clustering_mclust(new_coords, NULL)
  }
  # if the user did not choose the "default" output
  # but specified a desired number of clusters (in mclust_k)
  else {
    mclust_k <- as.numeric(cfg$mclust_k)
    clustering_mclust(new_coords, mclust_k)
  }
  print(proc.time() - ptm)
}
# end of model based clustering
