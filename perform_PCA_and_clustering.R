#!/usr/bin/env Rscript

# reads output from gene_info and produces PCA_and_clustering
# author: Simonida Zehr
# date: 18 June 2019


library(factoextra) # for pca
library(ggplot2) # for plotting
library(viridis) # for colours
library(stats) # for hclust and kmeans clustering


library(yaml) # for reading config file
cfg <- yaml.load_file("config.yml")

plot_pdf_andor_png <- function(theplot, thepath, spec) {
  # spec: if dimensions and resolution should be specified
  if (spec) {
    if (cfg$output_png) {
      png(paste(c(thepath,".png"), collapse=""), units="in", width=3.5, height=3.5, res=500)
      print(theplot)
      dev.off()
    }
    if (cfg$output_pdf) {
      pdf(paste(c(thepath,".pdf"), collapse=""), width=3.5, height=3.5)
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



# ::::::::::::::::: DATA PROCESSING ::::::::::::::::::::
rawdata <- read.table(paste(c(cfg$output_path,'gene_info/imputed_gene_table.txt'), collapse=''), row.names=1, header=TRUE)
# WITHOUT GC AND ONLY WITH PEARSON R NUCLEOTIDE FREQS

myvars <- unlist(strsplit(cfg$input_variables,",")) # variable to be used for PCA specified by user in config file
if (cfg$include_coverage){
  cov_vars <- c("c_cov","c_covsd","g_cov","g_covsd","g_covdev_c") # coverage variables
  header_vars <- colnames(rawdata) # variables in input file
  my_cov_vars <- grep(paste(cov_vars,collapse="|"),myvars, value=TRUE) # extract which of the coverage variables shall be used (selected by user in config)
  all_my_cov_vars <- grep(paste(my_cov_vars,collapse="_*[0-9]|"),header_vars, value=TRUE) #extract variables matching my_cov_vars from the header of rawdata -> '_X' appended
  myvars <- myvars[!myvars %in% my_cov_vars] # remove un-indexed cov variables from the user selection in config
  myvars <- c(myvars,all_my_cov_vars) # merge in config file selected variables with extracted coverage variables
}else {
  cov_vars <- c("c_cov","c_covsd","g_cov","g_covsd","g_covdev_c","g_cov_z_bool")
  myvars <- myvars[!myvars %in% cov_vars]
}
print(myvars)
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
unlink(paste(c(cfg$output_path, "PCA_and_clustering/PCA_results"), collapse=""), recursive=TRUE)
dir.create(paste(c(cfg$output_path, "PCA_and_clustering/PCA_results"), collapse=""), recursive=TRUE)



# CLEAN DATA (FROM NaNS) -- COLUMNWISE
# if there are any columns with at least 30% NaN values, save them to a new data frame
columns_with_nans <- data_columnfiltered1[colSums(is.na(data_columnfiltered1))/nrow(data_columnfiltered1) >= .3]
# and print them to a file
write.table(colnames(columns_with_nans), file=paste(c(cfg$output_path, "PCA_and_clustering/variables_excluded_from_PCA_and_clustering.txt"), collapse=""),col.names=FALSE, row.names=FALSE, quote=FALSE)
# keep only those columns with less than 30% NaN values
data_columnfiltered2 <- data_columnfiltered1[colSums(is.na(data_columnfiltered1))/nrow(data_columnfiltered1) < .3]


if (dim(table(subsetted_data["c_name"])) == 1) {
    cat("\nexcluded due to a variance of 0 (since only one contig is in the input set):\n", file=paste(c(cfg$output_path, "PCA_and_clustering/variables_excluded_from_PCA_and_clustering.txt"), collapse=""), append=TRUE)
    write.table(colnames(subsetted_data[, (names(subsetted_data) %in% contig_columns)]), file=paste(c(cfg$output_path, "PCA_and_clustering/variables_excluded_from_PCA_and_clustering.txt"), collapse=""), col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)
}


# CLEAN DATA FROM NaNS -- LISTWISE (i.e. by row, i.e. for gene)
# save the genes containing NaN values to a specific data frame
mydata_nans <- data_columnfiltered2[!complete.cases(data_columnfiltered2),]
# and print them to a file
write.csv(mydata_nans, file=paste(c(cfg$output_path, "PCA_and_clustering/genes_excluded_from_PCA_and_clustering.csv"), collapse=""), row.names=TRUE, quote=FALSE)
# keep working only with the genes without NaNs (complete rows / cases)
complete_data <- data_columnfiltered2[complete.cases(data_columnfiltered2),]
# output which genes have been used as an input for the PCA
# write.csv(mydata, file=paste(c(cfg$output_path, "PCA_and_clustering/genes_used_for_PCA_and_clustering.csv"), collapse=""), row.names=TRUE, quote=FALSE)



# FILTER GENES BASED ON COVERAGE
if (cfg$include_coverage){
  if (cfg$coverage_cutoff_mode=="transposons") {
    g_cov_median <- mean(complete_data[,"g_cov_0"])
    mydata <- complete_data[complete_data$g_cov_0 >= (g_cov_median),]
    #c_cov_median <- median(complete_data[,"c_cov"])
    #mydata <- g_filtered_data[g_filtered_data$c_cov >= (c_cov_median),]
  } else if (cfg$coverage_cutoff_mode=="contamination") {
    g_cov_median <- mean(complete_data[,"g_cov_0"])
    mydata <- complete_data[complete_data$g_cov_0 <= (g_cov_median),]
    #c_cov_median <- median(complete_data[,"c_cov"])
    #mydata <- g_filtered_data[g_filtered_data$c_cov <= (c_cov_median),]
  } else {
    mydata <- complete_data
  }
  if (dim(table(mydata["g_cov_z_bool"])) == 1) {
      #cat("\nexcluded due to a variance of 0 (no genes with conspicuos coverage profile):\n", file=paste(c(cfg$output_path, "PCA_and_clustering/variables_excluded_from_PCA_and_clustering.txt"), collapse=""), append=TRUE)
      #write.table(colnames(subsetted_data[, (names(subsetted_data) %in% contig_columns)]), file=paste(c(cfg$output_path, "PCA_and_clustering/variables_excluded_from_PCA_and_clustering.txt"), collapse=""), col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)
      mydata <- subset(mydata, select = -g_cov_z_bool)
}


# :::::::::::::::::::: PCA ::::::::::::::::::::::::::::::
print("PCA: ")
ptm <- proc.time()
mypca <- prcomp(mydata[2:length(mydata)], scale = TRUE)
print(proc.time() - ptm)
# get new coordinates of genes on ALL principal components
#explor(mypca) #FH
genecoords <- get_pca_ind(mypca)$coord
genecoordsDf <- data.frame(genecoords)

# write summary (standard deviation, explained variance, ...) and loadings to the
write.csv(summary(mypca)$importance, quote=FALSE, file=paste(c(cfg$output_path, "PCA_and_clustering/PCA_results/pca_summary.csv"), collapse=""))
write.csv(mypca$rotation, quote=FALSE, file=paste(c(cfg$output_path, "PCA_and_clustering/PCA_results/pca_loadings.csv"), collapse=""))


# scree plot
myscree <- fviz_eig(mypca) + labs(title="Scree Plot")
plot_pdf_andor_png(myscree, paste(c(cfg$output_path, "PCA_and_clustering/PCA_results/scree_plot"), collapse=""), FALSE)


# contributions of variables
mycontrib <- fviz_pca_var(mypca, col.var = "contrib", gradient.cols = viridis(length(mydata)), repel = TRUE) + labs(title="Contributions of Variables")
plot_pdf_andor_png(mycontrib, paste(c(cfg$output_path, "PCA_and_clustering/PCA_results/contribution_of_variables"), collapse=""), FALSE)


# genes and variables
genes_and_vars <- fviz_pca_biplot(mypca, label="var", col.var="black",  alpha.var="contrib", col.ind="skyblue2", alpha.ind="contrib", addEllipses=FALSE)  + labs(title="Genes and Variables in PCA")
plot_pdf_andor_png(genes_and_vars, paste(c(cfg$output_path, "PCA_and_clustering/PCA_results/genes_and_variables"), collapse=""), FALSE)



# HORNS PARALLEL ANALYSIS TO EXTRACT NUMBER OF COMPONENTS THAT (TOGETHER) EXPLAIN MORE VARIANCE THAN RANDOM
# perform only if requested
if (cfg$perform_parallel_analysis) {
  print("parallel analysis: ")
  ptm <- proc.time()
  library(paran) # for parallel analysis
  # open documents to store graph of parallel analysis in it
  png(paste(c(cfg$output_path, "PCA_and_clustering/PCA_results/parallel_analysis.png"), collapse=""))
  a<-dev.cur()
  pdf(paste(c(cfg$output_path, "PCA_and_clustering/PCA_results/parallel_analysis.pdf"), collapse=""))
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
raw_gene_table <- read.table(paste(c(cfg$output_path,'gene_info/raw_gene_table.txt'), collapse=''), row.names=1, header=TRUE)

# txt output of raw_gene_table with PCA coordinates
gene_table_coords <- transform(merge(raw_gene_table, new_coords, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
# add gene names as extra column for csv output (csv has per default no header for frist column)
gene_table_coords_w_g_names <- cbind(g_name = rownames(gene_table_coords), gene_table_coords)
write.csv(gene_table_coords_w_g_names, file=paste(c(cfg$output_path, "PCA_and_clustering/gene_table_coords.csv"), collapse=""), row.names=FALSE, quote=FALSE)



# _______________K-MEANS CLUSTERING_____________________
# (if requested by the user)
if (cfg$perform_kmeans) {
  print("k means clustering: ")
  ptm <- proc.time()
  # create directory to store the results in
  dir.create(paste(c(cfg$output_path, "PCA_and_clustering/k-means_clustering"), collapse=""))
  if (cfg$kmeans_k=="default") {

    # ---- 2 CLUSTERS ----
    kmeans_2clusters <- kmeans(new_coords, 2)
    # store cluster assignments in a data frame
    cluster_2 <- as.factor(kmeans_2clusters$cluster)
    plot_2 <- ggplot(new_coords, aes(x=Dim.1, y=Dim.2))+ scale_color_viridis(discrete=TRUE, begin=0.1, end=0.9) + geom_point(aes(shape=cluster_2, color=cluster_2)) + labs(color="cluster", shape="cluster") + theme_bw() + ggtitle("k-means clustering, k=2")
    # visual output
    plot_pdf_andor_png(plot_2, paste(c(cfg$output_path, "PCA_and_clustering/k-means_clustering/k-means_2groups"), collapse=""), TRUE)


    # text output
    # store cluster assignments in data frame
    cluster_2 <- as.data.frame(kmeans_2clusters$cluster)
    colnames(cluster_2) <- "cluster"
    # merge PCA coordinates and cluster assignments into one data frame
    genecoords_kmeans2 <- transform(merge(new_coords, cluster_2, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # merge original gene data and cluster assignments into one data frame
    raw_gene_table_cluster_kmeans2 <- transform(merge(raw_gene_table, genecoords_kmeans2, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)


    # output 1:
    # store g_names as an individual column (for output)
    raw_gene_table_cluster_with_g_names_kmeans2 <- cbind(g_name = rownames(raw_gene_table_cluster_kmeans2), raw_gene_table_cluster_kmeans2)
    # and write this info to a .csv file
    write.csv(raw_gene_table_cluster_with_g_names_kmeans2, file=paste(c(cfg$output_path, "PCA_and_clustering/k-means_clustering/k-means_2groups.MILTS.csv"), collapse=""), row.names=FALSE, quote=FALSE)

    # output 2:
    # create separate directory to store the files holding all genes by cluster assignment
    dir.create(paste(c(cfg$output_path, "PCA_and_clustering/k-means_clustering/genes_by_cluster/k-means_2"), collapse=""), recursive=TRUE)
    # split the genes by cluster membership
    cluster_lists <- split(rownames(raw_gene_table_cluster_kmeans2), raw_gene_table_cluster_kmeans2$cluster) # --> this results in separate lists, stored in a superlist (cluster_lists)
    # for each list in the super list: print the gene names (row names) to a separate file
    lapply(names(cluster_lists), function(x){write.table(cluster_lists[[x]], file =paste(c(cfg$output_path, "PCA_and_clustering/k-means_clustering/genes_by_cluster/k-means_2/cluster_", x, ".txt"), collapse=""),col.names=FALSE, row.names=FALSE, quote=FALSE)})


    # ---- 3 CLUSTERS ----
    kmeans_3clusters <- kmeans(new_coords, 3)
    cluster_3 <- as.factor(kmeans_3clusters$cluster)
    plot_3 <- ggplot(new_coords, aes(x=Dim.1, y=Dim.2))+ scale_color_viridis(discrete=TRUE, begin=0.1, end=0.9) + geom_point(aes(shape=cluster_3, color=cluster_3)) + labs(color="cluster", shape="cluster") + theme_bw() + ggtitle("k-means clustering, k=3")

    # visual output
    plot_pdf_andor_png(plot_3, paste(c(cfg$output_path, "PCA_and_clustering/k-means_clustering/k-means_3groups"), collapse=""), TRUE)

    # text output
    # store cluster assignments in data frame
    cluster_3 <- as.data.frame(kmeans_3clusters$cluster)
    colnames(cluster_3) <- "cluster"
    # merge PCA coordinates and cluster assignments into one data frame
    genecoords_kmeans3 <- transform(merge(new_coords, cluster_3, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # merge original gene data and cluster assignments into one data frame
    raw_gene_table_cluster_kmeans3 <- transform(merge(raw_gene_table, genecoords_kmeans3, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)

    # output 1:
    # store g_names as an individual column (for output)
    raw_gene_table_cluster_with_g_names_kmeans3 <- cbind(g_name = rownames(raw_gene_table_cluster_kmeans3), raw_gene_table_cluster_kmeans3)
    # and write this info to a .csv file
    write.csv(raw_gene_table_cluster_with_g_names_kmeans3, file=paste(c(cfg$output_path, "PCA_and_clustering/k-means_clustering/k-means_3groups.MILTS.csv"), collapse=""), row.names=FALSE, quote=FALSE)

    # output 2:
    # create separate directory to store the files holding all genes by cluster assignment
    dir.create(paste(c(cfg$output_path, "PCA_and_clustering/k-means_clustering/genes_by_cluster/k-means_3"), collapse=""), recursive=TRUE)
    # split the genes by cluster membership
    cluster_lists <- split(rownames(raw_gene_table_cluster_kmeans3), raw_gene_table_cluster_kmeans3$cluster) # --> this results in separate lists, stored in a superlist (cluster_lists)
    # for each list in the super list: print the gene names (row names) to a separate file
    lapply(names(cluster_lists), function(x){write.table(cluster_lists[[x]], file =paste(c(cfg$output_path, "PCA_and_clustering/k-means_clustering/genes_by_cluster/k-means_3/cluster_", x, ".txt"), collapse=""),col.names=FALSE, row.names=FALSE, quote=FALSE)})

    # ---- 4 CLUSTERS ----
    kmeans_4clusters <- kmeans(new_coords, 4)
    cluster_4 <- as.factor(kmeans_4clusters$cluster)
    plot_4 <- ggplot(new_coords, aes(x=Dim.1, y=Dim.2))+ scale_color_viridis(discrete=TRUE, begin=0.1, end=0.9) + geom_point(aes(shape=cluster_4, color=cluster_4)) + labs(color="cluster", shape="cluster") + theme_bw() + ggtitle("k-means clustering, k=4")

    # visual output
    plot_pdf_andor_png(plot_4, paste(c(cfg$output_path, "PCA_and_clustering/k-means_clustering/k-means_4groups"), collapse=""), TRUE)

    # text output
    # store cluster assignments in data frame
    cluster_4 <- as.data.frame(kmeans_4clusters$cluster)
    colnames(cluster_4) <- "cluster"
    # merge PCA coordinates and cluster assignments into one data frame
    genecoords_kmeans4 <- transform(merge(new_coords, cluster_4, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # merge original gene data and cluster assignments into one data frame
    raw_gene_table_cluster_kmeans4 <- transform(merge(raw_gene_table, genecoords_kmeans4, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)

    # output 1:
    # store g_names as an individual column (for output)
    raw_gene_table_cluster_with_g_names_kmeans4 <- cbind(g_name = rownames(raw_gene_table_cluster_kmeans4), raw_gene_table_cluster_kmeans4)
    # and write this info to a .csv file
    write.csv(raw_gene_table_cluster_with_g_names_kmeans4, file=paste(c(cfg$output_path, "PCA_and_clustering/k-means_clustering/k-means_4groups.MILTS.csv"), collapse=""), row.names=FALSE, quote=FALSE)

    # output 2:
    # create separate directory to store the files holding all genes by cluster assignment
    dir.create(paste(c(cfg$output_path, "PCA_and_clustering/k-means_clustering/genes_by_cluster/k-means_4"), collapse=""), recursive=TRUE)
    # split the genes by cluster membership
    cluster_lists <- split(rownames(raw_gene_table_cluster_kmeans4), raw_gene_table_cluster_kmeans4$cluster) # --> this results in separate lists, stored in a superlist (cluster_lists)
    # for each list in the super list: print the gene names (row names) to a separate file
    lapply(names(cluster_lists), function(x){write.table(cluster_lists[[x]], file =paste(c(cfg$output_path, "PCA_and_clustering/k-means_clustering/genes_by_cluster/k-means_4/cluster_", x, ".txt"), collapse=""),col.names=FALSE, row.names=FALSE, quote=FALSE)})


  }

  # if the user did not choose the "default" output
  # but specified a desired number of clusters (in kmeans_k)
  else {
    kmeans_k <- as.numeric(cfg$kmeans_k)

    kmeans_kclusters <- kmeans(new_coords, kmeans_k)
    cluster_k <- as.factor(kmeans_kclusters$cluster)
    plot_k <- ggplot(new_coords, aes(x=Dim.1, y=Dim.2))+ scale_color_viridis(discrete=TRUE, begin=0.1, end=0.9) + geom_point(aes(shape=cluster_k, color=cluster_k)) + labs(color="cluster", shape="cluster") + theme_bw() + ggtitle(paste(c("k-means clustering, k=",kmeans_k), collapse=""))

    # visual output
    plot_pdf_andor_png(plot_k, paste(c(cfg$output_path,"PCA_and_clustering/k-means_clustering/k-means_",kmeans_k, "groups"), collapse=""), TRUE)

    # text output
    # store cluster assignments in data frame
    cluster_k <- as.data.frame(kmeans_kclusters$cluster)
    colnames(cluster_k) <- "cluster"
    # merge PCA coordinates and cluster assignments into one data frame
    genecoords_kmeansk <- transform(merge(new_coords, cluster_k, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # merge original gene data and cluster assignments into one data frame
    raw_gene_table_cluster_kmeansk <- transform(merge(raw_gene_table, genecoords_kmeansk, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)

    # output 1:
    # store g_names as an individual column (for output)
    raw_gene_table_cluster_with_g_names_kmeansk <- cbind(g_name = rownames(raw_gene_table_cluster_kmeansk), raw_gene_table_cluster_kmeansk)
    # and write this info to a .csv file
    write.csv(raw_gene_table_cluster_with_g_names_kmeansk, file=paste(c(cfg$output_path, "PCA_and_clustering/k-means_clustering/k-means_", kmeans_k, "groups.MILTS.csv"), collapse=""), row.names=FALSE, quote=FALSE)

    # output 2:
    # create separate directory to store the files holding all genes by cluster assignment
    dir.create(paste(c(cfg$output_path, "PCA_and_clustering/k-means_clustering/genes_by_cluster/k-means_", kmeans_k), collapse=""), recursive=TRUE)
    # split the genes by cluster membership
    cluster_lists <- split(rownames(raw_gene_table_cluster_kmeansk), raw_gene_table_cluster_kmeansk$cluster) # --> this results in separate lists, stored in a superlist (cluster_lists)
    # for each list in the super list: print the gene names (row names) to a separate file
    lapply(names(cluster_lists), function(x){write.table(cluster_lists[[x]], file =paste(c(cfg$output_path, "PCA_and_clustering/k-means_clustering/genes_by_cluster/k-means_", kmeans_k, "/cluster_", x, ".txt"), collapse=""),col.names=FALSE, row.names=FALSE, quote=FALSE)})

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
  dir.create(paste(c(cfg$output_path, "PCA_and_clustering/hierarchical_clustering"), collapse=""))

  d <- dist(new_coords, method="euclidean")
  hclust_clustering <- hclust(d, method="ward.D2")

  if (cfg$hclust_k=="default") {

    # ---- 2 CLUSTERS ----
    cluster_2 <- as.factor(cutree(hclust_clustering, k=2))
    plot_2 <- ggplot(new_coords, aes(x=Dim.1, y=Dim.2))+ scale_color_viridis(discrete=TRUE, begin=0.1, end=0.9) + geom_point(aes(shape=cluster_2, color=cluster_2)) + labs(color="cluster", shape="cluster") + theme_bw() + ggtitle("hierarchical clustering, k=2")

    # visual output
    plot_pdf_andor_png(plot_2, paste(c(cfg$output_path, "PCA_and_clustering/hierarchical_clustering/hclust_2groups"), collapse=""), TRUE)
    # text output
    # store cluster assignments in data frame
    cluster_2 <- as.data.frame(cutree(hclust_clustering, k=2))
    colnames(cluster_2) <- "cluster"
    # merge PCA coordinates and cluster assignments into one data frame
    genecoords_hclust2 <- transform(merge(new_coords, cluster_2, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # merge original gene data and cluster assignments into one data frame
    raw_gene_table_cluster_hclust2 <- transform(merge(raw_gene_table, genecoords_hclust2, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # output 1:
    # store g_names as an individual column (for output)
    raw_gene_table_cluster_with_g_names_hclust2 <- cbind(g_name = rownames(raw_gene_table_cluster_hclust2), raw_gene_table_cluster_hclust2)
    # and write this info to a .csv file
    write.csv(raw_gene_table_cluster_with_g_names_hclust2, file=paste(c(cfg$output_path, "PCA_and_clustering/hierarchical_clustering/hclust_2groups.MILTS.csv"), collapse=""), row.names=FALSE, quote=FALSE)

    # output 2:
    # create separate directory to store the files holding all genes by cluster assignment
    dir.create(paste(c(cfg$output_path, "PCA_and_clustering/hierarchical_clustering/genes_by_cluster/hclust_2"), collapse=""), recursive=TRUE)
    # split the genes by cluster membership
    cluster_lists <- split(rownames(raw_gene_table_cluster_hclust2), raw_gene_table_cluster_hclust2$cluster) # --> this results in separate lists, stored in a superlist (cluster_lists)
    # for each list in the super list: print the gene names (row names) to a separate file
    lapply(names(cluster_lists), function(x){write.table(cluster_lists[[x]], file =paste(c(cfg$output_path, "PCA_and_clustering/hierarchical_clustering/genes_by_cluster/hclust_2/cluster_", x, ".txt"), collapse=""),col.names=FALSE, row.names=FALSE, quote=FALSE)})


    # ---- 3 CLUSTERS ----
    cluster_3 <- as.factor(cutree(hclust_clustering, k=3))
    plot_3 <- ggplot(new_coords, aes(x=Dim.1, y=Dim.2))+ scale_color_viridis(discrete=TRUE, begin=0.1, end=0.9) + geom_point(aes(shape=cluster_3, color=cluster_3)) + labs(color="cluster", shape="cluster") + theme_bw() + ggtitle("hierarchical clustering, k=3")

    # visual output
    plot_pdf_andor_png(plot_3, paste(c(cfg$output_path, "PCA_and_clustering/hierarchical_clustering/hclust_3groups"), collapse=""), TRUE)
    # text output
    # store cluster assignments in data frame
    cluster_3 <- as.data.frame(cutree(hclust_clustering, k=3))
    colnames(cluster_3) <- "cluster"
    # merge PCA coordinates and cluster assignments into one data frame
    genecoords_hclust3 <- transform(merge(new_coords, cluster_3, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)

    # merge original gene data and cluster assignments into one data frame
    raw_gene_table_cluster_hclust3 <- transform(merge(raw_gene_table, genecoords_hclust3, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # output 1:
    # store g_names as an individual column (for output)
    raw_gene_table_cluster_with_g_names_hclust3 <- cbind(g_name = rownames(raw_gene_table_cluster_hclust3), raw_gene_table_cluster_hclust3)
    # and write this info to a .csv file
    write.csv(raw_gene_table_cluster_with_g_names_hclust3, file=paste(c(cfg$output_path, "PCA_and_clustering/hierarchical_clustering/hclust_3groups.MILTS.csv"), collapse=""), row.names=FALSE, quote=FALSE)

    # output 2:
    # create separate directory to store the files holding all genes by cluster assignment
    dir.create(paste(c(cfg$output_path, "PCA_and_clustering/hierarchical_clustering/genes_by_cluster/hclust_3"), collapse=""), recursive=TRUE)
    # split the genes by cluster membership
    cluster_lists <- split(rownames(raw_gene_table_cluster_hclust3), raw_gene_table_cluster_hclust3$cluster) # --> this results in separate lists, stored in a superlist (cluster_lists)
    # for each list in the super list: print the gene names (row names) to a separate file
    lapply(names(cluster_lists), function(x){write.table(cluster_lists[[x]], file =paste(c(cfg$output_path, "PCA_and_clustering/hierarchical_clustering/genes_by_cluster/hclust_3/cluster_", x, ".txt"), collapse=""),col.names=FALSE, row.names=FALSE, quote=FALSE)})

    # ---- 4 CLUSTERS ----
    cluster_4 <- as.factor(cutree(hclust_clustering, k=4))
    plot_4 <- ggplot(new_coords, aes(x=Dim.1, y=Dim.2))+ scale_color_viridis(discrete=TRUE, begin=0.1, end=0.9) + geom_point(aes(shape=cluster_4, color=cluster_4)) + labs(color="cluster", shape="cluster") + theme_bw() + ggtitle("hierarchical clustering, k=4")

    # visual output
    plot_pdf_andor_png(plot_4, paste(c(cfg$output_path, "PCA_and_clustering/hierarchical_clustering/hclust_4groups"), collapse=""), TRUE)
    # text output
    # store cluster assignments in data frame
    cluster_4 <- as.data.frame(cutree(hclust_clustering, k=4))
    colnames(cluster_4) <- "cluster"
    # merge PCA coordinates and cluster assignments into one data frame
    genecoords_hclust4 <- transform(merge(new_coords, cluster_4, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)

    # merge original gene data and cluster assignments into one data frame
    raw_gene_table_cluster_hclust4 <- transform(merge(raw_gene_table, genecoords_hclust4, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # output 1:
    # store g_names as an individual column (for output)
    raw_gene_table_cluster_with_g_names_hclust4 <- cbind(g_name = rownames(raw_gene_table_cluster_hclust4), raw_gene_table_cluster_hclust4)
    # and write this info to a .csv file
    write.csv(raw_gene_table_cluster_with_g_names_hclust4, file=paste(c(cfg$output_path, "PCA_and_clustering/hierarchical_clustering/hclust_4groups.MILTS.csv"), collapse=""), row.names=FALSE, quote=FALSE)

    # output 2:
    # create separate directory to store the files holding all genes by cluster assignment
    dir.create(paste(c(cfg$output_path, "PCA_and_clustering/hierarchical_clustering/genes_by_cluster/hclust_4"), collapse=""), recursive=TRUE)
    # split the genes by cluster membership
    cluster_lists <- split(rownames(raw_gene_table_cluster_hclust4), raw_gene_table_cluster_hclust4$cluster) # --> this results in separate lists, stored in a superlist (cluster_lists)
    # for each list in the super list: print the gene names (row names) to a separate file
    lapply(names(cluster_lists), function(x){write.table(cluster_lists[[x]], file =paste(c(cfg$output_path, "PCA_and_clustering/hierarchical_clustering/genes_by_cluster/hclust_4/cluster_", x, ".txt"), collapse=""),col.names=FALSE, row.names=FALSE, quote=FALSE)})


  }

  # if the user did not choose the "default" output
  # but specified a desired number of clusters (in hclust_k)
  else {
    hclust_k <- as.numeric(cfg$hclust_k)

    cluster_k <- as.factor(cutree(hclust_clustering, k=hclust_k))
    plot_k <- ggplot(new_coords, aes(x=Dim.1, y=Dim.2))+ scale_color_viridis(discrete=TRUE, begin=0.1, end=0.9) + geom_point(aes(shape=cluster_k, color=cluster_k)) + labs(color="cluster", shape="cluster") + theme_bw() + ggtitle(paste(c("hierarchical clustering, k=",hclust_k), collapse=""))

    # visual output
    plot_pdf_andor_png(plot_k, paste(c(cfg$output_path,"PCA_and_clustering/hierarchical_clustering/hclust_",hclust_k, "groups"), collapse=""), TRUE)
    # text output
    # store cluster assignments in data frame
    cluster_k <- as.data.frame(cutree(hclust_clustering, k=hclust_k))
    colnames(cluster_k) <- "cluster"
    # merge PCA coordinates and cluster assignments into one data frame
    genecoords_hclustk <- transform(merge(new_coords, cluster_k, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # merge original gene data and cluster assignments into one data frame
    raw_gene_table_cluster_hclustk <- transform(merge(raw_gene_table, genecoords_hclustk, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # output 1:
    # store g_names as an individual column (for output)
    raw_gene_table_cluster_with_g_names_hclustk <- cbind(g_name = rownames(raw_gene_table_cluster_hclustk), raw_gene_table_cluster_hclustk)
    # and write this info to a .csv file
    write.csv(raw_gene_table_cluster_with_g_names_hclustk, file=paste(c(cfg$output_path, "PCA_and_clustering/hierarchical_clustering/hclust_", hclust_k, "groups.MILTS.csv"), collapse=""), row.names=FALSE, quote=FALSE)

    # output 2:
    # create separate directory to store the files holding all genes by cluster assignment
    dir.create(paste(c(cfg$output_path, "PCA_and_clustering/hierarchical_clustering/genes_by_cluster/hclust_", hclust_k), collapse=""), recursive=TRUE)
    # split the genes by cluster membership
    cluster_lists <- split(rownames(raw_gene_table_cluster_hclustk), raw_gene_table_cluster_hclustk$cluster) # --> this results in separate lists, stored in a superlist (cluster_lists)
    # for each list in the super list: print the gene names (row names) to a separate file
    lapply(names(cluster_lists), function(x){write.table(cluster_lists[[x]], file =paste(c(cfg$output_path, "PCA_and_clustering/hierarchical_clustering/genes_by_cluster/hclust_", hclust_k, "/cluster_", x, ".txt"), collapse=""),col.names=FALSE, row.names=FALSE, quote=FALSE)})

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
  dir.create(paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering"), collapse=""))
  library(dbscan)

  # IF DEFAULT
  if(cfg$dbscan_groups == "default") {
    run_a <- dbscan(new_coords, eps=0.15, minPts=5)
    cluster_a <- as.factor(run_a$cluster+1L)
    plot_a <- ggplot(new_coords, aes(x=Dim.1, y=Dim.2))+ geom_point(aes(shape=cluster_a, color=cluster_a)) + labs(color="cluster", shape="cluster") + theme_bw() + ggtitle("DB SCAN, eps=0.15, minPts=5")
    # visual output
    plot_pdf_andor_png(plot_a, paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/dbscan_eps0.15_minPts5"), collapse=""), TRUE)
    # text output
    # store cluster assignments in data frame
    cluster_a <- as.data.frame(run_a$cluster+1L)
    colnames(cluster_a) <- "cluster"
    rownames(cluster_a) <- rownames(new_coords)
    # merge PCA coordinates and cluster assignments into one data frame
    genecoords_dbscan_a <- transform(merge(new_coords, cluster_a, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # merge original gene data and cluster assignments into one data frame
    raw_gene_table_cluster_dbscan_a <- transform(merge(raw_gene_table, genecoords_dbscan_a, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # output 1:
    # store g_names as an individual column (for output)
    raw_gene_table_cluster_with_g_names_dbscan_a <- cbind(g_name = rownames(raw_gene_table_cluster_dbscan_a), raw_gene_table_cluster_dbscan_a)
    # and write this info to a .csv file
    write.csv(raw_gene_table_cluster_with_g_names_dbscan_a, file=paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/dbscan_eps0.15_minPts5.MILTS.csv"), collapse=""), row.names=FALSE, quote=FALSE)

    # output 2:
    # create separate directory to store the files holding all genes by cluster assignment
    dir.create(paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/genes_by_cluster/dbscan_eps0.15_minPts5"), collapse=""), recursive=TRUE)
    # split the genes by cluster membership
    cluster_lists <- split(rownames(raw_gene_table_cluster_dbscan_a), raw_gene_table_cluster_dbscan_a$cluster) # --> this results in separate lists, stored in a superlist (cluster_lists)
    # for each list in the super list: print the gene names (row names) to a separate file
    lapply(names(cluster_lists), function(x){write.table(cluster_lists[[x]], file =paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/genes_by_cluster/dbscan_eps0.15_minPts5/cluster_", x, ".txt"), collapse=""),col.names=FALSE, row.names=FALSE, quote=FALSE)})



    run_b <- dbscan(new_coords, eps=0.20, minPts=5)
    cluster_b <- as.factor(run_b$cluster+1L)
    plot_b <- ggplot(new_coords, aes(x=Dim.1, y=Dim.2))+ geom_point(aes(shape=cluster_b, color=cluster_b)) + labs(color="cluster", shape="cluster") + theme_bw() + ggtitle("DB SCAN, eps=0.20, minPts=5")
    # visual output
    plot_pdf_andor_png(plot_b, paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/dbscan_eps0.20_minPts5"), collapse=""), TRUE)
    # text output
    # store cluster assignments in data frame
    cluster_b <- as.data.frame(run_b$cluster+1L)
    colnames(cluster_b) <- "cluster"
    rownames(cluster_b) <- rownames(new_coords)
    # merge PCA coordinates and cluster assignments into one data frame
    genecoords_dbscan_b <- transform(merge(new_coords, cluster_b, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # merge original gene data and cluster assignments into one data frame
    raw_gene_table_cluster_dbscan_b <- transform(merge(raw_gene_table, genecoords_dbscan_b, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # output 1:
    # store g_names as an individual column (for output)
    raw_gene_table_cluster_with_g_names_dbscan_b <- cbind(g_name = rownames(raw_gene_table_cluster_dbscan_b), raw_gene_table_cluster_dbscan_b)
    # and write this info to a .csv file
    write.csv(raw_gene_table_cluster_with_g_names_dbscan_b, file=paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/dbscan_eps0.20_minPts5.MILTS.csv"), collapse=""), row.names=FALSE, quote=FALSE)

    # output 2:
    # create separate directory to store the files holding all genes by cluster assignment
    dir.create(paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/genes_by_cluster/dbscan_eps0.20_minPts5"), collapse=""), recursive=TRUE)
    # split the genes by cluster membership
    cluster_lists <- split(rownames(raw_gene_table_cluster_dbscan_b), raw_gene_table_cluster_dbscan_b$cluster) # --> this results in separate lists, stored in a superlist (cluster_lists)
    # for each list in the super list: print the gene names (row names) to a separate file
    lapply(names(cluster_lists), function(x){write.table(cluster_lists[[x]], file =paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/genes_by_cluster/dbscan_eps0.20_minPts5/cluster_", x, ".txt"), collapse=""),col.names=FALSE, row.names=FALSE, quote=FALSE)})



    run_c <- dbscan(new_coords, eps=0.20, minPts=10)
    cluster_c <- as.factor(run_c$cluster+1L)
    plot_c <- ggplot(new_coords, aes(x=Dim.1, y=Dim.2))+ geom_point(aes(shape=cluster_c, color=cluster_c)) + labs(color="cluster", shape="cluster") + theme_bw() + ggtitle("DB SCAN, eps=0.20, minPts=10")
    # visual output
    plot_pdf_andor_png(plot_c, paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/dbscan_eps0.20_minPts10"), collapse=""), TRUE)
    # text output
    # store cluster assignments in data frame
    cluster_c <- as.data.frame(run_c$cluster+1L)
    colnames(cluster_c) <- "cluster"
    rownames(cluster_c) <- rownames(new_coords)
    # merge PCA coordinates and cluster assignments into one data frame
    genecoords_dbscan_c <- transform(merge(new_coords, cluster_c, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # merge original gene data and cluster assignments into one data frame
    raw_gene_table_cluster_dbscan_c <- transform(merge(raw_gene_table, genecoords_dbscan_c, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # output 1:
    # store g_names as an individual column (for output)
    raw_gene_table_cluster_with_g_names_dbscan_c <- cbind(g_name = rownames(raw_gene_table_cluster_dbscan_c), raw_gene_table_cluster_dbscan_c)
    # and write this info to a .csv file
    write.csv(raw_gene_table_cluster_with_g_names_dbscan_c, file=paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/dbscan_eps0.20_minPts10.MILTS.csv"), collapse=""), row.names=FALSE, quote=FALSE)

    # output 2:
    # create separate directory to store the files holding all genes by cluster assignment
    dir.create(paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/genes_by_cluster/dbscan_eps0.20_minPts10"), collapse=""), recursive=TRUE)
    # split the genes by cluster membership
    cluster_lists <- split(rownames(raw_gene_table_cluster_dbscan_c), raw_gene_table_cluster_dbscan_c$cluster) # --> this results in separate lists, stored in a superlist (cluster_lists)
    # for each list in the super list: print the gene names (row names) to a separate file
    lapply(names(cluster_lists), function(x){write.table(cluster_lists[[x]], file =paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/genes_by_cluster/dbscan_eps0.20_minPts10/cluster_", x, ".txt"), collapse=""),col.names=FALSE, row.names=FALSE, quote=FALSE)})



    run_d <- dbscan(new_coords, eps=0.25, minPts=10)
    cluster_d <- as.factor(run_d$cluster+1L)
    plot_d <- ggplot(new_coords, aes(x=Dim.1, y=Dim.2))+ geom_point(aes(shape=cluster_d, color=cluster_d)) + labs(color="cluster", shape="cluster") + theme_bw() + ggtitle("DB SCAN, eps=0.25, minPts=10")
    # visual output
    plot_pdf_andor_png(plot_d, paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/dbscan_eps0.25_minPts10"), collapse=""), TRUE)
    # text output
    # store cluster assignments in data frame
    cluster_d <- as.data.frame(run_d$cluster+1L)
    colnames(cluster_d) <- "cluster"
    rownames(cluster_d) <- rownames(new_coords)
    # merge PCA coordinates and cluster assignments into one data frame
    genecoords_dbscan_d <- transform(merge(new_coords, cluster_d, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # merge original gene data and cluster assignments into one data frame
    raw_gene_table_cluster_dbscan_d <- transform(merge(raw_gene_table, genecoords_dbscan_d, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # output 1:
    # store g_names as an individual column (for output)
    raw_gene_table_cluster_with_g_names_dbscan_d <- cbind(g_name = rownames(raw_gene_table_cluster_dbscan_d), raw_gene_table_cluster_dbscan_d)
    # and write this info to a .csv file
    write.csv(raw_gene_table_cluster_with_g_names_dbscan_d, file=paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/dbscan_eps0.25_minPts10.MILTS.csv"), collapse=""), row.names=FALSE, quote=FALSE)

    # output 2:
    # create separate directory to store the files holding all genes by cluster assignment
    dir.create(paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/genes_by_cluster/dbscan_eps0.25_minPts10"), collapse=""), recursive=TRUE)
    # split the genes by cluster membership
    cluster_lists <- split(rownames(raw_gene_table_cluster_dbscan_d), raw_gene_table_cluster_dbscan_d$cluster) # --> this results in separate lists, stored in a superlist (cluster_lists)
    # for each list in the super list: print the gene names (row names) to a separate file
    lapply(names(cluster_lists), function(x){write.table(cluster_lists[[x]], file =paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/genes_by_cluster/dbscan_eps0.25_minPts10/cluster_", x, ".txt"), collapse=""),col.names=FALSE, row.names=FALSE, quote=FALSE)})



    run_e <- dbscan(new_coords, eps=0.25, minPts=5)
    cluster_e <- as.factor(run_e$cluster+1L)
    plot_e <- ggplot(new_coords, aes(x=Dim.1, y=Dim.2))+ geom_point(aes(shape=cluster_e, color=cluster_e)) + labs(color="cluster", shape="cluster") + theme_bw() + ggtitle("DB SCAN, eps=0.25, minPts=5")
    # visual output
    plot_pdf_andor_png(plot_e, paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/dbscan_eps0.25_minPts5"), collapse=""), TRUE)
    # text output
    # store cluster assignments in data frame
    cluster_e <- as.data.frame(run_e$cluster+1L)
    colnames(cluster_e) <- "cluster"
    rownames(cluster_e) <- rownames(new_coords)
    # merge PCA coordinates and cluster assignments into one data frame
    genecoords_dbscan_e <- transform(merge(new_coords, cluster_e, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # merge original gene data and cluster assignments into one data frame
    raw_gene_table_cluster_dbscan_e <- transform(merge(raw_gene_table, genecoords_dbscan_e, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # output 1:
    # store g_names as an individual column (for output)
    raw_gene_table_cluster_with_g_names_dbscan_e <- cbind(g_name = rownames(raw_gene_table_cluster_dbscan_e), raw_gene_table_cluster_dbscan_e)
    # and write this info to a .csv file
    write.csv(raw_gene_table_cluster_with_g_names_dbscan_e, file=paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/dbscan_eps0.25_minPts5.MILTS.csv"), collapse=""), row.names=FALSE, quote=FALSE)

    # output 2:
    # create separate directory to store the files holding all genes by cluster assignment
    dir.create(paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/genes_by_cluster/dbscan_eps0.25_minPts5"), collapse=""), recursive=TRUE)
    # split the genes by cluster membership
    cluster_lists <- split(rownames(raw_gene_table_cluster_dbscan_e), raw_gene_table_cluster_dbscan_e$cluster) # --> this results in separate lists, stored in a superlist (cluster_lists)
    # for each list in the super list: print the gene names (row names) to a separate file
    lapply(names(cluster_lists), function(x){write.table(cluster_lists[[x]], file =paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/genes_by_cluster/dbscan_eps0.25_minPts5/cluster_", x, ".txt"), collapse=""),col.names=FALSE, row.names=FALSE, quote=FALSE)})



    run_f <- dbscan(new_coords, eps=0.30, minPts=5)
    cluster_f <- as.factor(run_f$cluster+1L)
    plot_f <- ggplot(new_coords, aes(x=Dim.1, y=Dim.2))+ geom_point(aes(shape=cluster_f, color=cluster_f)) + labs(color="cluster", shape="cluster") + theme_bw() + ggtitle("DB SCAN, eps=0.30, minPts=5")
    # visual output
    plot_pdf_andor_png(plot_f, paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/dbscan_eps0.30_minPts5"), collapse=""), TRUE)
    # text output
    # store cluster assignments in data frame
    cluster_f <- as.data.frame(run_f$cluster+1L)
    colnames(cluster_f) <- "cluster"
    rownames(cluster_f) <- rownames(new_coords)
    # merge PCA coordinates and cluster assignments into one data frame
    genecoords_dbscan_f <- transform(merge(new_coords, cluster_f, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # merge original gene data and cluster assignments into one data frame
    raw_gene_table_cluster_dbscan_f <- transform(merge(raw_gene_table, genecoords_dbscan_f, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # output 1:
    # store g_names as an individual column (for output)
    raw_gene_table_cluster_with_g_names_dbscan_f <- cbind(g_name = rownames(raw_gene_table_cluster_dbscan_f), raw_gene_table_cluster_dbscan_f)
    # and write this info to a .csv file
    write.csv(raw_gene_table_cluster_with_g_names_dbscan_f, file=paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/dbscan_eps0.30_minPts5.MILTS.csv"), collapse=""), row.names=FALSE, quote=FALSE)

    # output 2:
    # create separate directory to store the files holding all genes by cluster assignment
    dir.create(paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/genes_by_cluster/dbscan_eps0.30_minPts5"), collapse=""), recursive=TRUE)
    # split the genes by cluster membership
    cluster_lists <- split(rownames(raw_gene_table_cluster_dbscan_f), raw_gene_table_cluster_dbscan_f$cluster) # --> this results in separate lists, stored in a superlist (cluster_lists)
    # for each list in the super list: print the gene names (row names) to a separate file
    lapply(names(cluster_lists), function(x){write.table(cluster_lists[[x]], file =paste(c(cfg$output_path, "PCA_and_clustering/DBSCAN_clustering/genes_by_cluster/dbscan_eps0.30_minPts5/cluster_", x, ".txt"), collapse=""),col.names=FALSE, row.names=FALSE, quote=FALSE)})
  }
  else if ((is.numeric(cfg$custom_eps))&&(is.numeric(cfg$custom_minPts))) {
  custom_run <- dbscan(new_coords, eps=cfg$custom_eps, minPts=cfg$custom_minPts)
  custom_cluster <- as.factor(custom_run$cluster+1L)
  custom_plot <- ggplot(new_coords, aes(x=Dim.1, y=Dim.2))+ geom_point(aes(shape=custom_cluster, color=custom_cluster)) + labs(color="cluster", shape="cluster") + theme_bw() + ggtitle(paste(c("DB SCAN, eps=", cfg$custom_eps, ", minPts=", cfg$custom_minPts), collapse=""))
  # visual output
  plot_pdf_andor_png(custom_plot, paste(c(cfg$output_path,"PCA_and_clustering/DBSCAN_clustering/dbscan_eps", cfg$custom_eps, "_minPts", cfg$custom_minPts), collapse=""), TRUE)
  # text output
  # store cluster assignments in data frame
  custom_cluster <- as.data.frame(custom_run$cluster+1L)
  colnames(custom_cluster) <- "cluster"
  rownames(custom_cluster) <- rownames(new_coords)
  # merge PCA coordinates and cluster assignments into one data frame
  genecoords_dbscan_custom <- transform(merge(new_coords, custom_cluster, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
  # merge original gene data and cluster assignments into one data frame
  raw_gene_table_cluster_dbscan_custom <- transform(merge(raw_gene_table, genecoords_dbscan_custom, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
  # output 1:
  # store g_names as an individual column (for output)
  raw_gene_table_cluster_with_g_names_dbscan_custom <- cbind(g_name = rownames(raw_gene_table_cluster_dbscan_custom), raw_gene_table_cluster_dbscan_custom)
  # and write this info to a .csv file
  write.csv(raw_gene_table_cluster_with_g_names_dbscan_custom, file=paste(c(cfg$output_path,"PCA_and_clustering/DBSCAN_clustering/dbscan_eps", cfg$custom_eps, "_minPts", cfg$custom_minPts, ".MILTS.csv"), collapse=""), row.names=FALSE, quote=FALSE)

  # output 2:
  # create separate directory to store the files holding all genes by cluster assignment
  dir.create(paste(c(cfg$output_path,"PCA_and_clustering/DBSCAN_clustering/dbscan_eps", cfg$custom_eps, "_minPts", cfg$custom_minPts), collapse=""), recursive=TRUE)
  # split the genes by cluster membership
  cluster_lists <- split(rownames(raw_gene_table_cluster_dbscan_custom), raw_gene_table_cluster_dbscan_custom$cluster) # --> this results in separate lists, stored in a superlist (cluster_lists)
  # for each list in the super list: print the gene names (row names) to a separate file
  lapply(names(cluster_lists), function(x){write.table(cluster_lists[[x]], file =paste(c(cfg$output_path,"PCA_and_clustering/DBSCAN_clustering/dbscan_eps", cfg$custom_eps, "_minPts", cfg$custom_minPts, "/cluster_", x, ".txt"), collapse=""),col.names=FALSE, row.names=FALSE, quote=FALSE)})
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
  dir.create(paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering"), collapse=""))
  library(mclust)

  if (cfg$mclust_k=="default") {
    # ---- 2 CLUSTERS ----
    model2 <- Mclust(new_coords, G=2)
    cluster_2 <- as.factor(model2$classification)
    plot_2 <- ggplot(new_coords, aes(x=Dim.1, y=Dim.2))+ scale_color_viridis(discrete=TRUE, begin=0.1, end=0.9) + geom_point(aes(shape=cluster_2, color=cluster_2)) + labs(color="cluster", shape="cluster") + theme_bw() + ggtitle("model-based clustering, k=2")

    # visual output
    plot_pdf_andor_png(plot_2, paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/mclust_2groups"), collapse=""), TRUE)
    # text output
    # store cluster assignments in data frame
    cluster_2 <- as.data.frame(model2$classification)
    colnames(cluster_2) <- "cluster"
    # merge PCA coordinates and cluster assignments into one data frame
    genecoords_mclust2 <- transform(merge(new_coords, cluster_2, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # merge original gene data and cluster assignments into one data frame
    raw_gene_table_cluster_mclust2 <- transform(merge(raw_gene_table, genecoords_mclust2, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # output 1:
    # store g_names as an individual column (for output)
    raw_gene_table_cluster_with_g_names_mclust2 <- cbind(g_name = rownames(raw_gene_table_cluster_mclust2), raw_gene_table_cluster_mclust2)
    # and write this info to a .csv file
    write.csv(raw_gene_table_cluster_with_g_names_mclust2, file=paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/mclust_2groups.MILTS.csv"), collapse=""), row.names=FALSE, quote=FALSE)

    # output 2:
    # create separate directory to store the files holding all genes by cluster assignment
    dir.create(paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/genes_by_cluster/mclust_2"), collapse=""), recursive=TRUE)
    # split the genes by cluster membership
    cluster_lists <- split(rownames(raw_gene_table_cluster_mclust2), raw_gene_table_cluster_mclust2$cluster) # --> this results in separate lists, stored in a superlist (cluster_lists)
    # for each list in the super list: print the gene names (row names) to a separate file
    lapply(names(cluster_lists), function(x){write.table(cluster_lists[[x]], file =paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/genes_by_cluster/mclust_2/cluster_", x, ".txt"), collapse=""),col.names=FALSE, row.names=FALSE, quote=FALSE)})


    # ---- 3 CLUSTERS ----
    model3 <- Mclust(new_coords, G=3)
    cluster_3 <- as.factor(model3$classification)
    plot_3 <- ggplot(new_coords, aes(x=Dim.1, y=Dim.2))+ scale_color_viridis(discrete=TRUE, begin=0.1, end=0.9) + geom_point(aes(shape=cluster_3, color=cluster_3)) + labs(color="cluster", shape="cluster") + theme_bw() + ggtitle("model-based clustering, k=3")

    # visual output
    plot_pdf_andor_png(plot_3, paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/mclust_3groups"), collapse=""), TRUE)
    # text output
    # store cluster assignments in data frame
    cluster_3 <- as.data.frame(model3$classification)
    colnames(cluster_3) <- "cluster"
    # merge PCA coordinates and cluster assignments into one data frame
    genecoords_mclust3 <- transform(merge(new_coords, cluster_3, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # merge original gene data and cluster assignments into one data frame
    raw_gene_table_cluster_mclust3 <- transform(merge(raw_gene_table, genecoords_mclust3, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # output 1:
    # store g_names as an individual column (for output)
    raw_gene_table_cluster_with_g_names_mclust3 <- cbind(g_name = rownames(raw_gene_table_cluster_mclust3), raw_gene_table_cluster_mclust3)
    # and write this info to a .csv file
    write.csv(raw_gene_table_cluster_with_g_names_mclust3, file=paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/mclust_3groups.MILTS.csv"), collapse=""), row.names=FALSE, quote=FALSE)

    # output 2:
    # create separate directory to store the files holding all genes by cluster assignment
    dir.create(paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/genes_by_cluster/mclust_3"), collapse=""), recursive=TRUE)
    # split the genes by cluster membership
    cluster_lists <- split(rownames(raw_gene_table_cluster_mclust3), raw_gene_table_cluster_mclust3$cluster) # --> this results in separate lists, stored in a superlist (cluster_lists)
    # for each list in the super list: print the gene names (row names) to a separate file
    lapply(names(cluster_lists), function(x){write.table(cluster_lists[[x]], file =paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/genes_by_cluster/mclust_3/cluster_", x, ".txt"), collapse=""),col.names=FALSE, row.names=FALSE, quote=FALSE)})

    # ---- 4 CLUSTERS ----
    model4 <- Mclust(new_coords, G=4)
    cluster_4 <- as.factor(model4$classification)
    plot_4 <- ggplot(new_coords, aes(x=Dim.1, y=Dim.2))+ scale_color_viridis(discrete=TRUE, begin=0.1, end=0.9) + geom_point(aes(shape=cluster_4, color=cluster_4)) + labs(color="cluster", shape="cluster") + theme_bw() + ggtitle("model-based clustering, k=4")

    # visual output
    plot_pdf_andor_png(plot_4, paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/mclust_4groups"), collapse=""), TRUE)
    # text output
    # store cluster assignments in data frame
    cluster_4 <- as.data.frame(model4$classification)
    colnames(cluster_4) <- "cluster"
    # merge PCA coordinates and cluster assignments into one data frame
    genecoords_mclust4 <- transform(merge(new_coords, cluster_4, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # merge original gene data and cluster assignments into one data frame
    raw_gene_table_cluster_mclust4 <- transform(merge(raw_gene_table, genecoords_mclust4, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # output 1:
    # store g_names as an individual column (for output)
    raw_gene_table_cluster_with_g_names_mclust4 <- cbind(g_name = rownames(raw_gene_table_cluster_mclust4), raw_gene_table_cluster_mclust4)
    # and write this info to a .csv file
    write.csv(raw_gene_table_cluster_with_g_names_mclust4, file=paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/mclust_4groups.MILTS.csv"), collapse=""), row.names=FALSE, quote=FALSE)

    # output 2:
    # create separate directory to store the files holding all genes by cluster assignment
    dir.create(paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/genes_by_cluster/mclust_4"), collapse=""), recursive=TRUE)
    # split the genes by cluster membership
    cluster_lists <- split(rownames(raw_gene_table_cluster_mclust4), raw_gene_table_cluster_mclust4$cluster) # --> this results in separate lists, stored in a superlist (cluster_lists)
    # for each list in the super list: print the gene names (row names) to a separate file
    lapply(names(cluster_lists), function(x){write.table(cluster_lists[[x]], file =paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/genes_by_cluster/mclust_4/cluster_", x, ".txt"), collapse=""),col.names=FALSE, row.names=FALSE, quote=FALSE)})

    # ---- OPTIMAL (number clusters defined by BIC) ----
    model_opt <- Mclust(new_coords)
    cluster_opt <- as.factor(model_opt$classification)
    plot_opt <- ggplot(new_coords, aes(x=Dim.1, y=Dim.2))+ scale_color_viridis(discrete=TRUE, begin=0.1, end=0.9) + geom_point(aes(shape=cluster_opt, color=cluster_opt)) + labs(color="cluster", shape="cluster") + theme_bw() + ggtitle("model-based clustering, optimal #clusters")

    # visual output
    plot_pdf_andor_png(plot_opt, paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/mclust_opt_groups"), collapse=""), TRUE)

    # text output
    # store cluster assignments in data frame
    cluster_opt <- as.data.frame(model_opt$classification)
    # add rownames (gene names) to the cluster assignments
    rownames(cluster_opt) <- rownames(new_coords)
    # add column title
    colnames(cluster_opt) <- "cluster"
    # merge PCA coordinates and cluster assignments into one data frame
    genecoords_mclustopt <- transform(merge(new_coords, cluster_opt, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # merge original gene data and cluster assignments into one data frame
    raw_gene_table_cluster_mclustopt <- transform(merge(raw_gene_table, genecoords_mclustopt, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # output 1:
    # store g_names as an individual column (for output)
    raw_gene_table_cluster_with_g_names_mclustopt <- cbind(g_name = rownames(raw_gene_table_cluster_mclustopt), raw_gene_table_cluster_mclustopt)
    # and write this info to a .csv file
    write.csv(raw_gene_table_cluster_with_g_names_mclustopt, file=paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/mclust_opt_groups.MILTS.csv"), collapse=""), row.names=FALSE, quote=FALSE)

    # output 2:
    # create separate directory to store the files holding all genes by cluster assignment
    dir.create(paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/genes_by_cluster/mclust_opt"), collapse=""), recursive=TRUE)
    # split the genes by cluster membership
    cluster_lists <- split(rownames(raw_gene_table_cluster_mclustopt), raw_gene_table_cluster_mclustopt$cluster) # --> this results in separate lists, stored in a superlist (cluster_lists)
    # for each list in the super list: print the gene names (row names) to a separate file
    lapply(names(cluster_lists), function(x){write.table(cluster_lists[[x]], file =paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/genes_by_cluster/mclust_opt/cluster_", x, ".txt"), collapse=""),col.names=FALSE, row.names=FALSE, quote=FALSE)})
  } else if (cfg$mclust_k=="BIC") {
    # ---- OPTIMAL (number clusters defined by BIC) ----
    model_opt <- Mclust(new_coords)
    cluster_opt <- as.factor(model_opt$classification)
    plot_opt <- ggplot(new_coords, aes(x=Dim.1, y=Dim.2))+ scale_color_viridis(discrete=TRUE, begin=0.1, end=0.9) + geom_point(aes(shape=cluster_opt, color=cluster_opt)) + labs(color="cluster", shape="cluster") + theme_bw() + ggtitle("model-based clustering, optimal #clusters")

    # visual output
    plot_pdf_andor_png(plot_opt, paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/mclust_opt_groups"), collapse=""), TRUE)

    # text output
    # store cluster assignments in data frame
    cluster_opt <- as.data.frame(model_opt$classification)
    # add rownames (gene names) to the cluster assignments
    rownames(cluster_opt) <- rownames(new_coords)
    # add column title
    colnames(cluster_opt) <- "cluster"
    # merge PCA coordinates and cluster assignments into one data frame
    genecoords_mclustopt <- transform(merge(new_coords, cluster_opt, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # merge original gene data and cluster assignments into one data frame
    raw_gene_table_cluster_mclustopt <- transform(merge(raw_gene_table, genecoords_mclustopt, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # output 1:
    # store g_names as an individual column (for output)
    raw_gene_table_cluster_with_g_names_mclustopt <- cbind(g_name = rownames(raw_gene_table_cluster_mclustopt), raw_gene_table_cluster_mclustopt)
    # and write this info to a .csv file
    write.csv(raw_gene_table_cluster_with_g_names_mclustopt, file=paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/mclust_opt_groups.MILTS.csv"), collapse=""), row.names=FALSE, quote=FALSE)

    # output 2:
    # create separate directory to store the files holding all genes by cluster assignment
    dir.create(paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/genes_by_cluster/mclust_opt"), collapse=""), recursive=TRUE)
    # split the genes by cluster membership
    cluster_lists <- split(rownames(raw_gene_table_cluster_mclustopt), raw_gene_table_cluster_mclustopt$cluster) # --> this results in separate lists, stored in a superlist (cluster_lists)
    # for each list in the super list: print the gene names (row names) to a separate file
    lapply(names(cluster_lists), function(x){write.table(cluster_lists[[x]], file =paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/genes_by_cluster/mclust_opt/cluster_", x, ".txt"), collapse=""),col.names=FALSE, row.names=FALSE, quote=FALSE)})
  }

  # if the user did not choose the "default" output
  # but specified a desired number of clusters (in mclust_k)
  else {
    mclust_k <- as.numeric(cfg$mclust_k)

    modelk <- Mclust(new_coords, G=mclust_k)
    cluster_k <- as.factor(modelk$classification)
    plot_k <- ggplot(new_coords, aes(x=Dim.1, y=Dim.2))+ scale_color_viridis(discrete=TRUE, begin=0.1, end=0.9) + geom_point(aes(shape=cluster_k, color=cluster_k)) + labs(color="cluster", shape="cluster") + theme_bw() + ggtitle(paste(c("model-based clustering, k=",mclust_k), collapse=""))

    # visual output
    plot_pdf_andor_png(plot_k, paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/mclust_",mclust_k, "groups"), collapse=""), TRUE)
    # text output
    # store cluster assignments in data frame
    cluster_k <- as.data.frame(modelk$classification)
    colnames(cluster_k) <- "cluster"
    # merge PCA coordinates and cluster assignments into one data frame
    genecoords_mclustk <- transform(merge(new_coords, cluster_k, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # merge original gene data and cluster assignments into one data frame
    raw_gene_table_cluster_mclustk <- transform(merge(raw_gene_table, genecoords_mclustk, by = "row.names", all = TRUE), row.names=Row.names, Row.names=NULL)
    # output 1:
    # store g_names as an individual column (for output)
    raw_gene_table_cluster_with_g_names_mclustk <- cbind(g_name = rownames(raw_gene_table_cluster_mclustk), raw_gene_table_cluster_mclustk)
    # and write this info to a .csv file
    write.csv(raw_gene_table_cluster_with_g_names_mclustk, file=paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/mclust_", mclust_k, "groups.MILTS.csv"), collapse=""), row.names=FALSE, quote=FALSE)

    # output 2:
    # create separate directory to store the files holding all genes by cluster assignment
    dir.create(paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/genes_by_cluster/mclust_", mclust_k), collapse=""), recursive=TRUE)
    # split the genes by cluster membership
    cluster_lists <- split(rownames(raw_gene_table_cluster_mclustk), raw_gene_table_cluster_mclustk$cluster) # --> this results in separate lists, stored in a superlist (cluster_lists)
    # for each list in the super list: print the gene names (row names) to a separate file
    lapply(names(cluster_lists), function(x){write.table(cluster_lists[[x]], file =paste(c(cfg$output_path, "PCA_and_clustering/model-based_clustering/genes_by_cluster/mclust_", mclust_k, "/cluster_", x, ".txt"), collapse=""),col.names=FALSE, row.names=FALSE, quote=FALSE)})

  }
  print(proc.time() - ptm)
}
# end of model based clustering
