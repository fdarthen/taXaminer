#!/usr/bin/env Rscript

library(ggplot2) # for plotting
library(viridis) # for colours
library(data.table)
library(plotly)
library(htmlwidgets)
library(yaml) # for reading config file

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[1]

cfg <- yaml.load_file(config_path)

plot_pdf_andor_png <- function(plot, plot_path, spec) {
    # spec: if dimensions and resolution should be specified
    if (spec) {
      if (cfg$output_png) {
        png(paste0(plot_path,".png"), units="in", width=6, height=6, res=500)
        print(plot)
        dev.off()
      }
      if (cfg$output_pdf) {
        pdf(paste0(plot_path,".pdf"), width=6, height=6)
        print(plot)
        dev.off()
      }
    }

    # if not (bc e.g. fviz plots may be impaired by this)
    else {
      if (cfg$output_png) {
        png(paste0(plot_path,".png"))
        print(plot)
        dev.off()
      }
      if (cfg$output_pdf) {
        pdf(paste0(plot_path,".pdf"))
        print(plot)
        dev.off()
      }
    }
  }

orca_static_3d <- function(plot, plot_filename, eye_x, eye_y, eye_z) {
    m <- list(l = 0, r = 0, b = 0)
    s <- list(camera = list(eye = list(x = eye_x, y = eye_y, z = eye_z)),
        xaxis=list(title="PC 1"),yaxis=list(title="PC 2"),zaxis=list(title="PC 3"),
        aspectmode='cube')
    plot <- plot %>% layout(title = 'taxonomic assignment', scene = s, margin = m)

    # save plots to json files to generate static plots exterenally with orca
    # orca function in R does not work and using system() does not work on cluster
    plotly_json(plot, FALSE) %>% cat(file = paste0(cfg$output_path,"tmp/", plot_filename, ".json"))

}

# CREATE MANUAL COLOR SCHEME FOR A LIST

#' Create qualitative colours
#' @param n number of colors
#' @param light light colors TRUE or FALSE
#' @return list of n different colors
#' @source Modified based on https://gist.github.com/peterk87/6011397
#' @examples
#' \dontrun{
#' qualitativeColours(5)
#' }

qualitativeColours <- function(n, light = FALSE) {
    # For more than 21 colours needed
    rich12equal <- c(
        "#000040", "#000093", "#0020E9", "#0076FF", "#00B8C2", "#04E466",
        "#49FB25", "#E7FD09", "#FEEA02", "#FFC200", "#FF8500", "#FF3300")
    # Qualitative colour schemes by Paul Tol
    ifelse (n >= 19 & n <= 21, return(grDevices::colorRampPalette(c(
        "#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD",
        "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA",
        "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77",
        "#771122", "#AA4455", "#DD7788"))(n)),
        ifelse (n >= 16 & n <= 18, return(grDevices::colorRampPalette(c(
            "#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD",
            "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77",
            "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455",
            "#DD7788"))(n)),
            ifelse (n == 15, return(grDevices::colorRampPalette(c(
                "#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA",
                "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711",
                "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77",
                "#771122", "#AA4455", "#DD7788"))(n)),
                ifelse(n > 12 & n <= 14, return(grDevices::colorRampPalette(c(
                    "#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7",
                    "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55",
                    "#F6C141", "#F1932D", "#E8601C", "#DC050C"))(n)),
                    ifelse(n > 9 & n <= 12,
                        ifelse(
                            light,
                            return(RColorBrewer::brewer.pal(
                                n = n, name = 'Set3')),
                            return(RColorBrewer::brewer.pal(
                                n = n, name = 'Paired'))),
                        ifelse(n <= 9,
                            ifelse(
                                light,
                                return(RColorBrewer::brewer.pal(
                                    n = n, name = 'Pastel1')),
                                return(RColorBrewer::brewer.pal(
                                    n = n, name = 'Set1'))),
                            return(grDevices::colorRampPalette(rich12equal)(n))
                        )
                    )
                )
            )
        )
    )
}

#' Get color for a list of items
#' @return list of colors for each element (same elements will have the same
#' color)
#' @param x input list
#' @export
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{qualitativeColours}}
#' @examples
#' items <- c("a", "b", "c")
#' getQualColForVector(items)

getQualColForVector <- function(x = NULL) {
    if (is.null(x)) stop("Input list is NULL!")
    types <- unique(x)
    types <- types[!is.na(types)]
    typeColors <- qualitativeColours(length(types))
    colorsTypes <- as.vector(x)
    countTypes <- countColors <- 1
    while (countTypes <= length(types)) {
        if (countColors > length(typeColors)) countColors <- 1
        colorsTypes[colorsTypes == types[countTypes]] <- typeColors[countColors]
        countColors <- countColors + 1
        countTypes <- countTypes + 1
    }
    return(unlist(colorsTypes))
}

#  _________________ READING DATA __________________ #

genes_coords_taxon <- data.table::fread(paste0(cfg$output_path,'taxonomic_assignment/gene_table_taxon_assignment.csv'), sep=",", header=TRUE)
query_label_name <- genes_coords_taxon$plot_label[genes_coords_taxon$plot_label != genes_coords_taxon$taxon_assignment & genes_coords_taxon$plot_label != "Unassigned" & genes_coords_taxon$plot_label != "Otherwise assigned"][1]
print(query_label_name)
if (is.null(query_label_name) | is.na(query_label_name)) {
    query_label_name <- "Undefined"
}

# ________________ PLOT PREPARATION _______________ #

# frequency information for labels in plots
label_count_table <- data.frame(table(genes_coords_taxon$plot_label)) # how often each label appears
label_count_table$plot_label_freq <- as.character(paste0(label_count_table$Var1," (",label_count_table$Freq,")"))
genes_coords_taxon <- merge(genes_coords_taxon, label_count_table[,c("Var1","plot_label_freq")], by.x = "plot_label", by.y = "Var1", all.x = TRUE)
label <- sort(unique(label_count_table$plot_label_freq))


# add sequence information
prot_fasta <- Biostrings::readAAStringSet(cfg$proteins_path)
protID <- names(prot_fasta)
query_seq <- paste(prot_fasta)
prot_sequences <- data.frame(protID, query_seq)
prot_sequences$protID <- gsub('\\ .*', '' , prot_sequences$protID)
prot_sequences$query_seq <- gsub("(.{70}?)", "\\1</br>", prot_sequences$query_seq)
genes_coords_taxon <- merge(genes_coords_taxon, prot_sequences, by="protID", all.x=TRUE)


# prepare data for hover window display
genes_coords_taxon$g_terminal[genes_coords_taxon$g_terminal == "0"] <- "no"
genes_coords_taxon$g_terminal[genes_coords_taxon$g_terminal == "1"] <- "yes"
g_cov_vars <- as.factor(grep("g_cov_[0-9]",colnames(genes_coords_taxon), value=TRUE))
genes_coords_taxon$g_coverages <- apply(subset(genes_coords_taxon, select=g_cov_vars),1,paste,collapse="; ")
g_covdev_vars <- as.factor(grep("g_covdev_c_[0-9]",colnames(genes_coords_taxon), value=TRUE))
genes_coords_taxon$g_covdeviations <- apply(subset(genes_coords_taxon, select=g_covdev_vars),1,paste,collapse="; ")

genes_coords_taxon$label_color[genes_coords_taxon$plot_label == query_label_name] <- "#404a4a"
genes_coords_taxon$label_color[genes_coords_taxon$plot_label == "Unassigned"] <- "#778899"
genes_coords_taxon$label_color[genes_coords_taxon$plot_label != query_label_name & genes_coords_taxon$plot_label != "Unassigned"] <- getQualColForVector(genes_coords_taxon$plot_label_freq[genes_coords_taxon$plot_label != query_label_name & genes_coords_taxon$plot_label != "Unassigned"])


genes_coords_taxon_query <- genes_coords_taxon[genes_coords_taxon$plot_label == query_label_name,]
genes_coords_taxon_unassigned <- genes_coords_taxon[genes_coords_taxon$plot_label == "Unassigned",]
genes_coords_taxon_rest <- genes_coords_taxon[genes_coords_taxon$plot_label != query_label_name & genes_coords_taxon$plot_label != "Unassigned",]

# 1D density of 1 element groups is not possible
# for 1D density only use taxons with abundancy of at least 2
rownames(label_count_table) <- label_count_table$Var1
genes_taxon_1dim <- genes_coords_taxon[genes_coords_taxon$plot_label %in% rownames(label_count_table[label_count_table$Freq > 1,]),]


# plot_label column as factor for marking in plots
taxon_factor <- as.factor(genes_coords_taxon$plot_label)
taxon_rest_factor <- as.factor(genes_coords_taxon_rest$plot_label)
taxon_query_factor <- as.factor(genes_coords_taxon_query$plot_label)
taxon_unassigned_factor <- as.factor(genes_coords_taxon_unassigned$plot_label)
taxons_1d <- as.factor(genes_taxon_1dim$plot_label)


# _______________PLOT GENERATION_____________________ #

# manual colorscale for plots with all taxons
myColors <- unique(genes_coords_taxon$label_color)
names(myColors) <- unique(genes_coords_taxon$plot_label_freq)
colScale <- scale_colour_manual(name = "species", labels=label, values = myColors)

if (nrow(genes_taxon_1dim) > 0) {
    # label for density plots
    label_1d <- c(sort(unique(genes_taxon_1dim$plot_label_freq)))

    # manual colorscale for 1D plots
    myColors <- unique(genes_taxon_1dim$label_color)
    names(myColors) <- unique(genes_taxon_1dim$plot_label_freq)
    colScale_1D <- scale_colour_manual(name = "species", labels=label_1d, values = myColors)

    plot_x <- ggplot(genes_taxon_1dim,aes(colour = plot_label_freq)) +
        colScale_1D +
        geom_density(data=genes_taxon_1dim, aes(x=Dim.1)) +
        theme_bw() +
        ggtitle("density of Dim.1")

    plot_pdf_andor_png(plot_x, paste(c(cfg$output_path, "taxonomic_assignment/density_x"), collapse=""), TRUE)

    plot_y <- ggplot(genes_taxon_1dim,aes(colour = plot_label_freq)) +
        colScale_1D +
        geom_density(data=genes_taxon_1dim, aes(x=Dim.2)) +
        theme_bw() +
        ggtitle("density of Dim.2")

    plot_pdf_andor_png(plot_y, paste(c(cfg$output_path, "taxonomic_assignment/density_y"), collapse=""), TRUE)
}


plot_2 <- ggplot(genes_coords_taxon,aes(colour = plot_label_freq)) +
    scale_shape_manual(name="species", labels=label, values=c(1,2,15,16,17,3,7,8,9,23,18,4,10,11,12,13,14,0)) +
    colScale
if (nrow(genes_coords_taxon_query) > 0) {
    plot_2 <- plot_2 +
        geom_point(data=genes_coords_taxon_query, aes(x=Dim.1, y=Dim.2)) + #, shape=taxon_query_factor
        stat_density_2d(data=genes_coords_taxon_query,aes(x=Dim.1, y=Dim.2), geom="polygon", alpha=0.05)
}
if (nrow(genes_coords_taxon_unassigned) > 0) {
    plot_2 <- plot_2 +
        geom_point(data=genes_coords_taxon_unassigned, aes(x=Dim.1, y=Dim.2)) #, shape=taxon_unassigned_factor
}
if (nrow(genes_coords_taxon_rest) > 0) {
    plot_2 <- plot_2 +
        geom_point(data=genes_coords_taxon_rest, aes(x=Dim.1, y=Dim.2)) #, shape=taxon_rest_factor
}
plot_2 <- plot_2 +
      theme_bw() +
      ggtitle("2D density")
plot_pdf_andor_png(plot_2, paste(c(cfg$output_path, "taxonomic_assignment/density_2d"), collapse=""), TRUE)

if (length(grep("Dim.",colnames(genes_coords_taxon),value=TRUE)) >= 3){

    if (nrow(genes_taxon_1dim) > 0) {
        plot_z <- ggplot(genes_taxon_1dim,aes(colour = plot_label_freq)) +
            colScale_1D +
            geom_density(data=genes_taxon_1dim, aes(x=Dim.3)) +
            theme_bw() +
            ggtitle("density of Dim.3")

        plot_pdf_andor_png(plot_z, paste(c(cfg$output_path, "taxonomic_assignment/density_z"), collapse=""), TRUE)
    }

    fig <- plot_ly(type="scatter3d", mode="markers", symbols=c('circle'))

    fig <- fig %>% add_markers(fig, data=genes_coords_taxon_query, x = ~Dim.1, y = ~Dim.2, z = ~Dim.3,
        hoverinfo="text",
        text = ~paste('</br>ID:',g_name, '| Taxonomic assignment:',taxon_assignment,'| Label:',plot_label,
            '</br>Coverage:',g_coverages, '(SD from contig mean:',g_covdeviations,')',
            '</br>Terminal:',g_terminal,'(Genes on contig:',c_num_of_genes,')',
            '</br>LCA:',lca_hit,
            '</br>Best hit:',best_hit,'(e-value:',evalue,')',
            '</br>Seq:',query_seq),
        textposition="bottom left",
        opacity=0.75,
        size=~I(25),
        color=~I(label_color),
        name=~plot_label_freq
    )

    fig <- fig %>% add_markers(fig, data=genes_coords_taxon_unassigned, x = ~Dim.1, y = ~Dim.2, z = ~Dim.3,
        hoverinfo="text",
        text = ~paste('</br>ID:',g_name, '| Taxonomic assignment:',taxon_assignment,'| Label:',plot_label,
            '</br>Coverage:',g_coverages, '(SD from contig mean:',g_covdeviations,')',
            '</br>Terminal:',g_terminal,'(Genes on contig:',c_num_of_genes,')',
            '</br>LCA:',lca_hit,
            '</br>Best hit:',best_hit,'(e-value:',evalue,')',
            '</br>Seq:',query_seq),
        textposition="bottom left",
        opacity=0.5,
        size=~I(25),
        color=~I(label_color),
        name=~plot_label_freq
    )

    fig <- fig %>% add_markers(fig, data=genes_coords_taxon_rest, x = ~Dim.1, y = ~Dim.2, z = ~Dim.3,
        hoverinfo="text",
        text = ~paste('</br>ID:',g_name, '| Taxonomic assignment:',taxon_assignment,'| Label:',plot_label,
            '</br>Coverage:',g_coverages, '(SD from contig mean:',g_covdeviations,')',
            '</br>Terminal:',g_terminal,'(Genes on contig:',c_num_of_genes,')',
            '</br>LCA:',lca_hit,
            '</br>Best hit:',best_hit,'(e-value:',evalue,')',
            '</br>Seq:',query_seq),
        textposition="bottom left",
        size=~I(25),
        color=~I(label_color),
        name=~plot_label_freq
    )

    orca_static_3d(fig, "2D_plot_1_2", 0, 0, 2.25)
    orca_static_3d(fig, "2D_plot_1_3", 0, 2.25, 0)
    orca_static_3d(fig, "2D_plot_2_3", 2.25, 0, 0)
    orca_static_3d(fig, "3D_plot", 1.55, 1.55, 1.55)

    fig <- fig %>% layout(title = 'taxonomic assignment', scene = list(xaxis=list(title="PC 1"),yaxis=list(title="PC 2"),zaxis=list(title="PC 3"),aspectmode='cube'))
    withr::with_dir(paste0(cfg$output_path, "taxonomic_assignment/"), saveWidget(as_widget(fig), file="3D_plot.html", selfcontained=FALSE))
}
