# SCRIPT TO ADD COMMA SEPARATED taXaminer LINE AS AN ATTRIBUTE TO THE 9TH COLUMN OF A GIVEN GFF FILE
# author: Simonida Zehr
# date: 14 June 2019

# requires Python 3.6.1 or higher
# call:
# python3 add_taXaminer_to_gff.py <path.to.gff> <path_to.my_clustering.taXaminer.csv>
# e.g.
# python3 add_taXaminer_to_gff.py my_annotations.gff hclust_2groups.taXaminer.csv 
# the file will create a my_annotations.with_taXaminer.<my_clustering>.gff file in the directory the GFF file is located in

import sys
gff_path=sys.argv[1] # 1st parameter passed to this script is the GFF file path
taxaminer_path=sys.argv[2] # 2nd parameter passed is the taXaminer path (csv)

# hash in which each gene (key) will be mapped to its taXaminer (value)
gene_taxaminer_mapping = {}
header=[]


# read taXaminer file to retrieve all genes taXaminer information is available on
with open(taxaminer_path, "r") as taxaminer:
    for line in taxaminer:
        # remove newline from line
        line = line.strip()

        # get the gene name (first entry in comma seaprated line)
        taxaminer_info_table = line.split(",")
        # remove g_name info from table and store it to variable
        g_name = taxaminer_info_table.pop(0)
        # remove c_name info from table
        taxaminer_info_table.pop(0)

        # for the first line, the g_name value will be "g_name" -> this is the header
        if g_name=="g_name":
            header = taxaminer_info_table # store header items to header table

        # for every gene info line:
        else:
           # generate a taXaminer info string in "tag=value" (separated by commas) format
            taxaminer_string = ";taXaminer="
            for (tag, value) in zip(header, taxaminer_info_table):
                taxaminer_string = taxaminer_string + tag + "=" + value + ","
            # remove last comma
            taxaminer_string = taxaminer_string[:-1]
            # assign the taXaminer line to the gene
            gene_taxaminer_mapping[g_name] = taxaminer_string


with open(gff_path, "r") as gff: # open GFF file
    # define output path name
    gff_with_taxaminer_path = gff_path[:-4] + ".with_taXaminer." + taxaminer_path[:-10] + ".gff"

    # open output file to write to
    with open(gff_with_taxaminer_path, "w") as gff_with_taxaminer:

        # read every line of the original GFF file
        for line in gff:
            # remove newline from line
            line = line.strip()

            # if anntations have been reached:
            if not line.startswith("#"):
                # split GFF line to access the individual column entries
                gff_array = line.split()

                # only for gene feature annotations:
                if (gff_array[2]=="gene" or gff_array[2]=="pseudogene"):
                    # get gene name from attributes entry in 9th column
                    attributes = gff_array[8].split(';') # get attributes entry
                    id_attribute = attributes[0].split('=') # something like ID=gene1
                    gene_name = id_attribute[1] # get only the name after the equals sign

                    # if there is taXaminer info on this gene
                    if gene_name in gene_taxaminer_mapping:
                        # add a taXaminer attribute to the list of attributes in the 9th column of the GFF file
#                        line = line + ";taXaminer=" + gene_taxaminer_mapping[gene_name]
                        line = line + gene_taxaminer_mapping[gene_name]
            # whether the line has been modified or not, write it to the output GFF
            gff_with_taxaminer.write(line + "\n")
