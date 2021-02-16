# reads GFF file and adds an annotation line for contigs

# information based on directives/pragma (##sequence-region scaffold...) or FASTA of contigs
# type is stated as "region"

import sys

gff_path=sys.argv[1] # input GFF
out_path=sys.argv[2] # output GFF
if len(sys.argv) > 3:
    fasta_path=sys.argv[3] # input GFF




def directives(gff_path,out_path):
    """ uses lines with pragmas to retrieve information for contig entry line """
    with open(gff_path, 'r') as gff: # open GFF file
        with open(out_path, 'w') as out: # open output file
            for line in gff:
                if line.startswith("##sequence-region "): # this line is a directive specifying a contig
                    out.write(line)
                    spline = line.strip().split(" ")
                    contig_name = spline[1]
                    contig_start = spline[2]
                    contig_end = spline[3]
                    out.write(contig_name+"\t"+"."+"\t"+"region"+"\t"+contig_start+"\t"+contig_end+"\t"+"."+"\t"+"+"+"\t"+"."+"\t"+"."+"\n")
                else:
                    out.write(line)


def no_directives_fasta(gff_path,out_path,fasta_path):
    """ gathers information about contigs from the FASTA
    uses header (everything after '>') and len of sequence for contig stats """
    with open(fasta_path, 'r') as fasta:
        scaffolds = {}
        size = 0
        scaffold = None
        for line in fasta:
            if line.startswith('>'):
                scaffolds[scaffold] = size
                scaffold = line.lstrip('>').strip()
                size = 0
            else:
                size += len(line.strip())
        scaffolds[scaffold] = size


    with open(gff_path, 'r') as gff: # open GFF file
        with open(out_path, 'w') as out: # open output file
            scaffold = None
            for line in gff:
                if line.startswith("scaffold"):
                    if line.split("\t")[0] == scaffold:
                        out.write(line)
                    else:
                        scaffold = line.split("\t")[0]
                        if int(scaffold.split("size")[1]) != int(scaffolds.get(scaffold)):
                            print(scaffold, scaffold.split("size")[1], scaffolds.get(scaffold))
                        out.write(scaffold+"\t"+"."+"\t"+"region"+"\t"+"1"+"\t"+str(scaffolds.get(scaffold))+"\t"+"."+"\t"+"+"+"\t"+"."+"\t"+"."+"\n")
                        out.write(line)
                        del scaffolds[scaffold]
                else:
                    out.write(line)
            if len(scaffolds) != 0:
                print(scaffolds.keys()) # featureless contigs
                for scaffold in scaffolds.keys():
                    if not scaffold is None:
                        out.write(scaffold+"\t"+"."+"\t"+"region"+"\t"+"1"+"\t"+str(scaffolds.get(scaffold))+"\t"+"."+"\t"+"+"+"\t"+"."+"\t"+"."+"\n")

if fasta_path:
    no_directives_fasta(gff_path,out_path,fasta_path)
else:
    directives(gff_path,out_path)
