# write file with longest transcript from gffread proteins file with header = >"transcript ID"\t"gene ID"\t"CDS length"

import sys
import yaml



fun=sys.argv[1]
config_path = sys.argv[2]


config_obj=yaml.safe_load(open(config_path,'r'))
proteins_path=config_obj['proteins_path'] # GFF file path
output_path=config_obj['output_path'] # GFF file path


def gffread_prot_fasta(proteins_path,output_path):
    """ the protein FASTA file was created within MILTS with gffread,
    so it contains in its header the length of the cds
    This information is used to retrieve the longest transcript """

    path_out = output_path+"tmp/tmp.longest_cds.protein.fasta"

    with open(proteins_path, 'r') as file_prot:
        gene_dict = {} #key=geneID, value=sequence
        seq = ''
        add = False
        for line in file_prot:
            if line.startswith('>'):
                if add:
                    gene_dict[gene] = seq
                    seq = '' # reset the sequence
                    add = False
                transcript,gene,length = line.split('\t')
                if gene in gene_dict.keys():
                    # if length of the new transcript is larger, set add bool to True
                    if ((len(gene_dict.get(gene))*3)+3) < int(length): #seq is in peptides (len*3) and gffread does not write final stop peptide (len*3+3)
                        add = True
                    else:
                        add = False
                else: # if no entry for gene is in dict, always add the transcript
                    add = True
            else:
                if add:
                    seq = seq + line.strip()

    # overwrites the input file
    with open(path_out, 'w') as file_out:
        for gene, seq in gene_dict.items():
            file_out.write('>'+gene+'\n')
            file_out.write(seq+'\n')

def regular_prot_fasta(proteins_path,output_path):
    """ the protein FASTA file was provided from outside,
    so it does not contain in its header the length of the cds
    the gene-protein ID matching table, generated with gffread is used for this table
    (it contains the length of the CDS as well) """

    matching_table_path = output_path+"tmp/tmp.prot_gene_matching.txt"
    path_out = output_path+"tmp/tmp.longest_cds.protein.fasta"

    with open(matching_table_path, 'r') as file_match:
        gene_dict = {} #key=geneID, value=(len,transcriptID)
        longest_transcripts = [] # list of transcripts to use (easy access)
        for line in file_match:
            spline = line.strip().split()
            transcriptID, geneID, cdslen = spline[0],spline[1],int(spline[2])
            if geneID in gene_dict.keys():
                if cdslen > gene_dict.get(geneID)[1]: # if new transcript is longer, replace
                    longest_transcripts.remove(gene_dict.get(geneID)[0])
                    gene_dict[geneID] = (transcriptID,cdslen)
                    longest_transcripts.append(transcriptID)
            else:
                gene_dict[geneID] = (transcriptID,cdslen)
                longest_transcripts.append(transcriptID)

    if longest_transcripts[0].startswith("transcript:"):
        prefix = "transcript:"
    else:
        prefix = ""

    # use the list of longest transcript IDs to write to the tmp prot fasta only those
    # proteins whose ID appears in the list
    with open(path_out, 'w') as file_out:
        with open(proteins_path, 'r') as file_prot:
            write_bool = False
            for line in file_prot:
                if line.startswith('>'):
                    transcript = line.split()[0].lstrip(">")
                    if not transcript.startswith(prefix):
                        transcript = prefix+transcript
                    if transcript in longest_transcripts:
                        write_bool = True
                        file_out.write(line)
                    else:
                        write_bool = False
                else:
                    if write_bool:
                        file_out.write(line)



if fun == "gffread":
    gffread_prot_fasta(proteins_path)
else:
    regular_prot_fasta(proteins_path,output_path)
