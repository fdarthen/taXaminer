# write file with longest transcript from gffread proteins file with header = >"transcript ID"\t"gene ID"\t"CDS length"

import yaml
config_obj=yaml.safe_load(open('./config.yml','r'))
proteins_path=config_obj['proteins_path'] # GFF file path


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
with open(proteins_path, 'w') as file_out:
    for gene, seq in gene_dict.items():
        file_out.write('>'+gene+'\n') 
        file_out.write(seq+'\n')
