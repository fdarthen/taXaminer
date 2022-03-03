#!/usr/bin/env python

"""Extract protein sequence for genes based on GFF and FASTA file

Compute protein FASTA by parsing gene sequences baed on coordinates from
GFF and FASTA with scaffold sequences (assembly FASTA) and converting
it to amino acid sequence

Expects processed config file
"""
__author__ = "Freya Arthen"
__version__ = "0.6.0"

import sys
from Bio import Seq
import logging

import prepare_and_check
import taxonomic_assignment

class Feature:
    """
    Objetc for GFF feature

    Attributes
    ----------
    contig : str

    id : str

    parent : str

    start : str

    end : str

    strand : str

    phase : str

    biotype : str

    transl_table : str

    children : str

    parsed_parent : str

    transcripts : str

    cdss : str

    coordinates : str

    """
    def __init__(self, contig, info_dict):

        self.contig = contig
        self.id = info_dict.get('ID')
        self.parent = info_dict.get('Parent')
        self.start = info_dict.get('start')
        self.end = info_dict.get('end')
        self.strand = info_dict.get('strand')
        self.phase = info_dict.get('phase')
        self.biotype = info_dict.get('biotype') if info_dict.get('biotype') else info_dict.get('gene_biotype')
        self.transl_table = info_dict.get('transl_table')

        self.children = {}
        self.parsed_parent = None

        # gene only
        self.transcripts = {} # key = transcript ID, value = CDS length
        self.cdss = {} # key = transcript ID, value = [CDS IDs]
        self.coordinates = [] # coordinates of CDS of currently longest transcript

def get_cds_coordinates(cfg):
    """ Parse GFF to find gene and their corresponding mRNA and/or CDS

    Parses the GFF and, while considering the special parsing options of
    taXamine, gathers information on genes, mRNA and CDS. Saves for each
    contig the list of genes IDs that are annotated on it. CDS and transcripts
    features are connected to genes. Transcript Feature object are stored
    together with genes in features dictionary; CDS features are saved in
    gene object.

    Args:
      cfg(obj): Config class object of config parameters

    Returns:
        contigs(dict): {scaffold ID : [gene IDs]}
        features(dict): {GFF feature ID : Feature object}
    """

    contigs = {} #key = contig, value = [gene ID]
    features = {} #key = feature ID, value = feature object


    with open(cfg.gff_path, "r") as gff_file:
        parent = None
        for line in gff_file:
            if not line.startswith("#"):
                spline = line.strip().split("\t")
                # only read lines that should be considered based on cfg.gff_source
                # if default is selected, read all lines
                if cfg.gff_source == "default"  or spline[1] == cfg.gff_source:
                    if spline[2] in [cfg.gff_gene_tag, cfg.gff_transcript_tag,
                                    cfg.gff_cds_tag]:
                        #gather information for GFF entries of type gene, mRNA and CDS
                        contig = spline[0]
                        info_dict={"start": int(spline[3]), "end": int(spline[4]), "strand": spline[6], "phase": spline[7]}
                        # save each attribute separately in dict
                        for elem in spline[8].split(";"):
                            # check if attributes are empty/faulty
                            if len(elem.split("=")) == 2:
                                key, value = elem.split("=")
                                info_dict[key] = value
                        feature = Feature(contig, info_dict)
                        if spline[2] == cfg.gff_gene_tag: # gene
                            # if biotype tag is available, use only genes with protein_coding tag
                            if (feature.biotype and feature.biotype == "protein_coding") or not feature.biotype:
                                features[feature.id] = feature
                                if contig in contigs.keys():
                                    contigs[contig].append(feature.id)
                                else:
                                    contigs[contig] = [feature.id]
                        else: # mRNA or CDS
                            # get Feature object of parent for current feature
                            if cfg.gff_gene_connection == 'parent/child' or cfg.gff_source in ['default', 'maker', 'augustus_masked']:
                                parent_attr = cfg.gff_parent_child_attr.get('parent').lower() if cfg.gff_parent_child_attr else 'parent'
                                parent = features.get(getattr(feature,parent_attr))
                            elif cfg.gff_gene_connection == 'inline':
                                parent = features.get(getattr(feature,cfg.gff_gene_attr.lower()))
                            else:
                                logging.error('Please check GFF parsing rule. Unclear how to connect gene and transcript')
                            if parent:
                                # which type is parent feature: gene or transcript?
                                feature.parsed_parent = parent.id
                                if spline[2] == cfg.gff_transcript_tag:
                                    # current feature -> transcript
                                    # parent -> gene
                                    features[feature.id] = feature
                                    parent.transcripts[feature.id] = 0
                                    parent.cdss[feature.id] = []
                                else:
                                    # current feature -> CDS
                                    # parent -> gene or transcript
                                    if parent.parsed_parent:
                                        # parsed_parent of parent is already assigned
                                        # parent -> transcript
                                        gene = features.get(parent.parsed_parent)
                                        transcript_id = parent.id
                                    else:
                                        # parent -> gene
                                        # thus transcripts are not specified and
                                        # no transcripts are available
                                        # current feature -> transcript and CDS
                                        gene = parent
                                        transcript_id = gene.id
                                        parent.transcripts[transcript_id] = 0
                                        if not parent.cdss.get(transcript_id):
                                            parent.cdss[transcript_id] = []
                                    gene.cdss[transcript_id].append(feature)
            elif "#FASTA" in line: # if FASTA block has been reached
                break
    return contigs, features

def get_longest_transcript(contigs, features):
    """Identify the transcript with the longest CDS for each gene.

    Computes length of each transcript associated with gene Features
    and stores the coordinates for the respective CDS. Also corrects
    for phase, adds strand and translational table information

    Args:
      contigs(dict): {scaffold ID : [gene IDs]}
      features(dict): {GFF feature ID : Feature object}
    """

    for contig, c_genes in contigs.items():
        cds = None
        for gene_id in c_genes:
            gene = features.get(gene_id)
            for transcript_id, cdss in gene.cdss.items():
                length = 0
                for cds in cdss:
                    length += cds.end-cds.start+1
                gene.transcripts[transcript_id] = length
            if not gene.transcripts:
                # no transcript was found for gene
                continue
            # identify transcript with longest CDS
            max_len_t = max(gene.transcripts, key=gene.transcripts.get)
            # get coordinates for that CDS
            for cds in gene.cdss.get(max_len_t):
                gene.coordinates.append((cds.start-1,cds.end,cds.phase))

            if not cds:
                continue

            if cds.strand == '+':
                gene.strand = False # CDS is on forward strand
            else:
                gene.strand = True
            phase = int(gene.coordinates[0][2]) if gene.coordinates[0][2] != '.' else 0
            # sort coordinates
            gene.coordinates = sorted(gene.coordinates, key=lambda k: k[0], reverse=gene.strand)
            # correct first coordinate for phase
            gene.coordinates[0] = (gene.coordinates[0][0]+phase,
                                gene.coordinates[0][1], gene.coordinates[0][2])

            # add translational table information if available
            gene.transl_table = cds.transl_table
            if not gene.transl_table:
                gene.transl_table = '1'


def set_seqs(proteins_file, contigs, features, current_contig, contig_seq):
    """Extract nuc sequences and translate to AA for genes on contig.

    Use coordinates of longest transcript in gene Feature object
    to extract the gene sequence and translate it into AA sequence

    Args:
      proteins_file(obj): file object for protein FASTA file
      contigs(dict): {scaffold ID : [gene IDs]}
      features(dict): {GFF feature ID : Feature object}
      current_contig(str): ID of currrent contig
      contig_seq(str): sequence for current_contig
    """

    for gene_id in contigs.get(current_contig):
        gene = features.get(gene_id)

        if not gene.transcripts:
            # no transcript was found for gene
            continue

        proteins_file.write(">"+gene_id+"\n")
        seq = ''
        for coordinate in gene.coordinates:
            # on reverse strand CDS needs to be reverse complemented before concatenated
            if gene.strand:
                seq += str(Seq.Seq(contig_seq[coordinate[0]:coordinate[1]]).reverse_complement())
            else:
                seq += str(Seq.Seq(contig_seq[coordinate[0]:coordinate[1]]))

        protein = str(Seq.Seq(seq).translate(table=int(gene.transl_table)))

        # if len(seq)%3 != 0:
        #     # prepend N to make sequence length a multiple of three
        #     seq = ("N"*(3-(len(seq)%3)))+seq
        # if "*" in protein_w_stop[:-1]:
            # print(gene.__dict__)
            # print(str(Seq.Seq(seq).translate(table=int(gene.transl_table))))
        # taking only the part to the first stop codon resulted in better matches
        # when using the augustus_masked features from GFF
        # protein = protein.split('*')[0]

        proteins_file.write(protein)
        proteins_file.write("\n")

def extract_seq(cfg, contigs, features):
    """Retrieve sequence for contig from FASTA and write protein FASTA.

    Read assembly FASTA, retrieve scaffold sequence and call set_seqs()
    to write AA sequences for genes to file

    Args:
      cfg(obj): Config class object of config parameters
      contigs(dict): {scaffold ID : [gene IDs]}
      features(dict): {GFF feature ID : Feature object}
    """

    current_contig = ""
    with open(cfg.proteins_path, "w") as proteins_file, \
        open(cfg.fasta_path, "r") as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                # False if no genes in contig dict (geneless contigs)
                if contigs.get(current_contig):
                    set_seqs(proteins_file, contigs, features, current_contig, contig_seq)
                # retrieve ID from current contig
                current_contig = line.strip().lstrip(">").split()[0]
                contig_seq = ''
            else:
                contig_seq += line.strip()
        # add genes from last contig
        if contigs.get(current_contig):
            set_seqs(proteins_file, contigs, features, current_contig, contig_seq)


def generate_fasta(cfg):
    """Generate FASTA file of protein sequences.

    Extract coordinates for genes and transcripts and their relation,
    determine the longest transcript for each gene,
    extract the nucleotide sequence and translate to amino acid sequence

    Args:
      cfg(obj): Config class object of config parameters
    """

    contigs, features = get_cds_coordinates(cfg)
    get_longest_transcript(contigs, features)
    extract_seq(cfg, contigs, features)


def main():
    """Call module directly with preprocessed config file"""
    config_path = sys.argv[1]
    # create class object with configuration parameters
    cfg = prepare_and_check.cfg2obj(config_path)

    generate_fasta(cfg)


if __name__ == '__main__':
    main()
