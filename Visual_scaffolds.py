#!/usr/bin/env python

import os
import sys
import random

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from Bio.Blast.Applications import NcbiblastxCommandline

import cairocffi as cairo



class Annotation(object):
    """
    This class is included to characterise annotation elements. For instance, each gene in the annotation
    will be included as an attribute of class Annotation. The class is run with the line from the annotation 
    file (GFF file) and it split the line into sequence ID, start and end position, chromosome, strain (plus 
    or minus) and type of annotation (ex. Gene, pseudogene, rRNA, rRNA_gene, snoRNA or others). 
    """
    def __init__(self,line):
        temporal_line = line.split('\n')[0].split('\t')
        if temporal_line[2] in included_regions:
            self.seqid = temporal_line[8].split(';')[1].split('=')[1]
            self.start = int(temporal_line[3])
            self.end = int(temporal_line[4])
            self.chromosome = temporal_line[0]
            self.strand = str(temporal_line[6])
            self.type_annotation = str(temporal_line[2])
    def extract_seq(self, genome_dictionary):
        """
        This function uses the start, end position and chromosome to extract the sequence from
        the reference sequence. Addicionally, if the annotation is in the reverse strain (minus)
        it gives the reverse complement sequence instead.
        """ 
        self.seq = genome_dictionary[self.chromosome][(int(self.start)-1):int(self.end)]
        if self.strand == "-":
            self.seq = self.seq.reverse_complement()

class Scaffold(object):
    """
    This class is used to store the de-novo sequence that it is aim to plot
    """
    def __init__(self, pat = None): 
        self.seq_data, self.names = format_seq(pat)
        self.fasta = pat
    def blast_scaffold(self, fasta_reference_seq):
        """
        This method produce a database from a fasta file and use it to blast the de-novo scaffold.
        It uses blastn for that.
        Home:
        /home/sitg/Downloads/ncbi-blast-2.5.0+/bin
        Office:
        /home/sergio/Downloads/ncbi-blast-2.5.0+/bin/
        pat_blast = "/home/sergio/Downloads/ncbi-blast-2.5.0+/bin/"
        """
        line_reference_db = 'makeblastdb -in ' + fasta_reference_seq + ' -dbtype nucl'
        os.system(line_reference_db)
        blastn_cline = NcbiblastxCommandline(cmd="blastn", \
            query=self.fasta, db=fasta_reference_seq, \
            evalue=0.00000001, outfmt='"6 qseqid sseqid pident length mismatch gaps qstart qend sstart send slen evalue bitscore"')
        self.blast_results = str(blastn_cline()[0]).split('\n')
        self.parse_blast_results(self.blast_results)
    def parse_blast_results(self, blast_results):
        """
        This method takes the blast results and parses it include only the best hit in each genetic region.
        The method starts at the begining of the scaffold and for each genetic regions it checks overlapping
        hits and selects the one with the highest alignment length. If two records have the same alignment length
        it will select the one with the lowest number of mismatched.
        """
        pos = 0
        cleanblast = []
        colour_panel = {}
        while pos <= (len(self.seq_data[self.names[0]])-5):
            values_selected_line = [0,0]
            selected_line = ""
            for line in blast_results:
                if line == '':
                    continue
                else:
                    temporal_line = line.split('\t')
                    alignmet_segment = int(temporal_line[3])
                    mismatches = int(temporal_line[4])
                    start = int(min([temporal_line[6],temporal_line[7]]))
                    end = int(max([temporal_line[6],temporal_line[7]]))
                    if (pos + 1) >= start and (pos + 1) < end:
                        if alignmet_segment > values_selected_line[0]:
                            print line
                            selected_line = line
                            values_selected_line = [alignmet_segment, mismatches]
                        elif alignmet_segment == values_selected_line[0] and mismatches < values_selected_line[1]:
                            selected_line = line
                            values_selected_line = [alignmet_segment, mismatches]
            if values_selected_line == [0,0]:
                pos += 100
            else:
                cleanblast.append(selected_line)
                temporal2_line = selected_line.split('\t')
                end_pos = int(max([temporal2_line[6],temporal2_line[7]]))
                pos = end_pos + 2
                id_annotation = temporal2_line[1]
                if not id_annotation in colour_panel:
                    colour_panel[id_annotation] = random_color()
        self.cleanblast = cleanblast
        self.colour_panel = colour_panel
    def plot_scaffold(self, output_file, scale):
        """
        This method takes the filteres blast results and creates a figure with the blasted annotations.
        It shows the direction of the hits, the name IDs and the proportion of the annotation aligning to the
        de-novo scaffold.
        It save the figure in the file "output_file" specified.
        """
        scaffold_len = float(len(self.seq_data[self.names[0]]))
        scaffold_unids = float(1000)/scaffold_len
        
        surface = cairo.SVGSurface(output_file, 1100, 100)
        cr = cairo.Context(surface)
        # Sets color to grey
        cr.set_source_rgb(0.6, 0.6, 0.6)
        # Draws two rectangles and fills them with grey
        cr.rectangle(50, 47.5, 1000, 5)
        cr.fill()
        cr.set_font_size(5)
        for line in self.cleanblast:
            temporal_line = line.split('\t')
            id_seq = temporal_line[1]
            if int(temporal_line[8]) < int(temporal_line[9]):
                start = int(float(temporal_line[6])*scaffold_unids)
                end = int(float(temporal_line[7])*scaffold_unids)
            else:
                start = int(float(temporal_line[7])*scaffold_unids)
                end = int(float(temporal_line[6])*scaffold_unids)        
            v1, v2, v3 = self.colour_panel[id_seq]
            cr.set_source_rgb(v1, v2, v3)
            cr.move_to((start+50),60)
            cr.line_to((end+50),55.5)
            cr.line_to((end+50),45.5)
            cr.line_to((start+50),40)
            cr.line_to((start+50),60)
            cr.fill()
            label_location = ((start+50) + (end+50))/2
            cr.move_to(label_location,35)
            cr.set_source_rgb(0, 0, 0)
            cr.rotate(200)
            cr.show_text(id_seq)
            cr.rotate(-200)
            label_coverage = "{0:.2f}".format(float(temporal_line[3])/float(temporal_line[10]))
            cr.move_to(label_location,80)
            cr.rotate(200)
            cr.show_text(label_coverage)
            cr.rotate(-200)
        cr.set_source_rgb(0, 0, 0)
        cr.rectangle(50, 85, int(scaffold_unids*scale), 1)
        cr.fill()
        cr.move_to(((50+(scaffold_unids*scale))/2),85)
        cr.show_text("10 kb")
        surface.finish()
    def blast_scaffold_multiple_seq(self, fasta_reference_seq):
        """
        This method produce a database from a fasta file and use it to blast the de-novo scaffold.
        It uses blastn for that.
        Home:
        /home/sitg/Downloads/ncbi-blast-2.5.0+/bin
        Office:
        /home/sergio/Downloads/ncbi-blast-2.5.0+/bin/
        pat_blast = "/home/sergio/Downloads/ncbi-blast-2.5.0+/bin/"
        """
        line_reference_db = 'makeblastdb -in ' + fasta_reference_seq + ' -dbtype nucl'
        os.system(line_reference_db)
        blastn_cline = NcbiblastxCommandline(cmd="blastn", \
            query=self.fasta, db=fasta_reference_seq, \
            evalue=0.00000001, outfmt='"6 qseqid sseqid pident length mismatch gaps qstart qend sstart send slen evalue bitscore"')
        self.blast_results = str(blastn_cline()[0]).split('\n')
        self.parse_blast_results_multiple_seq(self.blast_results)
    def parse_blast_results_multiple_seq(self, blast_results):
        """
        This method takes the blast results and parses it include only the best hit in each genetic region.
        The method starts at the begining of the scaffold and for each genetic regions it checks overlapping
        hits and selects the one with the highest alignment length. If two records have the same alignment length
        it will select the one with the lowest number of mismatched.
        """
        colour_panel = {}
        cleanblast = []
        for contig in self.names:
            pos = 0
            while pos <= (len(self.seq_data[contig])-5):
                values_selected_line = [0,0]
                selected_line = ""
                for line in blast_results:
                    if line == '':
                        continue
                    else:
                        temporal_line = line.split('\t')
                        if temporal_line[0] == contig:
                            alignmet_segment = int(temporal_line[3])
                            mismatches = int(temporal_line[4])
                            start = int(min([temporal_line[6],temporal_line[7]]))
                            end = int(max([temporal_line[6],temporal_line[7]]))
                            if (pos + 1) >= start and (pos + 1) < end:
                                if alignmet_segment > values_selected_line[0]:
                                    print line
                                    selected_line = line
                                    values_selected_line = [alignmet_segment, mismatches]
                                elif alignmet_segment == values_selected_line[0] and mismatches < values_selected_line[1]:
                                    selected_line = line
                                    values_selected_line = [alignmet_segment, mismatches]
                if values_selected_line == [0,0]:
                    pos += 100
                else:
                    cleanblast.append(selected_line)
                    temporal2_line = selected_line.split('\t')
                    end_pos = int(max([temporal2_line[6],temporal2_line[7]]))
                    pos = end_pos + 2
                    id_annotation = temporal2_line[1]
                    if not id_annotation in colour_panel:
                        colour_panel[id_annotation] = random_color()
        self.cleanblast = cleanblast
        self.colour_panel = colour_panel
    def plot_scaffold_multiple_seq(self, output_file):
        """
        This method takes the filteres blast results and creates a figure with the blasted annotations.
        It shows the direction of the hits, the name IDs and the proportion of the annotation aligning to the
        de-novo scaffold.
        It save the figure in the file "output_file" specified.
        """
        maxlength = max(len(s) for s in self.seq_data.values())
        scaffold_len = float(maxlength)
        scaffold_unids = float(1000)/scaffold_len
        number_contigs = int(len(self.names))
        
        surface = cairo.SVGSurface(output_file, 1100, number_contigs*100)
        cr = cairo.Context(surface)
        contig_number = 0
        for contig in self.names:
            base_line = (contig_number*100)+50
            contig_number += 1
            # Sets color to grey
            cr.set_source_rgb(0.6, 0.6, 0.6)
            # Draws two rectangles and fills them with grey
            cr.rectangle(50, base_line-2.5, float(len(self.seq_data[contig]))*scaffold_unids, 5)
            cr.fill()
            cr.set_font_size(5)
            for line in self.cleanblast:
                temporal_line = line.split('\t')
                if temporal_line[0] == contig:
                    id_seq = temporal_line[1]
                    if int(temporal_line[8]) < int(temporal_line[9]):
                        start = int(float(temporal_line[6])*scaffold_unids)
                        end = int(float(temporal_line[7])*scaffold_unids)
                    else:
                        start = int(float(temporal_line[7])*scaffold_unids)
                        end = int(float(temporal_line[6])*scaffold_unids)        
                    v1, v2, v3 = self.colour_panel[id_seq]
                    cr.set_source_rgb(v1, v2, v3)
                    cr.move_to((start+50),base_line+10)
                    cr.line_to((end+50),base_line+5.5)
                    cr.line_to((end+50),base_line-5.5)
                    cr.line_to((start+50),base_line-10)
                    cr.line_to((start+50),base_line+10)
                    cr.fill()
                    label_location = ((start+50) + (end+50))/2
                    cr.move_to(label_location,base_line-15)
                    cr.set_source_rgb(0, 0, 0)
                    cr.rotate(200)
                    cr.show_text(id_seq)
                    cr.rotate(-200)
                    label_coverage = "{0:.2f}".format(float(temporal_line[3])/float(temporal_line[10]))
                    cr.move_to(label_location,base_line+30)
                    cr.rotate(200)
                    cr.show_text(label_coverage)
                    cr.rotate(-200)
            cr.set_source_rgb(0, 0, 0)
            cr.rectangle(50, base_line+35, int(scaffold_unids*100000), 1)
            cr.fill()
            cr.move_to(((50+(scaffold_unids*100000))/2),base_line+35)
            cr.show_text("100 kb")
        surface.finish()

def format_seq(fasta_file):
    """ 
    This function is to process the fasta file. It opens the file and re-format all sequences.
    It produces a dictionary with sequences' names (as key) and their sequences (as values).
    Additionally, it produces a list with the sequences names. Although this list is not used 
    in this version of the tool, this list will be implemented in the future to be able to 
    index the sequences.
    """
    seq_dictionary = {}
    list_names = []
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        seq_dictionary[seq_record.id] = seq_record.seq
        list_names.append(seq_record.id)
    return seq_dictionary, list_names

def annotation_Parser(gff_file, reference_fa, included_regions,output_file = None):
    """
    This function takes a annotation file (GFF) and parsers it to create attributes of class
    Annotation. The output is a dictionary with all entries from the GFF file.
    If an output file is given, it creates additionally a fasta file with all parsed annotation
    and save it with that name.
    """
    annotation_gff = open(gff_file, 'r')
    reference_sequences, chromosome_names = format_seq(reference_fa)
    annotation = {}
    fasta_sequences = []
    for line in annotation_gff:
        if line[0] == str('#'):
            continue
        else:
            temporal_line = line.split('\n')[0].split('\t')
            if temporal_line[2] in included_regions:
                id_annotation_ref = temporal_line[8].split(';')[0].split('=')[1]
                annotation[id_annotation_ref] = Annotation(line)
                annotation[id_annotation_ref].extract_seq(reference_sequences)
                record = SeqRecord(Seq(str(annotation[id_annotation_ref].seq)), id=annotation[id_annotation_ref].seqid)
                fasta_sequences.append(record)
    if output_file != None:
        produce_fasta(fasta_sequences, output_file)
    return annotation

def produce_fasta(fasta_sequences, output):
    """
    this function is to write a fasta file with all sequences which will be used for the blast step.
    """
    output_handle = open(output, "w")
    SeqIO.write(fasta_sequences, output, "fasta")
    output_handle.close()

def random_color():
    """
    This function is to create random colours which will be used to plot the annotations
    """
    return [random.random() for _ in range(3)]
