# Visual_Scaffolds

Documents and scripts written by: Sergio Tusso  
Dept. of Evolutionary Biology  
Evolutionary Biology Center (EBC)  
Uppsala University  
Sweden  
2016  
  
Email: situssog@gmail.com  

# Description: 

This tool allows the visualization of genomic regions or scaffolds that were reconstructed de-novo. It is mainly aimed to process regions with repetitive sequences and/or with structural variances compared with a reference sequence. The tool takes the reference sequence and its annotation to extract the sequences of coding regions of sequences of interest.  Then it performs a blast with the de-novo scaffold to infer the location of different genes. It then produces a plot of the genomic region placing genes, the direction and proportion of the gene aligning to the scaffold.

# Dependecies:

- biopython # this is used to process sequences and fasta files
- blast - ncbi # used to create databases with fasta files (makeblastdb) and to blast sequences to scaffolds (blastn)
	- NcbiblastxCommandline
- cairocffi # visualisation package the create final figures

# Input: 
Fasta file with de-novo scaffold/s.  
Fasta file with reference genome and annotation (gff file). Alternatively this can be changes by a fasta file with names and sequences of interest. For instance a fasta file only with transposable sequences.  

  
# Instructions:
Copy the script Visual_scaffolds in the same forder with the other input files.  
The script is run as a module inside python/ipython.  
open python and load the module with the command:  

```python
import Visual_scaffolds
```

## Example using a fasta file with sequences of interes:

define fasta file with secuesces of interest:  

```python
file_reference_sequences = "fragments.fa"
```

open scaffold fasta:  

```python
scaffold_file = "full_MTR_h90.fasta"
scaffold_seq = Scaffold(scaffold_file)
```

blast sequences to scaffold  
example: scaffold_seq.blast_scaffold("fasta_file.fa")  
```python
scaffold_seq.blast_scaffold(file_reference_sequences)
```

plot scaffold with annotation:  
give output file between parentesis.  
example: scaffold_seq.plot_scaffold("output_figure_file.svg")  
```python
scaffold_seq.plot_scaffold("MTR_using_fasta.svg")
```

## Example using reference genome (whole genome in fasta file) and annotation file (GFF file)

Names of sequences in fasta file should be the same as chromosome names in annotation.
The tool only extract and plot 'gene', 'ncRNA_gene', 'pseudogene', 'rRNA', 'rRNA_gene', 
'snoRNA', 'snoRNA_gene', 'snRNA', 'snRNA_gene' and 'tRNA_gene'.  
  
If it is required to exclude some of them, delete them from the next list:  

```python
included_regions = ['gene', 'ncRNA_gene', 'pseudogene', 'rRNA', 'rRNA_gene', 'snoRNA', \
'snoRNA_gene', 'snRNA', 'snRNA_gene', 'tRNA_gene']
```

define annotation file and reference sequence:  
```python
annotation_file = "schizosaccharomyces_pombe.chr.gff3"
reference_genome_file = "PombeRef.fa"
```

create annotation dictionary:  
```python
annotation_dic = annotation_Parser(annotation_file, reference_genome_file, included_regions, "reference_parsed_gff.fa")
```

open scaffold fasta  
```python
scaffold_file = "full_MTR_h90.fasta"
scaffold_seq = Scaffold(scaffold_file)
```

blast sequences to scaffold  
example: scaffold_seq.blast_scaffold("fasta_file.fa")  
```python
scaffold_seq.blast_scaffold("reference_parsed_gff.fa")
```

plot scaffold with annotation:  
give output file between parentesis  
example: scaffold_seq.plot_scaffold("output_figure_file.svg")  
```python
scaffold_seq.plot_scaffold("MTR.svg")
```
