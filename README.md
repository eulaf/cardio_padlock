# cardio_padlock
Primer calling scripts for cardio padlock probes.

These scripts are designed to call padlock probe primers to cover all exons in a list of genes.

## Requirements
These scripts require the following programs to run:

* Perl 5
* [Ensembl API](http://www.ensembl.org/info/docs/api/index.html)
* [Primer3](http://primer3.sourceforge.net)
* [RepeatMasker](http://www.repeatmasker.org) (optional)
* [NCBI standalone BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (optional)

## Steps

#### 1. Download sequences for genes of interest
Create a file with the list of genes.  The genes can be designated using the 
gene name or Ensembl gene ID.  Run get_exons_fa.pl to get the gene sequences
and a list of SNPs overlapping the gene (this script requires the
Ensembl API and an internet connection).  For example,

```
get_exons_fa.pl --maf 0.1 genelist.txt
```

will retrieve sequences for all genes listed in genelist.txt and all SNPs
in the area of the sequences with a MAF >= 0.1.

Or if you only have a few genes:

```
get_exons_fa.pl --maf 0.1 GENE1 GENE2 ...
```

The script will create a series of files for each gene:

* GENE.exons.fa -- FASTA file of the merged exons for gene GENEID.
The merged exons should cover all exons from all transcripts of the
gene.

* GENE.exon_summary.txt -- tab-delimited file describing all gene
transcripts and exons and how they were merged.

* GENE.gene.fa -- FASTA file of the gene sequence with flanking 
sequence, merged exons and introns appearing as separate sequences.    

* GENE.snps.txt -- list of all SNPs in the area of the gene.

#### 2. Mask sequences (optional)
Run RepeatMasker on the \*.gene.fa files created previously to mask
common sequences from the primer caller.  

#### 3. Call primers
Call primers on the masked sequences (\*.masked) or gene sequences 
(\*.gene.fa) if not using RepeatMasker.  The script will attempt to call 
primers to cover all exons plus a 20bp flank around each exon.  Exons larger 
than the specified product size will be divided into smaller segments and 
primers called on each segment.

```
call_primers -label myRun *.gene.fa.masked
```

The script will create a series of files for each FASTA file:

* myRun.GENE.all_amplicons.fa -- FASTA file of amplicons for all 
primer pairs

* myRun.GENE.all_primers.fa -- FASTA file of primer sequences for
all primers called

* myRun.GENE.all_primers.txt -- tab-delimited file listing all primer pairs 
and their sequences, amplicon sequence, location, and other info.

* myRun.GENE.primers_found.txt -- tab-delimited file listing the number of 
primers found for each segment for all exons of GENE.  

* myRun.primers_found.txt -- tab-delimited file listing the number of primers 
found for each segment for all exons for all genes.

* myRun.primers_summary.txt -- tab-delimited file listing number of segments 
used to cover each exon and number of segments not covered by primers for all 
genes.

#### 4. BLAST primers against human genome (optional)
Run NCBI standalone BLAST of the primers against the human genome to enable
filtering of primers by number of BLAST hits.

For each \*.all_primers.fa file, run:

```
blastn -db human_genome.fasta -evalue 0.05 -max_target_seqs 5 
-outfmt 7 -task blastn -query myRun.GENE.all_primers.fa 
-out myRun.GENE.all_primers.blastn.txt
```

Then add the BLAST results to the \*.all_primers.txt files by running:

```
add_blast.pl *.all_primers.txt *.blastn.txt
```

The script adds two columns to each \*.all_primers.txt file: 

* num_blast_hits_left  -- number of BLAST hits of for the left primer
* num_blast_hits_right  -- number of BLAST hits of for the right primer

The new files have the name format myRun.GENE.all_primers.blast_counts.txt.


#### 5. Filter primers
Finally, to select the best primers for each segment and remove redundant
primers, run:

```
filter_primers.pl -label myRun *.all_primers.blast_counts.txt
```

or if BLAST was not run:

```
filter_primers.pl -label myRun *.all_primers.txt
```

The script creates the following files:

* myRun.GENE.filtered_amplicons.fa -- FASTA file of amplicon sequences for selected primer pairs
* myRun.GENE.filtered_primers.fa -- FASTA file of primer sequences for selected primer pairs
* myRun.GENE.filtered_primers.txt -- tab-delimited text file describing each primer pair chosen
* myRun.filtered_primers.txt, myRun.filtered_amplicons.fa, myRun.filtered_primers.fa -- compilation of all GENE files into one file

