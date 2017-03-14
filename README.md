---
layout: post2
permalink: /informatics_for_high-throughput_data_sequencing_2017_module5_lab/
title: Informatics for High-Throughput Sequencing Data 2017 Module 5 lab
header1: Informatics for High-Throughput Sequencing Data 2017
header2: Module 7 lab
image: CBW_High-throughput_icon.jpg
---


-----------------------

**This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.**

-----------------------

# CBW HT-seq Module 5 - Structural Variant Calling   

 
by Mathieu Bourgey, _Ph.D_



## Table of contents
1. [Introduction](#introduction)
2. [Original Setup](#setup)
2. [Align DNA with BWA-MEM](#align)
3. [Characterize the fragment size distribution](#fragments)
4. [Run DELLY to detect SVs](#delly)
5. [Setting up IGV for SV visualization](#IGV)
6. [Explore the SVs](#explore)
8. [(Optional) Look for other SVs](#otherSV)
9. [Acknowledgements](#ackno)



# CBW HT-seq Module 5 - Structural Variant Calling   

 
by Mathieu Bourgey, _Ph.D_




## Introduction
<a name="introduction"></a>

The goal of this practical session is to identify structural variants (SVs) in a human genome by identifying both discordant paired-end alignments and split-read alignments that. 

Discordant paired-end alignments conflict with the alignment patterns that we expect (i.e., concordant alignments) for the DNA library and sequencing technology we have used. 

For example, given a ~500bp paired-end Illumina library, we expect pairs to align in F/R orientation and we expect the ends of the pair to align roughly 500bp apart. Pairs that align too far apart suggest a potential deletion in the DNA sample's genome. As you may have guessed, the trick is how we define what "too far" is --- this depends on the fragment size distribution of the data. 


Split-read alignments contain SV breakpoints and consequently, then DNA sequences up- and down-stream of the breakpoint align to disjoint locations in the reference genome.


We will be working on a 1000 genome sample, NA12878. You can find the whole raw data on the 1000 genome website:
<http://www.1000genomes.org/data>

The dataset comes from the [Illumina Platinum Genomes Project](http://www.illumina.com/platinumgenomes/), which is a 50X-coverage dataset of the NA12891/NA12892/NA12878 trio. The raw data can be downloaded from the following [URL](http://www.ebi.ac.uk/ena/data/view/ERP001960).

NA12878 is the child of the trio while NA12891 and NA12892 are her parents.

![Pedigree](https://github.com/mbourgey/CBW_HTseq_module5/blob/master/img/Pedigree.png) 


For practical reasons we subsampled the reads from the sample because running the whole dataset would take way too much time and resources.
We're going to focus on the reads extracted from the chromosome 20.

## Original Setup
<a name="setup"></a>

### Amazon node

Read these [directions] (http://bioinformatics-ca.github.io/logging_into_the_Amazon_cloud/) for information on how to log in to your assigned Amazon node. 

### Software requirements

These are all already installed, but here are the original links.

  * [BWA-MEM](http://bio-bwa.sourceforge.net/)
  * [samtools](http://sourceforge.net/projects/samtools/)
  * [delly](https://github.com/tobiasrausch/delly)
  * [bcftools](https://samtools.github.io/bcftools/bcftools.html)
  * [LUMPY](https://github.com/arq5x/lumpy-sv)
  * [R](https://cran.r-project.org/)
  * [python](https://www.python.org/)
  
  
In this session, we will particularly focus on DELLY, a SV detection tool. DELLY is an integrated structural variant prediction method that can discover, genotype and visualize deletions, tandem duplications, inversions and translocations at single-nucleotide resolution in short-read massively parallel sequencing data. It uses paired-ends and split-reads to sensitively and accurately delineate genomic rearrangements throughout the genome. 

If you are interested in DELLY, you can read the full manuscript [here](http://bioinformatics.oxfordjournals.org/content/28/18/i333.abstract).

### Environment setup

```
#set up
export SOFT_DIR=/usr/local/
export WORK_DIR=~/workspace/HTseq/Module5/
export REF=$WORK_DIR/reference


rm -rf $WORK_DIR
mkdir -p $WORK_DIR/SVvariants
cd $WORK_DIR
ln -s ~/CourseData/HT_data/Module5/* .

```

***Note:***
    The `ln -s` command adds symbolic links of all of the files contained in the (read-only) `~/CourseData/HT_data/Module5` directory.
    
### Data files

The initial structure of your folders should look like this:   


```
ROOT
|-- bam/               # bam files (down sampled)
    `-- NA12878/             # Child sample directory
    `-- NA12891/             # Father sample directory
    `-- NA12892/             # Mother sample directory
`-- reference/               # hg19 reference and indexes
`-- scripts/                 # command lines scripts
`-- saved_results/           # precomputed final files
```

### Cheat sheets

* [Unix comand line cheat sheet](http://sites.tufts.edu/cbi/files/2013/01/linux_cheat_sheet.pdf)
* [commands file of this module](https://github.com/mbourgey/CBW_HTseq_module5/blob/master/scripts/commands.sh)



## Align DNA with BWA-MEM
<a name="align"></a>

This step has been done for you in the interest of time, but the commands are shown so that you can reproduce the results later. The advantage of using BWA-MEM in the context of SV discovery is that it produces both paired-end and split-read alignments in a single BAM output file. In contrast, prior to BWA-MEM, one typically had to use two different aligners in order to produce both high quality paired-end and split-read alignments.

In the alignment commands, note the use of the -M parameter to mark shorter split hits as secondary.



```
#################
#Align NA12878 #
#################
#bwa mem $REF/hg19.fa \
#fastq/NA12878_S1.chr20.20X.1.fq \
#fastq/NA12878_S1.chr20.20X.2.fq \
#-M -t 2| samtools view -S -b - \
#> bam/NA12878/NA12878_S1.chr20.20X.pairs.readSorted.bam


#################
#Align NA12891 #
#################
#bwa mem hg19.fa \
#fastq/NA12891_S1.chr20.20X.1.fq \
#fastq/NA12891_S1.chr20.20X.2.fq \
#-M -t 2 | samtools view -S -b - \
#> bam/NA12891/NA12891_S1.chr20.20X.pairs.readSorted.bam


#################
#Align NA12892 #
#################
#bwa mem hg19.fa \
#fastq/NA12892_S1.chr20.20X.1.fq \
#fastq/NA12892_S1.chr20.20X.2.fq \
#-M -t 2 | samtools view -S -b - \
#> bam/NA1289/NA12892_S1.chr20.20X.pairs.readSorted.bam
```


**Why should mark shorter split hits as secondary ?** [solution](https://github.com/mbourgey/CBW_HTseq_module5/blob/master/solutions/_aln1.md)


## Characterize the fragment size distribution
<a name="fragments"></a>

Before we can attempt to identify structural variants via discordant alignments, we must first characterize the insert size distribution 

**What means the insert size distribution ?** [solution](https://github.com/mbourgey/CBW_HTseq_module5/blob/master/solutions/_fragment1.md)

**How can we use the fragment size distribution in SV detection ?** [solution](https://github.com/mbourgey/CBW_HTseq_module5/blob/master/solutions/_fragment2.md)


The following script, taken from LUMPY extracts F/R pairs from a BAM file and computes the mean and stdev of the F/R alignments. It also generates a density plot of the fragment size distribution.

Let's calculate the fragment distribution for the three dataset: 

```
#NA12878
samtools view bam/NA12878/NA12878_S1.chr20.20X.pairs.readSorted.bam \
     | python scripts/pairend_distro.py \
     -r 101 \
     -X 4 \
     -N 10000 \
     -o SVvariants/NA12878_S1.chr20.20X.pairs.histo \
 > SVvariants/NA12878_S1.chr20.20X.pairs.params


#NA12891 
samtools view bam/NA12891/NA12891_S1.chr20.20X.pairs.readSorted.bam \
     | python scripts/pairend_distro.py \
     -r 101 \
     -X 4 \
     -N 10000 \
     -o SVvariants/NA12891_S1.chr20.20X.pairs.histo \
 > SVvariants/NA12891_S1.chr20.20X.pairs.params


#NA12892. 
samtools view bam/NA12892/NA12892_S1.chr20.20X.pairs.readSorted.bam \
     | python scripts/pairend_distro.py \
     -r 101 \
     -X 4 \
     -N 10000 \
     -o SVvariants/NA12892_S1.chr20.20X.pairs.histo \
 > SVvariants/NA12892_S1.chr20.20X.pairs.params
```



Let's take a peak at the first few lines of the histogram file that was produced:

```
head -n 10 SVvariants/NA12878_S1.chr20.20X.pairs.histo
```

Expected results:
>0  &nbsp;&nbsp;&nbsp;&nbsp; 0.0   
>1  &nbsp;&nbsp;&nbsp;&nbsp; 0.000200461060439   
>2  &nbsp;&nbsp;&nbsp;&nbsp; 0.000300691590659   
>3  &nbsp;&nbsp;&nbsp;&nbsp; 0.000300691590659   
>4  &nbsp;&nbsp;&nbsp;&nbsp; 0.000200461060439   
>5  &nbsp;&nbsp;&nbsp;&nbsp; 0.000200461060439   
>6  &nbsp;&nbsp;&nbsp;&nbsp; 0.000300691590659   
>7  &nbsp;&nbsp;&nbsp;&nbsp; 0.00010023053022   
>8  &nbsp;&nbsp;&nbsp;&nbsp; 0.00010023053022   
>9  &nbsp;&nbsp;&nbsp;&nbsp; 0.00010023053022   


Let's use R to plot the fragment size distribution. 



First, launch R from the command line.

```
R
```

Now, within R, execute the following commands:

```
size_dist <- read.table('SVvariants/NA12878_S1.chr20.20X.pairs.histo')
pdf(file = "SVvariants/fragment.hist.pdf") 
layout(matrix(1:3))
plot(size_dist[,1], size_dist[,2], type='h', main="NA12878 insert size") 
size_dist <- read.table('SVvariants/NA12891_S1.chr20.20X.pairs.histo') 
plot(size_dist[,1], size_dist[,2], type='h', main="NA12891 insert size") 
size_dist <- read.table('SVvariants/NA12892_S1.chr20.20X.pairs.histo') 
plot(size_dist[,1], size_dist[,2], type='h', main="NA12892 insert size") 
dev.off()
quit("no")

```

At this point, you should have the following files: 

Let's look at the 3 PDF files : 

open a web browser on your laptop, and navigate to `http://cbwXX.dyndns.info/`, where `XX` is the id of your node. You should be able to find there the directory hierarchy under `~/workspace/` on your node. open ```HTseq/Module5/SVvariants``` folder and open the pdf.



 * NA12878.fragment.hist.pdf
 * NA12891.fragment.hist.pdf
 * NA12892.fragment.hist.pdf


Spend some time thinking about what this plot means for identifying discordant alignments.

**What does the mean fragment size appear to be? Are all 3 graphs the same?** [solution](https://github.com/mbourgey/CBW_HTseq_module5/blob/master/solutions/_insert1.md)


## SV detection
<a name="delly"></a>

For germline SVs, calling is done by sample or in small batches to increase SV sensitivity & breakpoint precision. 

At this point we will only call ***deletions***. Here are the steps adapted from the DELLY [readme](https://github.com/tobiasrausch/delly/blob/master/README.md).



### First call

Let's call deletions:

```
#NA12878
delly call -t DEL -g $REF/hg19.fa -o SVvariants/NA12878.bcf -x $REF/hg19.excl \
  bam/NA12878/NA12878_S1.chr20.20X.pairs.posSorted.bam

#NA12891
delly call -t DEL -g $REF/hg19.fa -o SVvariants/NA12891.bcf -x $REF/hg19.excl \
  bam/NA12891/NA12891_S1.chr20.20X.pairs.posSorted.bam

#NA12892
delly call -t DEL -g $REF/hg19.fa -o SVvariants/NA12892.bcf -x $REF/hg19.excl \
  bam/NA12892/NA12892_S1.chr20.20X.pairs.posSorted.bam
```


`.bcf` files are compressed vcf files. To look at the output you can use bcftools


```
bcftools view SVvariants/NA12878.bcf | less -S
```



***Cheat:*** If these commands are taking too long, simply run the command `cp delly_call/* .`

**How many variant delly found in each sample ?** [solution](https://github.com/mbourgey/CBW_HTseq_module5/blob/master/solutions/_vcf1.md)

### Merge calls

We need to merge the SV sites into a unified site list:

```
delly merge -t DEL -m 500 -n 1000000 -o SVvariants/del.bcf -b 500 -r 0.5 \
  SVvariants/NA12878.bcf SVvariants/NA12891.bcf SVvariants/NA12892.bcf
```

Look at the output:

```
bcftools view SVvariants/del.bcf | less -S
```

**What can you notice different from the individual bcf ?** [solution](https://github.com/mbourgey/CBW_HTseq_module5/blob/master/solutions/_vcf2.md)


### Re-genotype in all samples

We need to re-genotype the merged SV site list across all samples. This can be run in parallel for each sample.

```
#NA12878
delly call -t DEL -g $REF/hg19.fa -v SVvariants/del.bcf -o SVvariants/NA12878.geno.bcf \
  -x $REF/hg19.excl bam/NA12878/NA12878_S1.chr20.20X.pairs.posSorted.bam

#NA12891
delly call -t DEL -g $REF/hg19.fa -v SVvariants/del.bcf -o SVvariants/NA12891.geno.bcf \
  -x $REF/hg19.excl bam/NA12891/NA12891_S1.chr20.20X.pairs.posSorted.bam

#NA12892
delly call -t DEL -g $REF/hg19.fa -v SVvariants/del.bcf -o SVvariants/NA12892.geno.bcf \
  -x $REF/hg19.excl bam/NA12892/NA12892_S1.chr20.20X.pairs.posSorted.bam
```

Look at the output:

```
bcftools view SVvariants/NA12878.geno.bcf | less -S
```

We now have our good candidate list genotype for each individual.


### Merge the new calls

Merge all re-genotyped samples to get a single VCF/BCF using bcftools merge. Also index the resulting file and create vcf file for visualization.

```
bcftools merge -O b -o SVvariants/merged.bcf SVvariants/NA12878.geno.bcf SVvariants/NA12891.geno.bcf SVvariants/NA12892.geno.bcf
bcftools index SVvariants/merged.bcf
bcftools view SVvariants/merged.bcf > SVvariants/merged.vcf
```


### Apply a filter for germline events

```
delly filter -t DEL -f germline -o SVvariants/germline.bcf SVvariants/merged.bcf
bcftools view SVvariants/germline.bcf > SVvariants/germline.vcf
```

**Do you know how to look at the resulting file?** [solution](https://github.com/mbourgey/CBW_HTseq_module5/blob/master/solutions/_vcf3.md)



## Setting up IGV for SV visualization
<a name="IGV"></a>

Launch IGV and load the merged calls and the germline calls using `File -> Load from URL` using:


 * http://cbwXX.dyndns.info/HTseq/Module5/SVvariants/merged.vcf
 * http://cbwXX.dyndns.info/HTseq/Module5/SVvariants/germline.vcf


***Note:*** Once again you will need to replace `XX` by your student number.

Navigate to the following location to see a deletion: `chr20:31,308,410-31,315,294`


Now load the bam files in the same way using:

 * http://cbwXX.dyndns.info/HTseq/Module5/bam/NA12878/NA12878_S1.chr20.20X.pairs.posSorted.bam
 * http://cbwXX.dyndns.info/HTseq/Module5/bam/NA12891/NA12891_S1.chr20.20X.pairs.posSorted.bam
 * http://cbwXX.dyndns.infoHTseq/Module5/bam/NA12892/NA12892_S1.chr20.20X.pairs.posSorted.bam


You should see something like this:

![Deletion](https://github.com/mbourgey/CBW_HTseq_module5/blob/master/img/deletion.png)

You can try to configure IGV such that we can more clearly see the alignments that support the SV prediction.

***Note*** - Do you remember how to make sure that the alignments are colored by insert size and orientation?*
 

## Explore the SVs
<a name="explore"></a>

**Is the variant at chr20:31,310,769-31,312,959 found in each member of the trio?** [solution](https://github.com/mbourgey/CBW_HTseq_module5/blob/master/solutions/_igv1.md)

**What are the genotypes for each member of the trio at this locus (e.g., hemizygous, homozygous)?** [solution](https://github.com/mbourgey/CBW_HTseq_module5/blob/master/solutions/_igv2.md) 

**What about the variant at chr20:37,054,372-37,056,562?** [solution](https://github.com/mbourgey/CBW_HTseq_module5/blob/master/solutions/_igv3.md)

**Does the evidence in the Moleculo track mimic the evidence in the Illumina track for NA12878?** [solution](https://github.com/mbourgey/CBW_HTseq_module5/blob/master/solutions/_igv4.md)

**What about chr20:42,269,896-42,278,072?** [solution](https://github.com/mbourgey/CBW_HTseq_module5/blob/master/solutions/_igv5.md)

Continue exploring the data!



## (Optional) Look for other SVs
<a name="otherSV"></a>

You can try using Delly to look for other types of SVs, for example using:

`delly -t DUP`

`delly -t INV`

`delly -t TRA`


## Acknowledgements
<a name="ackno"></a>

This module is heavily based on a previous module prepared by Aaron Quinlan and Guillaume Bourque. 

