Class 13: Genome informatics and high throughput sequencing
================
Yi Fu
5/14/2019

## 1\. Galaxy

The Galaxy web-based interface to a suite of bioinformatics tools for
genomic sequence analysis. Galaxy is free and comparatively easy to use.

Galaxy was originally written for genomic data analysis. However, the
set of available tools has been greatly expanded over the years and
Galaxy is now also used for gene expression, genome assembly,
epigenomics, transcriptomics and host of other sub-disciplines in
bioinformatics.

To begin our analysis of this data we will use [Galaxy on
Jetstream](https://galaxyproject.org/tutorials/g101/).

<img src="data/galaxy.png" width="1000"/>

## 2\. RNA-Seq pipeline

1.  Upload our fastqsanger sequences

<img src="data/upload.png" width="1000"/>

You should see this panel on the right telling you the submitted job has
finished.

<img src="data/status.png" width="300"/>

2.  Quality Control

<img src="data/qc.png" width="1000"/>

3.  Mapping RNA-Seq reads to genome

<img src="data/mapping.png" width="1000"/>

4.  Display at UCSC

<img src="data/ucsc1.png" width="1000"/>
<img src="data/ucsc2.png" width="1000"/>

5.  Calculate gene expression

<img src="data/count.png" width="1000"/>
