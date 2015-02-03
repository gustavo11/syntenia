# syntenia
Suite of scripts to generate figures depicting synteny and orthology between set of genomes


### Setup
You would need to download three repositories to use **syntenia**

```
$ github clone https://github.com/gustavo11/syntenia
$ github clone https://github.com/gustavo11/GFFLib
$ github clone https://github.com/gustavo11/Orthologia
```

### **gff2graph_ort_projections_no_contraction.pl**
```
usage:
gff2graph_ort_projections_no_contraction.pl <orts file> <list chrom/scaffolds> \
<list gff files> <list fasta files> <svg out>
```

**<orts file>** -  Output of OrthoMCL
<list chrom/scaffolds> - List of genomes of scaffolds that will depicted in the figure
FORMAT:


* **<list gff files>** - List of paths to GFF file containing the annotation
* **<list fasta files>** - List of paths to FASTA file containing the whole genomic sequence
* **<svg out>** - Name of the output
