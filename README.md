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

* **\<orts file\>** -  Output of OrthoMCL
* **\<list chrom/scaffolds\>** - List of genomes of scaffolds that will depicted in the figure. See format below
* **\<list gff files\>** - List of paths to GFF file containing the annotation
* **\<list fasta files\>** - List of paths to FASTA file containing the whole genomic sequence
* **\<svg out\>** - Name of the output


##### FORMAT of chrom/scaffold list
```
<org name on Ort cluster file> <scaffold/chrom>(start:end)orientation ...
```
(start:end)orientation" are optional. No need of them if the whole chrom/scaffold will be rendered or if it will be rendered in its current orientation

**Example:**
```
Loa_loa_V2 7000000145608817(1:200)- 7000000145609071 7000000145608793
C_elegans_WS224 7000000183869869-
```


