# syntenia
Suite of scripts to generate figures depicting synteny and orthology between set of genomes.

[Example of figure generated by **syntenia** scripts](http://www.google.com/imgres?imgurl=http://www.nature.com/ng/journal/v45/n5/images/ng.2585-F1.jpg&imgrefurl=http://www.nature.com/ng/journal/v45/n5/fig_tab/ng.2585_F1.html&h=350&w=946&tbnid=nc6NCfH-xQMlOM:&zoom=1&docid=426ljHgTbrBRQM&ei=aijRVMjxL4TgggT11IHABw&tbm=isch&ved=0CB4QMygAMAA)


### Setup
You will need to download three repositories to use **syntenia**

```
$ github clone https://github.com/gustavo11/syntenia
$ github clone https://github.com/gustavo11/GFFLib
$ github clone https://github.com/gustavo11/Orthologia
```

-----

# Scripts

### **gff2graph_ort_projections_no_contraction.pl**
```
usage:
gff2graph_ort_projections_no_contractions.pl [--RBH | --RBH_CALHOUN | --OMCL ] --orts <file> --chrom_list <file> --gff_list <file> --fasta_list --out <file> [--help]

```
* **--RBH,--RBH_Calhoun,--OMCL** - those flags indicate the format of the orthologous clusters file that will be used when rendering projections. 
See a description of each format below.

* **--orts** - orthologous clusters file that will be used when rendering projections. See format below.

* **--chrom_list** - list of chromosomes or scaffolds to be rendered. See format below.

* **--fasta_list** - file listing the path of GFF files (annotation) associated with each genome that will be rendered 

* **--gff_list** - file listing the path of FASTA files (sequence) associated with each genome that will be rendered 

* **--out** - output file in SVG format.

* **--help** - print this message **(Optional)**

<BR>
<BR>

##### FORMATS

* **RBH format:**

```
 828547707	prodigal	org1	transcript1	gene1	None	GAPDH
 828547707	prodigal	org2	transcript2	gene2	None	GAPDH
```

* **RBH_Calhoun format:**

```  
 828547707	org1	prodigal	transcript1	gene1	None	GAPDH
 828547707	org2	prodical	transcript2	gene2	None	GAPDH
```

* **OMCL format:**  

```
 ORTHOMCL0(3 genes,2 taxa): G001|gene1(G001) G001|gene2(G001) G002|gene3(G001)
 ORTHOMCL1(2 genes,1 taxa): G001|gene4(G001) G001|gene5(G001)
```

* **chrom_list file format:**

```
 <org name on ort. cluster file>\t<scaffold  or chromosome id>(<start>:<end>)<orientation>\t<scaffold  or chromosome id>(<start>:<end>)<orientation>
```

 "(start:end)orientation" are optional. There is no need of them if the whole chrom/scaffold needs to be rendered and if it will be rendered 
 in its current orientation. A start or end equals to -1 indicates to the program to render the chrom/scaffold from its first coordinates (start=-1)
 till its last coordinate (end=-1). A dash (minus sign) indicates sequences that should be render in the oposite orientation that they are 
 reported in the GFF and FASTA fiels  

**Ex.:**
```
 org1 chrom1(1:200)- chrom2(-1:400) chrom3
 org2 chrom1-
```




