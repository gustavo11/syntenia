# NAME

gff2graph\_ort\_projections\_no\_contractions.pl - Draws synteny graph

# SYNOPSIS

gff2graph\_ort\_projections\_no\_contractions.pl \[--RBH | --RBH\_CALHOUN | --OMCL \] --orts &lt;file> --chrom\_list &lt;file> --gff\_list &lt;file> --fasta\_list --out &lt;file> \[--help\]

# OPTIONS

**--RBH,--RBH\_Calhoun,--OMCL ** - those flags indicates the format of the orthologous clusters file that will be used when rendering projections. 
See  a description of each format with --help.

**--orts** - orthologous clusters file that will be used when rendering projections. See format with --help.

**--chrom\_list** - list of chromosomes or scaffolds to be rendered. See format with --help.

**--fasta\_list** - file listing the path of GFF files (annotation) associated with each genome that will be rendered 

**--gff\_list** - file listing the path of FASTA files (sequence) associated with each genome that will be rendered 

**--out** - output file in SVG format.

**--help** - print this message **(Optional)**

# DESCRIPTION

**\* RBH format:**

```
828547707      prodigal        org1    transcript1     gene1   None    GAPDH
828547707      prodigal        org2    transcript2     gene2   None    GAPDH
```

**\* RBH\_Calhoun format:**

```
828547707      org1    prodigal        transcript1     gene1   None    GAPDH
828547707      org2    prodical        transcript2     gene2   None    GAPDH
```

**\* OMCL format:**  

```
ORTHOMCL0(3 genes,2 taxa): G001|gene1(G001) G001|gene2(G001) G002|gene3(G001)
ORTHOMCL1(2 genes,1 taxa): G001|gene4(G001) G001|gene5(G001)
```

**\* chrom\_list file format:**

```
<org name on ort. cluster file>\t<scaffold  or chromosome id>(<start>:<end>)<orientation>\t<scaffold  or chromosome id>(<start>:<end>)<orientation>

"(start:end)orientation" are optional. There is no need of them if the whole chrom/scaffold needs to be rendered and if it will be rendered 
in its current orientation. A start or end equals to -1 indicates to the program to render the chrom/scaffold from its first coordinates (start=-1)
till its last coordinate (end=-1). A dash (minus sign) indicates sequences that should be render in the oposite orientation that they are 
reported in the GFF and FASTA fiels  
```

**Ex.:**
 org1 chrom1(1:200)- chrom2(-1:400) chrom3
 org2 chrom1-

# CONTACT

```
Gustavo C. Cerqueira (2015)
cerca11@gmail.com
gustavo@broadinstitute.org
```
