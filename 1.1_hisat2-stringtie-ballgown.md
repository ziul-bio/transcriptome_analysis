# Transcriptome analysis with Hisat2 - stringtie - ballgown

#### Autor: Luiz Carlos Vieira
#### Date: 25/03/2022

The aim of this analysis is to quantify the expression level of a long non-coding Rna, which
has been characterizaed by the lab.

Then with the lncRNA sequence in hands, we can try to find the expression level of this transcript, and verify 
the differential expression between groups.


## Data

Status:	Public on Sep 17, 2019

Title:	Phenotypic plasticity shapes genome architecture in the honeybee Apis mellifera (RNA-seq)

Organism:	Apis mellifera

Experiment type:	Expression profiling by high throughput sequencing

Summary: Female honeybees are specified as workers or queens based on diet during early development. Workers are essentially sterile with a reduced number of ovarioles and no spermatheca. In the presence of the queen (queen mandibular pheromone) and her brood, worker ovaries are kept in an inactive quiescent state. If the queen is removed, or lost, worker bees are able to sense this change in their environment and their ovaries undergo complete remodeling producing unfertilized haploid eggs that will produce male (drone bees). In this study we analyze gene expression in queen, worker, and laying worker ovaries using RNA-seq and explore differences in the chromatin landscape (focusing on H3K27me3).

Overall design:	RNA-seq to measure gene expression in queen, queen right worker and actively laying worker ovaries

Contributor(s):	Duncan EJ, Leask MP, Dearden PK

Citation(s): Leask M, Lovegrove M, Walker A, Duncan E et al. Evolution and genomic organization of the insect sHSP gene cluster and coordinate regulation in phenotypic plasticity. BMC Ecol Evol 2021 Aug 4;21(1):154. PMID: 34348652


Downloading the data from Series GSE120561 [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120561)
```bash
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR790/006/SRR7908186/SRR7908186.fastq.gz -o SRR7908186_W1_Worker_Pool_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR790/007/SRR7908187/SRR7908187.fastq.gz -o SRR7908187_W2_Worker_Pool_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR790/008/SRR7908188/SRR7908188.fastq.gz -o SRR7908188_AW1_Active_Pool_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR790/009/SRR7908189/SRR7908189.fastq.gz -o SRR7908189_AW2_Active_Pool_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR790/000/SRR7908190/SRR7908190.fastq.gz -o SRR7908190_Q1_Queen_Pool_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR790/001/SRR7908191/SRR7908191.fastq.gz -o SRR7908191_Q2_Queen_Pool_2.fastq.gz
```


## Find strandness of reads
```bash
# subset fastq files
zcat SRR7908186_W1_Worker_Pool_1.fastq.gz | head -n4000 > teste/test.fastq

# Creating a kallisto index
kallisto index --index="apis" ref_genome/genome_Apis.fa
```


## Running kallisto quant
```bash
KALLISTO_INDEX=/mnt/c/Users/luiz_/Downloads/transcriptomAnalysis/teste
OUT_DIR=/mnt/c/Users/luiz_/Downloads/transcriptomAnalysis/teste
SEQ_DIR=/mnt/c/Users/luiz_/Downloads/transcriptomAnalysis/teste

kallisto quant -i ${KALLISTO_INDEX} -o ${OUT_DIR}/test.un ${SEQ_DIR}/test_1.fg.gz ${SEQ_DIR}/test_2.fq.gz
kallisto quant -i ${KALLISTO_INDEX} -o ${OUT_DIR}/test.rf ${SEQ_DIR}/test_1.fg.gz ${SEQ_DIR}/test_2.fq.gz --rf-stranded
kallisto quant -i ${KALLISTO_INDEX} -o ${OUT_DIR}/test.fr ${SEQ_DIR}/test_1.fg.gz ${SEQ_DIR}/test_2.fq.gz --fr-stranded

paste ${OUT_DIR}/test.fr/abundance.tsv ${OUT_DIR}/test.rf/abundance.tsv ${OUT_DIR}/test.un/abundance.tsv | cut -f1,4,9,14  | awk 'BEGIN{sum1=0;sum2=0;sun3=0}{sum1+=$2;sum2+=$3;sum3+=$4}END{print sum1,sum2,sum3}' > ${OUT_DIR}/test.libtypetesting
less ${OUT_DIR}/test.libtype.txt | awk '{print $2/$1,$3/$1,$3/$2}' | awk '{if($1<0.3 && $3>3)print "stranded";else if($1>3 && $2>3)print "reverse";else print "unstranded"}' >> ${OUT_DIR}/test.libtype.txt

```


## Alignment with hisat2

First, using the python scripts included in the HISAT2 package, extract splice-site and exon information from the gene annotation file:
```bash
extract_splice_sites.py ref_genome/annotation_Apis.gtf > hisat2_index/apis.ss
extract_exons.py ref_genome/annotation_Apis.gtf > hisat2_index/apis.exon

# Second, build a HISAT2 index:
hisat2-build -p 4 --ss ref_genome/hisat2_index/apis.ss --exon ref_genome/hisat2_index/apis.exon ref_genome/genome/genome_apis.fa ref_genome/hisat2_index/apis

```


Align all RNA-seq data sets to the reference genome using hisat2. 
```bash
hisat2 --dta -p 4 -t -x ref_genome/hisat2_index/apis -U rawFastq/SRR7908186_W1_Worker_Pool_1.fastq.gz | samtools sort -o hisat2/w1.bam -

hisat2 --dta -p 4 -t -x ref_genome/hisat2_index/apis -U rawFastq/SRR7908187_W2_Worker_Pool_2.fastq.gz | samtools sort -o hisat2/w2.bam -

hisat2 --dta -p 4 -t -x ref_genome/hisat2_index/apis -U rawFastq/SRR7908188_AW1_Active_Pool_1.fastq.gz | samtools sort -o hisat2/aw1.bam -

hisat2 --dta -p 4 -t -x ref_genome/hisat2_index/apis -U rawFastq/SRR7908189_AW2_Active_Pool_2.fastq.gz | samtools sort -o hisat2/aw2.bam -

hisat2 --dta -p 4 -t -x ref_genome/hisat2_index/apis -U rawFastq/SRR7908190_Q1_Queen_Pool_1.fastq.gz | samtools sort -o hisat2/q1.bam -

hisat2 --dta -p 4 -t -x ref_genome/hisat2_index/apis -U rawFastq/SRR7908191_Q2_Queen_Pool_2.fastq.gz | samtools sort -o hisat2/q2.bam -

```



## Stringtie

StringTie automatically detects new genes and new isoforms if they are present in your experimental data, 
regardless of whether or not they appear in standard annotation files, and the protocol described here can 
discover differential expression affecting these genes.

stringtie [-p <cpus>] [-v] [-G <guide_gff>] [-l <prefix>] <in.bam ..> [-o <out.gtf>] [--rf] or [--fr]

```bash
stringtie -v -p 4 -G ref_genome/annotation_Apis.gff3 -l w1 hisat2/w1.bam -o stringtie/w1.gtf
stringtie -v -p 4 -G ref_genome/annotation_Apis.gff3 -l w2 hisat2/w2.bam -o stringtie/w2.gtf

stringtie -v -p 4 -G ref_genome/annotation_Apis.gff3 -l aw1 hisat2/aw1.bam -o stringtie/aw1.gtf
stringtie -v -p 4 -G ref_genome/annotation_Apis.gff3 -l aw2 hisat2/aw2.bam -o stringtie/aw2.gtf

stringtie -v -p 4 -G ref_genome/annotation_Apis.gff3 -l q1 hisat2/q1.bam -o stringtie/q1.gtf
stringtie -v -p 4 -G ref_genome/annotation_Apis.gff3 -l q2 hisat2/q2.bam -o stringtie/q2.gtf
```


## Merging the gtf files

Merge
```bash
stringtie --merge -p 4 -G ref_genome/annotation_Apis.gff3 merge_List.txt -o stringtie/stringtie_merged.gtf
```


## GffCompare

GffCompare examine how the transcripts compare to the reference annotation:

If GffCompare was run with the -r option (i.e. comparing with a reference annotation), tracking rows will contain a "class code" value showing the relationship between a transfrag and the closest reference transcript (where applicable). If the -r option was not used the rows will all contain “-” in their class code column. The same codes are also shown as the value of the attribute "class_code" in the output GTF file. The class codes are shown below in decreasing order of their priority.

Here we are interested in a predicted transcript which falls entirely within a reference intron, yhis will be marked as (class code "i").

```bash
gffcompare -G -r ref_genome/annotation_Apis.gff3 stringtie/stringtie_merged.gtf -o gffCompare/gffComm
```
The –r option is followed by the annotation file to use as reference, and the –G option tells gffcompare to compare all transcripts in the input transcripts.gtf file, even those that might be redundant.


## Estimate transcript abundances and create table counts for Ballgown:

### Expression estimation mode (-e)

When the -e option is used, the reference annotation file -G is a required input and StringTie will not attempt to assemble 
the input read alignments but instead it will only estimate the expression levels of the "reference" transcripts provided in the -G file.

With this option, no "novel" transcript assemblies (isoforms) will be produced, and read alignments not overlapping any of the given 
reference transcripts will be ignored, which may provide a considerable speed boost when the given set of reference transcripts is limited 
to a set of target genes for example.

Running stringtie with No Expression estimation mode
```bash
stringtie -v -B -p 4 -G stringtie/stringtie_merged.gtf -o ballgown/SRR7908186/w1.gtf hisat2/w1.bam
stringtie -v -B -p 4 -G stringtie/stringtie_merged.gtf -o ballgown/SRR7908187/w2.gtf hisat2/w2.bam

stringtie -v -B -p 4 -G stringtie/stringtie_merged.gtf -o ballgown/SRR7908188/aw1.gtf hisat2/aw1.bam
stringtie -v -B -p 4 -G stringtie/stringtie_merged.gtf -o ballgown/SRR7908189/aw2.gtf hisat2/aw2.bam

stringtie -v -B -p 4 -G stringtie/stringtie_merged.gtf -o ballgown/SRR7908190/q1.gtf hisat2/q1.bam
stringtie -v -B -p 4 -G stringtie/stringtie_merged.gtf -o ballgown/SRR7908191/q2.gtf hisat2/q2.bam
```


## Extracting transcript sequences

gffread can also be used to generate a FASTA file with the DNA sequences for all transcripts in a GFF file. 
For this operation a fasta file with the genomic sequences have to be provided as well. For example, one might
want to extract the sequence of all transfrags (defined as transcripts or transcript fragments that result from 
the assembly process) assembled from a StringTie or Cufflinks assembly session. 

The file genome.fa in this example would be a multi-fasta file with the genomic sequences of the target genome. 
This also requires that every contig or chromosome name found in the 1st column of the input GFF file (transcript.gtf in this example)
must have a corresponding sequence entry in chromosomes.fa. 

This should be the case in our example if genome.fa is the file corresponding to the same genome (index) that was used for mapping the reads.

Note that the retrieval of the transcript sequences this way is going to be much faster if a fasta index file (genome.fa.fai in this example)
is found in the same directory with the genomic fasta file. 

Instalation
```bash
conda install -c bioconda gffread 
```

Creating a index with the samtools utility prior to running gffread:
```bash
samtools faidx genome.fa
```

Extracting the sequences into all_transcripts.fasta file 
```bash
gffread stringtie/stringtie_merged.gtf -g ref_genome/genome_Apis.fa -w blast/all_transcripts.fasta
```
39397 transcripts were found


## Blast

Creating a database with all transcript
```bash
makeblastdb -in all_transcripts.fasta -dbtype nucl -input_type fasta -title "all_transcripts" -out db/all_transcripts_db
```
db info:  
39,397 sequences;  
123,754,570 total letters  

## Blast of lncRNA against all_transcripts db
```bash
blastn -db all_transcripts_db -query lncRNA.fna -out res_blast_lncRNA.txt
```
    Query= NC_037648.1:c11949345-11920888 lncRNA [organism=Apis mellifera] [GeneID=726407] [chromosome=LG11]
    Length=28458
                                                                        Score     E
    Sequences producing significant alignments:                         (Bits)  Value
    MSTRG.8421.2              (new)                                      2588    0.0  
    MSTRG.8421.3               (new)                                     1690    0.0  
    transcript:XM_001120691    (described)                               1681    0.0  
    transcript:XM_006566320                                              826     0.0  
    transcript:XM_006566319                                              826     0.0  
    transcript:XM_393800                                                 826     0.0 


## References:

1 - Pertea, Mihaela et al. “Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown.” Nature protocols vol. 11,9 (2016): 1650-67. doi:10.1038/nprot.2016.095  
2 - [Hisat2](http://daehwankimlab.github.io/hisat2/)  
3 - [Stringtie](http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)  
4 - [GFF Utilities: GffRead and GffCompare.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7222033/)  
5 - [BLAST® Command Line Applications User Manual](https://www.ncbi.nlm.nih.gov/books/NBK279690/)  

