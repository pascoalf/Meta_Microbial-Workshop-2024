# Meta_Microbial-Workshop-2024
Official repository for the Meta_Microbial workshop which aims at providing an overview of current tools to analyze next-generation sequence data and to make insightful interpretation towards the fields of microbial ecology and environmental microbiology. 
This repository provides material and instruction accordingly to the lessons provided during the workshop and the sections belows are organized to mimic the day by day timeline.
Access to a remote machine will be available to ensure the correct 

> [!IMPORTANT]
> It is recommended that participants arrive with all the required programs installed properly. Carefully reading this GitHub page is also recommended.

# Jupyter Notebook & Python

Check if you already have python installed. You  can do so by typing `python` in your terminal.

Otherwise, install python from here:
+ [Windows](https://www.python.org/downloads/windows/)
+ [Mac](https://www.python.org/downloads/macos/)
+ [Linux](https://www.python.org/downloads/source/)

After installing python, `pip` should be available as command.

+ Run `pip install notebook` to install jupyter notebook.
+ Run `jupyter notebook` to open jupyter notebook.

An additional package required for the HandsOn in Metagenomics is `BioPython`.
To install it, run `pip install biopython`.

# Access to the virtual machine - Instructions

In order to allow the participants to run the analysis, a server has been provided with 96 threads, 640 GB of RAM, and 2 TB of storage.
To access the machine, please do:

+ download the "my_ecdsa_key" file 
+ access using `ssh -i my_ecdsa_key alunoNUMBER@34.71.198.208`

`NUMBER` and the ssh key will be provided during the workshop.

# Tmux
> tmux is a terminal multiplexer. It lets you switch easily between several programs in one terminal, detach them (they keep running in the background) and reattach them to a different terminal.

`tmux` will be used to ensure your analysis will keep running in case of connection issues or other kind of technical problems.

To open a new session
```
tmux new-session -t name
```

To detach from the open session, ensuring it keeps running in the background press `CTRL+B` and then `D`.

To attach a pre-existing session
```
tmux attach-session -t name
```

To kill a pre-existing session
```
tmux kill-session -t name
```

To get a list of the currently existing session
```
tmux ls
```

# First Day
## Next-generation sequencing technologies and data generation
## Sampling protocols and standardization
[Folder with lecture content](https://github.com/Nik3939/Meta_Microbial-Workshop-2024/tree/main/Sampling_Protocols_and_Standardization)
## R tools for microbial ecology studies
[Folder with lecture content](https://github.com/Nik3939/Meta_Microbial-Workshop-2024/tree/main/R-tools-for-microbial-ecology)

## Python in microbial ecology studies
[Folder with lecture content](https://github.com/Nik3939/Meta_Microbial-Workshop-2024/tree/main/Python-microbial-ecology)
## Metabarcoding (amplicon sequencing) analysis workflow
[Folder with lecture content]([/Metbarcoding-lectures](https://github.com/Nik3939/Meta_Microbial-Workshop-2024/tree/main/Metabarcoding-lectures))

# Second Day
## Exploring Online Resources and Repositories I
## Insights on data visualization and analysis
[Folder with lecture content]([/Data-visualization](https://github.com/Nik3939/Meta_Microbial-Workshop-2024/tree/main/Data-visualization))

## Metagenomics (shotgun sequencing) analysis workflow
## Phylogenetic trees from high-throughput sequence data
## Exploring Online Resources and Repositories II
## Alternative statistical and data analyses

# Third Day
## Hands-on Metabarcoding</summary>
### From raw data to community structure insights
[Folder with lecture content](Hands-on-Metabarcoding/Downstream)

## Hands-on Metagenomics</summary>
## From raw data to taxonomic and functional insights 

### 1st Quality Check

The first step of our analysis is to check the quality of the metagenomic sequences coming out from the sequencing platform.
Such data are present in the folder `raw data` and `fastQC` is the program we are going to use the assess their quality. 

Let's create a directory to store the output of the program.
```
mkdir fastqc_r
```

Run fastqc.
```
fastqc -o fastqc_r --threads 2 raw_data/M19-81_METAG_R1.fastq raw_data/M19-81_METAG_R2.fastq raw_data/M19-84_METAG_R1.fastq raw_data/M19-84_METAG_R2.fastq raw_data/M19-88_METAG_R1.fastq raw_data/M19-88_METAG_R2.fastq
```

Now a html report is present for each sample in the folder `fastqc_r`.

### Trimming

After assessing the quality of our sequences we need to remove sequences with low quality. For that, we are going to use `trimmomatic` to trim our sequences.
```
trimmomatic PE -threads 96 -phred33 -trimlog trim_log.log -summary trim_sum.log raw_data/M19-81_METAG_R1.fastq raw_data/M19-81_METAG_R2.fastq M19-81_f_p.fastq M19-81_f_u.fastq M19-81_r_p.fastq M19-81_r_u.fastq LEADING:10 TRAILING:10 SLIDINGWINDOW:5:20 

trimmomatic PE -threads 96 -phred33 -trimlog trim_log.log -summary trim_sum.log raw_data/M19-84_METAG_R1.fastq raw_data/M19-84_METAG_R2.fastq M19-84_f_p.fastq M19-84_f_u.fastq M19-84_r_p.fastq M19-84_r_u.fastq LEADING:10 TRAILING:10 SLIDINGWINDOW:5:20 

trimmomatic PE -threads 96 -phred33 -trimlog trim_log.log -summary trim_sum.log raw_data/M19-88_METAG_R1.fastq raw_data/M19-88_METAG_R2.fastq M19-88_f_p.fastq M19-88_f_u.fastq M19-88_r_p.fastq M19-88_r_u.fastq LEADING:10 TRAILING:10 SLIDINGWINDOW:5:20 
```
### Assembly

Now that we completed the first quality step of our pipeline. It is time to assemble the reads in contigs using `megahit`.
```
megahit -t 2 -o mega -1 trimm/M19-81_f_p.fastq,trimm/M19-84_f_p.fastq,trimm/M19-88_f_p.fastq -2 trimm/M19-81_r_p.fastq,trimm/M19-84_r_p.fastq,trimm/M19-88_r_p.fastq
```
An output folder called `mega` will automatically created and it will store the results.

### 2nd Quality Check

Our second quality check in our pipelines it assesing the quality of the assembly using `quast` or, to be more precise, the version for metagenomics `metaquast`.
Let's run the command.
```
metaquast -t 96 -o quast mega/final.contigs.fa
```

### Binning

After ensuring the good quality of the assembly, we now need to prepare for the binning step. In order to do so, we need to map raw reads to our obtained assembly to use the depth information for grouping contigs in bins.
To do so, we first need to index the assembly using `bowtie2-build`. 
Let's run the command.
```
bowtie2-build --threads 2 -f ../mega/final.contigs.fa bw
```

Then we query the indexed assembly with the raw reads.
Let's run the command.
```
bowtie2 -x bw -1 ../raw_data/M19-81_METAG_R1.fastq -2 ../raw_data/M19-81_METAG_R2.fastq -q --phred33  --threads 2 > M19-81.sam 2> M19-81.log

bowtie2 -x bw -1 ../raw_data/M19-84_METAG_R1.fastq -2 ../raw_data/M19-84_METAG_R2.fastq -q --phred33  --threads 2 > M19-84.sam 2> M19-84.log 

bowtie2 -x bw -1 ../raw_data/M19-88_METAG_R1.fastq -2 ../raw_data/M19-88_METAG_R2.fastq -q --phred33  --threads 2 > M19-88.sam 2> M19-88.log
```

Now, we need to convert the `SAM` output to a sorted and binary file `BAM`. 
We do this by running `samtools`.
```
samtools sort M19-81.sam -o M19-81.bam && samtools sort M19-84.sam -o M19-84.bam && samtools sort M19-88.sam -o M19-88.bam
```

Let's make some cleaning. We move the results from `bowtie2` and `samtools` in a directory called `bowtie`. 
We create the directory.
```
mkdir bowtie
```
We move all the `SAM` and `BAM` file to the directory
```
mv *.sam bowtie/.
mv *.bam bowtie/.
```

Finally, we run the binning software `metabat` using as input the assembly and the alignment.
```
runMetaBat.sh mega/final.contigs.fa bowtie/M19-81.bam bowtie/M19-84.bam bowtie/M19-88.bam
```

##### 3rd Quality Check

Now it is time for the last quality check, where we assess the completeness and contamination of our bins, together with other informative statistcs.
Let's run the following command.
```
checkm lineage_wf -t 2 -x fa metabat/final.contigs.fa.metabat-bins-20240819_221948 checkm
```

##### Taxonomy Assignment

After obtaining our bins and having assesed their quality with `checkm`, we will assign the taxonomy using `gtdbt-tk` by runnning the following commands.
```
conda activate gtdbtk-2.4.0

conda env config vars set GTDBTK_DATA_PATH="../../mnt/disk_2TB/release220";

export GTDBTK_DATA_PATH="../../mnt/disk_2TB/release220"

gtdbtk classify_wf -x fa --cpus 2 --genome_dir metabat/final.contigs.fa.metabat-bins-20240819_221948 --out_dir gtdbtk --skip_ani_screen
```

##### Download results
Output from `checkm`, `gtdbtk`, and `metabat` can be further analyzed in `jupyter notebook` to produce some visualization and retrieve the abundance of the members of the microbial community.
We first download the files locally on our laptop.
```
scp -i Documents/CIIMAR/LEC_METAG_MICROBIAL/my_ecdsa_key nicola@34.71.198.208:gtdbtk/gtdbtk.bac120.summary.tsv Downloads/.

scp -i Documents/CIIMAR/LEC_METAG_MICROBIAL/my_ecdsa_key nicola@34.71.198.208:checkm/storage/bin_stats.analyze.tsv Downloads/.

scp -i Documents/CIIMAR/LEC_METAG_MICROBIAL/my_ecdsa_key nicola@34.71.198.208:metabat/final.contigs.fa.depth.txt Downloads/.

scp -r -i Documents/CIIMAR/LEC_METAG_MICROBIAL/my_ecdsa_key nicola@34.71.198.208:metabat/final.contigs.fa.metabat-bins-20240819_221948 Downloads/.
```
Then, we run `jupyter notebook`.

##### Functional Annotation

After exploring the previous results, focused on the MAG reconstruction and taxonomy composition of the community. We can investigate the functional potential of our putative organisms.

First, we predict the genes present in our MAGs using `prodigal`.
We create a folder to store our results called `prodig` with the following command.
```
mkdir prodig
```

We then run `prodigal`.
```
prodigal -i metabat/final.contigs.fa.metabat-bins-20240819_221948/bin.1.fa -o prodig/bin.1.genes -a prodig/bin1.proteins.faa
prodigal -i metabat/final.contigs.fa.metabat-bins-20240819_221948/bin.6.fa -o prodig/bin.6.genes -a prodig/bin6.proteins.faa
prodigal -i metabat/final.contigs.fa.metabat-bins-20240819_221948/bin.12.fa -o prodig/bin.12.genes -a prodig/bin12.proteins.faa
```

After obtaining the aminoacidic sequences of our predicted genes, we can query the `KEGG` database using `kofamscan`.
Let's run the following command.
```
exec_annotation --profile=../../mnt/disk_2TB/kofam_scan-master/profiles --ko-list=../../mnt/disk_2TB/kofam_scan-master/ko_list --cpu=96 -o kofam/bin.1.ko.txt prodig/bin1.proteins.faa
exec_annotation --profile=../../mnt/disk_2TB/kofam_scan-master/profiles --ko-list=../../mnt/disk_2TB/kofam_scan-master/ko_list --cpu=96 -o kofam/bin.6.ko.txt prodig/bin6.proteins.faa
exec_annotation --profile=../../mnt/disk_2TB/kofam_scan-master/profiles --ko-list=../../mnt/disk_2TB/kofam_scan-master/ko_list --cpu=96 -o kofam/bin.12.ko.txt prodig/bin12.proteins.faa
```

Now, we want to select only the matches above the pre-computed threshold. We do that by using the powerful utilities provided by `bash`. Namely, `cat` and `awk`.
Let's run the following commands.
```
cat bin.6.ko.txt | t='*' awk '$1==ENVIRON["t"]{print $2,$3}' > bin.6.mapper.txt
cat bin.1.ko.txt | t='*' awk '$1==ENVIRON["t"]{print $2, $3}' > bin.1.mapper.txt
cat bin.12.ko.txt | t='*' awk '$1==ENVIRON["t"]{print $2, $3}' > bin.12.mapper.txt
```

Now we can use the online tool [KEGG Mapper Reconstruct](https://www.genome.jp/kegg/mapper/reconstruct.html) to map the genes to KEGG database.

##### Functional Annotation - Biosynthetic gene clusters
We will now use antiSMASH for the identification and annotation of biosynthetic gene clusters. antiSMASH can be used by a public [web version(https://antismash.secondarymetabolites.org/#!/start) or using a local installation, as we will perform here.
As an example, we will run antiSMASH in the recovered high-quality MAG bin.12. 
```
conda activate antismash
```
```
/mnt/disk_2TB/anaconda3/bin/./antismash metabat/final.contigs.fa.metabat-bins-20240819_221948/bin.12.fa --cb-general --cb-subclusters --cb-knownclusters --genefinding-tool prodigal --output-dir arctic_bin12_antismash_output &
```

After finishing we can download the antiSMASH output folder, which includes a detailed HTML output report and the genbank files of identified BGCs.
```
scp -i my_ecdsa_key -r adriana@34.71.198.208:/home/adriana/arctic_bin12_antismash_output/ .
```

We can complement the antiSMASH analysis by running [BiG-SCAPE](https://github.com/medema-group/BiG-SCAPE) (Biosynthetic Gene Similarity Clustering and Prospecting Engine), which builds sequence similarity networks (SSNs) of BGCs and groups them into Gene Cluster Families (GCFs). With this we can visualize how closely related the recovered BGCs are from each other and from BGCs that codify for characterized compounds in the [MIBiG] (https://mibig.secondarymetabolites.org/) database. 
To run BiG-SCAPE we use an input the output genbank files from the antiSMASH analysis. 
```
conda activate bigscape
```
```
cd /mnt/disk_2TB/BiG-SCAPE-1.1.5/
```
```
python bigscape.py -i  /home/adriana/arctic_bin12_antismash_output --mix --mibig --include_singletons --cutoffs 0.3 0.7 -o /home/adriana/arctic_bin12_bigscape_output
```

After finishing we can download the BiG-SCAPE output folder, which includes a detailed HTML output report and the network files. 
```
scp -i my_ecdsa_key -r adriana@34.71.198.208:/home/adriana/arctic_bin12_bigscape_output/ .
```

