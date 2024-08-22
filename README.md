# Meta_Microbial-Workshop-2024
Official repository for the Meta_Microbial workshop which aims at providing an overview of current tools to analyze next-generation sequence data and to make insightful interpretation towards the fields of microbial ecology and environmental microbiology. 
This repository provides material and instruction accordingly to the lessons provided during the workshop and the sections belows are organized to mimic the day by day timeline.
Access to a remote machine will be available to ensure the correct 

# Jupyter Notebook & Python

Check if you already have python installed. You  can do so by typing `python` in your terminal.

Otherwise, install python from here:
+ [Windows](https://www.python.org/downloads/windows/)
+ [Mac](https://www.python.org/downloads/macos/)
+ [Linux](https://www.python.org/downloads/source/)

After installing python, `pip` should be available as command.

+ Run `pip install notebook` to install jupyter notebook.
+ Run `jupyter notebook` to open jupyter notebook.


# Access to the virtual machine - Instructions

In order to allow the participants to run the analysis, a server has been provided with 96 threads, 640 GB of RAM, and 2 TB of storage.
To access the machine, please do:

- [ ] download thhe "my_ecdsa_key" file 
- [ ] access using `ssh -i my_ecdsa_key alunoNUMBER@34.71.198.208`

Where `NUMBER` is a two-digit integer from 0 to 40.

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

# Local installation 
Necessary softwares:
+ Jupyter Notebook
+ Python 3
+ software n

> [!NOTE]
> Useful information that users should know, even when skimming content.

> [!TIP]
> It is .

> [!IMPORTANT]
> It is recommended that participants arrive with all the required programs installed properly. Carefully reading this GitHub page is also recommended.

> [!WARNING]
> Urgent info that needs immediate user attention to avoid problems.

> [!CAUTION]
> Advises about risks or negative outcomes of certain actions.

eventual footnotes [^1]

[^1]: reference

# First Day
## Next-generation sequencing technologies and data generation
## Sampling protocols and standardization
## R tools for microbial ecology studies
## Python in microbial ecology studies
## Metabarcoding (amplicon sequencing) analysis workflow

# Second Day
## Exploring Online Resources and Repositories I
## Insights on data visualization and analysis
## Metagenomics (shotgun sequencing) analysis workflow
## Phylogenetic trees from high-throughput sequence data
## Exploring Online Resources and Repositories II
## Alternative statistical and data analyses

# Third Day
## Hands-on Metabarcoding</summary>
### From raw data to community structure insights
## Hands-on Metagenomics</summary>
### From raw data to taxonomic and functional insights 

```
fastqc -o fastqc_r --threads 96 raw_data/M19-81_METAG_R1.fastq raw_data/M19-81_METAG_R2.fastq raw_data/M19-84_METAG_R1.fastq raw_data/M19-84_METAG_R2.fastq raw_data/M19-88_METAG_R1.fastq raw_data/M19-88_METAG_R2.fastq
```

```
trimmomatic PE -threads 96 -phred33 -trimlog trim_log.log -summary trim_sum.log raw_data/M19-81_METAG_R1.fastq raw_data/M19-81_METAG_R2.fastq M19-81_f_p.fastq M19-81_f_u.fastq M19-81_r_p.fastq M19-81_r_u.fastq LEADING:10 TRAILING:10 SLIDINGWINDOW:5:20 

trimmomatic PE -threads 96 -phred33 -trimlog trim_log.log -summary trim_sum.log raw_data/M19-84_METAG_R1.fastq raw_data/M19-84_METAG_R2.fastq M19-84_f_p.fastq M19-84_f_u.fastq M19-84_r_p.fastq M19-84_r_u.fastq LEADING:10 TRAILING:10 SLIDINGWINDOW:5:20 

trimmomatic PE -threads 96 -phred33 -trimlog trim_log.log -summary trim_sum.log raw_data/M19-88_METAG_R1.fastq raw_data/M19-88_METAG_R2.fastq M19-88_f_p.fastq M19-88_f_u.fastq M19-88_r_p.fastq M19-88_r_u.fastq LEADING:10 TRAILING:10 SLIDINGWINDOW:5:20 
```

```
megahit -t 96 -o mega -1 trimm/M19-81_f_p.fastq,trimm/M19-84_f_p.fastq,trimm/M19-88_f_p.fastq -2 trimm/M19-81_r_p.fastq,trimm/M19-84_r_p.fastq,trimm/M19-88_r_p.fastq
```

```
metaquast -t 96 -o quast mega/final.contigs.fa
```

```
bowtie2 -x bw -1 ../raw_data/M19-81_METAG_R1.fastq -2 ../raw_data/M19-81_METAG_R2.fastq -q --phred33  --threads 96 > M19-81.sam 2> M19-81.log

bowtie2 -x bw -1 ../raw_data/M19-84_METAG_R1.fastq -2 ../raw_data/M19-84_METAG_R2.fastq -q --phred33  --threads 96 > M19-84.sam 2> M19-84.log 

bowtie2 -x bw -1 ../raw_data/M19-88_METAG_R1.fastq -2 ../raw_data/M19-88_METAG_R2.fastq -q --phred33  --threads 96 > M19-88.sam 2> M19-88.log
```

```
samtools sort M19-81.sam -o M19-81.bam && samtools sort M19-84.sam -o M19-84.bam && samtools sort M19-88.sam -o M19-88.bam
```

```
runMetaBat.sh mega/final.contigs.fa bowtie/M19-81.bam bowtie/M19-84.bam bowtie/M19-88.bam
```

```
conda activate gtdbtk-2.4.0

conda env config vars set GTDBTK_DATA_PATH="../../mnt/disk_2TB/release220";

export GTDBTK_DATA_PATH="../../mnt/disk_2TB/release220"

gtdbtk classify_wf -x fa --cpus 96 --genome_dir metabat/final.contigs.fa.metabat-bins-20240819_221948 --out_dir gtdbtk --skip_ani_screen

```

```
scp -i Documents/CIIMAR/LEC_METAG_MICROBIAL/my_ecdsa_key nicola@34.71.198.208:gtdbtk/gtdbtk.bac120.summary.tsv Downloads/.

scp -i Documents/CIIMAR/LEC_METAG_MICROBIAL/my_ecdsa_key nicola@34.71.198.208:checkm/storage/bin_stats.analyze.tsv Downloads/.

scp -i Documents/CIIMAR/LEC_METAG_MICROBIAL/my_ecdsa_key nicola@34.71.198.208:metabat/final.contigs.fa.depth.txt Downloads/.

scp -r -i Documents/CIIMAR/LEC_METAG_MICROBIAL/my_ecdsa_key nicola@34.71.198.208:metabat/final.contigs.fa.metabat-bins-20240819_221948 Downloads/.
```
