# Meta_Microbial-Workshop-2024
Official repository for the Meta_Microbial workshop which aims at providing an overview of current tools to analyze next-generation sequence data and to make insightful interpretation towards the fields of microbial ecology and environmental microbiology. 
This repository provides material and instruction accordingly to the lessons provided during the workshop and the sections belows are organized to mimic the day by day timeline.
Access to a remote machine will be available to ensure the correct 

# Ubuntu Windows Subsystem

# Access to the virtual machine - Instructions
- [ ] step 1
- [ ] step 2
- [ ] step n

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
