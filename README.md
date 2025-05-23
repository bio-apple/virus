# Virus:An integrated data analysis pipeline for viruses(amplicon-based)

It can be used for various infectious disease, public health surveillance, 
and microbial research applications, including viral whole-genome sequencing, 
antimicrobial resistance marker analysis, bacterial and fungal identification, and more.

If you have any questions, feel free to email me.**Email:fanyucai3@gmail.com**

Keeping databases up to date: 

**pangolin verison:4.3.1 pangolin-data: 1.33**：export PATH=/opt/conda/envs/pangolin/bin:$PATH 

**nextclad version:3.13.2** nextclade

**Freyja version:1.5.3**.export PATH=/opt/conda/envs/Freyja/bin:$PATH

**Last Update:2025.05**


## Docker

<pre>docker pull fanyucai1/virus</pre> 

## Software

## Download the latest version of the Nextclade database
<pre>python3 modules/download_nextclade_db.py -d nextclade_db</pre> 

## Build reference
<pre>python3 modules/build_database.py -o ref/</pre>


