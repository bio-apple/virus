# virus

An integrated data analysis pipeline for viruses(amplicon-based), with a particular focus on keeping databases up to date: **pangolin**, **nextclade**, and **Freyja**.
It can be used for various infectious disease, public health surveillance, 
and microbial research applications, including viral whole-genome sequencing, 
antimicrobial resistance marker analysis, bacterial and fungal identification, and more.

If you have any questions, feel free to email me.
Email:fanyucai3@gmail.com

**Last Update:2025.05.19**


## Docker

<pre>docker pull fanyucai1/virus</pre> 

## Software

1.  Freyja version:1.5.3
<pre> docker run --rm fanyucai1/virus sh -c 'export PATH=/opt/conda/envs/Freyja/bin:$PATH && /opt/conda/envs/Freyja/bin/fryja --help'</pre>

2.  nextclade version:3.13.2
<pre> docker run --rm fanyucai1/virus sh -c 'nextclade -V'</pre> 

3.  pangolin verison:4.3.1 pangolin-data: 1.33
<pre>
docker run --rm fanyucai1/virus sh -c 'export PATH=/opt/conda/envs/pangolin/bin:$PATH && pangolin --all-versions'
</pre>

## Download the latest version of the Nextclade database
<pre>python3 modules/download_nextclade_db.py -d nextclade_db</pre> 

## Build reference
<pre>python3 modules/build_database.py -o ref/</pre>


