An integrated data analysis pipeline for viruses, with a particular focus on keeping databases up to date: **pangolin**, **nextclade**, and **Freyja**.

If you have any questions, feel free to email me.
Email:fanyucai3@gmail.com



**Last Update:2025.05.19**

# virus

## Docker

<pre>docker pull fanyucai1/virus</pre> 

## software

1.  Freyja version:1.5.3
<pre> docker run --rm fanyucai1/virus sh -c 'export PATH=/opt/conda/envs/Freyja/bin:$PATH && /opt/conda/envs/Freyja/bin/fryja --help'</pre>

2.  nextclade version:3.13.2
<pre> docker run --rm fanyucai1/virus sh -c 'nextclade -V'</pre> 

3.  pangolin verison:4.3.1 pangolin-data: 1.33
<pre>
docker run --rm fanyucai1/virus sh -c 'export PATH=/opt/conda/envs/pangolin/bin:$PATH && pangolin --all-versions'
</pre>

## nextclade_db

<pre>python3 sub/download_nextclade_db.py -d nextclade_db</pre> 




