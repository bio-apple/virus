## An integrated data analysis pipeline for viruses(amplicon-based)

### 1.Docker Last Update:2025.05

<pre>docker pull fanyucai1/virus</pre> 

databases and software list: 

    pangolin verison:4.3.1 
    pangolin-data: 1.33 
    export PATH=/opt/conda/envs/pangolin/bin:$PATH

    nextclad version:3.13.2
    nextclade

    Freyja version:1.5.3
    export PATH=/opt/conda/envs/Freyja/bin:$PATH

    

### 2.Download Nextclade database
<pre>python3 modules/download_nextclade_db.py -d nextclade_db</pre> 

### 3.Build reference
<pre>python3 modules/build_database.py -o ref/</pre>




