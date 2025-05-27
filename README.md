## An integrated data analysis pipeline for viruses(amplicon-based)

![flow-chart](./virus.png)

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

    

### 2.Database

**2-1:nextclade**
<pre>python3 modules/download_nextclade_db.py -d nextclade_db</pre> 

**2-2:virus genome and index**
<pre>python3 modules/build_database.py -o ref/</pre>

**2-3:ncbi nt virus**
<pre>
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
dnf install perl-Archive-Tar
dnf install perl-JSON-PP
perl ncbi-blast-2.16.0+/bin/update_blastdb.pl –showall
截止2025-05 NCBI blast本地化可下载数据库：
    Betacoronavirus
    ITS_RefSeq_Fungi
    28S_fungal_sequences
    18S_fungal_sequences
    ITS_eukaryote_sequences
    LSU_eukaryote_rRNA
    LSU_prokaryote_rRNA
    16S_ribosomal_RNA
    SSU_eukaryote_rRNA
    env_nt
    env_nr
    human_genome
    landmark
    mito
    mouse_genome
    nr
    nt_euk
    nt
    nt_others
    nt_prok
    nt_viruses
    pataa
    patnt
    pdbaa
    pdbnt
    ref_euk_rep_genomes
    ref_prok_rep_genomes
    ref_viroids_rep_genomes
    ref_viruses_rep_genomes
    refseq_select_rna
    refseq_select_prot
    refseq_protein
    refseq_rna
    swissprot
    tsa_nr
    tsa_nt
    taxdb
    core_nt
perl ncbi-blast-2.16.0+/bin/update_blastdb.pl nt_viruses --decompress
</pre>




