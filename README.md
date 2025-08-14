# An integrated data analysis pipeline for viruses

![flow-chart](./virus.png)

## Step1.Docker

<pre>docker pull fanyucai1/virus
docker tag fanyucai1/virus virus</pre>

## Step2.Prepare Database Last Update:2025.08.07

    mkdir -p /ref/

**2-1:update or downlaod nextclade**
<pre>
rm -rf /ref/nextclade_db/
python3 core/update_nextclade_db.py -d /ref/nextclade_db</pre> 

**2-2:virus genome and index**
<pre>
mkdir -p /ref/bowtie2/
python3 core/ref_index.py -o /ref/bowtie2/</pre>

The currently available list of reference genomes for viral species includes:

    Chikungunya_virus
    Dengue_virus_type_1
    Dengue_virus_type_2
    Dengue_virus_type_3
    Dengue_virus_type_4
    H10N4
    H1N1
    H3N2
    H5N1
    H5N6
    H5N8
    H6N1
    H6N2
    H7N9
    H9N2
    HIV-1
    Human_Metapneumovirus
    Human_adenovirus_B1
    Human_adenovirus_F
    Human_adenovirus_type_7
    Influenza_B_viruses_Victoria
    Marburg_Virus
    Measles_virus
    Monkeypox_virus
    Porcine_reproductive_and_respiratory_syndrome_virus_1
    RSV-A
    RSV-B
    SARS-CoV-2
    Yellow_fever_virus
    Zika_virus

**2-3:ncbi nt virus**

*Download NCBI database using BLAST*:https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
<pre>
dnf install perl-Archive-Tar
dnf install perl-JSON-PP
</pre>

*View currently available databases for download:*
<pre>perl ncbi-blast-2.16.0+/bin/update_blastdb.pl --showall

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
</pre>

Download **nt_viruses**:
<pre>perl ncbi-blast-2.16.0+/bin/update_blastdb.pl nt_viruses --decompress</pre>

**If above method fails, you can directly download the corresponding database files from the NCBI BLAST database using wget.** (https://ftp.ncbi.nlm.nih.gov/blast/db/)

**2-4:vsp database**
<pre>
mkdir -p /ref/VSP/
cd /ref/VSP/
ncbi-blast-2.14.1+/bin/blastdbcmd -db nt_viruses -entry_batch /ref/VSP/accession.list -outfmt "%f" > /ref/VSP/VSP.fasta
ncbi-blast-2.14.1+/bin/makeblastdb -dbtype nucl -in /ref/VSP/VSP.fasta
wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
tar xzvf taxdb.tar.gz
</pre>

**2-5:kraken2 database:https://benlangmead.github.io/aws-indexes/k2**
<pre>
mkdir -p /ref/kraken/
cd /ref/kraken/
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20250402.tar.gz
tar xvzf k2_pluspf_20250402.tar.gz
</pre>

**2-5:Download or build host(default:human) genome bowtie2:https://github.com/BenLangmead/bowtie-majref**
<pre>
mkdir -p /ref/host/human/
cd /ref/host/human/
wget https://genome-idx.s3.amazonaws.com/bt/grch38_1kgmaj_snvindels_bt2.zip
unzip grch38_1kgmaj_snvindels_bt2.zip
</pre>

## Step3:run pipeline
<pre>
usage: Virus NGS pipeline.
Email:fanyucai3@gmail.com
 [-h] -p1 PE1 [PE1 ...] [-p2 PE2 [PE2 ...]] -p PREFIX [PREFIX ...] -o OUTDIR -c CONFIG -l {50,75,100,150,200,250,300} [-e BED] [-r REF] [-b BOWTIE2]

options:
  -h, --help            show this help message and exit
  -p1 PE1 [PE1 ...], --pe1 PE1 [PE1 ...]
                        R1 fastq
  -p2 PE2 [PE2 ...], --pe2 PE2 [PE2 ...]
                        R2 fastq
  -p PREFIX [PREFIX ...], --prefix PREFIX [PREFIX ...]
                        prefix of output
  -o OUTDIR, --outdir OUTDIR
                        diretory of output
  -c CONFIG, --config CONFIG
                        config file
  -l {50,75,100,150,200,250,300}, --length {50,75,100,150,200,250,300}
                        read length

Reference/Bowtie2 Index Options:
  -e BED, --bed BED     bed file
  -r REF, --ref REF     ref fasta
  -b BOWTIE2, --bowtie2 BOWTIE2
                        directory contains reference bowtie2 index
</pre>

Command-line example: 

    python3 pipeline.py \
    -p1 test_data/sample_S7_L001_R1_001.fastq.gz \
    -p2 test_data/sample_S7_L001_R2_001.fastq.gz \
    -o outdir/sample -l 150 -c config.ini -p sample

**Relevant external resources:**

Bacterial and Viral Bioinformatics Resource Center (BV-BRC):https://www.bv-brc.org