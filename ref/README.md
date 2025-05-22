# Download NCBI Entrez
<pre>
    wget https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh
    ./install-edirect.sh
</pre>

mkdir /ref/


# Download reference and build bowtie index

<pre>

# Chikungunya_virus

mkdir -p /ref/bowtie2/Chikungunya_virus
cd /ref/bowtie2/Chikungunya_virus
efetch -db nucleotide -id NC_004162 -format fasta >NC_004162.fasta
bowtie2-build NC_004162.fasta NC_004162.fasta

# SARS-CoV-2

mkdir -p /ref/bowtie2/SARS-CoV-2
cd /ref/bowtie2/SARS-CoV-2
efetch -db nucleotide -id NC_045512 -format fasta >NC_045512.fasta
bowtie2-build NC_045512.fasta NC_045512.fasta

# Dengue_virus_type_1

mkdir -p /ref/bowtie2/Dengue_virus_type_1
efetch -db nucleotide -id NC_001477 -format fasta >NC_001477.fasta
bowtie2-build NC_001477.fasta NC_001477.fasta

# Dengue_virus_type_2
/ref/bowtie2/Dengue_virus_type_2
efetch -db nucleotide -id NC_001474 -format fasta >NC_001474.fasta
bowtie2-build NC_001474.fasta NC_001474.fasta

# Dengue_virus_type_3
mkdir -p /ref/bowtie2/Dengue_virus_type_3
efetch -db nucleotide -id NC_001475 -format fasta >NC_001475.fasta
bowtie2-build NC_001475.fasta NC_001475.fasta

# Dengue_virus_type_4
mkdir -p /ref/bowtie2/Dengue_virus_type_4
efetch -db nucleotide -id NC_002640 -format fasta >NC_002640.fasta
bowtie2-build NC_002640.fasta NC_002640.fasta

# Monkeypox_virus
mkdir -p mkdir -p /ref/bowtie2/Monkeypox_virus
efetch -db nucleotide -id NC_063383 -format fasta >NC_063383.fasta
bowtie2-build NC_063383.fasta NC_063383.fasta

</pre>
