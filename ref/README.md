# Download NCBI Entrez
<pre>
    wget https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh
    ./install-edirect.sh
</pre>

# Download reference and build bowtie index

**Chikungunya_virus**
<pre>efetch -db nucleotide -id NC_004162 -format fasta >NC_004162.fasta
bowtie2-build NC_004162.fasta NC_004162.fasta</pre>

**SARS-CoV-2**
<pre>efetch -db nucleotide -id NC_045512 -format fasta >NC_045512.fasta
bowtie2-build NC_045512.fasta NC_045512.fasta</pre>

**Dengue_virus_type_1**
<pre>efetch -db nucleotide -id NC_001477 -format fasta >NC_001477.fasta
bowtie2-build NC_001477.fasta NC_001477.fasta</pre>

**Dengue_virus_type_2**
<pre>efetch -db nucleotide -id NC_001474 -format fasta >NC_001474.fasta
bowtie2-build NC_001474.fasta NC_001474.fasta</pre>

**Dengue_virus_type_3**
<pre>efetch -db nucleotide -id NC_001475 -format fasta >NC_001475.fasta
bowtie2-build NC_001475.fasta NC_001475.fasta</pre>

**Dengue_virus_type_4**
<pre>efetch -db nucleotide -id NC_002640 -format fasta >NC_002640.fasta
bowtie2-build NC_002640.fasta NC_002640.fasta</pre>

**Monkeypox_virus**
<pre>efetch -db nucleotide -id MT903345 -format fasta >MT903345.fasta</pre>