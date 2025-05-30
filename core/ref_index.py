import os
import subprocess
import argparse

parse=argparse.ArgumentParser()
parse.add_argument("-o","--outdir",help="directory of output",required=True)
args=parse.parse_args()
args.outdir=os.path.abspath(args.outdir)

docker='virus:latest'
accession=['NC_004162','NC_045512','NC_001477','NC_001474','NC_001475','NC_002640','NC_063383','NC_039199','NC_001802']
name=['Chikungunya_virus','SARS-CoV-2','Dengue_virus_type_1','Dengue_virus_type_2','Dengue_virus_type_3','Dengue_virus_type_4','Monkeypox_virus','Human_Metapneumovirus','HIV-1']

# https://genomes.atcc.org/genomes/4f980dee15b2432f
# https://www.ncbi.nlm.nih.gov/nuccore/KX087101
# Zika virus strain ZIKV/Homo sapiens/PRI/PRVABC59/2015, complete genome
name.append("Zika_virus")
accession.append("KX087101")

#RSV-A Human respiratory syncytial virus B cRNA https://www.ncbi.nlm.nih.gov/nuccore/LR699737

name.append("RSV-A")
accession.append("LR699737")

#RSV-B Human respiratory syncytial virus B isolate hRSV/B/Australia/VIC-RCH056/2019, complete genome https://www.ncbi.nlm.nih.gov/nuccore/OP975389.1
name.append("RSV-B")
accession.append('OP975389.1')

#Human adenovirus type 7 #https://www.frontiersin.org/journals/virology/articles/10.3389/fviro.2024.1462907/full and VSPv2
name.append("Human_adenovirus_type_7")
accession.append('AC_000018')

#Human mastadenovirus B (HAdV-B)#https://www.frontiersin.org/journals/virology/articles/10.3389/fviro.2024.1462907/full and ChatGPT
name.append("Human_adenovirus_B1")
accession.append('NC_011203')

#Human adenovirus F
name.append("Human_adenovirus_F")
accession.append("NC_001454")

#https://www.who.int/news-room/fact-sheets/detail/influenza-(seasonal)
#Influenza A virus

##H1N1
## https://www.ncbi.nlm.nih.gov/nuccore/?term=A/Wisconsin/588/2019
#####HA(segment 4 MW626062) strain:A/Wisconsin/588/2019
#####NA(segment 6 MW626056) strain:A/Wisconsin/588/2019
#####PB2 (segment 1 NC_026438) strain:A/California/07/2009
#####PB1 (segment 2 NC_026435) strain:A/California/07/2009
#####PA  (segment 3 NC_026437) strain:A/California/07/2009
#####NP  (segment 5 NC_026436) strain:A/California/07/2009
#####MP  (segment 7 NC_026431) strain:A/California/07/2009
#####NS  (segment 8 NC_026432) strain:A/California/07/2009
name.append("H1N1")
ref=['MW626062','MW626056','NC_026438','NC_026435','NC_026437','NC_026436','NC_026431','NC_026432']
accession.append(",".join(ref))

##H3N2
#####HA EPI1857216 A/Darwin/6/2021(maybe NCBI:OQ718999)
#####NA EPI1857215 A/Darwin/6/2021(maybe NCBI:OQ718998)
#####PB1 (segment 2) NC_007372 A/New York/392/2004
#####NP (segment 5) NC_007369 A/New York/392/2004
#####NS (segment 8) NC_007370 A/New York/392/2004
#####MP (segment 7) NC_007367 A/New York/392/2004
#####PA (segment 3) NC_007371 A/New York/392/2004
#####PB2 (segment 1) NC_007373 A/New York/392/2004
name.append("H3N2")
ref=['OQ718999','OQ718998','NC_007372','NC_007369',' NC_007370','NC_007367','NC_007371','NC_007373']
accession.append(",".join(ref))

##H5N1 A/Goose/Guangdong/1/96(H5N1)
#####HA  AF144305
#####NA  AF144304
#####NP  AF144303
#####PA  AF144302
#####PB1 AF144301
#####PB2 AF144300
#####NS  AF144307
#####    AF144306
name.append("H5N1")
ref=['AF144300','AF144301','AF144302','AF144303',' AF144304','AF144305','AF144306','AF144307']
accession.append(",".join(ref))

## H7N9 A/Anhui/DEWH72-01/2013(H7N9)
name.append("H7N9")
ref=['CY181520','CY181519','CY181518','CY181513','CY181516','CY181515','CY181514','CY181517']
accession.append(",".join(ref))

#Influenza B viruses
#B/Yamagata or B/Victoria lineage
#B/Yamagata viruses have not been observed since 2020.
#####(Vic) HA KX058884 B/Brisbane/60/2008
#####(Vic) NA CY073894 B/Brisbane/60/2008
#####(Vic) PA (segment 3) CY115156 B/Brisbane/60/2008
#####(Vic) PB1 (segment 1) CY115157 B/Brisbane/60/2008
#####(Vic) NP (segment 5) CY115154 B/Brisbane/60/2008
#####(Vic) MP (segment 7) CY115152 B/Brisbane/60/2008
#####(Vic) PB2 (segment 2) CY115158 B/Brisbane/60/2008
#####(Vic) NS (segment 8) CY115155 B/Brisbane/60/2008
name.append("Influenza_B_viruses_Victoria")
ref=['KX058884','CY073894','CY115156','CY115157',' CY115154','CY115152','CY115158','CY115155']
accession.append(",".join(ref))

#H10N4
name.append('H10N4')
ref=['OP572383','OP572384','OP572385','OP572386','OP572387','OP572388','OP572389','OP572390']
accession.append(",".join(ref))

print(f"Currently supported species list:{name}")

for a,b in zip(accession,name):
    subprocess.check_call(f'mkdir -p {args.outdir}/{b}',shell=True)
    subprocess.check_call(f'docker run --rm -v {args.outdir}/{b}:/ref/ {docker} sh -c \'export PATH=/opt/conda/envs/kraken2/bin:/opt/conda/bin:$PATH && '
                          f'cd /ref/ && efetch -db nucleotide -id {a} -format fasta >{b}.fasta && '
                          f'bowtie2-build {b}.fasta {b}.fasta && '
                          f'samtools faidx {b}.fasta && '
                          f'bwa index -a bwtsw {b}.fasta\'',shell=True)
