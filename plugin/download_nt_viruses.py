import subprocess
infile=open("download_ncbi_viruses.sh","w")

for i in range (0,22):
    if i<10:
        infile.write(f"nohup wget https://ftp.ncbi.nlm.nih.gov/blast/db/nt_viruses.0{i}.tar.gz&\n"
                 f"nohup wget https://ftp.ncbi.nlm.nih.gov/blast/db/nt_viruses.0{i}.tar.gz.md5&\n")
    else:
        infile.write(f"nohup wget https://ftp.ncbi.nlm.nih.gov/blast/db/nt_viruses.{i}.tar.gz&\n"
                     f"nohup wget https://ftp.ncbi.nlm.nih.gov/blast/db/nt_viruses.{i}.tar.gz.md5&\n")

# Download taxdb contains:taxdb.bti、taxdb.btd、taxonomy4blast.sqlite3
infile.write("nohup wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz &\n"
             "nohup wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz.md5&\n")
infile.close()

subprocess.check_call(f'bash download_ncbi_viruses.sh', shell=True)

















