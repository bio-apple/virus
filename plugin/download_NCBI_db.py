import subprocess
import sys

infile=open("download_ncbi_viruses.sh","w")

def run(num):
    for i in range (0,num):
        if i<10:
            infile.write(f"nohup wget https://ftp.ncbi.nlm.nih.gov/blast/db/nt_viruses.0{i}.tar.gz&\n"
                     f"nohup wget https://ftp.ncbi.nlm.nih.gov/blast/db/nt_viruses.0{i}.tar.gz.md5&\n")

        else:
            infile.write(f"nohup wget https://ftp.ncbi.nlm.nih.gov/blast/db/nt_viruses.{i}.tar.gz&\n"
                         f"nohup wget https://ftp.ncbi.nlm.nih.gov/blast/db/nt_viruses.{i}.tar.gz.md5&\n")
    infile.close()
    subprocess.check_call(f'bash download_ncbi_viruses.sh', shell=True)

    for i in range(0,num):
        if i<10:
            subprocess.check_call(f'tar -xvf nt_viruses.0{i}.tar.gz',shell=True)
        else:
            subprocess.check_call(f'tar -xvf nt_viruses.{i}.tar.gz', shell=True)


if __name__=="__main__":
    if len(sys.argv)<2:
        print(f"Usage: python {sys.argv[0]} <num>\n"
              f"number files:nt_viruses.<num>.tar.gz\n"
              f"https://ftp.ncbi.nlm.nih.gov/blast/db/")
    else:
        num=int(sys.argv[1])
        run(num)












