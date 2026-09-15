import subprocess
import sys

def run(num):
    with open("download_ncbi_viruses.sh", "w") as infile:
        infile.write("#!/bin/bash\n\n")

        for i in range(num + 1):
            if i < 10:
                infile.write(
                    f"wget https://ftp.ncbi.nlm.nih.gov/blast/db/nt_viruses.0{i}.tar.gz &\n"
                )
                infile.write(
                    f"wget https://ftp.ncbi.nlm.nih.gov/blast/db/nt_viruses.0{i}.tar.gz.md5 &\n"
                )
            else:
                infile.write(
                    f"wget https://ftp.ncbi.nlm.nih.gov/blast/db/nt_viruses.{i}.tar.gz &\n"
                )
                infile.write(
                    f"wget https://ftp.ncbi.nlm.nih.gov/blast/db/nt_viruses.{i}.tar.gz.md5 &\n"
                )

        # 等待所有 wget 完成
        infile.write("\nwait\n")

    # 执行下载脚本（会一直等到 wait 结束）
    subprocess.check_call("bash download_ncbi_viruses.sh", shell=True)

    # 下载完成后再解压
    for i in range(num + 1):
        if i < 10:
            subprocess.check_call(f"tar -xvf nt_viruses.0{i}.tar.gz", shell=True)
        else:
            subprocess.check_call(f"tar -xvf nt_viruses.{i}.tar.gz", shell=True)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(
            f"Usage: python {sys.argv[0]} <num>\n"
            "number files: nt_viruses.<num>.tar.gz"
        )
    else:
        run(int(sys.argv[1]))







