import os
import subprocess
import argparse

docker="virus:latest"

def run(R1,outdir,prefix,R2=None,max_sequences=4000000):

    R1=os.path.abspath(R1)
    if R1.endswith(".gz"):
        sequence_count = int(os.popen(f'zcat {R1}|wc -l').read().strip("\n").split(" ")[0]) // 4
    else:
        sequence_count =int(os.popen(f'wc -l {R1}').read().strip("\n").split(" ")[0])//4
    print(sequence_count)
    if sequence_count>max_sequences:
        os.makedirs(outdir, exist_ok=True)
        outdir = os.path.abspath(outdir)
        cmd=f"docker run --rm -v {outdir}:/outdir/ -v {R1}:/raw_data/{R1.split('/')[-1]} "
        cmd1 = cmd + f"{docker} sh -c \'export PATH=/opt/conda/bin:$PATH && seqtk sample -s100 /raw_data/{R1.split('/')[-1]} {max_sequences} | pigz > /outdir/{prefix}.sub.R1.fastq.gz\'"
        if R2 is not None:
            R2=os.path.abspath(R2)
            cmd2 = cmd + f"-v {R2}:/raw_data/{R2.split('/')[-1]} {docker} sh -c \'export PATH=/opt/conda/bin:$PATH &&  seqtk sample -s100 /raw_data/{R2.split('/')[-1]} {max_sequences} | pigz > /outdir/{prefix}.sub.R2.fastq.gz\'"
            process1 = subprocess.Popen(cmd1, shell=True)
            process2 = subprocess.Popen(cmd2, shell=True)
            process1.wait()
            process2.wait()
        else:
            subprocess.check_call(cmd1,shell=True)
    return sequence_count

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Subsample virus sequences')
    parser.add_argument('-p1',"--pe1", help='R1 fastq file',required=True)
    parser.add_argument('-p2',"--pe2", help='R2 fastq file',default=None)
    parser.add_argument('-p',"--prefix", help='prefix of output',required=True)
    parser.add_argument('-o','--outdir', help='directory of output',required=True)
    parser.add_argument('-m','--max', help='Subsample sequence number,default=4,000,000',type=int,default=4000000)
    args = parser.parse_args()
    run(args.pe1,args.outdir,args.prefix,args.pe2,args.max)