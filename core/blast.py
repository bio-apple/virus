import os
import sys
import subprocess,argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

docker="virus:latest"

def blastn(query_fasta, db, out_file):
    query_fasta=os.path.abspath(query_fasta)
    out_file=os.path.abspath(out_file)
    db_dir=os.path.abspath(os.path.dirname(db))
    db_name=db.split("/")[-1]
    if not os.path.exists(db_dir+"/taxdb.btd"):
        print("Blast database not contains taxdb.btd.Downlaod:ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz")
        exit()
    cmd=(f'docker run --rm -v {query_fasta}:/raw_data/{query_fasta.split("/")[-1]} '
         f'-v {os.path.dirname(out_file)}:/outdir/ '
         f'-v {db_dir}/:/ref/ ')

    cmd+=(f'{docker} '
          f'sh -c \"export PATH=/opt/conda/envs/kraken2/bin:$PATH && '
          f'export BLASTDB=/ref/ && '
          f'blastn -db /ref/{db_name} '
          f'-query /raw_data/{query_fasta.split("/")[-1]} '
          f'-out /outdir/{out_file.split("/")[-1]} '
          f'-outfmt \'6 qseqid sacc pident length mismatch qcovs evalue bitscore score sscinames stitle\' '
          f'-perc_identity 98 '
          f'-max_target_seqs 5 '
          f'-evalue 1e-10 '
          f'-num_threads 5\"')
    print(cmd)
    subprocess.check_call(cmd, shell=True)
    print(f"Run blast Done: {query_fasta}")

def run(input_fasta,db, output_dir,prefix,num_parts):
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: Count total number of sequences
    total_seqs = 0
    with open(input_fasta) as f:
        for line in f:
            if line.startswith(">"):
                total_seqs += 1

    # Step 2: Adjust num_parts if needed
    if total_seqs < num_parts:
        num_parts = total_seqs
        print(f"Only {total_seqs} sequences detected. Adjusting num_parts to {num_parts}")

    # Step 3: Create output file handles
    files = [open(f"{output_dir}/part_{i}.fasta", "w") for i in range(num_parts)]

    # Step 4: Distribute sequences in round-robin
    with open(input_fasta) as f:
        count = 0
        buffer = []
        for line in f:
            if line.startswith(">"):
                if buffer:
                    idx = count % num_parts
                    files[idx].writelines(buffer)
                    count += 1
                    buffer = []
            buffer.append(line)
        if buffer:
            idx = count % num_parts
            files[idx].writelines(buffer)

    for f in files:
        f.close()

    print(f"Split fasta done: {output_dir}/part_0.fasta ~ part_{num_parts - 1}.fasta")
    # 返回实际的 part 数用于后续 parallel
    subprocess.check_call(f"rm -rf {output_dir}/part_*.blast.txt", shell=True)
    queries = [f"{output_dir}/part_{i}.fasta" for i in range(num_parts)]
    outputs = [f"{output_dir}/part_{i}.blast.txt" for i in range(num_parts)]

    with ThreadPoolExecutor(max_workers=num_parts) as executor:
        futures = [
            executor.submit(blastn, queries[i], db, outputs[i])
            for i in range(num_parts)
        ]
        for future in as_completed(futures):
            pass

    merged_out = os.path.join(output_dir, f"{prefix}.blast_all.txt")
    with open(merged_out, "w") as outfile:
        outfile.write(
            "#Query\tAccesion\tPer.Ident\tAlignment_Length\tMismatch\tQuery_Coverage\tE_value\tMax_Score\tTotal_Score\tScientific_Name\tDescription\n")
        for f in outputs:
            with open(f) as infile:
                outfile.writelines(infile.readlines())
            subprocess.check_call(f'rm -rf {f}', shell=True)
    subprocess.check_call(f'rm -rf {output_dir}/part_*.fasta', shell=True)
    print(f"\nAll task finished: {merged_out}")

if __name__ == "__main__":
    parser=argparse.ArgumentParser("")
    parser.add_argument("-q",'--query',help="query fasta sequence",required=True)
    parser.add_argument("-d",'--db_name',help="blast database name",required=True)
    parser.add_argument("-o","--outdir",help="directory of output",default=os.getcwd())
    parser.add_argument("-p","--prefix",help="prefix of output",required=True)
    parser.add_argument("-n",'--num_parts',help="number split part",type=int,default=10)
    args=parser.parse_args()
    run(args.query,args.db_name,args.outdir,args.prefix,args.num_parts)