import os
import sys
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

def split_fasta_to_n_files(input_fasta, output_dir, num_parts=10):
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
    return num_parts  # 返回实际的 part 数用于后续 parallel

def run_blastn(query_fasta, db, out_file):
    cmd = [
        "/staging/fanyucai/ncbi-blast-2.16.0+/bin/blastn",
        "-db", db,
        "-query", query_fasta,
        "-out", out_file,
        "-outfmt", "6 qseqid sacc pident length mismatch qcovs evalue bitscore score sscinames stitle",
        "-perc_identity", "98",
        "-max_target_seqs", "5",
        "-evalue", "1e-10",
        "-num_threads", "5",
    ]
    subprocess.run(cmd, check=True)
    print(f"Run blast Done: {query_fasta}")

def main():
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print("Usage: python parallel_blast.py input.fasta nt_viruses output_dir [num_parts]")
        sys.exit(1)

    input_fasta = sys.argv[1]
    db = sys.argv[2]
    output_dir = sys.argv[3]
    num_parts = int(sys.argv[4]) if len(sys.argv) == 5 else 10

    split_fasta_to_n_files(input_fasta, output_dir, num_parts)

    queries = [f"{output_dir}/part_{i}.fasta" for i in range(num_parts)]
    outputs = [f"{output_dir}/part_{i}.blast.txt" for i in range(num_parts)]

    with ThreadPoolExecutor(max_workers=num_parts) as executor:
        futures = [
            executor.submit(run_blastn, queries[i], db, outputs[i])
            for i in range(num_parts)
        ]
        for future in as_completed(futures):
            pass

    merged_out = os.path.join(output_dir, "blast_all.txt")
    with open(merged_out, "w") as outfile:
        outfile.write("#Query\tAccesion\tPer.Ident\tAlignment_Length\tMismatch\tQuery_Coverage\tE_value\tMax_Score\tTotal_Score\tScientific_Name\tDescription\n")
        for f in outputs:
            with open(f) as infile:
                outfile.writelines(infile.readlines())
            subprocess.check_call(f'rm -rf {f}', shell=True)
    subprocess.check_call(f'rm -rf {output_dir}/part_*.fasta',shell=True)
    print(f"\nAll task finished: {merged_out}")

if __name__ == "__main__":
    main()
