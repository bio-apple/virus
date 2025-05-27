import os
import sys
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

def split_fasta_to_10_files(input_fasta, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    files = [open(f"{output_dir}/part_{i}.fasta", "w") for i in range(10)]

    with open(input_fasta) as f:
        count = 0
        buffer = []
        for line in f:
            if line.startswith(">"):
                if buffer:
                    idx = count % 10
                    files[idx].writelines(buffer)
                    count += 1
                    buffer = []
            buffer.append(line)
        # 写入最后一条
        if buffer:
            idx = count % 10
            files[idx].writelines(buffer)
    for f in files:
        f.close()
    print(f"split fasta has done: {output_dir}/part_0.fasta ~ part_9.fasta")

def run_blastn(query_fasta, db, out_file):
    cmd = [
        "blastn",
        "-db", db,
        "-query", query_fasta,
        "-out", out_file,
        "-outfmt", "6",
        "-max_target_seqs", "10",
        "-evalue", "1e-10",
        "-num_threads", "2"
    ]
    subprocess.run(cmd, check=True)
    print(f"Done: {query_fasta}")

def main():
    if len(sys.argv) != 4:
        print("usage: python parallel_blast.py input.fasta nt_viruses output_dir")
        sys.exit(1)

    input_fasta = sys.argv[1]
    db = sys.argv[2]
    output_dir = sys.argv[3]

    split_fasta_to_10_files(input_fasta, output_dir)

    queries = [f"{output_dir}/part_{i}.fasta" for i in range(10)]
    outputs = [f"{output_dir}/part_{i}.blast.txt" for i in range(10)]

    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = [
            executor.submit(run_blastn, queries[i], db, outputs[i])
            for i in range(10)
        ]
        for future in as_completed(futures):
            pass  # 可选：处理异常

    # 合并所有输出
    merged_out = os.path.join(output_dir, "blast_all_merged.txt")
    with open(merged_out, "w") as outfile:
        for f in outputs:
            with open(f) as infile:
                outfile.writelines(infile.readlines())

    print(f"\nAll task finished：{merged_out}")

if __name__ == "__main__":
    main()
