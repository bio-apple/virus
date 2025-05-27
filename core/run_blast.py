import subprocess

def run(query_fasta, db, out_file):
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
