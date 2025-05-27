import gzip

def run(fastq_file, max_reads=10):
    # 判断是否是压缩文件
    open_func = gzip.open if fastq_file.endswith(".gz") else open

    total_len = 0
    read_count = 0

    with open_func(fastq_file, "rt") as f:
        while read_count < max_reads:
            header = f.readline()
            if not header:
                break  # 文件结束
            seq = f.readline().strip()     # 序列
            f.readline()                   # 跳过 '+'
            f.readline()                   # 跳过质量值

            total_len += len(seq)
            read_count += 1

    if read_count == 0:
        return 0
    read_length=total_len / read_count
    if  read_length in(301,151,101,251,51,201):
        read_length=read_length-1
    return int(read_length)  # 返回平均长度
