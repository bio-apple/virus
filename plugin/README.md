cd /staging/fanyucai/virus/ref/nt_viruses 
/staging/fanyucai/ncbi-blast-2.14.1+/bin/blastdbcmd -db nt_viruses -entry_batch /staging/fanyucai/virus/ref/VSP/accession.list -outfmt "%f" > /staging/fanyucai/virus/ref/VSP/VSP.fasta
/staging/fanyucai/ncbi-blast-2.14.1+/bin/
