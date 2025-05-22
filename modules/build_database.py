infile=open("virus.tsv","w")

infile.write("#IMAP_Name\tIMAP_Accession\tNextclade_Accession\n")
infile.write("chikungunya virus\tNC_004162\t\n")
infile.write("Dengue 1 virus\tKM204119\tNC_001477\n")
infile.write("Dengue 2 virus\t\tNC_001474\n")
infile.write("Dengue 3 virus\t\tNC_001475\n")
infile.write("Dengue 4 virus\t\tNC_002640\n")
infile.write("Zika virus\t\t")
infile.write("Human respiratory syncytial virus A(RSV-A)\t\tEPI_ISL_412866")
infile.write("Human respiratory syncytial virus B(RSV-B)\t\tOP975389")
infile.write("Monkeypox virus\tMT903345\tNC_063383")
infile.close()

