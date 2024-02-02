
# ['study_accession', 'sample_accession', 'experiment_accession', 'run_accession', 'tax_id', 'scientific_name', 'fastq_md5', 'fastq_ftp', 'sample_alias']
dct = {}
st = True
with open("../metadata/sample_names_kuderna_2023.tsv") as file:
    for line in file:
        if st:
            st = False
            continue
        line = line.rstrip().split("\t")
        if line[1] not in dct: dct[line[1]] = [line[3]]
        else: dct[line[1]].append(line[3])

sample_list = []
for i in dct.values():
    sample_list = sample_list + i

print(sample_list)