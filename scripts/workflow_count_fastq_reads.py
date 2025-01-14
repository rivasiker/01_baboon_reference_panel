from gwf import Workflow, AnonymousTarget
import os
import sys
from workflow import *

gwf = Workflow()

# conda activate baboon_genotyping


study_name_list = ["kuderna_2023", "robinson_2019", "vilgalys_fogel_2022", "rogers_2019", "snyder-mackler_2016", "wall_2016"]

for study_name in study_name_list:
    path_to_samples = f"/mnt/primevo/shared_data/sequencing_data/baboon_other_studies/{study_name}"
    for filename in os.listdir(path_to_samples):
        if filename.endswith('_1.fastq.gz'):
            sample_name = filename.split(".")[0]
            gwf.target_from_template(
                f"{sample_name}_count_fastq_gz_reads", 
                count_fastq_gz_reads(
                    infile = path_to_samples+"/"+filename,
                    outfile = f"../steps/17_count_fastq_gz_reads/{filename}.txt",
                    done = f"../done/17_count_fastq_gz_reads/{sample_name}"
                    )
                )

# Afterwards...
# sum=0; for file in *.txt; do [ -s "$file" ] && sum=$((sum + $(cat "$file"))); done; echo "The total sum is: $sum"