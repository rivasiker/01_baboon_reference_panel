from gwf import Workflow, AnonymousTarget
import os
from workflow import *

gwf = Workflow()

# conda activate baboon_genotyping


# chr_lst = [
#     "NC_018152.2", "NC_018153.2", "NC_018154.2", "NC_018155.2", 
#     "NC_018156.2", "NC_018157.2", "NC_018158.2", "NC_018159.2", 
#     "NC_018160.2", "NC_018161.2", "NC_018162.2", "NC_018163.2", 
#     "NC_018164.2", "NC_018165.2", "NC_018166.2", "NC_018167.2", 
#     "NC_018168.2", "NC_018169.2", "NC_018170.2", "NC_018171.2", 
#     "NC_018172.2", "NC_020006.2"
# ]

chr_lst = [
    "NC_044976.1",
    "NC_044977.1", "NC_044978.1", "NC_044979.1", "NC_044980.1", 
    "NC_044981.1", "NC_044982.1", "NC_044983.1", "NC_044984.1", 
    "NC_044985.1", "NC_044986.1", "NC_044987.1", "NC_044988.1", 
    "NC_044989.1", "NC_044990.1", "NC_044991.1", "NC_044992.1", 
    "NC_044993.1", "NC_044994.1", "NC_044995.1", "NC_044996.1", 
    "NC_044997.1", "NC_020006.2"
]

# chr_lst = ["NC_018169.2", "NC_018170.2"]
# sample_list = ['ERR10942715', 'ERR10942717', 'ERR10942719', 'ERR10942720', 'ERR10942721']
# path_to_samples = "../../../../shared_data/sequencing_data/baboon_other_studies/kuderna_2023/"

study_name_list = ["kuderna_2023", "robinson_2019", "vilgalys_fogel_2022"]
study_name_list = ["vilgalys_fogel_2022"]
# acc_lim = 500

for study_name in study_name_list:

    individual_list = ["SAMN10524564", "SAMN10524546", "SAMN10524552", "SAMN10524565", "SAMN10524554", "SAMN10524539", "SAMN10524543", "SAMN10524555", "SAMN10524547", "SAMN10524550", "SAMN10524563", "SAMN10524560", "SAMN10524551", "SAMN10524559", "SAMN11119509", "SAMN10524549", "SAMN10524562", "SAMN10524548", "SAMN11119508", "SAMN10524566", "SAMN11119507", "SAMN10524544", "SAMN09761236", "SAMN10524567", "SAMN10524561", "SAMN10524540", "SAMN10524545", "SAMN10524541", "SAMN10524556", "SAMN10524558", "SAMN10524542"]
    path_to_samples = f"/mnt/primevo/shared_data/sequencing_data/baboon_other_studies/{study_name}"

    dct = {}
    acc = 0
    with open(f"../metadata/sample_names_{study_name}.tsv") as file:
        for line in file:
            if acc == 0: 
                acc = 1
                continue
            line = line.rstrip().split("\t")
            if (study_name == "robinson_2019") and (line[1] not in individual_list):
                continue
            if line[1] not in dct: 
                dct[line[1]] = [line[3]]
            else: 
                dct[line[1]].append(line[3])
            
    # print(len(dct))

    # with open(f'../../../../shared_data/sequencing_data/baboon_other_studies/md5_{study_name}.md5', 'r') as file:
    #     sample_list = sorted(list(set([f.rstrip().split(' ')[1].split('/')[1].split('_')[0] for f in file])))

    # n_sample = 300
    # acc = 0
    for individual in list(dct.keys()):
    # for individual in individual_list:
        # print(individual)
        # if acc == n_sample: break
        # if os.path.exists(f"../done/07_get_coverage/{individual}"): continue
        # acc += 1
        for sample_name in dct[individual]:
            # if not os.path.exists(f'{path_to_samples}/{sample_name}_1.fastq.gz'): continue
            # if not os.path.exists(f'{path_to_samples}/{sample_name}_2.fastq.gz'): continue
            # if all(os.path.exists(f"../done/08_call_variants/{sample_name}_{chr}") for chr in chr_lst): continue
            gwf.target_from_template(
                f"{sample_name}_FASTQ_to_uBAM", 
                FASTQ_to_uBAM(
                    sample_name = sample_name,
                    fastq1 = f'{path_to_samples}/{sample_name}_1.fastq.gz',
                    fastq2 = f'{path_to_samples}/{sample_name}_2.fastq.gz',
                    outfile = f"../steps/01_uBAM_files/{sample_name}.bam",
                    done = f"../done/01_uBAM_files/{sample_name}",
                    tmp = f"../tmp/01_uBAM_files/"
                    )
                )
            gwf.target_from_template(
                f"{sample_name}_mark_adapters", 
                mark_adapters(
                    infile = f"../steps/01_uBAM_files/{sample_name}.bam",
                    outfile = f"../steps/02_mark_adapters/{sample_name}.bam",
                    tmp = f"../tmp/02_mark_adapters/",
                    done = f"../done/02_mark_adapters/{sample_name}",
                    done_prev = f"../done/01_uBAM_files/{sample_name}"
                    )
                )
            gwf.target_from_template(
                f"{sample_name}_map_reads", 
                map_reads(
                    infile = f"../steps/02_mark_adapters/{sample_name}.bam",
                    ref = "../data/reference_data/newref/GCF_008728515.1_Panubis1.0_genomic.fna",
                    outfile = f"../steps/03_mapped_reads/{sample_name}.mapped.bam",
                    tmp = f"../tmp/03_mapped_reads/",
                    done = f"../done/03_mapped_reads/{sample_name}",
                    done_prev = [f"../done/02_mark_adapters/{sample_name}", '../done/00_make_ref_index_and_dict/done']
                    )
                )

            gwf.target_from_template(
                f"{sample_name}_sort_coordinates", 
                sort_coordinates(
                    infile = f"../steps/03_mapped_reads/{sample_name}.mapped.bam",
                    outfile = f"../steps/04_sort_coordinates/{sample_name}.mapped.sorted.bam",
                    temp = "../tmp/04_sort_coordinates/",
                    done = f"../done/04_sort_coordinates/{sample_name}",
                    done_prev = f"../done/03_mapped_reads/{sample_name}"
                    )
                )
            
        gwf.target_from_template(
            f"{individual}_merge_bams", 
            merge_bams(
                infile = [f"../steps/04_sort_coordinates/{sample_name}.mapped.sorted.bam" for sample_name in dct[individual]],
                outfile = f"../steps/05_merge_bams/{individual}.mapped.sorted.bam",
                done = f"../done/05_merge_bams/{individual}",
                done_prev = [f"../done/04_sort_coordinates/{sample_name}" for sample_name in dct[individual]]
                )
            )

        gwf.target_from_template(
            f"{individual}_mark_and_remove_duplicates", 
            mark_and_remove_duplicates(
                infile = f"../steps/05_merge_bams/{individual}.mapped.sorted.bam",
                metrics = f"../steps/06_mark_and_remove_duplicates/{individual}.mapped.nodupes.bam.metrics",
                outfile = f"../steps/06_mark_and_remove_duplicates/{individual}.mapped.nodupes.bam",
                temp = "../tmp/06_mark_and_remove_duplicates/",
                done = f"../done/06_mark_and_remove_duplicates/{individual}",
                done_prev = f"../done/05_merge_bams/{individual}"
                )
            )
        
        gwf.target_from_template(
            f"{individual}_get_coverage", 
            get_coverage(
                infile = f"../steps/06_mark_and_remove_duplicates/{individual}.mapped.nodupes.bam",
                outfile = f"../steps/07_get_coverage/{individual}",
                done = f"../done/07_get_coverage/{individual}",
                done_prev = f"../done/06_mark_and_remove_duplicates/{individual}"
                )
            )
        
        gwf.target_from_template(
            f"{individual}_rename_merged_bam",
            rename_merged_bam(
                infile = f"../steps/06_mark_and_remove_duplicates/{individual}.mapped.nodupes.bam",
                outfile = f"../steps/08_rename_merged_bam/{individual}.mapped.nodupes.renamed.bam",
                done = f"../done/08_rename_merged_bam/{individual}",
                done_prev = f"../done/06_mark_and_remove_duplicates/{individual}"
                )
            )
        
        # if acc == acc_lim: break



# acc_lim = 500
for study_name in study_name_list:

    individual_list = ["SAMN10524564", "SAMN10524546", "SAMN10524552", "SAMN10524565", "SAMN10524554", "SAMN10524539", "SAMN10524543", "SAMN10524555", "SAMN10524547", "SAMN10524550", "SAMN10524563", "SAMN10524560", "SAMN10524551", "SAMN10524559", "SAMN11119509", "SAMN10524549", "SAMN10524562", "SAMN10524548", "SAMN11119508", "SAMN10524566", "SAMN11119507", "SAMN10524544", "SAMN09761236", "SAMN10524567", "SAMN10524561", "SAMN10524540", "SAMN10524545", "SAMN10524541", "SAMN10524556", "SAMN10524558", "SAMN10524542"]
    path_to_samples = f"../../../../shared_data/sequencing_data/baboon_other_studies/{study_name}"
    dct = {}
    acc = 0
    with open(f"../metadata/sample_names_{study_name}.tsv") as file:
        for line in file:
            if acc == 0: 
                acc = 1
                continue
            line = line.rstrip().split("\t")
            if (study_name == "robinson_2019") and (line[1] not in individual_list):
                continue
            if line[1] not in dct: 
                dct[line[1]] = [line[3]]
            else: 
                dct[line[1]].append(line[3])
    # acc = 0
    for individual in list(dct.keys()):
        # if not os.path.exists(f"../done/08_rename_merged_bam/{individual}"): continue
        # acc += 1
        for chr in chr_lst:
            # if os.path.exists(f"../done/09_call_variants/{individual}_{chr}"): continue
            gwf.target_from_template(
                f"{individual}_{chr}_call_variants", 
                call_variants_per_chromosome(
                    chr = chr,
                    ref = "../data/reference_data/newref/GCF_008728515.1_Panubis1.0_genomic.fna",
                    infile = f"../steps/08_rename_merged_bam/{individual}.mapped.nodupes.renamed.bam",
                    outfile = f"../steps/09_call_variants/{individual}_{chr}.g.vcf.gz",
                    done = f"../done/09_call_variants/{individual}_{chr}",
                    done_prev = [f"../done/07_get_coverage/{individual}", f"../done/08_rename_merged_bam/{individual}"]
                    )
                )
        # if acc == acc_lim: break
            
