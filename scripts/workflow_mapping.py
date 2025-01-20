from gwf import Workflow, AnonymousTarget
import os
from workflow import *

gwf = Workflow()

# conda activate baboon_genotyping


chr_lst = [
    "NC_044976.1",
    "NC_044977.1", "NC_044978.1", "NC_044979.1", "NC_044980.1", 
    "NC_044981.1", "NC_044982.1", "NC_044983.1", "NC_044984.1", 
    "NC_044985.1", "NC_044986.1", "NC_044987.1", "NC_044988.1", 
    "NC_044989.1", "NC_044990.1", "NC_044991.1", "NC_044992.1", 
    "NC_044993.1", "NC_044994.1", "NC_044995.1",
    "NC_044996.1", 
    "NC_044997.1", 
    "NC_020006.2"
]

# chr_lst = [
#     # "NC_044976.1",
#     # "NC_044977.1", "NC_044978.1", "NC_044979.1", "NC_044980.1", 
#     # "NC_044981.1", "NC_044982.1", "NC_044983.1", "NC_044984.1", 
#     # "NC_044985.1", "NC_044986.1", "NC_044987.1", "NC_044988.1", 
#     # "NC_044989.1", "NC_044990.1", "NC_044991.1", "NC_044992.1", 
#     # "NC_044993.1", "NC_044994.1", 
#     "NC_044995.1"
#     # "NC_044996.1", 
#     # "NC_044997.1", 
#     # "NC_020006.2"
# ]

# chr_lst = ["NC_044995.1"]

# chr_lst = ["NC_018169.2", "NC_018170.2"]
# sample_list = ['ERR10942715', 'ERR10942717', 'ERR10942719', 'ERR10942720', 'ERR10942721']
# path_to_samples = "../../../../shared_data/sequencing_data/baboon_other_studies/kuderna_2023/"

study_name_list = ["kuderna_2023", "robinson_2019", "vilgalys_fogel_2022", "snyder-mackler_2016", "wall_2016"]
study_name_list = ["kuderna_2023", "robinson_2019"]
# study_name_list = ["vilgalys_fogel_2022"]
# acc_lim = 500

for study_name in study_name_list:

    # individual_list = ["SAMN10524564", "SAMN10524546", "SAMN10524552", "SAMN10524565", "SAMN10524554", "SAMN10524539", "SAMN10524543", "SAMN10524555", "SAMN10524547", "SAMN10524550", "SAMN10524563", "SAMN10524560", "SAMN10524551", "SAMN10524559", "SAMN11119509", "SAMN10524549", "SAMN10524562", "SAMN10524548", "SAMN11119508", "SAMN10524566", "SAMN11119507", "SAMN10524544", "SAMN09761236", "SAMN10524567", "SAMN10524561", "SAMN10524540", "SAMN10524545", "SAMN10524541", "SAMN10524556", "SAMN10524558", "SAMN10524542"]
    # individual_list = ['SAMN11119506', 'SAMN09742522', 'SAMN09742523', 'SAMN09742524', 'SAMN09710567', 'SAMN11119513', 'SAMN09710569', 'SAMN08519680', 'SAMN09742526', 'SAMN09742527', 'SAMN09710571', 'SAMN09742528', 'SAMN09710572', 'SAMN09710573', 'SAMN09742529', 'SAMN09742530', 'SAMN09710574', 'SAMN09710575', 'SAMN09710576', 'SAMN09710577', 'SAMN09742531', 'SAMN09710578', 'SAMN09742532', 'SAMN09710579', 'SAMN09761229', 'SAMN12703524', 'SAMN09710580', 'SAMN09710581', 'SAMN09742533', 'SAMN09742534', 'SAMN09761230', 
    #                    'SAMN12703525', 'SAMN09761231', 'SAMN09761232', 'SAMN11119510', 'SAMN09742535', 'SAMN12703528', 'SAMN09761233', 'SAMN09742536', 'SAMN09761234', 'SAMN10524552', 'SAMN10524565', 'SAMN10524554', 'SAMN10524539', 'SAMN10524543', 'SAMN10524555', 'SAMN10524547', 'SAMN10524550', 'SAMN10524563', 'SAMN10524553', 'SAMN10524560', 'SAMN10524551', 'SAMN10524559', 'SAMN11119509', 'SAMN10524549', 'SAMN10524562', 'SAMN10524548', 'SAMN11119508', 'SAMN10524566', 'SAMN11119507', 'SAMN10524544', 'SAMN09761236', 
    #                    'SAMN10524567', 'SAMN10524561', 'SAMN10524540', 'SAMN10524545', 'SAMN09742537', 'SAMN10524541', 'SAMN09710582', 'SAMN10524556', 'SAMN10524558', 'SAMN10524542', 'SAMN09710583', 'SAMN09742538', 'SAMN09710584', 'SAMN09710585', 'SAMN09742539', 'SAMN09742540', 'SAMN10524557', 'SAMN09710586', 'SAMN09742541', 'SAMN10524564', 'SAMN10524546', 'SAMN09742542', 'SAMN09742543', 'SAMN12703534', 'SAMN11119511', 'SAMN09742544', 'SAMN11119512', 'SAMN09742545', 'SAMN10524536', 'SAMN09742546', 'SAMN10524537', 
    #                    'SAMN09742547', 'SAMN09710587', 'SAMN09710588', 'SAMN09742548', 'SAMN09742549', 'SAMN09710589', 'SAMN09710590', 'SAMN09710591', 'SAMN09742550', 'SAMN09742551', 'SAMN09742552', 'SAMN11119514', 'SAMN09742553']
    
    individual_list = [
        'SAMN11119506', 'SAMN09742522', 'SAMN09742523', 'SAMN09742524', 'SAMN09710567', 'SAMN11119513', 'SAMN09710569', 
        'SAMN08519680', 'SAMN09742526', 'SAMN09742527', 'SAMN09710571', 'SAMN09742528', 'SAMN09710572', 'SAMN09710573', 
        'SAMN09742529', 'SAMN09742530', 'SAMN09710574', 'SAMN09710575', 'SAMN09710576', 'SAMN09710577', 'SAMN09742531', 
        'SAMN09710578', 'SAMN09742532', 'SAMN09710579', 'SAMN09761229', 'SAMN09710580', 'SAMN09710581', 'SAMN09742533', 
        'SAMN09742534', 'SAMN09761230', 'SAMN09761231', 'SAMN09761232', 'SAMN11119510', 'SAMN09742535', 'SAMN09761233', 
        'SAMN09742536', 'SAMN09761234', 'SAMN10524552', 'SAMN10524565', 'SAMN10524554', 'SAMN10524539', 'SAMN10524543', 
        'SAMN10524555', 'SAMN10524547', 'SAMN10524550', 'SAMN10524563', 'SAMN10524553', 'SAMN10524560', 'SAMN10524551', 
        'SAMN10524559', 'SAMN11119509', 'SAMN10524549', 'SAMN10524562', 'SAMN10524548', 'SAMN11119508', 'SAMN10524566', 
        'SAMN11119507', 'SAMN10524544', 'SAMN09761236', 'SAMN10524567', 'SAMN10524561', 'SAMN10524540', 'SAMN10524545', 
        'SAMN09742537', 'SAMN10524541', 'SAMN09710582', 'SAMN10524556', 'SAMN10524558', 'SAMN10524542', 'SAMN09710583', 
        'SAMN09742538', 'SAMN09710584', 'SAMN09710585', 'SAMN09742539', 'SAMN09742540', 'SAMN10524557', 'SAMN09710586', 
        'SAMN09742541', 'SAMN10524564', 'SAMN10524546', 'SAMN09742542', 'SAMN09742543', 'SAMN11119511', 'SAMN09742544', 
        'SAMN11119512', 'SAMN09742545', 'SAMN09742546', 'SAMN09742547', 'SAMN09710587', 'SAMN09710588', 'SAMN09742548', 
        'SAMN09742549', 'SAMN09710589', 'SAMN09710590', 'SAMN09710591', 'SAMN09742550', 'SAMN09742551', 'SAMN09742552', 
        'SAMN11119514', 'SAMN09742553']
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
        # if os.path.exists(f"../done/16_count_reads/{individual}.mapped.nodupes.renamed.bam.mapped_reads.txt"): continue
        # acc += 1
        for sample_name in dct[individual]:
            # if not os.path.exists(f'{path_to_samples}/{sample_name}_1.fastq.gz'): continue
            # if not os.path.exists(f'{path_to_samples}/{sample_name}_2.fastq.gz'): continue
            # if all(os.path.exists(f"../done/08_call_variants/{sample_name}_{chr}") for chr in chr_lst): continue
            if (not os.path.exists(f'{path_to_samples}/{sample_name}_1.fastq.gz')) and (not os.path.exists(f'{path_to_samples}/{sample_name}_2.fastq.gz')):
                gwf.target_from_template(
                    f"{sample_name}_FASTQ_to_uBAM", 
                    FASTQ_to_uBAM_single_read(
                        sample_name = sample_name,
                        fastq1 = f'{path_to_samples}/{sample_name}.fastq.gz',
                        outfile = f"../steps/01_uBAM_files/{sample_name}.bam",
                        done = f"../done/01_uBAM_files/{sample_name}",
                        tmp = f"../tmp/01_uBAM_files/"
                        )
                    )
            else:
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
                    ref = "/mnt/primevo/work/iker_rivas_gonzalez/01_baboon_reference_panel/data/reference_data/newref/GCF_008728515.1_Panubis1.0_genomic.fna",
                    outfile = f"../steps/03_mapped_reads/{sample_name}.mapped.bam",
                    tmp = f"../tmp/03_mapped_reads/",
                    done = f"../done/03_mapped_reads/{sample_name}",
                    done_prev = [f"../done/02_mark_adapters/{sample_name}", '../done/00_make_ref_index_and_dict/done'],
                    min_seed_length = 50
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
            f"{individual}_filter_reads", 
            filter_reads_quality(
                infile = f"../steps/06_mark_and_remove_duplicates/{individual}.mapped.nodupes.bam",
                outfile = f"../steps/06_mark_and_remove_duplicates/{individual}.mapped.nodupes.filtered.bam",
                quality = 10,
                done = f"../done/06_mark_and_remove_duplicates/{individual}_filtered",
                done_prev = f"../done/06_mark_and_remove_duplicates/{individual}"
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
            f"{individual}_get_coverage_filtered", 
            get_coverage(
                infile = f"../steps/06_mark_and_remove_duplicates/{individual}.mapped.nodupes.filtered.bam",
                outfile = f"../steps/07_get_coverage/{individual}_filtered",
                done = f"../done/07_get_coverage/{individual}_filtered",
                done_prev = f"../done/06_mark_and_remove_duplicates/{individual}_filtered"
                )
            )

        gwf.target_from_template(
            f"{individual}_get_coverage_filtered_all", 
            get_coverage_all(
                infile = f"../steps/06_mark_and_remove_duplicates/{individual}.mapped.nodupes.filtered.bam",
                outfile = f"../steps/07_get_coverage/{individual}_filtered_all",
                done = f"../done/07_get_coverage/{individual}_filtered_all",
                done_prev = f"../done/06_mark_and_remove_duplicates/{individual}_filtered"
                )
            )
        
        gwf.target_from_template(
            f"{individual}_rename_merged_bam",
            rename_merged_bam(
                infile = f"../steps/06_mark_and_remove_duplicates/{individual}.mapped.nodupes.filtered.bam",
                outfile = f"../steps/08_rename_merged_bam/{individual}.mapped.nodupes.filtered.renamed.bam",
                done = f"../done/08_rename_merged_bam/{individual}_filtered",
                done_prev = f"../done/06_mark_and_remove_duplicates/{individual}"
                )
            )

        intermediate_bams = [
            f"../steps/05_merge_bams/{individual}.mapped.sorted.bam",
            f"../steps/06_mark_and_remove_duplicates/{individual}.mapped.nodupes.filtered.bam",
            # f"../steps/06_mark_and_remove_duplicates/{individual}.mapped.nodupes.bam",
            # f"../steps/08_rename_merged_bam/{individual}.mapped.nodupes.renamed.bam"
            ]
        for bam in intermediate_bams:
            gwf.target_from_template(
                f"{bam.split("/")[-1]}_count_all_reads",
                count_all_reads(
                    infile = bam,
                    outfile = f"../steps/16_count_reads/{bam.split("/")[-1]}.all_reads.txt",
                    done = f"../done/16_count_reads/{bam.split("/")[-1]}.all_reads",
                    done_prev = f"{"/".join(bam.split("/")[0:-1]).replace("steps", "done")}/{individual}"
                    )
                )
            gwf.target_from_template(
                f"{bam.split("/")[-1]}_count_mapped_reads",
                count_mapped_reads(
                    infile = bam,
                    outfile = f"../steps/16_count_reads/{bam.split("/")[-1]}.mapped_reads.txt",
                    done = f"../done/16_count_reads/{bam.split("/")[-1]}.mapped_reads",
                    done_prev = f"{"/".join(bam.split("/")[0:-1]).replace("steps", "done")}/{individual}"
                    )
                )
        
        # if acc == acc_lim: break



# acc_lim = 500
for study_name in study_name_list:

    # individual_list = ["SAMN10524564", "SAMN10524546", "SAMN10524552", "SAMN10524565", "SAMN10524554", "SAMN10524539", "SAMN10524543", "SAMN10524555", "SAMN10524547", "SAMN10524550", "SAMN10524563", "SAMN10524560", "SAMN10524551", "SAMN10524559", "SAMN11119509", "SAMN10524549", "SAMN10524562", "SAMN10524548", "SAMN11119508", "SAMN10524566", "SAMN11119507", "SAMN10524544", "SAMN09761236", "SAMN10524567", "SAMN10524561", "SAMN10524540", "SAMN10524545", "SAMN10524541", "SAMN10524556", "SAMN10524558", "SAMN10524542"]
    # individual_list = ['SAMN11119506', 'SAMN09742522', 'SAMN09742523', 'SAMN09742524', 'SAMN09710567', 'SAMN11119513', 'SAMN09710569', 'SAMN08519680', 'SAMN09742526', 'SAMN09742527', 'SAMN09710571', 'SAMN09742528', 'SAMN09710572', 'SAMN09710573', 'SAMN09742529', 'SAMN09742530', 'SAMN09710574', 'SAMN09710575', 'SAMN09710576', 'SAMN09710577', 'SAMN09742531', 'SAMN09710578', 'SAMN09742532', 'SAMN09710579', 'SAMN09761229', 'SAMN12703524', 'SAMN09710580', 'SAMN09710581', 'SAMN09742533', 'SAMN09742534', 'SAMN09761230', 
    #                    'SAMN12703525', 'SAMN09761231', 'SAMN09761232', 'SAMN11119510', 'SAMN09742535', 'SAMN12703528', 'SAMN09761233', 'SAMN09742536', 'SAMN09761234', 'SAMN10524552', 'SAMN10524565', 'SAMN10524554', 'SAMN10524539', 'SAMN10524543', 'SAMN10524555', 'SAMN10524547', 'SAMN10524550', 'SAMN10524563', 'SAMN10524553', 'SAMN10524560', 'SAMN10524551', 'SAMN10524559', 'SAMN11119509', 'SAMN10524549', 'SAMN10524562', 'SAMN10524548', 'SAMN11119508', 'SAMN10524566', 'SAMN11119507', 'SAMN10524544', 'SAMN09761236', 
    #                    'SAMN10524567', 'SAMN10524561', 'SAMN10524540', 'SAMN10524545', 'SAMN09742537', 'SAMN10524541', 'SAMN09710582', 'SAMN10524556', 'SAMN10524558', 'SAMN10524542', 'SAMN09710583', 'SAMN09742538', 'SAMN09710584', 'SAMN09710585', 'SAMN09742539', 'SAMN09742540', 'SAMN10524557', 'SAMN09710586', 'SAMN09742541', 'SAMN10524564', 'SAMN10524546', 'SAMN09742542', 'SAMN09742543', 'SAMN12703534', 'SAMN11119511', 'SAMN09742544', 'SAMN11119512', 'SAMN09742545', 'SAMN10524536', 'SAMN09742546', 'SAMN10524537', 
    #                    'SAMN09742547', 'SAMN09710587', 'SAMN09710588', 'SAMN09742548', 'SAMN09742549', 'SAMN09710589', 'SAMN09710590', 'SAMN09710591', 'SAMN09742550', 'SAMN09742551', 'SAMN09742552', 'SAMN11119514', 'SAMN09742553']    
    individual_list = [
        'SAMN11119506', 'SAMN09742522', 'SAMN09742523', 'SAMN09742524', 'SAMN09710567', 'SAMN11119513', 'SAMN09710569', 
        'SAMN08519680', 'SAMN09742526', 'SAMN09742527', 'SAMN09710571', 'SAMN09742528', 'SAMN09710572', 'SAMN09710573', 
        'SAMN09742529', 'SAMN09742530', 'SAMN09710574', 'SAMN09710575', 'SAMN09710576', 'SAMN09710577', 'SAMN09742531', 
        'SAMN09710578', 'SAMN09742532', 'SAMN09710579', 'SAMN09761229', 'SAMN09710580', 'SAMN09710581', 'SAMN09742533', 
        'SAMN09742534', 'SAMN09761230', 'SAMN09761231', 'SAMN09761232', 'SAMN11119510', 'SAMN09742535', 'SAMN09761233', 
        'SAMN09742536', 'SAMN09761234', 'SAMN10524552', 'SAMN10524565', 'SAMN10524554', 'SAMN10524539', 'SAMN10524543', 
        'SAMN10524555', 'SAMN10524547', 'SAMN10524550', 'SAMN10524563', 'SAMN10524553', 'SAMN10524560', 'SAMN10524551', 
        'SAMN10524559', 'SAMN11119509', 'SAMN10524549', 'SAMN10524562', 'SAMN10524548', 'SAMN11119508', 'SAMN10524566', 
        'SAMN11119507', 'SAMN10524544', 'SAMN09761236', 'SAMN10524567', 'SAMN10524561', 'SAMN10524540', 'SAMN10524545', 
        'SAMN09742537', 'SAMN10524541', 'SAMN09710582', 'SAMN10524556', 'SAMN10524558', 'SAMN10524542', 'SAMN09710583', 
        'SAMN09742538', 'SAMN09710584', 'SAMN09710585', 'SAMN09742539', 'SAMN09742540', 'SAMN10524557', 'SAMN09710586', 
        'SAMN09742541', 'SAMN10524564', 'SAMN10524546', 'SAMN09742542', 'SAMN09742543', 'SAMN11119511', 'SAMN09742544', 
        'SAMN11119512', 'SAMN09742545', 'SAMN09742546', 'SAMN09742547', 'SAMN09710587', 'SAMN09710588', 'SAMN09742548', 
        'SAMN09742549', 'SAMN09710589', 'SAMN09710590', 'SAMN09710591', 'SAMN09742550', 'SAMN09742551', 'SAMN09742552', 
        'SAMN11119514', 'SAMN09742553']
    
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
            if chr in ["NC_044996.1", "NC_044997.1", "NC_020006.2"]: 
                ploidy = "1" 
            else:
                ploidy = "2"
            gwf.target_from_template(
                f"{individual}_{chr}_call_variants", 
                call_variants_per_chromosome(
                    chr = chr,
                    ref = "/mnt/primevo/work/iker_rivas_gonzalez/01_baboon_reference_panel/data/reference_data/newref/GCF_008728515.1_Panubis1.0_genomic.fna",
                    infile = f"../steps/08_rename_merged_bam/{individual}.mapped.nodupes.filtered.renamed.bam",
                    outfile = f"../steps/09_call_variants/{individual}_{chr}.g.vcf.gz",
                    done = f"../done/09_call_variants/{individual}_{chr}",
                    done_prev = [f"../done/07_get_coverage/{individual}", f"../done/08_rename_merged_bam/{individual}_filtered"],
                    ploidy = ploidy
                    )
                )
        # if acc == acc_lim: break
            
