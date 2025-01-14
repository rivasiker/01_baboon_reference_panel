from gwf import Workflow, AnonymousTarget
from workflow import *
from os import listdir
from os.path import isfile, join
import pandas as pd

gwf = Workflow()
interval = 1_000_000

# conda activate baboon_genotyping

        
study_name_list = ["robinson_2019", "kuderna_2023"]
chr_lst = [
    "NC_044976.1",
    "NC_044977.1", "NC_044978.1", "NC_044979.1", "NC_044980.1", 
    "NC_044981.1", "NC_044982.1", "NC_044983.1", "NC_044984.1", 
    "NC_044985.1", "NC_044986.1", "NC_044987.1", "NC_044988.1", 
    "NC_044989.1", "NC_044990.1", "NC_044991.1", "NC_044992.1", 
    "NC_044993.1", "NC_044994.1", "NC_044995.1", 
    # "NC_044996.1", "NC_044997.1", "NC_020006.2"
]

chr_lst = [
    "NC_044995.1"
]

chr_length = pd.read_csv("./chr_length_Panubis1.0.tsv", sep = "\t")
chr_length = chr_length[["RefSeq seq accession", "Seq length"]]
chr_length = chr_length[chr_length["RefSeq seq accession"].isin(chr_lst)]
chr_length = pd.Series(chr_length["Seq length"].values, index=chr_length["RefSeq seq accession"]).to_dict()

# dct = {}
# for study_name in study_name_list:

#     individual_list = ["SAMN10524564", "SAMN10524546", "SAMN10524552", "SAMN10524565", "SAMN10524554", "SAMN10524539", "SAMN10524543", "SAMN10524555", "SAMN10524547", "SAMN10524550", "SAMN10524563", "SAMN10524560", "SAMN10524551", "SAMN10524559", "SAMN11119509", "SAMN10524549", "SAMN10524562", "SAMN10524548", "SAMN11119508", "SAMN10524566", "SAMN11119507", "SAMN10524544", "SAMN09761236", "SAMN10524567", "SAMN10524561", "SAMN10524540", "SAMN10524545", "SAMN10524541", "SAMN10524556", "SAMN10524558", "SAMN10524542"]
#     path_to_samples = f"../../../../shared_data/sequencing_data/baboon_other_studies/{study_name}"

#     acc = 0
#     with open(f"../metadata/sample_names_{study_name}.tsv") as file:
#         for line in file:
#             if acc == 0:
#                 acc = 1
#                 continue
#             line = line.rstrip().split("\t")
#             if (study_name == "robinson_2019") and (line[1] not in individual_list):
#                 continue
#             if line[1] not in dct:
#                 dct[line[1]] = [line[3]]
#             else:
#                 dct[line[1]].append(line[3])


# study_name = "robinson_founders_and_kuderna"
# for chr in chr_lst:
#     input_lst = []
#     done_lst = []
#     for start in range(1, chr_length[chr]+1, interval):
#         end = min([start+interval, chr_length[chr]])
#         gwf.target_from_template(
#             f"{study_name}_{chr}_{start}_{end}_combine_gVCFs", 
#             combine_gvcfs(
#                 infiles = [f"../steps/09_call_variants/{individual}_{chr}.g.vcf.gz" for individual in list(dct.keys())],
#                 chr = f"{chr}:{start}-{end}",
#                 workpath = f'../steps/10_combine_gvcfs_chunked/{study_name}/{chr}-{start}-{end}/',
#                 done = f"../done/10_combine_gvcfs_chunked/done_{study_name}_{chr}-{start}-{end}",
#                 done_prev = [f"../done/09_call_variants/{individual}_{chr}" for individual in list(dct.keys())]
#                 )
#             )
#         gwf.target_from_template(
#             f"{study_name}_{chr}_{start}_{end}_genotype_gVCFs", 
#             genotype_gvcfs(
#                 infile =  f"gendb://../steps/10_combine_gvcfs_chunked/{study_name}/{chr}-{start}-{end}/",
#                 ref = "../data/reference_data/newref/GCF_008728515.1_Panubis1.0_genomic.fna",
#                 outfile = f"../steps/11_final_vcf_chunked/{study_name}/{chr}-{start}-{end}.combined.vcf.gz",
#                 done = f"../done/11_final_vcf_chunked/done_{study_name}_{chr}-{start}-{end}",
#                 done_prev = f"../done/10_combine_gvcfs_chunked/done_{study_name}_{chr}-{start}-{end}"
#                 )
#             )
#         input_lst.append(f"../steps/11_final_vcf_chunked/{study_name}/{chr}-{start}-{end}.combined.vcf.gz")
#         done_lst.append(f"../done/11_final_vcf_chunked/done_{study_name}_{chr}-{start}-{end}")
#     gwf.target_from_template(
#         f"{study_name}_{chr}_concatenate_VCFs", 
#         concatenate_vcfs(
#             infiles =  input_lst,
#             outfile = f"../steps/12_concatenated_vcf_chunked/{study_name}.{chr}.vcf.gz",
#             done = f"../done/12_concatenated_vcf_chunked/done_{study_name}_{chr}",
#             done_prev = done_lst
#             )
#         )
    






individual_list_robinson_100 = [
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
# Excluded because of high duplication according to Gencove
individual_list_robinson_exclude = ["SAMN09761232", "SAMN09761233"]

study_name_list = ["robinson_2019", "kuderna_2023"]
dct = {}
for study_name in study_name_list:
    path_to_samples = f"../../../../shared_data/sequencing_data/baboon_other_studies/{study_name}"
    acc = 0
    with open(f"../metadata/sample_names_{study_name}.tsv") as file:
        for line in file:
            if acc == 0:
                acc = 1
                continue
            line = line.rstrip().split("\t")
            if (study_name == "robinson_2019") and (line[1] not in individual_list_robinson_100):
                continue
            if (study_name == "robinson_2019") and (line[1] in individual_list_robinson_exclude):
                continue
            if (study_name == "kuderna_2023") and (line[5] not in ["Papio anubis", "Papio cynocephalus"]):
                continue
            if line[1] not in dct:
                dct[line[1]] = [line[3]]
            else:
                dct[line[1]].append(line[3])

print(line[5])

print(list(dct.keys()))

study_name = "robinson_100_and_kuderna_anubis_yellow"
if not os.path.exists(f"../steps/10_combine_gvcfs_chunked/{study_name}/"):
    os.makedirs(f"../steps/10_combine_gvcfs_chunked/{study_name}/")
if not os.path.exists(f"../steps/11_final_vcf_chunked/{study_name}/"):
    os.makedirs(f"../steps/11_final_vcf_chunked/{study_name}/")
for chr in chr_lst:
    input_lst = []
    done_lst = []
    for start in range(1, chr_length[chr]+1, interval):
        end = min([start+interval, chr_length[chr]])
        gwf.target_from_template(
            f"{study_name}_{chr}_{start}_{end}_combine_gVCFs", 
            combine_gvcfs(
                infiles = [f"../steps/09_call_variants/{individual}_{chr}.g.vcf.gz" for individual in list(dct.keys())],
                chr = f"{chr}:{start}-{end}",
                workpath = f'../steps/10_combine_gvcfs_chunked/{study_name}/{chr}-{start}-{end}/',
                done = f"../done/10_combine_gvcfs_chunked/done_{study_name}_{chr}-{start}-{end}",
                done_prev = [f"../done/09_call_variants/{individual}_{chr}" for individual in list(dct.keys())]
                )
            )
        gwf.target_from_template(
            f"{study_name}_{chr}_{start}_{end}_genotype_gVCFs", 
            genotype_gvcfs(
                infile =  f"gendb://../steps/10_combine_gvcfs_chunked/{study_name}/{chr}-{start}-{end}/",
                ref = "../data/reference_data/newref/GCF_008728515.1_Panubis1.0_genomic.fna",
                outfile = f"../steps/11_final_vcf_chunked/{study_name}/{chr}-{start}-{end}.combined.vcf.gz",
                done = f"../done/11_final_vcf_chunked/done_{study_name}_{chr}-{start}-{end}",
                done_prev = f"../done/10_combine_gvcfs_chunked/done_{study_name}_{chr}-{start}-{end}"
                )
            )
        input_lst.append(f"../steps/11_final_vcf_chunked/{study_name}/{chr}-{start}-{end}.combined.vcf.gz")
        done_lst.append(f"../done/11_final_vcf_chunked/done_{study_name}_{chr}-{start}-{end}")
    gwf.target_from_template(
        f"{study_name}_{chr}_concatenate_VCFs", 
        concatenate_vcfs(
            infiles =  input_lst,
            outfile = f"../steps/12_concatenated_vcf_chunked/{study_name}.{chr}.vcf.gz",
            done = f"../done/12_concatenated_vcf_chunked/done_{study_name}_{chr}",
            done_prev = done_lst
            )
        )
