from gwf import Workflow, AnonymousTarget
from workflow import *
from os import listdir
from os.path import isfile, join
import pandas as pd

gwf = Workflow()

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

    for chr in chr_lst:

        gwf.target_from_template(
            f"{study_name}_{chr}_combine_gVCFs", 
            combine_gvcfs(
                infiles = [f"../steps/09_call_variants/{individual}_{chr}.g.vcf.gz" for individual in list(dct.keys())],
                chr = chr,
                workpath = f'../steps/10_combine_gvcfs/{study_name}/{chr}/',
                done = f"../done/10_combine_gvcfs/done_{study_name}_{chr}",
                done_prev = [f"../done/09_call_variants/{individual}_{chr}" for individual in list(dct.keys())]
                )
            )
    #     gwf.target_from_template(
    #         f"{study_name}_{chr}_genotype_gVCFs", 
    #         genotype_gvcfs(
    #             infile =  f"gendb://../steps/10_combine_gvcfs/{study_name}/{chr}/",
    #             ref = "../data/reference_data/newref/GCF_008728515.1_Panubis1.0_genomic.fna",
    #             outfile = f"../steps/11_final_vcf/{study_name}/{chr}.combined.vcf.gz",
    #             done = f"../done/11_final_vcf/done_{study_name}_{chr}",
    #             done_prev = f"../done/10_combine_gvcfs/done_{study_name}_{chr}"
    #             )
    #         )
    
    # gwf.target_from_template(
    #         f"{study_name}_concatenate_VCFs", 
    #         concatenate_vcfs(
    #             infiles =  [f"../steps/11_final_vcf/{study_name}/{chr}.combined.vcf.gz" for chr in chr_lst],
    #             outfile = f"../steps/12_concatenated_vcf/{study_name}.vcf.gz",
    #             done = f"../done/12_concatenated_vcf/done_{study_name}",
    #             done_prev = [f"../done/11_final_vcf/done_{study_name}_{chr}" for chr in chr_lst]
    #             )
    #         )


#### Add Robinson founder samples to Kuderna GenomicsDB

study_name_new = "robinson_2019"
individual_list = ["SAMN10524564", "SAMN10524546", "SAMN10524552", "SAMN10524565", "SAMN10524554", "SAMN10524539", "SAMN10524543", "SAMN10524555", "SAMN10524547", "SAMN10524550", "SAMN10524563", "SAMN10524560", "SAMN10524551", "SAMN10524559", "SAMN11119509", "SAMN10524549", "SAMN10524562", "SAMN10524548", "SAMN11119508", "SAMN10524566", "SAMN11119507", "SAMN10524544", "SAMN09761236", "SAMN10524567", "SAMN10524561", "SAMN10524540", "SAMN10524545", "SAMN10524541", "SAMN10524556", "SAMN10524558", "SAMN10524542"]
path_to_samples = f"../../../../shared_data/sequencing_data/baboon_other_studies/{study_name_new}"
dct = {}
acc = 0
with open(f"../metadata/sample_names_{study_name_new}.tsv") as file:
    for line in file:
        if acc == 0: 
            acc = 1
            continue
        line = line.rstrip().split("\t")
        if (study_name_new == "robinson_2019") and (line[1] not in individual_list):
            continue
        if line[1] not in dct: 
            dct[line[1]] = [line[3]]
        else: 
            dct[line[1]].append(line[3])


study_name_old = "kuderna_2023"
for chr in chr_lst:
    gwf.target_from_template(
        f"{study_name_new}_to_{study_name_old}_{chr}_update_genomicsdb", 
        update_genomicsdb(
            infiles = [f"../steps/09_call_variants/{individual}_{chr}.g.vcf.gz" for individual in list(dct.keys())],
            workpath = f'../steps/10_combine_gvcfs/{study_name_old}/{chr}/',
            done = f"../done/10_combine_gvcfs/done_{study_name_new}_to_{study_name_old}_{chr}",
            done_prev = f"../done/10_combine_gvcfs/done_{study_name_old}_{chr}",
            )
        )
    gwf.target_from_template(
            f"{study_name_new}_and_{study_name_old}_{chr}_genotype_gVCFs", 
            genotype_gvcfs(
                infile =  f"gendb://../steps/10_combine_gvcfs/{study_name}/{chr}/",
                ref = "../data/reference_data/newref/GCF_008728515.1_Panubis1.0_genomic.fna",
                outfile = f"../steps/11_final_vcf/{study_name_new}_and_{study_name_old}/{chr}.combined.vcf.gz",
                done = f"../done/11_final_vcf/done_{study_name_new}_and_{study_name_old}_{chr}",
                done_prev = f"../done/10_combine_gvcfs/done_{study_name_new}_to_{study_name_old}_{chr}"
                )
            )
    # gwf.target_from_template(
    #         f"{study_name}_concatenate_VCFs", 
    #         concatenate_vcfs(
    #             infiles =  [f"../steps/11_final_vcf/{study_name_new}_and_{study_name_old}/{chr}.combined.vcf.gz" for chr in chr_lst],
    #             outfile = f"../steps/12_concatenated_vcf/{study_name_new}_and_{study_name_old}.vcf.gz",
    #             done = f"../done/12_concatenated_vcf/done_{study_name_new}_and_{study_name_old}",
    #             done_prev = [f"../done/11_final_vcf/done_{study_name_new}_and_{study_name_old}_{chr}" for chr in chr_lst]
    #             )
    #         )
    

study_name = "vilgalys_fogel_2022"

path_to_samples = f"../../../../shared_data/sequencing_data/baboon_other_studies/{study_name}"

# For the old samples

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

# For the newly generated

onlyfiles = []
filepaths = []
for filepath in ['/mnt/primevo/shared_data/sequencing_data/novaseq_baboons/flow_cell1/Tung_8475_230606B7', 
                 '/mnt/primevo/shared_data/sequencing_data/novaseq_baboons/flow_cell2/Tung_8475_230530A6', 
                 '/mnt/primevo/shared_data/sequencing_data/novaseq_baboons/flow_cell3/Tung_8475_230531A7']:
    tmp = [f for f in listdir(filepath) if isfile(join(filepath, f)) and 'fastq'  in f and 'R1' in f]
    onlyfiles += tmp
    filepaths += [filepath]*len(tmp)
# dct = {}
for i, run in enumerate(onlyfiles):
    ilist = run.split('_')
    if "_".join(ilist[:-4]) not in dct: dct["_".join(ilist[:-4])] = []
    dct["_".join(ilist[:-4])].append((ilist[-4]+'_'+ilist[-3], filepaths[i]))

# print(list(dct.keys()))

study_name_new = "amboseli_all"
for chr in chr_lst:
    gwf.target_from_template(
        f"{study_name_new}_{chr}_combine_gVCFs", 
        combine_gvcfs(
            infiles = [f"../steps/09_call_variants/{individual}_{chr}.g.vcf.gz" for individual in list(dct.keys())],
            chr = chr,
            workpath = f'../steps/10_combine_gvcfs/{study_name_new}/{chr}/',
            done = f"../done/10_combine_gvcfs/done_{study_name_new}_{chr}",
            done_prev = [f"../done/09_call_variants/{individual}_{chr}" for individual in list(dct.keys())]
            )
        )
    gwf.target_from_template(
        f"{study_name_new}_{chr}_genotype_gVCFs", 
        genotype_gvcfs(
            infile =  f"gendb://../steps/10_combine_gvcfs/{study_name_new}/{chr}/",
            ref = "../data/reference_data/newref/GCF_008728515.1_Panubis1.0_genomic.fna",
            outfile = f"../steps/11_final_vcf/{study_name_new}/{chr}.combined.vcf.gz",
            done = f"../done/11_final_vcf/done_{study_name_new}_{chr}",
            done_prev = f"../done/10_combine_gvcfs/done_{study_name_new}_{chr}"
            )
        )
gwf.target_from_template(
    f"{study_name_new}_concatenate_VCFs", 
    concatenate_vcfs(
        infiles =  [f"../steps/11_final_vcf/{study_name_new}/{chr}.combined.vcf.gz" for chr in chr_lst],
        outfile = f"../steps/12_concatenated_vcf/{study_name_new}.vcf.gz",
        done = f"../done/12_concatenated_vcf/done_{study_name_new}",
        done_prev = [f"../done/11_final_vcf/done_{study_name_new}_{chr}" for chr in chr_lst]
        )
    )