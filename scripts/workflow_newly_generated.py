from gwf import Workflow, AnonymousTarget
import os
from workflow import *
from os import listdir
from os.path import isfile, join


gwf = Workflow()

onlyfiles = []
filepaths = []
for filepath in ['/mnt/primevo/shared_data/sequencing_data/novaseq_baboons/flow_cell1/Tung_8475_230606B7', 
                 '/mnt/primevo/shared_data/sequencing_data/novaseq_baboons/flow_cell2/Tung_8475_230530A6', 
                 '/mnt/primevo/shared_data/sequencing_data/novaseq_baboons/flow_cell3/Tung_8475_230531A7']:
    tmp = [f for f in listdir(filepath) if isfile(join(filepath, f)) and 'fastq'  in f and 'R1' in f]
    onlyfiles += tmp
    filepaths += [filepath]*len(tmp)


dct = {}
for i, run in enumerate(onlyfiles):
    ilist = run.split('_')
    if "_".join(ilist[:-4]) not in dct: dct["_".join(ilist[:-4])] = []
    dct["_".join(ilist[:-4])].append((ilist[-4]+'_'+ilist[-3], filepaths[i]))

"R1_001.fastq.gz"

# conda activate baboon_genotyping

chr_lst = [
    "NC_044976.1",
    "NC_044977.1", "NC_044978.1", "NC_044979.1", "NC_044980.1", 
    "NC_044981.1", "NC_044982.1", "NC_044983.1", "NC_044984.1", 
    "NC_044985.1", "NC_044986.1", "NC_044987.1", "NC_044988.1", 
    "NC_044989.1", "NC_044990.1", "NC_044991.1", "NC_044992.1", 
    "NC_044993.1", "NC_044994.1", "NC_044995.1", "NC_044996.1", 
    "NC_044997.1", "NC_020006.2"
]

study_name = "amboseli_newly_generated_2023"

for individual in list(dct.keys()):
    for sample_name, path_to_samples in dct[individual]:
        gwf.target_from_template(
            f"{individual.replace('-', '_')}_{sample_name}_FASTQ_to_uBAM", 
            FASTQ_to_uBAM(
                sample_name = sample_name,
                fastq1 = f'{path_to_samples}/{individual}_{sample_name}_R1_001.fastq.gz',
                fastq2 = f'{path_to_samples}/{individual}_{sample_name}_R2_001.fastq.gz',
                outfile = f"../steps/01_uBAM_files/{individual}_{sample_name}.bam",
                done = f"../done/01_uBAM_files/{individual}_{sample_name}",
                tmp = f"../tmp/01_uBAM_files/"
                )
            )
        gwf.target_from_template(
            f"{individual.replace('-', '_')}_{sample_name}_mark_adapters", 
            mark_adapters(
                infile = f"../steps/01_uBAM_files/{individual}_{sample_name}.bam",
                outfile = f"../steps/02_mark_adapters/{individual}_{sample_name}.bam",
                tmp = f"../tmp/02_mark_adapters/",
                done = f"../done/02_mark_adapters/{individual}_{sample_name}",
                done_prev = f"../done/01_uBAM_files/{individual}_{sample_name}"
                )
            )
        gwf.target_from_template(
            f"{individual.replace('-', '_')}_{sample_name}_map_reads", 
            map_reads(
                infile = f"../steps/02_mark_adapters/{individual}_{sample_name}.bam",
                ref = "../data/reference_data/newref/GCF_008728515.1_Panubis1.0_genomic.fna",
                outfile = f"../steps/03_mapped_reads/{individual}_{sample_name}.mapped.bam",
                tmp = f"../tmp/03_mapped_reads/",
                done = f"../done/03_mapped_reads/{individual}_{sample_name}",
                done_prev = [f"../done/02_mark_adapters/{individual}_{sample_name}", '../done/00_make_ref_index_and_dict/done']
                )
            )

        gwf.target_from_template(
            f"{individual.replace('-', '_')}_{sample_name}_sort_coordinates", 
            sort_coordinates(
                infile = f"../steps/03_mapped_reads/{individual}_{sample_name}.mapped.bam",
                outfile = f"../steps/04_sort_coordinates/{individual}_{sample_name}.mapped.sorted.bam",
                temp = "../tmp/04_sort_coordinates/",
                done = f"../done/04_sort_coordinates/{individual}_{sample_name}",
                done_prev = f"../done/03_mapped_reads/{individual}_{sample_name}"
                )
            )
    gwf.target_from_template(
            f"{individual.replace('-', '_')}_merge_bams", 
            merge_bams(
                infile = [f"../steps/04_sort_coordinates/{individual}_{sample_name}.mapped.sorted.bam" for sample_name, path_to_samples in dct[individual]],
                outfile = f"../steps/05_merge_bams/{individual}.mapped.sorted.bam",
                done = f"../done/05_merge_bams/{individual}",
                done_prev = [f"../done/04_sort_coordinates/{individual}_{sample_name}" for sample_name, path_to_samples in dct[individual]]
                )
            )
    gwf.target_from_template(
        f"{individual.replace('-', '_')}_mark_and_remove_duplicates", 
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
        f"{individual.replace('-', '_')}_get_coverage", 
        get_coverage(
            infile = f"../steps/06_mark_and_remove_duplicates/{individual}.mapped.nodupes.bam",
            outfile = f"../steps/07_get_coverage/{individual}",
            done = f"../done/07_get_coverage/{individual}",
            done_prev = f"../done/06_mark_and_remove_duplicates/{individual}"
            )
        )
    
    gwf.target_from_template(
        f"{individual.replace('-', '_')}_rename_merged_bam",
        rename_merged_bam(
            infile = f"../steps/06_mark_and_remove_duplicates/{individual}.mapped.nodupes.bam",
            outfile = f"../steps/08_rename_merged_bam/{individual}.mapped.nodupes.renamed.bam",
            done = f"../done/08_rename_merged_bam/{individual}",
            done_prev = f"../done/06_mark_and_remove_duplicates/{individual}"
            )
        )
    
for individual in list(dct.keys()):
    for chr in chr_lst:
        # if os.path.exists(f"../done/09_call_variants/{individual.replace('-', '_')}_{chr}"): continue
        gwf.target_from_template(
            f"{individual.replace('-', '_')}_{chr}_call_variants", 
            call_variants_per_chromosome(
                chr = chr,
                ref = "../data/reference_data/newref/GCF_008728515.1_Panubis1.0_genomic.fna",
                infile = f"../steps/08_rename_merged_bam/{individual}.mapped.nodupes.renamed.bam",
                outfile = f"../steps/09_call_variants/{individual}_{chr}.g.vcf.gz",
                done = f"../done/09_call_variants/{individual}_{chr}",
                done_prev = [f"../done/07_get_coverage/{individual}", f"../done/08_rename_merged_bam/{individual}"]
                )
            )
        
# for chr in chr_lst:
#     gwf.target_from_template(
#         f"{study_name}_{chr}_combine_gVCFs", 
#         combine_gvcfs(
#             infiles = [f"../steps/09_call_variants/{individual}_{chr}.g.vcf.gz" for individual in list(dct.keys())],
#             chr = chr,
#             workpath = f'../steps/10_combine_gvcfs/{study_name}/{chr}/',
#             done = f"../done/10_combine_gvcfs/done_{study_name}_{chr}",
#             done_prev = [f"../done/09_call_variants/{individual}_{chr}" for individual in list(dct.keys())]
#             )
#         )
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