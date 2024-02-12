from gwf import Workflow, AnonymousTarget
import os

gwf = Workflow()

# conda activate baboon_genotyping


#################################################################################################################################
################################################# ---- DOWNLOAD REFERENCE ---- ##################################################
#################################################################################################################################


def download_reference(accession, name, outfile):
    """Download reference file."""
    inputs = []
    outputs = [outfile, f"../data/reference_data/newref/{accession}_{name}_genomic.fna"]
    options = {'cores': 1, 'memory': '20g', 'walltime': "01:00:00"}
    spec = f"""
        cd ../data/reference_data/newref/
        datasets download genome accession {accession} --include genome 
        unzip ncbi_dataset.zip
        mv ./ncbi_dataset/data/{accession}/{accession}_{name}_genomic.fna ./
        rm -r ncbi_dataset*
        rm README.md
        gzip -c {accession}_{name}_genomic.fna > {accession}_{name}_genomic.fna.gz
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# gwf.target_from_template(
#     f"download_reference", 
#     download_reference(
#         accession = 'GCF_000264685.3',
#         name = 'Panu_3.0',
#         outfile = "../data/reference_data/newref/GCF_000264685.3_Panu_3.0_genomic.fna.gz"
#         )
#     )

gwf.target_from_template(
    f"download_reference", 
    download_reference(
        accession = 'GCF_008728515.1',
        name = 'Panubis1.0',
        outfile = "../data/reference_data/newref/GCF_008728515.1_Panubis1.0_genomic.fna.gz"
        )
    )

#############################################################################################################################
################################################# ---- MAKE REF INDEX ---- ##################################################
#############################################################################################################################

# Check genome size per contig, sorted by chromosome size
# cut -f1,2 GCF_008728515.1_Panubis1.0_genomic.fna.fai | sort -r -n -k 2 | head -n 30
# Print all chromosomes starting with NC
# cut -f1,2 GCF_008728515.1_Panubis1.0_genomic.fna..fai |  grep NC_

def make_ref_index_and_dict(infile, done):
    """Make index and dictionary for reference file."""
    inputs = [infile]
    outputs = [done] + [infile.replace('.fna', '.dict')] + [infile+'.'+i for i in ['amb', 'ann', 'bwt', 'fai', 'pac', 'sa']]
    options = {'cores': 1, 'memory': '20g', 'walltime': "02:00:00"}
    spec = f"""
        /mnt/primevo/shared_data/software/gatk/gatk-4.3.0.0/gatk CreateSequenceDictionary -R {infile}
        bwa index -a bwtsw {infile}
        samtools faidx {infile}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

gwf.target_from_template(
    f"make_ref_index_and_dict", 
    make_ref_index_and_dict(
        infile = '../data/reference_data/newref/GCF_008728515.1_Panubis1.0_genomic.fna',
        done = f"../done/00_make_ref_index_and_dict/done"
        )
    )

########################################################################################################################
################################################# ---- MAKE uBAM ---- ##################################################
########################################################################################################################

def FASTQ_to_uBAM(sample_name, fastq1, fastq2, outfile, done):
    """Make uBAM file from FASTQ file"""
    inputs = [fastq1, fastq2]
    outputs = [done, outfile]
    options = {'cores': 1, 'memory': '20g', 'walltime': "01:00:00"}
    spec = f"""
        picard FastqToSam F1=$PWD/{fastq1} F2=$PWD/{fastq2} O=$PWD/{outfile} SAMPLE_NAME={sample_name}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

#####################################################################################################################################
################################################# ---- Mark Illumina adapters ---- ##################################################
#####################################################################################################################################

def mark_adapters(infile, outfile, tmp, done, done_prev):
    """Mark adapter sequences in uBAM file."""
    inputs = [infile, done_prev]
    outputs = [done, outfile, outfile+'_metrics']
    options = {'cores': 1, 'memory': '20g', 'walltime': "01:00:00"}
    spec = f"""
        picard MarkIlluminaAdapters I=$PWD/{infile} O=$PWD/{outfile} M=$PWD/{outfile+'_metrics'} TMP_DIR=$PWD/{tmp}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

###########################################################################################################################
################################################# ---- MAPPING PIPE ---- ##################################################
###########################################################################################################################


def map_reads(infile, ref, outfile, tmp, done, done_prev):
    """Mapping pipe."""
    inputs = done_prev + [ref, infile]
    outputs = [done, outfile]
    options = {'cores': 4, 'memory': '40g', 'walltime': "12:00:00"}
    spec = f"""
        set -o pipefail
        picard SamToFastq I=$PWD/{infile} FASTQ=/dev/stdout INTERLEAVE=true TMP_DIR=$PWD/{tmp} | \
        bwa mem -M -t 4 -p $PWD/{ref} /dev/stdin | \
        picard MergeBamAlignment ALIGNED_BAM=/dev/stdin UNMAPPED_BAM=$PWD/{infile} OUTPUT=$PWD/{outfile} R=$PWD/{ref} TMP_DIR=$PWD/{tmp}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

########################################################################################################################
#########################################---- SORT BAM FILE BY COORDINATE ----##########################################
########################################################################################################################

def sort_coordinates(infile, outfile, temp, done, done_prev):
    """Sort bam file by coordinate."""
    inputs = [infile, done_prev]
    outputs = [done, outfile, outfile.replace('bam', 'bai')]
    options = {'cores': 1, 'memory': "20g", 'walltime': "01:00:00"}
    spec = f"""
        picard SortSam I=$PWD/{infile} O=$PWD/{outfile} SO=coordinate CREATE_INDEX=true TMP_DIR=$PWD/{temp} 
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


################################################################################################################################
######################################## ---- MERGE BAM FILES OF SAME INDIVIDUAL ---- ##########################################
################################################################################################################################

def merge_bams(infile, outfile, done, done_prev):
    """Sort bam file by coordinate."""
    inputs = [infile, done_prev]
    outputs = [done, outfile]
    options = {'cores': 8, 'memory': "20g", 'walltime': "04:00:00"}
    spec = f"""
        samtools merge -@8 {outfile} {" ".join(infile)}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

########################################################################################################################
######################################## ---- MARK AND REMOVE DUPLICATES ---- ##########################################
########################################################################################################################

def mark_and_remove_duplicates_old(infile, intermediate, metrics, outfile, temp, done, done_prev):
    """Mark and remove duplicates."""
    inputs = [infile, done_prev]
    outputs = [done, outfile, metrics]
    options = {'cores': 4, 'memory': '200g', 'walltime': "24:00:00"}
    spec = f"""
        picard SortSam I=$PWD/{infile} O=$PWD/{intermediate} SO=queryname CREATE_INDEX=true TMP_DIR=$PWD/{temp}
        picard MarkDuplicates I=$PWD/{intermediate} M=$PWD/{metrics} O=$PWD/{outfile} REMOVE_DUPLICATES=true CREATE_INDEX=true TMP_DIR=$PWD/{temp}
        rm {intermediate}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def mark_and_remove_duplicates(infile, metrics, outfile, temp, done, done_prev):
    """Mark and remove duplicates."""
    inputs = [infile, done_prev]
    outputs = [done, outfile, metrics]
    options = {'cores': 6, 'memory': '200g', 'walltime': "24:00:00"}
    spec = f"""
        picard MarkDuplicates I=$PWD/{infile} M=$PWD/{metrics} O=$PWD/{outfile} REMOVE_DUPLICATES=true CREATE_INDEX=true TMP_DIR=$PWD/{temp}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

#############################################################################################
#################################---- GET AVG COVERAGE ----##################################
#############################################################################################
def get_coverage(infile, outfile, done, done_prev):
    """Get average coverage across chromosomes at covered sites, for files already coordinate-sorted."""
    inputs = [infile, done_prev]
    outputs = [done, outfile]
    options = {'cores': 1, 'memory': "20g", 'walltime': "02:00:00"}
    spec = f"""
        samtools depth {infile} | awk '{{sum += $3}} END {{print sum / NR}}' > {outfile}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


##########################################################################################
#################################---- CALL VARIANTS ----##################################
##########################################################################################

def call_variants_per_chromosome(chr, ref, infile, outfile, done, done_prev, ploidy="2"):
    """Call variants with the option to set ploidy to 1 or 2 (default is 2)."""
    inputs = [done_prev, ref, infile]
    outputs = [done, outfile, outfile+'.tbi']
    options = {'cores': 32, 'memory': "16g", 'walltime': "06:00:00"}
    spec = f"""
        /mnt/primevo/shared_data/software/gatk/gatk-4.3.0.0/gatk HaplotypeCaller -R {ref} -I {infile} -L {chr} -ploidy {ploidy} --native-pair-hmm-threads 32 -ERC BP_RESOLUTION -O {outfile}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



##########################################################################################
#################################---- COMBINE gVCFs ----##################################
##########################################################################################

# def combine_gvcfs(infiles, ref, outfile, chr, done, done_prev):
#     """Combine gVCF files."""
#     variants = ''.join([f'-V {i} ' for i in infiles])
#     inputs = [done_prev, ref] + [infiles]
#     outputs = [done, outfile]
#     options = {'cores': 1, 'memory': "16g", 'walltime': "48:00:00"}
#     spec = f"""
#         gatk CombineGVCFs -R {ref} {variants} -O {outfile} -L {chr}
#         touch {done}
#     """
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
def combine_gvcfs(infiles, chr, workpath, done, done_prev):
    """Combine gVCF files."""
    variants = ''.join([f'-V {i} ' for i in infiles])
    inputs = [done_prev] + [infiles]
    outputs = [done, workpath]
    options = {'cores': 1, 'memory': "16g", 'walltime': "48:00:00"}
    spec = f"""
        /mnt/primevo/shared_data/software/gatk/gatk-4.3.0.0/gatk --java-options "-Xmx12g -Xms12g" GenomicsDBImport {variants} --genomicsdb-workspace-path {workpath} -L {chr}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


##########################################################################################
#################################---- GENOTYPE gVCFs ----#################################
##########################################################################################

# def genotype_gvcfs(infile, ref, outfile, done, done_prev):
#     """Genotype gVCFs."""
#     inputs = [done_prev, ref, infile]
#     outputs = [done, outfile]
#     options = {'cores': 1, 'memory': "16g", 'walltime': "48:00:00"}
#     spec = f"""
#     gatk GenotypeGVCFs -R {ref} -V {infile} -O {outfile}
#     touch {done}
#     """
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def genotype_gvcfs(infile, ref, outfile, done, done_prev):
    """Genotype gVCFs."""
    inputs = [done_prev, ref, infile.replace('gendb://', '')]
    outputs = [done, outfile, outfile+'.tbi']
    options = {'cores': 1, 'memory': "16g", 'walltime': "48:00:00"}
    spec = f"""
    /mnt/primevo/shared_data/software/gatk/gatk-4.3.0.0/gatk GenotypeGVCFs -R {ref} -V {infile} -O {outfile}
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)








# chr_lst = [
#     "NC_018152.2", "NC_018153.2", "NC_018154.2", "NC_018155.2", 
#     "NC_018156.2", "NC_018157.2", "NC_018158.2", "NC_018159.2", 
#     "NC_018160.2", "NC_018161.2", "NC_018162.2", "NC_018163.2", 
#     "NC_018164.2", "NC_018165.2", "NC_018166.2", "NC_018167.2", 
#     "NC_018168.2", "NC_018169.2", "NC_018170.2", "NC_018171.2", 
#     "NC_018172.2", "NC_020006.2"
# ]

chr_lst = [
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

study_name = "kuderna_2023"
path_to_samples = f"../../../../shared_data/sequencing_data/baboon_other_studies/{study_name}"

dct = {}
n_sample = 3
acc = 0
with open(f"../metadata/sample_names_{study_name}.tsv") as file:
    for line in file:
        if acc == 0:
            acc += 1
            continue
        line = line.rstrip().split("\t")
        if line[1] not in dct: 
            acc += 1
            if (acc-2) == n_sample: break
            dct[line[1]] = [line[3]]
        else: dct[line[1]].append(line[3])
        
print(dct)

# with open(f'../../../../shared_data/sequencing_data/baboon_other_studies/md5_{study_name}.md5', 'r') as file:
#     sample_list = sorted(list(set([f.rstrip().split(' ')[1].split('/')[1].split('_')[0] for f in file])))

acc = 0
for individual in list(dct.keys()):
    for sample_name in dct[individual]:
        if not os.path.exists(f'{path_to_samples}/{sample_name}_1.fastq.gz'): continue
        if not os.path.exists(f'{path_to_samples}/{sample_name}_2.fastq.gz'): continue
        # if all(os.path.exists(f"../done/08_call_variants/{sample_name}_{chr}") for chr in chr_lst): continue
        acc += 1
        gwf.target_from_template(
            f"{sample_name}_FASTQ_to_uBAM", 
            FASTQ_to_uBAM(
                sample_name = sample_name,
                fastq1 = f'{path_to_samples}/{sample_name}_1.fastq.gz',
                fastq2 = f'{path_to_samples}/{sample_name}_2.fastq.gz',
                outfile = f"../steps/01_uBAM_files/{sample_name}.bam",
                done = f"../done/01_uBAM_files/{sample_name}"
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
    
    # for chr in chr_lst:
    #     gwf.target_from_template(
    #         f"{individual}_{chr}_call_variants", 
    #         call_variants_per_chromosome(
    #             chr = chr,
    #             ref = "../data/reference_data/newref/GCF_008728515.1_Panubis1.0_genomic.fna",
    #             infile = f"../steps/06_mark_and_remove_duplicates/{individual}.mapped.nodupes.bam",
    #             outfile = f"../steps/08_call_variants/{individual}_{chr}.g.vcf.gz",
    #             done = f"../done/08_call_variants/{individual}_{chr}",
    #             done_prev = f"../done/07_get_coverage/{individual}"
    #             )
    #         )
    # # if acc == 3: break
    

# for chr in chr_lst:
#     gwf.target_from_template(
#         f"{chr}_combine_gVCFs", 
#         combine_gvcfs(
#             infiles = [f"../steps/08_call_variants/{sample_name}_{chr}.g.vcf.gz" for sample_name in sample_list],
#             chr = chr,
#             workpath = f'../steps/09_combine_gvcfs/{chr}/',
#             done = f"../done/09_combine_gvcfs/done_{chr}",
#             done_prev = [f"../done/08_call_variants/{sample_name}_{chr}" for sample_name in sample_list]
#             )
#         )
#     gwf.target_from_template(
#         f"{chr}_genotype_gVCFs", 
#         genotype_gvcfs(
#             infile =  f"gendb://../steps/09_combine_gvcfs/{chr}/",
#             ref = "../data/reference_data/newref/GCF_008728515.1_Panubis1.0_genomic.fna",
#             outfile = f"../steps/10_final_vcf/{chr}.combined.vcf.gz",
#             done = f"../done/10_final_vcf/done_{chr}",
#             done_prev = f"../done/09_combine_gvcfs/done_{chr}"
#             )
#         )
