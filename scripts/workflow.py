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

if not os.path.exists("../data/reference_data/newref/GCF_008728515.1_Panubis1.0_genomic.fna.gz"):
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
    options = {'cores': 1, 'memory': '20g', 'walltime': "06:00:00"}
    spec = f"""
        /mnt/primevo/shared_data/software/gatk/gatk-4.3.0.0/gatk CreateSequenceDictionary -R {infile}
        bwa index -a bwtsw {infile}
        samtools faidx {infile}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

if not os.path.exists("../done/00_make_ref_index_and_dict/done"):
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

def FASTQ_to_uBAM(sample_name, fastq1, fastq2, outfile, done, tmp):
    """Make uBAM file from FASTQ file"""
    inputs = [fastq1, fastq2]
    outputs = [done, outfile]
    options = {'cores': 1, 'memory': '40g', 
            # 'walltime': "24:00:00"
            'walltime': "UNLIMITED"
            }
    spec = f"""
        picard FastqToSam F1={fastq1} F2={fastq2} O=$PWD/{outfile} SAMPLE_NAME={sample_name} TMP_DIR=$PWD/{tmp}
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
    options = {'cores': 1, 'memory': '32g', 
            # 'walltime': "08:00:00"
            'walltime': "UNLIMITED"
            }

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
    options = {'cores': 16, 'memory': '30g', 'walltime': "UNLIMITED"}
    spec = f"""
        set -o pipefail
        picard SamToFastq I=$PWD/{infile} FASTQ=/dev/stdout INTERLEAVE=true TMP_DIR=$PWD/{tmp} | \
        bwa mem -M -t 16 -p $PWD/{ref} /dev/stdin | \
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
    options = {'cores': 1, 'memory': "20g", 'walltime': "UNLIMITED"}
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
    options = {'cores': 8, 'memory': "20g", 'walltime': "08:00:00"}
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
    outputs = [done, outfile, metrics, outfile.replace(".bam", ".bai")]
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
    options = {'cores': 1, 'memory': "20g", 'walltime': "08:00:00"}
    spec = f"""
        samtools depth {infile} | awk '{{sum += $3}} END {{print sum / NR}}' > {outfile}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, protect=outputs)


################################################################################################
#################################---- RENAME MERGED FILES ----##################################
################################################################################################
def rename_merged_bam(infile, outfile, done, done_prev):
    """Replace read names by sample name"""
    sample_name = infile.split("/")[-1].split(".")[0]
    inputs = [infile, done_prev]
    outputs = [done, outfile, outfile+".bai"]
    options = {'cores': 1, 'memory': "1g", 'walltime': "UNLIMITED"}
    spec = f"""
        samtools view -H {infile}  | sed "s/SM:[^\t]*/SM:{sample_name}/g" | samtools reheader - {infile} > {outfile}
        samtools index {outfile}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, protect=outputs)



##########################################################################################
#################################---- CALL VARIANTS ----##################################
##########################################################################################

def call_variants_per_chromosome(chr, ref, infile, outfile, done, done_prev, ploidy="2"):
    """Call variants with the option to set ploidy to 1 or 2 (default is 2)."""
    inputs = done_prev+[ref, infile]
    outputs = [done, outfile, outfile+'.tbi']
    options = {'cores': 32, 'memory': "50g", 'walltime': "UNLIMITED"}
    spec = f"""
        /mnt/primevo/shared_data/software/gatk/gatk-4.3.0.0/gatk HaplotypeCaller -R {ref} -I {infile} -L {chr} -ploidy {ploidy} --native-pair-hmm-threads 32 -ERC BP_RESOLUTION -O {outfile}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, protect=outputs)



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
    inputs = done_prev + infiles
    outputs = [done, workpath]
    options = {'cores': 1, 'memory': "100g", 'walltime': "UNLIMITED"}
    spec = f"""
        /mnt/primevo/shared_data/software/gatk/gatk-4.3.0.0/gatk --java-options "-Xmx90g -Xms90g" GenomicsDBImport {variants} --genomicsdb-workspace-path {workpath} -L {chr}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def update_genomicsdb(infiles, workpath, done, done_prev):
    """Update GenomicsDB with new gVCF files."""
    variants = ''.join([f'-V {i} ' for i in infiles])
    inputs = [done_prev] + infiles
    outputs = [done]
    options = {'cores': 1, 'memory': "100g", 'walltime': "UNLIMITED"}
    spec = f"""
        /mnt/primevo/shared_data/software/gatk/gatk-4.3.0.0/gatk --java-options "-Xmx90g -Xms90g" GenomicsDBImport {variants} --genomicsdb-update-workspace-path {workpath}
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
    options = {'cores': 1, 'memory': "80g", 'walltime': "UNLIMITED"}
    spec = f"""
    /mnt/primevo/shared_data/software/gatk/gatk-4.3.0.0/gatk GenotypeGVCFs -R {ref} -V {infile} -O {outfile}
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


############################################################################################
#################################---- CONCATENATE VCFS ----#################################
############################################################################################

def concatenate_vcfs(infiles, outfile, done, done_prev):
    """Concatenate VCFs."""
    inputs = done_prev + infiles
    outputs = [done, outfile]
    options = {'cores': 8, 'memory': "16g", 'walltime': "UNLIMITED"}
    spec = f"""
    bcftools concat --threads 8 -o {outfile} {" ".join(infiles)}
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


##############################################################################################
#################################---- COUNT MAPPED READS ----#################################
##############################################################################################

def count_all_reads(infile, outfile, done, done_prev):
    """Count all reads."""
    inputs = [done_prev, infile]
    outputs = [done, outfile]
    options = {'cores': 1, 'memory': "16g", 'walltime': "UNLIMITED"}
    spec = f"""
    samtools view -c {infile} > {outfile}
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, protect=outputs)

def count_mapped_reads(infile, outfile, done, done_prev):
    """Count mapped (primary aligned) reads."""
    inputs = [done_prev, infile]
    outputs = [done, outfile]
    options = {'cores': 1, 'memory': "16g", 'walltime': "UNLIMITED"}
    spec = f"""
    samtools view -c -F 260 {infile} > {outfile}
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec, protect=outputs)