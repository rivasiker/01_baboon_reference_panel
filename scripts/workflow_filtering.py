from gwf import Workflow, AnonymousTarget
import os

gwf = Workflow()

# conda activate baboon_genotyping

##############################################################################################################################
################################################# ---- FILTERING STEPS ---- ##################################################
##############################################################################################################################

def allele_frequency_per_site(infile, outfile, done, done_prev):
    """Calculates allele frequency for sites with a maximum of 2 alleles."""
    inputs = [infile, done_prev]
    outputs = [outfile+".frq", done]
    options = {'cores': 1, 'memory': '20g', 'walltime': "12:00:00"}
    spec = f"""
        vcftools --gzvcf {infile} --freq2 --out {outfile} --max-alleles 2
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def depth_per_individual(infile, outfile, done, done_prev):
    """Calculates mean depth of coverage per individual."""
    inputs = [infile, done_prev]
    outputs = [outfile+".idepth", done]
    options = {'cores': 1, 'memory': '20g', 'walltime': "12:00:00"}
    spec = f"""
        vcftools --gzvcf {infile} --depth --out {outfile}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def depth_per_site(infile, outfile, done, done_prev):
    """Calculates mean depth of coverage per site."""
    inputs = [infile, done_prev]
    outputs = [outfile+".ldepth.mean", done]
    options = {'cores': 1, 'memory': '20g', 'walltime': "12:00:00"}
    spec = f"""
        vcftools --gzvcf {infile} --site-mean-depth --out {outfile}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def quality_per_site(infile, outfile, done, done_prev):
    """Extract quality score per site."""
    inputs = [infile, done_prev]
    outputs = [outfile+".lqual", done]
    options = {'cores': 1, 'memory': '20g', 'walltime': "12:00:00"}
    spec = f"""
        vcftools --gzvcf {infile} --site-quality --out {outfile}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def missingness_per_individual(infile, outfile, done, done_prev):
    """Calculate proportion of missing data per individual."""
    inputs = [infile, done_prev]
    outputs = [outfile+".imiss", done]
    options = {'cores': 1, 'memory': '20g', 'walltime': "12:00:00"}
    spec = f"""
        vcftools --gzvcf {infile} --missing-indv --out {outfile}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def missingness_per_site(infile, outfile, done, done_prev):
    """Calculate proportion of missing data per site."""
    inputs = [infile, done_prev]
    outputs = [outfile+".lmiss", done]
    options = {'cores': 1, 'memory': '20g', 'walltime': "12:00:00"}
    spec = f"""
        vcftools --gzvcf {infile} --missing-site --out {outfile}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def het_and_inbred_per_individual(infile, outfile, done, done_prev):
    """Calculate heterozygosity and inbreeding coefficient per individual."""
    inputs = [infile, done_prev]
    outputs = [outfile+'.het', done]
    options = {'cores': 1, 'memory': '20g', 'walltime': "12:00:00"}
    spec = f"""
        vcftools --gzvcf {infile} --het --out {outfile}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def hardy_weinberg(infile, outfile, done, done_prev):
    """Calculate heterozygosity and inbreeding coefficient per individual."""
    inputs = [infile, done_prev]
    outputs = [outfile+".hwe", done]
    options = {'cores': 1, 'memory': '20g', 'walltime': "12:00:00"}
    spec = f"""
        vcftools --gzvcf {infile} --hardy --out {outfile}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def relatedness2(infile, outfile, done, done_prev):
    """Calculate kinship matrix."""
    inputs = [infile, done_prev]
    outputs = [outfile+".relatedness2", done]
    options = {'cores': 1, 'memory': '20g', 'walltime': "12:00:00"}
    spec = f"""
        vcftools --gzvcf {infile} --relatedness2 --out {outfile}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def relatedness(infile, outfile, done, done_prev):
    """Calculate kinship matrix."""
    inputs = [infile, done_prev]
    outputs = [outfile+".relatedness", done]
    options = {'cores': 1, 'memory': '20g', 'walltime': "12:00:00"}
    spec = f"""
        vcftools --gzvcf {infile} --relatedness --out {outfile}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def perform_qc(outfile, done, infile, done_prev, jobname):
    gwf.target_from_template(
        f"{jobname}_allele_frequency",
        allele_frequency_per_site(
            infile = infile, 
            outfile = f"{outfile}", 
            done = f"{done}_frq", 
            done_prev = done_prev
            )
        )
    
    gwf.target_from_template(
        f"{jobname}_depth_per_individual",
        depth_per_individual(
            infile = infile, 
            outfile = f"{outfile}", 
            done = f"{done}_idepth.mean", 
            done_prev = done_prev
            )
        )
    
    gwf.target_from_template(
        f"{jobname}_depth_per_site",
        depth_per_site(
            infile = infile, 
            outfile = f"{outfile}", 
            done = f"{done}_ldepth.mean", 
            done_prev = done_prev
            )
        )
    
    gwf.target_from_template(
        f"{jobname}_quality_per_site",
        quality_per_site(
            infile = infile, 
            outfile = f"{outfile}", 
            done = f"{done}_lqual", 
            done_prev = done_prev
            )
        )
    
    gwf.target_from_template(
        f"{jobname}_missingness_per_individual",
        missingness_per_individual(
            infile = infile, 
            outfile = f"{outfile}", 
            done = f"{done}_imiss", 
            done_prev = done_prev
            )
        )
    
    gwf.target_from_template(
        f"{jobname}_missingness_per_site",
        missingness_per_site(
            infile = infile, 
            outfile = f"{outfile}", 
            done = f"{done}_lmiss", 
            done_prev = done_prev
            )
        )
    
    gwf.target_from_template(
        f"{jobname}_het_and_inbred_per_individual",
        het_and_inbred_per_individual(
            infile = infile, 
            outfile = f"{outfile}", 
            done = f"{done}_het", 
            done_prev = done_prev
            )
        )

    gwf.target_from_template(
        f"{jobname}_hardy_weinberg",
        hardy_weinberg(
            infile = infile, 
            outfile = f"{outfile}", 
            done = f"{done}_hwe", 
            done_prev = done_prev
            )
        )
    
    gwf.target_from_template(
        f"{jobname}_relatedness",
        relatedness(
            infile = infile, 
            outfile = f"{outfile}", 
            done = f"{done}_relatedness", 
            done_prev = done_prev
            )
        )

    gwf.target_from_template(
        f"{jobname}_relatedness2",
        relatedness2(
            infile = infile, 
            outfile = f"{outfile}", 
            done = f"{done}_relatedness2", 
            done_prev = done_prev
            )
        )




############################################################################################
#################################---- QUALITY CONTROL ----##################################
############################################################################################


# study_name_list = ["NC_044995.1.combined"]
study_name_list = ["robinson_2019"]

for study_name in study_name_list:

    outfile_prexif = "../steps/13_qc"
    done_prexif = "../done/13_qc"
    infile = f"../steps/12_concatenated_vcf/{study_name}.vcf.gz"
    done_prev = f"../done/12_concatenated_vcf/done_{study_name}"
    # infile = f"../steps/11_final_vcf/kuderna_2023/{study_name}.vcf.gz"
    # done_prev = f"../done/11_final_vcf/done_kuderna_2023_NC_044995.1"
    outfile = f"{outfile_prexif}/{study_name}"
    done = f"{done_prexif}/{study_name}"

    perform_qc(outfile, done, infile, done_prev, study_name)

study_name_list = ["robinson_2019_and_kuderna_2023", "amboseli_all"]
chr = 'NC_044995.1'    
for study_name in study_name_list:
    outfile_prexif = "../steps/13_qc"
    done_prexif = "../done/13_qc"
    infile = f"../steps/11_final_vcf/{study_name}/{chr}.combined.vcf.gz"
    done_prev = f"../done/11_final_vcf/done_{study_name}_{chr}"
    # infile = f"../steps/11_final_vcf/kuderna_2023/{study_name}.vcf.gz"
    # done_prev = f"../done/11_final_vcf/done_kuderna_2023_NC_044995.1"
    outfile = f"{outfile_prexif}/{study_name}_{chr}"
    done = f"{done_prexif}/{study_name}_{chr}"

    perform_qc(outfile, done, infile, done_prev, study_name)

######################################################################################
#################################---- FILTERING ----##################################
######################################################################################

def qc_filtering(infile, outfile, done, done_prev, maf, miss, qual, min_depth, max_depth):
    """Calculate heterozygosity and inbreeding coefficient per individual."""
    inputs = [infile] + done_prev
    outputs = [outfile+".filtered.vcf.gz", done]
    options = {'cores': 1, 'memory': '20g', 'walltime': "12:00:00"}
    spec = f"""
        vcftools --gzvcf {infile} \
        --remove-indels --maf {maf} --max-missing {miss} --minQ {qual} \
        --min-meanDP {min_depth} --max-meanDP {max_depth} \
        --minDP {min_depth} --maxDP {max_depth} --out {outfile} --recode --stdout | gzip -c > \
        {outfile+".filtered.vcf.gz"}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

maf=0.1
miss=0.9
qual=30
min_depth=10
max_depth=80

study_name_list = ["robinson_2019"]

for study_name in study_name_list:

    outfile_prexif = "../steps/14_qc_filtering"
    done_prexif = "../done/14_qc_filtering"
    done_prexif_prev = "../done/13_qc"
    infile = f"../steps/12_concatenated_vcf/{study_name}.vcf.gz"
    done_prev = [f"{done_prexif_prev}/{study_name}_{i}" for i in ["frq", "idepth.mean", "ldepth.mean", "lqual", "imiss", "lmiss", "het"]]

    gwf.target_from_template(
            f"{study_name}_qc_filtering",
            qc_filtering(
                infile = infile, 
                outfile = f"{outfile_prexif}/{study_name}", 
                done = f"{done_prexif}/{study_name}", 
                done_prev = done_prev,
                maf = maf,
                miss = miss,
                qual = qual,
                min_depth = min_depth,
                max_depth = max_depth
                )
            )
    
study_name_list = ["robinson_2019_and_kuderna_2023", "amboseli_all"]
chr = 'NC_044995.1'

for study_name in study_name_list:

    outfile_prexif = "../steps/14_qc_filtering"
    done_prexif = "../done/14_qc_filtering"
    done_prexif_prev = "../done/13_qc"
    infile = f"../steps/11_final_vcf/{study_name}/{chr}.combined.vcf.gz"
    done_prev = [f"{done_prexif_prev}/{study_name}_{chr}_{i}" for i in ["frq", "idepth.mean", "ldepth.mean", "lqual", "imiss", "lmiss", "het"]]

    gwf.target_from_template(
            f"{study_name}_qc_filtering",
            qc_filtering(
                infile = infile, 
                outfile = f"{outfile_prexif}/{study_name}_{chr}", 
                done = f"{done_prexif}/{study_name}_{chr}", 
                done_prev = done_prev,
                maf = maf,
                miss = miss,
                qual = qual,
                min_depth = min_depth,
                max_depth = max_depth
                )
            )


############################################################################################################
#################################---- QUALITY CONTROL AFTER FILTERING ----##################################
############################################################################################################


study_name_list = ["robinson_2019"]
for study_name in study_name_list:
    outfile_prexif = "../steps/13_qc"
    done_prexif = "../done/13_qc"
    infile = f"../steps/14_qc_filtering/{study_name}.filtered.vcf.gz"
    done_prev = f"../done/14_qc_filtering/{study_name}"
    outfile = f"{outfile_prexif}/{study_name}_filtered"
    done = f"{done_prexif}/{study_name}_filtered"
    perform_qc(outfile, done, infile, done_prev, study_name+"_filtered")

# study_name_list = ["robinson_2019_and_kuderna_2023", "amboseli_all"]
# chrom = "NC_044995.1"
# for study_name in study_name_list:
#     outfile_prexif = "../steps/13_qc"
#     done_prexif = "../done/13_qc"
#     infile = f"../steps/11_final_vcf/{study_name}/{chrom}.combined.vcf.gz"
#     done_prev = f"../done/11_final_vcf/done_{study_name}_{chrom}"
#     outfile = f"{outfile_prexif}/{study_name}_{chrom}_filtered"
#     done = f"{done_prexif}/{study_name}_{chrom}_filtered"
#     perform_qc(outfile, done, infile, done_prev, study_name+"_filtered")


################################################################################
#################################---- PCA ----##################################
################################################################################

def run_pca(infile, done, done_prev, study_name):
    """Calculate heterozygosity and inbreeding coefficient per individual."""
    inputs = [infile, done_prev]
    outputs = [study_name+".prune.in", study_name+".eigenval", study_name+".eigenvec", 
               study_name+".bed", study_name+".bim", study_name+".fam", 
               done]
    options = {'cores': 1, 'memory': '20g', 'walltime': "12:00:00"}
    spec = f"""
        plink --vcf {infile} --double-id --allow-extra-chr --set-missing-var-ids @:# \
            --indep-pairwise 50 10 0.1 --out {study_name}
        plink --vcf {infile} --double-id --allow-extra-chr --set-missing-var-ids @:# \
            --extract {study_name}.prune.in \
            --make-bed --pca --out {study_name}
        touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# study_name_list = ["NC_044995.1.combined"]
study_name_list = ["robinson_2019"]
for s in study_name_list:
    study_name = s+".filtered"
    gwf.target_from_template(
            f"{study_name}_filtered_run_pca",
            run_pca(
                infile = f"../steps/14_qc_filtering/{study_name}.vcf.gz", 
                done = f"../done/15_run_pca/{study_name}", 
                done_prev = f"../done/14_qc_filtering/{s}",
                study_name = f"../steps/15_run_pca/{study_name}"
                )
            )
    

study_name_list = ["robinson_2019_and_kuderna_2023"]
chr = 'NC_044995.1'
for s in study_name_list:
    study_name = f"{s}_{chr}.filtered"
    gwf.target_from_template(
            f"{study_name}_run_pca",
            run_pca(
                infile = f"../steps/14_qc_filtering/{study_name}.vcf.gz", 
                done = f"../done/15_run_pca/{study_name}", 
                done_prev = f"../done/14_qc_filtering/{s}_{chr}",
                study_name = f"../steps/15_run_pca/{study_name}"
                )
            )


study_name_list = ["amboseli_all"]
chr = 'NC_044995.1'
for study_name in study_name_list:
    gwf.target_from_template(
            f"{study_name}_run_pca_unfiltered",
            run_pca(
                infile = f"../steps/11_final_vcf/{study_name}/{chr}.combined.vcf.gz", 
                done = f"../done/15_run_pca/{study_name}_unfiltered", 
                done_prev = f"../done/14_qc_filtering/{study_name}_{chr}",
                study_name = f"../steps/15_run_pca/{study_name}_unfiltered"
                )
            )