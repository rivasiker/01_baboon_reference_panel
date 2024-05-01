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






# study_name_list = ["NC_044995.1.combined"]
study_name_list = ["robinson_2019"]

for study_name in study_name_list:

    outfile_prexif = "../steps/13_qc"
    done_prexif = "../done/13_qc"
    infile = f"../steps/12_concatenated_vcf/{study_name}.vcf.gz"
    done_prev = f"../done/12_concatenated_vcf/done_{study_name}"
    # infile = f"../steps/11_final_vcf/kuderna_2023/{study_name}.vcf.gz"
    # done_prev = f"../done/11_final_vcf/done_kuderna_2023_NC_044995.1"

    gwf.target_from_template(
        f"{study_name}_allele_frequency",
        allele_frequency_per_site(
            infile = infile, 
            outfile = f"{outfile_prexif}/{study_name}", 
            done = f"{done_prexif}/{study_name}_frq", 
            done_prev = done_prev
            )
        )
    
    gwf.target_from_template(
        f"{study_name}_depth_per_individual",
        depth_per_individual(
            infile = infile, 
            outfile = f"{outfile_prexif}/{study_name}", 
            done = f"{done_prexif}/{study_name}_idepth.mean", 
            done_prev = done_prev
            )
        )
    
    gwf.target_from_template(
        f"{study_name}_depth_per_site",
        depth_per_site(
            infile = infile, 
            outfile = f"{outfile_prexif}/{study_name}", 
            done = f"{done_prexif}/{study_name}_ldepth.mean", 
            done_prev = done_prev
            )
        )
    
    gwf.target_from_template(
        f"{study_name}_quality_per_site",
        quality_per_site(
            infile = infile, 
            outfile = f"{outfile_prexif}/{study_name}", 
            done = f"{done_prexif}/{study_name}_lqual", 
            done_prev = done_prev
            )
        )
    
    gwf.target_from_template(
        f"{study_name}_missingness_per_individual",
        missingness_per_individual(
            infile = infile, 
            outfile = f"{outfile_prexif}/{study_name}", 
            done = f"{done_prexif}/{study_name}_imiss", 
            done_prev = done_prev
            )
        )
    
    gwf.target_from_template(
        f"{study_name}_missingness_per_site",
        missingness_per_site(
            infile = infile, 
            outfile = f"{outfile_prexif}/{study_name}", 
            done = f"{done_prexif}/{study_name}_lmiss", 
            done_prev = done_prev
            )
        )
    
    gwf.target_from_template(
        f"{study_name}_het_and_inbred_per_individual",
        het_and_inbred_per_individual(
            infile = infile, 
            outfile = f"{outfile_prexif}/{study_name}", 
            done = f"{done_prexif}/{study_name}_het", 
            done_prev = done_prev
            )
        )

