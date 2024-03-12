process prepare_sim_gm  {

    
    
    memory 5.GB
    //container = 'biocontainers/bcftools:v1.9-1-deb_cv1'
    //container = 'andersenlab/nemascan:20220407173056db3227'
    ///Users/ryanmckeown/anaconda3/envs/gm_matrix_test
    executor 'local'

    //memory params.eigen_mem
    publishDir "${params.out}/${sp}/${strain_set}/Markers", mode: 'copy', pattern: "*Genotype_Matrix.tsv"  
    //publishDir "${params.out}/${sp}/${strain_set}/Markers", mode: 'copy', pattern: "*.vcf.gz"  


    input:
        tuple val(sp), val(strain_set), path(vcf), path(vcf_index), path(selected_snps)

    output:
        tuple val(sp), val(strain_set), path(vcf), path("${sp}_${strain_set}_Genotype_Matrix.tsv")


    """
    bcftools view -S ${selected_snps} -Ou ${vcf} |\\
    bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' |\\
    sed 's/[[# 0-9]*]//g' |\\
    sed 's/:GT//g' |\\
    sed 's/0|0/-1/g' |\\
    sed 's/1|1/1/g' |\\
    sed 's/0|1/NA/g' |\\
    sed 's/1|0/NA/g' |\\
    sed 's/.|./NA/g'  |\\
    sed 's/0\\/0/-1/g' |\\
    sed 's/1\\/1/1/g'  |\\
    sed 's/0\\/1/NA/g' |\\
    sed 's/1\\/0/NA/g' |\\
    sed 's/.\\/./NA/g' > ${sp}_${strain_set}_Genotype_Matrix.tsv

    """

}

process prepare_sim_plink {
    cpus 6
    time '5m'
    memory 5.GB
    container = 'andersenlab/nemascan:20220407173056db3227'
    executor 'local'

    //memory params.eigen_mem

    input:
        tuple val(sp), val(strain_set), file(vcf), file(vcf_index), file(selected_snps)

    output:
        tuple val(sp), val(strain_set), file("${sp}_${strain_set}_Genotype_Matrix.tsv")


    """
    bcftools view -S ${selected_snps} -Ou ${vcf} |\\
    bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' |\\
    sed 's/[[# 0-9]*]//g' |\\
    sed 's/:GT//g' |\\
    sed 's/0|0/-1/g' |\\
    sed 's/1|1/1/g' |\\
    sed 's/0|1/NA/g' |\\
    sed 's/1|0/NA/g' |\\
    sed 's/.|./NA/g'  |\\
    sed 's/0\\/0/-1/g' |\\
    sed 's/1\\/1/1/g'  |\\
    sed 's/0\\/1/NA/g' |\\
    sed 's/1\\/0/NA/g' |\\
    sed 's/.\\/./NA/g' > ${sp}_${strain_set}_Genotype_Matrix.tsv

    """


}

process chrom_eigen_variants_sims_repeated  {

    tag { CHROM }

    cpus 6
    time '5m'
    memory 5.GB
    container = 'andersenlab/nemascan:20220407173056db3227'

    //memory params.eigen_mem

    input:
        tuple val(CHROM), val(sp), val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno), val(MAF), file(get_genomatrix_eigen)

    output:
        tuple val(sp), val(strain_set), val(strains), val(MAF), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno), emit: sim_geno_meta
        tuple val(sp), val(strain_set), val(strains), val(MAF), file("${CHROM}_${sp}_${strain_set}_${MAF}_independent_snvs.csv"), emit: sim_geno_eigen_join


    """
        cat ${geno} |\\
        awk -v chrom="${CHROM}" '{if(\$1 == "CHROM" || \$1 == chrom) print}' > ${CHROM}_gm.tsv
        Rscript --vanilla ${get_genomatrix_eigen} ${CHROM}_gm.tsv ${CHROM}
        mv ${CHROM}_independent_snvs.csv ${CHROM}_${sp}_${strain_set}_${MAF}_independent_snvs.csv
    """

}

process simulate_orthogroup_effects {
    label 'causal_ogs'
    errorStrategy 'ignore'
    executor 'local'
    conda '/home/rjm6024/.conda/envs/vcf_stats_1.0'


    input:
        tuple val(OGS), val(SIMREP), val(sp), val(strain_set), file(bim), file(create_causal_qtls)

    output:
        tuple val(sp), val(strain_set), file("${sp}_${strain_set}_${MAF}_${SIMID}_${SIMREP}_causal_og_vars.txt")


    """
        python ${create_causal_qtls} ${OGS} ${bim} ${sp}
        cat causal_og_vars.txt > ${sp}_${strain_set}_${MAF}_${SIMID}_${SIMREP}_causal_og_vars.txt
    """
}