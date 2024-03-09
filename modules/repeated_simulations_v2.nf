process prepare_genomatrix  {

    tag { CHROM }

    cpus 6
    time '5m'
    memory 5.GB
    container = 'andersenlab/nemascan:20220407173056db3227'

    //memory params.eigen_mem

    input:
        tuple val(sp), val(strain_set), file(bim), file(vcf), file(index), file(selected_snps)

    output:
        tuple val(sp), val(strain_set), file("${sp}_${strain_set}_Genotype_Matrix.tsv")


    """
    bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' ${vcf} |\\
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
        tuple val(sp), val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(SIMID), val(OGS), val(SIMREP), file(create_causal_qtls), file(master_snps_dir)

    output:
        tuple val(sp), val(strain_set), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(SIMREP), val(MAF), file(n_indep_tests), val(SIMID), val(OGS), file("${sp}_${strain_set}_${MAF}_${SIMID}_${SIMREP}_causal_og_vars.txt"), emit: pheno_inputs
        //tuple val(sp), val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(SIMID), val(OGS), val(SIMREP), file(master_snps_dir), file("${sp}_${strain_set}_${MAF}_${SIMID}_${SIMREP}_causal_og_vars.txt")


    """
        python ${create_causal_qtls} ${OGS} ${bim} ${master_snps_dir} ${sp}
        cat causal_og_vars.txt > ${sp}_${strain_set}_${MAF}_${SIMID}_${SIMREP}_causal_og_vars.txt
    """
}