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
        tuple val(sp), val(strain_set), path("${sp}_${strain_set}_Genotype_Matrix.tsv")


    """
    bcftools view -T ${selected_snps} -Ou ${vcf} |\\
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
    
    cpus 3
    time '5m'
    memory 5.GB
    container = 'andersenlab/nemascan:20220407173056db3227'
    executor 'local'
    
    publishDir "${params.out}/${sp}/${strain_set}/Markers", mode: 'copy', pattern: "*bim"

    input:
        tuple val(sp), val(strain_set), path(vcf), path(vcf_index), path(selected_snps)

    output:
        tuple val(sp), val(strain_set), path("TO_SIMS.bed"), path("TO_SIMS.bim"), path("TO_SIMS.fam"), path("TO_SIMS.map"), path("TO_SIMS.nosex"), path("TO_SIMS.ped"), path("TO_SIMS.log")


    """
    bcftools view -T ${selected_snps} -Oz -o renamed_chroms.vcf.gz ${vcf}
    tabix -p vcf renamed_chroms.vcf.gz
    plink --vcf renamed_chroms.vcf.gz \\
    --make-bed \\
    --snps-only \\
    --biallelic-only \\
    --geno \\
    --recode \\
    --out TO_SIMS \\
    --allow-extra-chr
    """


}

process chrom_eigen_variants_sims_repeated  {

    tag { CHROM }

    cpus 3
    time '5m'
    memory 5.GB
    executor 'local'
    container = 'andersenlab/nemascan:20220407173056db3227'

    //memory params.eigen_mem
    publishDir "${params.out}/${sp}/${strain_set}/Markers/eigendecomp", mode: 'copy', pattern: "*_independent_snvs.csv"  

    input:
        tuple val(CHROM), val(sp), val(strain_set), file(gm), file(get_genomatrix_eigen)

    output:
        tuple val(sp), val(strain_set), file("${CHROM}_${sp}_${strain_set}_independent_snvs.csv")


    """
        cat ${gm} |\\
        awk -v chrom="${CHROM}" '{if(\$1 == "CHROM" || \$1 == chrom) print}' > ${CHROM}_gm.tsv
        Rscript --vanilla ${get_genomatrix_eigen} ${CHROM}_gm.tsv ${CHROM}
        mv ${CHROM}_independent_snvs.csv ${CHROM}_${sp}_${strain_set}_independent_snvs.csv
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