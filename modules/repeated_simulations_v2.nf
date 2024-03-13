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
process collect_eigen_variants_sims_repeated {

    executor 'local'

    //publishDir "${params.out}/Genotype_Matrix", mode: 'copy'
    publishDir "${params.out}/${sp}/${strain_set}/Markers/eigendecomp", mode: 'copy', pattern: "*_total_independent_tests.txt"  
    container = 'andersenlab/nemascan:20220407173056db3227'
    cpus 1
    time '5m'
    memory 5.GB



    input:
        tuple val(sp), val(strain_set), path(tests)


    output:
        tuple val(sp), val(strain_set), path("${sp}_${strain_set}_total_independent_tests.txt")

    """
        cat *independent_snvs.csv |\\
        grep -v inde |\\
        awk '{s+=\$1}END{print s}' > ${sp}_${strain_set}_total_independent_tests.txt
    """

}

process simulate_orthogroup_effects {
    label 'causal_ogs'
    //errorStrategy 'ignore'
    executor 'local'
    conda '/home/rjm6024/.conda/envs/vcf_stats_1.0'
    //conda  '/Users/ryanmckeown/anaconda3/envs/dev_env'

    publishDir "${params.out}/Simulations/${SIMREP}/causal_vars", pattern: "*txt", overwrite: true

    input:
        tuple val(sp), val(strain_set), path(all_pop_snps_anno_bim), path(create_causal_qtls), val(SIMREP), val(OGS) 

    output:
        tuple val(sp), val(strain_set), val(SIMREP), path("${sp}_${strain_set}_${SIMREP}_causal_og_vars.txt")


    """
        python ${create_causal_qtls} ${OGS} ${all_pop_snps_anno_bim} ${sp}
        cat causal_og_vars.txt > ${sp}_${strain_set}_${SIMREP}_causal_og_vars.txt
    """
}

process sim_phenos {
    label 'sim_map_phenos'
    
    executor 'local'
    //errorStrategy 'retry'
    container = 'andersenlab/nemascan:20220407173056db3227'
    //publishDir "${params.out}/Simulations/${sp}/${SIMID}/Mappings", pattern: "*fastGWA", overwrite: true
    //publishDir "${params.out}/Simulations/${sp}/${SIMID}/Mappings", pattern: "*loco.mlma", overwrite: true
    publishDir "${params.out}/Phenotypes/", pattern: "*.phen", overwrite: true
    publishDir "${params.out}/Phenotypes/", pattern: "*.par", overwrite: true

    cpus 5
    time '20m'
    memory 10.GB

    input:
        tuple val(sp), val(strain_set), path(bed), path(bim), path(fam), path(map), path(nosex), path(ped), path(log), path(gm), path(n_indep_tests), val(SIMREP), path(loci), val(H2), path(check_vp)

    output:
        tuple val(sp), val(strain_set), val(SIMREP), path("${SIMREP}}_${sp}_${strain_set}_sims.phen"), path("${SIMREP}_${sp}_${strain_set}_sims.par")
    """
        gcta64 --bfile TO_SIMS \\
         --simu-qt \\
         --simu-causal-loci ${loci} \\
         --simu-hsq ${H2} \\
         --simu-rep 1 \\
         --thread-num 5 \\
         --out ${SIMREP}_${sp}_${strain_set}_sims
    """
}

process simulate_map_phenotypes {

    label 'sim_map_phenos'
    //tag {"${SIMREP} - ${H2} - ${MAF}"}

    //errorStrategy 'retry'
    //container = 'andersenlab/nemascan:20220407173056db3227'
    publishDir "${params.out}/Simulations/${sp}/${SIMID}/Mappings", pattern: "*fastGWA", overwrite: true
    publishDir "${params.out}/Simulations/${sp}/${SIMID}/Mappings", pattern: "*loco.mlma", overwrite: true
    publishDir "${params.out}/Simulations/${sp}/${SIMID}/Phenotypes", pattern: "*.phen", overwrite: true
    publishDir "${params.out}/Simulations/${sp}/${SIMID}/Phenotypes", pattern: "*.par", overwrite: true

    cpus 5
    time '20m'
    memory 10.GB



    input:
        tuple val(sp), val(strain_set), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), file(n_indep_tests), val(SIMID), file(loci), val(H2), file(check_vp)

    output:
        tuple file("TO_SIMS_${SIMREP}_${MAF}_${SIMID}_${sp}_${strain_set}.bed"), file("TO_SIMS_${SIMREP}_${MAF}_${SIMID}_${sp}_${strain_set}.bim"), file("TO_SIMS_${SIMREP}_${MAF}_${SIMID}_${sp}_${strain_set}.fam"), file("TO_SIMS_${SIMREP}_${MAF}_${SIMID}_${sp}_${strain_set}.map"), file("TO_SIMS_${SIMREP}_${MAF}_${SIMID}_${sp}_${strain_set}.nosex"), file("TO_SIMS_${SIMREP}_${MAF}_${SIMID}_${sp}_${strain_set}.ped"), file("TO_SIMS_${SIMREP}_${MAF}_${SIMID}_${sp}_${strain_set}.log"), val(SIMREP), file(loci), file("${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sims.phen"), file("${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sims.par"), emit: sim_phen_output
        tuple file("${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_lmm-exact_inbred.fastGWA"), file("${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_lmm-exact.loco.mlma"), file("${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_lmm-exact.log"), file("${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_lmm-exact_inbred.log"), emit: sim_GCTA_mapping_results
        
        path "${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_lmm-exact_inbred.fastGWA", emit: lmm_exact_inbred_analyze_sims
        path "${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_lmm-exact_inbred_pca.fastGWA", emit: lmm_exact_inbred_pca_analyze_sims
        path "${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_lmm-exact.loco.mlma", emit: lmm_exact_loco_analyze_sims
        path "${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_lmm-exact_pca.loco.mlma", emit: lmm_exact_loco_pca_analyze_sims
        path "${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sims.phen", emit: simphen_analyze_sims
        path "${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sims.par", emit: simgen_analyze_sims
        tuple val(sp), val(strain_set), val(SIMREP), val(H2), file(loci), file(gm), file(n_indep_tests), val(MAF),val(SIMID), val(OGS), file("${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_lmm-exact_inbred.fastGWA"), file("${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_lmm-exact_inbred_pca.fastGWA"),file("${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_lmm-exact.loco.mlma"), file("${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_lmm-exact_pca.loco.mlma"), file("${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sims.phen"), file("${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sims.par"),emit: gcta_intervals

    """
    gcta64 --bfile TO_SIMS \\
         --simu-qt \\
         --simu-causal-loci ${loci} \\
         --simu-hsq ${H2} \\
         --simu-rep 1 \\
         --thread-num 5 \\
         --out ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sims
    plink --bfile TO_SIMS \\
        --make-bed \\
        --snps-only \\
        --biallelic-only \\
        --maf ${MAF} \\
        --set-missing-var-ids @:# \\
        --geno \\
        --recode \\
        --out TO_SIMS_${SIMREP}_${MAF}_${SIMID}_${sp}_${strain_set} \\
        --allow-extra-chr \\
        --pheno ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sims.phen
    gcta64 --bfile TO_SIMS_${SIMREP}_${MAF}_${SIMID}_${sp}_${strain_set} \\
            --autosome --maf ${MAF} --make-grm \\
            --out TO_SIMS_${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_gcta_grm \\
            --thread-num 5
    gcta64 --bfile TO_SIMS_${SIMREP}_${MAF}_${SIMID}_${sp}_${strain_set} \\
            --autosome --maf ${MAF} --make-grm-inbred \\
            --out TO_SIMS_${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_gcta_grm_inbred \\
            --thread-num 5
    gcta64 --grm TO_SIMS_${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_gcta_grm_inbred \\
            --pheno ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sims.phen \\
            --reml --out check_vp \\
            --thread-num 5
    
    python ${check_vp} --check_vp check_vp.hsq --simulated_phenos ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sims.phen 

    if [[ -f "new_phenos.temp" ]]
    then
        rm ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sims.phen
        mv new_phenos.temp ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sims.phen
    fi

    gcta64 --grm TO_SIMS_${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_gcta_grm \\
           --make-bK-sparse ${params.sparse_cut} \\
           --out ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sparse_grm \\
           --thread-num 5
    gcta64 --grm TO_SIMS_${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_gcta_grm \\
           --pca 1 \\
           --out ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sparse_grm \\
           --thread-num 5
    gcta64 --mlma-loco \\
           --bfile TO_SIMS_${SIMREP}_${MAF}_${SIMID}_${sp}_${strain_set} \\
           --grm ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sparse_grm \\
           --out ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_lmm-exact \\
           --pheno ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sims.phen \\
           --maf ${MAF} \\
           --thread-num 5
    gcta64 --mlma-loco \\
           --bfile TO_SIMS_${SIMREP}_${MAF}_${SIMID}_${sp}_${strain_set} \\
           --grm ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sparse_grm \\
           --qcovar ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sparse_grm.eigenvec \\
           --out ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_lmm-exact_pca \\
           --pheno ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sims.phen \\
           --maf ${MAF} \\
           --thread-num 5

    gcta64 --grm TO_SIMS_${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_gcta_grm_inbred \\
          --make-bK-sparse ${params.sparse_cut} \\
          --out ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sparse_grm_inbred \\
          --thread-num 5
    gcta64 --grm TO_SIMS_${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_gcta_grm_inbred \\
          --pca 1 \\
          --out ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sparse_grm_inbred \\
          --thread-num 5
    gcta64 --fastGWA-lmm-exact \\
          --grm-sparse ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sparse_grm_inbred \\
          --bfile TO_SIMS_${SIMREP}_${MAF}_${SIMID}_${sp}_${strain_set} \\
          --out ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_lmm-exact_inbred \\
          --pheno ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sims.phen \\
          --maf ${MAF} \\
          --thread-num 5
    gcta64 --fastGWA-lmm-exact \\
          --grm-sparse ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sparse_grm_inbred \\
          --qcovar ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sparse_grm_inbred.eigenvec \\
          --bfile TO_SIMS_${SIMREP}_${MAF}_${SIMID}_${sp}_${strain_set} \\
          --out ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_lmm-exact_inbred_pca \\
          --pheno ${SIMREP}_${H2}_${MAF}_${SIMID}_${sp}_${strain_set}_sims.phen \\
          --maf ${MAF} \\
          --thread-num 5
    """
}
