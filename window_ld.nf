if( !nextflow.version.matches('>20.0') ) {
    println "This workflow requires Nextflow version 20.0 or greater -- You are running version $nextflow.version"
    println "On QUEST, you can use `module load python/anaconda3.6; source activate /projects/b1059/software/conda_envs/nf20_env`"
    exit 1
}

nextflow.preview.dsl=2

process marker_ld {
    container 'andersenlab/nemascan:20220407173056db3227'
    cpus 4
    memory 30.GB
    time '30m'
    
    publishDir "${params.out}/${sp}/${strain_set}", pattern: "*vcf.gz", overwrite: true
    publishDir "${params.out}/${sp}/${strain_set}", pattern: "*vcf.gz.tbi", overwrite: true
    publishDir "${params.out}/${sp}/${strain_set}", pattern: "*bim", overwrite: true
    publishDir "${params.out}/${sp}/${strain_set}", pattern: "*log", overwrite: true
    publishDir "${params.out}/${sp}/${strain_set}", pattern: "*frq", overwrite: true


    input:
        tuple val(sp), val(strain_set), val(strains), file(vcf), file(index), file(num_chroms), val(MAF)

    output:
        tuple val(strain_set), val(strains), file("${strain_set}_${MAF}.bed"), file("${strain_set}_${MAF}.bim"), file("${strain_set}_${MAF}.fam"), file("${strain_set}_${MAF}.map"), file("${strain_set}_${MAF}.nosex"), file("${strain_set}_${MAF}.ped"), file("${strain_set}_${MAF}.log"), file("${strain_set}_${MAF}_Genotype_Matrix.tsv"),  val(MAF), file("${strain_set}_${MAF}.frq"),  emit: sim_geno
        tuple val(strain_set), val(strains), val(MAF), file("renamed_chroms.vcf.gz"), file("renamed_chroms.vcf.gz.tbi"), emit: renamed_chrom_vcf_to_ld // This output is no longer used. Should be removed in the future.


    """
    bcftools annotate --rename-chrs ${num_chroms} ${vcf} |\\
    bcftools view -s `echo ${strains} | tr -d '\\n'` |\\
    bcftools filter -i N_MISSING=0 -Oz -o renamed_chroms.vcf.gz
    tabix -p vcf renamed_chroms.vcf.gz

   
    plink --vcf renamed_chroms.vcf.gz \\
    --make-bed \\
    --snps-only \\
    --biallelic-only \\
    --set-missing-var-ids @:# \\
    --geno \\
    --recode \\
    --out ${strain_set}_${MAF} \\
    --allow-extra-chr

    plink --vcf renamed_chroms.vcf.gz \\
    --snps-only \\
    --biallelic-only \\
    --set-missing-var-ids @:# \\
    --geno \\
    --freq \\
    --out ${strain_set}_${MAF}

    awk -F":" '\$1=\$1' OFS="\\t" plink.prune.in | \\
    sort -k1,1d -k2,2n > markers.txt
    bcftools query -l renamed_chroms.vcf.gz |\\
    sort > sorted_samples.txt
    
    bcftools view -v snps \\
    -S sorted_samples.txt \\
    -R markers.txt \\
    renamed_chroms.vcf.gz |\\
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
    sed 's/.\\/./NA/g' > ${strain_set}_${MAF}_Genotype_Matrix.tsv
    """
}

