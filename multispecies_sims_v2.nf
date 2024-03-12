#! usr/bin/env nextflow
if( !nextflow.version.matches('>20.0') ) {
    println "This workflow requires Nextflow version 20.0 or greater -- You are running version $nextflow.version"
    println "On QUEST, you can use `module load python/anaconda3.6; source activate /projects/b1059/software/conda_envs/nf20_env`"
    exit 1
}

nextflow.enable.dsl=2
date = new Date().format( 'yyyyMMdd' )
params.out = "Analysis_Results-${date}"


params.bin_dir = "${workflow.projectDir}/bin" // this is different for gcp
params.sparse_cut = 0.05
params.group_qtl = 1000
params.ci_size = 150
sthresh= "BF"
params.maf = 0.05

include {prepare_sim_gm; prepare_sim_plink; chrom_eigen_variants_sims_repeated} from './modules/repeated_simulations_v2.nf'

// load the population data from the input folder

params.proc_data 
//ce_pop_id = params.ce_pop_id
//cb_pop_id = params.cb_pop_id
//ct_pop_id = params.ct_pop_id

ce_pop_id = "ce.fullpop"
cb_pop_id = "cb.fullpop"
ct_pop_id = "ct.fullpop"

// load the data to create genotype matrix


// ce_geno = ["c_elegans",
//     "${ce_pop_id}",
//     file("${params.proc_data}/${ce_pop_id}/${ce_pop_id}_0.00.bim"), // the bim file containing all SNPs 
//     file("${params.proc_data}/${ce_pop_id}/renamed_chroms.vcf.gz"), // the VCF file thats been filtered to only strains 
//     file("${params.proc_data}/${ce_pop_id}/renamed_chroms.vcf.gz.tbi"), // the index file for the VCF
//     file("${params.proc_data}/${ce_pop_id}/selected_snps.txt") // the list of SNPs selected as gene based markers
//     ]


// Channel.from( ["c_elegans", "${ce_pop_id}"], ["c_briggsae","${cb_pop_id}"], ["c_tropicalis", "${ct_pop_id}"] ) \
//     | map { sp, strain_set -> [sp, \
//                                 strain_set, \
//                                 file("${params.proc_data}/${sp}/${strain_set}/${strain_set}_0.00.bim.bed.annotated"), \
//                                 file("${params.proc_data}/${sp}/${strain_set}/renamed_chroms.vcf.gz"), \
//                                 file("${params.proc_data}/${sp}/${strain_set}/renamed_chroms.vcf.gz.tbi") \
//                                 ]} \
//     | view()

//input channels for GM with all species
// prep_gm_ins = Channel.from( ["c_elegans", "${ce_pop_id}"], ["c_briggsae","${cb_pop_id}"], ["c_tropicalis", "${ct_pop_id}"] ) \
//     | map { sp, strain_set -> [sp, \
//                                 strain_set, \
//                                 file("${params.proc_data}/${sp}/${strain_set}/renamed_chroms.vcf.gz"), \
//                                 file("${params.proc_data}/${sp}/${strain_set}/renamed_chroms.vcf.gz.tbi"), \
//                                 file("${params.proc_data}/${sp}/${strain_set}/selected_snps.txt") \
//                                 ]} \
//     | prepare_sim_gm

//input channel for GM with one species for testing

workflow{
prep_gm_ins = Channel.from( ["c_elegans", "${ce_pop_id}"], ["c_briggsae", "${cb_pop_id}"] ) \
    .map { sp, strain_set -> [sp, \
                                strain_set, \
                                "/projects/b1059/projects/Ryan/ortholog_sims/pipeline_dev/Caeno_Scan/test_data/test_gm/${sp}/${sp}.vcf.gz", \
                                "/projects/b1059/projects/Ryan/ortholog_sims/pipeline_dev/Caeno_Scan/test_data/test_gm/${sp}/${sp}.vcf.gz.tbi", \
                                "/projects/b1059/projects/Ryan/ortholog_sims/pipeline_dev/Caeno_Scan/test_data/test_gm/${sp}/${sp}_all_snps.txt" \
                                ]} \
    // create a tuple
    .map { sp, strain_set, vcf, vcf_tbi, snp_list -> [sp, strain_set, vcf, vcf_tbi, snp_list]} \
    | prepare_sim_plink
    // eigen
    //contigs = Channel.from("1")
    //contigs = Channel.from(["1", "2", "3", "4", "5", "6"]) //Parallelize by chrom
    //contigs.combine(prepare_sim_gm.out) // Combine with Plink files and Genotype matrix + Sim INFO
    //    .combine(Channel.fromPath("bin/Get_GenoMatrix_Eigen.R")) | chrom_eigen_variants_sims_repeated
    


}
// load the data to simulate phenotypes 
// sp_causal_snps_inputs = [
//     ["c_elegans",
//     "${ce_pop_id}",
//     file("${params.proc_data}/${ce_pop_id}/${ce_pop_id}_0.00.bim.bed.annotated") // the list of SNPs in the population annotated with gene id and OG status
//     ],
//     ["c_briggsae",
//     "${cb_pop_id}",
//     file("${params.proc_data}/${cb_pop_id}/${cb_pop_id}_0.00.bim.bed.annotated") // the list of SNPs in the population annotated with gene id and OG status
//     ]
// ]



// create an input channel from the sp_causal_snps_inputs

//sp_causal_snps_inputs = Channel.from(sp_causal_snps_inputs)

//load the simulation key files
//File sim_key_file = new File("test_data/repeated_sim_keys.txt") ;

//ce_anno_snp = file("${params.proc_data}/${ce_pop_id}/${ce_pop_id}_0.00.bim.bed.annotated")
//cb_anno_snp = file("${params.proc_data}/${cb_pop_id}/${cb_pop_id}_0.00.bim.bed.annotated")
//ct_anno_snp = file("${params.proc_data}/${ct_pop_id}/${ct_pop_id}_0.00.bim.bed.annotated")

//anno_snps = [ce_anno_snp, cb_anno_snp, ct_anno_snp]

// join the simulation data to the input file that specifies simulation ids and the orthogroups selected for each rep
//sims = Channel.from(sim_key_file.collect { it.tokenize( ' ' ) }).map {SIMID, OGS -> [SIMID, OGS]}

//sims.view()

//ce_causal_snps =["c_elegans", "${ce_pop_id}", file("${params.proc_data}/${ce_pop_id}/${ce_pop_id}_0.00.bim.bed.annotated")]

//sim_inputs = sims.map {SIMID, OGS -> [SIMID, OGS, ce_causal_snps]}

// create a new input tuple that includes the simulation id, the orthogroups, and the causal snps
//sim_inputs = sims.combine(anno_snps)


//sim_inputs.view()

//    ],
//    ["c_briggsae",
//    ],
//    ["c_briggsae", file("${params.proc_data}/${cb_pop_id}/${cb_pop}_0.00.bim")],
//    ["c_tropicalis", file("${params.proc_data}/${ct_pop_id}/${ct_pop}_0.00.bim")]
//    ]
