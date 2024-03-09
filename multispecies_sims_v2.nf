#! usr/bin/env nextflow
if( !nextflow.version.matches('>20.0') ) {
    println "This workflow requires Nextflow version 20.0 or greater -- You are running version $nextflow.version"
    println "On QUEST, you can use `module load python/anaconda3.6; source activate /projects/b1059/software/conda_envs/nf20_env`"
    exit 1
}

nextflow.enable.dsl=2

params.bin_dir = "${workflow.projectDir}/bin" // this is different for gcp
params.sparse_cut = 0.05
params.group_qtl = 1000
params.ci_size = 150
sthresh= "BF"
params.maf = 0.05

include {} from './modules/repeated_simulations.nf'

// load the population data from the input folder

params.proc_data 
ce_pop_id = params.ce_pop_id
cb_pop_id = params.cb_pop_id
ct_pop_id = params.ct_pop_id

// load the data to create genotype matrix
pop_geno_inputs = [
    ["c_elegans",
    ${ce_pop_id},
    file("${params.proc_data}/${ce_pop_id}/${ce_pop_id}_0.00.bim"), // the bim file containing all SNPs 
    file("${params.proc_data}/${ce_pop_id}/renamed_chroms.vcf.gz"), // the VCF file thats been filtered to only strains 
    file("${params.proc_data}/${ce_pop_id}/renamed_chroms.vcf.gz.tbi"), // the index file for the VCF
    file("${params.proc_data}/${ce_pop_id}/selected_snps.txt") // the list of SNPs selected as gene based markers
    ]
]

// load the data to simulate phenotypes 
sp_causal_snps_inputs = [
    ["c_elegans",
    ${ce_pop_id},
    file("${params.proc_data}/${ce_pop_id}/${ce_pop_id}_0.00.bim.bed.annotated") // the list of SNPs in the population annotated with gene id and OG status
    ],
    ["c_briggsae",
    ${cb_pop_id},
    file("${params.proc_data}/${cb_pop_id}/${cb_pop_id}_0.00.bim.bed.annotated") // the list of SNPs in the population annotated with gene id and OG status
    ]
]
//load the simulation key files
File sim_key_file = new File("test_data/repeated_sim_keys.txt") ;

// join the simulation data to the input file that specifies simulation ids and the orthogroups selected for each rep
sims = Channel.from(sim_key_file.collect { it.tokenize( ' ' ) }).map {SIMID, OGS -> [SIMID, OGS]}

// for each species in the sp_causl_snps_inputs combine with the sims channel
// to create a channel of inputs for the repeated_simulations module
sim_inputs = sp_causal_snps_inputs.collectMany { species, pop_id, snp_file ->
    sims.map { SIMID, OGS -> [species, pop_id, snp_file, SIMID, OGS] }
}

sim_inputs.view()

//    ],
//    ["c_briggsae",
//    ],
//    ["c_briggsae", file("${params.proc_data}/${cb_pop_id}/${cb_pop}_0.00.bim")],
//    ["c_tropicalis", file("${params.proc_data}/${ct_pop_id}/${ct_pop}_0.00.bim")]
//    ]
