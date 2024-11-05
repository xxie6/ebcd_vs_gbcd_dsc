# simulate samples from 4 population balanced tree
baltree_4pop: 4pop_balancedtree_sim.R
    # module input and variables
    pop_sizes: rep(40, 4)
    n_genes: 1000
    branch_sds: rep(2, 7)
    indiv_sd: 1
    seed: R{1:10}
    # module decoration
    @ALIAS: args = List()
    # module output
    $data: data
    $true_L: data$LL
    $dat: data$YYt

# simulate samples from 10 group overlapping binary structure
group_overlap: group_overlapping_sim.R
    # module input and variables
    n: 100
    p: 1000
    k: 10
    group_sd: 1 
    indiv_sd: 1
    seed: R{1:10}
    pi0: 0.1
    # module decoration
    @ALIAS: args = List()
    # module output
    $data: data
    $true_L: data$LL
    $dat: data$YYt

# simulate samples from non-overlapping binary structure
group_nonoverlap: group_nonoverlapping_sim.R
    # module input and variables
    pop_sizes: rep(40, 4), c(20,50,30,60)
    n_genes: 1000
    branch_sds: rep(2, 4)
    indiv_sd: 1
    seed: R{1:10}
    # module decoration
    @ALIAS: args = List()
    # module output
    $data: data
    $true_L: data$LL
    $dat: data$YYt

pca: runpca.R
    # module input and variables
    input: $data
    # module output
    $fit_obj: pca_data$pca_fit
    $est_L: pca_data$est_L
    $est_LLt: pca_data$est_LLt

ebcd: runebcd.R
    # module input and variables
    seed: 1
    ebnm_fn: ebnm::ebnm_generalized_binary
    # module decoration
    @ALIAS: args = List()
    input: $data
    # module output
    $fit_obj: ebcd_data$ebcd_fit
    $est_L: ebcd_data$ebcd_fit$EL
    $est_LLt: ebcd_data$est_LLt
    
ebmfcov: runebmfcov.R
    # module input and variables
    ebnm_fn: ebnm::ebnm_generalized_binary
    @ALIAS: args = List()
    input: $data
    # module output
    $fit_obj: flash_cov_data$flash_cov_fit
    $est_L: flash_cov_data$flash_cov_fit$L_pm
    $est_LLt: flash_cov_data$est_LLt

gbcd: rungbcd.R
    # module input and variables
    input: $data
    # module output
    $fit_obj: gbcd_data$fit_obj
    $est_L: gbcd_data$fit_obj$L
    $est_LLt: gbcd_data$est_LLt

#drift: rundrift.R
    # module input and variables
    #input: $data
    # module output
    #$fit_obj: drift_data$drift_fit
    #$est_L: drift_data$drift_fit$EL

sindclus: runsindclus.R
    # module input and variables
    init_method: "default"
    additive_term: TRUE, FALSE
    off_diagonal: FALSE
    # module decoration
    @ALIAS: args = List()
    input: $data
    # module output
    $fit_obj: sindclus_data$sindclus_fit
    $est_L: sindclus_data$sindclus_fit$P
    $est_LLt: sindclus_data$est_LLt

sympres: runsympres.R
    # module input and variables
    init_method: "default"
    additive_term: TRUE, FALSE
    off_diagonal: FALSE
    # module decoration
    @ALIAS: args = List()
    input: $data
    # module output
    $fit_obj: sympres_data$sympres_fit
    $est_L: sympres_data$sympres_fit$P
    $est_LLt: sympres_data$est_LLt

crossprod_similarity: score_crossprod_sim.R
    # module input and variables
    est: $est_L
    truth: $true_L
    # module output
    $result: result

cov_L2_fit: score_L2_fit.R
    # module input and variables
    est: $est_LLt
    dat: $dat
    # module output
    $result: result

DSC:
    # module ensembles
    define:
      simulate: baltree_4pop, group_overlap, group_nonoverlap
      analyze1: pca, ebcd, ebmfcov, gbcd, sindclus, sympres
      #analyze2: drift
      score: crossprod_similarity, cov_L2_fit
    # pipelines
    run: 
      pipeline_1: simulate * analyze1 * score
      #pipeline_2: simulate * drift * crossprod_similarity
    output: cov_decomps_comparison_v2

