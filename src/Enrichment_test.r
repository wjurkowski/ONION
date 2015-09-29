#####################################################################################################
# Function for FDR corrected hypergeometric enrichment test
#####################################################################################################

## Arguments of HyperGeomFDR function:
#	steps:		the rounds of simulations (a single number)
#	pool:		background genes (character vector)
#	select:		genes to investigate (character vector)
#	DB: 		the genes set used for enrichment analysis (character list)
#	nthreads:	number of threads to use (a single number)

## Description of the hypergeometric test
# phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
# Arguments:
#        q: vector of quantiles representing the number of white balls
#           drawn without replacement from an urn which contains both
#           black and white balls.
#        m: the number of white balls in the urn.
#        n: the number of black balls in the urn.
#        k: the number of balls drawn from the urn.
# 
# x=length(intersect(select,DB_i))    	#Number of common genes between DB and select
# m=length(intersect(pool,DB_i))        #Number of common genes between DB and pool
# n=length(pool)-length(intersect(pool,DB_i))     #Number of non-pool genes among DB (setdiff)
# k=length(select)                    	#Number of genes in select
# P_val=dhyper(length(intersect(select,DB_i)), length(intersect(pool,DB_i)), length(pool)-length(intersect(pool,DB_i)), length(select))
# 
# P_val=(choose(length(intersect(pool,DB_i)),length(intersect(select,DB_i)))*
#     choose(length(pool)-length(intersect(pool,DB_i)),length(select)-length(intersect(select,DB_i))))/
#     choose(length(pool),length(select))

HyperGeomFDR=function(steps, pool, select, DB, nthreads=4) {
    DB_names=names(DB)
    num_DB=length(DB)
    size_pool=length(pool)
    size_select=length(select)
    
    DB_in_select=integer(num_DB)
    DB_in_pool=integer(num_DB)
    Genes_in_DB=integer(num_DB)
    P_val=double(num_DB)
    R_obs=integer(num_DB)
    
    # for every DB entity in the DB list
    for (i in 1:num_DB) {
        # create a vector of genes connected to the ith DB category
        DB_i=DB[[i]]
        # hypergometric test
        DB_in_select[i]=length(intersect(select,DB_i))	#q: number of common genes between a DBterm and select
        DB_in_pool[i]=length(intersect(pool,DB_i))	#m: number of common genes between DBterm and BackGround
        Genes_in_DB[i]=length(DB_i)			
							#n:  number of non-pool genes among DB
							#k: number of genes in select
	# phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
        P_val[i]=1-phyper(DB_in_select[i]-1, DB_in_pool[i], size_pool-DB_in_pool[i], size_select)
    }
	P_val_round=round(P_val, digits=15) ## can change the digits, this is important for the precision of '0' is R
    for (i in 1:num_DB) {
    	 R_obs[i]=sum(P_val_round<=P_val_round[i])
    }
    P_val_df=data.frame(DB_names, DB_in_select, DB_in_pool, Genes_in_DB, P=P_val, P_adj_Bonf=p.adjust(P_val, method="bonferroni"), P_adj_BH=p.adjust(P_val, method="BH"), R_obs)
    
    ######
    # simualtion
    ######
    R_exp=integer(num_DB)
    # random sampling from pool (background genes)
    # The time consuming step. The simulation here can be parallelized
   	require(snow)
	require(rlecuyer)
	seeds=sample(seq(1e4,1e6),6) # max number of seeds for RNGstream is 6
	cl=makeCluster(nthreads, type="SOCK")
	clusterSetupRNG(cl, type='RNGstream', seed=seeds)
	P_Sim_vec=clusterApply(cl, rep(ceiling(steps/nthreads), nthreads), sim_hyperGeom, pool, select,DB) # return a list
	stopCluster(cl)
	P_Sim_vec=as.vector(unlist(P_Sim_vec))
	for (l in 1:length(P_val_df$P)) {
	    P_Sim_round=round(P_Sim_vec, digits=15)
	    R_exp[l]=sum(P_Sim_round<=P_val_round[l])
	}
    P_val_df$R_exp=R_exp/steps
    P_val_df$FDR=P_val_df$R_exp/R_obs
    return(P_val_df)
}

sim_hyperGeom=function(steps, pool, select, DB) {
    DB_names=names(DB)
    num_DB=length(DB_names)
    P_Sim_mat=matrix(numeric(num_DB*steps), ncol=steps)
    size_pool=length(pool)
    size_select=length(select)

    for (j in 1:steps) {
        Rand.select=sample(pool, size_select)       
        for (i in 1:num_DB) {
            DB_i=DB[[i]]
            # hypergometric test
            P_Sim_mat[i,j]=1-phyper(length(intersect(Rand.select, DB_i))-1, length(intersect(pool, DB_i)), size_pool-length(intersect(pool, DB_i)), size_select)
        }
    }
    return(as.vector(P_Sim_mat))
}


