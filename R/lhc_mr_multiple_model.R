#' Main trait pair analysis using LHC-MR
#'
#' @param SP_list List resulting from calculate_SP. Contains the filtered dataset (by every 'SNP_filter'th SNP), the starting points to be used in the pair trait optimisation, the traits' intercepts,
#' the traits' polygenicity if nStep = 2, as well as some extra parameters like the cross-trait intercept and bidirectional causal effect estimated by IVW
#' @param trait.names Vector containing the trait names in the order they were used in merge_sumstats(): Exposure, Outcome
#' @param partition String indicating the partition name to be used for the "rslurm" parallelisation - equivalent to '-p, --partition' in SLURM commands
#' @param account String indicating the account name to be used for the "rslurm" parallelisation - equivalent to '-A, --account' in SLURM commands
#' @param param String indicating which model the likelihood function will be optimised with, either "comp" by default or "U" for a no-confounder model
#' @param paral_method String indicating which method to parallelise the optimisation over the number of sets of starting points. "rslurm" will submit the calculation to a SLURM cluster using
#' a 'Slurm' workload manager, "lapply" will parallelise the optimisation using 'mclapply' over a set number of cores but will go sequentially over the sets of starting points and thus take more time.
#' @param nCores Numerical value indicating number of cores to be used in 'mclapply' to parallelise the analysis. If not set (default value = NA), then it will be calculated as 2/3 of the available cores
#' @param nBlock Numerical value indicating the number of blocks to create from the block jackknife analysis, where at each iteration one block is left out and the optimisation is ran again for a single starting
#' point to obtain eventually 'nBlock' estimates and calculate the SE of the parameter estimates
#' @param M Numerical value indicating the number of SNPs used to calculate the LD reported in the LD file (for genotyped SNPs). Default value = 1e7
#'
#' @return Prints out a summary of the results
#' @export
#'
#' @importFrom rslurm slurm_apply get_slurm_out cleanup_files
#' @importFrom stringr str_detect str_split str_squish str_which
#' @importFrom dplyr mutate
#' @importFrom MASS fitdistr
#' @importFrom mixtools normalmixEM
#' @importFrom stats cov median var
#' @importFrom parallel mclapply detectCores
#'
#' @examples



lhc_mr_multiple_model = function(SP_list,trait.names,partition=NA,account=NA,param="comp",paral_method="rslurm",nCores=NA,nBlock=200,M=1e7){
  
  input.df_filtered = SP_list$input.df_filtered
  SP_matrix = SP_list$sp_mat
  iX = SP_list$iX
  iY = SP_list$iY
  piX = SP_list$piX
  piY = SP_list$piY
  SP_pair = nrow(SP_matrix)
  if(ncol(SP_matrix)==9) {
    nStep = 1
  } else if(ncol(SP_matrix)==7) {
    nStep = 2
  } else{
    cat(print("Input error in number of steps"))
    stop()
  }
  
  if(is.na(nCores)){
    nCores = max(1,floor((parallel::detectCores())/3))
  } else {
    if(nCores > parallel::detectCores()){
      cat(print("Core number chosen is greater than cores available\n"))
      stop()
    }
  }
  
  # sp_h2X   sp_h2Y     sp_tX       sp_tY       sp_axy      sp_ayx   sp_iXY
  
  if(param=="only_causal"){
    SP_matrix = dplyr::select(as.data.frame(SP_matrix),-c(sp_tX,sp_tY))
  }
  if(param=="only_latent"){
    SP_matrix = dplyr::select(as.data.frame(SP_matrix),-c(sp_axy,sp_ayx))
  }
  if(param=="only_a"){
    SP_matrix = dplyr::select(as.data.frame(SP_matrix),-c(sp_tX,sp_tY,sp_ayx))
  }
  if(param=="only_b"){
    SP_matrix = dplyr::select(as.data.frame(SP_matrix),-c(sp_tX,sp_tY,sp_axy))
  }
  if(param=="no_a"){
    SP_matrix = dplyr::select(as.data.frame(SP_matrix),-c(sp_axy))
  }
  if(param=="no_b"){
    SP_matrix = dplyr::select(as.data.frame(SP_matrix),-c(sp_ayx))
  }
  if(param=="NULL"){
    SP_matrix = dplyr::select(as.data.frame(SP_matrix),-c(sp_tX,sp_tY,sp_axy,sp_ayx))
  }
  
  
  # Generate a dataframe of lists for each row of parameters - input for rslurm/lapply
  par.df = data.frame(par=I(apply(SP_matrix,1,as.list)))
  
  # Parameters set-up
  piU=0.1
  bn = 2^7
  bins = 10
  #param="comp" #possibility to develop nesting in later version
  EXP = trait.names[1]
  OUT = trait.names[2]
  
  # Data set up
  nX = mean(input.df_filtered$N.x)  #get sample size for trait X
  nY = mean(input.df_filtered$N.y)  #get sample size for trait Y
  m0 = mean(input.df_filtered$M_LOCAL) #get number of SNPs used to calculate the local LD
  
  bX = input.df_filtered$TSTAT.x/sqrt(nX)   #get standardised beta for trait X
  bY = input.df_filtered$TSTAT.y/sqrt(nY)   #get standardised beta for trait Y
  betXY = cbind(bX,bY)
  
  ld = input.df_filtered$LDSC
  w8s = input.df_filtered$WEIGHT
  pi1 = input.df_filtered$PIK
  sig1 = input.df_filtered$SIGK
  
  # Generate a data frame of positions for block JK - used for SP generation
  uniqPos = seq(1:nrow(input.df_filtered))
  nBlock = nBlock
  nSNP = nrow(betXY) %/% nBlock
  limit = nBlock * nSNP
  leftover = as.numeric(nrow(betXY) %% nSNP )
  start_ind = uniqPos[seq(1, limit, nSNP)]
  end_ind = start_ind+(nSNP-1)
  end_ind[length(end_ind)]=end_ind[length(end_ind)]+leftover #add any leftover SNPs to last block
  
  JK_index=cbind(start_ind,end_ind)
  
  # Run analysis based on parallelisation method//number of steps
  if(nStep == 2){
    if(param=="only_causal"){
      parscale2 = rep(1e-1,5)
    }else if(param=="only_latent"){
      parscale2 = rep(1e-1,5)
    }else if(param=="only_a"){
      parscale2 = rep(1e-1,4)
    }else if(param=="only_b"){
      parscale2 = rep(1e-1,4)
    }else if(param=="no_a"){
      parscale2 = rep(1e-1,6)
    }else if(param=="no_b"){
      parscale2 = rep(1e-1,6)
    }else if(param=="NULL"){
      parscale2 = rep(1e-1,3)
    }else{
      parscale2 = rep(1e-1,7)
    }
    assign(x="betXY", value=betXY, env=.GlobalEnv)
    assign(x="pi1", value=pi1, env=.GlobalEnv)
    assign(x="sig1", value=sig1, env=.GlobalEnv)
    assign(x="w8s", value=w8s, env=.GlobalEnv)
    assign(x="m0", value=m0, env=.GlobalEnv)
    assign(x="M", value=M, env=.GlobalEnv)
    assign(x="nX", value=nX, env=.GlobalEnv)
    assign(x="nY", value=nY, env=.GlobalEnv)
    assign(x="piU", value=piU, env=.GlobalEnv)
    assign(x="piX", value=piX, env=.GlobalEnv)
    assign(x="piY", value=piY, env=.GlobalEnv)
    assign(x="iX", value=iX, env=.GlobalEnv)
    assign(x="iY", value=iY, env=.GlobalEnv)
    assign(x="param", value=param, env=.GlobalEnv)
    assign(x="bn", value=bn, env=.GlobalEnv)
    assign(x="bins", value=bins, env=.GlobalEnv)
    assign(x="parscale2", value=parscale2, env=.GlobalEnv)
    
    # wget https://github.com/LizaDarrous/lhcMR/raw/main/R/pairTrait_twoStep_likelihood.R
    
    if(paral_method=="lapply"){
      cat(print("Running optimisation"))
      test.res <- parallel::mclapply(par.df[[1]], function(x) {
        theta = unlist(x)
        test = optim(theta, pairTrait_twoStep_likelihood_multiple_model,
                     betX=betXY, pi1=pi1, sig1=sig1, w8s=w8s, M=M,
                     m0=m0, nX=nX, nY=nY, pi_U=piU, pi_X=piX, pi_Y=piY, i_X=iX, i_Y=iY,
                     bn=bn, bins=bins, model=param,
                     method = "Nelder-Mead",
                     control = list(maxit = 50,
                                    parscale = parscale2))
        
        list("mLL"=test$value,"par"=test$par,"conv"=test$convergence)
      }, mc.cores = nCores)
      
      
      test.res_tmp <- as.data.frame(t(matrix(unlist(test.res), nrow=length(unlist(test.res[1])))))
      colnames(test.res_tmp)  <- c("mLL",gsub("^sp_","",names(test.res[[1]]$par)),"conv")
      test.res_mx <- data.frame(matrix(0,ncol=9,nrow=NROW(test.res_tmp)))
      colnames(test.res_mx) <- c("mLL", "h2X","h2Y","tX","tY","axy","ayx","iXY","conv")
      test.res_mx[,colnames(test.res_tmp)] <- test.res_tmp
      test.res = test.res_mx
      
      
      test.res %>%
          dplyr::mutate(h2X = abs(h2X),
                        h2Y = abs(h2Y),
                        tX = abs(tX)) -> res_values
      
      res_values = cbind("SP"=c(1:nrow(res_values)),"mLL"=res_values[,1],"piX"=piX,"piY"=piY,res_values[,-1],"iX"=iX,"iY"=iY)
      write.csv(res_values, paste0("FullRes_",EXP,"-",OUT,"__",param,".csv"), row.names = FALSE)
      
      # Running blockJK
      cat(print("Running block JK"))
      res_ordered = res_values[order(res_values$mLL, decreasing = F),]
      
      if(param=="only_causal"){
        sp_mat2 = dplyr::select(res_ordered,,h2X,h2Y,axy,ayx,iXY)
      }
      if(param=="only_latent"){
        sp_mat2 = dplyr::select(res_ordered,,h2X,h2Y,tX,tY,iXY)
      }
      if(param=="only_a"){
        sp_mat2 = dplyr::select(res_ordered,,h2X,h2Y,axy,iXY)
      }
      if(param=="only_b"){
        sp_mat2 = dplyr::select(res_ordered,,h2X,h2Y,ayx,iXY)
      }
      if(param=="no_a"){
        sp_mat2 = dplyr::select(res_ordered,h2X,h2Y,tX,tY,ayx,iXY)
      }
      if(param=="no_b"){
        sp_mat2 = dplyr::select(res_ordered,h2X,h2Y,tX,tY,axy,iXY)
      }
      if(param=="NULL"){
        sp_mat2 = dplyr::select(res_ordered,h2X,h2Y,iXY)
      }
      if(param=="full"){
        sp_mat2 = dplyr::select(res_ordered,h2X,h2Y,tX,tY,axy,ayx,iXY)
      }
      
      sp_mat2 = sp_mat2[1,]
      
      par.df2 = data.frame(par=I(apply(sp_mat2,1,as.list)))
      par.df2 = merge(par.df2,JK_index)
      
      test.res1 <- parallel::mclapply(split(par.df2,1:nrow(par.df2)), function(x) {
        theta = unlist(x[1])
        start_ind = as.numeric(x[2])
        end_ind = as.numeric(x[3])

        test1 = optim(theta, pairTrait_twoStep_likelihood_multiple_model,
                      betXY=betXY[-(start_ind:end_ind),], pi1=pi1[-(start_ind:end_ind)], sig1=sig1[-(start_ind:end_ind)],
                      w8s=w8s[-(start_ind:end_ind)], M=M, pi_U=piU,
                      pi_X=piX, pi_Y=piY, i_X=iX, i_Y=iY,
                      m0=m0, nX=nX, nY=nY, bn=bn, bins=bins, model=param,
                      method = "Nelder-Mead",
                      control = list(maxit = 10,
                                     parscale = parscale2))
        
        list("mLL"=test1$value,"par"=test1$par,"conv"=test1$convergence,"start_ind"=start_ind, "end_ind"=end_ind)
      }, mc.cores = nCores)
      
      test.res1_tmp <- as.data.frame(t(matrix(unlist(test.res1), nrow=length(unlist(test.res1[1])))))
      colnames(test.res1_tmp)  <- c("mLL",gsub("^par.[0-9]+.","",names(test.res1[[1]]$par)),"conv","start_ind","end_ind")
      test.res1_mx <- data.frame(matrix(0,ncol=11,nrow=NROW(test.res1_tmp)))
      colnames(test.res1_mx) <- c("mLL", "h2X","h2Y","tX","tY","axy","ayx","iXY","conv","start_ind", "end_ind")
      test.res1_mx[,colnames(test.res1_tmp)] <- test.res1_tmp
      test.res1 = test.res1_mx

      test.res1 %>%
        dplyr::mutate(h2X = abs(h2X),
                      h2Y = abs(h2Y),
                      tX = abs(tX)) -> res_values2
      
      write.csv(res_values2, paste0("FullRes_",EXP,"-",OUT,"__",param,"_",nBlock,"blockJK.csv"), row.names = FALSE)
      
    }
    
  }
  model <- param
  # BlockJK analysis regardless of method - one step or two step
  res_min2 = res_values2
  
  res_minFil = dplyr::select(as.data.frame(res_min2), -c(mLL,conv,start_ind,end_ind))
  JK_res = as.data.frame(matrix(data=NA, nrow=ncol(res_minFil),ncol=17))
  colnames(JK_res)=c("Parameter","mean","median","se","se_JK","var","loglik_1comp","AIC_1comp","convergance_2comp","mu1","mu2","sigma1","sigma2","lambda1","lambda2","loglik_2comp","AIC_2comp")
  JK_res$Parameter=colnames(res_minFil)
  JK_res$mean = colMeans(res_minFil)
  JK_res$median = apply(res_minFil, 2, median)
  JK_res$se = apply(res_minFil, 2, sd)
  JK_res$se_JK = (apply(res_minFil, 2, sd))*sqrt(nBlock-1)#*sqrt(10)
  JK_res$var = apply(res_minFil, 2, var)
  
  for (x in c(1:ncol(res_minFil))){
    param = res_minFil[,x]
    #print(colnames(res_minFil)[x])
    Xf1 = MASS::fitdistr(param, "normal")
    Xf2 = tryCatch({capture.output(mixtools::normalmixEM(param, k=2, maxit=1e8))}, warning = function(warning_condition) {
      return(NA)
    }, error = function(error_condition) {
      print("Error: 2-component normal mixture could not be fit on this data.")
      return(NA)
    })
    AIC1 = 2*2 - 2*(Xf1$loglik)
    AIC2 = 2*4 - 2*(as.numeric(str_split(str_squish(Xf2[str_which(Xf2, "loglik")[1]+1]), " " )[[1]][2]))
    JK_res$loglik_1comp[x] = Xf1$loglik
    JK_res$loglik_2comp[x] = as.numeric(str_split(str_squish(Xf2[str_which(Xf2, "loglik")[1]+1]), " " )[[1]][2])
    JK_res$AIC_1comp[x] = AIC1
    JK_res$AIC_2comp[x] = AIC2
    
    
    if(!is.na(Xf2[1])){
      #print(colnames(res_minFil)[x])
      JK_res$convergance_2comp[x] = !any(str_detect(Xf2, "Warning: not convergent!"))
      JK_res[x,10:11] = as.numeric(str_split(str_squish(Xf2[str_which(Xf2, "mu")+1]), " " )[[1]][2:3])
      JK_res[x,12:13] = (as.numeric(str_split(str_squish(Xf2[str_which(Xf2, "sigma")+1]), " " )[[1]][2:3]))
      JK_res[x,14:15] = as.numeric(str_split(str_squish(Xf2[str_which(Xf2, "lambda")+1]), " " )[[1]][2:3])
    }
  }
  
  tstat = (JK_res$mu1 - JK_res$mu2) / sqrt(JK_res$sigma1^2 + JK_res$sigma2^2)
  t.pval = 2*pnorm(-abs(tstat))
  
  JK_res$tstat = tstat
  JK_res$tstat_pval = t.pval
  
  JK_res$ci_lower_JK = JK_res$mean - (1.96*JK_res$se_JK)
  JK_res$ci_upper_JK = JK_res$mean + (1.96*JK_res$se_JK)
  
  bimo = which(JK_res$tstat_pval == 0)
  
  JK_res$bimod = "FALSE"
  JK_res$bimod[bimo] = "TRUE"
  
  write.csv(JK_res,paste0(nBlock,"blockJK_summary_",EXP,"-",OUT,"__",model,".csv"), row.names = F)
  
  cov_matrix = cov(res_minFil)
  write.csv(cov_matrix,paste0(nBlock,"blockJK_VarCovMatrix_",EXP,"-",OUT,"__",model,".csv"))
  print("Done!")
  
  # Print out results in table format
  res_est = res_values[which(res_values$mLL==min(res_values$mLL)),]
  res_est = unlist(dplyr::select(res_est, -c(SP,mLL,conv)))
  res_JKse = rep(NA,length(res_est))
  if(nStep == 1){res_JKse[1:length(JK_res$se_JK)] = JK_res$se_JK}
  if(nStep == 2){res_JKse[3:(2+length(JK_res$se_JK))] = JK_res$se_JK}
  res_pval = 2*pnorm(-abs(res_est/res_JKse))
  res_tab = rbind(res_est,res_JKse,res_pval)
  rownames(res_tab)=c("Parameter estimates","SE-JK","Pval")
  
  write.csv(res_tab,paste0("SummarisedResults_",EXP,"-",OUT,"__",model,".csv"),row.names = T)
  print(res_tab)
  
}

