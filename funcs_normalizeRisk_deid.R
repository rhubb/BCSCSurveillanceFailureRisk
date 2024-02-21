### =========================================================================================================== ###
### A function to normalize risk across patient sites and PBC diagnosis year
### =========================================================================================================== ###

normalizeRisk_siteDxyear = function(model_object, data, weights=NULL, weight_level="response", new_lambda=NULL, designMx_siteDxyear){
  
  if(is.null(new_lambda)) lambda_choice = "lambda.min"
  else lambda_choice = new_lambda

  ## Get needed combination of site and PBC dxyear based on the names of weight:
  if(is.null(dim(weights)) & length(weights)>0) combs = names(weights)
  else combs = colnames(weights)
  risk_pr_combs = sapply(1:length(combs),function(comb){
    library(dplyr)
    designMx_comb = designMx_siteDxyear[comb,]
    if(is.matrix(data)) data=as_data_frame(data)
    data_with_comb = data %>% mutate(PatientSite_cB=designMx_comb["PatientSite_cB"],
                                     PatientSite_cC=designMx_comb["PatientSite_cC"],
                                     PatientSite_cD=designMx_comb["PatientSite_cD"],
                                     PatientSite_cE=designMx_comb["PatientSite_cE"],
                                     PatientSite_cF=designMx_comb["PatientSite_cF"],
                                     PatientSite_cG=designMx_comb["PatientSite_cG"],
                                     dxyear_ns.1=designMx_comb["dxyear_ns.1"],
                                     dxyear_ns.2=designMx_comb["dxyear_ns.2"],
                                     dxyear_ns.3=designMx_comb["dxyear_ns.3"],
                                     dxyear_ns.4=designMx_comb["dxyear_ns.4"])
    data_with_comb = as.matrix(data_with_comb)
    if(weight_level=="link"){
      risk_pr_comb = data_with_comb %*% (model_object[colnames(data_with_comb),1]) + model_object["(Intercept)",1]
    }
    else{
      ### Define expit function
      expit = function(x) exp(x)/(1+exp(x))
      risk_pr_comb = expit(data_with_comb %*% model_object[colnames(data_with_comb),1] + model_object["(Intercept)",1])
    }
    return(risk_pr_comb)
  })

  if(is.null(weights)) risk_pr_norm = rowMeans(risk_pr_combs)
  else if(is.null(dim(weights)) & length(weights)>0) risk_pr_norm = matrixStats::rowWeightedMeans(x = risk_pr_combs,w = weights)
  else risk_pr_norm = rowSums(risk_pr_combs * weights) / rowSums(weights)
  
  if(weight_level=="link") risk_pr_norm = expit(risk_pr_norm)
  return(risk_pr_norm)
}