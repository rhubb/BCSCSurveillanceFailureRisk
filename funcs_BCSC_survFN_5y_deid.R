BCSC_survFN_5y = function(data, 
                          var_name = c("race_ethnicity" = "racenci_c",
                                       "family_history" = "famhx_c",
                                       "breast_density" = "density_c",
                                       "PBC_stage" = "stage1",
                                       "PBC_grade" = "grade1_c",
                                       "PBC_ERPR" = "erplus1",
                                       "PBC_adj_therapy" = "adjthrpy_1",
                                       "menopause" = "menopht_2cat",
                                       "mammo_type" = "mammtype",
                                       "PBC_MOD" = "mod1_c",
                                       "PBC_histology" = "hist1_c",
                                       "PBC_surgery_type" = "srg1_new",
                                       "PBC_radiation" = "rt_new",
                                       "PBC_dx_age" = "dxage",
                                       "BMI" = "bmi_c",
                                       "surveillance_interval" = "mossinceprvsurvmm_cat",
                                       "time_since_prv_mammo" = "prvmos_c",
                                       "exam_age" = "age_c",
                                       "month_since_PBC_dx" = "mossince1stdx"
                                       ),
                          path)
# "race_ethnicity" = "racenci_c", # 1: White, 2: Black, 3: Asian, 4: Pacific Islander, 5: Native American, 6: Mixed, 7: Others, 8: Hispanic
# "family_history" = "famhx_c", # 0: No, 1: Yes
# "breast_density" = "density_c", # 1: A, 2: B, 3: C, 4: D
# "PBC_stage" = "stage1", # 1: DCIS, 2: I NOS - IB, 3: II NOS and IIA, 4: IIB+
# "PBC_grade" = "grade1_c", # 1: grade 1, 2: grade 2, 3: grade 3
# "PBC_ERPR" = "erplus1", # 1: ER negative and PR negative, 2: ER negative and PR positive, 3: ER positive and PR negative, 4: ER positive and PR positive
# "PBC_adj_therapy" = "adjthrpy_1", # 0: None, 1: Chemotherapy only, 2: Hormonal therapy only, 3: Both Chemo amd hormonal therapy
# "menopause" = "menopht_2cat", # 0: Pre/peri-menopausal, 1: Postmenopausal
# "mammo_type" = "mammtype", # 1: film, 2: DM, 5: DBT
# "PBC_MOD" = "mod1_c", # 1: screening detected, 2: intervally detected, 3: clinically detected
# "PBC_histology" = "hist1_c", # 0: DCIS, 1: Invasive NOS, 2: Invasive ductal, 3: Invasive lobular, 4: Invasive mixed
# "PBC_surgery_type" = "srg1_new", # 1: mastectomy, 2: BCS
# "PBC_radiation" = "rt_new", # 1: with RT, 2: without RT
# "PBC_dx_age" = "dxage", # continuous value between 21 and 103
# "BMI" = "bmi_c", # continuous value between 15.05 and 89.43
# "surveillance_interval" = "mossinceprvsurvmm_cat", # 1: 9-14 months, 2: 1st surveillance, 3: 3-8 months, 4: 15-23 months, 5: 24+ months
# "time_since_prv_mammo" = "prvmos_c", # continuous value between 3 and 254. Values outside the window will be forced to either take the min or max of this range.
# "exam_age" = "age_c" # continuous value between 23 and 104
# "month_since_PBC_dx" = "mossince1stdx" # continuous value not smaller than 6


{
  library(glmnet)
  library(dplyr)
  library(splines)

  ### === Load all needed information === ###
  load(paste0(path,"/5yr_model_needed_materials_deid.RData"))
  
  ## Materials needed for calculating the values of the splines and the normalization process.
  log_prvmos_ns_at12 = ns(x=log(12), knots=ns_fit_parameter$log_prvmos_ns_fit$knots, 
                          Boundary.knots=ns_fit_parameter$log_prvmos_ns_fit$Boundary.knots, 
                          intercept=ns_fit_parameter$log_prvmos_ns_fit$intercept.ind)
  log_prvmos_ns.1_at12 = log_prvmos_ns_at12[1]
  log_prvmos_ns.2_at12 = log_prvmos_ns_at12[2]
  log_prvmos_ns.3_at12 = log_prvmos_ns_at12[3]
  
  ## Define expit function
  expit = function(x) exp(x)/(1+exp(x))
  
  ### === Process the input data to match what the model would take for generating prediction === ###
  
  ## Ensure the input data is a data.frame or tibble
  if(all(!class(data) %in% c("matrix", "numeric", "data.frame", "data_frame"))) stop("Input data need to be one of the following: data.frame, matrix or tibble")
  else{
    if(class(data)=="matrix") data_df = data.frame(data, stringsAsFactors = FALSE)
    else data_df = data
  }
  ## Replace the current predictor names with what match to the model coefficients
  library(dplyr)
  
  if(nrow(data_df)==1) data_df = bind_rows(data_df,data_df)
  
  data_df_rename_tmp = data_df %>% select(all_of(var_name))
  
  var_name_model = c("racenci_c" = "race_ethnicity", 
                     "famhx_c" = "family_history", 
                     "density_c" = "breast_density", 
                     "stage1" = "PBC_stage", 
                     "grade1_c" = "PBC_grade", 
                     "erplus1" = "PBC_ERPR", 
                     "adjthrpy_1" = "PBC_adj_therapy", 
                     "menopht_2cat" = "menopause", 
                     "mammtype" = "mammo_type", 
                     "mod1_c" = "PBC_MOD", 
                     "hist1_c" = "PBC_histology", 
                     "srg1_new" = "PBC_surgery_type", 
                     "rt_new" = "PBC_radiation",
                     "dxage" = "PBC_dx_age",
                     "bmi_c" = "BMI",
                     "mossinceprvsurvmm_cat" = "surveillance_interval",
                     "prvmos_c" = "time_since_prv_mammo",
                     "age_c" = "exam_age",
                     "mossince1stdx" = "month_since_PBC_dx")
  
  data_df_rename = data_df_rename_tmp %>% select(all_of(var_name_model))
  
  ## Data checking before compute needed variables for the model: check the range of continuous predictors to match that in the data used to derive splines.
  
  if(min(data_df_rename$prvmos_c)<3 | max(data_df_rename$prvmos_c)>254) warning("The predictor 'months since previous mammography' needs to be at least 3 months and no larger than 254 months. Any value of this predictor less than 3 months are replaced by 3 and larger than 254 months by 254.")
  if(min(data_df_rename$age_c)<23 | max(data_df_rename$age_c)>104) warning("The predictor 'age at surveillance exam' needs to be at least 23 years and no larger than 104 years. Any value of this predictor less than 23 years are replaced by 23 and larger than 104 years by 104.")
  if(min(data_df_rename$dxage)<21 | max(data_df_rename$dxage)>103) warning("The predictor 'age at PBC diagnosis' needs to be at least 21 years and no larger than 103 years. Any value of this predictor less than 21 years are replaced by 21 and larger than 103 years by 103.")
  if(min(data_df_rename$bmi_c)<15.05 | max(data_df_rename$bmi_c)>89.43) warning("The predictor 'BMI' needs to be at least 15.05kg/m2 and no larger than 89.43kg/m2. Any value of this predictor less than 15.05 are replaced by 15.05 and larger than 89.43 by 89.43.")
  if(min(data_df_rename$mossince1stdx)<6) warning("The predictor 'month since PBC diagnosis' needs to be at least 6 months. Any value of this predictor less than 6 are replaced by 6.")
  
  ## Derive needed variables for the risk model:
  
  data_df_new = data_df_rename %>%  
    mutate(racenci_new = factor(racenci_c, levels=c(1,2,3,4,8,5,6,7), labels=c("White","Black","Asian","Asian","Others","Others","Others","Hispanic")),
           famhx_c = factor(famhx_c, levels=c("0","1"), labels=c("no","yes")),
           density_c = factor(density_c, levels=c("1","2","3","4")),
           stage1 = factor(stage1, levels=c("2","1","3","4"), labels = c("I NOS - IB", "DCIS", "II NOS and IIA", "IIB+")),
           # stage1 = relevel(stage1, ref="I NOS - IB"),
           grade1_c = factor(grade1_c, levels=c("1","2","3")),
           erplus1 = factor(erplus1, levels=c("1","2","3","4"), labels = c("ERnPRn","ERnPRp","ERpPRn","ERpPRp")),
           adjthrpy_1 = factor(adjthrpy_1, levels=c("0","1","2","3"), labels=c("none","chemoOnly","hormonalOnly","both")),
           menopht_2cat = factor(menopht_2cat, levels=c("1","0"), labels=c("Post","Pre")),         
           mammtype_new = factor(mammtype, levels = c("2","1","5"), labels = c("mamm_2D","film","DBT")),
           mod1_new = factor(mod1_c, levels=c("1","2","3"), labels = c("screening_detect","interval","clinical")),
           hist1_new = factor(hist1_c, levels=c("0","1","2","3","4"), labels = c("ductal","ductal","ductal","non-ductal","non-ductal")),
           srg1_new = factor(srg1_new, levels=c("1","2"), labels = c("mastectomy","BCS")),
           rt_new = factor(rt_new, levels=c("1","2"), labels = c("with RT", "without RT")),
           dxage = pmin(pmax(dxage,21),103),
           dxage_ns.1 = ns(x=dxage, knots=ns_fit_parameter$dxage_ns_fit$knots,
                           Boundary.knots = ns_fit_parameter$dxage_ns_fit$Boundary.knots,
                           intercept = ns_fit_parameter$dxage_ns_fit$intercept.ind)[,1],
           dxage_ns.2 = ns(x=dxage, knots=ns_fit_parameter$dxage_ns_fit$knots,
                           Boundary.knots = ns_fit_parameter$dxage_ns_fit$Boundary.knots,
                           intercept = ns_fit_parameter$dxage_ns_fit$intercept.ind)[,2],
           dxage_ns.3 = ns(x=dxage, knots=ns_fit_parameter$dxage_ns_fit$knots,
                           Boundary.knots = ns_fit_parameter$dxage_ns_fit$Boundary.knots,
                           intercept = ns_fit_parameter$dxage_ns_fit$intercept.ind)[,3],
           dxage_ns.4 = ns(x=dxage, knots=ns_fit_parameter$dxage_ns_fit$knots,
                           Boundary.knots = ns_fit_parameter$dxage_ns_fit$Boundary.knots,
                           intercept = ns_fit_parameter$dxage_ns_fit$intercept.ind)[,4],
           bmi_c = pmin(pmax(bmi_c,15.05),89.43),
           bmi_ns.1 = ns(x=bmi_c, knots=ns_fit_parameter$bmi_ns_fit$knots,
                         Boundary.knots = ns_fit_parameter$bmi_ns_fit$Boundary.knots,
                         intercept = ns_fit_parameter$bmi_ns_fit$intercept.ind)[,1],
           bmi_ns.2 = ns(x=bmi_c, knots=ns_fit_parameter$bmi_ns_fit$knots,
                         Boundary.knots = ns_fit_parameter$bmi_ns_fit$Boundary.knots,
                         intercept = ns_fit_parameter$bmi_ns_fit$intercept.ind)[,2],
           bmi_ns.3 = ns(x=bmi_c, knots=ns_fit_parameter$bmi_ns_fit$knots,
                         Boundary.knots = ns_fit_parameter$bmi_ns_fit$Boundary.knots,
                         intercept = ns_fit_parameter$bmi_ns_fit$intercept.ind)[,3],
           bmi_ns.4 = ns(x=bmi_c, knots=ns_fit_parameter$bmi_ns_fit$knots,
                         Boundary.knots = ns_fit_parameter$bmi_ns_fit$Boundary.knots,
                         intercept = ns_fit_parameter$bmi_ns_fit$intercept.ind)[,4],
           mossinceprvsurvmm_cat = factor(mossinceprvsurvmm_cat, levels=c("1","2","3","4","5"), 
                                          labels=c("9-14","1st surveillance","3-8","15-23","24+")),
           prvmos_c = pmin(pmax(prvmos_c,3),254),
           log_prvmos = log(prvmos_c),
           log_prvmos_ns.1 = ns(x=log_prvmos, knots = ns_fit_parameter$log_prvmos_ns_fit$knots,
                                Boundary.knots = ns_fit_parameter$log_prvmos_ns_fit$Boundary.knots,
                                intercept = ns_fit_parameter$log_prvmos_ns_fit$intercept.ind)[,1],
           log_prvmos_ns.2 = ns(x=log_prvmos, knots = ns_fit_parameter$log_prvmos_ns_fit$knots,
                                Boundary.knots = ns_fit_parameter$log_prvmos_ns_fit$Boundary.knots,
                                intercept = ns_fit_parameter$log_prvmos_ns_fit$intercept.ind)[,2],
           log_prvmos_ns.3 = ns(x=log_prvmos, knots = ns_fit_parameter$log_prvmos_ns_fit$knots,
                                Boundary.knots = ns_fit_parameter$log_prvmos_ns_fit$Boundary.knots,
                                intercept = ns_fit_parameter$log_prvmos_ns_fit$intercept.ind)[,3],
           age_c = pmin(pmax(age_c,23),104),
           age_ns.1 = ns(x=age_c, knots = ns_fit_parameter$age_ns_fit$knots,
                         Boundary.knots = ns_fit_parameter$age_ns_fit$Boundary.knots,
                         intercept = ns_fit_parameter$age_ns_fit$intercept.ind)[,1],
           age_ns.2 = ns(x=age_c, knots = ns_fit_parameter$age_ns_fit$knots,
                         Boundary.knots = ns_fit_parameter$age_ns_fit$Boundary.knots,
                         intercept = ns_fit_parameter$age_ns_fit$intercept.ind)[,2],
           age_ns.3 = ns(x=age_c, knots = ns_fit_parameter$age_ns_fit$knots,
                         Boundary.knots = ns_fit_parameter$age_ns_fit$Boundary.knots,
                         intercept = ns_fit_parameter$age_ns_fit$intercept.ind)[,3],
           age_ns.4 = ns(x=age_c, knots = ns_fit_parameter$age_ns_fit$knots,
                         Boundary.knots = ns_fit_parameter$age_ns_fit$Boundary.knots,
                         intercept = ns_fit_parameter$age_ns_fit$intercept.ind)[,4],
           yrsincdx_cat = cut(mossince1stdx,
                              breaks=c(6,11,35,59,83,120,1200),
                              labels=c("<1","1-2","3-4","5-6","7-9",">10"),
                              include.lowest = TRUE),
           yrsincdx_cat_new = relevel(yrsincdx_cat,ref="1-2")
    )

  ### === Calculate predicted risk with 30 LASSO models === ###
  
  risk_pr_30lasso = sapply(1:30, simplify=FALSE,function(mi_round){
    print(paste0("LASSO model No. ", mi_round))
    imp_act_tmp = data_df_new
    w = do.call(rbind,weight_list_deid[imp_act_tmp$racenci_new])
    
    #### Assign 0 to create dxyear_ns variables and "A" to create PatientSite_c (factor)
    imp_act_tmp = imp_act_tmp %>% mutate(
      dxyear_ns.1 = 0,
      dxyear_ns.2 = 0,
      dxyear_ns.3 = 0,
      dxyear_ns.4 = 0,
      PatientSite_c = factor("A", levels = c("A","B","C","D","E","F","G"))
    )
    
    #### Change the reference group of stage1 to match to that in the LASSO model:
    if(levels(imp_act_tmp$stage1)[1]!="I NOS - IB") imp_act_tmp$stage1 = relevel(imp_act_tmp$stage1, ref="I NOS - IB")
    
    #### Define model formula for LASSO
    ## Model fitting and prediction on actively imputed data with transform-impute-tranform (TIT) approach for interaction terms
    formula_glmnet = as.formula("~PatientSite_c + racenci_new + famhx_c + density_c + stage1 + 
    grade1_c + erplus1 + adjthrpy_1 + menopht_2cat + mammtype_new + 
    mossinceprvsurvmm_cat + mod1_new + hist1_new + srg1_new + 
    rt_new + yrsincdx_cat_new + age_ns.1 + age_ns.2 + age_ns.3 + 
    age_ns.4 + dxage_ns.1 + dxage_ns.2 + dxage_ns.3 + dxage_ns.4 + 
    dxyear_ns.1 + dxyear_ns.2 + dxyear_ns.3 + dxyear_ns.4 + log_prvmos_ns.1 + 
    log_prvmos_ns.2 + log_prvmos_ns.3 + bmi_ns.1 + bmi_ns.2 + 
    bmi_ns.3 + bmi_ns.4 + racenci_new:bmi_ns.1 + racenci_new:bmi_ns.2 + 
    racenci_new:bmi_ns.3 + racenci_new:bmi_ns.4 + density_c:bmi_ns.1 + 
    density_c:bmi_ns.2 + density_c:bmi_ns.3 + density_c:bmi_ns.4 + 
    density_c:racenci_new + mossinceprvsurvmm_cat:racenci_new")
    
    #### Get 1yr prediction at index exam
    library(glmnet)
    source(paste0(here(),"/funcs_normalizeRisk_deid.R"))
    input_x_index = model.matrix(formula_glmnet,data=imp_act_tmp)[,-1]
    risk_pr_1yr_index_w = normalizeRisk_siteDxyear(model_object = coef_lasso_primary_deid[[mi_round]], 
                                                   data = input_x_index, 
                                                   weights=w,
                                                   designMx_siteDxyear = uniq_comb_glmnet_input_deid)
    risk_pr_1yr_compete_index_w = normalizeRisk_siteDxyear(model_object = coef_lasso_compete_deid[[mi_round]],
                                                           data = input_x_index,
                                                           weights=w,
                                                           designMx_siteDxyear = uniq_comb_glmnet_input_deid)
    print("risk at index exam done!")
    
    #### Get 1yr prediction at 4 follow-up exams
    ## Assuming an annual surveillance schedule
    pred_constant = c("PatientSite_c", "racenci_new", "famhx_c", "density_c", "stage1",
                      "grade1_c", "erplus1", "adjthrpy_1", "menopht_2cat", "mammtype_new",
                      "mod1_new", "hist1_new", "srg1_new", "rt_new",
                      "dxage_ns.1", "dxage_ns.2", "dxage_ns.3", "dxage_ns.4",
                      "dxyear_ns.1", "dxyear_ns.2", "dxyear_ns.3", "dxyear_ns.4",
                      "bmi_ns.1", "bmi_ns.2", "bmi_ns.3", "bmi_ns.4")
    pred_assign = c("mossinceprvsurvmm_cat","log_prvmos_ns.1", "log_prvmos_ns.2", "log_prvmos_ns.3")
    pred_varying = c("age_ns.1", "age_ns.2", "age_ns.3", "age_ns.4", "yrsincdx_cat_new")
    ## Create input_x_constant, input_x_assign which are invariant across 4 follow-up exams
    input_x_constant = imp_act_tmp[,colnames(imp_act_tmp)%in%pred_constant]
    input_x_assign = imp_act_tmp[,colnames(imp_act_tmp)%in%pred_assign]
    input_x_assign$mossinceprvsurvmm_cat[input_x_assign$mossinceprvsurvmm_cat!="9-14"] = "9-14"
    input_x_assign[,"log_prvmos_ns.1"] = log_prvmos_ns.1_at12
    input_x_assign[,"log_prvmos_ns.2"] = log_prvmos_ns.2_at12
    input_x_assign[,"log_prvmos_ns.3"] = log_prvmos_ns.3_at12
    
    ## Create input_x_varying which is varying across 4 follow-up exams as a list
    risk_pr_1yr_4fu_w = sapply(1:4,simplify=FALSE,function(fu_yr){
      input_x_varying = data.frame(age_ns = ns(x = data_df_new$age_c+fu_yr, knots = ns_fit_parameter$age_ns_fit$knots,
                                               Boundary.knots = ns_fit_parameter$age_ns_fit$Boundary.knots, 
                                               intercept = ns_fit_parameter$age_ns_fit$intercept.ind),
                                   yrsincdx_cat_new = relevel(cut(data_df_new$mossince1stdx+12*fu_yr,
                                                                  breaks=c(6,11,35,59,83,120,1200),
                                                                  labels=c("<1","1-2","3-4","5-6","7-9",">10"),
                                                                  include.lowest = TRUE),ref="1-2")
      )
      input_x_fu = model.matrix(formula_glmnet,data=cbind(input_x_constant,input_x_assign,input_x_varying))[,-1]
      risk_pr_1yr_fu = normalizeRisk_siteDxyear(model_object = coef_lasso_primary_deid[[mi_round]],
                                                data = input_x_fu,
                                                weights=w,
                                                designMx_siteDxyear = uniq_comb_glmnet_input_deid)
      return(risk_pr_1yr_fu)
    })
    risk_pr_1yr_compete_4fu_w = sapply(1:4,simplify=FALSE,function(fu_yr){
      input_x_varying = data.frame(age_ns = ns(x = data_df_new$age_c+fu_yr, knots = ns_fit_parameter$age_ns_fit$knots,
                                               Boundary.knots = ns_fit_parameter$age_ns_fit$Boundary.knots,
                                               intercept = ns_fit_parameter$age_ns_fit$intercept.ind),
                                   yrsincdx_cat_new = relevel(cut(data_df_new$mossince1stdx+12*fu_yr,
                                                                  breaks=c(6,11,35,59,83,120,1200),
                                                                  labels=c("<1","1-2","3-4","5-6","7-9",">10"),
                                                                  include.lowest = TRUE),ref="1-2")
      )
      input_x_fu = model.matrix(formula_glmnet,data=cbind(input_x_constant,input_x_assign,input_x_varying))[,-1]
      risk_pr_1yr_fu = normalizeRisk_siteDxyear(model_object = coef_lasso_compete_deid[[mi_round]],
                                                data = input_x_fu,
                                                weights=w,
                                                designMx_siteDxyear = uniq_comb_glmnet_input_deid)
      return(risk_pr_1yr_fu)
    })
    
    print("risk at 4 follow-up exams done!")

    return(list(primary_w = append(list(risk_pr_1yr_index_w),risk_pr_1yr_4fu_w),
                compete_w = append(list(risk_pr_1yr_compete_index_w),risk_pr_1yr_compete_4fu_w)
                ))
  })

  ### === Calculate the average 1-year predictions from 30 models === ###
  ## 1yr risk of primary event at 5 rounds, race-stratified weighted normalization ==
  risk_pr_mx_lasso_1yr_index_w = do.call(cbind,sapply(risk_pr_30lasso, function(res) res[["primary_w"]][1]))
  risk_pr_mx_lasso_1yr_fu1_w = do.call(cbind,sapply(risk_pr_30lasso, function(res) res[["primary_w"]][2]))
  risk_pr_mx_lasso_1yr_fu2_w = do.call(cbind,sapply(risk_pr_30lasso, function(res) res[["primary_w"]][3]))
  risk_pr_mx_lasso_1yr_fu3_w = do.call(cbind,sapply(risk_pr_30lasso, function(res) res[["primary_w"]][4]))
  risk_pr_mx_lasso_1yr_fu4_w = do.call(cbind,sapply(risk_pr_30lasso, function(res) res[["primary_w"]][5]))
  risk_pr_lasso_1yr_index_w = rowMeans(risk_pr_mx_lasso_1yr_index_w)
  risk_pr_lasso_1yr_fu1_w = rowMeans(risk_pr_mx_lasso_1yr_fu1_w)
  risk_pr_lasso_1yr_fu2_w = rowMeans(risk_pr_mx_lasso_1yr_fu2_w)
  risk_pr_lasso_1yr_fu3_w = rowMeans(risk_pr_mx_lasso_1yr_fu3_w)
  risk_pr_lasso_1yr_fu4_w = rowMeans(risk_pr_mx_lasso_1yr_fu4_w)
  ## 1yr risk of competing event at 5 rounds, race-stratified weighted normalization ==
  risk_pr_mx_lasso_1yr_compete_index_w = do.call(cbind,sapply(risk_pr_30lasso, function(res) res[["compete_w"]][1]))
  risk_pr_mx_lasso_1yr_compete_fu1_w = do.call(cbind,sapply(risk_pr_30lasso, function(res) res[["compete_w"]][2]))
  risk_pr_mx_lasso_1yr_compete_fu2_w = do.call(cbind,sapply(risk_pr_30lasso, function(res) res[["compete_w"]][3]))
  risk_pr_mx_lasso_1yr_compete_fu3_w = do.call(cbind,sapply(risk_pr_30lasso, function(res) res[["compete_w"]][4]))
  risk_pr_mx_lasso_1yr_compete_fu4_w = do.call(cbind,sapply(risk_pr_30lasso, function(res) res[["compete_w"]][5]))
  risk_pr_lasso_1yr_compete_index_w = rowMeans(risk_pr_mx_lasso_1yr_compete_index_w)
  risk_pr_lasso_1yr_compete_fu1_w = rowMeans(risk_pr_mx_lasso_1yr_compete_fu1_w)
  risk_pr_lasso_1yr_compete_fu2_w = rowMeans(risk_pr_mx_lasso_1yr_compete_fu2_w)
  risk_pr_lasso_1yr_compete_fu3_w = rowMeans(risk_pr_mx_lasso_1yr_compete_fu3_w)
  risk_pr_lasso_1yr_compete_fu4_w = rowMeans(risk_pr_mx_lasso_1yr_compete_fu4_w)
  
  ### === Integrate 1-year risk over 5 surveillance rounds === ###
  risk_pr_1yr_no_prv_failur_or_compete_w_mx = cbind(index=1,
                                                    fu1=(1-risk_pr_lasso_1yr_compete_index_w)*(1-risk_pr_lasso_1yr_index_w),
                                                    fu2=(1-risk_pr_lasso_1yr_compete_fu1_w)*(1-risk_pr_lasso_1yr_fu1_w),
                                                    fu3=(1-risk_pr_lasso_1yr_compete_fu2_w)*(1-risk_pr_lasso_1yr_fu2_w),
                                                    fu4=(1-risk_pr_lasso_1yr_compete_fu3_w)*(1-risk_pr_lasso_1yr_fu3_w))
  risk_pr_lasso_5yr_w_mx= cbind(index=risk_pr_lasso_1yr_index_w,
                                fu1=risk_pr_lasso_1yr_fu1_w*apply(risk_pr_1yr_no_prv_failur_or_compete_w_mx[,1:2],MARGIN=1,prod),
                                fu2=risk_pr_lasso_1yr_fu2_w*apply(risk_pr_1yr_no_prv_failur_or_compete_w_mx[,1:3],MARGIN=1,prod),
                                fu3=risk_pr_lasso_1yr_fu3_w*apply(risk_pr_1yr_no_prv_failur_or_compete_w_mx[,1:4],MARGIN=1,prod),
                                fu4=risk_pr_lasso_1yr_fu4_w*apply(risk_pr_1yr_no_prv_failur_or_compete_w_mx[,1:5],MARGIN=1,prod))
  risk_pr_lasso_5yr_w = rowSums(risk_pr_lasso_5yr_w_mx)
  
  ### === Combine the predicted risk with the derived data and output the results === ###
  if(nrow(data_df_new)==length(risk_pr_lasso_5yr_w)){
    data_output = data_df_new %>% mutate(risk_5yr = risk_pr_lasso_5yr_w)
    if(nrow(data_df)==1 & nrow(data_output)==2) data_output = data_output[-2,]
  }
  else stop("Error: the length of predicted risk values is not the same as the number of rows of input data.")

  return(predicted_5yr_risk = data_output)
}
