## ===================================== ##
## This is an example of using the function
## 'funcs_BCSC_survFN_5yr' to generate 
## predicted 5=year risk.
## The input data needs to be either a
## data frame or tibble.
## Users can either use the same predictor
## name as shown in the example, or use 
## their own naming fasion for the predictors,
## but a list for variable-name mapping is
## needed for the program to identify the
## corresponding column for each predictor
## if they use user-defined variable names.

## See below for an example:
# var_name = c("race_ethnicity" = "variable_name_for_race/ethn",
#              "family_history" = "variable_name_for_family_history",
#              "breast_density" = "xxxx",
#              "PBC_stage" = "xxxx",
#              "PBC_grade" = "xxxx",
#              "PBC_ERPR" = "xxxx",
#              "PBC_adj_therapy" = "xxxx",
#              "menopause" = "xxxx",
#              "mammo_type" = "xxxx",
#              "PBC_MOD" = "xxxx",
#              "PBC_histology" = "xxxx",
#              "PBC_surgery_type" = "xxxx",
#              "PBC_radiation" = "xxxx",
#              "PBC_dx_age" = "xxxx",
#              "BMI" = "xxxx",
#              "surveillance_interval" = "xxxx",
#              "time_since_prv_mammo" = "xxxx",
#              "exam_age" = "xxxx",
#              "month_since_PBC_dx" = "xxxx"
# )
## ===================================== ##

## How to use the function:
BCSC_survFN_5y(data,
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
               path = "G:/CTRHS/Scc/Projects/AB468DB/test_materials/")
### data: input predictors as a data frame or tibble.
### var_name: a vector of characters for mapping of variable names. Users need to supply the corresponding variable name on the right hand side of each equal sign for each predictor. If the user follows the same name of predictors as shown here in input data, no input for this parameter is needed.
### path: a character string showing the path to the folder where the needed materials are stored.

### === NOTE for users === 
### In a scenario when the input data only has one observation (i.e., data has 1 row only), the program automatically duplicates the record to make the data have 2 rows for computational reasons. In this case, users will see two rows in the output data frame with identical values.
### === === === === === ===

test_input = data.frame(
  racenci_c = c(1,1,1,1), # 1: White, 2: Black, 3: Asian and Pacific Islander, 4: Others, 5: Hispanic
  famhx_c = c("0","0","0","0"), # 0: No, 1: Yes
  density_c = c("3","3","2","2"), # 1: A, 2: B, 3: C, 4: D
  stage1 = c(3,3,2,2), # 1: DCIS, 2: I NOS - IB, 3: II NOS and IIA, 4: IIB+
  grade1_c = c("3", "3","1","1"), # 1: grade 1, 2: grade 2, 3: grade 3
  erplus1 = c("1","1","4","4"), # 1: ER negative and PR negative, 2: ER negative and PR positive, 3: ER positive and PR negative, 4: ER positive and PR positive
  adjthrpy_1 = c("1","1","2","2"), # 0: None, 1: Chemotherapy only, 2: Hormonal therapy only, 3: Both Chemo amd hormonal therapy
  menopht_2cat = c("1","1","1","1"), # 0: Pre/peri-menopausal, 1: Postmenopausal
  mammtype = c("1","1","1","2"), # 1: film, 2: DM, 5: DBT
  mod1_c = c("2","2","1","1"), # 1: screening detected, 2: intervally detected, 3: clinically detected
  hist1_c = c("2","2","3","3"), # 0: DCIS, 1: Invasive NOS, 2: Invasive ductal, 3: Invasive lobular, 4: Invasive mixed
  srg1_new = c("1","1","2","2"), # 1: mastectomy, 2: BCS
  rt_new = c("2","2","1","1"), # 1: with RT, 2: without RT
  dxage = c(59,59,74,74), # continuous value between 21 and 103
  bmi_c = c(21.96570, 20.80086, 25.60605, 25.24025), # continuous value between 15.05 and 89.43
  mossinceprvsurvmm_cat = c("2","4","2","1"), # 1: 9-14 months, 2: 1st surveillance, 3: 3-8 months, 4: 15-23 months, 5: 24+ months
  prvmos_c = c(12, 15, 16, 5), # continuous value between 3 and 254. Values outside the window will be forced to either take the min or max of this range.
  age_c = c(60,62,75,76), # continuous value between 23 and 104
  mossince1stdx = c(12,27,17,29) # continous value no smaller than 6
)

library(here)
source(paste0(here(),"/funcs_BCSC_survFN_5y_deid.R"))
time1 = proc.time()
test_result_2 = BCSC_survFN_5y(data = test_input, path = here())
time2 = proc.time()

print(time2 - time1)

test_result_2$risk_5yr

