#The code describes how to implement ICE multimorbidity analysis
# install.packages("asaur")
# install.packages("cmprsk")
library(tableone)
library(dplyr)
library(tidyverse)
library(asaur)
library(survival)
library(cmprsk)
library(mstate)
options(rstudio.help.showDataPreview = FALSE)
result_1 <-read.table('/share/home/end1/r_script_data/multi_status_addcovariate_addbase_diseases_20230718.txt',stringsAsFactors=FALSE,header=TRUE,sep = "\t")
dim(result_1)
# 499793     71
# delete the participants who were diagnosed with COPD/Emph or IHD without follow-up
test<-subset(result_1,xxg_days<0|mzf_days<0|is_common_days<0)
dim(test)
head(test)
# 4644   71
result_1<-subset(result_1,xxg_days>=0 & mzf_days>=0 & is_common_days>=0)
dim(result_1) #495149     71
#change the days to year
result_1$xxg_days<-result_1$xxg_days/365.25
result_1$mzf_days<-result_1$mzf_days/365.25
result_1$is_common_days <- result_1$is_common_days/365.25
result_1$death_days <- result_1$death_days/365.25
result_1$fev1_fvc <- result_1$FEV1_1/result_1$FVC_1
table(result_1$is_common)
#  0      1 
# 486940   8209 

#Group baseline variables
#somoking status
# -3	Prefer not to answer
# 0	Never
# 1	Previous
# 2	Current
result_1<- result_1 %>%mutate(smoking_category=case_when(smoking_status==0|smoking_status==1 ~1,
                                                         smoking_status==2 ~2))
table(result_1$smoking_category)
# 1      2 
# 440484  51770 

#alcohol consumption
# 1	Daily or almost daily
# 2	Three or four times a week
# 3	Once or twice a week
# 4	One to three times a month
# 5	Special occasions only
# 6	Never
# -3	Prefer not to answer
result_1<- result_1%>%mutate(alcohol_category_new = case_when(alcohol_category==1 ~1,
                                                              alcohol_category==2 ~2,
                                                              alcohol_category==3 ~3,
                                                              alcohol_category==4|alcohol_category==5 ~4,
                                                              alcohol_category==6 ~5
))

table(result_1$alcohol_category_new)
#       1      2      3      4      5 
# 100237 113904 127539 112161  39834  
#ethnicity category

result_1 <- result_1%>%mutate(ethnity_category=case_when(ethnity=="1"|ethnity=="1001" | ethnity=="1002" | ethnity=="1003" |ethnity=="2001" | ethnity=="2002" | ethnity=="2003"| ethnity=="2004" ~1,
                                                         ethnity=="3"|ethnity=="3001" | ethnity=="3002" | ethnity=="3003"| ethnity=="3004" ~2,
                                                         ethnity=="4"|ethnity=="4001" | ethnity=="4002" | ethnity=="4003" ~3,
                                                         ethnity=="5"|ethnity=="6" ~4
))
table(result_1$ethnity_category)
# 1      2      3      4 
# 468718   9726   7912   6021 
#physical activity status
table(result_1$physical_activity_status)
#     0      1      2 
# 74860 161726 159946 
#education level
# 1 High: college or university degree.
# 2 Intermediate: A/AS levels or equivalent, 0 levels/GCSEs or equivalent.
# 3 Low: none of the aforementioned
result_1 <- result_1%>%mutate(education_category=case_when(education=="1"~1,
                                                           education=="2"|education=="3" ~2,
                                                           education=="4"| education=="5"|education=="6"|education=="-7"~3
))
table(result_1$education_category)

#    1      2      3 
# 159044 158510 167621 
#multiple deprivation index
result_1 <-  result_1%>%mutate(multiple_deprivation = case_when(multiple_deprivation_england !="" ~multiple_deprivation_england,
                                                                multiple_deprivation_scotland!=""~multiple_deprivation_scotland,
                                                                multiple_deprivation_wales!="" ~multiple_deprivation_wales
))
result_1 <- result_1[order(result_1$multiple_deprivation),]
m = floor(nrow(result_1)/4)
result_1[1:m,(ncol(result_1)+1)]="Q1"
result_1[(m+1):(2*m),ncol(result_1)] = "Q2"
result_1[(2*m+1):(3*m), ncol(result_1)]="Q3"
result_1[(3*m+1):(nrow(result_1)),ncol(result_1)]="Q4"
colnames(result_1)[(ncol(result_1))]=c("multiple_deprivation_category")
table(result_1$multiple_deprivation_category)
sum(table(result_1$multiple_deprivation_category))

# Q1     Q2     Q3     Q4 
# 123787 123787 123787 123788 

#bmi category
result_1 <-result_1%>%mutate(bmi_category = case_when(bmi<18.5 ~2,
                                                      bmi>=18.5 & bmi<25 ~1,
                                                      bmi>=25 & bmi<30 ~3,
                                                      bmi>=30 ~4))
table(result_1$bmi_category)
# 1      2      3      4 
# 160318   2584 209144 120080 

result_1 <- result_1%>%mutate(category = case_when(is_common==1 ~4,
                                                   is_common==0&is_xxg==1 ~2,
                                                   is_common==0&is_mzf==1 ~3,
                                                   is_xxg==0&is_mzf==0 ~1))
table(result_1$category)
#  1      2      3      4 
# 416663  51707  18570   8209 
#family history of IHD and COPD/Emph
result_1 <- result_1 %>% mutate(family_history_chd_category = case_when(family_history_chd==1~1,
                                                                        TRUE ~0))

table(result_1$family_history_chd_category)
#   0      1 
# 281856 213293 
result_1 <- result_1 %>% mutate(family_history_copd_category = case_when(family_history_copd==1~1,
                                                                         TRUE ~0))
table(result_1$family_history_copd_category)
#  0      1 
# 417937  77212 


result_1$education_category <- as.factor(result_1$education_category)
result_1$smoking_category <- as.factor(result_1$smoking_category)
result_1$family_history_chd_category <- as.factor(result_1$family_history_chd_category)
result_1$family_history_copd_category <- as.factor(result_1$family_history_copd_category)
result_1$bmi_category <- as.factor(result_1$bmi_category)
result_1$sex <- as.factor(result_1$sex)
result_1$alcohol_category_new <-as.factor(result_1$alcohol_category_new)
result_1$physical_activity_status <- as.factor(result_1$physical_activity_status)
result_1$ethnity_category <- as.factor(result_1$ethnity_category)
head(result_1)

table(result_1$is_death)


result_1 <- result_1 %>%mutate(is_common_status = case_when(is_death==1&is_common==1&is_common_days==death_days ~1,
                                                          is_death==1&is_common==0&is_common_days>=death_days ~2,
                                                            TRUE ~is_common))
result_1 <- result_1 %>%mutate(is_common_statusy_days = case_when(is_death==1&is_common==1&is_common_days==death_days ~death_days,
                                                                  is_death==1&is_common==0&is_common_days>=death_days ~death_days,
                                                                  TRUE ~ is_common_days))


table(result_1$is_common_status)
table(result_1$is_common)

dim(result_1)

result_test2 <- select(result_1,c("eid","is_common_days","is_common","is_common_statusy_days","is_common_status","sex" , "age" , "smoking_category","alcohol_category_new","ethnity_category", "education_category"  ,  "multiple_deprivation_category" ,
                                               "bmi_category"      ,   "Apolipoprotein_B"          ,   "Cholesterol"                  , "C_reactive_protein"         ,
                                               "Triglycerides","is_base_xxg","is_base_mzf","is_cerebrovascular_disease","is_hypertension","is_asthma","is_cancer","is_dm","family_history_chd_category","family_history_copd_category"))
dim(result_test2) # 495149     25
missing_rate_2 <- colMeans(is.na(result_test2))
missing_rate_2

result_test2 <- na.omit(result_test2)
dim(result_test2)
# 448442     24
table(result_test2$is_common_status)
#    0      1      2 
# 439064   7167   2211 


result_sub <- merge(result_1, result_test2[,c("eid","sex")],by="eid")
dim(result_sub)
#448442     83

lipid_drug <- read.table('/share/home/end1/r_script_data/baseline_lipid_status.txt',stringsAsFactors=FALSE,header=TRUE,sep = "\t")
head(lipid_drug)
dim(lipid_drug)
result_sub <- merge(result_sub, lipid_drug[,c("eid","lipid_drug_status")],by="eid")
dim(result_sub)
#448442     86
res.cox <- coxph(Surv(is_common_days, is_common) ~ sex+age+smoking_category+alcohol_category_new+ethnity_category+education_category+multiple_deprivation_category+bmi_category+Apolipoprotein_B+C_reactive_protein+Triglycerides+is_base_xxg+is_base_mzf+is_cerebrovascular_disease+is_dm+is_hypertension+is_asthma+is_cancer+family_history_chd_category+family_history_copd_category,data=result_test2)
summary(res.cox)
table(result_test2$alcohol_category_new)


res.cox <- coxph(Surv(is_common_days, is_common) ~ Apolipoprotein_B,data=result_sub)
summary(res.cox)

res.cox <- coxph(Surv(is_common_days, is_common) ~ Cholesterol,data=result_sub)
summary(res.cox)

res.cox <- coxph(Surv(is_common_days, is_common) ~ Apolipoprotein_B+Cholesterol,data=result_sub)
summary(res.cox)


head(result_sub)

no_lipid_data <- subset(result_sub,lipid_drug_status==0)
dim(no_lipid_data)
#372734     86


res.cox <- coxph(Surv(is_common_days, is_common) ~ sex.x+age+smoking_category+alcohol_category_new+ethnity_category+education_category+multiple_deprivation_category+bmi_category+Apolipoprotein_B+C_reactive_protein+Triglycerides+is_base_xxg+is_base_mzf+is_cerebrovascular_disease+is_dm+is_hypertension+is_asthma+is_cancer+family_history_chd_category+family_history_copd_category,data=no_lipid_data)
summary(res.cox)

res.cox <- coxph(Surv(is_common_days, is_common) ~ Apolipoprotein_B,data=no_lipid_data)
summary(res.cox)

res.cox <- coxph(Surv(is_common_days, is_common) ~ Cholesterol,data=no_lipid_data)
summary(res.cox)

res.cox <- coxph(Surv(is_common_days, is_common) ~ Apolipoprotein_B+Cholesterol,data=no_lipid_data)
summary(res.cox)



#generate table one
vars <- c( "age" ,"sex.x" , "ethnity_category" ,  "education_category"  ,  "multiple_deprivation_category" , "bmi_category",  "smoking_category"  , "alcohol_category_new","physical_activity_status" , "family_history_chd_category","family_history_copd_category",
           "systolic_blood_pressure"    ,  "diastolic_blood_pressure"  , "Apolipoprotein_B" ,  "Cholesterol" , "C_reactive_protein"  , "Triglycerides" , "is_base_xxg","is_base_mzf","is_cerebrovascular_disease","is_dm","is_hypertension","is_asthma","is_cancer","is_common_days","lipid_drug_status")
factorvars <- c("sex.x" , "ethnity_category" ,  "education_category"  ,  "multiple_deprivation_category" ,   "smoking_status"  , "physical_activity_status" , "alcohol_category_new","bmi_category", "family_history_chd_category","family_history_copd_category","is_base_xxg","is_base_mzf","is_cerebrovascular_disease","is_dm","is_hypertension","is_asthma","is_cancer","lipid_drug_status")
tableone_groups <- CreateTableOne(vars = vars, 
                                  strata = 'category', 
                                  data = result_sub, 
                                  factorVars = factorvars,
                                  addOverall = TRUE)
a_test <-print(tableone_groups,contDigits = 2,catDigits = 2)
# a_test <-print(tableone_groups)
write.csv(a_test,file='/share/home/end1/r_script_data/tableone_ukb_common_mortality_20230801.csv')

no_take_lipid_drug <-subset(result_sub, lipid_drug_status==0)
table(no_take_lipid_drug$category)
round(table(no_take_lipid_drug$category)/table(result_sub$category)*100,2)
dim(no_take_lipid_drug)
dim(no_take_lipid_drug)/dim(result_sub)



vars <- c( "age" ,"sex.x" , "is_common_days")
factorvars <- c("sex.x")
tableone_groups <- CreateTableOne(vars = vars, 
                                  strata = 'category', 
                                  data = result_sub, 
                                  factorVars = factorvars,
                                  addOverall = TRUE) 
# a_test <-print(tableone_groups,contDigits = 2)
a_test <-print(tableone_groups)


# result_test2 <- select(result_1,c("eid","sex" , "age" , "smoking_category","alcohol_category_new","ethnity_category","physical_activity_status", "education_category"  ,  "multiple_deprivation_category" ,
#                                   "bmi_category"      ,   "Apolipoprotein_B"          ,   "Cholesterol"                  , "C_reactive_protein"         ,
#                                   "Triglycerides","is_base_xxg","is_base_mzf","is_cerebrovascular_disease","is_hypertension","is_asthma","is_cancer","is_dm","family_history_chd_category","family_history_copd_category"))



# data_surv <- Surv(result_test2$is_common_days,result_test2$is_common)
# coxph.subdata <- coxph(data_surv ~ sex+age+smoking_category+alcohol_category_new+ethnity_category+education_category+multiple_deprivation_category+bmi_category+Apolipoprotein_B+Cholesterol+C_reactive_protein+Triglycerides+is_base_xxg+is_base_mzf+is_cerebrovascular_disease+is_dm+is_hypertension+is_asthma+is_cancer+family_history_chd_category+family_history_copd_category,data=result_test2,na.action=na.omit, method="breslow")
# step(coxph.subdata)
#Competing risks
# install.packages("riskRegression")
# install.packages("pec")
library('cmprsk')
# library(riskRegression)
# library(pec)
# csc<-CSC(Hist(is_common_statusy_days,is_common_status)~sex+age+smoking_category+alcohol_category_new,data=result_test2)

# 

# dim(test)
# prostateSurvival.Mat<-model.matrix(object =~sex+age+smoking_category+alcohol_category_new, data = result_test2)
# prostateSurvival.Mat
# prostateSurvival.Mat<-prostateSurvival.Mat[, -1]
# 
# cov <- data.frame(age = result_test2$age,
#                   sex = result_test2$sex,
#                   smoking_category = result_test2$smoking_category,
#                   alcohol_category_new = result_test2$alcohol_category_new)

test<-head(result_test2,400000)
cov <- data.frame(
                  sex = result_test2$sex
                 )

# table(result_test2$is_common_status)
rr.PROSTATA.Mort<-crr(ftime = result_test2$is_common_statusy_days, fstatus = result_test2$is_common_status, cov, failcode = 1, cencode = 0)
summary(rr.PROSTATA.Mort)





res.cox <- coxph(Surv(is_common_days, is_common) ~ sex+age+smoking_category+alcohol_category_new+ethnity_category+education_category+multiple_deprivation_category+bmi_category+Apolipoprotein_B+Cholesterol+C_reactive_protein+Triglycerides+is_base_xxg+is_base_mzf+is_cerebrovascular_disease+is_dm+is_hypertension+is_asthma+is_cancer+family_history_chd_category+family_history_copd_category,data=result_test2)
summary(res.cox)
table(result_test2$alcohol_category_new)


table(result_test2$sex)

dim(result_test2)

head(result_test2)
write.csv(result_test2,file = "/share/home/end1/r_script_data/select_variable_competing_risk.csv",col.names = TRUE,row.names = FALSE)
write.csv(result_test2,file = "/share/home/end1/r_script_data/select_variable_competing_risk.csv",col.names = TRUE,row.names = FALSE)


table(result_test2$sex)


table(result_test2$smoking_category)
table(result_test2$alcohol_category_new)


#Analyze people with comorbidities

multimorbidity <- subset(result_sub, category==4)
dim(multimorbidity) #7167   83


multimorbidity <- multimorbidity%>%mutate(mortality_status=case_when(is_death==1 ~1,
                                                                     TRUE ~0))
table(multimorbidity$mortality_status)
# 0    1 
# 4668 2499 

vars <- c( "age" ,"sex.x" , "ethnity_category" ,  "education_category"  ,  "multiple_deprivation_category" , "bmi_category",  "smoking_category"  , "alcohol_category_new","physical_activity_status" , "family_history_chd_category","family_history_copd_category",
           "systolic_blood_pressure"    ,  "diastolic_blood_pressure"  , "Apolipoprotein_B" ,  "Cholesterol" , "C_reactive_protein"  , "Triglycerides" , "is_base_xxg","is_base_mzf","is_cerebrovascular_disease","is_dm","is_hypertension","is_asthma","is_cancer","death_days","lipid_drug_status")
factorvars <- c("sex.x" , "ethnity_category" ,  "education_category"  ,  "multiple_deprivation_category" ,   "smoking_status"  , "physical_activity_status" , "alcohol_category_new","bmi_category", "family_history_chd_category","family_history_copd_category","is_base_xxg","is_base_mzf","is_cerebrovascular_disease","is_dm","is_hypertension","is_asthma","is_cancer","lipid_drug_status")
tableone_groups <- CreateTableOne(vars = vars, 
                                  strata = 'mortality_status', 
                                  data = multimorbidity, 
                                  factorVars = factorvars,
                                  addOverall = TRUE) 
a_test <-print(tableone_groups,contDigits = 2)


no_take_lipid_drug <-subset(multimorbidity, lipid_drug_status==0)
table(no_take_lipid_drug$mortality_status)
round(table(no_take_lipid_drug$mortality_status)/table(multimorbidity$mortality_status)*100,2)
dim(no_take_lipid_drug)
dim(no_take_lipid_drug)/dim(multimorbidity)



res.cox <- coxph(Surv(death_days,is_death) ~ Apolipoprotein_B,data=multimorbidity_test)
summary(res.cox)

res.cox <- coxph(Surv(death_days,is_death) ~ C_reactive_protein,data=multimorbidity_test)
summary(res.cox)

res.cox <- coxph(Surv(death_days,is_death) ~ Apolipoprotein_B+C_reactive_protein,data=multimorbidity_test)
summary(res.cox)

cor(multimorbidity_test$Apolipoprotein_B,multimorbidity_test$C_reactive_protein)
cor(multimorbidity_test$Apolipoprotein_B,multimorbidity_test$Cholesterol)


res.cox <- coxph(Surv(death_days,is_death) ~ sex.x+age+smoking_category+alcohol_category_new+ethnity_category+multiple_deprivation_category+bmi_category+Apolipoprotein_B+is_base_xxg+is_base_mzf+is_cerebrovascular_disease+is_dm+is_cancer,data=multimorbidity_test)
summary(res.cox)

res.cox <- coxph(Surv(death_days,is_death) ~ sex.x+age+smoking_category+alcohol_category_new+ethnity_category+multiple_deprivation_category+bmi_category+C_reactive_protein+is_base_xxg+is_base_mzf+is_cerebrovascular_disease+is_dm+is_cancer,data=multimorbidity_test)
summary(res.cox)

res.cox <- coxph(Surv(death_days,is_death) ~ sex.x+age+smoking_category+alcohol_category_new+ethnity_category+multiple_deprivation_category+bmi_category+Apolipoprotein_B+C_reactive_protein+is_base_xxg+is_base_mzf+is_cerebrovascular_disease+is_dm+is_cancer,data=multimorbidity_test)
summary(res.cox)


test<-subset(multimorbidity,is_death==1)%>%select(c('eid'))
write.csv(test,file = "/share/home/end1/r_script_data/death_id.csv",row.names = FALSE)


#The following analyzes the population who did not use lipid-lowering drugs at baseline.
non_lipid_data <- subset(multimorbidity, lipid_drug_status==0)
dim(non_lipid_data)

res.cox <- coxph(Surv(death_days,is_death) ~ Apolipoprotein_B,data=non_lipid_data)
summary(res.cox)

res.cox <- coxph(Surv(death_days,is_death) ~ sex.x+age+smoking_category+alcohol_category_new+ethnity_category+multiple_deprivation_category+bmi_category+Apolipoprotein_B+C_reactive_protein+is_base_xxg+is_base_mzf+is_cerebrovascular_disease+is_dm+is_cancer,data=non_lipid_data)
summary(res.cox)

table(non_lipid_data$mortality_status)

vars <- c( "age" ,"sex.x" , "ethnity_category" ,  "education_category"  ,  "multiple_deprivation_category" , "bmi_category",  "smoking_category"  , "alcohol_category_new","physical_activity_status" , "family_history_chd_category","family_history_copd_category",
           "systolic_blood_pressure"    ,  "diastolic_blood_pressure"  , "Apolipoprotein_B" ,  "Cholesterol" , "C_reactive_protein"  , "Triglycerides" , "is_base_xxg","is_base_mzf","is_cerebrovascular_disease","is_dm","is_hypertension","is_asthma","is_cancer","death_days","lipid_drug_status")
factorvars <- c("sex.x" , "ethnity_category" ,  "education_category"  ,  "multiple_deprivation_category" ,   "smoking_status"  , "physical_activity_status" , "alcohol_category_new","bmi_category", "family_history_chd_category","family_history_copd_category","is_base_xxg","is_base_mzf","is_cerebrovascular_disease","is_dm","is_hypertension","is_asthma","is_cancer","lipid_drug_status")
tableone_groups <- CreateTableOne(vars = vars,
                                  strata = 'mortality_status', 
                                  data = non_lipid_data, 
                                  factorVars = factorvars,
                                  addOverall = TRUE) 
a_test <-print(tableone_groups,contDigits = 2)



#On the basis of the original analysis, the death population is analyzed according to whether it is ICE or non_ice and what are the differences between the two.
death_status <-read.table('/share/home/end1/r_script_data/death_id_cause_category.txt',stringsAsFactors=FALSE,header=TRUE,sep = "\t")
dim(multimorbidity)
death_multimoribidity <- subset(multimorbidity,is_death==1)
dim(death_multimoribidity)

combined_death_data <- merge(death_multimoribidity,death_status[,c("eid","death_cause_ice")],by="eid")
dim(combined_death_data)

vars <- c( "age" ,"sex.x" , "ethnity_category" ,  "education_category"  ,  "multiple_deprivation_category" , "bmi_category",  "smoking_category"  , "alcohol_category_new","physical_activity_status" , "family_history_chd_category","family_history_copd_category",
           "systolic_blood_pressure"    ,  "diastolic_blood_pressure"  , "Apolipoprotein_B" ,  "Cholesterol" , "C_reactive_protein"  , "Triglycerides" , "is_base_xxg","is_base_mzf","is_cerebrovascular_disease","is_dm","is_hypertension","is_asthma","is_cancer","death_days","lipid_drug_status")
factorvars <- c("sex.x" , "ethnity_category" ,  "education_category"  ,  "multiple_deprivation_category" ,   "smoking_status"  , "physical_activity_status" , "alcohol_category_new","bmi_category", "family_history_chd_category","family_history_copd_category","is_base_xxg","is_base_mzf","is_cerebrovascular_disease","is_dm","is_hypertension","is_asthma","is_cancer","lipid_drug_status")
tableone_groups <- CreateTableOne(vars = vars, 
                                  strata = 'death_cause_ice', 
                                  data = combined_death_data, 
                                  factorVars = factorvars,
                                  addOverall = TRUE) 
a_test <-print(tableone_groups,contDigits = 2)
write.csv(a_test,file='/share/home/end1/r_script_data/tableone_ukb_mortality_20230731.csv')
