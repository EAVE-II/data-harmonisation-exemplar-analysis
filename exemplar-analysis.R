library(tidyr)
library(dplyr)
library(ggplot2)
library(lme4)
library(phsstyles)
library(broom.mixed)


get_data <- function(df){
  df %>% group_by(days_since_vaccineG10) %>%
    mutate(days_since_vaccineG10 = as.integer(as.character(days_since_vaccineG10))) %>% 
    summarise(igg=mean(value_as_number),err=sd(value_as_number)/sqrt(n()),n=n()) %>%
    filter(n>1) %>% return
}
perform_analysis <- function(df,nlme=F,startvec=NULL){
  
  data <- get_data(df) 
  
  f <- ~ a + b*days_since_vaccine*exp(-1*lambda*days_since_vaccine)
  func <- deriv(f,
                namevec=c("a","b","lambda"),
                function.arg=c("days_since_vaccine","a","b","lambda"))
  
  if(is.null(startvec)){
    startvec <- c(a=1000,b=600,lambda=0.04)
  }

  
  m0 <- tryCatch({
    nls(value_as_number ~ func(days_since_vaccine,a,b,lambda),
              start=startvec,
              control=nls.control(maxiter = 1000),
              data=df)
  },error=function(e){
      message(e)
      return (NULL)
  })
  
  fixed_estimates <- NULL
  if(!is.null(m0)){
    fixed_estimates <- tryCatch({
      broom.mixed::tidy(m0, conf.int = TRUE)
      },error=function(e){
          message(e)
          return (NULL)
      })
  }
  if(is.null(fixed_estimates)){
    m0 <- NULL
  }
  
  
  if(nlme){
    startvec <- coef(m0)
    startvec <- c(a=1200,b=300,lambda=0.03)
    startvec <- c(a=1000,b=600,lambda=0.04)
    startvec <- coef(m0)
    startvec <- c(a=1200,b=300,lambda=0.03)
    
    m1 <- nlmer(
      #value_as_number ~ func(days_since_vaccine,a,b,lambda) ~ (a|ageG) + (b|ageG),
      value_as_number ~ func(days_since_vaccine,a,b,lambda) ~ 
        (a|ageG) + (b|ageG) ,
      start=startvec,
      verbose=1,
      #control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
      data=df)
    summary(m1)
    fixed_estimates <-  broom.mixed::tidy(m1, effects = "fixed", conf.int = TRUE)
  }
  else{
    print ('no nlme')
  }
  
  
  model <- NULL
  a <- 1
  b <- 1
  lambda <- 1
  
  a2 <- 1
  b2 <- 1
  lambda2 <- 1
  
  a1 <- 1
  b1 <- 1
  lambda1 <- 1
  
  if(!is.null(m0)){
    model <- m0
    
    a <- fixed_estimates$estimate[[1]]
    b <- fixed_estimates$estimate[[2]]
    lambda <- fixed_estimates$estimate[[3]]
    
    a2 <- fixed_estimates$conf.low[[1]]
    b2 <- fixed_estimates$conf.low[[2]]
    lambda2 <- fixed_estimates$conf.low[[3]]
    
    a1 <- fixed_estimates$conf.high[[1]]
    b1 <- fixed_estimates$conf.high[[2]]
    lambda1 <- fixed_estimates$conf.high[[3]]
  }
 
  
  t <- seq(0,250,0.1)
  
  
  prediction <- as.data.frame(x=t) %>%
    mutate(y=func(t,a,b,lambda),
           y100=func(t,a1,b,lambda),
           y010=func(t,a,b1,lambda),
           y001=func(t,a,b,lambda1),
           y110=func(t,a1,b1,lambda),
           y101=func(t,a1,b,lambda1),
           y111=func(t,a1,b1,lambda1),
           y200=func(t,a2,b,lambda),
           y020=func(t,a,b2,lambda),
           y002=func(t,a,b,lambda2),
           y220=func(t,a2,b2,lambda),
           y202=func(t,a2,b,lambda2),
           y222=func(t,a2,b2,lambda2),
           yup=pmax(y100,y010,y001,y110,y101,y111,y200,y020,y002,y220,y202,y222,na.rm=T),
           ydown=pmin(y100,y010,y001,y110,y101,y111,y200,y020,y002,y220,y202,y222,na.rm=T)
    )
  
  return (list(
    model=model,
    prediction=prediction,
    data=data))
}


df_condition_occurrence <- read.csv("data/condition_occurrence.csv") %>% as_tibble
conditions <- df_condition_occurrence %>% group_by(condition_concept_id) %>%
  summarise(n=n()) %>% arrange(desc(n))
ordered_conditions <- conditions$condition_concept_id %>% as.character()


#m_filter <- 'n_risks > 2'
#m_save <- '2RISKS'

#m_filter <- 'n_risks == 1'
#m_save <- '1RISK'

#m_filter <- paste0(ordered_conditions[9],'==1')


#m_filter <- NULL
#m_save <- 'FULL'

#extract arguments passed to the R script
args = commandArgs(trailingOnly=TRUE)
#first argument passed is the name of a QCOVID variable
m_save = as.character(args[1])
m_filter = as.character(args[2])


if(is.na(m_filter)){
  m_save <- 'FULL'
  m_filter <- NULL
}

print (m_save)
print (m_filter)

#m_filter <- 'n_risks>2'
#m_save <- '3plusRISK'


concepts <- c("37003436"="Pf",
              "37003518"="Md",
              "35894915"="Az")

df_person <- read.csv("data/person.csv")

df_person <-  df_person %>% select(person_id,gender_concept_id,birth_datetime) %>% 
              mutate(age = as.integer(as.numeric(as.Date('2023-03-30') - as.Date(birth_datetime),unit='days')/365.25))


df_drug_exposure <- read.csv("data/drug_exposure.csv")

df_drug_exposure <- df_drug_exposure %>% 
                    select(drug_exposure_id,person_id,drug_concept_id,drug_exposure_start_date) %>% 
                    mutate(drug_concept_name = concepts[as.character(drug_concept_id)],
                            drug_exposure_start_date = as.Date(drug_exposure_start_date)) %>% 
                    group_by(person_id) %>% arrange(person_id,drug_exposure_start_date) %>%
                    mutate(dose=row_number()) %>% ungroup


df_measurement <- read.csv("data/measurement.csv")
df_measurement <- df_measurement %>% select(measurement_id,person_id,measurement_date,value_as_number) %>% 
                   mutate_at(c('measurement_date'),as.Date)


df_condition_occurrence <- df_condition_occurrence %>%
                           select(person_id,condition_concept_id) %>% 
                           pivot_wider(names_from=condition_concept_id,values_from=condition_concept_id) %>% 
                           mutate_at(vars(contains('Q_DIAG')),~if_else(is.na(.), 0, 1)) %>% 
                           mutate(n_risks = rowSums(across(contains('Q_DIAG')))) 


df_measurement <- df_measurement %>% select(measurement_id,person_id,measurement_date,value_as_number) %>% 
  mutate_at(c('measurement_date'),as.Date)



df_ana <- df_measurement %>% left_join(df_drug_exposure) %>% left_join(df_person)
df_ana <- df_ana %>% group_by(person_id,measurement_id) %>% 
           filter(measurement_date > drug_exposure_start_date) %>% 
           ungroup %>% group_by(person_id,measurement_id) %>% 
           filter(row_number() == n())



df_ana <- df_ana %>% left_join(df_condition_occurrence) %>% mutate_at(c('n_risks'),~replace_na(., 0)) 
df_ana %>% group_by(n_risks) %>% summarise(n=n())

nrow(df_ana)
m_filter
if (!is.null(m_filter)){
  df_ana <- df_ana %>% filter(eval(rlang::parse_expr(m_filter)))
}

nrow(df_ana)

df_ana <- df_ana %>% ungroup %>% 
           select(person_id,age,gender_concept_id,value_as_number,measurement_date,
                  drug_exposure_start_date,drug_concept_name,dose) %>% 
           mutate(days_since_vaccine = as.numeric(measurement_date - drug_exposure_start_date)) %>%
           select(-measurement_date,-drug_exposure_start_date) %>% 
           mutate(days_since_vaccineG = cut(days_since_vaccine,
                                           breaks=seq(0,250,20),
                                           labels=seq(10,240,20)),
                  days_since_vaccineG10 = cut(days_since_vaccine,
                                            breaks=seq(0,250,5),
                                            labels=seq(2.5,247.5,5))
                 ) %>%
           mutate(ageG = cut(age,breaks=seq(0,101,10)))


results <- list()


df_ana_1_az <- df_ana %>% filter(dose==1 & !is.na(value_as_number) & drug_concept_name=='Az' & 
                                   days_since_vaccine<80) 
nrow(df_ana_1_az)

results_1_az <- perform_analysis(df_ana_1_az,startvec=c(a=0,b=100,lambda=0.5))

prediction <- results_1_az$prediction
data <- results_1_az$data

p <- prediction %>% ggplot(aes(x=t,y=y,ymax=yup,ymin=ydown)) +
  geom_ribbon(fill='purple',alpha=0.2) +# scale_fill_continuous_phs(palette='main') + 
  geom_line(linetype='dashed') + 
  geom_pointrange(aes(x=days_since_vaccineG10,y=igg,ymin=igg-err,ymax=igg+err),data=data) +
  theme_classic(base_size=25) +
  labs(title='Vaccine = 1st dose AstraZeneca',
       x='Days since vaccination',
       y='Mean IgG titre [U/ml]')

results_1_az$p <- p
results_1_az$p
results[['1st Dose AstraZeneca']] <- results_1_az



df_ana_1_pf <- df_ana %>% filter(dose==1 & !is.na(value_as_number) & drug_concept_name=='Pf' & 
                                   days_since_vaccine<50) 
nrow(df_ana_1_pf)

results_1_pf <- perform_analysis(df_ana_1_pf)#,startvec=c(a=-20,b=100,lambda=0.1))

prediction <- results_1_pf$prediction
data <- results_1_pf$data

p <- prediction %>% ggplot(aes(x=t,y=y,ymax=yup,ymin=ydown)) +
  geom_ribbon(fill='purple',alpha=0.2) +# scale_fill_continuous_phs(palette='main') + 
  geom_line(linetype='dashed') + 
  geom_pointrange(aes(x=days_since_vaccineG10,y=igg,ymin=igg-err,ymax=igg+err),data=data) +
  theme_classic(base_size=25) +
  labs(title='Vaccine = 1st dose Pfizer',
       x='Days since vaccination',
       y='Mean IgG titre [U/ml]')

results_1_pf$p <- p
results_1_pf$p
results[['1st Dose Pfizer']] <- results_1_pf


df_ana_1_md <- df_ana %>% filter(dose==1 & !is.na(value_as_number) & drug_concept_name=='Md' & 
                                   days_since_vaccine<70) 
nrow(df_ana_1_md)

results_1_md <- perform_analysis(df_ana_1_md,startvec=c(a=0,b=1000,lambda=0.01))

prediction <- results_1_md$prediction
data <- results_1_md$data

p <- prediction %>% ggplot(aes(x=t,y=y,ymax=yup,ymin=ydown)) +
  geom_ribbon(fill='purple',alpha=0.2) +# scale_fill_continuous_phs(palette='main') + 
  geom_line(linetype='dashed') + 
  geom_pointrange(aes(x=days_since_vaccineG10,y=igg,ymin=igg-err,ymax=igg+err),data=data) +
  theme_classic(base_size=25) +
  labs(title='Vaccine = 1st dose Moderna',
       x='Days since vaccination',
       y='Mean IgG titre [U/ml]') +
  ylim(c(0,2200))

results_1_md$p <- p
results_1_md$p

results[['1st Dose Moderna']] <- results_1_md


df_ana_2_az <- df_ana %>% filter(dose==2 & !is.na(value_as_number) & drug_concept_name=='Az' & 
                                days_since_vaccine<170) 

results_2_az <- perform_analysis(df_ana_2_az,startvec=c(a=80,b=10,lambda=0.03))
prediction <- results_2_az$prediction
data <- results_2_az$data

p <- prediction %>% ggplot(aes(x=t,y=y,ymax=yup,ymin=ydown)) +
    geom_ribbon(fill='purple',alpha=0.2) +# scale_fill_continuous_phs(palette='main') + 
    geom_line(linetype='dashed') + 
    geom_pointrange(aes(x=days_since_vaccineG10,y=igg,ymin=igg-err,ymax=igg+err),data=data) +
    theme_classic(base_size=25) +
    labs(title='Vaccine = 2nd dose AstraZeneca',
         x='Days since vaccination',
         y='Mean IgG titre [U/ml]')
  
results_2_az$p <- p
results_2_az$p

results[['2nd Dose AstraZeneca']] <- results_2_az



df_ana_2_pf <- df_ana %>% filter(dose==2 & !is.na(value_as_number) & drug_concept_name=='Pf') 
nrow(df_ana_2_pf)

results_2_pf <- perform_analysis(df_ana_2_pf)

prediction <- results_2_pf$prediction
data <- results_2_pf$data

p <- prediction %>% ggplot(aes(x=t,y=y,ymax=yup,ymin=ydown)) +
  geom_ribbon(fill='purple',alpha=0.2) +# scale_fill_continuous_phs(palette='main') + 
  geom_line(linetype='dashed') + 
  geom_pointrange(aes(x=days_since_vaccineG10,y=igg,ymin=igg-err,ymax=igg+err),data=data) +
  theme_classic(base_size=25) +
  labs(title='Vaccine = 2nd dose Pfizer',
       x='Days since vaccination',
       y='Mean IgG titre [U/ml]')


results_2_pf$p <- p
results_2_pf$p

results[['2nd Dose Pfizer']] <- results_2_pf


df_ana_2_md <- df_ana %>% filter(dose==2 & !is.na(value_as_number) & drug_concept_name=='Md') 
nrow(df_ana_2_md)
results_2_md <- perform_analysis(df_ana_2_md,startvec=c(a=100,b=2000,lambda=0.5))

prediction <- results_2_md$prediction
data <- results_2_md$data

p <- prediction %>% ggplot(aes(x=t,y=y,ymax=yup,ymin=ydown)) +
  geom_ribbon(fill='purple',alpha=0.2) +# scale_fill_continuous_phs(palette='main') + 
  geom_line(linetype='dashed') + 
  geom_pointrange(aes(x=days_since_vaccineG10,y=igg,ymin=igg-err,ymax=igg+err),data=data) +
  theme_classic(base_size=25) +
  labs(title='Vaccine = 2nd dose Moderna',
       x='Days since vaccination',
       y='Mean IgG titre [U/ml]')

results_2_md$p <- p
results_2_md$p

results[['2nd Dose Moderna']] <- results_2_md


df_ana_3_pf <- df_ana %>% filter(dose>2 & !is.na(value_as_number) & drug_concept_name=='Pf') 
#table(df_ana_3$drug_concept_name)
nrow(df_ana_3_pf)

results_3_pf <- perform_analysis(df_ana_3_pf)
prediction <- results_3_pf$prediction
data <- results_3_pf$data

p <- prediction %>% ggplot(aes(x=t,y=y,ymax=yup,ymin=ydown)) +
  geom_ribbon(fill='purple',alpha=0.2) +# scale_fill_continuous_phs(palette='main') + 
  geom_line(linetype='dashed') + 
  geom_pointrange(aes(x=days_since_vaccineG10,y=igg,ymin=igg-err,ymax=igg+err),data=data) +
  theme_classic(base_size=25) +
  labs(title='Vaccine = 3+ Pfizer',
       x='Days since vaccination',
       y='Mean IgG titre [U/ml]')

results_3_pf$p <- p
results_3_pf$p

results[['3rd Dose Pfizer']] <- results_3_pf


df_ana_3_md <- df_ana %>% filter(dose>2 & !is.na(value_as_number) & drug_concept_name=='Md') 
#table(df_ana_3$drug_concept_name)
nrow(df_ana_3_md)

results_3_md <- perform_analysis(df_ana_3_md)
prediction <- results_3_md$prediction
data <- results_3_md$data

p <- prediction %>% ggplot(aes(x=t,y=y,ymax=yup,ymin=ydown)) +
  geom_ribbon(fill='purple',alpha=0.2) +# scale_fill_continuous_phs(palette='main') + 
  geom_line(linetype='dashed') + 
  geom_pointrange(aes(x=days_since_vaccineG10,y=igg,ymin=igg-err,ymax=igg+err),data=data) +
  theme_classic(base_size=25) +
  labs(title='Vaccine = 3+ Moderna',
       x='Days since vaccination',
       y='Mean IgG titre [U/ml]')

results_3_md$p <- p
results_3_md$p

results[['3rd Dose Moderna']] <- results_3_md



demo <- table(df_ana$ageG,df_ana$gender_concept_id,df_ana$dose,df_ana$drug_concept_name)

results[['meta']] <- list(demo=demo)

results[['meta']]

fname <- paste0(getwd(),"/results/",m_save,".rds")
print (fname)
saveRDS(results, file = fname)


