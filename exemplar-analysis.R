library(tidyr)
library(dplyr)
library(ggplot2)
library(lme4)
library(phsstyles)
library(broom.mixed)


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


df_ana <- df_measurement %>% left_join(df_drug_exposure) %>% left_join(df_person)
df_ana <- df_ana %>% group_by(person_id,measurement_id) %>% 
           filter(measurement_date > drug_exposure_start_date) %>% 
           ungroup %>% group_by(person_id,measurement_id) %>% 
           filter(row_number() == n())


df_ana <- df_ana %>% ungroup %>% 
           select(person_id,age,gender_concept_id,value_as_number,measurement_date,
                  drug_exposure_start_date,drug_concept_name,dose) %>% 
           mutate(days_since_vaccine = as.numeric(measurement_date - drug_exposure_start_date)) %>%
           select(-measurement_date,-drug_exposure_start_date) %>% 
           mutate(days_since_vaccineG = cut(days_since_vaccine,
                                           breaks=seq(0,250,20),
                                           labels=seq(10,240,20))
                 ) %>%
           mutate(ageG = cut(age,breaks=seq(60,101,7)))


df_ana_2 <- df_ana  %>% na.omit %>% filter(dose==2 & drug_concept_name=='Pf') 


f <- ~ a + b*days_since_vaccine*exp(-1*lambda*days_since_vaccine)
func <- deriv(f,
              namevec=c("a","b","lambda"),
              function.arg=c("days_since_vaccine","a","b","lambda"))

startvec <- c(a=1000,b=600,lambda=0.04)

m0 <- nls(value_as_number ~ func(days_since_vaccine,a,b,lambda),
          start=startvec,
          data=df_ana_2)

summary(m0)

bda=0.04)
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
          data=df_ana_2)
summary(m1)



a <- fixed_estimates$estimate[[1]]
b <- fixed_estimates$estimate[[2]]
lambda <- fixed_estimates$estimate[[3]]

a2 <- fixed_estimates$conf.low[[1]]
b2 <- fixed_estimates$conf.low[[2]]
lambda2 <- fixed_estimates$conf.low[[3]]

a1 <- fixed_estimates$conf.high[[1]]
b1 <- fixed_estimates$conf.high[[2]]
lambda1 <- fixed_estimates$conf.high[[3]]

t <- seq(0,250,0.1)

prediction <- as.data.frame(x=t) %>%
              mutate(y=func(t,a,b,lambda),
                     yup=func(t,a1,b1,lambda1),
                     ydown=func(t,a2,b2,lambda2))

prediction %>% ggplot(aes(x=t,y=y,ymax=yup,ymin=ydown)) +
               geom_ribbon(fill='purple',alpha=0.2) +# scale_fill_continuous_phs(palette='main') + 
               geom_line(linetype='dashed') + 
               geom_pointrange(aes(x=days_since_vaccineG,y=igg,ymin=igg-err,ymax=igg+err),data=data) +
               theme_classic(base_size=25) +
               labs(title='Vaccine = 2nd dose Pfizer',
                    x='Days since vaccination',
                    y='Mean IgG titre (Spike Roche)')