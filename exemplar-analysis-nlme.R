

# NLME ----


df <- df_ana %>% filter(dose>0 & !is.na(value_as_number)  
                        & days_since_vaccine < 170 
                        #drug_concept_name=='Az' & #|drug_concept_name=='Md') & 
                        & !(dose==1 & days_since_vaccine>80)
) %>% 
  mutate(
    n_risks = as.factor(ifelse(n_risks>4,'5+',n_risks)),
    #n_risks = as.factor(ifelse(n_risks>1,'2+',n_risks)),
    drug_concept_nameG = as.factor(drug_concept_name),
    sexG = as.factor(gender_concept_id),
    doseG = as.factor(case_when( dose < 3 ~ as.character(dose),
                                 TRUE ~ '3+')),
    doseProductG = as.factor(
      paste0(dose,' - ',drug_concept_name)
    )
  ) 


df %>% colnames


grp_name <- c(
  'ageG'='Age',
  'n_risks'='Number of Risk Groups',
  'dose'='Dose',
  'drug_concept_name'='Vaccine Product',
  'dose_history'='Dose History',
  'dose_history2'='Prior Dose History'
)



fit_nlmer <- function(df,startvec=c(a=1,b=1,lambda=1)){
  nlmer(
    value_as_number ~ func(days_since_vaccine, a, b, lambda) ~ 
      #+ (lambda|sexG)
      #+ (a|dose)
    + (a|dose_history2)
    #+ (b|dose)
    + (a | ageG)
    + (b |ageG) 
      + (a|n_risks)
    + (b|n_risks)
    #+ (lambda|n_risks)
    #+ (a|Q_DIAG_DIABETES_2)
    #+ (a|Q_DIAG_CKD3)
    #+ (a|drug_concept_name)
    + (b|drug_concept_name)
    #+ (lambda|ageG)
    + (lambda|drug_concept_name) 
    #+ (b|drug_concept_name)
    #value_as_number ~ func(days_since_vaccine, a, b, lambda) ~ (lambda|ageG) + (a|ageG) + (b|ageG)
    #+ (a|drug_concept_name) + (b|drug_concept_name) + (lambda|drug_concept_name)
    #+ (lambda|gender_concept_id) + (a|gender_concept_id) 
    #+ (a|days_between_vaccinesG) 
    ,
    start=startvec,
    verbose=1,
    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10000)),
    data=df)
}



results <- perform_analysis(df,nlme=T,startvec=c(a=800,b=50,lambda=3))
#results <- perform_analysis(df,nlme=T,startvec=c(a=4,b=2,lambda=2))
prediction <- results$prediction
data <- results$data
model <- results$model
summary(model)

res <- get_nlmer_results(model,modify=T)

term_name <- c(
  'a'='\u03B1',
  'b'='\u03B2',
  'lambda'='\u03BB'
)

results <- res$results
intercepts <- res$intercepts


results <- results %>% mutate(group = factor(grp_name[group],levels=grp_name),
                              level = case_when(
                                #level=='Pf' ~ 'Pfizer',
                                #level=='Az' ~ 'AstraZeneca',
                                #level=='Md' ~ 'Moderna',
                                TRUE ~ level
                              ),
                              term= factor(term_name[term],levels=term_name))


intercepts <- intercepts %>% mutate(term=factor(term_name[term],levels=term_name))


results <- results %>% filter(!grepl('NA',level))



p <- plot_nlmer_results(results,intercepts) 
p

pdf('nlme_pc_v4.pdf',width=10,height=7)
print (p)
dev.off()






png('nlme_pc.png',width=800,height=500)
print (p)
dev.off()


p <- prediction %>% ggplot(aes(x=t,y=y,ymax=yup,ymin=ydown)) +
  geom_ribbon(fill='purple',alpha=0.2) +# scale_fill_continuous_phs(palette='main') + 
  geom_line(linetype='dashed') + 
  geom_pointrange(aes(x=days_since_vaccineG10,y=igg,ymin=igg-err,ymax=igg+err),data=data) +
  theme_classic(base_size=25) +
  labs(title='Vaccine = 2nd dose Pfizer',
       x='Days since vaccination',
       y='Mean IgG titre [U/ml]')
p




prediction
model
data


results_2_pf <- perform_analysis(df)#,nlme=F,startvec=c(a=800,b=30,lambda=4))
prediction_nls <- results_2_pf$prediction
prediction_nls

prediction_both <- prediction %>% mutate(type='NLME') %>%
  rbind(
    prediction_nls %>% mutate(type='NLS')
  )

p <- prediction_both %>% ggplot(aes(x=t,y=y,ymax=yup,ymin=ydown)) +
  geom_ribbon(aes(fill=as.factor(type)),alpha=0.2) +# scale_fill_continuous_phs(palette='main') + 
  geom_line(aes(fill=as.factor(type)),linetype='dashed') + 
  geom_pointrange(aes(x=days_since_vaccineG10,y=igg,ymin=igg-err,ymax=igg+err),data=data) +
  theme_classic(base_size=25) +
  labs(title='Vaccine = 2nd dose Pfizer',
       x='Days since vaccination',
       y='Mean IgG titre [U/ml]')
p


results_2_pf$p <- p
results_2_pf$p






