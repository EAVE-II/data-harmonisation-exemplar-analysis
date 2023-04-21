# GAM fits ----

fit_gam <- function(gamFit,var,cats = NULL, value = NULL) {
  term_list <- list()
  for (term in labels(gamFit$terms)){
    new_term <- gamFit[["var.summary"]][[term]][[1]]
    term_list <- append(term_list, list(new_term))
  }
  
  names(term_list) <-  labels(gamFit$terms)
  
  if (var =='age'){
    term_list[['age']] <- seq(0,80,1)
  }
  else if(var =='days_since_first_measurement'){
    term_list[['days_since_first_measurement']] <- seq(0,450,1)
  }
  else if(var =='days_since_last_vac'){
    term_list[['days_since_last_vac']] <- seq(0,300,1)
  }
  else if(var =='days_since_vaccine'){
    term_list[['days_since_vaccine']] <- seq(0,150,0.1)
  }
  
  else{
    term_list[[var]] <- seq(0,300,1)
  }
  
  if (!is.null(cats)){
    for (lab in labels(cats)){
      term_list[[lab]] <- cats[[lab]]
    }
  }
  if(!is.null(value)){
    term_list[[var]] <- value
  }
  
  print (term_list)
  
  new_data <- expand.grid(term_list)
  
  pred <- predict.gam(gamFit,new_data,se.fit = TRUE)
  pred <- cbind(new_data, pred)
  #pred <- pred %>% mutate(cat=eavehelpers::get_label(paste0('cat',cat)))
  
  return (pred);
}

plot_gam <- function(gamFit,var, var_ref=0, cats = NULL) {
  
  
  palette <- 'main'
  if (!is.null(cats)){
    palette <- 'all'
  }
  
  pred <- fit_gam(gamFit,var,cats)
  
  
  ref <- pred %>% filter(!!as.name(var)==var_ref) %>% select(fit,cat) %>% rename(ref_fit = fit)
  
  
  pred <- pred %>% left_join(ref, by='cat')
  
  
  
  func <- function(a,b){
    return (exp(a - b));
  }
  
  
  p <- ggplot(pred, aes(x=!!as.name(var), y=func(fit,ref_fit), color=cat)) +
    geom_line(size = 0.5, linetype='dashed') +
    geom_ribbon(aes(x = !!as.name(var), ymin = func(fit-se.fit,ref_fit), ymax = func(fit+se.fit,ref_fit), color=cat, fill=cat), alpha = 0.3) +
    labs(y='Odds Ratio',x=var,color='',fill='') +
    scale_colour_discrete_phs(palette=palette) +
    geom_vline(xintercept=var_ref,linetype='dotted') +
    geom_hline(yintercept=1,linetype='dotted') +
    scale_fill_discrete_phs(palette=palette) +
    #scale_y_log10() +
    theme_classic()
  
  if(palette=='none'){
    p <- p + theme(legend.position="none")
  }
  return(p)
}

get_pterm_or_from_gam <- function (fit,model) {
  #extract all parametric term labels
  pterms <- labels(fit$pterms)
  or <- NULL
  #loop over each term
  for(term in pterms){
    #extract all possible levels from the varible
    #this works best when all these terms are factors (which is what we do for this analysis)
    levs <- levels(model[[term]])
    if (is.null(levs)){
      levs <- unique(sort(model[[term]]))
    }
    else{
      levs <- levels(droplevels(model[[term]]))
    }
    #get the first level as the reference
    ref <- levs[1]
    #get all other levels to calculate ORs with respect to the reference level
    values <- levs[2:length(levs)]
    #loop over the values
    for (value in values){
      #use the or_gam function to get the ORs (and CIs) of this parametic term
      tor <- or_gam(data = model, model = fit, pred = term, values=c(ref,value))
      #book this variable
      or <- or %>% rbind (tor)
    }
  }
  #do some cleaning up of the mini dataframe used to save the ORs and return
  or <- or %>% unite('names',sep='',c(predictor,value2)) %>% select(-value1)
  names(or)[2] <- 'OR'
  names(or)[3] <- 'LCL'
  names(or)[4] <- 'UCL'
  return (or)
}




df <- df_ana %>% filter(dose>0 & !is.na(value_as_number) & 
                          #drug_concept_name=='Az' & #|drug_concept_name=='Md') & 
                          days_since_vaccine<170) %>% 
  mutate(n_risks = as.factor(ifelse(n_risks>4,'5+',n_risks)),
         drug_concept_nameG = as.factor(drug_concept_name),
         sexG = as.factor(gender_concept_id),
         doseG = as.factor(case_when( dose < 3 ~ as.character(dose),
                                      TRUE ~ '3+')),
         doseProductG = as.factor(
           paste0(dose,' - ',drug_concept_name)
         )
  ) 

gfit <- gam(value_as_number ~ doseG +  ageG + sexG +
              s(days_since_vaccine,by=drug_concept_nameG,k=10),
            #s(days_since_vaccine,by=doseProductG,k=10),
            family=gaussian,
            data=df)

#get_pterm_or_from_gam(gfit,df)


pred <- fit_gam(gfit,'days_since_vaccine',
                cats=list(
                  drug_concept_nameG=levels(df$drug_concept_nameG)
                  #doseProductG=levels(df$doseProductG)
                )
) %>%
  as_tibble

pred %>% ggplot(aes(x=days_since_vaccine, y=fit, 
                    color=as.factor(drug_concept_nameG))) +
  geom_line(size = 0.5, linetype='dashed') +
  geom_ribbon(aes(ymin = fit-se.fit, ymax=fit+se.fit, color=drug_concept_nameG, 
                  fill=drug_concept_nameG), alpha = 0.3) +
  labs(y='IgG',x='Days Since Vaccination',color='Product',fill='Product') +
  #geom_vline(xintercept=var_ref,linetype='dotted') +
  #geom_hline(yintercept=1,linetype='dotted') +
  #scale_fill_discrete_phs(palette=palette) +
  #scale_y_log10() +
  theme_classic()



nrow(df)

results_2_pf <- perform_analysis(df,nlme=T,startvec=c(a=800,b=30,lambda=4))
prediction <- results_2_pf$prediction
data <- results_2_pf$data
model <- results_2_pf$model
summary(model)

res <- get_nlmer_results(model,modify=T)

plot_nlmer_results(res$results,res$intercepts) 

p <- prediction %>% ggplot(aes(x=t,y=y,ymax=yup,ymin=ydown)) +
  geom_ribbon(fill='purple',alpha=0.2) +# scale_fill_continuous_phs(palette='main') + 
  geom_line(linetype='dashed') + 
  geom_pointrange(aes(x=days_since_vaccineG10,y=igg,ymin=igg-err,ymax=igg+err),data=data) +
  theme_classic(base_size=25) +
  labs(title='Vaccine = 2nd dose Pfizer',
       x='Days since vaccination',
       y='Mean IgG titre [U/ml]')
p

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

