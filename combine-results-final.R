library(dplyr)
library(tibble)
library(ggplot2)
library(phsstyles)
library(forcats)
library(ggtext)

f <- ~ a + b*days_since_vaccine*exp(-0.01*lambda*days_since_vaccine)
func <- deriv(f,
              namevec=c("a","b","lambda"),
              function.arg=c("days_since_vaccine","a","b","lambda"))


cname_map <-
  c(
    "FULL"="Full Cohort",
  "Q_DIAG_DIABETES_1"="Diabetes (Type-I)",
"Q_DIAG_DIABETES_2"="Diabetes (Type-II)",
"Q_DIAG_SICKLE_CELL"="Sickle Cell Disease",
"Q_DIAG_RESP_CANCER"="Respiratory Cancer",
"Q_DIAG_ASTHMA"="Asthma",
"Q_DIAG_BLOOD_CANCER"="Haematological Cancer",
"Q_DIAG_CHD"="Coronary Heart Disease ",
"Q_DIAG_COPD"="Chronic Obstructive Pulmonary Disease",
"Q_DIAG_CKD_LEVEL"="Chronic Kidney Disease",
"Q_DIAG_CKD"="Chronic Kidney Disease",
"Q_DIAG_CKD3"="Chronic Kidney Disease (Level 3)",
"Q_DIAG_CKD4"="Chronic Kidney Disease (Level 4)",
"Q_DIAG_CKD5"="Chronic Kidney Disease (Level 5)",
"Q_DIAG_AF"="Atrial Fibrillation",
"Q_DIAG_CCF"="Heart Failure",
"Q_DIAG_EPILEPSY"="Epilepsy",
"Q_DIAG_FRACTURE"="Fracture of hip, spine or humerus",
"Q_DIAG_IMMU"="Immune Deficiency",
"Q_DIAG_NEURO"="Rare Neurone Disease",
"Q_DIAG_PARKINSONS"="Parkinsons",
"Q_DIAG_PULM_RARE"="Cystic Fibrosis or Bronchiectasis",
"Q_DIAG_PVD"="Peripheral Vascular Disease",
"Q_DIAG_RA_SLE"="Rheumatoid Arthritis",
"Q_DIAG_SEV_MENT_ILL"="Severe Mental Health Illness",
"Q_DIAG_STROKE"="Stroke",
"Q_DIAG_SICKLE_CELL"="Sickle Cell Disease",
"Q_DIAG_VTE"="Thrombosis or Pulmonary Embolus",
"Q_DIAG_CEREBRALPALSY"="Cerebral Palsy",
"Q_DIAG_CIRRHOSIS"="Cirrhosis",
"Q_DIAG_CONGEN_HD"="Congenital Heart Disease",
"Q_DIAG_HIV_AIDS"="HIV/Aids",
"Q_DIAG_PULM_HYPER"="Pulmonary Hypertension",
"Q_DIAG_DEMENTIA"="Dementia"
)

get_max <- function(df){

  #df <- predictions %>% filter(dataset=='bd' & dose=='1st Dose Moderna')

  max <- df %>% arrange(desc(y)) %>% head(1)
  max.low <- df %>% arrange(desc(ydown)) %>% head(1)
  max.high <- df %>% arrange(desc(yup)) %>% head(1)
  max_y <- max$y[[1]]
  max_t <- max$t[[1]]
  print (max_y)

  half <- df %>% filter(t>max_t & y < 0.5*max_y) %>% head(1)
  #if(nrow(half)<1){
  #  print ('half bad')
  #  return (NULL);
  #}

  max_y.low <- max.low$y[[1]]
  max_t.low <- max.low$t[[1]]


  half.low <- df %>% filter(t>max_t.low & y < 0.5*max_y.low) %>% head(1)

  #if(nrow(half.low)<1){
  #  print ('half.low bad')
  #  return (NULL);
  #}

  max_y.high <- max.high$y[[1]]
  max_t.high <- max.high$t[[1]]


  half.high <- df %>% filter(t>max_t.high & y < 0.5*max_y.high) %>% head(1)
  #if(nrow(half.high)<1){
  ##  print ('half.high bad')
  #  return (NULL);
  #}

  max_t.low <- min((df %>% filter(yup>max_y))$t)
  max_t.high <- max((df %>% filter(yup>max_y))$t)


  half_t <- NA
  half_t.low <- NA
  half_t.high <- NA

  if(nrow(half)>0){
    half_t <- half$t[[1]] - max_t
  }
  if(nrow(half.low)>0){
    half_t.low <- half.low$t[[1]] - max_t.low
  }
  if(nrow(half.high)>0){
    half_t.high <- half.high$t[[1]] - max_t.high
  }


  max_y <- tibble(group='max_y',dose=max$dose[[1]],
                  estimate=max_y,
                  conf.low=max$ydown[[1]],
                  conf.high=max$yup[[1]])

  max_t <- tibble(group='max_t',dose=max$dose[[1]],
                  estimate=max_t,
                  conf.low=max_t.low,
                  conf.high=max_t.high)

  half_t <- tibble(group='half_t',dose=max$dose[[1]],
                   estimate=half_t,
                   conf.low=half_t.low,
                   conf.high=half_t.high)

  return (max_y %>% rbind(max_t) %>% rbind(half_t));
}

get_prediction <- function(df,modify_lambda=F){
  a <- df$estimate[[1]]
  b <- df$estimate[[2]]
  lambda <- df$estimate[[3]]

  a1 <- df$conf.low[[1]]
  b1 <- df$conf.low[[2]]
  lambda1 <- df$conf.low[[3]]

  a2 <- df$conf.high[[1]]
  b2 <- df$conf.high[[2]]
  lambda2 <- df$conf.high[[3]]

  if(modify_lambda){
    lambda <- lambda*100.
    lambda1 <- lambda1*100.
    lambda2 <- lambda2*100.
  }


  t <- seq(0,250,0.1)

  prediction <- as.data.frame(x=t) %>%
    mutate(
      y=func(t,a,b,lambda),
      y100=func(t,a1,b,lambda),
      y010=func(t,a,b1,lambda),
      y001=func(t,a,b,lambda1),
      y110=func(t,a1,b1,lambda),
      y101=func(t,a1,b,lambda1),
      y011=func(t,a,b1,lambda1),
      y111=func(t,a1,b1,lambda1),
      y200=func(t,a2,b,lambda),
      y020=func(t,a,b2,lambda),
      y002=func(t,a,b,lambda2),
      y220=func(t,a2,b2,lambda),
      y202=func(t,a2,b,lambda2),
      y022=func(t,a,b2,lambda2),
      y222=func(t,a2,b2,lambda2),
      yup=pmax(y100,y010,y001,y110,y101,y011,y111,y200,y020,y002,y220,y202,y022,y222,na.rm=T),
      ydown=pmin(y100,y010,y001,y110,y101,y011,y111,y200,y020,y002,y220,y202,y022,y222,na.rm=T)
    ) %>% return
}



files <- list.files(path="results", pattern="FULL.rds", full.names=TRUE, recursive=T)
files

get_all_params <- function(model,modify_lambda=F){

  #model <- temp$`1st Dose AstraZeneca`

  fixed <- model$fixed_estimates
  if(is.null(fixed)){
    print ('--> fixed_estimates is null')
    return (NULL)
  }
  pred <- get_prediction(fixed,modify_lambda)

  par <- pred %>% select(contains('y')) %>%
      select_if(~ !any(is.na(.))) %>%
      summarise_all(function(x) {which.max(x)})

  max_i <- max(par)
  min_i <- min(par)

  max_t <- pred$t[par$y]
  max_t.high <- pred$t[max_i]
  max_t.low <- pred$t[min_i]

  max_y <- pred$y[par$y]
  max_y.low <- pred$ydown[par$y]
  max_y.high <- pred$yup[par$y]

  y0 <- pred$y[1]
  y0.high <- pred$yup[1]
  y0.low <- pred$ydown[1]

  half_t <- (pred %>% filter(t>max_t & y<(y0+0.5*(max_y-y0))))$t[1]

  half_t.low <- (pred %>% filter(t>max_t & ydown<(y0+0.5*(max_y-y0))))$t[1]


  half_t.high <- (pred %>% filter(t>max_t & yup<(y0+0.5*(max_y-y0))))$t[1]


  if(is.na(half_t.high)){
    half_t.high <- 2*half_t - half_t.low
  }

  if(is.na(half_t.low)){
    half_t.low <- 2*half_t - half_t.high
  }

  #pred %>% ggplot() +
  #  geom_line(aes(x=t,y=y)) +
  #  geom_line(aes(x=t,y=yup),linetype='dashed',color='red') +
  #  geom_line(aes(x=t,y=ydown),linetype='dashed',color='blue') +
  #  geom_vline(xintercept=half_t,linetype='dotted') +
  #  geom_vline(xintercept=half_t.high,linetype='dotted',color='red') +
  #  geom_vline(xintercept=half_t.low,linetype='dotted',color='blue')

  max_y <- tibble(group='max_y',
                  #dose=max$dose[[1]],
                  estimate=max_y,
                  conf.low=max_y.low,
                  conf.high=max_y.high)

  max_t <- tibble(group='max_t',
                  #dose=max$dose[[1]],
                  estimate=max_t,
                  conf.low=max_t.low,
                  conf.high=max_t.high)

  half_t <- tibble(group='half_t',
                   #dose=max$dose[[1]],
                   estimate=half_t,
                   conf.low=half_t.low,
                   conf.high=half_t.high)

  return (max_y %>% rbind(max_t) %>% rbind(half_t));

}


#files <- list.files(path="results/pc", pattern="*.rds", full.names=TRUE, recursive=FALSE)
#files
#files <- c('results/pc/1RISK.rds')
#files <- c('results/pc/0RISK.rds',
#           'results/pc/Q_DIAG_PVD.rds',
#           'results/pc/Q_DIAG_DIABETES_2.rds',
#           'results/pc/Q_DIAG_CIRRHOSIS.rds',
#           'results/pc/Q_DIAG_CHD.rds')

files <- list.files(path="results", pattern="^Q.*\\.rds$", full.names=TRUE, recursive=T)
files

files <- list.files(path="results", pattern="^.*RISK.*\\.rds$", full.names=TRUE, recursive=T)
files

files <- list.files(path="results", pattern="^.*AGE.*\\.rds$", full.names=TRUE, recursive=T)
files

all_stats <- NULL
datas <- NULL
predictions <- NULL

#files <- c('results/pc/FULL.rds')
#files <- c('results/ACE/FULL.rds','results/Elderly/FULL.rds')

for(fname in files){#files[5:length(files)]){#files[1:4]){
  name <- sub('\\.rds$', '', basename(fname) )
  dataset <- strsplit(fname,'/')[[1]][2]
  pc <- readRDS(fname)
  stats <- NULL

  labs <- labels(pc)[1:length(pc)-1]
  #labs <- c('2nd Dose Pfizer','2nd Dose AstraZeneca',
  #          '3rd Dose Pfizer','3rd Dose Moderna')

  #labs <- c('2nd Dose Pfizer','2nd Dose AstraZeneca','2nd Dose Moderna')
  for (dose in labs){

    data <- pc[[dose]]$data %>% as_tibble %>% mutate(dose=dose,
                                                     dataset=dataset,
                                                     name=name)

    data <- data %>%
            filter(across(any_of("n"), ~.x > 1)) %>%
            select(days_since_vaccineG10,igg,err,dose,dataset,name)

    datas <- rbind(datas,data)

    prediction <- pc[[dose]]$prediction

    if(is.null(prediction)){
      next
    }

    prediction <- prediction %>% select(t,y,yup,ydown) %>% mutate(dose=dose,
                                                                  dataset=dataset,
                                                                  name=name)

    if(prediction$y[[1]] == 1 & prediction$yup[[1]] == 1 & prediction$ydown[[1]] == 1){
      prediction$y <- NA
      prediction$yup <- NA
      prediction$ydown <- NA
    }

    predictions <- rbind(predictions,prediction)
    stat <- get_all_params(pc[[dose]])

    #stat <- get_max(prediction)
    if(is.null(stat)){
      next
    }
    stat <- stat %>% mutate(dose=dose,dataset=dataset,name=name)
    stats <- rbind(stats,stat)
  }
  all_stats <- rbind(all_stats,stats)

  #stats %>% ggplot(aes(y=dose,x=estimate,xmin=conf.low,xmax=conf.high)) +
  #          geom_pointrange() +
  #          facet_grid(dose ~ group,scales='free')

}

table(all_stats$dataset,all_stats$dose)



labs <- c(
  "1st Dose AstraZeneca","1st Dose Pfizer","1st Dose Moderna","2nd Dose AstraZeneca",
  "2nd Dose Pfizer","2nd Dose Moderna","3rd Dose ","4th Dose ","5th Dose "
)
lab <- "2nd Dose Pfizer"
for(lab in labs){
  f <- paste0('results/ace/data_ACE_',lab,'.csv')
  if(!file.exists(f)){
    next
  }
  data <- read.csv(f)
  f <- paste0('results/ace/tab_ACE_',lab,'.csv')
  prediction <- NULL
  if(file.exists(f)){
    fit <- read.csv(f)
    prediction <- get_prediction(fit,modify_lambda = T) %>%
      mutate(yup = y + 0.6*(yup-y),ydown = y + 0.6*(ydown-y))
  }

  data <- data %>% mutate(
    dataset='ace',
    dose=lab,
    name='FULL'
  ) %>% select(-X)

  datas <- rbind(datas,data)

  if(is.null(prediction)){
    next
  }

  prediction <- prediction %>% mutate(
    dataset='ace',
    dose=lab,
    name='FULL'
  ) %>% select(predictions%>%colnames)
  predictions <- rbind(predictions,prediction)

  #stat <- get_max(prediction)
  stat <- get_all_params(list('fixed_estimates'=fit),modify_lambda=T)

  if(is.null(stat)){
    next
  }
  stat <- stat %>% mutate(dose=lab,dataset='ace',name='FULL')
  all_stats <- rbind(all_stats,stat)

}


labs <- c(
  "1st Dose AstraZeneca","1st Dose Pfizer","1st Dose Moderna","2nd Dose AstraZeneca",
  "2nd Dose Pfizer","2nd Dose Moderna","3rd Dose ","4th Dose ","5th Dose "
  )


for(lab in labs){
  f <- paste0('results/elderly/data_BHAM_',lab,'.csv')
  if(!file.exists(f)){
    print (f)
    next
  }
  data <- read.csv(f)
  f <- paste0('results/elderly/tab_BHAM_',lab,'.csv')
  prediction <- NULL
  if(file.exists(f)){
    fit <- read.csv(f)
    prediction <- get_prediction(fit,modify_lambda = T) %>%
      mutate(yup = y + 0.6*(yup-y),ydown = y + 0.6*(ydown-y))
  }

  data <- data %>% mutate(
    dataset='elderly',
    dose=lab,
    name='FULL'
  ) %>% select(-X)

  datas <- rbind(datas,data)

  if(is.null(prediction)){
    next
  }

  prediction <- prediction %>% mutate(
    dataset='elderly',
    dose=lab,
    name='FULL'
  ) %>% select(predictions%>%colnames)
  predictions <- rbind(predictions,prediction)

  #stat <- get_max(prediction)
  stat <- get_all_params(list('fixed_estimates'=fit),modify_lambda = T)
  if(is.null(stat)){
    next
  }
  stat <- stat %>% mutate(dose=lab,dataset='elderly',name='FULL')
  all_stats <- rbind(all_stats,stat)

}

nmap <- c(
  'pc'='PHS Primary Care',
  'bd'='PHS Blood Donors',
  'ace'='ACE',
  'elderly'='Bham Elderly'
)


#ddose <- '1st'
#ddose <- '2nd'
#ddose <- '3rd'
#ddose <- '3rd|4th|5th'

for (ddose in c('1st','2nd','3rd|4th|5th')){
  p <- predictions %>% filter(grepl(ddose,dose)) %>% filter(ydown>-100) %>%
    ggplot(aes(x=t,y=y,ymax=yup,ymin=ydown,fill=dose)) +
    geom_ribbon(alpha=0.2) +# scale_fill_continuous_phs(palette='main') +
    geom_line(linetype='dashed') +
    geom_pointrange(aes(x=days_since_vaccineG10,
                        y=igg,ymin=igg-err,ymax=igg+err,
                        #shape=dataset,
                        color=dose),data=datas %>% filter(grepl(ddose,dose))) +
    labs(title='',
         fill='Dose',
         color='Dose',
         x='Days since vaccination',
         y='Mean IgG Titre') +
    scale_color_discrete_phs(palette='all') +
    scale_fill_discrete_phs(palette='all') +
    facet_wrap( ~ dataset,
                labeller=labeller(dataset=nmap),
                scales='free',
              nrow=2) +
    theme_classic(base_size=10) +
    theme(
      #legend.position="none",
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "lightgrey", colour = NULL, size = 1)
    )
  print (p)
  print (ddose)

  pdf(paste0("all_",ddose,".pdf"),
      width = 7.3, height = 4.22)
  print (p)
  dev.off()

}

risk_map <- c(
  '0RISK'='0 Risks',
  '1RISK'='1 Risk',
  '2RISK'='2 Risks',
  '3plusRISK'='3+ Risks'
)
#'1st','2nd','3rd|4th|5th'
#for (ddose in c('1st')){
for (ddose in c('1st','2nd','3rd|4th|5th')){

  data <- datas %>% filter(grepl(ddose,dose))

  p <- predictions %>% filter(grepl(ddose,dose)) %>% filter(ydown>-100) %>%
    mutate(yup=ifelse(yup>2200,2200,yup)) %>%
    mutate(ydown=ifelse(ydown>2200,2200,ydown)) %>%
    mutate(ydown=ifelse(dataset=='bd' & ydown<0,0,ydown)) %>%
    mutate(yup=ifelse(dataset=='bd' & yup>10,10,yup)) %>%
    ggplot(aes(x=t,y=y,ymax=yup,ymin=ydown,fill=dose)) +
    geom_ribbon(alpha=0.2) +# scale_fill_continuous_phs(palette='main') +
    geom_line(linetype='dashed') +
    geom_pointrange(aes(x=days_since_vaccineG10,
                        y=igg,ymin=igg-err,ymax=igg+err,
                        #shape=dataset,
                        color=dose),data=data) +
    labs(title='',
         fill='Dose',
         color='Dose',
         x='Days since vaccination',
         y='Mean IgG Titre') +
    scale_color_discrete_phs(palette='main') +
    scale_fill_discrete_phs(palette='main') +
    facet_grid( dataset ~ name,
                labeller=labeller(dataset=nmap,
                                  name=risk_map),
                scales='free') +
    theme_classic(base_size=10) +
    theme(
      #legend.position="none",
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "lightgrey", colour = NULL, size = 1)
    )
  print (p)
  print (ddose)

  pdf(paste0("risk_",ddose,".pdf"),
      width = 11.2, height = 4.22)
  print (p)
  dev.off()


}


age_map <- c(
  #'AGE0_10'='0-10 year olds',
  'AGE10_20'='10-20 year olds',
  'AGE20_30'='20-30 year olds',
  'AGE30_40'='30-40 year olds',
  'AGE40_50'='40-50 year olds',
  'AGE50_60'='50-60 year olds',
  'AGE60_70'='60-70 year olds',
  'AGE70_80'='70-80 year olds',
  'AGE80_90'='80-90 year olds'#,
  #'AGE90_100'='90-100 year olds'
)
#'1st','2nd','3rd|4th|5th'
doses <- unique(datas$dose)
doses
for (ddose in doses){

  data <- datas %>% filter(grepl(ddose,dose))

  p <- predictions %>% filter(grepl(ddose,dose)) %>% filter(ydown>-100) %>%
    mutate(yup=ifelse(yup>2200,2200,yup)) %>%
    mutate(ydown=ifelse(ydown>2200,2200,ydown)) %>%
    mutate(ydown=ifelse(dataset=='bd' & ydown<0,0,ydown)) %>%
    mutate(yup=ifelse(dataset=='bd' & yup>10,10,yup)) %>%
    mutate(name=age_map[name]) %>% filter(!is.na(name)) %>%
    ggplot(aes(x=t,
               y=y,
               ymax=yup,
               ymin=ydown,
               color=name,
               fill=name)) +
    geom_ribbon(alpha=0.2) +# scale_fill_continuous_phs(palette='main') +
    geom_line(linetype='dashed') +
    #geom_pointrange(aes(x=days_since_vaccineG10,
    #                    y=igg,ymin=igg-err,ymax=igg+err,
    #                    #shape=dataset,
    #                    color=dose),data=data) +
    labs(title='',
         fill='',
         color='',
         x='Days since vaccination',
         y='Mean IgG Titre') +
    scale_color_brewer(palette='Spectral') +
    scale_fill_brewer(palette='Spectral') +
    facet_grid( dataset ~ .,
                labeller=labeller(dataset=nmap,
                                  name=age_map),
                scales='free') +
    theme_classic(base_size=10) +
    theme(
      #legend.position="none",
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "lightgrey", colour = NULL, size = 1)
    )
  p
  print (p)
  print (ddose)

  pdf(paste0("ages ",ddose,".pdf"),
      width = 11.2, height = 4.22)
  print (p)
  dev.off()


}




cnames <- unique(datas$name)

cnames <- c("Q_DIAG_ASTHMA","Q_DIAG_COPD","Q_DIAG_AF","Q_DIAG_CCF",
           "Q_DIAG_PULM_RARE","Q_DIAG_PVD","Q_DIAG_PULM_HYPER","Q_DIAG_VTE",
           "Q_DIAG_FRACTURE","Q_DIAG_RA_SLE","Q_DIAG_STROKE","Q_DIAG_CIRRHOSIS",
           "Q_DIAG_NEURO","Q_DIAG_DEMENTIA","Q_DIAG_PARKINSONS","Q_DIAG_SEV_MENT_ILL",
           "Q_DIAG_DIABETES_1","Q_DIAG_DIABETES_2","Q_DIAG_EPILEPSY","Q_DIAG_CEREBRALPALSY",
           "Q_DIAG_BLOOD_CANCER","Q_DIAG_CKD3","Q_DIAG_CKD4","Q_DIAG_CKD5",
           "Q_DIAG_RESP_CANCER","Q_DIAG_SICKLE_CELL","Q_DIAG_CHD","Q_DIAG_CONGEN_HD")


for( i in seq(1,length(cnames),4)){

  j <- i + 3
  conds <- cnames[i:j]

  ddose <- '2nd'
  temp <- predictions %>% #filter(grepl(ddose,dose)) %>%
    filter(name %in% conds) %>%
    filter(ydown>-100)

  temp <- temp %>% group_by(dataset,name,dose) %>%
          mutate(max_diff = max(yup-y), mean_y=mean(y)) %>%
          filter(max_diff < 2*mean_y) %>%
          ungroup


  dtemp <- datas %>% #filter(grepl(ddose,dose)) %>%
    filter(name %in% conds)
  dtemp %>% group_by(name,dataset) %>% summarise(n=n())


  p <- temp %>%
    ggplot(aes(x=t,y=y,ymax=yup,ymin=ydown,fill=dose,color=dose)) +
    geom_ribbon(alpha=0.4) +# scale_fill_continuous_phs(palette='main') +
    geom_line(linetype='dashed') +
    #geom_pointrange(aes(x=days_since_vaccineG10,
    #                    y=igg,ymin=igg-err,ymax=igg+err,
    #                    fill=dose,
    #                    color=dose),data=dtemp) +
    labs(title='',
         fill='',
         color='',
         x='Days since vaccination',
         y='Mean IgG Titre') +
    scale_color_discrete_phs(palette='all') +
    scale_fill_discrete_phs(palette='all') +
    facet_grid( dataset ~ name,
                labeller=labeller(dataset=nmap,
                                  name=cname_map),
                scales='free') +
    theme_classic(base_size=10) +
    theme(
      #legend.position="none",
      #panel.spacing = unit(0, "lines"),
      panel.border = element_rect(fill = NA, color = "black", linetype = "dashed"),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "lightgrey", colour = NULL, size = 1)
    )
  print (p)
  pdf(paste0("qcovid_",ddose,"_",i,".pdf"),
      width = 12.3, height = 4.22)
  print (p)
  dev.off()
}




all_stats %>% colnames

temp <- c('FULL',cnames)
temp

temp <- lapply(temp,function(x) cname_map[x])

p <- all_stats %>% #conf.high>100000 |
  mutate(estimate=ifelse(estimate<2 |conf.high>50000,NA,estimate)) %>%
  mutate(conf.high=ifelse(conf.high>20000,20000,conf.high)) %>%
  mutate(group = case_when(
    group=='max_y' & dataset=='bd' ~ 'max_y2',
    group=='max_y' & dataset=='ace' ~ 'max_y3',
    group=='max_y' & dataset=='elderly' ~ 'max_y4',
    TRUE ~ group)) %>%
  mutate(group = factor(group,level=c('max_y','max_y2',
                                      'max_y3','max_y4',
                                      'max_t','half_t'))) %>%
  #filter(name!='FULL') %>%
  filter(name=='FULL') %>%
  #mutate(name = fct_rev(factor(cname_map[name],levels=temp))) %>%
  mutate(ddose = case_when(
    grepl('1st',dose) ~ '1st',
    grepl('2nd',dose) ~ '2nd',
    grepl('3rd',dose) ~ '3rd',
    grepl('4th',dose) ~ '4th',
    TRUE ~ '5th'
    )) %>%
  mutate(dddose = case_when(
    grepl('Astra',dose) ~ 'AstraZeneca',
    grepl('Pfizer',dose) ~ 'Pfizer',
    TRUE ~ 'Moderna'
  )) %>%
  mutate(dataset=factor(nmap[dataset],levels=unique(nmap))) %>%
  #filter(group=='max_t') %>%
    #!(group=='max_y' & conf.low<0)) %>%
  mutate(conf.high = pmax(conf.high,2*estimate-conf.low),
         conf.low = pmin(conf.low,2*estimate-conf.high)) %>%
  mutate(conf.low=ifelse(conf.low<0,0,conf.low)) %>%
  ggplot(aes(y= dose,#name,#cname_map[name],
             x=estimate,xmin=conf.low,xmax=conf.high,
             #shape=dataset,
             color=dataset)) +
  geom_pointrange(position=position_dodge(width=0.5)) +
  #xlim(0,100) +
  scale_color_discrete_phs(palette = 'all') +
  theme_bw(base_size=10) +
  facet_grid(ddose ~ group,
              #ddose #fct_rev(as.factor(name))
             #~ factor(dataset,levels=c('pc','bd','ace','elderly'),
            #                                   labels=c('PHS Primary Care','PHS Blood Donors',
             #                                           'ACE','Bham Elderly')),
             scales='free',#fct_rev(as.factor(name)),
             switch='x',
             space='free_y',
             labeller=labeller(dataset=nmap,
                               group=c(
                                 'half_t'='Days From Peak to Half Peak',
                                 'max_t'='Days to Peak',
                                 'max_y'='Maximum IgG (Diasorin)',
                                 'max_y2'='Maximum IgG (ELISA)',
                                 'max_y3'='Maximum IgG (In-House)',
                                 'max_y4'='Maximum IgG (Roche)'
                               ),
                               name=cname_map)) +
  labs(x='',
       y='',
       shape='',
       color='') +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = NA, color = "black"),
    strip.background = element_blank(),
    #strip.placement = 'outside',
    axis.text.y = element_markdown(angle=0),
    strip.text.y.right = element_blank(),
    strip.placement = "outside",
    legend.position="top"
  )
p

pdf('ymax_full.pdf',
    width = 12, height = 5)
print (p)
dev.off()




p <- all_stats %>% #conf.high>100000 |
  mutate(estimate=ifelse(estimate<5 |conf.high>50000,NA,estimate)) %>%
  mutate(conf.high=ifelse(conf.high>20000,20000,conf.high)) %>%
  filter(group=='half_t') %>%
  mutate(group = case_when(
    group=='max_y' & dataset=='bd' ~ 'max_y2',
    group=='max_y' & dataset=='ace' ~ 'max_y3',
    group=='max_y' & dataset=='elderly' ~ 'max_y4',
    TRUE ~ group)) %>%
  mutate(group = factor(group,level=c('max_y','max_y2',
                                      'max_y3','max_y4',
                                      'max_t','half_t'))) %>%
  #filter(name!='FULL') %>%
  filter(grepl('Q_',name)) %>%
  #mutate(name = fct_rev(factor(cname_map[name],levels=temp))) %>%
  mutate(ddose = case_when(
    grepl('1st',dose) ~ '1st',
    grepl('2nd',dose) ~ '2nd',
    grepl('3rd',dose) ~ '3rd',
    grepl('4th',dose) ~ '4th',
    TRUE ~ '5th'
  )) %>%
  mutate(dddose = case_when(
    grepl('Astra',dose) ~ 'AstraZeneca',
    grepl('Pfizer',dose) ~ 'Pfizer',
    TRUE ~ 'Moderna'
  )) %>%
  mutate(dataset=factor(nmap[dataset],levels=unique(nmap))) %>%
  #filter(group=='max_t') %>%
  #!(group=='max_y' & conf.low<0)) %>%
  mutate(conf.high = pmax(conf.high,2*estimate-conf.low),
         conf.low = pmin(conf.low,2*estimate-conf.high)) %>%
  mutate(conf.low=ifelse(conf.low<0,0,conf.low)) %>%
  mutate(name=cname_map[name]) %>%
  ggplot(aes(y=name,#name,#cname_map[name],
             x=estimate,xmin=conf.low,xmax=conf.high,
             #shape=dataset,
             color=dose)) +
  geom_pointrange(position=position_dodge(width=0.5)) +
  #xlim(0,100) +
  scale_color_discrete_phs(palette = 'all') +
  theme_bw(base_size=10) +
  facet_grid(name ~ dataset,
             #ddose #fct_rev(as.factor(name))
             #~ factor(dataset,levels=c('pc','bd','ace','elderly'),
             #                                   labels=c('PHS Primary Care','PHS Blood Donors',
             #                                           'ACE','Bham Elderly')),
             scales='free',#fct_rev(as.factor(name)),
             #switch='x',
             space='free_y',
             labeller=labeller(#dataset=nmap,
                               group=c(
                                 'half_t'='Days From Peak to Half Peak',
                                 'max_t'='Days to Peak',
                                 'max_y'='Maximum IgG (Diasorin)',
                                 'max_y2'='Maximum IgG (ELISA)',
                                 'max_y3'='Maximum IgG (In-House)',
                                 'max_y4'='Maximum IgG (Roche)'
                               ),
                               name=cname_map)) +
  labs(#x='Maximum IgG',
       #x='Days to Max IgG',
       x='Days From Peak to Half Peak',
       y='',
       shape='',
       color='') +
  theme_classic() +
  xlim(0,200) +
  theme(
    panel.background = element_rect(fill = NA, color = "black"),
    strip.background = element_blank(),
    #strip.text.x = element_blank(),
    #strip.placement = 'outside',
    axis.text.y = element_markdown(angle=0),
    strip.text.y.right = element_blank(),
    strip.placement = "outside",
    legend.position="right"
  )
p

pdf('thalf_all.pdf',
    width = 8, height = 10)
print (p)
dev.off()

print(all_stats %>% filter(group=='max_t'),n=100)




t1_p1 <- all_stats %>% filter(group=='max_t') %>%
  mutate(cell=paste0(sprintf("%.1f",estimate)," (",sprintf("%.1f",conf.high),' - ',
                     sprintf("%.1f",conf.low),")")) %>%
  mutate(name=cname_map[name]) %>%
  select(dataset,dose,name,cell) %>%
  arrange(dataset,name,dose)

t1_p2 <- t1_p1 %>% filter(dataset=='bd')
t1_p3 <- t1_p1 %>% filter(dataset=='ace')
t1_p4 <- t1_p1 %>% filter(dataset=='elderly')
t1_p1 <- t1_p1 %>% filter(dataset=='pc')



t1 <- t1_p1 %>% full_join(t1_p2,by=c('name','dose')) %>%
  full_join(t1_p3,by=c('name','dose')) %>%
  full_join(t1_p4,by=c('name','dose')) %>%
          mutate_all(~ifelse(is.na(.),'-',.))

write.csv(t1,'table2_max_y.csv')





datas_norm <- datas %>% group_by(dataset) %>%
              mutate(igg_orig = igg,
                     err_orig = err,
                     err = err/max(igg),
                     igg = igg/max(igg)) %>%
              ungroup


predictions_norm <- predictions %>% group_by(dataset) %>%
  mutate(ymax = max(y),
         y = y/ymax,
         yup=yup/ymax,
         ydown=ydown/ymax) %>%
  ungroup

names <- unique(predictions$dataset)
predictions$dataset = factor(predictions$dataset, levels = names)
datas$dataset = factor(datas$dataset, levels = names)

doses <- sort(unique(datas$dose))
predictions$dose = factor(predictions$dose, levels = doses)
datas$dose = factor(datas$dose, levels = doses)






