library(dplyr)
library(tibble)
library(ggplot2)


get_max <- function(df){
  max <- df %>% arrange(desc(y)) %>% head(1)

  max.low <- df %>% arrange(desc(ydown)) %>% head(1)
  max.high <- df %>% arrange(desc(yup)) %>% head(1)

  max_y <- max$y[[1]]
  max_t <- max$t[[1]]
  half <- df %>% filter(t>max_t & y < 0.75*max_y) %>% head(1)
  if(nrow(half)<1){
    return (NULL);
  }

  max_y.low <- max.low$y[[1]]
  max_t.low <- max.low$t[[1]]
  half.low <- df %>% filter(t>max_t.low & y < 0.75*max_y.low) %>% head(1)

  if(nrow(half.low)<1){
    return (NULL);
  }

  max_y.high <- max.high$y[[1]]
  max_t.high <- max.high$t[[1]]
  half.high <- df %>% filter(t>max_t.high & y < 0.75*max_y.high) %>% head(1)
  if(nrow(half.high)<1){
    return (NULL);
  }

  half_t <- half$t[[1]] - max_t
  half_t.low <- half.low$t[[1]] - max_t.low
  half_t.high <- half.high$t[[1]] - max_t.high

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

f <- ~ a + b*days_since_vaccine*exp(-1*lambda*days_since_vaccine)
func <- deriv(f,
              namevec=c("a","b","lambda"),
              function.arg=c("days_since_vaccine","a","b","lambda"))

get_sheets <- function(fname) {

  # getting info about all excel sheets
  sheets <- readxl::excel_sheets(fname)
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x))
  data_frame <- lapply(tibble, as.data.frame)

  # assigning names to data frames
  names(data_frame) <- sheets

  # print data frame
  print(data_frame)
}


path <- "results/ACE/ace.xlsx"
ace <- get_sheets(path)

results <- list()
ace_data <- ace[['1 - AZ']] %>% mutate(name='1st Dose AstraZeneca') %>%
            rbind(
              ace[['1 - Pf']] %>% mutate(name='1st Dose Pfizer')
              ) %>%
            rbind(
              ace[['2 - AZ']]  %>% select(days_since_vaccineG10,igg,err) %>% mutate(name='2nd Dose AstraZeneca')
            ) %>%
            rbind(
              ace[['2 - Pf']] %>% select(days_since_vaccineG10,igg,err) %>%  mutate(name='2nd Dose Pfizer')
            ) %>%
            rbind(
              ace[['3']] %>% select(days_since_vaccineG10,igg,err) %>% mutate(name='3 Doses')
            )

results_1_az <- list()
results_1_az$prediction <- NULL
results_1_az$data <- ace[['1 - AZ']]
results[['1st Dose AstraZeneca']] <- results_1_az

results_1_pf <- list()
results_1_pf$prediction <- NULL
results_1_pf$data <- ace[['1 - Pf']]
results[['1st Dose Pfizer']] <- results_1_pf

results_2_az <- list()
results_2_az$prediction <- NULL
results_2_az$data <- ace[['2 - AZ']] %>% select(days_since_vaccineG10,igg,err)
results[['2nd Dose AstraZeneca']] <- results_2_az

results_2_pf <- list()
results_2_pf$data <- ace[['2 - Pf']] %>% select(days_since_vaccineG10,igg,err)


ace_predictions <- ace[['2 - Pf']][2:4,5:8]

a <- ace_predictions[[1,2]]
b <- ace_predictions[[2,2]]
lambda <- ace_predictions[[3,2]]

a1 <- ace_predictions[[1,3]]
b1 <- ace_predictions[[2,3]]
lambda1 <- ace_predictions[[3,3]]

a2 <- ace_predictions[[1,4]]
b2 <- ace_predictions[[2,4]]
lambda2 <- ace_predictions[[3,4]]

t <- seq(0,250,0.1)

ace_prediction <- as.data.frame(x=t) %>%
  mutate(
    name = '2nd Dose Pfizer',
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
  )
results_2_pf$prediction <- ace_prediction

results[['2nd Dose Pfizer']] <- results_2_pf

results_3_pf <- list()
results_3_pf$prediction <- NULL
results_3_pf$data <- ace[['3']] %>% select(days_since_vaccineG10,igg,err)
results[['3rd Dose']] <- results_3_pf


fname <- paste0(getwd(),"/results/ACE/FULL.rds")
saveRDS(results, file = fname)


p <- ace_prediction %>%
  ggplot(aes(x=t,y=y,ymax=yup,ymin=ydown,fill=as.factor(name))) +
  geom_ribbon(alpha=0.2) +# scale_fill_continuous_phs(palette='main') +
  geom_line(linetype='dashed') +
  geom_pointrange(aes(x=days_since_vaccineG10,
                      y=igg,ymin=igg-err,ymax=igg+err,
                      color=as.factor(name)),data=ace_data) +
  labs(title='ACE',
       x='Days since vaccination',
       y='Mean IgG titre [BAU/ml]') +
  facet_wrap(name ~ .,scales='free_y') +
  theme_classic(base_size=10) +
  theme(
    legend.position="none",
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "lightgrey", colour = NULL, size = 1)
  )
p

files <- list.files(path="results", pattern="FULL.rds", full.names=TRUE, recursive=T)
files


#files <- list.files(path="results/pc", pattern="*.rds", full.names=TRUE, recursive=FALSE)
#files
#files <- c('results/pc/1RISK.rds')
#files <- c('results/pc/0RISK.rds',
#           'results/pc/Q_DIAG_PVD.rds',
#           'results/pc/Q_DIAG_DIABETES_2.rds',
#           'results/pc/Q_DIAG_CIRRHOSIS.rds',
#           'results/pc/Q_DIAG_CHD.rds')

all_stats <- NULL

datas <- NULL
predictions <- NULL

#files <- c('results/pc/FULL.rds')
#files <- c('results/ACE/FULL.rds','results/Elderly/FULL.rds')

for(fname in files){#files[5:length(files)]){#files[1:4]){
  name <- sub('\\.rds$', '', basename(fname) )
  dataset <- strsplit(fname,'/')[[1]][2]
  pc <- readRDS(fname)
  labels(pc)
  stats <- NULL

  labs <- labels(pc)[1:length(pc)-1]
  #labs <- c('2nd Dose Pfizer','2nd Dose AstraZeneca',
  #          '3rd Dose Pfizer','3rd Dose Moderna')

  for (dose in labs){
    data <- pc[[dose]]$data %>% as_tibble %>% mutate(dose=dose,dataset=dataset)

    data <- data %>%
            filter(across(any_of("n"), ~.x > 1)) %>%
            select(days_since_vaccineG10,igg,err,dose,dataset)

    datas <- rbind(datas,data)

    prediction <- pc[[dose]]$prediction

    if(is.null(prediction)){
      next
    }

    prediction <- prediction %>% select(t,y,yup,ydown) %>% mutate(dose=dose,dataset=dataset)

    if(prediction$y[[1]] == 1 & prediction$yup[[1]] == 1 & prediction$ydown[[1]] == 1){
      prediction$y <- NA
      prediction$yup <- NA
      prediction$ydown <- NA
    }

    predictions <- rbind(predictions,prediction)
    stat <- get_max(prediction)
    if(is.null(stat)){
      next
    }
    stat <- stat %>% mutate(name=name)
    stats <- rbind(stats,stat)
  }
  all_stats <- rbind(all_stats,stats)

  #stats %>% ggplot(aes(y=dose,x=estimate,xmin=conf.low,xmax=conf.high)) +
  #          geom_pointrange() +
  #          facet_grid(dose ~ group,scales='free')

}


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


nmap <- c(
  'ace'='ACE',
  'bd'='PHS Blood Donors',
  'pc'='PHS Primary Care',
  'Elderly'='Bham Elderly'
)


ddose <- '1st'
ddose <- '2nd'
#ddose <- '3rd'
predictions %>% filter(grepl(ddose,dose)) %>% filter(ydown>-100) %>%
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
         y='Mean IgG titre [BAU/ml]') +
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

all_stats %>% mutate(estimate=ifelse(estimate<5,NA,estimate)) %>%
          mutate(conf.high = pmax(conf.high,2*estimate-conf.low),
                 conf.low = pmin(conf.low,2*estimate-conf.high)) %>%
          ggplot(aes(y=name,
                     x=estimate,xmin=conf.low,xmax=conf.high,
                     color=group)) +
          geom_pointrange() +
          theme_bw(base_size=10) +
          facet_grid(dose ~ group,scales='free')



