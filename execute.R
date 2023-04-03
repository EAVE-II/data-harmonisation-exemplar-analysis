library(rstudioapi)
library(dplyr)
library(tidyr)

df_condition_occurrence <- read.csv("data/condition_occurrence.csv") %>% as_tibble
conditions <- df_condition_occurrence %>% group_by(condition_concept_id) %>%
  summarise(n=n()) %>% arrange(desc(n))
ordered_conditions <- conditions$condition_concept_id %>% as.character()
ordered_conditions

jobs <- lapply(ordered_conditions, function(cname) {
  print (cname)
  filter <- paste0(cname,'==1')
  launcherSubmitJob(
    cname,
    cluster = "Kubernetes",
    tags = c('EAVE-II','Serology'),
    environment = c(R_LIBS_USER='/mnt/homes/calumm09/R/x86_64-pc-linux-gnu-library/3.6'),
    command = '/opt/R/3.6.1/lib/R/bin/R',
    args =  c("--slave", "--no-save", "--no-restore", 
              paste("-f", paste0(getwd(),"/exemplar-analysis.R")),
              paste("--args",cname,filter)),
    workingDirectory = getwd(),
  )
})


filter <- 'n_risks==0'
cname <- '0RISK'
launcherSubmitJob(
  cname,
  cluster = "Kubernetes",
  tags = c('EAVE-II','Serology'),
  environment = c(R_LIBS_USER='/mnt/homes/calumm09/R/x86_64-pc-linux-gnu-library/3.6'),
  command = '/opt/R/3.6.1/lib/R/bin/R',
  args =  c("--slave", "--no-save", "--no-restore", 
            paste("-f", paste0(getwd(),"/exemplar-analysis.R")),
            paste("--args",cname,filter)),
  workingDirectory = getwd(),
)

filter <- 'n_risks==1'
cname <- '1RISK'
launcherSubmitJob(
  cname,
  cluster = "Kubernetes",
  tags = c('EAVE-II','Serology'),
  environment = c(R_LIBS_USER='/mnt/homes/calumm09/R/x86_64-pc-linux-gnu-library/3.6'),
  command = '/opt/R/3.6.1/lib/R/bin/R',
  args =  c("--slave", "--no-save", "--no-restore", 
            paste("-f", paste0(getwd(),"/exemplar-analysis.R")),
            paste("--args",cname,filter)),
  workingDirectory = getwd(),
)



filter <- 'n_risks==2'
cname <- '2RISK'
launcherSubmitJob(
  cname,
  cluster = "Kubernetes",
  tags = c('EAVE-II','Serology'),
  environment = c(R_LIBS_USER='/mnt/homes/calumm09/R/x86_64-pc-linux-gnu-library/3.6'),
  command = '/opt/R/3.6.1/lib/R/bin/R',
  args =  c("--slave", "--no-save", "--no-restore", 
            paste("-f", paste0(getwd(),"/exemplar-analysis.R")),
            paste("--args",cname,filter)),
  workingDirectory = getwd(),
)



filter <- 'n_risks>2'
cname <- '3plusRISK'
launcherSubmitJob(
  cname,
  cluster = "Kubernetes",
  tags = c('EAVE-II','Serology'),
  environment = c(R_LIBS_USER='/mnt/homes/calumm09/R/x86_64-pc-linux-gnu-library/3.6'),
  command = '/opt/R/3.6.1/lib/R/bin/R',
  args =  c("--slave", "--no-save", "--no-restore", 
            paste("-f", paste0(getwd(),"/exemplar-analysis.R")),
            paste("--args",cname,filter)),
  workingDirectory = getwd(),
)

