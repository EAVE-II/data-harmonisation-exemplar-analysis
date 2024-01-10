library(forcats)
library(ggforce)
library(ggplot2)

files <- list.files(path="results", pattern=".rds", full.names=TRUE, recursive=T)
files

large_table <- NULL
for(fname in files){
  name <- sub('\\.rds$', '', basename(fname))
  dataset <- strsplit(fname,'/')[[1]][2]
  data <- readRDS(fname)


  table <- data.frame(data$meta$demo) #%>% as_tibble
  if(nrow(table)==0){
    next
  }
  table <- table %>% as_tibble
  names(table) <- c('age','sex','dose','product','freq')
  large_table <- rbind(large_table,table %>% mutate(name=name,dataset=dataset))

  next

  p <- table %>%
    mutate(
      freq=ifelse(freq<5,0,freq),
      product = factor(product,levels=c('Az','Pf','Md')),
      dose = fct_collapse(dose,'Dose 1'='1','Dose 2'='2','Dose 3+'=c('3','4')),
      sex=ifelse(sex==8507,'Male','Female'),
      freq=ifelse(sex=='Male',freq,-freq)) %>%
    ggplot(aes(x=freq,y=age,fill=as.factor(sex))) +
    geom_col() +
    geom_vline(xintercept=0,linetype='dashed') +
    scale_x_continuous(label = ~abs(.))+
    facet_grid(dose ~ product,
               scale='free',
               labeller = labeller(product=c('Az'='AstraZeneca',
                                             'Pf'='Pfizer',
                                             'Md'='Moderna'))) +
    scale_fill_brewer(palette = 'Set2') +
    labs(title=name,x='Number of people',y='Age',fill='Sex') +
    theme_classic()

  print (p)

}

table(large_table$name)

full <- large_table %>% filter(name=='FULL')
full  %>% mutate(
  age=as.numeric( sub("\\D*(\\d+).*", "\\1", age)) + 5
) %>%
  group_by(dataset) %>%
  summarise(mean=sum(age*freq)/sum(freq),
            sem=mean/sqrt(sum(freq)))

full %>%
  group_by(dataset,sex) %>%
  summarise(n=sum(freq)) %>%
  group_by(dataset) %>%
  mutate(p=100*n/sum(n))

full %>%
  group_by(dataset) %>%
  summarise(n=sum(freq)) %>%
  group_by(dataset) %>%
  mutate(p=100*n/sum(n))


risks <- large_table %>% filter(grepl('RISK',name))
risks %>%
  group_by(dataset,name) %>%
  summarise(n=sum(freq)) %>%
  group_by(dataset) %>%
  mutate(p=100*n/sum(n))

irisks <- large_table %>% filter(grepl('Q_',name))
t1 <- irisks %>%
  filter(dataset=='pc') %>%
  group_by(dataset,name) %>%
  summarise(n=sum(freq),p=100*n/25682) %>%
  mutate(x=paste0(n," (",round(p,1),")")) %>%
  mutate(y=cname_map[name]) %>% ungroup %>%
  select(y,x)

t2 <- irisks %>%
  filter(dataset=='bd') %>%
  group_by(dataset,name) %>%
  summarise(n=sum(freq),p=100*n/25682) %>%
  mutate(x2=paste0(n," (",round(p,1),")")) %>%
  mutate(y=cname_map[name]) %>% ungroup %>%
  select(y,x2)

t1 %>% left_join(t2) %>% View


p <- large_table %>%
    filter(dataset=='pc') %>%
    filter(grepl('RISK',name)) %>%
    mutate(
      freq=ifelse(freq<5,0,freq),
      product = factor(product,levels=c('Az','Pf','Md')),
      name=factor(name,levels=c('0RISK','1RISK','2RISK','3plusRISK'),
                  labels=c('0','1','2','3+')),
      dose = fct_collapse(dose,'Dose 1'='1','Dose 2'='2','Dose 3+'=c('3','4')),
      sex=ifelse(sex==8507,'Male','Female'),
      freq=ifelse(sex=='Male',freq,-freq)) %>%
    ggplot(aes(x=freq,y=age,
               fill=as.factor(sex),
               alpha=as.factor(name))) +
    geom_col() +
    geom_vline(xintercept=0,linetype='dashed') +
    scale_alpha_discrete(range = c(1., 0.2)) +
    scale_x_continuous(label = ~abs(.))+
    facet_grid(product ~ dose,
               scale='free',
               labeller = labeller(product=c('Az'='AstraZeneca',
                                             'Pf'='Pfizer',
                                             'Md'='Moderna'))) +
    scale_fill_brewer(palette = 'Set2') +
    labs(#title='PHS Blood Donors',
         title='PHS Primary Care',
         x='Number of people',y='Age',
         fill='Sex',alpha='Number of Risks') +
    theme_classic()
p

pdf(paste0("demo_pc.pdf"),
    width = 9.5, height = 4.22)
print (p)
dev.off()

p <- large_table %>% filter(grepl('DIAG',name)) %>%
      mutate(dose = fct_collapse(dose,'Dose 1'='1','Dose 2'='2','Dose 3+'=c('3','4'))) %>%
      mutate(product = factor(product,levels=c('Az','Pf','Md'),
                              labels=c('AstraZeneca','Pfizer','Moderna'))) %>%
      mutate(age=as.numeric(sub("\\D*(\\d+).*", "\\1",age))+5) %>%
      mutate(name = as.factor(cname_map[name])) %>%
      group_by(dataset,name,dose,product) %>%
      summarise(n=n(),
                mean=sum(freq*age)/sum(freq),
                err=mean/sqrt(n)) %>%
      ggplot(aes(x=mean,xmin=mean+err,
                 xmax=mean-err,
                 y=name,
                 color=as.factor(dose),
                 shape=as.factor(product))) +
      geom_pointrange(position=position_dodge(width=1)) +
      labs(x='Mean Age',y='',color='',shape='') +
      theme_classic() +
      facet_grid_paginate(name ~ dataset,scale='free',space='free',
                          page=1,
                          ncol=2,
                          nrow=8,
                 labeller=labeller(dataset=c('bd'='PHS Blood Donors',
                                             'pc'='PHS Primary Care'))) +
      theme(
        panel.border = element_rect(fill = NA, color = "black", linetype = "dashed"),
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = "lightgrey", colour = NULL, size = 1),
        strip.text.y.right = element_blank()
      )

p
pdf(paste0("demo_qcvoid_1.pdf"),
    width = 6.5, height = 6)
print (p)
dev.off()


sail <- read.csv('sail_demo.csv')
sail

sail %>% mutate(
  age=as.numeric( sub("\\D*(\\d+).*", "\\1", age)) + 5
) %>%
  group_by(dataset) %>%
  summarise(mean=sum(age*freq)/sum(freq),
            sem=mean/sqrt(sum(freq)))

sail %>%
  group_by(dataset,sex) %>%
  summarise(n=sum(freq)) %>%
  group_by(dataset) %>%
  mutate(p=100*n/sum(n))

p <- sail %>%
  mutate(freq=ifelse(sex=='Male',freq,-freq)) %>%
  ggplot(aes(x=freq,y=age,
             fill=as.factor(sex))) +
  geom_col() +
  geom_vline(xintercept=0,linetype='dashed') +
  scale_alpha_discrete(range = c(1., 0.2)) +
  scale_x_continuous(label = ~abs(.))+
  scale_fill_brewer(palette = 'Set2') +
  labs(#title='PHS Blood Donors',
    title='NCSI datasets',
    x='Number of people',y='Age',
    fill='Sex',alpha='Number of Risks') +
  facet_grid(. ~ dataset,
             scale='free',
             labeller = labeller(dataset=c('ace'='ACE',
                                           'bham'='Bham Elderly',
                                           'Md'='Moderna'))) +
  theme_classic() +
  theme(panel.spacing = unit(2, "lines"))
p

pdf(paste0("demo_sail.pdf"),
    width = 6.5, height = 3.22)
print (p)
dev.off()



