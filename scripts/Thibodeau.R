##date 27th jan 2021

library(readxl)
library(tidyverse)

##read supplementary table 1
sup_table1 <- read_excel("data/supp_mcs.a003681_Supplemental_Table_S1_cancer_types_MUTYH_status.xlsx")
View(sup_table1)
names(sup_table1) <- gsub(" ","_",names(sup_table1))

##count how many samples in total?
sup_table1 %>% summarise(sum(count))
#731

##count number of breast
sup_table1 %>% filter(grepl("Breast", Oncotree_tumour_type)) %>% summarise(sum(count))
#158