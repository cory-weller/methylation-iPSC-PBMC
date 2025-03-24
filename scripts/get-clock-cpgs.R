#!/usr/bin/env Rscript

library(ExperimentHub)
library(data.table)
eh <- ExperimentHub()

pData <- query(eh , "methylclockData")

dt <- mcols(pData)
dt$EHID <- rownames(dt)
dat <- as.data.table(dt)

# Look at data saved within this experiment object
dat[, .SD, .SDcols=c('EHID','title')]

#      EHID                                            title
#     <char>                                           <char>
#  1: EH3913 Datasets to estimate cell counts for EEAA method
#  2: EH6068                                   CpGs BNN clock
#  3: EH6069                      Coefficients Bohlin's clock
#  4: EH6070                      Coefficients Hannum's clock
#  5: EH6071                     Coefficients Hobarth's clock # ACTUALLY HORVATH... mis-named in package
#  6: EH6072                      Coefficients Knight's clock
#  7: EH6073                         Coefficients Lee's clock
#  8: EH6074                        Coefficients Levine clock
#  9: EH6075                       Coefficients Mayne's clock
# 10: EH6076                         Coefficients PedBE clock
# 11: EH6077          Coefficients Horvath’s skin+blood clock
# 12: EH6078             Coefficients Telomere Length’s clock
# 13: EH6079                          Coefficients Wu's clock
# 14: EH6080                         Methylation Data Example
# 15: EH6081                             probe Annotation 21k
# 16: EH6082                                     Test Dataset
# 17: EH6083                                       References
# 18: EH7367                           Coefficients BLUPclock
# 19: EH7368                            Coefficients EN clock
# 20: EH7369                          Coefficients EPIC clock


# Manually extract the relevant objects with CpG information

BN_cpgs <- pData[['EH6068']]
Bohlin_cpgs <- pData[['EH6069']]$CpGmarker[-c(1)]       # excluded first element which is '(Intercept)'
Hannum_cpgs <- pData[['EH6070']]$CpGmarker
Horvath_cpgs <- pData[['EH6071']]$CpGmarker[-c(1)]      # excluded first element which is '(Intercept)'
Knight_cpgs <- pData[['EH6072']]$CpGmarker[-c(1)]      # excluded first element which is '(Intercept)'
Lee_cpgs <- pData[['EH6073']]$CpGmarker[-c(1)]      # excluded first element which is '(Intercept)'
Levine_cpgs <- pData[['EH6074']]$CpGmarker[-c(1)]      # excluded first element which is '(Intercept)'
Mayne_cpgs <- pData[['EH6075']]$CpGmarker[-c(1)]      # excluded first element which is '(Intercept)'
PedBE_cpgs <- pData[['EH6076']]$CpGmarker[-c(1)]      # excluded first element which is '(Intercept)'
Horvath_skin_blood_cpgs <- pData[['EH6077']]$CpGmarker[-c(1)]      # excluded first element which is '(Intercept)'
Wu_cpgs <- pData[['EH6079']]$CpGmarker[-c(1)]      # excluded first element which is '(Intercept)'
BLUP_cpgs <- pData[['EH7367']]$CpGmarker[-c(1)]      # excluded first element which is '(Intercept)'
EN_cpgs <- pData[['EH7368']]$CpGmarker[-c(1)]      # excluded first element which is '(Intercept)'
EPIC_cpgs <- pData[['EH7369']]$CpGmarker[-c(1)]      # excluded first element which is '(Intercept)'

# Remove other objects before saving Rdata object
rm(dt)
rm(eh)
rm(dat)
rm(pData)

# Export Rdata object
save.image(file='CLOCKS/clock-cpgs.Rdata')

# Export text files
writeLines(BN_cpgs, con=paste0('CLOCKS/BN_cpgs.txt'))
writeLines(Bohlin_cpgs, con=paste0('CLOCKS/Bohlin_cpgs.txt'))
writeLines(Hannum_cpgs, con=paste0('CLOCKS/Hannum_cpgs.txt'))
writeLines(Horvath_cpgs, con=paste0('CLOCKS/Horvath_cpgs.txt'))
writeLines(Knight_cpgs, con=paste0('CLOCKS/Knight_cpgs.txt'))
writeLines(Lee_cpgs, con=paste0('CLOCKS/Lee_cpgs.txt'))
writeLines(Levine_cpgs, con=paste0('CLOCKS/Levine_cpgs.txt'))
writeLines(Mayne_cpgs, con=paste0('CLOCKS/Mayne_cpgs.txt'))
writeLines(PedBE_cpgs, con=paste0('CLOCKS/PedBE_cpgs.txt'))
writeLines(Horvath_skin_blood_cpgs, con=paste0('CLOCKS/Horvath_skin_blood_cpgs.txt'))
writeLines(Wu_cpgs, con=paste0('CLOCKS/Wu_cpgs.txt'))
writeLines(BLUP_cpgs, con=paste0('CLOCKS/BLUP_cpgs.txt'))
writeLines(EN_cpgs, con=paste0('CLOCKS/EN_cpgs.txt'))
writeLines(EPIC_cpgs, con=paste0('CLOCKS/EPIC_cpgs.txt'))