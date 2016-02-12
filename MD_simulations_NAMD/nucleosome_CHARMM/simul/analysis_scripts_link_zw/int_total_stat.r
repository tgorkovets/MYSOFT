#R-script to get general statistics about interaction in nucleosome.
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(ggplot2)
library(reshape2)
library(xtable)
library(plyr)
##############
###############
##Loading data frames
#DNA-protein interactions
dna_prot<-read.csv('../analysis_data/dna_prot_raw_df.csv')
dna_prot_cryst<-read.csv('../analysis_data/dna_prot_raw_df_cryst.csv')

#Protein-protein interactions
# prot_prot<-read.csv('../analysis_data/prot_prot_raw_df.csv')
# prot_prot_cryst<-read.csv('../analysis_data/prot_prot_raw_df_cryst.csv')

# #Nucleosome-ion interactions
# ion_int<-read.csv('../analysis_data/ions_int_raw_df.csv')
# ion_int_cryst<-read.csv('../analysis_data/ions_int_raw_df_cryst.csv')

# #NUcleosome-water interactions

# wat_int<-read.csv('../analysis_data/wat_int_raw_df.csv')
# wat_int_cryst<-read.csv('../analysis_data/wat_int_raw_df_cryst.csv')


# prot_prot_cryst=subset(prot_prot_cryst,PROT1_chain!=PROT2_chain)
# prot_prot=subset(prot_prot,PROT1_chain!=PROT2_chain)


#Here we need to filter out HB and IP from contacts, and recalcualte WM interactions.
############
#data frame SQL-like magic
a=split(dna_prot_cryst,dna_prot_cryst$type)
#get rid of duplicated in contacts via tricky merge and join
m=merge(a$SC,a$SB,by.x=c("DNA_atom","DNA_chain","DNA_resid","PROT_atom","PROT_chain","PROT_resid"),by.y=c("DNA_atom","DNA_chain","DNA_resid","PROT_atom","PROT_chain","PROT_resid"),all.x=TRUE,suffix=c('','y'))
m2=merge(m,a$HB,by.x=c("DNA_atom","DNA_chain","DNA_resid","PROT_atom","PROT_chain","PROT_resid"),by.y=c("DNA_atom","DNA_chain","DNA_resid","PROT_atom","PROT_chain","PROT_resid"),all.x=TRUE,suffix=c('','z'))
a$vdW=m2[is.na(m2$typez) & is.na(m2$typey),-(11:20)]
a$vdW$type<-factor('vdW')

#Let's simlify now the water mediated interactions
a$WM=a$WM[!duplicated(a$WM),]

dna_prot_cryst=do.call(rbind,a)

a=split(dna_prot,dna_prot$type)
#get rid of duplicated in contacts via tricky merge and join
m=merge(a$SC,a$SB,by.x=c("Time","DNA_atom","DNA_chain","DNA_resid","PROT_atom","PROT_chain","PROT_resid"),by.y=c("Time","DNA_atom","DNA_chain","DNA_resid","PROT_atom","PROT_chain","PROT_resid"),all.x=TRUE,suffix=c('','y'))
m2=merge(m,a$HB,by.x=c("Time","DNA_atom","DNA_chain","DNA_resid","PROT_atom","PROT_chain","PROT_resid"),by.y=c("Time","DNA_atom","DNA_chain","DNA_resid","PROT_atom","PROT_chain","PROT_resid"),all.x=TRUE,suffix=c('','z'))
a$vdW=m2[is.na(m2$typez) & is.na(m2$typey),-(12:20)]
a$vdW$type<-factor('vdW')

#Let's simlify now the water mediated interactions
a$WM=a$WM[!duplicated(a$WM),]

dna_prot=do.call(rbind,a)

#-----------

# a=split(prot_prot_cryst,prot_prot_cryst$type)
# #get rid of duplicated in contacts via tricky merge and join
# m=merge(a$SC,a$SB,by.x=c("PROT1_atom","PROT1_chain","PROT1_resid","PROT2_atom","PROT2_chain","PROT2_resid"),by.y=c("PROT1_atom","PROT1_chain","PROT1_resid","PROT2_atom","PROT2_chain","PROT2_resid"),all.x=TRUE,suffix=c('','y'))
# m2=merge(m,a$HB,by.x=c("PROT1_atom","PROT1_chain","PROT1_resid","PROT2_atom","PROT2_chain","PROT2_resid"),by.y=c("PROT1_atom","PROT1_chain","PROT1_resid","PROT2_atom","PROT2_chain","PROT2_resid"),all.x=TRUE,suffix=c('','z'))

# a$vdW=m2[is.na(m2$typez) & is.na(m2$typey),-(11:20)]
# a$vdW$type<-factor('vdW')

# # head(a$vdW)


# #Let's simlify now the water mediated interactions
# a$WM=a$WM[!duplicated(a$WM),]
# # print(a)
# prot_prot_cryst=do.call(rbind,a)
# # head(subset(prot_prot_cryst,type=='SC'))

# a=split(prot_prot,prot_prot$type)
# #get rid of duplicated in contacts via tricky merge and join
# m=merge(a$SC,a$SB,by.x=c("Time","PROT1_atom","PROT1_chain","PROT1_resid","PROT2_atom","PROT2_chain","PROT2_resid"),by.y=c("Time","PROT1_atom","PROT1_chain","PROT1_resid","PROT2_atom","PROT2_chain","PROT2_resid"),all.x=TRUE,suffix=c('','y'))
# m2=merge(m,a$HB,by.x=c("Time","PROT1_atom","PROT1_chain","PROT1_resid","PROT2_atom","PROT2_chain","PROT2_resid"),by.y=c("Time","PROT1_atom","PROT1_chain","PROT1_resid","PROT2_atom","PROT2_chain","PROT2_resid"),all.x=TRUE,suffix=c('','z'))
# a$vdW=m2[is.na(m2$typez) & is.na(m2$typey),-(12:20)]
# a$vdW$type<-factor('vdW')

# #Let's simlify now the water mediated interactions
# a$WM=a$WM[!duplicated(a$WM),]

# prot_prot=do.call(rbind,a)
#----
##################

cat('======Analyzing DNA-protein interactions all:=========\n')
cat('Total interactions in crystal:\n')
dna_prot_cryst_sdf=data.frame(t(summary(dna_prot_cryst$type)),row.names=c('DNA-protein Crystal'))
print(dna_prot_cryst_sdf)
cat('In dynamics averaged over ', length(table(dna_prot$Time)),'frames:\n')
dna_prot_sdf=data.frame(t(round(summary(dna_prot$type)/length(table(dna_prot$Time)),digits=2)),row.names=c('DNA-protein'))
print(dna_prot_sdf)
# xtable(summary(dna_prot['type']))

dna_prot=subset(dna_prot, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid > 43)&(PROT_resid < 132))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid > 23)&(PROT_resid < 99))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid > 15)&(PROT_resid < 118))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid > 29)&(PROT_resid < 124))))
dna_prot_cryst=subset(dna_prot_cryst, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid > 43)&(PROT_resid < 132))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid > 23)&(PROT_resid < 99))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid > 15)&(PROT_resid < 118))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid > 29)&(PROT_resid < 124))))

cat('======Analyzing DNA-protein interactions core:=========\n')
cat('Total interactions in crystal:\n')
dna_prot_cryst_sdf=data.frame(t(summary(dna_prot_cryst$type)),row.names=c('DNA-protein Crystal'))
print(dna_prot_cryst_sdf)
cat('In dynamics averaged over ', length(table(dna_prot$Time)),'frames:\n')
dna_prot_sdf=data.frame(t(round(summary(dna_prot$type)/length(table(dna_prot$Time)),digits=2)),row.names=c('DNA-protein'))
print(dna_prot_sdf)
# xtable(summary(dna_prot['type']))


quit()
cat('======Analyzing protein-protein interactions:========\n')

# prot_prot_cryst=subset(prot_prot_cryst,type!='SC')
# prot_prot=subset(prot_prot,type!='SC')

cat('Total interactions in crystal:\n')
prot_prot_cryst_sdf=data.frame(t(summary(prot_prot_cryst$type)),row.names=c('Protein-protein Crystal'))
print(prot_prot_cryst_sdf)
cat('In dynamics averaged over ', length(table(prot_prot$Time)),'frames:\n')
prot_prot_sdf=data.frame(t(round(summary(prot_prot$type)/length(table(prot_prot$Time)),digits=2)),row.names=c('Protein-protein'))
print(prot_prot_sdf)

quit()


cat('Analyzing nucleosome-ion interactions:\n')
cat('Total interactions in crystal:\n')
ion_int_cryst_sdf=data.frame(t(summary(ion_int_cryst$type)),row.names=c('Ion-nucleosome Crystal'))
print(ion_int_cryst_sdf)
cat('In dynamics averaged over ', length(table(ion_int$Time)),'frames:\n')
ion_int_sdf=data.frame(t(round(summary(ion_int$type)/length(table(ion_int$Time)),digits=2)),row.names=c('Ion-nucleosome'))
print(ion_int_sdf)

quit()

cat('Analyzing nucleosome-water interactions:\n')
cat('Total H-bonds with water in crystal:\n')
cat(sum(wat_int_cryst$num_HB),'\n')
wat_int_cryst_sdf=data.frame(HB=sum(wat_int_cryst$num_HB),row.names=c('Water-nucleosome Crystal'))
print(wat_int_cryst_sdf)
cat('In dynamics average\n')
wat_int_sdf=data.frame(HB=sum(wat_int$num_HB),row.names=c('Water-nucleosome'))
print(wat_int_sdf)

RBIND <- function(datalist) {
  require(plyr)
  temp <- rbind.fill(datalist)
  rownames(temp) <- unlist(lapply(datalist, row.names))
  temp
}

df=RBIND(list(dna_prot_cryst_sdf,dna_prot_sdf,prot_prot_cryst_sdf,prot_prot_sdf,ion_int_cryst_sdf,ion_int_sdf,wat_int_cryst_sdf,wat_int_sdf))
xtable(df)
