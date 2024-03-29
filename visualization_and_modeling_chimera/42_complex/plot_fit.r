#R-script analysis of dna protein interactions
#Additional analysis asked by Zhurkin
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(ggplot2)
library(reshape2)
library(xtable)
library(plyr)
library(gridExtra)
##############
###############
##Loading data frames
#DNA-protein interactions
# dna_prot<-read.csv('../analysis_data/dna_prot_raw_df.csv')
pf<-read.csv('fit.dat',na.strings=c('NA'))
print(pf)
pf$overlap=pf$overlap/1000


pf$overlap=pf$overlap-mean(pf$overlap)
pf$cor=pf$cor-mean(pf$cor)
pf$cor_mean=pf$cor_mean-mean(pf$cor_mean)

pf=melt(pf,id.var=c('str_post'))
p=ggplot(data=pf,aes(x=str_post,y=value, color=variable))+geom_line()+geom_point()+xlab('Last protected position in nucleosome from dyad')+ylab('Relative value')+
ggtitle('Fitting score of models with different DNA unwrapping modes from nucleosome')

# +scale_x_log10(breaks=c(1,2,4,8,16,32,64,128,256,512,1024))+scale_y_log10(breaks=c(2.5,5,10,20,40,80,160))+xlab('Number of nodes (16 cores each)')+ylab('Productivity, ns/day')+scale_color_manual(values=c('blue','dark green','orange','red'),name='',breaks=c('ns_per_day','ns_per_day_impi','linear_mvapich2','linear_impi'),labels=c('Stampede MVAPICH2','Stampede IMPI','Linear MVAPICH2','Linear IMPI'))+ggtitle('Scaling on Stampede: MD simulations of nucleosome, 360K atoms')
ggsave(filename="fit.png",plot=p,width=10,height=5)
quit()
p=ggplot(data=s,aes(x=Number.of.nodes,y=value, color=variable))+geom_line()+geom_point()+scale_x_continuous(limits=c(0,256),breaks=c(1,2,4,8,16,32,64,128,256,512,1024))+scale_y_continuous(limits=c(0,160),breaks=c(2.5,5,10,20,40,80,160))+xlab('Number of nodes (16 cores each)')+ylab('Productivity, ns/day')+scale_color_manual(values=c('blue','dark green','orange','red'),name='',breaks=c('ns_per_day','ns_per_day_impi','linear_mvapich2','linear_impi'),labels=c('Stampede MVAPICH2','Stampede IMPI','Linear MVAPICH2','Linear IMPI'))+ggtitle('Scaling on Stampede: MD simulations of nucleosome, 360K atoms')
ggsave(filename="stampede_scaling_nonlog.png",plot=p,width=10,height=5)

scaling$percent_mvapich2=scaling$ns_per_day/scaling$linear_mvapich2*100
scaling$percent_impi=scaling$ns_per_day_impi/scaling$linear_impi*100
print(scaling)
sp=melt(scaling[,c('Number.of.nodes','percent_mvapich2','percent_impi')], id=c('Number.of.nodes'))

p=ggplot(data=sp,aes(x=Number.of.nodes,y=value,color=variable))+geom_line()+geom_point()+
scale_x_log10(limits=c(1,128),breaks=c(1,2,4,8,16,32,64,128,256,512,1024))+scale_y_continuous(limits=c(0,130),breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130))+
xlab('Number of nodes (16 cores each)')+ylab('Scaling efficiency')+
# scale_color_manual(values=c('red'),name='',breaks=c(''),labels=c('Stampede','Linear'))+
ggtitle('Scaling on Stampede: MD simulations of nucleosome, 360K atoms')+
scale_color_manual(values=c('red2','cyan3'),name='MPI version',breaks=c('percent_mvapich2','percent_impi'),labels=c('MVAPICH2','IMPI'))
ggsave(filename="stampede_scaling_percent.png",plot=p,width=10,height=5)


scaling$speedup_mvapich2=scaling$ns_per_day/scaling$ns_per_day[1]
scaling$speedup_impi=scaling$ns_per_day_impi/scaling$ns_per_day_impi[1]

ss=melt(scaling[,c('Number.of.nodes','speedup_mvapich2','speedup_impi')], id=c('Number.of.nodes'))


p=ggplot(data=ss,aes(x=Number.of.nodes,y=value,color=variable))+geom_line()+geom_point()+
scale_x_log10(breaks=c(1,2,4,8,16,32,64,128,256,512,1024))+scale_y_continuous(limits=c(0,130),breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130))+
xlab('Number of nodes (16 cores each)')+ylab('Speedup')+
# scale_color_manual(values=c('red'),name='',breaks=c(''),labels=c('Stampede','Linear'))+
ggtitle('Scaling on Stampede: MD simulations of nucleosome, 360K atoms')+
scale_color_manual(values=c('red2','cyan3'),name='MPI version',breaks=c('speedup_mvapich2','speedup_impi'),labels=c('MVAPICH2','IMPI'))

ggsave(filename="stampede_scaling_speedup.png",plot=p,width=10,height=5)



quit()
# dna_prot<-subset(dna_prot,DNA_resid>-74&DNA_resid<74)
# dna_prot<-subset(dna_prot,(PROT_chain%in%c('CHA','CHE') & PROT_resid>36)|(PROT_chain%in%c('CHB','CHF') & PROT_resid>15)|(PROT_chain%in%c('CHC','CHG') & PROT_resid>11&PROT_resid<119)|(PROT_chain%in%c('CHD','CHH') & PROT_resid>20))
dna_prot$hist_part='Tail'
dna_prot[with(dna_prot,(PROT_chain%in%c('CHA','CHE') & PROT_resid>36)|(PROT_chain%in%c('CHB','CHF') & PROT_resid>15)|(PROT_chain%in%c('CHC','CHG') & PROT_resid>11&PROT_resid<119)|(PROT_chain%in%c('CHD','CHH') & PROT_resid>20)),'hist_part']='Core'
dna_prot$hist_part<-factor(dna_prot$hist_part,levels=c('Core','Tail'))

dna_prot_cryst$hist_part='Tail'
dna_prot_cryst[with(dna_prot_cryst,(PROT_chain%in%c('CHA','CHE') & PROT_resid>36)|(PROT_chain%in%c('CHB','CHF') & PROT_resid>15)|(PROT_chain%in%c('CHC','CHG') & PROT_resid>11&PROT_resid<119)|(PROT_chain%in%c('CHD','CHH') & PROT_resid>20)),'hist_part']='Core'
dna_prot_cryst$hist_part<-factor(dna_prot_cryst$hist_part,levels=c('Core','Tail'))


dna_prot<-subset(dna_prot,type=='Z')
dna_prot_cryst<-subset(dna_prot_cryst,type=='Z')
#Let's define levels order for better graphics
prot_rn_lev=c('ARG','LYS','THR','SER','TYR','HSE','GLN','ASN','VAL','ILE','ALA','GLY','PRO','PHE','LEU','GLU','ASP','MET','CYS')
# type_lev=c('IP','HB','VdW','IM','WM')

#Let's assign level order
# dna_prot_cryst$type<-factor(dna_prot_cryst$type,levels=type_lev)
# dna_prot$type<-factor(dna_prot$type,levels=type_lev)

dna_prot$DNA_part<-factor(dna_prot$DNA_part,levels=c('phosphate','sugar','base'))
dna_prot_cryst$DNA_part<-factor(dna_prot_cryst$DNA_part,levels=c('phosphate','sugar','base'))

dna_prot$groove<-factor(dna_prot$groove,levels=c('minor','major'))
dna_prot_cryst$groove<-factor(dna_prot_cryst$groove,levels=c('minor','major'))


dna_prot_cryst$PROT_resname<-factor(dna_prot_cryst$PROT_resname,levels=prot_rn_lev)
dna_prot$PROT_resname<-factor(dna_prot$PROT_resname,levels=prot_rn_lev)

dna_prot$PROT_part<-factor(dna_prot$PROT_part,levels=c('side chain','backbone'))
dna_prot_cryst$PROT_part<-factor(dna_prot_cryst$PROT_part,levels=c('side chain','backbone'))


#canonical 14 arginines
# resname ARG and segname CHA CHE and resid 83 63
# resname ARG and segname CHB CHF and resid 45
# resname ARG and segname CHC CHG and resid 42 77
# resname ARG and segname CHD CHH and resid 30
#(resname ARG and segname CHA CHE and resid 49)
#let's construct mask for them
arg_mask=with(dna_prot,((PROT_chain %in% c('CHA','CHE'))&(PROT_resid %in% c(83,63,49)))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid %in% c(45)))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid %in% c(42,77)))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid %in% c(30))))
arg_mask_cryst=with(dna_prot_cryst,((PROT_chain %in% c('CHA','CHE'))&(PROT_resid %in% c(83,63,49)))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid %in% c(45)))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid %in% c(42,77)))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid %in% c(30))))

####General statistics section:
theme_set(theme_gray(base_size = 15))

##########Histograms with DNA_part classification



####Histogram of contacts 
theme_set(theme_gray(base_size = 12))

c<-ggplot(data=dna_prot_cryst[dna_prot_cryst$DNA_part=='base',],aes(x=PROT_resname,fill=groove))+
geom_bar(aes(alpha=hist_part,y=..count..),position='stack')+
scale_fill_manual(values=c("blue", "purple",'grey', "dark green"),name='DNA groove')+
scale_alpha_manual(values=c(1.0,0.2.5,5,0.6),name='Histone part',breaks=c('Core','Tail'))+ylab('Number of contacts')+xlab('Protein residues')+
scale_color_manual(values=c("blue", "purple",'grey', "dark green"),legend=NULL)+scale_y_continuous(limits=c(0,100),breaks = round(seq(0,100, by = 25),1))
c<-c+ggtitle('Direct interactions with BASES of DNA, X-ray data')
ggsave(filename="../analysis_data/int_dna_base_prot_hist_Z_cryst_poster.png",plot=c,width=10,height=5)

# ##------MD

md<-ggplot(data=dna_prot[dna_prot$DNA_part=='base',],aes(x=PROT_resname,fill=groove,weight=av_num))+
geom_bar(aes(alpha=hist_part,y=..count..),position='stack')+
scale_fill_manual(values=c("blue", "purple",'grey', "dark green"),name='DNA groove')+
scale_alpha_manual(values=c(1.0,0.2.5,5,0.6),name='Histone part',breaks=c('Core','Tail'))+ylab('Number of contacts')+xlab('Protein residues')+
scale_color_manual(values=c("blue", "purple",'grey', "dark green"),legend=NULL)
md<-md+ggtitle('Direct interactions with BASES of DNA, SIMULATIONS')
# png(file="../analysis_data/int_dan_prot_cryst.png",width=2000,height=800)
# ggplot(data=dna_prot_cryst,aes(x=PROT_resname,fill=DNA_part,alpha=type))+geom_bar(aes(color=DNA_part,y=..count..),position='dodged')+scale_fill_manual(values=c("red", "blue", "green", "grey", "purple"))+scale_alpha_manual(values=c(1,0.2.5,5,0.2,0.1))+scale_color_manual(values=c("red", "blue", "green", "grey", "purple"))
q<-arrangeGrob(c,md)
ggsave(filename="../analysis_data/int_dna_base_prot_hist_Z_poster.png",plot=md,width=10,height=5)

ggsave(filename="../analysis_data/int_dna_base_prot_hist_Z_comb_poster.png",plot=q,width=10,height=5)



c<-ggplot(data=dna_prot_cryst,aes(x=PROT_resname))+
geom_bar(aes(alpha=hist_part,y=..count..),position='stack')+
scale_fill_manual(values=c("blue", "purple",'grey', "dark green"),name='DNA groove')+
scale_alpha_manual(values=c(1.0,0.2.5,5,0.6),name='Histone part',breaks=c('Core','Tail'))+ylab('Number of contacts')+xlab('Protein residues')+
scale_color_manual(values=c("blue", "purple",'grey', "dark green"),legend=NULL)+scale_y_continuous(limits=c(0,650),breaks = round(seq(0,600, by = 100),1))
c<-c+ggtitle('Interactions with DNA, X-ray data')
ggsave(filename="../analysis_data/int_dna_all_prot_hist_Z_cryst_poster.png",plot=c,width=10,height=5)

# ##------MD

md<-ggplot(data=dna_prot,aes(x=PROT_resname,weight=av_num))+
geom_bar(aes(alpha=hist_part,y=..count..),position='stack')+
scale_fill_manual(values=c("blue", "purple",'grey', "dark green"),name='DNA groove')+
scale_alpha_manual(values=c(1.0,0.2.5,5,0.6),name='Histone part',breaks=c('Core','Tail'))+ylab('Number of contacts')+xlab('Protein residues')+
scale_color_manual(values=c("blue", "purple",'grey', "dark green"),legend=NULL)
md<-md+ggtitle('Interactions with DNA, SIMULATIONS')
# png(file="../analysis_data/int_dan_prot_cryst.png",width=2000,height=800)
# ggplot(data=dna_prot_cryst,aes(x=PROT_resname,fill=DNA_part,alpha=type))+geom_bar(aes(color=DNA_part,y=..count..),position='dodged')+scale_fill_manual(values=c("red", "blue", "green", "grey", "purple"))+scale_alpha_manual(values=c(1,0.2.5,5,0.2,0.1))+scale_color_manual(values=c("red", "blue", "green", "grey", "purple"))
q<-arrangeGrob(c,md)
ggsave(filename="../analysis_data/int_dna_all_prot_hist_Z_poster.png",plot=md,width=10,height=5)

ggsave(filename="../analysis_data/int_dna_all_prot_hist_Z_comb_poster.png",plot=q,width=10,height=5)



quit()


##############################
#Now let's make profiles along DNA sequence

# theme_set(theme_grey(base_size=30))

# dnaseq<-"ATCAATATCCACCTGCAGATACTACCAAAAGTGTATTTGGAAACTGCTCCATCAAAAGGCATGTTCAGCTGGAATCCAGCTGAACATGCCTTTTGATGGAGCAGTTTCCAAATACACTTTTGGTAGTATCTGCAGGTGGATATTGAT"
# seqdf=data.frame(sequence=substring(dnaseq, seq(1,nchar(dnaseq)),seq(1,nchar(dnaseq))),X=seq(0,nchar(dnaseq)-1))
# seqdf_t<-seqdf
# seqdf_t$X<-seqdf_t$X-73.0
# seqdf_t$Y<-rep(0,length(seqdf_t$X))
# # #Let;s average interactions for every resid.

################Total number of interactions
q=seq(93,-93,-1)

# ###Crystal
# d1=ddply(dna_prot_cryst,c("DNA_chain","DNA_resid","type"),function(df) c(num=nrow(df)))
# d1i=d1[d1$DNA_chain=='CHI',]
# d1j=d1[d1$DNA_chain=='CHJ',]
# d1j$DNA_resid<-q[d1j$DNA_resid+74]
# #Let's add zeros also

# d1zt=data.frame(DNA_resid=seq(73,-73,-1),DNA_chain=c('CHI'),num=c(0))
# t=data.frame(type=c('VdW','WM','IP','IM','HB'))
# d1z=merge(d1zt,t)

# d1<-rbind(d1i,d1j,d1z)

# tot_cryst=ddply(d1,c("DNA_resid","type"),summarize,number=sum(num))

# c<-ggplot(data=tot_cryst[tot_cryst$type %in% c('VdW','WM'),]) + ggtitle("DNA-protein interactions in nucleosome: crystal") + 
# xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=type),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("blue", "dark green"))#+#+xlim(-70,70)
# # geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)#+ylim(12.5,5,50)
# o<-ggplot(data=tot_cryst[tot_cryst$type %in% c('IP','HB'),]) + ggtitle("DNA-protein interactions in nucleosome: crystal") + 
# xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=type),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("red", "purple", "grey","dark green"))

# # q<-arrangeGrob(c,o)
# # ggsave(filename="../analysis_data/int_dna_prot_prof.png",plot=q,width=12.5,5,height=5)


# ##MD----------
# d1=ddply(dna_prot,c("DNA_chain","DNA_resid","type"),function(df) c(num=sum(df$av_num)))
# d1i=d1[d1$DNA_chain=='CHI',]
# d1j=d1[d1$DNA_chain=='CHJ',]
# d1j$DNA_resid<-q[d1j$DNA_resid+74]
# d1<-rbind(d1i,d1j,d1z)
# tot_md=ddply(d1,c("DNA_resid","type"),summarize,number=sum(num))


# cmd<-ggplot(data=tot_md[tot_md$type %in% c('VdW','WM'),]) + ggtitle("DNA-protein interactions in nucleosome: MD") + 
# xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=type),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("blue", "dark green"))#+#+xlim(-70,70)
# # geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)#+ylim(12.5,5,50)
# omd<-ggplot(data=tot_md[tot_md$type %in% c('IP','HB','IM'),]) + ggtitle("DNA-protein interactions in nucleosome: MD") + 
# xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=type),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("red", "purple","grey", "dark green"))


# q<-arrangeGrob(c,cmd,o,omd,ncol=1)
# ggsave(filename="../analysis_data/int_dna_prot_gen_prof.png",plot=q,width=12.5,5,height=10)


################Profile of interactions with bases, sugars, phosphates from MD
theme_set(theme_bw(base_size = 18))


d1=ddply(dna_prot,c("DNA_chain","DNA_resid","DNA_part"),function(df) c(num=sum(df$av_num)))
#Add zero frames

#let's invert chain J
d1[d1$DNA_chain=='CHJ','DNA_resid']=(d1[d1$DNA_chain=='CHJ','DNA_resid'])*-1
t=data.frame(DNA_part=c('base','sugar','phosphate'))
d1z=data.frame(DNA_resid=seq(93,-93,-1),DNA_chain=c('CHI'),num=c(0))

d1zp=merge(d1z,t)
d2zp=d1zp
d2zp[,'DNA_chain']=factor(c('CHJ'))
d1<-rbind(d1,d1zp,d2zp)
int_md=ddply(d1,c("DNA_resid",'DNA_part','DNA_chain'),summarize,number=sum(num))
int_md_all=ddply(d1,c("DNA_resid",'DNA_chain'),summarize,number=sum(num))



int_md_vse=int_md
#CHAIN I
#By int type
chi_part<-ggplot(data=int_md) + ggtitle("DNA-protein contacts in nucleosome along DNA, SIMULATIONS")  +
xlab("Super helix location")+ylab("Number of contacts")+geom_line(aes(x=DNA_resid/10,y=number,color=DNA_chain),size=1)+scale_y_continuous(breaks = round(seq(0,20, by = 5),1))+
scale_x_continuous(breaks = round(seq(-10,10, by = 1),1))+scale_color_manual(values=c("red", "blue"),name='DNA strand',breaks=c("CHI",'CHJ'),labels=c("Stand I","Strand J"))+
facet_grid(DNA_part~.,scales='free',space='free')



ggsave(filename="../analysis_data/int_dna_prot_Z_poster.png",plot=chi_part,width=17,height=7)


#The same for crystal


d1=ddply(dna_prot_cryst,c("DNA_chain","DNA_resid","DNA_part"),function(df) c(num=length(df$param1)))
#Add zero frames

#let's invert chain J
d1[d1$DNA_chain=='CHJ','DNA_resid']=(d1[d1$DNA_chain=='CHJ','DNA_resid'])*-1
t=data.frame(DNA_part=c('base','sugar','phosphate'))
d1z=data.frame(DNA_resid=seq(93,-93,-1),DNA_chain=c('CHI'),num=c(0))

d1zp=merge(d1z,t)
d2zp=d1zp
d2zp[,'DNA_chain']=factor(c('CHJ'))
d1<-rbind(d1,d1zp,d2zp)
int_md=ddply(d1,c("DNA_resid",'DNA_part','DNA_chain'),summarize,number=sum(num))
int_md_all=ddply(d1,c("DNA_resid",'DNA_chain'),summarize,number=sum(num))



int_md_vse=int_md
#CHAIN I
#By int type
chi_part<-ggplot(data=int_md) + ggtitle("DNA-protein contacts in nucleosome along DNA, NCP X-ray structure")  +
xlab("Super helix location")+ylab("Number of contacts")+geom_line(aes(x=DNA_resid/10,y=number,color=DNA_chain),size=1)+scale_y_continuous(breaks = round(seq(0,20, by = 5),1))+
scale_x_continuous(breaks = round(seq(-10,10, by = 1),1))+scale_color_manual(values=c("red", "blue"),name='DNA strand',breaks=c("CHI",'CHJ'),labels=c("Stand I","Strand J"))+
facet_grid(DNA_part~.,scales='free',space='free')


ggsave(filename="../analysis_data/int_dna_prot_Z_cryst_poster.png",plot=chi_part,width=17,height=7)


quit()

#Comparative
int_c$data=factor('crystal')
int_c_all$data=factor('crystal')
int_md$data=factor('MD')
int_md_all$data=factor('MD')

int=rbind(int_c,int_md)
int_all=rbind(int_c_all,int_md_all)

chi_part<-ggplot(data=int[ (int$DNA_chain=='CHI'),]) + ggtitle("DNA-protein interactions in nucleosome, by DNA_part, cryst+MD, chain I, Z-contacts") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,color=data,y=number),size=1)+
scale_x_continuous(breaks = round(seq(-100,100, by = 10),1))+scale_color_manual(values=c("blue", "dark green"))+
facet_grid(DNA_part~.,scales='free')

chi_all<-ggplot(data=int_all[ (int_all$DNA_chain=='CHI'),]) + ggtitle("DNA-protein interactions in nucleosome, all, cryst+MD, chain I, Z-contacts") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,color=data,y=number),size=1)+
scale_x_continuous(breaks = round(seq(-100,100, by = 10),1))+scale_color_manual(values=c("blue", "dark green"))

q<-arrangeGrob(chi_part,chi_all,ncol=1)

ggsave(filename="../analysis_data/int_dna_prot_Z_comb.png",plot=q,width=12.5,5,height=10)


#Let's do all the same for arginines


d1=ddply(subset(dna_prot,PROT_resname=='ARG'),c("DNA_chain","DNA_resid","DNA_part"),function(df) c(num=sum(df$av_num)))
#Add zero frames
t=data.frame(DNA_part=c('base','sugar','phosphate'))
d1z=data.frame(DNA_resid=seq(93,-93,-1),DNA_chain=c('CHI'),num=c(0))

d1zp=merge(d1z,t)
d2zp=d1zp
d2zp[,'DNA_chain']=factor(c('CHJ'))
d1<-rbind(d1,d1zp,d2zp)
int_md=ddply(d1,c("DNA_resid",'DNA_part','DNA_chain'),summarize,number=sum(num))
int_md_all=ddply(d1,c("DNA_resid",'DNA_chain'),summarize,number=sum(num))
int_md_ARG=int_md

#The same for crystal

d1=ddply(subset(dna_prot_cryst,PROT_resname=='ARG'),c("DNA_chain","DNA_resid","DNA_part"),function(df) c(num=length(df$param1)))
#Add zero frames
t=data.frame(DNA_part=c('base','sugar','phosphate'))
d1z=data.frame(DNA_resid=seq(93,-93,-1),DNA_chain=c('CHI'),num=c(0))

d1zp=merge(d1z,t)
d2zp=d1zp
d2zp[,'DNA_chain']=factor(c('CHJ'))
d1<-rbind(d1,d1zp,d2zp)
int_c=ddply(d1,c("DNA_resid",'DNA_part','DNA_chain'),summarize,number=sum(num))
int_c_all=ddply(d1,c("DNA_resid",'DNA_chain'),summarize,number=sum(num))
int_c_ARG=int_c

#Comparative
int_c$data=factor('crystal')
int_c_all$data=factor('crystal')
int_md$data=factor('MD')
int_md_all$data=factor('MD')

int=rbind(int_c,int_md)
int_all=rbind(int_c_all,int_md_all)

chi_part<-ggplot(data=int[ (int$DNA_chain=='CHI'),]) + ggtitle("DNA-ARG interactions in nucleosome, by DNA_part, cryst+MD, chain I, Z-contacts") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,color=data,y=number),size=1)+
scale_x_continuous(breaks = round(seq(-100,100, by = 10),1))+scale_color_manual(values=c("blue", "dark green"))+
facet_grid(DNA_part~.,scales='free')

chi_all<-ggplot(data=int_all[ (int_all$DNA_chain=='CHI'),]) + ggtitle("DNA-ARG interactions in nucleosome, all, cryst+MD, chain I, Z-contacts") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,color=data,y=number),size=1)+
scale_x_continuous(breaks = round(seq(-100,100, by = 10),1))+scale_color_manual(values=c("blue", "dark green"))

q<-arrangeGrob(chi_part,chi_all,ncol=1)

ggsave(filename="../analysis_data/int_dna_prot_ARG_Z_comb.png",plot=q,width=12.5,5,height=10)

#Let's do all the same for lys


d1=ddply(subset(dna_prot,PROT_resname=='LYS'),c("DNA_chain","DNA_resid","DNA_part"),function(df) c(num=sum(df$av_num)))
#Add zero frames
t=data.frame(DNA_part=c('base','sugar','phosphate'))
d1z=data.frame(DNA_resid=seq(93,-93,-1),DNA_chain=c('CHI'),num=c(0))

d1zp=merge(d1z,t)
d2zp=d1zp
d2zp[,'DNA_chain']=factor(c('CHJ'))
d1<-rbind(d1,d1zp,d2zp)
int_md=ddply(d1,c("DNA_resid",'DNA_part','DNA_chain'),summarize,number=sum(num))
int_md_all=ddply(d1,c("DNA_resid",'DNA_chain'),summarize,number=sum(num))

#The same for crystal

d1=ddply(subset(dna_prot_cryst,PROT_resname=='LYS'),c("DNA_chain","DNA_resid","DNA_part"),function(df) c(num=length(df$param1)))
#Add zero frames
t=data.frame(DNA_part=c('base','sugar','phosphate'))
d1z=data.frame(DNA_resid=seq(93,-93,-1),DNA_chain=c('CHI'),num=c(0))

d1zp=merge(d1z,t)
d2zp=d1zp
d2zp[,'DNA_chain']=factor(c('CHJ'))
d1<-rbind(d1,d1zp,d2zp)
int_c=ddply(d1,c("DNA_resid",'DNA_part','DNA_chain'),summarize,number=sum(num))
int_c_all=ddply(d1,c("DNA_resid",'DNA_chain'),summarize,number=sum(num))

#Comparative
int_c$data=factor('crystal')
int_c_all$data=factor('crystal')
int_md$data=factor('MD')
int_md_all$data=factor('MD')

int=rbind(int_c,int_md)
int_all=rbind(int_c_all,int_md_all)

chi_part<-ggplot(data=int) + ggtitle("DNA-LYS interactions in nucleosome, by DNA_part, cryst+MD, chain I, Z-contacts") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,color=data,y=number,linetype=DNA_chain),size=1)+
scale_x_continuous(breaks = round(seq(-100,100, by = 10),1))+scale_color_manual(values=c("blue", "dark green"))+
facet_grid(DNA_part~.,scales='free')

chi_all<-ggplot(data=int_all[ (int_all$DNA_chain=='CHI'),]) + ggtitle("DNA-LYS interactions in nucleosome, all, cryst+MD, chain I, Z-contacts") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,color=data,y=number),size=1)+
scale_x_continuous(breaks = round(seq(-100,100, by = 10),1))+scale_color_manual(values=c("blue", "dark green"))

q<-arrangeGrob(chi_part,chi_all,ncol=1)

ggsave(filename="../analysis_data/int_dna_prot_LYS_Z_comb.png",plot=q,width=12.5,5,height=10)


####Let's subdivide all arg and canonical arg


d1=ddply(subset(dna_prot,arg_mask),c("DNA_chain","DNA_resid","DNA_part"),function(df) c(num=sum(df$av_num)))
#Add zero frames
t=data.frame(DNA_part=c('base','sugar','phosphate'))
d1z=data.frame(DNA_resid=seq(93,-93,-1),DNA_chain=c('CHI'),num=c(0))

d1zp=merge(d1z,t)
d2zp=d1zp
d2zp[,'DNA_chain']=factor(c('CHJ'))
d1<-rbind(d1,d1zp,d2zp)
int_md=ddply(d1,c("DNA_resid",'DNA_part','DNA_chain'),summarize,number=sum(num))
int_md_all=ddply(d1,c("DNA_resid",'DNA_chain'),summarize,number=sum(num))
int_md_ARG14=int_md

int_md_vse$data=factor('all')
int_md_ARG$data=factor('ARG')
int_md_ARG14$data=factor('ARG14')
int=rbind(int_md_vse,int_md_ARG,int_md_ARG14)

q<-ggplot(data=int[(int$DNA_chain %in% c('CHI')),]) + ggtitle("MD: DNA-protein interactions, chain I, Z-contacts") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,color=data,y=number,linetype=DNA_chain),size=1)+
scale_x_continuous(breaks = round(seq(-100,100, by = 10),1))+scale_color_manual(values=c("blue", "green",'red'))+
facet_grid(DNA_part~.,scales='free')+theme(panel.grid.major = element_line(colour = "black"),panel.grid.minor = element_blank())

int$CF=factor(paste0(int$DNA_chain,sep="_",int$DNA_part,sep="_",int$data))
icsv=dcast(DNA_resid~CF,data=int,value.var='number')

write.csv(icsv,file="../analysis_data/int_dna_prot_Z_arg_comb.csv")
ggsave(filename="../analysis_data/int_dna_prot_Z_arg_comb.png",plot=q,width=12.5,5,height=10)

#The same for crystal



d1=ddply(subset(dna_prot_cryst,arg_mask_cryst),c("DNA_chain","DNA_resid","DNA_part"),function(df) c(num=length(df$param1)))
#Add zero frames
t=data.frame(DNA_part=c('base','sugar','phosphate'))
d1z=data.frame(DNA_resid=seq(93,-93,-1),DNA_chain=c('CHI'),num=c(0))

d1zp=merge(d1z,t)
d2zp=d1zp
d2zp[,'DNA_chain']=factor(c('CHJ'))
d1<-rbind(d1,d1zp,d2zp)
int_c=ddply(d1,c("DNA_resid",'DNA_part','DNA_chain'),summarize,number=sum(num))
int_c_all=ddply(d1,c("DNA_resid",'DNA_chain'),summarize,number=sum(num))
int_c_ARG14=int_c

int_c_vse$data=factor('all')
int_c_ARG$data=factor('ARG')
int_c_ARG14$data=factor('ARG14')
int=rbind(int_c_vse,int_c_ARG,int_c_ARG14)

q<-ggplot(data=int[ (int$DNA_chain %in% c('CHI')),]) + ggtitle("X-ray: DNA-protein interactions, chain I, Z-contacts") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,color=data,y=number,linetype=DNA_chain),size=1)+
scale_x_continuous(breaks = round(seq(-100,100, by = 10),1))+scale_color_manual(values=c("blue", "green",'red'))+
facet_grid(DNA_part~.,scales='free')+theme(panel.grid.major = element_line(colour = "black"),panel.grid.minor = element_blank())

int$CF=factor(paste0(int$DNA_chain,sep="_",int$DNA_part,sep="_",int$data))
icsv=dcast(DNA_resid~CF,data=int,value.var='number')
write.csv(icsv,file="../analysis_data/int_dna_prot_Z_arg_comb_cryst.csv")

ggsave(filename="../analysis_data/int_dna_prot_Z_arg_comb_cryst.png",plot=q,width=12.5,5,height=10)



quit()
#By contact type

d1=ddply(subset(dna_prot,type=='VdW'),c("DNA_chain","DNA_resid","C_type","DNA_part"),function(df) c(num=sum(df$av_num)))
#Add zero frames
d1zt=data.frame(DNA_resid=seq(73,-73,-1),DNA_chain=c('CHI'),num=c(0))
t=data.frame(DNA_part=c('base','sugar','phosphate'))
d1zp=merge(d1zt,t)
t2=data.frame(C_type=c('PN','PP','NN'))
d1z=merge(d1zp,t2)

d1<-rbind(d1,d1z)
int_md=ddply(d1,c("DNA_resid","DNA_chain",'DNA_part','C_type'),summarize,number=sum(num))

chi_dna_part_c_type<-ggplot(data=int_md[(int_md$DNA_chain=='CHI'),]) + ggtitle("DNA-protein interactions in nucleosome, MD, chain I, contacts by type") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=C_type),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
facet_grid(DNA_part~.,scales='free')
ggsave(filename="../analysis_data/int_dna_prot_dp_chi_c_type.png",plot=chi_dna_part_c_type,width=12.5,5,height=10)

quit()

####Let's specifically study ARGinines

#MD
d1=ddply(subset(dna_prot,PROT_resname=='ARG' & type %in% c('IP','HB','VdW') ),c('DNA_chain',"DNA_resid","DNA_part"),function(df) c(num=sum(df$av_num)))
#Add zero frames
d1zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHI'))
d2zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHJ'))
t=data.frame(DNA_part=c('base','sugar','phosphate'))
dzt<-rbind(d1zt,d2zt)
dz=merge(dzt,t)
d1=rbind(d1,dz)

int_md=ddply(d1,c('DNA_chain',"DNA_resid",'DNA_part'),summarize,number=sum(num))

chi_dna_part_arg_int<-ggplot(data=int_md) + ggtitle("DNA-protein interactions in nucleosome, MD, interactions with ARG (VdW+IP+HB)") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=DNA_chain),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
facet_grid(DNA_part~.,scales='free')
ggsave(filename="../analysis_data/int_dna_prot_dp_arg.png",plot=chi_dna_part_arg_int,width=12.5,5,height=10)

#Crystal
d1=ddply(subset(dna_prot_cryst,PROT_resname=='ARG' & type %in% c('IP','HB','VdW') ),c('DNA_chain',"DNA_resid","DNA_part"),function(df) c(num=length(df$param1)))
#Add zero frames
d1zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHI'))
d2zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHJ'))
t=data.frame(DNA_part=c('base','sugar','phosphate'))
dzt<-rbind(d1zt,d2zt)
dz=merge(dzt,t)
d1=rbind(d1,dz)

int_md=ddply(d1,c('DNA_chain',"DNA_resid",'DNA_part'),summarize,number=sum(num))

chi_dna_part_arg_int_cryst<-ggplot(data=int_md) + ggtitle("DNA-protein interactions in nucleosome, Crystal, interactions with ARG (VdW+IP+HB)") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=DNA_chain),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
facet_grid(DNA_part~.,scales='free')
ggsave(filename="../analysis_data/int_dna_prot_dp_arg_cryst.png",plot=chi_dna_part_arg_int_cryst,width=12.5,5,height=10)


####Let's specifically study LYSinines


d1=ddply(subset(dna_prot,PROT_resname=='LYS' & type %in% c('IP','HB','VdW') ),c('DNA_chain',"DNA_resid","DNA_part"),function(df) c(num=sum(df$av_num)))
#Add zero frames
d1zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHI'))
d2zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHJ'))
t=data.frame(DNA_part=c('base','sugar','phosphate'))
dzt<-rbind(d1zt,d2zt)
dz=merge(dzt,t)
d1=rbind(d1,dz)

int_md=ddply(d1,c('DNA_chain',"DNA_resid",'DNA_part'),summarize,number=sum(num))

chi_dna_part_lys_int<-ggplot(data=int_md) + ggtitle("DNA-protein interactions in nucleosome, MD, interactions with LYS (VdW+IP+HB)") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=DNA_chain),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
facet_grid(DNA_part~.,scales='free')
ggsave(filename="../analysis_data/int_dna_prot_dp_lys.png",plot=chi_dna_part_lys_int,width=12.5,5,height=10)



####Misc


####Let's specifically study ARGinines

#MD
d1=ddply(subset(dna_prot,PROT_resname=='ARG' & type %in% c('VdW') ),c('DNA_chain',"DNA_resid","DNA_part"),function(df) c(num=sum(df$av_num)))
#Add zero frames
d1zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHI'))
d2zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHJ'))
t=data.frame(DNA_part=c('base','sugar','phosphate'))
dzt<-rbind(d1zt,d2zt)
dz=merge(dzt,t)
d1=rbind(d1,dz)

int_md=ddply(d1,c('DNA_chain',"DNA_resid",'DNA_part'),summarize,number=sum(num))

chi_dna_part_arg_int<-ggplot(data=int_md) + ggtitle("DNA-protein interactions in nucleosome, MD, interactions with ARG (contacts only)") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=DNA_chain),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
facet_grid(DNA_part~.,scales='free')

#Crystal
d1=ddply(subset(dna_prot_cryst,PROT_resname=='ARG' & type %in% c('VdW') ),c('DNA_chain',"DNA_resid","DNA_part"),function(df) c(num=length(df$param1)))
#Add zero frames
d1zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHI'))
d2zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHJ'))
t=data.frame(DNA_part=c('base','sugar','phosphate'))
dzt<-rbind(d1zt,d2zt)
dz=merge(dzt,t)
d1=rbind(d1,dz)

int_md=ddply(d1,c('DNA_chain',"DNA_resid",'DNA_part'),summarize,number=sum(num))

chi_dna_part_arg_int_cryst<-ggplot(data=int_md) + ggtitle("DNA-protein interactions in nucleosome, Crystal, interactions with ARG (contacts only)") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=DNA_chain),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
facet_grid(DNA_part~.,scales='free')
q<-arrangeGrob(chi_dna_part_arg_int,chi_dna_part_arg_int_cryst,ncol=1)

ggsave(filename="../analysis_data/int_dna_prot_dp_arg_md_vs_cryst.png",plot=q,width=12.5,5,height=20)


