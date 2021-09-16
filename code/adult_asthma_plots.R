##########################################################
# Name of file: adult_asthma_plots.R
# Data release (if applicable):
# Original author(s): Chris Robertson chrisobertson@nhs.net
# Original date: 06 August 2020
# Latest update author (if not using version control) - Jiafeng Pan  jiafeng.pan@phs.scot
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: Run adult_asthma_cohort.R to get the cohort 
#                         runs through graphs
# Approximate run time: Unknown
##########################################################

# 01 Setup ####
#Libraries
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(survival)
library(survminer)
#Load data
#Location <- "/conf/"  # Server
Location <- "//isdsf00d03/"  # Desktop

#use df.org2 from 08_01a_children_cohort_asthma.r

setwd("Z:\\JP\\vaccine safety\\covid risk children results")



#z.xlab <- paste0("Days from ", format.Date(a_begin, "%d-%b-%Y"))
z.xlab <- "Date"
cutpt1 <- as.numeric(as.Date("2020-07-31")-as.Date("2020-03-01"))
cutpt2 <- as.numeric(as.Date("2021-04-30")-as.Date("2020-03-01"))
cutpt3 <- as.numeric(as.Date("2021-01-04")-as.Date("2020-03-01"))
cutpt4 <- as.numeric(as.Date("2021-05-17")-as.Date("2020-03-01"))
cutpt5 <- as.numeric(as.Date("2020-12-08")-as.Date("2020-03-01"))
#you can use this to speed up the model fitting if needed
#select all the events and a random sample of 100 controls per event
#change df to z.df in the coxph fitting below
#z.case <- filter(df, hosp_covid==1)
#z.control <- slice_sample(filter(df, hosp_covid==0), n=nrow(z.case)*100)
#z.df <- rbind(z.case,z.control)

#survival curve + hazard ratio

df <- df.org2


#splots <- list()
#k=1
#z.var="asthma_hosp_2yr"
#z.var="asthma_pres_2yrgp1"
#z.var="asthma_hosp_1yr"
#z.var="asthma_pres_1yrgp1"
z.timetoevent.list <- c("Time.To.Hosp", "Time.To.ICU.Death","Time.To.Death")
z.outcome.list <- c("hosp_covid","icu_death","death_covid")
z.label.outcome.list <- c("Covid Hospitalisation","Covid ICU or Death","Covid Death")
z.var.list <- c("asthma_hosp_2yr","asthma_pres_2yrgp1","asthma_hosp_1yr","asthma_pres_1yrgp1")
for(j in 1:3){
  z.timetoevent <- z.timetoevent.list[j]
  z.outcome <- z.outcome.list[j]
  z.label.outcome <- z.label.outcome.list[j]
for(i in 1:4){
  z.var <- z.var.list[i]
z <- survfit(Surv(get(z.timetoevent), get(z.outcome)) ~ get(z.var),weights = eave_weight, data=df)
z.ymax <- 1-min(z$surv)
#z.time <- Sys.time()
#z.fit <- coxph(Surv(Time.To.Hosp, hosp_covid) ~  pspline(ageYear) + Sex + simd2020_sc_quintile + ur6_2016_name + get(z.var),weights = eave_weight,data=df)
#print(Sys.time() - z.time)
#z.tt <- summary(z.fit)$conf.int
#z.t <- round(z.tt[nrow(z.tt),c(1,3,4)],2)
#z.t2 <- round(z.tt[nrow(z.tt)-1,c(1,3,4)],2)
#z.t2 <- rbind(c(1,NA,NA), z.t)
#dimnames(z.t2) <- list(NULL, c("HR.Death","LCL.Death","UCL.Death"))
#z.tab <- cbind(z.tab,z.t2)
#ym[k] <- z.ymax

g1 <-ggsurvplot(z, data = df, fun="event", risk.table = FALSE, legend.labs=levels(as.factor(df[,z.var])), censor=FALSE) +
  labs(x=z.xlab, y=paste0("Proportion ",z.label.outcome), title=paste0("Time to ", z.label.outcome," - ",z.var)) 
#g1$plot <- g1$plot +ggplot2::annotate("text",x=50,y=z.ymax*0.9, label = paste0("HR = ", z.t[1], "\n(",z.t[2],", ",z.t[3],")" ))
# g1$plot <- g1$plot +ggplot2::annotate("text",x=50,y=z.ymax*0.9, label = paste0("HR1 = ", z.t[1], "(",z.t[2],", ",z.t[3],")" ,"\nHR2 = ", z.t2[1], "(",z.t2[2],", ",z.t2[3],")"))
g1$plot <- g1$plot +ggplot2::annotate("text",x=50,y=z.ymax*0.5, label = "1st wave")
g1$plot <- g1$plot +ggplot2::annotate("text",x=cutpt1+120,y=z.ymax*0.5, label = "2nd wave")
g1$plot <- g1$plot +ggplot2::annotate("text",x=cutpt2+50,y=z.ymax*0.5, label = "3rd wave")
g1$plot <-g1$plot + geom_vline(xintercept = c(cutpt1,cutpt2),colour="grey", linetype = "longdash")
g1$plot <-g1$plot + geom_rect(aes(xmin=cutpt3, xmax=cutpt4, ymin=0, ymax=z.ymax),alpha=0.005,fill=5)
g1$plot <-g1$plot + geom_rect(aes(xmin=cutpt4, xmax=510, ymin=0, ymax=z.ymax),alpha=0.005,fill=7)
g1$plot <- g1$plot +ggplot2::annotate("text",x=cutpt3+50,y=z.ymax*0.9, label = "alpha\n dominant",col="dark green")
g1$plot <- g1$plot +ggplot2::annotate("text",x=cutpt4+30,y=z.ymax*0.9, label = "delta\n dominant", col="dark orange")
g1$plot <-g1$plot + geom_vline(xintercept = cutpt5,colour="red", linetype = "longdash")
g1$plot <- g1$plot +ggplot2::annotate("text",x=cutpt5-30,y=z.ymax*0.9, label = "vaccination\n starts", col="red")
g1$plot <- g1$plot + scale_x_continuous(breaks = seq(0,450,by=150), label=as.character(as.Date("2020-03-01")+seq(0,450,by=150)))
#g1
pdf(paste(z.outcome,"cumprop_events_over_time",z.var,".pdf",sep="_"), width = 8, height = 6)
#pdf("cumprop events over time 5to17years.pdf", width = 8, height = 6)
print(g1)
dev.off()
}
}


