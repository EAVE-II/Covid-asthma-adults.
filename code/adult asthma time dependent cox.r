##########################################################
# Name of file: dataset for time dependent cox.R
# Data release (if applicable):
# Original author(s): Chris Robertson chrisobertson@nhs.net
# Original date: 06 August 2020
# Latest update author (if not using version control) - Jiafeng Pan jiafeng.pan@phs.scot
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: reads in the cohort df.org2 and set up dataset for time dependent cox model.
#run cox model
# Approximate run time: Unknown
##########################################################

library(survival)

df <- df.org2
a.start <- as.Date("2020-03-01")
#time to variable already adjusted for death
z.rv <- "hosp_covid"
z.rv.time <- "Time.To.Hosp"

z.rv <- "icu_death"
z.rv.time <- "Time.To.ICU.Death"

z.rv <- "death_covid"
z.rv.time <- "Time.To.Death"

event_ids <- filter(df, df[,z.rv]==1) %>% pull(EAVE_LINKNO) %>% unique()
set.seed(21021854)
no_event_ids <- filter(df, df[,z.rv]==0) %>% slice_sample(n=100*length(event_ids)) %>% pull(EAVE_LINKNO) %>% unique()
z_df <- df %>% filter(EAVE_LINKNO %in% c(event_ids, no_event_ids))
wt <- nrow(filter(df, df[,z.rv]==0))/length(no_event_ids)
z_df$wt2 <- ifelse(z_df[,z.rv]==1,1,wt)
z_df$wt <- z_df$eave_weight*z_df$wt2

z_df$end <- a.start+z_df[,z.rv.time]
z_df <- subset(z_df, z_df[,z.rv.time]>=0)

z_df$dv1 <- z_df$date_vacc_1
z_df$dv2 <- z_df$date_vacc_1+28
z_df$dv3 <- z_df$date_vacc_2
z_df$dv4 <- z_df$date_vacc_2+28

z_df$dv2 <- ifelse(z_df$dv2>z_df$dv3, as.character(z_df$dv3),as.character(z_df$dv2))
z_df$dv2 <- as.Date(z_df$dv2)
#subset(z_df, as.numeric(date_vacc_2-date_vacc_1)<28)
#subset(z_df, EAVE_LINKNO=="EAVE0572605")

  
z.name<-c(names(z_df)[grepl("dv",names(z_df))])

z <- z_df %>% 
  dplyr::select(EAVE_LINKNO, all_of(z.name), "end", all_of(z.rv)) %>%
  mutate(start=a.start)
z$event = z[,z.rv]

#z <- subset(z, start<end)
z <- z %>% mutate_at(vars(start,end,dv1,dv2,dv3,dv4), funs(as.numeric))
  
z1 <- tmerge(z,z,id=EAVE_LINKNO, endpt = event(end, event), tstart=start-1, tstop=end)
z1 <- tmerge(z1,z, id=EAVE_LINKNO, per1=tdc(dv1))
z1 <- tmerge(z1,z, id=EAVE_LINKNO, per2=tdc(dv2))
z1 <- tmerge(z1,z, id=EAVE_LINKNO, per3=tdc(dv3))
z1 <- tmerge(z1,z, id=EAVE_LINKNO, per4=tdc(dv4)) 
tfun <- function(x)as.Date(x,origin="1970-01-01")
z1 <- z1 %>%  
  mutate_at(vars(start,end,dv1,dv2,dv3,dv4), funs(tfun))
z1 <- z1 %>%  
  mutate_at(vars(tstart,tstop), funs(tfun))

  
z1 <- z1 %>% mutate(vacc_status=per1+per2+per3+per4) 
table(z1$vacc_status,exclude=NULL)
z1$vacc_status <- factor(z1$vacc_status, levels=0:4, labels=c("Unvacc","Dose1_0_27","Dose1_28","Dose2_0_27","Dose2_28"))




zz <- left_join(z1, z_df[,c("EAVE_LINKNO","Sex", "ageYear", "simd2020_sc_quintile", "no_riskgp_gp" ,"no_hosp_2yrgp",
                        "bmigp","wt","asthma_hosp_1yr","asthma_hosp_2yr","asthma_los_1yrgp","asthma_los_2yrgp",
                        "asthma_pres_1yrgp1","asthma_pres_2yrgp1", "asthma")])

zz$start <- as.numeric(zz$tstart-a.start)+2
zz$stop <- as.numeric(zz$tstop-a.start)+2
levels(zz$simd2020_sc_quintile)[6] <- NA
#cox model

fun.cox <- function(z.var,z) {
  #z is a cox model fitted
  z.coefs <- summary(z)$coefficients
  z.est <- z.coefs[,"coef"]
  z.se <- z.coefs[,"se(coef)"]
  z.p <- z.coefs[,"p"]
  z.out <- cbind.data.frame(levels=dimnames(z.coefs)[[1]],HR=exp(z.est),LCL=exp(z.est-1.96*z.se),
                            UCL=exp(z.est+1.96*z.se),p=z.p,outcome=z.var)
  z.out$hrci <- paste(round(z.out$HR,2), " (", round(z.out$LCL,2), "-",round(z.out$UCL,2),")", sep="")
  
  return(z.out)
}

#z_exp.list <- c("asthma_hosp_1yr","asthma_hosp_2yr","asthma_los_1yrgp", "asthma_los_2yrgp",
#                "asthma_pres_1yrgp1","asthma_pres_2yrgp1")
z_exp.list <-"asthma"

i=1
for(z_exp in z_exp.list){
z.fmla <- as.formula(paste("Surv(start,stop,",z.rv,") ~  pspline(ageYear) + Sex + simd2020_sc_quintile + no_hosp_2yrgp + no_riskgp_gp + bmigp + vacc_status +",
                           z_exp))

z.fit <- coxph(z.fmla , data=zz, weights = wt)
#z.fit <- coxph(z.fmla , data=df, weights = eave_weight)
#summary(z.fit)

z.r <- fun.cox(z.rv,z.fit)
z.r$exp <- z_exp
if(i==1) z.rr <- z.r else
  z.rr <- rbind(z.rr,z.r)
i=i+1
print(i)
}

write.csv(z.rr,paste(z.rv,"asthma.adults2.csv",sep=""))
  
  





