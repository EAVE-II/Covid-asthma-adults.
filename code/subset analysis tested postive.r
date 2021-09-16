##########################################################
# Name of file: subset analysis tested positive.R
# Data release (if applicable):
# Original author(s): Jiafeng Pan jiafeng.pan@phs.scot
# Original date: 30 August 2020
# Latest update author (if not using version control) - Jiafeng Pan jiafeng.pan@phs.scot
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: model the risk for those with severe asthma among those tested postive for covid
# Approximate run time: Unknown
##########################################################
#for asthma severity analysis
#only for thoes tested positive
#look back 2 years from the date of test
#cox model using date of test as starting point
#add a spline of time between the date of test and 01-03-2020
#measure the vaccine status at the date of test


#use df.org from adult_asthma_cohort.r
#For rhose with no positive test give them an impute specimen date 7 days before the date of admission
#Some people are in hospital before the day of the test and get a negative time
#reset the time to hospitalisation to zero – remove the ones that are very negative


zz <- subset(df.org, result==1 |hosp_covid==1|icu_death==1|death_covid==1)
table(zz$result)
zz[is.na(zz$SpecimenDate), "SpecimenDate"] <- zz[is.na(zz$SpecimenDate), "date_hosp_covid"]-7
zz[is.na(zz$SpecimenDate), "SpecimenDate"] <- zz[is.na(zz$SpecimenDate), "date_icu_death"]-14
summary(zz$SpecimenDate)

#need to account for death here
zz$time.test.to.hosp <- zz$Time.To.Hosp + as.numeric(zz$SpecimenDate-as.Date("2020-03-01"))
summary(zz$time.test.to.hosp)
zz[zz$time.test.to.hosp <0,"time.test.to.hosp"] <- 0

zz$time.test.to.icudeath <- zz$Time.To.ICU.Death + as.numeric(zz$SpecimenDate-as.Date("2020-03-01"))
summary(zz$time.test.to.icudeath)
#zz[zz$time.test.to.icudeath <0,"time.test.to.icudeath"] <- 0

zz$time.test.to.death<- zz$Time.To.Death + as.numeric(zz$SpecimenDate-as.Date("2020-03-01"))
summary(zz$time.test.to.death)
#zz[zz$time.test.to.death <0,"time.test.to.death"] <- 0

df.org3 <- zz

#check vaccination status at the date of postive reuslts

zz$vacc_status <- ifelse(is.na(zz$date_vacc_1) | zz$SpecimenDate<zz$date_vacc_1, "Unvacc", NA)
zz$vacc_status <- ifelse(!is.na(zz$date_vacc_1) & zz$SpecimenDate>=zz$date_vacc_1 & zz$SpecimenDate<zz$date_vacc_1+28, "Dose1_0_27", zz$vacc_status)
zz$vacc_status <- ifelse(!is.na(zz$date_vacc_1) & zz$SpecimenDate>=zz$date_vacc_1+28 , "Dose1_28", zz$vacc_status)
zz$vacc_status <- ifelse(!is.na(zz$date_vacc_2) & zz$SpecimenDate>=zz$date_vacc_2 & zz$SpecimenDate<zz$date_vacc_2+28, "Dose2_0_27", zz$vacc_status)
zz$vacc_status <- ifelse(!is.na(zz$date_vacc_2) & zz$SpecimenDate>=zz$date_vacc_2+28 , "Dose2_28", zz$vacc_status)

zz$vacc_status <- factor(zz$vacc_status, levels=c("Unvacc","Dose1_0_27","Dose1_28","Dose2_0_27","Dose2_28"))

table(zz$vacc_status, exclude=NULL)
df.org3 <- zz




smr01_asthma <- readRDS(paste0(Location,"EAVE/GPanalysis/data/smr01_adult_asthma_admits.rds"))
z <- as.data.frame(smr01_asthma)
z$admit_date <- as.Date(substr(z$admit_date,1,10))
z$dis_date <- as.Date(substr(z$dis_date,1,10))
#table(z$main_diag_admit)
#all children with prior asthma hosp
dim(z)
z <- merge(z, df.org3[,c("EAVE_LINKNO", "SpecimenDate")])
zs <- subset(z,admit_date>=SpecimenDate-365.25*2 & admit_date<SpecimenDate)
zs$los_2yr <- as.numeric(zs$dis_date-zs$admit_date)
z.agg <- zs %>% group_by(EAVE_LINKNO) %>% 
  dplyr::summarise(admit_date_2yr=max(admit_date), los_2yr=sum(los_2yr))
zz <- left_join(df.org3, z.agg)
dim(zz)

zz$asthma_hosp_2yr <- ifelse(is.na(zz$admit_date_2yr)==FALSE, "prior hosp 2yr", as.character(zz$asthma))
table(zz$asthma_hosp_2yr)
zz$asthma_hosp_2yr <- factor(zz$asthma_hosp_2yr, level=c("No", "Yes", "prior hosp 2yr"))
levels(zz$asthma_hosp_2yr)[1:2] <- c("no asthma", "mild asthma")


zz$los_2yrgp <- cut(zz$los_2yr, breaks=c(-1,0,1,2,1000),label=c("0","1","2","3+"))
table(zz$los_2yrgp)
zz$asthma_los_2yr <- ifelse(is.na(zz$los_2yr)==FALSE, as.character(zz$los_2yrgp), as.character(zz$asthma))
table(zz$asthma_los_2yr)
zz$asthma_los_2yr <- factor(zz$asthma_los_2yr, level=c("No", "Yes", "0","1","2","3+"))
levels(zz$asthma_los_2yr)[1:2] <- c("no asthma", "mild asthma")


zz$asthma_los_2yrgp <- zz$asthma_los_2yr
levels(zz$asthma_los_2yrgp)[3:5] <- "0-2"


df.org3 <- zz

z <- readRDS("/conf/EAVE/GPanalysis/data/EAVEII_PIS_ASTHMA_2021-08-25.rds")
z <- as.data.frame(z)
dim(z)#145449
summary(z$prescribed_full_date)#
#z <- subset(z, `PI Approved Name` %in% c("PREDNISOLONE"))#
z <- merge(z, df.org3[,c("EAVE_LINKNO", "SpecimenDate")])
zs <-subset(z,prescribed_full_date>=SpecimenDate-365.25*2 & prescribed_full_date<SpecimenDate)
z.agg <- z %>% group_by(EAVE_LINKNO) %>% 
  dplyr::summarise( no_pres_2yr=length(EAVE_LINKNO))
table(z.agg$no_pres_2yr)

zz <- left_join(df.org3, z.agg)
dim(zz)
table(zz$no_pres)


zz$no_pres_2yrgp1 <- ifelse(zz$no_pres_2yr>=3,"3+",zz$no_pres_2yr)
#zz$no_pres_2yrgp2 <- ifelse(zz$no_pres_2yr>=2,"2+",zz$no_pres_2yr)

zz$asthma_pres_2yrgp1 <- ifelse(is.na(zz$no_pres_2yrgp1)==FALSE, zz$no_pres_2yrgp1, as.character(zz$asthma))
zz$asthma_pres_2yrgp1 <- factor(zz$asthma_pres_2yrgp1, level=c("No", "Yes", "1" ,"2" , "3+" ))
levels(zz$asthma_pres_2yrgp1)[1:2] <- c("no asthma", "0")
table(zz$asthma_pres_2yrgp1)

zz$time <- as.numeric(zz$SpecimenDate - as.Date("2020-03-01"))
summary(zz$time)
df.org3 <- zz



df <- df.org3

z_vars_use2 <- c("asthma_hosp_2yr","asthma_los_2yrgp", "asthma_pres_2yrgp1")
z_resp_vars2 <- c("result" ,"hosp_covid", "icu_death", "death_covid")
z <- df %>% dplyr::select_at(c(z_vars_use2, z_resp_vars2, "eave_weight"))
z <- z %>% pivot_longer(cols=all_of(z_vars_use2))
z.df <- z %>% group_by(name, value) %>% 
  dplyr::summarise(across(all_of(z_resp_vars2), ~ sum(.))) 
z1 <- z.df %>%  mutate(rate_hosp_covid = round(hosp_covid/result*1000,1),
                       rate_icu_death = round(icu_death/result*1000,1),
                       rate_death_covid = round(death_covid/result*1000,1)) %>% ungroup() %>% as.data.frame()

z1





z.rv <- "hosp_covid"
z.rv.time <- "time.test.to.hosp"

z.rv <- "icu_death"
z.rv.time <- "time.test.to.icudeath"

z.rv <- "death_covid"
z.rv.time <- "time.test.to.death"

z_exp.list <- c("asthma_hosp_2yr","zz$asthma_los_2yrgp", "asthma_pres_2yrgp1")

i=1
for(z_exp in z_exp.list){
  z.fmla <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  pspline(ageYear) + Sex + simd2020_sc_quintile + no_hosp_2yrgp + no_riskgp_gp + bmigp + vacc_status +",
                             z_exp))
  
  z.fit <- coxph(z.fmla , data=df, weights = eave_weight)

  z.r <- fun.cox(z.rv,z.fit)
  z.r$exp <- z_exp
  if(i==1) z.rr <- z.r else
    z.rr <- rbind(z.rr,z.r)
  i=i+1
  print(i)
}

write.csv(z.rr,paste(z.rv,"asthma.adults.sub.csv",sep=""))

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









