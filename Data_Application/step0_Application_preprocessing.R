rm(list = ls())


# --------------
#  library
# --------------
# library(ggplot2)
# library(patchwork) # for arranging multiple plots



# --------------
#   Path
# --------------
dir <- "path/Codes/"
dir_fun <- paste0(dir, "R_func/")
dir_save <- paste0("path/Results/")
dir1 <- "path/vaccine/p703/"
dir2 <- "path/vaccine/p704/"


# --------------
#   Data prep
# --------------
# survival
s703 <- read.csv(paste0(dir1, "v703_survival.csv"))
s704 <- read.csv(paste0(dir2, "v704_survival.csv"))

# covariates
covar703 <- read.csv(paste0(dir2, "v703_subject_master.csv"))
covar704 <- read.csv(paste0(dir2, "v704_subject_master.csv"))

# concentration
con703 <- read.csv(paste0(dir1, "nonmem.csv"))
con704 <- read.csv(paste0(dir2, "nonmem.csv"))

# eligible controls
dat_eliCtl <- read.csv(paste0(dir2, "amp_eligible_ctrls.csv"))

# mapping id
map703 <- read.csv(paste0(dir1, "id_pubid_ptid_703.csv"))
map704 <- read.csv(paste0(dir2, "id_pubid.csv"))

# risk_score
risk703 <- read.csv(paste0(dir1, "vaccine_ptids_with_riskscores.csv"))
risk704 <- read.csv(paste0(dir2, "vaccine_ptids_with_riskscores.csv"))

#----------------------------------------
# prepare HVTN703
# subset D51 concentration
subcon <- con703[con703$AVISITN == 501, ]
# initiate dat703 with IDs
id703 <- map703[match(subcon$ID, map703$ID),]
id703$ptidx <- as.numeric(gsub("-", "", id703$ptid))

# subset three other datasets
subcovar <- covar703[match(id703$pub_id, covar703$pub_id),]
subs <- s703[match(id703$ptid, s703$ptid),]
subrisk <- risk703[match(id703$ptidx, risk703$usubjid),]


# prepare geography info
# 1: South Africa -- 703
# 2: not South Africa -- 703
# 3: Brazil & Peru  --704
# 4: not Brazil & Peru (Switzerland & USA) --704
geo <- rep(2, length = nrow(id703))
geo[subcovar$country %in% c("South Africa")] <- 1

# prepare geography info
risk_score <- subrisk$standardized_risk_score


# final dataset including
# (1) survival outcome
# (2) D61 concentration
# (3) covariates: weights, risk score, geo, age group
dat703 <- data.frame(Trial = "703",
                     id = id703$ID,
                     pub_id = id703$pub_id,
                     ptid = id703$ptid,
                     time = subs$fudayswk104,
                     status = subs$statuswk104,
                     D61 = subcon$DV,
                     TRT = subcon$TRT,
                     weight = subcon$weight,
                     risk_score = risk_score,
                     geo = factor(geo),
                     agegrp = factor(subcovar$agegrp))




#----------------------------------------
# prepare HVTN704
# subset D51 concentration
subcon <- con704[con704$AVISITN == 501, ]
# initiate dat704 with IDs
id704 <- map704[match(subcon$ID, map704$ID),]
id704$ptid <- covar704$ptid[match(id704$PUB_ID, covar704$pub_id)]
id704$ptidx <- as.numeric(gsub("-", "", id704$ptid))

# subset three other datasets
subcovar <- covar704[match(id704$PUB_ID, covar704$pub_id),]
subs <- s704[match(id704$ptid, s704$ptid),]
subrisk <- risk704[match(id704$ptidx, risk704$usubjid),]


# prepare geography info
# 1: South Africa -- 703
# 2: not South Africa -- 703
# 3: Brazil & Peru  --704
# 4: not Brazil & Peru (Switzerland & USA) --704
geo <- rep(4, length = nrow(id704))
geo[subcovar$country %in% c("Brazil", "Peru")] <- 3

# prepare geography info
risk_score <- subrisk$standardized_risk_score


# final dataset including
# (1) Trial indicator
# (2) patient identification information: id, pub_id, ptid
# (3) survival outcome (time & status)
# (4) D61 concentration
# (5) treatment dosing group
# (6) covariates: weights, risk score, geo, age group
dat704 <- data.frame(Trial = "704",
                     id = id704$ID,
                     pub_id = id704$PUB_ID,
                     ptid = id704$ptid,
                     time = subs$fudayswk104,
                     status = subs$statuswk104,
                     D61 = subcon$DV,
                     TRT = subcon$TRT,
                     weight = subcon$weight,
                     risk_score = risk_score,
                     geo = factor(geo),
                     agegrp = factor(subcovar$agegrp))

# combine 703 and 704 data
dat <- rbind(dat703, dat704)


#----------------------------------------
# preparing sampling weights -- case-cohort sampling design
dat_eliCtl_noPDI <- filter(dat_eliCtl,is.na(pdistatus),delta2==1)

# eligible controls (numerator)
wt_N <- dat_eliCtl_noPDI %>% 
  group_by(Protocol,rx_code) %>% 
  summarise(n=n()) %>% 
  mutate(f=n,
         rx_code=ifelse(rx_code=="T1","10","30"),
         Protocol=ifelse(Protocol=="HVTN 703","703","704"),
         n=NULL) %>% 
  rename(study=Protocol,dose=rx_code,N=f)

# sampled set (denominator)
wt_n <- dat[dat$status==0,] %>% 
  group_by(Trial,TRT) %>% 
  summarise(n=n()) %>% 
  mutate(f=n,
         n=NULL) %>% 
  rename(study=Trial,dose=TRT,n=f)

# sampling weights for ctrl: stratified by Dose(10/30) and Trial(703/704)
wt_N$n <- wt_n$n
wt_N$wt <- wt_N$N / wt_n$n
# infection rate for sub-groups
# inf_rate <- c(26/642, 17/645, 31/895, 28/894)
inf_rate <- c(1,1,1,1)
names(inf_rate) <- c("703_d10", "703_d30", "704_d10", "704_d30")


#----------------------------------------
# add weights to dataset (sampling weights * infection rates)
dat <- dat %>% 
  group_by(Trial, TRT, status) %>% 
  mutate(wt = case_when(Trial == 703 & TRT == 10 & status == 1 ~ 1 / inf_rate[1],
                        Trial == 703 & TRT == 30 & status == 1 ~ 1 / inf_rate[2],
                        Trial == 704 & TRT == 10 & status == 1 ~ 1 / inf_rate[3],
                        Trial == 704 & TRT == 30 & status == 1 ~ 1 / inf_rate[4],
                        Trial == 703 & TRT == 10 & status == 0 ~ wt_N$wt[1] / inf_rate[1],
                        Trial == 703 & TRT == 30 & status == 0 ~ wt_N$wt[2] / inf_rate[2],
                        Trial == 704 & TRT == 10 & status == 0 ~ wt_N$wt[3] / inf_rate[3],
                        Trial == 704 & TRT == 30 & status == 0 ~ wt_N$wt[4] / inf_rate[4]))

saveRDS(dat, file = paste0(dir_save, "dat_pooled.RDS"))


#----------------------------------------
# one-hot encoding
#----------------------------------------
# geography
dummy_geo <- sapply(sort(unique(dat$geo)), function(x){ifelse(dat$geo == x, 1, 0)})
colnames(dummy_geo) <- paste0("geo_", 1:4)

# age
# 1: "18 - 20" 
# 2: "21 - 30" 
# 3: "31 - 40" 
# 4: "41 - 50"
dummy_agegrp <- sapply(sort(unique(dat$agegrp)), function(x){ifelse(dat$agegrp == x, 1, 0)})
colnames(dummy_agegrp) <- paste0("agegrp_", 1:4)


dat <- cbind(dat, dummy_geo, dummy_agegrp)
saveRDS(dat, file = paste0(dir_save, "dat_pooled_onehot.RDS"))





#-----------------
# PS fit
#-----------------
set.seed(241121)
fit_aw <- EstPropScore(a = dat$D61,
                       w = Ws_g,
                       no.folds = 10,
                       no.eval = 20)
# generate prediction of E(A=a | W)
gn_As <- fit_aw$cond.dens
saveRDS(gn_As, file = paste0(dir_save, "res/gn_As.rds"))
