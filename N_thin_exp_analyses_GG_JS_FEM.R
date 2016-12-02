################################################################################################
# R code to reproduce results in 
# Granath and Strengbom - Nitrogen fertilization reduces wild berry production in boreal forests
# Submitted to Journal of Forest ecology and Management
# Contact: Gustaf.Granath@gmail.com
#################################################################################################

# Load required packages
library(MCMCglmm)
library(ggplot2)
library(dplyr)

# load data
dat <- read.csv("N_thin_exp_Granath_Strengbom_FEM.csv")

# Calulate sums for subplots within plots. I.e., the total for 3.5 m2
sums <- aggregate(cbind(bilberry_fruits, bilberry_fungi, bilberry_herbivory, bilberry_cover,
                        cowberry_fruits, cowberry_fungi, cowberry_herbivory, cowberry_cover) ~ 
                    exp_name + treatment_name + latitude + year + thin+N+P, dat, sum)
# centralize latitude
sums$latitude.s <- scale(sums$latitude, scale=FALSE)

# Year is factor
sums$year <- factor(sums$year)

# General stats on subplot level (0.25 m2)
# mean, min, max
dd <- dat %>% 
  group_by(year) %>% 
  summarise_at(vars(bilberry_cover, bilberry_fungi, bilberry_herbivory,
                    cowberry_cover, cowberry_fungi, cowberry_herbivory),
               funs(mean, min, max)) %>%
  as.data.frame()
dd
# At plot level and per m2 
per.sqm <- sums
per.sqm$bilberry_fruits.ave <- sums$bilberry_fruits/3.5
per.sqm$cowberry_fruits.ave <- sums$cowberry_fruits/3.5
per.sqm$bilberry_cover.ave <- sums$bilberry_cover/14
per.sqm$cowberry_cover.ave <- sums$cowberry_cover/14
per.sqm$bilberry_fungi.ave <- sums$bilberry_fungi/14
per.sqm$cowberry_fungi.ave <- sums$cowberry_fungi/3.5
per.sqm$bilberry_herbivory.ave <- sums$bilberry_herbivory/14
per.sqm$cowberry_herbivory.ave <- sums$cowberry_herbivory/14

dd <- per.sqm %>% 
  group_by(year) %>% 
  summarise_at(vars(bilberry_fruits.ave, bilberry_cover.ave, bilberry_fungi.ave, bilberry_herbivory.ave,
                    cowberry_fruits.ave, cowberry_cover.ave, cowberry_fungi.ave, cowberry_herbivory.ave),
               funs(mean, min, max)) %>%
  as.data.frame()
dd

#Overview of treatments and sample size
table(sums[,c("thin","N", "P")])

##### Make map over experimental sites ####
# Figure 1

library(ggmap)

# creating a sample data.frame with lat/lon points
df <- as.data.frame(cbind(dat$longitude,dat$latitude))
df <- unique(df)
# getting the map
map.gg <- get_map(location = c(lon = mean(df$V1), lat = mean(df$V2)), zoom = 4,
                  maptype = "satellite", scale = 2)

# plotting the map with some points on it
ggg.map <- ggmap(map.gg) +
  geom_point(data = df, aes(x = V1, y = V2, fill = "red", alpha = 0.8), size = 3, shape = 21) +
  guides(fill=FALSE, alpha=FALSE, size=FALSE)
ggsave("gggMap.png")


#### Statistical analyses ####

##### Plant cover
# Gaussian models give same results as binomial and for illustrative purposes we
# present gaussian models where the effect is on the %-scale

# Make id variable to include repeated measurement random structure
sums$id <- factor(paste(sums$treatment_name,sums$exp_name, sep="_"))

prior = list(R=list(V=1, nu=0.002), 
             G=list(G1=list(V=diag(2), nu=0.002, alpha.mu=c(0,0), alpha.V=diag(2)*625^2), 
                    G2=list(V=1, nu=1,alpha.mu=0, alpha.V=625^2))) 
# average cover per plot
sums$bilberry_cover.ave <- sums$bilberry_cover/14

# Model with year*treatments did not indicate an interaction effect so 
# a simpler model is run.
bil.cov.mc <- MCMCglmm(bilberry_cover.ave ~ year + N*thin + P + latitude.s, 
                        random = ~idh(year):exp_name + id, data=sums, 
                        family = "gaussian",  nitt = 80000, burnin = 15000, 
                        thin=25, pr = TRUE, pl = TRUE, saveX = TRUE,  saveZ = TRUE, 
                        prior=prior, singular.ok=TRUE )
summary(bil.cov.mc)
# No treatment effects and very small effect sizes (<5%)

# Plot bilberry effect
# 1) the model matrix, 2) make model predictions and
# 3) fix contrasts to be plotted
raw.means <- aggregate(bilberry_cover ~ year+thin*N+P, sums, mean)
newDat <- raw.means[,1:4] 
newDat$latitude.s <- 0
X <- model.matrix(formula(bil.cov.mc$Fixed$formula,fixed.only=TRUE)[-2],
                  newDat)
res <- apply(bil.cov.mc$Sol[,3:7],1, function (x) x %*% t(X[,-c(1:2)]))
raw.means$eff <- rowMeans(res)
cis <- t(apply(res, 1, function (x) HPDinterval(as.mcmc(x))))
raw.means <- cbind(raw.means, cis)
colnames(raw.means)[7:8] <- c("lo_95", "up_95")
lat.eff <- summary(bil.cov.mc)$solutions[6,1:3] #latitude effect
names(lat.eff) <- colnames(raw.means)[6:8]
raw.means <- rbind(raw.means, c(raw.means[10,1:5], lat.eff) )
raw.means <- raw.means[c(1,3,5,7,9,11),]
raw.means$vars <- c("control", "Thin","N", "Thin+N", "Thin+N+P", "Latitude")

bil.cov.p <- ggplot(raw.means[-1,], aes(x=vars, y=eff)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=lo_95, ymax=up_95), 
                lwd=0.7, width=0) +
  geom_point(size=4, pch=21, fill="black") +
  ylim(c(-30, 30)) +
  xlab("") +
  ylab("Absolute change (%)") +
  theme(axis.text.x  = element_text(size=14),
        axis.text.y  = element_text(size=14),
        axis.title = element_text(size=14)) +
  annotate("text", label = "a)", y = -23, x = 5.3, size = 10) +
  coord_flip()
# done bilberry

# Cowberry

# MCMC version for normal errors and transformation.
# Results close to the more "accurate" binomial model which is more
# complicated to interpret and comapre, with the bilberry model.
prior=list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=625^2),
                                         G2=list(V=1, nu=1, alpha.mu=0, alpha.V=625^2)))
sums$cowberry_cover.ave <- sums$cowberry_cover/14
cow.cov.mc <- MCMCglmm(sqrt(cowberry_cover.ave) ~ year+N*thin + P + latitude.s, random = ~exp_name + id, 
                       data=sums, family = "gaussian",  nitt = 110000, burnin = 15000,
                       thin=50, pr = TRUE, pl = TRUE, saveX = TRUE,  saveZ = TRUE, prior=prior)
summary(cow.cov.mc)

#plot effect
raw.means <- aggregate(cowberry_cover.ave ~ year+thin*N+P, sums, mean)
newDat <- raw.means[,1:4] 
newDat$latitude.s <- 0
X <- model.matrix(formula(cow.cov.mc$Fixed$formula,fixed.only=TRUE)[-2],
                  newDat)
res <- apply(cow.cov.mc$Sol[,1:7],1, function (x) x %*% t(X))

# Sqrt tranformed response so for comparison with bilberry model, which is untransformed,
# effect sizes are recalculated on the response scale.
# Effects are calculated by taking the difference beyween treatment groups
# while contrasts are kept the same.
res_eff1 <- (apply(res, 2, function (x) x[c(3,5,7,9)]^2 - x[1]^2 ))
res <- res_eff1[c(1,2,3,4),] 
raw.means <- raw.means[c(3,5,7,9),]
raw.means$eff <- rowMeans(res)
cis <- t(apply(res, 1, function (x) HPDinterval(as.mcmc(x))))
raw.means <- cbind(raw.means, cis)
colnames(raw.means)[7:8] <- c("lo_95", "up_95")
lat.eff <- summary(cow.cov.mc)$solutions[6,1:3] #latitude effect
names(lat.eff) <- colnames(raw.means)[6:8]
raw.means <- rbind(raw.means, c(raw.means[10,1:5], lat.eff) )
raw.means$vars <- c("Thin", "N", "Thin+N", "Thin+N+P", "Latitude")

cow.cov.p <- ggplot(raw.means, aes(x=vars, y=eff)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=lo_95, ymax=up_95), 
                lwd=0.7, width=0) +
  geom_point(size=4, pch=21, fill="black") +
  ylim(c(-30, 30)) +
  xlab("") +
  ylab("Absolute change (%)") +
  theme(axis.text.x  = element_text(size=14),
        axis.text.y  = element_text(size=14),
        axis.title = element_text(size=14)) +
  annotate("text", label = "b)", y = -23, x = 5.3, size = 10) +
  coord_flip()

# Merge plots to produce Figure 2
library(gridExtra)
png("figure2_cover_plot.png", width=15, height=30, units="cm", res=300)
grid.arrange(bil.cov.p , cow.cov.p , ncol=1, nrow =2)
dev.off()

#### Fungal incidence ####

# Make id variable to include repeated measurement random structure
sums$id <- factor(paste(sums$treatment_name,sums$exp_name, sep="_"))
sums$bilberry_cover.ave <- sums$bilberry_cover/14 # calc mean cover

# Priors with offset to account for availabillity of plants
B= list(mu = matrix(c(0,0,0,0,0,0,1,0,0,0,0,0),12),V = diag(12)*(10))
prior = list(B=B, R=list(V=1, nu=0.002), 
             G=list(G1=list(V=diag(2), nu=0.002, alpha.mu=c(0,0), alpha.V=diag(2)*1000), 
                    G2=list(V=1, nu=0.002,alpha.mu=0, alpha.V=625^2))) 
diag(prior$B$V)[7]<-1e-9

sums$quad <- 1400 # number of total quadrats
bil.fu.mc <- MCMCglmm(cbind(bilberry_fungi, quad-bilberry_fungi)  ~ 
                        year*N*thin + year*P + latitude.s + log(bilberry_cover.ave), 
                       random = ~idh(year):exp_name + id, data=sums, family = "multinomial2",  
                       nitt = 80000, burnin = 15000, thin=25, 
                       pr = TRUE, pl = TRUE, saveX = TRUE,  saveZ = TRUE, prior=prior)

summary(bil.fu.mc) 

# plot effect
# see general explanation for cover
raw.means <- aggregate(bilberry_fungi ~ year*thin*N+year*P, sums, mean)
newDat <- raw.means[,1:4] 
newDat$latitude.s <- 0
newDat$bilberry_cover.ave <- 1
X <- model.matrix(formula(bil.fu.mc$Fixed$formula,fixed.only=TRUE)[-2],
                  newDat)
res <- apply(bil.fu.mc$Sol[,3:12],1, function (x) x %*% t(X[,-c(1:2)]))
raw.means$eff <- rowMeans(res)
cis <- t(apply(res, 1, function (x) HPDinterval(as.mcmc(x))))
raw.means <- cbind(raw.means, cis)
colnames(raw.means)[7:8] <- c("lo_95", "up_95")
lat.eff <- summary(bil.fu.mc)$solutions[6,1:3] #latitude effect
names(lat.eff) <- colnames(raw.means)[6:8]
raw.means <- rbind(raw.means, c(raw.means[10,1:5], lat.eff) )
raw.means$vars <- c("control", "control", "Thin-2014", "Thin-2015", "N-2014", "N-2015",
                    "Thin+N-2014", "Thin+N-2015", "Thin+N+P-2014", "Thin+N+P-2015", "latitude")
raw.means[,6:8] <- (exp(raw.means[,6:8]))

bil.fun.p <- ggplot(raw.means[-c(1:2),], aes(x=vars, y=eff)) + 
  geom_hline(yintercept=1, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=lo_95, ymax=up_95), 
                lwd=0.7, width=0) +
  geom_point(size=4, pch=21, fill="black") +
  ylim(c(0, 6)) +
  xlab("") +
  ylab("Odds ratio") +  # (treatment:control)
  theme(axis.text.x  = element_text(size=14),
        axis.text.y  = element_text(size=14),
        axis.title = element_text(size=14)) +
  annotate("text", label = "a)", y = 0, x = 9, size = 10) +
  coord_flip()

# cowberry
# Response is number of infected shoots. 
# Hence Poisson errors are used. 
sums$cowberry_cover.ave <- sums$cowberry_cover/14 # calc mean cover

# Priors with offset term (+0.5 in the model because there are zeros in the data)
B= list(mu = matrix(c(0,0,0,0,0,0,1,0,0,0,0,0),12),V = diag(12)*(10))
prior = list(B=B, R=list(V=1, nu=0.002), 
             G=list(G1=list(V=1, nu=0.002, alpha.mu=0, alpha.V=625^2), 
                    G2=list(V=1, nu=0.002,alpha.mu=0, alpha.V=625^2))) 
diag(prior$B$V)[7]<-1e-9

cow.fu.mc <- MCMCglmm(cowberry_fungi  ~ year*N*thin + year*P + latitude.s + log(cowberry_cover.ave+0.5), 
                      random = ~exp_name + id, data=sums, family = "poisson",  nitt = 80000, burnin = 15000,
                      thin=25, pr = TRUE, pl = TRUE, saveX = TRUE,  saveZ = TRUE, prior=prior)

summary(cow.fu.mc)

# plot effect
raw.means <- aggregate(cowberry_fungi ~ year*thin*N+year*P, sums, mean)
newDat <- raw.means[,1:4] 
newDat$latitude.s <- 0
newDat$cowberry_cover.ave <- 0.5
X <- model.matrix(formula(cow.fu.mc$Fixed$formula,fixed.only=TRUE)[-2],
                  newDat)
res <- apply(cow.fu.mc$Sol[,3:12],1, function (x) x %*% t(X[,-c(1:2)]))
raw.means$eff <- rowMeans(res)
cis <- t(apply(res, 1, function (x) HPDinterval(as.mcmc(x))))
raw.means <- cbind(raw.means, cis)
colnames(raw.means)[7:8] <- c("lo_95", "up_95")
lat.eff <- summary(cow.fu.mc)$solutions[6,1:3] #latitude effect
names(lat.eff) <- colnames(raw.means)[6:8]
raw.means <- rbind(raw.means, c(raw.means[10,1:5], lat.eff) )
raw.means$vars <- c("control", "control", "Thin-2014", "Thin-2015", "N-2014", "N-2015",
                    "Thin+N-2014", "Thin+N-2015", "Thin+N+P-2014", "Thin+N+P-2015", "latitude")
raw.means[,6:8] <- (exp(raw.means[,6:8]))

cow.fun.p <- ggplot(raw.means[-c(1:2),], aes(x=vars, y=eff)) + 
  geom_hline(yintercept=1, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=lo_95, ymax=up_95), 
                lwd=0.7, width=0) +
  geom_point(size=4, pch=21, fill="black") +
  ylim(c(0, 2)) +
  xlab("") +
  ylab("Rate ratio") +
  theme(axis.text.x  = element_text(size=14),
        axis.text.y  = element_text(size=14),
        axis.title = element_text(size=14)) +
  annotate("text", label = "b)", y = 0, x = 9, size = 10) +
  coord_flip()

# Merge plots to get figure 3
library(gridExtra)
png("figure3_fungi.png", width=30, height=15, units="cm", res=300)
grid.arrange(bil.fun.p , cow.fun.p , ncol=2, nrow =1)
dev.off()

#### Herbivory ####

# bilberry
# Fix data (see above)
sums$id <- factor(paste(sums$real_treat,sums$exp_name, sep="_"))
sums$billberry_cover.ave <- sums$bilberry_cover/14
sums$quad <- 1400

# Priors with an offset term
B= list(mu = matrix(c(0,0,0,0,0,0,1,0,0,0,0,0),12),V = diag(12)*(10))
prior = list(B=B, R=list(V=1, nu=0.002), 
             G=list(G1=list(V=diag(2), nu=0.002, alpha.mu=c(0,0), alpha.V=diag(2)*1000), 
                    G2=list(V=1, nu=0.002,alpha.mu=0, alpha.V=625^2))) 
diag(prior$B$V)[7]<-1e-9

bil.he.mc <- MCMCglmm(cbind(bilberry_herbivory, quad-bilberry_herbivory)  ~ 
                        year*N*thin + year*P + latitude.s + log(bilberry_cover.ave), 
                       random = ~idh(year):exp_name + id, data=sums, family = "multinomial2",  
                       nitt = 80000, burnin = 15000, thin=25, 
                       pr = TRUE, pl = TRUE, saveX = TRUE,  saveZ = TRUE, prior=prior)
summary(bil.he.mc) 
# A neg N effect on herbivory the first year, nothing else.
# Note! if site "929Nysund" is removed,
# there is an effect of latitude

# cowberry

# Fix data
sums$cowberry_cover.ave <- sums$cowberry_cover/14 # calc mean cover

# Priors with an offset term
B= list(mu = matrix(c(0,0,0,0,0,0,1,0,0,0,0,0),12),V = diag(12)*(10))
prior = list(B=B, R=list(V=1, nu=0.002), 
             G=list(G1=list(V=1, nu=0.002, alpha.mu=0, alpha.V=625^2), 
                    G2=list(V=1, nu=0.002,alpha.mu=0, alpha.V=625^2)))
diag(prior$B$V)[7]<-1e-9

cow.he.mc <- MCMCglmm(cbind(cowberry_herbivory, quad-cowberry_herbivory)  ~ 
                        year*N*thin + year*P + latitude.s + log(cowberry_cover.ave+0.5), 
                      random = ~exp_name + id, data=sums, family = "multinomial2",  
                      nitt = 80000, burnin = 15000, thin=25, 
                      pr = TRUE, pl = TRUE, saveX = TRUE,  saveZ = TRUE, prior=prior)
summary(cow.he.mc) # NO significant effects


#### Fruit production ####

# Make id variable to include repeated measurement random structure
sums$id <- factor(paste(sums$treatment_name, sums$exp_name, sep="_"))

# First model not controlling for precent cover (ie. no offset term)
# Prior 
prior = list( R=list(V=1, nu=0.002), 
              G=list(G1=list(V=diag(2), nu=0.002, alpha.mu=c(0,0), alpha.V=diag(2)*625^2), 
              G2=list(V=1, nu=0.002,alpha.mu=0, alpha.V=625^2))) 

#sums$thin <- relevel(sums$thin, ref="no_thin") # thinning as reference level instead
bil.mc <- MCMCglmm(bilberry_fruits ~ year*N*thin + year*P + latitude.s, 
                    random = ~idh(year):exp_name + id, 
                    data=sums, family = "poisson",  nitt = 80000, 
                    burnin = 15000, thin=25, pr = TRUE, pl = TRUE,
                    saveX = TRUE, saveZ = TRUE,  singular.ok=TRUE )
summary(bil.mc)
# random slope model with covariance between intercept and slope only 2 DIC smaller, 
# so a simpler model without covariance is OK.

# Estimate an overall year effect
# sums$year <- relevel(sums$year, ref="2015") # thinning as reference level instead
bil.mc.year <- MCMCglmm(bilberry_fruits ~ year + latitude.s, 
                        random = ~idh(year):exp_name + id, 
                         data=sums, family = "poisson",  nitt = 80000, 
                        burnin = 15000, thin=25, pr = TRUE, pl = TRUE,
                         saveX = TRUE, saveZ = TRUE,  singular.ok=TRUE )
summary(bill.mc.year)

# In the next model we controll for percent plant cover, i.e. add an offset for plant cover. 
# With the offset term we get fruits per percent cover.

# Prior 
B= list(mu = matrix(c(0,0,0,0,0,0,1,0,0,0,0,0),12),V = diag(12)*(10))
prior = list(B=B, R=list(V=1, nu=0.002), 
             G=list(G1=list(V=diag(2), nu=0.002, alpha.mu=c(0,0), alpha.V=diag(2)*1000), 
                    G2=list(V=1, nu=0.002,alpha.mu=0, alpha.V=625^2))) 
diag(prior$B$V)[7]<-1e-9

#sums2 <- sums[!(sums$exp_name=="929Nysund"),] #if removed an effect of latitude
bil.con.mc <- MCMCglmm(bilberry_fruits ~ year*N*thin + year*P + latitude.s + log(bilberry_cover), 
                       random = ~idh(year):exp_name + id, data=sums, family = "poisson",  
                       nitt = 80000, burnin = 15000, thin=25, 
                       pr = TRUE, pl = TRUE, saveX = TRUE,  saveZ = TRUE, 
                       prior=prior, singular.ok=TRUE )

summary(bil.con.mc)

# Plot effect
# first without offset term (ie not cotrolling for plant cover)
# See earlier explanation of these effect plots
raw.means <- aggregate(bilberry_fruits ~ year*thin*N+year*P, sums, mean)
newDat <- raw.means[,1:4] 
newDat$latitude.s <- 0
X <- model.matrix(formula(bil.mc$Fixed$formula,fixed.only=TRUE)[-2],
                  newDat)
res <- apply(bil.mc$Sol[,3:11],1, function (x) x %*% t(X[,-c(1:2)]))
raw.means$eff <- rowMeans(res)
cis <- t(apply(res, 1, function (x) HPDinterval(as.mcmc(x))))
raw.means <- cbind(raw.means, cis)
colnames(raw.means)[7:8] <- c("lo_95", "up_95")
lat.eff <- summary(bil.mc)$solutions[6,1:3] #latitude effect
names(lat.eff) <- colnames(raw.means)[6:8]
raw.means <- rbind(raw.means, c(raw.means[10,1:5], lat.eff) )
raw.means$vars <- c("control", "control", "Thin-2014", "Thin-2015", "N-2014", "N-2015",
                    "Thin+N-2014", "Thin+N-2015", "Thin+N+P-2014", "Thin+N+P-2015", "latitude")
raw.means[,6:8] <- (exp(raw.means[,6:8]))

# then with offset
raw.means.off <- aggregate(bilberry_fruits ~ year*thin*N+year*P, sums, mean)
newDat.off <- raw.means.off[,1:4] 
newDat.off$latitude.s <- 0
newDat.off$bilberry_cover <- 1 # per m2 if offset for area is in the model
X.off <- model.matrix(formula(bil.con.mc$Fixed$formula,fixed.only=TRUE)[-2],
                      newDat.off)
res.off <- apply(bil.con.mc$Sol[,3:12],1, function (x) x %*% t(X.off[,-c(1:2)]))
raw.means.off$eff <- rowMeans(res.off)
cis.off <- t(apply(res.off, 1, function (x) HPDinterval(as.mcmc(x))))
raw.means.off <- cbind(raw.means.off, cis.off)
colnames(raw.means.off)[7:8] <- c("lo_95", "up_95")
lat.eff.off <- summary(bil.con.mc)$solutions[6,1:3] #latitude effect
names(lat.eff.off) <- colnames(raw.means.off)[6:8]
raw.means.off <- rbind(raw.means.off, c(raw.means.off[10,1:5], lat.eff.off) )
raw.means.off$vars <- c("control", "control", "Thin-2014", "Thin-2015", "N-2014", "N-2015",
                        "Thin+N-2014", "Thin+N-2015", "Thin+N+P-2014", "Thin+N+P-2015", "latitude")
raw.means.off[,6:8] <- (exp(raw.means.off[,6:8]))

# Make plot
bil.p <- ggplot(raw.means[-c(1:2),], aes(x=vars, y=eff)) + 
  geom_hline(yintercept=1, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=lo_95, ymax=up_95), 
                lwd=0.7, colour="red", width=0, , position=position_nudge(x=0.1)) +
  geom_point(size=2, pch=21, aes(fill="yellow"), position=position_nudge(x=0.1)) +
  geom_point(data=raw.means.off[-c(1:2),], aes(x=vars, y=eff, fill="blue"), position=position_nudge(x=-0.1), 
             size=2, pch=21) +
  geom_errorbar(data=raw.means.off[-c(1:2),], aes(ymin=lo_95, ymax=up_95), 
                lwd=0.7, colour="black", width=0, , position=position_nudge(x=-0.1)) +
  ylim(c(0, 4)) +
  xlab("") +
  ylab("Rate ratio") +
  theme(axis.text.x  = element_text(size=14),
        axis.text.y  = element_text(size=14)) +
  coord_flip() +
  scale_fill_manual(name="Model",values=c("blue", "yellow"), breaks = c("yellow", "blue"), 
                    labels= c("Overall effect", "Controlling for plant cover")) +
  theme(legend.justification=c(1,1), legend.position=c(1,1))
ggsave("figure4_fruit_cbilberry.png", bil.p) 

# Cowberry

# Priors, no offset
prior = list(R=list(V=1, nu=0.002), 
             G=list(G1=list(V=diag(1), nu=0.002, alpha.mu=c(0), alpha.V=diag(1)*625^2), 
                    G2=list(V=1, nu=0.002,alpha.mu=0, alpha.V=625^2))) 

cow.mc <- MCMCglmm(cowberry_fruits ~ year*N*thin + year*P + latitude.s, 
                   random = ~ exp_name + id, data=sums, singular.ok=TRUE, 
                   family = "poisson",  nitt = 110000,
                   burnin = 15000, thin=50, pr = TRUE, 
                   pl = TRUE, saveX = TRUE,  saveZ = TRUE, prior=prior)
summary(cow.mc)
# DIC 6 higher without year-specific variances (~idh(year):exp_name ), 
# so not required

# Overall year effect
cow.mc.year <- MCMCglmm(cowberry_fruits ~ year + latitude.s, 
                        random = ~ exp_name + id, data=sums, singular.ok=TRUE, 
                        family = "poisson",  nitt = 110000,
                        burnin = 15000, thin=50, pr = TRUE, 
                        pl = TRUE, saveX = TRUE,  saveZ = TRUE, prior=prior)
summary(cow.mc.year)


# Model with offset, estimating the effect per percent plant cover
# Priors
B= list(mu = matrix(c(0,0,0,0,0,0,1,0,0,0,0,0),12),V = diag(12)*(10))
prior = list(B=B,R=list(V=1, nu=0.002), 
             G=list(G1=list(V=1, nu=0.002, alpha.mu=c(0), alpha.V=diag(1)*625^2), 
                    G2=list(V=1, nu=0.002,alpha.mu=0, alpha.V=625^2))) 
diag(prior$B$V)[7]<-1e-9

cow.con.mc <- MCMCglmm(cowberry_fruits ~ year*N*thin + year*P + latitude.s + log(cowberry_cover+0.5), 
                       random = ~exp_name + id, data=sums, singular.ok=TRUE, family = "poisson",  
                       nitt = 110000, burnin = 15000, thin=50, pr = TRUE, 
                       pl = TRUE, saveX = TRUE,  saveZ = TRUE, prior=prior)

summary(cow.con.mc)

# Plot effect
raw.means <- aggregate(cowberry_fruits ~ year*thin*N+year*P, sums, mean)
newDat <- raw.means[,1:4] 
newDat$latitude.s <- 0
X <- model.matrix(formula(cow.mc$Fixed$formula,fixed.only=TRUE)[-2],
                  newDat)
res <- apply(cow.mc$Sol[,3:11],1, function (x) x %*% t(X[,-c(1:2)]))
raw.means$eff <- rowMeans(res)
cis <- t(apply(res, 1, function (x) HPDinterval(as.mcmc(x))))
raw.means <- cbind(raw.means, cis)
colnames(raw.means)[7:8] <- c("lo_95", "up_95")
lat.eff <- summary(cow.mc)$solutions[6,1:3] #latitude effect
names(lat.eff) <- colnames(raw.means)[6:8]
raw.means <- rbind(raw.means, c(raw.means[10,1:5], lat.eff) )
raw.means$vars <- c("control", "control", "Thin-2014", "Thin-2015", "N-2014", "N-2015",
                    "Thin+N-2014", "Thin+N-2015", "Thin+N+P-2014", "Thin+N+P-2015", "latitude")

# then with offset term
raw.means.off <- aggregate(cowberry_fruits ~ year*thin*N+year*P, sums, mean)
newDat.off <- raw.means.off[,1:4] 
newDat.off$latitude.s <- 0
newDat.off$cowberry_cover <- 0.5 # per m2 if offset for area is in the model
X.off <- model.matrix(formula(cow.con.mc$Fixed$formula,fixed.only=TRUE)[-2],
                      newDat.off)
res.off <- apply(cow.con.mc$Sol[,3:12],1, function (x) x %*% t(X.off[,-c(1:2)]))
raw.means.off$eff <- rowMeans(res.off)
cis.off <- t(apply(res.off, 1, function (x) HPDinterval(as.mcmc(x))))
raw.means.off <- cbind(raw.means.off, cis.off)
colnames(raw.means.off)[7:8] <- c("lo_95", "up_95")
lat.eff.off <- summary(cow.con.mc)$solutions[6,1:3] #latitude effect
names(lat.eff.off) <- colnames(raw.means.off)[6:8]
raw.means.off <- rbind(raw.means.off, c(raw.means.off[10,1:5], lat.eff.off) )
raw.means.off$vars <- c("control", "control", "Thin-2014", "Thin-2015", "N-2014", "N-2015",
                        "Thin+N-2014", "Thin+N-2015", "Thin+N+P-2014", "Thin+N+P-2015", "latitude")

cow.p <- ggplot(raw.means[-c(1:2),], aes(x=vars, y=eff)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=lo_95, ymax=up_95), 
                lwd=0.7, colour="red", width=0, , position=position_nudge(x=0.1)) +
  geom_point(size=2, pch=21, aes(fill="yellow"), position=position_nudge(x=0.1)) +
  geom_point(data=raw.means.off[-c(1:2),], aes(x=vars, y=eff, fill="blue"), position=position_nudge(x=-0.1), 
             size=2, pch=21) +
  geom_errorbar(data=raw.means.off[-c(1:2),], aes(ymin=lo_95, ymax=up_95), 
                lwd=0.7, colour="black", width=0, , position=position_nudge(x=-0.1)) +
  ylim(c(-8, 7)) +
  xlab("") +
  ylab("Log rate ratio") +
  theme(axis.text.x  = element_text(size=14),
        axis.text.y  = element_text(size=14)) +
  coord_flip() +
  scale_fill_manual(name="Model",values=c("blue", "yellow"), breaks = c("yellow", "blue"), 
                    labels= c("Overall effect", "Controlling for plant cover")) +
  theme(legend.justification=c(0,1), legend.position=c(0,1))
ggsave("figure5_fruit_cowberry.png", cow.p) 


#### SEM model ####

# SEM mcmc way

# Fix variables
sums$id <- factor(paste(sums$treatment_name, sums$exp_name, sep="_"))
sums$quad <- 1400

# bilberry
sums$bilberry_cover.ave <- sums$bilberry_cover/14
sums$bilberry_fungi.ave <- sums$bilberry_fungi/14

# Priors
prior = list(R=list(V=1, nu=0.002), 
             G=list(G1=list(V=diag(2), nu=1.002, alpha.mu=c(0,0), alpha.V=diag(2)*625^2), 
                    G2=list(V=1, nu=0.002,alpha.mu=0, alpha.V=625^2)))

bil.cov.mc.sem <- MCMCglmm(bilberry_cover.ave  ~ N + thin + latitude.s + year, 
                           random = ~idh(year):exp_name + id, data=sums, 
                           family = "gaussian",  nitt = 80000, burnin = 15000, thin=25, 
                           pr = TRUE, pl = TRUE, saveX = TRUE,  saveZ = TRUE, prior=prior)

bil.fun.mc.sem <- MCMCglmm(cbind(bilberry_fungi, quad - bilberry_fungi) ~ 
                             bilberry_cover.ave + latitude.s + N + thin + year, 
                           random = ~idh(year):exp_name + id, data=sums, 
                           family = "multinomial2",  nitt = 80000, burnin = 15000, 
                           thin=25, pr = TRUE, pl = TRUE, saveX = TRUE,  saveZ = TRUE, prior=prior)

bil.fruit.mc.sem <- MCMCglmm(bilberry_fruits  ~ bilberry_fungi.ave + bilberry_cover.ave + 
                               latitude.s + N + thin + year, 
                             random = ~idh(year):exp_name + id, data=sums, family = "poisson",  
                             nitt = 80000, burnin = 15000, thin=25, pr = TRUE, pl = TRUE, 
                             saveX = TRUE,  saveZ = TRUE, prior=prior)

summary(bil.cov.mc.sem)
summary(bil.fun.mc.sem)
summary(bil.fruit.mc.sem)

# Cowberry
sums$cowberry_cover.ave <- sums$cowberry_cover/14

# priors
prior = list(R=list(V=1, nu=0.002), 
             G=list(G1=list(V=diag(1), nu=0.002, alpha.mu=c(0), alpha.V=diag(1)*625^2), 
                    G2=list(V=1, nu=0.002,alpha.mu=0, alpha.V=625^2)))

# Random slope not needed
cow.cov.mc.sem <- MCMCglmm(cbind(cowberry_cover, quad-cowberry_cover)  ~ 
                             N + thin + latitude.s + year, random = ~exp_name + id, 
                           data=sums, family = "multinomial2",  nitt = 80000, 
                           burnin = 15000, thin=25, pr = TRUE, pl = TRUE, 
                           saveX = TRUE,  saveZ = TRUE, prior=prior)

cow.fun.mc.sem <- MCMCglmm(cowberry_fungi  ~ cowberry_cover.ave + latitude.s + N + thin + year, 
                           random = ~exp_name + id, data=sums, family = "poisson",  
                           nitt = 80000, burnin = 15000, thin=25, pr = TRUE, 
                           pl = TRUE, saveX = TRUE,  saveZ = TRUE, prior=prior)

cow.fruit.mc.sem <- MCMCglmm(cowberry_fruits  ~ cowberry_fungi + cowberry_cover.ave + 
                               latitude.s + N + thin + year, 
                             random = ~exp_name + id, data=sums, family = "poisson",  
                             nitt = 80000, burnin = 15000, thin=25, pr = TRUE, 
                             pl = TRUE, saveX = TRUE,  saveZ = TRUE, prior=prior)

summary(cow.cov.mc.sem)
summary(cow.fun.mc.sem)
summary(cow.fruit.mc.sem)

