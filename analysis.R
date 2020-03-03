
##--------------------------------------------------------------------------------
## Authors: FJ Heather, DZ Childs
## Date: 22-05-18
## Title: Seabream growth IPM
## Description: For the paper - 

## Heather FJ, Childs DZ, Darnaude AM, Blanchard JL (2018) 
## Using an integral projection model to assess the effect of temperature on the
## growth of gilthead seabream Sparus aurata. PLoS ONE 13(5): e0196092. 

## https://doi.org/10.1371/journal.pone.0196092
##--------------------------------------------------------------------------------

# required packages
list.of.packages <- c("tidyverse", "nlme", "data.table")
new.packages     <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# loading required packages
lapply(list.of.packages, require, character.only = TRUE)

##--------------------------------------------------------------------------------
## Importing data
##--------------------------------------------------------------------------------

useData  <- 
  fread("https://datadryad.org/stash/downloads/file_stream/43742") %>%
  data.frame() %>% 
  select(-1) %>%
  mutate(covariate = summerT) %>% 
  filter(annuli_num > 0, OL.t1<2.7) %>%
  na.omit # this also removes fish_id == M101 (anomylous)

al_dat <- 
  fread("https://datadryad.org/stash/downloads/file_stream/43741") %>%
  data.frame() %>%
  select(-1) %>% 
  filter(rad4>1) %>%
  mutate(match(month_capture,month.name) + (age*12))

gr_dat <- 
  useData %>%
  select(fish.id, birth, month.capture, 
         year.capture, annuli_num, OL.t1, OL.t, age, covariate) %>%
  rename(fish_id = fish.id, birth_year = birth, 
         month_capture = month.capture, year_capture = year.capture)

lagoonData <- 
  fread("https://datadryad.org/stash/downloads/file_stream/43744") %>%
  data.frame()

summerT <- 
  fread("https://datadryad.org/stash/downloads/file_stream/43743")  %>% 
  data.frame() %>%
  select(-c(1,season))

last_obs <- gr_dat %>% 
  select(-birth_year, -month_capture, -year_capture) %>%
  group_by(fish_id) %>%
  mutate(last_obs = annuli_num == max(annuli_num)) %>%
  filter(last_obs)  %>%
  left_join(al_dat, "fish_id") %>% 
  mutate(month_capture = 
           factor(month_capture, levels = format(ISOdate(2000, 1:12, 1), "%B"))
  )

##--------------------------------------------------------------------------------


##--------------------------------------------------------------------------------
## Bootstrapping
##--------------------------------------------------------------------------------

in_dat <- 
  gr_dat %>% 
  filter(annuli_num == 1) %>%
  select(fish_id, OL.t, covariate) %>%
  na.omit

gr_dat <- 
  gr_dat %>% 
  select(fish_id, OL.t1, OL.t, covariate) %>%
  filter(fish_id %in% al_dat$fish_id) %>%
  na.omit

al_dat <- 
  al_dat %>% 
  select(fish_id, fish_length_mm, rad4,age, month_capture) %>%
  na.omit

al_dat_1 <- filter(al_dat,  fish_id %in% gr_dat$fish_id)
al_dat_2 <- filter(al_dat, !fish_id %in% gr_dat$fish_id)

n_boots <- 10000
in_ests <- gr_ests <- al_ests <- list()
for (i in 1:n_boots) {
  
  # resample the growth data + fit and store the init size model coefs
  in_dat_now <- sample_n(in_dat, size = nrow(in_dat), replace = TRUE)
  im <- lm(OL.t ~ covariate, data=in_dat_now)
  in_ests[[i]] <- as.list(c(imod=coef(im), imod.sd=summary(im)$sigma))
  
  # resample the growth data + fit and store the growth model coefs
  gr_dat_now <- sample_n(gr_dat, size = nrow(gr_dat), replace = TRUE)
  gm <- lm(OL.t1 ~ OL.t + covariate, data=gr_dat_now)
  gr_ests[[i]] <- as.list(c(gmod=coef(gm), gmod.sd=summary(gm)$sigma))
  
  # fit and store the allometry model coefficients
  al_dat_now <- sample_n(al_dat, size = nrow(al_dat), replace = TRUE)
  am <- gls(fish_length_mm ~ rad4, weights=varExp(form= ~ rad4), data=al_dat_now)
  al_ests[[i]] <- as.list(c(smod=coef(am), 
                            smod.sd.p0=summary(am)$sigma,
                            smod.sd.p=as.numeric(am$modelStruct$varStruct)))
}

mpars_bs <- bind_cols(bind_rows(gr_ests),
                      bind_rows(in_ests),
                      bind_rows(al_ests))
mpars_all <- apply(mpars_bs, 2, mean)

#models
im <- lm(OL.t ~ covariate, data=in_dat)
gm <- lm(OL.t1 ~ OL.t + covariate, data=gr_dat)
am <- gls(fish_length_mm ~ rad4, weights=varExp(form= ~ rad4), data=al_dat)
summary(im)
summary(gm)

##--------------------------------------------------------------------------------


##--------------------------------------------------------------------------------
## Projection functions
##--------------------------------------------------------------------------------

szkern <- function(TL, z, pars) {
  TLmean <- pars["smod.(Intercept)"] + pars["smod.rad4"]*z
  TLsd <- pars["smod.sd.p0"]*sqrt(exp(2*pars["smod.sd.p"]*z))
  return( dnorm(TL, mean = TLmean, sd = TLsd))
}

gkern <- function(z1, z, T, pars) {
  z1mean <- pars["gmod.(Intercept)"] + pars["gmod.OL.t"]*z + pars["gmod.covariate"]*T
  z1sd <- pars["gmod.sd"]
  return( dnorm(z1, mean = z1mean, sd = z1sd) )
}

idist <- function(z1, T, pars) {
  z1mean <- pars["imod.(Intercept)"] + pars["imod.covariate"]*T
  z1sd <- pars["imod.sd"]
  return( dnorm(z1, mean = z1mean, sd = z1sd))
}

mkdelta <- function(ipars) {
  (ipars["u"]-ipars["l"])/ipars["n"]
}

mkmesh <- function(ipars) {
  ipars["l"] + mkdelta(ipars) * (1:ipars["n"]-1/2)
}

get_quantile <- function(pdens, dlta, mesh, quantile = 0.95) {
  cdf <- cumsum(pdens) / sum(pdens) 
  range(mesh[cdf > (1-quantile)/2 & cdf < 1/2 + quantile/2])
}

# Otolith projection function  -----------------------

projectO <- function(iparsO, mpars, max_a, ref_temp) {
  mesh <- mkmesh(iparsO)
  dlta <- mkdelta(iparsO)
  gK <- outer(mesh, mesh, gkern, T=ref_temp, pars=mpars) * dlta
  n <- matrix(NA, nrow=length(mesh), ncol=max_a)
  n[,1] <- idist(mesh, T=ref_temp, pars=mpars)
  for (A in 2:max_a) {
    n[,A] <- (gK %*% n[,A-1])[,,drop=TRUE]
  }
  colnames(n) <- paste0("age_", 1:max_a - 1)
  return(list(n = n, # <- projections 
              ref_temp = ref_temp, # <- environment
              max_a = max_a, mesh = mesh, dlta = dlta)) # <- integr. params
}


# Body size projection function  ---------------------

projectB <- function(iparsO, iparsB, mpars, max_a, ref_temp) {
  meshB <- mkmesh(iparsB)
  dltaB <- mkdelta(iparsB)
  meshO <- mkmesh(iparsO)
  dltaO <- mkdelta(iparsO)
  
  # making the growth (# 1) and size (# 2) kernels 
  gK <- outer(meshO, meshO,  gkern, T=ref_temp, pars=mpars) *  dltaO # 1
  bK <- outer(meshB, meshO, szkern, pars=mpars) *  dltaO # 2
  
  # creating empty matrix to put output into
  nO <- matrix(NA, nrow=length(meshO), ncol=max_a)
  nB <- matrix(NA, nrow=length(meshB), ncol=max_a)
  
  # otolith size (nO) and fish body size (nB) distributions of age0 individuals:
  nO[,1] <- idist(meshO, T=ref_temp, pars=mpars)
  nB[,1] <- (bK %*% nO[,1])[,,drop=TRUE]
  
  # extend this distribution to all ages(using growth kernel and body kernel)
  for (A in 2:max_a) {
    nO[,A] <- (gK %*% nO[,A-1])[,,drop=TRUE]
    nB[,A] <- (bK %*% nO[,A]  )[,,drop=TRUE]
  }
  colnames(nB) <- paste0("age_", 1:max_a - 1)
  colnames(nO) <- paste0("age_", 1:max_a - 1)
  return(list(nB = nB, # <- projections 
              ref_temp = ref_temp, # <- environment
              max_a = max_a, meshB = meshB, dltaB = dltaB)) # <- integr. params
}


##--------------------------------------------------------------------------------


##--------------------------------------------------------------------------------
## Projections - Otoliths only, using mean paramater values 
##--------------------------------------------------------------------------------

# Simple mean projection  -----------------------------

proj_mean <-
  projectO(iparsO   = c(l=0.6, u=4, n=500), 
           mpars    = apply(mpars_bs, 2, mean),
           ref_temp = mean(unique(useData$covariate)),
           max_a    = 8)

## plot the expected age specific density functions 
par(mfrow=c(1,1))
with(proj_mean, {
  plot(mesh, n[,1], type="l", ylim=c(0,4), 
       xlab="Otolith Size", ylab="Density")
  for (A in 2:max_a) lines(mesh, n[,A]) 
})

## plot the size vs age relationship...
par(mfrow=c(1,1))
with(proj_mean, {
  mean_sz <- apply(n, 2, function(pdens) sum(pdens * mesh * dlta))
  plot(1:max_a-1, mean_sz, 
       xlab="Age", ylab="Otolith Size", ylim=c(1.2,3.0))
  ## ...and add the 90% quantiles
  quantile90 <- apply(n, 2, get_quantile, dlta, mesh)
  for (A in 1:max_a-1) lines(rep(A, length(quantile90[[A+1]])), quantile90[[A+1]])
})

## repeat the mean calculations for different temperatures (+/- 1 degree)
proj_mean_m1 <-
  projectO(iparsO   = c(l=0.6, u=4, n=500), 
           mpars    = apply(mpars_bs, 2, mean),
           ref_temp = mean(unique(useData$covariate)) - 1,
           max_a    = 8)
proj_mean_p1 <-
  projectO(iparsO   = c(l=0.6, u=4, n=500), 
           mpars    = apply(mpars_bs, 2, mean),
           ref_temp = mean(unique(useData$covariate)) + 1,
           max_a    = 8)

## plot the size vs age relationship...
par(mfrow=c(1,1))
with(proj_mean, {
  mean_sz <- apply(n , 2, function(pdens) sum(pdens*mesh*dlta))
  plot(1:max_a-1, mean_sz,  ylim=c(1.2,3.0),
       xlab="Age", ylab="Otolith Size", type="b")
})
with(proj_mean_m1, {
  mean_sz <- apply(n, 2, function(pdens) sum(pdens*mesh*dlta))
  lines(1:max_a-1, mean_sz, type="b", col="blue")
})
with(proj_mean_p1, {
  mean_sz <- apply(n, 2, function(pdens) sum(pdens*mesh*dlta))
  lines(1:max_a-1, mean_sz, type="b", col="red")
})

##--------------------------------------------------------------------------------
## Projections - Otoliths only, full bootstrap
##--------------------------------------------------------------------------------

# otolith projection  --------------------------------

mpars_bs_mat <- as.matrix(mpars_bs)

# because the second argument in the function is 1, this is basically running the projection (projectO) n.boot*n times
# for each row (boot) of mpars
all_proj <- apply(mpars_bs_mat, 1, function(mpars) {
  projectO(iparsO   = c(l=0.1, u=4, n=500), 
           mpars    = mpars,
           ref_temp = mean(unique(useData$covariate)),
           max_a    = 8)
})


# get the bootstrapped set of mean otolith sizes
# because the second argument is 2, this is applied n.boot times
mean_sz <- lapply(all_proj, function(proj) {
  with(proj, {
    apply(n, 2, function(pdens) sum(pdens * mesh * dlta)) %>%
      as.list %>% as_data_frame
  })
}) %>% bind_rows

# get the prediction intervals for each age class
pred_dist <- lapply(all_proj, function(x) x$n)
pred_dist <- Reduce("+", pred_dist)
pred_interv <- with(all_proj[[1]], {
  apply(pred_dist, 2, get_quantile, dlta, mesh)
})


# Age-specific density plot --------------------------
max_a <- 8
mesh <- mkmesh(c(l = 0.6, u = 4.0 , n=500 )) 
## plot the age specific density functions
pred_dist <- lapply(all_proj, function(x) x$n)
pred_dist2 <- Reduce("+", pred_dist)/ length(pred_dist) # mean over all nB
plot(mesh, pred_dist2[,1], type="l", ylim=c(0,4), xlim=c(1,3.5),
     xlab="Otolith radius (mm)", ylab="Density",
     cex.lab=1.5, cex.axis=1.5)
for (A in 2:max_a) lines(mesh, pred_dist2[,A])
# ----------------------------------------------------


# Otolith vs. Age plot -------------------------------

# empty plot
par(mar=c(5,5,1,1))
plot(1, 1, type = "n", pch = 20,
     xlab="Age (Years)", ylab="Otolith Size (mm)",
     xlim = c(0, 8), ylim=c(1.0, 3.5),
     cex.lab=1.5, cex.axis=1.5)

# add oberved data
with(useData, 
     points(jitter(annuli_num - 1), OL.t, 
            pch = 20, cex = 0.3, col = grey(0.3, alpha = 0.4)))

# add the prediction intervals
for (i in 1:ncol(pred_interv)) 
  lines(rep(i, 2)-1, pred_interv[,i], 
        col = rgb(1, 0, 0, alpha = 0.2), lwd = 6)

# now add the average size +/- CI
y_vals <- apply(mean_sz, 2, mean)
points(1:length(y_vals)-1, y_vals, pch = 20, cex = 1.1, col = "red")
y_vals <- apply(mean_sz, 2, quantile, probs = c(0.025, 0.975))
for (i in 1:ncol(y_vals)) lines(rep(i, 2)-1, y_vals[,i], lwd = 1.5, col = "red")
# ----------------------------------------------------

##--------------------------------------------------------------------------------

##--------------------------------------------------------------------------------
## Projections - Otoliths and Body size
##--------------------------------------------------------------------------------


# Full projection ------------------------------------

mpars_bs_mat <- as.matrix(mpars_bs)
iparsO <- c(l = 0.6, u = 4.0 , n=500 ) # otolith 
iparsB <- c(l = 0.0, u = 800, n=1000) # body
meshO <- mkmesh(iparsO) # otolith mesh
meshB <- mkmesh(iparsB) # body mesh
max_a <- 8

full_proj <- apply(mpars_bs_mat, 1, function(mpars) {
  projectB(iparsB   = iparsB,
           iparsO   = iparsO,
           mpars    = mpars,
           ref_temp = mean(unique(useData$covariate)),
           max_a    = max_a)
})

# ----------------------------------------------------

# Body size vs. Age plot -----------------------------

proj_plot <- function(full_proj){
  
  
  # empty plot
  par(mar=c(5,5,1,1))
  plot(1, 1, type = "n", pch = 20,
       xlab="Age (Years)", ylab="Fish Total Length (mm)",
       xlim = c(0, 8), ylim=c(0, 600),
       cex.lab=1.5, cex.axis=1.5)
  
  # add observed data points 
  with(al_dat, 
       points(jitter(age), fish_length_mm, 
              pch = 20, cex = 0.3, col = grey(0.3, alpha = 0.4)))
  
  # add the prediction intervals
  for (i in 1:ncol(pred_intervB)) 
    lines(rep(i, 2)-1, pred_intervB[,i], 
          col = rgb(1, 0, 0, alpha = 0.2), lwd = 6)
  
  # now add the average size +/- CI
  y_vals <- apply(mean_szB, 2, mean)
  points(1:length(y_vals)-1, y_vals, pch = 20, cex = 1.1, col = "red")
  y_vals <- apply(mean_szB, 2, quantile, probs = c(0.025, 0.975))
  for (i in 1:ncol(y_vals)) lines(rep(i, 2)-1, y_vals[,i], lwd = 1.5, col = "red")
  
}

# get the bootstrapped set of mean body sizes
mean_szB <- lapply(full_proj, function(proj) {
  with(proj, {
    apply(nB, 2, function(pdens) sum(pdens * meshB * dltaB)) %>%
      as.list %>% as_data_frame
  })
}) %>% bind_rows

# get the prediction intervals for each age class
pred_distB <- lapply(full_proj, function(x) x$nB)
pred_distB <- Reduce("+", pred_distB)
pred_intervB <- with(full_proj[[1]], {
  apply(pred_distB, 2, get_quantile, dltaB, meshB)
})

# ----------------------------------------------------

# Age-specific density plot --------------------------

meshB <- mkmesh(c(l = 0.0, u = 800, n=1000)) 
## plot the age specific density functions
pred_distB <- lapply(full_proj, function(x) x$nB)
pred_distB2 <- Reduce("+", pred_distB)/ length(pred_distB) # mean over all nB
plot(meshB, pred_distB2[,1], type="l", ylim=c(0,0.02), xlab="Fish total length (mm)", ylab="Density",
     cex.lab=1.5, cex.axis=1.5)
for (A in 2:max_a) lines(meshB, pred_distB2[,A])

# ----------------------------------------------------

# TL vs. Age (temp scenarios) plots ------------------

mpars_bs_mat <- as.matrix(mpars_bs)
iparsO <- c(l = 0.6, u = 4.0 , n=500 ) # otolith 
iparsB <- c(l = 0.0, u = 800, n=1000) # body
meshO <- mkmesh(iparsO) # otolith mesh
meshB <- mkmesh(iparsB) # body mesh
max_a <- 8

full_proj_p1 <- apply(mpars_bs_mat, 1, function(mpars) {
  projectB(iparsB   = iparsB,
           iparsO   = iparsO,
           mpars    = mpars,
           ref_temp = mean(unique(useData$covariate)+1),
           max_a    = max_a)
})

full_proj_m1 <- apply(mpars_bs_mat, 1, function(mpars) {
  projectB(iparsB   = iparsB,
           iparsO   = iparsO,
           mpars    = mpars,
           ref_temp = mean(unique(useData$covariate)-1),
           max_a    = max_a)
})

proj_plot(full_proj_m1)
proj_plot(full_proj_p1)


mean_szB_p1 <- lapply(full_proj_p1, function(proj) {
  with(proj, {
    apply(nB, 2, function(pdens) sum(pdens * meshB * dltaB)) %>%
      as.list %>% as_data_frame
  })
}) %>% bind_rows

mean_szB_m1 <- lapply(full_proj_m1, function(proj) {
  with(proj, {
    apply(nB, 2, function(pdens) sum(pdens * meshB * dltaB)) %>%
      as.list %>% as_data_frame
  })
}) %>% bind_rows


# get the prediction intervals for each age class
pred_distB_p1 <- lapply(full_proj_p1, function(x) x$nB)
pred_distB_p1 <- Reduce("+", pred_distB_p1)
pred_distB_m1 <- lapply(full_proj_m1, function(x) x$nB)
pred_distB_m1 <- Reduce("+", pred_distB_m1)

pred_intervB_m1 <- with(full_proj_m1[[1]], {
  apply(pred_distB_m1, 2, get_quantile, dltaB, meshB)
})
pred_intervB_p1 <- with(full_proj_p1[[1]], {
  apply(pred_distB_p1, 2, get_quantile, dltaB, meshB)
})


# empty plot
par(mar=c(5,5,1,1))
plot(1, 1, type = "n", pch = 20,
     xlab="Age (Years)", ylab="Fish Total Length (mm)",
     xlim = c(0, 8), ylim=c(0, 600),
     cex.lab=1.5, cex.axis=1.5)

# add observed data points 
with(al_dat, 
     points(jitter(age), fish_length_mm, 
            pch = 20, cex = 0.3, col = grey(0.3, alpha = 0.4)))

# add the prediction intervals
for (i in 1:ncol(pred_intervB)) 
  lines(rep(i, 2)-1, pred_intervB[,i], 
        col = rgb(0, 1, 0, alpha = 0.2), lwd = 6)
for (i in 1:ncol(pred_intervB_p1)) 
  lines(rep(i+0.1, 2)-1, pred_intervB_p1[,i], 
        col = rgb(1, 0, 0, alpha = 0.2), lwd = 6)
for (i in 1:ncol(pred_intervB_m1)) 
  lines(rep(i-0.1, 2)-1, pred_intervB_m1[,i], 
        col = rgb(0, 0, 1, alpha = 0.2), lwd = 6)


# now add the average size +/- CI
y_vals <- apply(mean_szB, 2, mean)
points(1:length(y_vals)-1, y_vals, pch = 20, cex = 1.1, col = "dark green")
y_vals <- apply(mean_szB, 2, quantile, probs = c(0.025, 0.975))
for (i in 1:ncol(y_vals)) lines(rep(i, 2)-1, y_vals[,i], lwd = 1.5, col = "dark green")

y_vals <- apply(mean_szB_p1, 2, mean)
points((1:length(y_vals)-1)+0.1, y_vals, pch = 20, cex = 1.1, col = "red")
y_vals <- apply(mean_szB_p1, 2, quantile, probs = c(0.025, 0.975))
for (i in 1:ncol(y_vals)) lines(rep(i+0.1, 2)-1, y_vals[,i], lwd = 1.5, col = "red")

y_vals <- apply(mean_szB_m1, 2, mean)
points((1:length(y_vals)-1)-0.1, y_vals, pch = 20, cex = 1.1, col = "blue")
y_vals <- apply(mean_szB_m1, 2, quantile, probs = c(0.025, 0.975))
for (i in 1:ncol(y_vals)) lines(rep(i-0.1, 2)-1, y_vals[,i], lwd = 1.5, col = "blue")

useData$CaptureM <- match(useData$month.capture,month.name)
useData$AgeM <- ((useData$annuli_num)*12) + useData$CaptureM - 1

al_dat$CaptureM <- match(al_dat$month_capture,month.name)
al_dat$AgeM <- (al_dat$age*12) + al_dat$CaptureM - 1

# --------------------------------------------------------------------------------
