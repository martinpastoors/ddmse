# ------------------------------------------------------------------------------
# report.R
#
# Generating DDMSE results to be included into Rmd files
#
# Laurie Kell and Martin Pastoors
#
# 09/02/2023 First coding
# ------------------------------------------------------------------------------

rm(list=ls())

iFig=0
iTab=0

library(FLCore)    # install.packages("FLCore", repos="http://flr-project.org/R")
library(FLBRP)     # install.packages("FLBRP", repos="http://flr-project.org/R")
library(FLasher)   # install.packages("FLasher", repos="http://flr-project.org/R")
library(FLAssess)  # install.packages("FLAssess", repos="http://flr-project.org/R")
library(FLife)     # install.packages("FLife", repos="http://flr-project.org/R")
library(ggplotFL)

library(popbio)

library(GGally)
library(ggpubr)

library(rjson) # install.packages("rjson"); used for dropbox function
library(RJSONIO)

library(plyr)
library(dplyr)
library(reshape)

library(readxl)
  
library(mgcv)
library(mgcViz)

library(tidyverse)

theme_set(theme_bw(16))

source("R/dd.R")
source("R/get_dropbox.r")
source("R/pm.R")

dropboxdir<-try(file.path(get_dropbox(), "pelagics"))

if ("try-error" %in% is(dropboxdir))
  dropboxdir="~/Dropbox/pelagics"

myparams <- data.frame(
  stock = c("mac"            ,  "whb"            ),
  m1    = c(0.15/(0.2^-0.288),  0.2/(300^-0.288) ),
  m2    = c(-0.28            ,  -0.28            ),
  m3    = c(NA               ,  NA               ),
  blim  = c(2000000          , 1500000           ),
  btrig = c(2580000          , 2250000           ),
  bmsy  = c(3500000          ,  NA               ),
  b     = c(0.3              ,  NA               ),
  k     = c(30.1             ,  NA               ), 
  w50   = c(0.2              ,  NA               ),
  wt1   = c(0.01             ,  NA               ))

mystk     <- "mac"; mystkname <- "Northeast Atlantic mackerel"; mylatin <- "Scomber scombrus" 
# mystk     <- "whb"; mystkname <- "Blue whiting";                mylatin <- "Micromesistius poutassou" 
maxyear <- 2050


# ------------------------------------------------------------------------------
# Operating model
# ------------------------------------------------------------------------------
resultsdir <- file.path(dropboxdir, "results", mystk, "figures")

btrig_blim <- myparams[myparams$stock==mystk,"btrig"]/myparams[myparams$stock==mystk,"blim"] #mac 1.29; whb 1.5
  
load(file.path(dropboxdir,paste0("data/om/",mystk, ".RData")))

combY <- as.data.frame(pmYear(window(get(mystk), start=2001)),drop=T) %>% mutate(stk=mystk, type="SAM non-DD") 

par=FLPar(c(m1  = myparams[myparams$stock==mystk,"m1"],
            m2  = myparams[myparams$stock==mystk,"m2"],
            bmsy= myparams[myparams$stock==mystk,"bmsy"],
            b   = myparams[myparams$stock==mystk,"b"],
            k   = myparams[myparams$stock==mystk,"k"], 
            w50 = myparams[myparams$stock==mystk,"m1"])) 

# if (is.numeric(myparams[myparams$stock==mystk,"wt1"])) {
#   stock.wt(   get(mystk))[1] = myparams[myparams$stock==mystk,"wt1"]
#   catch.wt(   get(mystk))[1] = myparams[myparams$stock==mystk,"wt1"]
#   landings.wt(get(mystk))[1] = myparams[myparams$stock==mystk,"wt1"]
#   discards.wt(get(mystk))[1] = myparams[myparams$stock==mystk,"wt1"]
# }

if (is.numeric(myparams[myparams$stock==mystk,"wt1"])) {
  stock.wt(   mac)[1] = myparams[myparams$stock==mystk,"wt1"]
  catch.wt(   mac)[1] = myparams[myparams$stock==mystk,"wt1"]
  landings.wt(mac)[1] = myparams[myparams$stock==mystk,"wt1"]
  discards.wt(mac)[1] = myparams[myparams$stock==mystk,"wt1"]
}

# Density dependent SAM
om_orig <- ddFn2(get(mystk), par); 
combY   <- bind_rows(combY,
                     as.data.frame(pmYear(window(om_orig, start=2001)),drop=T) %>% mutate(stk=mystk, type="SAM-DD")) 

# VPA assessment; not density dependent
om    <- om_orig+VPA(om_orig)
combY <- bind_rows(combY, 
                   as.data.frame(pmYear(window(om, start=2001)),drop=T) %>% mutate(stk=mystk, type="VPA non-DD") )

om    <- ddFn2(om,par)
combY <- bind_rows(combY, 
                   as.data.frame(pmYear(window(om, start=2001)),drop=T) %>% mutate(stk=mystk, type="VPA DD") )

p <-
  combY %>% 
  filter(qname %in% c("catch","ssb","rec","fbar")) %>% 
  mutate(dd = ifelse(grepl("non-DD", type), "non-DD", "DD")) %>%
  mutate(dd = factor(dd, levels=c("non-DD", "DD"))) %>% 
  ggplot(aes(x=year, y=median)) +
  theme_bw() +
  geom_ribbon(aes(fill=type, ymin=lowq, ymax=uppq, xmin=year, xmax=year), alpha=0.3) +
  geom_line(aes(colour=type)) +
  labs(title=paste(mystk, mystkname)) +
  facet_grid(qname~dd, scales="free_y") 

ggsave(p, 
       filename=file.path(resultsdir, paste(mystk, "om.jpg", sep="_")), 
       device="jpeg", 
       width=10, height=10, units="in")

# ------------------------------------------------------------------------------
# Density dependence: assumptions and estimation
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Management procedures
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Simulation results
# ------------------------------------------------------------------------------

# Forward projections with F=fmsy (fwd), 0.8xfmsy?, 1.0xfmsy, 1.2xfmsy?

sr=fmle(as.FLSR(om,model="bevholt"),control=list(silent=TRUE))
eq=FLBRP(om,sr=sr)

om=fwdWindow(om,end=maxyear,eq)
f =fbar(om)[,ac(2021:maxyear)]%=%refpts(eq)["msy","harvest"]
om=fwd(om,fbar=f,sr=eq)

om2=om
om3=om
om4=om

for (iYear in ac(2021:maxyear)){
  om2=ddFn(iYear,om2,par)
  om2=fwd(om2,fbar=f[,iYear]*0.8,sr=eq)
  
  om3=ddFn(iYear,om3,par)
  om3=fwd(om3,fbar=f[,iYear],sr=eq)
  
  om4=ddFn(iYear,om4,par)
  om4=fwd(om4,fbar=f[,iYear]*1.2,sr=eq)
}

plot(FLStocks("fwd"=om,"0.8"=om2,"1.0"=om3,"1.2"=om4))


