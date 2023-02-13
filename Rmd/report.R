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
  wt1   = c(0.01             ,  NA               ),
  maxage= c(12               ,  10               ))

save(myparams, file=file.path(dropboxdir, "data", "inputs", "myparams.RData"))

mystk     <- "mac"; mystkname <- "Northeast Atlantic mackerel"; mylatin <- "Scomber scombrus" 
# mystk     <- "whb"; mystkname <- "Blue whiting";                mylatin <- "Micromesistius poutassou" 
maxyear <- 2050

resultsdir <- file.path(dropboxdir, "results", mystk, "figures")

# ------------------------------------------------------------------------------
# Operating model
# ------------------------------------------------------------------------------

maxage <-  myparams[myparams$stock==mystk,"maxage"]

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
# Reference points
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Density dependence: assumptions and estimation
# ------------------------------------------------------------------------------

# HOW DO THESE PARAMETERS RELATE TO MYPARAMS; e.g. b??

# Get theoretical values from literature
par=FLife::teleost[c("linf","k","t0","l50","a","b"),dimnames(teleost)$iter==mylatin]

# Uses life history theory to derive parameters for biological relationships
par=FLife::lhPar(par,
                 m1=myparams[myparams$stock==mystk,"m1"],
                 m2=myparams[myparams$stock==mystk,"m2"],
                 m3=myparams[myparams$stock==mystk,"m3"])

# WHAT IS THE ASSUMED SELECTIVITY AND WHERE DOES IT COME FROM?

# derive a FLBRP object with Lorenzen natural mortality
eq =FLife::lhEql(par,m="lorenzen")

# plot the BRP object
ggsave(ggplotFL::plot(eq) + 
         theme_bw() +
         labs(title=paste(mystk, mystkname, "FLBRP object with DD M")), 
       filename=file.path(resultsdir, paste(mystk, "brp.jpg", sep="_")), 
       device="jpeg", 
       width=10, height=10, units="in")

# convert to FLStock and project
om=as(eq,"FLStock")
om=fwd(om,fbar=fbar(om)[,-1],sr=eq)
ggsave(ggplot(model.frame(FLQuants(Biomass=biomass(eq),Yield=catch(eq))))+
         geom_line(aes(Biomass,Yield)), 
       filename=file.path(resultsdir, paste(mystk, "eq.jpg", sep="_")), 
       device="jpeg", 
       width=10, height=10, units="in")


# Simulation of density dependence in mass-at-age for different levels of biomass relative to Bmsy
dat=ddWt(stock.wt(eq),biomass(eq),refpts(eq)["msy","biomass"]) 
dat=merge(as.data.frame(dat,drop=T),
          as.data.frame(biomass(eq)%/%refpts(eq)["msy","biomass"],drop=T),
          by="year")
names(dat)[3:4]=c("mass","biomass")

ggsave(ggplot(dat)+    
         theme_bw() +
         geom_line(aes(age,mass,col=biomass,group=year))+
         scale_x_continuous(limits=c(0,maxage), breaks=(seq(0,maxage,1)))+
         # scale_y_continuous(limits=c(0,800))+
         geom_line(aes(age,data),data=as.data.frame(stock.wt(eq)),col="red",size=2)+
         labs(x="Age", y="Mass-at-age", colour="Biomass/Bmsy") , 
       filename=file.path(resultsdir, paste(mystk, "sim_dd_mass.jpg", sep="_")), 
       device="jpeg", 
       width=10, height=10, units="in")

# Simulation of density dependence in maturity-at-age for different levels of biomass relative to Bmsy
# WHAT IS K HERE?
parMat=FLPar(c(k=0.1, w5=par["a"]*par["l50"]^par["b"]))
dat=transform(dat,mat=1/(1+exp(-parMat[1]*(mass-parMat[2])))) 

ggsave(ggplot(dat)+           
         theme_bw() +
         geom_line(aes(age,mat,col=biomass,group=year))+
         scale_x_continuous(limits=c(0,maxage), breaks=(seq(0,maxage,1)))+
         geom_line(aes(age,data),data=as.data.frame(mat(eq)),col="red",size=2)+
         labs(x="Age", y="Maturity-at-age", colour="Biomass/Bmsy"), 
       filename=file.path(resultsdir, paste(mystk, "sim_dd_maturity.jpg", sep="_")), 
       device="jpeg", 
       width=10, height=10, units="in")

# Simulation of density dependence in M-at-age for different levels of biomass relative to Bmsy
dat=transform(dat,m=c(par["m1"])*(mass^c(par["m2"])))

ggsave(ggplot(dat)+           
         geom_line(aes(age,m,col=biomass,group=year))+
         scale_x_continuous(limits=c(0,maxage), breaks=(seq(0,maxage,1)))+
         geom_line(aes(age,data),data=as.data.frame(m(eq)),col="red",size=2)+
         labs(x="Age", y="M-at-age", colour="Biomass/Bmsy"), 
       filename=file.path(resultsdir, paste(mystk, "sim_dd_m.jpg", sep="_")), 
       device="jpeg", 
       width=10, height=10, units="in")


# Comparison of equilibrium yields curves without DD and with DD mass

x         = propagate(eq,length(dimnames(eq)$year))
fbar(x)   = fbar(x)[,1]
fbar(x)[] = c(fbar(eq))*2

# landing weight at age and discard weight at age relative to stock weight
lsr       = landings.wt(eq)/stock.wt(eq)
dsr       = discards.wt(eq)/stock.wt(eq)

# iterate
for (i in seq(10)){
  stock.wt(   x)=ddWt(stock.wt(eq),biomass(x),refpts(eq)["msy","biomass"],alpha=-0.2)
  landings.wt(x)=stock.wt(x)*lsr
  discards.wt(x)=stock.wt(x)*dsr
  x=brp(x)}

df <- 
  bind_rows(
    model.frame(FLQuants(eq,biomass=function(x) biomass(x), catch=function(x) catch(x)),drop=T) %>% mutate(scen="no-DD"),
    model.frame(FLQuants(x, biomass=function(x) biomass(x), catch=function(x) catch(x)),drop=T) %>% mutate(scen="DD mass")
  )

# Comparison of equilibrium yields curves with DD mass+maturity

x=propagate(eq,length(dimnames(eq)$year))
fbar(x)[,1]=c(seq(0,max(fbar(eq)),length.out=51),seq(1,50,length.out=51)*max(fbar(eq)))[-52]
fbar(x)=fbar(x)[,1]

lsr=landings.wt(eq)/stock.wt(eq)
dsr=discards.wt(eq)/stock.wt(eq)
# iterate
for (i in seq(10)){
  stock.wt(   x)=ddWt(stock.wt(eq),biomass(x),refpts(eq)["msy","biomass"],alpha=-0.2)
  landings.wt(x)=stock.wt(x)*lsr
  discards.wt(x)=stock.wt(x)*dsr
  
  mat(     x)=1/(1+exp(-parMat[1]*(stock.wt(x)-parMat[2])))
  x=brp(x)
}

df <- 
  bind_rows(
    df,
    model.frame(FLQuants(x, biomass=function(x) biomass(x), catch=function(x) catch(x)),drop=T) %>% mutate(scen="DD mass+mat")
  )

# Comparison of equilibrium yields curves with DD mass+mat+M 

x=propagate(eq,length(dimnames(eq)$year))
fbar(x)[,1]=c(seq(0,max(fbar(eq)),length.out=51),seq(1,50,length.out=51)*max(fbar(eq)))[-52]
fbar(x)=fbar(x)[,1]

lsr=landings.wt(eq)/stock.wt(eq)
dsr=discards.wt(eq)/stock.wt(eq)
# iterate
for (i in seq(10)){
  stock.wt(   x)=ddWt(stock.wt(eq),biomass(x),refpts(eq)["msy","biomass"],alpha=-0.2)
  landings.wt(x)=stock.wt(x)*lsr
  discards.wt(x)=stock.wt(x)*dsr
  
  mat(     x)=1/(1+exp(-parMat[1]*(stock.wt(x)-parMat[2])))
  x=brp(x)
}

df <- 
  bind_rows(
    df,
    model.frame(FLQuants(x, biomass=function(x) biomass(x), catch=function(x) catch(x)),drop=T) %>% mutate(scen="DD mass+mat+M")
  )
ggplot(data=df, aes(biomass, catch))+
  theme_bw() +
  geom_line(aes(colour=scen)) +
  labs(x="Total biomass", y="Catch", colour="") +
  scale_colour_manual(values =c("no-DD"         = "black",
                                "DD mass"       = "red", 
                                "DD mass+mat"   = "blue",
                                "DD mass+mat+M" = "darkgreen")) +
  facet_wrap(~scen)


ggsave(ggplot(data=df, aes(biomass, catch))+
         theme_bw() +
         geom_line(aes(colour=scen)) +
         labs(x="Total biomass", y="Catch", colour="") +
         scale_colour_manual(values =c("no-DD"         = "black",
                                       "DD mass"       = "red", 
                                       "DD mass+mat"   = "blue",
                                       "DD mass+mat+M" = "darkgreen")), 
       filename=file.path(resultsdir, paste(mystk, "eq_yield.jpg", sep="_")), 
       device="jpeg", 
       width=10, height=10, units="in")



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


