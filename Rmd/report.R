 # ------------------------------------------------------------------------------
# report.R
#
# Generating DDMSE results to be included into Rmd files
#
# Laurie Kell and Martin Pastoors
#
# 09/02/2023 First coding
# 13/02/2023 Implemented all codes from cond, DD and OM
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
library(ggrepel)   # labelling points on a graph

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
  stock    = c("mac"                       , "whb"            ),
  mystkname=c("Northeast Atlantic mackerel", "Blue whiting"   ),
  mylatin = c("Scomber scombrus"           ,"Micromesistius poutassou"),
  m1      = c(0.15/(300^-0.288)            , 0.2/(150^-0.288) ), # 300 g MAC, 150 g WHB 
  m2      = c(-0.288                        , -0.288          ),
  m3      = c(NA                           , NA               ),
  m       = c("lorenzen"                   , "lorenzen"       ),
  blim    = c(2000000                      , 1500000          ),
  btrig   = c(2580000                      , 2250000          ),
  bmsy    = c(3500000                      , 4000000          ), # guesstimates
  alpha   = c(-0.2                         , -0.2             ),
  delta   = c(-0.2                         , -0.2             ),
  matk    = c(30                           , 0.1              ), # source: conditioning
  wt1     = c(0.01                         , NA               ), # for resetting weight at age 1
  minage  = c(0                            , 1                ),
  maxage  = c(12                           , 10               ),
  maxyear = c(2050                         , 2050             ))

save(myparams, file=file.path(dropboxdir, "data", "inputs", "myparams.RData"))

mystk     <- "mac";  
# mystk     <- "whb";  
# for (mystk in c("mac","whb")) {
  
  mystkname  <- myparams[myparams$stock==mystk,"mystkname"]
  mylatin    <- myparams[myparams$stock==mystk,"mylatin"]
  maxyear    <- myparams[myparams$stock==mystk,"maxyear"]
  minage     <- myparams[myparams$stock==mystk,"minage"]
  maxage     <- myparams[myparams$stock==mystk,"maxage"]
  nages      <- maxage - minage + 1  
  
  figuresdir <- file.path(dropboxdir, "results", mystk, "figures")
  tablesdir  <- file.path(dropboxdir, "results", mystk, "tables")
  
  # ------------------------------------------------------------------------------
  # 1. Introduction
  # ------------------------------------------------------------------------------
  
  # ------------------------------------------------------------------------------
  # 2. Conditioning
  # ------------------------------------------------------------------------------
  
  # Get FLStock and FLBRP object  ----------------------------------------------
  
  load(file.path(dropboxdir,"data","om",paste0(mystk,".RData"))) # only used for FLBRP part; comes from 
  stk   = get(mystk)
  stkR  = get(paste0(mystk,"R"))
  

  if (mystk == "mac") {
    ts    = as.data.frame(readxl::read_excel(file.path(dropboxdir,"data/inputs/dat.xlsx"),sheet="ts"))
    
    stock.wt(   stk)[1] = myparams[myparams$stock==mystk,"wt1"]
    catch.wt(   stk)[1] = myparams[myparams$stock==mystk,"wt1"]
    landings.wt(stk)[1] = myparams[myparams$stock==mystk,"wt1"]
    discards.wt(stk)[1] = myparams[myparams$stock==mystk,"wt1"]
    mat(        stk)[1] = 0
  } else if (mystk == "whb") {
    
    ts=model.frame(FLQuants(stk,tb=stock,ssb=ssb),drop=T)
    ts$tb=ts$tb/1e6
  }
  
  wt     = transform(as.data.frame(stock.wt(stk),drop=TRUE),cohort=year-as.numeric(age))
  wt     = merge(wt,ts,by="year")
  wt     = ddply(wt,.(age), transform, rsdl =data/mean(data))
  wt     = ddply(wt,.(age), transform, rsdl_=rsdl/var(data)^0.5)
  wt$age = factor(wt$age,levels=sort(unique(wt$age)))
  wt     = arrange(wt, year, age)
  
  mat    = transform(as.data.frame(mat(stk),drop=TRUE),cohort=year-as.numeric(age))
  mat    = merge(mat,ts)
  mat    = ddply(mat,.(age), transform, rsdl =data/mean(data))
  mat    = ddply(mat,.(age), transform, rsdl_=rsdl/var(data)^0.5) 
  mat$age= factor(mat$age,levels=sort(unique(mat$age)))
  mat    = arrange(mat, year, age)
  
  m      = transform(as.data.frame(m(stk),drop=TRUE),cohort=year-as.numeric(age))  
  m      = merge(m,ts)
  m      = ddply(m,.(age), transform, rsdl =data/mean(data))
  m      = ddply(m,.(age), transform, rsdl_=rsdl/var(data)^0.5) 
  m$age  = factor(m$age,levels=sort(unique(m$age))) 
  m      = arrange(m, year, age)
  
  # plot weight at age ---------------------------------------------------------
  p <-
    ggplot(wt)+  
    geom_line(aes(year,data,group=cohort),linetype=3,colour="grey10")+
    geom_line(aes(year,data,group=age,col=age))+
    geom_point(aes(year,data,col=age))+
    ggrepel::geom_text_repel(data = tail(wt, length(unique(wt$age))), 
                             aes(year, data,label=age, colour=age),
                             hjust = "outward", vjust=0.5, direction = "y", min.segment.length = 100) +
    xlab("Year")+ylab("Mass-at-age")+
    theme(legend.position="none")+
    # guides(colour = guide_legend(nrow = 1)) +
    scale_color_manual("Age",values=rainbow(length(unique(wt$age))))
  
  ggsave(p, 
         filename=file.path(figuresdir, paste(mystk, "waa.jpg", sep="_")), 
         device="jpeg", 
         width=10, height=6, units="in")

  # plot maturity at age -------------------------------------------------------
  p <-
    ggplot(mat)+  
    geom_line(aes(year,data,group=cohort),linetype=3,colour="grey10")+
    geom_line(aes(year,data,group=age,col=age))+
    geom_point(aes(year,data,col=age))+
    ggrepel::geom_text_repel(data = tail(mat, length(unique(mat$age))), 
                             aes(year, data,label=age, colour=age),
                             hjust = "outward", vjust=0.5, direction = "y", min.segment.length = 100) +
    xlab("Year")+ylab("Maturity-at-age")+
    theme(legend.position="none")+
    # guides(colour = guide_legend(nrow = 1)) +
    scale_color_manual("Age",values=rainbow(length(unique(mat$age))))
  
  ggsave(p, 
         filename=file.path(figuresdir, paste(mystk, "mat.jpg", sep="_")), 
         device="jpeg", 
         width=10, height=6, units="in")


  # plot m at age --------------------------------------------------------------
  p <-
    ggplot(m)+  
    geom_line(aes(year,data,group=cohort),linetype=3,colour="grey10")+
    geom_line(aes(year,data,group=age,col=age))+
    geom_point(aes(year,data,col=age))+
    ggrepel::geom_text_repel(data = tail(m, length(unique(m$age))), 
                             aes(year, data,label=age, colour=age),
                             hjust = "outward", vjust=0.5, direction = "y", min.segment.length = 100) +
    xlab("Year")+ylab("M-at-age")+
    theme(legend.position="none")+
    # guides(colour = guide_legend(nrow = 1)) +
    scale_color_manual("Age",values=rainbow(length(unique(m$age))))
  
  ggsave(p, 
         filename=file.path(figuresdir, paste(mystk, "m.jpg", sep="_")), 
         device="jpeg", 
         width=10, height=6, units="in")
  
  # Density dependence in weights at age ---------------------------------------
  
  a=ggplot(wt)+
    geom_boxplot(aes(age,data))+
    xlab("Mass-at-age")+ylab("Age")
  b=ggplot(subset(wt,rsdl<3&rsdl>0))+
    geom_point(aes(tb,rsdl))+
    geom_smooth(aes(tb,rsdl),method="lm")+
    facet_wrap(age~.)+
    xlab("Total biomass")+ylab("Residual")
  p <-
    ggarrange(a,b ,label.y=1,  
            widths = c(4,6), heights = c(1,1,1),
            nrow   = 1,          ncol    = 2)+
    theme_bw(16)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())
  ggsave(p, 
         filename=file.path(figuresdir, paste(mystk, "dd_waa.jpg", sep="_")), 
         device="jpeg", 
         width=10, height=8, units="in")
  
  
  # GAM model and diagnostics of weight vs Total biomass -----------------------
  
  # drop age 0 if existing
  wt2 = subset(wt,age!=0)
  nn  = length(unique(wt2$age))
  
  # run gammV  
  gmr=mgcViz::gammV(log(data)~s(tb/mean(tb), bs="tp")+age, data=wt2)
  
  # plot gammV
  jpeg(filename=file.path(figuresdir, paste(mystk, "gamm.jpg", sep="_")),
       width=10, height=8, units="in", res=300)
  plot(gmr) + labs(x="Total biomass", y="Weight deviation")
  dev.off()  
  
  # summary table
  fileConn <-file(file.path(tablesdir, paste(mystk, "gamm summary.txt", sep="_")))
  summary(gmr) %>% capture.output() %>% writeLines(., con=fileConn)
  close(fileConn)
  
  # summary coefficients
  fileConn <-file(file.path(tablesdir, paste(mystk, "gamm coefficients.txt", sep="_")))
  as.data.frame(coefficients(gmr)) %>% capture.output() %>% writeLines(., con=fileConn)
  close(fileConn)
  
  # plot(exp(coefficients(gmr)[2:nn]))  # NOT NEEDED?
  
  # connect prediction to wt dataset  
  wt2=cbind(wt2,hat=predict(gmr))
  wt2=merge(wt2,data.frame(age=1:nn,
                           mean=exp(c(0,coefficients(gmr)[2:nn]))))
  
  # Linear   
  p <-
    ggplot(wt2)+
    geom_line(aes(tb,exp(hat)/mean,col=age))+
    geom_point(aes(tb,data/mean,col=age))+
    xlab("Total biomass")+ylab("Mass-Mean")
  
  ggsave(p, 
         filename=file.path(figuresdir, paste(mystk, "lm_waa.jpg", sep="_")), 
         device="jpeg", 
         width=10, height=8, units="in")
  
  # Linear model summary  
  fileConn <-file(file.path(tablesdir, paste(mystk, "lm summary.txt", sep="_")))
  summary(lm(hat~tb,data=transform(wt2,tb=tb/mean(tb)))) %>% capture.output() %>% writeLines(., con=fileConn)
  close(fileConn)
  
  # Density dependence in maturity at age --------------------------------------
  
  a=ggplot(mat)+
    geom_boxplot(aes(age,data))+
    xlab("Maturity-at-age")+ylab("Age")
  b=ggplot(mat)+
    geom_point(aes(tb,rsdl))+
    geom_smooth(aes(tb,rsdl),method="lm")+
    facet_wrap(age~.)+
    xlab("Biomass")+ylab("Residual")
  p <-
    ggarrange(a,b ,label.y=1,  
            widths = c(4,6), heights = c(1,1,1),
            nrow   = 1,          ncol    = 2)+
    theme_bw(16)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())

  # Estimating maturity ogive --------------------------------------------------
  
  matFn <-function(wt,x) 1/(1+exp(-x[1]*(wt-x[2])))
  fn    <-function(x,dat) sum((dat$mat-matFn(dat$wt,x))^2)
  
  dat=model.frame(FLQuants(stk,"wt"=stock.wt,"mat"=function(x) mat(x)))
  dat$age=factor(dat$age)

  matPar=optim(par=c(k=0.5,w50=.2), fn, dat=dat)
  save(matPar,file=file.path(dropboxdir,paste0("data/om/",mystk,"matPar.RData")))
 
  p <-
    ggplot(dat,aes(wt,mat))+
    geom_point(aes(wt,mat,col=age))+
    geom_smooth(se=F,span=0.2)+
    xlab("Mass-at-age")+ylab("Matutrity-at-age")+
    geom_line(aes(wt,mat),col="red",data=data.frame(mat=matFn(dat$wt,matPar$par),wt=dat$wt))

  # print(paste("k=", matPar$par[1]))
  # print(paste("W50=", matPar$par[2]))

  ggsave(p, 
         filename=file.path(figuresdir, paste(mystk, "mat_ogive.jpg", sep="_")), 
         device="jpeg", 
         width=10, height=8, units="in")


  # COND = OK
  
  # ----------------------------------------------------------------------------
  # 3. Density dependence
  # ----------------------------------------------------------------------------
  
  # Get theoretical values from literature -------------------------------------
  par=FLife::teleost[c("linf","k","t0","l50","a","b"),dimnames(teleost)$iter==mylatin]
  
  # Uses life history theory to derive parameters for biological relationships
  par=FLife::lhPar(par,
                   m1=myparams[myparams$stock==mystk,"m1"],
                   m2=myparams[myparams$stock==mystk,"m2"],
                   m3=myparams[myparams$stock==mystk,"m3"])
  
  # WHAT IS THE ASSUMED SELECTIVITY AND WHERE DOES IT COME FROM?
  
  # derive a FLBRP object with Lorenzen natural mortality
  eq =FLife::lhEql(par,m=myparams[myparams$stock==mystk,"m"])
  
  
  # plot the BRP object --------------------------------------------------------
  p <-ggplotFL::plot(eq) + 
           theme_bw() +
           labs(title=paste("FLBRP object"))
  
  ggsave(p, 
         filename=file.path(figuresdir, paste(mystk, "brp.jpg", sep="_")), 
         device="jpeg", 
         width=10, height=8, units="in")
  
  # convert to FLStock and project
  om=as(eq,"FLStock")
  om=fwd(om,fbar=fbar(om)[,-1],sr=eq)

  # update the par object
  par=rbind(par,
            FLPar(c(bref =c(refpts(eq)["msy","biomass"]),
                    delta= myparams[myparams$stock==mystk,"delta"],
                    matk = myparams[myparams$stock==mystk,"matk"],
                    w50  =par["a"]*par["l50"]^par["b"],
                    alpha=myparams[myparams$stock==mystk,"alpha"])))

  fileConn <-file(file.path(tablesdir, paste(mystk, "par.txt", sep="_")))
  as.data.frame(par) %>% dplyr::select(-iter) %>% pander::pandoc.table(style="simple") %>% capture.output() %>% writeLines(., con=fileConn)
  close(fileConn)

  # save par file
  save(par,file=file.path(dropboxdir,paste0("data/om/par",mystk,".RData")))
  
  # Simulation of density dependence in mass-at-age ----------------------------
  # for different levels of biomass relative to Bmsy
  
  dat=ddWt(stock.wt(eq),biomass(eq),refpts(eq)["msy","biomass"]) 
  dat=merge(as.data.frame(dat,drop=T),
            as.data.frame(biomass(eq)%/%refpts(eq)["msy","biomass"],drop=T),
            by="year")
  names(dat)[3:4]=c("mass","biomass")
  
  p <- ggplot(dat)+    
           theme_bw() +
           geom_line(aes(age,mass,col=biomass,group=year))+
           scale_x_continuous(limits=c(0,maxage), breaks=(seq(0,maxage,1)))+
           # scale_y_continuous(limits=c(0,800))+
           geom_line(aes(age,data),data=as.data.frame(stock.wt(eq)),col="red",size=2)+
           labs(x="Age", y="Mass-at-age", colour="Biomass/Bmsy") 
  
  ggsave(p, 
         filename=file.path(figuresdir, paste(mystk, "sim_dd_mass.jpg", sep="_")), 
         device="jpeg", 
         width=10, height=10, units="in")
  
  # Simulation of density dependence in maturity-at-age ------------------------ 
  # for different levels of biomass relative to Bmsy
  
  dat=transform(dat,mat=1/(1+exp(-par["matk"]*(mass-par["w50"]))))
  
  p <- 
    ggplot(dat)+           
    theme_bw() +
    geom_line(aes(age,mat,col=biomass,group=year))+
    scale_x_continuous(limits=c(0,6), breaks=(seq(0,maxage,1)))+
    geom_line(aes(age,data),data=as.data.frame(mat(eq)),col="red",size=2)+
    labs(x="Age", y="Maturity-at-age", colour="Biomass/Bmsy")
  
  ggsave(p, 
         filename=file.path(figuresdir, paste(mystk, "sim_dd_maturity.jpg", sep="_")), 
         device="jpeg", 
         width=10, height=10, units="in")
  

  # Simulation of density dependence in M-at-age -------------------------------
  # for different levels of biomass relative to Bmsy
  
  dat=transform(dat,m=c(par["m1"])*(mass^c(par["m2"])))

  p <-
    ggplot(dat)+           
    geom_line(aes(age,m,col=biomass,group=year))+
    scale_x_continuous(limits=c(0,maxage), breaks=(seq(0,maxage,1)))+
    geom_line(aes(age,data),data=as.data.frame(m(eq)),col="red",size=2)+
    labs(x="Age", y="M-at-age", colour="Biomass/Bmsy")
  
  ggsave(p, 
       filename=file.path(figuresdir, paste(mystk, "sim_dd_m.jpg", sep="_")), 
       device="jpeg", 
       width=10, height=10, units="in")

  # Comparison of equilibrium yields curves without DD and with DD mass --------
  
  x         = propagate(eq,length(dimnames(eq)$year))
  fbar(x)   = fbar(x)[,1]
  fbar(x)[] = c(fbar(eq))*2
  
  # landing weight at age and discard weight at age relative to stock weight
  lsr       = landings.wt(eq)/stock.wt(eq)
  dsr       = discards.wt(eq)/stock.wt(eq)
  
  # iterate
  for (i in seq(10)){
    stock.wt(   x)=ddWt(stock.wt(eq),
                        biomass(x),
                        refpts(eq)["msy","biomass"],
                        alpha=-0.2)
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
    stock.wt(   x)=ddWt(stock.wt(eq),
                        biomass(x),
                        refpts(eq)["msy","biomass"],
                        alpha=-0.2)
    landings.wt(x)=stock.wt(x)*lsr
    discards.wt(x)=stock.wt(x)*dsr
    
    mat(x) =1/(1+exp(-par["matk"]*(stock.wt(x)-par["w50"])))
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
    
    mat(x) = 1/(1+exp(-par["matk"]*(stock.wt(x)-par["w50"])))
    m(x)   = par["m1"]%*%(stock.wt(x)%^%par["m2"])
    x      = brp(x)
  }
  
  df <- 
    bind_rows(
      df,
      model.frame(FLQuants(x, biomass=function(x) biomass(x), catch=function(x) catch(x)),drop=T) %>% mutate(scen="DD mass+mat+M")
    )
  
  p <-
    ggplot(data=df, aes(biomass, catch))+
    theme_bw() +
    geom_line(aes(colour=scen)) +
    labs(x="Total biomass", y="Catch", colour="") +
    scale_colour_manual(values =c("no-DD"         = "black",
                                  "DD mass"       = "red", 
                                  "DD mass+mat"   = "blue",
                                  "DD mass+mat+M" = "darkgreen")) 
  
  
  ggsave(p, 
         filename=file.path(figuresdir, paste(mystk, "eq_yield.jpg", sep="_")), 
         device="jpeg", 
         width=10, height=10, units="in")

  # Check 1 --------------------------------------------------------------------
  
  # generate FLquants with different F values
  f = FLQuant(c(rep(seq(0, 1,length.out=51)    ,each=101),
                rep(seq(1,25,length.out=61)[-1],each=101)),
              dimnames=list(year=dimnames(om)$year,
                            iter=seq(101)))  %*%
    refpts(eq)["msy","harvest"]
  
  # Add F's to eq object 
  fbar(eq)=f
  eq=brp(eq)
  
  om=as(eq,"FLStock")
  om=fwd(om,f=f[,-1],sr=eq)
  
  # simulate with DD for 100 years into the future (THIS TAKES ABOUT 10 MINUTES)
  om2=om
  
  ptm <- proc.time()
  
  for (year in dimnames(om2)$year[-1]){
    om2=ddFn(year,om2,par)
    om2=fwd(om2,f=f[,year],sr=eq)
  }
  catch(om2)=computeCatch(om2,slot="all")
  
  # Stop the clock
  proc.time() - ptm
  
  # plot(window(FLStocks("OM"=om,DD=om2),start=40))  # DIFFICULT TO INTERPRET

  t <-
    bind_rows(
      model.frame(FLQuants(om2[,dim(m(om2))[2]],"biomass"=biomass,"catch"=catch)) %>% mutate(scen="DD sim"),
      model.frame(FLQuants(om[,dim(m(om))[2]],"biomass"=biomass,"catch"=catch)) %>% mutate(scen="no-DD sim"),
      filter(df, scen=="DD mass+mat+M") %>% mutate(scen="DD equilibrium")
    ) 
  
  p <-
    ggplot(data=t) +  
    geom_line(aes(biomass,catch, col=scen))+ 
    labs(title="Check on equilibrium and simulation results") +
    scale_colour_manual(values =c("no-DD sim"         = "black",
                                  "DD sim"       = "red", 
                                  "DD equilibrium"   = "blue")) 
  
  ggsave(p, 
         filename=file.path(figuresdir, paste(mystk, "check_sim.jpg", sep="_")), 
         device="jpeg", 
         width=10, height=10, units="in")
  
  
  # ------------------------------------------------------------------------------
  # 4. OM
  # ------------------------------------------------------------------------------

  # load(file.path(dropboxdir,"data/inputs/ices.RData"))
  # ices=ices[["mac.27.nea"]]
  # ices=window(ices,start=1991)
  
  ices=window(stk,start=1991)
  name(ices) <- mystkname

  # Start with the SAM assessment in forward projection
  sr  = as.FLSR(ices,model=bevholtSV)
  sr  = fmle(sr,fixed=list(s=0.8,spr0=mean(spr0(ices))),control=list(silent=TRUE))
  eq  = FLBRP(ices,ab(sr))
  sam = fwdWindow(ices,eq,end=2050) 
  
  F  =propagate(window(fbar(sam),start=2020),101)
  F[]=rep(c(seq(0,                               c(refpts(eq)["msy",  "harvest"]),length.out=51),
            seq(c(refpts(eq)["msy",  "harvest"]),c(refpts(eq)["crash","harvest"])*1.2,length.out=51)[-1]),each=dim(F)[2])
  
  control=as(FLQuants("f"=F),"fwdControl")
  sam   =fwd(sam,control=control,sr=eq)
  
  df_om  =model.frame(FLQuants(sam[,"2050"], "Biomass"=function(x) biomass(x),
                                    "SSB"=function(x) ssb(x),
                                    "F"=function(x) fbar(x),
                                    "Yield"=function(x) catch(x)),drop=T) %>% mutate(scen="sam")
                

  # VPA assessment based on SAM assessment -------------------------------------
  
  vpa =ices+FLAssess:::VPA(ices)
  
  sr  =as.FLSR(vpa,model=bevholtSV)
  sr  =fmle(sr,fixed=list(s=0.8,spr0=mean(spr0(vpa))),control=list(silent=TRUE))
  eq  =FLBRP(vpa,ab(sr))
  vpa =fwdWindow(vpa,eq,end=2050) 
  
  F  =propagate(window(fbar(vpa),start=2020),101)
  F[]=rep(c(seq(0,                               c(refpts(eq)["msy",  "harvest"]),length.out=51),
            seq(c(refpts(eq)["msy",  "harvest"]),c(refpts(eq)["crash","harvest"])*1.2,length.out=51)[-1]),each=dim(F)[2])
  
  control=as(FLQuants("f"=F),"fwdControl")
  vpa    =fwd(vpa,control=control,sr=eq)
  
  df_om  =bind_rows(
    df_om,
    model.frame(FLQuants(vpa[,"2050"], "Biomass"=function(x) biomass(x),
                                       "SSB"=function(x) ssb(x),
                                       "F"=function(x) fbar(x),
                                       "Yield"=function(x) catch(x)),drop=T) %>% mutate(scen="vpa"))

  # VPA assessment with age varying M -----------------------------------------------------

  load(file.path(dropboxdir, paste0("data/om/par",mystk,".RData")))
  par["m1"]  =par["m1"]/10      # This is a bit awkward; every time you rerun it will become smaller
  par["w50"] =par["w50"]/1000
  par["matk"]=30                # WHAT DOES MATK DO ?? 
  
  fileConn <-file(file.path(tablesdir, paste(mystk, "par_om.txt", sep="_")))
  as.data.frame(par) %>% dplyr::select(-iter) %>% pander::pandoc.table(style="simple") %>% capture.output() %>% writeLines(., con=fileConn)
  close(fileConn)
  
  # load(file.path(dropboxdir,"data/inputs/ices.RData"))
  # vpaM=ices[["mac.27.nea"]]
  # vpaM=window(vpaM,start=1991)
  vpaM   = ices
  m(vpaM)=par["m1"]%*%(stock.wt(vpaM)%^%par["m2"])
  vpaM=vpaM+FLAssess:::VPA(vpaM)
  
  sr  =as.FLSR(vpaM,model=bevholtSV)
  sr  =fmle(sr,fixed=list(s=0.8,spr0=mean(spr0(vpa))),control=list(silent=TRUE))
  eq  =FLBRP(vpaM,ab(sr))
  vpaM=fwdWindow(vpaM,eq,end=2050) 
  
  F  =propagate(window(fbar(vpaM),start=2020),101)
  F[]=rep(c(seq(0,                               c(refpts(eq)["msy",  "harvest"]),length.out=51),
            seq(c(refpts(eq)["msy",  "harvest"]),c(refpts(eq)["crash","harvest"])*1.5,length.out=51)[-1]),each=dim(F)[2])
  
  control=as(FLQuants("f"=F),"fwdControl")
  vpaM  =fwd(vpaM,control=control,sr=eq)
  
  d1=model.frame(FLQuants(vpaM[,"2050"],"Biomass"=function(x) biomass(x),"SSB"=function(x) ssb(x),
                          "F"=function(x) fbar(x),"Yield"=function(x) catch(x)),drop=T)
  
  df_om  =bind_rows(
    df_om,
    d1 %>% mutate(scen="vpa-M"))
  
  df_helper <-
    df_om %>% 
    group_by(scen) %>% 
    summarise(Yield = max(Yield, na.rm=TRUE)) %>% 
    left_join(df_om)
  
  # Plot of yield vs biomass
  p <-
    df_om %>% 
    ggplot(aes(x=Biomass, y=Yield)) +
    theme_bw() +
    geom_line(aes(colour=scen)) +
    geom_segment(data=df_helper,
                 aes(x=Biomass, xend=Biomass, y=0, yend=Yield, colour=scen)) +
    geom_segment(data=df_helper,
                 aes(x=0, xend=Biomass, y=Yield, yend=Yield, colour=scen))
  
  ggsave(p,
         filename=file.path(figuresdir, paste(mystk, "om_biomass_yield.jpg", sep="_")),
         device="jpeg",
         width=10, height=10, units="in")
  
  # plot of F vs biomass
  p <-
    df_om %>% 
    ggplot(aes(x=Biomass, y=F)) +
    theme_bw() +
    geom_line(aes(colour=scen)) +
    geom_segment(data=df_helper,
                 aes(x=0, xend=Biomass, y=F, yend=F, colour=scen)) +
    geom_segment(data=df_helper,
                 aes(x=Biomass, xend=Biomass, y=0, yend=F, colour=scen))
  
  ggsave(p,
         filename=file.path(figuresdir, paste(mystk, "om_biomass_F.jpg", sep="_")),
         device="jpeg",
         width=10, height=10, units="in")
  
  
  ##############################################################################
  # Setting Bref !!
  
  par["bref"]=subset(d1, Yield==max(Yield))[,"Biomass"]
  
  ##############################################################################

  # OM with DD in mass ---------------------------------------------------------
  
  # Start the clock
  ptm <- proc.time()
  
  vpaDDM=vpaM
  
  for (i in ac(2020:2050)) {
    vpaDDM =ddFn(i,vpaDDM,par,massFlag=TRUE,matFlag=FALSE,mFlag=FALSE)
    control=as(FLQuants("f"=F[,i]),"fwdControl")
    vpaDDM =fwd(vpaDDM,control=control,sr=eq)}
  
  p1=ggplot(stock.wt(vpaDDM[,"2030"]))+geom_line(aes(age,data,group=iter))+
    geom_line(aes(age,data),data=as.data.frame(stock.wt(vpaM[,"2031"])),col="red")+xlab("Age")+ylab("Mass-at-age")
  p2=ggplot(mat(vpaDDM[,"2030"]))+geom_line(aes(age,data,group=iter))+
    geom_line(aes(age,data),data=as.data.frame(mat(vpaM[,"2031"])),col="red")+xlab("Age")+ylab("Mat-at-age")
  p3=ggplot(m(vpaDDM[,"2030"]))+geom_line(aes(age,data,group=iter))+
    geom_line(aes(age,data),data=as.data.frame(m(vpaM[,"2031"])),col="red")+xlab("Age")+ylab("M-at-age")
  
  d2=model.frame(FLQuants(vpaDDM[,"2050"],"Biomass"=function(x) biomass(x),"SSB"=function(x) ssb(x),
                          "F"=function(x) fbar(x),"Yield"=function(x) catch(x)),drop=T)
  
  df_om  =bind_rows(
    df_om,
    d2 %>% mutate(scen="vpa-DDM"))
  
  p4=ggplot(d2)+
    geom_line(aes(SSB,Yield))
  
  p <- ggarrange(p1,p2,p3,p4, ncol=2, nrow=2)

  ggsave(p,
         filename=file.path(figuresdir, paste(mystk, "VPADDM_4panels.jpg", sep="_")),
         device="jpeg",
         width=10, height=10, units="in")
  
  # Stop the clock
  proc.time() - ptm
  
  # OM with DD in mass and maturity --------------------------------------------
  
  # Start the clock
  ptm <- proc.time()
  
  vpaDDMM =vpaM
  
  for (i in ac(2020:2050)) {
    vpaDDMM =ddFn(i,vpaDDMM,par,massFlag=TRUE,matFlag=TRUE,mFlag=FALSE)
    control=as(FLQuants("f"=F[,i]),"fwdControl")
    vpaDDMM =fwd(vpaDDMM,control=control,sr=eq)}
  
  p1=ggplot(stock.wt(vpaDDMM[,"2030"]))+geom_line(aes(age,data,group=iter))+
    geom_line(aes(age,data),data=as.data.frame(stock.wt(vpaM[,"2031"])),col="red")+xlab("Age")+ylab("Mass-at-age")
  p2=ggplot(mat(vpaDDMM[,"2030"]))+geom_line(aes(age,data,group=iter))+
    geom_line(aes(age,data),data=as.data.frame(mat(vpaM[,"2031"])),col="red")+xlab("Age")+ylab("Mat-at-age")
  p3=ggplot(m(vpaDDMM[,"2030"]))+geom_line(aes(age,data,group=iter))+
    geom_line(aes(age,data),data=as.data.frame(m(vpaM[,"2031"])),col="red")+xlab("Age")+ylab("M-at-age")
  
  d3=model.frame(FLQuants(vpaDDMM[,"2050"],"Biomass"=function(x) biomass(x),"SSB"=function(x) ssb(x),
                          "F"=function(x) fbar(x),"Yield"=function(x) catch(x)),drop=T)
  df_om  =bind_rows(
    df_om,
    d3 %>% mutate(scen="vpa-DDMM"))
  
  p4=ggplot(d3)+
    geom_line(aes(SSB,Yield))+
    geom_line(aes(SSB,Yield),col="red",data=d1)
  
  p <- ggarrange(p1,p2,p3,p4, ncol=2, nrow=2)
  
  ggsave(p,
         filename=file.path(figuresdir, paste(mystk, "VPADDMM_4panels.jpg", sep="_")),
         device="jpeg",
         width=10, height=10, units="in")
  
  # Stop the clock
  proc.time() - ptm
  
  # OM with DD in mass and maturity and M --------------------------------------
  
  # Start the clock
  ptm <- proc.time()
  
  vpaDDMMM=vpaM
  
  for (i in ac(2020:2050)) {
    vpaDDMMM =ddFn(i,vpaDDMMM,par,massFlag=TRUE,matFlag=TRUE,mFlag=TRUE)
    control=as(FLQuants("f"=F[,i]),"fwdControl")
    vpaDDMMM =fwd(vpaDDMMM,control=control,sr=eq)}
  
  p1=ggplot(stock.wt(vpaDDMMM[,"2030"]))+geom_line(aes(age,data,group=iter))+
    geom_line(aes(age,data),data=as.data.frame(stock.wt(vpaM[,"2031"])),col="red")+xlab("Age")+ylab("Mass-at-age")
  p2=ggplot(mat(vpaDDMMM[,"2030"]))+geom_line(aes(age,data,group=iter))+
    geom_line(aes(age,data),data=as.data.frame(mat(vpaM[,"2031"])),col="red")+xlab("Age")+ylab("Mat-at-age")
  p3=ggplot(m(vpaDDMMM[,"2030"]))+geom_line(aes(age,data,group=iter))+
    geom_line(aes(age,data),data=as.data.frame(m(vpaM[,"2031"])),col="red")+xlab("Age")+ylab("M-at-age")
  
  d4=model.frame(FLQuants(vpaDDMMM[,"2050"],"Biomass"=function(x) biomass(x),"SSB"=function(x) ssb(x),
                          "F"=function(x) fbar(x),"Yield"=function(x) catch(x)),drop=T)
  df_om  =bind_rows(
    df_om,
    d4 %>% mutate(scen="vpa-DDMMM"))
  
  p4=ggplot(d4)+
    geom_line(aes(SSB,Yield))+
    geom_line(aes(SSB,Yield),col="red",data=d1)
  
  p <- ggarrange(p1,p2,p3,p4, ncol=2, nrow=2)
  
  ggsave(p,
         filename=file.path(figuresdir, paste(mystk, "VPADDMMM_4panels.jpg", sep="_")),
         device="jpeg",
         width=10, height=10, units="in")
  
  # Stop the clock
  proc.time() - ptm
  

  d=rbind(cbind("What"="M",d1),
          cbind("What"="DD-M",d2),
          cbind("What"="DD-MM",d3),
          cbind("What"="DD-MMM",d4))
  p <- ggplot(d)+
    geom_line(aes(SSB,Yield,col=What))

  ggsave(p,
         filename=file.path(figuresdir, paste(mystk, "comparing_DD.jpg", sep="_")),
         device="jpeg",
         width=10, height=10, units="in")
  
  # p <- ggplot(df_om)+
  #   geom_line(aes(SSB,Yield,col=scen))
  
  # ------------------------------------------------------------------------------
  # 5. MSE
  # ------------------------------------------------------------------------------
  
  # Forward projections with F=fmsy (fwd), 0.8xfmsy?, 1.0xfmsy, 1.2xfmsy?
  
  # sr=fmle(as.FLSR(om,model="bevholt"),control=list(silent=TRUE))
  # eq=FLBRP(om,sr=sr)
  # 
  # om=fwdWindow(om,end=maxyear,eq)
  # f =fbar(om)[,ac(2021:maxyear)]%=%refpts(eq)["msy","harvest"]
  # om=fwd(om,fbar=f,sr=eq)
  # 
  # om2=om
  # om3=om
  # om4=om
  # 
  # for (iYear in ac(2021:maxyear)){
  #   om2=ddFn(iYear,om2,par)
  #   om2=fwd(om2,fbar=f[,iYear]*0.8,sr=eq)
  #   
  #   om3=ddFn(iYear,om3,par)
  #   om3=fwd(om3,fbar=f[,iYear],sr=eq)
  #   
  #   om4=ddFn(iYear,om4,par)
  #   om4=fwd(om4,fbar=f[,iYear]*1.2,sr=eq)
  # }
  # 
  # plot(FLStocks("fwd"=om,"0.8"=om2,"1.0"=om3,"1.2"=om4))
  
  
  
# } # end of stk loop


