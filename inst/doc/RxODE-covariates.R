## ---- echo=FALSE---------------------------------------------------------
options(cli.unicode=FALSE, crayon.enabled=FALSE);
options(knitr.table.format = "html")
htmltools::img(src = knitr::image_uri("logo.png"), 
               alt = 'RxODE', 
               style = 'position:absolute; top:0; right:0; padding:10px; border: 0;')

## ------------------------------------------------------------------------
library(RxODE)
mod3 <- RxODE({
    KA=2.94E-01;
    CL=1.86E+01; 
    V2=4.02E+01; 
    Q=1.05E+01;
    V3=2.97E+02; 
    Kin0=1;
    Kout=1;
    EC50=200;
    ## The linCmt() picks up the variables from above
    C2   = linCmt();
    Tz= 8
    amp=0.1
    eff(0) = 1  ## This specifies that the effect compartment starts at 1.
    ## Kin changes based on time of day (like cortosol)
    Kin =   Kin0 +amp *cos(2*pi*(ctime-Tz)/24)
    d/dt(eff) =  Kin - Kout*(1-C2/(EC50+C2))*eff;
})


ev <- eventTable(amount.units="mg", time.units="hours") %>%
    add.dosing(dose=10000, nbr.doses=1, dosing.to=2) %>%
    add.sampling(seq(0,48,length.out=100));


 ## Create data frame of  8 am dosing for the first dose
cov.df  <- data.frame(ctime =(seq(0,48,length.out=100)+8) %% 24);


## ----results="asis"------------------------------------------------------

r1 <- solve(mod3, ev, covs=cov.df,covs_interpolation="linear")
rxHtml(r1)

## ----out.width="100%"----------------------------------------------------
par(mfrow=c(1,2))
matplot(r1[,"C2"], type="l", ylab="Central Concentration")
matplot(r1[,"eff"], type="l", ylab = "Effect")

## ----results="asis"------------------------------------------------------
r2 <- solve(mod3, ev, covs=cov.df,covs_interpolation="constant")
rxHtml(r2)

## ----out.width="100%"----------------------------------------------------
par(mfrow=c(1,2))
matplot(r2[,"C2"], type="l", ylab="Central Concentration")
matplot(r2[,"eff"], type="l", ylab = "Effect")

