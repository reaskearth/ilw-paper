rm(list=ls())

# LIBRARIES
library(data.table)
library(stats)
library(stringr)
library(quadprog)
library(ggplot2)
library(ggpattern)
library(Matrix)
library(fields)
library(arrow)

# PARAMETERS
catalog="STD"
baseline="1951-2020"
forecast="JUNE"
years=seq(1985,2024)
yini=years[1]
yend=years[length(years)]

# IO PATHS
path_project="./" # set project path if needed
path_verisk=paste0(path_project,"Data/Restricted_Data/ylt_Verisk_",catalog,"_cat.csv")
path_reask=paste0(path_project,"CBRA/Output/","STD","_",baseline,"_",forecast,"/")
path_pcs=paste0(path_project,"Data/Restricted_Data/PCS_RegionalSplit.csv")
path_landfall=paste0(path_project,"Data/Landfall.csv")
path_index=paste0(path_project,"Data/IndexationFactors.csv")
path_verisk_split=paste0(path_project,"Data/Restricted_Data/Verisk_",catalog,"_RegionalSplit.csv")
path_ilw=paste0(path_project,"Data/ILW_prices.csv")
path_climate=paste0(path_project,"Data/Reask_ClimateIndices.parquet")
path_out=paste0(path_project,"Output/Figures_",yini,"-",yend,"/",catalog,"/",forecast,"/")
path_save=paste0(path_project,"Output/",baseline,"_",yini,"-",yend,"_",catalog,"_",forecast,".RData")
if (!dir.exists(path_out)) dir.create(path_out,recursive=T)

# FUNCTIONS

read_bulk_ylt=function(year,path) {
  file_ylt=paste0(path,year,"/ylt-adjusted.csv")
  ylt=fread(file_ylt)
  ylt[,Year:=year]
  return(ylt)
} 

calc_percentile=function(vals,smp) {
  cdf=ecdf(smp)
  percs=cdf(vals)
  return(percs)
}

calc_tvar=function(losses,allocations,prob) {
  losses=t(allocations*t(losses))
  loss_total=rowSums(losses)
  ox=order(loss_total,decreasing=T)
  loss_total=sort(loss_total,decreasing=T)
  freq=seq(1,nrow(losses))/nrow(losses)
  ix=which(freq<=prob)
  ox=ox[ix]
  tvar_contribs=colMeans(losses[ox,])
  return(tvar_contribs)
}

calc_asset_returns=function(notional,prices,years,pcs) {
  
  ny=length(years)
  ns=nrow(prices)

  ILW=data.table(Strike=rep(prices$Strike,ny),
                 Zone=rep(prices$Zone,ny),
                 Year=rep(years,each=ns),
                 isStrike=FALSE,
                 Return_Realized=as.numeric(NA))
  
  for (y in years) {
    
    for (is in 1:ns){
      z=prices[is,Zone]
      s=prices[is,Strike]
      rol=prices[is,Price]
      
      if (z=="Nationwide") {
        pcs_loss=as.numeric(any(pcs[Year==y]$Loss>s))
      } else if (z=="Northeast") {
        pcs_loss=as.numeric(any(pcs[Year==y]$LossNe>s))
      } else if (z=="Southeast") {
        pcs_loss=as.numeric(any(pcs[Year==y]$LossSe>s))
      } else if (z=="Florida") {
        pcs_loss=as.numeric(any(pcs[Year==y]$LossFl>s))
      } else if (z=="Gulf") {
        pcs_loss=as.numeric(any(pcs[Year==y]$LossGl>s))
      } else {stop("Zone specified incorrect")}
    
      ILW[Year==y & Zone==z & Strike==s,isStrike:=pcs_loss>0]
      ILW[Year==y & Zone==z & Strike==s,Return_Realized:=notional*(rol-pcs_loss)]
    }
  }
  return(ILW)
}

optimize_portfolio=function(aversion,notional,prices,years,verisk,reask,pcs,reinvest=F) {

  ns=length(prices$Strike)
  ny=length(years)
  
  ILW_PORTF=data.table(Year=rep(years,each=ns),
                       Strike=rep(prices$Strike,ny),
                       Zone=rep(prices$Zone,ny),
                       Aversion=aversion,
                       Allocation_Verisk=as.numeric(NA),
                       Allocation_Reask=as.numeric(NA),
                       Return_Expected_Verisk=as.numeric(NA),
                       Return_Expected_Reask=as.numeric(NA),
                       Return_Realized_Verisk=as.numeric(NA),
                       Return_Realized_Reask=as.numeric(NA))
  
  notional_reask=notional
  notional_verisk=notional
  
  for (y in years) {
    reask_y=reask[Year==y]
    reask_loss=array(0.0,dim=c(10000,ns))
    verisk_loss=array(0.0,dim=c(10000,ns))
    rol=rep(NA,ns)
    pcs_loss=rep(0,ns)
    
    for (is in 1:ns) {
      s=prices[is,Strike]
      z=prices[is,Zone]
      if (z=="Nationwide") {
        t=reask_y[,.(Loss=any(Loss>s)),by="Sample"]
        reask_loss[t$Sample,is]=as.numeric(t$Loss)
        t=verisk[,.(Loss=any(Loss>s)),by="Sample"]
        verisk_loss[t$Sample,is]=as.numeric(t$Loss)
        pcs_loss[is]=as.numeric(any(pcs[Year==y]$Loss>s))
      } else if (z=="Northeast") {
        t=reask_y[,.(Loss=any(LossNe>s)),by="Sample"]
        reask_loss[t$Sample,is]=as.numeric(t$Loss)
        t=verisk[,.(Loss=any(LossNe>s)),by="Sample"]
        verisk_loss[t$Sample,is]=as.numeric(t$Loss)
        pcs_loss[is]=as.numeric(any(pcs[Year==y]$LossNe>s))
      } else if (z=="Southeast") {
        t=reask_y[,.(Loss=any(LossSe>s)),by="Sample"]
        reask_loss[t$Sample,is]=as.numeric(t$Loss)
        t=verisk[,.(Loss=any(LossSe>s)),by="Sample"]
        verisk_loss[t$Sample,is]=as.numeric(t$Loss)
        pcs_loss[is]=as.numeric(any(pcs[Year==y]$LossSe>s))
      } else if (z=="Florida") {
        t=reask_y[,.(Loss=any(LossFl>s)),by="Sample"]
        reask_loss[t$Sample,is]=as.numeric(t$Loss)
        t=verisk[,.(Loss=any(LossFl>s)),by="Sample"]
        verisk_loss[t$Sample,is]=as.numeric(t$Loss)
        pcs_loss[is]=as.numeric(any(pcs[Year==y]$LossFl>s))
      } else if (z=="Gulf") {
        t=reask_y[,.(Loss=any(LossGl>s)),by="Sample"]
        reask_loss[t$Sample,is]=as.numeric(t$Loss)
        t=verisk[,.(Loss=any(LossGl>s)),by="Sample"]
        verisk_loss[t$Sample,is]=as.numeric(t$Loss)
        pcs_loss[is]=as.numeric(any(pcs[Year==y]$LossGl>s))
      } else {stop("Zone specified incorrect")}
      rol[is]=prices[is,Price]
    }
    
    reask_cov=cov(reask_loss)
    verisk_cov=cov(verisk_loss)
    
    #reask_cov=nearPD(reask_cov,ensureSymmetry=T,maxit=1000)
    #verisk_cov=nearPD(verisk_cov,ensureSymmetry=T,maxit=1000)
    
    reask_er=rol-colMeans(reask_loss)
    verisk_er=rol-colMeans(verisk_loss)
    
    # NO SHORT POSITIONS, NO LEVERAGE (0<w<1, sum(w)=1)
    A=cbind(rep(1,ns), diag(1,ns,ns), diag(-1,ns,ns))
    b=c(1, rep(0,ns), rep(-1,ns))
    meq=1
    
    verisk_w=round(solve.QP(Dmat=aversion*verisk_cov, dvec=verisk_er, Amat=A, bvec=b, meq=meq)$solution,ns)
    reask_w=round(solve.QP(Dmat=aversion*reask_cov, dvec=reask_er, Amat=A, bvec=b, meq=meq)$solution,ns)
    
    verisk_alloc=verisk_w*notional_verisk
    reask_alloc=reask_w*notional_reask
    
    verisk_ret_real=verisk_alloc*(rol-pcs_loss)
    reask_ret_real=reask_alloc*(rol-pcs_loss)
    
    verisk_ret_exp=verisk_alloc*verisk_er
    reask_ret_exp=reask_alloc*reask_er
    
    verisk_var=verisk_alloc * verisk_cov%*%verisk_alloc 
    reask_var=reask_alloc * reask_cov%*%reask_alloc 
    
    verisk_tvar=calc_tvar(verisk_loss,verisk_alloc,0.1)
    reask_tvar=calc_tvar(reask_loss,reask_alloc,0.1)
      
    ILW_PORTF[Year==y,c("Allocation_Verisk",
                        "Allocation_Reask",
                        "Return_Expected_Verisk",
                        "Return_Expected_Reask",
                        "Variance_Verisk",
                        "Variance_Reask",
                        "TVaR10_Verisk",
                        "TVaR10_Reask",
                        "Return_Realized_Verisk",
                        "Return_Realized_Reask"):=list(verisk_alloc,reask_alloc,verisk_ret_exp,reask_ret_exp,verisk_var,reask_var,verisk_tvar,reask_tvar,verisk_ret_real,reask_ret_real)]
    
    if (reinvest==T) {
      notional_reask=notional_reask+sum(reask_ret_real)
      notional_verisk=notional_verisk+sum(verisk_ret_real)
    }

  }
  return(ILW_PORTF)
}

# MAIN

# READ VERISK YLT, APPLY REGIONAL SPLIT, AGGREGATE ANNUALLY, FILL-IN MISSING YEARS, CALCULATE RETURN PERIODS 
verisk=fread(path_verisk)
verisk_split=fread(path_verisk_split)[,.(event_id=EventId,PctSe,PctNe,PctFl,PctGl)]
verisk_occ=verisk_split[verisk,on=c("event_id")]
verisk_occ=verisk_occ[,.(Sample=year,
                         EventId=event_id,
                         Loss=loss/1e9,
                         LossSe=loss*PctSe/1e9,
                         LossNe=loss*PctNe/1e9,
                         LossFl=loss*PctFl/1e9,
                         LossGl=loss*PctGl/1e9)]
verisk_agg=verisk_occ[,.(Freq=.N,Loss=sum(Loss),
                         LossSe=sum(LossSe),
                         LossNe=sum(LossNe),
                         LossFl=sum(LossFl),
                         LossGl=sum(LossGl)),by="Sample"]
verisk_agg=verisk_agg[data.table(Sample=seq(1,10000)),on="Sample"]
verisk_agg[is.na(Freq),Freq:=0]
verisk_agg[is.na(Loss),Loss:=0]
verisk_agg[is.na(LossSe),LossSe:=0]
verisk_agg[is.na(LossNe),LossNe:=0]
verisk_agg[is.na(LossFl),LossFl:=0]
verisk_agg[is.na(LossGl),LossGl:=0]
setorder(verisk_agg,-LossSe)
verisk_agg[,RpSe:=10000/.I]
setorder(verisk_agg,-LossNe)
verisk_agg[,RpNe:=10000/.I]
setorder(verisk_agg,-LossGl)
verisk_agg[,RpGl:=10000/.I]
setorder(verisk_agg,-LossFl)
verisk_agg[,RpFl:=10000/.I]
setorder(verisk_agg,-Loss)
verisk_agg[,Rp:=10000/.I]

# CALCULATE STATISTICS ON VERISK
verisk_stat=verisk_agg[,.(q1Freq=quantile(Freq,0.1),q25Freq=quantile(Freq,0.25),q5Freq=quantile(Freq,0.5),q75Freq=quantile(Freq,0.75),q9Freq=quantile(Freq,0.9),
                      q1Loss=quantile(Loss,0.1),q25Loss=quantile(Loss,0.25),q5Loss=quantile(Loss,0.5),q75Loss=quantile(Loss,0.75),q9Loss=quantile(Loss,0.9),
                      q1LossSe=quantile(LossSe,0.1),q25LossSe=quantile(LossSe,0.25),q5LossSe=quantile(LossSe,0.5),q75LossSe=quantile(LossSe,0.75),q9LossSe=quantile(LossSe,0.9),
                      q1LossNe=quantile(LossNe,0.1),q25LossNe=quantile(LossNe,0.25),q5LossNe=quantile(LossNe,0.5),q75LossNe=quantile(LossNe,0.75),q9LossNe=quantile(LossNe,0.9),
                      meanLoss=sum(Loss)/10000,meanFreq=sum(Freq)/10000,meanLossSe=sum(LossSe)/10000, meanLossNe=sum(LossNe)/10000,
                      stdLoss=sd(Loss),stdFreq=sd(Freq),stdLossSe=sd(LossSe),stdLossNe=sd(LossNe))]

# READ REASK YLT FOR ALL HISTORICAL YEARS, APPLY REGIONAL SPLIT, AGGREGATE ANNUALLY, FILL-IN MISSING YEARS, CALCULATE RETURN PERIODS 
reask=rbindlist(lapply(years,read_bulk_ylt,path_reask))
reask[,event_id:=as.integer(str_replace_all(event_id,"reask_",""))]
reask_occ=verisk_split[reask,on="event_id"]
reask_occ=reask_occ[,.(Year,
                       Sample=sample,
                       EventId=event_id,
                       Loss=loss/1e9,
                       LossSe=loss*PctSe/1e9,
                       LossNe=loss*PctNe/1e9,
                       LossFl=loss*PctFl/1e9,
                       LossGl=loss*PctGl/1e9)]
reask_agg=reask_occ[,.(Freq=.N,Loss=sum(Loss),
                       LossSe=sum(LossSe),
                       LossNe=sum(LossNe),
                       LossFl=sum(LossFl),
                       LossGl=sum(LossGl)), by=c("Year","Sample")]
           
reask_agg=reask_agg[data.table(Year=rep(years,each=10000),Sample=rep(1:10000,length(years))),on=c("Year","Sample")]
reask_agg[is.na(Freq),Freq:=0]
reask_agg[is.na(Loss),Loss:=0]
reask_agg[is.na(LossSe),LossSe:=0]
reask_agg[is.na(LossNe),LossNe:=0]
reask_agg[is.na(LossFl),LossFl:=0]
reask_agg[is.na(LossGl),LossGl:=0]
setorder(reask_agg,Year,-LossSe)
reask_agg[,RpSe:=.SD[,10000/.I],by="Year"]
setorder(reask_agg,Year,-LossNe)
reask_agg[,RpNe:=.SD[,10000/.I],by="Year"]
setorder(reask_agg,Year,-LossFl)
reask_agg[,RpFl:=.SD[,10000/.I],by="Year"]
setorder(reask_agg,Year,-LossGl)
reask_agg[,RpGl:=.SD[,10000/.I],by="Year"]
setorder(reask_agg,Year,-Loss)
reask_agg[,Rp:=.SD[,10000/.I],by="Year"]

# CALCULATE STATISTICS ON REASK
reask_stat=reask_agg[,.(q1Freq=quantile(Freq,0.1),q25Freq=quantile(Freq,0.25),q5Freq=quantile(Freq,0.5),q75Freq=quantile(Freq,0.75),q9Freq=quantile(Freq,0.9),
                    q1Loss=quantile(Loss,0.1),q25Loss=quantile(Loss,0.25),q5Loss=quantile(Loss,0.5),q75Loss=quantile(Loss,0.75),q9Loss=quantile(Loss,0.9),
                    q1LossSe=quantile(LossSe,0.1),q25LossSe=quantile(LossSe,0.25),q5LossSe=quantile(LossSe,0.5),q75LossSe=quantile(LossSe,0.75),q9LossSe=quantile(LossSe,0.9),
                    q1LossNe=quantile(LossNe,0.1),q25LossNe=quantile(LossNe,0.25),q5LossNe=quantile(LossNe,0.5),q75LossNe=quantile(LossNe,0.75),q9LossNe=quantile(LossNe,0.9),
                    meanLoss=sum(Loss)/10000,meanFreq=sum(Freq)/10000,meanLossSe=sum(LossSe)/10000, meanLossNe=sum(LossNe)/10000,
                    stdLoss=sd(Loss),stdFreq=sd(Freq),stdLossSe=sd(LossSe),stdLossNe=sd(LossNe)),by="Year"]

# READ PCS, CALCULATE ZONE LOSSES (NE, SE), REMOVE LOSSES FROM CARIBBEAN AND PACIFIC, 
pcs=fread(path_pcs)
pcs=pcs[Year>=yini & Year<=yend & grepl("Hurricane",EventName)]
pcs=pcs[,.(ZoneLoss=sum(StateTotal)),by=c("PCSId","Cat#","Year","DateFrom","DateTo","EventName","EventTotal","Zone")]
pcs=pcs[Zone%in%c("Northeast","Florida","Gulf","Southeast_Other","US_Other")]
pcs[,UsTotal:=sum(ZoneLoss),by=c("PCSId","Cat#","Year","DateFrom","DateTo","EventName","EventTotal")]
pcs=dcast(pcs,formula= PCSId+Year+EventName+EventTotal+UsTotal~Zone,value.var = "ZoneLoss",fill=0.0)
pcs=pcs[,.(PCSId,Year,EventName,EventTotal,
           PctNe=Northeast/UsTotal,
           PctSe=(Florida+Gulf+Southeast_Other)/UsTotal,
           PctFl=Florida/UsTotal,
           PctGl=Gulf/UsTotal,
           PcsLoss=UsTotal)]
setorder(pcs,-Year)
pcs=pcs[PcsLoss>0]

landfall=fread(path_landfall)[,.(PCSId,Landfall)]
pcs=landfall[pcs,on=c("PCSId")]

# READ INDEXATION FACTORS AND INDEX PCS LOSSES
index=fread(path_index)
pcs=index[pcs,on="Year"][,.(Year,EventName,Landfall,PctNe,PctSe,PctFl,PctGl,PcsLoss,Index)]
pcs[,IndexedLoss:=PcsLoss*(1+Index)^(2024-Year+1)]
pcs_occ=pcs[,.(Year,EventName,Landfall,Loss=IndexedLoss/1e9,LossNe=IndexedLoss*PctNe/1e9,LossSe=IndexedLoss*PctSe/1e9,LossFl=IndexedLoss*PctFl/1e9,LossGl=IndexedLoss*PctGl/1e9)]

# READ CLIMATE INDICES
climate=data.table(read_parquet(path_climate))
climate=climate[init_timing=="MIDJUNE",.(ONI=mean(SST_NINO34_ONI),RONI=mean(SST_NINO34_RONI),AMM=mean(SST_AMM)),by="season"]
climate=climate[season>=1985 & season<=2024,.(Year=season,ONI,RONI,AMM)]

# AGGREGATE PCS LOSSES ANNUALY AND FILL-IN MISSING YEARS
pcs_agg=pcs_occ[,.(Freq=.N,
           FreqNe=.SD[Landfall=="Northeast",.N],
           FreqSe=.SD[grepl("Southeast",Landfall),.N],
           Loss=sum(Loss),
           LossNe=sum(LossNe),
           LossSe=sum(LossSe),
           LossFl=sum(LossFl),
           LossGl=sum(LossGl)),by="Year"]
pcs_agg=pcs_agg[data.table(Year=years),on="Year"]
pcs_agg[is.na(Freq),c("Freq","FreqNe","FreqSe","Loss","LossNe","LossSe","LossFl","LossGl"):=0]

# CALCULATE RETURN PERIODS AND PERCENTILES OF PCS LOSSES ACCORDING TO VERISK AND REASK
pcs_agg[,Verisk_LossPc:=calc_percentile(Loss,verisk_agg$Loss)]
pcs_agg[,Verisk_FreqPc:=calc_percentile(Freq,verisk_agg$Freq)]
for (y in years) {
  pcs_agg[Year==y,Reask_LossPc:=calc_percentile(Loss,reask_agg[Year==y]$Loss)]
  pcs_agg[Year==y,Reask_FreqPc:=calc_percentile(Freq,reask_agg[Year==y]$Freq)]
  
} 

## ILW

setorder(reask_occ,Year,Sample)
setorder(reask_agg,Year,Sample)
setorder(verisk_occ,Sample)
setorder(verisk_agg,Sample)
setorder(pcs_occ,Year)
setorder(pcs_agg,Year)

risk_aversion=exp(seq(log(0.1),log(1000),length.out=51))
ilw_prices=fread(path_ilw)


# SINGLE ILWs

ILW=calc_asset_returns(notional=1,ilw_prices,years,pcs_occ)

ILW_stat=ILW[,.(Return_Realized=mean(Return_Realized),
                Sharpe_Realized=mean(Return_Realized)/sd(Return_Realized)),
             by=c("Strike","Zone")]


# PORTFOLIO OF ALL ILWs

PORTF_ALL=rbindlist(lapply(risk_aversion,optimize_portfolio,notional=1,ilw_prices,years,verisk_occ,reask_occ,pcs_occ))

ALLOCATIONS_ALL=PORTF_ALL[,.(Allocation_Verisk_mean=mean(Allocation_Verisk),
                           Allocation_Reask_mean=mean(Allocation_Reask)),by=c("Aversion","Strike","Zone")]

TIMESERIES_ALL=PORTF_ALL[,.(Return_Realized_Verisk=sum(Return_Realized_Verisk),
                            Return_Realized_Reask=sum(Return_Realized_Reask),
                            Return_Expected_Verisk=sum(Return_Expected_Verisk),
                            Return_Expected_Reask=sum(Return_Expected_Reask),
                            Variance_Verisk=sum(Variance_Verisk),
                            Variance_Reask=sum(Variance_Reask),
                            TVaR10_Verisk=sum(TVaR10_Verisk),
                            TVaR10_Reask=sum(TVaR10_Reask)),by=c("Aversion","Year")]

RETURNS_ALL=TIMESERIES_ALL[,.(Return_Realized_Verisk_mean=mean(Return_Realized_Verisk),
                              Return_Realized_Reask_mean=mean(Return_Realized_Reask),
                              Return_Realized_Verisk_comp=prod(1+Return_Realized_Verisk)^(1/(yend-yini+1))-1,
                              Return_Realized_Reask_comp=prod(1+Return_Realized_Reask)^(1/(yend-yini+1))-1,
                              Return_Realized_Verisk_std=sd(Return_Realized_Verisk),
                              Return_Realized_Reask_std=sd(Return_Realized_Reask),
                              Average_Drawdown_Verisk=sum(Return_Realized_Verisk[Return_Realized_Verisk<0])/length(years),
                              Average_Drawdown_Reask=sum(Return_Realized_Reask[Return_Realized_Reask<0])/length(years),
                              Prob_Drawdown_Verisk=sum(Return_Realized_Verisk<0)/length(years),
                              Prob_Drawdown_Reask=sum(Return_Realized_Reask<0)/length(years), 
                              Return_Expected_Verisk_mean=mean(Return_Expected_Verisk),
                              Return_Expected_Reask_mean=mean(Return_Expected_Reask),
                              Variance_Verisk_mean=mean(Variance_Verisk),
                              Variance_Reask_mean=mean(Variance_Reask),
                              TVaR10_Verisk_mean=mean(TVaR10_Verisk),
                              TVaR10_Reask_mean=mean(TVaR10_Reask)),by="Aversion"]

# REINVESTING STRATEGY

aversion_as=3
aversion_cs=30

PORTF_AS=optimize_portfolio(aversion=aversion_as,1,ilw_prices,years,verisk_occ,reask_occ,pcs_occ,reinvest=T)
TIMESERIES_AS=PORTF_AS[,.(Allocation_Verisk=sum(Allocation_Verisk),
                          Allocation_Reask=sum(Allocation_Reask),
                          Return_Realized_Verisk=sum(Return_Realized_Verisk),
                          Return_Realized_Reask=sum(Return_Realized_Reask)),by="Year"]
TIMESERIES_AS[,Return_Normalized_Verisk:=Return_Realized_Verisk/Allocation_Verisk]
TIMESERIES_AS[,Return_Normalized_Reask:=Return_Realized_Reask/Allocation_Reask]
TIMESERIES_AS[,Return_Normalized_Differnce:=Return_Normalized_Reask - Return_Normalized_Verisk]

PORTF_CS=optimize_portfolio(aversion=aversion_cs,1,ilw_prices,years,verisk_occ,reask_occ,pcs_occ,reinvest=T)
TIMESERIES_CS=PORTF_CS[,.(Allocation_Verisk=sum(Allocation_Verisk),
                            Allocation_Reask=sum(Allocation_Reask),
                            Return_Realized_Verisk=sum(Return_Realized_Verisk),
                            Return_Realized_Reask=sum(Return_Realized_Reask)),by="Year"]
TIMESERIES_CS[,Return_Normalized_Verisk:=Return_Realized_Verisk/Allocation_Verisk]
TIMESERIES_CS[,Return_Normalized_Reask:=Return_Realized_Reask/Allocation_Reask]
TIMESERIES_CS[,Return_Normalized_Differnce:=Return_Normalized_Reask - Return_Normalized_Verisk]


### SAVE THE RDATA SO IF CAN BE QUICKLY LOADED FOR PLOTTING

save(ALLOCATIONS_ALL,
     ILW,ilw_prices,ILW_stat,climate,
     index,landfall,pcs,pcs_agg,pcs_occ,
     PORTF_ALL,PORTF_AS,PORTF_CS,
     reask,reask_agg,reask_occ,reask_stat,
     RETURNS_ALL,TIMESERIES_ALL,TIMESERIES_AS,TIMESERIES_CS,
     verisk,verisk_agg,verisk_occ,verisk_split,verisk_stat,
     baseline,catalog,forecast,years,aversion_as,aversion_cs,
     path_ilw,path_index,path_landfall,path_out,path_pcs,path_project,path_reask,path_save,path_verisk,path_verisk_split,
     file=path_save)

# PLOTS

# LOAD OUTPUT RDATA TO AVOID RUNNING THE SIMULATION
load(path_save) 

png(filename=paste0(path_out,"01_",catalog,"_",baseline,"_",forecast,"_LOSS_TS.png"),width=12,height=8,units="cm",res=300)
par(mar=c(4,4,1,1))
plot(c(),c(),col="black",xlim=c(1985,2025),ylim=c(0,200),xlab="Year",ylab="Annual Industry Loss (USD billion)",xaxt="n",cex.axis=0.9)
points(pcs_agg$Year,pcs_agg$Loss,pch=16,cex=0.7,col="black")
lines(pcs_agg$Year,pcs_agg$Loss,type="h",lty=1,lwd=0.5,col="black")
lines(reask_stat$Year,reask_stat$meanLoss,col="red",lty=2)
polygon(c(1900,2050,2050,1900),c(verisk_stat$q1Loss,verisk_stat$q1Loss,verisk_stat$q9Loss,verisk_stat$q9Loss),col=rgb(0,0,1,0.2),border=NA)
lines(c(1900,2050),c(verisk_stat$meanLoss,verisk_stat$meanLoss),col="blue",lty=2)
polygon(c(years,rev(years)),c(reask_stat$q1Loss,rev(reask_stat$q9Loss)),col=rgb(1,0,0,0.2),border=NA)
axis(1,at=seq(1985,2025,length.out=9),las=3)
axis(1,at=seq(1985,2025),labels=rep("",41))
legend("topright",legend=c("Historical","Seasonal Risk Model","Long-term Risk Model"),col=c("black","red","blue"),lty=c(NA,2,2),pch=c(16,NA,NA),cex=0.8,bg=NA)
legend("topright",legend=c("","",""),x.intersp=12.75,col=c(NA,rgb(1,0,0,0.2),rgb(0,0,1,0.2)),pch=c(NA,15,15),pt.cex=c(NA,1.5,1.5),bty="n",cex=0.8)
text(x=1985,y=200,labels="(a)",adj=c(0,1))
dev.off()

png(filename=paste0(path_out,"03_",catalog,"_",baseline,"_",forecast,"_AEP_ONI.png"),width=12,height=8,units="cm",res=300)
par(mar=c(4,4,1,1))
setorder(verisk_agg,-Rp)
setorder(reask_agg,Year,-Rp)
plot(c(),c(),xlim=c(0,800),ylim=c(1000,1),log="y",xlab="Annual Modeled Industry Loss (USD billion)", ylab="Return Period (years)",yaxt="n")
for (y in years) {
  k=reask_agg$Year==y
  c=rgb(0.5,0.5,0.5,0.5)
  if (climate[Year==y]$RONI >= 0.5) c=rgb(0.13,0.55,0.13,0.5)
  if (climate[Year==y]$RONI <= -0.5) c=rgb(0.5,0,1,0.5)
  lines(reask_agg[k]$Loss,reask_agg[k]$Rp,col=c,lwd=0.75)
}
lines(verisk_agg$Loss,verisk_agg$Rp,col="black",type="l",lwd=2,lty=2)
axis(2,at=c(1,10,100,1000),labels=c("1","10","100","1,000"))
axis(2,at=c(seq(2,9)%o%c(1,10,100)),labels=rep("",8*3),lwd.ticks=0.5,tck=-0.025)
#image.plot(legend.only=TRUE, zlim= range(clr$val), col = clr$clr, legend.lab = "ONI", horizontal = T,smallplot = c(0.17,0.95,0.95,0.99), legend.cex=0.8, legend.shrik=0.8)
legend("topright",legend=c("Long-term Risk Model","Seasonal Risk Model Neutral","Seasonal Risk Model El Niño","Seasonal Risk Model La Niña"),
       col=c("black","grey","forest green","purple"),lty=c(2,1,1,1),lwd=c(2,.75,.75,.75),cex=0.8)
text(x=50,y=1,label="(c)",adj=c(0,1))
dev.off()


png(filename=paste0(path_out,"03_",catalog,"_",baseline,"_",forecast,"_AEP_AMM.png"),width=12,height=8,units="cm",res=300)
par(mar=c(4,4,1,1))
setorder(verisk_agg,-Rp)
setorder(reask_agg,Year,-Rp)
plot(c(),c(),xlim=c(0,800),ylim=c(1000,1),log="y",xlab="Annual Modeled Industry Loss (USD billion)", ylab="Return Period (years)",yaxt="n")
for (y in years) {
  k=reask_agg$Year==y
  if (climate[Year==y]$AMM > 0) c=rgb(1,0,0,0.5)
  if (climate[Year==y]$AMM <= 0) c=rgb(0,0,1,0.5)
  lines(reask_agg[k]$Loss,reask_agg[k]$Rp,col=c,lwd=0.75)
}
lines(verisk_agg$Loss,verisk_agg$Rp,col="black",type="l",lwd=2,lty=2)
axis(2,at=c(1,10,100,1000),labels=c("1","10","100","1,000"))
axis(2,at=c(seq(2,9)%o%c(1,10,100)),labels=rep("",8*3),lwd.ticks=0.5,tck=-0.025)
#image.plot(legend.only=TRUE, zlim= range(clr$val), col = clr$clr, legend.lab = "ONI", horizontal = T,smallplot = c(0.17,0.95,0.95,0.99), legend.cex=0.8, legend.shrik=0.8)
legend("topright",legend=c("Long-term Risk Model","Seasonal Risk Model AMM+","Seasonal Risk Model AMM-"),
       col=c("black","red","blue"),lty=c(2,1,1),lwd=c(2,.75,.75),cex=0.8)
text(x=50,y=1,label="(d)",adj=c(0,1))
dev.off()

png(filename=paste0(path_out,"03_",catalog,"_",baseline,"_",forecast,"_AAL_CLIMO.png"),width=12,height=8,units="cm",res=300)
par(mar=c(4,4,1,4))
setorder(climate,Year)
setorder(reask_stat,Year)
plot(reask_stat$meanLoss,climate$RONI,col="blue",pch=1,xlab="Average Annual Modeled Loss (USD billion)",ylab="ENSO Index Forecast",xlim=c(0,100),ylim=c(-4,4))
lin=lm(y~x, data=data.table(x=reask_stat$meanLoss,y=climate$RONI))
r=cor(reask_stat$meanLoss,climate$RONI)
abline(lin,col="blue",lty=2)
text(x=70,y=-1.6,label=paste0("r = ",sprintf("%.1f",r)),adj=c(0,1),col="blue",cex=0.8)
text(x=0,y=4,label="(b)",adj=c(0,1))
par(new=T)
plot(reask_stat$meanLoss,climate$AMM,col="red",pch=4,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,100),ylim=c(-5.5,5.5))
lin=lm(y~x, data=data.table(x=reask_stat$meanLoss,y=climate$AMM))
r=cor(reask_stat$meanLoss,climate$AMM)
abline(lin,col="red",lty=2)
text(x=70,y=4.2,label=paste0("r = ",sprintf("%.1f",r)),adj=c(0,1),col="red",cex=0.8)
axis(4,at=c(-5,-2.5,0,2.5,5),labels=c("-5","-2.5","0","2.5","5"))
mtext("AMM Index Forecast", side=4, line=2.9)
legend("bottomright",legend=c("RONI Index (linear fit)", "AMM Index (linear fit)"),pch=c(1,4),col=c("blue","red"),lty=2,cex=0.8)
dev.off()

png(filename=paste0(path_out,"04_",catalog,"_",baseline,"_",forecast,"_ILW_RETURN.png"),width=12,height=12,units="cm",res=300)
par(mar=c(5,5,1,1),xpd=TRUE)
ILW_stat[,ZoneStrike:=paste0(Zone," ",Strike)]
s_nw=paste("Nationwide",sort(ilw_prices[Zone=="Nationwide",Strike]))
s_ne=paste("Northeast",sort(ilw_prices[Zone=="Northeast",Strike]))
s_fl=paste("Florida",sort(ilw_prices[Zone=="Florida",Strike]))
s_gl=paste("Gulf",sort(ilw_prices[Zone=="Gulf",Strike]))
cols=c(colorRampPalette(c("DimGrey","LightGrey"))(length(s_nw)),
       colorRampPalette(c("DarkBlue","LightBlue"))(length(s_ne)),
       colorRampPalette(c("DarkRed","LightCoral"))(length(s_fl)),
       colorRampPalette(c("DarkGreen","LightGreen"))(length(s_gl))
       )
plot(c(),c(),xlim=c(0,20),ylim=c(36,1),xlab="Average Return (%)",ylab="",yaxt="n",yaxs="i")
axis(2,at=seq(1,36),label=unique(ILW_stat$ZoneStrike),las=2,cex.axis=0.6)
i=0
for (zs in ILW_stat$ZoneStrike) {
  i=i+1
  lines(c(0,ILW_stat[ZoneStrike==zs]$Return_Realized*100),c(i,i),col=cols[i],lty=1,lwd=4,lend=2)
} 
dev.off()

png(filename=paste0(path_out,"04_",catalog,"_",baseline,"_",forecast,"_ILW_COVARIANCE.png"),width=15,height=12,units="cm",res=300)
par(mar=c(5,5,1,7))
ILW[,ZoneStrike:=paste0(Zone," ",Strike)]
tmp=dcast(ILW,Year ~ factor(ZoneStrike,levels=unique(ILW$ZoneStrike)),value.var = "Return_Realized")
tmp=as.matrix(tmp,rownames="Year")
image(cov(tmp),ylim=c(1,0),zlim=c(0,0.2),col=colorRampPalette(c("white","red"))(50),xaxt="n",yaxt="n")
image.plot(cov(tmp),zlim=c(0,0.2),col=colorRampPalette(c("white","red"))(50),legend.only=T,legend.lab="Return Covariance",legend.line=3,smallplot=c(0.8,0.83,0.21,0.96))
axis(1,at=seq(0,1,length.out=nrow(ilw_prices)),label=unique(ILW$ZoneStrike),las=2,cex.axis=0.6)
axis(2,at=seq(0,1,length.out=nrow(ilw_prices)),label=unique(ILW$ZoneStrike),las=2,cex.axis=0.6)
dev.off()

png(filename=paste0(path_out,"04_",catalog,"_",baseline,"_",forecast,"_ILW_PRICES.png"),width=8,height=10,units="cm",res=300)
par(mar=c(4,4,1,1))
plot(c(),c(),xlab="ILW Strike (USD billion)",ylab="Risk Premium (%)",xlim=c(0,100),ylim=c(0,50))
points(ilw_prices[Zone=="Nationwide"]$Strike,ilw_prices[Zone=="Nationwide"]$Price_Howden*100,col=rgb(0.5,0.5,0.5,0.3),pch=1,cex=0.75)
points(ilw_prices[Zone=="Nationwide"]$Strike,ilw_prices[Zone=="Nationwide"]$Price_Aon*100,col=rgb(0.5,0.5,0.5,0.3),pch=4,cex=0.75)
lines(ilw_prices[Zone=="Nationwide"]$Strike,ilw_prices[Zone=="Nationwide"]$Price*100,col="dark grey",lty=2)
points(ilw_prices[Zone=="Florida"]$Strike,ilw_prices[Zone=="Florida"]$Price_Howden*100,col=rgb(1,0,0,0.3),pch=1,cex=0.75)
points(ilw_prices[Zone=="Florida"]$Strike,ilw_prices[Zone=="Florida"]$Price_Aon*100,col=rgb(1,0,0,0.3),pch=4,cex=0.75)
lines(ilw_prices[Zone=="Florida"]$Strike,ilw_prices[Zone=="Florida"]$Price*100,col="red",lty=2)
points(ilw_prices[Zone=="Northeast"]$Strike,ilw_prices[Zone=="Northeast"]$Price_Howden*100,col=rgb(0,0,1,0.3),pch=1,cex=0.75)
points(ilw_prices[Zone=="Northeast"]$Strike,ilw_prices[Zone=="Northeast"]$Price_Aon*100,col=rgb(0,0,1,0.3),pch=4,cex=0.75)
lines(ilw_prices[Zone=="Northeast"]$Strike,ilw_prices[Zone=="Northeast"]$Price*100,col="blue",lty=2)
points(ilw_prices[Zone=="Gulf"]$Strike,ilw_prices[Zone=="Gulf"]$Price_Howden*100,col=rgb(0,1,0,0.3),pch=1,cex=0.75)
points(ilw_prices[Zone=="Gulf"]$Strike,ilw_prices[Zone=="Gulf"]$Price_Aon*100,col=rgb(0,1,0,0.3),pch=4,cex=0.75)
lines(ilw_prices[Zone=="Gulf"]$Strike,ilw_prices[Zone=="Gulf"]$Price*100,col="green",lty=2)
legend("topright",legend=c("Nationwide","Florida","Gulf","Northeast"),
       lty=c(2,2,2,2),
       x.intersp=1,
       pch=NA,
       col=c("dark grey","red","green","blue"),cex=0.8,pt.cex=0.75)
legend("topright",legend=c("","","",""),pch=1,x.intersp=7.8,col=c(rgb(0.5,0.5,0.5,0.3),rgb(1,0,0,0.3),rgb(0,1,0,0.3),rgb(0,0,1,0.3)),pt.cex=0.75,cex=0.8,bty="n")
legend("topright",legend=c("","","",""),pch=4,x.intersp=7.0,col=c(rgb(0.5,0.5,0.5,0.3),rgb(1,0,0,0.3),rgb(0,1,0,0.3),rgb(0,0,1,0.3)),pt.cex=0.75,cex=0.8,bty="n")
text(x=0,y=50,labels="(b)",adj=c(0,1))
dev.off()


png(filename=paste0(path_out,"06_",catalog,"_",baseline,"_",forecast,"_PORTF_RETURN_REALIZED.png"),width=12,height=8,units="cm",res=300)
par(mar=c(4,4,1,4))
plot(c(),c(),xlab=expression("Risk Aversion,"~lambda),ylab="",xlim=c(0.1,1000),ylim=c(0,25),xaxt="n",log="x")
title(ylab=expression("Average Realized Return,"~bar(italic(r))~"(%)"),line=2.8)
axis(side=1,at=c(0.1,1,10,100,1000),labels=c("0.1","1","10","100","1,000"))
axis(side=1,at=c(seq(2,9)%o%c(0.1,1,10,100)),labels=rep("",8*4),lwd.ticks=0.5,tck=-0.025)
points(RETURNS_ALL$Aversion,RETURNS_ALL$Return_Realized_Verisk_mean*100,col="blue",type="p",cex=0.8,pch=1,lty=2)
points(RETURNS_ALL$Aversion,RETURNS_ALL$Return_Realized_Reask_mean*100,col="red",type="p",cex=0.8,pch=4,lty=2)
text(x=0.1,y=25,labels="(a)",adj=c(0,1))
par(new=T)
plot(c(),c(),xlab="",ylab="",xlim=c(0.1,1000),ylim=c(-40,40),xaxt="n",yaxt="n",log="x")
tmp=(RETURNS_ALL$Return_Realized_Reask_mean-RETURNS_ALL$Return_Realized_Verisk_mean)/RETURNS_ALL$Return_Realized_Verisk_mean*100
polygon(c(RETURNS_ALL$Aversion,rev(RETURNS_ALL$Aversion)),c(tmp,rep(0,length(tmp))),col=rgb(0.5,0.5,0.5,0.2),border=NA)
lines(c(0.01,10000),c(0,0),lty=2,lwd=0.5)
axis(side=4,at=c(-40,-20,0,20,40),labels=c("-40","-20","0","20","40"))
mtext("Relative Difference (%)", side=4, line=2.9)
legend("bottomleft",
       legend=c("Long-term risk model","Seasonal risk model"),
       col=c("blue","red"),pch=c(1,4),cex=0.8)
dev.off()

png(filename=paste0(path_out,"06_",catalog,"_",baseline,"_",forecast,"_PORTF_RETURN_COMPOUND.png"),width=12,height=8,units="cm",res=300)
par(mar=c(4,4,1,4))
plot(c(),c(),xlab=expression("Risk Aversion,"~lambda),ylab="",xlim=c(0.1,1000),ylim=c(0,20),xaxt="n",log="x")
title(ylab=expression("Compound growth rate,"~italic(r)[c]~"(%)"),line=2.8)
axis(side=1,at=c(0.1,1,10,100,1000),labels=c("0.1","1","10","100","1,000"))
axis(side=1,at=c(seq(2,9)%o%c(0.1,1,10,100)),labels=rep("",8*4),lwd.ticks=0.5,tck=-0.025)
points(RETURNS_ALL$Aversion,RETURNS_ALL$Return_Realized_Verisk_comp*100,col="blue",type="p",cex=0.8,pch=1,lty=2)
points(RETURNS_ALL$Aversion,RETURNS_ALL$Return_Realized_Reask_comp*100,col="red",type="p",cex=0.8,pch=4,lty=2)
text(x=0.1,y=20,labels="(b)",adj=c(0,1))
par(new=T)
plot(c(),c(),xlab="",ylab="",xlim=c(0.1,1000),ylim=c(-40,40),xaxt="n",yaxt="n",log="x")
tmp=(RETURNS_ALL$Return_Realized_Reask_comp-RETURNS_ALL$Return_Realized_Verisk_comp)/RETURNS_ALL$Return_Realized_Verisk_comp*100
polygon(c(RETURNS_ALL$Aversion,rev(RETURNS_ALL$Aversion)),c(tmp,rep(0,length(tmp))),col=rgb(0.5,0.5,0.5,0.2),border=NA)
lines(c(0.01,10000),c(0,0),lty=2,lwd=0.5)
axis(side=4,at=c(-40,-20,0,20,40),labels=c("-40","-20","0","20","40"))
mtext("Relative Difference (%)", side=4, line=2.9)
legend("bottomleft",
       legend=c("Long-term risk model","Seasonal risk model"),
       col=c("blue","red"),pch=c(1,4),cex=0.8)
dev.off()

png(filename=paste0(path_out,"06_",catalog,"_",baseline,"_",forecast,"_PORTF_SHARPE.png"),width=12,height=8,units="cm",res=300)
par(mar=c(4,4,1,4))
plot(c(),c(),xlab=expression("Risk Aversion,"~lambda),ylab=expression("Ex-Post Sharpe Ratio,"~italic(S)~"(%)"),xlim=c(0.1,1000),ylim=c(40,100),xaxt="n",yaxt="n",log="x")
axis(side=1,at=c(0.1,1,10,100,1000),labels=c("0.1","1","10","100","1,000"))
axis(side=1,at=c(seq(2,9)%o%c(0.1,1,10,100)),labels=rep("",8*4),lwd.ticks=0.5,tck=-0.025)
axis(side=2,at=c(40,60,80,100),labels=c("40","60","80","100"))
axis(side=2,at=c(50,70,90),labels=c("","",""))
points(RETURNS_ALL$Aversion,RETURNS_ALL$Return_Realized_Verisk_mean/RETURNS_ALL$Return_Realized_Verisk_std*100,col="blue",type="p",cex=1,pch=1,lty=2)
points(RETURNS_ALL$Aversion,RETURNS_ALL$Return_Realized_Reask_mean/RETURNS_ALL$Return_Realized_Reask_std*100,col="red",type="p",cex=1,pch=4,lty=2)
text(x=0.1,y=100,labels="(c)",adj=c(0,1))
par(new=T)
plot(c(),c(),xlab="",ylab="",xlim=c(0.1,1000),ylim=c(-50,50),xaxt="n",yaxt="n",log="x")
tmp=(RETURNS_ALL$Return_Realized_Reask_mean/RETURNS_ALL$Return_Realized_Reask_std-RETURNS_ALL$Return_Realized_Verisk_mean/RETURNS_ALL$Return_Realized_Verisk_std)/(RETURNS_ALL$Return_Realized_Verisk_mean/RETURNS_ALL$Return_Realized_Verisk_std)*100
polygon(c(RETURNS_ALL$Aversion,rev(RETURNS_ALL$Aversion)),c(tmp,rep(0,length(tmp))),col=rgb(0.5,0.5,0.5,0.2),border=NA)
lines(c(0.01,10000),c(0,0),lty=2,lwd=0.5)
axis(side=4,at=c(-50,-25,0,25,50),labels=c("-50","-25","0","25","50"))
mtext("Relative Difference (%)", side=4, line=2.9)
legend("bottomright",
       legend=c("Long-term Risk Model",
                "Seasonal Risk Model"),
       col=c("blue","red"),pch=c(1,4),cex=0.8)
dev.off()

png(filename=paste0(path_out,"06_",catalog,"_",baseline,"_",forecast,"_PORTF_AVERAGE_DRAWDOWN.png"),width=12,height=8,units="cm",res=300)
par(mar=c(4,4,1,4))
plot(c(),c(),xlab=expression("Risk Aversion,"~lambda),ylab="",xlim=c(0.1,1000),ylim=c(-15,0),xaxt="n",log="x")
title(ylab=expression("Average Drawdown,"~bar(italic(d))~"(%)"),line=2.8)
axis(side=1,at=c(0.1,1,10,100,1000),labels=c("0.1","1","10","100","1,000"))
axis(side=1,at=c(seq(2,9)%o%c(0.1,1,10,100)),labels=rep("",8*4),lwd.ticks=0.5,tck=-0.025)
points(RETURNS_ALL$Aversion,RETURNS_ALL$Average_Drawdown_Verisk*100,col="blue",type="b",cex=0.8,pch=1,lty=2)
points(RETURNS_ALL$Aversion,RETURNS_ALL$Average_Drawdown_Reask*100,col="red",type="b",cex=0.8,pch=4,lty=2)
text(x=0.1,y=0,labels="(d)",adj=c(0,1))
par(new=T)
plot(c(),c(),xlab="",ylab="",xlim=c(0.1,1000),ylim=c(-50,50),xaxt="n",yaxt="n",log="x")
tmp=(RETURNS_ALL$Average_Drawdown_Reask-RETURNS_ALL$Average_Drawdown_Verisk)/RETURNS_ALL$Average_Drawdown_Verisk*100
polygon(c(RETURNS_ALL$Aversion,rev(RETURNS_ALL$Aversion)),c(tmp,rep(0,length(tmp))),col=rgb(0.5,0.5,0.5,0.2),border=NA)
lines(c(0.01,10000),c(0,0),lty=2,lwd=0.5)
axis(side=4,at=c(-50,-25,0,25,50),labels=c("-50","-25","0","25","50"))
mtext("Relative Difference (%)", side=4, line=2.9)
legend("bottomright",
       legend=c("Long-term risk model","Seasonal risk model"),
       col=c("blue","red"),pch=c(1,4),cex=0.8)
dev.off()

y=2017
png(filename=paste0(path_out,"06_",catalog,"_",baseline,"_",forecast,"_PORTF_ALLOC_REASK_",y,".png"),width=16,height=8,units="cm",res=300)
s_nw=paste("Nationwide",sort(ilw_prices[Zone=="Nationwide",Strike]))
s_ne=paste("Northeast",sort(ilw_prices[Zone=="Northeast",Strike]))
s_fl=paste("Florida",sort(ilw_prices[Zone=="Florida",Strike]))
s_gl=paste("Gulf",sort(ilw_prices[Zone=="Gulf",Strike]))
PORTF_ALL[,`Zone Strike`:=factor(paste(Zone,Strike),levels=c(s_ne,s_fl,s_gl,s_nw))]
style=data.table(`Zone Strike`=c(s_ne,s_fl,s_gl,s_nw),
                 Color=c(colorRampPalette(c("DarkBlue","LightBlue"))(length(s_ne)),
                         colorRampPalette(c("DarkRed","LightCoral"))(length(s_fl)),
                         colorRampPalette(c("DarkGreen","LightGreen"))(length(s_gl)),
                         colorRampPalette(c("DimGrey","LightGrey"))(length(s_nw))))
ILW[,`Zone Strike`:=factor(paste(Zone,Strike),levels=c(s_ne,s_fl,s_gl,s_nw))]
ILW[isStrike==TRUE,Pattern:="stripe"]
ILW[isStrike==FALSE,Pattern:="none"]
pats=ILW[Year==y,.(`Zone Strike`,Pattern)]
style=pats[style,on="Zone Strike"]
ggplot(PORTF_ALL[Year==y],aes(x=Aversion,y=Allocation_Reask,group=`Zone Strike`))+
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.x=element_text(size=12,color="black"),
        axis.text.y=element_text(size=12,color="black"),
        legend.text=element_text(size=8,color="black"),
        legend.box.spacing = unit(15, "pt"),
        legend.key.size = unit(10, "pt"),
        plot.title = element_text(hjust = 0))+
  geom_area_pattern(aes(pattern=`Zone Strike`,fill=`Zone Strike`),pattern_density=0.005,pattern_spacing=0.03,pattern_fill="white",pattern_color="white")+
  scale_fill_manual(values=style$Color)+
  scale_x_log10(name=expression("Risk Aversion,"~lambda),breaks=c(0.1,1,10,100,1000),labels=c("0.1","1","10","100","1,000"),expand=c(0, 0))+
  #scale_x_log10(name="Risk Aversion",breaks=c(0.1,1,10,100,1000),labels=c("0.1","1","10","100","1,000"),expand=c(0, 0))+
  scale_y_continuous(name="Seasonal Allocation (%)",breaks=seq(0,1,by=0.2),labels=c("0","20","40","60","80","100"),expand=c(0, 0))+
  scale_pattern_manual(values=style$Pattern)+
  guides(x=guide_axis_logticks(short=0.75,mid=0.75,long=1.25))+
  labs(pattern=NULL,fill=NULL,title=paste("(f)",y))
dev.off()

y=2012
png(filename=paste0(path_out,"06_",catalog,"_",baseline,"_",forecast,"_PORTF_ALLOC_VERISK_",y,".png"),width=16,height=8,units="cm",res=300)
s_nw=paste("Nationwide",sort(ilw_prices[Zone=="Nationwide",Strike]))
s_ne=paste("Northeast",sort(ilw_prices[Zone=="Northeast",Strike]))
s_fl=paste("Florida",sort(ilw_prices[Zone=="Florida",Strike]))
s_gl=paste("Gulf",sort(ilw_prices[Zone=="Gulf",Strike]))
PORTF_ALL[,`Zone Strike`:=factor(paste(Zone,Strike),levels=c(s_ne,s_fl,s_gl,s_nw))]
style=data.table(`Zone Strike`=c(s_ne,s_fl,s_gl,s_nw),
                 Color=c(colorRampPalette(c("DarkBlue","LightBlue"))(length(s_ne)),
                         colorRampPalette(c("DarkRed","LightCoral"))(length(s_fl)),
                         colorRampPalette(c("DarkGreen","LightGreen"))(length(s_gl)),
                         colorRampPalette(c("DimGrey","LightGrey"))(length(s_nw))))
ILW[,`Zone Strike`:=factor(paste(Zone,Strike),levels=c(s_ne,s_fl,s_gl,s_nw))]
ILW[isStrike==TRUE,Pattern:="stripe"]
ILW[isStrike==FALSE,Pattern:="none"]
pats=ILW[Year==y,.(`Zone Strike`,Pattern)]
style=pats[style,on="Zone Strike"]
ggplot(PORTF_ALL[Year==y],aes(x=Aversion,y=Allocation_Verisk,group=`Zone Strike`))+
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.x=element_text(size=12,color="black"),
        axis.text.y=element_text(size=12,color="black"),
        legend.text=element_text(size=8,color="black"),
        legend.box.spacing = unit(15, "pt"),
        legend.key.size = unit(10, "pt"),
        plot.title = element_text(hjust = 0))+
  geom_area_pattern(aes(pattern=`Zone Strike`,fill=`Zone Strike`),pattern_density=0.05,pattern_spacing=0.03,pattern_fill="white",pattern_color="white")+
  scale_fill_manual(values=style$Color)+
  scale_x_log10(name="Risk Aversion",breaks=c(0.1,1,10,100,1000),labels=c("0.1","1","10","100","1,000"),expand=c(0, 0))+
  scale_y_continuous(name="Long-term Allocation (%)",breaks=seq(0,1,by=0.2),labels=c("0","20","40","60","80","100"),expand=c(0, 0))+
  scale_pattern_manual(values=style$Pattern)+
  guides(x=guide_axis_logticks(short=0.75,mid=0.75,long=1.25))+
  labs(pattern=NULL,fill=NULL,title=paste("(c)",y))
dev.off()

png(filename=paste0(path_out,"07_",catalog,"_",baseline,"_",forecast,"_INVESTMENT_AGG_TS.png"),width=12,height=8,units="cm",res=300)
par(mar=c(4,4,1,4))
plot(c(),c(),xlab="Year",ylab="",xlim=c(1985,2025),ylim=c(0,150),xaxt="n")
title(ylab=expression("Invested Capital,"~italic(C)[y]),line=2.8)
axis(1,at=seq(1985,2025,length.out=9),las=3)
axis(1,at=seq(1985,2025),labels=rep("",41))
points(TIMESERIES_AS$Year,TIMESERIES_AS$Allocation_Reask+TIMESERIES_AS$Return_Realized_Reask,col="red",type="b",pch=4,lty=2)
points(TIMESERIES_AS$Year,TIMESERIES_AS$Allocation_Verisk+TIMESERIES_AS$Return_Realized_Verisk,col="blue",type="b",pch=1,lty=2)
text(x=2025,y=150,labels="(a) Aggressive",adj=c(1,1))
par(new=T)
plot(c(),c(),xlab="",ylab="",xlim=c(1985,2025),ylim=c(-250,250),xaxt="n",yaxt="n",log="x")
tmp=(TIMESERIES_AS$Allocation_Reask+TIMESERIES_AS$Return_Realized_Reask-TIMESERIES_AS$Allocation_Verisk-TIMESERIES_AS$Return_Realized_Verisk)/(TIMESERIES_AS$Allocation_Verisk+TIMESERIES_AS$Return_Realized_Verisk)*100
polygon(c(TIMESERIES_AS$Year,rev(TIMESERIES_AS$Year)),c(tmp,rep(0,length(tmp))),col=rgb(0.5,0.5,0.5,0.2),border=NA)
lines(c(1985,2025),c(0,0),lty=2,lwd=0.5)
axis(side=4,at=c(-250,-125,0,125,250),labels=c("-250","-125","0","125","250"))
mtext("Relative Difference (%)", side=4, line=2.9)
legend("topleft",
       legend=c("Long-term Risk Model","Seasonal Risk Model"),
       col=c("blue","red"),pch=c(1,4),cex=0.8)
dev.off()

png(filename=paste0(path_out,"07_",catalog,"_",baseline,"_",forecast,"_INVESTMENT_CON_TS.png"),width=12,height=8,units="cm",res=300)
par(mar=c(4,4,1,4))
plot(c(),c(),xlab="Year",ylab="",xlim=c(1985,2025),ylim=c(0,10),xaxt="n")
title(ylab=expression("Invested Capital,"~italic(C)[y]),line=2.8)
axis(1,at=seq(1985,2025,length.out=9),las=3)
axis(1,at=seq(1985,2025),labels=rep("",41))
points(TIMESERIES_CS$Year,TIMESERIES_CS$Allocation_Reask+TIMESERIES_CS$Return_Realized_Reask,col="red",type="b",pch=4,lty=2)
points(TIMESERIES_CS$Year,TIMESERIES_CS$Allocation_Verisk+TIMESERIES_CS$Return_Realized_Verisk,col="blue",type="b",pch=1,lty=2)
text(x=2025,y=10,labels="(b) Conservative",adj=c(1,1))
par(new=T)
plot(c(),c(),xlab="",ylab="",xlim=c(1985,2025),ylim=c(-50,50),xaxt="n",yaxt="n",log="x")
tmp=(TIMESERIES_CS$Allocation_Reask+TIMESERIES_CS$Return_Realized_Reask-TIMESERIES_CS$Allocation_Verisk-TIMESERIES_CS$Return_Realized_Verisk)/(TIMESERIES_CS$Allocation_Verisk+TIMESERIES_CS$Return_Realized_Verisk)*100
polygon(c(TIMESERIES_CS$Year,rev(TIMESERIES_CS$Year)),c(tmp,rep(0,length(tmp))),col=rgb(0.5,0.5,0.5,0.2),border=NA)
lines(c(1985,2025),c(0,0),lty=2,lwd=0.5)
axis(side=4,at=c(-50,-25,0,25,50),labels=c("-50","-25","0","25","50"))
mtext("Relative Difference (%)", side=4, line=2.9)
legend("topleft",
       legend=c("Long-term Risk Model","Seasonal Risk Model"),
       col=c("blue","red"),pch=c(1,4),cex=0.8)
dev.off()