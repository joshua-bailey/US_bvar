

#-----------------------------------
# Functions for out-of-sample forecasting
#-----------------------------------
backtesting=function(df,targets,start,end,by,h,var_names0,trans,p,n_draw,n_burn,start_tb,zerobound=FALSE,conditional=FALSE,conditional_vars=NULL){
  # This function generates one-step ahead backtesting results.
  start_bt=proc.time()
  
  # select variables for backtesting
  selected_vars=c('GDP','Consumption','Investment (nonres)',
                  'Exports','Imports',
                  'PCE Price Index','Core PCE Price Index','CPI-U',
                  'Payroll Employment','Labor Force','Unemployment Rate',
                  'TR10y Rate')
  
  inds=which(var_names0 %in% selected_vars)
  inds0=inds+1
  df=df[,append(inds0,1)]
  var_names0=var_names0[inds]
  trans=trans[inds]
  lag=4 # for yoy growth rates
  h_max=12 # for anchoring
  
  # forecast target
  tbs=list()
  for (target in targets){
    ind=which(var_names0==target)
    
    tb=data.frame(date=seq(from=start_tb, to=end, by=by),actual=NA,predicted=NA)
    tb[tb$date<start,'actual']=df[(df$yq>=start_tb)&(df$yq<start),(ind)]
    tb[tb$date<start,'predicted']=tb[tb$date<start,'actual']
    tbs[[target]]=tb
  }
  
  for (t in seq(from=start, to=end, by=by)){
    if (t==floor(t)){print(paste0('forecasting date:',t))}
    
    # prepare data
    data=df[df$yq<=(t-h/4),]
    data=data[,colnames(data)!='yq']
    data=transform_data_cmac(data,var_names0,trans)
    
    # forecasting
    tb0=backtesting_one_time(data,var_names0,p,n_draw,n_burn,h,h_max,zerobound)
    
    for (target in targets){
      ind=which(var_names0==target)
      f_value=tb0[ind]
      tb=tbs[[target]]
      if (trans[ind]=='log100'){tb[tb$date==t,'predicted']=exp(f_value/100)}
      else {tb[tb$date==t,'predicted']=f_value/100}
      tb[tb$date==t,'actual']=df[df$yq==t,(ind)]
      tbs[[target]]=tb
    }
  }
  
  # to growth rates
  for (target in targets){
    ind=which(var_names0==target)
    tb=tbs[[target]]
    if (trans[ind]=='log100'){
      tb$predicted=(tb$predicted/shift(tb$predicted,n=lag,fill=NA,type='lag')-1)*100
      tb$actual=(tb$actual/shift(tb$actual,n=lag,fill=NA,type='lag')-1)*100
    }
    tbs[[target]]=tb
  }
  
  # compare with others
  tbs=lapply(targets,backtesting_compare,tbs,h,start,end)
  names(tbs)=targets
  
  running_time=(proc.time()-start_bt)[3]
  print(paste('The elapsed time for backtesting: ',running_time/60,' minutes'))
  
  return(tbs)
}

backtesting_cbo=function(df,targets,start,end,by,h,var_names0,trans,p,n_draw,n_burn,start_tb,zerobound=FALSE,conditional=FALSE,conditional_vars=NULL){
  # This function generates one-step ahead backtesting results
  # to be compared with CBO forecasts
  
  start_bt=proc.time()
  
  # select variables for backtesting
  selected_vars=c('GDP','Consumption','Investment (nonres)',
                  'Exports','Imports',
                  'PCE Price Index','Core PCE Price Index','CPI-U',
                  'Payroll Employment','Labor Force','Unemployment Rate',
                  'TR10y Rate')
  
  inds=which(var_names0 %in% selected_vars)
  inds0=inds+1
  df=df[,append(inds0,1)]
  var_names0=var_names0[inds]
  trans=trans[inds]
  lag=4 # for annual growth rates
  h_max=12 # for anchoring
  by=1 # interval 1 year
  
  # prepare tables
  tbs=list()
  for (target in targets){
    ind=which(var_names0==target)
    tb=data.frame(date=seq(from=(start-lag/4), to=end, by=by),actual=NA,predicted=NA)
    average_quarterly=mean(df[(df$yq>=(start-lag/4-lag/4+1/4))&(df$yq<=(start-lag/4)),ind])
    tb[tb$date==(start-lag/4),'actual']=average_quarterly
    tb[tb$date<start,'predicted']=tb[tb$date<start,'actual']
    tbs[[target]]=tb
  }
  
  # cbo forecasts to be conditioned
  if (conditional==TRUE){
    if (grepl('cbo',conditional_vars)==TRUE){cbo=load_cbo_forecasts('1y')}
  }
  
  # one-step ahead forecasting
  for (t in seq(from=start, to=end, by=by)){
    if ((by==1/4) & ( t==floor(t))){print(paste0('forecasting date:',t))}
    if (by==1) {print(paste0('forecasting date:',t))}
    
    # prepare data
    data=df[df$yq<=(t-h/4),]
    data=data[,colnames(data)!='yq']
    data=transform_data_cmac(data,var_names0,trans)
    
    # forecasting
    tb0=backtesting_one_time(data,var_names0,p,n_draw,n_burn,h,h_max,
                             zerobound,average=TRUE,
                             conditional=conditional,conditional_vars=conditional_vars,
                             forecasting_date=t,cbo_forecasts=cbo)
    
    # update tables
    for (target in targets){
      ind=which(var_names0==target)
      f_value=tb0[ind]
      tb=tbs[[target]]
      if (trans[ind]=='log100'){tb[tb$date==t,'predicted']=exp(f_value/100)}
      else {tb[tb$date==t,'predicted']=f_value/100}
      average_quarterly=mean(df[(df$yq>=(t-lag/4+1/4))&(df$yq<=t),ind])
      tb[tb$date==t,'actual']=average_quarterly
      tbs[[target]]=tb
    }
  }
  
  # to growth rates (on averages)
  for (target in targets){
    ind=which(var_names0==target)
    tb=tbs[[target]]
    if (trans[ind]=='log100'){
      tb$predicted=(tb$predicted/shift(tb$actual,n=1,fill=NA,type='lag')-1)*100
      tb$actual=(tb$actual/shift(tb$actual,n=1,fill=NA,type='lag')-1)*100
    }
    tbs[[target]]=tb
  }
  
  # add other forecasts
  #tbs=lapply(targets,backtesting_compare,tbs,1,start,end) # use only last_value
  names(tbs)=targets
  
  running_time=(proc.time()-start_bt)[3]
  print(paste('The elapsed time for backtesting: ',running_time/60,' minutes'))
  
  return(tbs)
}

backtesting_one_time=function(data,var_names0,p,n_draw,n_burn,h,h_max,zerobound,average=FALSE,conditional=FALSE,conditional_vars=NULL,forecasting_date=NULL,cbo_forecasts=NULL){
  # This function generates forecasts for one time.
  res=estimate_bvar_bt(data,var_names0,p,n_draw,n_burn)
  Y=rbind(as.matrix(data),matrix(NA,h_max,ncol(data)))
  if (conditional==TRUE){
    if (grepl('cbo',conditional_vars)==TRUE){
      Y=conditional_cbo(Y,data,cbo_forecasts,forecasting_date,conditional_vars,var_names0,h)
    }
  }
  end_ind=nrow(data)
  n_sim=n_draw-n_burn
  base=NULL
  #scenario_name='base'
  scenario_name='PCE Price Index 6% higher in 3 years' # inflation targeting
  forecasts=generate_one_scenario(Y,end_ind,df,res,p,h_max,n_sim,
                                  var_names0,base,scenario_name,
                                  zerobound=zerobound,verbose=FALSE)
  
  #forecasts=entropic_tilting_pce(forecasts,var_names0)
  #inds_stable=check_stability(res,p)
  #fm=apply(forecasts[,,inds_stable==1],c(1,2),median)
  fm=apply(forecasts,c(1,2),median)
  #fm1=apply(forecasts,c(1,2),mean)
  fm_h=fm[h,]
  if (average==TRUE){fm_h=apply(fm[1:h,],2,mean)} # for average
  return(fm_h)
}

estimate_bvar_bt=function(data,var_names0,p,n_draw,n_burn){
  # This function estimates BVAR for backtesting.
  
  # prior
  soc0=bv_soc(mode = 0.2, sd = 1, min = 1e-04, max = 50) # sum-of-coefficients
  sur0=bv_soc(mode = 0.5, sd = 1, min = 1e-04, max = 50) # single-unit-root
  
  mode_psi=c(1,1,1,1,1,
             1,1,1,1,1,
             50,50)
  
  mn0=bv_mn(bv_lambda(mode=0.06),
            bv_alpha(mode=2),
            bv_psi(mode=mode_psi)) #scale = 1, shape = 1
  #bv_psi(mode=mode_psi))
  
  mh0=bv_metropolis(scale_hess = c(0.01))
  hyper0=c('lambda')
  priors0=bv_priors(hyper=hyper0,mn=mn0)
  
  # estimation  
  res=bvar(data=data,lags=p,n_draw = n_draw,n_burn = n_burn,n_thin = 1L,
           priors = priors0,mh = mh0,
           fcast = NULL,irf = NULL,verbose =FALSE)

  return(res)
}

backtesting_plotting=function(file_name,tb,target,start_bt,start_plot){
  # This function generates charts.
  
  pdf(file_name)
  fig=ggplot(data=tb[tb$date>=(start_plot+4),])+
    geom_line(aes(x=date,y=actual,color='actual'))+
    geom_point(aes(x=date,y=actual,color='actual'))+
    geom_line(data=tb[tb$date>=start_bt,],
              aes(x=date,y=predicted,color='predicted'))+
    geom_point(data=tb[tb$date>=start_bt,],
               aes(x=date,y=predicted,color='predicted'))+
    scale_colour_manual("",
                        breaks=c('actual','predicted'),
                        values=c('blue','red'))+
    ggtitle(target)+xlab(NULL)+ylab('%')+
    theme(legend.position = 'bottom')
  #print(fig)
  grid.arrange(grobs=list(fig,fig),nrow=2,ncol=1)
  dev.off()
}

backtesting_compare=function(target,tbs,h,start_bt,end_bt){
  # This function generates RMSEs of other forecasting schemes.
  tb=tbs[[target]]
  # one value
  tb$mean_tb=mean(tb[tb$date<start_bt,'actual'],na.rm=TRUE)
  
  # last value
  tb$last_value=shift(tb$actual,n=h,fill=NA,type='lag')
  
  # average of last 4 values
  tb$avg_lag=NA
  lag=4
  for (t in seq(from=start_bt,to=end_bt,by=1/4)){
    average_lag4=mean(tb[(tb$date<=(t-h/4)) & (tb$date>=(t-h/4-(lag-1)/4)),'actual'])
    tb[tb$date==t,'avg_lag']=average_lag4
    
  }
  
  tp=tb[tb$date>=start_bt,]
  rmse_mean_tb=sqrt(mean((tp$actual-tp$mean_tb)^2))
  rmse_last_value=sqrt(mean((tp$actual-tp$last_value)^2))
  rmse_avg_lag=sqrt(mean((tp$actual-tp$avg_lag)^2))
  
  rmse_tb=data.frame(mean=c(rmse_mean_tb),
                     last_value=c(rmse_last_value),
                     avg_lag=c(rmse_avg_lag))
  
  #return(list(tb=tb,rmse=rmse_tb))
  return(tb)
}

load_cbo_forecasts=function(sheet_name){
  # This function loads CBO forecasts.
  foldername="../data/"
  filename=paste0(foldername,'cbo_forecasts.xlsx')
  sheet_name=sheet_name
  cbo=read_excel(filename,sheet=sheet_name)
  colnames(cbo)=tolower(colnames(cbo))
  return(cbo)
}

conditional_cbo=function(Y,data,cbo_forecasts,forecasting_date,conditional_vars,var_names0,h){
  # This function imposes conditions based on CBO forecasts.
  
  fcst=cbo_forecasts[cbo_forecasts$year==floor(forecasting_date),][[conditional_vars]]
  var_name=ifelse(grepl('qgdp',conditional_vars),'GDP',NA)
  ind=which(var_names0 %in% var_name)
  end_ind=nrow(data)
  # target growth rates of average
  avg_prev=mean(exp(Y[(end_ind-3):end_ind,ind]/100))
  x=4*avg_prev*(1+fcst/100)/(1+(1+fcst/400)+(1+fcst/400)^2+(1+fcst/400)^3)
  Y[(end_ind+1),ind]=log(x)*100
  Y[(end_ind+2),ind]=log(x*(1+fcst/400))*100
  Y[(end_ind+3),ind]=log(x*(1+fcst/400)^2)*100
  Y[(end_ind+4),ind]=log(x*(1+fcst/400)^3)*100
  return(Y)
}
