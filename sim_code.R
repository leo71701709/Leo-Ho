rm(list = ls())
source('~/Leo/leo_v0.r')
library(lattice);
library(parallel)
library(rbenchmark)
library(stringr)


# simulation setting
nulltype=2
confounder=1
calibration=1
# nulltype=1;confounder=0;calibration=0
# nulltype=2;confounder=0;calibration=0
# nulltype=1;confounder=1;calibration=0
# nulltype=1;confounder=1;calibration=1
# nulltype=2;confounder=1;calibration=0
# nulltype=2;confounder=1;calibration=1

parallel_cal=1
Corenums= 40
PC=FALSE
# parallel_cal=1;Corenums= 6 ;PC=TRUE
# parallel_cal=1;Corenums= 40 ;PC=FALSE

alphaXcoeff=1
sample_size_list = c(100,200,500,1000)
repeat_size = 1000
##for nulltype=2, eta=0
eta=0.8
myseed = 9


# data generation setting
if(nulltype==2){
  zeta_list=c(0)
}else{
  zeta_list=seq(0, 0.56, by = 0.07)
  # zeta_list=c(0)
}

alphaX = alphaXcoeff * confounder

pvals.rec<-list()
power.rec<-list()
Type1error<-matrix(0L, nrow = 2*length(sample_size_list), ncol = 5)
data = vector(mode = 'list', length = repeat_size)

dat.gen<-function(sample_size,zeta){
  data = vector(mode = 'list', length = repeat_size)
  set.seed(myseed)
  
  X = c(rep(0, sample_size*0.49), rep(1, sample_size*0.51)) * confounder
  Z = ((X + rnorm(sample_size)) > 0.5) + 1
  for(counter in 1:repeat_size){
    set.seed(counter + myseed)
    
    if(nulltype==1){
      if(confounder==1){
        T1 = rweibull(sample_size, scale = exp(zeta*Z+ alphaX * X), shape = 1)
        T2 = eta*T1 + 0.5* rweibull(sample_size, scale = exp(0+ alphaX * X), shape = 1)
        
      }else{
        T1 = rweibull(sample_size, scale = exp(-1.6+zeta*Z+ alphaX * X), shape = 1)
        T2 = eta*T1 + 1* rweibull(sample_size, scale = exp(0+ alphaX * X), shape = 1)
      }
    }else{
      if(confounder==1){
        T1 = rweibull(sample_size, scale = exp(-1.2+ 0.56*Z+ alphaX * X), shape = 1)
        T2 = rweibull(sample_size, scale = exp(1.2+ 0*Z+ alphaX * X), shape = 1.5)
        
      }else{
        T1 = rweibull(sample_size, scale = exp(-1.6+ 0.56*Z+ alphaX * X), shape = 1)
        T2 = rweibull(sample_size, scale = exp(-1.6+ 0.56*Z+ alphaX * X), shape = 1.5)
      }
    }
    C = rweibull(sample_size, scale = 2, shape = 5)
    
    d2 = T2 < C
    T2 = pmin(T2, C)
    d1 = T1 < T2
    T1 = pmin(T1, T2)
    
    if(confounder+calibration<2){
      data[[counter]] = data.frame(T1 = T1, d1 = d1, T2 = T2, d2 = d2, Z = Z)
    }else{
      data[[counter]] = data.frame(T1 = T1, d1 = d1, T2 = T2, d2 = d2, Z = Z, X = X)
    }
  }
  return(data)
}


# Function_test
IE.test<-function(dat,parallel_cal){
  
  Calpvals<-function(counter){
    df<-dat[[counter]]
    result<-0
    try(result <-CASCR(df, plot_result = F, num_of_cores = 1, HO = TRUE, timer = FALSE, intervention = c(2, 1),threshold = 1e-15,variance_method = 'old'), silent = TRUE)
    if (length(result)>1){
      iehat=tcrossprod(t(result$IE$HO$weight),result$IE$HO$integrand)
      ie.vargamma<-result$IE$variance$alpha_variance
      ie.varbeta<-result$IE$variance$beta_variance
      ##pvalue
      pval.WLR<-pchisq(iehat^2/(ie.vargamma+ie.varbeta), df=1, lower.tail=FALSE)
      pval1<-pchisq(iehat^2/ie.vargamma, df=1, lower.tail=FALSE)
      pval2<-pchisq(iehat^2/ie.varbeta, df=1, lower.tail=FALSE)
      pval.IUT<-max(pval1, pval2) 
      rowData <- c(pval.WLR,pval.IUT, pval1,  pval2)
      return(rowData)
    }else{
      print(counter)
      return(0)
    }
  }
  
  if(parallel_cal>0){
    cl <- makeCluster(Corenums)
    clusterExport(cl, c("repeat_size","dat",'my_sort_mat',"CASCR",'downsample_func','get_WdL','my_rep_row',"data_preprocess","form_matrix","get_pd","get_beta_variance","get_alpha_variance","compute_variance","get_counterfactual_hazard","df_shift_to_cal_level","estimate_alpha","inv_coxinformation","estimate_effect","mycoxph","my_eva_fun","make_small","get_position","my_basehaz"))
    clusterEvalQ(cl, c(library(lattice)))
    output <- parLapply(cl, c(1:repeat_size), Calpvals)  
    stopCluster(cl)
    return(output)
  }else{
    error_counter_index<-0
    pval.WLR<-rep(0, length(dat))
    pval1<-rep(0, length(dat)) 
    pval2<-rep(0, length(dat))
    pval.IUT<-rep(0, length(dat))
    for(counter in 1:repeat_size){
      df<-dat[[counter]]
      result<-0
      try(result <-CASCR(df, plot_result = F, num_of_cores = 1, HO = TRUE, timer = FALSE, intervention = c(2, 1),threshold = 1e-6), silent = TRUE)
      if (length(result)>1){
        iehat=tcrossprod(t(result$IE$HO$weight),result$IE$HO$integrand)
        ie.vargamma<-result$IE$variance$alpha_variance
        ie.varbeta<-result$IE$variance$beta_variance
        ##pvalue
        pval.WLR[counter]<-pchisq(iehat^2/(ie.vargamma+ie.varbeta), df=1, lower.tail=FALSE)
        pval1[counter]<-pchisq(iehat^2/ie.vargamma, df=1, lower.tail=FALSE)
        pval2[counter]<-pchisq(iehat^2/ie.varbeta, df=1, lower.tail=FALSE)
        pval.IUT[counter]<-max(pval1[counter], pval2[counter])
      }else{ 
        error_counter_index<-c(error_counter_index,counter)
        # print(counter)
        next
      }
    }
    while(length(error_counter_index)>1){
      pval.WLR <- pval.WLR[-error_counter_index[length(error_counter_index)]]
      pval1 <- pval1[-error_counter_index[length(error_counter_index)]]
      pval2 <- pval2[-error_counter_index[length(error_counter_index)]]
      pval.IUT <- pval.IUT[-error_counter_index[length(error_counter_index)]]
      error_counter_index <- error_counter_index[-length(error_counter_index)]
    }
    return(cbind(pval.WLR, pval.IUT, pval1, pval2))
  }
}


# Setting for Record
Setting_table<-function(){
  #sample_size_list = c(100,200,500,1000)
  str1<-str_c("sample_size_list = ",sample_size_list,sep="")
  #repeat_size = 1000
  str2<-str_c("repeat_size = ",repeat_size,sep="")
  #eta=0.8
  str3<-str_c("eta = ",eta,sep="")
  #alphaXcoeff=0.3
  str4<-str_c("alphaXcoeff = ",alphaXcoeff,sep="")
  #parallel_cal=1
  str5<-str_c("parallel_cal = ",parallel_cal,sep="")
  #nulltype=1
  str6<-str_c("nulltype = ",nulltype,sep="")
  #confounder=1
  str7<-str_c("confounder = ",confounder,sep="")
  #calibration=0
  str8<-str_c("calibration = ",calibration,sep="")
  #myseed = 9
  str9<-str_c("myseed = ",myseed,sep="")
  table<-c(str1,str2,str3,str4,str5,str7,str7,str8,str9)
  write.table(table, file = paste(data_location,'setting.R',  sep = ''),sep = " ", quote = FALSE, na = "NA", row.names=FALSE)
}


# Type1error Table for Latex
Type1error_table<-function(Type1error){
  Type1error<-round(Type1error, digits = 1)
  str0<-"\\% & $p<0.2$ & $p<0.1$ & $p<0.05$ & $p<0.01$ & $p<0.005$ \\\\"
  str1<-"\\hline"
  str2<-"\\multicolumn{6}{c}{WLR} \\\\"
  str3<-"\\multicolumn{6}{c}{IUT} \\\\"
  #\hline
  #\% & $p<0.2$ & $p<0.1$ & $p<0.05$ & $p<0.01$ & $p<0.005$ \\
  table<-c(str1, str0)
  #\hline
  #\multicolumn{6}{c}{WLR} \\
  table<-c(table, str1, str2)
  #\hline
  #m=100 & 6.3(20) & 0.3(10) & 0(5) & 0(1) & 0(0.5) \\
  for(i in 1:length(sample_size_list)){
    table<-c(table, str1, str_c("m=",sample_size_list[i], " & ", as.character(Type1error[i,1]), "(20) & ", as.character(Type1error[i,2]), "(10) & ", as.character(Type1error[i,3]), "(5) & ", as.character(Type1error[i,4]), "(1) & ", as.character(Type1error[i,5]), "(0.5) \\\\", sep = ""))
  }
  table<-c(table, str1, str3)
  #\hline
  #m=100 & 6.3(20) & 0.3(10) & 0(5) & 0(1) & 0(0.5) \\
  for(i in (1+length(sample_size_list)):(2*length(sample_size_list))){
    table<-c(table, str1, str_c("m=",sample_size_list[i-4], " & ", as.character(Type1error[i,1]), "(20) & ", as.character(Type1error[i,2]), "(10) & ", as.character(Type1error[i,3]), "(5) & ", as.character(Type1error[i,4]), "(1) & ", as.character(Type1error[i,5]), "(0.5) \\\\", sep = ""))
  }
  #\hline
  table<-c(table, str1)
  write.table(table, file = paste(data_location,'type1error.txt',  sep = ''),sep = " ", quote = FALSE, na = "NA", row.names=FALSE)
}


# Data analysis
if(confounder+calibration==2){
  data_location<-paste('~/Leo/null', nulltype, eta, 'confounder', alphaXcoeff, 'calibration/',  sep = '')
}else if(confounder+calibration==1){
  data_location<-paste('~/Leo/null', nulltype, eta, 'confounder', alphaXcoeff,'/',  sep = '')
}else{
  data_location<-paste('~/Leo/null', nulltype, eta,'/',  sep = '')
}
if(PC){ data_location<-c('')}
Setting_table()
for(count.m in 1:length(sample_size_list)){
  power<-matrix(0L,nrow = 2, ncol = length(zeta_list))
  for(count.zeta in 1:length(zeta_list)){
    sample_size<-sample_size_list[count.m]
    zeta<-zeta_list[count.zeta]
    dat<-dat.gen(sample_size,zeta)
    pvals<-IE.test(dat,parallel_cal)
    
    if(parallel_cal>0){
      for(j in repeat_size:1){
        if(sum(pvals[[j]][1])=='NaN'){
          pvals<-pvals[-j]
        }else if(sum(pvals[[j]][1])==0){
          pvals<-pvals[-j]
        }
      }
      pvals<-do.call("rbind", pvals)
      colnames(pvals)<-c("WLR", "IUT", "S->M", "M->Y")
    }
    
    if(zeta==0){
      left_data_size <- dim(pvals)[1]
      tem<-colSums(pvals<0.2)/left_data_size*100
      Type1error[count.m,1]<-tem[1]
      Type1error[count.m+length(sample_size_list),1]<-tem[2]
      tem<-colSums(pvals<0.1)/left_data_size*100
      Type1error[count.m,2]<-tem[1]
      Type1error[count.m+length(sample_size_list),2]<-tem[2]
      tem<-colSums(pvals<0.05)/left_data_size*100
      Type1error[count.m,3]<-tem[1]
      Type1error[count.m+length(sample_size_list),3]<-tem[2]
      tem<-colSums(pvals<0.01)/left_data_size*100
      Type1error[count.m,4]<-tem[1]
      Type1error[count.m+length(sample_size_list),4]<-tem[2]
      tem<-colSums(pvals<0.005)/left_data_size*100
      Type1error[count.m,5]<-tem[1]
      Type1error[count.m+length(sample_size_list),5]<-tem[2]
      
      png(paste(data_location,'qqplot', sample_size,'.png', sep = ''),width=496, height=496)
      cut_off_list = (1:left_data_size)/left_data_size
      plot(-log10(cut_off_list), -log10(sort(pvals[, 1])), xlim = c(0, 4), ylim = c(0, 4), col = 'blue', cex = .5, xlab = 'Expected -log(p)', ylab = 'Observed -log(p)', main = paste('m = ', sample_size, sep = ''))
      points(-log10(cut_off_list), -log10(sort(pvals[, 2])), xlim = c(0, 4), ylim = c(0, 4), col = 'red', cex = .5)
      abline(0, 1, col = 'grey')
      dev.off()
      
      png(paste(data_location,'hist', sample_size,'.png', sep = ''),width=496, height=496)
      hist(pvals[, 1],xlab="p-value,WLR(b),IUT(r)",breaks = seq(0,1, by =0.1),ylim = c(0, repeat_size/14),ylab="Frequency",main = paste('m = ', sample_size, sep = ''),col = rgb(153, 204, 255, 100, maxColorValue = 255)) 
      hist(pvals[, 2],breaks = seq(0,1, by =0.1),col = rgb(255, 204, 229, 100, maxColorValue = 255),add=TRUE) 
      dev.off()
      
      pvals.rec<-c(pvals.rec,list(pvals))
    }
    
    tem<-colSums(pvals<0.05)/left_data_size*100
    power[1,count.zeta]<-tem[1]
    power[2,count.zeta]<-tem[2]
  }
  
  png(paste(data_location,'powerplot', sample_size,'.png', sep = ''),width=496, height=496)
  plot(zeta_list, power[1,], type = "b",cex = .5,xlim = c(0, max(zeta_list)), ylim = c(0, 100),col = 'blue',xlab='zeta',ylab='% of p<0.05',main = paste('m = ', sample_size, sep = ''))
  lines(zeta_list, power[2,],xlim = c(0, max(zeta_list)), ylim = c(0, 100),col='red',type="b",cex = .5)
  dev.off()
  
  power.rec<-c(power.rec,list(power))
}
save(Type1error,file=paste(data_location,'Type1error.rdata', sep = ''))
Type1error_table(Type1error)
save(pvals.rec,file=paste(data_location,'pvals_rec.rdata', sep = ''))
save(power.rec,file=paste(data_location,'power_rec.rdata', sep = ''))

if(2<1){
  sample_size_list = c(100,200,500,1000)
  for(i in 1:4){
    pvals<-pvals.rec[[i]]
    sample_size<-sample_size_list[i]
    png(paste('hist', sample_size,'.png', sep = ''),width=496, height=496)
    hist(pvals[, 1],xlab="p-value,WLR(b),IUT(r)",ylim = c(0, 10000/6),breaks = seq(0,1,0.1),ylab="Frequency",main = paste('m = ', sample_size, sep = ''),col = rgb(153, 204, 255, 100, maxColorValue = 255)) 
    hist(pvals[, 2],breaks = seq(0,1,0.1),col = rgb(255, 204, 229, 100, maxColorValue = 255),add=TRUE) 
    dev.off()
  }
}



