
print_with_time<-function(message, pad=100){
  #prints given message, leaving given buffer, and appends current time
  tms = max(pad-nchar(message), 3)
  cat(paste(message,paste(rep(".",tms),collapse=""),Sys.time(),"\n"))
}

get_platform_from_device<-function(y){
  x = tolower(y)
  ifelse((x=='doubledown' | x=='spa' | x=='webtoolbar'), 'Web', 
         ifelse(x %in% c("android","android_tablet"), "Android",
                ifelse(x %in% c("ipad","iphone"),"iOS",
                       ifelse(x %in% c("chromebar","mozbar"), "Toolbars",
                              ifelse(x=="metro",'Metro',"Other")))))    
}

get_type_from_device<-function(y){
  x = tolower(y)
  ifelse((x=='doubledown' | x=='spa' | x=='webtoolbar'), 'Web', 
         ifelse(x %in% c("ipad","android_tablet"), "Tablet",
                ifelse(x %in% c("android","iphone"),"Phone",
                       ifelse(x %in% c("chromebar","mozbar"), "Toolbars",
                              ifelse(x=="metro",'Metro',"Other")))))    
}

ncdg<-function(real, guess){
  #net cumulative discounted gain. real is a column vector, guess is a matricx (a single row matches a column element in real)
  print_with_time("Calculating score ...")
  n = length(real)
  scores = c()
  for (index in seq(1,n,1)){
    scores[index] = ncdg_row(real[index],guess[index,])
  }
  return(scores)
}


ncdg_row<-function(real, guess){   
  #real is a number, guess is a vector
  #DCG_k=∑i=(1..k) [2^rel_i−1] / log2(i+1),
  
  #nDCGk=DCG_k / IDCG_k,
  
  #where rel_i is the relevance of the result at position i (in this function, relevance is 0 or 1).
  
  #IDCG_k is the maximum possible (ideal) DCGDCG for a given set of queries. 
  #In this function, IDCG_k is 1
  #All NDCG calculations are relative values on the interval 0.0 to 1.0.
  
  m = match(real, guess)
  return(ifelse(is.na(m), 0, 1/log2(1+m)))
}


balance_nd<-function(data, field, weights=c()){
  #n-dimensional balancing, if weights not specified uses equal weight for all. 
  #else, users weights specified in weights. weights vector need to carry same names as factor levels of "field". 
  require(R.oo)#package for error handling
  require(plyr)#for count function
  n = dim(data)[1]
  print_with_time(sprintf("Total number of records: %.0f", n))
  lvls = unique(data[[field]])
  nl = length(lvls)
  print_with_time(sprintf("Target variable %s has %.0f levels", field, nl))  
  counts = count(data,field)
  counts$native_ratio = counts$freq/n
  counts$native_perc = sprintf("%.0f%%", counts$freq/n*100)
  
  if (length(weights)==0){
    #default: all levels have same ratio
    weights = rep(1, nl)
    names(weights)<-counts[[field]]
  } else {
    #normalize such that level with minimum positive weight has weight 1
    min_weight = min(weights[weights>0])
    if (is.null(names(weights))){
      throw ("weights needs to be a named vector")
    }
    if (any(weights<0)){
      throw("no weight can be negative.")
    }
    #now all weights are positive
    if (min_weight==0){
      print_with_time("at least one weight needs to be positive, resetting to uniform weights")
      weights = rep(1, nl)    
    } else if (min_weight!=1){
      #at least one weight positive, normalize since min_weight is not one
      print_with_time("normalizing such that min(positive weights) = 1. Original:")
      print(weights)
      weights = weights/min_weight
      print_with_time("new weights:")
      print(weights)
    }
  }
  counts$weight = 0
  for (index in seq(1, nl, 1)){
    lvl = counts[index, field]
    counts[index, "weight"] = ifelse(lvl %in% names(weights), weights[as.character(lvl)], 0)
  }
  counts$desired_ratio = counts$weight/sum(counts$weight)
  counts$bottlenecks = counts$freq/counts$desired_ratio
  bottleneck = floor(min(counts$bottlenecks)) #max possible total supported
  counts$num_instances = floor(bottleneck * counts$desired_ratio)
  
  data$dummy = 1
  sampled_data = subset(data, dummy==0)#initialize to empty data frame
  for (lvl in lvls){
    sub = data[data[[field]]==lvl,]
    num_instances= counts[counts[[field]]==lvl,"num_instances"]
    sampled_data = rbind(sampled_data, sample_df(sub, num_instances))
  }
  counts_sd = count(sampled_data, field)
  m = sum(counts_sd$freq)
  counts_sd$sampled_ratio = counts_sd$freq/m
  counts_sd = counts_sd[,-2]#remove freq
  counts = merge(counts,counts_sd, all.x = TRUE)
  print_with_time("Set up complete. Native prevalence and new ratios:")
  print(counts)
  return(sampled_data)
}


get_device<-function(x){
  devices<-c("NONE","MOZBAR","IEBAR","CHROMEBAR","FFLITEBAR","SAFARIBAR","WEBTOOLBAR",
             "STUMBLEVIDEO","WEBTB_MOBILE","API","IPHONE","ANDROID","IPAD","IPHONE_INFO","ANDROID_INFO","ANDROID_TABLET")
  names(devices)<-c("0","1","2","3","4","5","10","11","12","30","31","32","33","34","35","36")
  return(devices[as.character(x)])
}

sort_by<-function(df,field, is_decreasing=FALSE){
  #sorts a given data frame by a certain field
  #sort by multiple columns: dd[with(dd, order(-z, b)), ] where z and b are column names and minus indicates descending
  return(df[order(df[[field]], decreasing = is_decreasing),])
}

get_column_classes<-function(df){
  return(lapply(df, class))
}

impute_fields<-function(data_to_impute, exclude=c("userid")){
  # Fill in missing values. 
  # For numerical, use median
  # For factor, use most frequent.
  # For character, use empty string
  # Do not perform imputation for fields in excluded list
  original_num_records = dim(data_to_impute)[1]
  complete_num_records = dim(data_to_impute[complete.cases(data_to_impute),])[1]
  fields_to_inpute = names(which(sapply(data_to_impute,function(x) any(is.na(x)))))
  if (length(fields_to_inpute) > 0){
    print_with_time(paste((original_num_records-complete_num_records),
                          " records from ",original_num_records,"have missing values. Imputing columns :",paste(fields_to_inpute, collapse=",")))
  } else {
    print_with_time("Data does not have missing values. Returning itself.")
    return(data_to_impute)
  }
  
  for (field in fields_to_inpute){
    if (class(data_to_impute[[field]])=="numeric"){
      value = floor(median(data_to_impute[[field]], na.rm=TRUE)) #median
    } else if (class(data_to_impute[[field]])=="factor"){ #categorical
      value = names(summary(data_to_impute[[field]]))[1] #most frequent
    } else {#textual
      value = ""
    }
    data_to_impute[which(is.na(data_to_impute[[field]])),field] = value
  }
  
  final_complete_num_records = dim(data_to_impute[complete.cases(data_to_impute),])[1]
  if (original_num_records!=final_complete_num_records){
    stop(paste("Encountered error in imputing values. Original number of records : ",original_num_records,", final complete cases : ",final_complete_num_records))
  }
  print_with_time("Imputation completed")
  return(data_to_impute)
}

sample_df<-function(df, sample_size){
  #samples sample_size rows from data frame randomly
  print_with_time(paste("Taking a sample of size ",sample_size," from ",dim(df)[1]," records"))
  indices = c(1:length(df[,1]))
  sample_indices = sample(indices, sample_size, replace=FALSE)
  return(df[sample_indices, ])
}

read_delimited_data<-function(filename,
                              delimiter, #alternatives: ",", \t",...
                              has_header=FALSE, #if T, column names and types are ignored
                              colNames,#comma-delimited fields
                              colClasses #comma-delimited field types, such as numeric, factor, character
){
  #reads a delimited file. example call: demo_read_delimited_data()
  col_names = c(unlist(strsplit(colNames,",")))
  col_classes = c(unlist(strsplit(colClasses,",")))
  
  print_with_time(sprintf("Reading data from: %s",filename))
  if (has_header){
    data<-read.delim(filename, header = T, sep=delimiter) 
  } else {
    data<-read.delim(filename, header = F, sep=delimiter, col.names = col_names, colClasses = col_classes) 
  }
  num_rows = dim(data)[1]
  print_with_time(sprintf("Read %.0f records successfully. Top ten lines: ",num_rows))
  print(head(data, 10))
  return(data)
}

demo_read_delimited_data<-function(){
  #http://www.ssa.gov/oact/babynames/state/namesbystate.zip
  #http://www.cyclismo.org/tutorial/R/input.html#reading-a-csv-file
  
    read_delimited_data("/Users/oeskiyenenturk/Dropbox/Documents/Personal/Oytun/Programming/data/CA.TXT",
                        ",",
                        FALSE,
                        "state,gender,year,name,count",
                        "factor,factor,numeric,character,numeric")

}

pg_by_day<-function(ldata, input, output="/Users/oeskiyenenturk/su/datascience/PER-676/pg_by_day3.tsv", app=F){
  write.table(input, file=output, quote=F, sep="\t", row.names=F, col.names=F, append=app)
  agg = aggregate(cbind(stumbles,greats,mstumbles,mgreats) ~ date, ldata, sum)
  write.table(agg, file=output, quote=F, sep='\t', row.names=F, col.names=T, append=T)
}

pg_by_day_granular<-function(ldata, output, app=F, h=F){
  agg = aggregate(cbind(stumbles,greats,mstumbles,mgreats) ~ test + group + date + device + device_name + device_type + method, ldata, sum)
  write.table(agg, file=output, quote=F, sep='\t', row.names=F, col.names=h, append=app)
}

pg_by_user<-function(ldata_d, ldata_i, ldata_e, input=c('default3','implicit3','explicit3'), output="/Users/oeskiyenenturk/su/datascience/PER-676/pg_by_user3.tsv", app=F){
  agg_d = aggregate(cbind(mstumbles,mgreats) ~ userid, ldata_d, sum)
  agg_i = aggregate(cbind(mstumbles,mgreats) ~ userid, ldata_i, sum)
  agg_e = aggregate(cbind(mstumbles,mgreats) ~ userid, ldata_e, sum)
  agg_d$pg = ifelse(agg_d$mstumbles==0, 0, agg_d$mgreats/agg_d$mstumbles)
  agg_i$pg = ifelse(agg_i$mstumbles==0, 0, agg_i$mgreats/agg_i$mstumbles)
  agg_e$pg = ifelse(agg_e$mstumbles==0, 0, agg_e$mgreats/agg_e$mstumbles)
  num_quantiles=10
  qt = quantile(agg_d$mstumbles, seq((1/num_quantiles),1,(1/num_quantiles)))
  
  write_couple(output, input[1],input[2],agg_d,agg_i,"mstumbles","pg",F,10,qt) 
  write_couple(output, input[1],input[3],agg_d,agg_e,"mstumbles","pg",T,10,qt) 
  write_couple(output, input[2],input[3],agg_i,agg_e,"mstumbles","pg",T,10,qt) 
  
}

retention<-function(sdata_d, sdata_i, sdata_e,input=c('default3','implicit3','explicit3'), output, app=F, threshold=13){
  print_with_time(sprintf("Extracting retention for %s",input[1]))
  agg_d = extract_retention(sdata_d, threshold)
  print_with_time(sprintf("Extracting retention for %s",input[2]))
  agg_i = extract_retention(sdata_i, threshold)
  print_with_time(sprintf("Extracting retention for %s",input[3]))
  agg_e = extract_retention(sdata_e, threshold)
  
  names(agg_d)<-c("userid","mstumbles","retained","retained2")
  names(agg_i)<-c("userid","mstumbles","retained","retained2")
  names(agg_e)<-c("userid","mstumbles","retained","retained2")
  
  num_quantiles=10
  qt = quantile(agg_d$mstumbles, seq((1/num_quantiles),1,(1/num_quantiles)))
  
  write_couple(output, input[1],input[2],agg_d,agg_i,"mstumbles","retained",F,10,qt) 
  write_couple(output, input[1],input[2],agg_d,agg_i,"mstumbles","retained2",T,10,qt) 
  write_couple(output, input[1],input[3],agg_d,agg_e,"mstumbles","retained",T,10,qt) 
  write_couple(output, input[1],input[3],agg_d,agg_e,"mstumbles","retained2",T,10,qt) 
  write_couple(output, input[2],input[3],agg_i,agg_e,"mstumbles","retained",T,10,qt) 
  write_couple(output, input[2],input[3],agg_i,agg_e,"mstumbles","retained2",T,10,qt)   
}

write_couple<-function(output, d,i,agg_d,agg_i,quantile_field,metric_field, app,num_quantiles=10, qt){
  print_with_time(sprintf("%s vs. %s", d, i))
  write.table(sprintf("%s vs. %s", d, i), file=output, quote=F, sep="\t", row.names=F, col.names=F, append=app)
  columns = sprintf("\tsize(%s)\tmean(%s$%s)\tsize(%s)\tmean(%s$%s)\tt.test( mean(%s$%s) > mean(%s$%s) )",d,d,metric_field,i,i,metric_field,i,metric_field,d,metric_field)
  td = t.test(agg_i[[metric_field]], agg_d[[metric_field]], alternative="greater")
  values = sprintf("Overall\t%.0f\t%f\t%.0f\t%f\t%f", dim(agg_d)[1], mean(agg_d[[metric_field]]),dim(agg_i)[1], mean(agg_i[[metric_field]]), td$p.value)
  write.table(columns, file=output, quote=F, sep='\t', row.names=F, col.names=F, append=T)
  write.table(values, file=output, quote=F, sep='\t', row.names=F, col.names=F, append=T)
  #loop over quantiles
  
  for (index in seq(1,num_quantiles-1,1)){
    #higher = qt[index]
    #lower= ifelse(index==1, -1, qt[index-1])
    lower = qt[index]
    range = sprintf("Top %0.f %%", (num_quantiles-index)/num_quantiles*100)
    print_with_time(sprintf("Processing %s", range))    
    #label = sprintf("Quantile %.0f (%s in (%.0f,%.0f])", index, quantile_field, lower, higher)  
    label = sprintf("%s (%s > %.0f)", range, quantile_field, lower)  
    
    #agg_dp = agg_d[which(agg_d[[quantile_field]]>lower & agg_d[[quantile_field]]<=higher),] 
    #agg_ip = agg_i[which(agg_i[[quantile_field]]>lower & agg_i[[quantile_field]]<=higher),]
    agg_dp = agg_d[which(agg_d[[quantile_field]]>lower),] 
    agg_ip = agg_i[which(agg_i[[quantile_field]]>lower),]
    td_p = t.test(agg_ip[[metric_field]], agg_dp[[metric_field]], alternative="greater")
    values_p = sprintf("%s\t%.0f\t%f\t%.0f\t%f\t%f", label, dim(agg_dp)[1], mean(agg_dp[[metric_field]]),dim(agg_ip)[1], mean(agg_ip[[metric_field]]), td_p$p.value)
    write.table(values_p, file=output, quote=F, sep='\t', row.names=F, col.names=F, append=T)    
  } 
  print_with_time(sprintf("%s vs. %s completed", d, i))
}

extract_retention<-function(sdata_d, threshold=13){
  rdata_d = subset(sdata_d, requested>=join_timestamp & maxAttainableDayIndex>=threshold & dayIndex<=threshold)
  rdata_d$modified_mstumbles = ifelse(rdata_d$dayIndex<7,rdata_d$mstumbles,0)
  rdata_d$minDayIndex = ifelse(rdata_d$dayIndex %in% seq(7,threshold,1),rdata_d$dayIndex,threshold+1)
  rdata_d$maxDayIndex = ifelse(rdata_d$dayIndex %in% seq(7,threshold,1),rdata_d$dayIndex,0) 
  agg_d_mstumbles = aggregate(modified_mstumbles ~ userid, rdata_d, sum)
  agg_d_min = aggregate(minDayIndex ~ userid, rdata_d, min)
  agg_d_max = aggregate(maxDayIndex ~ userid, rdata_d, max)
  agg_d = merge(agg_d_mstumbles,agg_d_min)
  agg_d = merge(agg_d, agg_d_max)
  agg_d$retained = ifelse(agg_d$maxDayIndex==0,0,1)
  agg_d$retained2 = ifelse(agg_d$maxDayIndex>agg_d$minDayIndex,1,0)
  return(agg_d[,c("userid","modified_mstumbles","retained","retained2")])
}

extract_test_data<-function(fn, group, test ){
  cn = paste('user_id,join_timestamp,dayIndex,maxAttainableDayIndex,userid,device,method,topic,urlid,',
             'stumbles,incat_stumbles,mstumbles,requested,date,hour,greats,incat_greats,mgreats,bads,incat_bads,mbads,session_num,stumble_in_session', sep='')
  columnNames = unlist(strsplit(cn,','))
  df = data.frame(name=columnNames,type='numeric',stringsAsFactors=F)
  df[which(df$name=='date'),2] = 'character'
  data = read_delimited_data(filename=fn,
                                delimiter='\t', #alternatives: ",", \t",...
                                has_header=FALSE, #if T, column names and types are ignored
                                colNames=paste(df$name,collapse=','),      
                                colClasses=paste(df$type,collapse=',') #4,599,164(default), 4,574,690(implicit)
  )
  data$device_name = get_device(data$device)
  data$device_type = get_type_from_device(data$device_name)
  data$group = group
  data$test = test
  
  print_with_time(sprintf("Original # of stumbles from metrics: %d",dim(data)[1]))
  original_users = unique(data$user_id)
  print_with_time(sprintf("Original # of users: %d",length(original_users)))
  
  #users that stumbled at least once:
  sdata = subset(data, !is.na(join_timestamp))#4,522,084; 4,487,875
  print_with_time(sprintf("# of stumbles from users that stumbled at least once during test period: %d",dim(sdata)[1]))
  stumblers = unique(sdata$user_id)
  print_with_time(sprintf("# of users that stumbled at least once during test period: %d",length(stumblers)))
  
  #stumbles that happen within one week, in certain devices only
  ldata = subset(sdata, dayIndex<7 & requested>=join_timestamp & requested <= 1451127864 & device %in% c(1, 3, 10, 31, 32, 33, 36, 60))
  print_with_time(sprintf("# of stumbles from stumblers with correct device, within 7 days, after join time: %d",dim(ldata)[1]))
  final_users = unique(ldata$user_id)
  print_with_time(sprintf("# of users that stumbled at least once during test period with correct device, within 7 days, after join time: %d",length(final_users)))
  return(list("sdata"=sdata,"ldata"=ldata))
}

engagement<-function(sdata_d, sdata_i, sdata_e, output="/Users/oeskiyenenturk/su/datascience/PER-676/engagement3.tsv", app=F){
  #what we need:
  #median stumbles and seconds over all sessions, compared to two other groups
  #median stumbles and seconds over test sessions, compared to two other groups
  #median stumbles and seconds over after-test sessions, compared to two other groups
  #for each person: mean stumbles, seconds, median stumbles, seconds (overall, during and after the test). compare the means among three groups
  #finally, compare engagement on a daily level, among three groups. instead of date itself, use dayIndex
  d = sdata_d[1,"group"]
  print_with_time(sprintf("Processing %s", d))
  edata_d = extract_engagement_data(sdata_d)
  i = sdata_i[1,"group"]
  print_with_time(sprintf("Processing %s", i))
  edata_i = extract_engagement_data(sdata_i)
  e = sdata_e[1,"group"]
  print_with_time(sprintf("Processing %s", e))
  edata_e = extract_engagement_data(sdata_e)
  
  #median stumbles and seconds over all sessions, compared to two other groups
  print_with_time("Processing overall stumbles")
  num_quantiles = 10
  qt = seq(1/num_quantiles,1,1/num_quantiles)
  
  print_with_time("Processing overall stumbles")
  prefix = "overall_stumbles"
  print_quantiles(edata_d,edata_i,edata_e, "stumbles",
                  c(sprintf("%s_%s",prefix,d), sprintf("%s_%s",prefix,i),sprintf("%s_%s",prefix,e)),
                  output,qt,F)
  
  print_with_time("Processing overall seconds")
  prefix = "overall_seconds"
  print_quantiles(edata_d,edata_i,edata_e, "seconds",
                  c(sprintf("%s_%s",prefix,d), sprintf("%s_%s",prefix,i),sprintf("%s_%s",prefix,e)),
                  output,qt,T)
  
  print_with_time("Processing test stumbles")
  prefix = "test_stumbles"
  print_quantiles(subset(edata_d, dayIndex<7 & mstumbles>0),subset(edata_i, dayIndex<7 & mstumbles>0),subset(edata_e,dayIndex<7 & mstumbles>0), "stumbles",
                  c(sprintf("%s_%s",prefix,d), sprintf("%s_%s",prefix,i),sprintf("%s_%s",prefix,e)),
                  output,qt,T)
  print_with_time("Processing test seconds")
  prefix = "test_seconds"
  print_quantiles(subset(edata_d, dayIndex<7 & mstumbles>0),subset(edata_i, dayIndex<7 & mstumbles>0),subset(edata_e,dayIndex<7 & mstumbles>0), "seconds",
                  c(sprintf("%s_%s",prefix,d), sprintf("%s_%s",prefix,i),sprintf("%s_%s",prefix,e)),
                  output,qt,T)
  
  print_with_time("Processing after-test stumbles")
  prefix = "after_test_stumbles"
  print_quantiles(subset(edata_d, dayIndex>=7),subset(edata_i, dayIndex>=7),subset(edata_e,dayIndex>=7), "stumbles",
                  c(sprintf("%s_%s",prefix,d), sprintf("%s_%s",prefix,i),sprintf("%s_%s",prefix,e)),
                  output,qt,T)
  print_with_time("Processing after-test seconds")
  prefix = "after_test_seconds"
  print_quantiles(subset(edata_d, dayIndex>=7),subset(edata_i, dayIndex>=7),subset(edata_e,dayIndex>=7), "seconds",
                  c(sprintf("%s_%s",prefix,d), sprintf("%s_%s",prefix,i),sprintf("%s_%s",prefix,e)),
                  output,qt,T)
  
  sessions_by_dayIndex(edata_d, d, output)
  sessions_by_dayIndex(edata_i, i, output)
  sessions_by_dayIndex(edata_e, e, output)
  
  #aggregate by userid - mean stumbles and seconds for each person, during test period. Segment by mstumbles during the session
  print_with_time("Processing by userid during test")
  engagement_per_user(edata_d, edata_i,d,i,output,"stumbles")
  engagement_per_user(edata_d, edata_i,d,i,output,"seconds")
  engagement_per_user(edata_d, edata_e,d,e,output,"stumbles")
  engagement_per_user(edata_d, edata_e,d,e,output,"seconds")
  engagement_per_user(edata_i, edata_e,i,e,output,"stumbles")
  engagement_per_user(edata_i, edata_e,i,e,output,"seconds")
  
  #this section prints out engagement by user with device

  
  for (device in c("ANDROID","ANDROID_TABLET","CHROMEBAR","IPAD","IPHONE","MOZBAR","WEBTOOLBAR")){
    print_with_time(sprintf("Processing by userid during test (%s)", device))
    output = sprintf("/Users/oeskiyenenturk/su/datascience/PER-676/engagement34_peruser_%s.tsv", tolower(device))
    engagement_per_user(subset(edata_d, top_device_name==device), subset(edata_i, top_device_name==device),
                        sprintf("%s_%s",d, device),sprintf("%s_%s",i, device),output,"stumbles")
    engagement_per_user(subset(edata_d, top_device_name==device), subset(edata_i, top_device_name==device),
                        sprintf("%s_%s",d, device),sprintf("%s_%s",i, device),output,"seconds")
    engagement_per_user(subset(edata_d, top_device_name==device), subset(edata_e, top_device_name==device),
                        sprintf("%s_%s",d, device),sprintf("%s_%s",e, device),output,"stumbles")
    engagement_per_user(subset(edata_d, top_device_name==device), subset(edata_e, top_device_name==device),
                        sprintf("%s_%s",d, device),sprintf("%s_%s",e, device),output,"seconds")
    engagement_per_user(subset(edata_i, top_device_name==device), subset(edata_e, top_device_name==device),
                        sprintf("%s_%s",i, device),sprintf("%s_%s",e, device),output,"stumbles")
    engagement_per_user(subset(edata_i, top_device_name==device), subset(edata_e, top_device_name==device),
                        sprintf("%s_%s",i, device),sprintf("%s_%s",e, device),output,"seconds")
  }
  
  
  #==this section plots engagement by dayIndex  (for sessions with mstumbles>0) at device level  
  output = "/Users/oeskiyenenturk/su/datascience/PER-676/engagement3_devices.tsv"
  
  for (device in c("ANDROID","ANDROID_TABLET","CHROMEBAR","IPAD","IPHONE","MOZBAR","WEBTOOLBAR")){
    sessions_by_dayIndex(subset(edata_d, device_name==device), sprintf("%s_%s", d, device), output)
    sessions_by_dayIndex(subset(edata_i, device_name==device), sprintf("%s_%s", i, device), output)
    sessions_by_dayIndex(subset(edata_e, device_name==device), sprintf("%s_%s", e, device), output)
  }
  
  
  #aggregate by userid - mean stumbles and seconds for each person, AFTER test period. Segment by mstumbles during the test period.
  print_with_time("Processing by userid after test")
  engagement_aftertest_per_user(edata_d, e_data_i, d, i, "stumbles", "mstumbles", output, 5)    
  engagement_aftertest_per_user(edata_d, e_data_i, d, i, "seconds", "mstumbles", output, 10)  
  engagement_aftertest_per_user(edata_d, e_data_e, d, e, "stumbles", "mstumbles", output, 5)    
  engagement_aftertest_per_user(edata_d, e_data_e, d, e, "seconds", "mstumbles", output, 10)
  engagement_aftertest_per_user(edata_i, e_data_e, i, e, "stumbles", "mstumbles", output, 5)    
  engagement_aftertest_per_user(edata_i, e_data_e, i, e, "seconds", "mstumbles", output, 10)
}

engagement_aftertest_per_user<-function(edata1, e_data2, d1, d2, metric_field, quantile_field, output, num_quantiles){
  agg_by_user1 = aggregate(mstumbles ~ userid, edata1, sum)
  adata1 = subset(edata1, dayIndex>=7)
  agg_mean1 = aggregate(as.formula(sprintf("%s ~ userid", metric_field)), adata1, mean)
  agg_by_user1 = merge(agg_by_user1,agg_mean1) #here, mstumbles are stumbles in test period, while stumbles are mean of stumbles_per_session in after-test
  
  agg_by_user2 = aggregate(mstumbles ~ userid, edata2, sum)
  adata2 = subset(edata2, dayIndex>=7)
  agg_mean2 = aggregate(as.formula(sprintf("%s ~ userid", metric_field)), adata2, mean)
  agg_by_user2 = merge(agg_by_user2,agg_mean2) #here, mstumbles are stumbles in test period, while stumbles are mean of stumbles_per_session in after-test
  
  qt = quantile(agg_by_user1[[quantile_field]], seq(1/num_quantiles, 1, 1/num_quantiles)) 
  write_couple(output, d1,d2,agg_by_user1,agg_by_user2,"mstumbles",metric_field, T,num_quantiles, qt)
}

engagement_per_user<-function(edata_d, edata_i,d,i,output,metric_field){
  print_with_time(sprintf("T.test stumbles difference for test stumbles for %s vs. %s", d,i))
  qt = c(seq(1,6,1), seq(8,12,2),seq(16,32,4),seq(42,72,10)) # mstumbles
  #default vs. implicit
  #metric_field = "seconds"
  quantile_field = "mstumbles"
  
  columns = sprintf("\tsize(%s)\tmean(%s$%s)\tsize(%s)\tmean(%s$%s)\tt.test( mean(%s$%s) > mean(%s$%s) )",d,d,metric_field,i,i,metric_field,i,metric_field,d,metric_field)
  write.table(columns, file=output, quote=F, sep='\t', row.names=F, col.names=F, append=T)
  total = length(unique(subset(edata_d, dayIndex<7)[["userid"]]))
  for (index in seq(1,length(qt),1)){
    lower = ifelse(index==1, -1, qt[index-1])
    tdata_d = subset(edata_d, dayIndex<7 & mstumbles>lower)
    agg_mean_d = aggregate(as.formula(sprintf("%s ~ userid", metric_field)), tdata_d, mean)   
    tdata_i = subset(edata_i, dayIndex<7 & mstumbles>lower)
    agg_mean_i = aggregate(as.formula(sprintf("%s ~ userid", metric_field)), tdata_i, mean)
    perc = dim(agg_mean_d)[1]/total*100
    range = ifelse(index==1,"Overall",sprintf("Top %0.f %%", perc))
    print_with_time(sprintf("Processing %s", range)) 
    
    td_p = t.test(agg_mean_i[[metric_field]], agg_mean_d[[metric_field]], alternative="greater")
    label = sprintf("%s (%s > %.0f)", range, quantile_field, lower)  
    values_p = sprintf("%s\t%.0f\t%f\t%.0f\t%f\t%f", label, dim(agg_mean_d)[1], mean(agg_mean_d[[metric_field]]),dim(agg_mean_i)[1], mean(agg_mean_i[[metric_field]]), td_p$p.value)
    
    write.table(values_p, file=output, quote=F, sep='\t', row.names=F, col.names=F, append=T)
  }
}

engagement_by_device<-function(edata_d, edata_e, edata_i, input = c("default34", "implicit34", "explicit34"), plot_dir = "/Users/oeskiyenenturk/su/datascience/PER-676/plots"){
  print_with_time(sprintf("Plotting by-day-index during test period for %s", category))
  d = input[1]
  i = input[2]
  e = input[3]
  #we look at only sessions that have at least one test stumble
  for (device in c("ANDROID","ANDROID_TABLET","CHROMEBAR","IPAD","IPHONE","MOZBAR","WEBTOOLBAR")){
    seconds = data.frame(x = seq(0,6,1))
    seconds[[d]] = 0
    seconds[[i]] = 0
    seconds[[e]] = 0
    stumbles = seconds
    
    for (di in seconds$x){
      seconds[di+1,d] = quantile(subset(edata_d, dayIndex==di & mstumbles>0 & device_name==device)[["seconds"]], 0.5)
      stumbles[di+1,d] = quantile(subset(edata_d, dayIndex==di & mstumbles>0 & device_name==device)[["stumbles"]], 0.5)    
      seconds[di+1,i] = quantile(subset(edata_i, dayIndex==di & mstumbles>0 & device_name==device)[["seconds"]], 0.5)
      stumbles[di+1,i] = quantile(subset(edata_i, dayIndex==di & mstumbles>0 & device_name==device)[["stumbles"]], 0.5) 
      seconds[di+1,e] = quantile(subset(edata_e, dayIndex==di & mstumbles>0 & device_name==device)[["seconds"]], 0.5)
      stumbles[di+1,e] = quantile(subset(edata_e, dayIndex==di & mstumbles>0 & device_name==device)[["stumbles"]], 0.5) 
    }
    
    draw(seconds, measure = "seconds", output = file.path(plot_dir, 
                                                          sprintf("%s_test_seconds.png", tolower(device))), device_name=device)
    draw(stumbles, measure = "stumbles", output = file.path(plot_dir, 
                                                            sprintf("%s_test_stumbles.png", tolower(device))), device_name=device)
    
  }
  edata_d$platform  = get_platform_from_device(edata_d$device_name)
  edata_d$type  = get_type_from_device(edata_d$device_name)
  edata_i$platform  = get_platform_from_device(edata_i$device_name)
  edata_i$type  = get_type_from_device(edata_i$device_name)
  edata_e$platform  = get_platform_from_device(edata_e$device_name)
  edata_e$type  = get_type_from_device(edata_e$device_name)
  
  
  for (device in c("ANDROID","ANDROID_TABLET","CHROMEBAR","IPAD","IPHONE","MOZBAR","WEBTOOLBAR")){
    seconds = data.frame(x = seq(0,6,1))
    seconds[[d]] = 0
    seconds[[i]] = 0
    seconds[[e]] = 0
    stumbles = seconds
    
    for (di in seconds$x){
      seconds[di+1,d] = quantile(subset(edata_d, dayIndex==di & mstumbles>0 & device_name==device)[["seconds"]], 0.5)
      stumbles[di+1,d] = quantile(subset(edata_d, dayIndex==di & mstumbles>0 & device_name==device)[["stumbles"]], 0.5)    
      seconds[di+1,i] = quantile(subset(edata_i, dayIndex==di & mstumbles>0 & device_name==device)[["seconds"]], 0.5)
      stumbles[di+1,i] = quantile(subset(edata_i, dayIndex==di & mstumbles>0 & device_name==device)[["stumbles"]], 0.5) 
      seconds[di+1,e] = quantile(subset(edata_e, dayIndex==di & mstumbles>0 & device_name==device)[["seconds"]], 0.5)
      stumbles[di+1,e] = quantile(subset(edata_e, dayIndex==di & mstumbles>0 & device_name==device)[["stumbles"]], 0.5) 
    }
    
    draw(seconds, measure = "seconds", output = file.path(plot_dir, 
                                                          sprintf("%s_test_seconds.png", tolower(device))), device_name=device)
    draw(stumbles, measure = "stumbles", output = file.path(plot_dir, 
                                                            sprintf("%s_test_stumbles.png", tolower(device))), device_name=device)
    
  }
  
  for (device_platform in c("Toolbars","Web","Android","iOS")){
    seconds = data.frame(x = seq(0,6,1))
    seconds[[d]] = 0
    seconds[[i]] = 0
    seconds[[e]] = 0
    stumbles = seconds
    
    for (di in seconds$x){
      seconds[di+1,d] = quantile(subset(edata_d, dayIndex==di & mstumbles>0 & platform==device_platform)[["seconds"]], 0.5)
      stumbles[di+1,d] = quantile(subset(edata_d, dayIndex==di & mstumbles>0 & platform==device_platform)[["stumbles"]], 0.5)    
      seconds[di+1,i] = quantile(subset(edata_i, dayIndex==di & mstumbles>0 & platform==device_platform)[["seconds"]], 0.5)
      stumbles[di+1,i] = quantile(subset(edata_i, dayIndex==di & mstumbles>0 & platform==device_platform)[["stumbles"]], 0.5) 
      seconds[di+1,e] = quantile(subset(edata_e, dayIndex==di & mstumbles>0 & platform==device_platform)[["seconds"]], 0.5)
      stumbles[di+1,e] = quantile(subset(edata_e, dayIndex==di & mstumbles>0 & platform==device_platform)[["stumbles"]], 0.5) 
    }
    
    draw(seconds, measure = "seconds", output = file.path(plot_dir, 
                                                          sprintf("%s_test_seconds.png", tolower(device_platform))), device_name=device_platform)
    draw(stumbles, measure = "stumbles", output = file.path(plot_dir, 
                                                            sprintf("%s_test_stumbles.png", tolower(device_platform))), device_name=device_platform)
    
  }
  
  for (device_type in c( "Phone","Tablet")){
    seconds = data.frame(x = seq(0,6,1))
    seconds[[d]] = 0
    seconds[[i]] = 0
    seconds[[e]] = 0
    stumbles = seconds
    
    for (di in seconds$x){
      seconds[di+1,d] = quantile(subset(edata_d, dayIndex==di & mstumbles>0 & type==device_type)[["seconds"]], 0.5)
      stumbles[di+1,d] = quantile(subset(edata_d, dayIndex==di & mstumbles>0 & type==device_type)[["stumbles"]], 0.5)    
      seconds[di+1,i] = quantile(subset(edata_i, dayIndex==di & mstumbles>0 & type==device_type)[["seconds"]], 0.5)
      stumbles[di+1,i] = quantile(subset(edata_i, dayIndex==di & mstumbles>0 & type==device_type)[["stumbles"]], 0.5) 
      seconds[di+1,e] = quantile(subset(edata_e, dayIndex==di & mstumbles>0 & type==device_type)[["seconds"]], 0.5)
      stumbles[di+1,e] = quantile(subset(edata_e, dayIndex==di & mstumbles>0 & type==device_type)[["stumbles"]], 0.5) 
    }
    
    draw(seconds, measure = "seconds", output = file.path(plot_dir, 
                                                          sprintf("%s_test_seconds.png", tolower(device_type))), device_name=device_type)
    draw(stumbles, measure = "stumbles", output = file.path(plot_dir, 
                                                            sprintf("%s_test_stumbles.png", tolower(device_type))), device_name=device_type)
    
  }
  
}

draw<-function(data, output, measure="stumbles", device_name){
  myplot = multiple_lines_charts(data,                         
                        ticker_fill = c("white","#A1B5ED","#F1A5A6"),  #could be specified as a named vector, fields same as df
                        line_colors=c("black","#1F497D","#800000"),
                        xlabel = "dayIndex",
                        ylabel = sprintf("session length in %s", measure),
                        title = sprintf("session length in %s by dayIndex (%s)", measure, device_name),
                        legend_position = c(0.9, 0.2),
                        title_font_type = 'Calibri',
                        title_font_weight = "bold",
                        title_font_size = 14,
                        legend_font_type = 'Calibri',
                        legend_font_weight = "plain",
                        legend_font_size = 14,
                        x_axis_interval = 1
  )
  print(myplot)
  dev.copy(png,output,width=800,height=600)
  dev.off()
}

sessions_by_dayIndex<-function(edata_d, d, output){
  print_with_time(sprintf("Processing by-day-index during test period for %s", d))
  #we look at only sessions that have at least one test stumble
  di_test = data.frame(dayIndex = seq(0,6,1), size=0, median_stumbles=0, median_seconds=0, size_qualified=0, median_stumbles_qualified=0, median_seconds_qualified=0)
  for (di in di_test$dayIndex){
    data = subset(edata_d, dayIndex==di)
    qdata = subset(data, mstumbles>0)
    di_test[di+1,] = c(di, dim(data)[1],quantile(data$stumbles, 0.5), quantile(data$seconds, 0.5), dim(qdata)[1], quantile(qdata$stumbles, 0.5), quantile(qdata$seconds, 0.5))
  }
  names(di_test)<-paste(names(di_test),d,sep="_") 
  write.table(di_test, file=output, quote=F, sep='\t', row.names=F, col.names=T, append=T)
  

  
  require(plyr)
  di_after_test = subset(count(edata_d,"dayIndex"), dayIndex>=7)
  di_after_test$median_stumbles = 0
  di_after_test$median_seconds = 0  
  for (row in seq(1,dim(di_after_test)[1],1)){
    di = di_after_test[row,"dayIndex"]
    data = subset(edata_d, dayIndex==di)
    di_after_test[row,] = c(di, dim(data)[1],quantile(data$stumbles, 0.5), quantile(data$seconds, 0.5))    
  }
  names(di_after_test)<-paste(names(di_after_test),d,sep="_") 
  write.table(di_after_test, file=output, quote=F, sep='\t', row.names=F, col.names=T, append=T)
  
}

align_factor_levels<-function(train_subset, test_subset){
  #align levels
  column_classes = get_column_classes(train_subset)
  column_names = names(train_subset)
  factor_columns = column_names[column_classes=='factor']
  factor_columns = factor_columns[!(factor_columns %in% c("id","country_destination"))]
  for (column in factor_columns){
    test_levels = levels(test_subset[[column]])
    train_levels = levels(train_subset[[column]])
    missing_levels = setdiff(test_levels, train_levels)
    if (length(missing_levels)>0){
      print_with_time(sprintf("Creating unknown level for column %s", column))
      new_train_levels = c(train_levels,"-unknown-")
      new_test_levels = c(test_levels,"-unknown-")     
      levels(train_subset[[column]])<-new_train_levels
      levels(test_subset[[column]])<-new_test_levels    
      test_subset[!(test_subset[[column]] %in% train_levels), column] = "-unknown-" 
      test_subset[[column]] <- factor(test_subset[[column]])
      levels(test_subset[[column]])<-new_train_levels          
    } else if (length(train_levels) > length(test_levels)){
      print_with_time(sprintf("Adding test levels for column %s", column))      
      levels(test_subset[[column]]) <- train_levels
    }
    
    if (length(levels(train_subset[[column]]))>10){
      print_with_time(sprintf("Field %s in training data has the too many levels", column))
    }
    
    if (length(levels(test_subset[[column]]))>10){
      print_with_time(sprintf("Field %s in test data has the too many levels", column))
    }
  }
  return(list("train"=train_subset,"test"=test_subset))  
}

print_quantiles<-function(edata_d,edata_i,edata_e,field,headers,output,qt,app){
  data<-data.frame(overall_stumbles_d = quantile(edata_d[[field]], qt),
                              overall_stumbles_i = quantile(edata_i[[field]], qt),
                              overall_stumbles_e = quantile(edata_e[[field]], qt))
  names(data)<-headers
  write.table(data, file=output, quote=F, sep='\t', row.names=F, col.names=T, append=app)
}

extract_engagement_data<-function(sdata){
  #ldata = subset(sdata, dayIndex<7 & requested>=join_timestamp & requested <= 1451127864 & device %in% c(1, 3, 10, 31, 32, 33, 36, 60))  
  g = sdata[1,"group"]
  print_with_time(sprintf("subsetting to session after join %s", g))  
  edata = subset(sdata, requested>=join_timestamp)
  edata_compact = edata[,c("userid","dayIndex","session_num","stumble_in_session","requested","device","method","stumbles","mstumbles")]
  
  print_with_time(sprintf("aggregating session sum for %s", g))  
  agg_sum = aggregate(cbind(stumbles,mstumbles) ~ userid + session_num, edata_compact, sum)
  print_with_time(sprintf("aggregating session min for %s", g))  
  agg_min = aggregate(cbind(device, requested, dayIndex) ~ userid + session_num, edata_compact, min)
  print_with_time(sprintf("aggregating session max for %s", g))  
  agg_max = aggregate(cbind(device, requested, dayIndex) ~ userid + session_num, edata_compact, max)
  
  names(agg_min) <- c("userid","session_num","minDevice","minRequested","minDayIndex")
  names(agg_max) <- c("userid","session_num","maxDevice","maxRequested","maxDayIndex")
  
  print_with_time(sprintf("Merging min and max for %s", g))  
  agg_minmax = merge(agg_min, agg_max)#159,718 sessions
  
  print_with_time(sprintf("subsetting to sessions with same dayIndex and same device for %s", g))  
  agg_minmax_legit = subset(agg_minmax, (minDevice==maxDevice | minDevice==0) & minDayIndex==maxDayIndex)#147,571 sessions; around 7K lost to dayIndex, 1.5K to device discr
  agg_minmax_legit$seconds = agg_minmax_legit$maxRequested-agg_minmax_legit$minRequested+10
  agg_minmax_legit$device_name = get_device(agg_minmax_legit$maxDevice)
  agg_minmax_legit = agg_minmax_legit[,c("userid","session_num","device_name","maxDayIndex","seconds")]
  names(agg_minmax_legit)<-c("userid","session_num","device_name","dayIndex","seconds")
  print_with_time(sprintf("Merging sum and minmax for %s", g))   
  with_stumbles = merge(agg_sum, agg_minmax_legit)
  print_with_time(sprintf("completed %s", g))  
  return(with_stumbles)
}

tutorial1<-function(){
  #http://www.ceb-institute.org/bbs/wp-content/uploads/2011/09/handout_ggplot2.pdf 
  require(ggplot2)
  head(diamonds[,c("clarity","cut")])
  #diamonds is 54K rows, clarity and cut are fields
  #stacked bar chart, columns are clarity, stacks are cut, within each bar
  #x=clarity, y=cut
  qplot(clarity, data=diamonds, fill=cut, geom="bar")
  #same chart, using ggplot
  ggplot(diamonds, aes(clarity, fill=cut)) + geom_bar()
  
  #mtcars is 32 rows, contains columns mpg and wt (weight)
  head(mtcars[,c("wt","mpg")])
  # scatterplot, shows weight negatively correlated with mpg
  #x=wt, y=mpg. scatterplot seems to be the default
  qplot(wt, mpg, data=mtcars)
  #transformations on fields of data frame is acceptable
  qplot(log(wt), mpg - 10, data=mtcars)
  #in the following, the value from a third column is used to color-code data points
  qplot(wt, mpg, data=mtcars, color=qsec)
  qplot(wt, mpg, data=mtcars, color=cyl)
  #make data points bigger, but adds a 3 legend, which is not intended
  qplot(wt, mpg, data=mtcars, color=qsec, size=3)
  #removes funny 3 legend from previous plot
  qplot(wt, mpg, data=mtcars, colour=qsec, size=I(3))#this makes it constant, instead of mapping. 
  # use alpha blending
  qplot(wt, mpg, data=mtcars,alpha=qsec) #values between 0 (transparent) and 1 (opaque)
  # continuous scale vs. discrete scale
  head(mtcars)
  qplot(wt, mpg, data=mtcars, colour=cyl)#continuous color
  levels(as.factor(mtcars$cyl)) # 4,6,8
  qplot(wt, mpg, data=mtcars, colour=factor(cyl))#discrete color
  
  # use different aesthetic mappings
  qplot(wt, mpg, data=mtcars, shape=factor(cyl))#uses different tickers instad of color: circle (default), triangle, square
  #here, size of data points indicates different qsec values
  qplot(wt, mpg, data=mtcars, size=qsec)
  #here you combine aes and show yet another column (info overload)
  qplot(wt, mpg, data=mtcars, size=qsec, shape=factor(cyl), geom="point")
  
  #instead of stack, if you just want to plot bar graph for one field
  qplot(factor(cyl), data=mtcars,geom="bar")
  
  #turns bars 90 degrees, or more accurately, flips x and y axes
  qplot(factor(cyl), data=mtcars, geom="bar") + coord_flip()
  
  #fill changes color of bars, while color changes line color of bars
  qplot(factor(cyl), data=mtcars, geom="bar", fill=factor(cyl))
  qplot(factor(cyl), data=mtcars, geom="bar", colour=factor(cyl))
  
  #stacked bar chart, like in the beginning
  qplot(factor(cyl), data=mtcars, geom="bar", fill=factor(gear))
  
  #4 different variations of stacked bar chart
  qplot(clarity, data=diamonds, geom="bar", fill=cut, position="stack")#stacked bar chart
  qplot(clarity, data=diamonds, geom="bar", fill=cut, position="fill")#completes to 1 = 100%
  qplot(clarity, data=diamonds, geom="bar", fill=cut, position="dodge")#grouped bar chart
  qplot(clarity, data=diamonds, geom="bar", fill=cut, position="identity")#plots on top of each other, pretty useless
  
  #line charts, only identity and stack defined
  qplot(clarity, data=diamonds, geom="freqpoly", group=cut, colour=cut, position="identity")#as -is line chart
  qplot(clarity, data=diamonds, geom="freqpoly", group=cut, colour=cut, position="stack")#height is cumulative
}

tutorial2<-function(){
  #http://www.cookbook-r.com/Graphs/Bar_and_line_graphs_%28ggplot2%29/
  # A bar graph
  ggplot(data=dat1, aes(x=time, y=total_bill, fill=sex)) + 
    geom_bar(colour="black", stat="identity",
             position=position_dodge(),
             size=.3) +                        # Thinner lines
    scale_fill_hue(name="Sex of payer") +      # Set legend title
    xlab("Time of day") + ylab("Total bill") + # Set axis labels
    ggtitle("Average bill for 2 people") +     # Set title
    theme_bw()
  
  
  # A line graph
  ggplot(data=dat1, aes(x=time, y=total_bill, group=sex, shape=sex, colour=sex)) + 
    geom_line(aes(linetype=sex), size=1) +     # Set linetype by sex
    geom_point(size=3, fill="white") +         # Use larger points, fill with white
    expand_limits(y=0) +                       # Set y range to include 0
    scale_colour_hue(name="Sex of payer",      # Set legend title
                     l=30)  +                  # Use darker colors (lightness=30)
    scale_shape_manual(name="Sex of payer",
                       values=c(22,21)) +      # Use points with a fill color
    scale_linetype_discrete(name="Sex of payer") +
    xlab("Time of day") + ylab("Total bill") + # Set axis labels
    ggtitle("Average bill for 2 people") +     # Set title
    theme_bw() +
    theme(legend.position=c(.7, .4))           # Position legend inside
  # This must go after theme_bw
}

bar_vs_line_chart<-function(){
  random_normal = rnorm(1000)
  h = hist(random_normal, breaks=30, plot=FALSE)
  xy = data.frame(x_values=h$mids,Random=h$counts/sum(h$counts),Theoretical=calc_pnorm_bars(h$mids))
  
  library(ggplot2)
  library(reshape)
  library(plyr)
  library(scales)  
  
  #to see theme plot elemts, say ?theme
  x=h$mids;df=data.frame(x=x,Random=h$counts/sum(h$counts),Theoretical=calc_pnorm_bars(h$mids))
  ggplot(melt(df),aes(x,value,color=variable,fill=variable))+
    geom_bar(subset=.(variable=="Random"),stat="identity")+
    geom_line(subset=.(variable=="Theoretical"), size=1, linetype="dashed")+
    scale_fill_manual(values=c("#E69F00","blue"))+
    scale_color_manual(values=c("#E69F00","blue"))+
    xlab("bin") + ylab("count") +
    ggtitle("Random normal variable")+
    theme_bw() +
    annotate('text', x = 3, y = 0.07, 
             label = "mu==0~sigma^2==1",parse = TRUE,size=10)+ #\u03c3 is sigma
    theme(legend.position=c(.9, .9),
          title=element_text(family="Helvetica", face="bold", size=20, vjust=0.8),
          legend.title=element_blank(),
          legend.text=element_text(family="Helvetica", face="plain", size=12, hjust=0.8))
}

convert_to_vector<-function(df, ticker_shape){
  if (length(ticker_shape)==1){
    ticker_shapes = rep.int(ticker_shape,dim(df)[2]-1);
    names(ticker_shapes) = setdiff(names(df),"x")
  } else {
    ticker_shapes = ticker_shape
  }  
  return(ticker_shapes)
}



multiple_lines_charts<-function(df, 
                                #first column needs to be named x all other fields assumed lines to be plotted against x
                                ticker_size = 4,        #integer
                                ticker_shape = 21,      #round, could be specified as a named vector, fields same as df
                                ticker_fill = "black",  #could be specified as a named vector, fields same as df
                                line_size = 1,          #could be specified as a named vector, fields same as df
                                line_type = "solid",    #could be specified as a named vector, fields same as df
                                line_colors,
                                xlabel = "x axis please specify",
                                ylabel = "y axis please specify",
                                title = "graph title please specify",
                                legend_position = c(0.9, 0.9),
                                title_font_type = 'Helvetica',
                                title_font_weight = "bold",
                                title_font_size = 20,
                                legend_font_type = 'Helvetica',
                                legend_font_weight = "plain",
                                legend_font_size = 12,
                                x_axis_interval#could be interval or vector of values
                                ){

  ticker_shapes = convert_to_vector(df, ticker_shape)
  ticker_fills = convert_to_vector(df, ticker_fill)
  line_sizes = convert_to_vector(df, line_size)
  line_types = convert_to_vector(df, line_type)  

  
  library(ggplot2)
  library(reshape)
  library(plyr)
  library(scales)  
  
  #to see theme plot elemts, say ?theme
  mdf <- melt(df, 1, variable_name="Type")#keep column 1 as-is and melt others. 
  #which means you now have Type=column name, and third column (by default named value) is the value of that column, for that first column value
  
  if (length(x_axis_interval)==1){
    x_axis_values = round(seq(min(mdf$x), max(mdf$x), by = x_axis_interval),1)
  } else {
    x_axis_values = x_axis_interval
  }
  
#   if (length(y_axis_interval)==1){
#     y_axis_values = round(seq(min(mdf$value), max(mdf$value), by = y_axis_interval),1)
#   } else {
#     y_axis_values = y_axis_interval
#   }
  
  #shape indices: http://www.cookbook-r.com/Graphs/Shapes_and_line_types/
  
  myplot = ggplot(data=mdf, aes(x=x, y=value, group = Type, colour = Type, size=Type, linetype=Type, fill=Type)) +
    geom_line() +              
    geom_point(aes(shape=Type), size=4)+ #ticker size
    scale_color_manual(values=line_colors)+ #line color
    scale_fill_manual(values=ticker_fills)+ #ticker fill color
    scale_size_manual(values=line_sizes) +#line size
    scale_shape_manual(values=ticker_shapes) +                  # ticker shapes
    scale_linetype_manual(values=line_types)+ # line type
    scale_x_continuous(breaks = x_axis_values) +
#     scale_y_continuous(breaks = y_axis_values) + 
    xlab(xlabel) + ylab(ylabel) + #axis label
    ggtitle(title)+
    theme_bw() +
    theme(legend.position=legend_position,
          axis.line = element_line(colour = "black", size = 1),
          axis.text.x = element_text(size = 14),
          axis.title.x=element_text(vjust=0.1),
          axis.text.y = element_text(size = 14),
          axis.title.y=element_text(vjust=0.7),
          title=element_text(family=title_font_type, face=title_font_weight, size=title_font_size, vjust=0.8),
          legend.title=element_blank(),
          legend.text=element_text(family=legend_font_type, face=legend_font_weight, size=legend_font_size, hjust=0.8))
  return(myplot)
}



multiple_lines_chart_demo<-function(){
  
  random_normal = rnorm(1000)
  h = hist(random_normal, breaks=30, plot=FALSE)
  xy = data.frame(x_values=h$mids,Random=h$counts/sum(h$counts),Theoretical=calc_pnorm_bars(h$mids))
  
  library(ggplot2)
  library(reshape)
  library(plyr)
  library(scales)  
  
  #to see theme plot elemts, say ?theme
  x=h$mids;df=data.frame(x=x,Random=h$counts/sum(h$counts),Theoretical=calc_pnorm_bars(h$mids))
  
  mdf <- melt(df, 1, variable_name="Type")#keep column 1 as-is and melt others
  
  #shape indices: http://www.cookbook-r.com/Graphs/Shapes_and_line_types/
  
  ggplot(data=mdf, aes(x=x, y=value, group = Type, colour = Type, size=Type, linetype=Type, fill=Type)) +
    geom_line() +              
    geom_point(aes(shape=Type), size=4)+
    scale_color_manual(values=c("Theoretical" = "blue","Random" = "red"))+
    scale_fill_manual(values=c("Theoretical" = "white","Random" = "black"))+
    scale_size_manual(values=c("Theoretical" = 1,"Random" = 1)) +#line size
    scale_shape_manual(values=c("Theoretical" = 21,"Random" = 21)) +                  # Change shapes
    scale_linetype_manual(values=c("solid", "dotted"))+ # Change linetypes
    xlab("bin") + ylab("count") +
    ggtitle("Random normal variable")+
    theme_bw() +
    annotate('text', x = 3, y = 0.07, 
             label = "mu==0~sigma^2==1",parse = TRUE,size=10)+ #\u03c3 is sigma
    theme(legend.position=c(.9, .9),
          title=element_text(family="Helvetica", face="bold", size=20, vjust=0.8),
          legend.title=element_blank(),
          legend.text=element_text(family="Helvetica", face="plain", size=12, hjust=0.8))
}

misc<-function(){
  require(ggplot2)
  random_normal = rnorm(1000)
  h = hist(random_normal, breaks=30, plot=FALSE)
  xy = data.frame(x_values=h$mids,Random=h$counts/sum(h$counts),Theoretical=calc_pnorm_bars(h$mids))
  
  xy2 = data.frame(x_values=h$mids,y_values=calc_pnorm_bars(h$mids))
  ggplot(data=xy, aes(x=x_values, y=y_values, color="Random"))+
    geom_bar(colour="black", fill="#DD8888",stat="identity") + #identity means plot y value as-is. default is bin
    geom_line(aes(x=xy2$x_values, y=xy2$y_values, color="Theoretical"), size=1, colour="red", linetype="dashed")+
    #geom_point(colour="red", size=2, shape=21, fill="white")+this adds points to all series
    #expand_limits(y=0) + reset the max and min for scale
    #guides(fill=FALSE) + #turns off legend, as far as I can tell
    xlab("bin") + ylab("count") +
    ggtitle("Random normal variable")
  
  
  ggplot(data=xy,aes(x=x_values))+
    geom_bar(aes(y=y_values, fill=label),colour="black",stat="identity")+
    geom_line(aes(x=xy2$x_values, y=xy2$y_values, color="Theoretical"), size=1, linetype="dashed")+
    scale_fill_manual(values=c("#E69F00"))+ #applies to bar chart
    scale_color_manual(values=c("blue"))+ #applies to line
    scale_fill_discrete(breaks=c("Random","Theoretical"))
  
  
  
  #m<-ggplot(data=hst, aes(x=days))+geom_line(aes(y=retained, color="Retained"),size=1) 
  #m+geom_line(aes(y=left, color="Left"), size=1)+
  #labs(y = 'probability density', title=paste(attribute, retention_field, sep=" vs. "), color = 'Customers')
  #+scale_color_manual(values=c("Retained"="blue4", "Left"="red4"))
  #+theme(panel.background = element_rect(fill='white'), 
  #legend.text = element_text(size=12), 
  #legend.title = element_text(size=12))
}

calc_pnorm_bars<-function(x){
  #given mid points of histogram, calculates density of that bar by taking difference of distribution function
  require(binhf)#contains shift function
  y = pnorm(x)-pnorm(shift(x,1))
  y[1] = pnorm(x[1])
  return(y)
}

step<-function(k){
  if (k==0){
    return(1)
  } else {
    if (runif(1)<0.5){
      return(k-1)
    } else {
      return(k+1)
    }
  }
}

expected<-function(start,end,sample_size=1000){
  results = rep(NA, sample_size)
  for (iter in seq(1,sample_size,1)){
    #print(paste("iter: ",iter))
    steps = 0;
    current = start
    while (current<end){
      steps = steps+1
      current = step(current)
      #print(paste("steps:",steps,", current:", current))
    }
    results[iter] = steps
  }
  #print(results)
  return(mean(results))
}


process_retention<-function(threshold=0){
  filename="/Users/oeskiyenenturk/su/datascience/AB_testing/retention_2015-10-07.tsv"
  col_names = unlist(strsplit("group_name,user_id,days_user_came_back,ignore_if_lt13,mstumbles,mgreats,stumbles,greats,incat_stumbles,incat_greats",","))
  col_classes = unlist(strsplit("character,numeric,numeric,numeric,numeric,numeric,numeric,numeric,numeric,numeric",","))
  
  retention<-read.delim(filename, header = F,
                        col.names = col_names,
                        colClasses = col_classes
  )
  print(sprintf("Original dimension: %.0f",dim(retention)[1]))
  retention <- retention[complete.cases(retention),] #30,705
  print(sprintf("Complete cases: %.0f",dim(retention)[1]))
  fair = subset(retention, ignore_if_lt13>=13)       #18,042 
  print(sprintf("Users reaching day 13: %.0f",dim(fair)[1]))
  fair$isEngaged = ifelse(fair$days_user_came_back>threshold,1,0)
  
  default = subset(fair, group_name=="Default")
  random = subset(fair, group_name=="with_random_like_topic_experts")
  sorted = subset(fair, group_name=="with_sorted_like_topic_experts")
  qs = c(0, quantile(default$mstumbles,seq(0.1,1,0.1)))
  for (n in seq(1,10,1)){
    #default_n = subset(default, mstumbles>=n)
    #random_n = subset(random, mstumbles>=n)
    #sorted_n = subset(sorted, mstumbles>=n)
    
    default_n = subset(default, mstumbles>qs[n] & mstumbles<=qs[n+1])
    random_n = subset(random, mstumbles>qs[n] & mstumbles<=qs[n+1])
    sorted_n = subset(sorted, mstumbles>qs[n] & mstumbles<=qs[n+1])
    
    if (mean(default_n$isEngaged)<mean(random_n$isEngaged)){
      tt = t.test(default_n$isEngaged, random_n$isEngaged, alternative="less") 
      sign1 = "<"
    } else {
      tt = t.test(default_n$isEngaged, random_n$isEngaged, alternative="greater") 
      sign1 = ">"
    }
    if (mean(default_n$isEngaged)<mean(sorted_n$isEngaged)){
      tt2 = t.test(default_n$isEngaged, sorted_n$isEngaged, alternative="less")  
      sign2 = "<"
    } else {
      tt2 = t.test(default_n$isEngaged, sorted_n$isEngaged, alternative="greater") 
      sign2 = ">"
    }    
    str = sprintf("%.0f,%.0f,%.0f,%.0f,%.0f,%.1f%% %s %.1f%% (p-value = %f),%.1f%% %s %.1f%% (p-value = %f)",qs[n],qs[n+1],dim(default_n)[1], dim(random_n)[1], dim(sorted_n)[1],
                  tt$estimate[1]*100, sign1, tt$estimate[2]*100, tt$p.value,tt2$estimate[1]*100, sign2, tt2$estimate[2]*100, tt2$p.value)
    print(str)
  } 
  #write.table(subset(retention, group_name=="with_sorted_like_topic_experts"), file="/Users/oeskiyenenturk/su/datascience/AB_testing/sorted_mstumbles.tsv", quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
}

print_array<-function(name, value, sp=","){
  print(paste(name, " = [", paste(value,collapse=sp),"]", sep=""))
}

cube_step<-function(position){
  x=runif(1)
  if(x<(1/3)){
    position[1] = abs(1-position[1])
  } else if (x<(2/3)){
    position[2] = abs(1-position[2])
  } else{
    position[3] = abs(1-position[3])
  }
  return(position)
}

cube_steps<-function(current){
  iter = 0
  while(sum(current)<3){
    iter = iter+1
    current = cube_step(current)
    #print(current)
  }
  return(iter)
}

expected_steps<-function(start, num_iterations=100){
  results = rep(0,num_iterations)
  for (k in seq(1,num_iterations,1)){
    results[k] = cube_steps(start)
  }
  return(mean(results))
}

term<-function(n,k,i){
  return(choose(k-1+i, k-1)*factorial(k-1)*choose(n-k,i)*factorial(i)*factorial(n-k-i))
}

denominator<-function(n,k){
  sm = 0
  for (i in seq(0, (n-k),1)){
    t = term(n,k,i)
    sm = sm+t
  }
  return(sm)
}

numerator<-function(n,k){
  return(term(n,k,(n-k)))
}

ratios<-function(n){
  results = rep(0,n)
  for (k in seq(1,n,1)){
    results[k] = numerator(n,k)/denominator(n,k)
  }
  return(results)
}

read_mau<-function(filename='/Users/oeskiyenenturk/su/datascience/DATA-341/mau_201505.tsv'){
  #for edges formed within 201505, count distinct user_followed and count distinct user_following: 252462, 64277
  col_names = unlist(strsplit("userid,created,stumbles,month,user_follows,user_followed_by,first_time_followed_by,first_time_following",","))
  col_classes = unlist(strsplit("numeric,numeric,numeric,character,numeric,numeric",","))  
  mau<-read.delim(filename, header = F,
                  col.names = col_names,
                  colClasses = col_classes,
                  skipNul = T
  ) 
  total = dim(mau)[1]
  cat(sprintf("Read %0.f rows from %s\n",total,filename))
  mau = mau[complete.cases(mau),]
  cat(sprintf("complete.cases: %0.f rows\n",dim(mau)[1]))
  following_only = subset(mau, user_follows>0 & user_followed_by==0)
  followed_only = subset(mau, user_follows==0 & user_followed_by>0)
  both = subset(mau, user_follows>0 & user_followed_by>0)
  num_both = dim(both)[1]
  none = subset(mau, user_follows==0 & user_followed_by==0)
  exclude_none = subset(mau, user_follows>0 | user_followed_by>0)
  num_none = dim(none)[1]
  users_following = rbind(following_only,both)
  num_users_following = dim(users_following)[1]
  users_followed = rbind(followed_only, both)
  num_users_followed = dim(users_followed)[1]
  
  #cat(sprintf("%s,%0.f,%0.f,%0.f,%0.f\n",mau$month[1],dim(following_only)[1],dim(followed_only)[1],num_both,num_none))
  
  cat(sprintf("Only following: %0.f ( %0.f%% )\n", dim(following_only)[1], dim(following_only)[1]/total*100))
  cat(sprintf("Only followed: %0.f ( %0.f%% )\n", dim(followed_only)[1], dim(followed_only)[1]/total*100)) 
  cat(sprintf("Both: %0.f ( %0.f%% )\n", num_both, num_both/total*100)) 
  cat(sprintf("None: %0.f ( %0.f%% )\n", num_none, num_none/total*100)) 
  
  cat(sprintf("Number of users following someone: %0.f (%0.f%% of mau)\n", num_users_following, num_users_following/total*100))
  cat(sprintf("Number of users followed by someone: %0.f (%0.f%% of mau)\n", num_users_followed, num_users_followed/total*100))
  cat(sprintf("Average out-degree : %2.2f, average out-degree excluding None : %2.2f\n", 
              mean(mau$user_follows), mean(exclude_none$user_follows)))#3.65, 17.28
  cat(sprintf("Median out-degree : %2.2f, median out-degree excluding None : %2.2f\n", 
              quantile(mau$user_follows, 0,5), quantile(exclude_none$user_follows, 0.5)))#3.65, 17.28
  cat(sprintf("Average in-degree : %2.2f, average in-degree excluding None : %2.2f\n", 
              mean(mau$user_followed_by), mean(exclude_none$user_followed_by))) #3.63, 11.48
  cat(sprintf("Median in-degree : %2.2f, median in-degree excluding None : %2.2f\n", 
              quantile(mau$user_followed_by, 0.5), quantile(exclude_none$user_followed_by, 0.5))) #3.63, 11.48 
}

read_maus<-function(root = "/Users/oeskiyenenturk/su/datascience/DATA-341"){
  files = list.files(path=root,pattern="^mau.*tsv$")  
  for (file in files){
    #replace_nulls(root, file)
    read_mau(file.path(root, file))
  }
}

replace_nulls<-function(root="/Users/oeskiyenenturk/su/datascience/DATA-341", filename="mau_201505.tsv"){
  command = sprintf("cat %s | sed 's.NULL..g' > %s;mv %s %s", file.path(root, filename), file.path(root, "dummy"),
                    file.path(root, "dummy"), file.path(root, filename))
  print(command)
  system(command)
}

recip<-function(){
  col_names = unlist(strsplit("rank,is_followed_back,count",","))
  col_classes = unlist(strsplit("numeric,numeric,numeric",","))  
  recip<-read.delim("/Users/oeskiyenenturk/su/datascience/DATA-341/recip_by_rank_201509.tsv", header = F,
                    col.names = col_names,
                    colClasses = col_classes,
                    skipNul = T
  ) 
  agg = aggregate(count ~ rank, recip, FUN=sum)
  perc = quantile(rep(agg$rank, agg$count), seq(0.05, 0.95, 0.05))
}

maubots<-function(){
  col_names = unlist(strsplit("userid,level,dupe,prob",","))
  col_classes = unlist(strsplit("numeric,numeric,numeric,numeric",","))  
  recip<-read.delim("/Users/oeskiyenenturk/su/sucron/cron/analytics/UsersEnrichedNew/mau_combined.tsv", header = F,
                    col.names = col_names,
                    colClasses = col_classes,
                    skipNul = T
  ) 
  recip$dummy = 1
  recip$predicted = ifelse(recip$prob>0.8,1,0)
  recip$actual = ifelse(recip$level==0,1,0)
  aggregate(dummy ~ predicted + actual, recip, FUN=sum)  
  
}

