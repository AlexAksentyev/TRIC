source("Utilities.R")
library(dplyr)

## core functions ####
loadData <- function(file = "CosyAll_2016-06-24_11-57-53.dat", where = "~/Analysis/Jun16/DevReader/Data/", ...){
  require(data.table); require(tidyr)
  
  xtra = list(...)
  bound = NULL
  if("from" %in% names(xtra)) bound <- c(bound, start = timeToUTS(xtra$from))
  if("to" %in% names(xtra)) bound <- c(bound, finish = timeToUTS(xtra$to))
  
  ## loading all data
  Data = data.frame(fread(paste(where,file,sep="")))
  if("Fields" %in% names(xtra)) .flds = xtra$Fields
  else
    .flds = c("WD", "Month", "Day", "Time", "Year", as.vector(t(outer(c("BCT","FCT"), 1:2, paste, sep=""))), paste("SB", 1:4, sep=""))
  names(Data)<- .flds

  ## cleaning time data
  # Data$Date = paste(Data$Year, match(Data$Month, month.abb), Data$Day, sep="-")
  # Datetime = paste(Data$Date, Data$Time, sep=" ")
  # Data$UTS = timeToUTS(Datetime); rm(Datetime)
  Data %>% mutate(Month = match(Month, month.abb)) %>% 
    unite(Date, Year, Month, Day, sep="-") %>% unite(Clock, Date, Time, sep=" ") %>%
    mutate(Clock = as.POSIXct(Clock, origin="1970-1-1"), UTS = as.numeric(Clock)) -> Data
  
  ## picking relevant time frame
  if("start" %in% names(bound)) Data <- Data[Data$UTS >= bound["start"], ]
  if("finish" %in% names(bound)) Data <- Data[Data$UTS <= bound["finish"], ]
  rm(bound)
  
  ## rearranging columns
  Data <- Data%>%dplyr::select(-WD)
  # Data <- subset(Data, select = -c(1:3,5))
  # Data <- subset(Data, select = c(11,10,1:9))
  
  return(Data)
  
}

cleanData <- function(.Data){
  require(plyr); require(mosaic)
  
  ## cleaning spin data ##
  sb.i = which(substr(names(.Data),1,2) == "SB") %>% sort(decreasing=TRUE)
  
  sapply(sb.i, function(i, X) cut(X[,i], breaks = 2, labels = c("0", "1")), .Data) -> .Data[,sb.i]
  
  .Data %>% dplyr::select(sb.i) %>% (function(x) do.call(paste0, x)) %>% strtoi(base=2) %>%
    factor(levels=c(1,2,15), labels = c("U","D","N")) -> .Data$B.Spin
  # .Data <- subset(.Data, select = -sb.i); rm(sb.i)
  
  ## categorizing ##
  ## current
  .Data$FSgl = with(
    .Data, 
    cut(BCT2, breaks = c(min(BCT2)-100, findCut(BCT2), max(BCT2)+100), labels = c("F", "T"))
  )
  ## unit
  .Data %>% .unitize(.cond = c("F" = 1, "T" = 270)) -> .Data
   
  ## cycle units at this point include not only experiment 6, but also 5 instances of experiment 2; 
  ## the tail of a cycle ends abruptly enough, so we know where an exp6 cycle ends
  ## i remove those points more than 792 seconds away from the cycle end.
  ## those cycles only partially included (duration less than 792) are dropped from analysis
  cyc.dur = 792; # the duration of an exp6 cycle
  f = function(.sub){
    if(.sub$FSgl[1] == "F") return(.sub)
    
    rng = range(.sub$UTS); .sub.dur = diff(rng)
    if(.sub.dur>cyc.dur){
      which(.sub$UTS >= rng[2]-cyc.dur)[1]->i
      
      .sub <- .sub[-c(1:(i-1)),]
      
      if(rng[2] - .sub$UTS[1] < cyc.dur)
        data.frame(UTS = .sub$UTS[1]-1, Unit = .sub$Unit[1], FSgl = "T") %>% rbind.fill(.sub) -> .sub
        
      return(.sub)
    }
  }
  ddply(.Data, "Unit", f) -> .Data

  ## adding data ##
  ## target
  list.files(path = "./Data", pattern = "ABS") -> ABS.list
  llply(ABS.list, getABS, .frame = FALSE) %>% unlist() -> ABS
  which(.Data$UTS %in% as.numeric(names(ABS))) -> i
  if(length(i)>0){
    .Data$ABS <- NA
    .Data$ABS[i] <- ABS[as.character(.Data$UTS[i])]
    .Data$FABS <- cut(.Data$ABS, breaks = 2, labels = c("F","T"))
    rm(ABS)
  }else{
    f <- function(.sub){
      
      tg.sw = c("On" = 68, "->Off" = 430, "Off" = 430 + 90) + .sub$UTS[1]

      .sub %>% mutate(FABS = derivedFactor(
          "F" = (UTS > tg.sw[3]),
          "T" = (UTS > tg.sw[1] & UTS <= tg.sw[2]),
          .default = NA,
          .method = "last"
        )
      ) -> .sub
      
      return(.sub)
    }
    .Data <- ddply(.Data, "Unit", f)
  }
  ## cooler
  f <- function(.sub){

    ec.sw = c("On" = 1,"120" = 48, "45" = 78, "Off" = 770) + .sub$UTS[1]
    
    .sub %>% mutate(eCool = derivedFactor(
      "On" = (UTS > ec.sw[1]),
      "Inj"  = (UTS > ec.sw[2]),
      "Acc" = (UTS > ec.sw[3]),
      "Off" = (UTS > ec.sw[4]),
      .default = NA,
      .method = "last"
    )) -> .sub
    
    return(.sub)
    
  }
  .Data <- ddply(.Data, "Unit", f)
  
  ## transformations ##
  if(FALSE){
    filter(.Data, FSgl == "F") %>% dplyr::select(BCT1, BCT2) %>% as.matrix() %>%
      aaply(.margins = 2, .fun = median, na.rm = TRUE) -> ADC.offs; cat(ADC.offs[2]); cat("\n")
    .Data %>% mutate(BCT1 = BCT1 - ADC.offs[1], BCT2 = BCT2 - ADC.offs[2]) -> .Data
  }
  
  
  return(.Data)
}

getABS <- function(.file, .where = "./Data/", .frame = TRUE){

  con = file(paste(.where, .file, sep=""))
  
  ## get the initial time stamp
  readLines(con, n = 1) %>% (function(x) strsplit(x, split = c(" ", "\t"))[[1]]) %>% 
    (function(x) gsub("\t", "", x)) -> datetime
  datetime <- strptime(paste0(datetime[1], datetime[2], collapse = " "), format = "%Y%m%d %H%M%S")
  stt.uts = timeToUTS(datetime)
  
  ABS = read.table(con, sep = "\t", quote = "", skip = 14, col.names = c("UTS", "Flux")); 
  ABS <- mutate(ABS, UTS = UTS + stt.uts)
  
  if(!.frame){
    Flux = ABS$Flux; names(Flux) <- ABS$UTS
    return(Flux)
  }
  
  return(ABS)
}

## support functions ####
findModes<- function(x) {
  require(quantmod)
  
  unlist(sapply(2:(length(x)-1), FUN = function(i,x) if(x[i]>x[i-1] & x[i] > x[i+1]) i, x)) -> modes
  if ( length(modes) == 0 )
    modes = 'This is a monotonic distribution'
  
  return(modes)
}

.smartCut <- function(Data, which = "BCT2"){
  cut(
    Data[, which], 
    breaks = c(min(Data[, which])-100, findCut(Data[, which]), max(Data[, which])+100), 
    labels = c("F", "T")
  )
}
.unitize <- function(Data, .cond=c("F" = 1, "T" = 270)){
  rle(as.numeric(Data$FSgl)) -> x
  Data %>% mutate(Unit = rep(1:length(x$lengths), times = x$lengths)) %>%
    ddply("Unit", function(u) if(nrow(u)>=.cond[u$FSgl[1]]) u)
}

findCut <- function(.data){
  require(quantmod)
  density(.data)->d
  ind.v = findValleys(d$y); brk = d$x[ind.v[1]]
  
  return(brk)
}