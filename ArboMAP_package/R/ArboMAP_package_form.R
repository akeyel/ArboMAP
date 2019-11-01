#' ArboMAP: Arbovirus Modeling and Prediction to Forecast Mosquito-Borne Disease Outbreaks
#'
#' Predict mosquito-borne disease outbreaks using distributed lags approach.
#' 
#' @details Code developed by JKD and MCW, Geography and Environmental Sustainability, University of Oklahoma.
#' ACK reformatted the code to be in package format.
#'
#'@author Justin K. Davis, Michael C. Wimberly (justinkdavis@ou.edu, mcwimberly@ou.edu), Alexander C. Keyel (akeyel@albany.edu)
#' @name ArboMAP-package
NULL

#' Simplify Names (ArboMAP)
#'
#' The simplify names function converts district names to a more consistent format.
#' It does this by converting to lower case, removing spaces, removing "county" or "parish"
#' from the district name.
#'
#' For example "New Haven County" and "new haven" would both be converted to "newhaven"
#'
#' It does not correct for abbreviations. For example, it would not correctly
#' simplify Saint Lawrence, St. Lawrence and St Lawrence into a consistent name.
#'
#' It does not differentiate between repeat county names from different states.
#' E.g. it would have no way to distinguish Jefferson County in Colorado from
#' Jefferson County in New York.
#'
#' @details Written by J. Davis and M. Wimberly. Modified by A. Keyel. Original available at www.github.com/ecograph/ArboMAP. Used under GPL-3 license
#'
#' @param priornames A vector of names to be simplified
#'
#' @export simplifynames
simplifynames <- function(priornames=NULL) {

  # convert to lower case
  priornames <- tolower(priornames)
  
  # remove spaces
  priornames <- sub(pattern=" ", replacement="", x=priornames, fixed=TRUE)
  
  # remove district and parish
  priornames <- sub(pattern="county", replacement="", x=priornames, fixed=TRUE)
  priornames <- sub(pattern="parish", replacement="", x=priornames, fixed=TRUE)

  # return names
  return(priornames)
  
}

#packages <- c("reshape2", "ggplot2", "gridExtra",
#              "lme4","pracma","dplyr","maptools",
#              "raster","spdep","mgcv","sp","rgdal",
#              "GISTools","data.table","splines","maps",
#              "broom","mapproj", "Hmisc")

#' ArboMAP (main external function call)
#' 
#' ArboMAP is a set of software to be used in the RStudio envionment to model and predict vector-borne diseases, especially arboviruses transmitted by mosquitoes. ArboMAP was developed by the EcoGRAPH research group at the University of Oklahoma with support from NASA. We are happy to answer your questions, and we would appreciate feedback on how these tools could be improved. Please contact justinkdavis@ou.edu for technical issues about the original code, or mcwimberly@ou.edu for questions about the arbovirus modeling project. For issues related to the package/function modifications, please contact akeyel@albany.edu.
#'
#' @param humandatafile Data on human cases, by date of onset. For models that do
#' not consider date of onset, these data will be recompiled. This can be given as the name of a text file, or can be entered as an R data object.
#' @param mosqfile Data on district of collection, collection date,
#' west nile virus pool result, size of the mosquito pool, and species. This can be the name of a text file or as an R data object.
#' @param districtshapefile A shapefile containing the districts in polygon format.
#' @param stratafile A text file with data on how to stratifiy the district polygons
#' @param weatherpathstr A directory containing the weather data to be aggregated.
#' The weather data can contain multiple files, and these can overlap in date range,
#' as long as they all have the same columns.
#' The weather data files in this folder should have a district column, a doy (day of year) column
#' a year column, and then an entry for each independent variable of interest
#' @param weathersummaryfile If this file is in weatherpathstr, it will be ignored.
#' Note that no consolidated summary file is currently being created.
#' @param maxobservedhumandate The last date for which human data should be included. This should NOT include the year for which forecasts are made.
#' @param weekinquestion The week to be used in the analysis
#' @param var1name The first variable to be used in the model. These can be identified using the model selection code.
#' @param var2name The second variable to be used in the model.
#' @param compyear1 The first year to use for comparison to the current year
#' @param compyear2 The second comparison year
#' @param results.path The path for the model results
#' @param original.metrics 1 runs the code to correspond to the original version of ArboMAP. 0 makes departures from the code for the calculation of number of human cases predicted in the year
#'
#'@export ArboMAP
ArboMAP = function(humandatafile, mosqfile, districtshapefile, stratafile, weatherpathstr, weathersummaryfile, maxobservedhumandate, weekinquestion, var1name, var2name, compyear1, compyear2, results.path, original.metrics = 0){
# inputs coming in a future version of ArboMAP
#  #' @param modeltype The type of model to be used. Choices are 'fixeddfcr' for a fixed df cubic regression spline or 'thinplate' for a thin-plate spline. 'fixeddfcr' is the default.
#  #' @param use.cluster Indicator variable. If 1, the code will attempt to use parallel processing for the second human data step. Otherwise it will not. 0 is the default.
  
 
  # where do we want the outputs?
  graphicoutputdir <- sprintf("%s\\graphical outputs\\", results.path)
  fullcasematoutputdir <- sprintf("%s\\case matrix with estimates\\",   results.path)
  mosqmatoutputdir <- sprintf("%s\\mosquito matrix with estimates\\", results.path)
  
  # Create any missing directories
  if (!file.exists(graphicoutputdir)){dir.create(graphicoutputdir, recursive = TRUE)}
  if (!file.exists(fullcasematoutputdir)){dir.create(fullcasematoutputdir, recursive = TRUE)}
  if (!file.exists(mosqmatoutputdir)){dir.create(mosqmatoutputdir, recursive = TRUE)}
  
  # probably don't want to modify what follows, but you have some options if you're comfortable

  # makes sure we round to the previous Sunday, so that this week is included
  weekinquestionSun <- weekinquestion - (as.numeric(strftime(weekinquestion, '%u')) %% 7)
  weekinquestionSat <- weekinquestionSun + 6
  weekinquestionSunstr <- strftime(weekinquestionSun, '%A')
  weekinquestionSatstr <- strftime(weekinquestionSat, '%A')

  # figure out which year this is and begin all weeks on Sunday
  maxmosqyear <- as.numeric(format(weekinquestion, "%Y"))
  maxdesiredhumandate <- as.Date(paste(maxmosqyear,
                                     "-12-31",sep=""))
  maxdesiredhumandate <- maxdesiredhumandate - (as.numeric(strftime(maxdesiredhumandate, '%u')) %% 7)
  maxdesiredhumandatestr <- strftime(maxdesiredhumandate, '%A')

  # make sure the max desired human date is no earlier than the week in question
  # otherwise, there will be no predictions for the week in question
  maxdesiredhumandate <- max(weekinquestionSun+7, maxdesiredhumandate)

  # set up lag and regression information
  laglen   <- 181
  dlagdeg  <- 8

  results = call.ArboMAP(humandatafile, mosqfile, districtshapefile, stratafile, weatherpathstr, weathersummaryfile, maxobservedhumandate, weekinquestion, var1name, var2name, compyear1, compyear2, results.path, graphicoutputdir, fullcasematoutputdir, mosqmatoutputdir, weekinquestionSun, weekinquestionSat, weekinquestionSunstr, weekinquestionSatstr, maxmosqyear, maxdesiredhumandate, maxdesiredhumandatestr, laglen, dlagdeg, original.metrics)

  # Create a list of inputs that can be used later as function inputs
  # This does not include inputs that were fed into the function, as those will already be in memory
  inputs = list(graphicoutputdir, fullcasematoutputdir, mosqmatoutputdir, weekinquestionSun, weekinquestionSat, weekinquestionSunstr, weekinquestionSatstr, maxmosqyear, maxdesiredhumandate, maxdesiredhumandatestr, laglen, dlagdeg)
  
  return(list(results, inputs))
}

# Create a sub-function to contain all the individual function calls
call.ArboMAP = function(humandatafile, mosqfile, districtshapefile, stratafile, weatherpathstr, weathersummaryfile, maxobservedhumandate, weekinquestion, var1name, var2name, compyear1, compyear2, results.path, graphicoutputdir, fullcasematoutputdir, mosqmatoutputdir, weekinquestionSun, weekinquestionSat, weekinquestionSunstr, weekinquestionSatstr, maxmosqyear, maxdesiredhumandate, maxdesiredhumandatestr, laglen, dlagdeg, original.metrics){
  
  #**# Place plots that would be knit into a pdf into a standalone pdf?
  # Call each chunks as a function
  
  # Read in weather data
  if (typeof(weatherpathstr) == "character"){
    weather = read.weather.data(weatherpathstr, weathersummaryfile)
  # If it is already read in, just use the data object
  }else{  weather = weatherpathstr  }
  
  
  # Plot weather data and create daily extremes object
  dailyextr = plot.weather.data(weather, graphicoutputdir, weekinquestionSun)
  
  # Calculate vector infection data
  wnv.out = vector.infection.data(mosqfile, maxmosqyear)
    wnv = wnv.out[[1]]
    nrow2 = wnv.out[[2]]
    nrow3 = wnv.out[[3]]
    minmosqyear = wnv.out[[4]]
    numpos = wnv.out[[5]]
    perpos = wnv.out[[6]]
    
  # Additional mosquito data processing
  mdp2 = mosquito.data.process2(wnv, stratafile, mosqmatoutputdir, graphicoutputdir, weekinquestionSun, minmosqyear, maxmosqyear)
    wnv = mdp2[[1]]
    infectglm = mdp2[[2]]
    randeffs = mdp2[[3]]
    strata = mdp2[[4]]
    
  # Process human data 
  human.out = process.human.data(humandatafile, maxdesiredhumandate, maxobservedhumandate, minmosqyear, maxmosqyear, weekinquestionSun, districtshapefile, stratafile, graphicoutputdir)
    fullcasemat = human.out[[1]]
    projected_districts.df = human.out[[2]]
    totcase = human.out[[3]]
    observedfraction = human.out[[4]]
  
  # Estimate relationship between number of positive mosquito districts and expected number of human cases
  pdthc = positive.districts.to.human.cases(wnv, fullcasemat)
  slope = pdthc[[1]]
  intercept = pdthc[[2]]
    
  # Create mosquito plots by month
  mosq.by.month(wnv, infectglm, minmosqyear, maxmosqyear, compyear1, compyear2, graphicoutputdir, weekinquestionSun)
  
  # Second human data processing step
  phd2 = process.human.data2(dailyextr, fullcasemat, randeffs, var1name, var2name, strata, maxdesiredhumandate, minmosqyear, maxmosqyear, laglen, graphicoutputdir, compyear1, compyear2, weekinquestionSun)
    tempdf = phd2[[1]]
    tempdf2 = phd2[[2]]
    fullcasemat = phd2[[3]]

  # Write human case matrix
  cm.out = write.case.matrix.and.calc.cases(fullcasemat, fullcasematoutputdir, maxmosqyear)
    fullcasemat = cm.out[[1]]
    positivesthisyear = cm.out[[2]]
    multiplier = cm.out[[3]]
    casesthisyear = cm.out[[4]]

  # Create the proportion positive plot
  thisplot2 = proportion.positive.plot(tempdf,tempdf2, maxmosqyear, compyear1, compyear2, graphicoutputdir, weekinquestionSun)
  
  # Create the raw risk map
  make.raw.risk.map(districtshapefile, fullcasemat, graphicoutputdir, weekinquestionSun)

  # Create the relative risk map
  make.riskcalcs(fullcasemat, weekinquestion, projected_districts.df, graphicoutputdir, weekinquestionSun, thisplot2)

  # Calculate number of cases expected for the week
  n.districts = length(unique(wnv$district))
  Districts.With.Cases.Focal.Week = round(n.districts*tempdf2$est,1) # 66 was the number of counties in South Dakota. Replaced on 2019-07-17 with a variable
  #Districts.With.Cases.Focal.Week = round(multiplier*tempdf2$est,1) # formerly hardcoded to be 66, changed to multiplier which is the number of human cases expected, given a human case, and $est is the probability of a positive district #**# This option removed on 2019-07-17 after discussing with Justin and correcting 66 to be n.districts. This version does not multiply by number of counties, and tempdf2$est gives the average probability of a human case for a county. 
  if (original.metrics == 1){
    human.cases = casesthisyear # casesthisyear
  }
  if (original.metrics == 0){
    human.cases = round(positivesthisyear * slope + intercept, 1)
  }
  
  # compile the main results into the results object
  results = list(weeks.cases = Districts.With.Cases.Focal.Week, annual.positive.district.weeks = positivesthisyear, multiplier = multiplier, annual.human.cases = human.cases)

  return(results)
}


#' Read Weather Data (ArboMAP)
#'
#' Read in the ArboMAP formatted weather data. Data in county-year format for the RF1 model
#' should be added using the files.to.add entry in the rf1.inputs object.
#'
#' Data in this format can be downloaded from Google Earth Engine, see ArboMAP documentation
#' for instructions
#'
#' @details Written by J. Davis and M. Wimberly. Modified by A. Keyel. Original available at www.github.com/ecograph/ArboMAP. Used under GPL-3 license
#'
#' @param weatherpathstr A directory containing the weather data to be aggregated.
#' The weather data can contain multiple files, and these can overlap in date range,
#' as long as they all have the same columns.
#' The weather data files in this folder should have a district column, a doy (day of year) column
#' a year column, and then an entry for each independent variable of interest
#' @param weathersummaryfile If this file is in weatherpathstr, it will be ignored.
#' Note that no consolidated summary file is created.
#'
#' @export read.weather.data
read.weather.data = function(weatherpathstr, weathersummaryfile){
  # load and concat files
  weatherlist <- list.files(path=weatherpathstr, pattern="(.csv)", recursive=FALSE)
  weather <- data.frame()
  for (i in 1:length(weatherlist)) {
  
    if (weatherlist[[i]] != weathersummaryfile) {
  
      tempdf <- read.csv(paste(weatherpathstr, weatherlist[[i]], sep="")) 
      weather <- bind_rows(weather, tempdf)  
    
    }
  
  }
  weather$date <- as.Date(paste(weather$year,
                              weather$doy,
                              sep="-"),
                        "%Y-%j")

  # get rid of duplicated rows
  weather$districtdate <- paste(weather$district, weather$date)
  weather <- subset(weather, !duplicated(weather$districtdate))
  return(weather)
}


plot.weather.data = function(weather, graphicoutputdir, weekinquestionSun){
  # plot normals and this year
  weather <- group_by(weather,
                    doy)
  doymet  <- dplyr::summarize(weather,
                     med_tmeanc = quantile(tmeanc, probs=0.50, na.rm=TRUE),
                     med_rmean  = quantile(rmean,  probs=0.50, na.rm=TRUE),
                     med_vpd    = quantile(vpd,    probs=0.50, na.rm=TRUE),
                     max_tmeanc = max(tmeanc, na.rm=TRUE),
                     max_rmean  = max(rmean,  na.rm=TRUE),
                     max_vpd    = max(vpd,    na.rm=TRUE),
                     min_tmeanc = min(tmeanc, na.rm=TRUE),
                     min_rmean  = min(rmean,  na.rm=TRUE),
                     min_vpd    = min(vpd,    na.rm=TRUE))
  weather <- ungroup(weather)

  thisyear <- max(weather$year, na.rm=TRUE)
  thisyear <- subset(weather, year == thisyear)
  thisyear <- group_by(thisyear, doy)
  thisyear <- dplyr::summarize(thisyear,
                      med_tmeanc = quantile(tmeanc, probs=0.50, na.rm=TRUE),
                      med_rmean  = quantile(rmean,  probs=0.50, na.rm=TRUE),
                      med_vpd    = quantile(vpd,    probs=0.50, na.rm=TRUE),
                      max_tmeanc = max(tmeanc, na.rm=TRUE),
                      max_rmean  = max(rmean,  na.rm=TRUE),
                      max_vpd    = max(vpd,    na.rm=TRUE),
                      min_tmeanc = min(tmeanc, na.rm=TRUE),
                      min_rmean  = min(rmean,  na.rm=TRUE),
                      min_vpd    = min(vpd,    na.rm=TRUE))

  tempdf <- left_join(doymet, thisyear, by="doy")
  tempdf <- tempdf[!is.na(tempdf$med_tmeanc.x),]

  meantemp.x <- round(mean(tempdf$med_tmeanc.x, na.rm=TRUE), 1)
  meantemp.y <- round(mean(tempdf$med_tmeanc.y, na.rm=TRUE), 1)
  meanvapd.x <- round(mean(tempdf$med_vpd.x, na.rm=TRUE), 1)
  meanvapd.y <- round(mean(tempdf$med_vpd.y, na.rm=TRUE), 1)
  
  plot1 <- ggplot() + geom_line(data=doymet, aes(x=doy, y=med_tmeanc)) +
  geom_ribbon(data=doymet, aes(x=doy, ymin=min_tmeanc, ymax=max_tmeanc), alpha=0.3) +
  geom_line(data=thisyear, aes(x=doy, y=med_tmeanc), color="red", size=1) +
  xlab("Day of the year") + ylab(var1name) +
  ggtitle(paste(var1name, 
                max(weather$year, na.rm=TRUE),
                sep=" "))
  plot2 <- ggplot() + geom_line(data=doymet, aes(x=doy, y=log(med_vpd+1))) +
   geom_ribbon(data=doymet, aes(x=doy, ymin=log(min_vpd+1), ymax=log(max_vpd+1)), alpha=0.3) +
   geom_line(data=thisyear, aes(x=doy, y=log(med_vpd+1)), color="red", size=1) +
   ggtitle(paste(var2name,
                max(weather$year, na.rm=TRUE),
                sep=" ")) +
    xlab("Day of the year") + ylab(var2name)
  grid.arrange(plot1, plot2, nrow=2)

  # simplify the district names
  weather$district <- simplifynames(weather$district)

  # put aside this set for use in the regression
  dailyextr <- weather
  rm(weather)
  gc()

  save1 = ggsave(sprintf("%svar1_%s.png", graphicoutputdir, weekinquestionSun), plot1)
  save2 = ggsave(sprintf("%svar2_%s.png", graphicoutputdir, weekinquestionSun), plot2)

  return(dailyextr)
}

vector.infection.data = function(mosqfile, maxmosqyear){
  # Allow input of an R object instead of a .csv file
  if (typeof(mosqfile) == "character"){
    wnv <- read.csv(mosqfile, stringsAsFactors=FALSE)
  }else{  wnv = mosqfile  }
  wnv$col_date <- as.Date(wnv$col_date, "%m/%d/%Y")
  wnv$year <- as.numeric(format(wnv$col_date, "%Y"))
  
  # convert district to factor
  wnv$district <- simplifynames(wnv$district)
  wnv$district <- factor(wnv$district)

  # figure out how many rows we start with
  nrow1 <- nrow(wnv)

  # create some variables we can use to filter
  wnv$col_year <- as.numeric(format(wnv$col_date, "%Y"))
  wnv$doy      <- as.numeric(format(wnv$col_date, "%j"))
  wnv$weeknum  <- as.numeric(format(wnv$col_date, "%U"))
  wnv$species <- NULL
  wnv$district <- factor(wnv$district)
  wnv$district <- droplevels(wnv$district)

  # figure out which years we're modeling
  minmosqyear <- min(wnv$year, na.rm=TRUE)

  # get rid of those which don't have a result
  wnv <- wnv[which(!is.na(wnv$wnv_result)),]
  wnv <- wnv[which(!is.na(wnv$doy)),]

  # delete anything before a certain day
  wnv <- wnv[which(wnv$doy >= 100),]
  # delete anything after a certain day
  wnv <- wnv[which(wnv$doy <= 212),]

  # after cleaning, how many do we have?
  nrow2 <- nrow(wnv)
  nrow3 <- nrow(wnv[wnv$year == maxmosqyear,])

  tempdf <- wnv[wnv$year == maxmosqyear,]
  tempdf <- tempdf[!is.na(tempdf$wnv_result),]
  wnvdenominator <- nrow(tempdf)
  wnvnumerator <- nrow(tempdf[tempdf$wnv_result == 1,])

  numpos <- wnvnumerator
  perpos <- 100*round(wnvnumerator/wnvdenominator, 3)

  wnv.out = list(wnv, nrow2, nrow3, minmosqyear, numpos, perpos)  
}

mosquito.data.process2 = function(wnv, stratafile, mosqmatoutputdir, graphicoutputdir, weekinquestionSun, minmosqyear, maxmosqyear){
  # import district identifiers
  strata <- read.csv(stratafile)
  strata$district <- simplifynames(strata$district)
  strata <- strata[c("district", "strata")]
  wnv <- merge(x=wnv, y=strata,
             by.x="district",
             by.y="district",
             all.x=TRUE)

  # figure out how many distinct districts we have left
  districtlist           <- data.frame(district=unique(wnv$district))
  distinctdistricts      <- length(unique(wnv$district))
  districtlist$districtnum <- seq(from=1, to=distinctdistricts, by=1)
  wnv <- merge(x=wnv, y=districtlist,
             by="district",
             all=TRUE)

  # create a variable that at least has a little chance of being orthogonal to 1.
  wnv$dminus <- wnv$doy - mean(wnv$doy, na.rm=TRUE)

  # make sure all the observations have a stratum and year
  wnv <- wnv[!is.na(wnv$strata),]
  wnv <- wnv[!is.na(wnv$year),]

  # run a random effect model on orthogonalized data
  infectglm <- glmer(wnv_result ~ 1+dminus+
                    (0+1|year) +
                    (0+dminus|year) +
                    (0+1|strata:year) + 
                    (0+dminus|strata:year),
                  family=binomial(),
                  data=wnv)

  wnv$est <- predict(infectglm, newdata=wnv, type="response")
  write.csv(x=wnv, file=paste(mosqmatoutputdir,
                            "mosqmatrix.csv",
                            sep=""))

  # predict random effects for all years
  # dminus cannot be set to simply 0 - we have to make sure this isn't read in as a factor
  randeffs <- expand.grid(strata=unique(wnv$strata),
                        year=minmosqyear:maxmosqyear,
                        dminus=0.00001)
  randeffs <- randeffs[which(!is.na(randeffs$strata)),]
  # if you don't allow new levels, the most recent year might not have an estimate
  randeffs$mosqinfect <- predict(infectglm, newdata=randeffs, allow.new.levels = TRUE)
  randeffs$stratayear <- paste(randeffs$strata, randeffs$year, sep=":")
  randeffs$adjmosqinfect <- randeffs$mosqinfect - mean(randeffs$mosqinfect,na.rm=TRUE)
  randeffs$stratum <- factor(randeffs$strata)
  thisplot <- ggplot(randeffs) + geom_line(aes(x=year, y=adjmosqinfect, group=stratum, color=stratum)) +
    geom_abline(slope=0, intercept=0, linetype=2) +
    scale_x_continuous(breaks=minmosqyear:maxmosqyear) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
        panel.grid.minor=element_blank(),
        panel.grid.major.y=element_blank())+
    xlab("") + ylab("Relative risk due to\nvector infection growth rate")
  
  plot(thisplot)
  ggsave(sprintf("%svector infection rates_%s.png", graphicoutputdir, weekinquestionSun), thisplot)
  
  return(list(wnv, infectglm, randeffs, strata))
}

process.human.data = function(humandatafile, maxdesiredhumandate, maxobservedhumandate, minmosqyear, maxmosqyear, weekinquestionSun, districtshapefile, stratafile, graphicoutputdir){
  
  # import data
  if (typeof(humandatafile) == "character"){
    human <- read.csv(humandatafile)
  # If it is already read in, just use the data object
  }else{  human = humandatafile  }

  begrow <- nrow(human)
  human$chardate <- as.character(levels(human$date))[as.numeric(human$date)]
  human$date <- as.Date(human$date, "%m/%d/%Y")
  human$creationyear <- as.numeric(format(human$date, "%Y"))
  human$creationmonth <- as.numeric(format(human$date, "%m"))

  # simplify district names
  human$district <- simplifynames(human$district)
  human$district <- factor(human$district)
  human$doy <- as.numeric(format(human$date, "%j"))

  # retain only those in the right date range
  hWNVminyear <- minmosqyear
  hWNVmaxyear <- maxmosqyear-1
  human <- human[which((human$creationyear >= hWNVminyear) & (human$creationyear <= hWNVmaxyear)),]

  # retain only those in months 5-11
  human <- human[which((human$creationmonth >= 5) & (human$creationmonth <= 11)),]

  # create the full list of weeks
  minobservedhumandate <- min(human$date)
  
  # set up the data frame so that it ends at the maxdesiredhumandate 
  filledweeks <- seq(from=maxdesiredhumandate, to=minobservedhumandate, by=-7)
  
  fullcasemat <- expand.grid(sort(unique(human$district)), filledweeks)
  names(fullcasemat) <- c("district", "weekstartdate")
  head(fullcasemat)
  fullcasemat$anycases   <- rep(0, nrow(fullcasemat))
  fullcasemat$totalcases <- rep(0, nrow(fullcasemat))
  fullcasemat$observed   <- 1*((fullcasemat$weekstartdate >= minobservedhumandate)&
                               (fullcasemat$weekstartdate <= maxobservedhumandate))
  for (i in 1:nrow(fullcasemat)) {
    
    thisweekstartdate <- fullcasemat$weekstartdate[i]
    thisdistrict        <- fullcasemat$district[i]
    
    tempcases <- human[which(human$district == thisdistrict),]
    tempcases <- tempcases[which(tempcases$date >= thisweekstartdate),]
    tempcases <- tempcases[which(tempcases$date <= (thisweekstartdate + 6)),]
    
    fullcasemat$anycases[i] <- 1*(nrow(tempcases) > 0)
    if (nrow(tempcases) > 0) {
      
      fullcasemat$totalcases[i] <- nrow(tempcases)
      
    }
    
  }
  totcase <- sum(fullcasemat$totalcases, na.rm=TRUE)
  anypos  <- sum(fullcasemat$anycases, na.rm=TRUE)
  

  # figure out what percentage of cases we'll likely have seen before the start of this week
  weekinquestionSundoy <- as.numeric(format(weekinquestionSun, "%j"))
  fullcasemat$doy      <- as.numeric(format(fullcasemat$weekstartdate, "%j"))
  
  tempdf <- fullcasemat[fullcasemat$doy < (weekinquestionSundoy+7),]
  observedbefore <- sum(tempdf$totalcases, na.rm=TRUE)
  observedtotal  <- sum(fullcasemat$totalcases, na.rm=TRUE)
  weekinquestionreformat <- format(weekinquestionSun, "%m-%d")
  observedfraction <- 100*round(observedbefore / observedtotal, 2)

  district_shapes <- readShapePoly(districtshapefile)
  # simplify name
  district_shapes$district <- simplifynames(district_shapes$NAME)
  crs(district_shapes) <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80     +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  projected_districts <- spTransform(district_shapes, crs("+proj=longlat +datum=WGS84 +no_defs"))
  projected_districts@data$id = rownames(projected_districts@data)
  projected_districts.df <- tidy(projected_districts)
  projected_districts.df <- left_join(projected_districts.df, projected_districts@data, by="id")
  
  stratamapcsv <- read.csv(stratafile)
  stratamapcsv$district <- simplifynames(stratamapcsv$district)
  stratamapcsv <- stratamapcsv[c("district", "strata")]
  
  projected_districts.df <- left_join(projected_districts.df, stratamapcsv,
                                     by="district")
  
  if (length(projected_districts.df$strata) == 0){
    projected_districts.df$stratum <- factor(projected_districts.df$strata.x)
  }else{
    projected_districts.df$stratum = factor(projected_districts.df$strata)
  }

  thisplot <- ggplot(projected_districts.df) +
        aes(long,lat,fill=stratum,group=group,id=id,guides=FALSE) +
        geom_polygon() + xlab("") + ylab("") +
        geom_path(color="black") +
        theme(legend.position="bottom") +
        coord_map() + ggtitle("State stratification map") +
          theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              axis.ticks = element_blank(), axis.text.x = element_blank(),
              axis.text.y = element_blank(), axis.title.x=element_blank(),
              axis.title.y = element_blank(),
              legend.position="right",
              legend.key.width=unit(1, "cm"),
              legend.key.height=unit(0.5,"cm"))
  plot(thisplot)
  
  ggsave(sprintf("%sstrata_%s.png", graphicoutputdir, weekinquestionSun), thisplot)
  
  human.out = list(fullcasemat, projected_districts.df, totcase, observedfraction)
  return(human.out)
}

# Calculate the relationship between positive districts to human cases in that year
positive.districts.to.human.cases = function(wnv, fullcasemat){
  
  # Calculate number of WNV+ districts in historical data
  wnv.by.week.district.year = aggregate(wnv$wnv_result, by = list(wnv$weeknum, wnv$district, wnv$year), max, na.rm = TRUE)
  wnv.positives = wnv.by.week.district.year[wnv.by.week.district.year$x > 0, ] # Just keep those with a positive
  # Add together all positive weeks for each district and year
  positive.district.weeks = aggregate(wnv.positives$x, by = list(wnv.positives$Group.2, wnv.positives$Group.3), sum) # Using x as that will give the appropriate sums
  colnames(positive.district.weeks) = c("district", "year", "PositiveDistrictWeeks")
  # Group.1 is district, Group.2 is year
  positive.district.weeks$district_year = sprintf("%s_%s", positive.district.weeks$district, positive.district.weeks$year)
  
  # Calculate number of human cases for each year of historical data
  fullcasemat$year = strftime(fullcasemat$weekstartdate, "%Y")
  district.cases = aggregate(fullcasemat$totalcases, by = list(fullcasemat$district, fullcasemat$year), sum)
  colnames(district.cases) = c("district", "year", "HumanCases")
  district.cases$district_year = sprintf("%s_%s", district.cases$district, district.cases$year)

  # Merge based on district and year
  annual.positives = merge(positive.district.weeks, district.cases, by = "district_year") #, all = TRUE #**# No point in keeping all, to only drop them in the next step.
  
  #**# THIS MISSES HUMAN CASES THAT OCCURRED WITHOUT A POSITIVE DISTRICT WEEK BUT AVOIDS ZERO-INFLATION
  # Run simple regression to get relationship between district-week positives and human cases
  regression = lm(annual.positives$HumanCases~annual.positives$PositiveDistrictWeeks) #**# May want something more sophisticated: i.e. year and district as random effects
  
  slope = regression$coefficients[2]
  intercept = regression$coefficients[1]
  
  reg = summary(regression)
  intercept.p = reg$coefficients[1,4]
  slope.p = reg$coefficients[2,4]
  # Return slope estimate #**# Can add quantile regression or some other approach to better incorporate variation across years
  return(list(slope, intercept, slope.p, intercept.p))
  
}

mosq.by.month = function(wnv, infectglm, minmosqyear, maxmosqyear, compyear1, compyear2, graphicoutputdir, weekinquestionSun){
  
  mosqmopreds <- expand.grid(strata=unique(wnv$strata),
                             year=minmosqyear:maxmosqyear,
                             doy=seq(from=min(wnv$doy, na.rm=TRUE),
                                     to  =max(wnv$doy, na.rm=TRUE), by=1))
  mosqmopreds$dminus <- mosqmopreds$doy - mean(wnv$doy, na.rm=TRUE)
  mosqmopreds$preds <- predict(infectglm, newdata=mosqmopreds, type="response", allow.new.levels=TRUE)
  
  mosqmopreds <- group_by(mosqmopreds, dminus, year)
  mosqmopreds <- dplyr::summarize(mosqmopreds,
                           preds=mean(preds, na.rm=TRUE),
                           doy=mean(doy, na.rm=TRUE))
  thisyear1   <- mosqmopreds[mosqmopreds$year == maxmosqyear,]
  comparison1 <- mosqmopreds[mosqmopreds$year == compyear1,]
  comparison2 <- mosqmopreds[mosqmopreds$year == compyear2,]
  
  thisyeardot <- wnv[wnv$year == maxmosqyear,]
  
  if (sum(!is.na(thisyeardot$doy)) > 20) {
  
    thisyeardot$rounddoy <- cut2(thisyeardot$doy, g=6)
    thisyeardot <- group_by(thisyeardot, rounddoy)
    thisyeardot <- dplyr::summarize(thisyeardot,
                          meanpos = mean(wnv_result, na.rm=TRUE),
                          meandoy = mean(doy, na.rm=TRUE))
    thisplot <- ggplot(data=mosqmopreds) + geom_line(data=mosqmopreds, aes(x=doy, y=preds, group=year),
                                                     color="grey", alpha=0.5) +
      geom_line(data=comparison1, aes(x=doy, y=preds), color="blue") +
      geom_line(data=comparison2, aes(x=doy, y=preds), color="blue", linetype=2) +
      geom_line(data=thisyear1, aes(x=doy, y=preds), color="red") +
      geom_point(data=thisyeardot, aes(x=meandoy, y=meanpos), color="red") +
      xlab("day of the year") + ylab("Vector pool positive rate") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"))
    
  } else {
    
    thisplot <- ggplot(data=mosqmopreds) + geom_line(data=mosqmopreds, aes(x=doy, y=preds, group=year),
                                                     color="grey", alpha=0.5) +
      geom_line(data=comparison1, aes(x=doy, y=preds), color="blue") +
      geom_line(data=comparison2, aes(x=doy, y=preds), color="blue", linetype=2) +
      geom_line(data=thisyear1, aes(x=doy, y=preds), color="red") +
      xlab("day of the year") + ylab("Vector pool positive rate") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"))
    
  }
  
  
  plot(thisplot)
  ggsave(sprintf("%smosqinfectgrowthrates_%s.png", graphicoutputdir, weekinquestionSun), thisplot)

}

process.human.data2 = function(dailyextr, fullcasemat, randeffs, var1name, var2name, strata, maxdesiredhumandate, minmosqyear, maxmosqyear, laglen, graphicoutputdir, compyear1, compyear2, weekinquestionSun){
  
  # select which variables we're going to use!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dailyextr <- data.frame(dailyextr)
  dailyextr$var1 <- dailyextr[,which(colnames(dailyextr) == var1name)]
  dailyextr$var2 <- dailyextr[,which(colnames(dailyextr) == var2name)]
  dailyextr <- dailyextr[c("district","date","var1","var2")]
  dailyextr$district <- simplifynames(dailyextr$district)
  # fill in missing and future climatological data
  totaldailyextr <- expand.grid(district=unique(dailyextr$district),
                                date=seq(from=as.Date("2003-01-01", "%Y-%m-%d"),
                                         to=maxdesiredhumandate,
                                         by=1))
  dailyextr <- merge(x=totaldailyextr, y=dailyextr,
                     by.x=c("district","date"),
                     by.y=c("district","date"),
                     all.x=TRUE)
  dailyextr$doy <- as.numeric(format(dailyextr$date, "%j"))
  
  dailyextr <- group_by(dailyextr, district, doy)
  districtdoymean <- dplyr::summarize(dailyextr,
                             meanvar1=mean(var1, na.rm=TRUE),
                             meanvar2=mean(var2, na.rm=TRUE))
  dailyextr <- ungroup(dailyextr)                           
  
  dailyextr <- left_join(dailyextr, districtdoymean,
                     by=c("district","doy"))
  
  dailyextr$var1[is.na(dailyextr$var1)] <- dailyextr$meanvar1[is.na(dailyextr$var1)]
  dailyextr$var2[is.na(dailyextr$var2)] <- dailyextr$meanvar2[is.na(dailyextr$var2)]
  
  # garbage collection
  rm(districtdoymean)
  dailyextr$meanvar1 <- NULL
  dailyextr$meanvar2 <- NULL
  gc()
  
  # we need the districts to be in lower case to merge with the gridMET data
  fullcasemat$district <- simplifynames(fullcasemat$district)
  
  datalagger <- expand.grid(unique(fullcasemat$district),
                            unique(fullcasemat$weekstartdate),
                            seq(from=0, to=laglen-1, by=1))
  names(datalagger) <- c("district","date","lag")
  datalagger$laggeddate <- datalagger$date-datalagger$lag
  
  datalagger <- left_join(datalagger, dailyextr,
                          by=c("district"="district",
                               "laggeddate"="date"))
  
  # garbage collection
  rm(dailyextr)
  gc()
  
  # pivot
  mean1data <- dcast(datalagger, district + date ~ lag, value.var="var1")
  names(mean1data) <- paste("var1_",names(mean1data),sep="")
  head(mean1data)
  
  mean2data <- dcast(datalagger, district + date ~ lag, value.var="var2")
  names(mean2data) <- paste("var2_",names(mean2data),sep="")
  head(mean2data)
  
  # garbage collection
  rm(datalagger)
  gc()
  
  # and put all this lagged info back into the total case matrix
  fullcasemat <- merge(x=fullcasemat, y=mean1data,
                       by.x=c("district", "weekstartdate"),
                       by.y=c("var1_district","var1_date"),
                       all.x=TRUE)
  fullcasemat <- merge(x=fullcasemat, y=mean2data,
                       by.x=c("district", "weekstartdate"),
                       by.y=c("var2_district","var2_date"),
                       all.x=TRUE)
  
  # garbage collection
  rm(mean1data)
  rm(mean2data)
  gc()
  
  lagframe <- data.frame(x=seq(from=1, to=laglen, by=1))
  alpha <- 1/4
  distlagfunc <- ns(lagframe$x, intercept=TRUE,
                    knots=quantile(lagframe$x,
                                   probs=seq(from=alpha, to=1-alpha, by=alpha),
                                   na.rm=TRUE))
  #matplot(distlagfunc, type="l") #**# Commented out to suppress the plot
  dlagdeg <- size(distlagfunc)[2]
  
  # create actual distributed lag summaries
  band1summaries <- matrix(rep(0, nrow(fullcasemat)*dlagdeg), nrow(fullcasemat), dlagdeg)
  band2summaries <- matrix(rep(0, nrow(fullcasemat)*dlagdeg), nrow(fullcasemat), dlagdeg)
  
  names(fullcasemat)
  min1index <- which(names(fullcasemat) == "var1_0")
  min2index <- which(names(fullcasemat) == "var2_0")
  fullcasemat <- data.frame(fullcasemat)
  band1temp <- as.matrix(fullcasemat[,(min1index:(min1index+laglen-1))])
  band2temp <- as.matrix(fullcasemat[,(min2index:(min2index+laglen-1))])
  for (j in 1:dlagdeg) {
    
    band1summaries[,j] <- band1temp %*% distlagfunc[,j]
    band2summaries[,j] <- band2temp %*% distlagfunc[,j]
    
  }
  
  fullcasemat$band1summaries <- band1summaries
  fullcasemat$band2summaries <- band2summaries
  
  # garbage collection
  rm(band1summaries)
  rm(band2summaries)
  gc()
  
  # include the mosquito summary statistic
  stratayearmosq <- randeffs
  stratayearmosq <- stratayearmosq[c("stratayear","mosqinfect")]
  names(stratayearmosq) <- c("stratayear","MIRsummarystat")
  
  # import district identifiers
  fullcasemat <- merge(x=fullcasemat, y=strata,
                       by.x="district", by.y="district",
                       all.x=TRUE)
  fullcasemat$year <- as.numeric(format(fullcasemat$weekstartdate, "%Y"))
  fullcasemat$stratayear <- paste(fullcasemat$strata, fullcasemat$year, sep=":")
  fullcasemat$year <- NULL
  
  fullcasemat <- merge(x=fullcasemat, y=stratayearmosq,
                       by.x="stratayear", by.y="stratayear",
                       all.x=TRUE)
  
  # calculate missing MIRsummarystat
  stratayearmosq$strata <- substr(stratayearmosq$stratayear, 1, 3)
  tempdf <- aggregate(stratayearmosq$MIRsummarystat,
                      by=list(stratayearmosq$strata),
                      FUN=mean,
                      na.rm=TRUE)
  names(tempdf) <- c("strata","MIRsummarystat2")
  
  fullcasemat <- merge(x=fullcasemat, y=tempdf,
                       by.x="strata", by.y="strata",
                       all.x=TRUE)
  fullcasemat$MIRsummarystat[is.na(fullcasemat$MIRsummarystat)] <- fullcasemat$MIRsummarystat2[is.na(fullcasemat$MIRsummarystat)]
  
  fullcasemat$strata <- NULL
  fullcasemat$stratayear <- NULL

  # get rid of districts that never have a positive case in the time period
  fullcasemat$doy <- as.numeric(format(fullcasemat$weekstartdate, "%j"))
  firstreg  <- glm(anycases ~ 0+MIRsummarystat+district+
                     band1summaries+
                     band2summaries+
                     band1summaries:doy+
                     band2summaries:doy,
                   family=binomial(), data=fullcasemat,
                   subset=observed==1)
  
  summary(firstreg)
  
  
  # show predictions
  preds <- predict(firstreg, type="link", newdata=fullcasemat, se.fit=TRUE, newlevels=TRUE)
  fullcasemat$est <- predict(firstreg, type="response", newdata=fullcasemat, newlevels=TRUE)
  
  preds <- as.matrix(cbind(preds$fit, preds$se.fit))
  preds <- data.frame(preds)
  names(preds) <- c("est", "sd")
  preds$upper <- preds$est + 1.96*preds$sd
  preds$lower <- preds$est - 1.96*preds$sd
  preds$est   <- preds$est
  
  preds$upper <- 1 / (1 + exp(-preds$upper))
  preds$lower <- 1 / (1 + exp(-preds$lower))
  preds$est <- 1 / (1 + exp(-preds$est))
  
  
  tempfull <- fullcasemat
  tempfull$preds <- preds$est
  
  tempfull$month <- as.numeric(format(tempfull$weekstartdate, "%m"))
  tempfull <- tempfull[which(tempfull$month >= 6),]
  tempfull <- tempfull[which(tempfull$month <= 10),]
  
  tempdf <- data.frame(weekstartdate=fullcasemat$weekstartdate,
                       obs=fullcasemat$anycases,
                       est=preds$est)
  tempdf <- group_by(tempdf, weekstartdate)
  tempdf <- dplyr::summarize(tempdf,
                      obs=mean(obs, na.rm=TRUE),
                      est=mean(est, na.rm=TRUE))
  tempdf <- ungroup(tempdf)
  tempdf2 <- tempdf
  tempdf2$week <- as.numeric(format(tempdf2$weekstartdate, "%U"))
  lowe <- 12
  hiwe <- 52-12
  
  tempdf2 <- tempdf2[tempdf2$week >= lowe,]
  tempdf2 <- tempdf2[tempdf2$week <= hiwe,]
  
  tempdf2$newwe <- (tempdf2$week - lowe) / (hiwe - lowe + 1)
  tempdf2$year <- as.numeric(format(tempdf2$weekstartdate, "%Y"))
  tempdf2$newdate <- tempdf2$year + tempdf2$newwe
  thisplot <- ggplot(tempdf2) + geom_line(aes(x=newdate,
                                 y=obs)) +
    geom_line(aes(x=newdate, y=est), color="red") +
    ggtitle("Statewide model predictions") +
    xlab("") + ylab("") +
    scale_x_continuous(breaks=seq(from=minmosqyear,to=(maxmosqyear+1),by=1),
                       labels=c(seq(from=minmosqyear, to=maxmosqyear, by=1), ""),
                       limits=c((minmosqyear+0.125),(maxmosqyear+1))) +
    ylab("Proportion of districts positive") +
    theme(panel.grid.minor.x=element_blank()) +
    theme(axis.text.x = element_text(hjust=-0.225))
  plot(thisplot)
  ggsave(sprintf("%spredictions_%s.png", graphicoutputdir, weekinquestionSun), thisplot)
  
  tempdf <- fullcasemat
  tempdf <- group_by(tempdf, weekstartdate)
  tempdf <- dplyr::summarize(tempdf,
                      obs=mean(anycases, na.rm=TRUE),
                      est=mean(est, na.rm=TRUE))
  tempdf$year <- as.numeric(format(tempdf$weekstartdate, "%Y"))
  tempdf$month <- as.numeric(format(tempdf$weekstartdate, "%m"))
  
  tempdf <- tempdf[tempdf$month >= 5,]
  tempdf <- tempdf[tempdf$month <= 11,]
  
  thisyear <- tempdf[tempdf$year == maxmosqyear,]
  comparison1 <- tempdf[tempdf$year == compyear1,]
  comparison2 <- tempdf[tempdf$year == compyear2,]
  
  comparison1$weekstartdate <- comparison1$weekstartdate + (maxmosqyear - compyear1)*365
  comparison2$weekstartdate <- comparison2$weekstartdate + (maxmosqyear - compyear2)*365
  
  tempdf <- bind_rows(thisyear,
                      comparison1,
                      comparison2)
  
  tempdf$est[tempdf$year == compyear1] <- tempdf$obs[tempdf$year == compyear1]
  tempdf$est[tempdf$year == compyear2] <- tempdf$obs[tempdf$year == compyear2]
  tempdf$year <- factor(tempdf$year)
  
  tempdf2 <- thisyear[thisyear$weekstartdate >= weekinquestionSun,]
  tempdf2 <- tempdf2[1,]
  
  return(list(tempdf, tempdf2, fullcasemat))
}

write.case.matrix.and.calc.cases = function(fullcasemat, fullcasematoutputdir, maxmosqyear){
  write.csv(x=fullcasemat[c("district","weekstartdate",
                            "anycases","totalcases", "est")],
            file=paste(fullcasematoutputdir, "case matrix.csv", sep=""))
  fullcasemat$year <- as.numeric(format(fullcasemat$weekstartdate, "%Y"))
  thisyear <- fullcasemat[fullcasemat$year == maxmosqyear,]
  previousyears <- fullcasemat[fullcasemat$year < maxmosqyear,]
  sumtotal <- sum(previousyears$totalcases, na.rm=TRUE)
  sumany   <- sum(previousyears$anycases, na.rm=TRUE)
  
  multiplier <- round(sumtotal / sumany, 2)
  positivesthisyear <- round(sum(thisyear$est, na.rm=TRUE), 1)
  casesthisyear     <- round(positivesthisyear * multiplier, 1)
  
  return(list(fullcasemat, positivesthisyear, multiplier, casesthisyear))
}

proportion.positive.plot = function(tempdf, tempdf2, maxmosqyear, compyear1, compyear2, graphicoutputdir, weekinquestionSun){

  thisplot <- ggplot(tempdf) + geom_line(data=tempdf, aes(x=weekstartdate,
                             y=est,
                             color=year,
                             linetype=year)) +
    geom_point(data=tempdf2, aes(x=weekstartdate, y=est), color="red", size=4)+
    xlab("") + ylab("") +
    scale_x_date(date_labels="%b", date_breaks="1 month",
                 limits=c(as.Date(paste(maxmosqyear, "-05-01", sep=""), "%Y-%m-%d"),
                          as.Date(paste(maxmosqyear, "-11-15", sep=""), "%Y-%m-%d")),
                 labels=c("May","Jun","Jul","Aug","Sep","Oct","Nov","")) +
    ylab("Proportion of districts positive") +
    theme(panel.grid.minor.x=element_blank()) +
    scale_linetype_manual(values=c(3,2,1)) +
    scale_color_manual(values=c("black","black","red"))
  thisplot1 <- thisplot + ggtitle(paste("Estimates for ",
                  maxmosqyear, " compared to observations in ", compyear1, " and ", compyear2,
                  "\nwith week beginning ", weekinquestionSun, " highlighted", sep="")) +
    theme(axis.text.x = element_text(hjust=-1.25))
  plot(thisplot1)
  thisplot2 <- thisplot + theme(legend.position="bottom") +
    theme(axis.text.x = element_text(hjust=-0.50))
  ggsave(sprintf("%sestimates_%s.png", graphicoutputdir, weekinquestionSun), thisplot2,
         width=4,
         height=4)
  
  return(thisplot2)  
}

make.raw.risk.map = function(districtshapefile, fullcasemat, graphicoutputdir, weekinquestionSun){
  district_shapes <- readShapePoly(districtshapefile)
  district_shapes$NAME <- simplifynames(district_shapes$NAME)
  crs(district_shapes) <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80     +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  projected_districts <- spTransform(district_shapes, crs("+proj=longlat +datum=WGS84 +no_defs"))
  projected_districts@data$id = rownames(projected_districts@data)
  projected_districts.df <- tidy(projected_districts)
  projected_districts.df <- left_join(projected_districts.df, projected_districts@data, by="id")
  projected_districts.df$district <- simplifynames(projected_districts.df$NAME)
  
  # get this week's predictions
  thisweek <- fullcasemat[fullcasemat$weekstartdate >= weekinquestion,]
  thisweek <- thisweek[thisweek$weekstartdate == min(thisweek$weekstartdate, na.rm=TRUE),]
  thisweek <- thisweek[c("district", "est")]
  
  projected_districts.df <- left_join(projected_districts.df, thisweek,
                                     by="district")
  projected_districts.df$est[is.na(projected_districts.df$est)] <- min(projected_districts.df$est, na.rm=TRUE)
  #```
  #```{r shapefile2, include=TRUE, echo=FALSE, fig.width=7, fig.height=3.5}  
  
  thisplot <- ggplot(projected_districts.df) +
        aes(long,lat,fill=est,group=group,id=id,guides=FALSE) +
        geom_polygon() + xlab("") + ylab("") +
        geom_path(color="black") +
        theme(legend.position="bottom") +
        scale_fill_distiller(palette = "Spectral", limits=c(0,1),
                                                   breaks = c(0,1),
                                                   labels = c("Will definitely not\nreport any cases",
                                                              "Will definitely\nreport some cases"),
                                                       name = "") +
        coord_map() + ggtitle(paste("Estimate for week beginning ",
                              weekinquestionSun,
                              sep="")) +
          theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              axis.ticks = element_blank(), axis.text.x = element_blank(),
              axis.text.y = element_blank(), axis.title.x=element_blank(),
              axis.title.y = element_blank(),
              legend.position="bottom",
              legend.key.width=unit(1, "cm"),
              legend.key.height=unit(0.25,"cm"))
  plot(thisplot)
  #```
  #```{r silentsave89712, include=FALSE, echo=FALSE}
  ggsave(sprintf("%smap absolute_%s.png", graphicoutputdir, weekinquestionSun), thisplot)
}

make.riskcalcs = function(fullcasemat, weekinquestion, projected_districts.df, graphicoutputdir, weekinquestionSun, thisplot2){
  
  thisweek <- fullcasemat[fullcasemat$weekstartdate >= weekinquestion,]
  thisweek <- thisweek[thisweek$weekstartdate == min(thisweek$weekstartdate, na.rm=TRUE),]
  thisweek <- thisweek[c("district", "est", "doy")]
  
  approxdoy <- thisweek$doy[1]
  
  fullcasemat$year <- as.numeric(format(fullcasemat$weekstartdate, "%Y"))
  
  doypreds <- data.frame()
  for (curyear in unique(fullcasemat$year)) {
    
    for (curdistrict in unique(fullcasemat$district)) {
      
      thisdf <- fullcasemat[fullcasemat$year == curyear,]
      thisdf <- thisdf[thisdf$district == curdistrict,]
      
      if(sum(!is.na(thisdf$est)) > 1) {
      
        tempdf <- data.frame(district = curdistrict,
                             year = curyear,
                             est = approx(x=thisdf$doy,
                                          y=thisdf$est,
                                          xout=approxdoy)$y)
        
        doypreds <- bind_rows(doypreds, tempdf)
        
      }
      
    }
    
  }
  
  doypreds2 <- data.frame()
  for (curdistrict in unique(doypreds$district)) {
  
    tempdf <- doypreds[doypreds$district == curdistrict,]
    tempdf$percentile <- rank(tempdf$est, ties.method="random") / length(tempdf$est)
    
    doypreds2 <- bind_rows(doypreds2, tempdf)
    
  }
  doypreds <- doypreds2
  rm(doypreds2)
  
  riskalpha <- 0.25
  doypreds$riskcategory <- " About average    "
  doypreds$riskcategory[doypreds$percentile <= riskalpha/2] <- " Lower than usual    "
  doypreds$riskcategory[doypreds$percentile >= 1-(riskalpha/2)] <- " Higher than usual    "
  doypredscurrent <- doypreds[doypreds$year == as.numeric(format(weekinquestion, "%Y")),]
  projected_districts.df <- left_join(projected_districts.df, doypredscurrent,
                                     by="district")
  projected_districts.df$riskcategory[is.na(projected_districts.df$riskcategory)] <- " Not able to model    "
  thisplotrisk <- ggplot(projected_districts.df) +
        aes(long,lat,fill=riskcategory,group=group,guides=FALSE) +
        geom_polygon() + 
        geom_path(color="black") +
        coord_map() + 
        ggtitle("") +
        theme(text=element_text(size=15)) +
        theme(plot.title=element_text(size=15),
              legend.key=element_rect(fill="white")) +
        theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              axis.ticks = element_blank(), axis.text.x = element_blank(),
              axis.text.y = element_blank(), axis.title.x=element_blank(),
              axis.title.y = element_blank(),
              legend.position="bottom",
              legend.key.width=unit(1, "cm"),
              legend.key.height=unit(0.5,"cm")) +
        theme(legend.direction="vertical") +
        scale_fill_manual(name=paste("Risk for ", weekinquestionSun, " to ",
                                     format(weekinquestionSat, "%m-%d"), sep=""),
                          values=c(" Lower than usual    "="lightblue",
                                   " About average    "="yellow",
                                   " Higher than usual    "="red",
                                   " Not able to model    "="grey"))
  
  
  plot(thisplotrisk)
  #```
  #```{r silentsave1245, include=FALSE, echo=FALSE}
  ggsave(sprintf("%smap relative_%s.png", graphicoutputdir, weekinquestionSun), thisplotrisk,
         width=4,
         height=4)
  
  #**# COMMENTED OUT TO PREVENT sitegraphic from showing up in output. This also prevents the creation of this figure
  #**# Better to just make this a separate function in its own function with include = FALSE
  #**# But commenting it out was quicker.
  # sitegraphic <- grid.arrange(thisplotrisk +
  #                             theme(legend.position="right",
  #                                   legend.direction="vertical") +
  #                             scale_fill_manual(name="Risk this week",
  #                             values=c(" Lower than usual    "="lightblue",
  #                                      " About average    "="yellow",
  #                                      " Higher than usual    "="red",
  #                                      " Not able to model    "="grey")) +
  #                             #ylab("hidden axis") +
  #                             #theme(axis.title.y = element_text(colour = "white")) +
  #                             theme(legend.key.width=unit(0.25, "cm"),
  #                             legend.key.height=unit(0.50,"cm"),
  #                             text=element_text(size=10)) + 
  #                             theme(plot.margin = unit(c(0,0,0,0.4), "cm")),
  #   
  #                             # t r b l
  #                             
  #                             thisplot2 +
  #                             theme(legend.position="right") + ylab("") +
  #                             theme(axis.text.x=element_text(family="mono"))+
  #                             theme(panel.grid.minor=element_blank(),
  #                                   axis.text.y=element_blank(),
  #                                   axis.ticks=element_blank())+
  #                             theme(axis.text.x = element_text(hjust=-0.35)) +
  #                             theme(plot.margin = unit(c(0.1,0.2,0.1,-0.1), "cm")), nrow=2) 
  # 
  # ggsave(paste(graphicoutputdir,
  #              "site graphic.png",
  #              sep=""), sitegraphic,
  #        width=4,
  #        height=4)

}



