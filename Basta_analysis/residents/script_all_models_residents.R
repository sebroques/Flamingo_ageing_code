dir <- "D:/Senescence_analysis/Bird_data/Flamingo/Analyses_per_migration_strategy/residents"
setwd(dir)

#==================================================================================================================================================================================
#==================================================================================================================================================================================
#Age-dependent survival in Phoenicopterus roseus
#Data: capture-recapture, Tour-du-Valat, adult observations only
#Migration strategy: resident
# 
#==================================================================================================================================================================================
#==================================================================================================================================================================================

#==========================================================================================================================
# Load package

library(BaSTA)
library(snowfall)

#==========================================================================================================================
# Import dataset & check the data

Flamingo_adult<-read.csv2("Flamingo_adult_resident_1981.csv", dec=",", h=T)
head(Flamingo_adult)
datacheck<-DataCheck(Flamingo_adult,studyStart=1981, studyEnd=2020,autofix=rep(1,7), silent=FALSE)

#=========================================================================================================================
# Run Gompertz model, simple function / 1981-2020 / time-dependent recapture
Gompertz_simple<-basta(object=datacheck$newData,studyStart=1981,studyEnd=2020, model="GO",shape="simple", 
                              recaptTrans=c(1981, 1982, 1983, 1984, 1985, 1986, 1987,
                              1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998,
                              1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009,
                              2010, 2011,2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020), 
                              niter = 500000, burnin= 5000, thinn = 50, nsim = 4, 
                              parallel = TRUE, ncpus = 4, minAge = 1, lifeTable=TRUE)
# Generate some useful plots
plot(Gompertz_simple)
plot(Gompertz_simple, plot.trace = FALSE)
plot(Gompertz_simple, fancy = TRUE)

# Export outputs
write.csv2(Gompertz_simple$params, "00_Gompertz_simple_resident.csv")
write.csv2(Gompertz_simple$coefficients, "01_Gompertz_simple_resident.csv")
write.csv2(Gompertz_simple$mortQuant, "02_Gompertz_simple_resident.csv")

#=========================================================================================================================
# Run Gompertz model, bathtub function / 1981-2020 / time-dependent recapture
Gompertz_bathtub<-basta(object=datacheck$newData,studyStart=1981,studyEnd=2020, model="GO",shape="bathtub", 
                       recaptTrans=c(1981, 1982, 1983, 1984, 1985, 1986, 1987,
                                     1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998,
                                     1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009,
                                     2010, 2011,2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020), 
                       niter = 500000, burnin= 5000, thinn = 50, nsim = 4, 
                       parallel = TRUE, ncpus = 4, minAge = 1, lifeTable=TRUE)
# Generate some useful plots
plot(Gompertz_bathtub)
plot(Gompertz_bathtub, plot.trace = FALSE)
plot(Gompertz_bathtub, fancy = TRUE)

# Export outputs
write.csv2(Gompertz_bathtub$params, "00_Gompertz_bathtub_resident.csv")
write.csv2(Gompertz_bathtub$coefficients, "01_Gompertz_bathtub_resident.csv")
write.csv2(Gompertz_simple$mortQuant, "02_Gompertz_bathtub_resident.csv")

#=========================================================================================================================
# Run Gompertz model, Makeham function / 1981-2020 / time-dependent recapture
Gompertz_Makeham<-basta(object=datacheck$newData,studyStart=1981,studyEnd=2020, model="GO",shape="Makeham", 
                        recaptTrans=c(1981, 1982, 1983, 1984, 1985, 1986, 1987,
                                      1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998,
                                      1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009,
                                      2010, 2011,2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020), 
                        niter = 500000, burnin= 5000, thinn = 50, nsim = 4, 
                        parallel = TRUE, ncpus = 4, minAge = 1, lifeTable=TRUE)
# Generate some useful plots
plot(Gompertz_Makeham)
plot(Gompertz_Makeham, plot.trace = FALSE)
plot(Gompertz_Makeham, fancy = TRUE)

# Export outputs
write.csv2(Gompertz_Makeham$params, "00_Gompertz_Makeham_resident.csv")
write.csv2(Gompertz_Makeham$coefficients, "01_Gompertz_Makeham_resident.csv")
write.csv2(Gompertz_Makeham$mortQuant, "02_Gompertz_Makeham_resident.csv")

#=========================================================================================================================
# Run logistic model, simple function / 1981-2020 / time-dependent recapture
Logistic_simple<-basta(object=datacheck$newData,studyStart=1981,studyEnd=2020, model="LO",shape="simple", 
                        recaptTrans=c(1981, 1982, 1983, 1984, 1985, 1986, 1987,
                                      1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998,
                                      1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009,
                                      2010, 2011,2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020), 
                        niter = 500000, burnin= 5000, thinn = 50, nsim = 4, 
                        parallel = TRUE, ncpus = 4, minAge = 1, lifeTable=TRUE)
# Generate some useful plots
plot(Logistic_simple)
plot(Logistic_simple, plot.trace = FALSE)
plot(Logistic_simple, fancy = TRUE)

# Export outputs
write.csv2(Logistic_simple$params, "00_Logistic_simple_resident.csv")
write.csv2(Logistic_simple$coefficients, "01_Logistic_simple_resident.csv")
write.csv2(Logistic_simple$mortQuant, "02_Logistic_simple_resident.csv")

#=========================================================================================================================
# Run logistic model, bathtub function / 1981-2020 / time-dependent recapture
Logistic_bathtub<-basta(object=datacheck$newData,studyStart=1981,studyEnd=2020, model="LO",shape="bathtub", 
                       recaptTrans=c(1981, 1982, 1983, 1984, 1985, 1986, 1987,
                                     1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998,
                                     1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009,
                                     2010, 2011,2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020), 
                       niter = 500000, burnin= 5000, thinn = 50, nsim = 4, 
                       parallel = TRUE, ncpus = 4, minAge = 1, lifeTable=TRUE)
# Generate some useful plots
plot(Logistic_bathtub)
plot(Logistic_bathtub, plot.trace = FALSE)
plot(Logistic_bathtub, fancy = TRUE)

# Export outputs
write.csv2(Logistic_bathtub$params, "00_Logistic_bathtub_resident.csv")
write.csv2(Logistic_bathtub$coefficients, "01_Logistic_bathtub_resident.csv")
write.csv2(Logistic_bathtub$mortQuant, "02_Logistic_bathtub_resident.csv")

#=========================================================================================================================
# Run logistic model, Makeham function / 1981-2020 / time-dependent recapture
Logistic_Makeham<-basta(object=datacheck$newData,studyStart=1981,studyEnd=2020, model="LO",shape="Makeham", 
                        recaptTrans=c(1981, 1982, 1983, 1984, 1985, 1986, 1987,
                                      1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998,
                                      1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009,
                                      2010, 2011,2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020), 
                        niter = 500000, burnin= 5000, thinn = 50, nsim = 4, 
                        parallel = TRUE, ncpus = 4, minAge = 1, lifeTable=TRUE)
# Generate some useful plots
plot(Logistic_Makeham)
plot(Logistic_Makeham, plot.trace = FALSE)
plot(Logistic_Makeham, fancy = TRUE)

# Export outputs
write.csv2(Logistic_Makeham$params, "00_Logistic_Makeham_resident.csv")
write.csv2(Logistic_Makeham$coefficients, "01_Logistic_Makeham_resident.csv")
write.csv2(Logistic_Makeham$mortQuant, "02_Logistic_Makeham_resident.csv")

#=========================================================================================================================
# Run Weibull model, simple function / 1981-2020 / time-dependent recapture
Weibull_simple<-basta(object=datacheck$newData,studyStart=1981,studyEnd=2020, model="WE",shape="simple", 
                       recaptTrans=c(1981, 1982, 1983, 1984, 1985, 1986, 1987,
                                     1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998,
                                     1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009,
                                     2010, 2011,2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020), 
                       niter = 500000, burnin= 5000, thinn = 50, nsim = 4, 
                       parallel = TRUE, ncpus = 4, minAge = 1, lifeTable=TRUE)
# Generate some useful plots
plot(Weibull_simple)
plot(Weibull_simple, plot.trace = FALSE)
plot(Weibull_simple, fancy = TRUE)

# Export outputs
write.csv2(Weibull_simple$params, "00_Weibull_simple_resident.csv")
write.csv2(Weibull_simple$coefficients, "01_Weibull_simple_resident.csv")
write.csv2(Weibull_simple$mortQuant, "02_Weibull_simple_resident.csv")

#=========================================================================================================================
# Run logistic model, bathtub function / 1981-2020 / time-dependent recapture
Weibull_bathtub<-basta(object=datacheck$newData,studyStart=1981,studyEnd=2020, model="WE",shape="bathtub", 
                        recaptTrans=c(1981, 1982, 1983, 1984, 1985, 1986, 1987,
                                      1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998,
                                      1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009,
                                      2010, 2011,2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020), 
                        niter = 500000, burnin= 5000, thinn = 50, nsim = 4, 
                        parallel = TRUE, ncpus = 4, minAge = 1, lifeTable=TRUE)
# Generate some useful plots
plot(Weibull_bathtub)
plot(Weibull_bathtub, plot.trace = FALSE)
plot(Weibull_bathtub, fancy = TRUE)

# Export outputs
write.csv2(Weibull_bathtub$params, "00_Weibull_bathtub_resident.csv")
write.csv2(Weibull_bathtub$coefficients, "01_Weibull_bathtub_resident.csv")
write.csv2(Weibull_bathtub$mortQuant, "02_Weibull_bathtub_resident.csv")

#=========================================================================================================================
# Run logistic model, Makeham function / 1981-2020 / time-dependent recapture
Weibull_Makeham<-basta(object=datacheck$newData,studyStart=1981,studyEnd=2020, model="WE",shape="Makeham", 
                        recaptTrans=c(1981, 1982, 1983, 1984, 1985, 1986, 1987,
                                      1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998,
                                      1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009,
                                      2010, 2011,2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020), 
                        niter = 500000, burnin= 5000, thinn = 50, nsim = 4, 
                        parallel = TRUE, ncpus = 4, minAge = 1, lifeTable=TRUE)
# Generate some useful plots
plot(Weibull_Makeham)
plot(Weibull_Makeham, plot.trace = FALSE)
plot(Weibull_Makeham, fancy = TRUE)

# Export outputs
write.csv2(Weibull_Makeham$params, "00_Weibull_Makeham_resident.csv")
write.csv2(Weibull_Makeham$coefficients, "01_Weibull_Makeham_resident.csv")
write.csv2(Weibull_Makeham$mortQuant, "02_Weibull_Makeham_resident.csv")

#=========================================================================================================================
# Run exponential model
Exponential<-basta(object=datacheck$newData,studyStart=1981,studyEnd=2020, model="EX",shape="simple", 
                   recaptTrans=c(1981, 1982, 1983, 1984, 1985, 1986, 1987,
                                 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998,
                                 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009,
                                 2010, 2011,2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020),
                   niter = 500000, burnin= 5000, thinn = 50, nsim = 4, 
                   parallel = TRUE, ncpus = 4, minAge = 1, lifeTable=TRUE)
# Generate some useful plots
plot(Exponential)
plot(Exponential, plot.trace = FALSE)
plot(Exponential, fancy = TRUE)

# Export outputs
write.csv2(Exponential$params, "00_Exponential_resident.csv")
write.csv2(Exponential$coefficients, "01_Exponential_resident.csv")
write.csv2(Exponential$mortQuant, "02_Exponential_resident.csv")

#=========================================================================================================================
# Get DIC for all models

Gompertz_simple$DIC
Gompertz_bathtub$DIC
Gompertz_Makeham$DIC
Logistic_simple$DIC
Logistic_bathtub$DIC
Logistic_Makeham$DIC
Weibull_simple$DIC
Weibull_bathtub$DIC
Weibull_Makeham$DIC
Exponential$DIC

save.image("workspace_30112023.RData")



