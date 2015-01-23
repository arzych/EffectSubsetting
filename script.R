##Package Installation

if (!require(PerformanceAnalytics)){ 
  install.packages(PerformanceAnalytics) 
} 


if (!require(plyr)){ 
  install.packages(plyr) 
} 

if (!require(sfsmisc)){ 
  install.packages(sfsmisc) 
} 

if (!require(boot)){ 
  install.packages(boot) 
} 

library(PerformanceAnalytics)
library(plyr)
library(sfsmisc)
library(boot)


##Column Naming
colnames(DATAFILE)<-c( "studyname", "type", 
                       "R", "SEr", "Varr", "Z", 
                       "SEz", "Varz", "metaanalysis", 
                       "species", "n", "n2", "responsevar",
                       "venue", "treatment", "treatment_category"
) 


########################################################
####Effect Size vs. Sample Size, sort-of Funnel Plot####
########################################################

png("Pre-Subset/Pre-subset Effect Size vs. Sample Size.png", width=1200, height=1000)
par(mfrow=c(2,1))
correlatesample<-cbind(log10(DATAFILE$n),DATAFILE$R)
colnames(correlatesample)<-c("ln Sample", "Correlation")
plot(correlatesample, main="Correlation vs. Sample Size")
abline((mean(DATAFILE$R)), 0, col="red")
zsample<-cbind(log10(DATAFILE$n),DATAFILE$Z)
colnames(zsample)<-c("ln Sample", "Fishers_Z")
plot(zsample, main="Z vs. Sample Size")
abline((mean(DATAFILE$Z)), 0, col="red")


###################################################################
##~~~~Summary Data and Pre-treatment Counts Per Meta-Analysis~~~~##
###################################################################


sink("Pre-Subset/Pre-subset_summary_output_by_input_metaanalysis.txt")
Rmeta <- aggregate( R ~ metaanalysis,
                    DATAFILE,
                    summary )
print(Rmeta)
Rmetanum <- aggregate( R ~ metaanalysis,
                       DATAFILE,
                       length )
print(Rmetanum)
Zmeta <- aggregate( Z ~ metaanalysis,
                        DATAFILE,
                        summary )
Zmetanum <- aggregate( Z ~ metaanalysis,
                        DATAFILE,
                       length )
print(Zmetanum)
sink()

png("Pre-Subset/Pre-subsetRbymeta.png", width=2000, height=1400)
plot(Rmeta, main="Correlation by Source Meta-Analysis", cex.axis=0.6)
png("Pre-Subset/Pre-subsetZbymeta.png", width=2000, height=1400)
plot(Zmeta, main="Z by Source Meta-Analysis", cex.axis=0.6)

## CSV Output of counts
write.csv(cbind(Rmeta,Zmeta), "Pre-Subset/Pre-subset_datapermeta-analysis.csv")


############################
##~~~Outlier Extraction~~~##
############################

## Extract rows with Z >= 1.5
attach(DATAFILE)
outliers <- DATAFILE[abs(DATAFILE$Z) >= 1.5, ]
detach(DATAFILE)
write.csv(outliers, "Pre-Subset/Outliers/Pre-subsetExtractedOutliers.csv")


## Count outliers per meta-analysis
sink("Pre-Subset/Outliers/Pre-subset_Outlier_per_input_Meta-analysis_counts.txt")
count(outliers$metaanalysis)
sink()


#####################################################################
#####################################################################
#####################################################################
##########*****~~~~BEGINNING OF SUBSETTED ANALYSIS~~~~*****##########
#####################################################################
#####################################################################


#####################################################################
##  Calculates Median of subsets (species followed by sample size, ## 
##               author, then response variable)                   ##
#####################################################################


attach(DATAFILE)
Rdata <-aggregate(R, by=list(species, n, 
                             studyname, responsevar,
                             metaanalysis), 
                  FUN = median)                   

Zdata <-aggregate(Z, by=list(species, n, 
                             studyname, responsevar, 
                             metaanalysis),
                  FUN = median)
detach(DATAFILE)

##Stores Median values in output file
outputDATAFILE <- cbind(Rdata,Zdata)
write.csv(outputDATAFILE, "Post-Subset/SubsettedBymedSpeciesDatafile.csv")


## Extract vector of output effect size data from output file
medRdata <- Rdata[,6]
medZdata <- Zdata[,6]


#####################################################################
##          Meta-Analysis Counts Pre- and Post-Subsetting          ##
#####################################################################

sink("Post-Subset/Pre_and_post_subset_counts_by_meta_analysis.txt")
cat("pre-Median")
print(count(DATAFILE$metaanalysis))
cat("Total for Pre-Median")
print(nrow(DATAFILE))
cat("post-Median")
print(count(Rdata$Group.5))
cat("Total for Post-Median")
print(nrow(Rdata))
sink()

metacount<-cbind((count(DATAFILE$metaanalysis)),(count(Rdata$Group.5)))
row.names(metacount)<-metacount[,1]
metacount<-metacount[,c(-1,-3)]

png("Post-Subset/Meta-Analysis Counts Pre- and Post-Median.png", width=1400, height=1000)
u <- metacount$freq
i <- metacount$freq.1
height <- rbind(u, i)
mp <- barplot(height, beside = TRUE,
              ylim = c(0, 400.0), width=20,names.arg = row.names(metacount), 
              cex.names=0.6, cex.axis=0.75, 
              main="Meta-Analysis Counts Pre- and Post-Median")
text(mp, height, labels = format(height, 4),
     pos = 3, cex = 0.75)

write.csv(metacount, "Post-Subset/Meta-Analysis Counts Pre- and Post-Median.csv")



#######################################################
## ~~Make some plots, output is all exported to wd~~ ##
#######################################################

##Histograms and Density Plots
png("Post-Subset/SubsettedHistogram.png", width=1200, height=800)
par(mfrow=c(2,2))

hist(medRdata, prob=TRUE, breaks=100, main="Species Median Correlation Data", xlim=c(-1,1), ylim=c(0,3))
xfit<-seq(min(medRdata),max(medRdata),length=100) 
yfit<-dnorm(xfit,mean=mean(medRdata),sd=sd(medRdata)) 
lines(xfit, yfit, col="blue", lwd=2)


hist(medZdata, prob=TRUE, breaks=100, main="Species Median Species Median Z Data", xlim=c(-4,4), ylim=c(0,3))
xfit<-seq(min(medZdata),max(medZdata),length=100) 
yfit<-dnorm(xfit,mean=mean(medZdata),sd=sd(medZdata)) 
lines(xfit, yfit, col="blue", lwd=2)

R <- density(medRdata)
plot(R, main = "Density of Species Median Correlation Data", xlim=c(-1,1), ylim=c(0,2))
polygon(R, col="lightsteelblue1", border="white")
xfit<-seq(min(medRdata),max(medRdata),length=100) 
yfit<-dnorm(xfit,mean=mean(medRdata),sd=sd(medRdata)) 
lines(xfit, yfit, col="blue", lwd=2)

Z <- density(medZdata)
plot(Z, main = "Density of Species Median Z Data", xlim=c(-4,4), ylim=c(0,2))
polygon(Z, col="lightsteelblue1", border="white")
xfit<-seq(min(medZdata),max(medZdata),length=100) 
yfit<-dnorm(xfit,mean=mean(medZdata),sd=sd(medZdata)) 
lines(xfit, yfit, col="blue", lwd=2)

##Presentable Plot
png("Post-Subset/Zdistribution.png", width=1200, height=1000)
hist(medZdata, prob=FALSE, breaks=100, main="Median-Subsetted Z Data", xlab= "Fisher's Z", xlim=c(-2,3), ylim=c(0,50))
abline(v=median(medZdata), col="red")


## QQ Plot
png("Post-Subset/SubsettedQQandECDF.png", width=1200, height=800)
par(mfrow=c(2,2))
qqnorm(medRdata, main="R data Q-Q plot")
qqline(medRdata)

qqnorm(medZdata, main="Z data Q-Q plot")
qqline(medZdata)

## Empirical Cumulative Distribution Overlaid with Normal Curve 

chart.ECDF(medRdata, main="ECDF of Species Median Correlation Data")
chart.ECDF(medZdata, main="ECDF of Species Median Z Data")

## Boxplots
png("Post-Subset/SubsettedBoxplot.png", width=1200, height=800)
par(mfrow=c(1,2))
boxplot(medRdata, main="Boxplot of Species Median Correlation Data")
boxplot(medZdata, main="Boxplot of Species Median Z Data", ylim=c(-4,4))


#########################################
####~~SummaryData_and_NormalityTest~~####
#########################################

sink("Post-Subset/SubsettedSummaryData_and_NormalityTest.txt")

summary(medRdata)
shapiro.test(medRdata)

summary(medZdata)
shapiro.test(medZdata)

sink()



#############################
##~~~~~~ECDF Analysis~~~~~~##
#############################

x<-ecdf(medRdata)
y<-ecdf(medZdata)


sink("Post-Subset/SubsettedDataQuantilesMedian.csv")

##R data Quantiles
cat("R quantiles")
quantR<-quantile(x, probs=seq(0, 1, 0.01))
qR<-data.frame(quantR) 

print(qR)

#Binning

q <- data.frame(Categories = c("% negative","%  0 - 0.1","%  0.1 - 0.3",
                               "%  0.3 - 0.5","%  0.5 - 1"), 
                
                Values = c( (x(0)-x(-1)),
                            (x(0.1)-x(0)),
                            (x(0.3)-x(0.1)),
                            (x(0.5)-x(0.3)),
                            (x(1)-x(0.5))
                )
)

print(q)

##Z data Quantiles
cat("Z quantiles")
quantZ<-quantile(y, probs=seq(0, 1, 0.01))
qZ<-data.frame(quantR) 

print(qZ)

#Binning
w <- data.frame(Categories = c("% negative","%  0 - 0.1","%  0.1 - 0.310",
                               "%  0.310 - 0.549","%  0.549 - 1", "% > 1"), 
                
                Values = c( (y(0)-y(-2)),
                            (y(0.1)-y(0)),
                            (y(0.310)-y(0.1)),
                            (y(0.549)-y(0.310)),
                            (y(1)-y(0.549)),
                            (y(4)-y(1))
                )
)

print(w)

sink()

##Binned by small/medium/large
sink("Post-Subset/SubsettedBinnedProbability.txt")
cat("Rdata")
print(q)
cat("Zdata")
print(w)
sink()


##Barplot of distribution binned by small/medium/large etc.

png("Post-Subset/SubsettedQuantileplot.png", width=2000, height=1400)
par(mfrow=c(1,2))

barplot(q$Values, names.arg=q$Categories, border=NA, col="lightblue",
        ylim=c(0,0.5), main="Effect Size Distribution, correlation")

barplot(w$Values, names.arg=w$Categories, border=NA, col="lightblue",
        ylim=c(0,0.5), main="Effect Size Distribution, Fisher's Z")


#############################
##~~~~~~Conf Interval~~~~~~##
#############################

#plot with 95% bands
png("Post-Subset/ECDF with CI.png", width=2000, height=1400)
par(mfrow=c(1,2))
ecdf.ksCI(medRdata)
ecdf.ksCI(medZdata)

sink("Post-Subset/ConfidenceMeasures.txt")
##K-S statistic, basically the plus-or-minus value for distribution
KSd(length(medRdata))
KSd(length(medZdata))

##Bootstrapping for R
meanFunc <- function(x,i){mean(x[i])}
bootMean <- boot(medRdata ,meanFunc,10000)
#SE
sd(bootMean$t)
##bias
mean(bootMean$t) - bootMean$t0
##lower 95%
quantile(bootMean$t,.025)
##upper
quantile(bootMean$t, 0.975)

##BootStrapping for Z
meanFunc <- function(x,i){mean(x[i])}
bootMeanZ <- boot(medZdata ,meanFunc,10000)
#SE
sd(bootMeanZ$t)
##bias
mean(bootMeanZ$t) - bootMeanZ$t0
##lower 95%
quantile(bootMeanZ$t,.025)
##upper
quantile(bootMeanZ$t, 0.975)
sink()


