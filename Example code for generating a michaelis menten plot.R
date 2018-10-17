library(dplyr)
library(ggplot2)
library(tidyr)

############# Plotting Data ############# 

conc.uM <- c(0.5,1,2.5,3.5,5,7.5,10,15,25,50,70,75,100)
rate <- c(0.6,1.1,2.1,2.3,3.7,3.,4.3,4.8,5.3,6.0,5.1,5.7,5.8)

plot(conc.uM, rate, main="Plot Title", xlab="conc (uM)", ylab="rate (uM/min)")
lines(conc.uM, rate, lty="dotted", col="red")


############# Fitting the Michaelis-Menten equation ############# 

# Fit Michaelis-Menten equation
# input the data and plot
conc <- c(0.5,1,2.5,3.5,5,7.5,10,15,25,50,70,75,100)
rate <- c(0.6,1.1,2.1,2.3,3.7,3.,4.3,4.8,5.3,6.0,5.1,5.7,5.8)
test.df <- data.frame(conc,rate)
plot(test.df$conc,test.df$rate)

# perform the fitting
mm.nls <- nls(rate ~ (Vmax * conc / (Km + conc)), data=test.df, start=list(Km=5, Vmax=6))
summary(mm.nls)

# extract coefficients 
Km <- unname(coef(mm.nls)["Km"])
Vmax <- unname(coef(mm.nls)["Vmax"])

# plot data and plot line of best fit
x <- c(0:100)
y <- (Vmax*x/(Km+x))
lines(x,y, lty="dotted", col="blue")

# confidence intervals of parameters
confint(mm.nls)

# look at residuals and plot
mm.resid <- resid(mm.nls)
plot(test.df$conc, mm.resid)

# add weighting to fit
test.df$weight <- 1/test.df$conc^2
mm.weight.nls <- nls(rate ~ (Vmax * conc / (Km + conc)), data=test.df, start=list(Km=5, Vmax=6), weight=test.df$weight)
summary(mm.weight.nls)

############# Fitting the Michaelis-Menten equation with substrate inhibition ############# 

# Fit Michaelis-Menten equation with substrate inhibition
# input the data and plot
conc <- c(0.01,0.02,0.03,0.04,0.06,0.08,0.1,0.3,0.3,0.4,0.6,0.8,1,1,2,3,4,5,5,5,6)
rate <- c(3.1,5.2,5.6,5.9,7.2,7.6,8.4,9.2,10.4,10,11,10.9,10.3,10.4,10.1,9.6,9.4,9.6,8.6,8.7,8.5)
test.df <- data.frame(conc,rate)
plot(test.df$conc,test.df$rate)

# perform the fitting, look at plot to estimate start values
mminhib.nls <- nls(rate ~ (Vmax * conc / (Km + conc*(1+conc/Ks))), data=test.df, start=list(Km=0.2, Vmax=11, Ks=1))
summary(mminhib.nls)

# extract coefficients 
Km <- unname(coef(mminhib.nls)["Km"])
Vmax <- unname(coef(mminhib.nls)["Vmax"])
Ks <- unname(coef(mminhib.nls)["Ks"])

# plot data and plot line of best fit
x <- c(0:60)/10
y <- (Vmax*x/(Km+x*(1+x/Ks)))
lines(x,y, lty="dotted", col="blue")

# confidence intervals of parameters
confint(mminhib.nls)

# look at residuals and plot
mminhib.resid <- resid(mminhib.nls)
plot(test.df$conc, mminhib.resid)



############# Determining IC50 ############# 

# Calculate IC50
# enter the data into a data frame and plot
conc.uM <- c(300,150,75,38,19,5,2,1,0.6)
percent.activity <- c(2,7,12,22,36,53,67,83,85)
ic50.df <- data.frame(conc.uM, percent.activity)
ic50.df$conc.nM <- ic50.df$conc.uM * 1000
ic50.df$logconc.nM <- log10(ic50.df$conc.nM)
plot(ic50.df$logconc.nM,ic50.df$percent.activity)

# estimate initial values of the curve by examining the plot
# add a 4-parameter rodbard curve with these initial parameters to the plot 
# to check they are reasonable initial estimates
x <- c(2:12/2)
y <- 0+(80-0)/(1+(x/4)^10)
lines(x,y, col="red")

# fit the data using the nls() function to the 4-parameter logistic model
rodbard.fit <- nls(formula(percent.activity ~ bot+(top-bot)/(1+(logconc.nM/logic50)^slope)), algorithm="port", data=ic50.df, start=list(bot=0, top=80, logic50=4, slope=10), lower=c(bot=-Inf, top=-Inf, logic50=0, slope=-Inf) )

# generate a summary of the fit
summary(rodbard.fit)

# extract coefficients
top <- unname(coef(rodbard.fit)["top"])
bot <- unname(coef(rodbard.fit)["bot"])
logic50 <- unname(coef(rodbard.fit)["logic50"])
slope <- unname(coef(rodbard.fit)["slope"])

# calculate line of best fit and add to plot
y.fit <- bot+(top bot)/(1+(x/logic50)^slope)
lines(x,y.fit, col="green")

# log scale to linear scale to get IC50
# remember the results are in nM
10^logic50

############# Fitting competitive inhibition #############

# Fit competitive inhibition
# input the raw data into a data frame
conc <- c(3.9,1.9,0.7,0.5,0.4,0.3,3.9,1.9,0.7,0.5,0.4,0.3,3.9,1.9,0.7,0.5,0.4,0.3,3.9,1.9,0.7,0.5,0.4,3.9,1.9,0.7,0.5,0.4)
rate <- c(0.19,0.19,0.18,0.14,0.13,0.08,0.18,0.17,0.15,0.13,0.11,0.11,0.2,0.15,0.11,0.1,0.09,0.07,0.16,0.14,0.1,0.08,0.07,0.16,0.13,0.08,0.08,0.07)
inhibitor <- c(0,0,0,0,0,0,100,100,100,100,100,100,300,300,300,300,300,300,500,500,500,500,500,700,700,700,700,700)

kic.df <- data.frame(conc, rate, inhibitor)

# determine concentrations of inhibitor used in experiment
inhibitor.conc <- unique(kic.df$inhibitor)

# create colors for each inhibitor concentration
inhib.color <- rainbow(length(inhibitor.conc))

# generate a blank plot and then plot the raw data
plot(kic.df$conc,kic.df$rate, pch="")
for (i in 1:length(inhibitor.conc)) {
  points(subset(kic.df$conc, kic.df$inhibitor==inhibitor.conc[i]), subset(kic.df$rate, kic.df$inhibitor==inhibitor.conc[i]), col=inhib.color[i], pch=1)
}

# perform the fitting
kic.nls <- nls(rate ~ (Vmax * conc / (Km*(1+inhibitor/Kic) + conc)), data=kic.df, start=list(Km=0.5, Vmax=.2, Kic=300))

# generate a summary of the fit
summary(kic.nls)

# confidence intervals of parameters
confint(kic.nls)

# extract coefficients
Km <- unname(coef(kic.nls)["Km"])
Vmax <- unname(coef(kic.nls)["Vmax"])
Kic <- unname(coef(kic.nls)["Kic"])

# use the values directly in the equation to calculate line of best fit
fit.data <- expand.grid(x=(1:40)/10, inhib=inhibitor.conc)
fit.data$y <- Vmax*fit.data$x/(Km*(1+fit.data$inhib/Kic)+fit.data$x)

# plot lines of best fit
for (i in 1:length(inhibitor.conc)) {
  lines(subset(fit.data$x, fit.data$inhib==inhibitor.conc[i]), subset(fit.data$y, fit.data$inhib==inhibitor.conc[i]), col=inhib.color[i])
}

#######################################################################################
############# Script with Uncompetitive, Competitive and Mixed inhibition #############
#######################################################################################

# Fit multiple types of inhibition
# input the raw data into a data frame
conc <- c(3.9,1.9,0.7,0.5,0.4,0.3,3.9,1.9,0.7,0.5,0.4,0.3,3.9,1.9,0.7,0.5,0.4,0.3,3.9,1.9,0.7,0.5,0.4,3.9,1.9,0.7,0.5,0.4)
rate <- c(0.19,0.19,0.18,0.14,0.13,0.08,0.18,0.17,0.15,0.13,0.11,0.11,0.2,0.15,0.11,0.1,0.09,0.07,0.16,0.14,0.1,0.08,0.07,0.16,0.13,0.08,0.08,0.07)
inhibitor <- c(0,0,0,0,0,0,100,100,100,100,100,100,300,300,300,300,300,300,500,500,500,500,500,700,700,700,700,700)

ki.df <- data.frame(conc,rate,inhibitor)
ki.df$inv.conc <- 1/ki.df$conc
ki.df$inv.rate <- 1/ki.df$rate

# determine concentrations of inhibitor used in experiment
inhibitor.conc <- unique(ki.df$inhibitor)

# create colors for each inhibitor concentration
inhib.color <- rainbow(length(inhibitor.conc))

# perform the fittings
kic.nls <- nls(rate ~ (Vmax * conc / (Km*(1+inhibitor/Kic) + conc)), data=ki.df, start=list(Km=0.5, Vmax=.2, Kic=300))
mixed.nls <- nls(rate ~ (Vmax * conc / (Km*(1+inhibitor/Kic) + conc*(1+inhibitor/Kiu))), data=ki.df, start=list(Km=0.5, Vmax=.2, Kic=300, Kiu=100))
kiu.nls <- nls(rate ~ (Vmax * conc / (Km + conc*(1+inhibitor/Kiu))), data=ki.df, start=list(Km=0.5, Vmax=.2, Kiu=100))

# generate a summary of the fits
summary(kic.nls)
summary(mixed.nls)
summary(kiu.nls)

# extract coefficients - competitive inhibition
kic.Km <- unname(coef(kic.nls)["Km"])
kic.Vmax <- unname(coef(kic.nls)["Vmax"])
kic.Kic <- unname(coef(kic.nls)["Kic"])

# extract coefficients - mixed inhibition
mixed.Km <- unname(coef(mixed.nls)["Km"])
mixed.Vmax <- unname(coef(mixed.nls)["Vmax"])
mixed.Kic <- unname(coef(mixed.nls)["Kic"])
mixed.Kiu <- unname(coef(mixed.nls)["Kiu"])

# extract coefficients - uncompetitive inhibition
kiu.Km <- unname(coef(kiu.nls)["Km"])
kiu.Vmax <- unname(coef(kiu.nls)["Vmax"])
kiu.Kiu <- unname(coef(kiu.nls)["Kiu"])

# use the values directly in the equation to calculate line of best fit
fit.data <- expand.grid(x=(1:40)/10, inhib=inhibitor.conc)
fit.data$inv.x <- 1/fit.data$x
fit.data$kic.y <- kic.Vmax*fit.data$x/(kic.Km*(1+fit.data$inhib/kic.Kic)+fit.data$x)
fit.data$mixed.y <- mixed.Vmax*fit.data$x/(mixed.Km*(1+fit.data$inhib/mixed.Kic)+fit.data$x*(1+fit.data$inhib/mixed.Kiu))
fit.data$kiu.y <-kiu.Vmax*fit.data$x/(kiu.Km+fit.data$x*(1+fit.data$inhib/kiu.Kiu))
fit.data$inv.kic.y <- 1/fit.data$kic.y
fit.data$inv.mixed.y <- 1/fit.data$mixed.y
fit.data$inv.kiu.y <- 1/fit.data$kiu.y

############ Plot Data and Best Fit #############

# plot lines of best fit - competitive
# generate a blank plot and then plot the raw data
plot(ki.df$conc,ki.df$rate, pch="", main="Competitive")
for (i in 1:length(inhibitor.conc)) {
  points(subset(ki.df$conc, ki.df$inhibitor==inhibitor.conc[i]), subset(ki.df$rate, ki.df$inhibitor==inhibitor.conc[i]), col=inhib.color[i], pch=1)
}
for (i in 1:length(inhibitor.conc)) {
  lines(subset(fit.data$x, fit.data$inhib==inhibitor.conc[i]), subset(fit.data$kic.y, fit.data$inhib==inhibitor.conc[i]), col=inhib.color[i])
}

# plot lines of best fit - mixed
# generate a blank plot and then plot the raw data
plot(ki.df$conc,ki.df$rate, pch="", main="Mixed")
for (i in 1:length(inhibitor.conc)) {
  points(subset(ki.df$conc, ki.df$inhibitor==inhibitor.conc[i]), subset(ki.df$rate, ki.df$inhibitor==inhibitor.conc[i]), col=inhib.color[i], pch=1)
}
for (i in 1:length(inhibitor.conc)) {
  lines(subset(fit.data$x, fit.data$inhib==inhibitor.conc[i]), subset(fit.data$mixed.y, fit.data$inhib==inhibitor.conc[i]), col=inhib.color[i])
}

# plot lines of best fit - uncompetitive
# generate a blank plot and then plot the raw data
plot(ki.df$conc,ki.df$rate, pch="", main="Uncompetitive")
for (i in 1:length(inhibitor.conc)) {
  points(subset(ki.df$conc, ki.df$inhibitor==inhibitor.conc[i]), subset(ki.df$rate, ki.df$inhibitor==inhibitor.conc[i]), col=inhib.color[i], pch=1)
}
for (i in 1:length(inhibitor.conc)) {
  lines(subset(fit.data$x, fit.data$inhib==inhibitor.conc[i]), subset(fit.data$kiu.y, fit.data$inhib==inhibitor.conc[i]), col=inhib.color[i])
}

############ Plot Residuals #############

# look at residuals and plot
ki.df$kic.resid<- resid(kic.nls)
ki.df$mixed.resid<- resid(mixed.nls)
ki.df$kiu.resid<- resid(kiu.nls)

# plot Residuals - competitive
# generate a blank plot and then plot the raw data
plot(ki.df$conc,ki.df$kic.resid, pch="", main="Residuals - Competitive")
for (i in 1:length(inhibitor.conc)) {
  points(subset(ki.df$conc, ki.df$inhibitor==inhibitor.conc[i]), subset(ki.df$kic.resid, ki.df$inhibitor==inhibitor.conc[i]), col=inhib.color[i], pch=1)
}

# plot Residuals - mixed
# generate a blank plot and then plot the raw data
plot(ki.df$conc,ki.df$mixed.resid, pch="", main="Residuals - Mixed")
for (i in 1:length(inhibitor.conc)) {
  points(subset(ki.df$conc, ki.df$inhibitor==inhibitor.conc[i]), subset(ki.df$mixed.resid, ki.df$inhibitor==inhibitor.conc[i]), col=inhib.color[i], pch=1)
}

# plot Residuals - uncompetitive
# generate a blank plot and then plot the raw data
plot(ki.df$conc,ki.df$kiu.resid, pch="", main="Residuals - Uncompetitive")
for (i in 1:length(inhibitor.conc)) {
  points(subset(ki.df$conc, ki.df$inhibitor==inhibitor.conc[i]), subset(ki.df$kiu.resid, ki.df$inhibitor==inhibitor.conc[i]), col=inhib.color[i], pch=1)
}

############ Plot Lineweaver Burk #############

# plot lines of best fit - competitive
# generate a blank plot and then plot the raw data
plot(ki.df$inv.conc,ki.df$inv.rate, pch="", main="Lineweaver Burk - Competitive")
for (i in 1:length(inhibitor.conc)) {
  points(subset(ki.df$inv.conc, ki.df$inhibitor==inhibitor.conc[i]), subset(ki.df$inv.rate, ki.df$inhibitor==inhibitor.conc[i]), col=inhib.color[i], pch=1)
}
for (i in 1:length(inhibitor.conc)) {
  lines(subset(fit.data$inv.x, fit.data$inhib==inhibitor.conc[i]), subset(fit.data$inv.kic.y, fit.data$inhib==inhibitor.conc[i]), col=inhib.color[i])
}

# plot lines of best fit - mixed
# generate a blank plot and then plot the raw data
plot(ki.df$inv.conc,ki.df$inv.rate, pch="", main="Lineweaver Burk - Mixed")
for (i in 1:length(inhibitor.conc)) {
  points(subset(ki.df$inv.conc, ki.df$inhibitor==inhibitor.conc[i]), subset(ki.df$inv.rate, ki.df$inhibitor==inhibitor.conc[i]), col=inhib.color[i], pch=1)
}
for (i in 1:length(inhibitor.conc)) {
  lines(subset(fit.data$inv.x, fit.data$inhib==inhibitor.conc[i]), subset(fit.data$inv.mixed.y, fit.data$inhib==inhibitor.conc[i]), col=inhib.color[i])
}

# plot lines of best fit - uncompetitive
# generate a blank plot and then plot the raw data
plot(ki.df$inv.conc,ki.df$inv.rate, pch="", main="Lineweaver Burk - Uncompetitive")
for (i in 1:length(inhibitor.conc)) {
  points(subset(ki.df$inv.conc, ki.df$inhibitor==inhibitor.conc[i]), subset(ki.df$inv.rate, ki.df$inhibitor==inhibitor.conc[i]), col=inhib.color[i], pch=1)
}
for (i in 1:length(inhibitor.conc)) {
  lines(subset(fit.data$inv.x, fit.data$inhib==inhibitor.conc[i]), subset(fit.data$inv.kiu.y, fit.data$inhib==inhibitor.conc[i]), col=inhib.color[i])
}



