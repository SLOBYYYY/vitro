vitro = read.csv('data.csv', header=T)
# Variable explanation:
# Check it here if in doubt: http://rsbweb.nih.gov/ij/docs/guide/146-30.html
# Slice: Index/Identifier of the slice
# Count: Total cell count (calculated from?)
# Total Area: Occupied pixels of cells on the canvas
# Average Size: Total Area / Count
# % Area: Occupied pixels / Total pixels in percent
# Mean: Mean value of pixels
# Mode: Mode of pixels
# Perimeter: Perimeter of all the cells
# Major: The primary axis of the best fitting ellipse
# Minor: The secondary axis of the best fitting ellipse 
# Angle: The angle of the best fitting ellipse (between X axis and Major axis)
# Circularity: How 'round' the cells are. Between 0 and 1 (4*Pi*Area/Perimeter^2).
#           Values may not be valid for very small particles
# Solidity: Area/Convex area. Between 0 and 1.
# Feret: Feret's diameter measures the longest distance between 2 points along a selection boundary.
#        Feret holds the length itself while other feret variables control the angle and the position.
#        More info --> http://en.wikipedia.org/wiki/Feret_diameter
# FeretX: The X of the starting coordinate of the Feret diameter
# FeretY: The Y of the starting coordinate of the Feret diameter
# FeretAngle: Houses the angle of Feret diameter (between 0 and 180)
# MinFeret: Shortest feret diameter
# Integrated density: The sum of the 'values' of the pixels.
# Median: Median value of pixels
# Skewness: Skewness of distribution. Negative means left skewness, positive means right skewness
# Kurtosis: How 'peaky' the distribution is compared to normal distribution.
#           Negative means flatter distribution, positive means higher
attach(vitro)

symnum(cor(vitro))

plot(Total.Area ~ Slice, type='l')

# Circularity stops the slow growing rate after a while
plot(Circ. ~ Slice, type='l', col='blue', main="Vitro data analysis", xlab="Slices", ylab="Area")
lm.mean = lm(Circ. ~ poly(Slice,3))
fitted.mean = predict(lm.mean)
lines(fitted.mean, col='red')

# Examine FeretAngle
hist(FeretAngle)
qqnorm(FeretAngle)
qqline(FeretAngle)
# Not your typical normal distribution, too 'thick' in the middle

# Feret and MinFeret are more or less moving together but the rate is 
# not really significant
plot(Slice, MinFeret, type="l", ylim=c(5, 30))
lines(Slice, Feret, type='l', col='red')
lines(Slice, Feret-MinFeret, type='l', col='blue')

plot(Slice, log(Average.Size), type="l", col="green")

par(mfrow=c(1,1))
plot(Slice, X.Area, type="l", col="red")
lm.areaPerc = lm(X.Area ~ poly(Slice,3))
fitted.areaPerc = predict(lm.areaPerc)
lines(fitted.areaPerc, col="blue")

par(mfrow=c(3,1))
plot(Major ~ Slice, type="l", col="red")
plot(Minor ~ Slice, type="l", col="green")
plot(Angle ~ Slice, type="l", col="blue")

par(mfrow=c(2,1))
plot(Circ. ~ Slice, type="l", col="red")
plot(Solidity ~ Slice, type="l", col="green")


# Time series with regular, moving average with smoothing 3 and 5
if (!'TTR' %in% rownames(installed.packages())) {
    install.packages('TTR')
}
library(TTR)
par(mfrow=c(2,1))
areaTs = ts(Total.Area)
plot.ts(areaTs)
areaSMA3 = SMA(Total.Area, n=3)
plot.ts(areaSMA3)
areaSMA6 = SMA(Total.Area, n=6)
plot.ts(areaSMA6)
areaSMA9 = SMA(Total.Area, n=9)
plot.ts(areaSMA9)

par(mfrow=c(1,1))
areaTs.model = HoltWinters(areaTs, beta=F, gamma=F)
plot(areaTs.model)

library(forecast)
areaTs.prediction = forecast.HoltWinters(areaTs.model, h=50)
plot.forecast(areaTs.prediction)

# Use a correlogram to see if there's a correlation between forecast errors
# (stored in residuals). We use a lag value from 1 to 20
acf(areaTs.prediction$residuals, lag.max=20)
# It shows no outstanding values

# Test whether there's significant evidence for non-zero correlations
# we use Ljung-Box test
Box.test(areaTs.prediction$residuals, lag=20, type='Ljung-Box')
# p is very high (0.8) so there's no evidence

# We should also check that the forecast errors are normally distributed
plotForecastErrors = function (forecastErrors) {
    binSize = IQR(forecastErrors) / 4
    sd = sd(forecastErrors)
    min = min(forecastErrors) - sd * 5
    max = max(forecastErrors) + sd * 3
    norm = rnorm(10000, mean = 0, sd=sd)
    min2 = min(forecastErrors)
    max2 = max(forecastErrors)
    if (min2 < min) { min = min2 }
    if (max2 > max) { max = max2 }
    bins = seq(min, max, binSize)
    hist(forecastErrors, col='red', freq=F, breaks=bins)
    hist = hist(norm, plot=F, breaks=bins)
    points(hist$mids, hist$density, type='l', col='blue', lwd=2)
}

# Looks more or less normal
plotForecastErrors(areaTs.prediction$residuals)

# Now try to introduce some slope estimation for forecasts
areaTs.model2 = HoltWinters(areaTs, gamma=F, l.start=86500, b.start=10)
plot(areaTs.model2)

areaTs.prediction2 = forecast.HoltWinters(areaTs.model2, h=50)
plot(areaTs.prediction2)

acf(areaTs.prediction2$residuals, lag.max=20)

Box.test(areaTs.prediction2$residuals, lag=20, type='Ljung-Box')
# Very high p value, so no autocorrelation

plotForecastErrors(areaTs.prediction2$residuals)
# Looks pretty normal