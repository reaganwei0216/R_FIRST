############################################################
# Collated by Prof. Ching-Shih Tsou (Ph.D.) at the Institute of Information and Decision Sciences/Data Science Applications Research Center, NTUB(國立臺北商業大學資訊與決策科學研究所/資料科學應用研究中心); CARS(中華R軟體學會), DSBA(台灣資料科學與商業應用協會)
############################################################

##### Just a little bit review about data structures in R
### Supplement 1: List subsetting
# Create a list first
g <- "My First List"
h <- c(25, 26, 18, 39)
j <- matrix(1:10, nrow=5, byrow=T)
j
?matrix
k <- c("one", "two", "three")
mylist <- list(title=g, ages=h, j, k)
mylist

# Different list subsetting
mylist[[2]] # get the element(s) of a list
class(mylist[[2]]) # a one-dimensional numeric vector
mylist[["ages"]] # another way to subsetting a list
mylist$ages
mylist[[1:2]] # Error in mylist[[1:2]] : subscript out of bounds

mylist[2] # get a part of a list, same as below
class(mylist[2]) # a list
mylist[1:2]

### Supplement 2: Accessing other data structures (copied from CPU.R)
# vector indexing
x <- 1:5
x
x[c(1,2,2,3,3,3,4,4,4,4)]
names(x) <- c("one", "two", "three", "four", "five")
x[4]
x[-4]
x[1:4]
x[-(1:4)]
x[c(1,4,2)]
x["four"]
x>3
x[x>3]
x[x > 3 & x < 5]
x[x %in% c(1, 3, 5)]

# matrix indexing
x <- matrix(1:12, nrow=3, ncol=4)
x
dimnames(x) <- list(paste("row", 1:3, sep=''),paste("col", 1:4, sep=''))
x[3,4]
x[3,]
x[,4]
x[,c(1,3)]
x["row3",]
x[,"col4"]

# logical indexing
x <- 1:10
x[c(T,T,T,F,F,F,F,F,F,F)] # T: keep it; F: drop out
x[c(rep(T,3), rep(F,7))]
x <= 3
x[x <= 3]

### Supplement 3: Finding help (should be included in Basics course)
help.start() # Manuals, Reference, and Miscellaneous Material

help('plot') # function name already known (or help(plot))
?plot # an alias of help

help.search('plot') # function name unknown, search string from Vignette Names, Code Demonstrations, and Help Pages (titles and keywords)
??plot # an alias of help.search

apropos('plot') # function name with 'plot'

find('plot') # It is a different user interface to the same task. simple.words: logical; if TRUE, the what argument is only searched as whole word.
find('plot', simple.words = FALSE) # look for packages with function 'plot'

### Supplement 4: Concept of masking (should be included in Basics course)
library(psych)
describe # psych::describe
describe(iris)
Hmisc::describe(iris)

### Supplement 5: Detaching and removing package (should be included in Basics course)
search()
?detach
detach(package:psych)

remove.packages("psych")
library(psych)

install.packages('psych')

### Supplement 6: Attaching and detaching dataset (should be included in Basics course)
search()
attach(iris)
search()
Sepal.Length
detach(iris)
Sepal.Length
search()

### Supplement 7: S3 - object-oriented programming
set.seed(168)
layout(matrix(c(1,1,2:5,6,6),4,2, byrow=TRUE))
weight <- seq(50, 70, length=10) + rnorm(10,5,1)
height <- seq(150, 170, length=10) + rnorm(10,6,2)
test <- data.frame(weight, height) # type [tab]
test.lm <- lm(weight ~ height, data=test)
class(test.lm)
class(AirPassengers)
plot(test)
plot(test.lm)
plot(AirPassengers)
layout(c(1))

### Supplement 8: Data Objects in {UsingR} and {datasets} (copied from CPU.R)
library(UsingR) # for datasets bumpers, aid, firstchi, and crime
firstchi
help(firstchi) # ? is an alias of help
class(firstchi)
names(firstchi)
bumpers
help(bumpers)
class(bumpers)
names(bumpers)
head(crime)
help(crime)
class(crime)
names(crime)
dimnames(crime)
Harman23.cor
help(Harman23.cor)
class(Harman23.cor)
names(Harman23.cor)
Harman23.cor$cov
class(Harman23.cor$cov)
dimnames(Harman23.cor$cov) # NOT names !
Titanic
help(Titanic)
class(Titanic) # "table" means "array"
dimnames(Titanic)

### Supplement 9: Vectorization
x <- 1:10
x
x^2

##### Supplement 10: apply, lapply, and sapply function for replacing the for loops
(m <- matrix(1:10, nrow=5, byrow=T))
apply(m, 1, mean)
apply(m, 2, mean)

temp <- list(x = c(1,3,5), y = c(2,4,6), z = c("a","b"))
temp
lapply(temp, length)
sapply(temp, length)

### Slides Codes start from here !!!
# install.packages('DMwR')
library(DMwR)
library(help=DMwR)

### Data Understanding
data(algae)
algae$mxPH
algae[,4]

# 數據結構
str(algae) # compactly display the structure of an arbitrary R object

# 看前幾筆數據
head(algae) # return the first part of an object
# 看最末幾筆數據
tail(algae) # return the last part of an object

# 檢視說明文件
help(algae) # the primary interface to the help systems (R documentation)

# 檢視資料集屬性
attributes(algae) # access an object's attributes

# 檢視變數名稱
names(algae) # get or set the names of an object

# 檢視資料物件類型
class(algae) # return object classes

# 變數個數
length(algae) # get or set the length of vectors (including lists) and factors

# 資料框行列維度
dim(algae) # retrieve or set the dimension of an object

# 檢視摘要統計值
summary(algae) # a generic function used to produce object summaries

### Data Cleaning
is.na(algae$mxPH)
which(is.na(algae$mxPH)) # which one is true?
mxPH.na.omit <- na.omit(object=algae$mxPH) # try mxPH.na.omit
length(mxPH.na.omit)
attributes(mxPH.na.omit)
attr(mxPH.na.omit, "na.action")
na.fail(algae$mxPH)
complete.cases(algae) # a logical vector of n.obs
which(!complete.cases(algae)) # who are the 16 incomplete cases ?

algae[!complete.cases(algae),]

algae[which(!complete.cases(algae)),] # take a look at the 16 incomplete cases
algae1 <- algae[complete.cases(algae),] # remove them from algae

apply(matrix(1:10, ncol=2), 1, mean)

apply(algae, 1, function(x) sum(is.na(x))) # count the number of missing variables for each observation

options(digits=2)
(Student <- c("John Davis", "Angela Williams", "Bullwinkle Moose", "David Jones", "Janice Markhammer", "Cheryl Cushing", "Reuven Ytzrhak", "Greg Knox", "Joel England", "Mary Rayburn"))
(Math <- c(502, 600, 412, 358, 495, 512, 410, 625, 573, 522))
Science <- c(95, 99, 80, 82, 75, 85, 80, 95, 89, 86)
English <- c(25, 22, 18, 15, 20, 28, 15, 30, 27, 18)
(roster <- data.frame(Student, Math, Science, English, stringsAsFactors=FALSE))
z <- scale(roster[,2:4])
apply(z, 2, mean)
apply(z, 2, sd)
score <- apply(z, 1, mean)
(roster <- cbind(roster, score))
(y <- quantile(score, c(.8,.6,.4,.2))) # it's a named numeric vector
roster$grade[score >= y[1]] <- "A" # append a new column to roster
roster$grade[score < y[1] & score >= y[2]] <- "B"
roster$grade[score < y[2] & score >= y[3]] <- "C"
roster$grade[score < y[3] & score >= y[4]] <- "D"
roster$grade[score < y[4]] <- "F"
name <- strsplit((roster$Student), " ")
(lastname <- sapply(name, "[", 2))
(firstname <- sapply(name, "[", 1))
roster <- cbind(firstname,lastname, roster[,-1])
roster <- roster[order(lastname, firstname),]
roster

sort(lastname, firstname) # An error occurs
sort(lastname) # what is the difference between sort and order?
sort(lastname, index.return=T)

x <- data.frame(a=c(1,2,4,5,6), x=c(9,12,14,21,8))
y <- data.frame(a=c(1,3,4,6), y=c(8,14,19,2))
merge(x,y) # natural join
merge(x,y,all=TRUE) # full outer join
merge(x,y,all.x=TRUE) # left outer join
merge(x,y,all.y=TRUE) # right outer join

head(iris, 6)
str(iris)
summary(iris)

library(Hmisc)
describe(iris)

# Be naive by subset
setosa <- subset(iris, Species=='setosa')
versicolor <- subset(iris, Species=='versicolor')
virginica <- subset(iris, Species=='virginica')
# Use a loop
results <- data.frame() # for gathering results
for (species in unique(iris$Species)) { # How smart it is !
  tmp <- subset(iris, Species==species) # get subset here
  count <- nrow(tmp)
  mean <- mean(tmp$Sepal.Length)
  median <- median(tmp$Sepal.Length)
  results <- rbind(results, data.frame(species, count, mean, median))
} #
results
# tapply {base}, input a vector that is grouped by a list of one or more factors

tapply(iris$Sepal.Length, iris$Species, FUN=length)
tapply(iris$Sepal.Length, iris$Species, FUN=mean)
tapply(iris$Sepal.Length, iris$Species, FUN=median) # Can multiple functions be set in FUN? No!!!
cbind(count=tapply(iris$Sepal.Length, iris$Species, FUN=length), mean=tapply(iris$Sepal.Length, iris$Species, FUN=mean), median=tapply(iris$Sepal.Length, iris$Species, FUN=median))

# aggregate {stats}, input a vector, matrix, or data frame that is grouped by a list of one or more factors
aggregate(Sepal.Length~Species, data=iris, FUN='length') # '~' is grouped by 
aggregate(Sepal.Length~Species, data=iris, FUN='mean')
aggregate(Sepal.Length~Species, data=iris, FUN='median') # Can multiple functions be set in FUN?  No!!!

# summaryBy {doBy}, performe multiple operations (functions)
library('doBy')
summaryBy(Sepal.Length~Species, data=iris, FUN=function(x) {c(count=length(x), mean=mean(x), median=median(x))}) # the columns are automatically given

# ddply {plyr}, multiple-part grouping variables and functions
library('plyr')
ddply(iris, 'Species', function(x) c(count=nrow(x), mean=colMeans(x[-5])))
# 

### Data Reshaping
pop <- read.csv('http://www.census.gov/popest/data/national/totals/2012/files/NST-EST2012-popchg2010_2012.csv')
pop <- pop[,c("Name","POPESTIMATE2010","POPESTIMATE2011","POPESTIMATE2012")]
colnames(pop)
colnames(pop) <- c('state', seq(2010, 2012))
head(pop, 2)
library('doBy')
top <- orderBy(~-2010, pop) # decreasing by 2010
head(top,3)
library('reshape2') # for melt() and dcast()
mtop <- melt(top, id.vars='state', variable.name='year', value.name='population') # wide to long
tail(mtop) # easier for analysis, plotting, database storage, etc.
mtop[order(mtop$state),]

dcast(mtop, state~year, value.var='population')[1:5,] # long to wide (sometimes you really can put the toothpaste back into the tube)

##### Exploring Data
custdata <- read.table("/Users/vince/cstsouMac/RandS/Rexamples/custdata.tsv", header=T, sep='\t')
summary(custdata)
summary(custdata[,c("is.employed", "housing.type", "recent.move", "num.vehicles")])
summary(custdata$income)
summary(custdata$age)
summary(custdata$income)
Income = custdata$income/1000
summary(Income)
library(ggplot2)
ggplot(custdata) +
  geom_histogram(aes(x=age),
                 binwidth=5, fill="gray")
diff(range(custdata$age))/30 # 4.88934
library(scales)
ggplot(custdata) + geom_density(aes(x=income)) +
  scale_x_continuous(labels=dollar) # right-skewed
ggplot(custdata) + geom_density(aes(x=income)) +
  scale_x_log10(breaks=c(100,1000,10000,100000), labels=dollar) +
  annotation_logticks(sides="bt") # bottom and top
ggplot(custdata) + geom_bar(aes(x=marital.stat), fill="gray")
ggplot(custdata) +
  geom_bar(aes(x=state.of.res), fill="gray") +
  coord_flip() +
  theme(axis.text.y=element_text(size=rel(0.8)))
statesums <- table(custdata$state.of.res)   # aggregates the data by state of residence -- exactly the information that barchart plots.
statef <- as.data.frame(statesums)   # Convert the table object to a data frame using as.data.frame(). The default column names are "Var1" and "Freq".
colnames(statef)<-c("state.of.res", "count") 	# Rename the columns for readability.
summary(statef)  	# Notice that the default ordering for the state.of.res variable is alphabetical.
statef <- transform(statef,
                    state.of.res=reorder(state.of.res, count)) # Use the reorder() function to set the state.of.res variable to be count-ordered. Use the transform() function.
#   to apply the transformation to the statef data frame.
summary(statef) # The state.of.res variable is now count ordered.
ggplot(statef)+ geom_bar(aes(x=state.of.res,y=count),
                         stat="identity",	# use stat="identity" to plot the data exactly as given.
                         fill="gray") +
  coord_flip() +                                       	# Flip the axes and reduce the size of the label text as before.
  theme(axis.text.y=element_text(size=rel(0.8)))
custdata2 <- subset(custdata, (custdata$age > 0 & custdata$age < 100 & custdata$income > 0))
cor(custdata2$age, custdata2$income)
ggplot(custdata2, aes(x=age, y=income)) + geom_point() + ylim(0, 200000)
ggplot(custdata2, aes(x=age, y=income)) + geom_point() + stat_smooth(method="lm") + ylim(0, 200000)
ggplot(custdata2, aes(x=age, y=income)) +
  geom_point() + geom_smooth() +
  ylim(0, 200000)
ggplot(custdata2, aes(x=age, y=as.numeric(health.ins))) +
  geom_point(position=position_jitter(w=0.05, h=0.05)) +
  geom_smooth()
library(hexbin)
ggplot(custdata2, aes(x=age, y=income)) +
  geom_hex(binwidth=c(5, 10000)) +
  geom_smooth(color="white", se=F) +
  ylim(0,200000)
ggplot(custdata) + geom_bar(aes(x=marital.stat,
                                fill=health.ins))
ggplot(custdata) + geom_bar(aes(x=marital.stat,
                                fill=health.ins),
                            position="dodge")
ggplot(custdata, aes(x=marital.stat)) +
  geom_bar(aes(fill=health.ins), position="fill") +
  geom_point(aes(y=-0.05), size=0.75, alpha=0.3,
             position=position_jitter(h=0.01))
ggplot(custdata2) +
  geom_bar(aes(x=housing.type, fill=marital.stat ),
           position="dodge") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(custdata2) +
  geom_bar(aes(x=marital.stat), position="dodge",
           fill="darkgray") +
  facet_wrap(~housing.type, scales="free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# title  : R - basic-iii-day2
# date   : 2014.8.3
# author : Ming-Chang Lee
# email  : alan9956@gmail.com
# RWEPA  : http://rwepa.blogspot.tw/

# Simple linear regression
?lm
# my.lm <- lm(formula, data="xxx")
# formula: y ~ x1 + x2 + ... +xn, actually a multiple regression
# end

# women: Average Heights and Weights for American Women
# y: weight
# x: height
fit.lm <- lm(weight ~ height, data=women) # '~' stands for weight is modeled by height
class(fit.lm)
summary(fit.lm)
# weight = -87.52+3.45*height

# verify residuals
names(fit.lm)
women$weight   # actual
fitted(fit.lm) # predicted
residuals(fit.lm) # residual=actual-predicted
women$weight - fitted(fit.lm)

# plot data
plot(women$height,women$weight, xlab="Height (in inches)", ylab="Weight (in pounds)", main="Average Heights and Weights for American Women")
abline(fit.lm, col="red")
# end

# Polynomial regression
fit.poly.lm <- lm(weight ~ height + I(height^2), data=women)
summary(fit.poly.lm)
# weight = 261.88 - 7.35*height + 0.083*height^2

# plot data with polynomial regression
plot(women$height,women$weight, main="Polynomial regression", xlab="Height (in inches)", ylab="Weight (in lbs)")
lines(women$height, fitted(fit.poly.lm), col="blue") # A polynomial curve! NOT a line!!!
# end

# cubic polynomial regression
fit.cubic.lm <- lm(weight ~ height + I(height^2) +I(height^3), data=women)
summary(fit.cubic.lm)

# plot data with cubic polynomial regression
plot(women$height,women$weight, main="Cubic polynomial regression", xlab="Height (in inches)", ylab="Weight (in lbs)")
lines(women$height, fitted(fit.cubic.lm), col="blue")
# end

# scatterplot{car}
library(car)
scatterplot(weight ~ height, data=women, pch=19, spread=FALSE,
            lty=2, # lty=2, dashed line in linear model
            main="Women Age 30-39",
            xlab="Height (inches)",
            ylab="Weight (lbs.)")
# end

# comparisons among lm, resistant, robust
library(MASS) # for lqs() & rlm()
library(nutshell) # for data shiller.index
data(shiller) # nutshell package
hpi.lm <- lm(Real.Home.Price.Index~Year, data=shiller.index)
hpi.rlm <- rlm(Real.Home.Price.Index~Year, data=shiller.index)
hpi.lqs <- lqs(Real.Home.Price.Index~Year, data=shiller.index)
plot(Real.Home.Price.Index~Year, pch=19, cex=0.3, data=shiller.index)
abline(reg=hpi.lm, lty=1)
abline(reg=hpi.rlm, lty=2, col="red")
abline(reg=hpi.lqs, lty=3, col="green")
legend(x=1900, y=200, legend=c("lm", "rlm", "lqs"), lty=c(1, 2, 3), col=c(1,2,3))

# nonlinear regression
# example - population growth
library(car)
data(USPop)
attach(USPop)
plot(year, population)

# nls demo
time <- 0:21
pop.mod <- nls(population ~ beta1/(1 + exp(beta2 + beta3*time)),
               start=list(beta1 = 350, beta2 = 4.5, beta3 = -0.3),
               trace=TRUE)
summary(pop.mod)

lines(year, fitted.values(pop.mod), lwd=2, col="red")

plot(year, residuals(pop.mod), type="b") # type="b" means that plotting points and line both
abline(h=0, lty=2) # low-level plooting

# Local Polynomial Regression
library(locfit)

# OldFaithful <- read.csv("OldFaithful.csv")
faithful[1:3,]

## density histograms and smoothed density histograms
## time of eruption

hist(faithful$eruptions,freq=FALSE)
fit1 <- locfit(~lp(eruptions),data=faithful) # lp is a local polynomial model term for Locfit models. Usually, it will be the only term on the RHS of the model formula.

par(mfrow=c(1,2))
hist(faithful$eruptions,freq=FALSE)
fit1 <- locfit(~lp(eruptions),data=faithful)
plot(fit1)
par(mfrow=c(1,1))

## waiting time to next eruption
hist(faithful$waiting,freq=FALSE)
fit2 <- locfit(~lp(waiting),data=faithful)
plot(fit2)


## general cross validation (GCV) of smoothing constant
## for waiting time to next eruption
alpha <- seq(0.20,1,by=0.01)
n1 <- length(alpha)
g <- matrix(nrow=n1,ncol=4)
for (k in 1:length(alpha)) {
  g[k,] <- gcv(~lp(waiting,nn=alpha[k]),data=faithful) # nn: Nearest neighbor component of the smoothing parameter
}
g # nn=0.66
alpha[47]

plot(g[,4]~g[,3],ylab="GCV",xlab="degrees of freedom")
abline(h=min(g[,4]), lty=4, col="red")
## minimum at nn = 0.66
fit2 <- locfit(~lp(waiting,nn=0.66,deg=2), data=faithful)
plot(fit2)

# polynominlal vs.regression
fitreg <- lm(waiting~eruptions,data=faithful)
fit3 <- locfit(waiting~lp(eruptions),data=faithful) # deg=2
plot(fit3, get.data=T)
abline(fitreg, col="red")
legend("topleft",legend=c("polynominlal","regression"),lty=1,col=c(1,2),inset=0.05)

# anova
##### One-Way Analysis 0f Variance
### Understanding the F Distribution
xaxis <- seq(0,6,.1)
xaxis
y1 <- df(xaxis, 3, 10)
y2 <- df(xaxis, 4, 15)
y3 <- df(xaxis, 10, 29)
plot(xaxis, y3, type="l", main="Comparing F Distributions")
points(xaxis, y2, type="l", col="red")
points(xaxis, y1, type="l", col="blue")

### Using the F Distribution to Test the Equality of Two Variances
load("/Users/vince/cstsouMac/RandS/Rexamples/Beginning R Data/stackeddata.rda")
str(stackeddata)
stackeddata$Class <- as.factor(stackeddata$Class)
attach(stackeddata)
var.test(PostTest~Class)
tapply(PostTest, Class, var)
1/.3045669 # Ratio of variance greater than one
TwoTailedP <- 2*(1 - pf(3.283351, 15, 20)) #  We have to reverse the numerator and denominator degrees of freedom
TwoTailedP # same as that reported by R's var.test function
detach(stackeddata)

### One-Way ANOVA
##An Example of the One-Way ANOVA
mpg=c(34,35,34.3,35.5,35.8,35.3,36.5,36.4,37,37.6,33.3,34,34.7,33,34.9) 
brand=c("A","A","A","A","A","B","B","B","B","B","C","C","C","C","C")
mileage=data.frame(mpg=mpg,brand=brand)
attach(mileage)
factor(brand)
mileage

boxplot(mpg~brand)

group <- factor(brand)
group

results <- aov(mpg~group)
results # k-1=3-1=2 and n-k=15-3=12
summary(results)
summary(aov(mpg ~ group)) # same as above !

model.tables(results, type = "means")
?model.tables

## Tukey HSD Test
TukeyHSD(results)

## Bonferroni-Corrected Post Hoc Comparisons
# require(stats)
# pairwise.t.test(mpg, group, p.adjust.method = "bonferroni")

### Using the anova Function
anova(lm(mpg ~ group))

##### Advanced Analysis of ANOVA
###Two-Way ANOVA
##An Example of a Two-Way ANOVA
load("/Users/vince/cstsouMac/RandS/Rexamples/Beginning R Data/twowayexample.rda")
twowayexample # higher is better
attach(twowayexample)
Format <- factor(Format)
Subject <- factor(Subject)
Format
Subject
tapply(Satisfaction, Format, mean)
tapply(Satisfaction, Subject, mean)
results <- aov(Satisfaction~ Format *Subject , data=twowayexample)
summary(results)

?TukeyHSD
TukeyHSD(results, "Format") # Online-Classroom
TukeyHSD(results, "Subject") # All are significant

##Examining Interactions
interaction.plot(Format, Subject, Satisfaction) # x-axis, trace variable, response variable
interaction.plot(Subject, Format, Satisfaction)
deatch(twowayexample)

### Repeated-Measures ANOVA
load("/Users/vince/cstsouMac/RandS/Rexamples/Beginning R Data/repeated.rda")
repeated # six research subjects and four weeks
attach(repeated)

id <- factor(id)
time <- factor(time)
results <- aov(fitness ~ time + Error(id/time))
summary(results)

result <- tapply(fitness, time, mean)
plot(result, type='o', xlab='Time', ylab='Fitness Level') # the average trend is a positive one

### Mixed-Factorial ANOVA
load("/Users/vince/cstsouMac/RandS/Rexamples/Beginning R Data/mixed.rda")
mixed # age is the between-groups factor and distraction (closed eyes, simple distraction, and complex distraction) is the within-subject factor
with(mixed, {
  id <- factor(id)
  age <- factor(age)
  distr <- factor(distr)
}) # it seems no need to do it !
attach(mixed)
interaction.plot(age, distr, score)

library(ez)
?ezANOVA
ezANOVA(data=mixed, dv=score, wid=id, within=distr, between=age)
# http://personality-project.org/r/r.anova.html

##### Supplements of Regression - 1
## Example: Predicting Medical Expenses ----
## Step 2: Exploring and preparing the data ----
insurance <- read.csv("~/cstsouMac/RandS/Rexamples/MLwithR/insurance.csv", stringsAsFactors = TRUE)
str(insurance) # BMI provides a sense of how over or under-weight a person is relative to their height. BMI is equal to weight (in kilograms) divided by height (in meters) squared. An ideal BMI is within the range of 18.5 to 24.9.


# summarize the charges variable
summary(insurance$charges) # right-skewed

# histogram of insurance charges
hist(insurance$charges)

# table of region
table(insurance$region) # nearly evenly among four geographic regions

# exploring relationships among features: correlation matrix
cor(insurance[c("age", "bmi", "children", "charges")]) # age and bmi appear to have a moderate correlation, meaning that as age increases, so does bmi. There is also a moderate correlation between age and charges, bmi and charges, and children and charges. 

# visualing relationships among features: scatterplot matrix
pairs(insurance[c("age", "bmi", "children", "charges")]) # pairs {graphics}

# more informative scatterplot matrix
library(psych)
pairs.panels(insurance[c("age", "bmi", "children", "charges")]) # correlation ellipse, the dot at the center of the ellipse indicates the point of the mean value for the x axis variable and y axis variable.

## Step 3: Training a model on the data ----
ins_model <- lm(charges ~ age + children + bmi + sex + smoker + region, data = insurance)
ins_model <- lm(charges ~ ., data = insurance) # this is equivalent to above

# see the estimated beta coefficients
ins_model

## Step 4: Evaluating model performance ----
# see more detail about the estimated beta coefficients
summary(ins_model) # the majority of predictions were between $2,850 over the true value and $1,400 under the true value, the R-squared value is 0.7494, we know that nearly 75 percent of the variation in the dependent variable is explained by our model

## Step 5: Improving model performance ----

# add a higher-order "age" term
insurance$age2 <- insurance$age^2

# BMI may have zero impact on medical expenditures for individuals in the normal weight range, but it may be strongly related to higher costs for the obese (that is, BMI of 30 or above).
# add an indicator for BMI >= 30
insurance$bmi30 <- ifelse(insurance$bmi >= 30, 1, 0)

# To interact the obesity indicator (bmi30) with the smoking indicator (smoker), we would write a formula in the form charges ~ bmi30*smoker.
# create final model
ins_model2 <- lm(charges ~ age + age2 + children + bmi + sex +
                   bmi30*smoker + region, data = insurance)

summary(ins_model2)
# The R-squared value has improved from 0.75 to about 0.87. The higher-order age2 term is statistically significant, as is the obesity indicator, bmi30. The interaction between obesity and smoking suggests a massive effect; in addition to the increased costs of over $13,404 for smoking alone, obese smokers spend another $19,810 per year. 


##### Supplements for Regression - 2
library(nutshell)
data(team.batting.00to08)
help(team.batting.00to08) # ? is an alias of help()
str(team.batting.00to08)

library(nutshell)
data(team.batting.00to08)
runs.mdl <- lm(formula=runs ~ singles+doubles+triples+homeruns+walks+hitbypitch+sacrificeflies+ stolenbases+caughtstealing, data=team.batting.00to08)
class(runs.mdl)
summary(runs.mdl) # a generic function, actully calls the summary.lm(...)
names(runs.mdl)
runs.mdl$coefficients # 10 coefficients
runs.mdl$residuals # 270 residuals
runs.mdl$effects # 4*68+2=274
runs.mdl$rank # 1 (intercept) + 9 (indep. var.) = 10 秩
runs.mdl$fitted.values # 270 fitted values
runs.mdl$assign # 0 1 2 3 4 5 6 7 8 9
runs.mdl$qr # 270 * 10
runs.mdl$df.residual # 260
runs.mdl$xlevels # a record of the levels of the factors used in fitting
runs.mdl$call
runs.mdl$terms # (1+9) * 9 with 0 and 1 entries
runs.mdl$model # 270 * 10, the model frame used (raw data)
op <- par(mfrow=c(2,2))
plot(runs.mdl) # call plot.lm
par(op)
# model fit (residuals) diagnostics
op <- par(mfrow=c(2,3))
plot(runs.mdl, which=1:6) # a generic function, actully calls the plot.lm(...)
par(op)

influence(runs.mdl) # hat value and dfbetas

confint(runs.mdl)

# selecting the best regression variables
reduced.runs.mdl <- step(runs.mdl, direction='backward') # the smaller the AIC or BIC, the better the fit

min.runs.mdl <- lm(runs ~ 1, data=team.batting.00to08) # the minimum model
fwd.runs.mdl <- step(min.runs.mdl, direction='forward', scope=( ~ singles+doubles+triples+homeruns+walks+hitbypitch+sacrificeflies+stolenbases+caughtstealing), trace=0)
summary(fwd.runs.mdl)

# ANOVA statistics for the fitted model
anova(runs.mdl)
anova(fwd.runs.mdl, runs.mdl) # Akaike's An Information Criterion (-2*log-likehood + k*npar, k=2), k=log(n) for BIC or SBC (Schwarz's Bayesian criterion), the smaller the AIC or BIC, the better the fit, but there is not statistically significant ! 

# predicting
dreamTeam <- data.frame(homeruns=270, singles=1000, doubles=400, walks=600, triples=25, sacrificeflies=55, hitbypitch=44) # creat the new data
dreamTeam
predict(fwd.runs.mdl, newdata=dreamTeam, interval='confidence') # a generic function, actully calls the predict.lm(...)
predict(fwd.runs.mdl, newdata=dreamTeam, interval='prediction') # wider than confidence interval

# refining the model based on previous calculations
runs.mdl2 <- update(runs.mdl, formula=runs ~ singles + doubles + triples + homeruns + walks + hitbypitch + stolenbases + caughtstealing + 0) # or -1
summary(runs.mdl2)

##### Supplements for ANOVA
# Comparing means across more than two groups
data(mort06.smpl)
str(mort06.smpl);help(mort06.smpl)
tapply(mort06.smpl$age, INDEX=list(mort06.smpl$Cause), FUN=summary)

length(which(mort06.smpl$Cause=='Suicide')) # verify the valid N from Deducer and NA's from summary(...)

m <- aov(age~Cause, data=mort06.smpl)
class(m)
summary(m)

# By lm(...)
mort06.smpl.lm <- lm(age~Cause, data=mort06.smpl)
anova(mort06.smpl.lm) # samev as aov

# By oneway.test(...)
?oneway.test
oneway.test(age~Cause, data=mort06.smpl)
oneway.test(age~Cause, data=mort06.smpl, var.equal=TRUE) # there is an argument about var.equal, same as aov and anova

# Multiple comparisons
TukeyHSD(m)
op <- par(mar=c(5.1,8.1,4.1,0.1))
plot(TukeyHSD(m), las=1, cex.lab=.5)
par(op)
# Creating an Interaction Plot
# install.packages('faraway')
library(faraway)
data(rats)
interaction.plot(rats$poison, rats$treat, rats$time)
# Note the error bars produced by Rcmdr is not shown here

data(births2006.smpl)
births2006.cln <- births2006.smpl[births2006.smpl$WTGAIN<99 & !is.na(births2006.smpl$WTGAIN),]
tapply(X=births2006.cln$WTGAIN, INDEX=births2006.cln$DOB_MM, FUN=mean) # the weight gain increases slightly during winter months, but is this difference statistically significant? Let's check it now.
aov(WTGAIN~DOB_MM, births2006.cln)
summary(aov(WTGAIN~DOB_MM, births2006.cln))

oneway.test(WTGAIN~DOB_MM, births2006.cln)

oneway.test(WTGAIN~DOB_MM, births2006.cln, var.equal=TRUE)

################# Linear Model Selection and Regularization
# Best Subset Selection
library(ISLR)
?Hitters
names(Hitters) # 打數、安打數、全壘打數、得分數、打點(RBI)、保送、年數、前六項職涯統計、何聯盟(1986年底)、東區或西區(division series, league champion series, world series)、觸殺、助殺、失誤、薪資、何聯盟(1987年初)
rownames(Hitters)
dim(Hitters) # 322 players*above 20 variables
sum(is.na(Hitters$Salary)) # 59
Hitters=na.omit(Hitters)
sum(complete.cases(ISLR::Hitters)) # 263 complete cases   雙冒號運算子「::」：存取指定命名空間中的物件。http://ppt.cc/Rzqg
dim(Hitters) # 263 players*above 20 variables
sum(is.na(Hitters)) # No more NAs
library(leaps) # A package for regression subset selection including exhaustive search selection
regfit.full=regsubsets(Salary~.,Hitters)
summary(regfit.full)
regfit.full=regsubsets(Salary~.,data=Hitters,nvmax=19) # Default nvmax is 8
reg.summary=summary(regfit.full) # reg.summary$obj + reg.summary$outmat
names(reg.summary)
reg.summary$obj
reg.summary$outmat
reg.summary$rsq # The rsq increases from 32% to almost 55%, when all variables are included.

op <- par(mfrow=c(2,2))
plot(reg.summary$rss,xlab="Number of Variables",ylab="RSS",type="l")
plot(reg.summary$adjr2,xlab="Number of Variables",ylab="Adjusted RSq",type="l")
which.max(reg.summary$adjr2)
points(11,reg.summary$adjr2[11], col="red",cex=2,pch=20)
plot(reg.summary$cp,xlab="Number of Variables",ylab="Cp",type='l')
which.min(reg.summary$cp)
points(10,reg.summary$cp[10],col="red",cex=2,pch=20)
plot(reg.summary$bic,xlab="Number of Variables",ylab="BIC",type='l')
which.min(reg.summary$bic)
points(6,reg.summary$bic[6],col="red",cex=2,pch=20)
par(op)

op <- par(mfrow=c(2,2))
plot(regfit.full,scale="r2")
plot(regfit.full,scale="adjr2")
plot(regfit.full,scale="Cp")
plot(regfit.full,scale="bic")
par(op)

?plot.regsubsets
class(regfit.full) # class "regsubsets"
coef(regfit.full,6) # extract the coefficient of the best six-variable model

# Forward and Backward Stepwise Selection
regfit.fwd=regsubsets(Salary~.,data=Hitters,nvmax=19,method="forward")
summary(regfit.fwd)

regfit.bwd=regsubsets(Salary~.,data=Hitters,nvmax=19,method="backward")
summary(regfit.bwd)

# The best seven-variable model identified by three selection methods are different.
coef(regfit.full,7)
coef(regfit.fwd,7)
coef(regfit.bwd,7)

### Ridge Regression and the Lasso

x=model.matrix(Salary~.,Hitters)[,-1]
y=Hitters$Salary

# Ridge Regression

library(glmnet) # Lasso and elastic-net regularized linear models
grid=10^seq(10,-2,length=100)
ridge.mod=glmnet(x,y,alpha=0,lambda=grid) # full range of scenarios for the tuning parameter lambda and standardize = TRUE
dim(coef(ridge.mod)) # 20 (19 predictors plus an intercept) * 100 (one for each value of lambda)
ridge.mod$lambda[50]
coef(ridge.mod)[,50] # much smalled estimates
# sqrt(sum(coef(ridge.mod)[-1,50]^2))
ridge.mod$lambda[60]
coef(ridge.mod)[,60]
sqrt(sum(coef(ridge.mod)[-1,60]^2))
predict(ridge.mod,s=50,type="coefficients")[1:20,]

set.seed(1)
train=sample(1:nrow(x), nrow(x)/2)
test=(-train)
y.test=y[test]
ridge.mod=glmnet(x[train,],y[train],alpha=0,lambda=grid, thresh=1e-12) # convergence threshold is defaulted to 1E-7

ridge.pred=predict(ridge.mod,s=4,newx=x[test,]) # arbitrarily choose lambda to 4
mean((ridge.pred-y.test)^2) # much lower test MSE

mean((mean(y[train])-y.test)^2) # MSE of the base case (only the intercept)
ridge.pred=predict(ridge.mod,s=1e10,newx=x[test,]) # lambda is very large, so all coefficients are forced to be zeros
mean((ridge.pred-y.test)^2) # same as the base case

ridge.pred=predict(ridge.mod,s=0,newx=x[test,],exact=T) # attention to s=0 and exact=TRUE for the exact least squares
mean((ridge.pred-y.test)^2)
# lm(y~x, subset=train)
predict(ridge.mod,s=0,exact=T,type="coefficients")[1:20,] # none of the coefficients are zero - ridge regression does not perform variable selection !

# The Lasso

lasso.mod=glmnet(x[train,],y[train],alpha=1,lambda=grid) #alpha=1 is the lasso penalty, and alpha=0 the ridge penalty
plot(lasso.mod)

# end

