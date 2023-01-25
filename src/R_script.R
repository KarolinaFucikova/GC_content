#install.packages("ggplot2")
library("ggplot2")
library("phytools")

# setwd("C:/Users/k.fucikova/Desktop/WJT24_31")
data.name = "data/GC_content.csv"
gc_data <- read.csv(file=data.name,header=TRUE)
## changing column names as needed
colnames(gc_data)
## if some columns start with numbers, they should be renamed like so:
##names(gc_data)[names(gc_data) == "X18S_GC"] <- "SSU_GC"

#excluding Trochiscia and M. homosphaera which don't have any environmental data
gc_clean <- gc_data[-c(13,30),]


# simple scatter plot of gc as a function of average temp in the warmest month
# and other exploratory plots
plot(gc_data$max_temp,gc_data$SSU_GC,na.rm=TRUE)
plot(gc_data$cp_GC_total,gc_data$cp_GC_coding,na.rm=TRUE)
plot(gc_data$max_temp,gc_data$cp_GC_coding,na.rm=TRUE)
plot(gc_data$max_temp,gc_data$cp_GC_total,na.rm=TRUE)
plot(gc_data$cp_percent_coding,gc_data$cp_GC_coding,na.rm=TRUE)
plot(gc_data$cp_size,gc_data$cp_GC_coding,na.rm=TRUE)

## not a lot of obvious correlation, except for total GC to coding GC

## 18S GC as function of temperature, divided by habitat type
fancyplot<-ggplot(gc_clean, aes(x=max_temp, y=SSU_GC, color=habitat, shape=habitat)) + 
  geom_point(size=5) + geom_smooth(method=lm, aes(fill=habitat))
fancyplot 
# exploring some other options with color coding
fancyplot<-ggplot(gc_clean, aes(x=max_temp, y=SSU_GC, color=habitat, shape=habitat)) + 
  geom_point(size=5) + geom_smooth(method=lm, se=FALSE) +
  scale_color_manual(values=c('#56B4E9','#999999','#E69F00'))
fancyplot
# the above looks much nicer - note that we overwrote the earlier "fancyplot"

fancyplot<-ggplot(gc_clean, aes(x=max_temp, y=cp_GC_coding, color=habitat, shape=habitat)) + 
  geom_point(size=5) + geom_smooth(method=lm, aes(fill=habitat))
fancyplot 

## no correlation between max temp and cp GC coding; except maybe slight negative correlation in aquatic species


WJTtree <- "data/cp_nt_concatenated.con.tre"
WJT<-read.nexus(WJTtree)
#node number subtending OCC group happens to be 103 - this is not particularly easy to figure out
#I am rerooting the tree to have OCC as outgroup
# the current WJT can be viewed at anytime by simply running the line 
# plot(WJT)

#nodelabels(text=1:WJT$Nnode,node=1:WJT$Nnode+Ntip(WJT))
#getDescendants(WJT,103,curr=NULL)
#getSisters(WJT,106,mode=c("number","label"))
WJT<-reroot(WJT,103)
plot(WJT)

# mapping traits and inferring ancestral states
gc_data <- read.csv(file=data.name,header=TRUE, row.names=1)
# pick the column to map and store in a new object
gc_ssu_tomap <- subset(gc_data,select=SSU_GC)
## same as
## gc_ssu_tomap <- subset(gc_data,select=-c(habitat,max_temp,min_temp,precipitation,cp_size,cp_GC_total,cp_percent_coding,cp_GC_coding))
gc_cp_tomap <- subset(gc_data,select=cp_GC_total)
gc_cpcoding_tomap <- subset(gc_data,select=cp_GC_coding)

#turn the GC data into a vector
gc_ssu<-as.matrix(gc_ssu_tomap)[,1]
gc_cp<-as.matrix(gc_cp_tomap)[,1]
gc_cpcoding<-as.matrix(gc_cpcoding_tomap)[,1]

#estimate ancestral states, substitute objects from above
fit_ssu<-fastAnc(WJT,gc_ssu,vars=TRUE,CI=TRUE)
# anc.ML is the other method to infer ancestral states
# columns with missing data have to be handled separately
fit_cp<-fastAnc(WJT,gc_cp,vars=TRUE,CI=TRUE)
fit_cpcoding<-fastAnc(WJT,gc_cpcoding,vars=TRUE,CI=TRUE)
fit_cpcoding
# here we are mapping the GC content in chloroplast coding regions (fit_cpcoding) and overall chloroplast GC (fit_cp)

#pick variable to map, draw the tree
map <- contMap(WJT,gc_cpcoding,plot=FALSE)
plot(map,legend=0.5*max(nodeHeights(WJT)),fsize=c(0.7,0.7))
## some options to re-colorize
## what is the length of the current color ramp?
##n<-length(map$cols)
## change to grey scale
## map$cols[1:n]<-grey(0:(n-1)/(n-1))

## the default rainbow palette has high as blue and low as red
## this function will swap the colors to be more intuitive
setMap<-function(x,...){
  if(hasArg(invert)) invert<-list(...)$invert
  else invert<-FALSE
  n<-length(x$cols)
  if(invert) x$cols<-setNames(rev(x$cols),names(x$cols))
  else x$cols[1:n]<-colorRampPalette(...)(n)
  x
}

plot(setMap(map,invert=TRUE))
## more making pretty tree - ladderize
map$tree<-ladderize.simmap(map$tree)

plot(setMap(map,invert=TRUE),fsize=c(0.7,0.7))

#picking ribosomal GC to map, draw the tree
map <- contMap(WJT,gc_ssu,plot=FALSE)
plot(map,legend=0.5*max(nodeHeights(WJT)),fsize=c(0.7,0.7))

# recolorize using the invert function above
plot(setMap(map,invert=TRUE))
## more making pretty tree - ladderize
map$tree<-ladderize.simmap(map$tree)

plot(setMap(map,invert=TRUE),fsize=c(0.7,0.7))
# in this case, WJT stands out as extremely high in ribosomal GC content

# plotting in black and white
bw.contMap<-setMap(map,c("white","black"))
plot(bw.contMap, lwd=2, fsize=c(0.6,0.6))
##########################

## a couple of ways to correct for the effect of phylogeny in the temp-GC regression
## PIC is one
## gotta drop taxa with missing temp data from the tree
gc_data <- read.csv(file="data/GC_content.csv",header=TRUE, row.names=1)
WJTtree <- "data/cp_nt_concatenated.con.tre"
WJT<-read.nexus(WJTtree)
WJT<-reroot(WJT,103)

miss<-which(is.na(gc_data[,2]))
WJT_clean <- drop.tip(WJT,miss)
plot(WJT_clean)
# and also drop them from the table
data_clean <- gc_data[-miss,]
#max_temp<-setNames(data_clean[,"max_temp"],
#                        rownames(data_clean))
#SSU_GC<-setNames(data_clean[,"X18S_GC"],rownames(data_clean))
#pic.temp<-pic(max_temp,WJT_clean)
#pic.gc<-pic(SSU_GC,WJT_clean)
# the above are from tutorial but don't seem to work, cause warnings/errors
# below code seems to work ok

pic.temp<-pic(data_clean$max_temp,WJT_clean)
pic.gc<-pic(data_clean$SSU_GC,WJT_clean)

fit.pic<-lm(pic.gc~pic.temp -1)
# the -1 1 specifies that the regression is through the origin 
# (the intercept is set to zero) as recommended by Garland et al., 1992.
# https://www.r-phylo.org/wiki/HowTo/Phylogenetic_Independent_Contrasts
fit.pic


summary(fit.pic)
plot(pic.temp, pic.gc,
     xlab="Temperature PIC",
     ylab="18S GC Content PIC",bg="grey",
     cex=1.8,pch=21, cex.lab=1.3, cex.axis=1.2)
abline(fit.pic,lwd=2,lty="dashed",col="grey")


## labeling contrasts in the tree to see which point is which
x<-data_clean$max_temp
y<-data_clean$SSU_GC
names(x)<-row.names(data_clean)
names(y)<-row.names(data_clean)
plot(WJT_clean)
nodelabels(round(pic.temp,1), adj = c(0,0), frame="n")
nodelabels(round(pic.gc,1), adj = c(0,+0.5), frame="n")
# WJT/Spermatozopsis contrast is not the outliers; those are Chloromonas and Ankyra/Atractomorpha

### trying the same thing for cp GC coding
gc_data <- read.csv(file="data/GC_content.csv",header=TRUE, row.names=1)
WJTtree <- "data/cp_nt_concatenated.con.tre"
WJT<-read.nexus(WJTtree)
WJT<-reroot(WJT,103)

miss<-which(is.na(gc_data[,2]))
WJT_clean <- drop.tip(WJT,miss)
plot(WJT_clean)
# and also drop them from the table
data_clean <- gc_data[-miss,]
#max_temp<-setNames(data_clean[,"max_temp"],
#                        rownames(data_clean))
#SSU_GC<-setNames(data_clean[,"X18S_GC"],rownames(data_clean))
#pic.temp<-pic(max_temp,WJT_clean)
#pic.gc<-pic(SSU_GC,WJT_clean)
# the above are from tutorial but don't seem to work, cause warnings/errors
# below code seems to work ok

pic.temp<-pic(data_clean$max_temp,WJT_clean)
pic.gc<-pic(data_clean$cp_GC_coding,WJT_clean)

fit.pic<-lm(pic.gc~pic.temp+0)
fit.pic
## not sure what's going on with the labels; maybe it doesn't match the "Ankyra judayi" and Ankyra_judayi
## but checking tip labels and row names seems to match

summary(fit.pic)
plot(pic.temp,pic.gc,xlab="PICs for average temperature in warmest month",
     ylab="PICs for chloroplast coding GC content",bg="grey",
     cex=1.4,pch=21)
abline(fit.pic,lwd=2,lty="dashed",col="red")

