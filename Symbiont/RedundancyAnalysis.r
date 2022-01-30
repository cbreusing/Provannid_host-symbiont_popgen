library(vegan)
library(ggplot2)
library(ape)
library(tidyverse)
library(poppr)
library(pegas)
library(adegenet)
library(ade4)
library(psych)

setwd("/Users/Corinna/Documents/Work/Beinart_Lab/Snail_metagenomics/A_boucheti/Epsilon")

data <- t(read.table("Epsilon.Freebayes.FINAL.GT.FORMAT", header=T, row.names=1))
newdat <- apply(data, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))

env <- read.table("Epsilon.RDA.txt", header=T)

identical(rownames(newdat), env[,1]) 

pdf("Variable_correlation.pdf")
pairs.panels(env[,3:7], scale=T)
dev.off()

pred <- subset(env, select=-c(Depth, Longitude, Year))

rda <- capscale(newdat ~ Fluid + Condition(Latitude), data = pred, na.action = na.exclude, add = "cailliez")

RsquareAdj(rda)
vif.cca(rda) #Should be below 10

signif.full <- anova.cca(rda, parallel=getOption("mc.cores"))
signif.full
signif.axis <- anova.cca(rda, by="axis", parallel=getOption("mc.cores"))
signif.axis
signif.margin <- anova.cca(rda, by="margin", parallel=getOption("mc.cores"))
signif.margin

env$Vent = factor(env$Vent, levels = c("Kilo Moana", "Tow Cam", "ABE"))
col <- c("Kilo Moana" = "plum4", "Tow Cam" = "plum3", "ABE" = "lightsteelblue")

pdf("Epsilon_RDA.pdf")
plot(rda, type="n", scaling=3)
points(rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)
points(rda, display="sites", pch=22, cex=2, col="gray32", scaling=3, bg=col[env$Vent])
text(rda, scaling=3, display="bp", col="black", cex=1)
legend("bottomleft", legend=levels(env$Vent), bty="n", col="gray32", pch=23, cex=1, pt.bg=col)
dev.off()

load.rda <- scores(rda, choices = c(1), display = "species")
rownames(load.rda) <- colnames(newdat)

hist(load.rda[,1], main="Loadings on RDA1", xlab="RDA1")

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)    
  x[x < lims[1] | x > lims[2]] 
}
  
cand1 <- outliers(load.rda[,1],2.5)

ncand <- length(cand1)

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
colnames(cand1) <- c("axis","snp","loading")

cand <- rbind(cand1)
cand$snp <- as.character(cand$snp)

foo <- matrix(nrow=(ncand), ncol=1)
colnames(foo) <- c("Fluid")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- newdat[,nam]
  foo[i,] <- apply(pred[3], 2, function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  

length(cand$snp[duplicated(cand$snp)])
cand <- cand[!duplicated(cand$snp),]

for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,5] <- names(which.max(abs(bar[4])))
  cand[i,6] <- max(abs(bar[4]))
}

colnames(cand)[5] <- "predictor"
colnames(cand)[6] <- "correlation"

table(cand$predictor)

write.table(cand, file="Candidate_selected_SNPs_conditioned.txt")

sel <- cand$snp
col.cand <- rownames(rda$CCA$v)

for (i in 1:length(sel)) {
  foo <- match(sel[i],col.cand)
  col.cand[foo] <- "black"
}

col.cand[grep("black",col.cand,invert=TRUE)] <- "darkgray"

pdf("Epsilon_RDA2.pdf")
plot(rda, type="n", scaling=3)
points(rda, display="species", pch=20, cex=0.7, col=col.cand, scaling=3)
points(rda, display="sites", pch=22, cex=2, col="gray32", scaling=3, bg=col[env$Vent])
text(rda, scaling=3, display="bp", col="black", cex=1)
legend("bottomleft", legend=levels(env$Vent), bty="n", col="gray32", pch=23, cex=1, pt.bg=col)
dev.off()

