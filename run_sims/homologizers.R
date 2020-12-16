library(ggplot2)
library(gridExtra)
library(dplyr)
setwd("/scratch/joh97948/homologizer/homeolog_phasing/run_sims/")

simstats = read.csv("sim_summary.txt",header=TRUE,sep="\t")
simstats$nIndiv = as.factor(simstats$nIndiv)
summary(simstats)

noloESS = simstats[simstats$ESS > 200,]
noloPP = noloESS[noloESS$MinTetraPP > 0.9,]

p1 = ggplot(data=noloESS[noloESS$Type=="unlinked",],aes(x=nIndiv,y=TetraDist,col=TetraCorrect)) + geom_jitter(size=3,position=position_jitterdodge(0.4)) + ggtitle("Unlinked Gene Trees")
p2 = ggplot(data=noloESS[noloESS$Type=="shared",],aes(x=nIndiv,y=TetraDist,col=TetraCorrect)) + geom_jitter(size=3,position=position_jitterdodge(0.4)) + ggtitle("Shared Gene Trees")
g = arrangeGrob(p1,p2,nrow=2)
plot(g)
ggsave("TetraCorrectPhase.png",g)

noloESS %>% 
  group_by(Type,nIndiv) %>% count(TetraCorrect)
