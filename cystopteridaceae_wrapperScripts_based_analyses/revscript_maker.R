############################################## 
#   REV SCRIPT MAKER FOR HOMOLOGIZER 
#   by: Matt Johnson, Texas Tech University
#   December, 2020
##############################################

########## INPUTS ##########
# 1. A gene copy map that will set initial phase for copies by sample and subgenome
#       This should be a comma separated values file (CSV)
#       Gene columns should be in the same order as the gene list
############################

########## OUTPUTS #########
# 1. A file with homeolog phasing commands for RevBayes
# 2. A file with homeolog move commands for RevBayes
############################

options(error=traceback)
args = commandArgs(trailingOnly=TRUE)
if(length(args) < 1){
    print("Usage: Rscript --vanilla revscript_maker.R genecopymap.csv")
    quit()
}


#genelist = read.table(args[1])
#genelist = genelist[,1]

genecopymap = read.csv(args[1],header=T,stringsAsFactors=TRUE)

#Genes start in column 3
numLoci = length(genecopymap) - 2

missingtaxa.text = "data[%d].addMissingTaxa(\"%s\")\n"
missingtaxa.commands = c()
initial.phase.text = "data[%d].setHomeologPhase(\"%s\", \"%s\")\n"
initial.phase.commands = rep(NA,numLoci*nrow(genecopymap))
my.index = 1
for(gene in 3:length(genecopymap)){
    for(sample in 1:nrow(genecopymap)){
        if(grepl("BLANK",genecopymap[sample,gene])){
            missingtaxa.commands = append(missingtaxa.commands,
                                        sprintf(missingtaxa.text,
                                            gene-2,
                                            genecopymap[sample,gene]
                                            ))
        }
        
        initial.phase.commands[my.index] = sprintf(initial.phase.text,
                                                    gene-2,
                                                    genecopymap[sample,gene],
                                                    genecopymap[sample,2]
                                                    )
        my.index = my.index + 1
                
    }
}
cat(missingtaxa.commands,file="InitialPhase.rev")
cat(initial.phase.commands,file="InitialPhase.rev",append=TRUE)

########## Write mvHomeologPhase commands ##########

mv.text = "moves[++mvi] = mvHomeologPhase(ctmc[%d], \"%s\", \"%s\", weight=w)\n"
subgenomes.by.sample = split(genecopymap$Subgenome,genecopymap$Sample)

make.phase.comb = function(subgenome.samples,geneNum){
  #subgenome.combinations = combn(subgenome.samples,2)
  #apply(subgenome.combinations,2,function(x){
  #  sprintf(mv.text,
  #          geneNum,
  #          genecopymap[x[1],2],
  #          genecopymap[x[2],2]
  #  )
  #})
  out_text = ""
  for (i in 1:(length(subgenome.samples)-1)) {
    for (j in (i+1):length(subgenome.samples)) {
        out_text = paste0(out_text, sprintf(mv.text,
                                            geneNum,
                                            subgenome.samples[i],
                                            subgenome.samples[j]))
    }
  }
  return(out_text)
}

mv.commands = rep(NA,numLoci*length(subgenomes.by.sample))
my.index=1
for(gene in 1:numLoci){
  for(sample in 1:length(subgenomes.by.sample)){
    mv.commands[my.index] = make.phase.comb(subgenomes.by.sample[[sample]],gene)
    my.index = my.index + 1
  }
}

cat(mv.commands,file="PhaseMoves.rev")
