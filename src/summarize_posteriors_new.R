#
# Outputs a file containing all non-zero joint phase assignments per locus for each sample.
# filename: joint_phase_probs_SAMPLE_NAME.csv
# 
# example file for a hexaploid sample:
# locus     hp_tip1 hp_tip2 hp_tip3 joint_prob
# locus1    seq1    seq3    seq2    0.40
# locus1    seq1    seq2    seq3    0.30
# locus1    seq2    seq1    seq3    0.20
# locus1    seq2    seq3    seq1    0.05
# locus1    seq3    seq2    seq1    0.03
# locus1    seq3    seq1    seq2    0.02
# locus2    seq1    seq3    seq2    0.40
# locus2    seq1    seq2    seq3    0.30
# locus2    seq2    seq1    seq3    0.20
# locus2    seq2    seq3    seq1    0.05
# locus2    seq3    seq2    seq1    0.03
# locus2    seq3    seq1    seq2    0.02
# locus3    seq1    seq3    seq2    0.40
# locus3    seq1    seq2    seq3    0.30
# locus3    seq2    seq1    seq3    0.20
# locus3    seq2    seq3    seq1    0.05
# locus3    seq3    seq2    seq1    0.03
# locus3    seq3    seq1    seq2    0.02
# 
#
# Also output a file containing the non-zero per tip marginal posterior probabilities of each 
# phase assignment:
# filename: marginal_phase_probs.csv
# 
# example file:
# tip_name  locus   phase       marginal_prob
# hp_tip1   locus1  seq1        0.72
# hp_tip1   locus1  seq3        0.20
# hp_tip1   locus1  seq2        0.08
# hp_tip1   locus2  seq1        0.88
# hp_tip1   locus2  seq2        0.10
# hp_tip1   locus2  seq3        0.02
# hp_tip1   locus3  seq1        0.79
# hp_tip1   locus3  seq3        0.11
# hp_tip1   locus3  seq2        0.10
# hp_tip2   locus1  seq3        0.65
#

args=commandArgs(trailingOnly=TRUE)
if(length(args) < 2){
    print("Usage: Rscript --vanilla prefix genecopymap.csv")
    quit()
    }
# directory that contains the MCMC log files
#input_dir = 'output'
prefix = args[1]
genecopymapFn = args[2]
genecopymap = read.csv(genecopymapFn,header=T)

numLoci = length(genecopymap) - 2
samples = split(genecopymap$Subgenome,genecopymap$Sample)
loci = names(genecopymap)[3:length(genecopymap)]


# polyploid samples and their tips
#samples = list('hexaploid_1'  = c('hexa_X_1', 'hexa_Y_1', 'hexa_Z_1'),
#               'tetraploid_1' = c('tetra_X_1','tetra_Y_1'))

# what percentage of MCMC samples to exclude?
burnin = 0.2

# for each sample loop over each locus
marginal_results = data.frame()
for (sample in names(samples)) {
    joint_results = data.frame()
    for (i in 1:numLoci) {

        # read in file and exclude burnin
        f_in = paste0(prefix, '_gene', i, '_phase.log')
        d = read.csv(f_in, sep='\t')
        d = d[floor(nrow(d)*burnin):nrow(d),]

        # get joint phase assignments for this locus
        d1 = d[, samples[[sample]]]
        joint_results_locus = as.data.frame(table(d1))
        joint_results_locus$joint_prob = joint_results_locus$Freq / sum(joint_results_locus$Freq)
        joint_results_locus$locus = loci[i]
        joint_results = rbind(joint_results, joint_results_locus) 

        # get marginal posterior probs
        for (tip in samples[[sample]]) {
            m = as.data.frame(table(d[tip]))
            m$marginal_prob = m$Freq / sum(m$Freq)
            m$phase = m$Var1
            m = within(m, rm(Freq))
            m = within(m, rm(Var1))
            m$locus = loci[i]
            m$tip_name = tip
            marginal_results = rbind(marginal_results, m)
        }
    }
    joint_results = within(joint_results, rm(Freq))
    out_file = paste0(prefix, '_joint_phase_probs_', sample, '.csv')
    write.csv(joint_results, out_file, row.names=FALSE)
}
out_file = paste0(prefix, '_marginal_phase_probs.csv')
write.csv(marginal_results, out_file, row.names=FALSE)
