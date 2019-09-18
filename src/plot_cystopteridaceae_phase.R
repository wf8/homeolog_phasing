#
#       tip1_map, tip1_pp, tip2_map, tip2_pp, ....
# gene1,
# gene2,
# gene3,
#
#

library(ggplot2)
library(plyr)

genes = c('APP', 'GAP', 'IBR', 'PGI')
burnin = 0.1

results = data.frame()
for (i in 1:length(genes)) {

    f_in = paste0('output/cystopteridaceae_phase', i, '.log')

    # read in file and exclude burnin
    d = read.csv(f, sep='\t')
    d = d[floor(nrow(d)*burnin):nrow(d),]

    locus_results = data.frame()
    variables = names(d)
    for (j in 2:length(variables)) {
        
        variable = variables[j]
        c = count(d, vars=variable)
        prob = max(c['freq'])/sum(c['freq'])
        map = which(c['freq'] == max(c['freq']))
        map = c[variable][map,]

        locus_results[1,paste0(variable, '_map')] = map
        locus_results[1,paste0(variable, '_pp')] = prob

    }
    row.names(locus_results) = c(genes[i])
    results = rbind(results, locus_results)
}
results = t(results)
write.csv(results, 'output/cystopteridaceae_phase.csv')
