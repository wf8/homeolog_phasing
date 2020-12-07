library(coda)
library(ggplot2)

for (rep in 0:999) {
    line_out = paste0(rep)
    for (dir in c('output_w_dummy/', 'output_no_dummy/')) {
        line_out = paste0(line_out, ',')
        in_file = paste0(dir, rep, '/phasing.log')
        d = read.csv(in_file, sep='\t')
        ess = effectiveSize(d$Posterior[500:1999])[[1]]
        line_out = paste0(line_out, ess)
    }
    line_out = paste0(line_out, '\n')
    cat(line_out, file=paste0('ess.csv'), append=TRUE)
}

d = read.csv('ess.csv', header=FALSE, col.names=c('rep', 'ESS_with_dummy' ,'ESS_no_dummy'))
p = ggplot(d) + 
    geom_point(aes(x=ESS_no_dummy, y=ESS_with_dummy), alpha=0.5) + 
    geom_hline(yintercept=200, linetype='dashed') +
    geom_vline(xintercept=200, linetype='dashed') +
    annotate('text', x=400, y=750, label=paste0('proportion not converged'), size=3) +
    annotate('text', x=400, y=700, label=paste0('no dummy: ', sum(d$ESS_no_dummy < 200)/1000), size=3) +
    annotate('text', x=400, y=650, label=paste0('w/ dummy: ',sum(d$ESS_with_dummy < 200)/1000), size=3) +
    theme_classic()
ggsave('MCMC_convergence.pdf', p, width=5, height=4) 
