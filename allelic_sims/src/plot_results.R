library(ggplot2)

# For sim_num 0 - 499 there is *no* true allelic variation.
# For sim_num 500 - 999 there is allelic variation in loci 3 & 4.


d1 = read.csv('bf_output_no_dummy/marginals4.csv', header=FALSE, col.names=c('rep', 'marginal'))
d1 = d1[order(d1$rep),]
d2 = read.csv('bf_output_w_dummy/marginals4.csv', header=FALSE, col.names=c('rep', 'marginal'))
d2 = d2[order(d2$rep),]

bf = vector()
rep = vector()
haplotypes = vector()
x = 1
for (r in c(0:49,500:549)) {
    b = d2[d2$rep == r,][1,2] - d1[d1$rep == r,][1,2]
    bf = append(bf, b)
    rep = append(rep, x)
    x = x + 1
    if (r < 50) {
        haplotypes = append(haplotypes, '2')
    } else {
        haplotypes = append(haplotypes, '3')
    }
}
d = data.frame(bf=bf, rep=rep, haplotypes=haplotypes)
p = ggplot(d) + 
    geom_point(aes(x=rep, y=bf, color=haplotypes), alpha=0.7) + 
    geom_hline(yintercept=0, linetype='dashed') +
    geom_segment(x=50, xend=100, y=100, yend=100, linetype='dotted') +
    scale_color_manual(values=c('firebrick','grey34')) +
    ylab('Bayes factor') +
    xlab('Simulation replicate') +
    facet_wrap(vars(haplotypes), ncol=2, scales="free") +
    theme_classic() +
    theme(strip.background=element_blank(), strip.text=element_blank(), strip.placement = "outside") +
    guides(color=guide_legend(title='True number\nof mul-tree tips'))
ggsave('bayes_factors_sims2.pdf', p, width=7, height=3) 
