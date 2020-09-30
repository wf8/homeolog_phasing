#
#
#
library(ggplot2)

min_pp = c(0.0, 0.7, 0.8, 0.9, 0.95)
d = read.csv('new_results.tsv', sep='\t')

# plot by seq_length
final_plot_data = data.frame()
for (pp in min_pp) {
    plot_data = d[d$nIndiv == 0.0001,]
    plot_data = plot_data[plot_data$AvgTetraPP > pp,]
    plot_data$seq_size_bin = cut(plot_data$seq_length, breaks=c(0,200,400,600,800,1000))
    pd = aggregate(as.numeric(plot_data$TetraCorrect)-1, by=list(Category=plot_data$seq_size_bin), FUN=mean)
    pd$min_pp = pp
    final_plot_data = rbind(pd, final_plot_data)
}
p = ggplot(final_plot_data) + 
    geom_point(aes(x=Category, y=x, color=min_pp)) + 
    geom_line(aes(x=Category, y=x, color=min_pp)) + 
    theme_classic()
ggsave('results/seq_length.pdf', p, width=5, height=4) 

# plot by ILS
unlinked_RF = 32

# plot by TetraDist
# plot by min(TetraAMinDip, TetraBMinDip)

# plot by min(HexaABDist, HexaACDist, HexaBCDist)
# plot by min(HexaAMinDip, HexaBMinDip, HexaCMinDip)


