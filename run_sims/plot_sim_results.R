#
# Plot density ridge plots and proportion correct under all
# the simulated conditions.
# Will Freyman
#
library(ggplot2)
library(gridExtra)

d = read.csv('results/new_results.tsv', sep='\t')

# plot by seq_length
min_d = 25
plot_data = d[d$nIndiv == 0.0001,]
plot_data = plot_data[plot_data$TetraDist > min_d,]
plot_data = plot_data[plot_data$TetraAMinDip < min_d |
                      plot_data$TetraBMinDip < min_d,]
plot_data$correct = plot_data$TetraCorrect
plot_data$pp = plot_data$AvgTetraPP
plot_data1 = plot_data
plot_data = d[d$nIndiv == 0.0001,]
plot_data$correct = plot_data$HexaCorrect
plot_data = plot_data[plot_data$HexaABDist > min_d & 
                      plot_data$HexaACDist > min_d & 
                      plot_data$HexaBCDist > min_d,]
plot_data = plot_data[plot_data$HexaAMinDip < min_d |
                      plot_data$HexaBMinDip < min_d |
                      plot_data$HexaCMinDip < min_d,]
plot_data$pp = plot_data$AvgHexaPP
plot_data = rbind(plot_data, plot_data1)
plot_data$correct = as.logical(plot_data$correct)
plot_data$seq_size_bin = cut(plot_data$seq_length, breaks=seq(0, 1000, 200), include.lowest=TRUE,ordered_result=TRUE)

plot_data$bin = ''
plot_data$h = 0
plot_data$s = 0
plot_data$h2 = 0
h = 0
plot_data = plot_data[!(is.na(plot_data$seq_size_bin)),]
for (bin in levels(plot_data$seq_size_bin)) {
    if (nrow(plot_data[plot_data$seq_size_bin == bin & plot_data$correct == T,]) > 0) {
        plot_data[plot_data$seq_size_bin == bin & plot_data$correct == T,]$bin = paste0(as.character(bin), '_correct')
        plot_data[plot_data$seq_size_bin == bin & plot_data$correct == T,]$h = h
        plot_data[plot_data$seq_size_bin == bin & plot_data$correct == T,]$h2 = h + 0.25
        prop = nrow(plot_data[plot_data$seq_size_bin == bin & plot_data$correct == T,]) /
               nrow(plot_data[plot_data$seq_size_bin == bin,])
        plot_data[plot_data$seq_size_bin == bin & plot_data$correct == T,]$s = prop*1.5
    }
    if (nrow(plot_data[plot_data$seq_size_bin == bin & plot_data$correct == F,]) > 0) {
        plot_data[plot_data$seq_size_bin == bin & plot_data$correct == F,]$bin = paste0(as.character(bin), '_wrong')
        plot_data[plot_data$seq_size_bin == bin & plot_data$correct == F,]$h = h + 0.5
        plot_data[plot_data$seq_size_bin == bin & plot_data$correct == F,]$h2 = h + 0.25
        prop = nrow(plot_data[plot_data$seq_size_bin == bin & plot_data$correct == F,]) /
               nrow(plot_data[plot_data$seq_size_bin == bin,])
        plot_data[plot_data$seq_size_bin == bin & plot_data$correct == F,]$s = prop*1.5*2
    }
    h = h + 2
}
plot_data$bin = as.factor(plot_data$bin)
p1 = ggplot(plot_data) + 
    geom_density_ridges2(aes(x=pp, y=h, group=h, fill=correct, color=correct, scale=s), 
                        panel_scaling=TRUE, rel_min_height=0.005, alpha=0.6,
                        jittered_points = TRUE, #stat = "density",height = ..density..,
                        bandwidth = 0.03,
    position = position_points_jitter(width=0.03, height=0.1),
    point_size = 1, point_alpha = 0.4) + 
    ylab('Sequence size (bp)') +
    xlab('Posterior probability\n  ') +
    xlim(c(0.0,1.0)) +
    geom_segment(data=data.frame(y=c(0.25,2.25,4.25,6.25,8.25)),
                 aes(x=0.0, xend=1.0, y=y, yend=y), orientation='y', 
                 linetype='dotted', alpha=0.5, color='grey80') +
    scale_y_continuous(breaks=c(0.25,2.25,4.25,6.25,8.25), labels=c('1-200', '201-400', '401-600', '601-800', '801-1000'), limits=c(0,9.5)) +
    scale_fill_manual(values=c('firebrick','grey34')) +
    scale_color_manual(values=c('firebrick','grey34')) +
    theme_classic() +
    theme(legend.position = "none") 

plot_data$h2 = as.factor(plot_data$h2)
pd = aggregate(as.numeric(as.logical(plot_data$correct)), by=list(y=plot_data$h2), FUN=mean)
pd2 = aggregate(as.numeric(as.logical(plot_data$correct)), by=list(y=plot_data$h2), FUN=var)
pd$y = as.numeric(as.character(pd$y))
pd$z = as.numeric(as.character(pd2$x))
pd$zstart = pd$x - pd$z
pd$zend = pd$x + pd$z
if (nrow(pd[pd$zstart < 0.75,]) > 0)
    pd[pd$zstart < 0.75,]$zstart = 0.75
if (nrow(pd[pd$zend > 1,]) > 0)
    pd[pd$zend > 1,]$zend = 1
#pd[pd$y == 6.25,]$x = 0.98
p2 = ggplot(pd) + 
    geom_segment(aes(x=0.75, xend=1.0, y=y, yend=y), orientation='y', 
                 linetype='dotted', alpha=0.5, color='grey80') +
    geom_segment(aes(x=1.0, xend=1.0, y=0, yend=9), orientation='y', 
                 linetype='dotted', alpha=0.4, color='grey80') +
    geom_segment(aes(x=0.75, xend=0.75, y=0, yend=9), orientation='y', 
                 linetype='dotted', alpha=0.4, color='grey80') +
    geom_segment(aes(x=zstart, xend=zend, y=y, yend=y), 
                 orientation='y', alpha=0.5, size=1, color='darkorange2') +
    geom_point(aes(y=y, x=x), color='grey44') + 
    geom_line(aes(y=y, x=x), orientation='y', color='grey44') + 
    #geom_smooth(aes(y=y, x=x), method='loess',se=FALSE) + 
    xlab('Proportion\ncorrect') +
    ylim(c(0.0,9.5)) +
    scale_x_continuous(breaks=c(0.75,1.0), limits=c(0.75,1.01)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.line.y=element_blank(),
          axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 
#p = grid.arrange(p1, p2, nrow=1, widths=c(4,1))
#ggsave('results/seq_length.pdf', p, width=5, height=2.5) 



# plot by ILS
unlinked_RF = 32
min_d = 25
plot_data = d[d$nIndiv >= 0.001,]
plot_data = plot_data[plot_data$TetraDist > min_d,]
plot_data = plot_data[plot_data$TetraAMinDip < min_d | 
                      plot_data$TetraBMinDip < min_d,]
plot_data$correct = plot_data$TetraCorrect
plot_data$pp = plot_data$AvgTetraPP
plot_data1 = plot_data
plot_data = d[d$nIndiv >= 0.001,]
plot_data$correct = plot_data$HexaCorrect
plot_data = plot_data[plot_data$HexaABDist > min_d & 
                      plot_data$HexaACDist > min_d & 
                      plot_data$HexaBCDist > min_d,]
plot_data = plot_data[plot_data$HexaAMinDip < min_d |
                      plot_data$HexaBMinDip < min_d | 
                      plot_data$HexaCMinDip < min_d,]
plot_data$pp = plot_data$AvgHexaPP
plot_data = rbind(plot_data, plot_data1)
plot_data$correct = as.logical(plot_data$correct)
plot_data$ILS = plot_data$RF/unlinked_RF
plot_data$ILS_bin = cut(plot_data$ILS, breaks=c(0, 0.0001, 0.1, 0.2, 0.3, 1.0), include.lowest=TRUE, ordered_result=TRUE)

plot_data$bin = ''
plot_data$h = 0
plot_data$s = 0
plot_data$h2 = 0
h = 0
plot_data = plot_data[!(is.na(plot_data$ILS_bin)),]
for (bin in levels(plot_data$ILS_bin)) {
    if (nrow(plot_data[plot_data$ILS_bin == bin & plot_data$correct == T,]) > 0) {
        plot_data[plot_data$ILS_bin == bin & plot_data$correct == T,]$bin = paste0(as.character(bin), '_correct')
        plot_data[plot_data$ILS_bin == bin & plot_data$correct == T,]$h = h
        plot_data[plot_data$ILS_bin == bin & plot_data$correct == T,]$h2 = h + 0.25
        prop = nrow(plot_data[plot_data$ILS_bin == bin & plot_data$correct == T,]) /
               nrow(plot_data[plot_data$ILS_bin == bin,])
        plot_data[plot_data$ILS_bin == bin & plot_data$correct == T,]$s = prop*1.5
    }
    if (nrow(plot_data[plot_data$ILS_bin == bin & plot_data$correct == F,]) > 0) {
        plot_data[plot_data$ILS_bin == bin & plot_data$correct == F,]$bin = paste0(as.character(bin), '_wrong')
        plot_data[plot_data$ILS_bin == bin & plot_data$correct == F,]$h = h + 0.5
        plot_data[plot_data$ILS_bin == bin & plot_data$correct == F,]$h2 = h + 0.25
        prop = nrow(plot_data[plot_data$ILS_bin == bin & plot_data$correct == F,]) /
               nrow(plot_data[plot_data$ILS_bin == bin,])
        plot_data[plot_data$ILS_bin == bin & plot_data$correct == F,]$s = prop*1.5*2
    }
    h = h + 2
}
plot_data$bin = as.factor(plot_data$bin)
p3 = ggplot(plot_data) + 
    geom_density_ridges2(aes(x=pp, y=h, group=h, fill=correct, color=correct, scale=s), 
                        panel_scaling=TRUE, rel_min_height=0.005, alpha=0.6,
                        jittered_points = TRUE, #stat = "density",height = ..density..,
                        bandwidth = 0.03,
    position = position_points_jitter(width=0.03, height=0.1),
    point_size = 1, point_alpha = 0.4) + 
    ylab('ILS') +
    xlab('Posterior probability\n  ') +
    xlim(c(0.0,1.0)) +
    geom_segment(data=data.frame(y=c(0.25,2.25,4.25,6.25,8.25,10.25,12.25)),
                 aes(x=0.0, xend=1.0, y=y, yend=y), orientation='y', 
                 linetype='dotted', alpha=0.5, color='grey80') +
    scale_y_continuous(breaks=c(0.25,2.25,4.25,6.25,8.25,10.25,12.25), labels=c('None', '0.0-0.1', '0.1-0.2', '0.2-0.3', '0.3-0.5', '0.5-0.75', '0.75-1.0'), limits=c(0,9.5)) +
    scale_fill_manual(values=c('firebrick','grey34')) +
    scale_color_manual(values=c('firebrick','grey34')) +
    theme_classic() +
    theme(legend.position = "none") 

plot_data$h2 = as.factor(plot_data$h2)
pd = aggregate(as.numeric(as.logical(plot_data$correct)), by=list(y=plot_data$h2), FUN=mean)
pd2 = aggregate(as.numeric(as.logical(plot_data$correct)), by=list(y=plot_data$h2), FUN=var)
pd$y = as.numeric(as.character(pd$y))
pd$z = as.numeric(as.character(pd2$x))
pd$zstart = pd$x - pd$z
pd$zend = pd$x + pd$z
if (nrow(pd[pd$zstart < 0.75,]) > 0)
    pd[pd$zstart < 0.75,]$zstart = 0.75
if (nrow(pd[pd$zend > 1,]) > 0)
    pd[pd$zend > 1,]$zend = 1
p4 = ggplot(pd) + 
    geom_segment(aes(x=0.75, xend=1.0, y=y, yend=y), orientation='y', 
                 linetype='dotted', alpha=0.5, color='grey80') +
    geom_segment(aes(x=1.0, xend=1.0, y=0, yend=9), orientation='y', 
                 linetype='dotted', alpha=0.4, color='grey80') +
    geom_segment(aes(x=0.75, xend=0.75, y=0, yend=9), orientation='y', 
                 linetype='dotted', alpha=0.4, color='grey80') +
    geom_segment(aes(x=zstart, xend=zend, y=y, yend=y), 
                 orientation='y', alpha=0.5, size=1, color='darkorange2') +
    geom_point(aes(y=y, x=x), color='grey44') + 
    geom_line(aes(y=y, x=x), orientation='y', color='grey44') + 
    #geom_smooth(aes(y=y, x=x), method='loess',se=FALSE) + 
    xlab('Proportion\ncorrect') +
    ylim(c(0.0,9.5)) +
    scale_x_continuous(breaks=c(0.75,1.0), limits=c(0.75,1.01)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.line.y=element_blank(),
          axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 
#p = grid.arrange(p1, p2, nrow=1, widths=c(4,1))
#ggsave('results/ILS.pdf', p, width=5, height=2.5) 





# plot by Distance between subgenomes
min_d = 25
plot_data = d[d$nIndiv == 0.001,]
plot_data = plot_data[plot_data$TetraAMinDip < min_d | 
                      plot_data$TetraBMinDip < min_d,]
plot_data$correct = plot_data$TetraCorrect
plot_data$pp = plot_data$AvgTetraPP
plot_data$distance = plot_data$TetraDist
plot_data1 = plot_data
plot_data = d[d$nIndiv == 0.001,]
plot_data = plot_data[plot_data$HexaAMinDip < min_d |
                      plot_data$HexaBMinDip < min_d | 
                      plot_data$HexaCMinDip < min_d,]
plot_data = transform(plot_data, distance=pmin(HexaABDist, HexaACDist, HexaBCDist))
plot_data$pp = plot_data$AvgHexaPP
plot_data$correct = plot_data$HexaCorrect

plot_data = rbind(plot_data, plot_data1)
plot_data$correct = as.logical(plot_data$correct)

plot_data$distance_bin = cut(plot_data$distance, breaks=c(0, 5, 25, 50, 100), include.lowest=TRUE, ordered_result=TRUE)

plot_data$bin = ''
plot_data$h = 0
plot_data$s = 0
plot_data$h2 = 0
h = 0
plot_data = plot_data[!(is.na(plot_data$distance_bin)),]
for (bin in levels(plot_data$distance_bin)) {
    if (nrow(plot_data[plot_data$distance_bin == bin & plot_data$correct == T,]) > 0) {
        plot_data[plot_data$distance_bin == bin & plot_data$correct == T,]$bin = paste0(as.character(bin), '_correct')
        plot_data[plot_data$distance_bin == bin & plot_data$correct == T,]$h = h
        plot_data[plot_data$distance_bin == bin & plot_data$correct == T,]$h2 = h + 0.25
        prop = nrow(plot_data[plot_data$distance_bin == bin & plot_data$correct == T,]) /
               nrow(plot_data[plot_data$distance_bin == bin,])
        plot_data[plot_data$distance_bin == bin & plot_data$correct == T,]$s = prop*1.5
    }
    if (nrow(plot_data[plot_data$distance_bin == bin & plot_data$correct == F,]) > 0) {
        plot_data[plot_data$distance_bin == bin & plot_data$correct == F,]$bin = paste0(as.character(bin), '_wrong')
        plot_data[plot_data$distance_bin == bin & plot_data$correct == F,]$h = h + 0.5
        plot_data[plot_data$distance_bin == bin & plot_data$correct == F,]$h2 = h + 0.25
        prop = nrow(plot_data[plot_data$distance_bin == bin & plot_data$correct == F,]) /
               nrow(plot_data[plot_data$distance_bin == bin,])
        plot_data[plot_data$distance_bin == bin & plot_data$correct == F,]$s = prop*1.5
    }
    h = h + 2
}
plot_data$bin = as.factor(plot_data$bin)
p5 = ggplot(plot_data) + 
    geom_density_ridges2(aes(x=pp, y=h, group=h, fill=correct, color=correct, scale=s), 
                        panel_scaling=TRUE, rel_min_height=0.005, alpha=0.6,
                        jittered_points = TRUE, #stat = "density",height = ..density..,
                        bandwidth = 0.03,
    position = position_points_jitter(width=0.03, height=0.1),
    point_size = 1, point_alpha = 0.4) + 
    ylab('Distance between\nsubgenomes') +
    xlab('Posterior probability\n  ') +
    xlim(c(0.0,1.0)) +
    geom_segment(data=data.frame(y=c(0.25,2.25,4.25,6.25)),
                 aes(x=0.0, xend=1.0, y=y, yend=y), orientation='y', 
                 linetype='dotted', alpha=0.5, color='grey80') +
    scale_y_continuous(breaks=c(0.25,2.25,4.25,6.25), labels=c('0.0-0.05', '0.05-0.25', '0.25-0.50', '0.5-1.0'), limits=c(0,7.5)) +
    scale_fill_manual(values=c('firebrick','grey34')) +
    scale_color_manual(values=c('firebrick','grey34')) +
    theme_classic() +
    theme(legend.position = "none") 

plot_data$h2 = as.factor(plot_data$h2)
pd = aggregate(as.numeric(as.logical(plot_data$correct)), by=list(y=plot_data$h2), FUN=mean)
pd2 = aggregate(as.numeric(as.logical(plot_data$correct)), by=list(y=plot_data$h2), FUN=var)
pd$y = as.numeric(as.character(pd$y))
pd$z = as.numeric(as.character(pd2$x))
pd$zstart = pd$x - pd$z
pd$zend = pd$x + pd$z
if (nrow(pd[pd$zstart < 0.0,]) > 0)
    pd[pd$zstart < 0.0,]$zstart = 0.0
if (nrow(pd[pd$zend > 1,]) > 0)
    pd[pd$zend > 1,]$zend = 1
p6 = ggplot(pd) + 
    geom_segment(aes(x=0.0, xend=1.0, y=y, yend=y), orientation='y', 
                 linetype='dotted', alpha=0.5, color='grey80') +
    geom_segment(aes(x=1.0, xend=1.0, y=0, yend=7.5), orientation='y', 
                 linetype='dotted', alpha=0.4, color='grey80') +
    geom_segment(aes(x=0.0, xend=0.0, y=0, yend=7.5), orientation='y', 
                 linetype='dotted', alpha=0.4, color='grey80') +
    geom_segment(aes(x=zstart, xend=zend, y=y, yend=y), 
                 orientation='y', alpha=0.5, size=1, color='darkorange2') +
    geom_point(aes(y=y, x=x), color='grey44') + 
    geom_line(aes(y=y, x=x), orientation='y', color='grey44') + 
    #geom_smooth(aes(y=y, x=x), method='loess',se=FALSE) + 
    xlab('Proportion\ncorrect') +
    ylim(c(0.0,7.5)) +
    scale_x_continuous(breaks=c(0.0,1.0), limits=c(-0.01,1.01), labels=c('0.0','1.0')) +
    theme_classic() +
    theme(legend.position = "none",
          axis.line.y=element_blank(),
          axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 
#p = grid.arrange(p1, p2, nrow=1, widths=c(4,1))
#ggsave('results/genome_distance.pdf', p, width=5, height=2.5) 



# plot by Distance to nearest diploid
min_d = 25
plot_data = d[d$nIndiv == 0.001,]
plot_data = plot_data[plot_data$TetraDist > min_d,]
plot_data$correct = plot_data$TetraCorrect
plot_data$pp = plot_data$AvgTetraPP
plot_data$distance = pmin(plot_data$TetraAMinDip, plot_data$TetraBMinDip)
plot_data1 = plot_data
plot_data = d[d$nIndiv == 0.001,]
plot_data = plot_data[plot_data$HexaABDist > min_d & 
                      plot_data$HexaACDist > min_d & 
                      plot_data$HexaBCDist > min_d,]
plot_data$distance = pmin(plot_data$HexaAMinDip, plot_data$HexaBMinDip, plot_data$HexaCMinDip)
plot_data$pp = plot_data$AvgHexaPP
plot_data$correct = plot_data$HexaCorrect

plot_data = rbind(plot_data, plot_data1)
plot_data$correct = as.logical(plot_data$correct)

#plot_data$distance_bin = cut(plot_data$distance, breaks=c(0, 10, 25, 50, 75, 100), include.lowest=TRUE, ordered_result=TRUE)
plot_data$distance_bin = cut(plot_data$distance, breaks=c(0, 5, 25, 50, 100), include.lowest=TRUE, ordered_result=TRUE)

plot_data$bin = ''
plot_data$h = 0
plot_data$s = 0
plot_data$h2 = 0
h = 0
plot_data = plot_data[!(is.na(plot_data$distance_bin)),]
for (bin in levels(plot_data$distance_bin)) {
    if (nrow(plot_data[plot_data$distance_bin == bin & plot_data$correct == T,]) > 0) {
        plot_data[plot_data$distance_bin == bin & plot_data$correct == T,]$bin = paste0(as.character(bin), '_correct')
        plot_data[plot_data$distance_bin == bin & plot_data$correct == T,]$h = h
        plot_data[plot_data$distance_bin == bin & plot_data$correct == T,]$h2 = h + 0.25
        prop = nrow(plot_data[plot_data$distance_bin == bin & plot_data$correct == T,]) /
               nrow(plot_data[plot_data$distance_bin == bin,])
        plot_data[plot_data$distance_bin == bin & plot_data$correct == T,]$s = prop*1.5
    }
    if (nrow(plot_data[plot_data$distance_bin == bin & plot_data$correct == F,]) > 0) {
        plot_data[plot_data$distance_bin == bin & plot_data$correct == F,]$bin = paste0(as.character(bin), '_wrong')
        plot_data[plot_data$distance_bin == bin & plot_data$correct == F,]$h = h + 0.5
        plot_data[plot_data$distance_bin == bin & plot_data$correct == F,]$h2 = h + 0.25
        prop = nrow(plot_data[plot_data$distance_bin == bin & plot_data$correct == F,]) /
               nrow(plot_data[plot_data$distance_bin == bin,])
        plot_data[plot_data$distance_bin == bin & plot_data$correct == F,]$s = prop*1.5
    }
    h = h + 2
}
plot_data$bin = as.factor(plot_data$bin)
p7 = ggplot(plot_data) + 
    geom_density_ridges2(aes(x=pp, y=h, group=h, fill=correct, color=correct, scale=s), 
                        panel_scaling=TRUE, rel_min_height=0.005, alpha=0.6,
                        jittered_points = TRUE, #stat = "density",height = ..density..,
                        bandwidth = 0.03,
    position = position_points_jitter(width=0.03, height=0.1),
    point_size = 1, point_alpha = 0.4) + 
    ylab('Distance to\nnearest diploid') +
    xlab('Posterior probability\n  ') +
    xlim(c(0.0,1.0)) +
    geom_segment(data=data.frame(y=c(0.25,2.25,4.25,6.25)),
                 aes(x=0.0, xend=1.0, y=y, yend=y), orientation='y', 
                 linetype='dotted', alpha=0.5, color='grey80') +
    scale_y_continuous(breaks=c(0.25,2.25,4.25,6.25), labels=c('0.0-0.05', '0.05-0.25', '0.25-0.5', '0.5-1.0'), limits=c(0,7.5)) +
    scale_fill_manual(values=c('firebrick','grey34')) +
    scale_color_manual(values=c('firebrick','grey34')) +
    theme_classic() +
    theme(legend.position = "none") 

plot_data$h2 = as.factor(plot_data$h2)
pd = aggregate(as.numeric(as.logical(plot_data$correct)), by=list(y=plot_data$h2), FUN=mean)
pd2 = aggregate(as.numeric(as.logical(plot_data$correct)), by=list(y=plot_data$h2), FUN=var)
pd$y = as.numeric(as.character(pd$y))
pd$z = as.numeric(as.character(pd2$x))
pd$zstart = pd$x - pd$z
pd$zend = pd$x + pd$z
if (nrow(pd[pd$zstart < 0.0,]) > 0)
    pd[pd$zstart < 0.0,]$zstart = 0.0
if (nrow(pd[pd$zend > 1,]) > 0)
    pd[pd$zend > 1,]$zend = 1
p8 = ggplot(pd) + 
    geom_segment(aes(x=0.0, xend=1.0, y=y, yend=y), orientation='y', 
                 linetype='dotted', alpha=0.5, color='grey80') +
    geom_segment(aes(x=1.0, xend=1.0, y=0, yend=7.5), orientation='y', 
                 linetype='dotted', alpha=0.4, color='grey80') +
    geom_segment(aes(x=0.0, xend=0.0, y=0, yend=7.5), orientation='y', 
                 linetype='dotted', alpha=0.4, color='grey80') +
    geom_segment(aes(x=zstart, xend=zend, y=y, yend=y), 
                 orientation='y', alpha=0.5, size=1, color='darkorange2') +
    geom_point(aes(y=y, x=x), color='grey44') + 
    geom_line(aes(y=y, x=x), orientation='y', color='grey44') + 
    #geom_smooth(aes(y=y, x=x), method='loess',se=FALSE) + 
    xlab('Proportion\ncorrect') +
    ylim(c(0.0,7.5)) +
    scale_x_continuous(breaks=c(0.0,1.0), limits=c(-0.01,1.01), labels=c('0.0','1.0')) +
    theme_classic() +
    theme(legend.position = "none",
          axis.line.y=element_blank(),
          axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 
#p = grid.arrange(p1, p2, nrow=1, widths=c(4,1))
#ggsave('results/diploid_distance.pdf', p, width=5, height=2.5) 

p = grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow=4, widths=c(4,1))
ggsave('results/combined_plots.pdf', p, width=5, height=10) 
