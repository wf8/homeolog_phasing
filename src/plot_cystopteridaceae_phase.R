#
#       tip1_map, tip1_pp, tip2_map, tip2_pp, ....
# gene1,
# gene2,
# gene3,
#
#

library(ggplot2)
library(plyr)
library(magrittr)
library(tidyr)
library(dplyr)
library(ggtree)


tree_file = 'src/output/cystopteridaceae_v2_map_rooted.tre'
input_dir = 'src/output/cystopteridaceae_v2__'
output_dir = 'src/output/'

# names of the loci in the log file
loci = c('APP', 'GAP', 'IBR', 'PGI')

# what percentage of MCMC samples to exclude?
burnin = 0.6

# modified from ggtree gheatmap
homologized = function (p, data, data_labels, offset = 0, width = 1, low = "green", mid, high = "red",
    color = "white", colnames = TRUE, colnames_position = "bottom",
    colnames_angle = 0, colnames_level = NULL, colnames_offset_x = 0,
    colnames_offset_y = 0, font.size = 4, family = "", hjust = 0.5,
    legend_title = "value")
{
    colnames_position %<>% match.arg(c("bottom", "top"))
    variable <- value <- lab <- y <- NULL
    width <- width * (p$data$x %>% range(na.rm = TRUE) %>% diff)/ncol(data)
    isTip <- x <- y <- variable <- value <- from <- to <- NULL
    df <- p$data
    nodeCo <- intersect(df %>% filter(is.na(x)) %>% select(.data$parent,
        .data$node) %>% unlist(), df %>% filter(!is.na(x)) %>%
        select(.data$parent, .data$node) %>% unlist())
    labCo <- df %>% filter(.data$node %in% nodeCo) %>% select(.data$label) %>%
        unlist()
    selCo <- intersect(labCo, rownames(data))
    isSel <- df$label %in% selCo
    df <- df[df$isTip | isSel, ]
    start <- max(df$x, na.rm = TRUE) + offset
    dd <- as.data.frame(data)
    dd2 <- as.data.frame(data_labels)
    i <- order(df$y)
    i <- i[!is.na(df$y[i])]
    lab <- df$label[i]
    dd <- dd[match(lab, rownames(dd)), , drop = FALSE]
    dd2 <- dd2[match(lab, rownames(dd2)), , drop = FALSE]
    dd$y <- sort(df$y)
    dd2$y <- sort(df$y)
    dd$lab <- lab
    dd2$lab <- lab
    dd <- gather(dd, variable, value, -c(lab, y))
    dd2 <- gather(dd2, variable, value, -c(lab, y))
    i <- which(dd$value == "")
    if (length(i) > 0) {
        dd$value[i] <- NA
        dd2$value[i] <- NA
    }
    if (is.null(colnames_level)) {
         dd$variable <- factor(dd$variable, levels = colnames(data))
    }
    else {
        dd$variable <- factor(dd$variable, levels = colnames_level)
    }
    V2 <- start + as.numeric(dd$variable) * width
    mapping <- data.frame(from = dd$variable, to = V2)
    mapping <- unique(mapping)
    dd$x <- V2
    dd2$x <- V2
    dd$width <- width
    dd2$width <- width
    dd[[".panel"]] <- factor("Tree")
    dd2[[".panel"]] <- factor("Tree")
    if (is.null(color)) {
        p2 <- p + geom_tile(data = dd, aes(x, y, fill = value),
            width = width, inherit.aes = FALSE)
    }
    else {
        p2 <- p + geom_tile(data = dd, aes(x, y, fill = value),
            width = width, color = color, inherit.aes = FALSE)
        p2 <- p2 + geom_text(data = dd2, aes(x, y, label=value), size=1, inherit.aes = FALSE)
        
        # TODO
        #print(dd)
        dd3 = data.frame()
        start_x = max(dd$x)
        height = max(dd$y)
        margin = 0.006
        for (y in unique(dd$y)) {
            pp = mean(dd[dd$y == y, 'value'], na.rm=TRUE)
            dd4 = data.frame(pp = pp, x = pp/200 + margin + start_x, y=y)
            dd3 = rbind(dd3, dd4)
        }
        #print(dd3)
        p2 <- p2 + geom_segment(aes(x=start_x+margin, xend=1/200 + start_x + margin, y=0.2, yend=0.2), size=0.5, inherit.aes = FALSE)
        p2 <- p2 + geom_segment(aes(x=1/200+start_x+margin, xend=1/200+start_x+margin, y=0.5, yend=height), color='grey85', linetype='dotted', size=0.35, inherit.aes = FALSE)
        p2 <- p2 + geom_segment(aes(x=start_x+margin, xend=start_x+margin, y=0.5, yend=height), color='grey85', linetype='dotted', size=0.35, inherit.aes = FALSE)
        p2 <- p2 + geom_point(data = dd3, aes(x, y, color=pp), size=1.25, inherit.aes = FALSE, show.legend=FALSE)
        p2 <- p2 + geom_text(label='0.0', x=start_x+margin, y=-0.2, size=1.25, color='grey50')
        p2 <- p2 + geom_text(label='1.0', x=start_x+margin+1/200, y=-0.2, size=1.25, color='grey50')
    }
    if (is(dd$value, "numeric")) {
        midpoint = max(dd$value, na.rm=TRUE) - min(dd$value, na.rm=TRUE)
        midpoint = midpoint/2 + min(dd$value, na.rm=TRUE) 
        p2 <- p2 + scale_fill_gradient2(low = low, mid=mid, high = high, midpoint=midpoint,
            na.value = "white", name = legend_title)
        p2 <- p2 + scale_color_gradient2(low = low, mid=mid, high = high, midpoint=midpoint,
            na.value = "white", name = legend_title)
            #na.value = NA, name = legend_title)
    }
    else {
        p2 <- p2 + scale_fill_discrete(na.value = NA, name = legend_title)
    }
    if (colnames) {
        if (colnames_position == "bottom") {
            y <- 0
        }
        else {
            y <- max(p$data$y) + 1
        }
        mapping$y <- y
        mapping[[".panel"]] <- factor("Tree")
        p2 <- p2 + geom_text(data = mapping, aes(x = to, y = y,
            label = from), size = font.size, family = family,
            inherit.aes = FALSE, angle = colnames_angle, nudge_x = colnames_offset_x,
            nudge_y = colnames_offset_y, hjust = hjust)
    }
    p2 <- p2 + theme(legend.position = "right")
    if (!colnames) {
        p2 <- p2 + scale_y_continuous(expand = c(0, 0))
    }
    attr(p2, "mapping") <- mapping
    return(p2)
}



# polyploid samples and their tips
samples = list('xCystG_7974' = c('xCystG_7974_A', 'xCystG_7974_B'),
               'xCystC_7974' = c('xCystC_7974_A', 'xCystC_7974_B'),
               'G_rob_7945' = c('G_rob_7945_A', 'G_rob_7945_B'),
               'G_rem_4862' = c('G_rem_4862_A', 'G_rem_4862_B'),
               'G_oyaC_8739' = c('G_oyaC_8739_A', 'G_oyaC_8739_B'),
               'G_dry_7981' = c("G_dry_7981a_A","G_dry_7981a_B","G_dry_7981_A","G_dry_7981_B"),
               'G_dis_7751' = c("G_dis_7751a_A","G_dis_7751_A"),
               'G_con_6979' = c("G_con_6979_A","G_con_6979_B"),
               'C_uta_6848' = c("C_uta_6848_A","C_uta_6848_B"),
               'C_tenu_6387' = c("C_tenu_6387_A","C_tenu_6387_B"),
               'C_tas_6379' = c("C_tas_6379_C","C_tas_6379_B","C_tas_6379_A"),
               'C_sudB_8674' = c("C_sudB_8674_A","C_sudB_8674_B"),
               'C_pelA_6055' = c("C_pelA_6055_A","C_pelA_6055_B"),
               'C_mon_7943' = c("C_mon_7943_B","C_mon_7943_A"),
               'C_fraB_7248' = c("C_fraB_7248_A","C_fraB_7248_B"),
               'C_fraA_7009' = c("C_fraA_7009_A","C_fraA_7009_B"),
               'C_dia_6380' = c("C_dia_6380_A","C_dia_6380_B"),
               'A_tenC_8745' = c("A_tenC_8745_A","A_tenC_8745_B"),
               'A_tenB_8704' = c("A_tenB_8704_A","A_tenB_8704_B"),
               'A_tai_6137' = c("A_tai_6137_A","A_tai_6137_B"))

# populate empty dataframes to hold results
map_prob_results = data.frame()
joint_map_phase_results = data.frame()
for (sample in names(samples)) {
    sample_joint_map_prob = data.frame()
    sample_joint_map_phase = data.frame()
    for (i in 1:length(loci)) {
        sample_joint_map_prob[1, loci[i]] = 0.0
        sample_joint_map_phase[1, loci[i]] = ''
        row.names(sample_joint_map_prob) = c(sample)
        row.names(sample_joint_map_phase) = c(sample)
    }
    map_prob_results = rbind(map_prob_results, sample_joint_map_prob)
    joint_map_phase_results = rbind(joint_map_phase_results, sample_joint_map_phase)
}

# for each sample loop over each locus
marginal_results = data.frame()
for (sample in names(samples)) {
    joint_results = data.frame()
    for (i in 1:length(loci)) {

        # read in file and exclude burnin
        f_in = paste0(input_dir, 'phase', i, '.log')
        d = read.csv(f_in, sep='\t')
        d = d[floor(nrow(d)*burnin):nrow(d),]

        # get joint phase assignments for this locus
        d1 = d[, samples[[sample]]]
        joint_results_locus = as.data.frame(table(d1))
        joint_results_locus$joint_prob = joint_results_locus$Freq / sum(joint_results_locus$Freq)
        joint_results_locus$locus = loci[i]
        joint_results = rbind(joint_results, joint_results_locus) 
    
        # get the MAP joint phase for the plot
        map = which(joint_results_locus['joint_prob'] == max(joint_results_locus['joint_prob']))
        for (tip in samples[[sample]]) {
            #map_prob_results[tip,loci[i]] = joint_results_locus[map, 'joint_prob']
            joint_map_phase_results[tip,loci[i]] = as.character(joint_results_locus[map, tip])
        }

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
    out_file = paste0(output_dir, '/joint_phase_probs_', sample, '.csv')
    write.csv(joint_results, out_file, row.names=FALSE)
}
out_file = paste0(output_dir, '/marginal_phase_probs.csv')
write.csv(marginal_results, out_file, row.names=FALSE)

# get marginal probs for the joint MAP phase
for (sample in names(samples)) {
    for (tip in samples[[sample]]) {
        for (i in 1:length(loci)) {
            m = marginal_results[marginal_results$phase == joint_map_phase_results[tip,loci[i]] &
                                 marginal_results$locus == loci[i] &
                                 marginal_results$tip == tip, 'marginal_prob']
            map_prob_results[tip,loci[i]] = m
        }
    }
}

tree = treeio::read.beast(tree_file)
p = ggtree(tree) 
p = p + geom_tiplab(size=2, align=T, linesize=0.25, offset=0.0005)  
p = homologized(p, map_prob_results, joint_map_phase_results, 
                offset=0.012, low="#EE0000", mid="#FF0099", high="#DDDDFF", 
             colnames_position="top", font.size=2, width=0.5,
             legend_title="Posterior\nProbability") 
p = p + theme(legend.text=element_text(size=6),
              legend.title=element_text(size=8))
#p = p + scale_fill_viridis_c(option="D", name="Posterior\nProbability", na.value="white")
#p = p + scale_color_viridis_c(option="D", name="Posterior\nProbability", na.value="white")
ggsave('homologized_joint_MAP_v2.pdf')

