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
library(ggtree)

tree_file = 'output_good/cystopteridaceae_rooted.tree'
genes = c('APP', 'GAP', 'IBR', 'PGI')
burnin = 0.1

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
    }
    if (is(dd$value, "numeric")) {
        midpoint = max(dd$value, na.rm=TRUE) - min(dd$value, na.rm=TRUE)
        midpoint = midpoint/2 + min(dd$value, na.rm=TRUE) 
        p2 <- p2 + scale_fill_gradient2(low = low, mid=mid, high = high, midpoint=midpoint,
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


prob_results = data.frame()
seq_results = data.frame()
for (i in 1:length(genes)) {

    f_in = paste0('output_good/cystopteridaceae_phase', i, '.log')
    
    # read in file and exclude burnin
    d = read.csv(f_in, sep='\t')
    d = d[floor(nrow(d)*burnin):nrow(d),]

    prob_locus_results = data.frame()
    seq_locus_results = data.frame()
    variables = names(d)
    for (j in 2:length(variables)) {
        
        variable = variables[j]
        c = plyr::count(d, vars=variable)
        prob = max(c['freq'])/sum(c['freq'])
        map = which(c['freq'] == max(c['freq']))
        map = c[variable][map,]

        #seq_locus_results[1,paste0(variable, '_map')] = map
        #prob_locus_results[1,paste0(variable, '_pp')] = prob
        seq_locus_results[1, variable] = map
        prob_locus_results[1, variable] = prob

    }
    row.names(prob_locus_results) = c(genes[i])
    row.names(seq_locus_results) = c(genes[i])
    prob_results = rbind(prob_results, prob_locus_results)
    seq_results = rbind(seq_results, seq_locus_results)
}
prob_results = t(prob_results)
seq_results = t(seq_results)
#write.csv(results, 'output/cystopteridaceae_phase1.csv')

tree = treeio::read.beast(tree_file)
p = ggtree(tree) 
p = p + geom_tiplab(size=2, align=T, linesize=0.25, offset=0.0005)  
p = homologized(p, prob_results, seq_results, offset=0.009, low="#EE0000", mid="#FF0099", high="#DDDDFF", 
             colnames_position="top", font.size=2, width=0.5,
             legend_title="Posterior\nProbability") 
p = p + theme(legend.text=element_text(size=6),
              legend.title=element_text(size=8))
#+ scale_fill_viridis_c(option="D", name="Posterior\nProbability", na.value="white")
ggsave('homologized.pdf')

