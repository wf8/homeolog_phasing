#
# Plot distances between gene trees simulated under the MSC.
#
# Tree distance is calculated as both RF and as in Kuhner and Felsenstein (1994);
# similar to RF distance but taking branch lengths into account.
#
library(ape)
library(ggplot2)

num_loci = 4
num_species = 20
num_reps = 100
pop_sizes = c(0.1, 1.0, 10, 100, 1000)

plot_data = data.frame()
for (sim_num in 1:num_reps) {
    for (pop_size in pop_sizes) {
        species_tree = read.tree(paste0(pop_size, '/', sim_num, '/species_tree.tree'))
        species_tree_unrooted = unroot(species_tree)
        for (locus in 1:num_loci) {
            gene_tree_path = paste0(pop_size, '/', sim_num, '/gene_tree', locus, '.tree')
            gene_tree = read.tree(gene_tree_path)
            d_KF = dist.topo(x=species_tree, y=gene_tree, method="score")[[1]]
            gene_tree_unrooted = unroot(gene_tree)
            d_RF = dist.topo(x=gene_tree_unrooted, y=species_tree_unrooted)[[1]]
            d = data.frame(pop_size=pop_size,
                           d_KF=d_KF,
                           d_RF=d_RF)
            plot_data = rbind(plot_data, d)
        }
    }
}
for (sim_num in 1:(num_reps-1)) {
    for (pop_size in pop_sizes) {
        species_tree1 = read.tree(paste0(pop_size, '/', sim_num, '/species_tree.tree'))
        species_tree_unrooted1 = unroot(species_tree1)
        species_tree2 = read.tree(paste0(pop_size, '/', sim_num + 1, '/species_tree.tree'))
        species_tree_unrooted2 = unroot(species_tree2)
        d_KF = dist.topo(x=species_tree1, y=species_tree2, method="score")[[1]]
        d_RF = dist.topo(x=species_tree_unrooted1, y=species_tree_unrooted2)[[1]]
        d = data.frame(pop_size='unlinked',
                       d_KF=d_KF,
                       d_RF=d_RF)
        plot_data = rbind(plot_data, d)
    }
}
plot_data$pop_size = as.factor(plot_data$pop_size)
p = ggplot(plot_data) + 
    geom_point(aes(x=d_RF, y=log(d_KF), color=pop_size), alpha=0.5) + 
    theme_classic()
ggsave('pop_size_and_tree_discordance.pdf', p, width=5, height=4) 
