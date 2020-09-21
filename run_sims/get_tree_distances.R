#
# Tree distance is calculated as both RF and as in Kuhner and Felsenstein (1994);
# similar to RF distance but taking branch lengths into account.
#
library(ape)

species_tree = read.tree('species_tree.tree')
species_tree_unrooted = unroot(species_tree)
KF = vector()
RF = vector()
num_loci = 4
for (locus in 1:num_loci) {
    gene_tree_path = paste0('gene_tree', locus, '.tree')
    gene_tree = read.tree(gene_tree_path)
    for (i in 1:length(gene_tree$tip.label)) {
        gene_tree$tip.label[i] = substring(gene_tree$tip.label[i], 1, nchar(gene_tree$tip.label[i])-2)
    }
    d_KF = dist.topo(x=species_tree, y=gene_tree, method="score")[[1]]
    KF = c(KF, d_KF)
    gene_tree_unrooted = unroot(gene_tree)
    d_RF = dist.topo(x=gene_tree_unrooted, y=species_tree_unrooted)[[1]]
    RF = c(RF, d_RF)
}
cat(paste0(mean(KF), ',', mean(RF)))
