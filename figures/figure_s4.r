#!/usr/bin/env Rscript

library(ape)
library(phangorn)
library(magrittr)

create_tree_for_plotting = function(dataset, tree_files, rooted_tree_files, keep=NULL) {

    ## Read the tree that was rooted by using FigTree
    rtree = read.tree(rooted_tree_files[[dataset]])
    rtree[['tip.label']] = gsub('\\.fasta', '', rtree[['tip.label']])

    ## Read the tree directly from parsnp, which contains the bootstrap
    ## values as labels
    btree = read.tree(tree_files[[dataset]])
    btree[['tip.label']] = gsub('\\.fasta', '', btree[['tip.label']])

    ## Transfer bootstraps to rooted tree
    rbtree = phangorn::addConfidences(rtree, btree)

    ## Extract the desired tips
    if(!is.null(keep)) {
        drop.tip(rbtree, rbtree[['tip.label']][!(rbtree[['tip.label']] %in% keep)],
                 trim.internal=TRUE, collapse.singles=TRUE)
    } else {
        rbtree
    }
}

tree_dir = '/home/tiffinp/epste051/project/jgi_sequencing/results/ortho/medmel/sibeliaz/2020-04-09_split/parsnp-tree' 
tree_files = sapply(c('full'='output/parsnp.tree',
                      'chromosome'='output-chromosome/parsnp.tree',
                      'psyma'='output-pSymA/parsnp.tree',
                      'psymb'='output-pSymB/parsnp.tree'),
                    function(f) file.path(tree_dir, f))
rooted_tree_files = sapply(c('full'='rooted_by_figtree_2020-12-26.nw',
                             'chromosome'='chromosome-rooted_by_figtree_2021-07-27.nw',
                             'psyma'='pSymA-rooted_by_figtree_2021-07-27.nw',
                             'psymb'='pSymB-rooted_by_figtree_2021-07-27.nw'),
                           function(f) file.path(tree_dir, f))

rpcs = read.csv('/home/tiffinp/epste051/project/select_reseq/results/phenotypes/gwas_phenotypes/spanfran2/2020-11-13_pca/strain_pca_scores.tsv',
                sep='\t')


## Color palette for plotting phenotypes on the tree
pal = c("#3a4d6b", "#3d6da2", "#799a96", "#ccbe6a", "#ffec99")
cr = colorRamp(pal)

## Color palette for the community assignments
sp2_tree_colors = c('BG'='gray',
                    'SP2'='darkred')
all_tree_colors = c('BG'='gray80',
                    'SP2'='#c44245',
                    'C68'='#00b6f1',
                    'C101'='#d9bf0d',
                    'Both'='#6a28c7')
pchs = c('SP2'=19,
         'C101'=15,
         'C68'=17,
         'Both'=8)
tip_offset = 0.0015

## Use color palette to assign a color to each strain
x = rpcs[, 'rhizobia_score_PC1']
pc1 = structure((x - min(x)) / (max(x) - min(x)), names=rpcs[, 'strain'])
#pc1 = structure(rank(x) / length(x), names=rpcs[, 'strain'])
cl1 = structure(apply(cr(pc1), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255)),
                names=names(pc1))
x = rpcs[, 'rhizobia_score_PC2']
pc2 = structure((x - min(x)) / (max(x) - min(x)), names=rpcs[, 'strain'])
#pc2 = structure(rank(x) / length(x), names=rpcs[, 'strain'])
cl2 = structure(apply(cr(pc2), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255)),
                names=names(pc2))

keep = scan('/home/tiffinp/epste051/project/jgi_sequencing/results/ortho/medmel/sibeliaz/2020-04-09_split/meliloti/strains.txt',
            what='character') %>%
    .[grepl('^MAG', .)] %>%
    gsub('__.+', '', .)
keep2 = scan('/home/tiffinp/epste051/project/jgi_sequencing/results/ortho/medmel/sibeliaz/2020-04-09_split/meliloti/strains.txt',
            what='character') %>%
    gsub('__.+', '', .) %>%
    .[!grepl('USDA1025', .)]
gwa_panel = scan('/home/tiffinp/epste051/project/select_reseq/data/strain/lists/spanfran2_meliloti.txt',
                 what='character')

panels = list('SP2'=gwa_panel,
              'C68'=scan('/home/tiffinp/epste051/project/select_reseq/data/strain/lists/68_strains.txt', what='character'),
              'C101'=scan('/home/tiffinp/epste051/project/select_reseq/data/strain/lists/101Strains.txt', what='character')
              )
panels[['C101']] = ifelse(panels[['C101']] %in% c('3011', '2019', '3085'),
                          panels[['C101']], paste0('USDA', panels[['C101']]))
panels[['C68']] = ifelse(panels[['C68']] %in% c('3011', '2019', '3085') | is.na(as.numeric(panels[['C68']])),
                          panels[['C68']], paste0('USDA', panels[['C68']]))

sp2_tree = create_tree_for_plotting('full', tree_files, rooted_tree_files, keep)
replicon_trees = lapply(c('chromosome', 'psyma', 'psymb'),
                        function(d) create_tree_for_plotting(d, tree_files, rooted_tree_files, keep))
all_tree = create_tree_for_plotting('full', tree_files, rooted_tree_files, keep2)
all_replicon_trees = lapply(c('chromosome', 'psyma', 'psymb'),
                            function(d) create_tree_for_plotting(d, tree_files, rooted_tree_files, keep2))
names(all_replicon_trees) = c('chromosome', 'psyma', 'psymb')
names(replicon_trees) = c('chromosome', 'psyma', 'psymb')


## File with SpanFran2 community, just the trees
pdf('SP2_tree.2021-07-23.pdf', width=8, height=8)

plot(sp2_tree, type='fan', 
     tip.color=ifelse(sp2_tree[['tip.label']] %in% gwa_panel, sp2_tree_colors['SP2'], sp2_tree_colors['BG']),
     no.margin=TRUE, cex=0.6, label.offset=0.0005)
bts = as.numeric(sp2_tree[['node.label']])
cl = ifelse(bts > 0.95, 'grey20', ifelse(bts > 0.80, 'gray70', NA))
nodelabels(pch=19, cex=0.75, col=cl)
add.scale.bar(cex=1.5)

for(r in names(replicon_trees)) {
    tr = replicon_trees[[r]]
    plot(tr, type='fan', 
         tip.color=ifelse(tr[['tip.label']] %in% gwa_panel, sp2_tree_colors['SP2'], sp2_tree_colors['BG']),
         no.margin=TRUE, cex=0.6, label.offset=0.0005)
    bts = as.numeric(tr[['node.label']])
    cl = ifelse(bts > 0.95, 'grey20', ifelse(bts > 0.80, 'gray70', NA))
    nodelabels(pch=19, cex=0.75, col=cl)
    add.scale.bar(cex=1.5)
}

dev.off()


#-----------------------------------------------------------------------

## File with phenotypes plotted onto the SpanFran2 trees
pdf('SP2_tree_PC1.2021-07-23.pdf', width=8, height=8)

plot(sp2_tree, type='fan', root.edge=TRUE,
     show.tip.label=FALSE,
     no.margin=TRUE, cex=0.6, label.offset=0.0005)
tiplabels(pch=19, col=cl1[sp2_tree[['tip.label']]])

for(r in names(replicon_trees)) {
    tr = replicon_trees[[r]]
    #tr2 = drop.tip(tr, tr[['tip.label']][!(tr[['tip.label']] %in% gwa_panel)],
    #               trim.internal=TRUE, collapse.singles=TRUE)
    tr2 = tr
    plot(tr2, type='fan', root.edge=TRUE,
         show.tip.label=FALSE,
         no.margin=TRUE, cex=0.6, label.offset=0.0005)
    tiplabels(pch=19, col=cl1[tr2[['tip.label']]])
}
dev.off()


pdf('SP2_tree_PC2.2021-07-23.pdf', width=8, height=8)

plot(sp2_tree, type='fan', root.edge=TRUE,
     show.tip.label=FALSE,
     no.margin=TRUE, cex=0.6, label.offset=0.0005)
tiplabels(pch=19, col=cl2[sp2_tree[['tip.label']]])

for(r in names(replicon_trees)) {
    tr = replicon_trees[[r]]
    #tr2 = drop.tip(tr, tr[['tip.label']][!(tr[['tip.label']] %in% gwa_panel)],
    #               trim.internal=TRUE, collapse.singles=TRUE)
    tr2 = tr
    plot(tr2, type='fan', root.edge=TRUE,
         show.tip.label=FALSE,
         no.margin=TRUE, cex=0.6, label.offset=0.0005)
    tiplabels(pch=19, col=cl2[tr2[['tip.label']]])
}

dev.off()


#-----------------------------------------------------------------------


## File with all the strains
pdf('full_tree.2021-07-23.pdf', width=15, height=15)

plot(all_tree, type='fan', 
     tip.color=ifelse(all_tree[['tip.label']] %in% gwa_panel, sp2_tree_colors['SP2'], sp2_tree_colors['BG']),
     no.margin=TRUE, cex=0.6, label.offset=0.0005,
     main='All replicons combined')

#bts = as.numeric(all_tree[['node.label']])
#cl = ifelse(bts > 0.95, 'grey20', ifelse(bts > 0.80, 'gray70', NA))
#nodelabels(pch=19, cex=0.75, col=cl)
add.scale.bar(cex=1.5)

i = 0
for(community in names(panels)) {
    tiplabels(pch=pchs[community],
              col=ifelse(all_tree[['tip.label']] %in% panels[[community]],
                         all_tree_colors[community], all_tree_colors['BG']),
              offset=i*tip_offset)
    i = i + 1
}


for(r in names(all_replicon_trees)) {
    tr = all_replicon_trees[[r]]
    o = max(cophenetic(tr))


    color = structure(character(length(tr[['tip.label']])), names=tr[['tip.label']])
    pch = structure(numeric(length(tr[['tip.label']])), names=tr[['tip.label']])
    for(i in seq_along(color)) {
        s = names(color)[i]
        if(s %in% panels[['SP2']]) {
            color[i] = sp2_tree_colors[['SP2']]
            pch[i] = NA
        } else if(s %in% panels[['C101']] && s %in% panels[['C68']]) {
            pch[i]  = pchs['Both']
            color[i] = all_tree_colors[['Both']]
        } else if(s %in% panels[['C101']]) {
            pch[i]  = pchs['C101']
            color[i] = all_tree_colors[['C101']]
        } else if(s %in% panels[['C68']]) {
            pch[i]  = pchs['C68']
            color[i] = all_tree_colors[['C68']]
        } else {
            color[i] = sp2_tree_colors[['BG']]
            pch[i] = NA
        }
    }

    plot(tr, type='fan', 
         tip.color=color,
         no.margin=TRUE, cex=0.6,
         main=r, label.offset=o/300)
    bts = as.numeric(tr[['node.label']])
    cl = ifelse(bts > 0.95, 'grey20', ifelse(bts > 0.80, 'gray70', NA))
    nodelabels(pch=19, cex=0.75, col=cl)
    tiplabels(offset=strwidth('USDA0000') * 0.75, pch=pch, col=color)
    add.scale.bar(cex=1.5)

}

dev.off()
