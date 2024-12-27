setwd("~/cdkwan/alex_kwan/paper_figures/figures/top12_no11/figure_5/pseudobulk_clusters_counts/")

librarian::shelf(
  kwanlibr,
  edgeR,
  ggplot2,
  ggfortify,
  dplyr,
  DescTools,
  latex2exp,
)

samples = read.csv('pseudobulk_metadata.csv')
samples = samples %>%
  rename(sampleName = label, files = filepath, cluster = id)

clusters = c('8','10','12')
cluster.names = c('CIN-PV', 'CIN-SST', 'CIN-VIP')
genes_to_label = c('Eif4ebp2','Alkbh5','Hnrnpa1','Hnrnpa3','Eif4a1','Ppia')

for (this.cluster.index in 1:3) {
  for (this.time in c('1h','2h','4h','24h','72h')) {
    for (this.drug in c('Ket','Psilo')) {
      these.samples = samples %>% filter(Timepoint %in% c('0h', this.time) &
                                           Drug %in% c('none', this.drug) &
                                           cluster == clusters[this.cluster.index])
      these.samples$group = case_match(these.samples$Drug,
                                       'none' ~ 'ctrl',
                                       this.drug ~ 'trt')
      these.samples$group = factor(these.samples$group, levels=c('ctrl','trt'))

      dge = edgeR::readDGE(these.samples, header=TRUE)
      dge$counts = dge$counts %>% head(-5)

      high_quality = edgeR::filterByExpr(dge, min.count=10)
      dge = dge[ high_quality , , keep.lib.sizes=FALSE]
      dge = edgeR::calcNormFactors(dge)

      design.matrix = model.matrix(~ droplevels(these.samples[["group"]]))
      print(design.matrix)
      dge <- edgeR::estimateDisp(dge, design.matrix)

      fit = edgeR::glmFit(dge, design.matrix)
      lrt = edgeR::glmLRT(fit, coef = 2)
      lrt$table$ctrl.logCPM = cpm(dge, log=TRUE)[,dge$design[,2] == 0] %>%
        rowMeans() %>% as.vector()
      lrt$table$cKO.logCPM = cpm(dge, log=TRUE)[,dge$design[,2] == 1] %>%
        rowMeans() %>% as.vector()
      lrt$table$gene_name = rownames(lrt$table)

      final.table = as.data.frame(edgeR::topTags(lrt, n=nrow(dge$counts)))
      final.table$cluster = clusters[this.cluster.index]
      final.table$time = this.time
      final.table$drug = this.drug
      write.csv(final.table, file=paste0('DEGs_', cluster.names[this.cluster.index], '_', this.time, '_', this.drug, '.csv'))

      kwanlibr::make_volcano(
        lrt,
        figure_title = paste0(this.time, ' ', this.drug, ' ', cluster.names[this.cluster.index], ' cell DEGs'),
        filename = paste0(cluster.names[this.cluster.index], '_', this.time, '_', this.drug),
        figure_dir = '.',
        fdr = 0.05,
        xdiff = 3,
        ymax = 15,
        label_genes = genes_to_label
      )
    }
  }
}

cin_de = do.call('rbind', lapply(Sys.glob('DEGs_CIN-*_*h_Psilo.csv'), FUN=function(f) {
  df = read.csv(f)
  df = df %>% filter(gene_name %in% c(genes_to_label))
}))
cin_de %>%
  ggplot(aes(time, logFC, color=cluster)) +
  geom_line() +
  facet_wrap(~ gene_name, ncol=3)
write_csv(cin_de, 'cin_logFC_pseudobulk.csv')

#####################################################
#         Lumped Samples
#####################################################

lump.samples = read.csv('lump_pseudobulk_metadata.csv')
lump.samples = lump.samples %>%
  rename(sampleName = label, files = filepath)

for (this.time in c('1h','2h','4h','24h','72h')) {
  for (this.drug in c('Ket','Psilo')) {
    these.samples = lump.samples %>% filter(Timepoint %in% c('0h', this.time) &
                                         Drug %in% c('none', this.drug))
    these.samples$group = case_match(these.samples$Drug,
                                     'none' ~ 'ctrl',
                                     this.drug ~ 'trt')
    these.samples$group = factor(these.samples$group, levels=c('ctrl','trt'))

    dge = edgeR::readDGE(these.samples, header=TRUE)
    dge$counts = dge$counts %>% head(-5)

    high_quality = edgeR::filterByExpr(dge, min.count=10)
    dge = dge[ high_quality , , keep.lib.sizes=FALSE]
    dge = edgeR::calcNormFactors(dge)

    design.matrix = model.matrix(~ droplevels(these.samples[["group"]]))
    print(design.matrix)
    dge <- edgeR::estimateDisp(dge, design.matrix)

    fit = edgeR::glmFit(dge, design.matrix)
    lrt = edgeR::glmLRT(fit, coef = 2)
    lrt$table$ctrl.logCPM = cpm(dge, log=TRUE)[,dge$design[,2] == 0] %>%
      rowMeans() %>% as.vector()
    lrt$table$cKO.logCPM = cpm(dge, log=TRUE)[,dge$design[,2] == 1] %>%
      rowMeans() %>% as.vector()
    lrt$table$gene_name = rownames(lrt$table)

    final.table = as.data.frame(edgeR::topTags(lrt, n=nrow(dge$counts)))
    write.csv(final.table, file=paste0('DEGs_CIN_', this.time, '_', this.drug, '.csv'))

    kwanlibr::make_volcano(
      lrt,
      figure_title = paste0(this.time, ' ', this.drug, ' CIN cell DEGs'),
      filename = paste0('CIN_lump_', this.time, '_', this.drug),
      figure_dir = '.',
      fdr = 0.05,
      xdiff = 3,
      ymax = 15,
      label_genes = genes_to_label
    )

    p = kwanlibr::draw_volcano(
      lrt,
      figure_title = paste0(this.time, ' ', this.drug, ' CIN cell DEGs'),
      fdr = 0.05,
      xdiff = 3,
      ymax = 15
    )
    p = label_volcano_plot_genes(p=p, lrt=lrt, label_genes = genes_to_label)

    dir.create(file.path('.', 'void_volcanos'))

    p = p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    kwanlibr::ggsave_vector_raster(
      filename = file.path('.', 'void_volcanos', paste0('volcano_', this.time, '_', this.drug, '_aspect')),
      width = 6, height = 6, dpi=600,
      plot = p
    )

    p = p + theme_void()
    kwanlibr::ggsave_vector_raster(
      filename = file.path('.', 'void_volcanos', paste0('volcano_', this.time, '_', this.drug, '_void')),
      width = 6, height = 6, dpi=600,
      plot = p
    )
  }
}
