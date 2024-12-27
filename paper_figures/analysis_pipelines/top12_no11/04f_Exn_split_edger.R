setwd("~/cdkwan/alex_kwan/paper_figures/figures/top12_no11/figure_4_exn_split_pseudobulk/")

librarian::shelf(
  kwanlibr,
  edgeR,
  ggplot2,
  ggfortify,
  dplyr,
  DescTools,
  latex2exp,
)

samples = read.csv('pseudobulk_counts/pseudobulk_metadata.csv')
samples = samples %>%
  rename(sampleName = label, files = filepath, cluster = id)

clusters = c('7','2','3','0','5')
cluster.names = c('L5-PT', 'L4-5-IT', 'L4-IT', 'L2-3-IT', 'L5-6-IT')
genes_to_label = c('Camk1d', 'Grik3', 'Grin2a', 'Grin2b', 'Shisa7', 'Shisa9')

for (this.cluster.index in 1:length(clusters)) {
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

      p = kwanlibr::draw_volcano(
        lrt,
        figure_title = paste0(this.time,
                              ' ', this.drug,
                              ' ', cluster.names[this.cluster.index],
                              ' cell DEGs'),
        fdr = 0.05,
        xdiff = 3,
        ymax = 15
      )
      p = label_volcano_plot_genes(p=p, lrt=lrt, label_genes = genes_to_label)

      dir.create(file.path('.', 'void_volcanos'))

      p = p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      kwanlibr::ggsave_vector_raster(
        filename = file.path('.', 'void_volcanos', paste0('volcano_', this.drug,
                                                          '_', this.time,
                                                          '_', cluster.names[this.cluster.index],
                                                          '_aspect')),
        width = 6, height = 6, dpi=600,
        plot = p
      )

      p = p + theme_void()
      kwanlibr::ggsave_vector_raster(
        filename = file.path('.', 'void_volcanos', paste0('volcano_', this.drug,
                                                          '_', this.time,
                                                          '_', cluster.names[this.cluster.index],
                                                          '_void')),
        width = 6, height = 6, dpi=600,
        plot = p
      )
    }
  }
}

excitatory_de = do.call('rbind', lapply(Sys.glob('*Psilo.csv'), FUN=function(f) {
  df = read.csv(f)
  df = df %>% filter(gene_name %in% c(genes_to_label))
}))
excitatory_de %>%
  ggplot(aes(time, logFC, color=cluster)) +
  geom_line() +
  facet_wrap(~ gene_name, ncol=3)
write_csv(excitatory_de, 'excitatory__logFC_pseudobulk.csv')
