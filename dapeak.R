library(Seurat)
library(Signac)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(GenomicFeatures)
library(rtracklayer)

int <- seurat_obj
table(int@active.ident)
DefaultAssay(int) <- "ATAC"

int$ASIP_pos <- ifelse(as.numeric(int@assays$RNA['ASIP',])>0,'ASIP_pos','ASIP_neg')
int$anno_comb <- paste(int$stage,int$annotation,int$ASIP_pos,sep='_')
peaks_kept <- rownames(int)[(sub('-.*$','',rownames(slot(seu[[seu@active.assay]],name = 'counts'))) %in% seqnames(BSgenome.MySpecies.MyLab.Custom3))]
int <- subset(int, features=peaks_kept)
# int <- int[(sub('-.*$','',rownames(slot(seu[[seu@active.assay]],name = 'counts'))) %in% seqnames(BSgenome.MySpecies.MyLab.Custom3)),]

# Load GTF
gtf <- import("../Coturnix_japonica.Coturnix_japonica_2.0.113.gtf")
# Keep only gene-level entries
genes <- gtf#[gtf$type == "gene"]
genes$gene_name[is.na(genes$gene_name)] <- genes$gene_id[is.na(genes$gene_name)]
# Check metadata columns
head(mcols(genes))
Annotation(int) <- genes
AnnotationPlot(object = int, region = c("1-29554-39554"))

logFC_thresh <- 0.5; pct_thresh <- 10
logFC_thresh <- 0.1; pct_thresh <- 5


da_peaks <- presto:::wilcoxauc.Seurat(int, group_by = "anno_comb", assay = 'data', seurat_assay = "ATAC", verbose= TRUE)
da_peaks <-  dplyr::filter(da_peaks, padj < 0.01, logFC > logFC_thresh, pct_in > pct_thresh) %>% arrange(group, -auc) 
ClosestFeature_da_peaks <- ClosestFeature(int, regions = StringToGRanges(da_peaks$feature, sep = c("-", "-")))
head(ClosestFeature_da_peaks)
da_peaks <- cbind(da_peaks,ClosestFeature_da_peaks[match(da_peaks$feature,ClosestFeature_da_peaks$query_region),c('gene_name','gene_biotype','type','distance')])
names(da_peaks)[which(names(da_peaks)=='gene_name')] <- 'nearest_gene'
da_peaks$in_peak_list <- ifelse(da_peaks$feature %in% linked_feature_name, 'TRUE','FALSE')
write.csv(da_peaks,file=paste('da_peaks','logFC',logFC_thresh, 'pct', pct_thresh,'.csv',sep = '_'))
da_peaks[da_peaks$feature %in% linked_feature_name,c('feature','group','logFC','auc','padj','pct_in','pct_out','nearest_gene','gene_biotype','type','distance')]

da_peaks$stage <- sub('_.*$','',da_peaks$group)
da_peaks$cluster <- sub('_.*$','',sub('^[^_]*_','',da_peaks$group))
da_peaks$ASIP <-  sub('^[^_]*_','',sub('^[^_]*_','',da_peaks$group))
da_peaks$value <- 1 
da_peaks_summary <- da_peaks %>%
  group_by(stage, cluster, ASIP) %>%
  summarise(Freq = n(), .groups = "drop")
p_dap <- ggplot(da_peaks_summary, aes(x = cluster, y = Freq, fill = ASIP)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ stage) + ylab('# of differential peaks') + 
  cowplot::theme_cowplot()
ggsave(p_dap,file=paste('number_of_da_peaks','logFC',logFC_thresh, 'pct', pct_thresh,'.pdf',sep = '_'),width = 6, height = 4.5)


de_genes <- presto:::wilcoxauc.Seurat(int, group_by = "annotation", assay = 'data', seurat_assay = "RNA", verbose= TRUE)
de_genes <-  dplyr::filter(de_genes, padj < 0.01, logFC > 0.1, pct_in > 5) %>% arrange(group, -auc)

#### peak-gene relation
library(data.table);library(GenomicRanges)
library(Rcpp)
sourceCpp("~/Documents/GitHub/CCI2TF2CCI/src/row_correlation.cpp") # download from https://github.com/GreenleafLab/ArchR/blob/master/src/Correlation.cpp
sourceCpp("~/Documents/GitHub/CCI2TF2CCI/src/knn.cpp") # download from https://github.com/GreenleafLab/ArchR/blob/master/src/KNN_Utils.cpp

CollapseToLongestTranscript <- function(ranges) {
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]],
        gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c(
    "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
  )
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}
gene.coords <- CollapseToLongestTranscript(ranges = Annotation(object = int[['ATAC']]))
peak.coords <- StringToGRanges(rownames(int@assays$ATAC@counts), sep = c("-", "-"))
geneStart <- GRanges(gene.coords@seqnames, IRanges(gene.coords@ranges@start, width = 1), name = gene.coords$gene_name)
genome.info <- data.frame(Chrom=as.character(gene.coords@seqnames), Starts= gene.coords@ranges@start, 
                          Ends=gene.coords@ranges@start+1, genes= gene.coords$gene_name)
genome.info <- genome.info[!duplicated(genome.info$genes),]

# define functions
#Get Peak Set
DefaultAssay(int) <- 'ATAC' #'shared_ATAC'
peak.matrix <- slot(int[[int@active.assay]],name = 'counts') # counts slot
idx.keep <- (rowSums(x = peak.matrix) > 0)
rna.matrix <-slot(int[['RNA']], name= 'data')
idx.keep.rna <- which(rowSums(x = slot(seu[['RNA']], name= 'counts')) > 100)

myArchR_peak2geneLinkage <- function(seu,peak.assay='ATAC',peak.slot='counts',rna.assay='RNA',rna.slot='data',
                                     corCutOff = 0.45,
                                     FDRCutOff = 1e-04){
  #Get Peak Set
  DefaultAssay(seu) <- peak.assay #'shared_ATAC'
  peak.matrix <- slot(seu[[seu@active.assay]],name = peak.slot) # counts slot
  #idx.keep <- (rowSums(x = peak.matrix) > 0)
  peak.matrix <- peak.matrix[idx.keep, , drop = FALSE]
  peak.ranges <- StringToGRanges(rownames(peak.matrix), sep = c("-", "-"))
  mcols(peak.ranges)$peak_id <- rownames(peak.matrix)
  peakSet <- peak.ranges
  
  #Gene Info
  RNA_assay_name <- rna.assay #'RNA'
  rna.matrix <-slot(seu[[RNA_assay_name]], name= rna.slot)
  #rna.matrix <- rna.matrix[rowSums(x = slot(seu[[RNA_assay_name]], name= 'counts')) > 100, , drop = FALSE]
  rna.matrix <- rna.matrix[idx.keep.rna,, drop = FALSE]
  geneSet <- genome.info[which(genome.info$genes %in% rownames(rna.matrix)),]
  rna.matrix <- rna.matrix[geneSet$genes,]
  geneStart <- GRanges(geneSet$Chrom, IRanges(geneSet$Starts, width = 1), name = geneSet$genes, idx = 1:dim(geneSet)[1])
  
  #Check Chromosomes
  chri <- gtools::mixedsort(unique(paste0(seqnames(peakSet))))
  chrj <- gtools::mixedsort(unique(paste0(seqnames(geneStart))))
  chrij <- intersect(chri, chrj)
  
  #Get Reduced Dims
  rD <- seu@reductions$integrated_lsi@cell.embeddings
  #Subsample
  idx <- sample(seq_len(nrow(rD)), 500, replace = !nrow(rD) >= 500)
  #KNN Matrix
  data = rD; query = rD[idx,]; k = min(100,round(1/10*dim(rD)[1]))
  library(nabor)
  knnObj <- nabor::knn(data = data, query = query, k = k)$nn.idx
  overlapCutoff = 0.8
  
  
  keepKnn <- determineOverlapCpp(knnObj, floor(overlapCutoff * k))
  knnObj <- knnObj[keepKnn==0,]
  #Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
    rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList
  
  #Group Matrix RNA
  rna_mask <- sapply(seq_len(length(knnObj)), function(x) colnames(rna.matrix) %in%
                       knnObj[[x]])
  rna_mask <- Matrix::Matrix(rna_mask)
  rna_new <- as.matrix(rna.matrix %*% rna_mask)
  
  #Group Matrix ATAC
  atac_mask <- sapply(seq_len(length(knnObj)), function(x) colnames(peak.matrix) %in%
                        knnObj[[x]])
  atac_mask <- Matrix::Matrix(atac_mask)
  atac_new <- as.matrix(peak.matrix %*% atac_mask)
  
  rna_new <- t(t(rna_new) / colSums(rna_new)) * 1e6
  atac_new <- t(t(atac_new) / colSums(atac_new)) * 1e6
  
  rna_new <- log(rna_new + 1)
  atac_new <- log(atac_new + 1)
  
  #Overlaps
  o <- DataFrame(findOverlaps( resize(geneStart, 2 * 250000 + 1, "center"), peakSet, ignore.strand = TRUE))
  
  #Get Distance from Fixed point A B
  o$distance <- distance(geneStart[o[,1]] , peakSet[o[,2]] )
  colnames(o) <- c("B", "A", "distance")
  
  # Computing Correlations
  o$Correlation <- rowCorCpp(as.integer(o$A), as.integer(o$B), atac_new, rna_new)
  o$VarAssayA <- ArchR:::.getQuantiles(matrixStats::rowVars(atac_new))[o$A]
  o$VarAssayB <- ArchR:::.getQuantiles(matrixStats::rowVars(rna_new))[o$B]
  o$TStat <- (o$Correlation / sqrt((pmax(1-o$Correlation^2, 0.00000000000000001, na.rm = TRUE))/(ncol(atac_new)-2))) #T-statistic P-value
  o$Pval <- 2*pt(-abs(o$TStat), ncol(atac_new) - 2)
  o$FDR <- p.adjust(o$Pval, method = "fdr")
  out <- o[, c("A", "B", "Correlation", "FDR", "VarAssayA", "VarAssayB")]
  colnames(out) <- c("idxATAC", "idxRNA", "Correlation", "FDR", "VarQATAC", "VarQRNA")
  out$gene <- rownames(rna_new)[out$idxRNA]
  out$peak <- rownames(atac_new)[out$idxATAC]
  out <- out[,c(1,8,2,7,3:6)]
  out <- out[,c(2,4,5,6)]
  # out <- out[which(out$FDR <= FDRCutOff),] 
  out <- as.data.frame(out)
  # corCutOff = 0.45; FDRCutOff = 1e-04
  #out <- out[which(out$Correlation>=corCutOff & out$FDR <= FDRCutOff & out$VarQATAC >=varCutOffATAC & out$VarQRNA >= varCutOffRNA),]
  # out <- out[which(out$Correlation>=corCutOff & out$FDR <= FDRCutOff),]
  return(out)}
out <- myArchR_peak2geneLinkage(int)
linked_feature_name <- c("20-458334-459126", "20-532878-533632", "20-947449-948352", "20-1045751-1046617", "20-1127074-1128014", "20-1214806-1215670", "20-1401520-1402466", "20-1416306-1417189", "20-1417445-1418133", "20-1423597-1424499", "20-1426900-1427824", "20-1433773-1434635", "20-1461721-1462625", "20-1491896-1492806", "20-1496227-1497052", "20-1498594-1499385", "20-1503570-1504475", "20-1529550-1530453", "20-1584956-1585842", "20-1589229-1590088", "20-1892023-1892958", "20-2149966-2150859", "20-2282694-2283585", "20-2333951-2334727", "20-793174-794058", "20-1439538-1440464", "20-1452364-1453271", "20-1547661-1548536")


### build custom genome using BSgenome
### not run 
# fa_file <- "/Volumes/T7 Shield/Quail_scMultiome/custom_genome/Coturnix_japonica.Coturnix_japonica_2.0.dna.toplevel.regchr.fa"
# twobit_file <- "/Volumes/T7 Shield/Quail_scMultiome/custom_genome/custom_genome.2bit"
# fa <- Biostrings::readDNAStringSet(fa_file)
# rtracklayer::export.2bit(fa, twobit_file)
# forgeBSgenomeDataPkg("/Volumes/T7 Shield/Quail_scMultiome/custom_genome/MyCustomGenome.seed",seqs_srcdir = '/Volumes/T7 Shield/Quail_scMultiome/custom_genome/')
# # R CMD build BSgenome.MySpecies.MyLab.Custom3 # from the 
# # R CMD check BSgenome.MySpecies.MyLab.Custom3_1.0.0.tar.gz
# # R CMD INSTALL  BSgenome.MySpecies.MyLab.Custom3_1.0.0.tar.gz

library(JASPAR2020)
library(TFBSTools)
library(BSgenome.MySpecies.MyLab.Custom3)

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
int <- AddMotifs(object = int, genome = BSgenome.MySpecies.MyLab.Custom3, pfm = pfm)
int <- RegionStats(int, genome = BSgenome.MySpecies.MyLab.Custom3)
metafeatures <- Seurat::GetAssayData(
  object = int[['ATAC']], slot = 'meta.features'
)
metafeatures <- metafeatures[!is.na(metafeatures$GC.percent),]
da_peaks_TF <- do.call(rbind,lapply(unique(da_peaks$group), function(x){
  peaks <- da_peaks$feature[da_peaks$group==x]
  # peaks <- out$peak[out$gene == 'ASIP' & out$Correlation < -0.2]
  idx_match <- match(peaks,rownames(metafeatures))
  idx_match <- idx_match[!is.na(idx_match)]
  query.feature <- metafeatures[idx_match,]
  features.choose <- metafeatures
  bg_peaks <- MatchRegionStats(
    meta.feature = features.choose,
    query.feature = query.feature,
    features.match = c("GC.percent"),n = max(dim(query.feature)[1],40000)
      
  )
  enriched.motifs <- FindMotifs(
    object = int,
    background = bg_peaks,
    features = rownames(query.feature))
  enriched.motifs$group <- x
  return(enriched.motifs)
}))
# enriched.motifs$padj <- p.adjust(enriched.motifs$pvalue, method = "BH")
# enriched.motifs <- enriched.motifs %>% dplyr::filter(padj < 0.05)
da_peaks_TF$padj <- p.adjust(da_peaks_TF$pvalue, method = "BH")
da_peaks_TF_top <- da_peaks_TF %>% dplyr::filter(padj < 0.05) %>% 
  group_by(group) %>% 
  slice_max(order_by = -padj, n=25, with_ties = F)
print(da_peaks_TF_top, n=100)

peak_split <- do.call(rbind, strsplit(peaks, "-", fixed = TRUE))
peak_gr <- GRanges(
  seqnames = peak_split[, 1],
  ranges = IRanges(
    start = as.integer(peak_split[, 2]),
    end = as.integer(peak_split[, 3])
  )
)
motif_ix <- matchMotifs(pfm, subject=peak_gr, 
            genome = BSgenome.MySpecies.MyLab.Custom3,
            out = "scores")
motif_ix@assays@data$motifMatches


#### DirectNet analysis
# devtools::install_github("zhanglhbioinfor/DIRECT-NET")
# devtools::install_github("immunogenomics/presto")
set.seed(1234)
library(DIRECTNET)
library(xgboost)
library(cicero)
library(AnnotationDbi)
library(Seurat)
library(Signac)
library(patchwork)
library(dplyr)
library(ggplot2)
options(stringsAsFactors = FALSE)
dir.create("DirectNet_results") # put direct-net results

#### DIRECT-NET refer https://htmlpreview.github.io/?https://github.com/zhanglhbioinfor/DIRECT-NET/blob/main/tutorial/demo_DIRECTNET_PBMC.html
# preparing genome info
genome.info <- data.frame(
  Chrom = as.character(seqnames(genes)),
  Starts = start(genes),
  Ends = start(genes) + 1,   # mimic your code using TXSEQSTART + 1
  genes = mcols(genes)$gene_name  # or gene_id if you prefer
)
# run DIRECT-Net with modified `my_Run_DIRECT_NET`
source('~/Documents/GitHub/CCI2TF2CCI/R/modeling.R') # modified `Run_DIRECT_NET`
Idents(int) <- 'anno_comb'
table(int@active.ident)
pos_idents <- as.character(unique(int@active.ident))
pos_idents <- pos_idents[grepl("_pos", pos_idents)]
pos_idents <- setdiff(pos_idents, 'E4_M1_ASIP_pos') 
start_time <- Sys.time()
linked_feature_name <- c("20-458334-459126", "20-532878-533632", "20-947449-948352", "20-1045751-1046617", "20-1127074-1128014", "20-1214806-1215670", "20-1401520-1402466", "20-1416306-1417189", "20-1417445-1418133", "20-1423597-1424499", "20-1426900-1427824", "20-1433773-1434635", "20-1461721-1462625", "20-1491896-1492806", "20-1496227-1497052", "20-1498594-1499385", "20-1503570-1504475", "20-1529550-1530453", "20-1584956-1585842", "20-1589229-1590088", "20-1892023-1892958", "20-2149966-2150859", "20-2282694-2283585", "20-2333951-2334727", "20-793174-794058", "20-1439538-1440464", "20-1452364-1453271", "20-1547661-1548536")
DIRECT_NET_Result_list <- lapply(pos_idents, 
function(x) {
int_x <- subset(int, subset=anno_comb==x)
DIRECT_NET_Result_x <- my_Run_DIRECT_NET(int_x, peakcalling = FALSE, k_neigh = min(floor(dim(int_x)[2]/20),30), atacbinary = TRUE,
     max_overlap=0.5, size_factor_normalize = TRUE,
     genome.info = genome.info, nthread = 1, reduction.name = NULL,
     focus_markers = 'ASIP') # unique(deg_list$all$gene)
DIRECT_NET_Result_x$Peak1 <- paste(DIRECT_NET_Result_x$Chr, DIRECT_NET_Result_x$Starts, DIRECT_NET_Result_x$Ends,sep='_')
DIRECT_NET_Result_x$Peak22 <- gsub('_','-',DIRECT_NET_Result_x$Peak2)
DIRECT_NET_Result_x$function_type1 <- DIRECT_NET_Result_x$function_type
DIRECT_NET_Result_x$function_type <- 'HC'
DIRECT_NET_Result_x$overlap <- DIRECT_NET_Result_x$Peak2 %in% gsub('-','_',linked_feature_name)
DIRECT_NET_Result_x})
Sys.time() - start_time
names(DIRECT_NET_Result_list) <- pos_idents
DIRECT_NET_Result_list <- lapply(1:length(DIRECT_NET_Result_list), function(x) {
  DIRECT_NET_Result_list[[x]]$condition <-  names(DIRECT_NET_Result_list)[x]
  return(DIRECT_NET_Result_list[[x]])
})
names(DIRECT_NET_Result_list) <- pos_idents
results_full <- do.call(rbind, DIRECT_NET_Result_list)
results_overlap_only <- do.call(rbind, lapply(DIRECT_NET_Result_list, function(x) {x[x$overlap==T,]}))
write.csv(results_full, file = 'DirectNet_results_ASIP_pos_full.csv')
write.csv(results_overlap_only, file = 'DirectNet_results_ASIP_pos_overlap.csv')

results_full_motif <- lapply(names(DIRECT_NET_Result_list),function(x){
  peaks <- results_full$Peak22[results_full$condition==x]
  # peaks <- out$peak[out$gene == 'ASIP' & out$Correlation < -0.2]
  idx_match <- match(peaks,rownames(metafeatures))
  idx_match <- idx_match[!is.na(idx_match)]
  query.feature <- metafeatures[idx_match,]
  features.choose <- metafeatures
  bg_peaks <- MatchRegionStats(
    meta.feature = features.choose,
    query.feature = query.feature,
    features.match = c("GC.percent"),n = max(dim(query.feature)[1],40000)
    
  )
  enriched.motifs <- FindMotifs(
    object = int,# subset(int, subset=  anno_comb==),
    background = bg_peaks,
    features = rownames(query.feature))
  enriched.motifs$group <- x
  return(enriched.motifs)
})
names(results_full_motif) <- pos_idents
write.csv(do.call(rbind, results_full_motif), file = 'DirectNet_results_ASIP_pos_motif.csv')
motif.obj <- Motifs(int)

lapply(results_full_motif,function(x) {
  p_motif <- MotifPlot(int, motifs =  rownames(x)[1:10])
  # ggsave(filename = paste(x$group[1],'top10_motif.pdf',sep='_'),plot = wrap_plots(p_motif, ncol=4),
  #         width = 20, height=12)
  # return(NULL)
  })

p_vln <- VlnPlot(int, features = results_overlap_only$Peak22,group.by = 'anno_comb',ncol = 5)
ggsave(filename = 'vln_peaks.pdf',plot = p_vln,width = 50,height = 40, units = 'cm')

### coverage plot
# Load GTF
gene_anno <- rtracklayer::readGFF("../Coturnix_japonica.Coturnix_japonica_2.0.113.gtf")
# rename some columns to match requirements
gene_anno$chromosome <- paste0("", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_name
gene_anno$gene[is.na(gene_anno$gene)] <- gene_anno$gene_id[is.na(gene_anno$gene)]
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene
marker <- 'ASIP'

### 1. Arch Link plots
for(j in 1:length(DIRECT_NET_Result_list)){
  DIRECT_NET_Result_tmp <- DIRECT_NET_Result_list[[j]][which(DIRECT_NET_Result_list[[j]]$gene %in% 
                                                               marker), , drop = FALSE]
  conns <- DIRECT_NET_Result_tmp[, 5:7]
  names(conns) <- c("Peak1", "Peak2", "coaccess")
  rownames(conns) <- NULL
  conns$coaccess <- as.numeric(conns$coaccess) 
  pdf(paste('ArchLink_',names(DIRECT_NET_Result_list)[j],'.pdf',sep=''),width = 10,height= 7.5)
  plot_connections(conns, DIRECT_NET_Result_tmp$Chr[1], as.numeric(DIRECT_NET_Result_tmp$Starts[1]) - 
                     250000, as.numeric(DIRECT_NET_Result_tmp$Starts[1]) + 
                     250000, gene_model = gene_anno, coaccess_cutoff = 0, 
                   connection_width = 0.5, connection_color = "#99000d", 
                   peak_color = "#666666", collapseTranscripts = "longest")
  dev.off()
}
### 2. peak activation/inhibition classification
results_full$corr_global <- out$Correlation[match(paste(results_full$gene,results_full$Peak22),
                                                  paste(out$gene,out$peak))] 
out_by_cluster <- do.call(rbind, lapply(1:length(pos_idents),function(j){
  out_j <- myArchR_peak2geneLinkage(subset(int, subset = anno_comb==pos_idents[j]))
  out_j$condition <- pos_idents[j]
  return(out_j)
}))
results_full$corr_in_group <- out_by_cluster$Correlation[match(paste(results_full$gene,results_full$Peak22,results_full$condition),
                                                  paste(out_by_cluster$gene,out_by_cluster$peak,out_by_cluster$condition))] 
da_peaks <- presto:::wilcoxauc.Seurat(int, group_by = "anno_comb", assay = 'data', seurat_assay = "ATAC", verbose= TRUE)
results_full$dap_logFC <- da_peaks$logFC[match(paste(results_full$condition,results_full$Peak22),
                                               paste(da_peaks$group,da_peaks$feature))]
results_full$dap_padj <- da_peaks$padj[match(paste(results_full$condition,results_full$Peak22),
                                               paste(da_peaks$group,da_peaks$feature))]
#results_full <- subset(results_full, select = -correlation)
write.csv(results_full, file = 'DirectNet_results_ASIP_pos_full_corr.csv')

### 2.1 activator, repressor
results_full_motif <- lapply(names(DIRECT_NET_Result_list),function(x){
  peaks <- results_full$Peak22[results_full$condition==x & results_full$corr_in_group>0]
  # peaks <- out$peak[out$gene == 'ASIP' & out$Correlation < -0.2]
  idx_match <- match(peaks,rownames(metafeatures))
  idx_match <- idx_match[!is.na(idx_match)]
  query.feature <- metafeatures[idx_match,]
  features.choose <- metafeatures
  bg_peaks <- MatchRegionStats(
    meta.feature = features.choose,
    query.feature = query.feature,
    features.match = c("GC.percent"),n = max(dim(query.feature)[1],40000)
    
  )
  enriched.motifs1 <- FindMotifs(
    object = int,# subset(int, subset=  anno_comb==),
    background = bg_peaks,
    features = rownames(query.feature))
  enriched.motifs1$group <- x
  enriched.motifs1$activator_repressor <- 'activator'
  
  peaks <- results_full$Peak22[results_full$condition==x & results_full$corr_in_group<0]
  # peaks <- out$peak[out$gene == 'ASIP' & out$Correlation < -0.2]
  idx_match <- match(peaks,rownames(metafeatures))
  idx_match <- idx_match[!is.na(idx_match)]
  query.feature <- metafeatures[idx_match,]
  features.choose <- metafeatures
  bg_peaks <- MatchRegionStats(
    meta.feature = features.choose,
    query.feature = query.feature,
    features.match = c("GC.percent"),n = max(dim(query.feature)[1],40000)
    
  )
  enriched.motifs2 <- FindMotifs(
    object = int,# subset(int, subset=  anno_comb==),
    background = bg_peaks,
    features = rownames(query.feature))
  enriched.motifs2$group <- x
  enriched.motifs2$activator_repressor <- 'repressor'
  
  enriched.motifs <- rbind(enriched.motifs1,enriched.motifs2)
  return(enriched.motifs)
})
names(results_full_motif) <- pos_idents
results_full_motif_df <- do.call(rbind, results_full_motif)
results_full_motif_df <- results_full_motif_df[results_full_motif_df$pvalue<0.05,]

write.csv(results_full_motif_df, file = 'DirectNet_results_ASIP_pos_motif_activator_repressor.csv')


### 3. plot peaks
fragments1 <-  UpdatePath(
  int@assays$ATAC@fragments[[1]],
  new.path = "../E4/atac_fragments.tsv.gz"
)
fragments8 <- UpdatePath(
  int@assays$ATAC@fragments[[8]],
  new.path = "../E5/atac_fragments.tsv.gz"
)
Fragments(int) <- NULL
Fragments(int)<- list(fragments1, fragments8)

int_subset <- subset(int, subset=anno_comb %in% c(
  "E4_M1_ASIP_pos","E4_M1_ASIP_neg",
  "E5_ML_ASIP_pos","E5_ML_ASIP_neg",
  "E5_M1_ASIP_pos", "E5_M1_ASIP_neg", 
  "E5_M2_ASIP_pos", "E5_M2_ASIP_neg", 
  "E5_M3_ASIP_pos", "E5_M3_ASIP_neg"
))
int_subset$anno_comb <- factor(int_subset$anno_comb, levels = c(
  "E4_M1_ASIP_pos","E4_M1_ASIP_neg",
  "E5_ML_ASIP_pos","E5_ML_ASIP_neg",
  "E5_M1_ASIP_pos", "E5_M1_ASIP_neg", 
  "E5_M2_ASIP_pos", "E5_M2_ASIP_neg", 
  "E5_M3_ASIP_pos", "E5_M3_ASIP_neg"
))
## coverage plot
unique_peaks <- (unique(results_full$Peak22))
for(j in 1:length(unique_peaks)){
  cov_plot_j <- CoveragePlot(
    object = int_subset, group.by = 'anno_comb',
    region = unique_peaks[j],
    annotation = FALSE,
    peaks = FALSE
  )
  ggsave(filename = paste(unique_peaks[j],'.pdf',sep=''), cov_plot_j,width = 10, height = 16, units = 'cm')
}


## heatmap
results_full_raw <- results_full
results_full <- results_full[results_full$corr_in_group>0,]
results_full <- results_full[results_full$corr_in_group<0,]

int_subset <- RegionMatrix(
  int_subset,
  StringToGRanges(results_full$Peak22, sep = c("-", "-")),
  'heatmap',
  group.by = 'anno_comb',
  upstream = 3000,
  downstream = 3000,
  verbose = TRUE
)
p_heatmap <- RegionHeatmap(int_subset, 'heatmap',nrow = 5)
ggsave(filename = 'regionheatmap_DirectNet_peaks.pdf',p_heatmap,width = 10, height = 16, units = 'cm')

x = GetAssayData(object = int_subset, slot = "positionEnrichment")
par(mfrow=c(3,2))
plot(rowSums(x$heatmap$E4_M1_ASIP_neg)/1585,rowSums(x$heatmap$E4_M1_ASIP_pos)/13);abline(a=0,b=1)
plot(rowSums(x$heatmap$E5_ML_ASIP_neg)/1665,rowSums(x$heatmap$E5_ML_ASIP_pos)/138);abline(a=0,b=1)
plot(rowSums(x$heatmap$E5_M1_ASIP_neg)/319,rowSums(x$heatmap$E5_M1_ASIP_pos)/66);abline(a=0,b=1)
plot(rowSums(x$heatmap$E5_M2_ASIP_neg)/1539,rowSums(x$heatmap$E5_M2_ASIP_pos)/68);abline(a=0,b=1)
plot(rowSums(x$heatmap$E5_M3_ASIP_neg)/226,rowSums(x$heatmap$E5_M3_ASIP_pos)/49);abline(a=0,b=1)

int_subset <- RegionMatrix(
  int_subset, StringToGRanges(results_full$Peak22[results_full$condition=='E5_ML_ASIP_pos'], sep = c("-", "-")),
  'heatmap_E5_ML', group.by = 'anno_comb',upstream = 3000, downstream = 3000,verbose = TRUE
)
p_heatmap <- RegionHeatmap(int_subset, 'heatmap_E5_ML',nrow = 5)
ggsave(filename = 'regionheatmap_DirectNet_E5_ML.pdf',p_heatmap,width = 10, height = 16, units = 'cm')

int_subset <- RegionMatrix(
  int_subset, StringToGRanges(results_full$Peak22[results_full$condition=='E5_M1_ASIP_pos'], sep = c("-", "-")),
  'heatmap_E5_M1', group.by = 'anno_comb',upstream = 3000, downstream = 3000,verbose = TRUE
)
p_heatmap <- RegionHeatmap(int_subset, 'heatmap_E5_M1',nrow = 5)
ggsave(filename = 'regionheatmap_DirectNet_E5_M1.pdf',p_heatmap,width = 10, height = 16, units = 'cm')

int_subset <- RegionMatrix(
  int_subset, StringToGRanges(results_full$Peak22[results_full$condition=='E5_M2_ASIP_pos'], sep = c("-", "-")),
  'heatmap_E5_M2', group.by = 'anno_comb',upstream = 3000, downstream = 3000,verbose = TRUE
)
p_heatmap <- RegionHeatmap(int_subset, 'heatmap_E5_M2',nrow = 5)
ggsave(filename = 'regionheatmap_DirectNet_E5_M2.pdf',p_heatmap,width = 10, height = 16, units = 'cm')

int_subset <- RegionMatrix(
  int_subset, StringToGRanges(results_full$Peak22[results_full$condition=='E5_M3_ASIP_pos'], sep = c("-", "-")),
  'heatmap_E5_M3', group.by = 'anno_comb',upstream = 3000, downstream = 3000,verbose = TRUE
)
p_heatmap <- RegionHeatmap(int_subset, 'heatmap_E5_M3',nrow = 5)
ggsave(filename = 'regionheatmap_DirectNet_E5_M3.pdf',p_heatmap,width = 10, height = 16, units = 'cm')

### quantify the difference 
heatmap_mode <- c('heatmap','heatmap_E5_ML','heatmap_E5_M1', 'heatmap_E5_M2','heatmap_E5_M3')
pop_pairs <- t(matrix(c(
  "E4_M1_ASIP_pos","E4_M1_ASIP_neg",
  "E5_ML_ASIP_pos","E5_ML_ASIP_neg",
  "E5_M1_ASIP_pos", "E5_M1_ASIP_neg", 
  "E5_M2_ASIP_pos", "E5_M2_ASIP_neg", 
  "E5_M3_ASIP_pos", "E5_M3_ASIP_neg"
),nrow=2))

pvalue_mtx <- matrix(NA,length(heatmap_mode),dim(pop_pairs)[1])
dimnames(pvalue_mtx) <- list(heatmap_mode,paste(pop_pairs[,1],pop_pairs[,2],sep='_vs_'))
slope_mtx <- pvalue_mtx
for(i in 1:length(heatmap_mode)){
  heatmap_mode_i <- heatmap_mode[i]
  data_i <- int_subset@assays$ATAC@positionEnrichment[[heatmap_mode_i]][c(pop_pairs)]
  cell_i <- int_subset@assays$ATAC@positionEnrichment[[heatmap_mode_i]]$function.parameters$cells[c(pop_pairs)]
  for(j in 1:dim(pop_pairs)[1]){
    z1 <- pop_pairs[j,1]
    z2 <- pop_pairs[j,2]
    x <- rowSums(data_i[[z1]])/cell_i[z1]
    y <- rowSums(data_i[[z2]])/cell_i[z2]
    test_j <- t.test(x, y, paired = TRUE, alternative = "two.sided")
    pvalue_mtx[i,j] <- test_j$p.value
    lm_j <- lm(x~y+0)
    slope_mtx[i,j] <-  lm_j$coefficients
  }
}
write.csv(pvalue_mtx,'regionheatmap_pvalue.csv')
write.csv(slope_mtx,'regionheatmap_slope.csv')

