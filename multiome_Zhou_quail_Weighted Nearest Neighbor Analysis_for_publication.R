library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(SeuratWrappers)

library(SoupX)
library(DropletUtils)
library(DoubletFinder)

setwd("/home/miller/Desktop/Zhou_scAnalysis/Quail")
getwd()

#############################################
###soupX ambient RNA removal###
#Define sample directories (modify this list with actual paths)
sample_dirs <- c("E4", "E5")  # Add all your sample folder names

#Loop through each sample
for (sample in sample_dirs) {
  
  # Construct paths to h5 files
  filtered_path <- file.path(sample, "filtered_feature_bc_matrix.h5")
  raw_path <- file.path(sample, "raw_feature_bc_matrix.h5")
  
  # ad in h5 files
  filt.matrix <- Read10X_h5(filtered_path, use.names = TRUE)$`Gene Expression`
  raw.matrix  <- Read10X_h5(raw_path, use.names = TRUE)$`Gene Expression`
  
  # reate Seurat object
  srat <- CreateSeuratObject(counts = filt.matrix)
  print(paste("Processing:", sample))
  
  # reate SoupX object
  soup.channel <- SoupChannel(raw.matrix, filt.matrix)
  
  # Remove unnecessary variables to free memory
  rm(filt.matrix)
  rm(raw.matrix)
  
  # Perform SCTransform, PCA, UMAP, and Clustering
  srat <- SCTransform(srat, verbose = FALSE) %>% 
    RunPCA() %>% 
    RunUMAP(dims = 1:40) %>% 
    FindNeighbors(dims = 1:40) %>% 
    FindClusters(verbose = TRUE)
  
  # Extract metadata and clustering info
  meta <- srat@meta.data
  umap <- srat@reductions$umap@cell.embeddings
  soup.channel <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup.channel <- setDR(soup.channel, umap)
  
  # Remove unnecessary objects
  rm(srat)
  rm(meta)
  
  # stimate ambient RNA contamination
  soup.channel <- autoEstCont(soup.channel)
  
  # Check high background genes (Optional)
  print(head(soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = TRUE), ], n = 20))
  
  # Generate adjusted count matrix
  adj.matrix <- adjustCounts(soup.channel, roundToInt = TRUE)
  
  # Save corrected count matrix
  output_dir <- paste0("soupX_", sample)  # Unique output directory for each sample
  DropletUtils:::write10xCounts(output_dir, adj.matrix)
  
  print(paste("Finished processing:", sample))
}


###create seurat obj###
#https://matthieuxmoreau.github.io/EarlyPallialNeurogenesis/html-Reports/Quality_Control.html#Use_Scrublet_to_detect_obvious_doublets
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(GenomicRanges)
library(gintools)
library(future)

setwd("/home/miller/Desktop/Zhou_scAnalysis/Quail")
getwd()

## Import GFF file as a GRanges object (do not make unique(gr)!).
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)

file <- "./Coturnix_japonica.Coturnix_japonica_2.0.113.gtf"

gr <- rtracklayer::import(file)


#make seqinfo object as genome of CreateChromatinAssay
#https://rdrr.io/bioc/GenomeInfoDb/man/Seqinfo-class.html
#seqlength was calculate by samtools faidx
x <- Seqinfo(seqnames=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","Z","MT"),
             seqlengths=c(175656249,134823679,100941466,82189922,54019043,31416776,33419797,26894524,21326773,18483700,17799862,17221944,15706747,12824963,11697260,344938,8971062,9547476,8760791,12446475,5945000,4126574,5051087,5496777,2511919,4932196,4779900,4019625,67000979,16697),
             isCircular=c(FALSE, FALSE, FALSE, FALSE,FALSE, FALSE,FALSE, FALSE, FALSE, FALSE,FALSE, FALSE,FALSE, FALSE, FALSE, FALSE,FALSE, FALSE,FALSE, FALSE, FALSE, FALSE,FALSE, FALSE,FALSE, FALSE, FALSE, FALSE,FALSE, TRUE),
             genome="Coturnix_japonica_2.0")


# Define sample names (modify based on your dataset)
samples <- c("E4", "E5")  # Replace with actual sample names

# Loop through each sample
for (sample in samples) {
  
  # Construct paths for RNA (SoupX corrected) and ATAC data
  rna_path <- paste0("./soupX_", sample)
  atac_h5_path <- file.path(sample, "filtered_feature_bc_matrix.h5")
  atac_fragments <- file.path(sample, "atac_fragments.tsv.gz")
  
  # Read RNA (SoupX corrected) data
  adj.matrix <- Read10X(rna_path)
  
  # Create Seurat object for RNA
  srt.obj <- CreateSeuratObject(counts = adj.matrix, assay = "RNA", project = paste0(sample, "_RNA"))
  print(paste("Processing RNA for:", sample))
  
  # Read ATAC data from the HDF5 file
  inputdata.10x <- Read10X_h5(atac_h5_path)
  atac_counts <- inputdata.10x$Peaks
  
  # Create Chromatin Assay for ATAC and add it to the Seurat object
  srt.obj[["ATAC"]] <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = x,   # Ensure 'x' is properly defined (e.g., "hg38" or "mm10")
    fragments = atac_fragments,
    annotation = gr  # Ensure 'gr' contains genome annotations (e.g., EnsDb object)
  )
  
  print(paste("Processing ATAC for:", sample))
  
  # Save the processed Seurat object
  saveRDS(srt.obj, file = paste0("Multiome_Zhou_quail_", sample, "_soupX.obj.RDS"))
  
  print(paste("Finished processing:", sample))
}


#####We perform basic QC based on the number of detected molecules for each modality as well as mitochondrial percentage.
#srt.obj[["percent.mt"]] <- PercentageFeatureSet(srt.obj, pattern = "^MT-")

library(Seurat)
library(dplyr)
library(ggplot2)
library(SeuratWrappers)


seurat_list <- list(
  "E4" = readRDS("Multiome_Zhou_quail_E4_soupX.obj.RDS"),
  "E5" = readRDS("Multiome_Zhou_quail_E5_soupX.obj.RDS")
)


for (sample in names(seurat_list)) {
  
  srt.obj <- seurat_list[[sample]]
  print(paste("Processing:", sample))
  
 
  
  # Calculate mitochondrial percentage
  srt.obj[["percent.mt"]] <- PercentageFeatureSet(srt.obj, pattern = '(^ND1$|^MT-|^ND3$|^ND4$|^ND4L$|^ND5$|^ND6$|^CYTB$|^COX3$|^ATP6$|^ATP8$)')
  
  # Calculate ribosomal percentage
  srt.obj[["percent.ribo"]] <- PercentageFeatureSet(srt.obj, pattern = '(^RPL|^RPS|^MRP)')
  
  # Visualize QC metrics
  print(VlnPlot(srt.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), 
                ncol = 2, pt.size = 0))
  
  # Save the modified Seurat object back into the list
  seurat_list[[sample]] <- srt.obj
}


##filter out percent.mt one by one
seurat_list$E4 <- subset(seurat_list$E4, subset = percent.mt < 40)
seurat_list$E5 <- subset(seurat_list$E5, subset = percent.mt < 30)


#### next QC
for (sample in names(seurat_list)) {
  
  srt.obj <- seurat_list[[sample]]
  print(paste("Processing:", sample))
  
  # Compute RNA and ATAC thresholds
  RNA.max <- round(mean(srt.obj$nCount_RNA) + 2 * sd(srt.obj$nCount_RNA), digits = -2)
  RNA.min <- round(mean(srt.obj$nCount_RNA) - 2 * sd(srt.obj$nCount_RNA), digits = -2)
  ATAC.max <- round(mean(srt.obj$nCount_ATAC) + 2 * sd(srt.obj$nCount_ATAC), digits = -2)
  ATAC.min <- round(mean(srt.obj$nCount_ATAC) - 2 * sd(srt.obj$nCount_ATAC), digits = -2)
  
  # Set minimum parameters to 0 if negative
  RNA.min <- ifelse(RNA.min < 0, 0, RNA.min)
  ATAC.min <- ifelse(ATAC.min < 0, 0, ATAC.min)
  
  # Filter cells based on thresholds
  srt.obj_subset <- subset(srt.obj, subset = nCount_RNA < RNA.max & 
                             nCount_RNA > RNA.min &
                             nCount_ATAC < ATAC.max & 
                             nCount_ATAC > ATAC.min )
  
  # Assign new metadata for filtered object
  srt.obj_subset$stage <- sample
  
  # Save the filtered Seurat object
  saveRDS(srt.obj_subset, file = paste0("Multiome_Zhou_quail_",sample,"_filtered.soupX.obj.RDS"))
  
  # Generate and display violin plots
  print(
    VlnPlot(srt.obj_subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), 
            ncol = 3, log = TRUE, pt.size = 0) + NoLegend()
  )
  
  print(paste("Finished processing:", sample))
}


####doublet removal####
library(Seurat)
library(DoubletFinder)
library(dplyr)

# Define a list of loaded Multiome Seurat objects
seurat_dirs <- list(
  "E4" = readRDS("Multiome_Zhou_quail_E4_filtered.soupX.obj.RDS"),
  "E5" = readRDS("Multiome_Zhou_quail_E5_filtered.soupX.obj.RDS")
)


##start here##
# Initialize an empty list to store Seurat objects with doublets removed
doublet_objects <- list()

# Iterate through the list of Seurat objects
for (i in names(seurat_dirs)) {
  #Access the current Seurat object
  seurat_object <- seurat_dirs[[i]]
  
  seurat_object <- NormalizeData(seurat_object) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA(npcs = 10) %>%
    RunUMAP(dims = 1:10) %>%
    FindNeighbors(dims = 1:10, verbose = F) %>%
    FindClusters(verbose = F)
  
  #Calculate expected number of doublets
  n_cells <- ncol(seurat_object)
  homotypic_prop <- modelHomotypic(seurat_object@meta.data$seurat_clusters)  # Adjust if clusters are not present
  n_expected_doublets <- round(0.085 * n_cells)  # Assuming 8.5% doublet rate
  
  #Run DoubletFinder
  seurat_object <- doubletFinder(
    seurat_object, 
    PCs = 1:10, 
    pN = 0.25, 
    pK = 0.26,  # Adjust `pK` based on optimization
    nExp = n_expected_doublets, 
    reuse.pANN = FALSE
  )
  
  gc()
  
  #Identify and remove doublets
  # change name of metadata column with Singlet/Doublet information
  colnames(seurat_object@meta.data)[grepl('DF.classifications.*', colnames(seurat_object@meta.data))] <- "doublet_col"
  seurat_object <- subset(seurat_object, subset = doublet_col == "Singlet")
  
  #Add the filtered object to the list
  doublet_objects[[i]] <- seurat_object
}


#(optional) If you haven’t optimized the pK parameter, add the following steps before running doubletFinder (need runPCA in advance)
sweep_res <- paramSweep(doublet_objects$E4, PCs = 1:10, sct = FALSE)
sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
test <- find.pK(sweep_stats)
View(test)
optimal_pK <- find.pK(sweep_stats)$pK[which.max(find.pK(sweep_stats)$BCmetric)]


####merge & wnn integration (not done yet!!)
# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)


multiome_list <-readRDS("Multiome_Zhou_quail_ALLlist_scrublet.soupX.obj.RDS")

# Ensure consistent chromosome naming for ATAC assay
common.chroms <- Reduce(intersect, lapply(multiome_list, function(x) seqlevels(x[["ATAC"]])))

# Apply chromosome filtering to all objects
multiome_list <- lapply(multiome_list, function(obj) {
  seqlevels(obj[["ATAC"]]) <- common.chroms
  return(obj)
})

multiome_E4<-multiome_list$E4
multiome_E4$stage<- "E4"

multiome_E5<-multiome_list$E5
multiome_E5$stage<- "E5"

###RNA integration
#for merge
merged.obj <- merge(multiome_E4, y=multiome_E5, add.cell.ids = c("E4", "E5"), project = "Zhou_multiome")


# split the dataset into a list of two seurat objects (stim and CTRL)
srt.list <- SplitObject(merged.obj, split.by = "stage")


# normalize and identify variable features for each dataset independently
srt.list <- lapply(X = srt.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = srt.list)
srt.list <- lapply(X = srt.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})


# Find integration anchors using CCA
anchors <- FindIntegrationAnchors(object.list = srt.list, 
                                  anchor.features = features,
                                  k.anchor = 20,
                                  reduction = "cca")

# Integrate data using CCA
integrated_rna <- IntegrateData(anchorset = anchors)


# Scale the integrated data
integrated_rna <- ScaleData(integrated_rna, verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>% 
  RunUMAP(dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.3, verbose = FALSE)

# Visualization
p1 <- DimPlot(integrated_rna, reduction = "umap", group.by = "stage")
p2 <- DimPlot(integrated_rna, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
              repel = TRUE)
p1 + p2

DimPlot(integrated_rna, reduction = "umap", split.by = "stage")

###ATAC integration
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(GenomicRanges)
library(ggplot2)
library(dplyr)

#https://github.com/satijalab/seurat/issues/6537
integrated_rna <- readRDS("Multiome_Zhou_quail_ALLlist_RNAintegrated_scrublet.soupX.obj.RDS")
integrated_rna

srt.list <- SplitObject(integrated_rna, split.by = "stage")

multiome_E4<-srt.list$E4
DefaultAssay(multiome_E4) <- "ATAC"

multiome_E5<-srt.list$E5
DefaultAssay(multiome_E5) <- "ATAC"


# quantify multiome peaks in the scATAC-seq dataset
counts <- FeatureMatrix(
  fragments = Fragments(multiome_E4),
  features = granges(multiome_E5),
  cells = colnames(multiome_E4)
)

# add new assay with multiome peaks
multiome_E4[['ATAC']] <- CreateChromatinAssay(
  counts = counts,
  fragments = Fragments(multiome_E4)
)

# compute LSI
DefaultAssay(multiome_E4) <- "ATAC"
multiome_E4 <- FindTopFeatures(multiome_E4, min.cutoff = 10) %>%
  RunTFIDF() %>%
  RunSVD()

DefaultAssay(multiome_E5) <- "ATAC"
multiome_E5 <- FindTopFeatures(multiome_E5, min.cutoff = 10) %>%
  RunTFIDF() %>%
  RunSVD()

#Next we can merge the multiome and scATAC datasets together and observe that there is a difference between them that appears to be due to the batch (experiment and technology-specific variation).
# first add dataset-identifying metadata
multiome_E4$dataset <- "ATAC_E4"
multiome_E5$dataset <- "ATAC_E5"

# merge
atac.combined <- merge(multiome_E4, multiome_E5)

# process the combined dataset
atac.combined <- FindTopFeatures(atac.combined, min.cutoff = 10)
atac.combined <- RunTFIDF(atac.combined)
atac.combined <- RunSVD(atac.combined)
atac.combined <- RunUMAP(atac.combined, reduction = "lsi", dims = 2:30)

DimPlot(atac.combined, group.by = "dataset")

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(multiome_E4, multiome_E5),
  anchor.features = rownames(multiome_E4),
  reduction = "rlsi",
  k.anchor = 30,
  dims = 2:30
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = atac.combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", reduction.name = "ATAC_UMAP", dims = 2:30)
DimPlot(integrated, group.by = "stage")
integrated


# Combine RNA and ATAC
integreated_multiome <- integrated_rna
integreated_multiome[["ATAC"]]<- atac.combined[['ATAC']]


# Two reductions: integrated_lsi and umap.atac are added
integreated_multiome@reductions<- append(integrated_rna@reductions,integrated@reductions)
integreated_multiome

# WNN
integreated_multiome <- FindMultiModalNeighbors(integreated_multiome, reduction.list = list("pca", "integrated_lsi"), dims.list = list(1:30, 2:30))
integreated_multiome <- RunUMAP(integreated_multiome, 
                                nn.name = "weighted.nn", 
                                reduction.name = "wnn.umap", 
                                reduction.key = "wnnUMAP_",
                                return.model=TRUE,
                                seed.use=1234)



# Print all graph names stored in the Seurat object
print(names(integreated_multiome@graphs))

integreated_multiome <- FindClusters(integreated_multiome, graph.name = "wsnn", algorithm = 3,resolution = seq(0.1,0.9,0.2), verbose = FALSE)

sapply(grep("res",colnames(integreated_multiome@meta.data),value = TRUE),
       function(x) length(unique(integreated_multiome@meta.data[,x])))

# Plot integrated clusters
Idents(integreated_multiome) <- "wsnn_res.0.1"
integreated_multiome

DimPlot(integreated_multiome, reduction = "wnn.umap", split.by = "stage")


#export seurat obj to loupe
library(loupeR)

create_loupe_from_seurat(
  integreated_multiome,
  output_dir = "/home/miller/Desktop/Zhou_scAnalysis/Quail",
  output_name = "Multiome_Zhou_quail_WNNintegrated_scrublet.soupX.obj",
  dedup_clusters = FALSE,
  feature_ids = NULL,
  executable_path = NULL,
  force = FALSE
)

DimPlot(integreated_multiome, reduction = "umap", split.by = "stage")



#######CellChat analysis####
library(Seurat)
library(dplyr)
library(ggplot2)
library(SeuratDisk)
library(future)
library(SeuratData)
library(CellChat)
library(patchwork)

library(NMF)
library(ggalluvial)

#Using Existing Site Packages in Miniconda with reticulate in R
library(reticulate)

#check current conda envs
conda_list(conda = "auto")

use_condaenv("miniconda3", required = TRUE) #chage the env as needed

# Check the current Python version
reticulate::py_config()

CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

setwd("/home/miller/Desktop/Zhou_scAnalysis/Quail")
getwd()

integreated_multiome<-readRDS(file = "Sub_Multiome_Zhou_quail_motif_ASIP_linkage_WNNintegrated_scrublet.soupX.obj.RDS")

Idents(integreated_multiome) <- "wsnn_res.0.1"

head(integreated_multiome@active.ident)

#change cluster ID for cellchat ("0" is not an allowed cluster name)

integreated_multiome <- RenameIdents(object = integreated_multiome, c(`0` = "Me_L",`1` = "Me_A1",
                                                                      `2` = "Me_A2",`3` = "Me_A3",
                                                                      `4` = "Epi",`7` = "Mel",
                                                                      `8` = "Diff_Epi"))
                                     
DimPlot(integreated_multiome, reduction = "umap", split.by = "stage")


#use RNA assay for cellchat!!! avoid using integrated assay
DefaultAssay(integreated_multiome) <- "RNA"

integreated_multiome <- DietSeurat(
  integreated_multiome,
  assays = "RNA")

#add active identity to metadata for Cellchat 
integreated_multiome$annotated_cell <- Idents(integreated_multiome)

#feed integrated seurat obj to make integrated cellchat obj
object.list <- SplitObject(integreated_multiome, split.by = "stage")

##'E4'
object.list$'E4'$samples <- "E4"

cellChat.E4 <- createCellChat(object = object.list$'E4', group.by = "annotated_cell", assay = "RNA")
#Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellChat.E4@DB <- CellChatDB
cellChat.E4 <- subsetData(cellChat.E4) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellChat.E4 <- identifyOverExpressedGenes(cellChat.E4) %>% identifyOverExpressedInteractions() %>%
  computeCommunProb(type="truncatedMean",trim=0.01) %>% filterCommunication(min.cells = 10) %>%
  computeCommunProbPathway() %>% aggregateNet() %>% netAnalysis_computeCentrality()


##'E5'
object.list$'E5'$samples <- "E5"

cellChat.E5 <- createCellChat(object = object.list$'E5', group.by = "annotated_cell", assay = "RNA")
#Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellChat.E5@DB <- CellChatDB
cellChat.E5 <- subsetData(cellChat.E5) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellChat.E5 <- identifyOverExpressedGenes(cellChat.E5) %>% identifyOverExpressedInteractions() %>%
  computeCommunProb(type="truncatedMean",trim=0.01) %>% filterCommunication(min.cells = 10) %>%
  computeCommunProbPathway() %>% aggregateNet() %>% netAnalysis_computeCentrality()


############################################################

#E4+E5
object.list <- list("E4" = cellChat.E4, "E5"=cellChat.E5)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


###start here###
# Users can now export the merged CellChat object and the list of the two separate objects for later use
save(cellchat, file = "cellchat_E4E5merged_humanDB_multiome_quail_Zhou.RData")

ptm = Sys.time()
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#sep the circle plot by treatment
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:
     length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

#Heatmap showing differential number of interactions or interaction strength among different cell populations across two datasets
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

#Circle plot showing differential number of interactions or interaction strength among different cell populations across two datasets
#The differential number of interactions or interaction strength in the cell-cell communication network between two datasets can be visualized using circle plot, where red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
par(mfrow = c(1,2), xpd=T)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

#Identify cell populations with significant changes in sending or receiving signals
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#have to adjust plot window size!!!
patchwork::wrap_plots(plots = gg)

#Identify the signaling "changes" of specific cell populations
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(1, 2), label.size = 3,idents.use = "Me_A1")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(1, 2), label.size = 3,idents.use = "Me_A2")
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(1, 2), label.size = 3,idents.use = "Me_A3")
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(1, 2), label.size = 3,idents.use = "Me_L")
gg5 <- netAnalysis_signalingChanges_scatter(cellchat, comparison = c(1, 2), label.size = 3,idents.use = "Epi")
patchwork::wrap_plots(plots = list(gg1,gg2))
patchwork::wrap_plots(plots = list(gg3,gg4))
patchwork::wrap_plots(plots = list(gg4,gg5))



#Identify signaling groups based on their "functional" similarity
ptm = Sys.time()

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
#have to adjust plot window size!!!
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

#Identify signaling groups based on "structure" similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)

#Compute and visualize the pathway distance in the learned joint manifold
rankSimilarity(cellchat, comparison2 = c(1, 2), type = "functional")

##Identify altered signaling with distinct interaction strength##
#(A) Compare the overall information flow of each signaling pathway or ligand-receptor pair
gg1 <- rankNet(cellchat, comparison = c(1, 2),mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, comparison = c(1, 2),mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
gg1 + gg2

#(B) Compare outgoing (or incoming) signaling patterns associated with each cell population
#outgoing heatmap
library(ComplexHeatmap)

#outgoing
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, union(object.list[[i+1]]@netP$pathways, object.list[[i+1]]@netP$pathways))
ht1 <- netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 12)
ht2 <- netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 12)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
#incoming heatmap
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 12, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 12, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


##Identify the up-gulated and down-regulated signaling ligand-receptor pairs##
#Identify dysfunctional signaling by comparing the "communication probabities"
ptm = Sys.time()

netVisual_bubble(cellchat,  comparison = c(1, 2), angle.x = 45)

gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased in E5 than E4", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased in E5 than E4", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2

####Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
pathways.show <- c("WNT") 

#loop for pathways in all the 3 samples
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

table(object.list$E5@netP$pathways)

#individual sample
netVisual_aggregate(object.list[[1]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[1]))

pathways.show <- c("WNT") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]] , ht_gap = unit(0.5, "cm"))

# Chord diagram
pathways.show <- c("WNT") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

###pairwise comparisons###
##Identify dysfunctional signaling by using "differential expression" analysis##
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "E5"
neg.dataset = "E4"

# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
#> Use the joint cell labels from the merged CellChat object

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = pos.dataset,ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = neg.dataset,ligand.logFC = -0.05, receptor.logFC = NULL)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

#Visualize the identified up-regulated and down-regulated signaling ligand-receptor pairs
#(A) Bubble plot
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

#(C) Wordcloud plot
# visualize the enriched ligands in the first condition
install.packages("wordcloud")
library(wordcloud)
computeEnrichmentScore(net.up, species = 'human', variable.both = TRUE)
computeEnrichmentScore(net.down, species = 'human', variable.both = TRUE)

##Compare the signaling gene expression distribution between different datasets##
cellchat@meta$samples = factor(cellchat@meta$samples, levels = c("E4","E5")) # set factor level
plotGeneExpression(cellchat, signaling = "WNT", split.by = "samples", colors.ggplot = T, type = "violin")


