# Load libraries
library(Seurat)
library(tidyverse)  
library(ggplot2)
library(CellChat)
library(future)

# Load object
obj <- readRDS('/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/cellchat/obj_cellchat.rds')
obj
obj@meta.data |> head() 
obj$myeloid_subtype |> unique()

# Change assay
DefaultAssay(obj) <- 'SCT'

# # Join Layers (if needed)
# obj <- JoinLayers(obj)
# # backup <- JoinLayers(backup)
# obj

# Prepare for cell chat 
data.input <- obj[["SCT"]]$data
labels <- obj$myeloid_subtype
meta <- data.frame(labels = labels, row.names = names(labels))

# Create cellchat obj
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

# cellchat <- createCellChat(object = obj, group.by = "myeloid_subtype", assay = "SCT")

# Call database 
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

# Use a subset of CellChatDB 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # Use secreted signaling
showDatabaseCategory(CellChatDB.use)
# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB)

# Set the used database in the object 
cellchat@DB <- CellChatDB.use

# Subset expression data of signaling genes
cellchat <- subsetData(cellchat)

library(future)
plan("multisession", workers = 4) # multicore, sequential
options(future.globals.maxSize = 450 * 1024^3)  
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Time estimation
# execution.time = Sys.time() - ptm
# print(as.numeric(execution.time, units = "secs"))

# Compute the communication probability
# ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean")

# Compute infer cellular communication network
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Extract the inferred cellular communication network
df.net <- subsetCommunication(cellchat, slot.name = "netP")
df.net$pathway_name |> unique()

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate the cell-cell communication network
cellchat <- aggregateNet(cellchat)
# execution.time = Sys.time() - ptm
# print(as.numeric(execution.time, units = "secs"))
saveRDS(cellchat, "/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/cellchat/obj_cellchat_done_v2.rds")