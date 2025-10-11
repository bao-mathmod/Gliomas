# SCRIPT: 01_load_data.R
# PURPOSE: Load all ATAC-seq samples based on the metadata.csv file.
# VERSION: 2.0 - Handles PRJNA941288 RData objects correctly.

# --- Section 1: Setup ---
library(Seurat)
library(Signac)
library(tidyverse)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(biovizBase)

cat("--- Step 1: Setup Complete ---\n")

# --- Section 2: Load Metadata and Annotations ---
metadata_path <- "/mnt/18T/chibao/gliomas/data/upstream/scATAC/metadata.csv"
metadata <- read_csv(metadata_path)

# Extract gene annotations from EnsDb
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

print(head(metadata))
cat("\n--- Step 2: Metadata Loaded ---\n")

# --- Section 3: Iterative Object Creation ---
seurat_objects <- list()

for (i in 1:nrow(metadata)) { # start from 21 to skip already processed samples
    
    sample_info <- metadata[i, ]
    sample_id <- sample_info$sample_id
    project_id <- sample_info$project_id
    data_type <- sample_info$data_type
    
    cat(paste0("\nProcessing sample ", i, "/", nrow(metadata), ": ", sample_id, " (", data_type, ")\n"))
    
    seurat_obj <- NULL
    
    if (data_type == 'multiome_fragments') {
        
        frag_path <- sample_info$fragments_path
        
        # --- MODIFICATION START ---
        
        # 1. Create a Fragment object from the file path first.
        cat("   -> Creating Fragment object...\n")
        frag_obj <- CreateFragmentObject(path = frag_path)
        
        # 2. Create the gene activity matrix, passing the Fragment object *inside a list*.
        cat("   -> Calculating gene activity matrix from fragments...\n")
        gene_activity_counts <- FeatureMatrix(
          fragments = list(frag_obj),  # Note the list() wrapper
          features = annotation,
          cells = NULL
        )
        
        # 3. Create the ChromatinAssay. It's best practice to pass the Fragment object here too.
        chrom_assay <- CreateChromatinAssay(
            counts = gene_activity_counts,
            fragments = frag_obj, # Use the object instead of the path
            genome = "hg38",
            annotation = annotation
        )
        # --- MODIFICATION END ---
        
        seurat_obj <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC")
        
    } else if (data_type == 'snatac_matrix') {

        matrix_path <- sample_info$matrix_path
        peaks_path <- sample_info$peaks_path
        barcodes_path <- sample_info$barcodes_path
        
        # --- NEW LOGIC (v5) ---
        # Add a special handler for the lifted-over samples from PRJNA578617
        
        if (project_id == 'PRJNA578617') {
            cat("   -> Detected Lifted-Over sample. Performing matrix correction...\n")
            
            # --- This block will fix the dimension mismatch for lifted-over data ---
            
            # Infer paths to the original hg19 peaks and the unmapped log
            base_name <- gsub("_peaks.bed.gz", "", basename(peaks_path))
            original_dir <- file.path("/mnt/18T/chibao/gliomas/data/output_cell/snATAC/PRJNA578617")
            peaks_hg19_path <- file.path(original_dir, paste0(base_name, "_peaks.bed.gz"))
            unmapped_log_path <- file.path(dirname(peaks_path), paste0(base_name, "_unmapped.log"))
            
            # Read all necessary files
            mat_triplet <- Matrix::readMM(file = matrix_path)
            barcodes <- read.table(barcodes_path, col.names = c("barcode"))
            peaks_hg38 <- read.table(peaks_path, col.names=c("chr","start","end"), comment.char="#")
            peaks_hg19 <- read.table(peaks_hg19_path, col.names=c("chr","start","end"))
            unmapped_peaks <- read.table(unmapped_log_path, col.names=c("chr","start","end"))
            
            # Create unique names for original and unmapped peaks to find which rows to discard
            peaks_hg19$name <- paste(peaks_hg19$chr, peaks_hg19$start, peaks_hg19$end, sep="-")
            unmapped_peaks$name <- paste(unmapped_peaks$chr, unmapped_peaks$start, unmapped_peaks$end, sep="-")
            
            # Get the numerical indices of the rows to discard from the original matrix
            rows_to_discard <- which(peaks_hg19$name %in% unmapped_peaks$name)
            
            # Filter the matrix by keeping rows that are NOT in the discard list
            mat_summary <- summary(mat_triplet)
            rows_to_keep <- mat_summary$i[!mat_summary$i %in% rows_to_discard]
            
            # Since the matrix is now filtered, we need to create a new one. 
            # This is complex, so for simplicity we will convert to a dense matrix, filter, then go back to sparse.
            # This may be memory intensive for very large matrices but is robust.
            mat_dense <- as.matrix(mat_triplet)
            mat_corrected_dense <- mat_dense[-rows_to_discard, ]
            mat_corrected <- as(mat_corrected_dense, "sparseMatrix")
            rm(mat_dense, mat_corrected_dense) # Clean up memory

            # --- Now create the Seurat object with the corrected matrix ---
            feature.names <- paste(peaks_hg38$chr, peaks_hg38$start, peaks_hg38$end, sep = "-")
            rownames(mat_corrected) <- feature.names
            colnames(mat_corrected) <- barcodes$barcode
            
            chrom_assay <- CreateChromatinAssay(
              counts = mat_corrected,
              genome = 'hg38',
              annotation = annotation
            )

        } else {
            # --- This is the robust code for all other standard matrix files ---
            
            mat <- Matrix::readMM(file = matrix_path)
            barcodes <- read.table(file = barcodes_path, col.names = c("barcode"))
            
            first_line_peaks <- readLines(peaks_path, n = 1)
            
            if (grepl(":", first_line_peaks) && grepl("-", first_line_peaks)) {
                cat("   -> Detected 1-column peak format. Parsing...\n")
                peaks_df <- read.table(file = peaks_path, col.names = c("peak"))
                feature.names <- peaks_df$peak
            } else {
                cat("   -> Detected 3-column BED format. Parsing...\n")
                peaks_df <- read.table(file = peaks_path, col.names = c("chr", "start", "end"), comment.char = "#")
                feature.names <- paste(peaks_df$chr, peaks_df$start, peaks_df$end, sep = "-")
            }
            
            rownames(mat) <- feature.names
            colnames(mat) <- barcodes[, 1]

            chrom_assay <- CreateChromatinAssay(
                counts = mat,
                genome = 'hg38',
                annotation = annotation
            )
        }
        
        seurat_obj <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC")
        
    } else if (data_type == 'snatac_rdata') {
        cat("   -> Special RData case. Loading gene activity and chromVAR data.\n")
        
        # Load gene activities
        load(sample_info$rdata_path) # Assumes 'rdata_path' points to the _rna file
        
        # Load chromVAR data (assuming file name can be inferred)
        chromvar_path <- gsub("_rna.Rdata", "_chromvar.Rdata", sample_info$rdata_path)
        chromvar_path <- gsub("/5.._rna/", paste0("/", sample_info$sample_id, "_chromvar/"), chromvar_path)
        load(chromvar_path)
        
        # Create a Seurat object. Note we are creating an "RNA" assay, not "ATAC".
        seurat_obj <- CreateSeuratObject(counts = gene.activities, assay = "RNA")
        
        # Add the chromVAR data as a new Assay
        seurat_obj[["chromvar"]] <- CreateAssayObject(counts = chromvar_obj)
        
        # Clean up loaded objects from the environment before the next loop
        rm(gene.activities, chromvar_obj)
    }
    
    # Add metadata to the Seurat object
    seurat_obj$project_id <- project_id
    seurat_obj$sample_id <- sample_id
    
    # Add the created object to our list
    seurat_objects[[sample_id]] <- seurat_obj
}

cat("\n--- Step 3: All samples loaded ---\n")

# --- Section 4: Save Results ---
output_file <- "01_seurat_object_list.rds"
saveRDS(seurat_objects, file = output_file)

cat(paste0("\nProcessing complete. Object list saved to: ", output_file, "\n"))
print(names(seurat_objects))