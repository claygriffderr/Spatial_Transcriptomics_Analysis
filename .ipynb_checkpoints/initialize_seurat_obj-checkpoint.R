#' initialize the Seurat object for future use
#' @author Clay Griffin-Derr

#' @description Save a .rds file to your local directory with SCTransform(), RunPCA(), RunUMAP(), FindNeighbors(), and FindClusters(). Parameters for PCA, UMAP and Neighbors are set to 30, 1:30, 1:30

#' @param home_location directory which directly holds uncompressed Xenium data of interest
#' @param save_location Directory to store output files
#' @param fov_name Assign unique name to each segment/run/slice (function only runs one at a time)
#' @param resolution The resolution for FindClusters()
#' @param output_name File name of the output, make sure to include .rds

#' @return A Seurat object already loaded with the data -- can immediately assign varName <- create_Seurat_Obj()

#' @import Seurat
#' @import SeuratObject



create_Seurat_Obj <- function(home_location, save_location, fov_name, resolution, output_name) {
      
    # -- 1 setup directories
    home <- home_location
    save.location <- save_location
        #Validate directory
    if (!dir.exists(home)) {
        stop("The 'home' directory does not exist. Please check the path.")
      }
    if(!dir.exists(save_location)) {
        stop("The 'output_name' directory does not exist. Please check the path.")
        }
    
    obj <- LoadXenium(data.dir = home, fov = fov_name)
    obj <- subset(obj, subset = nCount_Xenium > 0)
    
    print(obj)
    print(names(obj@assays))


    # -- 2 Normalize using SCTransform
    obj <- SCTransform(obj, assay = "Xenium")
    
    # -- 3 PCA, UMAP, Neighbors, Clusters
    obj <- RunPCA(obj, npcs = 30) #compresses expression data into #principal components (major axes of variation) 
    
    obj <- RunUMAP(obj, dims = 1:30)
    obj <- FindNeighbors(obj, dims = 1:30)

    obj <- FindClusters(obj, resolution = resolution)

    full_path <- file.path(save_location, output_name)
    saveRDS(obj, file = full_path)


    # -- 4 Return the object so it can be immediately assigned
    return(obj)
}