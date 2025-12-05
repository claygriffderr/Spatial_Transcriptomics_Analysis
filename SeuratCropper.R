library(Seurat)
library(data.table)
library(sf)

PolygonCropFOV_MemSafe <- function(fov, polygon, pyxel_size = 0.2125){
  
  # --- 1. Setup & Pre-calculations ---
  polygon <- st_set_crs(polygon, NA)
  bbox <- st_bbox(polygon)
  
  # Correct BBox Scaling (Divide, don't multiply!)
  min_x <- bbox["xmin"] * pyxel_size
  max_x <- bbox["xmax"] * pyxel_size
  min_y <- bbox["ymin"] * pyxel_size
  max_y <- bbox["ymax"] * pyxel_size
  
  # --- 2. Process Centroids ---
  message("Processing centroids...")
  coords <- GetTissueCoordinates(fov, which = "centroids")
  
  # Get the OFFICIAL names from the FOV object
  official_cell_names <- Cells(fov)
  
  # Sanity Check: Ensure lengths match before we assume order is safe
  if(nrow(coords) != length(official_cell_names)){
    # If they don't match, try to find the 'cell' column which might hold barcodes
    if("cell" %in% colnames(coords)){
        official_cell_names <- coords$cell
    } else {
        stop("Mismatch: Coordinate rows (", nrow(coords), ") != Cell names (", length(official_cell_names), ").")
    }
  }

  # Quick Pre-filter: Only keep cells inside the BBox
  keep_mask <- coords[,1] >= min_x & coords[,1] <= max_x &
               coords[,2] >= min_y & coords[,2] <= max_y
  
  # --- [FIX IS HERE] ---
  # Use the mask on the OFFICIAL names, ignoring the integer rownames of 'coords'
  cells_to_check <- official_cell_names[keep_mask]
  coords_subset  <- coords[keep_mask, ]
  
  if (nrow(coords_subset) == 0){
    warning("No cells found in search area. Returning NULL.")
    return(NULL)
  }
  
  # Precise Polygon Check
  coords_scaled <- coords_subset
  coords_scaled[,1] <- coords_scaled[,1] / pyxel_size # Scale px -> microns
  coords_scaled[,2] <- coords_scaled[,2] / pyxel_size
  
  centroids_sf <- st_as_sf(as.data.frame(coords_scaled), coords = c("x", "y"))
  contained <- st_contains(polygon, centroids_sf, sparse = FALSE)[1,]
  cells_to_keep <- cells_to_check[contained]
  
  message("Found ", length(cells_to_keep), " cells.")
  
  if (length(cells_to_keep) == 0){
      warning("Cells were in BBox but none inside the exact Polygon. Returning NULL.")
      return(NULL)
  }
  
  # --- 3. Process Molecules ---
  molecules <- NULL
  if (!is.null(fov$molecule)){
    message("Processing molecules (Chunked Filter)...")
    molecules_dt <- rbindlist(lapply(names(fov$molecule), function(nm){
      x <- fov$molecule[[nm]]
      
      # Try standard access, fallback to function if needed
      raw_coords <- tryCatch(x@coords, error = function(e) GetTissueCoordinates(x))
      
      idx <- raw_coords[,1] >= min_x & raw_coords[,1] <= max_x &
             raw_coords[,2] >= min_y & raw_coords[,2] <= max_y
      
      if (!any(idx)) return(NULL)
      
      subset_coords <- raw_coords[idx, , drop = FALSE]
      dt <- as.data.table(subset_coords)
      dt[, gene := nm]
      return(dt)
    }), use.names = FALSE)
    
    if (nrow(molecules_dt) > 0) {
      molecules_dt[, x_scaled := x / pyxel_size]
      molecules_dt[, y_scaled := y / pyxel_size]
      mols_sf <- st_as_sf(molecules_dt, coords = c("x_scaled", "y_scaled"), crs = NA)
      is_inside <- st_contains(polygon, mols_sf, sparse = FALSE)[1,]
      molecules_dt <- molecules_dt[is_inside]
      molecules_dt[, c("x_scaled", "y_scaled") := NULL]
      
      if (nrow(molecules_dt) > 0){
        molecules <- CreateMolecules(molecules_dt)
      }
    }
    message("Done processing molecules.")
  }
  
  # --- 4. Final Subset ---
  message("Creating final object...")
  
  new_fov <- subset(x = fov, cells = cells_to_keep)
  
  if (!is.null(molecules)){
     new_fov[["molecules"]] <- molecules
  }
  
  return(new_fov)
}


PolygonCropSeurat <- function(object, polygons, pyxel_size = 0.2125){

    message("Start cropping fovs...")
    cropped_fovs <- lapply(names(polygons), function(name){
        fov <- object@images[[name]]
        # CALL THE NEW MEMORY SAFE FUNCTION
        fov <- PolygonCropFOV_MemSafe(fov, polygons[[name]], pyxel_size)
        message("Done: ", name)
        return(fov)
    })
    # ... (Rest of wrapper remains the same)
    names(cropped_fovs) <- names(polygons)
    
    # ... Merge and return logic from previous snippet ...
    untouched_fovs <- object@images[!(names(object@images) %in% names(cropped_fovs))]
    all_fovs <- c(untouched_fovs, cropped_fovs)
    cells_to_keep <- unlist(sapply(all_fovs, Cells))
    new_object <- subset(object, cells = cells_to_keep)
    new_object@images <- all_fovs 
    
    return(new_object)
}