# scDNS  class
scDNSobjClass <- function(){
  setClassUnion("matrixORdgCMatrix", c("matrix", "dgCMatrix"))
  setClassUnion("matrixORarray", c("matrix", "array"))
  setClass("scDNS", representation(
    # attribute
    counts =  "matrixORdgCMatrix",
    data =  "matrixORdgCMatrix",
    Network = 'data.frame',
    GroupLabel = 'character',
    GeneVariability= 'numeric',
    Div.Parameters = 'list',
    NEA.Parameters = 'list',
    NEAModel = 'list',
    JDensity_A = 'matrixORarray',
    JDensity_B = 'matrixORarray',
    uniCase = 'character',
    Zscore = 'data.frame',
    scZscore = 'matrixORdgCMatrix',
    Other = 'list'
  ))
}
scDNSobjClass()


seurat2scDNSObj <- function(sob, imputedAssay = "MAGIC_RNA", GroupBy = NULL, ...) {
  # Check if GroupBy argument is provided
  if (is.null(GroupBy)) {
    stop("Please provide cell group information via the 'GroupBy' parameter.")
  }

  # Detect Seurat major version (4 or 5)

  # Extract raw counts matrix
  counts_mat <- Seurat::GetAssayData(sob, assay = "RNA", slot = "counts")

  # Extract imputed expression matrix
  if (imputedAssay %in% names(sob@assays)) {
    data_mat <- Seurat::GetAssayData(sob, assay = imputedAssay, slot = "data")
  } else {
    stop(paste0("Assay '", imputedAssay, "' not found in the Seurat object."))
  }

  # Extract group labels from meta.data
  if (!GroupBy %in% colnames(sob@meta.data)) {
    stop(paste0("Column '", GroupBy, "' not found in sob@meta.data."))
  }
  group_label <- as.character(sob@meta.data[[GroupBy]])

  # Construct scDNS object
  scDNSob <- CreatScDNSobject(
    counts = counts_mat,
    data = data_mat,
    Network = scDNSBioNet,
    GroupLabel = group_label,
    uniCase = unique(group_label),
    ...
  )

  return(scDNSob)
}



setMethod("show", "scDNS",
          function(object) {
            cat("Gene expression matrix size (gene × cell):", dim(object@data), "\n")
            cat("The number of edges in the network:", nrow(object@Network), "\n")
            cat("Group information:", object@uniCase, "\n")
          }
)

setMethod("colnames", signature(x = "scDNS"), function(x) {
  if (!is.null(x@data)) {
    return(colnames(x@data))
  } else {
    warning("Data attribute not available.")
    return(NULL)
  }
})
setMethod("rownames", signature(x = "scDNS"), function(x) {

  if (!is.null(x@data)) {
    return(rownames(x@data))
  } else {
    warning("Data attribute not available.")
    return(NULL)
  }
})






