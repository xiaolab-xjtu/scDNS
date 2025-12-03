#' CreatScDNSobject
#'
#' creat scDNS object
#'
#' @param counts a matrix represents gene count expression, with rows being genes and columns being cells
#' @param data a matrix represents gene normalized expression, with rows being genes and columns being cells.
#' @param Network a two-column data frame represents the biological network
#' @param GroupLabel a string vector representing the grouping information of the cells
#' @param k a numerical value represents the kth (def: 10) nearest neighbor point
#' @param n.grid An integer value representing the number of grid (def:60) used to calculate the joint probability density
#' @param n.coarse An integer value representing the number of grid (def:20) after coarse-graining
#' @param loop.size An integer value representing the network size (def:3000) allocated to each compute core
#' @param parallel.sz An integer representing the number of parallel cores (def:1)
#' @param verbose Logical value, whether progression messages should be printed in the terminal.
#' @param exclude.zero Logical value indicating whether to exclude 0 values when calculating the joint probability density (def: FALSE)
#' @param NoiseRemove Logical value, indicating whether to remove noise that may lead to extremely low probability density (def:TRUE)
#' @param CoarseGrain Logical value, whether to coarse-grain the grid(def: TRUE).
#' @param returnCDS A logical value indicating whether to return the conditional probability density
#' @param ds.method Method to calculate joint probability density, optional knn and kde (def:knn)
#' @param divergence Method for calculating divergence, optional jsd (Jensen-Shannon divergence) and kld (Kullback–Leibler divergence)(def:jsd)
#' @param h A parameter used for kernel density estimation, which only takes effect when the density estimation method is kde (def:1).
#' @param n.dropGene Integer, number of cells sampled for model building (def:3000)
#' @param n.randNet Integer, number of cells sampled for model building (def:3000)
#' @param sdBias Numeric value, bias coefficient used to penalize low degree genes(>=1) (def:1.1)
#' @param parllelModel parllelModel=c('foreach','bplapply')[1]
#' @param rb.jsd
#' @param Div_weight
#' @param uniCase unique cases of groups
#' @param noiseSd
#'
#' @return
#' @export
#'
#' @examples
CreatScDNSobject <- function(counts,
                             data = NULL,
                             Network = NULL,
                             GroupLabel = NULL,
                             uniCase = unique(GroupLabel),
                             k= 10,
                             n.grid = 60,
                             NoiseRemove = TRUE,
                             CoarseGrain = TRUE,
                             n.coarse=20,
                             ds.method = c('knn','kde')[1],
                             divergence = c('jsd','kld')[1],
                             noiseSd=0.01,
                             h = 1,
                             rb.jsd=F,
                             Div_weight=sqrt(as.vector(matrix(1:n.coarse,nrow = n.coarse,ncol = n.coarse)*t(matrix(1:n.coarse,nrow = n.coarse,ncol = n.coarse)))),
                             loop.size=3000,
                             parallel.sz = 3,
                             verbose = TRUE,
                             exclude.zero = FALSE,
                             returnCDS = TRUE,
                             n.dropGene = 3000,
                             n.randNet = 20000,
                             sdBias=1.1,
                             parllelModel=c('foreach','bplapply')[1]){
  #01
  if(is.null(Network)){
    data(scDNSBioNet)
    message('No network is provided, we load the internal network.')
    Network <- scDNSBioNet
  }
  if(is.null(data)){
    message('Data is not provided, we will normalize with NormalizeData functions.')
    data <- NormalizeData(counts)
  }
  #02
  Div.Parameters <- list(k= k,
                         n.grid = n.grid,
                         NoiseRemove = NoiseRemove,
                         CoarseGrain = CoarseGrain,
                         n.coarse = n.coarse,
                         ds.method = ds.method,
                         divergence = divergence,
                         h = h,
                         noiseSd=noiseSd,
                         Div_weight = Div_weight,
                         rb.jsd=rb.jsd,
                         loop.size=loop.size,
                         parallel.sz = parallel.sz,
                         verbose = verbose,
                         exclude.zero = exclude.zero,
                         parllelModel=parllelModel)

  NEA.Parameters <- list(do.impute = FALSE,
                         n.dropGene = n.dropGene,
                         n.randNet = n.randNet,
                         sdBias=sdBias)
  #03
  message('get dropout for each gene pairs.')
  Network = getDoubleDropout(Network = Network,counts = counts)
  # smooth OutLier value
  # message('Smooth OutLier value.')
  # data = SmoothOutLier_UP(data)
  #
  GeneVariability <- CalGeneInfAmount(NormalizeData(counts),
                   probs = 0.99,
                   rescale = FALSE,
                   plot = FALSE)
  gc()

  scDNSobject <- new(Class = "scDNS",
                     counts = counts,
                     data = data,
                     Network = Network,
                     GroupLabel = GroupLabel,
                     Zscore = data.frame(),
                     scZscore = matrix(),
                     GeneVariability= GeneVariability,
                     Div.Parameters = Div.Parameters,
                     NEA.Parameters = NEA.Parameters,
                     NEAModel=list(),
                     JDensity_A = matrix(),
                     JDensity_B = matrix(),
                     uniCase = uniCase,
                     Other=list())
  scDNSobject


}



#' Network with divergence scores for each edge
#'
#' Computes Jensen-Shannon or Kullback-Leibler divergence between conditional density distributions of network edges (gene pairs) in two biological contexts.
#'
#' @param scDNSobject scDNS object
#' @param k a numerical value represents the kth (def: 10) nearest neighbor point
#' @param n.grid An integer value representing the number of grid (def:60) used to calculate the joint probability density
#' @param n.coarse An integer value representing the number of grid (def:20) after coarse-graining
#' @param loop.size An integer value representing the network size (def:3000) allocated to each compute core
#' @param parallel.sz An integer representing the number of parallel cores (def:1)
#' @param verbose Logical value, whether progression messages should be printed in the terminal.
#' @param exclude.zero Logical value indicating whether to exclude 0 values when calculating the joint probability density (def: FALSE)
#' @param NoiseRemove Logical value, indicating whether to remove noise that may lead to extremely low probability density (def:TRUE)
#' @param CoarseGrain Logical value, whether to coarse-grain the grid(def: TRUE).
#' @param returnCDS A logical value indicating whether to return the conditional probability density
#' @param ds.method Method to calculate joint probability density, optional knn and kde (def:knn)
#' @param divergence Method for calculating divergence, optional jsd (Jensen-Shannon divergence) and kld (Kullback–Leibler divergence)(def:jsd)
#' @param h A parameter used for kernel density estimation, which only takes effect when the density estimation method is kde (def:1).
#' @param returnCDS Logical value, whether to return conditional density matrix (def: FALSE)
#' @param parllelModel
#'
#' @return
#' scDNSobject
#' @export
#'
#' @examples
scDNS_1_CalDivs <- function(scDNSobject,
                            k = NULL,
                            n.grid = NULL,
                            NoiseRemove = NULL,
                            CoarseGrain = NULL,
                            n.coarse = NULL,
                            ds.method = NULL,
                            divergence = NULL,
                            h = NULL,
                            noiseSd=NULL,
                            loop.size = NULL,
                            parallel.sz = NULL,
                            verbose = NULL,
                            exclude.zero = NULL,
                            returnCDS = NULL,
                            parllelModel = NULL){
  if(!is.null(k)){scDNSobject@Div.Parameters$k=k}
  if(!is.null(n.grid)){scDNSobject@Div.Parameters$n.grid=n.grid}
  if(!is.null(NoiseRemove)){scDNSobject@Div.Parameters$NoiseRemove=NoiseRemove}
  if(!is.null(CoarseGrain)){scDNSobject@Div.Parameters$CoarseGrain=CoarseGrain}
  if(!is.null(n.coarse)){scDNSobject@Div.Parameters$n.coarse=n.coarse}
  if(!is.null(h)){scDNSobject@Div.Parameters$h=h}
  if(!is.null(loop.size)){scDNSobject@Div.Parameters$loop.size=loop.size}
  if(!is.null(parallel.sz)){scDNSobject@Div.Parameters$parallel.sz=parallel.sz}
  if(!is.null(verbose)){scDNSobject@Div.Parameters$verbose=verbose}
  if(!is.null(exclude.zero)){scDNSobject@Div.Parameters$exclude.zero=exclude.zero}
  if(!is.null(returnCDS)){scDNSobject@Div.Parameters$returnCDS=returnCDS}
  if(!is.null(parllelModel)){scDNSobject@Div.Parameters$parllelModel=parllelModel}

  # Div.Parameters <- list(k= k,
  #                        n.grid = n.grid,
  #                        NoiseRemove = NoiseRemove,
  #                        CoarseGrain = CoarseGrain,
  #                        n.coarse = n.coarse,
  #                        ds.method = ds.method,
  #                        divergence = divergence,
  #                        h = h,
  #                        loop.size=loop.size,
  #                        parallel.sz = parallel.sz,
  #                        verbose = verbose,
  #                        exclude.zero = exclude.zero,
  #                        parllelModel = parllelModel)

  # scDNSobject@Div.Parameters <- Div.Parameters

  NetDivs <- getKLD_cKLDnetwork(scDNSobject@data,
                     Network = scDNSobject@Network,
                     GroupLabel = scDNSobject@GroupLabel,
                     k = scDNSobject@Div.Parameters$k,
                     n.grid = scDNSobject@Div.Parameters$n.grid,
                     NoiseRemove = scDNSobject@Div.Parameters$NoiseRemove,
                     CoarseGrain = scDNSobject@Div.Parameters$CoarseGrain,
                     n.coarse = scDNSobject@Div.Parameters$n.coarse,
                     ds.method = scDNSobject@Div.Parameters$ds.method,
                     divergence = scDNSobject@Div.Parameters$divergence,
                     h = scDNSobject@Div.Parameters$h,
                     noiseSd = scDNSobject@Div.Parameters$noiseSd,
                     loop.size=scDNSobject@Div.Parameters$loop.size,
                     parallel.sz = scDNSobject@Div.Parameters$parallel.sz,
                     verbose = scDNSobject@Div.Parameters$verbose,
                     exclude.zero = scDNSobject@Div.Parameters$exclude.zero,
                     returnCDS = FALSE,
                     Div_weight = scDNSobject@Div.Parameters$Div_weight,
                     parllelModel = scDNSobject@Div.Parameters$parllelModel)
  scDNSobject@Network <- NetDivs$Network
  scDNSobject@JDensity_A <- NetDivs$ContextA_DS
  scDNSobject@JDensity_B <- NetDivs$ContextB_DS
  scDNSobject@uniCase <- NetDivs$uniCase
  rm(NetDivs)
  invisible(gc())
  scDNSobject

}

#' scDNS_2_creatNEAModel
#'
#' Creates a null model for network enrichment analysis by random sampling of genes and networks.
#'
#' @param scDNSobject scDNSobject
#' @param n.dropGene Integer, number of cells sampled for model building (def:3000)
#' @param n.randNet Integer, number of cells sampled for model building (def:20000)
#' @param sdBias Numeric value, bias coefficient used to penalize low degree genes(>=1) (def:1.1)
#'
#' @return
#' scDNSobject
#' @export
#'
#' @examples
scDNS_2_creatNEAModel <- function(scDNSobject,
                                  n.dropGene = NULL,
                                  n.randNet = NULL,
                                  sdBias=NULL){

  if(!is.null(n.dropGene)){scDNSobject@NEA.Parameters$n.dropGene=n.dropGene}
  if(!is.null(n.randNet)){scDNSobject@NEA.Parameters$n.randNet=n.randNet}
  if(!is.null(sdBias)){scDNSobject@NEA.Parameters$sdBias=sdBias}

  NEAModel <- creatNEAModel(counts=scDNSobject@counts,
                            ExpData=scDNSobject@data,
                            do.impute = scDNSobject@NEA.Parameters$do.impute,
                            n.dropGene = scDNSobject@NEA.Parameters$n.dropGene,
                            n.randNet = scDNSobject@NEA.Parameters$n.randNet,
                            k = scDNSobject@Div.Parameters$k,
                            Likelihood=scDNSobject@GeneVariability,
                            GroupLabel=scDNSobject@GroupLabel,
                            n.grid = scDNSobject@Div.Parameters$n.grid,
                            n.coarse=scDNSobject@Div.Parameters$n.coarse,
                            noiseSd = scDNSobject@Div.Parameters$noiseSd,
                            Div_weight = scDNSobject@Div.Parameters$Div_weight,
                            CoarseGrain=scDNSobject@Div.Parameters$CoarseGrain,
                            loop.size=scDNSobject@Div.Parameters$loop.size,
                            parallel.sz = scDNSobject@Div.Parameters$parallel.sz,
                            verbose = scDNSobject@Div.Parameters$verbose,
                            exclude.zero = scDNSobject@Div.Parameters$exclude.zero,
                            NoiseRemove = scDNSobject@Div.Parameters$NoiseRemove,
                            divergence = scDNSobject@Div.Parameters$divergence,
                            ds.method = scDNSobject@Div.Parameters$ds.method,
                            h=scDNSobject@Div.Parameters$h,
                            sdBias=scDNSobject@NEA.Parameters$sdBias,
                            rb.jsd = scDNSobject@Div.Parameters$rb.jsd,
                            parllelModel =  scDNSobject@Div.Parameters$parllelModel)
  scDNSobject@NEAModel <- NEAModel
  scDNSobject
}


#' scDNS_3_GeneZscore
#'
#' @param scDNSobject scDNSobject
#'
#' @return
#' scDNSobject
#' @export
#'
#' @examples
scDNS_3_GeneZscore <- function(scDNSobject){
  scDNSobject <- getZscore(EdgeScore = scDNSobject,
                            NEAModel = scDNSobject@NEAModel,
                            Likelihood = scDNSobject@GeneVariability )
  scDNSobject
}

#' scDNS_4_scContribution
#'
#' @param scDNSobject scDNSobject
#' @param topGene Integer value, representing the number of top significant genes(def:100).
#' When sigGene is not NULL, topGene is invalid.
#' @param sigGene A vector of strings representing the set of genes of interest (def:NULL)
#' @param q.th Adjusted p-value threshold
#'
#' @return
#' scDNSobject
#' @export
#'
#' @examples
scDNS_4_scContribution <- function(scDNSobject,
                                   topGene=NULL,
                                   Pval_col = 'p_comb',
                                   sigGene=NULL,
                                   q.th=0.01,...){
  if(is.null(topGene)){
    if(is.null(sigGene)){
      Zscore <- scDNSobject@Zscore
      p_adj <- p.adjust(Zscore[,Pval_col],method = 'BH')
      sigGene <- Zscore$Gene[p_adj<q.th]
      if(length(sigGene)==0){
        message('The number of sigGene is 0. We used top genes (100)')
      }else{
        message('The number of sigGene is ',length(sigGene),'.')
      }
    }
  }



  nx <-  sqrt(ncol(scDNSobject@JDensity_A))
  # scCon_res <- scContribution(EdgeScore = scDNSobject@Network,
  #                    Zscores = scDNSobject@Zscore,
  #                    ExpData = scDNSobject@data,
  #                    nx = nx,
  #                    rmZero = scDNSobject@Div.Parameters$exclude.zero,
  #                    topGene = topGene,
  #                    sigGene = sigGene,
  #                    NEAModel = scDNSobject@NEAModel)
  # scCon_res <- scContribution(scDNSobject,
  #                          nx=nx,
  #                          rmZero = scDNSobject@Div.Parameters$exclude.zero,
  #                          topGene=topGene,
  #                          sigGene=sigGene)
  # scCon_res <- scContribution_v2(scDNSobject,
  #                                nx=nx,
  #                                rmZero = scDNSobject@Div.Parameters$exclude.zero,
  #                                topGene=topGene,
  #                                sigGene=sigGene,
  #                                rb = scDNSobject@Div.Parameters$rb.jsd)
  # scCon_res <- scContribution_v3(scDNSobject,
  #                             nx=nx,
  #                             # rmZero = scDNSobject@Div.Parameters$exclude.zero,
  #                             topGene=topGene,
  #                             sigGene=sigGene,
  #                             rb = scDNSobject@Div.Parameters$rb.jsd,...)
  scCon_res <- scContribution_v5(scDNSobject,
                                 nx=nx,
                                 # rmZero = scDNSobject@Div.Parameters$exclude.zero,
                                 topGene=topGene,
                                 sigGene=sigGene,
                                 rb = scDNSobject@Div.Parameters$rb.jsd,...)

  scDNSobject@scZscore <- scCon_res
  scDNSobject
}


#' scDNS_5_cluster
#'
#' @param scDNSobj scDNS object
#' @param Sobj seurat object
#' @param biasToscDNS The weight of scZscore during clustering, a floating point number (default 1).
#' When equal to 1, scZscore and PCA's principal component 1 have the same L2 norm.
#' @param resolution The resolution parameter in Seurat's FindClusters function controls the number of clusters.
#' A higher resolution will lead to more clusters, while a lower resolution will lead to fewer clusters.
#' @param merge.red Single-cell transcriptomics dimensionality reduction object (default:pca) used for integration with scZScore.
#' @param red.NewName New scDNS dimensionality reduction name (default: scDNS)
#'
#' @return
#' seurat object
#' @export
#'
#' @examples
scDNS_5_cluster <- function(scDNSobj,
                          Sobj,biasToscDNS=1,resolution=0.5,merge.red='pca',red.NewName='scDNS'){
  scCC <- Zscore <- scDNSobj@scZscore
  scCC_assay <- SeuratObject::CreateAssayObject(data = scCC)
  Sobj[['scCC']] <- scCC_assay
  scDNSLoading <- matrixStats::colSums2(scCC)
  names(scDNSLoading) <- colnames(scCC)
  Sobj$scDNSLoading <- scDNSLoading[colnames(Sobj)]
  Sobj <- scDNS_Cluster(Sobj,
                        scCC = scCC,
                        reduction = merge.red,
                        red.name = red.NewName,
                        biasToscDNS = biasToscDNS,
                        resolution = resolution)
  Sobj
}





