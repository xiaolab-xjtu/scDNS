#' getKLD_cKLDnetwork
#'
#' @param ExpData a matrix represents gene expression, with rows being genes and columns being cells.
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
#' @param parllelModel parllelModel c('foreach','bplapply'). def:foreach
#' @param Div_weight sqrt(as.vector(matrix(1:n.coarse,nrow = n.coarse,ncol = n.coarse)*t(matrix(1:n.coarse,nrow = n.coarse,ncol = n.coarse)))),
#' @param rb.jsd (def:FALSE)
#'
#' @return
#'
#' @keywords internal
#'
#' @examples
#' @keywords internal
getKLD_cKLDnetwork<-function(ExpData,
                             Network = NULL,
                             GroupLabel = NULL,
                             k = 10,
                             n.grid = 60,
                             NoiseRemove = TRUE,
                             CoarseGrain = TRUE,
                             n.coarse=20,
                             ds.method = c('knn','kde'),
                             divergence = c('jsd','kld'),
                             h = 1,
                             noiseSd=0.01,
                             loop.size=3000,
                             parallel.sz = 10,
                             verbose = TRUE,
                             exclude.zero = FALSE,
                             Div_weight = sqrt(as.vector(matrix(1:n.coarse,nrow = n.coarse,ncol = n.coarse)*t(matrix(1:n.coarse,nrow = n.coarse,ncol = n.coarse)))),
                             returnCDS = TRUE,
                             parllelModel=c('foreach','bplapply')[1],
                             rb.jsd=FALSE){
  sqrtForNeg <- function(m){
    sign(m)*sqrt(abs(m))
  }
  # information about grouping
  if(!is.null(GroupLabel)){
    uniCase <- unique(GroupLabel)%>%sort()
    GA_id <-  GroupLabel==uniCase[1] # groups A
    GB_id <-  GroupLabel==uniCase[2] # group B
  }else{
    stop('Parameter GroupLabel is required.')
  }
  # prepare biological network
  geneSds <- MatrixGenerics::rowSds(ExpData)
  # print(dim(ExpData))
  ExpData = ExpData[!(geneSds==0|is.na(geneSds)),]
  # print(dim(ExpData))
  # Filtering network.
  Network <- suppressMessages(preFilterNet(Network,ExpData))
  ExpData <- ExpData[rownames(ExpData)%chin%unlist(Network[,1:2]),]

  RawExp_noise <- ExpData + rnorm(length(ExpData), mean = 0, sd = noiseSd)

  # map expression values to between 1 and n.grid
  # ExpDataDT <- mappingMinMaxRow(ExpData,minD = 1,maxD = n.grid,rmZero = exclude.zero)
  ExpDataDT <- (RawExp_noise - rowMins(ExpData))/(rowMaxs(ExpData) - rowMins(ExpData)) * (n.grid-1) + 1 # V2
  ExpDataDT <- sqrtForNeg(ExpDataDT) # The expression values were squared
  geneSds <- MatrixGenerics::rowSds(ExpDataDT)
  # print(dim(ExpData))
  ExpDataDT <- ExpDataDT[!(geneSds==0|is.na(geneSds)),]
  # print(dim(ExpData))
  # message('Bio net:')
  # print(dim(Network))
  Network <- invisible(preFilterNet(Network,ExpData))
  # GridData <- matlab::repmat(1:n.grid,1,nrow(ExpData))
  zeorData <- ExpDataDT==0
  if(exclude.zero){
    ExpDataDT[as.matrix(zeorData)] <-  10^5 # to avoid bias to zero
  }
  ExpDataDT = as.data.table(t(ExpDataDT))
  if(nrow(Network)<(parallel.sz*loop.size)&nrow(Network)>loop.size){
    loop.size = nrow(Network)/parallel.sz
  }
  ## @@@ context A----
  message('context A')
  ContextA_R = get_DS_cDS_MI_DREMI(ExpDataDT=ExpDataDT[GA_id,], #
                                   Network = Network,
                                   k = k,
                                   h= h,
                                   parallel.sz = parallel.sz,
                                   n.grid = n.grid,
                                   n.coarse = n.coarse,
                                   loop.size = loop.size,
                                   verbose = verbose,
                                   label=uniCase[1],
                                   parllelModel=parllelModel,
                                   CoarseGrain = CoarseGrain,
                                   kde.method = ds.method[1])
  ## @@@ context B----
  message('context B')
  ContextB_R = get_DS_cDS_MI_DREMI(ExpDataDT=ExpDataDT[GB_id,],
                                   Network = Network,
                                   k = k,
                                   h= h,
                                   parallel.sz = parallel.sz,
                                   n.grid = n.grid,
                                   n.coarse = n.coarse,
                                   loop.size = loop.size,
                                   verbose = verbose,
                                   label=uniCase[2],
                                   parllelModel=parllelModel,
                                   CoarseGrain = CoarseGrain,
                                   kde.method = ds.method[1])
  # DensityShowPoint(ExpDataDT[GB_id,]$SLC3A2,ExpDataDT[GB_id,]$`MT-CYB`)|DensityShowPoint(ExpDataDT[GA_id,]$SLC3A2,ExpDataDT[GA_id,]$`MT-CYB`)
  NetNcol = ncol(ContextA_R$Network)
  Network_MI_DREMI = cbind(ContextA_R$Network[,(NetNcol-2):NetNcol],ContextB_R$Network[,(NetNcol-2):NetNcol])
  MI_types = colnames(Network_MI_DREMI)
  Network_MI_DREMI = Network_MI_DREMI[,c(MI_types[stringr::str_detect(MI_types,'^MI')],sort(MI_types[stringr::str_detect(MI_types,'^DREMI')]))]
  Network = cbind(Network,Network_MI_DREMI)
  #
  if(NoiseRemove&ds.method[1]=='knn'){ # Whether to remove some rasters with extremely low density
    if(CoarseGrain){
      message('CoarseGrain')
      NoNoiseCDS = removeNoiseInCDS_towDS(DS1 = ContextA_R$Density,
                                          DS2 = ContextB_R$Density,
                                          cDSList_1 = ContextA_R$cDS,
                                          cDSList_2 = ContextB_R$cDS,
                                          nx = n.coarse,
                                          min.pt = 0.5,
                                          min.density =k/((k)^2*pi)*(n.grid^2/n.coarse^2))
      ContextA_R$cDS=NoNoiseCDS$cDSList_1
      ContextB_R$cDS=NoNoiseCDS$cDSList_2
      # calculate DREMI
      Network[,paste('DREMI_1',uniCase[1],sep = '_')] = MI_from_ked2d_v(x = ContextA_R$cDS$cDS_1,logbase = 2,nx = n.coarse)
      Network[,paste('DREMI_2',uniCase[1],sep = '_')] = MI_from_ked2d_v(x = ContextA_R$cDS$cDS_2,logbase = 2,nx = n.coarse)
      Network[,paste('DREMI_1',uniCase[2],sep = '_')] = MI_from_ked2d_v(x = ContextB_R$cDS$cDS_1,logbase = 2,nx = n.coarse)
      Network[,paste('DREMI_2',uniCase[2],sep = '_')] = MI_from_ked2d_v(x = ContextB_R$cDS$cDS_2,logbase = 2,nx = n.coarse)
    }else{
      NoNoiseCDS = removeNoiseInCDS_towDS(DS1 = ContextA_R$Density,
                                          DS2 = ContextB_R$Density,
                                          cDSList_1 = ContextA_R$cDS,
                                          cDSList_2 = ContextB_R$cDS,
                                          nx = n.grid,
                                          min.pt = 0.5,
                                          min.density =k/((k)^2*pi))
      ContextA_R$cDS=NoNoiseCDS$cDSList_1
      ContextB_R$cDS=NoNoiseCDS$cDSList_2
      # calculate DREMI
      Network[,paste('DREMI_1',uniCase[1],sep = '_')] = MI_from_ked2d_v(x = ContextA_R$cDS$cDS_1,logbase = 2,nx = n.grid)
      Network[,paste('DREMI_2',uniCase[1],sep = '_')] = MI_from_ked2d_v(x = ContextA_R$cDS$cDS_2,logbase = 2,nx = n.grid)
      Network[,paste('DREMI_1',uniCase[2],sep = '_')] = MI_from_ked2d_v(x = ContextB_R$cDS$cDS_1,logbase = 2,nx = n.grid)
      Network[,paste('DREMI_2',uniCase[2],sep = '_')] = MI_from_ked2d_v(x = ContextB_R$cDS$cDS_2,logbase = 2,nx = n.grid)
    }


  }
  invisible(gc())
  # @@@ calculate KLD, cKLD----
  # Div_weight = sqrt(as.vector(matrix(1:20,nrow = 20,ncol = 20)*t(matrix(1:20,nrow = 20,ncol = 20))))
  message('calculate div')
  if(divergence[1]=='kld'){
    Network[,'Div'] = KLD_batch_matrix(Ctl = ContextA_R$cDS$DS,
                                       Pert = ContextB_R$cDS$DS,
                                       symmetrical = TRUE)

    Network[,'cDiv_D1'] = KLD_batch_matrix(Ctl = ContextA_R$cDS$cDS_1,
                                           Pert = ContextB_R$cDS$cDS_1,
                                           symmetrical = TRUE)
    Network[,'cDiv_D2'] = KLD_batch_matrix(Ctl = ContextA_R$cDS$cDS_2,
                                           Pert = ContextB_R$cDS$cDS_2,
                                           symmetrical = TRUE)
  }else{
    if(rb.jsd){
      Network[,'Div'] = JSD_batch_matrix_rawMt_robust(rawDensity1  = ContextA_R$cDS$DS,
                                                      rawDensity2 = ContextB_R$cDS$DS)

      Network[,'cDiv_D1'] = JSD_batch_matrix_rawMt_robust(rawDensity1 = ContextA_R$cDS$cDS_1,
                                                          rawDensity2 = ContextB_R$cDS$cDS_1)
      Network[,'cDiv_D2'] = JSD_batch_matrix_rawMt_robust(rawDensity1 = ContextA_R$cDS$cDS_2,
                                                          rawDensity2 = ContextB_R$cDS$cDS_2)
    }else{
      print('weight jsd')
      Network[,'Div'] = colSums2(t(JSD_batch_matrix_rawMt(Ctl = ContextA_R$cDS$DS,
                                         Pert = ContextB_R$cDS$DS))*Div_weight)

      Network[,'cDiv_D1'] = colSums2(t(JSD_batch_matrix_rawMt(Ctl = ContextA_R$cDS$cDS_1,
                                                              Pert = ContextB_R$cDS$cDS_1))*Div_weight)
      Network[,'cDiv_D2'] = colSums2(t(JSD_batch_matrix_rawMt(Ctl = ContextA_R$cDS$cDS_2,
                                                              Pert = ContextB_R$cDS$cDS_2))*Div_weight)
    }

  }

  invisible(gc())
  # # i=5 # test
  # p1<-DensityPlot_raster(ContextA_R$CGDensity[i,])|DensityPlot_raster(ContextA_R$cDS$cDS_1[i,],addSideBar = F,transpose = F)|
  #   DensityPlot_raster(ContextA_R$cDS$cDS_2[i,],addSideBar = F,transpose = F)
  # p2<-DensityPlot_raster(ContextB_R$CGDensity[i,])|DensityPlot_raster(ContextB_R$cDS$cDS_1[i,],addSideBar = F,transpose = F)|
  #   DensityPlot_raster(ContextB_R$cDS$cDS_2[i,],addSideBar = F,transpose = F)
  # CombinePlots(list(p1,p2),ncol = 1)
  if(returnCDS){
    res  = list(Network=Network,
                ContextA_DS = ContextA_R$Density,
                ContextB_DS = ContextB_R$Density,
                ContextA_cDS = ContextA_R$cDS,
                ContextB_cDS = ContextB_R$cDS,
                uniCase=uniCase)

  }else{
    res  = list(Network=Network,
                ContextA_DS = ContextA_R$Density,
                ContextB_DS = ContextB_R$Density,
                uniCase=uniCase)
  }


  return(res)
}
#' @keywords internal
setNoiseToZero<-function(x,noiseIndex){
  x[noiseIndex]=0
  x
}

#' setNoiseToNA
#'
#' @param x
#' @param noiseIndex
#'
#' @return
#' @keywords internal
#'
#' @examples
setNoiseToNA<-function(x,noiseIndex){
  x[noiseIndex]=NA
  x
}

#' getKNNRadiusParallel
#'
#' @param ExpDataDT
#' @param Network
#' @param k
#' @param parallel.sz
#' @param GridData
#' @param loop.size
#' @param verbose
#'
#' @return
#' @keywords internal
#' @export
#' @noRd
#'
#' @examples
getKNNRadiusParallel <- function(ExpDataDT,
                                 Network,
                                 k=k,
                                 parallel.sz=1,
                                 GridData,
                                 parllelModel=c('foreach','bplapply')[1],
                                 loop.size=3000,
                                 verbose){
  loop.time <- ceiling(nrow(Network)/loop.size)
  # Calculate radius
  if(loop.time>1){
    # parallel
    loop.classter = infotheo::discretize(1:nrow(Network), disc="equalwidth", nbins=loop.time)[,1]
    if(parllelModel!='foreach'){
      #----------------------------------#
      # bplapply parallel----------------
      #----------------------------------#
      BPPARAM = BiocParallel::SerialParam(progressbar=verbose)
      if (parallel.sz > 1L && class(BPPARAM) == "SerialParam") {
        BPPARAM=BiocParallel::MulticoreParam(progressbar=verbose, workers=parallel.sz, tasks=100)
      } else if (parallel.sz == 1L && class(BPPARAM) != "SerialParam") {
        parallel.sz <- BiocParallel::bpnworkers(BPPARAM)
      } else if (parallel.sz > 1L && class(BPPARAM) != "SerialParam") {
        bpworkers(BPPARAM) <- parallel.sz
      }
      Radius_res <- bplapply(as.list(1:loop.time), function(j,
                                                            ExpDataDT,
                                                            Network,
                                                            loop.classter,
                                                            GridData,
                                                            k) {
        tempNet <-  Network[loop.classter%in%j,]
        Radius_md <- getKNNRadius(ExpDataDT=ExpDataDT,
                                  Network=tempNet,
                                  GridData=GridData,
                                  k=k) # Radius_md: Matrix of radius,
        #
        return(Radius_md)
      },ExpDataDT,Network,loop.classter,GridData,k,BPPARAM=BPPARAM)
      Radius_res = do.call('rbind',Radius_res)
      invisible(gc())
    }else{
      #--------------------------------#
      # foreach parallel-------------
      #-------------------------------#
      # Calculate radius

      cl <- parallel::makeCluster(parallel.sz)
      doSNOW::registerDoSNOW(cl)

      # # allowing progress bar to be used in foreach
      pb <- progress::progress_bar$new(
        format = "[:bar] :elapsed | eta: :eta",
        total = loop.time, # used for the foreach loop
        clear=FALSE)
      progress <- function(n){
        pb$tick()
      }

      opts <- list(progress = progress)

      Radius_res <- foreach(j = 1:loop.time, .combine = 'rbind',
                            .export=c("getKNNRadius"),
                            .packages=c('data.table','magrittr'),.options.snow = opts,.verbose=verbose) %dopar% {
                              tempNet <-  Network[loop.classter%in%j,]
                              Radius_md <- getKNNRadius(ExpDataDT=ExpDataDT,
                                                        Network=tempNet,
                                                        GridData=GridData,
                                                        k=k)
                              return(Radius_md)
                            }
      parallel::stopCluster(cl)

      invisible(gc())
    }
  }else{
    # Not parallel
    Radius_md <- getKNNRadius(ExpDataDT=ExpDataDT,
                              Network=Network,
                              GridData=GridData,
                              k=k) # Radius_md: Matrix of radius,
    Radius_res <- Radius_md
  }
  Radius_res
}


getKDENetParallel <- function(ExpDataDT,
                              Network,
                              h=h,
                              parallel.sz=1,
                              GridData,
                              parllelModel=c('foreach','bplapply')[1],
                              loop.size=3000,
                              verbose){
  loop.time <- ceiling(nrow(Network)/loop.size)
  # Calculate radius
  if(loop.time>1){

    # parallel
    loop.classter = infotheo::discretize(1:nrow(Network),
                                         disc="equalwidth",
                                         nbins=loop.time)[,1]
    if(parllelModel!='foreach'){
      #----------------------------------#
      # bplapply parallel----------------
      #----------------------------------#
      BPPARAM = BiocParallel::SerialParam(progressbar=verbose)
      if (parallel.sz > 1L && class(BPPARAM) == "SerialParam") {
        BPPARAM=BiocParallel::MulticoreParam(progressbar=verbose, workers=parallel.sz, tasks=100)
      } else if (parallel.sz == 1L && class(BPPARAM) != "SerialParam") {
        parallel.sz <- BiocParallel::bpnworkers(BPPARAM)
      } else if (parallel.sz > 1L && class(BPPARAM) != "SerialParam") {
        bpworkers(BPPARAM) <- parallel.sz
      }
      ds <- bplapply(as.list(1:loop.time), function(j,
                                                    ExpDataDT,
                                                    Network,
                                                    loop.classter,
                                                    GridData,
                                                    h) {
        tempNet <-  Network[loop.classter%in%j,]
        ds <- getKDEnet(ExpDataDT=ExpDataDT,
                        Network=tempNet,
                        GridData=GridData,
                        h=h) #
        #
        return(ds)
      },ExpDataDT,Network,loop.classter,GridData,h,BPPARAM=BPPARAM)
      ds = do.call('rbind',ds)
      invisible(gc())
    }else{
      #--------------------------------#
      # foreach parallel-------------
      #-------------------------------#
      # Calculate radius

      cl <- parallel::makeCluster(parallel.sz)
      doSNOW::registerDoSNOW(cl)

      # # allowing progress bar to be used in foreach
      pb <- progress::progress_bar$new(
        format = "[:bar] :elapsed | eta: :eta",
        total = loop.time, # used for the foreach loop
        clear=FALSE)
      progress <- function(n){
        pb$tick()
      }

      opts <- list(progress = progress)

      ds <- foreach(j = 1:loop.time, .combine = 'rbind',
                    .export=c("getKNNRadius"),
                    .packages=c('data.table','magrittr'),.options.snow = opts,.verbose=verbose) %dopar% {
                      tempNet <-  Network[loop.classter%in%j,]
                      ds <- getKDEnet(ExpDataDT=ExpDataDT,
                                      Network=tempNet,
                                      GridData=GridData,
                                      h=h) #
                      #
                      return(ds)
                    }
      stopCluster(cl)

      invisible(gc())
    }

  }else{
    # Not parallel
    ds <- getKDEnet(ExpDataDT=ExpDataDT,
                    Network=Network,
                    GridData=GridData,
                    h=h) #
  }
  ds
}


getKDEnet<-function(ExpDataDT,
                    Network,
                    GridData,h=h){

  # tictoc::tic()
  # bg=expand.grid(seq(1,60,10),seq(1,60,10))%>%as.matrix() # background point
  pb <- progress::progress_bar$new(
    format = "[:bar] :percent ETA: :eta",
    total = nrow(Network)  # 设置总的迭代次数
  )

  if(!is.matrix(GridData)){ # The process of normalizing data
    ds <-  matrix(0,nrow=nrow(Network),ncol=length(GridData)^2) # rad is radius.
    # Extracting the expression matrix of genes included in a biological network.
    ExpDataDT =  ExpDataDT[,c(unlist(Network[,1:2])%>%unique()),with=F]
    # Removing unused genes.
    queryData = expand.grid(GridData,GridData)
    Network = as.matrix(Network[,1:2])

    for(j in 1:nrow(Network)){
      d1 <-  as.matrix(ExpDataDT[,Network[j,],with=F])
      ds[j,]=TDA::kde(d1,queryData,h = h)
      pb$tick()
    }
  }else{
    ds <-  matrix(0,nrow=nrow(Network),ncol=nrow(GridData)^2) # rad is radius.
    # Extracting the expression matrix of genes included in a biological network.
    ExpDataDT =  ExpDataDT[,c(unlist(Network[,1:2])%>%unique()),with=F]
    # Removing unused genes.
    GridData = GridData[,colnames(ExpDataDT)]
    Network = as.matrix(Network)

    for(j in 1:nrow(Network)){
      d1 <-  as.matrix(ExpDataDT[,Network[j,],with=F])
      # d1[d1[,1]*d1[,2]==1,]=10^5
      # # d1 = rbind(d1,bg)
      # x <-  FNN::get.knnx(data=d1,
      #                     query=expand.grid(GridData[,Network[j,1]],GridData[,Network[j,2]]),
      #                     k = k, algorithm = "kd_tree")$nn.dist
      # rad[j,] <- x[,k]
      ds[j,]=TAD::kde(d1,expand.grid(GridData[,Network[j,1]],GridData[,Network[j,2]]),h = h)
      pb$tick()
    }
  }



  # tictoc::toc()
  return(ds)
}


#' CoarseGrainingDensity
#'
#' @param rawDensity a matrix,
#' @param n.coarse
#' @param rm.noise
#' @param gridPoint
#'
#' @return
#' @keywords internal
#'
#' @examples
CoarseGrainingDensity<-function(rawDensity,
                                n.coarse,
                                n.grid=NULL){
  if(is.null(n.grid)){
    gridPoint = expand.grid(1:sqrt(ncol(rawDensity)),1:sqrt(ncol(rawDensity))) # position on density plot
  }else{
    gridPoint = expand.grid(1:n.grid,1:n.grid) # position on density plot
  }
  coarse_1 <- infotheo::discretize(gridPoint[,1],
                                   disc="equalwidth",
                                   nbins=n.coarse)[,1]
  coarse_2 <- infotheo::discretize(gridPoint[,2],
                                   disc="equalwidth",
                                   nbins=n.coarse)[,1]
  coarse_label <- paste(coarse_1,coarse_2,sep = '_')
  CGDensity <-  base::rowsum(t(rawDensity),coarse_label,reorder = F)%>%t() #/table(coarse_label)[1] # coarse
  CGDensity
}

#' rmTechnologyNoise
#'
#' @param DensityD
#' @param min.pt
#' @param min.density
#'
#' @return
#' @keywords internal
#'
#' @examples
rmTechnologyNoise<-function(DensityD,
                            min.pt=0.7,
                            min.density=10/((10)^2*pi)*9){
  densityThershold <- pmin(apply(DensityD, 1, function(x)quantile(x,min.pt)),min.density)
  DensityD[(DensityD-densityThershold)<0]=0
  DensityD
}


getTechnologyNoiseThreshold<-function(DensityD,
                                      min.pt=0.7,
                                      min.density=10/((10)^2*pi)*9){
  densityThershold <- pmin(apply(DensityD, 1, function(x)quantile(x,min.pt)),min.density)
  densityThershold
}


KLD_batch_matrix<-function(Ctl,Pert,symmetrical=TRUE,base=exp(1)){
  Ctl = Ctl+1e-200
  Pert = Pert+1e-200
  Ctl = Ctl/base::rowSums(Ctl)
  Pert = Pert/base::rowSums(Pert)
  if(symmetrical){
    kld <- (rowSums(Ctl*log(Ctl/Pert,base =base ))+rowSums(Pert*log(Pert/Ctl,base =base )))/2
  }else{
    kld <- rowSums(Ctl*log(Ctl/Pert,base =base ))
  }
  kld
}


get_DS_cDS_MI_DREMI<-function(ExpDataDT,
                              Network,
                              k,
                              h=3,
                              parallel.sz,
                              n.grid,
                              n.coarse,
                              loop.size,
                              verbose,
                              label='Group',
                              CoarseGrain=TRUE,
                              kde.method = c('kde','knn'),
                              parllelModel=c('foreach','bplapply')[1]){
  if(kde.method[1]=='kde'){
    RawDS <- getKDENetParallel(ExpDataDT =ExpDataDT,
                               Network = Network,
                               h = h,
                               parallel.sz = parallel.sz,
                               GridData=1:n.grid,
                               parllelModel=parllelModel,
                               loop.size =loop.size,
                               verbose = verbose)
  }else{
    # (1)Calculate radius
    Radius_res <- getKNNRadiusParallel(ExpDataDT =ExpDataDT,
                                       Network = Network,
                                       k = k,
                                       parallel.sz = parallel.sz,
                                       # GridData=1:n.grid, # ❤
                                       GridData=sqrt(1:n.grid), # ❤ # The expression values were squared
                                       parllelModel=parllelModel,
                                       loop.size =loop.size,
                                       verbose = verbose)
    Radius_resNA = Radius_res
    # print(submt(Radius_res))
    # Radius_res[1,]=10^5 # to remove bias to
    Radius_resNA[Radius_resNA==0]=NA
    Radius_res = pmax(Radius_res,MatrixGenerics::rowMins(Radius_resNA,na.rm = T))
    # Radius_res[Radius_res<1]=1
    # Avoid infinite density when the radius is 0
    # Radius_res[Radius_res==0]=10
    # GridData = expand.grid(1:n.grid,1:n.grid)
    # BonderGrid = c(which(GridData[,1]==1),which(GridData[,2]==1))
    # Radius_res[,BonderGrid] = Radius_res[,BonderGrid]*1.5 # The radius of boundary points is increased to 1.5 times
    RawDS <- getRawDensityKNN_byRadius(Radius_res,k = k)
    # DensityPlot_raster(RawDS[1,])|
    # DensityPlot_raster(Radius_res[1,])
  }


  if(CoarseGrain){
    # (2) Calculate coarse-grained density
    CGDS <- CoarseGrainingDensity(RawDS,n.coarse,n.grid = n.grid) # Coarse Graining Density(CGDS)
    rm(RawDS)
    cDS <- getConditionalDenstiy(DS = CGDS,nx = n.coarse)
    # (4) Calculate MI
    Network[,paste('MI',label,sep = '_')] = MI_from_ked2d_v(x = cDS$DS,logbase = 2,nx = n.coarse)
    # (4) Calculate DREMI
    Network[,paste('DREMI_1',label,sep = '_')] = MI_from_ked2d_v(x = cDS$cDS_1,logbase = 2,nx = n.coarse)
    Network[,paste('DREMI_2',label,sep = '_')] = MI_from_ked2d_v(x = cDS$cDS_2,logbase = 2,nx = n.coarse)
  }else{
    CGDS = RawDS
    cDS <- getConditionalDenstiy(DS = RawDS,nx = n.grid)
    # (4) Calculate MI
    Network[,paste('MI',label,sep = '_')] = MI_from_ked2d_v(x = cDS$DS,logbase = 2,nx = n.grid)
    # (4) Calculate DREMI
    Network[,paste('DREMI_1',label,sep = '_')] = MI_from_ked2d_v(x = cDS$cDS_1,logbase = 2,nx = n.grid)
    Network[,paste('DREMI_2',label,sep = '_')] = MI_from_ked2d_v(x = cDS$cDS_2,logbase = 2,nx = n.grid)
  }


  return(list(Network=Network,Density=CGDS,cDS=cDS))
}

get_cDS_MI_DREMI_fromRawDS<-function(RawDS,
                                     Network,
                                     n.grid,
                                     n.coarse,
                                     label='Group',
                                     CoarseGrain=TRUE){
  if(ncol(RawDS)==n.coarse^2){
    CGDS = RawDS
    cDS <- getConditionalDenstiy(DS = RawDS,nx = n.coarse)
    # (1) Calculate MI
    Network[,paste('MI',label,sep = '_')] = MI_from_ked2d_v(x = cDS$DS,logbase = 2,nx = n.coarse)
    # (2) Calculate DREMI
    Network[,paste('DREMI_1',label,sep = '_')] = MI_from_ked2d_v(x = cDS$cDS_1,logbase = 2,nx = n.coarse)
    Network[,paste('DREMI_2',label,sep = '_')] = MI_from_ked2d_v(x = cDS$cDS_2,logbase = 2,nx = n.coarse)
  }
  if(ncol(RawDS)==n.grid^2){
    if(CoarseGrain){
      # (2)Calculate coarse-grained density
      CGDS <- CoarseGrainingDensity(RawDS,n.coarse,n.grid = n.grid) # Coarse Graining Density(CGDS)
      cDS <- getConditionalDenstiy(DS = CGDS,nx = n.coarse)
      # (3) Calculate MI
      Network[,paste('MI',label,sep = '_')] = MI_from_ked2d_v(x = cDS$DS,logbase = 2,nx = n.coarse)
      # (4) Calculate DREMI
      Network[,paste('DREMI_1',label,sep = '_')] = MI_from_ked2d_v(x = cDS$cDS_1,logbase = 2,nx = n.coarse)
      Network[,paste('DREMI_2',label,sep = '_')] = MI_from_ked2d_v(x = cDS$cDS_2,logbase = 2,nx = n.coarse)
    }else{
      CGDS = RawDS
      cDS <- getConditionalDenstiy(DS = RawDS,nx = n.grid)
      # (3) Calculate MI
      Network[,paste('MI',label,sep = '_')] = MI_from_ked2d_v(x = cDS$DS,logbase = 2,nx = n.grid)
      # (4) Calculate DREMI
      Network[,paste('DREMI_1',label,sep = '_')] = MI_from_ked2d_v(x = cDS$cDS_1,logbase = 2,nx = n.grid)
      Network[,paste('DREMI_2',label,sep = '_')] = MI_from_ked2d_v(x = cDS$cDS_2,logbase = 2,nx = n.grid)
    }
  }
  return(list(Network=Network,Density=CGDS,cDS=cDS))
}


getConditionalDenstiy<-function(DS,nx=NULL){
  if(is.null(nx)){
    nx = sqrt(ncol(DS))
  }
  gridCorasePoint = expand.grid(1:nx,1:nx)
  DS = DS/rowSums(DS) # to p value
  # (3) Calculate conditional density
  # G1->G2
  cds_G1 = NormalizeToOne_bygorup(DS,
                                  gridCorasePoint[,1]) # normalize to one by condition G1
  cds_G1[is.nan(cds_G1)]=0 # set nan to zeros
  cds_G1 = cds_G1/rowSums(cds_G1)  # to p value
  # G2->G1
  cds_G2 = NormalizeToOne_bygorup(DS,
                                  gridCorasePoint[,2])# normalize to one by condition G2
  cds_G2[is.nan(cds_G2)]=0 #DensityPlot_raster(cds_G2_CA[1,],addSideBar = F,transpose = T) # test function
  cds_G2 = cds_G2/rowSums(cds_G2)  # to p value
  return(list(DS=DS,cDS_1=cds_G1,cDS_2=cds_G2))
}

getConditionalDenstiy_NoiseRemove<-function(DS,nx=NULL,threshold_v){
  if(is.null(nx)){
    nx = sqrt(ncol(DS))
  }
  gridCorasePoint = expand.grid(1:nx,1:nx)
  # DS = DS/rowSums(DS) # to p value
  # (3) Calculate conditional density
  # G1->G2
  cds_G1 = NormalizeToOne_bygorup_withNoiseRemove(DS,
                                                  gridCorasePoint[,1],threshold_v) # normalize to one by condition G1
  cds_G1[is.na(cds_G1)]=0 # set nan to zeros
  cds_G1 = cds_G1/rowSums(cds_G1)  # to p value
  # G2->G1
  cds_G2 = NormalizeToOne_bygorup_withNoiseRemove(DS,
                                                  gridCorasePoint[,2],threshold_v)# normalize to one by condition G2
  cds_G2[is.na(cds_G2)]=0 #DensityPlot_raster(cds_G2_CA[1,],addSideBar = F,transpose = T) # test function
  cds_G2 = cds_G2/rowSums(cds_G2)  # to p value
  return(list(DS=DS,cDS_1=cds_G1,cDS_2=cds_G2))
}


getConditionalDenstiy_NoiseRemove_twoDS<-function(DS1,DS2,min.pt,nx=NULL,min.density=0.286789){


  if(is.null(nx)){
    nx = sqrt(ncol(DS1))
  }
  gridCorasePoint = expand.grid(1:nx,1:nx)
  threshold_r = getTechnologyNoiseThreshold_byColValue(DS1 = DS1,
                                                       DS2 = DS2,
                                                       min.pt = 0.5,
                                                       nx = nx,
                                                       min.density = min.density)




  getCDSList <- function(DS,GroupLabel,ThresoldList){
    cds_G1 = NormalizeToOne_bygorup_withNoiseRemove(DS,
                                                    GroupLabel[,1],threshold_v = ThresoldList$rowTh) # normalize to one by condition G1
    cds_G1[is.na(cds_G1)]=0 # set nan to zeros
    cds_G1 = cds_G1/rowSums(cds_G1)  # to p value
    # G2->G1
    cds_G2 = NormalizeToOne_bygorup_withNoiseRemove(DS,
                                                    GroupLabel[,2],threshold_v=ThresoldList$colTh)# normalize to one by condition G2
    cds_G2[is.na(cds_G2)]=0 #DensityPlot_raster(cds_G2_CA[1,],addSideBar = F,transpose = T) # test function
    cds_G2 = cds_G2/rowSums(cds_G2)  # to p value
    return(list(cDS_1=cds_G1,cDS_2=cds_G2))
  }
  DS1_cds =getCDSList(DS = DS1,GroupLabel = gridCorasePoint,ThresoldList = threshold_r$DS1_th)
  DS2_cds =getCDSList(DS = DS2,GroupLabel = gridCorasePoint,ThresoldList = threshold_r$DS2_th)



  return(list(DS1_cds=DS1_cds,DS2_cds=DS2_cds))
}



NormalizeToOne_bygorup_withNoiseRemove<-function(ColData,label,threshold_v) {

  if(!is.matrix(ColData)){
    ColData = matrix(ColData,nrow = 1)
  }
  uni_label = unique(label)
  for(i in 1:length(uni_label)){
    RmaxData=rowMaxs(ColData[,label%in%uni_label[i],drop=F])
    RmaxData[(RmaxData-threshold_v)<0]=NA
    ColData[,label%in%uni_label[i]] = ColData[,label%in%uni_label[i],drop=F]/RmaxData
  }
  ColData
}




rowMaxsByColGroup<-function(m,colGroup,na.rm=FALSE){

  uni_label = unique(colGroup)
  maxMatrix = matrix(NA,nrow = nrow(m),ncol = length(uni_label))
  for(i in 1:length(uni_label)){
    RmaxData=rowMaxs(m[,colGroup%in%uni_label[i],drop=F],na.rm = na.rm)
    maxMatrix[,i]=RmaxData
  }
  colnames(maxMatrix)=uni_label
  rownames(maxMatrix)=rownames(m)
  maxMatrix
}

getTechnologyNoiseThreshold_byColValue <- function(DS1,
                                                   DS2,
                                                   min.pt=0.5,
                                                   nx=20,
                                                   min.density=k/((k)^2*pi)*(n.grid^2/n.coarse^2)){
  GroupLabel = expand.grid(1:nx,1:nx)
  th1 =  getTechnologyNoiseThreshold(DS1,
                                     min.pt = min.pt,
                                     min.density = min.density)
  th2 =  getTechnologyNoiseThreshold(DS2,
                                     min.pt = min.pt,
                                     min.density = min.density)
  NoiseGridData = (DS1-th1)<0&(DS2-th2)<0

  rowMaxValueDS1 = rowMaxsByColGroup(setNoiseToNA(DS1,NoiseGridData),GroupLabel[,1],na.rm = TRUE)
  rowMaxValueDS1[is.infinite(rowMaxValueDS1)]=NA
  rowTh_minMax_DS1 = rowMins( rowMaxValueDS1,na.rm = TRUE)
  colMaxValueDS1 = rowMaxsByColGroup(setNoiseToNA(DS1,NoiseGridData),GroupLabel[,2],na.rm = TRUE)
  colMaxValueDS1[is.infinite(colMaxValueDS1)]=NA
  colTh_minMax_DS1 = rowMins( colMaxValueDS1,na.rm = TRUE)
  #
  rowMaxValueDS2 = rowMaxsByColGroup(setNoiseToNA(DS2,NoiseGridData),GroupLabel[,1],na.rm = TRUE)
  rowMaxValueDS2[is.infinite(rowMaxValueDS2)]=NA
  rowTh_minMax_DS2 = rowMins( rowMaxValueDS2,na.rm = TRUE)
  colMaxValueDS2 = rowMaxsByColGroup(setNoiseToNA(DS2,NoiseGridData),GroupLabel[,2],na.rm = TRUE)
  colMaxValueDS2[is.infinite(colMaxValueDS2)]=NA
  colTh_minMax_DS2 = rowMins( colMaxValueDS2,na.rm = TRUE)
  # MaxValueDS1 = cbind(rowMaxsByColGroup(setNoiseToNA(DS1,NoiseGridData),GroupLabel[,1],na.rm = TRUE),
  #       rowMaxsByColGroup(setNoiseToNA(DS1,NoiseGridData),GroupLabel[,2],na.rm = TRUE))
  # MaxValueDS1[is.infinite(MaxValueDS1)]=NA
  # th_minMax_DS1 = rowMins( MaxValueDS1,na.rm = TRUE)
  # #
  # MaxValueDS2 = cbind(rowMaxsByColGroup(setNoiseToNA(DS2,NoiseGridData),GroupLabel[,1],na.rm = TRUE),
  #                     rowMaxsByColGroup(setNoiseToNA(DS2,NoiseGridData),GroupLabel[,2],na.rm = TRUE))
  # MaxValueDS2[is.infinite(MaxValueDS2)]=NA
  # th_minMax_DS2 = rowMins( MaxValueDS2,na.rm = TRUE)
  DS1_th = data.frame(rowTh=rowTh_minMax_DS1,colTh=colTh_minMax_DS1)
  DS2_th = data.frame(rowTh=rowTh_minMax_DS2,colTh=colTh_minMax_DS2)
  return(list(DS1_th=DS1_th,DS2_th=DS2_th))
}

removeNoiseInCDS<-function(cDS,NoiseGridData,ColGroupLabel,nx=NULL){
  if(is.null(nx)){
    nx = sqrt(ncol(cDS))
  }

  NoiseRowOrCol=base::rowsum(t(NoiseGridData+0),ColGroupLabel,reorder = F)%>%t()==nx

  uni_label = unique(ColGroupLabel)
  for(i in 1:length(uni_label)){
    cDS[NoiseRowOrCol[,uni_label[i]],ColGroupLabel%in%uni_label[i]] = 0
  }
  cDS

}

getTechnologyNoiseData_towDS <- function(DS1,
                                         DS2,
                                         min.pt=0.5,
                                         min.density=k/((k)^2*pi)*(n.grid^2/n.coarse^2)){

  th1 =  getTechnologyNoiseThreshold(DS1,
                                     min.pt = min.pt,
                                     min.density = min.density)
  th2 =  getTechnologyNoiseThreshold(DS2,
                                     min.pt = min.pt,
                                     min.density = min.density)
  NoiseGridData = (DS1-th1)<0&(DS2-th2)<0
  return(NoiseGridData)
}

#' removeNoiseInCDS_towDS
#'
#' @param DS1
#' @param DS2
#' @param cDSList_1 data.frame(cDS_1=cDS_1,cDS_2=cDS_2)
#' @param cDSList_2 data.frame(cDS_1=cDS_1,cDS_2=cDS_2)
#' @param nx
#'
#' @return
#' @keywords internal
#'
#' @examples
#'
removeNoiseInCDS_towDS<-function(DS1,
                                 DS2,
                                 cDSList_1,
                                 cDSList_2,
                                 nx=NULL,
                                 min.pt=0.5,
                                 min.density=0.2864789){
  if(is.null(nx)){
    nx = sqrt(ncol(DS1))
  }
  GroupLabel = expand.grid(1:nx,1:nx)
  NoiseGridData<- getTechnologyNoiseData_towDS(DS1,
                                               DS2,
                                               min.pt=min.pt,
                                               min.density=min.density)
  # Context A
  cDSList_1$cDS_1 = removeNoiseInCDS(cDSList_1$cDS_1,NoiseGridData,ColGroupLabel = GroupLabel[,1],nx = nx)
  cDSList_1$cDS_2 = removeNoiseInCDS(cDSList_1$cDS_2,NoiseGridData,ColGroupLabel = GroupLabel[,2],nx = nx)
  # Context B
  cDSList_2$cDS_1 = removeNoiseInCDS(cDSList_2$cDS_1,NoiseGridData,ColGroupLabel = GroupLabel[,1],nx = nx)
  cDSList_2$cDS_2 = removeNoiseInCDS(cDSList_2$cDS_2,NoiseGridData,ColGroupLabel = GroupLabel[,2],nx = nx)
  return(list(cDSList_1=cDSList_1,cDSList_2=cDSList_2))

}


#' getKNNRadius
#' Calculation of distance radius based on KNN method.
#' @param ExpData a data.table representing the gene expression matrix. The columns represent genes, while the rows represent samples or cells.
#' @param Network  a two-column data.frame representing a biological network.
#' @param GridData  a query data matrix. The rows are genes, and the columns are query points.
#' @param k  a positive integer, which is the k parameter in the KNN method.
#'
#' @return
#' @keywords internal
#'
#' @examples
getKNNRadius <- function(ExpDataDT,
                         Network,
                         GridData,
                         k=10){

  # tictoc::tic()
  pb <- progress::progress_bar$new(
    format = "[:bar] :percent ETA: :eta",
    total = nrow(Network)  # 设置总的迭代次数
  )
  # bg=expand.grid(seq(1,60,10),seq(1,60,10))%>%as.matrix() # background point
  if(!is.matrix(GridData)){ # normalized gene expression
    rad <-  matrix(0,nrow=nrow(Network),ncol=length(GridData)^2) # rad is radius.
    # Extracting the expression matrix of genes included in a biological network.
    ExpDataDT =  ExpDataDT[,c(unlist(Network[,1:2])%>%unique()),with=F]
    # Removing unused genes.
    queryData = expand.grid(GridData,GridData)
    Network = as.matrix(Network)

    for(j in 1:nrow(Network)){
      d1 <-  as.matrix(ExpDataDT[,Network[j,1:2]%>%as.character(),with=F])
      # d1 = rbind(d1,bg)
      d1[d1[,1]*d1[,2]==1,]=10^5
      x <-  FNN::get.knnx(data=d1,
                          query=queryData,
                          k = k, algorithm = "kd_tree")$nn.dist
      rad[j,] <- x[,k]
      pb$tick()
    }
  }else{
    rad <-  matrix(0,nrow=nrow(Network),ncol=nrow(GridData)^2) # rad is radius.
    # Extracting the expression matrix of genes included in a biological network.
    ExpDataDT =  ExpDataDT[,c(unlist(Network[,1:2])%>%unique()),with=F]
    # Removing unused genes.
    GridData = GridData[,colnames(ExpDataDT)]
    Network = as.matrix(Network)

    for(j in 1:nrow(Network)){
      d1 <-  as.matrix(ExpDataDT[,Network[j,1:2]%>%as.character(),with=F])
      d1[d1[,1]*d1[,2]==1,]=10^5
      # d1 = rbind(d1,bg)
      x <-  FNN::get.knnx(data=d1,
                          query=expand.grid(GridData[,Network[j,1]],GridData[,Network[j,2]]),
                          k = k, algorithm = "kd_tree")$nn.dist
      rad[j,] <- x[,k]
      pb$tick()
    }
  }



  # tictoc::toc()
  return(rad)
}

#' getRawDensityKNN_byRadius
#' k/(π*r^2)
#' @param Radius_md a matrix,
#'
#' @return
#' @keywords internal
#'
#' @examples
getRawDensityKNN_byRadius<-function(Radius_md,k){
  # calculating radius
  # Radius_md[Radius_md>=k]=k #
  # bg = expand.grid(seq(1,n.grid,n.grid/6),seq(1,n.grid,n.grid/6))%>%as.matrix()
  # bg.r=FNN::get.knnx(data=bg,query=expand.grid(1:n.grid,1:n.grid),
  #                    k = k, algorithm = "kd_tree")$nn.dist[,k]
  # bg.ds = k/((bg.r+10^-64)^2*pi)
  rawDS <- k/((Radius_md+10^-64)^2*pi) # calculating circular area
  rawDS
  # rawDS+mean(bg.ds)
}

#' NormalizeToOne_bygorup
#' Rescale each value of the conditional density estimate by its column-maximum
#' @param ColData
#' @param label
#'
#' @return
#' @keywords internal
#'
#' @examples
NormalizeToOne_bygorup <- function(ColData,label) {
  if(!is.matrix(ColData)){
    ColData = matrix(ColData,nrow = 1)
  }
  uni_label = unique(label)
  for(i in 1:length(uni_label)){
    ColData[,label%in%uni_label[i]] = ColData[,label%in%uni_label[i],drop=F]/MatrixGenerics::rowMaxs(ColData[,label%in%uni_label[i],drop=F])
  }
  ColData
}


#' mappingMinMaxRow
#' maping each row of data to between minD and maxD
#' @param ma a matrix
#' @param minD minimal value. The default is 1
#' @param maxD maximum value. The default is 60
#' @param rmZero Whether to exclude 0. The default is TRUE
#'
#' @return a mapped matrix
#' @keywords internal
#'
#' @examples
#' mappingMinMaxRow(matrix(0:24,5),minD=1,maxD=60,rmZero = TRUE)
mappingMinMaxRow<-function(ma,
                           minD=1,
                           maxD=60,
                           rmZero = TRUE){
  mappedma=apply(ma,1,function(x){
    if(rmZero){
      x1 = x[x!=0]
      x_range = range(x1)
      x2 = (x-x_range[1])/(x_range[2]-x_range[1])*(maxD-minD)+minD
      x2[x==0]=0
      x2

    }else{
      x_range = range(x)
      (x-x_range[1])/(x_range[2]-x_range[1])*(maxD-minD)+minD
    }

  })
  t(mappedma)
}

JSD_batch_matrix <- function(Ctl,Pert,base=exp(1)){
  Ctl = Ctl+1e-200
  Pert = Pert+1e-200
  Ctl = Ctl/base::rowSums(Ctl)
  Pert = Pert/base::rowSums(Pert)
  M <- (Pert + Ctl) / 2

  jsd <- (rowSums(Ctl * log(Ctl / M,base = base))+rowSums(Pert * log(Pert / M,base = base)))/2
  jsd
}

JSD_batch_matrix_rawMt_robust <- function(rawDensity1,
                                          rawDensity2){
  #Enhance robustness and reduce sensitivity to subtle variations

  Indenx_r1  <- matrix(1:(22*22),22,22)
  Indenx_cen <- Indenx_r1[2:21,2:21]

  Indenx_up  <- Indenx_r1[1:20,2:21]
  Indenx_dn  <- Indenx_r1[3:22,2:21]
  Indenx_right  <- Indenx_r1[2:21,3:22]
  Indenx_left  <- Indenx_r1[2:21,1:20]
  Indenx_up_left  <- Indenx_r1[1:20,1:20]
  Indenx_up_right  <- Indenx_r1[1:20,3:22]
  Indenx_dn_left  <- Indenx_r1[3:22,1:20]
  Indenx_dn_right  <- Indenx_r1[3:22,3:22]

  TepMatrix1 <- matrix(0,nrow(rawDensity1),ncol = 22*22)
  TepMatrix2 <- matrix(0,nrow(rawDensity1),ncol = 22*22)
  TepMatrix3 <- matrix(0,nrow(rawDensity1),ncol = 22*22)
  TepMatrix4 <- matrix(0,nrow(rawDensity1),ncol = 22*22)
  TepMatrix5 <- matrix(0,nrow(rawDensity1),ncol = 22*22)
  TepMatrix6 <- matrix(0,nrow(rawDensity1),ncol = 22*22)
  TepMatrix7 <- matrix(0,nrow(rawDensity1),ncol = 22*22)
  TepMatrix8 <- matrix(0,nrow(rawDensity1),ncol = 22*22)


  TepMatrix1[,Indenx_up] <- rawDensity1
  TepMatrix2[,Indenx_dn] <- rawDensity1
  TepMatrix3[,Indenx_right] <- rawDensity1
  TepMatrix4[,Indenx_left] <- rawDensity1
  TepMatrix5[,Indenx_up_left] <- rawDensity1
  TepMatrix6[,Indenx_up_right] <- rawDensity1
  TepMatrix7[,Indenx_dn_left] <- rawDensity1
  TepMatrix8[,Indenx_dn_right] <- rawDensity1

  jsd1 <- JSD_batch_matrix_rawMt(TepMatrix1[,Indenx_cen,drop=F],rawDensity2)
  jsd2 <- JSD_batch_matrix_rawMt(TepMatrix2[,Indenx_cen,drop=F],rawDensity2)
  jsd3 <- JSD_batch_matrix_rawMt(TepMatrix3[,Indenx_cen,drop=F],rawDensity2)
  jsd4 <- JSD_batch_matrix_rawMt(TepMatrix4[,Indenx_cen,drop=F],rawDensity2)
  jsd5 <- JSD_batch_matrix_rawMt(TepMatrix5[,Indenx_cen,drop=F],rawDensity2)
  jsd6 <- JSD_batch_matrix_rawMt(TepMatrix6[,Indenx_cen,drop=F],rawDensity2)
  jsd7 <- JSD_batch_matrix_rawMt(TepMatrix7[,Indenx_cen,drop=F],rawDensity2)
  jsd8 <- JSD_batch_matrix_rawMt(TepMatrix8[,Indenx_cen,drop=F],rawDensity2)

  jsd0 <- JSD_batch_matrix_rawMt(rawDensity1,rawDensity2)

  jsd_min <- pmin(jsd1,jsd2,jsd3,jsd4,jsd5,jsd6,jsd7,jsd8,jsd0)

  # jsd_rb <- (jsd0+jsd_min)/2
  jsd_rb <- jsd_min
  jsd_rb
  # rawDensity
}



MI_from_ked2d_v <- function(x,logbase=2,nx=20){
  #z_0 = z
  xyLalbel=expand.grid(1:nx,1:nx)

  x = pmax(x,1e-68)
  P_xy = x/rowSums(x)
  P_x = base::rowsum(t(P_xy),xyLalbel[,1],reorder = F)#/table(coarse_label)[1] # coarse
  P_x = t(P_x[xyLalbel[,1],])
  P_y = base::rowsum(t(P_xy),xyLalbel[,2],reorder = F)#/table(coarse_label)[1] # coarse
  P_y = t(P_y[xyLalbel[,2],])
  MI = rowSums(P_xy*log(P_xy/(P_x*P_y),logbase))
  MI
}

JSD_batch_matrix_rawMt <- function(Ctl,Pert,base=exp(1)){
  Ctl = Ctl+1e-200
  Pert = Pert+1e-200
  Ctl = Ctl/base::rowSums(Ctl)
  Pert = Pert/base::rowSums(Pert)
  M <- (Pert + Ctl) / 2

  jsd_mt <- (Ctl * log(Ctl / M,base = base)+Pert * log(Pert / M,base = base))/2
  jsd_mt
}

