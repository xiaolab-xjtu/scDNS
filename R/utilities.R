#' rmNotExpressedGene
#' Remove genes with low expression ratio in cells
#' @param sob seurat object
#' @param min.cellN The minimum number of cells in which a gene can be expressed (def:10)
#'
#' @return sob seurat object
#' @export
#'
#' @examples
rmNotExpressedGene<-function(sob,min.cellN=10){

  rmSomeGene_Seurate <-   function(object, pattern = NULL, features = NULL, col.name = NULL,
                                   assay = NULL){
    assay <- assay %||% DefaultAssay(object = object)
    if (!is.null(x = features) && !is.null(x = pattern)) {
      warning("Both pattern and features provided. Pattern is being ignored.")
    }
    features <- features %||% grep(pattern = pattern, x = rownames(x = object[[assay]]),
                                   value = TRUE)
    RemianFeatures=rownames(object[[assay]])[!rownames(object[[assay]])%in%features]
    message('before:')
    print(dim(object))
    object = subset(x = object,features =RemianFeatures)
    message('after:')
    print(dim(object))

    return(object)
  }
  # Output a logical vector for every gene on whether the more than zero counts per cell
  counts <- GetAssayData(object = sob, slot = "counts")
  nonzero <- counts > 0

  # Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
  rm_genes <- Matrix::rowSums(nonzero) < min.cellN
  sob=rmSomeGene_Seurate(sob,features = names(rm_genes)[rm_genes])
  sob
}

#' Magic4MultipleData
#' Imputing Gene Expression with Rmagic for Different Conditions
#' @param sob seurat object
#' @param split.by group name to split data
#'
#' @return seurat object
#' @export
#'
#' @examples
Magic4MultipleData <- function(sob,split.by=NULL){
  sob_list <- SplitObject(sob,split.by =split.by )
  sob_list <- lapply(sob_list, function(x)Rmagic::magic(x))
  sob <- merge(sob_list[[1]],unlist(sob_list[2:length(sob_list)]))
  sob
}

#' CalGeneInfAmount
#' Evaluate the ability of genes to respond to valid information based on the standard deviation of gene expression
#'
#' @param ExpMat Normalized gene expression matrix, n genes × m samples
#' @param probs Quantile of non-zero gene expression values, default is 0.99
#'
#' @return
#' Returns the informative ability of each gene, in the range 0-1. The closer to 1, the greater the ability.
#' @keywords internal
#'
#' @examples
CalGeneInfAmount <- function(ExpMat,
                             probs = 0.99,
                             rescale=FALSE,
                             plot = FALSE){
  # Likelihoodx = CalGeneInfAmount(ExpMat = iAT2_data@assays$RNA@data)
  # iAT2_data@assays$RNA@meta.features$PercentExp = rowSums2(iAT2_data@assays$RNA@counts!=0)/ncol(iAT2_data)
  # iAT2_data@assays$RNA@meta.features$Likelihood = Likelihoodx
  # fit_sigmoid(iAT2_data@assays$RNA@meta.features$PercentExp, iAT2_data@assays$RNA@meta.features$Likelihoodx)$p +
  #   xlab('Percent Expressed')+ ylab('Likelihoodx')+AddGrid()
  if(rescale){
    # normalize to one
    ExpMat[ExpMat==0] = NA
    maxr = rowQuantiles(ExpMat,probs = probs,na.rm = T)
    maxr[maxr==0]=Inf
    ExpMat <- ExpMat/maxr
    ExpMat[ExpMat>1] = 1
    ExpMat[is.na(ExpMat)] = 0
  }

  # calculte sd
  sdx = matrixStats::rowSds(as.matrix(ExpMat))
  LikelihoodModel <- likelihoodFromMixEM_2d(log10(sdx),plot = plot)
  Likelihoodx <- LikelihoodModel$Cal_likelihood(log10(sdx),
                                                fitModel = LikelihoodModel$fitModel,
                                                likelihoodToLower = F)
  names(Likelihoodx) <- rownames(ExpMat)
  Likelihoodx

}






#' SmoothOutLier_UP
#'
#' @param mt gene expression matrix, n genes × m samples
#' @param zero.rm logical, whether to exclude zero values
#' @param percentile
#'
#' @return
#' @keywords internal
#'
#' @examples
SmoothOutLier_UP <- function(mt,zero.rm=TRUE,percentile=0.99){
  mt2 = mt
  if(zero.rm){
    mt[mt==0]=NA
  }
  Percentile_q1 = matrixStats::rowQuantiles(mt,probs =  percentile,na.rm = zero.rm)
  Percentile_q2 = matrixStats::rowQuantiles(mt,probs =  percentile*2-1,na.rm = zero.rm)
  OutL = DetectOutlierBy3kBoxplot_UP(mt2,zero.rm = zero.rm)
  OutL[which.na(OutL)]= FALSE
  # OutL[is.na(OutL)]=FALSE
  mt[OutL]=NA
  RMax_2 = matrixStats::rowMaxs(mt,na.rm = TRUE) # 非outlier的每一行的最大值
  mt3 = mt2
  mt3[!OutL]=Inf
  mt3 = matrixStats::rowRanks(mt3)*as.vector(Percentile_q1-Percentile_q2)/matrixStats::rowSums2(OutL)+RMax_2 # outlier的排序*间隔+最大值
  mt3[is.infinite(RMax_2),]=0
  mt2[OutL]=mt3[as.matrix(OutL)]
  mt2

}





#' DetectOutlierBy3kBoxplot_UP
#' Detect outliers
#' @param m  gene expression matrix, n genes × m samples
#' @param zero.rm logical, whether to exclude zero values
#'
#' @return
#' @keywords internal
#'
#' @examples
DetectOutlierBy3kBoxplot_UP <- function(m,zero.rm=TRUE){
  if(!zero.rm){
    eachPercent = sparseMatrixStats::rowQuantiles(m,na.rm = TRUE,probs = c(0.25,0.75),drop = F)
    gab = (eachPercent[,2]-eachPercent[,1])*3
    up_lim = eachPercent[,2]+gab
    # dn_lim = eachPercent[,1]-gab
    # (m-up_lim)>0|(m-dn_lim)<0
    (m-up_lim)>0
  }else{
    m2 =m
    m[m==0]=NA
    eachPercent = sparseMatrixStats::rowQuantiles(m,na.rm = TRUE,probs = c(0.25,0.75),drop = F)
    gab = (eachPercent[,2]-eachPercent[,1])*3
    up_lim = eachPercent[,2]+gab
    # dn_lim = eachPercent[,1]-gab
    # (m2-up_lim)>0|(m2-dn_lim)<0
    (m2-up_lim)>0
  }


}




#' getDoubleDropout
#'
#' @param Network
#' @param counts gene counts expression matrix, n genes × m samples
#'
#' @return
#' @keywords internal
#'
#' @examples
getDoubleDropout <- function(Network,counts){
  Network = preFilterNet(BioNet = Network,CountData = counts)
  counts = counts[rownames(counts)%in%(unlist(Network[,1:2])%>%unique()),]
  # dropoutMatrix = ncol(counts)-Matrix::tcrossprod(counts==0 )
  dropoutMatrix = Matrix::tcrossprod(counts!=0 ) # n
  # dropoutMatrix = Matrix::tcrossprod(counts!=0 ) # number of cells expressing both G1 and G2
  Network$npoint = dropoutMatrix[sub2ind(match(Network[,1],rownames(dropoutMatrix)),
                                         match(Network[,2],rownames(dropoutMatrix)),
                                         nrow = nrow(dropoutMatrix),ncol = ncol(dropoutMatrix))]
  Network$npoint_log = log10(Network$npoint+1)
  Network$npoint_ds = round(Network$npoint_log,1)
  Network
}

#' preFilterNet
#'
#' @param BioNet BioNet
#' @param CountData gene expression counts matrix ,genes × cells
#'
#' @return
#' @keywords internal
#'
#' @examples
preFilterNet <- function(BioNet,CountData){
  #@@@@@@@@@@@@@
  message('Step1: remove edges that are not in the expression data')
  reamin_label <- BioNet[,1]%in%rownames(CountData)&BioNet[,2]%in%rownames(CountData)
  BioNet_filterd = BioNet[reamin_label,]
  message('remain  edges (farction):')
  print(nrow(BioNet_filterd)/nrow(BioNet))
  #@@@@@@
  message('Step2: remove duplicated edges')
  BioNet_filterd = detect_Duplicate_edge4net(BioNet_filterd,ingnoreDireaction = TRUE,returnLogical = FALSE)
  message('remain non-dupicated edges (farction):')
  print(nrow(BioNet_filterd)/nrow(BioNet))

  #@@@@@@@@@@@@@
  message('Step3: selef connected interactions')
  reamin_label <- BioNet_filterd[,1]!=BioNet_filterd[,2]
  BioNet_filterd = BioNet_filterd[reamin_label,]
  message('remain  edges (farction):')
  print(nrow(BioNet_filterd)/nrow(BioNet))
  BioNet_filterd
}




#' Detect and handle duplicate edges in a network matrix or data frame.
#'
#' This function identifies duplicate edges in a two-column network representation.
#' It can optionally treat edges as undirected (ignoring direction) and return either
#' a logical vector indicating non-duplicate rows or a deduplicated network.

#'
#' @param net A two-column matrix or data frame representing edges (e.g., from-to pairs).
#' @param ingnoreDireaction Logical. If TRUE, edges are treated as undirected,
#' so (A, B) and (B, A) are considered duplicates. Default is TRUE.
#' @param returnLogical Logical. If TRUE, returns a logical vector where
#' TRUE indicates non-duplicate rows; otherwise, returns the deduplicated network.
#'
#' @return If returnLogical is TRUE, a logical vector of the same length as the number of rows in net,
#' with TRUE for unique (non-duplicate) edges. Otherwise, a matrix or data frame
#' containing only the unique edges.
#' @keywords internal

#'
#' @examples
detect_Duplicate_edge4net <- function(net,
         ingnoreDireaction = TRUE,
         returnLogical = FALSE){

  #
  if(ingnoreDireaction){
    net2 <- net[,1:2]
    net2[net[,1]>net[,2],] <- net[net[,1]>net[,2],c(2,1)]
    lable_1 <-  paste(net2[,1],net2[,2])
    Logical_1 <- duplicated(lable_1)
  }else{
    lable_1 <-  paste(net[,1],net[,2])
    # remove 1
    Logical_1 <- duplicated(lable_1)
  }
  if(returnLogical){
    return(!Logical_1)
  }else{
    net = net[!Logical_1,]
    return(net)
  }
}


#' Convert 2D subscripts (row and column indices) to linear indices.
#'
#' This function converts matrix subscripts (row and column positions) into a single linear index,
#' assuming column-major order (as used in R). It is useful for indexing elements in a matrix
#' when working with flattened (vectorized) representations.
#'
#' @param r Row index (1-based).
#' @param c Column index (1-based).
#' @param nrow Number of rows in the matrix.
#' @param ncol Number of columns in the matrix (not used in computation but included for clarity).
#'
#' @return A numeric scalar representing the linear index corresponding to the given row and column.
#' @export
#'
#' @examples
sub2ind <- function(r=2,c=3,nrow=5,ncol=5){
  ind = (c-1)*nrow + r
  return(ind)
}


#' Find indices of missing values (NA) in a vector or array.

#'
#' This function returns the positions of NA values in the input object, by wrapping
#' which(is.na(x), ...) for convenience.
#'
#' @param x An object (typically a vector or array) to check for missing values.
#' @param ... Additional arguments passed to which(), such as arr.ind = TRUE
#' for array indexing.
#'
#' @return An integer vector (or array index if arr.ind = TRUE is used) giving the
#' positions of NA values in x.
#' @export
#'
#' @examples
which.na <- function(x,...){
  which(is.na(x),...)
}


likelihoodFromMixEM_2d <- function(x,plot=TRUE){
  n.distrubtion=2 # Mixture of three normal distributions.
  muStart=seq(min(x),max(x),(max(x)-min(x))/(n.distrubtion+1))
  muStart = muStart[2:(length(muStart)-1)]
  fitModel <- mixtools::normalmixEM(x, lambda = 1/n.distrubtion, mu = muStart, sigma = sd(x)/3)
  if(plot){
    plot(fitModel, density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8,
         main2="", xlab2="x")
  }

  Cal_likelihood<-function(x,fitModel,likelihoodToLower=TRUE){
    if(likelihoodToLower){
      D1=pnorm(x,fitModel$mu[1],sd = fitModel$sigma[1],lower.tail = F) # perturbed
      D2=pnorm(x,fitModel$mu[2],sd = fitModel$sigma[2],lower.tail = T) # commn
      LR = D1/(D1+D2)
      LR
    }else{
      D1=pnorm(x,fitModel$mu[1],sd = fitModel$sigma[1],lower.tail = F) # perturbed
      D2=pnorm(x,fitModel$mu[2],sd = fitModel$sigma[2],lower.tail = T) # commn
      LR = D2/(D1+D2)
      LR
    }
    # D1=pnorm(x,fitModel$mu[3],sd = fitModel$sigma[3]) # perturbed
    # D2=pnorm(x,fitModel$mu[2],sd = fitModel$sigma[2],lower.tail = F) # commn
    # D3=pnorm(x,fitModel$mu[1],sd = fitModel$sigma[2],lower.tail = F) # control
    # LR = D1/(D1+D2+D3)
    # LR
  }

  return(list(fitModel=fitModel,Cal_likelihood=Cal_likelihood))
}


getSubNetByNode <- function(Net,
                            Node,
                            exclude=FALSE,
                            Model=c('or','and')[1]){
  # getSubNetByNode
  if(exclude){
    if(Model[1]=='or'){
      Net[!(Net[,1]%in%Node|Net[,2]%in%Node),]
    }else{
      Net[!(Net[,1]%in%Node&Net[,2]%in%Node),]
    }

  }else{
    if(Model[1]=='or'){
      Net[Net[,1]%in%Node|Net[,2]%in%Node,]
    }else{
      Net[Net[,1]%in%Node&Net[,2]%in%Node,]
    }
  }
}


#' theme_cowplot
#'
#' @param font_size
#' @param font_family
#' @param line_size
#' @param rel_small
#' @param rel_tiny
#' @param rel_large
#'
#' @return
#' @export
#'
#' @examples
theme_cowplot_i <- function (font_size = 14, font_family = "", line_size = 0.5,
                             rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)
{
  half_line <- font_size/2
  small_size <- rel_small * font_size
  theme_grey(base_size = font_size, base_family = font_family) %+replace%
    theme(line = element_line(color = "black", linewidth = line_size,
                              linetype = 1, lineend = "butt"), rect = element_rect(fill = NA,
                                                                                   color = NA, size = line_size, linetype = 1), text = element_text(family = font_family,
                                                                                                                                                    face = "plain", color = "black", size = font_size,
                                                                                                                                                    hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                                                                                                                                    margin = margin(), debug = FALSE), axis.line = element_line(color = "black",
                                                                                                                                                                                                                linewidth = line_size, lineend = "square"), axis.line.x = NULL,
          axis.line.y = NULL, axis.text = element_text(color = "black",
                                                       size = small_size), axis.text.x = element_text(margin = margin(t = small_size/4),
                                                                                                      vjust = 1), axis.text.x.top = element_text(margin = margin(b = small_size/4),
                                                                                                                                                 vjust = 0), axis.text.y = element_text(margin = margin(r = small_size/4),
                                                                                                                                                                                        hjust = 1), axis.text.y.right = element_text(margin = margin(l = small_size/4),
                                                                                                                                                                                                                                     hjust = 0), axis.ticks = element_line(color = "black",
                                                                                                                                                                                                                                                                           size = line_size), axis.ticks.length = unit(half_line/2,
                                                                                                                                                                                                                                                                                                                       "pt"), axis.title.x = element_text(margin = margin(t = half_line/2),
                                                                                                                                                                                                                                                                                                                                                          vjust = 1), axis.title.x.top = element_text(margin = margin(b = half_line/2),
                                                                                                                                                                                                                                                                                                                                                                                                      vjust = 0), axis.title.y = element_text(angle = 90,
                                                                                                                                                                                                                                                                                                                                                                                                                                              margin = margin(r = half_line/2), vjust = 1),
          axis.title.y.right = element_text(angle = -90, margin = margin(l = half_line/2),
                                            vjust = 0), legend.background = element_blank(),
          legend.spacing = unit(font_size, "pt"), legend.spacing.x = NULL,
          legend.spacing.y = NULL, legend.margin = margin(0,
                                                          0, 0, 0), legend.key = element_blank(), legend.key.size = unit(1.1 *
                                                                                                                           font_size, "pt"), legend.key.height = NULL,
          legend.key.width = NULL, legend.text = element_text(size = rel(rel_small)),
          legend.text.align = NULL, legend.title = element_text(hjust = 0),
          legend.title.align = NULL, legend.position = "right",
          legend.direction = NULL, legend.justification = c("left",
                                                            "center"), legend.box = NULL, legend.box.margin = margin(0,
                                                                                                                     0, 0, 0), legend.box.background = element_blank(),
          legend.box.spacing = unit(font_size, "pt"),
          panel.background = element_blank(), panel.border = element_blank(),
          panel.grid = element_blank(), panel.grid.major = NULL,
          panel.grid.minor = NULL, panel.grid.major.x = NULL,
          panel.grid.major.y = NULL, panel.grid.minor.x = NULL,
          panel.grid.minor.y = NULL, panel.spacing = unit(half_line,
                                                          "pt"), panel.spacing.x = NULL, panel.spacing.y = NULL,
          panel.ontop = FALSE, strip.background = element_rect(fill = "grey80"),
          strip.text = element_text(size = rel(rel_small),
                                    margin = margin(half_line/2, half_line/2, half_line/2,
                                                    half_line/2)), strip.text.x = NULL, strip.text.y = element_text(angle = -90),
          strip.placement = "inside", strip.placement.x = NULL,
          strip.placement.y = NULL, strip.switch.pad.grid = unit(half_line/2,
                                                                 "pt"), strip.switch.pad.wrap = unit(half_line/2,
                                                                                                     "pt"), plot.background = element_blank(),
          plot.title = element_text(face = "bold", size = rel(rel_large),
                                    hjust = 0, vjust = 1, margin = margin(b = half_line)),
          plot.subtitle = element_text(size = rel(rel_small),
                                       hjust = 0, vjust = 1, margin = margin(b = half_line)),
          plot.caption = element_text(size = rel(rel_tiny),
                                      hjust = 1, vjust = 1, margin = margin(t = half_line)),
          plot.tag = element_text(face = "bold", hjust = 0,
                                  vjust = 0.7), plot.tag.position = c(0, 1), plot.margin = margin(half_line,
                                                                                                  half_line, half_line, half_line), complete = TRUE)
}


#' AddBox
#' add box to ggplot
#' @param size
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
AddBox <- function(size=1.1,...){
  Boxed <-  theme(axis.line = element_blank(),panel.background = element_rect(fill = "white", colour = NA),
                  panel.border = element_rect(fill = NA, colour = "black",size=size,...))
  return(Boxed)
}


List2dataFrame <- function (l) {
  rId <- unlist(l, use.names = FALSE)
  lId <- rep(names(l), sapply(l, length))
  d<-data.frame(C1 = lId,C2=rId,stringsAsFactors = F)
  return(d)
}


getCloseseData <- function(data,query,returnIndex=FALSE){
  if(returnIndex){
    FNN::get.knnx(data,query,k=1)$nn.index[,1]
  }else{
    data[FNN::get.knnx(data,query,k=1)$nn.index[,1]]
  }
}


#' DensityPlotDF_withPoint
#'
#' @param x
#' @param y
#' @param colors
#' @param ncolor
#' @param pt.alpha
#' @param fill.alpha
#' @param pt.size
#' @param pt.color
#' @param fill.alpha.cut
#' @param showAllpoint
#'
#' @return
#' @export
#'
#' @examples
DensityPlotDF_withPoint <- function(x,y,colors=c('green','yellow','red','darkred'),
         ncolor=20,pt.alpha=0.5,fill.alpha=0.5,pt.size=0.5,
         pt.color='grey10',fill.alpha.cut=0.15,showAllpoint=TRUE){
  Pdata = data.frame(x=x,y=y)
  # ggplot()+geom_point(data =Pdata,mapping = aes(x,y),alpha=pt.alpha,color=pt.color,size=pt.size)+
  #   geom_density2d_filled(data =Pdata,mapping = aes(x,y), n=100,alpha=fill.alpha)+
  #   # geom_raster(mapping = aes(x = x ,y =  y,fill = ..density..), data = Pdata,interpolate = TRUE)+
  #   scale_fill_manual(values = colorRampPalette(c('transparent',colors), interpolate = c("linear"))(ncolor))+
  #   theme_cowplot()+AddBox()+AddGrid()


  if(showAllpoint){
    ggplot()+geom_point(data =Pdata,mapping = aes(x,y),alpha=pt.alpha,color=pt.color,size=pt.size)+
      stat_density_2d(data =Pdata,mapping = aes(x,y,fill = after_stat(density)),
                      geom = "raster",
                      contour = FALSE,interpolate =TRUE,alpha=fill.alpha)+scale_fill_gradientn(colours =  colorRampPalette(c('transparent',colors),
                                                                                                                           interpolate = c("linear"))(ncolor))+
      theme_cowplot_i()+AddBox()+scale_alpha_continuous(range = c(0, 1), guide = "none")
  }else{
    ggplot()+geom_point(data =Pdata,mapping = aes(x,y),alpha=pt.alpha,color=pt.color,size=pt.size)+
      stat_density_2d(data =Pdata,mapping = aes(x,y,fill = after_stat(density), alpha=ifelse( after_stat(density) < fill.alpha.cut, 0, 1)),
                      geom = "raster",
                      contour = FALSE,interpolate =TRUE)+scale_fill_gradientn(colours =  colorRampPalette(c('transparent',colors),
                                                                                                          interpolate = c("linear"))(ncolor))+
      theme_cowplot_i()+AddBox()+scale_alpha_continuous(range = c(0, 1), guide = "none")
  }


}
#' NoAxes2
#' NoAxes for ggplot
#' @param ...
#' @param keep.axis.text
#' @param keep.axis.title
#' @param keep.ticks
#'
#' @return
#' @export
#'
#' @examples
NoAxes2 <- function (..., keep.axis.text = FALSE,keep.axis.title = FALSE,  keep.ticks = FALSE) {
  blank <- element_blank()
  no.axes.theme <- theme(axis.line.x = blank, axis.line.y = blank,
                         validate = TRUE, ...)
  if (!keep.axis.text) {
    no.axes.theme <- no.axes.theme + theme(axis.text.x = blank,
                                           axis.text.y = blank,
                                           validate = TRUE, ...)
  }
  if (!keep.axis.title) {
    no.axes.theme <- no.axes.theme + theme(axis.title.x = blank, axis.title.y = blank,
                                           validate = TRUE, ...)
  }
  if (!keep.ticks) {
    no.axes.theme <- no.axes.theme + theme(axis.ticks.x = blank,
                                           axis.ticks.y = blank, validate = TRUE, ...)
  }
  return(no.axes.theme)
}




fit_poly <- function(x,y,degree = 1, raw = FALSE,...){
  InputData <- data.frame(x=as.numeric(x),y=as.numeric(y))
  modelx  <- lm(y ~ poly(x, degree = degree, raw = raw,...), data = InputData)
  InputData$y_p <- predict(modelx,data.frame(x=InputData$x))
  Linedata <- data.frame(x=seq(min(InputData$x),max(InputData$x),(max(InputData$x)-min(InputData$x))/30))
  Linedata$y <-  predict(modelx,data.frame(x=Linedata$x))

  InputData$x_ds <- getCloseseData(Linedata$x,InputData$x)
  Linedata = Linedata[Linedata$x%in%unique(InputData$x_ds),]
  Linedata$y_m <- aggregate(InputData$y,list(InputData$x_ds),mean)[,2]
  p <- ggplot()+geom_point(data = InputData,mapping = aes(x,y),size = 1,color='brown',alpha=0.1)+
    # geom_boxplot(data = InputData,mapping = aes(x_ds,y))+
    geom_line(data = Linedata,mapping = aes(x,y),size=1.5,color='grey80')+
    geom_point(data = Linedata,mapping = aes(x,y_m),size=2,color='green',alpha=1)+
    theme_cowplot_i()+AddBox()+ggtitle('y ~ x')+
    xlab('x')+ylab('y') # Size ~ mean
  return(list(predictP=modelx,p=p))
}

GeneInteraction <- function(scDNSobject,
         Nodes,
         EdgeID = NULL,
         subEdgeID = 1,fillColorData=NULL){
  Network <- scDNSobject@Network
  if (is.null(EdgeID)) {
    if (length(Nodes) == 1) {
      EdgeID = Network[, 1] %in% Nodes | Network[, 2] %in%
        Nodes
      if (sum(EdgeID) > 1) {
        print(sum(EdgeID))
        EdgeID = which(EdgeID)[subEdgeID]
      }
    }
    else {
      EdgeID = Network[, 1] %in% Nodes & Network[, 2] %in%
        Nodes
    }
  }
  if (sum(EdgeID + 0) != 1 & length(EdgeID) != 1) {
    stop("EdgeID is mutiple or dose not exit")
  }
  # subnet <- Network[EdgeID, ]
  # print(EdgeID)
  if(!is.null(fillColorData)){
    Pdata <- data.frame(G1=scDNSobject@data[Network[EdgeID,1],],
                        G2=scDNSobject@data[Network[EdgeID,2],],
                        scContubions=fillColorData,
                        label=scDNSobject@GroupLabel)
    ggplot(Pdata,aes(G1,G2,color=label,fill=scContubions))+geom_point(shape=21,alpha=0.9,size=2,stroke=0.1)+labs(x=Network[EdgeID,1],y=Network[EdgeID,2])+theme_cowplot_i()
  }else{
    Pdata <- data.frame(G1=scDNSobject@data[Network[EdgeID,1],],
                        G2=scDNSobject@data[Network[EdgeID,2],],
                        # scContubions=fillColorData,
                        label=scDNSobject@GroupLabel)
    ggplot(Pdata,aes(G1,G2,color=label))+geom_point()+labs(x=Network[EdgeID,1],y=Network[EdgeID,2])+theme_cowplot_i()
  }

  # Pdata <- data.frame(G1=ExpData['AP3B2',],G2=ExpData['FTL',],label=randLabel)
  # Pdata <- data.frame(G1=counts[CadiateNet2[7659,1],],G2=counts[CadiateNet2[7659,2],],label=randLabel)
  # ggplot(Pdata,aes(G1,G2,color=label))+geom_point()+labs(x=Network[EdgeID,1],y=Network[EdgeID,2])+theme_cowplot_i()
  # DensityPlotDF_withPoint(Pdata$G1,Pdata$G2)+labs(x=Network[EdgeID,1],y=Network[EdgeID,2])
}


replace2 <- function(x,RawData,RepData){
  x=as.character(x)
  RawData=as.character(RawData)
  RepData=as.character(RepData)
  matchID=match(x,RawData);
  x[!is.na(matchID)]=RepData[matchID[!is.na(matchID)]];
  x
}


#' theme_pretty
#'
#' @param fontsize font size. def:10
#' @param font Helvetica
#'
#' @return
#' @export
#'
#' @examples
theme_pretty <- function(fontsize = 10, font = "Helvetica"){
  nl <- theme_bw(base_size = fontsize) + theme(panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank(),
                                               axis.text = element_text(colour = "black", family = font),
                                               legend.key = element_blank(),
                                               strip.background = element_rect(colour = "black",  fill = "white"))
  return(nl)
}

#' theme_pretty_NoBox
#'
#'  theme_pretty with no box
#'
#' @param fontsize font size.def:10
#' @param font Helvetica
#'
#' @return
#' @export
#'
#' @examples
theme_pretty_NoBox <- function (fontsize = 10, font = "Helvetica") {
  nl <- theme_bw(base_size = fontsize) + theme(panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank(),
                                               axis.text = element_text(colour = "black", family = font),
                                               legend.key = element_blank(),
                                               strip.background = element_rect(colour = "black",  fill = "white"),
                                               axis.line.x.top = element_blank(),
                                               axis.line.y.right = element_blank(),
                                               rect = element_blank(),
                                               axis.line = element_line(color = "black", lineend = "square",linewidth = rel(0.5)))
  return(nl)
}

copyRowColname <- function(pasteM,copyM){
  rownames(pasteM) = rownames(copyM)
  colnames(pasteM) = colnames(copyM)
  pasteM
}

getAccByScore <- function(x){
  rank(x,ties.method = 'random')/length(x)
}
