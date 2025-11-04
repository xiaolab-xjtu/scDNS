#' densityCompare
#' visualize the joint and conditional probability densities for a single gene-pair (an "edge") across two different biological contexts,
#' @param scDNSobject scDNS object
#' @param Nodes The gene(s) of interest. If two gene names are provided (e.g., c("GeneA", "GeneB")), the function plots the edge between them. If only one gene name is provided, the function attempts to plot the edge corresponding to the index specified by subEdgeID.
#' @param EdgeID (Optional) A pre-calculated index or logical vector specifying the exact edge row in the network to plot. If NULL (default), the function determines the edge using the Nodes parameter. Must select only one edge.
#' @param interpolate Specifies whether to use interpolation when rendering the density raster plots. Default is FALSE.
#' @param filp If TRUE, it transposes the axes (flips the X and Y variables) for the Joint Probability a
#' @param subEdgeID (Only used when Nodes has a length of 1) An index to select one specific edge from the list of all edges connected to the single gene provided in Nodes. Default is 1 (the first edge found).
#' @param titleSzie The base font size to be used for plot titles and subtitles, passed to the internal plotting function FontSize. Default is 10.
#'
#' @return
#' @export
#'
#' @examples
densityCompare <- function(scDNSobject,Nodes,EdgeID=NULL,interpolate=FALSE,filp=FALSE,subEdgeID=1,titleSzie=10){
  Network <- scDNSobject@Network
  if(is.null(EdgeID)){
    if(length(Nodes)==1){
      EdgeID=Network[,1]%in%Nodes|Network[,2]%in%Nodes
      if(sum(EdgeID)>1){
        print(sum(EdgeID))
        EdgeID = which(EdgeID)[subEdgeID]
      }
    }else{
      EdgeID=Network[,1]%in%Nodes&Network[,2]%in%Nodes
    }

  }
  if(sum(EdgeID+0)!=1&length(EdgeID)!=1){
    stop('EdgeID is mutiple or dose not exit')
  }
  uniCase <- scDNSobject@uniCase
  subnet <-  Network[EdgeID,]
  ContextA_R <- getConditionalDenstiy(DS = scDNSobject@JDensity_A[EdgeID,,drop = F])
  ContextB_R <- getConditionalDenstiy(DS = scDNSobject@JDensity_B[EdgeID,,drop = F])
  if(filp){
    p1<-DensityPlot_raster(ContextA_R$DS[1,],interpolate = interpolate,transpose = T)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+NoLegend()+labs(subtitle='Joint probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextA_R$cDS_1[1,],addSideBar = T,transpose = F,interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(subtitle='Conditional probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextA_R$cDS_2[1,],addSideBar = T,transpose = T,interpolate = interpolate)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+labs(subtitle='Conditional probability')+plot_annotation('Conetxt A')+FontSize(main=titleSzie)
    p2<-DensityPlot_raster(ContextB_R$DS[1,],interpolate = interpolate,transpose = T)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+NoLegend()+labs(subtitle='Joint probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextB_R$cDS_1[1,],addSideBar = T,transpose = F,interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(subtitle='Conditional probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextB_R$cDS_2[1,],addSideBar = T,transpose = T,interpolate = interpolate)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+labs(subtitle='Conditional probability')+
      plot_annotation(title = 'Conetxt B')+FontSize(main=titleSzie)
    p1/p2+plot_layout(guides='collect')
  }else{
    p1<-DensityPlot_raster(ContextA_R$DS[1,],interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(title='Joint probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextA_R$cDS_1[1,],addSideBar = T,transpose = F,interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(title='Conditional probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextA_R$cDS_2[1,],addSideBar = T,transpose = T,interpolate = interpolate)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+labs(title='Conditional probability')+plot_annotation('Conetxt A')+FontSize(main=titleSzie)
    p2<-DensityPlot_raster(ContextB_R$DS[1,],interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(title='Joint probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextB_R$cDS_1[1,],addSideBar = T,transpose = F,interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(title='Conditional probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextB_R$cDS_2[1,],addSideBar = T,transpose = T,interpolate = interpolate)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+labs(title='Conditional probability')+plot_annotation('Conetxt B')+FontSize(main=titleSzie)
    p1/p2+plot_layout(guides='collect')
  }
  if(filp){
    p1<-DensityPlot_raster(ContextA_R$DS[1,],interpolate = interpolate,transpose = T)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+NoLegend()+labs(subtitle='Joint probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextA_R$cDS_1[1,],addSideBar = T,transpose = F,interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(subtitle='Conditional probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextA_R$cDS_2[1,],addSideBar = T,transpose = T,interpolate = interpolate)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+labs(subtitle='Conditional probability')+plot_annotation(uniCase[1])+FontSize(main=titleSzie)
    p2<-DensityPlot_raster(ContextB_R$DS[1,],interpolate = interpolate,transpose = T)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+NoLegend()+labs(subtitle='Joint probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextB_R$cDS_1[1,],addSideBar = T,transpose = F,interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(subtitle='Conditional probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextB_R$cDS_2[1,],addSideBar = T,transpose = T,interpolate = interpolate)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+labs(subtitle='Conditional probability')+
      plot_annotation(title = uniCase[2])+FontSize(main=titleSzie)
    p1/p2+plot_layout(guides='collect')
  }else{
    p1<-DensityPlot_raster(ContextA_R$DS[1,],interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(title='Joint probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextA_R$cDS_1[1,],addSideBar = T,transpose = F,interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(title='Conditional probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextA_R$cDS_2[1,],addSideBar = T,transpose = T,interpolate = interpolate)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+labs(title='Conditional probability')+plot_annotation(uniCase[1])+FontSize(main=titleSzie)
    p2<-DensityPlot_raster(ContextB_R$DS[1,],interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(title='Joint probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextB_R$cDS_1[1,],addSideBar = T,transpose = F,interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(title='Conditional probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextB_R$cDS_2[1,],addSideBar = T,transpose = T,interpolate = interpolate)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+labs(title='Conditional probability')+plot_annotation(uniCase[2])+FontSize(main=titleSzie)
    p1/p2+plot_layout(guides='collect')
  }

}




#' densityCompare2
#'
#' @param netRes
#' @param Nodes
#' @param EdgeID
#' @param interpolate
#' @param filp
#' @param subEdgeID
#' @param titleSzie
#'
#' @return
#' @keywords internal
#'
#' @examples
densityCompare2 <- function(netRes,Nodes,EdgeID=NULL,interpolate=FALSE,filp=FALSE,subEdgeID=1,titleSzie=10){
  Network <- netRes$Network
  if(is.null(EdgeID)){
    if(length(Nodes)==1){
      EdgeID=Network[,1]%in%Nodes|Network[,2]%in%Nodes
      if(sum(EdgeID)>1){
        print(sum(EdgeID))
        EdgeID = which(EdgeID)[subEdgeID]
      }
    }else{
      EdgeID=Network[,1]%in%Nodes&Network[,2]%in%Nodes
    }

  }
  if(sum(EdgeID+0)!=1&length(EdgeID)!=1){
    stop('EdgeID is mutiple or dose not exit')
  }
  uniCase <- netRes$uniCase
  subnet <-  Network[EdgeID,]
  ContextA_R <- getConditionalDenstiy(DS = netRes$ContextA_DS[EdgeID,,drop = F])
  ContextB_R <- getConditionalDenstiy(DS = netRes$ContextB_DS[EdgeID,,drop = F])
  if(filp){
    p1<-DensityPlot_raster(ContextA_R$DS[1,],interpolate = interpolate,transpose = T)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+NoLegend()+labs(subtitle='Joint probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextA_R$cDS_1[1,],addSideBar = T,transpose = F,interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(subtitle='Conditional probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextA_R$cDS_2[1,],addSideBar = T,transpose = T,interpolate = interpolate)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+labs(subtitle='Conditional probability')+plot_annotation('Conetxt A')+FontSize(main=titleSzie)
    p2<-DensityPlot_raster(ContextB_R$DS[1,],interpolate = interpolate,transpose = T)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+NoLegend()+labs(subtitle='Joint probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextB_R$cDS_1[1,],addSideBar = T,transpose = F,interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(subtitle='Conditional probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextB_R$cDS_2[1,],addSideBar = T,transpose = T,interpolate = interpolate)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+labs(subtitle='Conditional probability')+
      plot_annotation(title = 'Conetxt B')+FontSize(main=titleSzie)
    p1/p2+plot_layout(guides='collect')
  }else{
    p1<-DensityPlot_raster(ContextA_R$DS[1,],interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(title='Joint probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextA_R$cDS_1[1,],addSideBar = T,transpose = F,interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(title='Conditional probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextA_R$cDS_2[1,],addSideBar = T,transpose = T,interpolate = interpolate)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+labs(title='Conditional probability')+plot_annotation('Conetxt A')+FontSize(main=titleSzie)
    p2<-DensityPlot_raster(ContextB_R$DS[1,],interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(title='Joint probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextB_R$cDS_1[1,],addSideBar = T,transpose = F,interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(title='Conditional probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextB_R$cDS_2[1,],addSideBar = T,transpose = T,interpolate = interpolate)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+labs(title='Conditional probability')+plot_annotation('Conetxt B')+FontSize(main=titleSzie)
    p1/p2+plot_layout(guides='collect')
  }

}


#' Visualizing Conditional Density Distributions Across Biological Contexts
#'
#' visualize and compare the conditional density distributions of gene pairs (edges) between two different biological contexts,
#' likely representing a control (Context A) and a stimulation or perturbed state (Context B)
#' @param scDNSobject object of scDNS
#' @param Nodes A character vector specifying the gene(s) of interest.
#' @param topEdge The number of top-ranked edges (based on context divergence) to visualize (default is 10).
#' @param EdgeID An optional pre-computed logical vector to directly select the edges of interest in the network. If NULL (default), it's calculated using Nodes.
#' @param ncol The number of columns to use when arranging the resulting plots (default is 5).
#'
#' @return ggplot object
#' @export
#'
#' @examples
#'
densityCompare3 <- function (scDNSobject, Nodes, topEdge=10,
                             EdgeID = NULL,
                             ncol=5)
{
  Network <- scDNSobject@Network
  if (is.null(EdgeID)) {
    if (length(Nodes) == 1) {
      EdgeID = Network[, 1] %in% Nodes | Network[, 2] %in%
        Nodes

    }
    else {
      EdgeID = Network[, 1] %in% Nodes & Network[, 2] %in%
        Nodes
    }
  }
  message(scDNSobject@uniCase,'red','blue')
  uniCase <- scDNSobject@uniCase
  subnet <- Network[EdgeID, ]
  subnet$id <- 1:nrow(subnet)
  ContextA_R <- getConditionalDenstiy(DS = scDNSobject@JDensity_A[EdgeID,
                                                                  , drop = F])
  ContextB_R <- getConditionalDenstiy(DS = scDNSobject@JDensity_B[EdgeID,
                                                                  , drop = F])
  tbtheme <- gridExtra::ttheme_default( base_size = 5,
                                        core = list(bg_params = list(fill = NA, col = NA,alpha=0.5),  # 单元格背景透明
                                                    padding = unit(c(1, 1), "mm")) ,
                                        colhead = list(bg_params = list(fill = NA, col = NA,alpha=0.5),padding = unit(c(1, 1), "mm"))       # 表头背景透明
  )
  subnet$cDiv_D1.target <- subnet$cDiv_D1
  subnet$cDiv_D1.target[subnet[,2]%in%Nodes] <- subnet$cDiv_D2[subnet[,2]%in%Nodes]
  subnet <- arrange(subnet,desc(cDiv_D1.target))
  topx <- min(c(nrow(subnet),topEdge))
  ncolor=100
  plist <- vector("list", length = nrow(subnet[1:topx,]))
  for (i in 1:nrow(subnet[1:topx,])){
    idx <- subnet[i,]$id
    plist[[i]] <- local({
      if(subnet[i,1]%in%Nodes){
        CTL_Pdata <- getDenstiyPlotDREMI(z = matrix(ContextA_R$cDS_1[idx,],sqrt(length(ContextA_R$cDS_1[idx,]))))
        STIM_Pdata <- getDenstiyPlotDREMI(z = matrix(ContextB_R$cDS_1[idx,],sqrt(length(ContextA_R$cDS_1[idx,]))))
        Gene1 <- subnet[i,1]
        Gene2 <- subnet[i,2]
        MIDF <- subnet[i,str_detect(colnames(subnet),'DREMI_1')|colnames(subnet)%in%'cDiv_D1']
        colnames(MIDF) <- str_remove(colnames(MIDF),'_1')
        colnames(MIDF) <- str_remove(colnames(MIDF),'_D1')
      }else{
        CTL_Pdata <- getDenstiyPlotDREMI(z = matrix(ContextA_R$cDS_2[idx,],sqrt(length(ContextA_R$cDS_1[idx,])))%>%t())
        STIM_Pdata <- getDenstiyPlotDREMI(z = matrix(ContextB_R$cDS_2[idx,],sqrt(length(ContextA_R$cDS_1[idx,])))%>%t())
        Gene1 <- subnet[i,2]
        Gene2 <- subnet[i,1]
        MIDF <- subnet[i,str_detect(colnames(subnet),'DREMI_2')|colnames(subnet)%in%'cDiv_D2']
        colnames(MIDF) <- str_remove(colnames(MIDF),'_2')
        colnames(MIDF) <- str_remove(colnames(MIDF),'_D2')
      }
      MIDF <- round(MIDF,2)
      df <- data.frame(x = 0, y = 19, tb = I(list(MIDF)))
      CTL_fit <- fit_conditional_curve(peaks = CTL_Pdata$Pdata_line-1)
      r2_v <- sapply(CTL_fit$models, function(x)x$r2)
      others <- r2_v[names(r2_v) != "loess"]

      if (max(others) > 0.8&length(others)>1) {
        CTL_best_name <- names(others)[which.max(others)]
      } else {
        CTL_best_name <- "loess"
      }

      STIM_fit <- fit_conditional_curve(peaks = STIM_Pdata$Pdata_line-1)
      r2_v <- sapply(STIM_fit$models, function(x)x$r2)
      others <- r2_v[names(r2_v) != "loess"]

      if (max(others) > 0.8&length(others)>1) {
        STIM_best_name <- names(others)[which.max(others)]
      } else {
        STIM_best_name <- "loess"
      }

      p <- ggplot()+geom_raster(data = CTL_Pdata$Pdata,mapping = aes(gridx-1,gridy-1,fill = z),interpolate = T,alpha=0.5)+
        scale_fill_gradientn(colours = c('transparent','red'))+

        stat_function(fun = function(x){predict(CTL_fit$models[[CTL_best_name]]$fit,newdata = data.frame(x=x))},
                      color = "darkred", size = 1,xlim = c(0, 19)) +
        ggnewscale::new_scale_fill()+
        geom_raster(data = STIM_Pdata$Pdata,mapping = aes(gridx-1,gridy-1,fill = z),interpolate = T,alpha=0.5)+
        scale_fill_gradientn(colours = c('transparent','darkblue'))+
        stat_function(fun = function(x){predict(STIM_fit$models[[STIM_best_name]]$fit,newdata = data.frame(x=x))},
                      color = "darkblue", size = 1,xlim = c(0, 19)) +theme_pretty()+labs(x=Gene1,y=Gene2)+NoLegend()+NoAxes2(keep.axis.title = T)+
        ggpmisc::geom_table(data = df,
                   aes(x = x, y = y, label = tb),size =2,table.theme = tbtheme)
      print(i)
      p
    })

  }
  plist_patch <- wrap_plots(plist,ncol = ncol)
  plist_patch
}

#' Volcano Plot of Log2 Fold-Change versus scDNS Z-score
#'
#' plot_diffFC_scDNS_Zscore is designed to generate a Volcano-like plot that visualizes the relationship between gene expression fold-change (Log2 FC) and a network stability score (Z-score) derived from single-cell Differential Network Stability (scDNS) analysis.
#' @param Zscores A data frame or object containing the gene-level network stability Z-scores (calculated by scDNS)
#' @param scDNSob scDNS object
#' @param sob Seurat Object This is used to calculate Fold Change if it's not already present in Zscores.
#' @param col_Zs The name of the column in the Zscores data frame that contains the Z-scores (network stability metric). Default is 'Zscores.ZsPlus'
#' @param col_log2Fc The name of the column in the Zscores data frame that contains the avg_log2FC values. Default is 'avg_log2FC'.
#' @param fc_th The threshold for avg_log2FC (fold-change).  Default is 0.25.
#' @param Zs_th The threshold for the Z-score (network stability score). Genes with a Z-score greater than this value are considered to have significantly altered network stability. Default is 3.
#' @param TextSize The font size for the gene labels plotted using geom_text_repel. Default is 2.5.
#' @param pt.size A vector of two values specifying the point sizes: pt.size 1 for unlabeled/other genes and pt.size 2 for highlighted/labeled genes. Default is c(0.25, 0.5).
#' @param group.by The metadata column in the sob (Seurat object) to use for grouping and calculating fold-change (e.g., condition).
#' @param ident.up The identity class (group) to compare against the rest
#' @param interstingGene A vector of gene names to be explicitly labeled on the plot. These genes will be highlighted regardless of their score. Default is NULL
#' @param highlightGene A vector of gene names to be explicitly highlighted (enlarged point size) on the plot. Defaults to the value of interstingGene
#' @param TopGene The number of top genes based on the Zscores column (col_Zs) to be automatically labeled on the plot. Default is NULL.
#' @param darkOtherGene If TRUE, all genes that are not labeled or highlighted (label is NA and not a highlightGene).
#'
#' @return
#' @export
#'
#' @examples
plot_diffFC_scDNS_Zscore <- function(Zscores,
                                     scDNSob,
                                     sob,
                                     col_Zs=c('Zscores.ZsPlus'),
                                     col_log2Fc='avg_log2FC',
                                     fc_th=0.25,
                                     Zs_th=3,
                                     TextSize=2.5,
                                     pt.size=c(0.25,0.5),
                                     group.by,
                                     ident.up,
                                     interstingGene=NULL,
                                     highlightGene=interstingGene,
                                     TopGene=NULL,
                                     darkOtherGene=FALSE){
  if(is.null(Zscores)){
    Zscores <- Zscores@Zscore
  }
  if(!col_log2Fc%in%colnames(Zscores)){
    geneFC <- FoldChange(sob,group.by = group.by,ident.1 = ident.up)
    Zscores$avg_log2FC <- geneFC[Zscores$Gene,]$avg_log2FC
    message('claculate avg_log2FC')
  }else{
    message('use pre-calculated avg_log2FC')
  }
  Zscores$Zscore <- Zscores[,col_Zs]
  Zscores$avg_log2FC <- Zscores[,col_log2Fc]
  Zscores$label=NA
  if(is.null(interstingGene)&!is.null(TopGene)){
    # Zscores$label=NA
    # Zscores$label[Zscores$Gene%in%interstingGene] <- Zscores$Gene[Zscores$Gene%in%interstingGene]
    Zscores$label[top_n_vector(Zscores[,col_Zs],TopGene,returnLogical = T)] <- Zscores$Gene[top_n_vector(Zscores[,col_Zs],TopGene,returnLogical = T)]
  }
  if(!is.null(interstingGene)&is.null(TopGene)){
    # Zscores$label=NA
    Zscores$label[Zscores$Gene%in%interstingGene] <- Zscores$Gene[Zscores$Gene%in%interstingGene]
  }
  if(!is.null(interstingGene)&!is.null(TopGene)){
    # Zscores$label=NA
    Zscores$label[Zscores$Gene%in%interstingGene] <- Zscores$Gene[Zscores$Gene%in%interstingGene]
    Zscores$label[top_n_vector(Zscores[,col_Zs],TopGene,returnLogical = T)] <- Zscores$Gene[top_n_vector(Zscores[,col_Zs],TopGene,returnLogical = T)]
  }

  Zscores$label0 <- pt.size[1]
  Zscores$label0[!is.na(Zscores$label)] <- pt.size[2]
  Zscores$label0[Zscores$Gene%in%highlightGene] <- pt.size[2]
  Zscores <- arrange(Zscores,label0)

  Zscores <- Zscores %>%
    mutate(color = case_when(
      Zscore <= Zs_th            ~ "grey90",
      avg_log2FC >  fc_th             ~ "#E23E3E",
      avg_log2FC < -fc_th             ~ "#2E91E5",
      TRUE                        ~ "#F7B731"
    ))
  if(sum(is.na(Zscores$label))==nrow(Zscores)){
    ggplot(Zscores,aes(avg_log2FC,Zscore,color=color))+geom_point(size=pt.size[2])+theme_pretty(12)+
      geom_hline(yintercept = Zs_th,lty=2,color='grey')+geom_vline(xintercept = c(-fc_th,fc_th),lty=2,color='grey')+
      scale_color_identity()+labs(x='Log2 FC',y='Zscore')
  }else{
    if(darkOtherGene){
      Zscores$color[is.na(Zscores$label)&!Zscores$color=='grey90'] <- alpha(Zscores$color[is.na(Zscores$label)&!Zscores$color=='grey90'],0.5)
      Zscores$color[!is.na(Zscores$label)|Zscores$Gene%in%highlightGene] <- 'black'
      ggplot(Zscores,aes(avg_log2FC,Zscore,color=color,size=label0))+geom_point()+theme_pretty(12)+
        geom_hline(yintercept = Zs_th,lty=2,color='grey')+geom_vline(xintercept = c(-fc_th,fc_th),lty=2,color='grey')+
        geom_text_repel(aes(label=label),color='black', min.segment.length = 0,      # Draw leader lines for all points.
                        force_pull = 2,max.overlaps = 100,size=TextSize)+scale_color_identity()+labs(x='Log2 FC',y='Zscore')+scale_size_identity()
    }else{
      ggplot(Zscores,aes(avg_log2FC,Zscore,color=color))+geom_point(size=0.5)+theme_pretty(12)+
        geom_hline(yintercept = Zs_th,lty=2,color='grey')+geom_vline(xintercept = c(-fc_th,fc_th),lty=2,color='grey')+
        geom_text_repel(aes(label=label),color='black', min.segment.length = 0,      # Draw leader lines for all points.
                        force_pull = 2,max.overlaps = 100,size=TextSize)+scale_color_identity()+labs(x='Log2 FC',y='Zscore')
    }

  }

}


#' fit_conditional_curve
#'
#' @param peaks
#'
#' @return
#' @keywords internal
#'
#' @examples
fit_conditional_curve <- function(peaks) {
  library(dplyr)

  # Step 1: peak value of density
  # peaks <- df %>%
  #   group_by(x) %>%
  #   slice_max(order_by = z, n = 1, with_ties = FALSE) %>%
  #   ungroup()

  # R2 functions
  calc_r2 <- function(obs, pred) {
    ss_res <- sum((obs - pred)^2)
    ss_tot <- sum((obs - mean(obs))^2)
    return(1 - ss_res/ss_tot)
  }

  results <- list()

  # ---- Linear ----
  fit_linear <- lm(y ~ x, data = peaks)
  pred_linear <- predict(fit_linear, newdata = peaks)
  results$linear <- list(
    fit = fit_linear,
    r2 = calc_r2(peaks$y, pred_linear),
    predict = function(newx) predict(fit_linear, newdata = data.frame(x=newx))
  )

  # ---- Sigmoid ----
  try({
    # start_vals <- list(A = diff(range(peaks$y)),
    #                    k = 2.5,
    #                    x0 = median(peaks$x),
    #                    B = min(peaks$y))
    start_vals <- estimate_sigmoid_start(peaks)
    fit_sigmoid <- minpack.lm::nlsLM(y ~ A / (1 + exp(-k * (x - x0))) + B,
                                     data = peaks,
                                     start = start_vals,
                                     control = list(maxiter = 1000))
    pred_sigmoid <- predict(fit_sigmoid, newdata = peaks)
    results$sigmoid <- list(
      fit = fit_sigmoid,
      r2 = calc_r2(peaks$y, pred_sigmoid),
      predict = function(newx) predict(fit_sigmoid, newdata = data.frame(x=newx))
    )
  }, silent = TRUE)

  # ---- Double-sigmoid ----
  try({
    start_vals <- list(A1 = diff(range(peaks$y))/2,
                       k1 = 1,
                       x01 = quantile(peaks$x, 0.3),
                       A2 = diff(range(peaks$y))/2,
                       k2 = 1,
                       x02 = quantile(peaks$x, 0.7),
                       B = min(peaks$y))
    # start_vals <- estimate_double_sigmoid_start(peaks)
    fit_dsig <- minpack.lm::nlsLM(y ~ A1 / (1 + exp(-k1 * (x - x01))) +
                                    A2 / (1 + exp(-k2 * (x - x02))) + B,
                                  data = peaks,
                                  start = start_vals,
                                  control = list(maxiter = 1000))
    pred_dsig <- predict(fit_dsig, newdata = peaks)
    results$double_sigmoid <- list(
      fit = fit_dsig,
      r2 = calc_r2(peaks$y, pred_dsig),
      predict = function(newx) predict(fit_dsig, newdata = data.frame(x=newx))
    )
  }, silent = TRUE)

  # ---- LOESS ----
  fit_loess <- loess(y ~ x, data = peaks, span = 0.3)
  pred_loess <- predict(fit_loess, newdata = peaks)
  results$loess <- list(
    fit = fit_loess,
    r2 = calc_r2(peaks$y, pred_loess),
    predict = function(newx) predict(fit_loess, newdata = data.frame(x=newx))
  )

  # ---- find the best model ----
  r2_values <- sapply(results, function(m) m$r2)
  best_model <- names(which.max(r2_values))

  return(list(
    peaks = peaks,
    models = results,
    best_model = best_model,
    best_r2 = max(r2_values)
  ))
}

#' estimate_sigmoid_start
#'
#' @param peaks
#'
#' @return
#' @keywords internal
#'
#' @examples
estimate_sigmoid_start <- function(peaks){
  y_min <- min(peaks$y)
  y_max <- max(peaks$y)
  A <- y_max - y_min
  B <- y_min

  #
  half_y <- B + A/2
  x0 <- peaks$x[which.min(abs(peaks$y - half_y))]

  # estimation of k
  dy <- diff(peaks$y)
  dx <- diff(peaks$x)
  slopes <- dy/dx
  k <- max(slopes) / A  # 归一化斜率

  # avoid k=0
  if(k == 0) k <- 0.1

  list(A=A, B=B, x0=x0, k=k)
}

#' estimate_double_sigmoid_start
#'
#' @param peaks
#'
#' @return
#' @keywords internal
#'
#' @examples
estimate_double_sigmoid_start <- function(peaks){
  x <- peaks$x
  y <- peaks$y

  B <- min(y, na.rm=TRUE)

  dy <- diff(y)
  pos_peaks <- which(dy > 0)

  if(length(pos_peaks) < 2){
    x01 <- x[which.min(abs(y - (B + (max(y, na.rm=TRUE)-B)/2)))]
    x02 <- x01 + 1
  } else {
    x01 <- x[pos_peaks[1]]
    x02 <- x[pos_peaks[ceiling(length(pos_peaks)/2)]]
  }

  #
  x01 <- min(max(x01, min(x)), max(x))
  x02 <- min(max(x02, min(x)), max(x))

  A1 <- y[which.min(abs(x - x01))] - B
  A2 <- y[which.min(abs(x - x02))] - B
  A1 <- ifelse(is.finite(A1) & A1 != 0, A1, 1e-3)
  A2 <- ifelse(is.finite(A2) & A2 != 0, A2, 1e-3)

  k1 <- max(dy, na.rm=TRUE)/A1
  k2 <- max(dy, na.rm=TRUE)/A2
  k1 <- ifelse(is.finite(k1) & k1 != 0, k1, 0.01)
  k2 <- ifelse(is.finite(k2) & k2 != 0, k2, 0.01)

  list(A1=A1, k1=k1, x01=x01, A2=A2, k2=k2, x02=x02, B=B)

  list(A1=A1, k1=k1, x01=x01, A2=A2, k2=k2, x02=x02, B=B)
}


getDenstiyPlotDREMI <- function(z){
  grid = expand.grid(1:nrow(z),1:ncol(z))
  Pdata = data.frame(gridx = grid[,1],gridy=grid[,2],z = matrix(z,ncol = 1))
  max_d = apply(z, 1,function(x)which(x==max(x)))


  # z2 = (z[1:19,1:19]+z[2:20,2:20])/2
  # grid2 = expand.grid(1:nrow(z2)+0.5,1:ncol(z2)+0.5)
  # Pdata2 = data.frame(gridx = grid2[,1],gridy=grid2[,2],z = matrix(z2,ncol = 1))
  # ggplot(Pdata,aes(gridx,gridy,z = z))+geom_contour_filled(bins =ncolor)+scale_fill_manual(
  #   values = colorRampPalette(colors, interpolate = c("linear"))(ncolor))+
  #   theme_cowplot()+xlab('G1')+ylab('G2')+geom_point()+NoLegend()+AddBox()
  if(is.matrix(max_d)){
    Pdata_line=data.frame(x=rep(nrow(z)/2,nrow(z)),y=rep(nrow(z)/2,nrow(z)),stringsAsFactors = F)
  }else{
    for(i in 1:(length(max_d)-2)){
      allcombinations = expand.grid(max_d[[i]],max_d[[i+1]],max_d[[i+2]])
      rowsd = matrixStats::rowSds(allcombinations%>%as.matrix())
      minID = which(rowsd==min(rowsd))
      minID = minID[which.min(abs((allcombinations[,3]-allcombinations[,2]))[minID])]
      max_d[[i]]=allcombinations[minID,1]
      max_d[[i+1]]=allcombinations[minID,2]
      max_d[[i+2]]=allcombinations[minID,3]
    }
    Pdata_line = data.frame(x=1:nrow(z),y=unlist(max_d),stringsAsFactors = F)
  }
  list(Pdata=Pdata,Pdata_line=Pdata_line)
  # ggplot()+geom_contour_filled(data = Pdata,mapping = aes(gridx,gridy,z = z),bins =ncolor)+
  #   scale_fill_manual(values = colorRampPalette(colors, interpolate = c("linear"))(ncolor))+
  #   geom_point(data = Pdata,mapping = aes(gridx,gridy),fill='grey50',color='transparent',alpha=0.1,shape=21)+
  #   geom_line(data = Pdata_line,mapping = aes(x,y),size=2,color='white')+
  #   theme_cowplot()+xlab('G1')+ylab('G2')+NoLegend()+AddBox()
  #
  #
  # Pdata3 = rbind(Pdata,Pdata2)
  # ggplot(Pdata3,aes(gridx,gridy,z = z))+geom_contour_filled(bins =ncolor)+scale_fill_manual(
  #   values = colorRampPalette(colors, interpolate = c("linear"))(ncolor))+
  #   theme_cowplot()+xlab('G1')+ylab('G2')+geom_point()
}
