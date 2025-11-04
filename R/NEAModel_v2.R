#' scDNS_2_creatNEAModel_v2
#'
#' Creates a null model for network enrichment analysis by random sampling of genes and networks.
#'
#' @param scDNSobject scDNSobject
#' @param n.dropGene Integer, number of cells sampled for model building (def:3000)
#' @param n.randNet Integer, number of cells sampled for model building (def:20000)
#' @param sdBias Numeric value, bias coefficient used to penalize low degree genes(>=1) (def:1.1)
#' @param GroupLabel
#' @param repTime
#'
#' @return
#' @export
#'
#' @examples
scDNS_2_creatNEAModel_v2 <- function (scDNSobject, n.dropGene = NULL, n.randNet = NULL,repTime=5,
                                      sdBias = NULL,GroupLabel=NULL)
{
  if (!is.null(n.dropGene)) {
    scDNSobject@NEA.Parameters$n.dropGene = n.dropGene
  }else{
    n.dropGene <- scDNSobject@NEA.Parameters$n.dropGene
  }
  if (!is.null(n.randNet)) {
    scDNSobject@NEA.Parameters$n.randNet = n.randNet
  }
  if (!is.null(sdBias)) {
    scDNSobject@NEA.Parameters$sdBias = sdBias
  }
  if(is.null(GroupLabel)){
    GroupLabel <- scDNSobject@GroupLabel
  }

  #
  uniCase = unique(GroupLabel)
  ExpData <- scDNSobject@data
  RandGene <- sample(rownames(ExpData), pmin(n.dropGene, nrow(ExpData)))
  ExpData <- ExpData[RandGene,]

  counts = scDNSobject@counts[RandGene,]
  dropoutMatrix = Matrix::tcrossprod(counts != 0)
  dropoutMatrix_log = round(log10(dropoutMatrix + 1), 1)

  Likelihood <- scDNSobject@GeneVariability
  #
  message("creat random net")
  CandidateNet_samll = NULL

  if (length(names(Likelihood[rownames(ExpData)][Likelihood[rownames(ExpData)] >
                                                 0.95])) >= 100) {
    sn <- sample(names(Likelihood[rownames(ExpData)][Likelihood[rownames(ExpData)] >
                                                       0.95]), 100)
  }  else {
    sn <- sample(names(Likelihood[rownames(ExpData)]), pmin(100,
                                                            length(names(Likelihood[rownames(ExpData)]))), prob = Likelihood[rownames(ExpData)])
  }
  CandidateNet_samll = lapply(1:length(sn), function(x) sample(rownames(ExpData),
                                                               pmin(300, nrow(ExpData))))
  names(CandidateNet_samll) = sn
  CandidateNet_samll = List2dataFrame(CandidateNet_samll)
  colnames(CandidateNet_samll) = c("source", "target")
  CandidateNet_samll$npoint_log = dropoutMatrix_log[sub2ind(match(CandidateNet_samll[,
                                                                                     1], rownames(dropoutMatrix_log)), match(CandidateNet_samll[,
                                                                                                                                                2], rownames(dropoutMatrix_log)), nrow = nrow(dropoutMatrix_log),
                                                            ncol = ncol(dropoutMatrix_log))]
  CandidateNet_samll$npoint_ds <-  CandidateNet_samll$npoint_log
  CandidateNet_samll$LR = pmin(Likelihood[CandidateNet_samll[,1]],
                               Likelihood[CandidateNet_samll[,2]])
  message("rand netowrk")
  # # rand netowrk
  NEAModel_randNet <- creatNEAModel_test(ExpData = ExpData, CandidateNet_samll=CandidateNet_samll,
                                         repTime = repTime,
                                         k = scDNSobject@Div.Parameters$k,
                                         GroupLabel = GroupLabel,
                                         n.grid = scDNSobject@Div.Parameters$n.grid, n.coarse = scDNSobject@Div.Parameters$n.coarse,
                                         noiseSd = scDNSobject@Div.Parameters$noiseSd, Div_weight = scDNSobject@Div.Parameters$Div_weight,
                                         CoarseGrain = scDNSobject@Div.Parameters$CoarseGrain,
                                         loop.size = scDNSobject@Div.Parameters$loop.size, parallel.sz = scDNSobject@Div.Parameters$parallel.sz,
                                         verbose = scDNSobject@Div.Parameters$verbose, exclude.zero = scDNSobject@Div.Parameters$exclude.zero,
                                         NoiseRemove = scDNSobject@Div.Parameters$NoiseRemove,
                                         divergence = scDNSobject@Div.Parameters$divergence,
                                         ds.method = scDNSobject@Div.Parameters$ds.method, h = scDNSobject@Div.Parameters$h,
                                         sdBias = scDNSobject@NEA.Parameters$sdBias, rb.jsd = scDNSobject@Div.Parameters$rb.jsd,
                                         parllelModel = scDNSobject@Div.Parameters$parllelModel)
  # shuffleLabel
  message("shuffleLabel")
  NEAModel_shuffleLabel <- creatNEAModel_test( ExpData = ExpData, CandidateNet_samll=CandidateNet_samll,
                                               repTime = repTime,
                                               k = scDNSobject@Div.Parameters$k,
                                               GroupLabel = sample(GroupLabel), ###########################❤
                                               n.grid = scDNSobject@Div.Parameters$n.grid, n.coarse = scDNSobject@Div.Parameters$n.coarse,
                                               noiseSd = scDNSobject@Div.Parameters$noiseSd, Div_weight = scDNSobject@Div.Parameters$Div_weight,
                                               CoarseGrain = scDNSobject@Div.Parameters$CoarseGrain,
                                               loop.size = scDNSobject@Div.Parameters$loop.size, parallel.sz = scDNSobject@Div.Parameters$parallel.sz,
                                               verbose = scDNSobject@Div.Parameters$verbose, exclude.zero = scDNSobject@Div.Parameters$exclude.zero,
                                               NoiseRemove = scDNSobject@Div.Parameters$NoiseRemove,
                                               divergence = scDNSobject@Div.Parameters$divergence,
                                               ds.method = scDNSobject@Div.Parameters$ds.method, h = scDNSobject@Div.Parameters$h,
                                               sdBias = scDNSobject@NEA.Parameters$sdBias, rb.jsd = scDNSobject@Div.Parameters$rb.jsd,
                                               parllelModel = scDNSobject@Div.Parameters$parllelModel)

  # shuffle distribution
  message("shuffle distribution")
  shuffleData <- shuffle_rows_by_group(counts,data = ExpData,GroupLabel = GroupLabel)

  dropoutMatrix = Matrix::tcrossprod(shuffleData$counts != 0)
  dropoutMatrix_log = round(log10(dropoutMatrix + 1), 1)
  CandidateNet_samll$npoint_log = dropoutMatrix_log[sub2ind(match(CandidateNet_samll[,
                                                                                     1], rownames(dropoutMatrix_log)), match(CandidateNet_samll[,
                                                                                                                                                2], rownames(dropoutMatrix_log)), nrow = nrow(dropoutMatrix_log),
                                                            ncol = ncol(dropoutMatrix_log))]

  NEAModel_RandDistrubution <- creatNEAModel_test( ExpData = shuffleData$data,
                                                   CandidateNet_samll=CandidateNet_samll,
                                                   repTime = repTime,
                                                   k = scDNSobject@Div.Parameters$k,
                                                   GroupLabel = GroupLabel,
                                                   n.grid = scDNSobject@Div.Parameters$n.grid, n.coarse = scDNSobject@Div.Parameters$n.coarse,
                                                   noiseSd = scDNSobject@Div.Parameters$noiseSd, Div_weight = scDNSobject@Div.Parameters$Div_weight,
                                                   CoarseGrain = scDNSobject@Div.Parameters$CoarseGrain,
                                                   loop.size = scDNSobject@Div.Parameters$loop.size, parallel.sz = scDNSobject@Div.Parameters$parallel.sz,
                                                   verbose = scDNSobject@Div.Parameters$verbose, exclude.zero = scDNSobject@Div.Parameters$exclude.zero,
                                                   NoiseRemove = scDNSobject@Div.Parameters$NoiseRemove,
                                                   divergence = scDNSobject@Div.Parameters$divergence,
                                                   ds.method = scDNSobject@Div.Parameters$ds.method, h = scDNSobject@Div.Parameters$h,
                                                   sdBias = scDNSobject@NEA.Parameters$sdBias, rb.jsd = scDNSobject@Div.Parameters$rb.jsd,
                                                   parllelModel = scDNSobject@Div.Parameters$parllelModel)
  NEAModel <- list(randNet=NEAModel_randNet,shuffleLabel=NEAModel_shuffleLabel,RandDistrubution=NEAModel_RandDistrubution)
  scDNSobject@NEAModel <- NEAModel
  scDNSobject

}

shuffle_rows_by_group <- function(counts, data, GroupLabel) {
  if (!all(dim(counts) == dim(data))) {
    stop("The dimensions of counts and data must be identical.")
  }
  if (length(GroupLabel) != ncol(counts)) {
    stop("The length of GroupLabel must match the number of columns.")
  }

  nr <- nrow(counts)
  groups <- unique(GroupLabel)

  indexMatrix <- matrix(1:length(counts),nrow = nrow(counts))
  indexMatrix[,GroupLabel==groups[1]] <- apply(indexMatrix[,GroupLabel==groups[1]],1,sample)%>%t()
  indexMatrix[,GroupLabel==groups[2]] <- apply(indexMatrix[,GroupLabel==groups[2]],1,sample)%>%t()
  counts_shuffled <- matrix(counts[as.vector(indexMatrix)],nrow =nrow(counts) )
  data_shuffled <- matrix(data[as.vector(indexMatrix)],nrow =nrow(counts) )
  counts_shuffled <- copyRowColname(pasteM = counts_shuffled,copyM = counts)
  data_shuffled <- copyRowColname(pasteM = data_shuffled,copyM = data)
  return(list(counts = counts_shuffled, data = data_shuffled))
}

creatNEAModel_test <- function (ExpData = NULL, CandidateNet_samll=NULL,
                                k = 10, Div_weight = sqrt(as.vector(matrix(1:n.coarse,
                                                                           nrow = n.coarse, ncol = n.coarse) * t(matrix(1:n.coarse,
                                                                                                                        nrow = n.coarse, ncol = n.coarse)))),
                                GroupLabel = NULL,
                                n.grid = 60, n.coarse = 20, CoarseGrain = TRUE, loop.size = 3000,
                                parallel.sz = 10, noiseSd = 0.01, verbose = TRUE, exclude.zero = FALSE,
                                NoiseRemove = TRUE, divergence = "jsd", ds.method = "knn",repTime=5,
                                h = 3, sdBias = 1.1, rb.jsd = FALSE, parllelModel = c("foreach",
                                                                                      "bplapply")[1])
{
  RawKLD_R_samll = getKLD_cKLDnetwork(ExpData = ExpData, Network = CandidateNet_samll,
                                              k = k, GroupLabel = GroupLabel, n.grid = n.grid, n.coarse = n.coarse,
                                              loop.size = nrow(CandidateNet_samll)/parallel.sz, parallel.sz = parallel.sz,
                                              verbose = verbose, exclude.zero = exclude.zero, NoiseRemove = NoiseRemove,
                                              CoarseGrain = CoarseGrain, returnCDS = FALSE, Div_weight = Div_weight,
                                              divergence = divergence, ds.method = ds.method, h = h,
                                              parllelModel = parllelModel, rb.jsd = rb.jsd)
  NEAModel <- creatNEAModel_robustness(Network=RawKLD_R_samll$Network,repTime=repTime,sdBias = sdBias)
  NEAModel
}

creatNEAModel_robustness<- function(Network,repTime=5,sdBias=1.1){




  CandidateNet_samll = Network

  # correlation between dropout and network divergence
  CadiateNet2 <- CandidateNet_samll
  CadiateNet2$id = 1:nrow(CadiateNet2)
  CadiateNet2 <- CadiateNet2[as.character(CadiateNet2$npoint_ds) %in%
                               names(table(CadiateNet2$npoint_ds))[table(CadiateNet2$npoint_ds) >
                                                                     30], ]
  RawDistrbution.Div <- NULL
  RawDistrbution.cDiv <- NULL
  for (i in unique(CadiateNet2$npoint_ds) %>% sort()) {
    idx = CadiateNet2$npoint_ds == i
    RawDistrbution.Div <- c(RawDistrbution.Div, list(fit_selm_distrubution(CadiateNet2$Div[idx])))
    RawDistrbution.cDiv <- c(RawDistrbution.cDiv, list(fit_selm_distrubution(c(CadiateNet2$cDiv_D1[idx],
                                                                                       CadiateNet2$cDiv_D2[idx]))))
  }
  names(RawDistrbution.Div) <- unique(CadiateNet2$npoint_ds) %>%
    sort()
  names(RawDistrbution.cDiv) <- unique(CadiateNet2$npoint_ds) %>%
    sort()

  #
  CandidateNet_samll = ad.NetRScore3(netScore = CandidateNet_samll,
                                     RawDistrbution.Div, RawDistrbution.cDiv)
  ad.dropModel = list(AdModelList_1 = RawDistrbution.Div,
                      AdModelList_2 = RawDistrbution.cDiv)
  CandidateNet_samll = CandidateNet_samll[CandidateNet_samll$LR >
                                            0.99, ]
  NetList = tapply(c(1:nrow(CandidateNet_samll), 1:nrow(CandidateNet_samll)),
                   c(CandidateNet_samll$source, CandidateNet_samll$target),
                   list)
  DegreeD = table(unlist(CandidateNet_samll[, 1:2]))
  degreeRange = floor(10^seq(log10(1), log10(max(DegreeD)),
                             0.05)) %>% unique()


  Meancoef.Div <- NULL
  sdcoef.Div <- NULL
  Meandata.Div <- NULL
  Meancoef.cDiv <- NULL
  sdcoef.cDiv <- NULL
  Meandata.cDiv <- NULL
  RandRawScore_Mtime <- NULL
  for(i in 1:repTime){
    message(i)
    RandRawScore = NULL
    for (i in degreeRange) {
      Netid = lapply(NetList, function(x) if (length(x) >=
                                              i) {
        sample(x, i)
      }    else {
        NA
      })
      Netid = Netid[sapply(Netid, length) == i]
      CandidateGenes = names(Netid)
      if (length(CandidateGenes) >= 50) {
        if (i < 20) {
          CandidateGenes = sample(CandidateGenes, 50)
          xtime = 1:50
        }      else {
          CandidateGenes = sample(CandidateGenes, pmin(length(CandidateGenes),
                                                       200))
          xtime = 1:length(CandidateGenes)
        }
        for (j in xtime) {
          TempNet = CandidateNet_samll[Netid[[CandidateGenes[j]]],
          ]
          cDiv.data = TempNet[, c("cDiv_D1.chiSquare",
                                  "cDiv_D2.chiSquare")] %>% as.matrix()
          cDiv.ind = which(TempNet[, 1:2] == CandidateGenes[j],
                           arr.ind = TRUE)
          Raw.Div = sum(TempNet$Div.chiSquare * TempNet$LR)
          Raw.cDiv = sum(cDiv.data[cDiv.ind] * TempNet$LR)
          Raw.Or = sum(pmax(TempNet$Div.chiSquare * TempNet$LR,
                            cDiv.data[cDiv.ind] * TempNet$LR))
          Raw.Plus = sum(TempNet$Div.chiSquare * TempNet$LR +
                           cDiv.data[cDiv.ind] * TempNet$LR)
          Raw.Mpl = sum(TempNet$Div.chiSquare * TempNet$LR *
                          cDiv.data[cDiv.ind] * TempNet$LR)
          Ras_temp = data.frame(gene = CandidateGenes[j], degree = i,
                                Raw.Div = Raw.Div, Raw.cDiv = Raw.cDiv, Raw.Or = Raw.Or,
                                Raw.Plus = Raw.Plus, Raw.Mpl = Raw.Mpl)
          RandRawScore = rbind(RandRawScore, Ras_temp)
        }
      }
    }

    RandRawScore$degree = ceiling(RandRawScore$degree * 0.99)

    Div_Mean_SD_coef_temp <- sea_fit_Mean_sd_coef(Dsize = ceiling(RandRawScore$degree),
                                                  RawScores = RandRawScore$Raw.Div, coarse = F,
                                                  mean_1 = mean(RandRawScore$Raw.Div[RandRawScore$degree ==  1]),
                                                  sd_1 = sd(RandRawScore$Raw.Div[RandRawScore$degree ==  1]) * sdBias)


    cDiv_Mean_SD_coef_temp <- sea_fit_Mean_sd_coef(Dsize = RandRawScore$degree,
                                                   RawScores = RandRawScore$Raw.cDiv,
                                                   coarse = F,
                                                   mean_1 = mean(RandRawScore$Raw.cDiv[RandRawScore$degree == 1]),
                                                   sd_1 = sd(RandRawScore$Raw.cDiv[RandRawScore$degree ==  1]) * sdBias)



    Meancoef.Div <- rbind(Meancoef.Div,Div_Mean_SD_coef_temp$Meancoef)
    sdcoef.Div <- rbind(sdcoef.Div,Div_Mean_SD_coef_temp$sdcoef)
    Meandata.Div <- rbind(Meandata.Div,Div_Mean_SD_coef_temp$Meandata)

    Meancoef.cDiv <- rbind(Meancoef.cDiv,cDiv_Mean_SD_coef_temp$Meancoef)
    sdcoef.cDiv <- rbind(sdcoef.cDiv,cDiv_Mean_SD_coef_temp$sdcoef)
    Meandata.cDiv <- rbind(Meandata.cDiv,cDiv_Mean_SD_coef_temp$Meandata)

    RandRawScore_Mtime <- rbind(RandRawScore_Mtime,RandRawScore)
  }


  model_avg.Div <- make_model(colMeans2(Meancoef.Div)[1],colMeans2(Meancoef.Div)[2],colMeans2(Meancoef.Div)[3])
  model_sd.Div <- make_model(colMeans2(sdcoef.Div)[1],colMeans2(sdcoef.Div)[2],colMeans2(sdcoef.Div)[3])

  RandRawScore_Mtime$Zscore.Div <- (RandRawScore_Mtime$Raw.Div-predict(model_avg.Div,newdata = data.frame(x=RandRawScore_Mtime$degree)))/
    predict(model_sd.Div,newdata = data.frame(x=RandRawScore_Mtime$degree))


  model_avg.cDiv <- make_model(colMeans2(Meancoef.cDiv)[1],colMeans2(Meancoef.cDiv)[2],colMeans2(Meancoef.cDiv)[3])
  model_sd.cDiv <- make_model(colMeans2(sdcoef.cDiv)[1],colMeans2(sdcoef.cDiv)[2],colMeans2(sdcoef.cDiv)[3])

  RandRawScore_Mtime$Zscore.cDiv <- (RandRawScore_Mtime$Raw.cDiv-predict(model_avg.cDiv,newdata = data.frame(x=RandRawScore_Mtime$degree)))/
    predict(model_sd.cDiv,newdata = data.frame(x=RandRawScore_Mtime$degree))
  # hist(RandRawScore_Mtime$Zscore.cDiv)


  # Div


  Meandata.Div$PredictMeanM <- predict(model_avg.Div,newdata = data.frame(x=Meandata.Div$MeaS+1))
  Meandata.Div$PredictStdM <- predict(model_sd.Div,newdata = data.frame(x=Meandata.Div$MeaS+1))
  p2 <- ggplot(Meandata.Div, aes(MeaS + 1, StdM)) + geom_point(size = 5) +
    geom_line(aes(y = PredictStdM), size = 2, color = "red") +
    theme_cowplot_i() + AddBox() + ggtitle("Degree ~ Sd") +
    xlab("Degree") + ylab("Std")
  p1 <- ggplot(Meandata.Div, aes(MeaS + 1, MeaM)) + geom_point(size = 5) +
    geom_line(aes(y = PredictMeanM), size = 2, color = "red") +
    theme_cowplot_i() + AddBox() + ggtitle("Degree ~ Mean") +
    xlab("Degree") + ylab("Mean")

  y = RandRawScore_Mtime$Zscore.Div
  fit_ghyp <- ghyp::fit.ghypuv(y,silent = T)
  pgh    <- ghyp::pghyp(y, object = fit_ghyp, lower.tail = FALSE)
  p3 <- ggplot(data = data.frame(Zscore=y), aes(Zscore)) + # geom_line()
    geom_histogram(mapping = aes(y = after_stat(density)),colour = "white", fill = "black",bins = 50) +
    stat_function(fun = ghyp::dghyp, args = list(object = fit_ghyp),col = "red",size = 2)+
    theme_cowplot_i() +AddBox()+ggtitle('Z scores')
  p_Div <- p1|p2|p3
  ZscoreFit.Div <- list(CombinePs = p_Div,
                        ParmetersInput = list(Meanfit=model_avg.Div,Stdfit=model_sd.Div,ParmetersZsD=fit_ghyp),
                        cal_Zscore_Pvalue = function (RawScore, Dsize, ParmetersInput){
                          Zscores = (RawScore - predict(ParmetersInput$Meanfit, data.frame(x = Dsize)))/predict(ParmetersInput$Stdfit,
                                                                                                                data.frame(x = Dsize))
                          Zscores = as.numeric(Zscores)
                          Pvalues    <- ghyp::pghyp(Zscores, object = ParmetersInput$ParmetersZsD , lower.tail = FALSE, rel.tol=1e-6)
                          Pvalues2    <- ghyp::dghyp(Zscores, object = ParmetersInput$ParmetersZsD )
                          Pvalues[Pvalues == 0] = Pvalues2[Pvalues == 0]
                          data.frame(Zscores = Zscores, Pvalues = Pvalues)}
  )
  # cDiv


  Meandata.cDiv$PredictMeanM <- predict(model_avg.cDiv,newdata = data.frame(x=Meandata.cDiv$MeaS+1))

  Meandata.cDiv$PredictStdM <- predict(model_sd.cDiv,newdata = data.frame(x=Meandata.cDiv$MeaS+1))
  p2 <- ggplot(Meandata.cDiv, aes(MeaS + 1, StdM)) + geom_point(size = 5) +
    geom_line(aes(y = PredictStdM), size = 2, color = "red") +
    theme_cowplot_i() + AddBox() + ggtitle("Degree ~ Sd") +
    xlab("Degree") + ylab("Std")
  p1 <- ggplot(Meandata.cDiv, aes(MeaS + 1, MeaM)) + geom_point(size = 5) +
    geom_line(aes(y = PredictMeanM), size = 2, color = "red") +
    theme_cowplot_i() + AddBox() + ggtitle("Degree ~ Mean") +
    xlab("Degree") + ylab("Mean")

  y = RandRawScore_Mtime$Zscore.cDiv
  fit_ghyp <- ghyp::fit.ghypuv(y,silent = T)
  pgh    <- ghyp::pghyp(y, object = fit_ghyp, lower.tail = FALSE)
  p3 <- ggplot(data = data.frame(Zscore=y), aes(Zscore)) + # geom_line()
    geom_histogram(mapping = aes(y = after_stat(density)),colour = "white", fill = "black",bins = 50) +
    stat_function(fun = ghyp::dghyp, args = list(object = fit_ghyp),col = "red",size = 2)+
    theme_cowplot_i() +AddBox()+ggtitle('Z scores')
  p_cDiv <- p1|p2|p3

  ZscoreFit.cDiv<- list(CombinePs = p_cDiv,
                        ParmetersInput = list(Meanfit=model_avg.cDiv,Stdfit=model_sd.cDiv,ParmetersZsD=fit_ghyp),
                        cal_Zscore_Pvalue = function (RawScore, Dsize, ParmetersInput){
                          Zscores = (RawScore - predict(ParmetersInput$Meanfit, data.frame(x = Dsize)))/predict(ParmetersInput$Stdfit,
                                                                                                                data.frame(x = Dsize))
                          Zscores = as.numeric(Zscores)
                          Pvalues    <- ghyp::pghyp(Zscores, object = ParmetersInput$ParmetersZsD , lower.tail = FALSE, rel.tol=1e-6)
                          Pvalues2    <- ghyp::dghyp(Zscores, object = ParmetersInput$ParmetersZsD )
                          Pvalues[Pvalues == 0] = Pvalues2[Pvalues == 0]
                          data.frame(Zscores = Zscores, Pvalues = Pvalues)}
  )

  # Zscore_plus

  y = RandRawScore_Mtime$Zscore.cDiv+RandRawScore_Mtime$Zscore.Div
  fit_ghyp <- ghyp::fit.ghypuv(y,silent = T)
  pgh    <- ghyp::pghyp(y, object = fit_ghyp, lower.tail = FALSE)
  p_ZsPlus <- ggplot(data = data.frame(Zscore=y), aes(Zscore)) + # geom_line()
    geom_histogram(mapping = aes(y = after_stat(density)),colour = "white", fill = "black",bins = 50) +
    stat_function(fun = ghyp::dghyp, args = list(object = fit_ghyp),col = "red",size = 2)+
    theme_cowplot_i() +AddBox()+ggtitle('Z scores')
  p_ZsPlus
  ZscoreFit.ZsPlus <- list(reslutP=p_ZsPlus,
                           ParmetersInput = fit_ghyp,
                           cal_Pvalue = function (x, ParmetersInput){

                             Pvalues   <- ghyp::pghyp(x, object = ParmetersInput , lower.tail = FALSE, rel.tol=1e-6)
                             Pvalues2  <- ghyp::dghyp(x, object = ParmetersInput )
                             Pvalues[Pvalues == 0] = Pvalues2[Pvalues == 0]
                             data.frame(x = x, Pvalues = Pvalues)
                           }
  )
  #
  NEAModel <- list(DegreeData=Network,
                   ad.dropModel=ad.dropModel,
                   ZscoreFit.Div=ZscoreFit.Div,
                   ZscoreFit.cDiv=ZscoreFit.cDiv,
                   ZscoreFit.ZsPlus=ZscoreFit.ZsPlus)
}



ad.NetRScore3 <- function (netScore, DivFit, cDivFit)
{
  Pchi.Div = toChiSquareX(netScore$Div, netScore$npoint_ds,
                          DivFit)
  Pchi.cDiv_D1 = toChiSquareX(netScore$cDiv_D1, netScore$npoint_ds,
                              cDivFit)
  Pchi.cDiv_D2 = toChiSquareX(netScore$cDiv_D2, netScore$npoint_ds,
                              cDivFit)
  Pchi.Div = addLabel2Colnames(Pchi.Div, label = "Div.", before = TRUE)
  Pchi.cDiv_D1 = addLabel2Colnames(Pchi.cDiv_D1, label = "cDiv_D1.",
                                   before = TRUE)
  Pchi.cDiv_D2 = addLabel2Colnames(Pchi.cDiv_D2, label = "cDiv_D2.",
                                   before = TRUE)
  Pchi = cbind(Pchi.Div, Pchi.cDiv_D1, Pchi.cDiv_D2)
  netScore[, colnames(Pchi)] <- Pchi
  netScore
}


toChiSquareX <- function(rs, bias, DistrbutionList) {
  Pvalues <- rs
  uniBias <- unique(bias)
  Num_Model <- names(DistrbutionList)
  for (i in 1:length(uniBias)) {
    ind <- getCloseseData(data = as.numeric(Num_Model),
                                 uniBias[i], returnIndex = T)
    Pvalues[bias == uniBias[i]]  <-  DistrbutionList[[ind]]$cal_Pvalue(rs[bias ==
                                                                         uniBias[i]], ParmetersInput = DistrbutionList[[ind]]$ParmetersInput)$Pvalues
  }
  Pvalues[Pvalues == 0] <-  9.99999999998465e-313
  chiSquare = sqrt(-log10(Pvalues))
  data.frame(Pvalues = Pvalues, chiSquare = chiSquare)
}

make_model <- function(a,b,c) {
  model <- list(
    a =  a,
    b =  b,
    c =  c,
    fsd_mean = function(x) a * (x^b) + c
  )
  class(model) <- "sd_mean_model"
  model
}

#  predict
predict.sd_mean_model <- function(object, newdata) {
  object$fsd_mean(newdata$x)
}



sea_fit_Mean_sd_coef <- function (Dsize, RawScores, Dstd = NULL, winLength = 10, gap = 0.6,
                                  coarse = FALSE, sd_1 = NULL, mean_1 = NULL)
{
  if (is.null(sd_1)) {
    sd_1 = sd(RawScores/Dsize)
  }
  if (is.null(mean_1)) {
    mean_1 = mean(RawScores/Dsize)
  }
  Dsize = Dsize - 1
  if (is.null(Dstd)) {
    if (min(table(Dsize)) > 20) {
      winLength = min(table(Dsize))
      gap = 1
      coarse = FALSE
    }
    Meandata = MeanSD_By_slidingWindows(ranks = Dsize, scores = RawScores,
                                        winLength = winLength, gap = gap, coarse = coarse)
    colnames(Meandata) = c("MeaS", "MeaM", "StdM")
  }  else {
    Meandata = data.frame(MeaS = Dsize, MeaM = RawScores,
                          StdM = Dstd)
  }

  fmean = function(x, a, b) {
    a * I(x^b) + mean_1
  }
  fsd = function(x, a, b) {
    a * I(x^b) + sd_1
  }
  inital_coef <- lm(Meandata$StdM ~ Meandata$MeaS)
  resultSDPa = minpack.lm::nlsLM(StdM ~ fsd(MeaS, a, b), data = Meandata,
                                 start = list(a = as.numeric(coef(inital_coef)[2]), b = 1))
  for (i in 1:150) {
    resultSDPa = minpack.lm::nlsLM(StdM ~ fsd(MeaS, a, b),
                                   data = Meandata, start = coef(resultSDPa))
  }
  resultSDPa = minpack.lm::nlsLM(StdM ~ fsd(MeaS, a, b), data = Meandata,
                                 start = coef(resultSDPa))

  inital_coef <- lm(Meandata$MeaM ~ Meandata$MeaS)
  resultMePa = minpack.lm::nlsLM(MeaM ~ fmean(MeaS, a, b),
                                 data = Meandata, start = list(a = as.numeric(coef(inital_coef)[2]),
                                                               b = 1))
  for (i in 1:150) {
    resultMePa = minpack.lm::nlsLM(MeaM ~ fmean(MeaS, a,
                                                b), data = Meandata, start = coef(resultMePa))
  }
  resultMePa = minpack.lm::nlsLM(MeaM ~ fmean(MeaS, a, b),
                                 data = Meandata, start = coef(resultMePa))

  sdcoef <- c(coef(resultSDPa),sd_1=sd_1)
  Meancoef <- c(coef(resultMePa),mean_1=mean_1)
  return(list(Meancoef=Meancoef,sdcoef=sdcoef,Meandata=Meandata))


}

getZscore_v2 <- function (EdgeScore, NEAModel, GeneLikelihood,EdgeDataSpecific=NULL)
{
  if ("data.frame" %in% class(EdgeScore)) {
    message("input is a datafrme")
    Network <-  EdgeScore
  }
  else {
    message("input is a scNDS object")
    Network <- EdgeScore@Network
  }
  Network <- ad.NetRScore3(netScore = Network, NEAModel$ad.dropModel$AdModelList_1,
                           NEAModel$ad.dropModel$AdModelList_2)
  LR <- pmin(GeneLikelihood[Network[, 1]], GeneLikelihood[Network[,
                                                                  2]])
  if(is.null(EdgeDataSpecific)){
    EdgeDataSpecific <- 1
  }else{
    EdgeDataSpecific[EdgeDataSpecific>1] <- 1
    EdgeDataSpecific[EdgeDataSpecific<0] <- 0
    LR = pmin(LR, tissueLR)
  }

  chiSquare.LR <- Network[, c("Div.chiSquare", "cDiv_D1.chiSquare",
                              "cDiv_D2.chiSquare")] * as.numeric(LR)
  chiSquare.LR <- addLabel2Colnames(chiSquare.LR, ".LR")
  Network[, colnames(chiSquare.LR)] <- chiSquare.LR
  GeneDegree <- tapply(c(LR, LR), list(c(Network[, 1], Network[,
                                                              2])), function(x) sum(x/max(x)))

  RawS.Div <- getRawScore.Div(net = Network[, 1:3], Div = Network$Div.chiSquare.LR)
  RawS.cDiv  <-  getRawScore.cDiv(net = Network[, 1:3], cDiv_D1 = Network$cDiv_D1.chiSquare.LR,
                               cDiv_D2 = Network$cDiv_D2.chiSquare.LR)

  Zscores <- cbind(RawS.Div[, c(1, 3)],
                  addLabel2Colnames(RawS.Div[, 2, drop = F], label = ".Div"),
                  addLabel2Colnames(RawS.cDiv[, 2, drop = F], label = ".cDiv"))
  Zscores$degree.LR <- GeneDegree[Zscores$Gene]
  ZS.Div <- NEAModel$ZscoreFit.Div$cal_Zscore_Pvalue(RawScore = Zscores$RawScore.Div[1],
                                                    Dsize = Zscores$degree.LR[1], ParmetersInput = NEAModel$ZscoreFit.Div$ParmetersInput)
  ZS.cDiv <- NEAModel$ZscoreFit.cDiv$cal_Zscore_Pvalue(Zscores$RawScore.cDiv,
                                                      Dsize = Zscores$degree.LR, ParmetersInput = NEAModel$ZscoreFit.cDiv$ParmetersInput)
  ZS.ZsPlus <- NEAModel$ZscoreFit.ZsPlus$cal_Pvalue(ZS.Div[,
                                                           1] + ZS.cDiv[, 1], NEAModel$ZscoreFit.ZsPlus$ParmetersInput)
  colnames(ZS.ZsPlus) <- c("Zscores", "Pvalues")
  Zscores  <- cbind(Zscores, addLabel2Colnames(ZS.Div, label = ".Div"),
                  addLabel2Colnames(ZS.cDiv, label = ".cDiv"),
                  addLabel2Colnames(ZS.ZsPlus, label = ".ZsPlus"))
  RankData <- apply(Zscores[, stringr::str_detect(colnames(Zscores),
                                                 "^Zscores.")], 2, function(x) rank(-x, ties.method = "random"))
  colnames(RankData) <- stringr::str_replace(colnames(RankData),
                                            pattern = "^Zscores.", "Rank.")
  Acc = round((1 - (RankData - 1)/nrow(RankData)) * 100, 2)
  colnames(Acc) <- stringr::str_replace(colnames(Acc), pattern = "^Rank.",
                                       "Acc.")
  Zscores <- cbind(Zscores, RankData, Acc)
  if ("data.frame" %in% class(EdgeScore)) {
    list(Zscores=Zscores,Network=Network)
  }  else {
    EdgeScore@Network = Network
    EdgeScore@Zscore = Zscores
    EdgeScore
  }
}

#' Compute gene-level Z-scores v2
#'
#' This function calculates gene-level Z-scores from a scDNS object using
#' three types of null network models: random edge rewiring, label shuffling,
#' and random distribution.
#'
#' @param scDNSobject A scDNS object containing the gene interaction network
#' @param reCreatNEA Whether to re-generate NEA models with stochastic
#' @param PositiveGene Optional character vector specifying a set of positive
#' @param repTime Number of repetitions used when creating each NEA
#' @param testTime Integer. Number of independent stochastic tests performed to
#' identify the most robust seed with the best recovery accuracy of
#' PositiveGene.
#' @param FDRmethods Character. Multiple testing correction method for p-values,
#' such as "BH", "BY", or "bonferroni".
#' @param corr_adjust Whether to apply correlation adjustment when
#' combining multiple Z-scores using the Stouffer method.
#'
#' @return A scDNS object with the following updated slots:
#' - NEAModel$ZscoreList: List of Z-score results for all NEA models
#' - Zscore: Data frame containing combined Z-scores, p-values, adjusted p-values, and accuracy metrics
#' - Network: Updated network information used for computation
#'
#' @export
#'
#' @examples
scDNS_3_GeneZscore_v2 <- function(scDNSobject,reCreatNEA=FALSE,PositiveGene=NULL,repTime=1,testTime=5,FDRmethods='BH',corr_adjust = T)
{
  NEAModel <- scDNSobject@NEAModel
  if(reCreatNEA&!is.null(PositiveGene)){
    Acc_record <- NULL
    for(i in 1:testTime){
      #
      set.seed(123+i)
      Network_rand <- NEAModel$shuffleLabel$DegreeData
      NEAModel_shuffleLabel <- creatNEAModel_robustness(Network=Network_rand,repTime=repTime)

      #
      Network_rand <- NEAModel$RandDistrubution$DegreeData
      NEAModel_RandDistrubution <- creatNEAModel_robustness(Network=Network_rand,repTime=repTime)

      #
      Network_rand <- NEAModel$randNet$DegreeData
      NEAModel_randNet <- creatNEAModel_robustness(Network=Network_rand,repTime=repTime)

      #
      NEAModel_temp <- list(randNet=NEAModel_randNet,shuffleLabel=NEAModel_shuffleLabel,RandDistrubution=NEAModel_RandDistrubution)

      Res_randNet <- getZscore_v2(EdgeScore = scDNSobject@Network, NEAModel = NEAModel_temp$randNet,
                                  GeneLikelihood = scDNSobject@GeneVariability)
      # Res_randNet$Zscores[PositiveGene,]
      Res_shuffleLabel <- getZscore_v2(EdgeScore = scDNSobject@Network, NEAModel = NEAModel_temp$shuffleLabel,
                                       GeneLikelihood = scDNSobject@GeneVariability)
      # Res_shuffleLabel$Zscores[PositiveGene,]
      Res_RandDistrubution <- getZscore_v2(EdgeScore = scDNSobject@Network, NEAModel = NEAModel_temp$RandDistrubution,
                                           GeneLikelihood = scDNSobject@GeneVariability)

      # Res_RandDistrubution$Zscores[PositiveGene,]
      Zs_data <- data.frame(randNet=Res_randNet$Zscores$Zscores.ZsPlus,
                            shuffleLabel=Res_shuffleLabel$Zscores$Zscores.ZsPlus,
                            RandDistrubution=Res_RandDistrubution$Zscores$Zscores.ZsPlus)
      rownames(Zs_data) <- Res_randNet$Zscores$Gene
      p_combine <- combine_stouffer(df = Zs_data[,1:3], robust_center = F, corr_adjust = corr_adjust,FDRmethods = FDRmethods)
      p_combine$Acc <- getAccByScore(p_combine$combined_z)
      print(123+i)
      print(p_combine[rownames(Zs_data)%in%PositiveGene,])
      # sum(p_combine$p_adj<0.1)
      Acc_record <- c(Acc_record,mean(p_combine[rownames(Zs_data)%in%PositiveGene,]$Acc))
    }


    set.seed(123+which.max(Acc_record))
    Network_rand <- NEAModel$shuffleLabel$DegreeData
    NEAModel_shuffleLabel <- creatNEAModel_robustness(Network=Network_rand,repTime=repTime)

    #
    Network_rand <- NEAModel$RandDistrubution$DegreeData
    NEAModel_RandDistrubution <- creatNEAModel_robustness(Network=Network_rand,repTime=repTime)

    #
    Network_rand <- NEAModel$randNet$DegreeData
    NEAModel_randNet <- creatNEAModel_robustness(Network=Network_rand,repTime=repTime)

    NEAModel_temp <- list(randNet=NEAModel_randNet,shuffleLabel=NEAModel_shuffleLabel,RandDistrubution=NEAModel_RandDistrubution)

    Res_randNet <- getZscore_v2(EdgeScore = scDNSobject@Network, NEAModel = NEAModel_temp$randNet,
                                GeneLikelihood = scDNSobject@GeneVariability)
    # Res_randNet$Zscores[c('TP53','TERT','PDGFRA'),]
    Res_shuffleLabel <- getZscore_v2(EdgeScore = scDNSobject@Network, NEAModel = NEAModel_temp$shuffleLabel,
                                     GeneLikelihood = scDNSobject@GeneVariability)
    # Res_shuffleLabel$Zscores[c('TP53','TERT','PDGFRA'),]
    Res_RandDistrubution <- getZscore_v2(EdgeScore = scDNSobject@Network, NEAModel = NEAModel_temp$RandDistrubution,
                                         GeneLikelihood = scDNSobject@GeneVariability)

    # Res_RandDistrubution$Zscores[c('TP53','TERT','PDGFRA'),]
    Zs_data <- data.frame(randNet=Res_randNet$Zscores$Zscores.ZsPlus,
                          shuffleLabel=Res_shuffleLabel$Zscores$Zscores.ZsPlus,
                          RandDistrubution=Res_RandDistrubution$Zscores$Zscores.ZsPlus)
    rownames(Zs_data) <- Res_randNet$Zscores$Gene
    p_combine <- combine_stouffer(df = Zs_data[,1:3], robust_center = F, corr_adjust = F,FDRmethods = FDRmethods)
    p_combine$Acc <- getAccByScore(p_combine$combined_z)
    p_combine$seed <- 123+which.max(Acc_record)
    Zscore <- list(randNet=Res_randNet,
                   shuffleLabel=Res_shuffleLabel,
                   RandDistrubution=Res_RandDistrubution,
                   p_combine=p_combine)
    Zscore
  }else{
    #
    Res_randNet <- getZscore_v2(EdgeScore = scDNSobject@Network, NEAModel = NEAModel$randNet,
                                GeneLikelihood = scDNSobject@GeneVariability)
    #
    Res_shuffleLabel <- getZscore_v2(EdgeScore = scDNSobject@Network, NEAModel = NEAModel$shuffleLabel,
                                     GeneLikelihood = scDNSobject@GeneVariability)
    #
    #
    Res_RandDistrubution <- getZscore_v2(EdgeScore = scDNSobject@Network, NEAModel = NEAModel$RandDistrubution,
                                         GeneLikelihood = scDNSobject@GeneVariability)

    #
    Zs_data <- data.frame(randNet=Res_randNet$Zscores$Zscores.ZsPlus,
                          shuffleLabel=Res_shuffleLabel$Zscores$Zscores.ZsPlus,
                          RandDistrubution=Res_RandDistrubution$Zscores$Zscores.ZsPlus)
    rownames(Zs_data) <- Res_randNet$Zscores$Gene
    p_combine <- combine_stouffer(df = Zs_data[,1:3], robust_center = F, corr_adjust = F,FDRmethods = FDRmethods)
    p_combine$Acc <- getAccByScore(p_combine$combined_z)

    Zscore <- list(randNet=Res_randNet,
                   shuffleLabel=Res_shuffleLabel,
                   RandDistrubution=Res_RandDistrubution,
                   p_combine=p_combine)
    Zscore
  }
  scDNSobject@NEAModel$ZscoreList <- Zscore
  Zscore_res <- Res_RandDistrubution$Zscores
  Zscore_res[,colnames(p_combine)] <- p_combine
  scDNSobject@Zscore <- Zscore_res
  scDNSobject@Network <- Res_RandDistrubution$Network
  scDNSobject
}




combine_stouffer <- function(df,CenterCut=c(0.3,0.7),robust_center = TRUE,corr_adjust=TRUE,FDRmethods = 'BH' ){

  Z <- apply(as.matrix(df), 2, function(x){
    qx <- quantile(x,CenterCut)
    cx <- x[x>qx[1]&x<qx[2]]
    mx <- mean(cx)
    sdx <- sd(cx)
    x <- (x-mx)/sdx
  })

  if (corr_adjust) {
    Sigma <- cov(scale(Z))
  } else {
    Sigma <- diag(ncol(Z))
  }
  w <- rep(1, ncol(Z))
  den <- sqrt( t(w) %*% Sigma %*% w )
  num <- as.numeric(Z %*% w)

  combined_z <- num / as.numeric(den)

  p_comb <- 1 - pnorm(combined_z) # one side
  p_comb[p_comb==0] <- dnorm(combined_z)[p_comb==0]
  print(sum(p.adjust(p_comb,'BH')<0.01))
  print(sum(p.adjust(p_comb,'BH')<0.05))
  print(sum(p.adjust(p_comb,'BH')<0.1))
  data.frame(df, combined_z, p_comb,p_adj=p.adjust(p_comb,FDRmethods),acc=getAccByScore(combined_z)*100)
}

scDNS_cluster_4 <- function(sob,scZscore,DoUmap=TRUE,resolution=c(0.2,0.5,0.8),zsFeatureSel=FALSE){
  sob <- FindVariableFeatures(sob)
  data1 <- sob@assays$RNA@data[sob@assays$RNA@var.features,]
  data2 <- scZscore[,colnames(sob)]
  data1 <- addLabel2Rowname(data1,label = 'raw',sep = '_',before = T)
  if(zsFeatureSel){
    slecteGene <- select_hvg_activity(scZscore,method = 'mean_abs')
    data2 <- addLabel2Rowname(data2[slecteGene$selected,],label = 'scDNS',sep = '_',before = T)
  }else{
    data2 <- addLabel2Rowname(data2,label = 'scDNS',sep = '_',before = T)
  }


  data1 <- rowScale2(data1)
  data2 <- rowScale2(data2)
  data_merge <- rbind(data1,data2)
  sobx <- CreateSeuratObject(counts = data_merge)

  sobx@assays$RNA@data <- data_merge%>%as.matrix()
  sobx@assays$RNA@scale.data <- as.matrix(data_merge)
  sobx@assays$RNA
  sobx <- RunPCA(sobx,features = rownames(sobx))
  if(DoUmap){
    sobx <- RunUMAP(sobx,dims = 1:50)
    sob@reductions[['UMAP_scDNS']] <- sob@reductions[['umap']]
  }
  #
  sobx <- FindNeighbors(sobx)
  for(i in resolution){
    sobx <- FindClusters(sobx,resolution = resolution)
  }
  clusters <- sobx@meta.data[,str_detect(colnames(sobx@meta.data),'RNA_snn_res'),drop=F]
  colnames(clusters) <- str_replace(colnames(clusters),'RNA_snn_res',replacement = 'scDNS_snn')
  # sobx@meta.data[,setdiff(colnames(sob),colnames(sobx))] <- sob@meta.data[,setdiff(colnames(sobx),colnames(sob))]
  sob@meta.data[,colnames(clusters)] <- clusters
  sob@graphs[['scDNS_snn']] <- sobx@graphs$RNA_snn
  sob$scDNS_cluster <- Idents(sobx)
  sob@reductions[['PAC_scDNS']] <- sob@reductions[['pca']]

  sob
}
scDNS_cluster_5 <- function(sob,scZscore,
                            reduction='pca',
                            biasToscDNS=1,
                            red.name='scDNSAll',
                            resolution=c(0.2),
                            doUMAP=TRUE,
                            doTSNE=FALSE,
                            doCluster =TRUE,
                            DoResiduals=FALSE){
  Loading <- rowSums2(scZscore[,colnames(sob)])
  if(DoResiduals){
    message('DoResiduals')
    model <- lm(Loading ~ oldRed)
    Loading <- model$residuals
  }
  Loading <- Loading*(norm(sob@reductions[[reduction]]@cell.embeddings[,1,drop=F],type = "2")/
                        norm(matrix(Loading,ncol = 1),type = "2"))*biasToscDNS
  NewEmbeding <- sob@reductions[[reduction]]@cell.embeddings
  NewEmbeding <- cbind(NewEmbeding,`TempX`=Loading)
  colnames(NewEmbeding) <- paste(red.name,1:ncol(NewEmbeding),sep = '_')
  NewEmbeding <- CreateDimReducObject(embeddings = NewEmbeding)

  sob[[red.name]] <- NewEmbeding

  sob <- ClusterByAnyReduction(sob,
                               RedName = red.name,
                               resolution=resolution,
                               doUMAP=TRUE,
                               doTSNE=TRUE)
  sob

}

scDNS_cluster_6 <- function(sob,
                            scZscore,
                            graph.name=NULL,
                            group.by=NULL,
                            cluster=NULL,
                            DoUmap=TRUE,
                            max_reso=2,
                            min_reso=0.01,
                            target_k=c(2),
                            subcluster.name='scDNS_sub',
                            max.inter=12,
                            minCellinSub=20,
                            zsFeatureSel=FALSE){
  if(!is.null(group.by)){
    Idents(sob) <- sob@meta.data[,group.by]
  }
  if(is.null(cluster)){
    cluster = Idents(sob)%>%unique()%>%as.character()
  }
  if(is.null(graph.name)){
    sob <- FindVariableFeatures(sob)


    sob[['pca_scDNS']] <- RunPCA(rowScale2(scZscore))
    # rownames(sob@reductions$ica_scDNS@cell.embeddings) <- colnames(sob)
    sob <- RunPCA(sob)
    # rownames(sob@reductions$ica@cell.embeddings) <- colnames(sob)
    # combined$RNA_snn_res.0.5
    sob[['pca_scDNS']]@cell.embeddings <- norm(sob[['pca']]@cell.embeddings[,1],type = '2')/norm(sob[['pca_scDNS']]@cell.embeddings[,1],type = '2')*
      sob[['pca_scDNS']]@cell.embeddings


    sob <- FindMultiModalNeighbors(
      object = sob,
      reduction.list = list("pca","pca_scDNS"),
      dims.list = list(1:50, 1:50),l2.norm = F,
      modality.weight.name = "scDNS.weight"
    )

    graph.name <- 'wsnn'
    if(DoUmap){
      sob <- RunUMAP(sob, nn.name = "weighted.nn", reduction.name = "scDNS_umap_c6")
    }

  }
  cellcluster <- Idents(sob)
  cellcluster <- setNames(as.character(cellcluster), names(cellcluster))
  table(cellcluster)
  for(ck in 1:length(cluster)){
    cluster_i <- cluster[ck]

    sub_cells <- WhichCells(object = sob, idents = cluster_i)
    sub.graph <- as.Graph(x = sob[[graph.name]][sub_cells,
                                                sub_cells])
    max_reso_i  <- max_reso
    min_reso_i <- min_reso

    for(i in 1:max.inter){
      sub.clusters <- FindClusters(object = sub.graph,
                                   resolution = mean(c(max_reso_i,min_reso_i)),
                                   algorithm = 1)
      if(length(levels(sub.clusters[[1]]))==target_k){
        break()
      }
      if(length(levels(sub.clusters[[1]]))>target_k){
        max_reso_i <-  mean(c(max_reso_i,min_reso_i))
      }
      if(length(levels(sub.clusters[[1]]))<target_k){
        min_reso_i <-  mean(c(max_reso_i,min_reso_i))
      }

    }
    message(cluster[ck])
    if(min(table( sub.clusters[, 1]))<minCellinSub){
      message(c('The minimum number of cells in cluster of ',cluster_i, ' is less than ',minCellinSub))
    }else{
      sub.clusters[, 1] <- paste(cluster_i, sub.clusters[, 1], sep = "_")
      cellcluster[rownames(sub.clusters)] <- sub.clusters[, 1]
    }
  }
  table(cellcluster)
  sob[[subcluster.name]] <- cellcluster
  sob
}
