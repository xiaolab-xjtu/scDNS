#' filterNetowrkFromGAT
#'
#' @param scDNSob scDNSob
#' @param GAT_outNet
#'
#' @return scDNSob
#' @export
#'
#' @examples
filterNetowrkFromGAT <- function(scDNSob,GAT_outNet='./data/attention_layer1_epoch1000_wide_confidence.csv'){
  message('raw network will be sorted in scDNSob@Other$Network')
  scDNSob@Other$Network <- scDNSob@Network
  GATNet <- read.csv(GAT_outNet)
  NetworkLabel <- scDNSob@Network

  Hcon_index <- paste(NetworkLabel[,1],NetworkLabel[,2],sep = '_')%in%GATNet$edge[GATNet$confidence =='High_confidence']
  scDNSob@Network <- scDNSob@Network[Hcon_index,]
  scDNSob
}

#' filterNetowrkFromDREMI
#'
#' @param scDNSob scDNSob
#' @param rmP the threshold of significance(def:0.01)
#'
#' @return
#' @export
#'
#' @examples
filterNetowrkFromDREMI <- function(scDNSob,rmP=0.01){
  message('raw network will be sorted in scDNSob@Other$Network')
  scDNSob@Other$Network <- scDNSob@Network
  DRMEI_th <- quantile(scDNSob@NEAModel$RandDistrubution$DegreeData[,
                                                                    stringr::str_detect(colnames(scDNSob@NEAModel$RandDistrubution$DegreeData),
                                                                               pattern = 'DREMI')]%>%as.matrix(),1-rmP)
  maxDREMI <- MatrixGenerics::rowMaxs(scDNSob@Network[,stringr::str_detect(colnames(scDNSob@Network),pattern = 'DREMI')]%>%as.matrix())
  message('Non-significant interactions (based on DREMI) are removed from the network.')
  message(sum(maxDREMI<DRMEI_th))
  scDNSob@Network <- scDNSob@Network[maxDREMI>=DRMEI_th,]
  scDNSob
}
#' getZscore
#'
#' @param EdgeScore a list or dataframe from getKLD_cKLDnetwork
#' @param NEAModel NEA Mdoel from creatNEAModel
#' @param Likelihood a vector from CalGeneInfAmount
#'
#' @return
#' @keywords internal
#'
#' @examples
getZscore<-function(EdgeScore,NEAModel,Likelihood){


  if('data.frame'%in%class(EdgeScore)){
    message('input is a datafrme')
    Network = EdgeScore
  }else{
    message('input is a scNDS object')
    Network = EdgeScore@Network
  }

  Network = ad.NetRScore3(netScore = Network,NEAModel$ad.dropModel$AdModelList_1,NEAModel$ad.dropModel$AdModelList_2)

  #
  LR = pmin(Likelihood[Network[,1]],Likelihood[Network[,2]])
  # chiSquare.LR  <- Network[,stringr::str_detect(colnames(Network),'chiSquare')]*as.numeric(LR)
  chiSquare.LR  <- Network[,c('Div.chiSquare','cDiv_D1.chiSquare','cDiv_D2.chiSquare')]*as.numeric(LR)
  chiSquare.LR <- addLabel2Colnames(chiSquare.LR,'.LR')
  Network[,colnames(chiSquare.LR)] <- chiSquare.LR
  # Network$Div.chiSquare.LR = Network$Div.chiSquare *LR
  GeneDegree = tapply(c(LR,LR), list(c(Network[,1],Network[,2])), function(x)sum(x/max(x)))




  RawS.Div = getRawScore.Div(net = Network[,1:3],
                             Div = Network$Div.chiSquare.LR)
  RawS.cDiv = getRawScore.cDiv(net = Network[,1:3],
                               cDiv_D1 = Network$cDiv_D1.chiSquare.LR,
                               cDiv_D2 = Network$cDiv_D2.chiSquare.LR)

  RawS.OR = getRawScore.Or(net = Network[,1:3],
                           Div = Network$Div.chiSquare,
                           cDiv_D1 = Network$cDiv_D1.chiSquare.LR,
                           cDiv_D2 = Network$cDiv_D2.chiSquare.LR)

  RawS.Plus = getRawScore.Plus(net = Network[,1:3],
                               Div = Network$Div.chiSquare.LR,
                               cDiv_D1 = Network$cDiv_D1.chiSquare.LR,
                               cDiv_D2 = Network$cDiv_D2.chiSquare.LR)
  RawS.Mpl = getRawScore.Mpl(net = Network[,1:3],
                             Div = Network$Div.chiSquare.LR,
                             cDiv_D1 = Network$cDiv_D1.chiSquare.LR,
                             cDiv_D2 = Network$cDiv_D2.chiSquare.LR)

  Zscores = cbind(RawS.Div[,c(1,3)],
                  addLabel2Colnames(RawS.Div[,2,drop=F],label = '.Div'),
                  addLabel2Colnames(RawS.cDiv[,2,drop=F],label = '.cDiv'),
                  addLabel2Colnames(RawS.OR[,2,drop=F],label = '.OR'),
                  addLabel2Colnames(RawS.Plus[,2,drop=F],label = '.Plus'),
                  addLabel2Colnames(RawS.Mpl[,2,drop=F],label = '.Mpl'))
  Zscores$degree.LR = GeneDegree[Zscores$Gene]



  ZS.Div = NEAModel$ZscoreFit.Div$cal_Zscore_Pvalue(Zscores$RawScore.Div,
                                                    Dsize = Zscores$degree.LR,
                                                    ParmetersInput = NEAModel$ZscoreFit.Div$ParmetersInput)
  # ZS.Div <- cbind(RawS.Div,ZS.Div)
  # ZS.Div = ZscoreFit.Div$cal_Zscore_Pvalue(RawS.Div$RawScore,
  #                                                     Dsize = RawS.Div$degree,
  #                                                     ParmetersInput = ZscoreFit.Div$ParmetersInput)

  ZS.cDiv = NEAModel$ZscoreFit.cDiv$cal_Zscore_Pvalue(Zscores$RawScore.cDiv,
                                                      Dsize = Zscores$degree.LR,
                                                      ParmetersInput = NEAModel$ZscoreFit.cDiv$ParmetersInput)
  ZS.OR = NEAModel$ZscoreFit.OR$cal_Zscore_Pvalue(Zscores$RawScore.OR,
                                                  Dsize = Zscores$degree.LR,
                                                  ParmetersInput = NEAModel$ZscoreFit.OR$ParmetersInput)
  ZS.Plus = NEAModel$ZscoreFit.Plus$cal_Zscore_Pvalue(Zscores$RawScore.Plus,
                                                      Dsize = Zscores$degree.LR,
                                                      ParmetersInput = NEAModel$ZscoreFit.Plus$ParmetersInput)

  # ZS.Plus = ZscoreFit.Plus$cal_Zscore_Pvalue(RawS.Plus$RawScore,
  #                                                     Dsize = RawS.Plus$degree,
  #                                                     ParmetersInput = NEAModel$ZscoreFit.Plus$ParmetersInput)

  ZS.Mpl = NEAModel$ZscoreFit.Mpl$cal_Zscore_Pvalue(Zscores$RawScore.Mpl,
                                                    Dsize = Zscores$degree.LR,
                                                    ParmetersInput = NEAModel$ZscoreFit.Mpl$ParmetersInput)
  ZS.ZsPlus <- NEAModel$ZscoreFit.ZsPlus$cal_Pvalue(ZS.Div[,1]+ZS.cDiv[,1],
                                                         NEAModel$ZscoreFit.ZsPlus$ParmetersInput)
  colnames(ZS.ZsPlus) <- c('Zscores','Pvalues')


  Zscores = cbind(Zscores,
                  addLabel2Colnames(ZS.Div,label = '.Div'),
                  addLabel2Colnames(ZS.cDiv,label = '.cDiv'),
                  addLabel2Colnames(ZS.OR,label = '.OR'),
                  addLabel2Colnames(ZS.Plus,label = '.Plus'),
                  addLabel2Colnames(ZS.Mpl,label = '.Mpl'),
                  addLabel2Colnames(ZS.ZsPlus,label = '.ZsPlus'))
  RankData = apply(Zscores[,stringr::str_detect(colnames(Zscores),'^Zscores.')], 2, function(x)rank(-x,ties.method = 'random'))
  colnames(RankData) = stringr::str_replace(colnames(RankData),pattern = '^Zscores.','Rank.')
  Acc = round((1-(RankData-1)/nrow(RankData))*100,2)
  colnames(Acc) = stringr::str_replace(colnames(Acc),pattern = '^Rank.','Acc.')
  Zscores = cbind(Zscores,RankData,Acc)
  if('data.frame'%in%class(EdgeScore)){
    # Zscores[Zscores$Gene%in%'DSP',]
    # fit_poly(log10(Zscores$degree.Div),Zscores$Zscores.Plus)
    Zscores
  }else{
    EdgeScore@Network = Network
    EdgeScore@Zscore = Zscores
    # Zscores[Zscores$Gene%in%'DSP',]
    # fit_poly(log10(Zscores$degree.Div),Zscores$Zscores.Plus)
    EdgeScore
  }

}


addLabel2Colnames <-function(x,label,sep = '',before=FALSE){
  if(before){
    colnames(x) = paste(label,colnames(x),sep = sep)
  }else{
    colnames(x) = paste(colnames(x),label,sep = sep)
  }

  x
}
getRawScore.Div<- function(net,Div){
  degree = table(c(net[,1],net[,2]))
  Rawdata=data.frame(Gene=names(degree),RawScore=0,degree=as.numeric(degree))
  rownames(Rawdata)=Rawdata$Gene
  Div.r  <- tapply(c(Div,Div), c(net[,1],net[,2]), sum)

  Rawdata[names(Div.r),'RawScore']=Div.r
  Rawdata
}


getRawScore.cDiv<-function(net,cDiv_D1,cDiv_D2){
  degree = table(c(net[,1],net[,2]))
  Rawdata=data.frame(Gene=names(degree),RawScore=0,degree=as.numeric(degree))
  rownames(Rawdata)=Rawdata$Gene

  raw.1=tapply(cDiv_D1, net[,1], sum)
  raw.2=tapply(cDiv_D2, net[,2], sum)
  Rawdata[names(raw.1),'RawScore']=Rawdata[names(raw.1),'RawScore']+raw.1
  Rawdata[names(raw.2),'RawScore']=Rawdata[names(raw.2),'RawScore']+raw.2
  Rawdata
}
getRawScore.Or<-function(net,Div,cDiv_D1,cDiv_D2){
  degree = table(c(net[,1],net[,2]))
  Rawdata=data.frame(Gene=names(degree),RawScore=0,degree=as.numeric(degree))
  rownames(Rawdata)=Rawdata$Gene

  raw.1=tapply(pmax(cDiv_D1,Div), net[,1], sum)
  raw.2=tapply(pmax(cDiv_D2,Div), net[,2], sum)
  Rawdata[names(raw.1),'RawScore']=Rawdata[names(raw.1),'RawScore']+raw.1
  Rawdata[names(raw.2),'RawScore']=Rawdata[names(raw.2),'RawScore']+raw.2
  Rawdata
}
getRawScore.Plus<-function(net,Div,cDiv_D1,cDiv_D2){
  degree = table(c(net[,1],net[,2]))
  Rawdata=data.frame(Gene=names(degree),RawScore=0,degree=as.numeric(degree))
  rownames(Rawdata)=Rawdata$Gene

  raw.1=tapply(cDiv_D1+Div, net[,1], sum)
  raw.2=tapply(cDiv_D2+Div, net[,2], sum)
  Rawdata[names(raw.1),'RawScore']=Rawdata[names(raw.1),'RawScore']+raw.1
  Rawdata[names(raw.2),'RawScore']=Rawdata[names(raw.2),'RawScore']+raw.2
  Rawdata
}
getRawScore.Mpl<-function(net,Div,cDiv_D1,cDiv_D2){
  degree = table(c(net[,1],net[,2]))
  Rawdata=data.frame(Gene=names(degree),RawScore=0,degree=as.numeric(degree))
  rownames(Rawdata)=Rawdata$Gene

  raw.1=tapply(cDiv_D1*Div, net[,1], sum)
  raw.2=tapply(cDiv_D2*Div, net[,2], sum)
  Rawdata[names(raw.1),'RawScore']=Rawdata[names(raw.1),'RawScore']+raw.1
  Rawdata[names(raw.2),'RawScore']=Rawdata[names(raw.2),'RawScore']+raw.2
  Rawdata
}


# toChiSquareX <- function(rs,bias,DistrbutionList){
#   Pvalues <-  rs
#   uniBias <- unique(bias)
#   Num_Model <- names(DistrbutionList)
#   for(i in 1:length(uniBias)){
#     ind=getCloseseData(data = as.numeric(Num_Model),uniBias[i],returnIndex = T)
#     Pvalues[bias==uniBias[i]] = DistrbutionList[[ind]]$cal_Pvalue(rs[bias==uniBias[i]],
#                                                                   ParmetersInput = DistrbutionList[[ind]]$ParmetersInput)$Pvalues
#   }
#
#   Pvalues[Pvalues==0]=1e-312
#   chiSquare = sqrt(-log10(Pvalues))
#   data.frame(Pvalues=Pvalues,chiSquare=chiSquare)
# }

