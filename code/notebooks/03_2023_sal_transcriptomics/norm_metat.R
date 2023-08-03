##test.main
#example script to show how to format the input count data for taxon-specific scaling 

# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

#BiocManager::install("compcodeR")
require(DESeq2)
require(tidyverse)
generate.org.mat <- function(nSamples = 6, nOrg = 5, nVars = 1000, nDiffFunc = 100, CountVec = c(1e7,5e6,1e6,5e5,1e5), Sampling.Rate.Mat)
{
  ##use to generate simulated data
  #nSamples = number of samples per condition
  #nOrg = number of organisms to simulate
  #nVars = number of features to generate
  #nDiffFunc = number of differentially expressed functions
  #CountVec = base library size for each organism
  #Sampling.Rate.Mat: library size parameters; rows = samples (i.e. 12 for 6 samples per condition), columns = organisms 

  require(compcodeR)
  
  #error checking    
  if (ncol(Sampling.Rate.Mat) != nOrg | nrow(Sampling.Rate.Mat) != nSamples*2 | length(CountVec) != nOrg)
  {
    cat('dimension mismatch for samples or organisms, please check.','\n')
    return(-1)
  }
  
  B.gen <- list()
  for (k in 1:nOrg)
  {
    B.gen[[k]] <- generateSyntheticData(dataset = "B", n.vars = nVars,
                                        samples.per.cond = nSamples, n.diffexp = nDiffFunc,
                                        repl.id = 1, seqdepth = CountVec[k],
                                        fraction.upregulated = 0.5,
                                        filter.threshold.total = -1,
                                        filter.threshold.mediancpm = -1,
                                        minfact = Sampling.Rate.Mat[,k],
                                        maxfact = Sampling.Rate.Mat[,k]
    )
  }
  return(B.gen)
}


DESeq2.norm.mat <- function(Xmat,cond,type)
{
  
  #Xmat = raw count data for one species
  Xmat.col = ncol(Xmat)
  Xmat.row = nrow(Xmat)
  colData <- data.frame(condition = cond, type = type)
  Xmat = round(Xmat)
  storage.mode(Xmat) <- 'integer'
  dds <- DESeqDataSetFromMatrix(countData = Xmat, colData = colData, design = ~condition)
  colData(dds)$condition = factor(colData(dds)$condition, levels=unique(cond))
  dds <- DESeq(dds, quiet = TRUE)
  #normalize the data
  YMat <- Xmat/rep(dds@colData@listData$sizeFactor, each = (Xmat.row))
  return(YMat)
}

DESeq2.result <- function(Xmat,cond,type)
{
  Xmat.col = ncol(Xmat)
  Xmat.row = nrow(Xmat)
  Xmat = round(Xmat)
  storage.mode(Xmat) <- 'integer'
  colData <- data.frame(condition = cond, type = type)
  dds <- DESeqDataSetFromMatrix(countData = Xmat, colData = colData, design = ~condition)
  
  colData(dds)$condition = factor(colData(dds)$condition,
                                  levels=unique(cond))
  
  #stop DESeq2 from performing additional normalization
  normFactors <- matrix(1,ncol = Xmat.col, nrow = Xmat.row)
  normalizationFactors(dds) <- normFactors
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds)
  return(res)
}

abind.matrix <- function(Xarray=NULL,XMat)
{
  require(abind)
  #Xmat = raw count data for one species
  #Xmat rows = features, coloums = samples
  #Xarray = array containing data for already added species
  Xarray <- abind(Xarray, XMat, along = 3)
  return(Xarray)
}

DESeq2.tax.specific <- function(Xarray,cond.vec,type.vec)
{
  require(DESeq2)
  #Xarray array[feature,sample,organism]
  #cond.vec = the different condition types of the sample
  #type.vec = the type of sequence, for example 'single-read'
  
  nOrg = dim(Xarray)[3]
  nSamples = dim(Xarray)[2]
  #check input 
  if (length(cond.vec) != nSamples | length(type.vec) != nSamples)
  {
    return(-1)
  }
  
  Scaled.Mat <- matrix(0,ncol = ncol(Xarray), nrow = nrow(Xarray))
  cat("scaling ",nOrg," different matrices.\n")
  for (i in 1:nOrg)
  {
    Scaled.Mat <- Scaled.Mat + 
      DESeq2.norm.mat(Xarray[,,i], cond.vec,type.vec)
  }
  return(DESeq2.result(Scaled.Mat,cond.vec,type.vec))
}


test.main <- function()
{
  #returns the taxon-specific scaling normalized DESeq2 results of a metatranscriptome
  nOrg = 5                                                      #Simulation: number of organism to simulate
  nSamples = 6                                                  #Simulation: number of samples per condition
  CountVec = c(1e7,5e6,1e6,5e5,1e5)                             #Simulation: base library size (LS)
  Sampling.Rate.Mat = matrix(0,ncol = nOrg, nrow = nSamples*2)  #Simulation: LS variation for each organism and sample
  cond.vec <- c(rep('A',nSamples),rep('B',nSamples))            #conditions for each sample
  type.vec <- rep('1',2*nSamples)                               #type of read
  
  for (k in 1:nOrg)
  {
    Sampling.Rate.Mat[,k] = sample(seq(0.5,2,0.1), size = nSamples*2, replace = TRUE)
  }
  #replace this part with your count data##
  res = generate.org.mat(nSamples = nSamples, nOrg = nOrg, CountVec = CountVec, Sampling.Rate.Mat = Sampling.Rate.Mat)
  
  #create array for each organism count data matrix
  array = NULL
    for (k in 1:nOrg)
    {
      array = abind.matrix(array,res[[k]]@count.matrix)
    }  
  
  
  #########################################

  return(DESeq2.tax.specific(Xarray = array, cond.vec = cond.vec, type.vec = type.vec))
}

