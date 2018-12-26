#'@title load_FCtable
#'
#'@description Loading a large table of fold change values for 20,639 human genes in 5,158 cell conditions.
#'@return A fold change matrix of 5,158 rows and 20,639 columns.
#'@export
load_FCtable = function()
{
  cat('FC matrix is now loading ...\n')
  f = system.file('extdata', 'FC1.rds', package = 'bimir')
  fc = readRDS(f)
  cat("1/10\n")
  for(i in 2:10)
  {
  	f = system.file('extdata', paste('FC',i,'.rds',sep=""), package = 'bimir')
  	fc_temp = readRDS(f)
  	fc = cbind(fc,fc_temp)
 	cat(i,"/10\n",sep="")
  }
  cat('Loaded!\n')
  return(fc)
}

#'@title Load human gene list
#'@description 20,639 human genes
#'@return A vector of 20,639 human genes
#'@export
#'@import utils
getGenelist = function()
{
  f = system.file('extdata','genelist.rds', package='bimir')
  genelist = readRDS(f)
  genelist = gsub('-','.',genelist,fixed = T)
  return(genelist)
}

#'@title Load human miRNA list
#'@description human mature miRNA list
#'@return A vector of 2,632 human mature miRNAs
#'@export
#'@import utils
getmiRNAlist = function()
{
  f = system.file('extdata', 'miRNAlist.rds', package='bimir')
  miRNAlist = readRDS(f)
  return(miRNAlist)
}

#'@title Load human miRNA target index
#'@description human miRNA target index
#'@return A vector of human miRNA target index
#'@import utils
getAllmiRNATargetIndex = function()
{
  f = system.file('extdata', 'miRNATargets.rds', package='bimir')
  index = readRDS(f)
  return(index)
}

#'Get sequence-based miRNA targets
#'
#'@description Returns miRNA targets predicted from three or more algorithms.
#'
#'@param miRNA miRNA name. One of miRNAs listed in "miRNAlist" data.
#'@return A vector of sequence-based miRNA targets
#'@examples get_miR_target('hsa-miR-1-3p')
#'@export
get_miR_target = function(miRNA)
{
  genelist = getGenelist()
  miRNAlist = getmiRNAlist()
  miRNATargets = getAllmiRNATargetIndex()
  mir_index = which(miRNAlist%in%miRNA)
  if(length(mir_index)==0){stop(paste(miRNA,'is not in the miRNAlist. Type getmiRNAlist() to refer to the available miRNAs'))}
  index_target = miRNATargets[mir_index]
  index_target = unlist(strsplit(index_target,split="\t"))
  index_target = as.numeric(index_target)
  index_target = unique(index_target)
  targets = genelist[index_target]
  return(targets)
}

#'Get MIRprofile
#'
#'@description Returns MIR profile for input miRNA.
#'
#'@param miRNA miRNA name. One of miRNAs listed in "miRNAlist" data.
#'@param FCtable FC matrix of 20,639 genes for 5,158 conditions. If NULL, it will be loaded from the local.
#'@param FCcutoff Fold change cutoff in log2 scale. Default = log2(1.3). To extract MIR profile for target down-regulation in test condition, take negative value.
#'@return A matrix of fold change value of sequence-based miRNA targets under selected cell conditions.
#'@examples fc=load_FCtable(); MIRprofile = getMIRprofile('hsa-miR-1-3p',fc,log2(1.3))
#'@import stats utils
#'@export
getMIRprofile = function(miRNA, FCtable=NULL, FCcutoff=log2(1.3))
{
  if(is.null(miRNA)){stop('human miRNA name must be given. Refer to the list of miRNAs from "miRNAlist" data in the package.')}
  if(is.null(FCtable)){ FCtable = load_FCtable() }
  if(FCcutoff==0){stop('FCcutoff must be larger or lesser than 0.')}
  miRNAlist = getmiRNAlist()
  genelist = getGenelist()
  mir_index = which(miRNAlist%in%miRNA)
  if(length(mir_index)==0){stop('There is no such microRNA in the DB. Refer to the list of miRNAs by typing getmiRNAlist().')}
  ConditionOverlapP = function(X, seqTar, FCcut)
  {
    seqTar = intersect(genelist, seqTar)
    geneidx = if(FCcutoff>0){which(X>FCcut)}else{which(X<FCcut)}
    genes = genelist[geneidx]
    L1 = length(intersect(genes, seqTar)) - 1
    L2 = length(genes)
    L3 = 20639 - L2
    L4 = length(seqTar)
    p = phyper(L1, L2, L3, L4, lower.tail=F)
    p = signif(p, 4)
    return(p)
  }
  MIRtargets = get_miR_target(miRNA)
  overlapP = apply(FCtable, 1, ConditionOverlapP, seqTar = MIRtargets, FCcut = FCcutoff)
  overlapFDR = p.adjust(overlapP, method = 'fdr')
  signifCondition = which(overlapFDR<0.05)

  result_profile = FCtable[signifCondition, MIRtargets]
  return(result_profile)
}

