#' Get experimantal condition information
#'
#' @return Experimental condition table
expCondition = function()
{
  f = system.file('extdata','Experimental_conditions.txt', package = 'bimir')
  data = read.delim(f, sep="\t")
  return(data)
}

#' Get matrix size
#'
#' @param target a matrix
#' @return The size of input matrix
size = function(target)
{
  if (is.null(dim(target))) { return(0) }
  return(dim(target)[1] * dim(target)[2])
}

#' Get zero rate
#'
#' @param target a matrix
#' @return The ratio of zero in a binarized matrix
#' @export
get_zero_rate=function(target) { return(length(target[target == 0])/length(target)) }

#' Bicluster extension process
#'
#' @param MIR_profile A binarized microRNA profile
#' @param seed_row List of a seed row index
#' @param seed_col List of a seed column index
#' @param zero_ratio Zero ratio to be allowed.
#' @return extended matrix
#' @export
bicluster_extension_process = function( MIR_profile, seed_row, seed_col, zero_ratio )
{
  SEED = MIR_profile[seed_row,seed_col]
  biclust_row = rownames(SEED)
  biclust_col = colnames(SEED)

  while(get_zero_rate(SEED) < zero_ratio )
  {
    biclust_row = which(rownames(MIR_profile)%in%rownames(SEED))
    biclust_col = which(colnames(MIR_profile)%in%colnames(SEED))
    n1 = length(biclust_row)/length(biclust_col)
    n2 = 1/n1
    if(n1 >= 2){
      n=n1
      modified_zero_rate = function(target){return((length(target[target == 0])+ceiling(n))/length(target))}
      row_zero = apply(X=MIR_profile[-biclust_row, biclust_col], 1, modified_zero_rate)
      col_zero = apply(X=MIR_profile[biclust_row, -biclust_col], 2, get_zero_rate)
    }
    if(n2 >= 2){
      n=n2
      modified_zero_rate = function(target){return((length(target[target == 0])+ceiling(n))/length(target))}
      row_zero = apply(X = MIR_profile[-biclust_row, biclust_col], 1, get_zero_rate)
      col_zero = apply(X = MIR_profile[biclust_row, -biclust_col], 2, modified_zero_rate)
    }
    if(n1 < 2 & n2 < 2){
      row_zero = apply(X = MIR_profile[-biclust_row, biclust_col], 1, get_zero_rate)
      col_zero = apply(X = MIR_profile[biclust_row, -biclust_col], 2, get_zero_rate)
    }

    min_row = min(row_zero)
    min_col = min(col_zero)

    l1 = length(which(row_zero == min_row))
    l2 = length(which(col_zero == min_col))

    if(min_row < min_col | (min_row == min_col & l1 >= l2))
    {
      change = which(row_zero == min_row)
      tt = MIR_profile[-biclust_row,biclust_col][change,]
      if(is.null(dim(tt)))
      {
        biclust_row = c(biclust_row,as.integer(change))
        UPDATED = rbind(SEED,tt)
        rownames(UPDATED)[dim(UPDATED)[1]] = names(change)
      }
      else
      {
        UPDATED=rbind(SEED,tt)
        biclust_row = sort(c(biclust_row,as.integer(change)))
      }

      if(get_zero_rate(UPDATED)<zero_ratio) {
        SEED = UPDATED
        rowind = which(rownames(MIR_profile)%in%rownames(SEED))
        colind = which(colnames(MIR_profile)%in%colnames(SEED))
        SEED=MIR_profile[sort(rowind),sort(colind)]}
      else { return(SEED) }
    }
    if(min_row > min_col | (min_row == min_col & l1 < l2))
    {
      change = which(col_zero == min_col)
      tt = MIR_profile[biclust_row,-biclust_col][,change]
      if(is.null(dim(tt)))
      {
        biclust_col = c(biclust_col,as.integer(change))
        UPDATED=cbind(SEED,tt)
        colnames(UPDATED)[dim(UPDATED)[2]] = names(change)
      }
      else
      {
        UPDATED = cbind(SEED,tt)
        biclust_col = sort(c(biclust_col,as.integer(change)))
      }

      if(get_zero_rate(UPDATED)<zero_ratio) {
        SEED = UPDATED
        rowind = which(rownames(MIR_profile)%in%rownames(SEED))
        colind = which(colnames(MIR_profile)%in%colnames(SEED))
        SEED=MIR_profile[sort(rowind),sort(colind)]
      }
      else { return(SEED) }
    }

    biclust_row=which(rownames(MIR_profile)%in%rownames(SEED))
    biclust_col=which(colnames(MIR_profile)%in%colnames(SEED))
  }
}

#' Bicluster trimming process
#'
#' @param table A matrix that will be reduced
#' @param zero_ratio Zero ratio allowed.
#' @return Reduced matrix
#' @export
bicluster_trimming_process=function( table, zero_ratio )
{
  zero_rate_rows = apply(table,1,get_zero_rate); max_row_zero_rate = max(zero_rate_rows)
  zero_rate_cols = apply(table,2,get_zero_rate); max_col_zero_rate = max(zero_rate_cols)

  while( max_row_zero_rate > zero_ratio || max_col_zero_rate > zero_ratio )
  {
    if( max_row_zero_rate >= max_col_zero_rate )
    {
      delete_target = which(zero_rate_rows == max_row_zero_rate , arr.ind = T)
      table=table[-delete_target,]
    }

    else
    {
      delete_target = which(zero_rate_cols == max_col_zero_rate , arr.ind = T)
      table=table[,-delete_target]
    }
    zero_rate_rows = apply(table,1,get_zero_rate); max_row_zero_rate = max(zero_rate_rows)
    zero_rate_cols = apply(table,2,get_zero_rate); max_col_zero_rate = max(zero_rate_cols)
  }
  return(table)
}

#'Progressive bicluster extension
#'
#' @param MIR_profile Binarized microRNA profile
#' @param biclust_row List of index of row in seed biclust
#' @param biclust_col List of index of columns in seed biclust
#' @param step_number The number of extension process
#' @param finalZR final zero rate allowed.
#' @return list of row and column symbols of extended bicluster
#' @export
PBE = function( MIR_profile, biclust_row, biclust_col, step_number, finalZR)
{
  nth_step = 1:step_number
  intZRs = finalZR * nth_step / step_number
  size_back = -1
  for ( intZR in intZRs )
  {
    extended_seed = bicluster_extension_process(MIR_profile, biclust_row, biclust_col, intZR)
    reducted_seed = bicluster_trimming_process(extended_seed, intZR)

    if(size_back > size(reducted_seed)){break}
    else{size_back = size(reducted_seed)}
    biclust_row = match( rownames( reducted_seed ), rownames( MIR_profile ) )
    biclust_col = match( colnames( reducted_seed ), colnames( MIR_profile ) )
  }
  return( list(rows = rownames( reducted_seed ),cols = colnames( reducted_seed )))
}

#' Make all biclusters
#'
#' @param MIR_PROFILE A matrix of MIR profile.
#' @param FCcutoff Binarization fold change cutoff in log2-scale for MIR_PROFILE. Target up-regulation biclusters will be created for positive FC cutoff, and down-regulation biclusters will be generated with negative FC cutoff.
#' @param REPETITION The number of running ensemble bicluster function.
#' @param STEP_NUMBER The number of extension process
#' @param ZERORATE Final zero rate allowed.
#' @import biclust
#' @return list of seed and biclusters
#' @export
make_biclust = function(MIR_PROFILE, FCcutoff, REPETITION, STEP_NUMBER, ZERORATE)
{
  if(FCcutoff==0){stop('FCcutoff must be larger or lesser than zero.')}
  if(FCcutoff>0){PROFILE_POS = TRUE; PROFILE_NEG=FALSE}
  if(FCcutoff<0){PROFILE_NEG = TRUE; PROFILE_POS=FALSE}
  ### set tab to make biclust
  biclusts = list()
  seedlist = list()
  ##### DO ##### option set negative / positive
  for(i in 1:REPETITION)
  {
    print(proc.time())
    if(PROFILE_NEG)
    {
      t=MIR_PROFILE
      t=data.matrix(t)
      t[is.na(t)] = 0

      #Binarize the MIR profile
      profile = t
      index0 = which(t > FCcutoff)
      index1 = which(t <= FCcutoff)
      profile[index0] = 0
      profile[index1] = 1
      profile.seed = profile

      # Generation of seed biclusters using ensemble BIMAX algorithm
      rn = runif(n = 1, min = 1, max = 1E8)
      set.seed(rn)
      ensembled = try(ensemble(profile.seed, bimax.grid(), bootstrap = FALSE),silent = T)
      if(class(ensembled) == "try-error") {print("sorry. A seed cannot be found from the profile")
      }else
      {
        for(iterator in 1:ensembled@Number)
        {
          seed = biclusternumber(ensembled, iterator)
          biclust_row = seed[[1]]$Rows
          biclust_col = seed[[1]]$Cols

          try({biclust_step = PBE(profile, biclust_row, biclust_col,STEP_NUMBER,ZERORATE);
          seedlist[[length(seedlist)+1]] = list(rows = biclust_row, cols = biclust_col)
          biclusts = append(biclusts,list(biclust_step))
          } )
        }
      }
    }
    if(PROFILE_POS)
    {
      t=MIR_PROFILE
      t=data.matrix(t)
      t[is.na(t)] = 0

      # Binarize the MIR profile
      profile = t
      index0 = which(t < FCcutoff)
      index1 = which(t >= FCcutoff)
      profile[index0] = 0
      profile[index1] = 1
      profile.seed = profile

      # Generation of the seed biclusters
      rn = runif(n = 1, min = 1, max = 1E8)
      set.seed(rn)
      ensembled = try(ensemble(profile.seed,bimax.grid(), bootstrap = FALSE),silent = T)
      if(class(ensembled) == "try-error") {print("sorry. A seed cannot be found from the profile")
      }else
      {
        for(iterator in 1:ensembled@Number)
        {
          seed = biclusternumber(ensembled, iterator)
          biclust_row = seed[[1]]$Rows
          biclust_col = seed[[1]]$Cols
          seed = profile[biclust_row,biclust_col]

          try({biclust_step = PBE(profile, biclust_row, biclust_col,STEP_NUMBER,ZERORATE);
          seedlist[[length(seedlist)+1]] = list(rows=biclust_row, cols = biclust_col)
          biclusts = append(biclusts,list(biclust_step))
          })
        }
      }
    }
  }
  return(list(biclusts = biclusts,seedlist = seedlist))
}

#' Calculate similarity of two biclusters
#'
#' @param biclust1 First bicluster matrix
#' @param biclust2 Second bicluster matrix
#' @return Similarity of two biclusters
similarity = function(biclust1,biclust2) # Pairwise distance of two biclusters
{
  introw = intersect(rownames(biclust1),rownames(biclust2))
  intcol = intersect(colnames(biclust1),colnames(biclust2))
  nr = length(introw)
  nc = length(intcol)
  minsize = min(size(biclust1), size(biclust2))
  if(!(length(introw)*length(intcol))) {return(1)}
  return(1- (nr*nc / minsize ))
}

#' Pariwise distance matrix for the entire list of biclusters
#'
#' @param biclusts Lists of biclusters.
#' @param MIR_profile binarized microRNA profile
#' @return pariwise distance matrix
comatrix = function(biclusts, MIR_profile) # Pairwise distance matrix for the entire list of biclusters
{
  profile = MIR_profile
  cotab = matrix(0,length(biclusts),length(biclusts))

  for(i in 1:length(biclusts))
  {
    bicluster_1 = profile[which(rownames(profile)%in%biclusts[[i]]$rows),which(colnames(profile)%in%biclusts[[i]]$cols)]
    for(j in i:length(biclusts)) { cotab[i,j] = similarity(bicluster_1, profile[ which(rownames(profile)%in%biclusts[[j]]$rows) , which(colnames(profile)%in%biclusts[[j]]$cols)]) }
  }
  ind=lower.tri(cotab)
  cotab[ind]= t(cotab)[ind]
  return(cotab)
}

#' Merging similar biclusters
#'
#' @param cut A vector of bicluster classificatino
#' @param biclusts List of Biclusters to be merged
#' @param i iteration
#' @param MIR_profile Binarized microRNA profile
#' @param ZERORATE Final zero ratio allowed.
#' @return merged biclusters
merge = function(cut,biclusts,i, MIR_profile, ZERORATE)
{
  group_index = which(cut == i,arr.ind = T)
  rowin = c()
  colin = c()
  for(i in group_index)
  {
    rowin = append(rowin,biclusts[[i]]$rows)
    colin = append(colin,biclusts[[i]]$cols)
  }
  rowin = unique(rowin)
  colin = unique(colin)

  rough = MIR_profile[which(rownames(MIR_profile)%in%rowin),which(colnames(MIR_profile)%in%colin)]
  clear = bicluster_trimming_process(rough,ZERORATE)
  return(list(rows = rownames(clear),cols = colnames(clear)))
}

#' Merging similar biclusters
#'
#' @param MIR_profile binarized microRNA profile
#' @param tree_cutoff Similarity cutoff
#' @param biclusts list of extended biclusters
#' @param seedlist list of seeds
#' @param FCcutoff Fold change cutoff (in log2 sclale) for binarization of MIR profile
#' @param mir Name of miRNA
#' @param ZERORATE Final zero ratio allowed
#' @return Extended bicluster list
#' @export
merge_bicluster=function(MIR_profile, tree_cutoff = 0.5,biclusts, seedlist, FCcutoff=log2(1.3), mir='testmiRNA', ZERORATE)
{
  if(FCcutoff>0){PROFILE_POS = TRUE; PROFILE_NEG=FALSE}
  if(FCcutoff<0){PROFILE_NEG = TRUE; PROFILE_POS=FALSE}

  biclusts_o = biclusts
  if(length(biclusts)==0){stop('no biclusters to merge')}
  dis = as.dist(comatrix(biclusts, MIR_profile))

  if(max(dis)==0){
    biclusts = list()
    biclusts[[1]] = list(rows = biclusts_o[[1]]$rows, cols = biclusts_o[[1]]$cols)
    temp = rep(1, length(biclusts_o))
  }else if(min(dis) > tree_cutoff){
    biclusts = biclusts_o
    temp = 1:length(biclusts)
  }else{
    temp = list()
    while(min(dis) < tree_cutoff)
    {
      cut = cutree(hclust(dis,method = "average"),h=tree_cutoff)
      temp = append(temp, list(cut))

      biclusts_new = list()
      for(i in 1:length(unique(cut)))
      {
        biobject = merge(cut,biclusts,i, MIR_profile, ZERORATE)
        biclusts_new = append(biclusts_new,list(biobject))
      }
      biclusts = biclusts_new; rm(biclusts_new)
      dis = as.dist(comatrix(biclusts,MIR_profile))
    }
  }

  print(length(biclusts))

  level = length(temp)-1

  if(level>0)
  {
    for(i in level:1)
    {
      assign = temp[[(i+1)]][temp[[i]]]
    }
  }else if(level == 0){assign = unlist(temp)}

  return(biclusts)
  # for(i in 1:length(biclusts))
  # {
  #   merged_row = biclusts[[i]]$rows
  #   merged_col = biclusts[[i]]$cols
  #
  #   rowind = which(rownames(MIR_profile)%in%merged_row)
  #   colind = which(colnames(MIR_profile)%in%merged_col)
  #
  #   bicluster_name = if(PROFILE_NEG){paste(mir, "_biclust_down_",i,'.txt',sep="")
  #   }else if(PROFILE_POS){paste(mir, "_biclust_up_",i,'.txt',sep="")}
  #
  #   write.table(MIR_profile[rowind,colind], file = bicluster_name,quote = F)
  # }

}

#' Progressive extension and merging of biclusters
#'
#' @param MIR_profile A matrix of microRNA profile where rows are cell conditions and columns are sequence-based miRNA targets.
#' @param biclust.path Directory where biclusts will be saved
#' @param mir Name of microRNA
#' @param FCcutoff Log2-fold change binarization cutoff for seed. Default = log2(1.3)
#' @param REPETITION The number of repetition that the ensemble function runs to extract seed biclusters. Default = 10
#' @param STEP_NUMBER The number of extension process. Default = 10
#' @param ZERORATE Final zero rate allowed. Default = 0.1
#' @param tree_cutoff Similarity cutoff. Default = 0.5
#' @return Extended and merged biclusters in the directory assigned to ‘biclust.path’.
#' @import biclust
#' @export
PBE_MERGE = function(MIR_profile, mir, biclust.path='./', FCcutoff=log2(1.3), REPETITION=10, STEP_NUMBER=10, ZERORATE=0.1, tree_cutoff=0.5)
{
  Expcond = expCondition()
  if(FCcutoff>0){PROFILE_POS=TRUE; PROFILE_NEG=FALSE}
  if(FCcutoff<0){PROFILE_POS=FALSE; PROFILE_NEG=TRUE}
  MIR_profile = data.matrix(MIR_profile)
  if(nrow(MIR_profile)==0){stop('Empty profile.')}

  biclusts_and_seeds = make_biclust(MIR_profile, FCcutoff, REPETITION, STEP_NUMBER, ZERORATE)

  if(length(biclusts_and_seeds$biclusts)==0){stop("No biclusters exnteded")
  }else{
    # Merge the extended biclusters
    setwd(biclust.path)
    dir.create(mir, showWarnings = F)
    MIR_profile_b = MIR_profile
    if(PROFILE_POS){
      index0 = which(MIR_profile < FCcutoff)
      index1 = which(MIR_profile >= FCcutoff)
      MIR_profile_b[index0] = 0
      MIR_profile_b[index1] = 1
    }
    if(PROFILE_NEG){
      index0 = which(MIR_profile > FCcutoff)
      index1 = which(MIR_profile <= FCcutoff)
      MIR_profile_b[index0] = 0
      MIR_profile_b[index1] = 1
    }
    setwd(mir)
    bclist = merge_bicluster(MIR_profile_b, tree_cutoff, biclusts=biclusts_and_seeds$biclusts, seedlist=biclusts_and_seeds$seedlist, FCcutoff, mir,ZERORATE)


    for(i in 1:length(bclist))
    {
      merged_row = bclist[[i]]$rows
      merged_col = bclist[[i]]$cols

      rowind = which(rownames(MIR_profile)%in%merged_row)
      colind = which(colnames(MIR_profile)%in%merged_col)

      bicluster_file = if(PROFILE_NEG){paste("biclust_down_",i,'.txt',sep="")
      }else if(PROFILE_POS){paste("biclust_up_",i,'.txt',sep="")}

      condition_file = if(PROFILE_NEG){paste("Experimental_condition_down_",i,'.txt',sep="")
      }else if(PROFILE_POS){paste("Experimental_condition_up_",i,'.txt',sep="")}

      genelist_file = if(PROFILE_NEG){paste("Targetlist_down_",i,'.txt',sep="")
      }else if(PROFILE_POS){paste("Targetlist_up_",i,'.txt',sep="")}

      result_biclust = MIR_profile[rowind,colind]
      biclust_exp_cond = Expcond[which(Expcond$Symbol%in%rownames(result_biclust)),]
      biclust_genelist = colnames(result_biclust)

      write.table(result_biclust, file =bicluster_file,quote = F)
      write.table(biclust_exp_cond, file = condition_file, quote = F, sep = "\t", row.names = F)
      write(biclust_genelist, file = genelist_file, ncolumns = 1)
    }
    cat(paste('Finished. Check your folder:', biclust.path))
  }
}
