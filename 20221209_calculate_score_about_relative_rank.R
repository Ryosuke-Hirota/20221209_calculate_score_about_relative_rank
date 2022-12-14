# This script is to calculate relative rank, p value and score based on rank of correlation coefficient.
# made 2022/12/09

# make new directory
setwd("C:/Rdata")
dir.create("20221209_calculate_score_about_relative_rank")

# import resutl
# this result is located at ""
setwd("C:/Rdata/20221108_correlaiton_between_RBP_level_and_miRNA_level_with_excluding_outlier")
cor.result <-read.table("summary_of_correlaiton_between_RBP_level_and_miRNA_level_with_excluding_outlier_and_considering_no_expression.txt",sep="\t",header = T,stringsAsFactors = F)
cor.result[,1] <-paste0(cor.result[,1],"_vs_",cor.result[,2])
cor.result <-cor.result[,c(1,5)]
colnames(cor.result)[1] <-"combination"

# import table about correlation coefficient
# this table is located at ""
setwd("C:/Rdata/20221206_correlaiton_analysis_between_RBP_and_miRNA_in_each_site")
cor.table <-read.table("table_of_correlation_coefficient_by_each_site.txt",sep="\t",header = T,stringsAsFactors = F)

# merge two tables
cor.result.table <-merge(cor.result,cor.table,by="combination")

# get name of site primary
site <-colnames(cor.table)[-1]

# set cutoff
cutoff <-seq(0,900,50)

setwd("C:/Rdata/20221209_calculate_score_about_relative_rank")

#calculate relative ranking of positive or negative
for (i in 13:length(cutoff)) {
  for (k in 3:27) {
      
  # set cutoff
  data <-cor.result.table[cor.result.table[,2]>=cutoff[i],c(1,k)]
  data[,2] <-ifelse(is.na(data[,2]),0,data[,2])
  
  # make table
  if(k==3){  
    p.rank.df <-as.data.frame(matrix(nrow = nrow(data),ncol = 26))
    p.rank.df[,1] <-data[,1]
    colnames(p.rank.df) <-colnames(cor.table)
    
    n.rank.df <-as.data.frame(matrix(nrow = nrow(data),ncol = 26))
    n.rank.df[,1] <-data[,1]
    colnames(n.rank.df) <-colnames(cor.table)
  }  
  
  # order by rank about positive correlation and calculate relative rank
  p.rank <-data[order(data[,2],decreasing = T),]
  p.rank[,2] <-1:nrow(p.rank)
  p.rank[,2] <-p.rank[,2]/nrow(p.rank)-1/2/nrow(p.rank)
  p.rank <-p.rank[order(p.rank[,1]),]
  p.rank.df[,k-1] <-p.rank[,2]
  
  # order by rank about negative correlation and calculate relative rank
  n.rank <-data[order(data[,2]),]
  n.rank[,2] <-1:nrow(n.rank)
  n.rank[,2] <-n.rank[,2]/nrow(n.rank)-1/2/nrow(n.rank)
  n.rank <-n.rank[order(n.rank[,1]),]
  n.rank.df[,k-1] <-n.rank[,2]
  
  # write table about rank
  if(k==27){
  write.table(p.rank.df,paste0("relative_rank_of_positive_correlation_cutoff_",cutoff[i],".txt"),sep="\t",row.names = F,quote = F)
  write.table(n.rank.df,paste0("relative_rank_of_negative_correlation_cutoff_",cutoff[i],".txt"),sep="\t",row.names = F,quote = F)
  }
  }}

# make list of relative rank
rank <-list.files(path="C:/Rdata/20221209_calculate_score_about_relative_rank",pattern = "relative_rank")
rank.list <-as.list(NULL)
for (i in 1:length(cutoff)) {
  l <-grep(paste0("_",cutoff[i],".txt"),rank)
  tables <-rank[l]
  rank.list[[i]] <-tables
  }


#cutoff sample > 25
#setwd("C:/Rdata")
#cell.df <-read.table("CCLE_cell_lines_whose_have_over_25_samples.txt",header = T,sep="\t",stringsAsFactors = F)
#cell <-cell.df[,1]

# cutoff sample > 30
setwd("C:/Rdata/rho_for_RBP_vs_miRNA_by_site")
site.cutoff <-read.table("CCLE_cell_lines_whose_have_over_30_samples.txt",sep="\t",header = T,stringsAsFactors = F)
cell <-site.cutoff[,1]

# set function to calculate pvalue of chi square test
cal.pvalue <-function(x,y){
  r <-nrow(x)
  for (i in 1:r) {
  value <-sum(log(x[i,]))*-2
  p <-pchisq(value,2*y,lower.tail = F)
  if(i==1){
    p.values <-p
  }else{
  p.values <-append(p.values,p)
  }}
  return(p.values)
  }

setwd("C:/Rdata/20221209_calculate_score_about_relative_rank")
for (i in 1:length(cutoff)) {
  # read each list of relative rank
  rank1 <-rank.list[[i]][[2]]
  rank2 <-rank.list[[i]][[1]]
  
  p.rank.table <-read.table(rank1,sep="\t",header = T,stringsAsFactors = F)
  n.rank.table <-read.table(rank2,sep="\t",header = T,stringsAsFactors = F)
  
  #keep tissue that satisfy 25 samples over
  c <-match(cell,colnames(p.rank.table))
  p.rank.table <-p.rank.table[,c(1,c)]
  n.rank.table <-n.rank.table[,c(1,c)]
  
  # make each table to calculate p.value (positive and negative)
  pp.df <-p.rank.table
  rownames(pp.df) <-p.rank.table[,1]
  pp.df <-pp.df[,-1]

  np.df <-n.rank.table
  rownames(np.df) <-n.rank.table[,1]
  np.df <-np.df[,-1]

  #calculate p.value based on chi-squared distribution with 2*10 degress of freedom
  pp <-cal.pvalue(pp.df,length(cell))
  np <-cal.pvalue(np.df,length(cell))

  # write each table about p.value
  pp.data <-p.rank.table[,c(1,2)]
  pp.data[,2] <-pp
  colnames(pp.data)[2] <-"positive_p.value"
  
  np.data <-n.rank.table[,c(1,2)]
  np.data[,2] <-np
  colnames(np.data)[2] <-"negative_p.value"
  
  # write each table about p.value
  write.table(pp.data,paste0("pvalue_of_positive_rank_",cutoff[i],".txt"),sep="\t",quote = F,row.names = F)
  write.table(np.data,paste0("pvalue_of_negative_rank_",cutoff[i],".txt"),sep="\t",quote = F,row.names = F)
  }

# make each list about p.value
pp.lists <-list.files(path="C:/Rdata/20221209_calculate_score_about_relative_rank",pattern = "pvalue_of_positive")
np.lists <-list.files(path = "C:/Rdata/20221209_calculate_score_about_relative_rank",pattern = "pvalue_of_negative")

# calculate score
for (i in 1:length(cutoff)) {
  # import lists of pvalue
  p.list <-read.table(pp.lists[i],sep="\t",header = T,stringsAsFactors = F)
  np.list <-read.table(np.lists[i],sep="\t",header = T,stringsAsFactors = F)
  
  # merge 
  score.df <-merge(p.list,np.list,by="combination")
  score.df[,4] <-NA
  colnames(score.df)[4] <-"score"
  
  #calculate score by size of each p.value
  ps <-score.df[,2]<score.df[,3]
  ns <-score.df[,2]>score.df[,3]
  score.df[ps,4] <-log10(2*score.df[ps,2])*-1
  score.df[ns,4] <-log10(2*score.df[ns,3])
  
  # write table of score
  write.table(score.df,paste0("table_of_score_",cutoff[i],".txt"),sep="\t",quote = F,row.names = F)
}