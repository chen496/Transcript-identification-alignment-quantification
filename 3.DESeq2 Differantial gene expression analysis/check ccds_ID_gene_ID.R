

setwd("/media/root/10AF08F010AF08F03/root/WeiLi_cDNA_analysis")

###install Biostrings to read fasta file
## BiocManager::install("Biostrings")
library("Biostrings")
fastFile<-readDNAStringSet("CCDS_nucleotide.current.fna")
seq_name<-names(fastFile)
length(seq_name)
length(unique((seq_name)))


sequence<-paste(fastFile)
df<-data.frame(seq_name,seq_name)

library(stringr)
# temp<-matrix(unlist(str_split(seq_name,"\\|")),ncol=3,byrow = TRUE)
# ccds_ID<-temp[,1]
# ccds_ID<-unique(ccds_ID)
# length(ccds_ID)


CCDS.current<-read.delim(file="CCDS.current.txt")
dim(CCDS.current)
which(table(CCDS.current$gene_id)==5)
CCDS.current[CCDS.current$gene_id==CCDS.current$gene_id[11787],]




CCDS.current$gene_id[11787]

head(CCDS.current)
CCDS.current$ccds_id
ccds_ID<-paste0(CCDS.current$ccds_id,"|Mm108|chr",CCDS.current$X.chromosome)
gene_id<-matrix(CCDS.current$gene_id,ncol=1)
rownames(gene_id)<-ccds_ID
setdiff(seq_name,intersect(seq_name,ccds_ID))
gene_id[seq_name,]





####check whether all ccds_IDs have been mapped to gene_IDs
length(intersect(ccds_ID,CCDS.current$ccds_id))
length(setdiff(intersect(ccds_ID,CCDS.current$ccds_id),ccds_ID))


## mouse_ref_mapping.txt is a mapping from ccds_ID to genes. As shown in the snippet below, 
## each line of the mapping file contains a gene identifier and a ccds_ID.
CCDS.current_ref_mapping<-cbind.data.frame(gene_id[seq_name,],seq_name)
head(CCDS.current_ref_mapping)
write.table(CCDS.current_ref_mapping,file="CCDS_nucleotide.current_ref_mapping.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)





### install RSEM
# tar -xzf RSEM-1.2.25.tar.gz
# cd RSEM-1.2.25
# make -j 8
# make ebseq
# cd ..
# cd ..




##### check the alignment results
##### the reanscript_IDs actually are ccds_IDs
NGS128_genes_results<-read.table("RSEM_result/NGS128_S21.genes.results",header = TRUE)
dim(NGS128_genes_results)
head(NGS128_genes_results)

NGS128_genes_results[1:200,]

length(which(NGS128_genes_results$TPM!=0))

NGS128_isoforms_results<-read.table("RSEM_result/NGS128_S21.isoforms.results",header = TRUE,stringsAsFactors=F)
dim(NGS128_isoforms_results)
head(NGS128_isoforms_results)
NGS128_isoforms_results[1:200,]
#### the nunber of isoforms with nonzero TPM
length(which(NGS128_isoforms_results$TPM!=0))

index<-order(NGS128_isoforms_results[,"TPM"],decreasing = TRUE)
NGS128_isoforms_results[index,]



NGS130_genes_results<-read.table("RSEM_result/NGS130_S23.genes.results",header = TRUE)
dim(NGS130_genes_results)
head(NGS130_genes_results)

NGS130_genes_results[1:200,]

length(which(NGS130_genes_results$TPM!=0))



#### differential analysis
# Each line describes a gene and contains 7 fields: 
# the gene name, posterior probability of being equally expressed (PPEE), 
# posterior probability of being differentially expressed (PPDE), 
# posterior fold change of condition 1 over condition 2 (PostFC), 
# real fold change of condition 1 over condition 2 (RealFC), 
# mean count of condition 1 (C1Mean) and mean count of condition 2 (C2Mean). 
# For fold changes, PostFC is recommended over the RealFC and you can find the definition of these two fold changes 
# in the description of PostFC function of EBSeq vignette. 
diff_res<-read.table("RSEM_result/differential_analysis/GeneMat.de.txt",header = TRUE,row.names = 1)
dim(diff_res)
head(diff_res)
diff_gene<-rownames(diff_res)
# diff_res[order(diff_res[,"PostFC"],decreasing = TRUE),]
range(diff_res$PostFC)
diff_index<-vector()
for(i in 1:length(diff_gene)){
  tmp<-which(CCDS.current$gene_id==diff_gene[i])
  if(length(tmp)>1){
    s<-CCDS.current$gene[tmp]
    if(length(unique(s))>1){
      print(i)
    }
    diff_index[i]<-tmp[1]
  }else{
    diff_index[i]<-tmp
  }
    
  
}
diff_gene[6]

write.csv(cbind.data.frame(CCDS.current[diff_index,],diff_res),file="differential_analysis_res.csv",row.names = FALSE)


####
temp<-cbind.data.frame(CCDS.current[diff_index,],diff_res)
head(temp)
#### 
seq2current<-read.delim("CCDS2Sequence.current.txt")
seq2current<-seq2current[seq2current$source=="NCBI",]
head(seq2current)
### find the protein_id according to the ccds_id
protein_id<-list()
for(i in 1:length(temp$ccds_id)){

  tmp<-which(seq2current$X.ccds==as.character(temp$ccds_id[i]))
  s<-as.character(seq2current$protein_ID[tmp])
  protein_id[[i]]<-s

    if(length(unique(s))>1){
      print(i)
  }
}

library(stringr)
diff_protein<-vector()
for(i in 1:length(protein_id)){
  diff_protein[i]<-str_c(protein_id[[i]],collapse = "\\")
}
write.csv(cbind.data.frame(CCDS.current[diff_index,],diff_res,diff_protein),file="differential_analysis_include_proteinName_res.csv",row.names = FALSE)





############## differential analysis results of isoforms
CCDS.current<-read.delim(file="CCDS.current.txt")
dim(CCDS.current)
head(CCDS.current)
CCDS.current$ccds_id
ccds_ID<-paste0(CCDS.current$ccds_id,"|Mm108|chr",CCDS.current$X.chromosome)

diff_res<-read.table("RSEM_result/differential_analysis/GeneMat.isoforms.de_fdr0.05.txt",header = TRUE,row.names = 1)
dim(diff_res)
head(diff_res)
diff_gene<-rownames(diff_res)
# diff_res[order(diff_res[,"PostFC"],decreasing = TRUE),]
range(diff_res$PostFC)
diff_index<-vector()
for(i in 1:length(diff_gene)){
  tmp<-which(ccds_ID==diff_gene[i])
  if(length(tmp)>1){
    s<-CCDS.current$gene[tmp]
    if(length(unique(s))>1){
      print(i)
    }
    diff_index[i]<-tmp[1]
  }else{
    diff_index[i]<-tmp
  }
  
  
}
diff_gene[6]

write.csv(cbind.data.frame(CCDS.current[diff_index,],diff_res),file="differential_analysis_isoforms_fdr0.05_res.csv",row.names = FALSE)



which(CCDS.current$gene_id)

####
temp<-cbind.data.frame(CCDS.current[diff_index,],diff_res)
head(temp)
#### 
seq2current<-read.delim("CCDS2Sequence.current.txt")
seq2current<-seq2current[seq2current$source=="NCBI",]
head(seq2current)
### find the protein_id according to the ccds_id
protein_id<-list()
for(i in 1:length(temp$ccds_id)){
  
  tmp<-which(seq2current$X.ccds==as.character(temp$ccds_id[i]))
  s<-as.character(seq2current$protein_ID[tmp])
  protein_id[[i]]<-s
  
  if(length(unique(s))>1){
    print(i)
  }
}

library(stringr)
diff_protein<-vector()
for(i in 1:length(protein_id)){
  diff_protein[i]<-str_c(protein_id[[i]],collapse = "\\")
}
write.csv(cbind.data.frame(CCDS.current[diff_index,],diff_res,diff_protein),file="differential_analysis_isoforms_include_proteinName_res.csv",row.names = FALSE)







NGS130_genes_results<-read.table("RSEM_result/NGS130_S23.genes.results",header = TRUE)
colum_names<-colnames(NGS130_genes_results)
#### replace the transcript_id with ccds_id
colum_names[2]<-"ccds.id"
colum_names
NGS130_isoforms_results<-read.table("RSEM_result/NGS130_S23.isoforms.results",header = TRUE)

NGS128_genes_results<-read.table("RSEM_result/NGS128_S21.genes.results",header = TRUE)
NGS128_isoforms_results<-read.table("RSEM_result/NGS128_S21.isoforms.results",header = TRUE)

colnames(NGS130_genes_results)<-colum_names
colnames(NGS130_isoforms_results)<-colum_names
colnames(NGS128_genes_results)<-colum_names
colnames(NGS128_isoforms_results)<-colum_names


length(which(NGS130_genes_results$TPM!=0))
length(which(NGS128_genes_results$TPM!=0))

length(which(NGS130_isoforms_results$TPM!=0))
length(which(NGS128_isoforms_results$TPM!=0))

write.csv(NGS130_genes_results[order(NGS130_genes_results[,"TPM"],decreasing = TRUE),],file="NGS130_genes.csv",row.names = FALSE)
write.csv(NGS130_isoforms_results[order(NGS130_isoforms_results[,"TPM"],decreasing = TRUE),],file="NGS130_isoforms.csv",row.names = FALSE)
write.csv(NGS128_genes_results[order(NGS128_genes_results[,"TPM"],decreasing = TRUE),],file="NGS128_genes.csv",row.names = FALSE)
write.csv(NGS128_isoforms_results[order(NGS128_isoforms_results[,"TPM"],decreasing = TRUE),],file="NGS128_isoforms.csv",row.names = FALSE)









setwd("/media/root/10AF08F010AF08F03/root/WeiLi_cDNA_analysis")
NGS128<-read.table("RSEM_result/NGS128_S21.genes.results",header = TRUE)

NGS128[order(NGS128[,"TPM"],decreasing = TRUE),]
NGS128[order(NGS128[,"TPM"],decreasing = TRUE),][1:3,]

NGS130<-read.table("RSEM_result/NGS130_S23.genes.results",header = TRUE)
NGS130[order(NGS130[,"TPM"],decreasing = TRUE),]
NGS130[order(NGS130[,"TPM"],decreasing = TRUE),][1:3,]

for(i in 1:nrow(NGS128)){
  if(NGS128$gene_id[i]!=NGS130$gene_id[i])
    print(i)
}

library(DESeq2)
all((NGS128$gene_id) == (NGS130$gene_id))
cts1<-as.matrix(cbind(as.integer(NGS128$expected_count),
                      as.integer(NGS130$expected_count)))
rownames(cts1)<-NGS128$gene_id
colnames(cts1)<-c("NGS128","NGS130")
head(cts1)
coldata1<-cbind.data.frame(condition=c("NGS128","NGS130"),type=c("paired-end",
                                                                 "paired-end"))

coldata1$condition <- factor(coldata1$condition)
coldata1$type <- factor(coldata1$type)
coldata1
rownames(coldata1)<-colnames(cts1)


all(rownames(coldata1) == colnames(cts1))


library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts1,
                              colData = coldata1,
                              design = ~ 1)
dds
dim(dds)
keep <- rowSums(counts(dds)) > 10
dds <- dds[keep,]
dds
dds <- DESeq(dds)
res <- results(dds)
res

index<-vector()
t<-rownames(res)
for(i in 1:nrow(res)){
  I<-str_split(t[i],"\\.")[[1]][1]
  index[i]<-which(CCDS.current$gene_id==I)[1]
  
}
head(CCDS.current[index,])
temp<-CCDS.current[index,]



### find the protein_id according to the ccds_id
protein_id<-list()
for(i in 1:length(temp$ccds_id)){
  
  tmp<-which(seq2current$X.ccds==as.character(temp$ccds_id[i]))
  s<-as.character(seq2current$protein_ID[tmp])
  protein_id[[i]]<-s
  
  if(length(unique(s))>1){
    print(i)
  }
}
library(stringr)
diff_protein<-vector()
for(i in 1:length(protein_id)){
  diff_protein[i]<-str_c(protein_id[[i]],collapse = "\\")
}
res_temp<-cbind.data.frame(CCDS.current[index,c("gene", "gene_id","ccds_id")],res$log2FoldChange,diff_protein)
res_temp<-res_temp[order(res_temp[,"res$log2FoldChange"],decreasing = TRUE),]
colnames(res_temp)<-c("gene","gene_id","ccds_id","log2FoldChange","protein_NCBI")
dim(res_temp)

NGS128[keep,]

NGS128_RPK_count<-NGS128[keep,]
NGS130_RPK_count<-NGS130[keep,]
all(rownames(NGS128_RPK_count) == rownames(NGS130_RPK_count))
A<-cbind(NGS128_RPK_count$expected_count/NGS128_RPK_count$effective_length*1000,
         NGS128_RPK_count$effective_length,
         NGS128_RPK_count$expected_count,
         NGS130_RPK_count$expected_count/NGS130_RPK_count$effective_length*1000,
         NGS130_RPK_count$effective_length,
         NGS130_RPK_count$expected_count)
rownames(A)<-NGS128_RPK_count$gene_id
### RPK: divide the read counts by the length of each gene in kilobases.
colnames(A)<-c("NGS128_RPK","NGS128_effective_length","NGS128_expected_count",
               "NGS130_RPK","NGS130_effective_length","NGS130_expected_count")
all(res_temp$gene_id %in% rownames(A))

A[as.character(res_temp$gene_id),]
res_temp<-cbind(res_temp,A[as.character(res_temp$gene_id),])
colnames(res_temp)<-c("gene","gene_id","ccds_id","log2FoldChange","protein_NCBI",
                      "NGS128_RPK","NGS128_effective_length","NGS128_expected_count",
                      "NGS130_RPK","NGS130_effective_length","NGS130_expected_count")
write.csv(res_temp,file="Deseq2_differential_analysis_gene_include_proteinName_res.csv",row.names = FALSE)

res_temp1<-res_temp




NGS128<-read.table("RSEM_result/NGS128_S21.isoforms.results",header = TRUE)
which(table(as.character(NGS128$gene_id))==3)[1]

NGS128[which(NGS128$gene_id==100434),]


NGS128[order(NGS128[,"TPM"],decreasing = TRUE),]
NGS128[order(NGS128[,"TPM"],decreasing = TRUE),][1:3,]

NGS130<-read.table("RSEM_result/NGS130_S23.isoforms.results",header = TRUE)
NGS130[order(NGS130[,"TPM"],decreasing = TRUE),]
NGS130[order(NGS130[,"TPM"],decreasing = TRUE),][1:3,]

for(i in 1:nrow(NGS128)){
  if(NGS128$gene_id[i]!=NGS130$gene_id[i])
    print(i)
}

library(DESeq2)
all((NGS128$gene_id) == (NGS130$gene_id))
cts1<-as.matrix(cbind(as.integer(NGS128$expected_count),
                      as.integer(NGS130$expected_count)))
rownames(cts1)<-NGS128$gene_id
colnames(cts1)<-c("NGS128","NGS130")
head(cts1)
coldata1<-cbind.data.frame(condition=c("NGS128","NGS130"),type=c("paired-end",
                                                                 "paired-end"))

coldata1$condition <- factor(coldata1$condition)
coldata1$type <- factor(coldata1$type)
coldata1
rownames(coldata1)<-colnames(cts1)


all(rownames(coldata1) == colnames(cts1))


library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts1,
                              colData = coldata1,
                              design = ~ 1)
dds
dim(dds)
keep <- rowSums(counts(dds)) > 10
dds <- dds[keep,]
dds
dds <- DESeq(dds)
res <- results(dds)
res

str_split(rownames(res),"\\.")

which(rownames(res)=="81910.1")
which(NGS128$gene_id=="81910.1")
which(NGS128$gene_id=="81910")
NGS128[c(26614,26615),]


index<-vector()
t<-rownames(res)
for(i in 1:nrow(res)){
  I<-str_split(t[i],"\\.")[[1]][1]
  index[i]<-which(CCDS.current$gene_id==I)[1]

}

head(CCDS.current[index,])
temp<-CCDS.current[index,]



### find the protein_id according to the ccds_id
protein_id<-list()
for(i in 1:length(temp$ccds_id)){
  
  tmp<-which(seq2current$X.ccds==as.character(temp$ccds_id[i]))
  s<-as.character(seq2current$protein_ID[tmp])
  protein_id[[i]]<-s
  
  if(length(unique(s))>1){
    print(i)
  }
}
library(stringr)
diff_protein<-vector()
for(i in 1:length(protein_id)){
  diff_protein[i]<-str_c(protein_id[[i]],collapse = "\\")
}
res_temp<-cbind.data.frame(CCDS.current[index,c("gene", "gene_id","ccds_id")],res$log2FoldChange,diff_protein)
res_temp<-res_temp[order(res_temp[,"res$log2FoldChange"],decreasing = TRUE),]
colnames(res_temp)<-c("gene","gene_id","ccds_id","log2FoldChange","protein_NCBI")
dim(res_temp)

NGS128[keep,]

A<-cbind(NGS128_RPK_count$expected_count/NGS128_RPK_count$effective_length*1000,
         NGS128_RPK_count$effective_length,
         NGS128_RPK_count$expected_count,
         NGS130_RPK_count$expected_count/NGS130_RPK_count$effective_length*1000,
         NGS130_RPK_count$effective_length,
         NGS130_RPK_count$expected_count)
rownames(A)<-NGS128_RPK_count$gene_id
### RPK: divide the read counts by the length of each gene in kilobases.
colnames(A)<-c("NGS128_RPK","NGS128_effective_length","NGS128_expected_count",
               "NGS130_RPK","NGS130_effective_length","NGS130_expected_count")
all(res_temp$gene_id %in% rownames(A))

A[as.character(res_temp$gene_id),]
res_temp<-cbind(res_temp,A[as.character(res_temp$gene_id),])
colnames(res_temp)<-c("gene","gene_id","ccds_id","log2FoldChange","protein_NCBI",
                      "NGS128_RPK","NGS128_effective_length","NGS128_expected_count",
                      "NGS130_RPK","NGS130_effective_length","NGS130_expected_count")
write.csv(res_temp,file="Deseq2_differential_analysis_isoforms_include_proteinName_res.csv",row.names = FALSE)

write.csv(CCDS.current,file="CCDS.current.csv")

res_temp$gene
res_temp1$gene
length(intersect(res_temp$gene,res_temp1$gene))

setdiff(res_temp1$gene,res_temp$gene)

table(as.character(res_temp$gene))


CCDS.current$gene[table(CCDS.current$gene)==3]
