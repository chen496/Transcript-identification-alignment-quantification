setwd("/media/root/10AF08F010AF08F03/root/WeiLi_cDNA_analysis")
CCDS.current<-read.delim(file="CCDS.current.txt")
dim(CCDS.current)
head(CCDS.current)
ccds_ID<-paste0(CCDS.current$ccds_id,"|Mm108|chr",CCDS.current$X.chromosome)
rownames(CCDS.current)<-ccds_ID
head(CCDS.current)



RSEM_NGS128<-read.table("RSEM_result/NGS128_S21.isoforms.results",header = TRUE)
RSEM_NGS130<-read.table("RSEM_result/NGS130_S23.isoforms.results",header = TRUE)
dim(RSEM_NGS130)
head(RSEM_NGS130)
RSEM_NGS130_nonzero<-RSEM_NGS130[RSEM_NGS130$expected_count!=0,]
RSEM_NGS128_nonzero<-RSEM_NGS128[RSEM_NGS128$expected_count!=0,]
dim(RSEM_NGS130_nonzero)
dim(RSEM_NGS128_nonzero)

l<-union(RSEM_NGS130_nonzero$transcript_id,RSEM_NGS128_nonzero$transcript_id)
temp<-RSEM_NGS128
rownames(temp)<-RSEM_NGS128$transcript_id
write.csv(cbind(temp[l,],CCDS.current[l,c(3,4,5)]),file="total_result/RSEM_NGS128.csv",row.names = FALSE)

temp<-RSEM_NGS130
rownames(temp)<-RSEM_NGS130$transcript_id
write.csv(cbind(temp[l,],CCDS.current[l,c(3,4,5)]),file="total_result/RSEM_NGS130.csv",row.names = FALSE)


NGS128<-read.csv("total_result/RSEM_NGS128.csv")
NGS130<-read.csv("total_result/RSEM_NGS130.csv")
all((NGS128$transcript_id ) == (NGS130$transcript_id ))
cts1<-as.matrix(cbind(as.integer(NGS128$expected_count),
                      as.integer(NGS130$expected_count)))
rownames(cts1)<-NGS128$transcript_id 
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
write.csv(cbind(res,CCDS.current[rownames(res),c(3,4,5)]),file="total_result/Deseq_RSEM.csv")









Kallisto_NGS128<-read.table("Kallisto_result/output_NGS128/abundance.tsv",header = TRUE)
Kallisto_NGS130<-read.table("Kallisto_result/output_NGS130/abundance.tsv",header = TRUE)
dim(Kallisto_NGS130)
head(Kallisto_NGS130)
Kallisto_NGS130_nonzero<-Kallisto_NGS130[Kallisto_NGS130$est_counts!=0,]
Kallisto_NGS128_nonzero<-Kallisto_NGS128[Kallisto_NGS128$est_counts!=0,]

length(union(Kallisto_NGS130_nonzero$target_id,Kallisto_NGS128_nonzero$target_id))
l<-union(Kallisto_NGS130_nonzero$target_id,Kallisto_NGS128_nonzero$target_id)
temp<-Kallisto_NGS128
rownames(temp)<-Kallisto_NGS128$target_id
write.csv(cbind(temp[l,],CCDS.current[l,c(3,4,5)]),file="total_result/Kallisto_NGS128.csv",row.names = FALSE)




temp<-Kallisto_NGS130
rownames(temp)<-Kallisto_NGS130$target_id
write.csv(cbind(temp[l,],CCDS.current[l,c(3,4,5)]),file="total_result/Kallisto_NGS130.csv",row.names = FALSE)

NGS128<-read.csv("total_result/Kallisto_NGS128.csv")
NGS130<-read.csv("total_result/Kallisto_NGS130.csv")
head(NGS128)
all((NGS128$target_id ) == (NGS130$target_id))
cts1<-as.matrix(cbind(as.integer(NGS128$est_counts),
                      as.integer(NGS130$est_counts)))
rownames(cts1)<-NGS128$target_id 
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
dim(res)
write.csv(cbind(res,CCDS.current[rownames(res),c(3,4,5)]),file="total_result/Deseq_Kallisto.csv")









Salmon_NGS128<-read.table("Salmon_result/quants_NGS128/quant.sf",header = TRUE)
Salmon_NGS130<-read.table("Salmon_result/quants_NGS130/quant.sf",header = TRUE)
dim(Salmon_NGS130)
head(Salmon_NGS130)
Salmon_NGS130_nonzero<-Salmon_NGS130[Salmon_NGS130$NumReads!=0,]
Salmon_NGS128_nonzero<-Salmon_NGS128[Salmon_NGS128$NumReads!=0,]
dim(Salmon_NGS130_nonzero)
dim(Salmon_NGS128_nonzero)

length(union(Salmon_NGS130_nonzero$Name,Salmon_NGS128_nonzero$Name))
l<-union(Salmon_NGS130_nonzero$Name,Salmon_NGS128_nonzero$Name)
temp<-Salmon_NGS128
rownames(temp)<-Salmon_NGS128$Name
write.csv(cbind(temp[l,],CCDS.current[l,c(3,4,5)]),file="total_result/Salmon_NGS128.csv",row.names = FALSE)

temp<-Salmon_NGS130
rownames(temp)<-Salmon_NGS130$Name
write.csv(cbind(temp[l,],CCDS.current[l,c(3,4,5)]),file="total_result/Salmon_NGS130.csv",row.names = FALSE)


NGS128<-read.csv("total_result/Salmon_NGS128.csv")
NGS130<-read.csv("total_result/Salmon_NGS130.csv")
head(NGS128)
all((NGS128$Name) == (NGS130$Name))
cts1<-as.matrix(cbind(as.integer(NGS128$NumReads),
                      as.integer(NGS130$NumReads)))
rownames(cts1)<-NGS128$Name
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
dim(res)
write.csv(cbind(res,CCDS.current[rownames(res),c(3,4,5)]),file="total_result/Deseq_Salmon.csv")





length(intersect(RSEM_NGS130_nonzero$transcript_id,Salmon_NGS130_nonzero$Name))
length(intersect(Kallisto_NGS130_nonzero$target_id,Salmon_NGS130_nonzero$Name))


length(intersect(RSEM_NGS130_nonzero$transcript_id,Salmon_NGS130_nonzero$Name))
length(intersect(intersect(Kallisto_NGS130_nonzero$target_id,Salmon_NGS130_nonzero$Name),RSEM_NGS130_nonzero$transcript_id))



RSEM<-read.csv("total_result/Deseq_RSEM.csv",row.names = 1)
Kallisto<-read.csv("total_result/Deseq_Kallisto.csv",row.names = 1)
Salmon<-read.csv("total_result/Deseq_Salmon.csv",row.names = 1)

l<-intersect(intersect(rownames(RSEM),rownames(Kallisto)),rownames(Salmon))

temp<-cbind(RSEM[l,]$log2FoldChange,Kallisto[l,]$log2FoldChange,Salmon[l,]$log2FoldChange)
rownames(temp)<-l
colnames(temp)<-c("RSEM_log2FoldChange","Kallisto_log2FoldChange","Salmon_log2FoldChange")
write.csv(cbind(temp,CCDS.current[l,c(3,4,5)]),"total_result/DESeq_intersect_three_methods.csv")

