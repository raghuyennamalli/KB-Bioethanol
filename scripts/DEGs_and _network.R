######################################################
#boxplot getting outliers in dataset
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

######################################################
#Sample-Sample distance heatmap
library("RColorBrewer")

vsd <- vst(dds,blind=FALSE)

sample_dists <- dist(t(assay(vsd)))
dist_matrix <- as.matrix(sample_dists)

rownames(dist_matrix) <- paste(vsd$condition, sep="-")
colnames(dist_matrix) <- NULL

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(dist_matrix,
         clustering_distance_rows=sample_dists,
         clustering_distance_cols=sample_dists,
         col=colors)
         
######################################################
#rawcounts to tpms conversion
geneids = row.names(data)
library(mygene)
gene.info = getGenes(geneids, fields = "all")
len = (gene.info$genomic_pos.end - gene.info$genomic_pos.start)

tpm3 <- function(data,len) {
  x <- data/len
  return(t(t(x)*1e6/colSums(x)))
}

tpm = tpm3(data,len)
rownames(tpm) = NULL

data1 = data
tpm = data.frame(tpm)
data1$ara_1_tpm = tpm$ara_1
data1$ara_2_tpm = tpm$ara_2
data1$glu_1_tpm = tpm$glu_1
data1$glu_2_tpm = tpm$glu_2
data1$xyl_tpm = tpm$xyl

write.csv(data1, "tpm_counts.csv")

######################################################
#Running DESeq2
library("DESeq2")

#loading count matrix
data = read.csv("data.csv", header=T, row.counts=1)
info = read.table("info.txt", header=T, sep="\t")
dds = DESeqDataSetFromMatrix(data, info, ~condition)

#removing low-expressed genes
keep = rowSums(counts(dds))>= 10
dds = dds[keep,]

#running the DESeq main command
ddsDE = DESeq(dds)

#normalised read counts
normCounts = counts(ddsDE, normalized=T)
write.csv(normCounts,"ddsDE.csv")

#DESeq results
res = results(ddsDE, alpha = 0.01)

#ordering the DESeq results by adjusted p-value
resOrdered = res[order(res$padj),]
write.csv(resOrdered,"DESeq_results_ordered.csv")

#DEGs-method1
'''res$sig = ifelse(res$padj<=0.01,"yes","no")
res = na.omit(res)
result = res[res$sig == 'yes',]'''

#DEGs-method2
resSig = res[ res$padj < 0.01,]

#DESeq analysis plots

#MAplot
ggplot(res,aes(x = log10(basemean), y=log2FoldChange, color=sig)) + geom_point()

#volcanoplot
ggplot(res,aes(x=log2FoldChange, y = -log10(padj), color = sig)) + geom_point()

#heatmap
normCounts <- read.csv("normalized_counts.csv",row.names=1)
sig <- subset(res, padj <= 0.05)
allsig <- merge(normCounts,sig,by = 0)
sigcounts <- allsig[,2:7]
pheatmap(log2(sigcounts+1),scale = 'row',show_rownames = F, treeheight_row = 0, treeheight_col = 0)

######################################################
#Heatmap of the counts
ntd <- normTransform(dds)
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df)

######################################################
#log2() trans on counts data
sig_genes <- rownames(result)
counts <- read.csv("data.csv",header=T,row.names=1)
log_counts <- log2(counts+1)
log_counts <- log_counts[rownames(log_counts) %in% sig_genes,]

#######################################################
#correlation matrix
cordist <- function(dat){
	cor_matrix <- cor(t(dat))
	dist_matrix <- as.matrix(dist(dat,diag=TRUE,upper=TRUE))
	dist_matrix <- log1p(dist_matrix)
	dist_matrix <- 1- (dist_matrix/max(dist_matrix))

	sign(cor_matrix) * ((abs(cor_matrix)+dist_matrix)/2)
}

#######################################################
#similarity matrix
sim_matrix <- cordist(log_counts)
heatmap.2(t(sim_matrix),col =redgreen(75),labRow=NA,labcol=NA,trace='none',dendrogram='row',xlab='Gene',ylab='Gene',main='Similarity matrix', density.info='none',revC=T)

#######################################################
#Adjacency matrix
adj_matrix <- adjacency.fromSimilarity(sim_matrix,power=6,type='signed')
gene_ids = row_names(adj_matrix)
adj_matrix <- matrix(adj_matrix,nrow=nrow(adj_matrix))
rownames(adj_matrix) = gene_ids
colnames(adj_matrix) = gene_ids
heatmap.2(t(adj_matrix),col =redgreen(75),labRow=NA,labcol=NA,trace='none',dendrogram='row',xlab='Gene',ylab='Gene',main='Adjecency matrix', density.info='none',revC=T)

#######################################################
#converting na to geneid for genesymbol
symbol = info$symbol
symbol[is.na(symbol)] = genes$g

#Network generation and export
gene_tree <- hclust(as.dist(1 - adj_matrix),method="average")
module_labels <- cutreeDynamicTree(dendro=gene_tree,minModuleSize=15,deepSplit=TRUE)
module_colors <- labels2colors(module_labels)
gene_info = data.frame(gene_ids)
gene_info <- cbind(gene_info,module=module_colors)
gene_info$color_rgb < col2hex(gene_info$module)

export_network_to_graphml <- 
function (adj_mat, filename=NULL, weighted=TRUE, threshold=0.5, max_edge_ratio=3, nodeAttr=NULL, nodeAttrDataFrame=NULL, edgeAttributes=NULL, verbose=FALSE) 
{
  library('igraph')
  
  if (is.null(filename)) 
  	{
    filename='network.graphml'
  	}
  
  max_edges <- max_edge_ratio * nrow(adj_mat)
  
  edge_to_total_ratio <- max_edges / length(adj_mat)
  edge_limit_cutoff <- as.numeric(quantile(abs(adj_mat), 1 - edge_to_total_ratio))
  
  min_threshold <- as.numeric(quantile(abs(adj_mat), 0.9999))
  
  threshold <- min(min_threshold, max(threshold, edge_limit_cutoff))
  
  adj_mat[abs(adj_mat) < threshold] <- 0
  
  orphaned <- (colSums(adj_mat) == 0)
  adj_mat <- adj_mat[!orphaned, !orphaned]
  
  if (!is.null(nodeAttr)) 
  	{
    nodeAttr <- nodeAttr[!orphaned]
  	}
  if (!is.null(nodeAttrDataFrame)) 
  	{
    nodeAttrDataFrame <- nodeAttrDataFrame[!orphaned,]
  	}
  
  is_zero     <- adj_mat == 0
  is_negative <- adj_mat < 0
  
  adj_mat <- (abs(adj_mat) - threshold) / (max(adj_mat) - threshold)
  adj_mat[is_zero] <- 0
  adj_mat[is_negative] <- -adj_mat[is_negative]
  
  if (verbose) 
  	{
    message(sprintf("Outputting matrix with %d nodes and %d edges", 
                    nrow(adj_mat), sum(adj_mat > 0)))
  	}
  

  if (weighted) 
  	{
    g <- graph.adjacency(adj_mat, mode='undirected', weighted=TRUE, diag=FALSE)
  	} 
  else 
  	{
    adj_mat[adj_mat != 0] <- 1
    g <- graph.adjacency(adj_mat, mode='undirected', diag=FALSE)
  	}
  

  if (!is.null(nodeAttr)) 
  	{
    g <- set.vertex.attribute(g, "attr", value=nodeAttr)
  	}
  

  if (!is.null(nodeAttrDataFrame)) 
  	{
    for (colname in colnames(nodeAttrDataFrame)) 
    		{
      g <- set.vertex.attribute(g, colname, value=nodeAttrDataFrame[,colname])
    		}
  	}
  
  edge_correlation_negative <- c()
  
  edge_list <- get.edgelist(g)
  
  for (i in 1:nrow(edge_list)) 
  	{
    from <- edge_list[i, 1]    
    to   <- edge_list[i, 2]    
  	}
  
  write.graph(g, filename, format='graphml')
  return(g)
}

net <- export_network_to_graphml(adj_matrix,filename='5PercentNetwork.graphml',threshold=0.4,nodeAttrDataFrame=gene_info)
##########################################