# A wrapper around prcomp to ensure the sign of the 1st PC matches the raw data and to provide principal components in a more accessible way
PCA <- function(x, scale=T){
	means <- apply(x, 1, mean, na.rm=T)
	result <- prcomp(x, scale=scale)
	result$x <- result$x * sign(cor(result$x[, "PC1"], means)) #Keep sensible orientation
	return(c(result, as.data.frame(result$x), "means"=list(as.numeric(means))))
}

# Correlation based on similarity of t-statistics
tcor <- function(variable, groups){
	tmatrix = sapply(unique(groups), function(group.y){
		sapply(unique(groups), function(group.x){
			data.x = variable[which(groups==group.x)]
			data.y = variable[which(groups==group.y)]
			if(var(data.x) < 1e-20) data.x = jitter(data.x, amount=1)
			if(var(data.y) < 1e-20) data.y = jitter(data.y, amount=1)
			t.test(data.x, data.y)$statistic
		})
	})
	colnames(tmatrix) = unique(groups)
	rownames(tmatrix) = unique(groups)
	return(tmatrix)
}

# Correlation based on the number of common elements in two lists
vcor <- function(list.x, list.y=list.x, normalise=F){
  vmatrix = sapply(list.y, function(group.y){
    sapply(list.x, function(group.x){
      score = length(intersect(group.x, group.y))
      if(normalise) score = score / min(length(group.x), length(group.y))
      if(length(group.x)==0 | length(group.y)==0) score = NA
      return(score)
    })
  })
  return(vmatrix)
}

# Convert mm into inches
mm2in <- function(mm){return(mm / 25.4)}

# Vectorised version of variance filter
colVars <- function(x, na.rm=F) {
	colSums(t(t(x) - colMeans(x, na.rm=na.rm))^2, na.rm=na.rm) / (nrow(x) - 1)
}

# Filter raw RNA sequencing data to exclude very low-expressed (absent) genes
cpm.presence.filter <- function(counts, groups, nGroupsToPass = 1, fractionSamplesPerGroupToPass = 0.75, threshold = 1){
	group.indices <- lapply(unique(groups), function(group) which(groups==group) )
	rowSums(do.call("cbind", lapply(group.indices, function(i) (rowSums((counts[, i] / colSums(counts[, i], na.rm=T)[col(counts[, i])] * 1e6) >= threshold, na.rm=T) / length(i)) >= fractionSamplesPerGroupToPass )), na.rm=T) >= nGroupsToPass  
}

# A wrapper around the topGO package to just do enrichment (without the differential expression interface of the original code)
GOForIt <- function(reference.set, query.set, id.type="ensembl", species="mouse", min.node.size=10, use.ontology="BP", algorithm="weight01", filter="none"){
	reference <- rep(0, length(reference.set))
	names(reference) <- reference.set
	reference[query.set] <- 1
	getSet <- function(set){return(set==1)}
	GOdata = NULL
	if(species=="mouse"){
	GOdata <- new(
		"topGOdata",
		description = "Simple session", 
		ontology = use.ontology,
		allGenes = reference, 
		geneSelectionFun = getSet,
		nodeSize = min.node.size,
		annotationFun = annFUN.org, 
		mapping = "org.Mm.eg.db",
		ID = id.type	)
	}
	if(species=="human"){
	GOdata <- new(
		"topGOdata",
		description = "Simple session", 
		ontology = use.ontology,
		allGenes = reference, 
		geneSelectionFun = getSet,
		nodeSize = min.node.size,
		annotationFun = annFUN.org, 
		mapping = "org.Hs.eg.db",
		ID = id.type	)
	}
	resultFisher <- topGO::runTest(GOdata, algorithm=algorithm, statistic="fisher")
	filter.score <- p.adjust(score(resultFisher), method=filter)
	adj.score <- p.adjust(score(resultFisher), method="fdr")
	top.nodes <- which(filter.score < 0.05)
	result <- "Enrichment failed"
	if(length(top.nodes > 0)){
		result <- cbind(topGO::GenTable(GOdata, "Fisher"=resultFisher, topNodes=length(top.nodes), numChar=1000), "adjP"=signif(adj.score[top.nodes][order(adj.score[top.nodes])],2))
	}
	return(list("GOdata"=GOdata, "score"=filter.score, "result.table"=result))
}
