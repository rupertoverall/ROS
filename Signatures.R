source('Extras.R')

calculate.signatures = function(counts, groups, p=0.05, top=NA, filter=NA){
  # Counts must be raw reads
  counts = as.matrix(counts)
  if(!is.factor(groups)) groups = as.factor(groups) 
  # The parameter 'filter' MUST be provided (T/F) [to protect against old code which should fail if not aware of this option]
  if(is.na(filter)){ stop("The logical parameter 'filter' is missing")}
  else if(filter){
    # Keep only features which have more than 1 cpm in at least 3/4 of the samples for at least one group
    present = cpm.presence.filter(counts, groups, nGroupsToPass=1, fractionSamplesPerGroupToPass=0.75, threshold=1)
    counts = counts[present, ]
  }
  cpm = edgeR::cpm(counts, prior.count=1, log=T)
  
  # Filter to select genes that are regulated at all
  design = model.matrix(~groups)
  fit = edgeR::glmFit(estimateDisp(calcNormFactors(edgeR::DGEList(counts=counts, group=groups)), design), design)
  top.tags = NULL
  if(!is.na(top)){
  	top.tags = edgeR::topTags(edgeR::glmLRT(fit, coef=2:length(levels(groups))), n=top, p.value=1)$table
  }else{
  	top.tags = edgeR::topTags(edgeR::glmLRT(fit, coef=2:length(levels(groups))), n=nrow(counts), p.value=p)$table
  }
  
  # Enriched if _significantly_ greater than mean
  # The idea here is to protect against cases where most of the groups have high expression and only one is low. Such genes would seem to be poor candidates for being classed as 'enriched' in all the high celltypes.
  all.means = rowMeans(cpm[rownames(top.tags), ])
  enriched = lapply(levels(groups), function(group){
    group.indices = which(groups==group)
    tests = sapply(rownames(top.tags), function(i){
      if(var(cpm[i, group.indices]-all.means[i])<1e-20) 1 else unname(t.test(cpm[i, group.indices], mu=all.means[i], alternative="greater")$p.value)
    })
    names(which(tests < 0.05))
  })
  names(enriched) = levels(groups)

  # Signature gene if both enriched AND in a cluster of its own
  # Clusters are calculated by absolute t-statistic distances between groups
  signatures = lapply(levels(groups), function(group){
		is.signature = sapply(enriched[[group]], function(v){
			tmatrix = tcor(cpm[v, ], groups) # tcor() from include()
			tmatrix[is.na(tmatrix)] = 0 # Try to catch invariant groups (i.e. non-expressed)
			clusters = cutree(hclust(as.dist(abs(tmatrix))), k=2)
			cluster.sizes = sapply(clusters,function(cluster) length(which(clusters==cluster)) )
			cluster.sizes[group] == 1 # The only one in its group?
		})
		enriched[[group]][unlist(is.signature)]
  })
  names(signatures) = levels(groups)

  return(
    list("enriched"=enriched, "signatures"=signatures, "top.tags"=top.tags, "cpm"=cpm)
  )
}

collapse.replicate.names = function(rowname.ids){
  # When converting between IDs or species, 1:many mapping means there are replicate rownames - identified by trailing subindices
  # This function strips the vector of names and just returns the unique base ids
  stripped.ids = unique(sub("\\..*", "", rowname.ids))
  return(stripped.ids)
}

clean.replicate.names = function(rowname.ids, collapse=T){
  # When converting between IDs or species, 1:many mapping means there are replicate rownames - identified by trailing subindices
  # This function strips the vector of names and just returns the unique base ids
  stripped.ids = sub("\\..*", "", rowname.ids)
  if(collapse) stripped.ids = unique(stripped.ids)
  return(stripped.ids)
}

write.results.to.excel = function(input, file, id="ENSEMBL", background=NA, species="mouse"){
  # input is a list of data.frames (or matrices)
  # The first column must be a gene ID. 
  # Remaining columns will be written as metadata columns
  # 'id' is the geneid used. Either "ENSEMBL" or "ENTREZID"
  # 'species is either "mouse" or "human".
  group.names = names(input)
  workbook = XLConnect::loadWorkbook(file, create=T)
  XLConnect::createSheet(workbook, group.names) #Create empty sheets (which will be overwritten) and...
  XLConnect::removeSheet(workbook, setdiff(XLConnect::getSheets(workbook), group.names)) # ...delete other existing sheets.
  db = NULL
  if(species=="mouse") {db = org.Mm.eg.db::org.Mm.eg.db}
  if(species=="human") {db = org.Hs.eg.db::org.Hs.eg.db}
  if(id=="ENSEMBL"){
    # Gene lists
    . = lapply(group.names, function(group){
      XLConnect::createSheet(workbook, group)
      XLConnect::clearSheet(workbook, group)
      data = select(db, keys=as.matrix(input[[group]])[, 1], columns=c("SYMBOL", "ENTREZID"), "ENSEMBL")
      meta = as.data.frame(input[[group]])
      colnames(meta) = colnames(input[[group]])
      data = merge(data, meta, by=1)
      XLConnect::writeWorksheet(workbook, data, group)
      XLConnect::setColumnWidth(workbook, group, 1:ncol(data), -1)
    })
    # GO
    require(topGO) # This does not work unless namespace is loaded
    if(length(background)<2) background = keys(db, keytype="ENSEMBL")
    background.genes = select(db, keys=background, columns=c("SYMBOL", "ENTREZID"), "ENSEMBL")
    background.ENSEMBL = unique(background.genes$ENSEMBL[!is.na(background.genes$ENSEMBL)])
    go.results = lapply(group.names, function(group){
      go = GOForIt(background.ENSEMBL, collapse.replicate.names(as.matrix(input[[group]])[, 1]), id.type='ensembl', species=species)
    })
    names(go.results) = group.names
    . = lapply(group.names, function(group){
      XLConnect::createSheet(workbook, paste0(group, "_GO"))
      XLConnect::clearSheet(workbook, paste0(group, "_GO"))
      data = go.results[[group]]$result.table
      XLConnect::writeWorksheet(workbook, data, paste0(group, "_GO"))
      XLConnect::setColumnWidth(workbook, paste0(group, "_GO"), 1:ncol(data), -1)
    })
    # Reactome
    background.ENTREZID = unique(background.genes$ENTREZID[!is.na(background.genes$ENTREZID)])
    reactome.results = sapply(group.names, function(group){
      query = select(db, keys=collapse.replicate.names(as.matrix(input[[group]])[, 1]), columns=c("ENTREZID"), "ENSEMBL")
      query = unique(query$ENTREZID[!is.na(query$ENTREZID)])
      reactome.results = ReactomePA::enrichPathway(query, universe=background.ENTREZID, organism="mouse", pvalueCutoff=0.05, pAdjustMethod = "BH")
    })
    names(reactome.results) = group.names
    . = lapply(group.names, function(group){
      XLConnect::createSheet(workbook, paste0(group, "_Reactome"))
      XLConnect::clearSheet(workbook, paste0(group, "_Reactome"))
      data = as.data.frame(reactome.results[[group]])
      XLConnect::writeWorksheet(workbook, data, paste0(group, "_Reactome"))
      XLConnect::setColumnWidth(workbook, paste0(group, "_Reactome"), 1:ncol(data), -1)
    })
    XLConnect::saveWorkbook(workbook)
  }
  if(id=="ENTREZID"){
    # Gene lists
    . = lapply(group.names, function(group){
      createSheet(workbook, group)
      clearSheet(workbook, group)
      data = select(db, keys=as.matrix(input[[group]])[, 1], columns=c("SYMBOL", "ENSEMBL"), "ENTREZID")
      meta = as.data.frame(input[[group]])
      colnames(meta) = colnames(input[[group]])
      data = merge(data, meta, by=1)
      writeWorksheet(workbook, data, group)
      setColumnWidth(workbook, group, 1:ncol(data), -1)
    })
    # GO
    require(topGO) # This does not work unless namespace is loaded
    if(length(background)<2) background = keys(db, keytype="ENTREZID")
    background.ENTREZID = unique(background[!is.na(background)])
    go.results = lapply(group.names, function(group){
      go = GOForIt(background.ENTREZID, collapse.replicate.names(as.matrix(input[[group]])[, 1]), id.type='entrez', species=species)
    })
    names(go.results) = group.names
    . = lapply(group.names, function(group){
      createSheet(workbook, paste0(group, "_GO"))
      clearSheet(workbook, paste0(group, "_GO"))
      data = go.results[[group]]$result.table
      writeWorksheet(workbook, data, paste0(group, "_GO"))
      setColumnWidth(workbook, paste0(group, "_GO"), 1:ncol(data), -1)
    })
    # Reactome
    reactome.results = sapply(group.names, function(group){
      query = collapse.replicate.names(as.matrix(input[[group]])[, 1])
      query = unique(query[!is.na(query)])
      reactome.results = ReactomePA::enrichPathway(query, universe=background, organism="mouse", pvalueCutoff=0.05, pAdjustMethod = "BH")
    })
    names(reactome.results) = group.names
    . = lapply(group.names, function(group){
      createSheet(workbook, paste0(group, "_Reactome"))
      clearSheet(workbook, paste0(group, "_Reactome"))
      data = as.data.frame(reactome.results[[group]])
      writeWorksheet(workbook, data, paste0(group, "_Reactome"))
      setColumnWidth(workbook, paste0(group, "_Reactome"), 1:ncol(data), -1)
    })
    saveWorkbook(workbook)
  }
}

##########
