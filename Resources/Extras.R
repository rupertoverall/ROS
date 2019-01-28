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

# This function plots only the heatmap as a raster (bitmap) image and allows normal plot labelling
plot.heatmap.raster <- function(plot.data, x.groups=rep(1, ncol(plot.data)), y.groups=rep(1, nrow(plot.data)), x.labels=colnames(plot.data), y.labels=rownames(plot.data), x.group.colours=colourise.factor(x.groups, scheme = "auto"), title="", title.height=2, x.label.height=2, y.label.width=4, y.group.width = 0.1, cex=1, title.cex=1, x.label.cex=cex, y.label.cex=cex, title.adjust=NA, x.label.adjust=NA, y.label.adjust=NA, las=1, main="", xleft=0, ybottom=0, xright=1, ytop=1, colour.palette=colorRampPalette(c("#339900", "yellow", "#d40000"))(256), na.colour=NA, scale=FALSE, row.scale=FALSE, absolute=FALSE, overlay=F, ...){
  def.par <- par(no.readonly = TRUE)
  k <- length(unique(x.groups))
  d <- plot.data
  min.d <- min(d, na.rm=T)
  if(!absolute) d <- d - min.d
  max.d <- max(d, na.rm=T)
  if(!absolute) d <- d / max.d
  if(is.numeric(scale)) scale <- ((scale - min.d) / max.d)
  if(row.scale){
    if(is.numeric(scale)){
      scale <- sapply(1:nrow(d), function(n){
        row <- d[n, ]
        min.row <- suppressWarnings(min(row, na.rm=T))
        row <- row - min.row
        max.row <- suppressWarnings(max(row, na.rm=T))
        row <- row / max.row
        return(((scale[n] - min.row) / max.row))
      })
    }
    d <- t(apply(d, 1, function(row){
      row <- row - suppressWarnings(min(row, na.rm=T)) # If all are NA, then warning (with result = Inf)
      row <- row / suppressWarnings(max(row, na.rm=T)) # If all are NA, then warning (with result = Inf)
      # If all are NA / Inf = NA so the above code will just return NA anyway
      row[is.nan(row)] <- 0.5 # If only one value is present, set it to be the midpoint
      return(row) 
    }))
  }
  if(is.numeric(scale)){
    d <- d - scale
    d <- d / suppressWarnings(max(abs(d), na.rm=T))
    d <- (d + 1) / 2
  }
  middles <- ((1:ncol(plot.data) / ncol(plot.data)) - (1 / ncol(plot.data) / 2))
  centres <- (1 - (1:nrow(plot.data) / nrow(plot.data)) + (1 / nrow(plot.data) / 2))
  par(mar = c(0, y.label.width, x.label.height + title.height, 0), las=las)
  plot(0, 0, type="n", xlim=c(0 - y.group.width, 1), ylim=c(0, 1), axes=F, xlab="", ylab="", new=(!overlay))#, ...)
  rasterImage(x.group.colours, xleft, ybottom, 0 - y.group.width, ytop, interpolate=F)
  rasterImage(colour.matrix(d, colour.palette, na.colour=na.colour), xleft, ybottom, xright, ytop, interpolate=F)
  axis(3, at=middles, labels=x.labels, cex.axis=x.label.cex, line=x.label.adjust, tick=F)
  axis(2, at=centres, labels=y.labels, cex.axis=y.label.cex, line=y.label.adjust, tick=F)
  title(main=title, line=title.adjust, cex.main=title.cex)
  try(par(def.par), silent=TRUE)
}

# Returns a vector of colours corresponding to factor levels
colourise.factor = function(factor, scheme = "auto"){
	if(!is.factor(factor)) factor = factor(factor, levels = unique(factor)) # Ensure levels retain their order 
	indices = as.numeric(factor)
	return(colour.schemes$distinct(length(levels(factor)), scheme = scheme)[indices])
}

# A legend for the heapmaps
colour.scale.legend <- function(range = seq(0, 1, by=0.1), colour.palette = colour.schemes$greenred.3, cex.axis = 0.6){
	colour.scale = as.character(colour.matrix(seq(0, 1, length.out = length(range)), colour.palette=colour.palette))
	def.par <- par(no.readonly = TRUE)
	par(mar=c(2, 0, 0, 0), lwd=.4)
	plot(rep(0, length(colour.scale)), type="n", xlim=c(0, length(colour.scale)), ylim=c(-1, 1), axes=F, xlab="", ylab="")
	for(n in 1:length(colour.scale)) rect(n-1, -1, n, 1, col=colour.scale[n])
	axis(1, at=1:length(colour.scale) - .5, labels=range, tick=F, line=-1, cex.axis=cex.axis)
	try(par(def.par), silent=TRUE)
}

# Some custom colour schemes
colour.schemes <- list(
	"warmred.2"=colorRampPalette(c("#FEF8E0", "#d40000"))(256), 
	"blue.2"=colorRampPalette(c("#FEFEFC", "#1040C0"))(256), 
	"greenred.3"=colorRampPalette(c("#339900", "yellow", "#d40000"))(256), 
	"greenwhitered.3"=colorRampPalette(c("#339900", "white", "#d40000"))(256), 
	"bluewhitered.3"=colorRampPalette(c("#2956a5", "white", "#e04006"))(256), 
	"bluetored.dark.3"=colorRampPalette(c("#003ea3", "white", "#bc3301"))(256), 
	"bluewhiteorange.3"=colorRampPalette(c("#2956a5", "#f7f7f7", "#e27324"))(256), 
	"brewer.3"=colorRampPalette(c("#998ec3", "#f7f7f7", "#f1a340"))(256), 
	"rainbow"=colorRampPalette(c("#D40000", "#FEFA00", "#339900", "#1040C0"))(256),
	"distinct" = function(n, scheme = "default"){
		distinct.colours = NULL
		if(scheme == "default"){
			distinct.colours = c("#4444FF","#FF4444","#44FF44","#DDDDDD","#CCCCFF","#FFCCCC","#FFCCFF","#88FFFF","#FFFF88","#FF88FF","#2222FF","#FF2222","#22FF22","#AAAAFF","#FFAAAA","#AAFFAA","#6666FF","#FF6666","#66FF66","#EEEEFF","#FFEEEE","#EEFFEE")[1:n]
			return(distinct.colours[1:n])
		}else if(scheme == "auto"){
			max.n = 12
			i = do.call("c", lapply((2:max.n)^2, function(x){ c(matrix(subset.even(1:(x-1), (x/2)), nrow=2, byrow = T)) / x * sum((1:max.n)^2) }))
			distinct.colours = colorRampPalette(c("#FF2000", "#2000FF", "#20FF00", "#FF2000"))(sum((1:max.n)^2))[order(i)]
			return(distinct.colours[1:n])
		}else if(scheme == "pretty"){
			max.n = 8
			i = do.call("c", lapply((2:max.n)^2, function(x){ c(matrix(subset.even(1:(x-1), (x/2)), nrow=2, byrow = T)) / x * sum((1:max.n)^2) }))
			distinct.colours = c( "FF4000", "#CC0C00", "#CC1800", "#CC2500", "#CC3100", "#CC3D00", "#CC4900", "#CC5600", "#CC6200", "#CC6E00", "#CC7A00", "#CC8700", "#CC9300", "#CC9F00", "#CCAB00", "#CCB800", "#CCC400", "#C8CC00", "#BCCC00", "#AFCC00", "#A3CC00", "#97CC00", "#8BCC00", "#7ECC00", "#72CC00", "#66CC00", "#5ACC00", "#4ECC00", "#41CC00", "#35CC00", "#29CC00", "#1DCC00", "#10CC00", "#04CC00", "#00CC08", "#00CC14", "#00CC21", "#00CC2D", "#00CC39", "#00CC45", "#00CC52", "#00CC5E", "#00CC6A", "#00CC76", "#00CC83", "#00CC8F", "#00CC9B", "#00CCA7", "#00CCB4", "#00CCC0", "#00CCCC", "#00C0CC", "#00B4CC", "#00A7CC", "#009BCC", "#008FCC", "#0083CC", "#0076CC", "#006ACC", "#005ECC", "#0052CC", "#0045CC", "#0039CC", "#002DCC", "#0021CC", "#0014CC", "#0008CC", "#0400CC", "#1000CC", "#1D00CC", "#2900CC", "#3500CC", "#4100CC", "#4E00CC", "#5A00CC", "#6600CC", "#7200CC", "#7E00CC", "#8B00CC", "#9700CC", "#A300CC", "#AF00CC", "#BC00CC", "#C800CC", "#CC00C4", "#CC00B8", "#CC00AB", "#CC009F", "#CC0093", "#CC0087", "#CC007A", "#CC006E", "#CC0062", "#CC0056", "#CC0049", "#CC003D", "#CC0031", "#CC0025", "#CC0018", "#CC000C", "#CC0000")[order(i)]
			return(distinct.colours[1:n])
		}else if(scheme == "brewer"){
			distinct.colours = c( "#a6cee3", "#ff7f00", "#b2df8a", "#6a3d9a", "#ffff99", "#1f78b4", "#fb9a99", "#b15928", "#33a02c", "#e31a1c", "#cab2d6", "#fdbf6f")
			return(distinct.colours[1:n])
		}else if(scheme == "annette"){
			distinct.colours = c( "#4675A2", "#A49967", "#EA8004", "#745F8D", "#4F8E00", "#939393", "#CCCC33", "#931100", "#66A0DC", "#1E1D1D", "#FFFFFF", "#5C5637", "#FF9200", "#005492", "#A64A44")
			return(distinct.colours[1:n])
		}else if(scheme == "rainbow"){
			distinct.colours = subset.even(colour.schemes$rainbow, n)
		}else{
			
		}
		return(distinct.colours)
	}
)

# Converts a matrix of values from 0 to 1 into colours on a suitable scale
colour.matrix <- function(m, colour.palette=colorRampPalette(c("black", "white"))(256), rev=F, na.colour=NA){
	m <- as.matrix(m)
	cm <- matrix(colour.palette[round(m * (length(colour.palette) - 1)) + 1], nrow=nrow(m))
	cm[is.na(cm)] <- na.colour
	if(rev) cm <- as.matrix(cm[rev(1:nrow(cm)), ])
	return(cm)
}

# Get n evenly-spaced values from a vector
subset.even <- function(v, n) v[round(seq(1, length(v), (length(v)-1)/(n-1)))]
