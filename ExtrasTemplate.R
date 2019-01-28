#A set of standard functions for common use
options(digits=15)

include <- function(){
	#List of all the functions and variables in this source file
	print("version 2.124")
	print("(function) read.tab")
	print("(function) write.tab")
	print("(function) write.list")
	print("(function) write.textable")
	print("(function) count")
	print("(function) se")
	print("(function) t.statistic")
	print("(function) first")
	print("(function) cohen")
	print("(function) cohen.simple")
	print("(function) hedges")
	print("(function) oneway.dunnett")
	print("(function) oneway.anova")
	print("(function) statistics")
	print("(function) group.statistics")
	print("(function) barplot.errors")
	print("(list) std.colours")
	print("(list) nice.colours")
	print("(list) colour.schemes $distinct(n, scheme) is a function generating n distict colours following a scheme in 'default', 'auto', 'pretty'") 
	print("(function) colourise.factor(factor, scheme) Turns a factor into a vector of colours, one for each level (using colour.schemes$distinct())")
#print("(vector) colour.palette")
	print("(function) pcor")
	print("(function) ndc")
	print("(function) get.entrez.symbol")
	print("(function) corner")
	print("(function) plot.text ['write.text' is deprecated]") 
	print("(function) as.clean.data.frame")
	print("(function) swarmplot")
	print("(function) mango.query")
	print("(function) roman")
	print("(function) ROMAN")
	print("(function) aba.diff.structure")	
	print("(function) GOForIt(reference.set, gene.sets, id.type='ensembl') Requires package 'topGO'")
	print("(function) explodeCluster")
	print("(function) gse")
	print("(function) tf.analysis")
	print("(function) raw.presence.filter")
	print("(function) cpm.presence.filter")
	print("(function) na.presence.filter")
	print("(function) report.anova.html")
	print("(function) tosentence")
	print("(function) df2table")
	print("(function) plot.heatmap")
	print("(function) plot.heatmap2")
	print("(function) plot.clustering")
	print("(function) colour.matrix(matrix, colour.palette)")
	print("(function) balance.range.matrix(matrix) Converts a matrix with positve and negative values into the range -1,1 retaining the balance around 0")
	print("(function) colour.scale.legend(matrix, colour.palette)")
	print("(function) extreme")
	print("(function) PCA")
	print("(function) tcor")
	print("(function) vcor")
	print("(function) ecor")
	print("(function) eqcor")
	print("(function) GO.reduce")
	print("(function) rangify")
	print("(function) rowVars")
	print("(function) colVars")
	print("(function) t.intervals")
	print("(function) variance.filter")
	print("(function) subset.even")
	print("(function) split.vector")
	print("(function) collapse.cols")
	print("(function) in2mm")
	print("(function) mm2in")
	print("(function) hit(x, y)")
	print("(function) alpha(col, alpha Colour in any format, alpha as fraction)")	
}

read.tab <- function(file,row.names=NULL, header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="", ...){
	read.table(file=file, row.names=row.names, header=header, stringsAsFactors=stringsAsFactors, sep=sep, quote=quote, ...)
}

write.tab <- function(data,file,col.names=TRUE, quote=FALSE, sep="\t", row.names=FALSE, ...){
	if(row.names!=FALSE){ 
		data <- cbind(row.names(data), data)
		colnames(data)[1] <- row.names
	}
	write.table(data, file=file, quote=quote, sep=sep, row.names=F, col.names=col.names, ...)             
}

write.list <- function(data,file,col.names=TRUE, quote=FALSE, sep="\t", row.names=FALSE, ...){
	write.table(x=apply(data,2,as.character), file=file, quote=FALSE, sep="\t", row.names=FALSE, col.names=F, ...)             
}

count <- function(x){
	sum(!is.na(x))
}

se <- function(x){
	sd(x, na.rm=TRUE) / sqrt(count(x))
}

t.statistic <- function(values, groups){
	diff(by(values, groups, mean)) / sqrt(sum( by(values, groups, var) / by(values, groups, length) ))
}

cohen.simple <- function(x, y){
	#x = data for the test group
	#y = data for the control group
	(mean(x, na.rm=T)-mean(y, na.rm=T)) / sqrt( ( var(x, na.rm=T) + var(y, na.rm=T) ) /  2 )
}

first <- function(x){
	as.character(x[1])
}

cohen <- function(x, y){
	#x = data for the test group
	#y = data for the control group
	mean.x <- mean(x, na.rm=T)
	ss.x <- sum((x - mean.x) ^ 2, na.rm=T)   
	mean.y <- mean(y, na.rm=T)
	ss.y <- sum((y - mean.y) ^ 2, na.rm=T) ^ 2  
	(mean.x - mean.y) / sqrt( (ss.x + ss.y) / (count(x)+count(y)-2) )
}

hedges <- function(x, y){
	#x = data for the test group
	#y = data for the control group
	(mean(x, na.rm=T)-mean(y, na.rm=T)) / sqrt( ( (count(x)-1)*var(x, na.rm=T) + (count(y)-1)*var(y, na.rm=T) ) /  count(x)+count(y)-2 )
}

oneway.dunnett <- function(data, factor){
	data <- as.numeric(data)
	factor <- as.factor(factor)
	m <- lm(data~factor)
	summary(multcomp::glht(m, linfct = multcomp::mcp(factor = "Dunnett")))$test
}

oneway.anova <- function(data, factor){
	m <- aov(data~factor)
	s <- summary(m)[[1]]
	return(list("statistics"=s[,"F value"][-nrow(s)], "pvalues"=s[,"Pr(>F)"][-nrow(s)], "df"=s[,"Df"]))
}

statistics<-function(values, factor, rounding){
    Mean<-round(c(by(values, factor, mean, na.rm=T)), rounding)
	SEM<-round(c(by(values, factor, se)), rounding)
	n<-round(c(by(values, factor, function(x)count(x))), rounding)
	test<-t.test(values~factor, var.equal=FALSE)
	t_test<-data.frame(df=round(test$parameter, rounding), t=round(test$statistic, 2), p=signif(test$p.value, 2))
	return(list(Group.Data=cbind(Mean, SEM, n), Welch.t=t_test))
}

group.statistics <- function(x, factor, ...){	
	result <- by(x, factor, stats, ...)
	return(data.frame(t(sapply(result, c))))
}

barplot.errors <- function(values, factors, xlab="", ylab="", las=2, error.bars="se", error.range="both", error.colour="black", error.width=0.08, error.lwd=2, col=rep("#bebebeff",length(levels(factors))), ylim=c( min(0,(min(means,na.rm=T)-(max(errors,na.rm=T)*2))),max(0,(max(means,na.rm=T)+(max(errors,na.rm=T)*2))) ), metrics=F, ...){
	if(!is.list(factors)) factors <- list(factors)
	means <- as.numeric(by(values, factors, mean, na.rm=TRUE))
	errors <- numeric()
	if(is.character(error.bars)){
		if(error.bars[1]=="sd"){errors <- by(values, factors, sd, na.rm=TRUE)
		}else if(error.bars[1]=="se"){errors <- by(values, factors, se)}
	}else if(is.numeric(error.bars)){errors = error.bars}
	levels<-list()
	for(n in 1:length(factors))	levels[[n]]<-levels(as.factor(factors[[n]]))
	level.names<-apply(rev(expand.grid(levels)), 1, paste, collapse=".")
	if(error.width=="NA") y.limits <- c(0, max(means))
	bp <- barplot(
		means,
		beside=TRUE,
		space=c(0.2,1),
		ylab=ylab,
		xlab=xlab,
		ylim=ylim,
		col=col,
		names.arg=level.names,
		las=las,
		...
	)
	if(!is.na(error.width) & length(errors)>0){
		if(error.range=="upper") arrows(bp, means+errors, bp, means, angle=90, code=1, length=error.width, col=error.colour, lwd=error.lwd)
		else if(error.range=="lower") arrows(bp, means, bp, means-errors, angle=90, code=2, length=error.width, col=error.colour, lwd=error.lwd)
		else arrows(bp, means+errors, bp, means-errors, angle=90, code=3, length=error.width, col=error.colour, lwd=error.lwd)
	}
	centre <- bp
	maxy <- means+errors
	miny <- means-errors
	if(metrics) {return(list(centre=centre, maxy=maxy, miny=miny, ylim=ylim))}
}

# As in NICE project...
nice.colours <- list("Stem cells"="#a70b02", "Type-2"="#df7b31", "Type-3"="#f1af22", "Neurons"="#3e617e", "Astrocytes"="#5016cd", "Oligodendrocytes"="#c90aac", "Microglia"="#a0a0a0")

std.colours <- list("d2"="#0169c9ff", "b6"="#d40000ff", "f1"="#484586ff", "type1"="#a60a01ff", "type2a"="#df7b31ff", "type2b"="#f19c22ff", "type3"="#f1af22ff", "earlyimmatureneuron"="#8e9834ff", "lateimmatureneuron"="#557e3eff", "neuron"="#3e617eff", "astrocyte"="#5016cdff", "oligodendrocyte"="#c90aacff", "microglia"="#a0a0a0ff", "lightred"="#f44800ff", "lightblue"="#0360f1ff", "red"="#c83b00ff", "blue"="#024fc6ff", "darkred"="#a02f00ff", "darkblue"="#023f9eff", "warmred"="#bf471bff", "warmblue"="#0c53c2ff")

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

#
#
#
#
#
#

#
#
#
#
colourise.factor = function(factor, scheme = "auto"){
	if(!is.factor(factor)) factor = factor(factor, levels = unique(factor)) # Ensure levels retain their order 
	indices = as.numeric(factor)
	return(colour.schemes$distinct(length(levels(factor)), scheme = scheme)[indices])
}

#colour.palette <- c( "#CC0C00", "#CC1800", "#CC2500", "#CC3100", "#CC3D00", "#CC4900", "#CC5600", "#CC6200", "#CC6E00", "#CC7A00", "#CC8700", "#CC9300", "#CC9F00", "#CCAB00", "#CCB800", "#CCC400", "#C8CC00", "#BCCC00", "#AFCC00", "#A3CC00", "#97CC00", "#8BCC00", "#7ECC00", "#72CC00", "#66CC00", "#5ACC00", "#4ECC00", "#41CC00", "#35CC00", "#29CC00", "#1DCC00", "#10CC00", "#04CC00", "#00CC08", "#00CC14", "#00CC21", "#00CC2D", "#00CC39", "#00CC45", "#00CC52", "#00CC5E", "#00CC6A", "#00CC76", "#00CC83", "#00CC8F", "#00CC9B", "#00CCA7", "#00CCB4", "#00CCC0", "#00CCCC", "#00C0CC", "#00B4CC", "#00A7CC", "#009BCC", "#008FCC", "#0083CC", "#0076CC", "#006ACC", "#005ECC", "#0052CC", "#0045CC", "#0039CC", "#002DCC", "#0021CC", "#0014CC", "#0008CC", "#0400CC", "#1000CC", "#1D00CC", "#2900CC", "#3500CC", "#4100CC", "#4E00CC", "#5A00CC", "#6600CC", "#7200CC", "#7E00CC", "#8B00CC", "#9700CC", "#A300CC", "#AF00CC", "#BC00CC", "#C800CC", "#CC00C4", "#CC00B8", "#CC00AB", "#CC009F", "#CC0093", "#CC0087", "#CC007A", "#CC006E", "#CC0062", "#CC0056", "#CC0049", "#CC003D", "#CC0031", "#CC0025", "#CC0018", "#CC000C", "#CC0000");

# Partial correlation coefficient
# From formulas in Sheskin, 3e
# This code from http://ww2.coastal.edu/kingw/statistics/R-tutorials/multregr.html
# a,b=variables to be correlated, c=variable to be partialled out of both
pcor <- function(a,b,c){
     (cor(a,b)-cor(a,c)*cor(b,c))/sqrt((1-cor(a,c)^2)*(1-cor(b,c)^2))
}

# Normalised difference from control
ndc <- function(x,control=null){
	return((x - control) / control)
}

#Submit an Entrez GeneID to the NCBI server to retrieve the latest Gene Symbol annotation
get.entrez.symbol<-function(geneid){
	node<-xmlRoot(xmlTreeParse(paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&retmode=xml&id=",geneid,sep="")));
	elements<-c("Entrezgene", "Entrezgene_gene", "Gene-ref", "Gene-ref_locus");
	for(element in elements){
		result<-xmlElementsByTagName(node, element);
		node<-result[[1]];
	}
	symbol<-xmlValue(node);
	print(symbol); #It takes a while, this provides feedback that the function is still working
	return(c(geneid,symbol));
}

#Get the 'corner' of a large matrix. Like head, but better.
corner<-function(m, corner="tl", n=10){
	n.rows<-min(n, nrow(m))
	n.cols<-min(n, ncol(m))	
	if(corner=="tl") return(m[1:n.rows, 1:n.cols]) #Top left
	if(corner=="tr") return(m[1:n.rows, (ncol(m)-n.cols):ncol(m)]) #Top right
	if(corner=="bl") return(m[(nrow(m)-n.rows):nrow(m), 1:n.cols]) #Bottom left
	if(corner=="br") return(m[(nrow(m)-n.rows):nrow(m), (ncol(m)-n.cols):ncol(m)]) #Bottom right
}

#Creates a plot that only contains text. Designed for use in multi-plot layouts
plot.text<-function(text, valign="centre", halign="centre", v=0, h=0, position=0, ...){
	if(v==0){
		if(valign=="top") v <- 0.9
		if(valign=="centre") v <- 0
		if(valign=="bottom") v <- -0.9
	}
	if(h==0){
		if(halign=="left"){
			h <- -1
			if(position==0) position = 4 #Default positions based on horizontal text alignment
		}
		if(halign=="centre"){
			h <- 0
			if(position==0) position = NULL
		}
		if(halign=="right"){
			h <- 1
			if(position==0) position = 2
		}
	}
	old.mar<-par("mar")
	old.cex<-par("cex")
	par( mar=c(0, 0, 0, 0), ...)
	plot(0, xlim=c(-1, 1), ylim=c(-1, 1), ann=F, xaxt="n", yaxt="n", bty="n", type="n")
	text(h, v, text, pos=position)
	par(mar=old.mar)
	par(cex=old.cex)
}
#Deprecated 
write.text <- plot.text 

#Casts a matrix into a data frame with control over the column types
as.clean.data.frame<-function(matrix, classes, row.names=rownames(matrix), col.names=colnames(matrix)){
	df<-as.data.frame(matrix)
	for(column in 1:ncol(matrix)){
		cast<-match.fun(paste("as", classes[column], sep="."))
		df[, column]<-cast(as.character(matrix[, column])) #as.character ensures factors are decoded
	}
	rownames(df) <- row.names
	colnames(df) <- col.names
	return(df)
}

swarmplot <- function(values, factors, mean.width=0.5, error.lwd=3, errors="se", pch=21, col=std.colours$red, spacing=1.5, xlab="", ...){
	if(!is.list(factors)) factors <- list(factors)
	values <- as.numeric(values)
	formula <- eval(parse(text =paste0("values~",paste(paste0("as.factor(factors[[",1:length(factors),"]])"), collapse="+"))))
	beeswarm::beeswarm(formula, pch=pch, col=col, spacing=spacing, xlab=xlab, ...)
	means <- by(values, factors, mean, na.rm=T)
	error.values <- NULL
	if(errors=="se") error.values <- as.numeric(by(values, factors, se))
	if(errors=="sd") error.values <- as.numeric(by(values, factors, sd))
	if(errors=="none") errors <- ""
	for(n in 1:length(means)){
		offset <- mean.width / length(means)
		x <- n
		y <- means[n]
		segments(x-offset, y, x+offset, y, col="black", lwd=error.lwd)
		if(errors!=""){
			segments(x, y-error.values[n], x, y+error.values[n], col="black", lwd=error.lwd)
			segments(x-offset/3, y+error.values[n], x+offset/3, y+error.values[n], col="black", lwd=error.lwd)
			segments(x-offset/3, y-error.values[n], x+offset/3, y-error.values[n], col="black", lwd=error.lwd)
		}
	}
}

mango.query <- function(user="", password="", database="mango", pmid="", geneid="", process="", cellstage="", effect="", evidence="", species="", type="", expression="true", negatives="true", print.url=FALSE){
	login <- paste(user, password, sep=":")
	if(login==":") login <- ""
	url <- paste0("http://", login, database, ".adult-neurogenesis.de/xml/annotations?pmid=", pmid, "&geneid=", geneid, "&process=", process, "&cellstage=", cellstage, "&effect=", effect, "&evidence=", evidence,  "&species=", species,  "&type=", type,  "&expression=", expression,  "&negatives=", negatives)
	if(print.url) print(url)
	xml.result = XML::xmlParse(url)   
	query.result <- XML::xmlToList(xml.result)$`query-result`
	query.result.df <- as.data.frame(t(as.data.frame(lapply(query.result, unlist))), stringsAsFactors=F)
	return(query.result.df)
}

roman <- function(numbers){
	result <- sapply(numbers, function(number){
	if(number>3999) return("NA")
	else{
		l <- rev(unlist(strsplit(as.character(number), "")))
		ones <- c("i", "ii", "iii", "iv", "v", "vi", "vii", "viii", "ix")[as.numeric(l[1])]
		tens <- c("x", "xx", "xxx", "xl", "l", "lx", "lxx", "lxxx", "xc")[as.numeric(l[2])]
		hundreds <- c("c", "cc", "ccc", "cd", "d", "dc", "dcc", "dccc", "cm")[as.numeric(l[3])]
		thousands <- c("m", "mm", "mmm")[as.numeric(l[4])]
	}
	roman <- c(thousands, hundreds, tens, ones)
	paste0(roman[!is.na(roman)], collapse="")
	})
	return(result)
}

ROMAN <- function(number) toupper(roman(number))

aba.diff.structure <- function(contrast, target, page=0){
	s1 <-  XML::xmlRoot(XML::xmlTreeParse(paste0("http://api.brain-map.org/api/v2/data/query.xml?criteria=model::Structure,rma::criteria,[acronym$il'", contrast, "'],ontology[abbreviation$eq'Mouse']")))
	s1.id <- XML::xmlValue(XML::xmlElementsByTagName(s1[[1]][[1]], "id")[[1]])
	s2 <-  XML::xmlRoot(XML::xmlTreeParse(paste0("http://api.brain-map.org/api/v2/data/query.xml?criteria=model::Structure,rma::criteria,[acronym$il'", target, "'],ontology[abbreviation$eq'Mouse']")))
	s2.id <- XML::xmlValue(XML::xmlElementsByTagName(s2[[1]][[1]], "id")[[1]])
	#In the following, threshold1 seems always to be set at 0, 50 in the ABA searches.
	aba <- XML::xmlRoot(XML::xmlTreeParse(paste0("http://api.brain-map.org/api/v2/data/query.xml?criteria=service::mouse_differential[set$eq'mouse'][structures1$eq", s1.id, "][structures2$eq", s2.id, "][threshold1$eq0,50][threshold2$eq1,50][num_rows$eq2000][start_row$eq", 2000 * page,"]")))
	genes <- XML::xmlElementsByTagName(aba[[1]], "object")
	aba.symbols <- sapply(genes, function(gene){XML::xmlValue(XML::xmlElementsByTagName(gene, "gene-symbol")[[1]])})
	aba.geneids <- sapply(genes, function(gene){XML::xmlValue(XML::xmlElementsByTagName(gene, "entrez-id")[[1]])})
	aba.fcs <- sapply(genes, function(gene){XML::xmlValue(XML::xmlElementsByTagName(gene, "fold-change")[[1]])})
	result <- cbind("symbol"=aba.symbols, "geneid"=aba.geneids, "fc"=aba.fcs)
	row.names(result) <- 1:nrow(result) 
	return(as.clean.data.frame(result, c("character", "character", "numeric")))
}

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

explodeCluster <- function(clustering, naming.scheme="numbers"){
	gene.sets <- lapply(names(table(clustering)), function(cluster){
		names(clustering)[which(clustering==cluster)] 
	})
	names(gene.sets) <- paste0("Cluster_", 1:length(gene.sets))
	if(naming.scheme=="letters"){
		names(gene.sets) <- paste0("Cluster_", letters[1:length(gene.sets)])
	}
	if(naming.scheme=="LETTERS"){
		names(gene.sets) <- paste0("Cluster_", LETTERS[1:length(gene.sets)])
	}
	if(naming.scheme=="roman"){
		names(gene.sets) <- paste0("Cluster_", roman(1:length(gene.sets)))
	}
	if(naming.scheme=="ROMAN"){
		names(gene.sets) <- paste0("Cluster_", ROMAN(1:length(gene.sets)))
	}
	return(gene.sets)
}

#MSigDB has been crunched and geneids converted (where possible) 
gse <- function(query.set, reference.set, category, p.sig=0.05){
	if(!exists("msigdb")){
		msigdb <<- read.tab('http://research.rupertoverall.net/MSIGDB.tab')
		print("msigdb loaded")
	}
	category.data <- NULL
	if(category=="all") category.data <- msigdb
	else category.data <- msigdb[which(msigdb$SubCategory==category), ]
	enrichment <- t(apply(category.data, 1, function(row){
		gene.set <- strsplit(row[["MouseGeneIDs"]], ", ", fixed=T)[[1]]		
		m <- length(intersect(gene.set, reference.set)) #Total number of possible hits in the reference population
		n <- length(reference.set) - m #The remainder of the reference population
		k <- length(query.set)
		hits <- intersect(query.set, intersect(gene.set, reference.set))
		q <- length(hits) #How many hits in the query
		p <- phyper(q, m, n, k, lower.tail=F) + dhyper(q, m, n, k)		
		return(c(row[["GeneSet"]], NA, p, row[["Category"]], row[["SubCategory"]], row[["Description"]], paste(hits, collapse=", ")))
	}))
	enrichment.df <- as.clean.data.frame(enrichment, c("character", "numeric", "numeric", "character", "character", "character", "character"), col.names=c("GeneSet", "AdjP", "RawP", "Category", "SubCategory", "Description", "Hits"))  
	enrichment.df$AdjP <- p.adjust(enrichment.df$RawP, method="fdr")
	enrichment.df <- enrichment.df[which(enrichment.df$AdjP < p.sig),]
	return(enrichment.df[order(enrichment.df$RawP), ])
}

#Transcription factor mappings for mouse promoters (mm10; +2000, -200) using motifs from MotifDb (90% similarity) for mouse, human and rat.
tf.analysis <- function(query.set, reference.set, species="all", p.sig=0.05){
	if(!exists("tf.mapping")){
		tf.mapping <<- read.tab('http://www.biotec.tu-dresden.de/~ruperto/rcourse/TF%20mapping.tab')
		print("tf.mapping loaded")
	}
	if(species=="all") category.data <- tf.mapping
	else category.data <- tf.mapping[which(tolower(tf.mapping$Species)==tolower(species)), ]	
	enrichment <- t(apply(category.data, 1, function(row){
		gene.set <- strsplit(row[["Targets"]], ", ", fixed=T)[[1]]		
		m <- length(intersect(gene.set, reference.set)) #Total number of possible hits in the reference population
		n <- length(reference.set) - m #The remainder of the reference population
		k <- length(query.set)
		hits <- intersect(query.set, intersect(gene.set, reference.set))
		q <- length(hits) #How many hits in the query
		p <- phyper(q, m, n, k, lower.tail=F) + dhyper(q, m, n, k)		
		return(c(row[["Species"]], row[["TF"]], NA, p, paste(hits, collapse=", ")))
	}))
	enrichment.df <- as.clean.data.frame(enrichment, c("character", "character", "numeric", "numeric", "character"), col.names=c("Species", "TF", "AdjP", "RawP", "Hits"))  
	enrichment.df$AdjP <- p.adjust(enrichment.df$RawP, method="fdr")
	enrichment.df <- enrichment.df[which(enrichment.df$AdjP < p.sig),]
	return(enrichment.df[order(enrichment.df$AdjP), ])
}

#Fast filtering of expression data by groups
raw.presence.filter <- function(counts, groups, nGroupsToPass = 1, fractionSamplesPerGroupToPass = 0.75, threshold = 1){
	group.indices <- lapply(unique(groups), function(group) which(groups==group) )
	rowSums(do.call("cbind", lapply(group.indices, function(i) (rowSums(counts[, i] >= threshold, na.rm=T) / length(i)) >= fractionSamplesPerGroupToPass )), na.rm=T) >= nGroupsToPass  
}

cpm.presence.filter <- function(counts, groups, nGroupsToPass = 1, fractionSamplesPerGroupToPass = 0.75, threshold = 1){
	group.indices <- lapply(unique(groups), function(group) which(groups==group) )
	rowSums(do.call("cbind", lapply(group.indices, function(i) (rowSums((counts[, i] / colSums(counts[, i], na.rm=T)[col(counts[, i])] * 1e6) >= threshold, na.rm=T) / length(i)) >= fractionSamplesPerGroupToPass )), na.rm=T) >= nGroupsToPass  
}

na.presence.filter <- function(counts, groups, nGroupsToPass = 1, fractionSamplesPerGroupToPass = 0.75){
	group.indices <- lapply(unique(groups), function(group) which(groups==group) )
	rowSums(do.call("cbind", lapply(group.indices, function(i) (rowSums(!is.na(counts[, i]), na.rm=T) / length(i)) >= fractionSamplesPerGroupToPass )), na.rm=T) >= nGroupsToPass  
}


report.anova.html <- function(x, term=1){
	a <- anova(lm(x))
	df <- a[[1]][c(term,length(a[[1]]))]
	f <- a[[4]][term]
	p <- a[[5]][term]	
	return(paste0("(<i>F</i><sub>",paste(df, collapse=", "),"</sub> = ", format(round(f,2), nsmall=2, scientific=F), "; <i>p</i> = ", sub("e", " &times 10<sup>", format(p, digits=3, scientific=T)), "</sup>)"))
}

tosentence <- function(s, split=NA) {
	s <- tolower(s)
		s <- as.character(sapply(s, function(ss){
			ss <- strsplit(ss, split)[[1]]
			substring(ss, 1, 1) <- toupper(substring(ss, 1, 1))
			if(!is.na(split)) ss <- paste(ss, collapse=split)
			return(ss)
		}))
	return(s)
}

df2table <- function(df, write.rows=F, row.header=""){
	dfl <- as.list(as.data.frame(df)) # Cast into list to avoid problems with numbers of rows (1 row -> vector)
	sanitise <- function(str) gsub('([#$%&~_\\^\\\\{}])', '\\\\\\1', str, perl = TRUE)
	names(dfl) <- sanitise(names(dfl))
	dfl <- lapply(dfl, sanitise)
	dfl <- (lapply(dfl, function(col) sub("([0-9.]+)(e)([\\+|\\-])([0-9]+)", "$\\1  \\\\times 10^{\\3\\4}$", col) )) # Scientific notation
	df <- as.data.frame(do.call("cbind", dfl))
	colnames(df) <- names(dfl)
	if(write.rows){
		row.names <- rep("", nrow(df))
		if(length(rownames(df))!=0) row.names <- sanitise(rownames(df))
		df <- cbind(row.names, df)
		colnames(df)[1] <- row.header
	}
	output <- paste0("\t\t\\begin{tabular}{", paste(rep('c', ncol(df)), collapse=''), "}\n",
						  "\t\t\t\\toprule\n",
						  paste0("\t\t\t", paste(colnames(df), collapse=" & "),"\\\\\n"),
						  "\t\t\t\\midrule\n",
						  paste0(apply(df, 1, function(row) paste0("\t\t\t", paste((row), collapse=" & "),"\\\\\n") ), collapse=""),
						  "\t\t\t\\bottomrule\n",
						  "\t\t\\end{tabular}\n")
	return(output)
}

write.textable <- function(df, file, write.rows=F){
	this.write.rows <- write.rows
	this.file <- file
	cat(df2table(df, write.rows=this.write.rows), file=this.file)
}

plot.heatmap <- function(plot.data, plot.groups=rep(1, nrow(plot.data)), label.width=1, xleft=0, ybottom=0, xright=1, ytop=1, colour.palette=colorRampPalette(c("#339900", "yellow", "#d40000"))(256), na.colour=NA, scale=FALSE, row.scale=FALSE, cex=1, overlay=F, ...){
	def.par <- par(no.readonly = TRUE)
	k <- length(unique(plot.groups))
	l <- NULL
	if(is.na(plot.groups[1]) & is.na(label.width[1])) l <- layout(mat=t(1), widths=c(10))
	else if(is.na(plot.groups[1])) l <- layout(mat=t(1:2), widths=c(label.width, 10))
	else if(is.na(label.width[1])) l <- layout(mat=t(1:2), widths=c(1, 10))
	else l <- layout(mat=t(1:3), widths=c(1, label.width, 10))
	par(mar = c(0,0,0,0))
	if(!is.na(plot.groups[1])){
		colours <-  colorRampPalette(c("white", "black"))(k)
		plot(0, 0, type="n", xlim=c(0, 1), ylim=c(0, 1), axes=F, xaxs="i", yaxs="i")
		rasterImage(colours[plot.groups], 0, 0, 1, 1, interpolate = F)
		centres <- 1 - sapply(1:k, function(n) mean(range(which(plot.groups==n)) / length(plot.groups)) - (1 / length(plot.groups) / 2) )
		text(0.5, centres, ROMAN(1:k), font=2, col=c(rep("black", floor(k/2)), rep("white", k - floor(k/2))))
	}
	if(!is.na(label.width[1])){
		plot(0, 0, type="n", xlim=c(0, 1), ylim=c(0, 1), axes=F, xaxs="i", yaxs="i")
		centres <- 1 - (1:nrow(plot.data) / nrow(plot.data)) + (1 / nrow(plot.data) / 2)
		text(0.5, centres, rownames(plot.data), cex=cex)
	}
	d <- plot.data
	min.d <- min(d, na.rm=T)
	d <- d - min.d
	max.d <- max(d, na.rm=T)
	d <- d / max.d
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
	par(mar = c(0,0,0,0))
	plot(0, 0, type="n", xlim=c(0, 1), ylim=c(0, 1), axes=F, xaxs="i", yaxs="i", new=(!overlay))
	rasterImage(colour.matrix(d, colour.palette, na.colour=na.colour), xleft, ybottom, xright, ytop, interpolate=F, ...)
	par(def.par)
}

plot.heatmap.raster <- function(plot.data, x.groups=rep(1, ncol(plot.data)), y.groups=rep(1, nrow(plot.data)), x.labels=colnames(plot.data), y.labels=rownames(plot.data), x.group.colours=colourise.factor(x.groups, scheme = "auto"), title="", title.height=2, x.label.height=2, y.label.width=4, y.group.width = 0.1, cex=1, title.cex=1, x.label.cex=cex, y.label.cex=cex, title.adjust=NA, x.label.adjust=NA, y.label.adjust=NA, las=1, main="", xleft=0, ybottom=0, xright=1, ytop=1, colour.palette=colorRampPalette(c("#339900", "yellow", "#d40000"))(256), na.colour=NA, scale=FALSE, row.scale=FALSE, absolute=FALSE, overlay=F, ...){
  # This version plots only the heatmap as a raster image and allows normal plot labelling
  # It does not yet fully support the group highlighting
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

# For backwards compatibility
plot.heatmap2 <- plot.heatmap.raster

plot.clustering <- function(dist, k, plot.data, headers=NULL, colour.palette=colorRampPalette(c("#339900", "yellow", "#d40000"))(256), na.colour=NA, labels=FALSE, label.adjust=0, scale=FALSE, row.scale=FALSE, cluster.method="complete"){
	def.par <- par(no.readonly = TRUE)
	clust <- flashClust::hclust(dist, method=cluster.method)
	dendrogram <- as.dendrogram(clust, hang=-1)
	l <- layout(mat=matrix(1:6, ncol=3, byrow=T), widths=c(6, 1, 10), heights=c(1, nrow(plot.data)))
	#Header row (first two empty)
	par(mar = c(0,0,0,0))
	plot(0, type="n", axes=F)
	plot(0, type="n", axes=F)
	
	if(length(headers)==0 | length(unique(headers))<=1){
		par(mar = c(0,0,0,0))
		plot(0, type="n", axes=F)
	}else{
		colours <- colour.matrix((headers - 1)  / (length(unique(headers)) - 1), colorRampPalette(c("white", "black"))(256))
		par(mar = c(0,0,0,0))
		plot(0, 0, type="n", xlim=c(0, 1), ylim=c(0, 1), axes=F, xaxs="i", yaxs="i")
		rasterImage(t(colours), 0, 0, 1, 1, interpolate = F)
	}
	#Plots
	if(labels){
		label.width <- (2.2 * max(nchar(rownames(plot.data))) / (72/par()$ps)) + label.adjust
		par(mar = c(0,0,0,label.width))
		plot(dendrogram, horiz=T, axes=F, yaxs="i")
	}
	if(!labels){
		par(mar = c(0,0,0,0))
		plot(rev(dendrogram), horiz=T, axes=F,yaxs="i", leaflab="none") # Reversed as horizontal orientation ends up having cluters bottom-to-top (I think this is counterintuitive)
	}
	par(mar = c(0,0,0,0))
	cluster.ids <- cutree(clust, k=k)[clust$order]
	if(length(unique(cluster.ids))==1){
		plot(0, type="n", axes=F)
	}else{
		colours <- colour.matrix((cluster.ids - 1) / (k - 1), rev=F, colorRampPalette(c("white", "black"))(256))
		plot(0, 0, type="n", xlim=c(0, 1), ylim=c(0, 1), axes=F, xaxs="i", yaxs="i")
		rasterImage(colours, 0, 0, 1, 1, interpolate = F)
		centres <- sapply(1:k, function(n) mean(range(which(cluster.ids==n)) / length(cluster.ids)) - (1 / length(cluster.ids) / 2) )
		text(0.5, centres, ROMAN(1:k), font=2, col=c(rep("black", floor(k/2)), rep("white", k - floor(k/2))))
	}
	d <- plot.data[clust$order,]
	min.d <- min(d)
	d <- d - min.d
	max.d <- max(d)
	d <- d / max.d
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
	par(mar = c(0,0,0,0))
	plot(0, 0, type="n", xlim=c(0, 1), ylim=c(0, 1), axes=F, xaxs="i", yaxs="i")
	rasterImage(colour.matrix(d, rev=F, colour.palette, na.colour=na.colour), 0, 0, 1, 1, interpolate=F)
	par(def.par)
	return(cutree(clust, k=k))
}

colour.matrix <- function(m, colour.palette=colorRampPalette(c("black", "white"))(256), rev=F, na.colour=NA){
	m <- as.matrix(m)
	cm <- matrix(colour.palette[round(m * (length(colour.palette) - 1)) + 1], nrow=nrow(m))
	cm[is.na(cm)] <- na.colour
	if(rev) cm <- as.matrix(cm[rev(1:nrow(cm)), ])
	return(cm)
}

# Converts a matrix with positve and negative values into the range -1,1 retaining the balance around 0
balance.range.matrix = function(m){
	m <- as.matrix(m)
	extreme.m = max(abs(min(m, na.rm = T)), abs(max(m, na.rm = T)), na.rm = T)
	m / extreme.m
}

extreme = function(m){
	m <- as.matrix(m)
	max(abs(min(m, na.rm = T)), abs(max(m, na.rm = T)), na.rm = T)
}

colour.scale.legend <- function(range = seq(0, 1, by=0.1), colour.palette = colour.schemes$greenred.3, cex.axis = 0.6){
	colour.scale = as.character(colour.matrix(seq(0, 1, length.out = length(range)), colour.palette=colour.palette))
	def.par <- par(no.readonly = TRUE)
	par(mar=c(2, 0, 0, 0), lwd=.4)
	plot(rep(0, length(colour.scale)), type="n", xlim=c(0, length(colour.scale)), ylim=c(-1, 1), axes=F, xlab="", ylab="")
	for(n in 1:length(colour.scale)) rect(n-1, -1, n, 1, col=colour.scale[n])
	axis(1, at=1:length(colour.scale) - .5, labels=range, tick=F, line=-1, cex.axis=cex.axis)
	try(par(def.par), silent=TRUE)
}

# PCA <- function(x){
# 	x <- as.data.frame(x)
# 	means <- apply(x, 1, mean, na.rm=T)
# 	result <- prcomp(formula(paste0("~", paste(colnames(x), collapse="+"))), data=x, na.action=na.exclude, scale=T)
# 	return(c(result, as.data.frame(result$x), "means"=list(as.numeric(means))))
# }
PCA <- function(x, scale=T){
	means <- apply(x, 1, mean, na.rm=T)
	result <- prcomp(x, scale=scale)
	result$x <- result$x * sign(cor(result$x[, "PC1"], means)) #Keep sensible orientation
	return(c(result, as.data.frame(result$x), "means"=list(as.numeric(means))))
}

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

ecor <- function(x, y=x){
	result = apply(as.matrix(x), 2, function(col.x){
		sqrt( colSums((y - col.x) ^ 2) ) 
	})
	return(result)
}

eqcor <- function(list.x, list.y=list.x, use.na=F){
  eqmatrix = sapply(list.y, function(group.y){
    sapply(list.x, function(group.x){
			score = sum(group.x==group.y, na.rm = T)
			if(use.na) score = score + sum(is.na(group.x)&is.na(group.y))
      return(score)
    })
  })
  return(eqmatrix)
}

GO.reduce <- function(go.results, n=20){
	k <- min(n, nrow(go.results))
	rownames(go.results) <- go.results[,1]
	go.cot <- sapply(rownames(go.results), function(go) GOSemSim::mgoSim(go, rownames(go.results), ont="BP", organism="mouse", measure="Wang", combine=NULL ) )
	d <- as.dist(1 - go.cot)
	cluster.index <- cutree(hclust(d), k)
	cluster.representatives <- sapply(1:k, function(cluster.id) names(which(cluster.index==cluster.id))[1] ) #First is most significant in ordered list 
	reduced.table <- as.data.frame(go.results[cluster.representatives, ])
	clustered.table <- go.results[order(cluster.index, go.results$adjP), ]
	adjP <- reduced.table[, "adjP"]
	n <- reduced.table[, "Significant"]
	term <- reduced.table[, "Term"]
	return(list("reduced.table"=reduced.table, "clustered.table"=clustered.table, "adjP"=adjP, "n"=n, "Term"=term))
}

rangify = function(x){
  # Consider the values in the vector x as lengths...
  # ...and split into a list of ranges such that a vector of the length sum(x) is split into pieces with lengths defined by x
  # e.g. x = c(2, 4, 7)
  # rangify(x) = list(1:2, 3:6, 7:13)
  starts = c(0, cumsum(x)[-length(cumsum(x))])+ 1
  ends = cumsum(x)
  apply(rbind(starts, ends), 2, function(col) col[1]:col[2])
}

rowVars <- function(x, na.rm=F) {
	# Vectorised version of variance filter
	rowSums((x - rowMeans(x, na.rm=na.rm))^2, na.rm=na.rm) / (ncol(x) - 1)
}

colVars <- function(x, na.rm=F) {
	# Vectorised version of variance filter
	colSums(t(t(x) - colMeans(x, na.rm=na.rm))^2, na.rm=na.rm) / (nrow(x) - 1)
}

t.intervals <- function(data, factor, na.rm=FALSE){
	# Vectorised calculation of t.statistic between intervals
	# The intervals will be calculated based on the ordering of the factor (which refers to how the data columns are arranged).
	if(!is.factor(factor)) stop("The 2nd parameter must be a factor.")
	if(is.vector(data)) data = t(data) # Convert a vector (a single row of data) to a 1-row matrix
	interval.results = lapply(2:(length(levels(factor))), function(i){
		A = data[, factor==levels(factor)[i-1], drop = F]
		B = data[, factor==levels(factor)[i], drop = F]
		(rowMeans(B, na.rm=na.rm) - rowMeans(A, na.rm=na.rm)) / sqrt((rowVars(A, na.rm=na.rm) / ncol(A)) + (rowVars(B, na.rm=na.rm) / ncol(B)) )
	})
	interval.results = do.call("cbind", interval.results)
	colnames(interval.results) = sapply(2:length(levels(factor)), function(i) paste(levels(factor)[i:(i-1)], collapse="-"))
	return(interval.results)
}

variance.filter <- function(m, fraction=0.05){
	# Removes the rows with the lowest variance (below the threshold given)
	m[which(rowVars(m) > sort(rowVars(m))[max(1, floor(nrow(m)*fraction))]), ]
}

seq.tri <- function(x){
	# Where x is the number rows/columns of a matrix.
	# Returns the coordinates of the lower triangle in a 2-column matrix.
	cbind("Source"=rep(1:x, (1:x)-1), "Target"=unlist(lapply(1:(x-1), seq)))
}

subset.even <- function(v, n) v[round(seq(1, length(v), (length(v)-1)/(n-1)))]

split.vector <- function(v, n){
  ncols = floor(length(v)/n)
  c(split(v[1:(n * ncols)], rep(1:ncols, each=n)), list(v[(n * ncols + 1):length(v)]))
}

collapse.cols = function(data, factor, fn){
	factor.list = factor
	if(!is.list(factor.list)) factor.list = list(factor)
	aggregated = t(aggregate(t(data), factor.list, fn))
	groups = aggregated[1:length(factor.list), , drop=FALSE]
	values = apply(as.matrix(aggregated[-(1:length(factor.list)), ]), 2, as.numeric)
	colnames(values) = apply(groups, 2, paste, collapse=".")
	return(values)
}

in2mm <- function(inches){return(25.4 * inches)}

mm2in <- function(mm){return(mm / 25.4)}

hit <- function(x, y){return(as.data.frame(do.call("rbind", lapply(1:length(x), function(xx) return(cbind(queryHits = xx, subjectHits = which(y==x[xx]))) )), stringsAsFactors = F))}
