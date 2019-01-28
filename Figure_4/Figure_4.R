source('Extras.R')

print(load("Data/SVZ.RData")) # "svz"
print(load("Data/ROS.RData")) # "ros"

# Additional plot colours
green = "#5E9B3A"
blue = "#3346A8"
red = "#A70B02"

## PCA plot
pca = PCA(t(ros$cpm))
var.explained = round((pca$sdev)^2 / sum(pca$sdev^2) * 100)
colour.factor = sapply(ros$Group, function(name) ros$group.colours[name] )
pdf(file="Figure_4/ROS QCPrincipalComponentPlots.pdf", width=5, height=5)
plot(pca$PC1, pca$PC2, pch=20, cex=1.5, xlim=c(-130, 130), ylim=c(-130, 130), col=colour.factor, xlab=paste0("PC1 (", var.explained[1], " %)"), ylab=paste0("PC2 (", var.explained[2], " %)"))
legend("topright", c("hiROS", "midROS", "loROS"), col=ros$group.colours, ncol=1, pch=19, lty=0, text.col="black")
title(main="Principal component clustering of samples")
dev.off()

## Venn diagramme
pdf(file="Figure_4/ROS Venn.pdf", width=mm2in(85), height=mm2in(85))
plot(0, type="n", bty="n", axes=F, xlab="", ylab="") # (Clear plot for stupid Venn code)
# Enriched overlap
VennDiagram::draw.triple.venn(
  area1=length(ros$enriched$`+++`), 
  area2=length(ros$enriched$`++`), 
  area3=length(ros$enriched$`+`), 
  n12=length(intersect(ros$enriched$`+++`, ros$enriched$`++`)), 
  n23=length(intersect(ros$enriched$`++`, ros$enriched$`+`)), 
  n13=length(intersect(ros$enriched$`+++`, ros$enriched$`+`)), 
  n123=length(intersect(ros$enriched$`+++`, intersect(ros$enriched$`++`, ros$enriched$`+`))), 
  category=names(ros$enriched), cex=0.7, cat.cex=0.7, cat.fontfamily=rep("sans", 3), fontfamily="sans", lty="blank", fill=ros$group.colours)
dev.off()

venn.numbers = function(){
	hi=length(ros$enriched$`+++`) 
	mid=length(ros$enriched$`++`) 
	lo=length(ros$enriched$`+`)
	hi.mid=length(intersect(ros$enriched$`+++`, ros$enriched$`++`))
	lo.mid=length(intersect(ros$enriched$`++`, ros$enriched$`+`))
	hi.lo=length(intersect(ros$enriched$`+++`, ros$enriched$`+`))
	hi.mid.lo=length(intersect(ros$enriched$`+++`, intersect(ros$enriched$`++`, ros$enriched$`+`)))
	return(t(t(c(
		hi.only = hi - (hi.mid + hi.lo + hi.mid.lo), 
		mid.only = mid - (hi.mid + lo.mid + hi.mid.lo),
		lo.only = mid - (hi.lo + lo.mid + hi.mid.lo),
		hi.mid.only = hi.mid - hi.mid.lo,
		lo.mid.only = lo.mid - hi.mid.lo,
		hi.lo.only = hi.lo - hi.mid.lo,
		hi.mid.lo.only = hi.mid.lo
	))))
}
venn.numbers()

## Functional annotation (run once)
#$ source('Signatures.R')
#$ write.results.to.excel(ros$signatures, file="Figure_4/SignatureEnrichment.xlsx")

## Transcript overlap
pdf("Figure_4/SVZ overlaps.pdf", width=mm2in(85), height=mm2in(85))
par(mfrow=c(1, 1), cex=0.5)
intersections.ros = vcor(ros$signatures, svz$signatures, normalise=F)
intersections.ros = intersections.ros / unlist(lapply(ros$signatures, length)) # Normalised to SVZ region totals
barplot(intersections.ros, col=ros$group.colours, beside=TRUE, main="Overlap of genes of each ROS class in hippocampus regions")
dev.off()

## Cellular process markers
# Hand-curated lists
curated.markers = list(
	cell.cycle.markers = unique(data.frame(readxl::read_excel('Resources/Gene_lists.xlsx', sheet="cell_cycle_markers"), stringsAsFactors=F)[, "Ensembl.ID"]),
	ng.markers = unique(data.frame(readxl::read_excel('Resources/Gene_lists.xlsx', sheet="NG_markers"), stringsAsFactors=F)[, "Ensembl.ID"]),
	cell.cycle = unique(data.frame(readxl::read_excel('Resources/Gene_lists.xlsx', sheet="cell_cycle"), stringsAsFactors=F)[, "Ensembl.ID"]),
	ng.trajectory = unique(data.frame(readxl::read_excel('Resources/Gene_lists.xlsx', sheet="NG_trajectory"), stringsAsFactors=F)[, "Ensembl.ID"]),
	ros.markers = unique(data.frame(readxl::read_excel('Resources/Gene_lists.xlsx', sheet="ROS_markers"), stringsAsFactors=F)[, "Ensembl.ID"]),
	other.cell.cycle = unique(data.frame(readxl::read_excel('Resources/Gene_lists.xlsx', sheet="other_cell_cycle"), stringsAsFactors=F)[, "Ensembl.ID"]),
	forkhead = unique(data.frame(readxl::read_excel('Resources/Gene_lists.xlsx', sheet="forkhead"), stringsAsFactors=F)[, "Ensembl.ID"]),
	tf.genes.q = unique(data.frame(readxl::read_excel('Resources/Gene_lists.xlsx', sheet="TFs_waterfall_upinQ"), stringsAsFactors=F)[, "Ensembl.ID"]),
	tf.genes.a = unique(data.frame(readxl::read_excel('Resources/Gene_lists.xlsx', sheet="TFs_waterfall_upinA"), stringsAsFactors=F)[, "Ensembl.ID"])
)
curated.markers.names = c(cell.cycle.markers="Cell cycle markers", ng.markers="Neurogenesis markers", cell.cycle="Cell cycle genes", ng.trajectory="Neurogenic trajectory", ros.markers="ROS_markers", other.cell.cycle="Other cell cycle genes", forkhead="Forkhead family", tf.genes.q="TFs_waterfall_upinQ", tf.genes.a="TFs_waterfall_upinA")
curated.markers.pc1 = lapply(names(curated.markers), function(name) PCA(t(ros$cpm[intersect(rownames(ros$cpm), curated.markers[[name]]), ]))$PC1 )
names(curated.markers.pc1) = names(curated.markers)
#
pdf("Figure_4/Curated_lists.pdf", width=mm2in(85), height=mm2in(85))
select.markers = curated.markers[c("ros.markers", "cell.cycle.markers", "ng.markers")]
select.markers.pc1 = curated.markers.pc1[names(select.markers)]
plot(1:length(levels(ros$Group)), type="n", ylim=range(do.call("rbind", select.markers.pc1)), ylab = "arb. unit")
select.markers.colours = c(ros.markers=red, cell.cycle.markers=green, ng.markers=blue)
for(name in names(select.markers)) lines(by(select.markers.pc1[[name]], ros$Group, mean), col=select.markers.colours[name], lwd=2, cex=.6)
dev.off()

## Transcription factors from Shin et al.
pdf("Figure_4/TF_lists.pdf", width=mm2in(85), height=mm2in(85))
tf.markers = curated.markers[grep("tf\\.genes", names(curated.markers))]
tf.markers.pc1 = curated.markers.pc1[names(tf.markers)]
plot(1:length(levels(ros$Group)), type="n", ylim=range(do.call("rbind", tf.markers.pc1)), ylab = "arb. unit")
tf.markers.colours = c(tf.genes.q=red, tf.genes.a=green)
for(name in names(tf.markers)) lines(by(tf.markers.pc1[[name]], ros$Group, mean), col=tf.markers.colours[name], lwd=2, cex=.6)
dev.off()

## Expression data from Shin et al. (pseudotime data downloaded directly from the Supplementary Data provided with the Cell paper)
if(!file.exists("Resources/Shin.xlsx")) download.file('https://www.cell.com/cms/10.1016/j.stem.2015.07.013/attachment/e662afe8-70fe-4b31-87bc-d10b75ad0d95/mmc7.xlsx', "Resources/Shin.xlsx")
shin.expression = data.frame(readxl::read_excel("Resources/Shin.xlsx", skip = 2), row.names = 1, stringsAsFactors=F)
shin.pseudotime = shin.expression["Pseudotime", ] # This is the x-position in pseudotime
shin.expression = shin.expression[-1, ] # Remaining rows are expression data
shin.expression = shin.expression[rowSums(shin.expression) > 0, ] # Exclude those genes not expressed in any cell

# Convert gene symbols to ENSEMBL IDs to match our ROS lists
shin.mapping = AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=rownames(shin.expression), columns=c("ENSEMBL"), "SYMBOL")
shin.mapping = shin.mapping[match(rownames(shin.expression), shin.mapping$SYMBOL), ]
shin.mapping = shin.mapping[!is.na(shin.mapping$ENSEMBL), ]
shin.common = shin.expression[shin.mapping$SYMBOL, ]

# ROS clusters projected onto Shin
ros.pc1 = lapply(levels(ros$Group), function(group){
	common = shin.mapping$ENSEMBL %in% ros$signatures[[group]]
	data = t(shin.expression[common, ])
	data = data[, which(colVars(data, na.rm = T) > 0)] # Invariant variables will crash prcomp
	print(ncol(data))
	PCA(data)$PC1
})
names(ros.pc1) = levels(ros$Group)
#
pdf(paste0("Figure_4/Pseudotime.pdf"), width=mm2in(85), height=mm2in(85))
plot(NA, xlim=c(0, 1), ylim=c(-8, 8), type="n", main="Expression in Shin et al.", ylab="PC1 expression\n(smooth spline interpolation)", xlab="pseudotime")
for(group in levels(ros$Group)){
	smooth = smooth.spline(shin.pseudotime, ros.pc1[[group]], spar=0.7)
	lines(smooth, col=ros$group.colours[[group]], lwd=2)
}
dev.off()

writeLines(capture.output(sessionInfo()), "Figure_4/sessionInfo.txt")
##########
