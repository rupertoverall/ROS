source('Extras.R')

print(load("Data/SVZ.RData")) # "svz"

# Figure 1 A
pca = PCA(t(svz$cpm))
var.explained = round((pca$sdev)^2 / sum(pca$sdev^2) * 100)

colour.factor = sapply(svz$Group, function(name) svz$group.colours[name] )

pdf(file="Figure_2/SVZ QCPrincipalComponentPlots.pdf", width=5, height=5)
plot(pca$PC1, pca$PC2, pch=20, cex=1.5, xlim=c(-130, 130), ylim=c(-130, 130), col=colour.factor, xlab=paste0("PC1 (", var.explained[1], " %)"), ylab=paste0("PC2 (", var.explained[2], " %)"))
legend("topright", names(svz$group.colours), col=svz$group.colours, ncol=1, pch=19, lty=0, text.col="black")
title(main="Principal component clustering of samples")
dev.off()

pdf(file="Figure_2/Venn.pdf", width=mm2in(85), height=mm2in(85))
# Enriched 'overlap' (no overlap by definition - differential expression with 2 groups makes them disjoint)
# Make overlap equal to the number of expressed genes that are not enriched in either
overlap = nrow(svz$cpm) - (length(svz$enriched$DG) + length(svz$enriched$SVZ))
VennDiagram::draw.pairwise.venn(
  area1 = length(svz$enriched$DG) + overlap, 
  area2 = length(svz$enriched$SVZ) + overlap, 
  cross.area = overlap, 
  category=names(svz$enriched), cex=0.7, cat.cex=0.7, cat.fontfamily=rep("sans", 2), fontfamily="sans", lty="blank", fill=svz$group.colours)
dev.off()

## Functional annotation (this takes a while to run)
source('Signatures.R')
write.results.to.excel(svz$signatures, file="Figure_2/SignatureEnrichment.xlsx", background=rownames(svz$counts))

# GO bar graphs
pdf(paste0("Figure_2/GO bar graphs.pdf"), width=12, height=5)
par(mar=c(4, 27, 3, 1), mfrow=c(1, 1))
. = sapply(c("DG", "SVZ"), function(go.group){
  go.set = as.data.frame(readxl::read_excel(paste0("Figure_2/SignatureEnrichment.xlsx"), sheet=paste0(go.group, "_GO")))
  bp = barplot(rev(-log10(go.set[1:10, "adjP"])), las=1, horiz=T, col=svz$group.colours[[go.group]], xlim=c(0, 30))
  text(par("usr")[1]-0.5, bp, labels=rev(go.set[1:10, "Term"]), srt=0, pos=2, offset=0, xpd=T)
  title(xlab=expression(paste(-log[10], "(", italic("p"), "-value)")), line=2.5)
  text(par("usr")[1], par("usr")[4]+0.5, labels=paste0("Top 10 GO terms for ", go.group), cex=1.2, font=2, srt=0, pos=3, offset=0, xpd=T)
})
dev.off()

writeLines(capture.output(sessionInfo()), "Figure_2/sessionInfo.txt")

##########
