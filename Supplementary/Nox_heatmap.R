source('Resources/Extras.R')

print(load("Data/ROS.RData")) # "ros"

# ROS drops
# Re-run EdgeR comparing hi and lo to mid. 
centerGroup = relevel(ros$Group, "++") # Intercept is midROS
design = model.matrix(~centerGroup)
fit = edgeR::glmFit(edgeR::estimateDisp(edgeR::calcNormFactors(edgeR::DGEList(counts=ros$counts, group=centerGroup)), design), design)
drop.1 = edgeR::topTags(edgeR::glmLRT(fit, coef=2), n=nrow(ros$counts), p.value=1)$table
drop.2 = edgeR::topTags(edgeR::glmLRT(fit, coef=3), n=nrow(ros$counts), p.value=1)$table

nox.markers = c("Cybb", "Cyba", "Ncf1", "Ncf2", "Ncf4", "Noxo1", "Noxa1", "Rac1", "Rac2", "Enox1", "Enox2", "Nox1", "Nox3", "Nox4", "Duox1", "Duox2")
mapping = AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=nox.markers, columns=c("ENSEMBL"), "SYMBOL")
rownames(mapping) = mapping$SYMBOL # Also ensures uniqueness
mapping = mapping[nox.markers, ] # Ensure correct ordering
mapping = mapping[mapping$ENSEMBL %in% rownames(ros$counts), ] # Some of the genes have no read coverage in this dataset

# Filtered by significance
significance.data = data.frame("drop1" = drop.1[mapping$ENSEMBL, "FDR"] < 0.05, "drop2" = drop.2[mapping$ENSEMBL, "FDR"] < 0.05)
fc.data = data.frame("drop1" = -drop.1[mapping$ENSEMBL, "logFC"], "drop2" = drop.2[mapping$ENSEMBL, "logFC"]) # Signs of drop1 reversed due to level ordering
for(i in 1:ncol(fc.data)) fc.data[!significance.data[, i], i] = NA
rownames(fc.data) = mapping$ENSEMBL

# Balance plot data
extreme.value = ceiling(max(abs(range(fc.data, na.rm = T))))
plot.data = (fc.data / extreme.value + 1) / 2


png("Supplementary/Nox_heatmap.png", width = 85, height = (nrow(plot.data) * 3 + 25), units = "mm", res = 600)
plot.heatmap.raster(
	plot.data, 
	x.groups=colnames(plot.data), 
	x.labels=colnames(plot.data), 
	y.labels=mapping$SYMBOL,
	title="Changes in expression at each 'ROS drop'",
	title.cex=0.7,
	title.height=1,
	title.adjust=1,
	x.label.height=1, 
	y.label.width=2, 
	x.label.cex=0.7, 
	y.label.cex=0.4, 
	x.label.adjust=-1,
	y.label.adjust=-1.1,
	y.group.width = 0,
	colour.palette=colour.schemes$bluetored.dark.3,
	absolute=T
)
dev.off()

png(file="Supplementary/Nox_heatmap_ColourScale.png", width=85, height=18, units="mm", res=600)
colour.scale.legend(range = seq(-extreme.value, extreme.value, by = .5), colour.schemes$bluetored.dark.3)
dev.off()


writeLines(capture.output(sessionInfo()), "Supplementary/sessionInfo_nox_heatmap.txt")
##########



##########
