source('Resources/Extras.R')

print(load("Data/ROS.RData")) # "ros"

## Heatmap of expression
# Get an extended list of genes associated with adult neurogenesis
ng.markers = data.frame(readxl::read_excel('Resources/Gene_lists.xlsx', sheet="NG_markers"), row.names = "Ensembl.ID", stringsAsFactors=F)
ng.markers = ng.markers[intersect(rownames(ng.markers), rownames(ros$cpm)), ]
plot.data = ros$cpm[rownames(ng.markers), order(ros$Group)]
colnames(plot.data) = ros$Group[order(ros$Group)]

png("Supplementary/Expression_heatmap.png", width = 85, height = (nrow(plot.data) * 3 + 15), units = "mm", res = 600)
plot.heatmap.raster(
	plot.data, 
	x.groups=colnames(plot.data), 
	x.labels=colnames(plot.data), 
	y.labels=ng.markers$Symbol,
	title="Expression of neurogenesis genes",
	title.cex=0.7,
	title.height=0,
	title.adjust=0,
	x.label.height=1,
	y.group.width = 0,
	y.label.width=2, 
	x.label.cex=0.7, 
	y.label.cex=0.4, 
	x.label.adjust=-2,
	y.label.adjust=-1.1,
	colour.palette=colour.schemes$bluetored.dark.3,
	absolute=F
)
dev.off()

png(file="Supplementary/Expression_heatmap_ColourScale.png", width=85, height=13, units="mm", res=600)
extreme.value = ceiling(max(abs(range(plot.data, na.rm = T))))
colour.scale.legend(range = seq(-extreme.value, extreme.value, by = 2), colour.schemes$bluetored.dark.3, cex.axis = 0.5)
dev.off()


writeLines(capture.output(sessionInfo()), "Figure_4/sessionInfo_Expression_heatmap.txt")
##########
