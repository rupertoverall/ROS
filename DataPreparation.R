# This script creates the canonical data files (used for analysis and GEO submission)
# The RData archives contain output from the 'calculate.signatures' workflow.
source('http://research.rupertoverall.net/Signatures.R') # All diferential analysis now done by 'calculate.signatures'

##
# SVZ

gse = GEOquery::getGEO("GSE124095")
raw.counts = read.tab('~/genreg_ruperto/Sequencing_Data/SVZ/BAM/GeneCounts.tab', row.names=1)
n.counts.raw = dim(raw.counts)
raw.headers = read.tab('../Additional Data/SVZ_headers.tab', stringsAsFactors=F, row.names=1)
raw.headers = raw.headers[colnames(raw.counts), ]
headers = subset(raw.headers, Environment=="STD" & Sort=="POS")
raw.counts = raw.counts[, rownames(headers)]
counts = raw.counts[cpm.presence.filter(raw.counts, headers$Region, threshold=1), ] # _One CPM_ as threshold
n.counts.present = dim(counts)
Group  = factor(headers$Region, levels=c("DG", "SVZ"))
svz = calculate.signatures(counts, Group, filter=F) # Counts has already been thresholded
pca = PCA(t(svz$cpm))
outlier = colnames(svz$cpm)[which(pca$PC2 > 100)]
not.outlier = colnames(svz$cpm)[which(pca$PC2 < 100)]
headers = headers[not.outlier, ]
counts = counts[, not.outlier]
n.counts.cleaned = dim(counts)
# Re-run base analysis
Group  = factor(headers$Region, levels=c("DG", "SVZ"))
group.colours = c("DG"="#f03900", "SVZ"="#99cc00")
svz = calculate.signatures(counts, Group, filter=F) # Counts has already been thresholded
sessionInfo = sessionInfo()
report = list(n.counts.raw=n.counts.raw, n.counts.present=n.counts.present, n.counts.cleaned=n.counts.cleaned, sessionInfo=sessionInfo)
svz = c(svz, list(Group=Group, group.colours=group.colours, headers=headers, counts=counts, report=report))
save(svz, file="SVZ.RData")

##
# ROS
raw.counts = read.tab('~/genreg_ruperto/Sequencing_Data/ROS/BAM/GeneCounts.tab', row.names=1)
n.counts.raw = dim(raw.counts)
raw.headers = data.frame(readxl::read_excel('../Additional Data/raw_data_Nes-GFP_ROSsorted_sequencing.xls'), stringsAsFactors=F, row.names=1)
raw.headers = raw.headers[grep("L", rownames(raw.headers)), ] # One sample not sequenced
headers = raw.headers[colnames(raw.counts), ]
counts = raw.counts[cpm.presence.filter(raw.counts, headers$ROSGroup, threshold=1), ] # _One CPM_ as threshold
n.counts.present = dim(counts)
Group  = factor(sapply(headers$ROSGroup, function(x) paste0(rep("+", x), collapse="") ), levels=c("+++", "++", "+"))
ros = calculate.signatures(counts, Group, filter=F) # Counts has already been thresholded
pca = PCA(t(ros$cpm))
# One outlier samples is present - which we shall remove.
outlier = colnames(ros$cpm)[which(pca$PC2 < -100)]
not.outlier = colnames(ros$cpm)[which(pca$PC2 > -100)]
headers = headers[not.outlier, ]
counts = counts[, not.outlier]
n.counts.cleaned = dim(counts)
# Re-run base analysis
Group  = factor(sapply(headers$ROSGroup, function(x) paste0(rep("+", x), collapse="") ), levels=c("+++", "++", "+"))
group.colours = c("+++"="#a70b02", "++"="#df7b31", "+"="#f1af22")
ros = calculate.signatures(counts, Group, filter=F) # Counts has already been thresholded
sessionInfo = sessionInfo()
report = list(n.counts.raw=n.counts.raw, n.counts.present=n.counts.present, n.counts.cleaned=n.counts.cleaned, sessionInfo=sessionInfo)
ros = c(ros, list(Group=Group, group.colours=group.colours, headers=headers, counts=counts, report=report))
save(ros, file="ROS.RData")

