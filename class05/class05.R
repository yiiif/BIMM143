#' ---
#' title: "Class 5: R Graphics"
#' author: "Yi Fu"
#' date: "Apr 16th, 2019"
#' ---

# Class 5 R graphics

# 2A. Line plot
baby.weight=read.table("bimm143_05_rstats/weight_chart.txt",header=T)
plot(baby.weight$Age,baby.weight$Weight,xlab="Age (months)",ylab="Weight (kg)",
     pch=15,cex=1.5,lwd=2,type="o",ylim=c(2,10),
     main="Baby Weight with Age")

# 2B. Boxplot
feature.count=read.table("bimm143_05_rstats/feature_counts.txt",header=T,sep="\t")

## below, left, above, right
par(mar=c(c(3.1, 11.1, 4.1, 2)))

barplot(feature.count$Count,names.arg=feature.count$Feature,horiz=T,
        las=1,xlim=c(0,80000),
        main="Number of features in the mouse GRCm38 genome")

# 2C. Histogram
x=c(rnorm(10000),rnorm(10000)+4)
hist(x,breaks=80)

# 3A. Providing Color Vector
gender.count=read.table("bimm143_05_rstats/male_female_counts.txt",header=T,sep="\t")
barplot(gender.count$Count,names.arg=gender.count$Sample,
        las=2,col=rainbow(nrow(gender.count)))

## Try different plots
barplot(gender.count$Count,names.arg=gender.count$Sample,
        las=2,col=c("red","blue"))

# 3B. Coloring by value
genes=read.table("bimm143_05_rstats/up_down_expression.txt",header=T,sep="\t")
table(genes$State)
plot(genes$Condition1, genes$Condition2, col=genes$State, 
     xlab="Expression condition 1", ylab="Expression condition 2")

palette(c("blue","gray","red"))
plot(genes$Condition1, genes$Condition2, col=genes$State,
     xlab="Expression condition 1", ylab="Expression condition 2")

# 3C. Dynamic use of color
# Lets plot expresion vs gene regulation
meth <- read.delim("bimm143_05_rstats/expression_methylation.txt")
plot(meth$gene.meth, meth$expression)

# Plot changing the plot character ('pch') to a solid circle
dcols <- densCols(meth$gene.meth, meth$expression)
plot(meth$gene.meth, meth$expression, col = dcols, pch = 20)

# Find the indices of genes with above 0 expresion and Plot just these genes
inds <- meth$expression > 0
plot(meth$gene.meth[inds], meth$expression[inds])

## Make a denisty color vector for these genes and plot
dcols <- densCols(meth$gene.meth[inds], meth$expression[inds])
plot(meth$gene.meth[inds], meth$expression[inds], col = dcols, pch = 20)

## Make a custom denisty color 
dcols.custom <- densCols(meth$gene.meth[inds], meth$expression[inds],
                         colramp = colorRampPalette(c("blue2",
                                                      "green2",
                                                      "red2",
                                                      "yellow")) )
plot(meth$gene.meth[inds], meth$expression[inds], col = dcols.custom, pch = 20)

## Plot the promoter.meth column against the gene.meth column.
plot(meth$promoter.meth, meth$gene.meth)
dcols.custom <- densCols(meth$promoter.meth[inds], meth$gene.meth[inds],
                         colramp = colorRampPalette(c("blue","red")))
plot(meth$promoter.meth, meth$gene.meth, 
     ylab="Gene Methylation", xlab="Promoter Methylation",
     col = dcols.custom)
