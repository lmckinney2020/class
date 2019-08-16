#' ---
#' title: "R basics - afternoon session"
#' author: "FAES-AABI"
#' date: "August 8, 2018"
#' output:
#'   word_document:
#'     toc: yes
#' ---
#' 
## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(root.dir = 'C://Users/michaloa/Documents')
knitr::opts_chunk$set(echo=TRUE,tidy=TRUE, dpi=300)
knitr::opts_chunk$set(warning=FALSE,message=FALSE)

#' 
#' August 8, 2018
#' 8:30 to 9:30 AM : R and RStudio - basic concepts (Presentation)
#' 9.30 to 12 Noon : R programming - variables, data structures, objects, functions, operators,
#' packages, debugging, scripts, bioconductor (Hands-on)
#' 1 to 4.30 PM: R programming - reading, writing, generating, manipulating, visualizing
#' data and markups (Hands-on)
#' 
#' ### 12 Packages
#' R packages (libraries) are a collection of R functions and sample data. They are stored under a directory called "library" in the R environment. Only a subset of available packages are installed by default, other libraries are added by the user as needed.
#' 
## ------------------------------------------------------------------------
# check where libraries are stored
.libPaths()

# the list of all installed packages
#library()

# load a library to your current working environment
library(survival)

# the list of all currently loaded packages
#search()


#' 
#' #### 12.1 R repositories
#' - CRAN https://cran.r-project.org/web/packages/available_packages_by_name.html
#' - Bioconductor https://www.bioconductor.org/packages/release/BiocViews.html#___Software
#' - Github https://github.com/
#' 
## ----fi.width=2,fig.height=1---------------------------------------------

## install a package from CRAN
#install.packages('heatmaply')

# check if loads correctly; it might be necessary to install additional required packages
library(heatmaply)

## install a package from Bioconductor

#. source installing function
source('https://bioconductor.org/biocLite.R')

#. install Bioconductor and default set of packages
#biocLite() 

#. install Bioconductor packages by name
#biocLite(c('DESeq2'))
library(DESeq2)

#. install Bioconductor workflow package
#biocLite(c('arrays'))
library(arrays)
#vignette('arrays')

## install packages from github

#.first get installer from CRAN
#install.packages('githubinstall')
library(githubinstall)
#. install packages
#githubinstall("heatmaply")

## detach package from current workspace
detach("package:heatmaply")

## call function from installed, but not loaded package
heatmaply::heatmaply(matrix(sample(15),3,5))

#' 
#' #### 12.1 R datasets
#' 
## ------------------------------------------------------------------------
# the list of sample datasets in current workspace
data()

# the list of datasets from a particular package
data(package = "ggplot2")
data(package = "heatmaply")

# the list of all available datasets from installed packages
data(package = .packages(all.available = TRUE))


#' 
#' ### 13 Data import
#' Data can be read into R from different sources and file formats.
#' 
#' #### 13.1 txt/csv/tsv
#' 
## ------------------------------------------------------------------------
# set directory where files are locally stored
folder = "C://Users/michaloa/Documents"

# set file name with extension (list.files could be helpful)
# list.files(folder, patt='.txt')
filename = "rData_survivalPackage_lungCancer.txt"

# read.table
path = file.path(folder,filename)
dat = read.table(path, sep='\t', header=TRUE, as.is=TRUE) #as.is=TRUE does not transform categorical variables to factors
dat_lung = dat

# read.delim (the same as read.table, but with sep='\t' and header=TRUE by default)
dat = read.delim(path, as.is=FALSE)

# find end of comment lines and skip lines from top
filename = "rData_survivalPackage_lungCancer_withComment.txt"
path = file.path(folder,filename)
readLines(path, n=5)
dat = read.delim(path, skip=2)

# set file name with extension (list.files could be helpful)
# list.files(folder, patt='.csv')
filename = "rData_survivalPackage_lungCancer.csv"

# read.csv
path = file.path(folder,filename)
dat = read.csv(path)


#' 
#' #### 13.1 xls/xlsx files
#' 
## ------------------------------------------------------------------------
# set directory where files are locally stored
folder = "C://Users/michaloa/Documents"
filename = "rData_survivalPackage_lungCancer.xlsx"
path = file.path(folder,filename)
library(XLConnect)
workbook = loadWorkbook(path)
sheets = getSheets(workbook)
sheets
dat_list = vector(mode='list', length=length(sheets))
names(dat_list) = make.names(sheets)
for(i in 1:length(sheets)){
  dat_list[[i]] = readWorksheet(workbook, sheet=sheets[i])
}
str(dat_list)
dat = do.call(cbind,dat_list)
dat = data.frame(dat)

#' 
#' #### 13.2 read URL
#' 
## ------------------------------------------------------------------------
library(data.table)
url = "http://assets.datacamp.com/blog_assets/chol.txt"
dat = as.data.frame(fread(url))
dim(dat)
names(dat)
dat[1:3,]
dat_chol = dat

#' 
#' ##### 13.3 compressed files
#' can read .gz, .bz2, .xz, or .zip files
#' 
## ------------------------------------------------------------------------
library(readr)
folder = "C://Users/michaloa/Documents"
filename = "rData_survivalPackage_lungCancer.zip"
path = file.path(folder,filename)
datt = read_delim(path, delim='\t')
dat = as.data.frame(datt)

#' 
#' ### 14 Data export
#' 
#' #### 14.1 tsv/csv/xlsx
#' 
## ------------------------------------------------------------------------
out = dat_chol
folder = "C://Users/michaloa/Documents"

# write to tab delimited file
filename = "dat_chol.txt"
path = file.path(folder,filename)
write.table(out, file=path, sep='\t', row=FALSE, quote=FALSE)

# write to comma separated file
filename = "dat_chol.csv"
path=file.path(folder,filename)
write.csv(out, file=path)

# write to excel file
filename = "dat_chol.xlsx"
path=file.path(folder,filename)

#. create new workbok
workbook <- loadWorkbook(path, create = TRUE)

#.create a worksheet
createSheet(workbook, name = "cholesterol")

#.write dataset to worksheet
writeWorksheet(workbook, out, sheet = "cholesterol", header=TRUE, startRow = 1, startCol = 1)

#.save workbook
saveWorkbook(workbook)

#' 
#' #### 14.2 save RData
#' 
#' The RData format (usually with extension .rdata or .rda) is a format designed for storing a complete R workspace or selected objects from a workspace in a form that can be loaded back by R.
#' 
## ------------------------------------------------------------------------
save(dat_chol, dat_lung, folder, dat_list, file=file.path(folder, 'example_data.RData'))
load(file.path(folder,'example_data.RData'))

#' 
#' 
#' ### 15 Data manipulation
#' 
#' ##### NeveRM, "A collection of breast cancer cell lines for the study of functionally distinct cancer subtypes"; [PMC2730521](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2730521/)
#' 
#' - 15.1 load data
#' 
## ------------------------------------------------------------------------
load(file.path(folder, "NeveG.RData"))
dat <- NeveG_data
str(dat)
summary(dat)

#' 
#' 
#' #### 15.2 reshape (wide/long formats)
#' 
## ------------------------------------------------------------------------
library(reshape2)

# wide format (each cell line data are spread across columns in one row)
dat_wide = NeveG_data

# melt to long (cell line data are stacked by measurement variables)
dat_long = melt(dat_wide, measure.vars=5:13, id.vars=1:4, variable.name='Gene', value.name='GEX')

# cast to wide format (matrix-like data.frame only with measurement variables)
dat_wide2 = acast(dat_long, Gene~cellLine, value.var = "GEX")

save(dat_wide, dat_long, dat_wide2, file=file.path(folder,'data_shapes.rda'))


#' 
#' 
#' #### 15.2 row/column statitstics
#' 
## ------------------------------------------------------------------------
dat = dat_wide2

# apply() function
means <- apply(dat, 1, mean, na.rm = T)
medians <- apply(dat, 1, median, na.rm = T)
stdevs <- apply(dat, 1, sd, na.rm = T)

# faster implementation of a series of commonly used descriptive statistics; rowMeans/colMeans, rowSds/colSds, rowMedians/colMedians, etc.

# apply() does not allow functions to be nested, if needed a function needs to be written and called as one

calc.SEM = function(x){
  n <- sum(!is.na(x))
  sd(x, na.rm = T)/sqrt(n)
}
sems <- apply(dat, 1, calc.SEM)
size <- apply(dat, 1, function(x) {sum(!is.na(x))})
cv <- stdevs/means*100

## COMBINE stats with data.frame() function
stats_wide <- data.frame(N = size, Mean = means, Median = medians, SD = stdevs, SEM = sems, CV = cv)
round(stats_wide,2)

## SAVE as *.RData object
save(stats_wide, file = 'Stats.RData')

#' 
#' #### 15.3 operations on lists
#' 
## ------------------------------------------------------------------------
# apply t.test on each row and output results to a list (default type of R t.test function)

dat = dat_wide2
group = dat_wide$ER
group1 = which(group == '[+]')
group0 = which(group == '[-]')
fun.test = function(x, g0, g1){
  return(t.test(x[g1], x[g0], var.equal=TRUE))
}
test_list = apply(dat,1,fun.test, g1=group1, g0=group0)
length(test_list)

# inspect single result
test = test_list[[1]]
names(test)

# retrieve p-values from each test with lapply
pvalue = lapply(test_list, function(x) x$p.value)
pvalue = unlist(pvalue)

# retrieve t-values and group means from each test with sapply
tvalue = sapply(test_list, function(x) x$statistic)
mvalue = sapply(test_list, function(x) x$estimate)


#' 
#' ### 16 Data visualizations
#' 
#' #### 16.1 Basic R graphics
#' 
## ---- fig.width=8, fig.height=8------------------------------------------
# Boxplots (wide format) 
boxplot(dat_wide[, 5:13], col = 'grey', main='Box plot', ylab = 'log2', boxwex = 0.3, las=2)

# Boxplot + dotplot with grouping (long format)
# combine box plot and  dot plot
par(oma = c(2,0,0,0))
boxplot(GEX ~ ER+Gene, data = dat_long,
        outline = FALSE, varwidth = TRUE,
        col = c('coral','lightblue'),
        main='Box plot + Dot plot',ylab = 'log2',
        las = 2)
stripchart(GEX ~ ER+Gene, data = dat_long,
           method = 'jitter', vertical = TRUE,
           pch=16, add=T)

# Scatterplot (select a pair of columns from wide format)

# set variables
x <- dat_wide$ABAT
y <- dat_wide$ESR1
limit <- range(c(x,y))
r <- round(cor(x,y),2)

# make plot
par(pty='s')
plot(x, y, ylim = limit, xlim = limit, xlab = 'ABAT', ylab = 'ESR1')
abline(0,1)
fit <- lm(y ~ x)
abline(fit, col='blue')
legend('bottomright', col = c('black','blue'), lty = 1, lwd = 2, legend = c('identity line', 'regression line'), cex=0.7)
mtext(paste('r', r, sep=' = '))
title('Scatter plot and correlation')

## Barplot with error bars
par(mfrow = c(1,2))

# plot means
d <- stats_wide
name <- rownames(stats_wide)
bar_position=barplot( d$Mean, names = name, las=2, density = 20, col='lightblue', border='darkblue', ylim=c(0,8))
arrows(bar_position, d$Mean, bar_position, d$Mean+d$SEM, mode=1, angle=90, length=0.02)
title(ylab = 'Geometric Mean', main = 'Bar plot')




#' 
#' #### 16.2 Grammar of Graphics
#' 
#' {ggplot2} is a high-level R package developed by Hadley Wickham based on Leland Wilkinson's concept of the grammar of graphics ("The Grammar of Graphics", 1999/2005).
#' 
#' In practical sense, with ggplot2 building blocks of a plot are independently specified and then combined to create complex multi-layered graphics.
#' 
#' Main building blocks  
#' 
#' - data
#' - aesthetic mapping 
#' - geometric object 
#' - statistical transformations
#' - coordinate system
#' - scales
#' - themes
#' - faceting
#' 
#' Basic plotting functions  
#' 
#' - qplot() - makes quick plots, resembles standard graphics
#' - ggplot() - the main plot function, allows detailed control of the graphical elements
#'  
#' ##### 16.2.1 qplot()  
#' 
## ----fig.width = 6, fig.height = 6---------------------------------------
# prepare histograms
f1A <- qplot(data = dat_wide, ageY , geom=c('histogram'))
f1B <- qplot(data = dat_wide, MEGF9, geom=c('histogram'))
# prepare plot with data points and smoothed line
f2 <- qplot(data = dat_wide, ageY, MEGF9, geom=c('point','smooth'))
# prepare boxplots
f3  <- qplot(data = dat_wide, ER, ageY, geom = c('boxplot'))
# arrange plots
library(gridExtra)
grid.arrange(f1A, f1B, f2, f3, ncol=2)

#' 
#' ##### 16.2.2 ggplot()
#' 
#' - mapping data to aesthetics (choosing axis, color/fill, shape/linetype, size)
## ------------------------------------------------------------------------
plot.data <- ggplot( data = dat_long, aes(x = geneCluster, y = GEX, fill = geneCluster) )
plot.blanc <- ggplot( data = dat_long, aes(x = geneCluster, y = GEX, color = geneCluster) )


#' 
#' - geoms (the actual marks plotted)
## ------------------------------------------------------------------------
violin <- geom_violin()
dot <- geom_jitter(shape = 16, position = position_jitter(0.2), size=1)
box <- geom_boxplot(outlier.size = -1, width = 0.25, fill = 'antiquewhite', color = 'antiquewhite', alpha = 0.5)

#' 
#' - stat 
## ------------------------------------------------------------------------
stat.median <- stat_summary(fun.data=median, geom="point", shape=23, size=2)
#stat.error <- stat_summary(fun.data = mean_err, geom="pointrange")

#' 
#' - theme and plot attributes
## ------------------------------------------------------------------------
theme <- theme_bw()
leg <- theme(legend.position = 'none')
titles <- labs(title = "Violin Plots", x = "", y = "log2 gene expression")

#' 
#' - scales
## ------------------------------------------------------------------------
paint.fill <- scale_fill_manual(values = c("red","orange", "blue"))
paint.col <- scale_color_manual(values = c("red","orange", "blue"))
y.axis <- scale_y_continuous(breaks = seq(2,14,2))

#' 
#' - combine layers
## ------------------------------------------------------------------------
plot.data + violin + dot + box + titles + leg + y.axis + paint.fill + coord_flip()
plot.blanc + geom_violin(trim = FALSE)  + titles + theme + leg + paint.col

#' 
#' 
#' - Faceting: define subsets as the levels of a single grouping variable with facet_wrap()
## ------------------------------------------------------------------------
dat.plot <- ggplot(data = dat_long, aes(x = GEX, fill = ER))
histo <- geom_histogram(binwidth=1, position = 'identity', alpha = 0.6)
panel.wrap <- facet_wrap(~Gene)

dat.plot + histo + panel.wrap + labs(x = 'log2 Expression') + theme_light()

#' 
#' - Faceting: define subsets as the crossing of two grouping variables with facet.grid()
## ------------------------------------------------------------------------
dat.plot <- ggplot(data = dat_long, aes(x = GEX, fill = ER))
histo <- geom_histogram(binwidth=1, position = 'identity', alpha = 0.6)
panel.grid <- facet_grid(ER ~ geneCluster)
theme.bg <-  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.plot + histo + panel.grid + theme.bg


#' 
#' - Good ggplot2 summary from R Studio
#' [ggplot2-cheatsheet](https://www.rstudio.com/wp-content/uploads/2015/03/ggplot2-cheatsheet.pdf)
#' 
#' 
#' ### 17 RMarkdown
#' *RMarkdown* allows to create documents that serve as a neat and reproducible documentation of an analysis. *RMarkdown* is a very simple 'markup' language which provides methods for creating documents with headers, images, links etc. from plain text files, while keeping the code alongside its output (graphs, tables, etc.) with conventional text to explain it, a bit like a notebook. Markdown documents can be converted to other file types like .html , .pdf, or .doc
#' 
#' **To create a new RMarkdown file (.Rmd), select File -> New File -> R Markdown in RStudio, then choose the file type (HTML here)**
#' 
#' - Table
## ------------------------------------------------------------------------
load("C://Users/michaloa/Documents/example_data.RData")
tab = dat_chol

library(DT)
datatable(tab, colnames=names(tab), options = list(scrollX = FALSE, keys = TRUE, pageLength = 10),  caption=paste("Table 1"), rownames = TRUE)


#' 
#' - Figure
## ----fig.width=2, fig.height=2, dpi=300----------------------------------
dat = dat_wide2
cormat = cor(dat)
library(heatmaply)
heatmaply(cormat, cexRow = 0.5, cexCol=0.5)

#' 
#' #### 17.1 Code chunk instructions
#' 
#' #### 17.2 Formatting text
#' 
#' *Italic*
#' 
#' **Bold**
#' 
#' * Unordered list item
#' 
#' 1. Ordered list item
#' 
#' [Link to RStudio RMarkdown](https://rmarkdown.rstudio.com/index.html)
#' 
#' 
