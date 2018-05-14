# source("https://bioconductor.org/biocLite.R")
# biocLite("iSEE")
library("iSEE")
library("scRNAseq")
data("allen")
library("scater")
sce <- as(allen, "SingleCellExperiment")
counts(sce) <- assay(sce, "tophat_counts")
sce <- normalize(sce)
sce <- runPCA(sce)
sce <- runTSNE(sce)
# iSEE(sce)
# couple of genes to check: (Zeisel, Science 2015; Tasic, Nature Neuroscience 2016)
# Tbr1 (TF required for the final differentiation of cortical projection neurons);
# Snap25 (pan-neuronal); 
# Rorb (mostly L4 and L5a); 
# Foxp2 (L6)


## Got these from... "Display panel settings"
################################################################################
# Settings for reduced dimension plots
################################################################################

redDimPlotArgs <- new('DataFrame', nrows=5L, rownames=paste0('redDimPlot', seq_len(5)))
redDimPlotArgs[['Type']] <- c(1L, 1L, 1L, 1L, 1L)
redDimPlotArgs[['XAxis']] <- c(1L, 1L, 1L, 1L, 1L)
redDimPlotArgs[['YAxis']] <- c(2L, 2L, 2L, 2L, 2L)
redDimPlotArgs[['DataBoxOpen']] <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
redDimPlotArgs[['VisualBoxOpen']] <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
redDimPlotArgs[['SelectBoxOpen']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
redDimPlotArgs[['SelectByPlot']] <- c("---", "---", "---", "---", "---")
redDimPlotArgs[['SelectEffect']] <- c("Transparent", "Transparent", "Transparent", "Transparent", "Transparent")
redDimPlotArgs[['SelectAlpha']] <- c(0.1, 0.1, 0.1, 0.1, 0.1)
redDimPlotArgs[['SelectColor']] <- c("#FF0000", "red", "red", "red", "red")

tmp <- vector('list', 5)
redDimPlotArgs[['BrushData']] <- tmp

tmp <- vector('list', 5)
tmp[[1]] <- "Color"
tmp[[2]] <- "Color"
tmp[[3]] <- "Color"
tmp[[4]] <- "Color"
tmp[[5]] <- "Color"
redDimPlotArgs[['VisualChoices']] <- tmp
redDimPlotArgs[['PointSize']] <- c(1, 1, 1, 1, 1)
redDimPlotArgs[['PointAlpha']] <- c(1, 1, 1, 1, 1)
redDimPlotArgs[['Downsample']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
redDimPlotArgs[['SampleRes']] <- c(200, 200, 200, 200, 200)
redDimPlotArgs[['FontSize']] <- c(1, 1, 1, 1, 1)
redDimPlotArgs[['LegendPosition']] <- c("Bottom", "Bottom", "Bottom", "Bottom", "Bottom")

tmp <- vector('list', 5)
redDimPlotArgs[['ZoomData']] <- tmp

tmp <- vector('list', 5)
tmp[[1]] <- structure(c(-1.89613850834849, -4.23568865593679, 4.42064689013992, 8.74881466317827, 
                        10.5034772738695, 10.3864997664901, -1.89613850834849, 15.691945504741, 9.20530407803766, 
                        2.05525614178515, 1.31813779784159, 4.34032300801018, 8.76303307167152, 15.691945504741
), .Dim = c(7L, 2L), closed = TRUE, flipped = FALSE)
redDimPlotArgs[['LassoData']] <- tmp
redDimPlotArgs[['ColorBy']] <- c("Feature name", "None", "None", "None", "None")
redDimPlotArgs[['ColorByDefaultColor']] <- c("#000000", "black", "black", "black", "black")
redDimPlotArgs[['ColorByColData']] <- c("NREADS", "NREADS", "NREADS", "NREADS", "NREADS")
redDimPlotArgs[['ColorByRowTable']] <- c("Row statistics table 1", "---", "---", "---", "---")
redDimPlotArgs[['ColorByFeatName']] <- c(6307L, 1L, 1L, 1L, 1L)
redDimPlotArgs[['ColorByFeatNameAssay']] <- c(6L, 6L, 6L, 6L, 6L)

################################################################################
# Settings for column data plots
################################################################################

colDataPlotArgs <- new('DataFrame', nrows=5L, rownames=paste0('colDataPlot', seq_len(5)))
colDataPlotArgs[['YAxis']] <- c("NREADS", "NREADS", "NREADS", "NREADS", "NREADS")
colDataPlotArgs[['XAxis']] <- c("None", "None", "None", "None", "None")
colDataPlotArgs[['XAxisColData']] <- c("NALIGNED", "NALIGNED", "NALIGNED", "NALIGNED", "NALIGNED")
colDataPlotArgs[['DataBoxOpen']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
colDataPlotArgs[['VisualBoxOpen']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
colDataPlotArgs[['SelectBoxOpen']] <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
colDataPlotArgs[['SelectByPlot']] <- c("Feature assay plot 3", "---", "---", "---", "---")
colDataPlotArgs[['SelectEffect']] <- c("Transparent", "Transparent", "Transparent", "Transparent", "Transparent")
colDataPlotArgs[['SelectAlpha']] <- c(0.1, 0.1, 0.1, 0.1, 0.1)
colDataPlotArgs[['SelectColor']] <- c("#FF0000", "red", "red", "red", "red")

tmp <- vector('list', 5)
colDataPlotArgs[['BrushData']] <- tmp

tmp <- vector('list', 5)
tmp[[1]] <- "Color"
tmp[[2]] <- "Color"
tmp[[3]] <- "Color"
tmp[[4]] <- "Color"
tmp[[5]] <- "Color"
colDataPlotArgs[['VisualChoices']] <- tmp
colDataPlotArgs[['PointSize']] <- c(1, 1, 1, 1, 1)
colDataPlotArgs[['PointAlpha']] <- c(1, 1, 1, 1, 1)
colDataPlotArgs[['Downsample']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
colDataPlotArgs[['SampleRes']] <- c(200, 200, 200, 200, 200)
colDataPlotArgs[['FontSize']] <- c(1, 1, 1, 1, 1)
colDataPlotArgs[['LegendPosition']] <- c("Bottom", "Bottom", "Bottom", "Bottom", "Bottom")

tmp <- vector('list', 5)
colDataPlotArgs[['ZoomData']] <- tmp

tmp <- vector('list', 5)
colDataPlotArgs[['LassoData']] <- tmp
colDataPlotArgs[['ColorBy']] <- c("None", "None", "None", "None", "None")
colDataPlotArgs[['ColorByDefaultColor']] <- c("#000000", "black", "black", "black", "black")
colDataPlotArgs[['ColorByColData']] <- c("NREADS", "NREADS", "NREADS", "NREADS", "NREADS")
colDataPlotArgs[['ColorByRowTable']] <- c("---", "---", "---", "---", "---")
colDataPlotArgs[['ColorByFeatName']] <- c(1L, 1L, 1L, 1L, 1L)
colDataPlotArgs[['ColorByFeatNameAssay']] <- c(6L, 6L, 6L, 6L, 6L)

################################################################################
# Settings for feature assay plots
################################################################################

featAssayPlotArgs <- new('DataFrame', nrows=5L, rownames=paste0('featAssayPlot', seq_len(5)))
featAssayPlotArgs[['Assay']] <- c(6L, 6L, 6L, 6L, 6L)
featAssayPlotArgs[['XAxis']] <- c("None", "None", "Feature name", "None", "None")
featAssayPlotArgs[['XAxisColData']] <- c("NREADS", "NREADS", "NREADS", "NREADS", "NREADS")
featAssayPlotArgs[['XAxisFeatName']] <- c(1L, 1L, 15644L, 1L, 1L)
featAssayPlotArgs[['XAxisRowTable']] <- c("---", "---", "---", "---", "---")
featAssayPlotArgs[['YAxisFeatName']] <- c(15644L, 6307L, 6307L, 1L, 1L)
featAssayPlotArgs[['YAxisRowTable']] <- c("---", "---", "---", "---", "---")
featAssayPlotArgs[['DataBoxOpen']] <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
featAssayPlotArgs[['VisualBoxOpen']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
featAssayPlotArgs[['SelectBoxOpen']] <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
featAssayPlotArgs[['SelectByPlot']] <- c("Reduced dimension plot 1", "Reduced dimension plot 1", "Reduced dimension plot 1", 
                                         "---", "---")
featAssayPlotArgs[['SelectEffect']] <- c("Transparent", "Transparent", "Color", "Transparent", "Transparent")
featAssayPlotArgs[['SelectAlpha']] <- c(0.1, 0.1, 0.1, 0.1, 0.1)
featAssayPlotArgs[['SelectColor']] <- c("#FF0000", "#FF0000", "#CC00FF", "red", "red")

tmp <- vector('list', 5)
tmp[[3]] <- list(xmin = -0.28574301818977, xmax = 4.104681314094, ymin = 7.1341867267412, ymax = 12.60569466396, 
                 mapping = list(x = "X", y = "Y"), domain = list(left = -0.713907040920163, right = 14.9920478593234, 
                                                                 bottom = -0.61977934161258, top = 13.0153661738642), range = list(left = 38.1286922089041, 
                                                                                                                                   right = 449.520547945205, bottom = 465.836151541096, top = 24.7473723724842), 
                 log = list(x = NULL, y = NULL), direction = "xy", brushId = "featAssayPlot3_Brush", 
                 outputId = "featAssayPlot3")
featAssayPlotArgs[['BrushData']] <- tmp

tmp <- vector('list', 5)
tmp[[1]] <- "Color"
tmp[[2]] <- "Color"
tmp[[3]] <- "Color"
tmp[[4]] <- "Color"
tmp[[5]] <- "Color"
featAssayPlotArgs[['VisualChoices']] <- tmp
featAssayPlotArgs[['PointSize']] <- c(1, 1, 1, 1, 1)
featAssayPlotArgs[['PointAlpha']] <- c(1, 1, 1, 1, 1)
featAssayPlotArgs[['Downsample']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
featAssayPlotArgs[['SampleRes']] <- c(200, 200, 200, 200, 200)
featAssayPlotArgs[['FontSize']] <- c(1, 1, 1, 1, 1)
featAssayPlotArgs[['LegendPosition']] <- c("Bottom", "Bottom", "Bottom", "Bottom", "Bottom")

tmp <- vector('list', 5)
featAssayPlotArgs[['ZoomData']] <- tmp

tmp <- vector('list', 5)
featAssayPlotArgs[['LassoData']] <- tmp
featAssayPlotArgs[['ColorBy']] <- c("None", "None", "None", "None", "None")
featAssayPlotArgs[['ColorByDefaultColor']] <- c("#000000", "#000000", "#000000", "black", "black")
featAssayPlotArgs[['ColorByColData']] <- c("NREADS", "NREADS", "NREADS", "NREADS", "NREADS")
featAssayPlotArgs[['ColorByRowTable']] <- c("---", "---", "---", "---", "---")
featAssayPlotArgs[['ColorByFeatName']] <- c(1L, 1L, 1L, 1L, 1L)
featAssayPlotArgs[['ColorByFeatNameAssay']] <- c(6L, 6L, 6L, 6L, 6L)

################################################################################
# Settings for row statistics tables
################################################################################

rowStatTableArgs <- new('DataFrame', nrows=5L, rownames=paste0('rowStatTable', seq_len(5)))
rowStatTableArgs[['Selected']] <- c(6307L, 1L, 1L, 1L, 1L)
rowStatTableArgs[['Search']] <- c("Foxp2", "", "", "", "")

tmp <- vector('list', 5)
tmp[[1]] <- ""
tmp[[2]] <- character(0)
tmp[[3]] <- character(0)
tmp[[4]] <- character(0)
tmp[[5]] <- character(0)
rowStatTableArgs[['SearchColumns']] <- tmp
rowStatTableArgs[['SelectBoxOpen']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
rowStatTableArgs[['SelectByPlot']] <- c("---", "---", "---", "---", "---")

################################################################################
# Settings for row data plots
################################################################################

rowDataPlotArgs <- new('DataFrame', nrows=0L, rownames=paste0('rowDataPlot', seq_len(0)))
rowDataPlotArgs[['YAxis']] <- character(0)
rowDataPlotArgs[['XAxis']] <- character(0)
rowDataPlotArgs[['XAxisRowData']] <- character(0)
rowDataPlotArgs[['DataBoxOpen']] <- logical(0)
rowDataPlotArgs[['VisualBoxOpen']] <- logical(0)
rowDataPlotArgs[['SelectBoxOpen']] <- logical(0)
rowDataPlotArgs[['SelectByPlot']] <- character(0)
rowDataPlotArgs[['SelectEffect']] <- character(0)
rowDataPlotArgs[['SelectAlpha']] <- numeric(0)
rowDataPlotArgs[['SelectColor']] <- character(0)

tmp <- vector('list', 0)
rowDataPlotArgs[['BrushData']] <- tmp

tmp <- vector('list', 0)
rowDataPlotArgs[['VisualChoices']] <- tmp
rowDataPlotArgs[['PointSize']] <- numeric(0)
rowDataPlotArgs[['PointAlpha']] <- numeric(0)
rowDataPlotArgs[['Downsample']] <- logical(0)
rowDataPlotArgs[['SampleRes']] <- numeric(0)
rowDataPlotArgs[['FontSize']] <- numeric(0)
rowDataPlotArgs[['LegendPosition']] <- character(0)

tmp <- vector('list', 0)
rowDataPlotArgs[['ZoomData']] <- tmp

tmp <- vector('list', 0)
rowDataPlotArgs[['LassoData']] <- tmp
rowDataPlotArgs[['ColorBy']] <- character(0)
rowDataPlotArgs[['ColorByDefaultColor']] <- character(0)
rowDataPlotArgs[['ColorByRowData']] <- character(0)
rowDataPlotArgs[['ColorByRowTable']] <- character(0)
rowDataPlotArgs[['ColorByFeatName']] <- integer(0)
rowDataPlotArgs[['ColorByFeatNameColor']] <- character(0)

################################################################################
# Settings for custom column plots
################################################################################

customColPlotArgs <- new('DataFrame', nrows=0L, rownames=paste0('customColPlot', seq_len(0)))
customColPlotArgs[['Function']] <- character(0)
customColPlotArgs[['DataBoxOpen']] <- logical(0)
customColPlotArgs[['VisualBoxOpen']] <- logical(0)
customColPlotArgs[['SelectBoxOpen']] <- logical(0)
customColPlotArgs[['SelectByPlot']] <- character(0)
customColPlotArgs[['SelectEffect']] <- character(0)
customColPlotArgs[['SelectAlpha']] <- numeric(0)
customColPlotArgs[['SelectColor']] <- character(0)

tmp <- vector('list', 0)
customColPlotArgs[['BrushData']] <- tmp

tmp <- vector('list', 0)
customColPlotArgs[['VisualChoices']] <- tmp
customColPlotArgs[['PointSize']] <- numeric(0)
customColPlotArgs[['PointAlpha']] <- numeric(0)
customColPlotArgs[['Downsample']] <- logical(0)
customColPlotArgs[['SampleRes']] <- numeric(0)
customColPlotArgs[['FontSize']] <- numeric(0)
customColPlotArgs[['LegendPosition']] <- character(0)

tmp <- vector('list', 0)
customColPlotArgs[['ZoomData']] <- tmp

tmp <- vector('list', 0)
customColPlotArgs[['LassoData']] <- tmp
customColPlotArgs[['ColorBy']] <- character(0)
customColPlotArgs[['ColorByDefaultColor']] <- character(0)
customColPlotArgs[['ColorByColData']] <- character(0)
customColPlotArgs[['ColorByRowTable']] <- character(0)
customColPlotArgs[['ColorByFeatName']] <- integer(0)
customColPlotArgs[['ColorByFeatNameAssay']] <- integer(0)

################################################################################
# Settings for heat maps
################################################################################

heatMapPlotArgs <- new('DataFrame', nrows=5L, rownames=paste0('heatMapPlot', seq_len(5)))
heatMapPlotArgs[['Assay']] <- c(6L, 6L, 6L, 6L, 6L)
heatMapPlotArgs[['FeatNameBoxOpen']] <- c(TRUE, FALSE, FALSE, FALSE, FALSE)

tmp <- vector('list', 5)
tmp[[1]] <- c(15644L, 6307L, 17992L, 6306L, 6309L)
tmp[[2]] <- 1L
tmp[[3]] <- 1L
tmp[[4]] <- 1L
tmp[[5]] <- 1L
heatMapPlotArgs[['FeatName']] <- tmp
heatMapPlotArgs[['ColDataBoxOpen']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)

tmp <- vector('list', 5)
tmp[[1]] <- "NREADS"
tmp[[2]] <- "NREADS"
tmp[[3]] <- "NREADS"
tmp[[4]] <- "NREADS"
tmp[[5]] <- "NREADS"
heatMapPlotArgs[['ColData']] <- tmp
heatMapPlotArgs[['FeatNameSource']] <- c("Row statistics table 1", "---", "---", "---", "---")

tmp <- vector('list', 5)
tmp[[1]] <- "Centered"
tmp[[2]] <- "Centered"
tmp[[3]] <- "Centered"
tmp[[4]] <- "Centered"
tmp[[5]] <- "Centered"
heatMapPlotArgs[['CenterScale']] <- tmp
heatMapPlotArgs[['Lower']] <- c(NA, -Inf, -Inf, -Inf, -Inf)
heatMapPlotArgs[['Upper']] <- c(NA, Inf, Inf, Inf, Inf)
heatMapPlotArgs[['ColorScale']] <- c("blue-white-orange", "purple-black-yellow", "purple-black-yellow", "purple-black-yellow", 
                                     "purple-black-yellow")

tmp <- vector('list', 5)
heatMapPlotArgs[['ZoomData']] <- tmp
heatMapPlotArgs[['SelectBoxOpen']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
heatMapPlotArgs[['SelectByPlot']] <- c("---", "---", "---", "---", "---")
heatMapPlotArgs[['SelectEffect']] <- c("Transparent", "Transparent", "Transparent", "Transparent", "Transparent")
heatMapPlotArgs[['SelectAlpha']] <- c(0.1, 0.1, 0.1, 0.1, 0.1)
heatMapPlotArgs[['SelectColor']] <- c("red", "red", "red", "red", "red")


################################################################################
# Initial panel settings
################################################################################

initialPanels <- DataFrame(
  Name=c("Reduced dimension plot 1", "Feature assay plot 1", "Row statistics table 1", "Heat map 1", 
         "Feature assay plot 2", "Feature assay plot 3", "Column data plot 1"),
  Width=c(4L, 4L, 4L, 4L, 4L, 4L, 4L),
  Height=c(500L, 500L, 500L, 500L, 500L, 500L, 500L)
)




### Launching!


iSEE(sce,
     redDimArgs = redDimPlotArgs, 
     colDataArgs = colDataPlotArgs,
     featAssayArgs = featAssayPlotArgs,
     rowStatArgs = rowStatTableArgs,
     rowDataArgs = NULL, 
     heatMapArgs = heatMapPlotArgs,
     initialPanels = initialPanels,
     appTitle = "eRum2018 iSEE demo!")





## Got these from... "Extract the R code"

## The following list of commands will generate the plots created using iSEE.
## Copy them into a script or an R session containing your SingleCellExperiment.
## All commands below refer to your SingleCellExperiment object as `se`.

se <- sce
colData(se)[,"sizeFactors(se)"] <- sizeFactors(se)
colormap <- ExperimentColorMap()
colormap <- synchronizeAssays(colormap, se)
all_coordinates <- list()
custom_col_fun <- NULL

################################################################################
# Defining brushes
################################################################################

all_brushes <- list()
all_brushes[['featAssayPlot3']] <- list(xmin = -0.28574301818977, xmax = 4.104681314094, ymin = 7.1341867267412, ymax = 12.60569466396, 
                                        mapping = list(x = "X", y = "Y"), domain = list(left = -0.713907040920163, right = 14.9920478593234, 
                                                                                        bottom = -0.61977934161258, top = 13.0153661738642), range = list(left = 38.1286922089041, 
                                                                                                                                                          right = 449.520547945205, bottom = 465.836151541096, top = 24.7473723724842), 
                                        log = list(x = NULL, y = NULL), direction = "xy", brushId = "featAssayPlot3_Brush", 
                                        outputId = "featAssayPlot3")

################################################################################
# Defining lassos
################################################################################

all_lassos <- list()
all_lassos[['redDimPlot1']] <- structure(c(-1.89613850834849, -4.23568865593679, 4.42064689013992, 8.74881466317827, 
                                           10.5034772738695, 10.3864997664901, -1.89613850834849, 15.691945504741, 9.20530407803766, 
                                           2.05525614178515, 1.31813779784159, 4.34032300801018, 8.76303307167152, 15.691945504741
), .Dim = c(7L, 2L), closed = TRUE, flipped = FALSE)

################################################################################
## Reduced dimension plot 1
################################################################################

red.dim <- reducedDim(se, 1);
plot.data <- data.frame(X = red.dim[, 1], Y = red.dim[, 2], row.names=colnames(se));
plot.data$ColorBy <- assay(se, 6)[6307,];
plot.data <- subset(plot.data, !is.na(X) & !is.na(Y));

# Saving data for transmission
all_coordinates[['redDimPlot1']] <- plot.data

# Creating the plot
ggplot() +
  geom_point(aes(x = X, y = Y, color = ColorBy), alpha = 1, plot.data, size=1) +
  labs(x = "Dimension 1", y = "Dimension 2", color = "Foxp2\n(logcounts)", title = "(1) PCA") +
  coord_cartesian(xlim = range(plot.data$X, na.rm = TRUE),
                  ylim = range(plot.data$Y, na.rm = TRUE), expand = TRUE) +
  scale_color_gradientn(colors=assayColorMap(colormap, 6L, discrete=FALSE)(21), na.value='grey50') +
  scale_fill_gradientn(colors=assayColorMap(colormap, 6L, discrete=FALSE)(21), na.value='grey50') +
  theme_bw() +
  theme(legend.position = 'bottom', legend.text=element_text(size=9), legend.title=element_text(size=11),
        axis.text=element_text(size=10), axis.title=element_text(size=12), title=element_text(size=12)) +
  geom_polygon(aes(x = x, y = y), alpha=0.25, color='#3c8dbc', 
               data=data.frame(x = all_lassos[['redDimPlot1']][,1], y = all_lassos[['redDimPlot1']][,2]), 
               inherit.aes=FALSE, fill = '#D8E8F1') +
  scale_fill_manual(values = c('TRUE' = '#3c8dbc', 'FALSE' = '#D8E8F1')) +
  guides(shape = 'none')

################################################################################
## Feature assay plot 1
################################################################################

plot.data <- data.frame(Y=assay(se, 6)[15644L,], row.names = colnames(se))
plot.data$X <- factor(character(ncol(se)))
plot.data <- subset(plot.data, !is.na(X) & !is.na(Y));

# Receiving point selection
selected_pts <- mgcv::in.out(all_lassos[['redDimPlot1']], cbind(as.numeric(all_coordinates[['redDimPlot1']]$X), as.numeric(all_coordinates[['redDimPlot1']]$Y)))
plot.data$SelectBy <- rownames(plot.data) %in% rownames(all_coordinates[['redDimPlot1']])[selected_pts]

# Saving data for transmission
all_coordinates[['featAssayPlot1']] <- plot.data

# Setting up plot coordinates
plot.data$GroupBy <- plot.data$X;
set.seed(100);
plot.data$jitteredX <- vipor::offsetX(plot.data$Y,
                                      x=plot.data$X, width=0.4, varwidth=FALSE, adjust=1,
                                      method='quasirandom', nbins=NULL) + as.integer(plot.data$X);

# Creating the plot
ggplot() +
  geom_violin(aes(x = X, y = Y, group = GroupBy), alpha = 0.2, data=plot.data, scale = 'width', width = 0.8) +
  geom_point(aes(y = Y, x = jitteredX), subset(plot.data, !SelectBy), alpha = 0.10, color='#000000', size=1) +
  geom_point(aes(y = Y, x = jitteredX), subset(plot.data, SelectBy), color='#000000', size=1) +
  labs(x = "", y = "Rorb (logcounts)", title = "Rorb ") +
  coord_cartesian(ylim = range(plot.data$Y, na.rm=TRUE), expand = TRUE) +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  theme(legend.position = 'bottom', legend.text=element_text(size=9), 
        legend.title=element_text(size=11), legend.box = 'vertical',
        axis.text.x = element_text(angle=90, size=10, hjust=1, vjust=0.5), 
        axis.text.y=element_text(size=10), 
        axis.title=element_text(size=12), title=element_text(size=12))

################################################################################
## Feature assay plot 2
################################################################################

plot.data <- data.frame(Y=assay(se, 6)[6307L,], row.names = colnames(se))
plot.data$X <- factor(character(ncol(se)))
plot.data <- subset(plot.data, !is.na(X) & !is.na(Y));

# Receiving point selection
selected_pts <- mgcv::in.out(all_lassos[['redDimPlot1']], cbind(as.numeric(all_coordinates[['redDimPlot1']]$X), as.numeric(all_coordinates[['redDimPlot1']]$Y)))
plot.data$SelectBy <- rownames(plot.data) %in% rownames(all_coordinates[['redDimPlot1']])[selected_pts]

# Saving data for transmission
all_coordinates[['featAssayPlot2']] <- plot.data

# Setting up plot coordinates
plot.data$GroupBy <- plot.data$X;
set.seed(100);
plot.data$jitteredX <- vipor::offsetX(plot.data$Y,
                                      x=plot.data$X, width=0.4, varwidth=FALSE, adjust=1,
                                      method='quasirandom', nbins=NULL) + as.integer(plot.data$X);

# Creating the plot
ggplot() +
  geom_violin(aes(x = X, y = Y, group = GroupBy), alpha = 0.2, data=plot.data, scale = 'width', width = 0.8) +
  geom_point(aes(y = Y, x = jitteredX), subset(plot.data, !SelectBy), alpha = 0.10, color='#000000', size=1) +
  geom_point(aes(y = Y, x = jitteredX), subset(plot.data, SelectBy), color='#000000', size=1) +
  labs(x = "", y = "Foxp2 (logcounts)", title = "Foxp2 ") +
  coord_cartesian(ylim = range(plot.data$Y, na.rm=TRUE), expand = TRUE) +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  theme(legend.position = 'bottom', legend.text=element_text(size=9), 
        legend.title=element_text(size=11), legend.box = 'vertical',
        axis.text.x = element_text(angle=90, size=10, hjust=1, vjust=0.5), 
        axis.text.y=element_text(size=10), 
        axis.title=element_text(size=12), title=element_text(size=12))

################################################################################
## Feature assay plot 3
################################################################################

plot.data <- data.frame(Y=assay(se, 6)[6307L,], row.names = colnames(se))
plot.data$X <- assay(se, 6)[15644L,];
plot.data <- subset(plot.data, !is.na(X) & !is.na(Y));

# Receiving point selection
selected_pts <- mgcv::in.out(all_lassos[['redDimPlot1']], cbind(as.numeric(all_coordinates[['redDimPlot1']]$X), as.numeric(all_coordinates[['redDimPlot1']]$Y)))
plot.data$SelectBy <- rownames(plot.data) %in% rownames(all_coordinates[['redDimPlot1']])[selected_pts]

# Saving data for transmission
all_coordinates[['featAssayPlot3']] <- plot.data

# Creating the plot
ggplot() +
  geom_point(aes(x = X, y = Y), alpha=1, data=subset(plot.data, !SelectBy), color='#000000', size=1) +
  geom_point(aes(x = X, y = Y), alpha=1, data=subset(plot.data, SelectBy), color="#CC00FF", size=1) +
  labs(x = "Rorb (logcounts)", y = "Foxp2 (logcounts)", title = "Foxp2 vs Rorb") +
  coord_cartesian(xlim = range(plot.data$X, na.rm = TRUE),
                  ylim = range(plot.data$Y, na.rm = TRUE), expand = TRUE) +
  theme_bw() +
  theme(legend.position = 'bottom', legend.text=element_text(size=9), legend.title=element_text(size=11),
        axis.text=element_text(size=10), axis.title=element_text(size=12), title=element_text(size=12)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), color='#00a65a', alpha=0, 
            data=do.call(data.frame, all_brushes[['featAssayPlot3']][c('xmin', 'xmax', 'ymin', 'ymax')]), inherit.aes=FALSE)

################################################################################
## Column data plot 1
################################################################################

plot.data <- data.frame(Y = colData(se)[,"NREADS"], row.names=colnames(se));
plot.data$X <- factor(character(ncol(se)))
plot.data <- subset(plot.data, !is.na(X) & !is.na(Y));

# Receiving point selection
selected_pts <- shiny::brushedPoints(all_coordinates[['featAssayPlot3']], all_brushes[['featAssayPlot3']])
plot.data$SelectBy <- rownames(plot.data) %in% rownames(selected_pts);

# Saving data for transmission
all_coordinates[['colDataPlot1']] <- plot.data

# Setting up plot coordinates
plot.data$GroupBy <- plot.data$X;
set.seed(100);
plot.data$jitteredX <- vipor::offsetX(plot.data$Y,
                                      x=plot.data$X, width=0.4, varwidth=FALSE, adjust=1,
                                      method='quasirandom', nbins=NULL) + as.integer(plot.data$X);

# Creating the plot
ggplot() +
  geom_violin(aes(x = X, y = Y, group = GroupBy), alpha = 0.2, data=plot.data, scale = 'width', width = 0.8) +
  geom_point(aes(y = Y, x = jitteredX), subset(plot.data, !SelectBy), alpha = 0.10, color='#000000', size=1) +
  geom_point(aes(y = Y, x = jitteredX), subset(plot.data, SelectBy), color='#000000', size=1) +
  labs(x = "", y = "NREADS", title = "NREADS ") +
  coord_cartesian(ylim = range(plot.data$Y, na.rm=TRUE), expand = TRUE) +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  theme(legend.position = 'bottom', legend.text=element_text(size=9), 
        legend.title=element_text(size=11), legend.box = 'vertical',
        axis.text.x = element_text(angle=90, size=10, hjust=1, vjust=0.5), 
        axis.text.y=element_text(size=10), 
        axis.title=element_text(size=12), title=element_text(size=12))

################################################################################
## Heat map 1
################################################################################

value.mat <- as.matrix(assay(se, 6)[c(15644L, 6307L, 17992L, 6306L, 6309L), , drop=FALSE]);
value.mat <- t(scale(t(value.mat), center = TRUE, scale = FALSE));
plot.data <- reshape2::melt(value.mat, varnames = c('Y', 'X'));

plot.data[['OrderBy1']] <- colData(se)[['NREADS']][match(plot.data$X, rownames(colData(se)))];
plot.data <- dplyr::arrange(plot.data, OrderBy1);
plot.data$X <- factor(plot.data$X, levels = unique(plot.data$X));

# Creating the heat map
p0 <- ggplot(plot.data, aes(x = X, y = Y)) +
  geom_raster(aes(fill = value)) +
  labs(x='', y='') +
  scale_fill_gradientn(colors=c('blue','white','orange'), 
                       values=c(0,0.519232194351143,1), 
                       limits=c(-10.4982053701138,9.72050503406837), na.value='grey50') +
  scale_y_discrete(expand=c(0, 0)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line=element_blank());
heatlegend <- cowplot::get_legend(p0 + theme(legend.position='bottom'));

# Adding annotations
legends <- list()

p1 <- ggplot(plot.data, aes(x = X, y = 1)) +
  geom_raster(aes(fill = OrderBy1)) +
  labs(x='', y='') +
  scale_y_continuous(breaks=1, labels='NREADS') +
  scale_fill_gradientn(colors=colDataColorMap(colormap, 'NREADS', discrete=FALSE)(21L), na.value='grey50', name='NREADS') +
  theme(axis.text.x=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(), 
        rect=element_blank(), line=element_blank(), axis.title.y=element_blank(), 
        plot.margin = unit(c(0,0,-0.5,0), 'lines'));
legends[[1]] <- cowplot::get_legend(p1 + theme(legend.position='bottom', plot.margin = unit(c(0,0,0,0), 'lines')))

# Laying out the grid
cowplot::plot_grid(
  cowplot::plot_grid(
    p1 + theme(legend.position='none'),
    p0 + theme(legend.position='none'),
    ncol=1, align='v', rel_heights=c(0.1, 1)), 
  heatlegend, ncol=1, rel_heights=c(0.9, 0.1))

################################################################################
## To guarantee the reproducibility of your code, you should also
## record the output of sessionInfo()
sessionInfo()




