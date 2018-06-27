######## User-defined variables ########

#dataPath <- "~/scClustViz_files/mouse/1and2_Jun7/MouseLiver_forViz.RData"
dataPath <- "~/scClustViz_files/mouse/Nov28_Dec6_Healthy_BRCAKnockout_aggregate/Mouse1_forViz.RData"
##  ^ Point this to the output file from PrepareInputs.R
##  If you set a default resolution in the Shiny app, it will save to the same directory.

vizScriptPath <- "~/scClustViz/" 
##  ^ Point this to the directory in which the "app.R" Shiny script resides

species <- "mouse" 
##  ^ Set species ("mouse"/"human").  
##  If other, add the annotation database from Bioconductor to the egDB <- switch() expression below.

#### List known cell-type markers ####

cellMarkers <- list("Other"=c(), "TNK"=c("Cd3d","Cd68","Nkg7"))
##  ^ If you have canonical marker genes for expected cell types, list them here (see example above).
##  The Shiny app will attempt to label clusters in the tSNE projection by highest median gene expression.
##  Otherwise leave the list blank.

########################################



######## Code to run the Shiny app ########
library(shiny)
library(cluster)
library(gplots)
library(scales)
library(viridis)
library(RColorBrewer)
library(TeachingDemos)

egDB <- switch(species,
               mouse={ library(org.Mm.eg.db); "org.Mm.eg.db" },
               human={ library(org.Hs.eg.db); "org.Hs.eg.db" },
               stop("
Set species please!  
If not mouse/human, add your species' annotation database from Bioconductor:  
source('https://bioconductor.org/biocLite.R')
biocLite('org.Xx.eg.db')
"))

rainbow2 <- function(n,a=1) {
  require(scales)
  hues = seq(15, 375, length = n + 1)
  alpha(hcl(h = hues, l = 60, c = 100)[1:n],a)
}

if (length(cellMarkers) < 1) {
  cellMarkersS <- cellMarkersU <- list()
} else {
  cellMarkersS <- apply(combn(seq_along(cellMarkers),2),2,function(X) do.call(intersect,unname(cellMarkers[X])))
  try(names(cellMarkersS) <- apply(combn(seq_along(cellMarkers),2),2,function(X) paste(X,collapse="&")),silent=T)
  cellMarkersS <- cellMarkersS[sapply(cellMarkersS,length) > 0]
  cellMarkersU <- lapply(cellMarkers,function(X) X[!X %in% unlist(cellMarkersS)])
}

load(dataPath)
temp_dataPath <- strsplit(dataPath,"/|\\\\")
dataPath <- sub(temp_dataPath[[1]][length(temp_dataPath[[1]])],"",dataPath)
if (dataPath == "") { dataPath <- "./" }
dataTitle <- sub("\\..+$|_forViz\\..+$","",temp_dataPath[[1]][length(temp_dataPath[[1]])])
rm(temp_dataPath)

if (file.exists(paste0(dataPath,dataTitle,"_savedRes.RData"))) {
  load(paste0(dataPath,dataTitle,"_savedRes.RData"))
} else {
  savedRes <- NULL
}

silDist <- dist(dr_clust,method="euclidean")  
##  ^ precalculating distances in reduced dimensionality space for the silhouette plot.

for (l in names(CGS)) {
  for (i in names(CGS[[l]])) {
    CGS[[l]][[i]]$MTCrank <- rank(CGS[[l]][[i]]$MTC,ties.method="min")/nrow(CGS[[l]][[i]])
    CGS[[l]][[i]]$cMu <- rownames(CGS[[l]][[i]]) %in% unlist(cellMarkersU)
    CGS[[l]][[i]]$cMs <- rownames(CGS[[l]][[i]]) %in% unlist(cellMarkersS)
    CGS[[l]][[i]]$overCut <- CGS[[l]][[i]]$MTC > mean(CGS[[l]][[i]]$MTC)
    CGS[[l]][[i]]$genes <- rownames(CGS[[l]][[i]])
  }
}

if (length(cellMarkers) < 1) {
  clusterID <- sapply(colnames(cl),function(X) rep("",nrow(cl)),simplify=F)
} else {
  clusterID <- sapply(CGS,function(Z) {
    temp <- names(cellMarkers)[sapply(Z,function(Y) 
      which.max(sapply(cellMarkers,function(X) median(Y$MTC[rownames(Y) %in% X]))))]
    names(temp) <- names(Z)
    return(temp)
  },simplify=F)
}

#### Run the Shiny App!  ####
runApp(vizScriptPath)
