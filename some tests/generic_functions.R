#!/usr/bin/env Rscript

# Author: Suhaib Mohammed

## all generic functions used for meta-analysis of gene expression

### average technical replicates blue prints
technical_replicate_average_bp<-function(exp,expAcc,eCounts) {
  
  cat("\n======================")
  cat ("\nProcessing",expAcc,"\n") 
  cat("======================")
  
  if (any(c("technical_replicate") %in% colnames(colData(exp)))) {
    allCounts <-vector()
    allCombCounts<-vector()
    # techical replicates
    assayGroups<- colData(exp)[,grep("AtlasAssayGroup", colnames(colData(exp)))]
    celltypes<- unique(colData(exp)[,grep("cell_type", colnames(colData(exp)))])
    
    # for each individual cell type
    for (cell_type in celltypes) {
      ## check across all the asssociated cell types that technical rep occurs
      cell_type.df=colData(exp)[colData(exp)[,"cell_type"]==cell_type,]
      avgCounts <- vector()
      avgTechRepCounts<-vector()
      combCounts <-vector()
      
      ## combining techincal reps
      if (any(grep("^g", cell_type.df$technical_replicate))) {
        
        # aggregating multiple technical tissyes for the same donor
        # for each individual tissue for a donor 
        TechRepGroups<-unique(grep("^g", cell_type.df$technical_replicate,value=TRUE))
        
        for (group in TechRepGroups) {
          ## averaging by taking median expression value for each replicate
          cat("\ntechnical rep found in cell_type - ", cell_type)
         #cat("\naveraging -",expAcc,"-", cell_type, "-", group) 
          avgTechRepCounts<- as.matrix(apply(eCounts[,rownames(cell_type.df[which(cell_type.df$technical_replicate==group),])],1,median))
          colnames(avgTechRepCounts)<- paste(expAcc, group, cell_type, sep="_")
          avgCounts<-cbind(avgCounts,avgTechRepCounts)
          #print(paste0("AvgCounts - colnames",colnames(avgCounts)))
        }
      }
      ## combining no tech reps
      if (any(cell_type.df$technical_replicate == "  ")) {
        ## no technical group each assay in an biological replicate
        noTechRepAssayGroups<-cell_type.df[which(cell_type.df$technical_replicate == "  "),"AtlasAssayGroup"]
        #cat("\ncombining experiments that has no technical replicates-",cell_type,"-", noTechRepAssayGroups)
        noTechRepAssayNames<-rownames(cell_type.df[which(cell_type.df$technical_replicate == "  "),])
        #cat("\ncombining assays that has no technical replicates-",  noTechRepAssayNames)
        combCounts<-as.matrix(eCounts[,noTechRepAssayNames])
        #cat("\ncombining -",expAcc,"-", noTechRepAssayNames, cell_type,"\n")
        colnames(combCounts)<- paste(expAcc, noTechRepAssayNames,  cell_type, sep="_")
        #print(paste0("ComCounts - colnames",colnames(combCounts)))
      }
      ## combine non-technical assays with averaged technical replicate assays. 
      ## if any tehnical replicates found
      # cat("\ncombining technical and non-tech replicates\n")
      allCombCounts <- cbind(avgCounts,combCounts)
      allCounts<-cbind(allCounts,allCombCounts)
      #print(paste0("dim =", dim(allCounts)))
      #print(paste0("AllCounts - colnames",colnames(allCounts)))
    }
  }
  
  ## combine non-technical assays with averaged technical replicate assays. 
  else if (any(c("individual") %in% colnames(colData(exp)))){
    cat("\nNo technical rep for", expAcc)
    #cat("\ncombining experiment that has no technical replicates")
    noTechRepAssayNames<-rownames(colData(exp))
    noTechRepAssayTissue<-colData(exp)[noTechRepAssayNames,]$cell_type
    allCounts<-as.matrix(eCounts[,noTechRepAssayNames])
    #cat("\ncombining -",expAcc,"-", noTechRepAssayNames, "-", noTechRepAssayTissue,"\n")
    colnames(allCounts)<- paste(expAcc, noTechRepAssayNames,  noTechRepAssayTissue, sep="_")
  }
  return(allCounts)
}  
#################################################

### average technical replicates
technical_replicate_average_gtex<-function(exp,expAcc) {
  
  cat("\n======================")
  cat ("\nProcessing",expAcc,"\n") 
  cat("======================")
  
  if (any(c("individual") %in% colnames(colData(exp)))) {
    allCounts <-vector()
    avgTechRepCounts<-vector()
    combCounts <-vector()
    avgCounts <- vector()
    
    # techical replicates
    assayGroups<- colData(exp)[,grep("AtlasAssayGroup", colnames(colData(exp)))]
    individuals<- unique(colData(exp)[,grep("individual", colnames(colData(exp)))])
    
    # for each individual donor
    for (donor in individuals) {
      ## check across all the donors which tisses have technical replicates and average them
      donor.tissues.df<-colData(exp)[colData(exp)[,"individual"]==donor,]
      donor.tissues<-donor.tissues.df[,"organism_part"]
      
      if (any(table(donor.tissues)>1)) {
        # aggregating multiple technical tissyes for the same donor
        # for each individual tissue for a donor 
        for (tissue in names(which(table(donor.tissues)>1))) {
          ## averaging by taking median expression value for each replicate
          #cat("\ntechnical rep found in donor - ", donor)
          #cat("\naveraging -",expAcc,"-", donor, "-", tissue) 
          avgTechRepCounts<- as.matrix(apply(eCounts[,rownames(donor.tissues.df[which(tissue==donor.tissues),])],1,median))
          colnames(avgTechRepCounts)<- paste(expAcc, donor,  tissue, sep="_")
          avgCounts<-cbind(avgCounts,avgTechRepCounts)
        }
      }
      
      ## combining assays that thas has no technical replicates
      # assay names for single technical rep tissues for the donor
      ## combining expriments that thas has no technical replicates
      if (any(table(donor.tissues)==1)) {
        noTechRepAssayGroups<-names(which(table(donor.tissues)==1))
        #cat("\ncombining experiments that has no technical replicates-", noTechRepAssayGroups)
        noTechRepAssayNames<-rownames(donor.tissues.df[match(names(which(table(donor.tissues)==1)), donor.tissues),])
        noTechRepAssayTissue<-donor.tissues.df[noTechRepAssayNames,]$organism_part
        combCounts<-as.matrix(eCounts[,noTechRepAssayNames])
        colnames(combCounts)<-noTechRepAssayNames
        #cat("\ncombining -",expAcc,"-", donor, "-", noTechRepAssayTissue,"\n")
        colnames(combCounts)<- paste(expAcc, noTechRepAssayNames,  noTechRepAssayTissue, sep="_")
      }
      
      
      ## combine non-technical assays with averaged technical replicate assays. 
      ## if any tehnical replicates found
      #cat("\ncombining technical replicates\n")
      allCounts <- cbind(avgCounts,combCounts)
    }
  }
  ## combine non-technical assays with averaged technical replicate assays. 
  else {
    cat("\nNo technical rep for", expAcc)
    #cat("\ncombining experiment that has no technical replicates")
    noTechRepAssayNames<-rownames(colData(exp))
    noTechRepAssayTissue<-colData(exp)[noTechRepAssayNames,]$organism_part
    allCounts<-as.matrix(eCounts[,noTechRepAssayNames])
    #cat("\ncombining -",expAcc,"-", noTechRepAssayNames, "-", noTechRepAssayTissue,"\n")
    colnames(allCounts)<- paste(expAcc, noTechRepAssayNames,  noTechRepAssayTissue, sep="_")
  }
  return(allCounts)
}  
#################################################
## function to merge two unsymmetrical matrices by row
mergeX<- function(a,b){
  if(!is.matrix(a)){
    print(paste("a is not a matrix"))
  }
  if(!is.matrix(b)){
    print(paste("b is not a matrix")) 
  }
  else 
    comb<-merge(a, b, all=TRUE, by=0)
    exp <- comb[,-1]
    rownames(exp) <- comb[,1]
  
    # make NA to "0"
    exp[is.na(exp)]<-0
  return(exp)
}

#################################################
## function for color consistency based on the tissue types
colorOrder <-function(x) {
  colors.tab<-cbind(names(table(x)),rainbow(length(table(x))))
  colors.order<-matrix(nrow=length(x), ncol=1, dimnames=list(x,"color"))
  
  for (i in 1:nrow(colors.tab)) {
    colors.order[grep(colors.tab[i,1],x)]<-colors.tab[i,2]
  }
  return(colors.order)
}   
#################################################
## functions to plot heatmap
plot_heatmap<-function(exp, name) {
   suppressMessages(require("RColorBrewer"))
  suppressMessages (require("ggplot2"))

  exp<-log10(exp + 1)
  
  prDatTall <- data.frame(sample = rep(colnames(exp), each = nrow(exp)),  
                          gene = rownames(exp), expression = as.vector(exp))
  
  jBuPuFun <- colorRampPalette(brewer.pal(n = 9, "BuPu"))
  paletteSize <- 256
  jBuPuPalette <- jBuPuFun(paletteSize)
  
  ggplot(prDatTall, aes(x = sample, y = gene, fill = expression )) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_blank()) + ylab("") + xlab("") +
    geom_tile() + 
    ggtitle(paste0(name," log expression"), subtitle = NULL) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    scale_fill_gradient2(low = jBuPuPalette[1],
                         mid = jBuPuPalette[paletteSize/2],
                         high = jBuPuPalette[paletteSize],
                         midpoint = (max(prDatTall$expression) + min(prDatTall$expression)) / 2,
                         name = "Expression")
}
#################################################
## summarising gene expression across tissues and cell types
summary_tissues<-function(all.norm){
  x<-sapply(strsplit(colnames(all.norm), split="_"),"[",3)
  names <- colnames(all.norm)

  # experiment detailed expression values
  expNames <- t(as.data.frame(sapply(names, function(x) strsplit(x, "_"))))
  expNormData <- all.norm
  colnames(expNormData) <- paste(expNames[,1], expNames[,3], sep="_")

  # summarisation (aggregation) of expression values by median across tissue/cell types
  expData <-expNormData
  tissue.names<-sapply(strsplit(colnames(expData),split="_"),"[",2)
  colnames(expData) <- tissue.names
  expData <- t(apply(expData, 1, function(x) tapply(x, colnames(expData), median)))

  return (expData)
}
#################################################
## histogram density plot
plot.hist.dens <- function(exp)
{
  exp<-log(exp+1)
  junk.x = NULL
  junk.x = c(junk.x, density(exp))

  xr <- range(junk.x$x)
  yr <- range(junk.x$y)
  
  plot(density(exp), xlim = xr, ylim=yr, main = "Density Log-RawCounts", col="red", lwd=2, xlab="log(expression counts)")
  legend("topright",c("raw-counts"), col=c("red"),lwd=2)
}