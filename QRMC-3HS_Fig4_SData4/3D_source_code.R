
#'
plotUpdown <- function(data2plot,contrast,plot.title=NULL,angle=0,
                       hjust=0.5,vjust=0.5,
                       col.up=NULL,
                       col.down=NULL){
  data2plot$contrast <- factor(data2plot$contrast,levels = unique(data2plot$contrast))
  data2plot$regulation <- factor(data2plot$regulation,levels = c('up-regulated','down-regulated'))
  g <- ggplot(data2plot,aes(x=contrast,y=number,group=regulation,
                            fill=regulation,label=number))+
    geom_bar(stat='identity')+
    geom_text(size = 3, position = position_stack(vjust = 0.5))+
    theme_bw()+
    labs(title=plot.title)
  g <- g+theme(axis.text.x = element_text(angle = angle, hjust = hjust,vjust = vjust))
  if(!is.null(col.up) & !is.null(col.down))
    g <- g+scale_fill_manual(values=c(col.up,col.down))
  g
}

######################################################################################################################
#' Plot gene and transcript expression profiles
#' @param data.exp a data.frame of transcript level expression, e.g. read counts and TPM.
#' @param gene a gene name to plot.
#' @param mapping a data.frame of transcript-gene mapping (first column is transcript list and second column is gene list).
#' @param genes.ann a data.frame of gene annotation with first column "target" and second column "annotation". Default is \code{NULL}.
#' @param trans.ann a data.frame of transcript annotation with first column "target" and second column "annotation".
#' Default is \code{NULL}.
#' @param trans.expressed a vector of expressed transcripts. If not \code{NULL}, only the profiles of provided transcripts
#' will be shown in the plot.
#' @param trans.highlight a vector of transcripts. If not \code{NULL}, all the other transcripts which are not in this vector will be hidden in the plot.
#' @param reps a vecotr of replicate labels, which provide the grouping information of biolocial replciates to
#' calculate the average expression in each conditions.
#' @param groups a vecotr of grouping labels for \code{reps}, which provide the grouping information to slice the profile plot into different segments, e.g. samples in "treatment"
#' vs "Control", samples in different development stages, etc. The vector length of \code{groups} is the same to the length of \code{reps}.
#' @param sliceProfile logical, whether to slice profile plot according to \code{groups} labels.
#' @param x.lab,y.lab characters for x-axis and y-axis labels, respectively.
#' @param marker.size a number passed to \code{geom_point} to control the point size in the plot.
#' @param legend.text.size a number to control the legend text size.
#' @param error.type error type to make the error bars. Options are: "stderr" for standard error and
#' "sd" for standard deviation.
#' @param error.bar logical, whether to show error bars on the profile plot.
#' @param show.annotation logical, whether to show the annotations provided in \code{genes.ann} and \code{trans.ann} on the plot.
#' @param plot.gene logical, whether to show the gene level expression on the plot. The gene expression is the total of all the
#' transcript expression.
#' @details \code{plotAbundance} is used to plot the gene and/or trnascript level abundance while \code{plotPS} is used to plot
#' the percent spliced (PS) of transcript. PS is defined as the ratio of transcript expression to the total of all transcript
#' expression (i.e. gene expression) in the same gene.
#'
#' @return a plot in \code{ggplot} format.
#' @export
#' @rdname plot.abundance
#' @examples
#' data(exp.data)
#' plotAbundance(data.exp = exp.data$trans_TPM,
#'               gene = 'AT1G01020',
#'               mapping = exp.data$mapping,
#'               genes.ann = exp.data$genes.ann,
#'               trans.ann = exp.data$trans.ann,
#'               reps = exp.data$samples$condition)
#'
#' plotPS(data.exp = exp.data$trans_TPM,
#'               gene = 'AT1G01020',
#'               mapping = exp.data$mapping,
#'               genes.ann = exp.data$genes.ann,
#'               trans.ann = exp.data$trans.ann,
#'               reps = exp.data$samples$condition)
#' @seealso \code{\link{data.error}}
plotAbundance<- function(
  data.exp,
  gene,
  mapping,
  genes.ann=NULL,
  trans.ann=NULL,
  trans.expressed=NULL,
  trans.highlight=NULL,
  reps,
  groups=NULL,
  sliceProfile=F,
  x.lab='Conditions',
  y.lab='TPM',
  marker.size=3,
  legend.text.size=11,
  error.type='stderr',
  error.bar=T,
  show.annotation=T,
  plot.gene=T
){
  #####################################################################
  ##prepare plot data
  #####################################################################
  rownames(mapping) <- mapping$TXNAME
  if(!is.null(trans.expressed)){
    data.exp <- data.exp[trans.expressed,]
    mapping <- mapping[trans.expressed,]
  }
  rownames(genes.ann) <- genes.ann$target
  rownames(trans.ann) <- trans.ann$target
  trans <- mapping$TXNAME[mapping$GENEID==gene]
  plot.title <- paste0('Gene: ',gene)
  # expression.sum <- if(length(trans.idx)==1) sum(data.exp[trans.idx,]) else rowSums(data.exp[trans.idx,])
  # trans.idx <- trans.idx[expression.sum>0]

  if(length(trans)==0){
    message(paste0(gene, ' is not in the dataset'))
    return(NULL)
  }

  data2plot<-data.exp[trans,,drop=F]

  if(length(trans)==1)
    data2plot <- rbind(data2plot,data2plot) else data2plot <- rbind(colSums(data2plot),data2plot)

  genes.idx <- gene
  trans.idx <- trans
  if(show.annotation){
    if(gene %in% genes.ann$target){
      idx <- genes.ann[gene,'annotation']
      genes.idx <- paste0(gene,' (',idx,')')
    }

    if(any(trans %in% trans.ann$target)){
      trans.idx <- trans
      idx <- which(trans %in% trans.ann$target)
      trans.idx[idx] <- paste0(trans.idx[idx],' (',trans.ann[trans.idx[idx],'annotation'],')')
    }
  }
  rownames(data2plot) <- c(genes.idx,trans.idx)

  ##--mean and error
  mean2plot <- t(rowmean(t(data2plot),group = reps))
  sd2plot <- by(t(data2plot),INDICES = reps,FUN = function(x){
    apply(x,2,function(y){
      data.error(y,error.type = error.type)
    })
  })
  sd2plot <- do.call(cbind,sd2plot)
  sd2plot <- sd2plot[rownames(mean2plot),]

  mean2plot <- reshape2::melt(mean2plot)
  sd2plot <- reshape2::melt(sd2plot)
  colnames(mean2plot) <- c('Targets','Conditions','mean')
  colnames(sd2plot) <- c('Targets','Conditions','error')
  data2plot <- merge(mean2plot,sd2plot)

  ##--plot gene or not
  if(!plot.gene){
    data2plot <- data2plot[-which(data2plot$Targets %in% genes.idx),]
    plot.title <- paste0('Gene: ',genes.idx)
  }

  data2plot$Conditions <- factor(data2plot$Conditions,levels = unique(reps))
  if(!is.null(groups)){
    group.idx <- data.frame(reps=reps,groups=groups)
    idx2check <- as.vector(unlist(by(data = group.idx,INDICES = group.idx$reps,FUN = function(x){
      length(unique(x$groups))
    })))
    if(any(idx2check>1)){
      groups <- NULL
    } else {
      group.idx <- group.idx[!duplicated(group.idx$reps),]
      if(length(unique(group.idx$groups))==1){
        groups <- NULL
      } else {
        rownames(group.idx) <- group.idx$reps
        data2plot$Group <- factor(group.idx[data2plot$Conditions,'groups'])
      }
    }
  }

  if(!is.null(trans.highlight)){
    if(plot.gene){
      data2plot <- droplevels(data2plot[data2plot$Targets %in% c(genes.idx,trans.highlight),])
    } else {
      data2plot <- droplevels(data2plot[data2plot$Targets %in% trans.highlight,])
    }
  }

  legend.ncol <- ceiling(length(unique(data2plot$Targets))/15)

  profiles <- ggplot(data2plot,aes(x=Conditions,y=mean,group=Targets,color=Targets))+
    geom_line(size=1)+
    geom_point(size=marker.size,aes(fill=Targets,shape=Targets))+
    scale_shape_manual(name="Targets",values=c(25:0,25:0,rep(25,500)))+
    labs(x=x.lab,y=y.lab,title=plot.title)+
    theme_bw()+
    theme(panel.grid = element_blank(),legend.text=element_text(color='black',size=legend.text.size))+
    guides(
      fill=guide_legend(ncol = legend.ncol),
      shape=guide_legend(ncol = legend.ncol),
      color=guide_legend(ncol = legend.ncol))
  if(error.bar)
    profiles <- profiles+
    geom_errorbar(aes(ymin=mean-error,ymax=mean+error),width=0.1,color='black',size=0.3)
  if(sliceProfile & !is.null(groups))
    profiles <- profiles+facet_grid(.~Group,scales = 'free_x')
  profiles
}

#' @rdname plot.abundance
#' @export
plotPS <- function(
  data.exp,
  gene,
  mapping,
  genes.ann=NULL,
  trans.ann=NULL,
  trans.expressed=NULL,
  reps,
  groups=NULL,
  sliceProfile=F,
  y.lab='TPM',
  x.lab='Conditions',
  marker.size=3,
  legend.text.size=11,
  show.annotation=T
){
  #####################################################################
  ##prepare plot data
  #####################################################################
  rownames(mapping) <- mapping$TXNAME
  if(!is.null(trans.expressed)){
    data.exp <- data.exp[trans.expressed,]
    mapping <- mapping[trans.expressed,]
  }
  rownames(genes.ann) <- genes.ann$target
  rownames(trans.ann) <- trans.ann$target
  trans <- mapping$TXNAME[mapping$GENEID==gene]
  # expression.sum <- if(length(trans.idx)==1) sum(data.exp[trans.idx,]) else rowSums(data.exp[trans.idx,])
  # trans.idx <- trans.idx[expression.sum>0]

  if(length(trans)==0){
    message(paste0(gene, ' is not in the dataset'))
    return(NULL)
  }
  data2plot<-data.exp[trans,,drop=F]
  data2plot <- t(rowmean(t(data2plot),reps))

  if(length(trans)==1)
    data2plot <- rbind(data2plot,data2plot) else data2plot <- rbind(colSums(data2plot),data2plot)

  genes.idx <- gene
  trans.idx <- trans
  if(show.annotation){
    if(gene %in% genes.ann$target){
      idx <- genes.ann[gene,'annotation']
      genes.idx <- paste0(gene,' (',idx,')')
    }

    if(any(trans %in% trans.ann$target)){
      trans.idx <- trans
      idx <- which(trans %in% trans.ann$target)
      trans.idx[idx] <- paste0(trans.idx[idx],' (',trans.ann[trans.idx[idx],'annotation'],')')
    }
  }
  rownames(data2plot) <- c(genes.idx,trans.idx)
  legend.ncol <- ceiling(length(trans.idx)/15)
  plot.title <- paste0('Gene: ',genes.idx)
  ##PS value
  data2plot[trans.idx,] <- t(t(data2plot[trans.idx,])/data2plot[genes.idx,])
  data2plot <- data2plot[-which(rownames(data2plot) %in% genes.idx),,drop=F]
  data2plot <- reshape2::melt(data2plot)
  colnames(data2plot) <- c('Targets','Conditions','PS')
  data2plot$Conditions <- factor(data2plot$Conditions,levels = unique(reps))

  if(!is.null(groups)){
    group.idx <- data.frame(reps=reps,groups=groups)
    idx2check <- as.vector(unlist(by(data = group.idx,INDICES = group.idx$reps,FUN = function(x){
      length(unique(x$groups))
    })))
    if(any(idx2check>1)){
      groups <- NULL
    } else {
      group.idx <- group.idx[!duplicated(group.idx$reps),]
      if(length(unique(group.idx$groups))==1){
        groups <- NULL
      } else {
        rownames(group.idx) <- group.idx$reps
        data2plot$Group <- factor(group.idx[data2plot$Conditions,'groups'])
      }
    }
  }
  profiles <- ggplot(data2plot,aes(x=Conditions,y=PS,group=Targets,color=Targets))+
    geom_bar(stat='identity',aes(fill=Targets))+
    labs(x=x.lab,y=y.lab,title=plot.title)+
    theme_bw()+
    theme(panel.grid = element_blank(),legend.text=element_text(color='black',size=legend.text.size))+
    guides(
      fill=guide_legend(ncol = legend.ncol),
      shape=guide_legend(ncol = legend.ncol),
      color=guide_legend(ncol = legend.ncol))

  if(sliceProfile & !is.null(groups))
    profiles <- profiles+facet_grid(.~Group,scales = 'free_x')
  profiles
}

######################################################################################################################
#' Bar plot of GO annotation
#' @param go.table a data.frame of GO annotation table, with first column "Category" (i.e. BP, CC and MF),
#' second column of Go "Term" and remaining columns with significant statistics, e.g. "Count", "-log10(FDR)", etc.
#' @param col.idx a character to indicate which column of statistics to use to report the significance, e.g.
#' \code{col.idx="-log10(FDR)"}. \code{col.idx} must match to the column names of \code{go.table}.
#' @param plot.title the titile to show on the plot.
#' @return a plot in \code{ggplot} format.
#' @export
#' @examples
#' data(go.table)
#' plotGO(go.table = go.table,col.idx = '-log10(FDR)',plot.title = 'GO annotation: DE genes')
#'
plotGO <- function(go.table,col.idx,plot.title='GO annotation',go.text.cut=50,legend.position='right'){
  class(go.table) <- c("GO", class(go.table))
  data2plot <- go.table[,c('Category','Term',colnames(go.table)[which(colnames(go.table)==col.idx)])]
  colnames(data2plot) <- c('Category','Term','Value')
  data2plot <- by(data2plot,data2plot$Category,function(x){
    x[order(x$Value,decreasing = T),]
  })
  data2plot <- do.call(rbind,data2plot)
  idx<-sapply(data2plot$Term,function(x){
    x<-as.vector(x)
    if(nchar(x)>go.text.cut)
      x<-paste0(substr(x,1,go.text.cut),'...')
    x
  })
  data2plot$Term <- factor(idx,levels = rev(idx))
  g <- ggplot(data2plot,aes(x=Term,y=Value))+
    geom_bar(stat='identity',aes(fill=Category))+
    coord_flip() +
    theme_bw()+
    theme(axis.title.y = element_blank(),legend.position=legend.position)+
    facet_grid(Category~.,scales = 'free_y',space='free_y')+
    labs(y=col.idx,title=plot.title)
  g
}

######################################################################################################################
#' Make principal component analysis (PCA) plot of individual samples
#' @param data2pca a matrix/data.frame to make the PCA plot of which the rows are defined individuals to plot.
#' @param dim1,dim2 characters to indicate the PCs to plot on x-axis (dim1) and y-axis (dim2), e.g. \code{dim1='PC1'},
#' \code{dim2='PC2'}; \code{dim1='PC2'}, \code{dim2='PC3'}; etc.
#' @param groups a vector of characters to specify the groups of conditions. The scatter points will be coloured
#' according to provided group information.
#' @param plot.title the title to show on the plot.
#' @param ellipse.type add colour shadow to the sample groups. Options are: "none","ellipse" and "polygon".
#' @param add.label logical, whether to add sample labels to the scatter points.
#' @param adj.label logical, whether to adjust the position of sample labels to avoid overlapping.
#' @return a plot in \code{ggplot} format.
#' @export
#' @examples
#' data(exp.data)
#' library(edgeR)
#' library(RUVSeq)
#' library(ggplot2)
#' trans_counts <- exp.data$trans_counts
#'
#' ##------> sample informationsamples <- exp.data$samples
#' samples <- exp.data$samples[,c('condition','brep','srep')]
#'
#' ##------> sum sequencing replicates
#' idx <- paste0(samples$condition,'_',samples$brep)
#' y <- sumarrays(trans_counts,group = idx)
#'
#' ##------> update sample information after sum
#' samples <- samples[samples$srep==samples$srep[1],]
#'
#' ##------> normalisation
#' dge <- DGEList(counts=y)
#' dge <- calcNormFactors(dge)
#' data2pca <- t(counts2CPM(obj = dge,Log = T))
#'
#' ##------> plot PCA
#' groups <- samples$brep
#' plotPCAind(data2pca = data2pca,dim1 = 'PC1',dim2 = 'PC2',
#'              groups = groups,ellipse.type='polygon')
#'
#' ##------> remove batch effects
#' trans_batch <- remove.batch(read.counts = y,
#'                       condition = samples$condition,
#'                       method = 'RUVr')
#'
#' ##------> normalisation and plot PCA again
#' dge <- DGEList(counts=trans_batch$normalizedCounts)
#' dge <- suppressWarnings(calcNormFactors(dge))
#' data2pca <- t(counts2CPM(obj = dge,Log = T))
#' plotPCAind(data2pca = data2pca,dim1 = 'PC1',dim2 = 'PC2',
#'              groups = groups,ellipse.type='polygon')

plotPCAind <- function(data2pca,
                         dim1='PC1',dim2='PC2',
                         groups,
                         plot.title='PCA plot',
                         ellipse.type=c('none','ellipse','polygon'),
                         add.label=T,adj.label=T){

  ellipse.type <- match.arg(ellipse.type,c('none','ellipse','polygon'))
  dim1 <- toupper(dim1)
  dim2 <- toupper(dim2)

  fit <- prcomp(data2pca,scale = T)
  fit.stat <- summary(fit)$importance
  dim1.p <- round(fit.stat[2,dim1]*100,2)
  dim2.p <- round(fit.stat[2,dim2]*100,2)
  data2plot <- data.frame(groups=groups,dim1=fit$x[,dim1],dim2=fit$x[,dim2])
  data2plot$groups <- factor(data2plot$groups,levels = unique(data2plot$groups))
  data2plot$labels <- rownames(data2plot)

  g <- ggplot(data2plot,aes(x=dim1,y=dim2),frame=T)+
    geom_point(aes(colour=groups,shape=groups))+
    geom_hline(yintercept = 0,linetype=2)+
    geom_vline(xintercept=0,linetype=2)+
    theme_bw()+
    theme(panel.border = element_blank())+
    labs(x=paste0(dim1,' (',dim1.p,'%)'),
         y=paste0(dim2,' (',dim2.p,'%)'),
         title=plot.title)+
    scale_shape_manual(values = rep(c(18:0,19:25),3)[1:length(unique(groups))])+
    coord_cartesian(xlim = c(min(data2plot$dim1)*1.1,max(data2plot$dim1)*1.1),
                    ylim = c(min(data2plot$dim2)*1.1,max(data2plot$dim2)*1.1))

  if(ellipse.type=='ellipse')
    g <- g+stat_ellipse(geom = "polygon",aes(fill=groups,colour=groups),alpha=0.2)
  if(ellipse.type=='polygon'){
    hulls <- plyr::ddply(data2plot, "groups", function(x){
      x[chull(x$dim1,x$dim2),]
    })
    g <- g+geom_polygon(data = hulls,aes(x=dim1,y=dim2,fill=groups),alpha = 0.3)
  }
  if(add.label & !adj.label)
    g <- g+geom_text(aes(label=labels,colour=groups),hjust=0, vjust=0)
  if(add.label & adj.label)
    g <- g+ggrepel::geom_text_repel(label=rownames(data2pca),aes(colour=groups))
  g

}


#
#' Boxplot of expression distribution in each sample.
#' @param data.before a numeric matrix of expression before normalisation, e.g. raw read counts.
#' @param data.after a numeric matrix of expression after normalisation, e.g. \eqn{\log_2}-CPM.
#' @param condition a vector of conditions, which is used to distinguish and colour the conditions of
#' samples.
#' @param sample.name a vector of sample names, which is used to name and order the samples in the boxplot.
#' @return a list with two elements: "g1" and "g2", which are \code{ggplot} format boxplots of the data
#' before and after normalisation.
#' @export
#'
boxplotNormalised <- function(data.before,data.after,condition,sample.name){
  data2plot <- data.before
  data2plot <- t(apply(data2plot,2,function(x) boxplot.stats(x)$stats))
  data2plot <- cbind(condition=condition,samples=sample.name,data.frame(data2plot))
  data2plot <- reshape2::melt(data = data2plot,id.vars = c('condition','samples'))
  data2plot$samples <- factor(data2plot$samples,levels = sample.name)
  g1 <- ggplot(data2plot,aes(x=samples,y=value))+
    stat_boxplot(geom = "errorbar", width = 0.3) +
    geom_boxplot(aes(fill=condition))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5,vjust = 0.5))+
    labs(x='samples',y='read counts',title='Data distribution before normalisation')

  data2plot <- data.after
  data2plot <- t(apply(data2plot,2,function(x) boxplot.stats(x)$stats))
  data2plot <- cbind(condition=condition,samples=sample.name,data.frame(data2plot))
  data2plot <- reshape2::melt(data = data2plot,id.vars = c('condition','samples'))
  data2plot$samples <- factor(data2plot$samples,levels = sample.name)
  g2 <- ggplot(data2plot,aes(x=samples,y=value))+
    stat_boxplot(geom = "errorbar", width = 0.3) +
    geom_boxplot(aes(fill=condition))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5,vjust = 0.5))+
    labs(x='samples',y='log2(CPM)',title='Data distribution after normalisation')

  if(length(sample.name)>20){
    g1 <- g1+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())
    g2 <- g2+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())
  }
  return(list(g1=g1,g2=g2))
}

#' Plot Euler diagram
#' @param x (1) a list of individuals in different sets for comparisons or (2) a vecotr of numbers of set relations.
#' @param fill a vector of colours to fill the diagram.
#' @param shape geometric shape used in the diagram.
#' @param ... arguments passed down to the \code{\link{eulerr::euler}} function.
#' @return a Euler diagram
#' @examples
#' x <- letters[1:10]
#' y <- letters[5:15]
#' z <- set2(x,y)
#' combo <- c(x=length(z$x.only),y=length(z$y.only),"x&y"=length(z$xy))
#' plotEulerDiagram(list(x=x,y=y))
#' plotEulerDiagram(combo)
#' @export
#' @seealso \code{\link{eulerr::euler}}.
#'
plotEulerDiagram <- function(x,
                             fill = NULL,
                             shape = c("ellipse", "circle"),scale.plot=1,
                             ...){
  if(any(is.null(fill)))
    fill <- gg.color.hue(length(x))
  shape <- match.arg(shape,c("ellipse", "circle"))
  fit <- eulerr::euler(x,...)
  # grid.newpage()
  g <- plot(fit,quantities = TRUE,
            labels = list(font =1),
            fill=fill,shape = shape)
  vp = viewport(height=unit(1*scale.plot, "npc"),
                width=unit(1*scale.plot, "npc"))
  g <- arrangeGrob(g,vp=vp)
  g
  # grid.newpage()
  # grid.rect(vp=vp,gp=gpar(lty=1, col="white", lwd=0))
  # grid.draw(g)
}

#' Volcano plot
#' @param data2plot a dataframe with a column of "target" names, a column of x coordinates, a column of y coordinates and
#' a column of "Significant" or "Not significant" to distinguish significance.
#' @param xlab the x-axis label.
#' @param ylab the y-axis label.
#' @param title the title of the plot.
#' @param col0 colour of not significant points.
#' @param col1 colour of significant points.
#' @param size the point size on the plot.
#' @param top.n show target labels of top.n distance to (0,0) coordinate.
#' @return a Volcano diagram
#' @export

#'
plotVolcano <- function(data2plot,
                        xlab,
                        ylab,
                        title='Volcano plot',
                        col0='black',col1='red',
                        size=1,
                        top.n=10){
  data2plot$distance <- sqrt((data2plot$x)^2+(data2plot$y)^2)
  if('Significant' %in% data2plot$significance){
    data2label <- data2plot[data2plot$significance=='Significant',]
    data2label <- droplevels(data2label[order(data2label$distance,decreasing = T),][1:top.n,])
    g <- ggplot(data = data2plot,aes(x=x,y=y))+
      geom_point(aes(colour=significance),size=size)+
      scale_color_manual(values=c('Not significant'=col0,'Significant'=col1))+
      theme_bw()+
      labs(x=xlab,y=ylab,title = title)+
      scale_x_continuous(breaks = pretty(data2plot$x, n = 10))
    q <- g+
      geom_label_repel(data=data2label,
                       aes(x=x,y=y,label=target,fill=contrast),
                       alpha=0.6,
                       nudge_y      = -0.35,
                       direction    = "both",
                       segment.color = "black",
                       hjust        = 1,
                       segment.size = 0.2)
  } else {
    g <- ggplot(data = data2plot,aes(x=x,y=y))+
      geom_point(aes(colour=significance),size=size)+
      scale_color_manual(values=c('Not significant'=col0))+
      theme_bw()+
      labs(x=xlab,y=ylab,title = title)+
      scale_x_continuous(breaks = pretty(data2plot$x, n = 10))
    q <- g
  }
  q
}

rowmean<-function(x,group,reorder=F,na.rm=T){
  order.idx<-as.character(unique(group))
  if (reorder)
    order.idx<-gtools::mixedsort(order.idx)

  counts <- table(group)[order.idx]
  sums <- rowsum(x, group = group)[order.idx,]
  means <- (diag(1/counts)) %*% as.matrix(sums)
  rownames(means) <- order.idx
  if (na.rm)
    means[is.na(means)] <- 0
  return(means)
}


remove.batch<-function(read.counts,
                       condition,
                       design=NULL,
                       contrast=NULL,
                       group=NULL,
                       method=c('RUVr','RUVg','RUVs'),
                       cIdx=NULL,
                       k=1,...){
  start.time <- Sys.time()
  method <- match.arg(method,c('RUVr','RUVg','RUVs'))
  read.counts<-as.matrix(round(read.counts,0))

  ##########################################################################
  #---get the negative control
  if(is.null(cIdx) | method=='RUVr'){
    if(is.null(design)){
      condition<-factor(condition,levels = unique(condition))
      design<-model.matrix(~0+condition)
      colnames(design)<-gsub('condition','',colnames(design))
    }
    if(is.null(contrast)){
      contrast <- unique(condition)
      contrast <- paste0(contrast[-1],'-',contrast[1])
    }

    contrast <- makeContrasts(contrasts = contrast, levels=design)
    message('Estimate norm factor...')
    y <- DGEList(counts=read.counts, group=group)
    y <- calcNormFactors(y)
    message('Estimate Common Dispersion for Negative Binomial GLMs ...')
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    message('Fit genewise Negative Binomial Generalized Linear Models and calculate residuals...')
    fit <- glmQLFit(y, design)

    lrt <- glmQLFTest(fit, contrast = contrast)
    top <- topTags(lrt, n=nrow(read.counts))$table
    empirical <- top[order(top$PValue,decreasing = T),]
    empirical <- rownames(empirical)[empirical$PValue>0.1]
  }

  switch(method,
         RUVr = {
           empirical <- rep(TRUE,dim(read.counts)[1])
           message('Remove Unwanted Variation Using RUVr...')
           res <- residuals(fit, type="deviance")
           results <-  RUVr(x = read.counts,cIdx = empirical,residuals = res,k=k)
         },
         RUVg = {
           message('Remove Unwanted Variation Using RUVg...')
           results <- RUVg(x = read.counts,cIdx = empirical,k=k)
         },
         RUVs = {
           message('Remove Unwanted Variation Using RUVs...')
           results <- RUVs(x = read.counts,cIdx = empirical,k=k,scIdx = makeGroups(condition))
         }
  )

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message(paste('Time for analysis:',round(time.taken,3),attributes(time.taken)$units))
  message('Done!!!')
  results$method <- method
  return(results)
}

condition2design <- function(condition,batch.effect=NULL){
  design.data <- data.frame(condition=condition)
  if(!is.null(batch.effect)){
    colnames(batch.effect) <- paste0('batch',1:ncol(batch.effect))
    design.data <- data.frame(design.data,batch.effect)
  }
  design.fomula <- as.formula(paste0('~0+',paste0(colnames(design.data),collapse = '+')))
  design <- model.matrix(design.fomula, data=design.data)
  idx <- grep('condition',colnames(design))
  design.idx <- colnames(design)
  design.idx <- substr(design.idx,nchar('condition')+1,nchar(design.idx))
  colnames(design)[idx] <- design.idx[idx]
  design
}


#limma

limma.pipeline <- function(dge,
                           contrast,
                           span=0.5,
                           design,
                           deltaPS=NULL,
                           diffAS=F,
                           adjust.method='BH',
                           block=NULL,
                           voomWeights=F,...){
  start.time <- Sys.time()
  results <- list()
  if(is.null(dge$genes)){
    diffAS <- F
  } else {
    if(max(table(dge$genes$GENEID)) < 2)
      diffAS <- F
  }

  ##########################################################
  ##--limma voom
  cat('> Limma-voon to estimate mean-vriance trend ...','\n')
  if(voomWeights){
    voom.object<-voomWithQualityWeights(dge,design,plot=F,span = span)
  } else {
    voom.object<-voom(dge,design,plot=F,span = span,...)
  }

  results$voom.object <- voom.object
  targets <- rownames(voom.object$E)

  ##########################################################
  ##--Fit block
  if(!is.null(block)){
    cat('> Fit block information','\n')
    corfit <- duplicateCorrelation(voom.object,design,block=block)
    correlation <- corfit$consensus
    results$corfit <- corfit
  } else {
    correlation <- NULL
  }
  results$block <- block
  results$correlation <- correlation

  ##########################################################
  ##--Fit a basic linear model
  cat('> Fit a basic linear model ...','\n')
  fit.lmFit <- lmFit(voom.object, design,block = block,correlation = correlation)
  results$fit.lmFit <- fit.lmFit

  ##########################################################
  ##--Fit a basic linear model
  cat('> Fit the contrast model ...','\n')
  contrast.matrix <- makeContrasts(contrasts = contrast, levels=design)
  print(paste0('Contrast groups: ',paste0(contrast,collapse = '; ')))
  fit.contrast<-contrasts.fit(fit.lmFit, contrast.matrix)
  results$fit.contrast<-fit.contrast

  ##########################################################
  ##--Fit a eBayes model
  cat('> Fit a eBayes model ...','\n')
  fit.eBayes<-eBayes(fit.contrast)
  results$fit.eBayes<-fit.eBayes

  ##########################################################
  ##--Testing statistics for each contrast group
  cat('> Testing for each contrast group ...','\n')
  DE.pval.list <- lapply(contrast,function(i){
    x <- topTable(fit.eBayes,
                  coef=i,
                  adjust.method =adjust.method,
                  number = Inf)
    x <- x[targets,]
    x
  })
  names(DE.pval.list) <- contrast
  results$DE.pval.list<-DE.pval.list

  ###---DE pval and lfc
  DE.pval <- do.call(cbind,lapply(DE.pval.list,FUN = function(x) x$adj.P.Val))
  DE.lfc <- do.call(cbind,lapply(DE.pval.list,FUN = function(x) x$logFC))
  rownames(DE.pval) <- rownames(DE.lfc) <- targets
  colnames(DE.pval) <- colnames(DE.lfc) <- contrast

  DE.stat <- summaryStat(x = DE.pval,y = DE.lfc,
                          target = rownames(DE.pval),
                          contrast = contrast,
                          stat.type = c('adj.pval','log2FC'))

  # DE.stat <- cbind(DE.pval,DE.lfc)
  # colnames(DE.stat) <- c(paste0('pval:',contrast),paste0('lfc:',contrast))
  results$DE.pval<-DE.pval
  results$DE.lfc<-DE.lfc
  results$DE.stat<-DE.stat

  ##########################################################
  ##--Testing statistics for across all contrast groups
  cat('> Testing across all contrast groups ...','\n')
  DE.stat.overalltest<-topTable(fit.eBayes,number = Inf,
                            coef = contrast,adjust.method =adjust.method )
  # DE.stat.overalltest <- DE.stat.overalltest[targets,]
  col.idx <- gsub('-','.',contrast)
  col.idx <- grep(paste0(col.idx,collapse = '|'),colnames(DE.stat.overalltest))
  DE.stat.overalltest <- data.frame(target=rownames(DE.stat.overalltest),
                                    contrast='overall',DE.stat.overalltest,row.names = NULL)
  colnames(DE.stat.overalltest)[col.idx] <- gsub('[.]','-',colnames(DE.stat.overalltest)[col.idx])
  results$DE.stat.overalltest<-DE.stat.overalltest

  if(diffAS){
    cat('> Fit a splicing model ...','\n')
    if(is.null(deltaPS))
      stop('Please provide deltaPS for DAS analysis...')
    fit.splice<-diffSplice(fit.contrast, geneid = 'GENEID')
    results$fit.splice<-fit.splice

    # ##########################################################
    ##---DTU transcripts
    ##DTU transcript pval list
    genes.idx <- unique(fit.splice$genes$GENEID)
    trans.idx <- unique(fit.splice$genes$TXNAME)

    DTU.pval.list<-lapply(contrast,function(i){
      y<-topSplice(fit.splice, coef=i, test="t", number=Inf, FDR=10000)
      rownames(y)<-y$TXNAME
      z <- y[trans.idx,]
      z
    })
    names(DTU.pval.list) <- contrast

    ##DTU transcript pvals
    DTU.pval <- lapply(contrast,function(i){
      x <- DTU.pval.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })
    DTU.pval <- do.call(cbind,DTU.pval)
    colnames(DTU.pval) <- contrast

    ##DTU transcript deltaPS
    DTU.deltaPS <- deltaPS[rownames(DTU.pval),,drop=F]

    DTU.stat <- summaryStat (x = DTU.pval,y = DTU.deltaPS,
                            target = rownames(DTU.pval),
                            contrast = contrast,
                            stat.type = c('adj.pval','deltaPS'))

    # DTU.stat <- cbind(DTU.pval,DTU.deltaPS)
    # colnames(DTU.stat) <- c(paste0('pval:',contrast),paste0('deltaPS:',contrast))

    results$DTU.pval.list<-DTU.pval.list
    results$DTU.pval<-DTU.pval
    results$DTU.deltaPS<-DTU.deltaPS
    results$DTU.stat<-DTU.stat

    # ##########################################################
    ##---DAS genes
    ##---max deltaPS
    maxdeltaPS <- by(deltaPS[fit.splice$genes$TXNAME,,drop=F],
                     INDICES = fit.splice$genes$GENEID,
                     function(x){
                       apply(x,2,function(i) i[abs(i)==max(abs(i))][1])
                     },simplify = F)
    maxdeltaPS <- do.call(rbind,maxdeltaPS)
    maxdeltaPS <- maxdeltaPS[genes.idx,,drop=F]
    results$maxdeltaPS<-maxdeltaPS

    ##---F test
    DAS.pval.F.list<-lapply(contrast,function(i){
      y<-topSplice(fit.splice, coef=i, test="F", number=Inf, FDR=10000)
      rownames(y)<-y$GENEID
      y[genes.idx,]
    })
    names(DAS.pval.F.list) <- contrast

    DAS.pval.F <- lapply(contrast,function(i){
      x <- DAS.pval.F.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })

    DAS.pval.F <- do.call(cbind,DAS.pval.F)
    DAS.pval.F <- DAS.pval.F[genes.idx,,drop=F]
    colnames(DAS.pval.F) <- contrast

    DAS.F.stat <- summaryStat (x = DAS.pval.F,y = maxdeltaPS,
                             target = rownames(DAS.pval.F),
                             contrast = contrast,
                             stat.type = c('adj.pval','maxdeltaPS'))

    # DAS.F.stat <- cbind(DAS.pval.F,maxdeltaPS)
    # colnames(DAS.F.stat) <- c(paste0('pval:',contrast),paste0('MaxdeltaPS:',contrast))
    #
    results$DAS.pval.F.list<-DAS.pval.F.list
    results$DAS.pval.F<-DAS.pval.F
    results$DAS.F.stat<-DAS.F.stat

    ##---simes test
    DAS.pval.simes.list<-lapply(contrast,function(i){
      y<-topSplice(fit.splice, coef=i, test="simes", number=Inf, FDR=10000)
      rownames(y)<-y$GENEID
      y[genes.idx,]
    })
    names(DAS.pval.simes.list) <- contrast

    DAS.pval.simes <- lapply(contrast,function(i){
      x <- DAS.pval.simes.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })

    DAS.pval.simes <- do.call(cbind,DAS.pval.simes)
    DAS.pval.simes <- DAS.pval.simes[genes.idx,,drop=F]
    colnames(DAS.pval.simes) <- contrast

    DAS.simes.stat <- summaryStat (x = DAS.pval.simes,y = maxdeltaPS,
                               target = rownames(DAS.pval.simes),
                               contrast = contrast,
                               stat.type = c('adj.pval','maxdeltaPS'))
    #
    # DAS.simes.stat <- cbind(DAS.pval.simes,maxdeltaPS)
    # colnames(DAS.simes.stat) <- c(paste0('pval:',contrast),paste0('MaxdeltaPS:',contrast))

    results$DAS.pval.simes.list<-DAS.pval.simes.list
    results$DAS.pval.simes<-DAS.pval.simes
    results$DAS.simes.stat<-DAS.simes.stat
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat(paste0('Time for analysis: ',round(time.taken,3)))
  cat('\n','> Done!!! ','\n')
  return(results)
}

summaryDEtarget <- function(stat,cutoff=c(adj.pval=0.01,log2FC=1)){
  names(cutoff) <- c('adj.pval','log2FC')
  stat$up.down <- 'up-regulated'
  stat$up.down[stat[,'log2FC']<0] <- 'down-regulated'
  idx <- (abs(stat[,'adj.pval'])<=cutoff['adj.pval']) & (abs(stat[,'log2FC'])>=cutoff['log2FC'])
  stat <- stat[idx,]
  # split(stat , f = stat$contrast)
}

#' Summary DAS genes and DTU transcripts
#' @param stat a data.frame object with first column of "target", second column of "contrast",
#' third and fourth columns of DAS gene/DTU transcript statistics (i.e. "adj.pval" and "maxdeltaPS"/"deltaPS", respectively).
#' @param cutoff a numeric vector of cut-offs for the statistics.
#' @return a data.frame object, which is a subset of input \code{stat} after applying the \code{cutoff}.
#' @export
#'
summaryDAStarget <- function(stat,lfc,cutoff=c(adj.pval=0.01,deltaPS=0.1)){
  names(cutoff) <- c('adj.pval','deltaPS')
  lfc <- lfc[which(lfc$target %in% stat$target),]
  stat <- merge(stat,lfc)
  stat$up.down <- 'up-regulated'
  stat$up.down[stat[,'log2FC']<0] <- 'down-regulated'
  idx <- (abs(stat[,'adj.pval'])<=cutoff['adj.pval']) & (abs(stat[,grep('deltaPS',colnames(stat))])>=cutoff['deltaPS'])
  stat <- stat[idx,]
  # split(stat , f = stat$contrast)
}

#' Compare DE and DAS results
#' @param DE_genes,DAS_genes,DE_trans,DTU_trans data.frame objects of DE genes, DAS genes, DE transcripts and DTU transcripts,
#' which are the outputs of \code{\link{summaryDEtarget}} and \code{\link{summaryDAStarget}}.
#' @param contrast a vector of contrast groups, e.g. \code{contrast = c('B-A','C-A')}, which compares condition B and C to
#' condition A.
#' @details In each contrast group, \code{DEvsDAS} compares the DE and DAS genes; \code{DEvsDTU} compares the DE and DTU
#' transcripts; and \code{summary3Dnumber} summarises the DE/DAS/DTU gene/transcript numbers.
#'
#' @return a data.frame with first column of contrast goups, second column of DE only genes/transcripts, third column of
#' DE&DAS genes or DE&DTU transcripts and fourth column of DAS only genes or DTU only transcripts.
#' @export
#' @rdname DEvsDAS
DEvsDAS <- function(DE_genes,DAS_genes,contrast,consensus=F){
  idx <- c('DEonly','DE&DAS','DASonly')
  x <-  split(DE_genes$target,DE_genes$contrast)
  y <-  split(DAS_genes$target,DAS_genes$contrast)

  num <- lapply(contrast,function(i){
    x <- set2(x[[i]],y[[i]])
    names(x) <- idx
    x
  })
  names(num) <- contrast

  n1 <- lapply(contrast,function(i){
    data.frame(contrast=i,t(sapply(num[[i]],length)))
  })
  n1 <- do.call(rbind,n1)
  colnames(n1) <- c('Contrast',idx)
  if(consensus){
    if(length(contrast)==1){
      n1
    } else {
      num <- unlist(num,recursive=FALSE)
      n2 <- sapply(idx,function(i){
        length(Reduce(intersect,num[grep(paste0('\\.',i,'$'),names(num))]))
      })
      n2 <- data.frame(Contrast='Intersection',t(n2))
      names(n2) <- c('Contrast',idx)
      rbind(n1,n2)
    }
  } else {
    n1
  }
}

#' @export
#' @rdname DEvsDAS
DEvsDTU <- function(DE_trans,DTU_trans,contrast,consensus=F){
  idx <-  c('DEonly','DE&DTU','DTUonly')
  x <-  split(DE_trans$target,DE_trans$contrast)
  y <-  split(DTU_trans$target,DTU_trans$contrast)

  num <- lapply(contrast,function(i){
    x <- set2(x[[i]],y[[i]])
    names(x) <- idx
    x
  })
  names(num) <- contrast

  n1 <- lapply(contrast,function(i){
    data.frame(contrast=i,t(sapply(num[[i]],length)))
  })
  n1 <- do.call(rbind,n1)
  colnames(n1) <- c('Contrast',idx)

  if(consensus){
    if(length(contrast)==1){
      n1
    } else {
      num <- unlist(num,recursive=FALSE)
      n2 <- sapply(idx,function(i){
        length(Reduce(intersect,num[grep(paste0('\\.',i,'$'),names(num))]))
      })
      n2 <- data.frame(Contrast='Intersection',t(n2))
      names(n2) <- c('Contrast',idx)
      rbind(n1,n2)
    }
  } else {
    n1
  }
}

#' @export
#' @rdname DEvsDAS
summary3Dnumber <- function(DE_genes,DAS_genes,DE_trans,DTU_trans,contrast,consensus=F){
  # n1 <- lapply(DE_genes,function(x) x$target)
  idx <- factor(DE_genes$contrast,levels = contrast)
  n1 <- split(DE_genes$target,idx)
  n2 <- length(Reduce(intersect,n1))
  x<- data.frame(`DE genes`=c(sapply(n1,length),`Intersection`=n2))

  # n1 <- lapply(DAS_genes,function(x) x$target)
  idx <- factor(DAS_genes$contrast,levels = contrast)
  n1 <- split(DAS_genes$target,idx)
  n2 <- length(Reduce(intersect,n1))
  x <- cbind(x,data.frame(`DAS genes`=c(sapply(n1,length),`Intersection`=n2)))

  # n1 <- lapply(DE_trans,function(x) x$target)
  idx <- factor(DE_trans$contrast,levels = contrast)
  n1 <- split(DE_trans$target,idx)
  n2 <- length(Reduce(intersect,n1))
  x <- cbind(x,data.frame(`DE transcripts`=c(sapply(n1,length),`Intersection`=n2)))


  # n1 <- lapply(DTU_trans,function(x) x$target)
  idx <- factor(DTU_trans$contrast,levels = contrast)
  n1 <- split(DTU_trans$target,idx)
  n2 <- length(Reduce(intersect,n1))
  x <- cbind(x,data.frame(`DTU transcripts`=c(sapply(n1,length),`Intersection`=n2)))

  x <- data.frame(contrast=rownames(x),x,row.names = NULL)
  colnames(x) <- c('contrast','DE genes','DAS genes','DE transcripts','DTU transcripts')
  if(length(contrast)==1 | consensus==F)
    x <- x[-nrow(x),,drop=F]
  x
}


#' Merge two testing statistics
#' @param x,y data.frame object of statistics (e.g. p-values, log2-fold changes, deltaPS), with rows of targets and
#' columns statistics in contrast groups.
#' @param contrast a vector of contrast groups, e.g. \code{contrast = c('B-A','C-A')}. If \code{NULL}, the column names of
#' \code{x} are used as \code{contrast}.
#' @param stat.type a vector of statistic types of \code{x} and \code{y}, respectively, e.g.
#' \code{stat.type = c('adj.pval','lfc')}. \code{stat.type} is passed to third and fourth column names of the output data.frame.
#' @param srot.by a column name to sort the output data.frame.
#' @return \code{summaryStat} returns a data.frame object with first column of target, second column of contrast groups, third column name of
#' statistics \code{x} and fourth column of statistic \code{y}.
#' @export
summaryStat  <- function(x,y,target,
                         contrast = NULL,
                         stat.type = c('adj.pval','lfc'),
                         srot.by = stat.type[1]){
  x <- x[target,,drop=F]
  y <- y[target,,drop=F]
  if(is.null(contrast))
    contrast <- colnames(x)
  x <- x[,contrast,drop=F]
  y <- y[,contrast,drop=F]

  stat <- lapply(contrast,function(i){
    z <- data.frame(targets =target,contrast=i, x[,i],y[,i],row.names = NULL)
    colnames(z) <- c('target','contrast',stat.type)
    z <- z[order(z[,srot.by]),]
  })
  names(stat) <- contrast
  stat <- do.call(rbind,stat)
  rownames(stat) <- NULL
  stat
}

transAbundance2PS <- function(transAbundance=NULL,
                              PS=NULL,
                              contrast,
                              condition,
                              mapping){
  colnames(mapping) <- c('TXNAME','GENEID')
  if(!is.null(transAbundance))
    rownames(transAbundance) <- mapping$TXNAME

  if(is.null(PS)){
    transAbundance.mean <- t(rowmean(t(transAbundance),group = condition,reorder = F))
    ##---PS
    PS <- rowratio(x = transAbundance.mean[as.vector(mapping$TXNAME),],
                   group = as.vector(mapping$GENEID), reorder = F)
    PS <- PS[mapping$TXNAME,]
  }

  contrast.matrix <- makeContrasts(contrasts = contrast,levels = factor(condition,levels = unique(condition)))
  deltaPS <- data.frame(as.matrix(PS)%*%as.matrix(contrast.matrix),check.names = F)
  ##---deltaPSI
  # deltaPS <- lapply(strsplit(contrast,'-'),function(x){
  #   PS[,x[1]]-PS[,x[2]]
  # })
  # deltaPS <- do.call(cbind,deltaPS)
  colnames(deltaPS) <- contrast
  rownames(PS) <- rownames(deltaPS) <- mapping$TXNAME
  list(PS=PS,deltaPS=deltaPS)
}


rowratio <- function (x, group, reorder = F, na.rm = T){
  y <- rowsum(x, group = group)
  y <- y[group, ]
  ratio = x/y
  rownames(ratio) <- rownames(x)
  if (na.rm)
    ratio[is.na(ratio)] <- 0
  if (reorder & !is.null(rownames(ratio)))
    ratio <- ratio[gtools::mixedsort(rownames(ratio)), ]
  return(ratio)
}


set2 <- function(x,y){
  x.only <- setdiff(x,y)
  xy <- intersect(x,y)
  y.only <- setdiff(y,x)
  results <- list(x.only=x.only,xy=xy,y.only=y.only)
  attributes(results) <- list(
    x.only='x-x&y',
    xy='x&y',
    y.only='y-x&y'
  )
  names(results) <- c('x.only','xy','y.only')
  return(results)
}

#############################################################################################
##-----------------------------------------------------------------------------------------
#' @rdname set2
#' @export

set3 <- function(x,y,z){
  a5 <- Reduce(intersect, list(x,y,z))
  a2 <- setdiff(intersect(x,y),a5)
  a4 <- setdiff(intersect(x,z),a5)
  a6 <- setdiff(intersect(y,z),a5)
  a1 <- setdiff(x,c(a2,a4,a5))
  a3 <- setdiff(y,c(a2,a5,a6))
  a7 <- setdiff(z,c(a4,a5,a6))
  results <- list(a1=a1,a2=a2,a3=a3,a4=a4,a5=a5,a6=a6,a7=a7)
  # attributes(results) <- list(
  #   a1='x-(y&z)',
  #   a2='(x&y)-z',
  #   a3='y-(x&z)',
  #   a4='(x&z)-y',
  #   a5='x&y&z',
  #   a6='(y&z)-x',
  #   a7='z-(x&y)'
  # )
  name.idx <- c('x-(y&z)','(x&y)-z','y-(x&z)','(x&z)-y','x&y&z','(y&z)-x','z-(x&y)')
  name.idx <- paste0(paste0('a',1:7,':'),name.idx)
  names(results) <- name.idx
  return(results)
}

gg.color.hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

trend.mean.variance <- function(obj, design, lib.size = NULL, span = 0.5, ...){
  if (is(obj, "DGEList")) {
    counts <- obj$counts
    if (is.null(lib.size))
      lib.size <- with(obj$samples, lib.size * norm.factors)
  } else {
    counts <- as.matrix(obj)
    if(is.null(lib.size))
      lib.size <- colSums(counts)
  }
  y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
  fit <- lmFit(y, design, ...)
  if (is.null(fit$Amean))
    fit$Amean <- rowMeans(y, na.rm = TRUE)
  sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)
  sy <- sqrt(fit$sigma)
  allzero <- rowSums(counts) == 0
  if (any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  l <- lowess(sx, sy, f = span)
  list(fit=fit,sx=sx,sy=sy,l=l)
}


#' Compare expression mean-variance trends before and after low expression filters.
#' @param counts.raw a matrix/data.frame of raw read counts before low expression filters.
#' @param counts.filtered a matrix/data.frame of read counts after low expression filters.
#' @param condition a vector of characters to distinguish conditions of samples (e.g. c('A','A','B','B')), which is used to make the design
#' matrix to fit linear regression model.
#' @param span a numeric value passed to \code{\link{lowess}} smoothing window as a proportion to fit the mean-variance trend.
#' @param ... additional arguments passed to \code{\link{trend.mean.variance}} function.
#' @return a list object with elements "fit.raw" and "fit.flitered", which are the \code{\link{trend.mean.variance}} fit results
#' for \code{counts.raw} and \code{counts.filtered}, respectively.
#' @export
#' @seealso \code{\link{trend.mean.variance}} and \code{\link[limma]{voom}}.
check.mean.variance <- function(counts.raw,
                                counts.filtered,
                                condition,
                                span = 0.5,...){
  ###---design matrix
  condition <- factor(condition,levels = unique(condition))
  design<-model.matrix(~0+condition)
  colnames(design) <- gsub('condition','',colnames(design))

  ###---generate DGEList object
  message('=> Generate DGEList object')
  dge.raw <- DGEList(counts=counts.raw)
  dge.raw <- calcNormFactors(dge.raw)
  dge.filtered<- DGEList(counts=counts.filtered)
  dge.filtered <- calcNormFactors(dge.filtered)

  ###---fit mean-variance trend
  message('=> Fit mean-variance trend')
  #before filter
  fit.raw <- trend.mean.variance(obj = dge.raw,design = design,span = span,...)
  #after filter
  fit.filtered <- trend.mean.variance(obj = dge.filtered,design = design,span = span,...)

  message('Done!!!')
  return(list(fit.raw=fit.raw,fit.filtered=fit.filtered))
}

#' Plot the mean-variance trend
#' @param x a numeric vector of mean value from the results of \code{\link{trend.mean.variance}}.
#' @param y a numeric vector of variance value from the results of \code{\link{trend.mean.variance}}.
#' @param x.lab,y.lab a string of x-axis label and y-axis label, respectively.
#' @param main plot title.
#' @param ... additional arguments passed to \code{\link{plot}}.
#' @return a plot
#' @export
plotMeanVariance <- function(x,y,l,fit.line.col='red',
                               x.lab = "log2( count size + 0.5 )",
                               y.lab = "Sqrt(standard deviation)",
                               main="Mean-variance trend",lwd=1.5,...){
  plot(x, y, pch = 16, cex = 0.25,xlab=x.lab,ylab=y.lab,main=main,...)
  lines(l,col=fit.line.col,lwd=lwd)
}

edgeR.pipeline <- function(dge,
                           method=c('glm','glmQL'),
                           design,
                           contrast,
                           diffAS=F,
                           deltaPS=NULL,
                           adjust.method='BH'){

  start.time <- Sys.time()
  if(is.null(dge$genes)){
    diffAS <- F
  } else {
    if(max(table(dge$genes$GENEID)) < 2)
      diffAS <- F
  }

  results <- list()
  method <- match.arg(method,c('glm','glmQL'))
  targets <- rownames(dge$counts)
  cat('> Estimate dispersion ...','\n')
  Disp <- estimateGLMCommonDisp(dge, design)
  Disp <- estimateGLMTagwiseDisp(Disp, design)
  contrast.matrix <- makeContrasts(contrasts = contrast, levels=design)

  switch(method,
         glm = {
           cat('> Fit Genewise Negative Binomial Generalized Linear Models...','\n')
           fit <- glmFit(Disp, design = design)
           results$fit.glm <- fit

           ###individual test
           cat('> Fit the contrast model ...','\n')
           DE.pval.list <- lapply(contrast,function(i){
             x <- glmLRT(fit, contrast=contrast.matrix[,i,drop=F])
             y <- topTags(x,n = Inf,adjust.method = adjust.method)
             y <- y[targets,]
           })
           ###overall test
           x <- glmLRT(fit, contrast=contrast.matrix)
           DE.stat.overalltest <- topTags(x,n = Inf,adjust.method = adjust.method)
           # DE.stat.overalltest <- DE.stat.overalltest[targets,]
           DE.stat.overalltest <- data.frame(target=rownames(DE.stat.overalltest),
                                             contrast='overall',DE.stat.overalltest,row.names = NULL)
           results$DE.stat.overalltest<-DE.stat.overalltest
         },
         glmQL = {
           cat('> Genewise Negative Binomial Generalized Linear Models with Quasi-likelihood Tests...','\n')
           fit <- glmQLFit(Disp,design = design)
           results$fit.glmQL <- fit

           ###individual test
           cat('> Fit the contrast model ...','\n')
           DE.pval.list <- lapply(contrast,function(i){
             x <- glmQLFTest(fit, contrast=contrast.matrix[,i,drop=F])
             y <- topTags(x,n = Inf,adjust.method = adjust.method)
             y <- y[targets,]
           })

           ###overall test
           x <- glmQLFTest(fit, contrast=contrast.matrix)
           DE.stat.overalltest <- topTags(x,n = Inf,adjust.method = adjust.method)
           # DE.stat.overalltest <- DE.stat.overalltest[targets,]
           DE.stat.overalltest <- data.frame(target=rownames(DE.stat.overalltest),
                                             contrast='overall',DE.stat.overalltest,row.names = NULL)
           results$DE.stat.overalltest<-DE.stat.overalltest
         }
  )
  names(DE.pval.list) <- contrast
  results$DE.pval.list <- DE.pval.list

  DE.pval <- do.call(cbind,lapply(DE.pval.list,FUN = function(x) x$table$FDR))
  DE.lfc <- do.call(cbind,lapply(DE.pval.list,FUN = function(x) x$table$logFC))
  DE.rawpval <- do.call(cbind,lapply(DE.pval.list,FUN = function(x) x$table$PValue))#Added by Max Coulter 17/05/21
  rownames(DE.pval) <- rownames(DE.lfc) <- rownames(DE.rawpval) <- targets
  colnames(DE.pval) <- colnames(DE.lfc) <- colnames(DE.rawpval) <- contrast

  DE.stat <- summaryStat (x = DE.pval,y = DE.lfc,
                          target = rownames(DE.pval),
                          contrast = contrast,
                          stat.type = c('adj.pval','log2FC'))
  results$DE.pval<-DE.pval
  results$DE.lfc<-DE.lfc
  results$DE.stat<-DE.stat
  results$DE.rawpval<-DE.rawpval#Added by Max Coulter 17/05/21

  ##################################################################

  ###---DAS
  if(diffAS){
    fit.splice <- lapply(contrast,function(i){
      cat(paste0('\ndiffSpliceDGE of contrast: ', i))
      diffSpliceDGE(fit, contrast=contrast.matrix[,i,drop=F], geneid="GENEID")
    })
    names(fit.splice) <- contrast
    results$fit.splice<-fit.splice

    genes.idx <- unique(fit.splice[[1]]$genes$GENEID)
    trans.idx <- unique(fit.splice[[1]]$genes$TXNAME)

    DTU.pval.list<-lapply(contrast,function(i){
      fit.splice.i <- fit.splice[[i]]
      y<-topSpliceDGE(fit.splice.i, test="exon", number=Inf, FDR=10000)
      rownames(y)<-y$TXNAME
      z <- y[trans.idx,]
      z
    })
    names(DTU.pval.list) <- contrast

    ##DTU transcript pvals
    DTU.pval <- lapply(contrast,function(i){
      x <- DTU.pval.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })
    DTU.pval <- do.call(cbind,DTU.pval)
    colnames(DTU.pval) <- contrast

    DTU.deltaPS <- deltaPS[rownames(DTU.pval),,drop=F]
    DTU.stat <- summaryStat (x = DTU.pval,y = DTU.deltaPS,
                             target = rownames(DTU.pval),
                             contrast = contrast,
                             stat.type = c('adj.pval','deltaPS'))

    results$DTU.pval.list<-DTU.pval.list
    results$DTU.pval<-DTU.pval
    results$DTU.deltaPS<-DTU.deltaPS
    results$DTU.stat<-DTU.stat

    # ##########################################################
    ##---DAS genes
    ##---max deltaPS
    maxdeltaPS <- by(deltaPS[fit.splice[[1]]$genes$TXNAME,,drop=F],
                     INDICES = fit.splice[[1]]$genes$GENEID,
                     function(x){
                       apply(x,2,function(i) i[abs(i)==max(abs(i))][1])
                     },simplify = F)
    maxdeltaPS <- do.call(rbind,maxdeltaPS)
    maxdeltaPS <- maxdeltaPS[genes.idx,,drop=F]
    results$maxdeltaPS<-maxdeltaPS

    ##---F test
    DAS.pval.F.list<-lapply(contrast,function(i){
      fit.splice.i <- fit.splice[[i]]
      y<-topSpliceDGE(fit.splice.i, test="gene", number=Inf, FDR=10000)
      rownames(y)<-y$GENEID
      z <- y[genes.idx,]
      z
    })
    names(DAS.pval.F.list) <- contrast

    DAS.pval.F <- lapply(contrast,function(i){
      x <- DAS.pval.F.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })

    DAS.pval.F <- do.call(cbind,DAS.pval.F)
    DAS.pval.F <- DAS.pval.F[genes.idx,,drop=F]
    colnames(DAS.pval.F) <- contrast

    DAS.F.stat <- summaryStat (x = DAS.pval.F,y = maxdeltaPS,
                               target = rownames(DAS.pval.F),
                               contrast = contrast,
                               stat.type = c('adj.pval','maxdeltaPS'))

    results$DAS.pval.F.list<-DAS.pval.F.list
    results$DAS.pval.F<-DAS.pval.F
    results$DAS.F.stat<-DAS.F.stat

    ##---Simes test
    DAS.pval.simes.list<-lapply(contrast,function(i){
      fit.splice.i <- fit.splice[[i]]
      y<-topSpliceDGE(fit.splice.i, test="Simes", number=Inf, FDR=10000)
      rownames(y)<-y$GENEID
      z <- y[genes.idx,]
      z
    })
    names(DAS.pval.simes.list) <- contrast

    DAS.pval.simes <- lapply(contrast,function(i){
      x <- DAS.pval.simes.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })

    DAS.pval.simes <- do.call(cbind,DAS.pval.simes)
    DAS.pval.simes <- DAS.pval.simes[genes.idx,,drop=F]
    colnames(DAS.pval.simes) <- contrast

    DAS.simes.stat <- summaryStat (x = DAS.pval.simes,y = maxdeltaPS,
                                   target = rownames(DAS.pval.simes),
                                   contrast = contrast,
                                   stat.type = c('adj.pval','maxdeltaPS'))

    results$DAS.pval.simes.list<-DAS.pval.simes.list
    results$DAS.pval.simes<-DAS.pval.simes
    results$DAS.simes.stat<-DAS.simes.stat
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat(paste0('Time for analysis: ',round(time.taken,3)))
  cat('\n','> Done!!! ','\n')
  return(results)
}

sumarrays<-function(x,group=NULL){
  if(is.null(group))
    group<-colnames(x)
  colnames(x)<-group

  x<-rowsum(t(x),group = group)
  #order the columns as the input group odering
  x<-data.frame(t(x)[,unique(group)])
  colnames(x) <- unique(group)
  return(x)
}


low.expression.filter <- function(abundance,
                                  mapping,
                                  abundance.cut=1,
                                  sample.n=3,
                                  Log=F,
                                  unit=c('counts','TPM')){
  unit <- match.arg(unit,c('counts','TPM'))
  colnames(mapping) <- c('TXNAME','GENEID')
  rownames(mapping) <- mapping$TXNAME
  if(unit == 'counts'){
    filter.idx<-CPM.filter(counts = abundance,sample.n = sample.n,cpm.cut = abundance.cut,Log = Log)
  } else {
    filter.idx<-TPM.filter(TPM = abundance,sample.n = sample.n,tpm.cut = abundance.cut)
  }
  trans_high <- names(filter.idx)[filter.idx==T]
  genes_high <- unique(mapping[trans_high,]$GENEID)
  list(trans_high=trans_high,genes_high=genes_high,mapping_high=mapping[trans_high,])
}

CPM.filter<-function(counts,sample.n=3,cpm.cut=1,...){
  rowSums(counts2CPM(counts,...)>=cpm.cut)>=sample.n
}

counts2CPM <- function (obj,lib.size = NULL,Log=F){
  if (is(obj, "DGEList")) {
    counts <- obj$counts
    if (is.null(lib.size))
      lib.size <- with(obj$samples, lib.size * norm.factors)
  } else {
    counts <- as.matrix(obj)
    if(is.null(lib.size))
      lib.size <- colSums(counts)
  }
  if(Log){
    t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
  } else {
    t(t(counts)/(lib.size) * 1e+06)
  }
}
