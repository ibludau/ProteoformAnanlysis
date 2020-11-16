plot.tracesList <- function(traces,
                            design_matrix = NULL,
                            collapse_conditions = FALSE,
                            aggregateReplicates=FALSE,
                            log=FALSE,
                            legend = TRUE,
                            PDF=FALSE,
                            name="Traces",
                            plot = TRUE,
                            isoformAnnotation = FALSE,
                            colour_by = "id",
                            error_bar=TRUE,
                            highlight=NULL,
                            highlight_col=NULL,
                            colorMap=NULL,
                            monomer_MW=TRUE) {
  #.tracesListTest(traces)
  if(!is.null(design_matrix)){
    if(!all(design_matrix$Sample_name %in% names(traces))){
      stop("Invalid design matrix")
    }
  }else{
    design_matrix <- data.table(Sample_name = names(traces),
                                Condition = "",
                                Replicate = 1:length(traces))
  }
  tracesList <- lapply(names(traces), function(tr){
    res <- toLongFormat(traces[[tr]]$traces)
    res$Condition <- design_matrix[Sample_name == tr, Condition]
    res$Replicate <- design_matrix[Sample_name == tr, Replicate]
    res
  })
  traces_long <- do.call("rbind", tracesList)
  if(colour_by!="id") {
    if(!colour_by %in% names(traces[[1]]$trace_annotation)){
      stop("colour_by is not availbale in trace_annotation.")
    }
    isoform_annotation <- lapply(names(traces), function(tr){subset(traces[[tr]]$trace_annotation,select=c("id",colour_by))})
    isoform_annotation <- unique(do.call("rbind", isoform_annotation))
    traces_long <- merge(traces_long,isoform_annotation, by.x="id",by.y="id")
    traces_long[,line:=paste0(get(colour_by),id)]
  }
  ## Create a common fraction annotation
  traces_frac <- unique(do.call("rbind", lapply(traces, "[[", "fraction_annotation")))
  traces_frac <- unique(subset(traces_frac, select = names(traces_frac) %in% c("id","molecular_weight")))
  traces_long <- merge(traces_long,traces_frac,by.x="fraction",by.y="id")
  
  if(!is.null(highlight)){
    traces_long$outlier <- gsub("\\(.*?\\)","",traces_long$id) %in% gsub("\\(.*?\\)","",highlight)
    if(!any(traces_long$outlier)) highlight <- NULL
  }
  
  geom.text.size = 3
  theme.size = (14/5) * geom.text.size
  
  ## Create a reproducible coloring for the peptides plotted
  if(!is.null(colorMap)){
    if(!all(unique(traces_long[[colour_by]]) %in% names(colorMap))){
      stop("Invalid colorMap specified. Not all traces to be plotted are contained in the colorMap")
    }
  }else{
    cbPalette <- c("#56B4E9","#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999")
    ids <- sort(unique(traces_long[[colour_by]]))
    if (length(ids) <= length(cbPalette)) {
      colorMap <- cbPalette[1:length(unique(traces_long[[colour_by]]))]
      names(colorMap) <- ids
    } else {
      colorMap <- createGGplotColMap(ids)
    }
  }
  if (aggregateReplicates){
    traces_long[,meanIntensity := mean(intensity), by=c("Condition","fraction","id")]
    traces_long[,sdIntensity := sd(intensity), by=c("Condition","fraction","id")]
    if(colour_by!="id") {
      traces_long <- unique(subset(traces_long,select=c("id",colour_by, "line","fraction","Condition","molecular_weight","meanIntensity","sdIntensity")))
    } else {
      traces_long <- unique(subset(traces_long,select=c("id","fraction","Condition","molecular_weight","meanIntensity","sdIntensity")))
    }
    traces_long[,Replicate := "average"]
    setnames(traces_long, "meanIntensity", "intensity")
  }

  p <- ggplot(traces_long) +
    xlab('fraction') +
    ylab('intensity') +
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    theme(plot.margin = unit(c(1,.5,.5,.5),"cm")) +
    ggtitle(name) +
    scale_color_manual(values=colorMap) #+
  #theme(plot.title = element_text(vjust=19,size=10))
  
  if(collapse_conditions){
    if(colour_by == "id") {
      p <- p + facet_grid(~ Replicate) +
        geom_line(aes_string(x='fraction', y='intensity', color='id', lty = 'Condition'))
      if (aggregateReplicates){
        if(error_bar){
          p <- p + geom_errorbar(aes(x=fraction,ymin=ifelse(intensity-sdIntensity < 0, 0, intensity-sdIntensity), ymax=intensity+sdIntensity, color=id, lty = Condition), width=0.2, size=0.3, position=position_dodge(0.05))
        }
      }
    } else {
      message("Collapsing of conditions is not jet compatible with colouring by your selected id type.
      Plot conditions separately instead.")
      p <- p + facet_grid(Condition ~ Replicate) +
        geom_line(aes_string(x='fraction', y='intensity', color=colour_by, group='line'))
      if (aggregateReplicates){
        if(error_bar){
          p <- p + geom_errorbar(aes(x=fraction,ymin=ifelse(intensity-sdIntensity < 0, 0, intensity-sdIntensity), ymax=intensity+sdIntensity, color=colour_by), width=0.2, size=0.3, position=position_dodge(0.05))
        }
      }
    }
  }else{
    if (length(unique(traces_long$Replicate)) > 1) {
      if(colour_by == "id") {
        p <- p + facet_grid(Condition ~ Replicate) +
          geom_line(aes_string(x='fraction', y='intensity', color='id'))
        if (aggregateReplicates){
          if(error_bar){
            p <- p + geom_errorbar(aes(x=fraction,ymin=ifelse(intensity-sdIntensity < 0, 0, intensity-sdIntensity), ymax=intensity+sdIntensity, color=id), width=0.2, size=0.3, position=position_dodge(0.05))
          }
        }
      } else {
        p <- p + facet_grid(Condition ~ Replicate) +
          geom_line(aes_string(x='fraction', y='intensity', color=colour_by, group='line'))
      }
      if (aggregateReplicates){
        if(error_bar){
          p <- p + geom_errorbar(aes(x=fraction,ymin=ifelse(intensity-sdIntensity < 0, 0, intensity-sdIntensity), ymax=intensity+sdIntensity, color=colour_by), width=0.2, size=0.3, position=position_dodge(0.05))
        }
      }
    } else {
      if(colour_by == "id") {
        p <- p + facet_grid(Condition ~ .) +
          geom_line(aes_string(x='fraction', y='intensity', color='id'))
        if (aggregateReplicates){
          if(error_bar){
            p <- p + geom_errorbar(aes(x=fraction,ymin=ifelse(intensity-sdIntensity < 0, 0, intensity-sdIntensity), ymax=intensity+sdIntensity, color=id), width=0.2, size=0.3, position=position_dodge(0.05))
          }
        }
      } else {
        p <- p + facet_grid(Condition ~ .) +
          geom_line(aes_string(x='fraction', y='intensity', color=colour_by, group='line'))
        if (aggregateReplicates){
          if(error_bar){
            p <- p + geom_errorbar(aes(x=fraction,ymin=ifelse(intensity-sdIntensity < 0, 0, intensity-sdIntensity), ymax=intensity+sdIntensity, color=colour_by), width=0.2, size=0.3, position=position_dodge(0.05))
          }
        }
      }
    }
  }
  
  if(!is.null(highlight)){
    legend_peps <- unique(traces_long[outlier == TRUE, id])
    if(is.null(highlight_col)){
      if(collapse_conditions){
        p <- p +
          geom_line(data = traces_long[outlier == TRUE], aes_string(x='fraction', y='intensity', color='id', lty = 'Condition'), lwd=2) +
          scale_color_manual(values = colorMap, breaks = legend_peps)
      }else{
        p <- p +
          geom_line(data = traces_long[outlier == TRUE], aes_string(x='fraction', y='intensity', color='id'), lwd=2) +
          scale_color_manual(values = colorMap, breaks = legend_peps)
      }
      
    }else{
      ## legend_map <- unique(ggplot_build(p)$data[[1]]$colour)
      ## names(legend_map) <- unique(p$data$id)
      ## legend_map[legend_peps] <- highlight_col
      ## legend_vals <- rep(highlight_col, ceiling(length(legend_peps)/ length(highlight_col)))[1:length(legend_peps)]
      if(collapse_conditions){
        p <- p +
          geom_line(data = traces_long[outlier == TRUE],
                    aes(x=fraction, y=intensity, lty = Condition, group = interaction(Condition, id), color = id),
                    lwd=2) +
          # scale_color_discrete(guide = F)
          scale_color_manual(values = colorMap, breaks = legend_peps)
        # guides(lty = FALSE)
        # scale_color_manual(limits = legend_peps, values = rep(highlight_col, length(legend_peps))) +
        # geom_line(aes_string(x='fraction', y='intensity', color='id'))
      }else{
        p <- p +
          geom_line(data = traces_long[outlier == TRUE], aes_string(x='fraction', y='intensity', color = 'id'),
                    lwd=2) +
          # scale_color_discrete(guide = F)
          scale_color_manual(values = colorMap, breaks = legend_peps)
        ## scale_color_manual(values = legend_map, limits = legend_peps)
        # guides(lty = FALSE)
        # scale_color_manual(limits = legend_peps, values = rep(highlight_col, length(legend_peps))) +
        # geom_line(aes_string(x='fraction', y='intensity', color='id'))
      }
      
    }
  }
  
  if ("molecular_weight" %in% names(traces_frac)) {
    fraction_ann <- traces_frac
    tr <- lm(log(fraction_ann$molecular_weight) ~ fraction_ann$id)
    intercept <- as.numeric(tr$coefficients[1])
    slope <- as.numeric(tr$coefficients[2])
    mwtransform <- function(x){exp(slope*x + intercept)}
    MWtoFraction <- function(x){round((log(x)-intercept)/(slope), digits = 0)}
    mw <- round(fraction_ann$molecular_weight, digits = 0)
    
    breaks_frac <- seq(min(fraction_ann$id),max(fraction_ann$id),10)
    breaks_MW <- round(fraction_ann[id %in% seq(min(fraction_ann$id),max(fraction_ann$id),10)]$molecular_weight)
    p <- p + scale_x_continuous(name="fraction",
                                breaks=breaks_frac,
                                labels=breaks_frac,
                                sec.axis = dup_axis(trans = ~.,
                                                    breaks=breaks_frac,
                                                    labels = breaks_MW,
                                                    name = "MW (kDa)"))
    
    if (monomer_MW==TRUE){
      if ("protein_mw" %in% names(traces[[1]]$trace_annotation)) {
        ann_tab <- lapply(traces, function(t){subset(t$trace_annotation, select=c(colour_by,"protein_mw"))})
        ann_tab <- unique(do.call(rbind,ann_tab))
        subunitMW.dt <- data.table(id=ann_tab[[colour_by]],mw=ann_tab$protein_mw)
        subunitMW.dt$fraction <- MWtoFraction(subunitMW.dt$mw)
        subunitMW.dt[,boundary:=MWtoFraction(2*mw)]
        if (length(unique(subunitMW.dt$mw)) > 1) {
          p <- p + geom_point(data = subunitMW.dt, mapping = aes(x = fraction, y = Inf, colour=id), shape=18,size=5,alpha=.5)
        } else {
          p <- p + geom_vline(data = unique(subunitMW.dt), aes(xintercept = fraction), colour="red", linetype="dashed", size=.5)
          p <- p + geom_vline(data = unique(subunitMW.dt), aes(xintercept = boundary), colour="red", linetype="dashed", size=.5, alpha=0.5)
        }
      } else {
        message("No molecular weight annotation of the traces. Cannot plot monomer molecular weight.")
      }
    }
  } else {
    p <- p + scale_x_continuous(name="fraction",
                                breaks=seq(min(traces$fraction_annotation$id),
                                           max(traces$fraction_annotation$id),10),
                                labels=seq(min(traces$fraction_annotation$id),
                                           max(traces$fraction_annotation$id),10))
  }
  
  if (log) {
    p <- p + scale_y_log10('log(intensity)')
  }
  
  p <- p + theme(axis.text = element_text(size = theme.size, colour="black"))
  
  p <- p + theme(legend.position="bottom") +
    theme(legend.text=element_text(size=theme.size)) +
    theme(legend.title=element_blank()) +
    guides(col = guide_legend(ncol = 4))
  
  if (!legend) {
    p <- p + theme(legend.position="none")
  }
  
  if (length(ids) > 40) {
    p <- p + theme(legend.position="none")
  }
  
  if(PDF){
    pdf(paste0(name,".pdf"),width=5,height=4)
  }
  if(plot){
    plot(p)
  }else{
    return(p)
  }
  if(PDF){
    dev.off()
  }
}
