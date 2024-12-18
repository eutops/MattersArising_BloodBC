plot_single_roc <- function(type1, index1, col1 = "black", title1 = "all", style = "default",
                            direction1 = "<", generalTitle="", labelSize=2.9){
  
  require(pROC)
  require(ggsci)
  cols <- pal_lancet(palette="lanonc", alpha = 0.8)(8) 
  
  #--------------------------------#  
  # make annotation labels
  auc <- round(as.numeric(roc(type1,index1, quiet=T, ci=T, direction = direction1)$auc),digits=2)
  cil <- round(as.numeric(roc(type1,index1, quiet=T, ci=T, direction = direction1)$ci[1]),digits=2)
  ciu <- round(as.numeric(roc(type1,index1, quiet=T, ci=T, direction = direction1)$ci[3]),digits=2)
  anno1 <- paste('AUC (', title1, ') = ',auc,'\n(95% CI: ', cil,'-',ciu,')',sep='')
  
  
  
  if(style == "lancet"){
    anno1 <- gsub("[.]", "·", anno1)
    anno2 <- gsub("[.]", "·", anno2)
  }
  
  #--------------------------------#  
  
  roc1 <- roc(type1, index1, direction = direction1)
  title1 <- as.character(title1)
  
  
  ggplot() +
    geom_path(aes(x=1-roc1$specificities,
                  y=(roc1$sensitivities)),
              colour = col1,
              size = 0.7) +
    annotate("segment", x = 0, y = 0,
             xend = 1, yend = 1,
             colour = "gray60") +
    xlab("1 - Specificity") +
    ylab("Sensitivity") +
    theme_minimal() +
    theme(plot.title = element_text(size=10),
          panel.grid = element_blank(),
          axis.line = element_line()) +
    annotate(geom='label',
             x=0.65,
             y=0.1,
             label=anno1,
             fill='black',
             colour='white',
             size=labelSize) +
    ggtitle(generalTitle)  
  
}
