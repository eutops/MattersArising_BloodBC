plot_double_roc <- function(type1, index1, type2, index2, col1 = "black", col2 = "blue", title1 = "all", title2 = "below 30", style = "default",
                            direction1 = "<",direction2 = "<", generalTitle="", labelSize=2.9, textcol=F){
  
  require(pROC)
  require(ggsci)
  cols <- pal_lancet(palette="lanonc", alpha = 0.8)(8) 
  
  #--------------------------------#  
  # make annotation labels
  auc <- round(as.numeric(roc(type1,index1, quiet=T, ci=T, direction = direction1)$auc),digits=2)
  cil <- round(as.numeric(roc(type1,index1, quiet=T, ci=T, direction = direction1)$ci[1]),digits=2)
  ciu <- round(as.numeric(roc(type1,index1, quiet=T, ci=T, direction = direction1)$ci[3]),digits=2)
  anno1 <- paste('AUC (', title1, ') = ',auc,'\n(95% CI: ', cil,'-',ciu,')',sep='')
  
  # make annotation labels
  auc2 <- round(as.numeric(roc(type2,index2, quiet=T, ci=T, direction = direction2)$auc),digits=2)
  cil2 <- round(as.numeric(roc(type2,index2, quiet=T, ci=T, direction = direction2)$ci[1]),digits=2)
  ciu2 <- round(as.numeric(roc(type2,index2, quiet=T, ci=T, direction = direction2)$ci[3]),digits=2)
  anno2 <- paste('AUC (', title2, ') = ',auc2,'\n(95% CI: ', cil2,'-',ciu2,')',sep='')
  
  if(style == "lancet"){
    anno1 <- gsub("[.]", "·", anno1)
    anno2 <- gsub("[.]", "·", anno2)
  }
  
  #--------------------------------#  
  
  roc1 <- roc(type1, index1, direction = direction1)
  roc2 <- roc(type2, index2, direction = direction2)
  title1 <- as.character(title1)
  title2 <- as.character(title2)
  
  if (class(textcol)!="character") {
  # Determine text colour of annotation labels based on input hex code of col1 and col2
  get_contrast <- function(hex_code) {
    # Calculate luminance from hex code
    luminance_from_hex <- function(hex) {
      hex <- gsub("#", "", hex)
      if (nchar(hex) != 6) {
        stop("Invalid hex code length.")
      }
      r <- strtoi(substr(hex, 1, 2), 16L)
      g <- strtoi(substr(hex, 3, 4), 16L)
      b <- strtoi(substr(hex, 5, 6), 16L)
      lum <- 0.2126 * r + 0.7152 * g + 0.0722 * b
      return(lum / 255)  # Normalize to range [0, 1]
    }
    
    # Get luminance from hex code
    lum <- tryCatch(
      {
        luminance_from_hex(hex_code)
      },
      error = function(e) {
        warning("Invalid hex code.")
        return(NA)
      }
    )
    
    # Assign black or white based on contrast
    if (!is.na(lum)) {
      if (lum > 0.4) {
        return("black")  # if background is light, text should be black
      } else {
        return("white")  # if background is dark, text should be white
      }
    } else {
      return(NULL)
    }
  }
  
  textcol1 <- get_contrast(col1)
  textcol2 <- get_contrast(col2)
  
  }else{
  
    textcol1 <- textcol[1]
    textcol2 <- textcol[2]
    
    
}
  
  ggplot() +
    geom_path(aes(x=1-roc1$specificities,
                  y=(roc1$sensitivities)),
              colour = col1,
              size = 0.7) +
    geom_path(aes(x=1-roc2$specificities,
                  y=(roc2$sensitivities)),
              colour = col2,
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
             x=0.55,
             y=0.4,
             label=anno1,
             fill=col1,
             colour=textcol1,
             size=labelSize)  +
    annotate(geom='label',
             x=0.55,
             y=0.1,
             label=anno2,
             fill=col2,
             colour=textcol2,
             size=labelSize)  +
    ggtitle(generalTitle)  
  
}
