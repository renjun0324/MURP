
#' MURP_plot
#'
#' @description
#' Scatter plot of the relationship between k and pseudo-BIC
#'
#' @param murpResult: MURP Result
#'
#' @import ggplot2
#'
#' @export
#'

KBicPlot <- function(murpResult = NULL){

  murp = murpResult
  df = data.frame(k = murp$k,
                  bic = murp$BIC,
                  stringsAsFactors = FALSE)
  df$recommended_k = murp$Recommended_K

  p = ggplot(df,aes(x = k, y = bic)) +
    geom_point(size = 3, alpha = 0.6, colour = "gray20") +
    geom_vline(xintercept = df$recommended_k, color = "darkgreen",
               alpha = 0.8, size = 0.6, linetype = "dashed") +
    labs(x = "Number of MURPs / Cell Number", y = "pseudo-BIC") +
    theme(panel.background = element_rect(fill='transparent', color="black"),
          strip.text = element_text(size = 10),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          panel.border = element_rect(fill='transparent', color='black'),
          plot.title = element_text(size = 12.5, hjust = 0), # title
          plot.subtitle = element_text(size = 11.5, hjust = 0), # subtitle
          legend.key = element_rect( fill = "white"),
          axis.title.x = element_text(vjust = -1.5, size = 15, colour = 'black'), # face = "bold"
          axis.title.y = element_text(vjust = 1.5, size = 15, colour = 'black'), # face = "bold"
          #axis.line = element_line(colour = 'black'),
          axis.ticks = element_blank(),
          #axis.text.y = element_blank(),
          axis.text.x = element_text(vjust = -0.5, size = 12, colour = 'black'),
          axis.text.y = element_text(vjust = 0.5, size = 12, colour = 'black'),
          legend.text = element_text(vjust = 0.4, size = 15, colour = 'black'),
          legend.title = element_text(vjust = 0.4, size = 15, colour = 'black'),
          legend.key.size = unit(0.9, "cm") )

  return(p)
}

#' MURPNestedGridPlot
#'
#' @Description:
#' Partial zoom of KBicPlot
#'
#' @import dplyr
#' @import ggplot2
#'
#' @param murp form MURP
#'
#' @export
MURPNestedGridPlot <- function(murpResult = NULL){

  #require(dplyr)
  murp = murpResult
  rj.ftheme <-   theme(panel.background = element_rect(fill='transparent', color='black'),
                       panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank(),
                       #panel.border = element_rect(fill='transparent', color='transparent'),
                       plot.title = element_text(size = 14), # centers, hjust = 0.5
                       plot.subtitle= element_text(size = 10),
                       #plot.caption = element_text()，
                       legend.key = element_rect( fill = "white"),
                       axis.title.x = element_text(vjust = -1.5, size = 14, colour = 'black'), # face = "bold"
                       axis.title.y = element_text(vjust = 1.5, size = 14, colour = 'black'), # face = "bold"
                       #axis.title=element_text(size = 12),
                       #axis.ticks=element_blank(),
                       #axis.line = element_line(colour = 'black'),
                       axis.ticks = element_blank(),
                       #axis.text.y = element_blank(),
                       axis.text.x = element_text(vjust = -0.5, size = 11, colour = 'black'),
                       axis.text.y = element_text(vjust = 0.5, size = 11, colour = 'black') )
  rj.ftheme.small <-   theme(panel.background = element_rect(fill='transparent', color='black'),
                             panel.grid.major=element_blank(),
                             panel.grid.minor=element_blank(),
                             #panel.border = element_rect(fill='transparent', color='transparent'),
                             plot.title = element_text(size = 8), # centers, hjust = 0.5
                             plot.subtitle= element_text(size = 10),
                             #plot.caption = element_text()，
                             legend.key = element_rect( fill = "white"),
                             #axis.title=element_text(size = 12),
                             #axis.ticks=element_blank(),
                             #axis.line = element_line(colour = 'black'),
                             axis.ticks = element_blank(),
                             #axis.text.y = element_blank(),
                             axis.text.x = element_text(vjust=0.5, size = 8)   )

  bi <- data.frame(K = murp$k,
                   BIC = murp$BIC)

  bi_local <- bi[order(bi$K,decreasing = FALSE),][1:10,]
  highlight_bi <- bi_local %>% filter(BIC==min(bi_local$BIC))

  p1 <- ggplot(bi, aes(x = K, y = BIC)) +
    geom_point(color = '#2F4F4F', alpha = 0.6, size = 3) +
    geom_vline(aes(xintercept = murp$Recommended_K), colour="#990000", linetype="dashed") +
    labs(title="Likelihood Function Value Corresponding to each K",
         x = "Number of MURPs / Cell Number", y = "pseudo-BIC") +
    theme(axis.text.y = element_text(hjust=0.5, angle = 90)) +
    rj.ftheme

  p2 <- ggplot(bi_local, aes(x = K, y = BIC)) +
    geom_point(color = 'sienna1', alpha = 0.6, size = 1.8) +
    geom_point(data = highlight_bi, aes(x = K, y = BIC), color = 'blue', alpha = 0.7) +
    geom_vline(aes(xintercept = murp$Recommended_K), colour="#990000", linetype="dashed") +
    labs(x='', y='') +
    # scale_y_continuous(labels = scientific,
    #                    limits = c(min(bi_local$loglike),
    #                               max(bi_local$loglike)),
    #                    breaks = seq(min(bi_local$loglike),
    #                                 max(bi_local$loglike),
    #                                 (max(bi_local$loglike)-min(bi_local$loglike)))) +
    # theme(axis.text.y = element_text(hjust=0.5, angle = 90)) +
    scale_y_continuous(breaks = NULL) +
    rj.ftheme.small

  g <- ggplotGrob(p2)
  p3 <- p1 + annotation_custom(g, xmin = max(bi$K)*0.52, xmax = max(bi$K)*0.92,
                               ymin = abs(max(bi$BIC))-( abs(max(bi$BIC))-abs(min(bi$BIC)) )/4*2.05,
                               ymax = max(bi$BIC) )

  return(p3)

}
