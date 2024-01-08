# Functions BM_Seq
# Axel Rosendahl Huber 10-1-2018
# Aim: Generate barplots for the contributions of signatures
# Contributions are derived from MutationalPatterns nmf_res function
#=========================  FUNCTION PLOT CONTRIBUTION ============================#

contribution_bar_plot <- function(signature_contribution, times_sd = 1, selection = (0:length(signature_contribution))) {
  contribution_select <- signature_contribution[,selection]
  contribution_corr <- contribution_select/sum(contribution_select)
  h <- data.frame(rowSums(contribution_corr))
  colnames(h) <- "contribution"
  h$Signatures <- rownames(h)
  sd <- apply(contribution_corr, 1, sd)*times_sd
  h$sdmax <- h$contribution + sd
  h$sdmin <- h$contribution - sd
  print(h)
  ggplot(h,  aes(x = reorder(Signatures, contribution),  y= contribution, ymax = sdmax, ymin = sdmin)) + 
    geom_bar(stat = "identity", fill = c("#7CAE00","#F8766D","#00BFC4", "#C77CFF"), color = "black") + 
    coord_flip() + 
    geom_errorbar(width = 0.2, color = "black") + 
    theme_classic() 
}
