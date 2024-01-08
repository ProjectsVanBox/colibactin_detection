# Function: Age_plot
# Axel Rosendahl Huber 11-1-2018
age_plot <- function(data_frame, path) {
      plot <- ggplot(data_frame,
                aes(x = Age, y = norm_muts )) +
                geom_point(shape=21, col = "black", fill = "red" ,size = 2) +
                geom_smooth(method = "lm", size = .6, color = 'red') +
                scale_color_manual(values = "#FF0000") +
                scale_fill_manual(values = "#FF0000") +
                ylab("SNVs per genome") +
                xlab("Age (years)") +
                scale_y_continuous(limits = c(0, 2000), breaks = c(0,1000,2000)) +
                theme_BM()
      ggsave(path, plot, width = 5, height = 5)
      return(plot)
}
age_plot_0 <- function(data_frame,  xlimit) {     
              plot <- ggplot(data_frame, aes(x = Age, y = norm_muts)) +
                   geom_point(shape=21,  size = 2.3, col = "black", fill = "red", alpha = 0.7) +
                   geom_smooth(method = "lm", size = .6, color = "red", se = F) +
                   ylab("SNVs per genome") +
                   xlab("Age (years)") +
                   scale_x_continuous(limits = c(0, xlimit), breaks = seq(0, xlimit, by = 20)) +
                   expand_limits(y = 0) +  scale_y_continuous(expand = c(0, 0), limits = c(0, 1500), breaks = seq(0, 2000, 500)) +
                   axis.line.x
                   theme_BM()
            
  return(plot)
}
