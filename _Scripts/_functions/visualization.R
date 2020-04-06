### Code for plotting and visualizing TraceLab figures ###


### Import required packages ###

library(tidyr)
library(ggplot2)



### Plotting frames as points ###

plot_figure_pts <- function(px, py, res, col = NULL) {
  dat <- data.frame(x = px, y = py)
  if (is.null(col)) {
    plt <- ggplot(dat, aes(x = x, y = y))
  } else {
    dat$color <- col
    plt <- ggplot(dat, aes(x = x, y = y, color = color))
  }
  plt +
    geom_point(alpha = 0.5) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, res[1])) +
    scale_y_reverse(expand = c(0, 0), limits = c(res[2], 0)) +
    coord_fixed()
}


plot_tracing_pts <- function(px, py, tx, ty, res) {
  dat <- data.frame(
    n = 1:length(px),
    fig.x = px, fig.y = py,
    trace.x = tx, trace.y = ty
  )
  dat <- as_tibble(dat) %>%
    gather(2:5, key = "name", value = "val") %>%
    separate(name, c("type", "coord"), sep = "\\.") %>%
    spread(coord, val)
  ggplot(dat, aes(x = x, y = y, color = type)) +
    geom_point(alpha = 0.5) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, res[1])) +
    scale_y_reverse(expand = c(0, 0), limits = c(res[2], 0)) +
    coord_fixed()
}
