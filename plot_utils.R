library(alphahull)
library(scales)

#' Plot focus, position and show strategy regions based on downward null
#' heuristic.
#'
#' @param data A data frame containing one row with columns `focus` and `position`
#' @param blur_data A data frame that is created by `dom_make_blur_data()`
#'
#' @return
#' @export
#'
#' @examples
adj_dom_plot_strategy <- function(data, blur_data, show_data_ci = FALSE, show_legend = FALSE) {
  # TODO: Check input data is in the right format (option to not have ci in data if show_data_ci == FALSE)
  # Get convex hull of blur data with error bars
  ahuld.conv <- convex_hull(blur_data)
  
  polygons <- make_polygons(blur_data)
  bully.poly <- polygons[["bully"]]
  clcomps.poly <- polygons[["clcomps"]]
  undef.poly <- polygons[["undef"]]
  
  ########## plot showing strategy polygons
  # plot real data
  with(data, plot(focus, position,
                  ylim = c(-0.1, 1), xlim = c(-0.1, 1), las = 1,
                  col = "blue", pch = 5,
                  xlab = "focus", ylab = "position", cex = 2
              
  ))
  
  # Add strategy polygons to bottom plotted layer
  # undefined (light grey)
  polygon(undef.poly, col = scales::alpha("grey", 0.4),
          border = scales::alpha("black", 0.5), lty = 2, lwd = 1.5)
  
  # bully (blue)
  polygon(bully.poly, col = scales::alpha("blue", 0.2),
          border = scales::alpha("black", 0.5), lty = 2, lwd = 1.5)
  
  # close competitors (red)
  polygon(clcomps.poly, col = scales::alpha("red", 0.2),
          border = scales::alpha("black", 0.5), lty = 2, lwd = 1.5)
  
  # downward heuristic (black, plotted over white)
  polygon(ahuld.conv, col = scales::alpha("white", 1), lty = 2, lwd = 1.5) # plot white without border
  polygon(ahuld.conv, col = scales::alpha("black", 0), # plot black border only
          border = scales::alpha("black", 0.4), lwd = 1.5)
  
  
  if(show_legend){
    legend("bottomright",
           c(
             "Downward heuristic",
             "Bullying",
             "Close competitors",
             "Undefined"
           ),
           # bty="n",
           col = "black",
           pt.bg = c(
             scales::alpha("black", 0.2),
             scales::alpha("blue", 0.2),
             scales::alpha("red", 0.2),
             scales::alpha("grey", 0.2)
           ),
           pch = 22, cex = 1, pt.cex = 1.5
    )
  }
  
  # add points and segments to upper layer
  
  # plot line for increasingly blurred downward heuristic, text above top of high position errorbars
  lines(x = blur_data$focus, y = blur_data$position)
  points(x = blur_data$focus, y = blur_data$position)
  text(x = blur_data$focus, y = blur_data$position_ci_hi + 0.025,
       label = blur_data$blur, cex = 0.6)
  
  # horizontal focus error, at blurred position
  graphics::segments(
    x0 = blur_data$focus_ci_lo, x1 = blur_data$focus_ci_hi,
    y0 = blur_data$position, y1 = blur_data$position,
    col = scales::alpha("black", 0.7), lwd = 1.5
  )
  
  # vertical position error, at blurred focus
  graphics::segments(
    x0 = blur_data$focus, x1 = blur_data$focus,
    y0 = blur_data$position_ci_lo, y1 = blur_data$position_ci_hi,
    col = scales::alpha("black", 0.7), lwd = 1.5
  )
  
  # add real data
  if (show_data_ci) {
    # horizontal focus error, at real position
    graphics::segments(
      x0 = data$focus_ci_lo, x1 = data$focus_ci_hi,
      y0 = data$position, y1 = data$position,
      col = scales::alpha("blue", 0.7), lwd = 2
    )
    
    # vertical position error, at real focus
    graphics::segments(
      x0 = data$focus, x1 = data$focus,
      y0 = data$position_ci_lo, y1 = data$position_ci_hi,
      col = scales::alpha("blue", 0.7), lwd = 2
    )
  }
  with(data, points(focus, position, col = "white", bg = "blue",
                    pch = 23, cex = 1.75))
  
}
