library(hexSticker)
library(tidyverse)
library(sysfonts)
library(geomtextpath)

font_add_google(name = "Lora", regular.wt = "600")
font_add_google(name = "Poppins", bold.wt = "800")

make_weave <- function(strand_col = "black", background_col = "white") {
  toplot <- 0:2 |>
    map(\(i) {
      offset <- i * 2 * pi/3

      tibble(
        x = seq(1, 9.5, 0.01),
        y = sin(x + offset),
        z = (x + offset) %/% (pi / 2)
      ) |>
        mutate(
          z = factor(z),
          z = as.numeric(z) + i * 100,
          order = (x + offset) %% (pi) < (pi/2)
        )
    }) |>
    list_rbind(names_to = "braid") |>
    arrange(order) |>
    mutate(z = factor(z, levels = unique(z)))

  toadd <- tibble(
    x = seq(1, 5, 0.01),
    y = sin(x +  4 * pi/3),
    text = "omics"
  )

  toplot |>
    ggplot(aes(x, y)) +
    ## Background sections
    geom_line(
      mapping = aes(group = z),
      data = filter(toplot, order),
      linewidth = 3,
      colour = background_col
    ) +
    geom_line(
      mapping = aes(group = z),
      data = filter(toplot, order),
      linewidth = 2,
      colour = strand_col,
      lineend = "round"
    ) +
    ## Foreground sections
    geom_line(
      mapping = aes(group = z),
      data = filter(toplot, !order),
      linewidth = 3,
      colour = background_col
    ) +
    geom_line(
      mapping = aes(group = z),
      data = filter(toplot, !order),
      linewidth = 2,
      colour = strand_col,
      lineend = "round"
    ) +
    # geom_textpath(
    #   mapping = aes(label = text),
    #   data = toadd,
    #   #family = "Poppins",
    #   fontface = "bold",
    #   colour = background_col,
    #   size = 3,
    #   text_only = TRUE,
    #   hjust = 0,
    #   spacing = -200,
    #   text_smoothing = 10
    # ) +
    scale_y_continuous(limits = c(-3, 3)) +
    theme_void()
}

col_border <- "#003399"
col_content <- "#ffffff"
col_background <- "#003399"
s <- sticker(
  make_weave(col_content, col_background),
  package = "moiraine",
  p_size = 20,
  p_color = col_content,
  p_x = 1,
  p_y = 1.2,
  p_family = "Lora",
  p_fontface = "italic",
  s_x = 1,
  s_y = 0.8,
  s_width = 1.3,
  s_height = 0.7,
  h_color = col_border,
  h_fill = col_background,
  filename = "man/figures/logo.png"
)
plot(s)
