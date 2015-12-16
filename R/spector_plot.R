#' Plots for SpECtOR
#'
#' @param id_file
#' @param res_df
#' @param res_p
#' @param f_head
#' @param plot_type
#' @param plot_var
#' @param out_F
#'
#' @return ggplot variables
#' @export
#'
spector_plot <- function(id_file = NULL,
                      res_df = NULL,
                      res_p = NULL,
                      f_head = FALSE,
                      plot_type = "boxplot",
                      plot_var =  "rms",
                      out_F = FALSE) {
  # Read id assignment file file
  if (is.null(id_file)) {
    message("No file provided for sample grouping")
  } else {
    id_df <- .read_id_assign(id_path = id_file, f_head = f_head)
  }

  if (is.null(res_df) && is.null(res_p)) {
    stop("Neither results data-frame nor results path provided")
  }

  if (is.null(res_df) && !is.null(res_p)) {
    res_df <- .spector_plot_read(res_p)
  }

  plot_df <- .spector_id_res(res_df = res_df, id_df = id_df)

  if (plot_var == "rms") {
    plot_df <- plot_df %>%
      dplyr::mutate(R_rms = replace(R_rms, which(R_rms == 0), NA)) %>%
      dplyr::mutate(res_y = 1 / R_rms) %>%
      dplyr::select(-c(R_a, R_rms))
  } else if (plot_var == "mean") {
    plot_df <- plot_df %>%
      dplyr::mutate(R_a = replace(R_a, which(R_a == 0), NA)) %>%
      dplyr::mutate(res_y = 1 / R_a) %>%
      dplyr::select(-c(R_a, R_rms))
  } else if (plot_var == "Df") {
    plot_df <- plot_df %>%
      dplyr::mutate(Df = replace(Df, which(Df == 0), NA)) %>%
      dplyr::mutate(res_y = Df) %>%
      dplyr::select(-c(Df))
  } else {
    stop("Please choose a valid plotting variable.
          Valid options are: c('rms', 'mean', 'Df')")
  }

  print(plot_type)
  spector_pl <- list()

  if ("boxplot" %in% plot_type) {
    spector_boxplot <- plot_box(plot_df, id_df)
    spector_pl$spector_boxplot <- spector_boxplot
    if (out_F != FALSE) {
      ggsave(filename = paste(out_F, "spector_box.pdf", sep =""), spector_boxplot)
      plot_path <- paste(out_F, "spector_box.Rdata", sep = "")
      save(spector_boxplot, file = plot_path)
    } else {
      message("Plot not saved no output folder provided")
    }
  }

  if ("circle" %in% plot_type) {
    spector_circplot <- plot_circ(plot_df)
    spector_pl$spector_circplot <- spector_circplot
    if (out_F != FALSE){
      ggsave(filename = paste(out_F, "spector_circ.pdf", sep =""), spector_circplot)
      plot_path <- paste(out_F, "spector_circ.Rdata", sep = "")
      save(spector_circplot, file = plot_path)
    } else {
      message("Plot not saved no output folder provided")
    }
  }

  if ("dot" %in% plot_type) {
    spector_dotplot <- plot_dot(plot_df)
    spector_pl$spector_dotplot <- spector_dotplot
    if (out_F != FALSE) {
      if (!is.null(spector_dotplot)) {
        ggsave(filename = paste(out_F, "spector_dot.pdf", sep =""),
          spector_dotplot)
      }
      plot_path <- paste(out_F, "spector_dot.Rdata", sep = "")
      save(spector_dotplot, file = plot_path)
    }
  }

  if (sum(c("boxplot", "circle", "dot") %in% plot_type) == 0) {
    stop("Please choose a valid plot_type")
  }
return(spector_pl)
}

plot_box <- function(plot_df, id_df) {

  if (!is.null(plot_df$s_prep)) {
    plotgg_box <- plot_df %>%
    group_by(id, s_prep, base_id, gr_bam)
  } else if (!is.null(plot_df$base_id)) {
    plotgg_box <- plot_df %>%
    group_by(id, base_id, gr_bam)
  } else if (!is.null(plot_df$gr_bam)) {
    plotgg_box <- plot_df %>%
    group_by(id, gr_bam)
  } else {
    plotgg_box <- plot_df %>%
    group_by(id)
  }


    plotgg_box <- plotgg_box %>%
    summarise(M = mean(res_y, na.rm = TRUE),
              Q0 = boxplot.stats(res_y)$stats[1],
              Q1 = boxplot.stats(res_y)$stats[2],
              Q2 = boxplot.stats(res_y)$stats[3],
              Q3 = boxplot.stats(res_y)$stats[4],
              Q4 = boxplot.stats(res_y)$stats[5])

  spector_boxplot <- ggplot2::ggplot(plotgg_box, aes(x = id, ymin = Q0, lower = Q1,
                            middle = Q2, upper = Q3, ymax = Q4)) +
    spector_theme() +
    xlab("Sample ID") +
    ylab("Roughness metric")


  if (!is.null(plotgg_box$s_prep)) {
    spector_boxplot <- spector_boxplot +
      scale_fill_manual(name = "Sample type",
                        breaks = .sample_col(id_df, plot_df$s_prep),
                        values = c(RColorBrewer::brewer.pal(
                              length(unique(plot_df$s_pre)) - 1, "Set1"), "gray")) +
        geom_boxplot(stat = 'identity', aes(fill = s_prep))
  } else if (!is.null(plotgg_box$gr_bam)) {
    spector_boxplot <- spector_boxplot +
        geom_boxplot(stat = 'identity', aes(fill = gr_bam)) +
        guides(fill = FALSE)
  } else {
    spector_boxplot <- spector_boxplot +
      geom_boxplot(stat = 'identity', aes(fill = id)) +
      guides(fill = FALSE)
  }

  if (!is.null(plotgg_box$gr_bam) && !is.na(plotgg_box$gr_bam)) {
    spector_boxplot <- spector_boxplot + facet_grid(.~gr_bam, scales = "free_x")
  }

  return(spector_boxplot)
}

plot_circ <- function(plot_df) {
  circ_df <- plot_df %>%
    dplyr::group_by(id, gr_bam) %>%
    dplyr::summarise(mu = mean(res_y, na.rm = TRUE),
              range = boxplot.stats(res_y)$stats[5] -
                      boxplot.stats(res_y)$stats[1])

  if (!is.null(plot_df$base_id)) {
    circ_base <- plot_df %>%
      dplyr::filter(base_id == id) %>%
      dplyr::group_by(id, base_id, gr_bam) %>%
      dplyr::summarise(mu = mean(res_y, na.rm = TRUE),
            range = boxplot.stats(res_y)$stats[5] - boxplot.stats(res_y)$stats[1])

    circ_df <- circ_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(mu_d = mu - circ_base$mu[circ_base$base_id == base_id],
                    range_d = range -
                    circ_base$range[circ_base$base_id == base_id])
    }

  spector_circplot <- ggplot2::ggplot(circ_df, aes(x = 1, y = 1)) +
    geom_point(aes(size = mu + range), fill = "dodgerblue", shape = 21) +
    geom_point(aes(size = mu + range), colour = "dodgerblue") +
    geom_point(aes(size = mu), fill = "darkseagreen", shape = 21) +
    scale_size_continuous(range = c(20, 50), guide = FALSE) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_blank(),
          panel.background = element_rect(colour="seashell", fill='seashell'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.key = element_blank(),
          strip.text.x = element_text(face = "plain", size=14),
          strip.background = element_rect(colour="seashell", fill='seashell')) +
    facet_wrap(~id)

    if (!is.null(plot_df$base_id)) {
      spector_circplot <- spector_circplot +
        geom_point(aes(size = mu + range, fill = range_d), shape = 21) +
        geom_point(aes(size = mu + range, colour = range_d)) +
        geom_point(aes(size = mu, fill = mu_d), shape = 21) +
        scale_colour_gradient2(name = "RM Quartile range",
          space = "Lab", high = "dodgerblue", low = "darkseagreen") +
        scale_fill_gradient2(name = "RM", space = "Lab", high ="darkorange")
      }


  return(spector_circplot)
}

plot_dot <- function(plot_df) {
  if (is.null(plot_df$base_id)) {
    warning("No baseline provided, therefore no 'dotplot' produced.
            \nThis plot type requires a baseline.\n", call. = FALSE)
    return()
  }

  dot_base <- plot_df %>%
    dplyr::filter(base_id == id) %>%
    dplyr::group_by(id, s_prep, base_id, gr_bam) %>%
    dplyr::summarise(mu = mean(res_y, na.rm = TRUE))

  dot_df <- plot_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(mu_d = res_y - dot_base$mu[dot_base$base_id == base_id])

  spector_dotplot <- ggplot2::ggplot(dot_df, aes(x = chr, y = res_y)) +
    geom_point(aes(col = mu_d)) +
    scale_colour_gradient2(name = "spector - spector(baseline)",
      space = "Lab", low ="gray", mid = "gray") +
    theme_bw() +
    xlab("Chromosome") +
    ylab("RM genomic region") +
    theme(strip.background = element_rect(colour = NA),
    axis.text.x = element_text(size=5, , angle = 45, hjust = 1,
      colour = "darkgray"),
    axis.text.y = element_text(size=11),
    axis.title.x = element_text(face='bold', size=22),
    axis.title.y = element_text(face='bold', angle=90, size=22),
    plot.title = element_text(face='bold', size=22),
    legend.key = element_blank(),
    panel.margin = grid::unit(0, "lines"),
    strip.text.x = element_text(size=12),
    strip.background = element_rect(colour="#DBDBDB", fill='#DBDBDB'))

  if (!is.null(dot_df$s_prep) && !is.na(dot_df$s_prep)) {
    spector_dotplot <- spector_dotplot + facet_grid(gr_bam ~ s_prep)
  } else if (!is.null(dot_df$gr_bam) && !is.na(dot_df$gr_bam)) {
    spector_dotplot <- spector_dotplot + facet_grid(id ~.)
  }

  return(spector_dotplot)
}

#==================================================================================
### Auxiliary function used for plotting
.spector_plot_read <- function(res_path) {
  res_file <- paste(res_path, list.files(path = res_path,
                                pattern = "*out_merged.csv"), sep = "")

  tmp <- dplyr::tbl_df(read.csv(res_file, header = TRUE,
                          stringsAsFactors = FALSE, strip.white = TRUE))
}

.spector_id_res <- function(res_df, id_df) {
  id_df <- .id_base(id_df)
  dplyr::left_join(x = res_df, y = id_df, by = c("id" = "id_bam"))
}

.id_base <- function(id_df) {
  id_df <- id_df %>%
    dplyr::mutate(fs_bam = gsub(".bam", "", basename(as.character(fs_bam))))

  if (!is.null(id_df$baseline)) {
    if (id_df$baseline[1] %in% id_df$fs_bam) {
      id_df <- id_df %>%
        dplyr::rowwise() %>%
        dplyr::mutate(base_id = id_df$id_bam[id_df$fs_bam == baseline])
    } else if (id_df$baseline[1] %in% id_df$id_bam) {
      id_df <- id_df %>%
        dplyr::rowwise() %>%
        dplyr::mutate(base_id = id_df$id_bam[id_df$id_bam == baseline])
    }
  }
  return(id_df)
}

.sample_col <- function(id_df, s_prep) {
  s_b <- .id_base(id_df) %>%
  dplyr::filter(fs_bam == baseline | id_bam == baseline) %>%
  dplyr::select(s_prep) %>%
  unique() %>%
  .$s_prep

  s_n <- unique(s_prep)

  c(s_n[s_n == s_b], s_n[!(s_n == s_b)])
}

spector_theme <- function(base_size = 11, base_family = "Helvetica") {
theme_bw() +
theme(
  text = element_text(family = base_family, face = "plain",
                      colour = "black", size = base_size,
                      hjust = 0.5, vjust = 0.5, angle = 0,
                      lineheight = 0.9),
  strip.background = element_rect(colour = NA),
  axis.text.x = element_text(size=10, angle = 45, hjust = 1, colour = "darkgray"),
  axis.text.y = element_text(size=14),
  axis.title.x = element_text(face='bold', size=22),
  axis.title.y = element_text(face='bold', angle=90, size=22),
  axis.ticks.x = element_line(colour = "darkgray"),
  plot.title = element_text(face='bold', size=22),
  legend.key = element_blank(),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 14),
  strip.text.x = element_text(size=12),
  strip.background = element_rect(colour="#DBDBDB", fill='#DBDBDB')
  )
}
