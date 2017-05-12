#' GIM
#'
#' Compute GIM for samples or groups of samples
#'
#' @param res_df results data frame
#' @param grp string. Variable that groups results into samples.
#' @param r_size integer. Region size used in computation. The defualt is `NULL`
#'        in which case it is automatically calculated form `res_df`.
#'
#' @return
#' @export
#'
#' @importFrom dplyr filter select_ mutate full_join select_ group_by_
#' @importFrom purrr map map_dbl
#' @importFrom tidyr nest
#' @importFrom stringr str_c
calc_gim <- function(res_df, grp = "id_bam", r_size = NULL) {

  if (is.null(r_size)) {
    r_size <- regionSize(res_df)
  }

  r_size <- as.integer(r_size)

  base_pdf <- base_gim %>%
    filter(region_size == r_size)

  res_df <- res_df %>%
    select_(grp, "las") %>%
    group_by_(grp)

  res_n <- res_df %>%
    summarise(n = n())

  res_flt <-  res_df %>%
    filter(las >= (min(base_pdf$x) - 1) & las <= max(base_pdf$x))

  flt_smr <- res_flt %>%
    summarise(n_flt = n()) %>%
    full_join(res_n, by = "id_bam") %>%
    mutate(reg_discarded = if_else(is.na(n_flt), "100%",
                   str_c(((n - n_flt) / n) * 100, "%")))


  res_flt %>%
    nest(las, .key = las) %>%
    mutate(
      pdf = map(las, samplePurrr, base_pdf$x),
      D_JS = map_dbl(pdf, djsPurrr, base_pdf$pdf)) %>%
    full_join(flt_smr, by = "id_bam") %>%
    select_(grp, quote(D_JS), quote(reg_discarded))
}

samplePurrr <- function(las, x_v) {
  samplePDF(las$las, breaks = c(x_v[1] - 1, x_v)) %>%
    with(pdf)
}


samplePDF <- function(data_v, breaks = NULL, min_brk = 50, max_brk = 300,
                       n_brks = 1) {

  if (!("numeric" %in% class(data_v))) {
    stop("'data_v' needs to be a vector")
  }
  if (is.null(breaks)) {
    bin_brks <- seq(min_brk, max_brk, by = n_brks)
  } else if (length(breaks) == 2) {
    bin_brks <- seq(breaks[1], breaks[2], by = n_brks)
  } else if (length(breaks) > 2) {
    bin_brks <- breaks
  }

  # Check if changes in max bin needed
  if (max(bin_brks) < max(data_v)) {
    bin_brks[length(bin_brks)] <- max(data_v)
  }

  tmp.freq <- hist(data_v, breaks = bin_brks, plot = F)
  tmp.freq <- tmp.freq$count / (length(data_v) * diff(bin_brks))


  data_frame(
    x = bin_brks[-n_brks],
    pdf = tmp.freq / sum(tmp.freq)
    )

}

djsPurrr <- function(pdf, base) {
  jsd(pdf, base)
}

jsd <- function(p_pdf, q_pdf) {
  p_bar <- 0.5 * (p_pdf + q_pdf)

  sqrt(0.5 * (KL_pdf(p_pdf, p_bar, 0) + KL_pdf(q_pdf, p_bar, 0)))

}

KL_pdf <- function(p_pdf, q_pdf, eps = 10^-10) {

  # for numeric stability numbers smaller than `eps` are replaced by `eps`
  if (!is.null(eps)) {
    p_pdf[p_pdf < eps] <- eps
    q_pdf[q_pdf < eps] <- eps
  }
  # ensure pdf not frequency

  tmp.D <- p_pdf * log(p_pdf / q_pdf)
  sum(if_else(is.nan(tmp.D), 0, tmp.D))
}
