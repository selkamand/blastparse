
blast_summary_table <- function(blast){
  assert_is_blast_object(blast)


  dplyr::tibble(
    `Query Sequences` = as.numeric(blast@blast_config$num_seqs_blasted),
    `Queries with a strong hits` = dplyr::n_distinct(blast@blast_df[["qaccver"]][blast@blast_df[["strong_hit"]]]),
    `% Queries with a strong hit` = `Queries with a strong hits`/ `Query Sequences` * 100,
  )
}


blast_summary_donut <- function(blast){
  df = dplyr::tibble(
    Property = c("Strong hit", "No strong hits"),
    Counts = c(
      dplyr::n_distinct(blast@blast_df[["qaccver"]][blast@blast_df[["strong_hit"]]]),
      as.numeric(blast@blast_config$num_seqs_blasted) - dplyr::n_distinct(blast@blast_df[["qaccver"]][blast@blast_df[["strong_hit"]]])
      )
    )

  ggiraphExtra::ggDonut(df, ggplot2::aes(donuts=Property, count=Counts), interactive = TRUE, addDonutLabel = TRUE, title = "Proportion of sequences with at Least 1 strong hit")
}

#' Donut Pie of Most Common Microbes
#'
#' @param blast blast object. See [blast_parse()]
#'
#' @return ggiraph object
#' @export
blast_top_microbes_donut <- function(blast){
  assert_is_blast_object(blast)
  tophits_df = blast_top_hits(blast)

 ggiraphExtra::ggPieDonut(
   data = tophits_df,
   ggplot2::aes(pies = scientific_name_genus, donuts = scientific_name_species),
   title = "Distribution of Top Microbial Hits",
   interactive = TRUE
   )
}

blast_species_level_summary <- function(tophits_df){
  species_level_summary = tophits_df |>
    dplyr::group_by(.data[["scientific_name_species"]]) |>
    dplyr::summarise(
      n_queries_tophit = dplyr::n_distinct(qaccver),
      n_queries_uniquely_tophit = dplyr::n_distinct(qaccver[!multiple_species_in_tophits]),
      n_queries_non_uniquely_tophit = dplyr::n_distinct(qaccver[multiple_species_in_tophits])
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(scientific_name_species = forcats::fct_reorder(.data[["scientific_name_species"]], n_queries_tophit))
  #browser()

  # Add 'Unclassified'
  # species_level_summary <- rbind(
  #   species_level_summary,
  #   data.frame(
  #     scientific_name_species = "Unclassified",
  #     )
  # )
}


#' Blast Top Microbes
#' Takes a blast object and finds the most appopriate top hit for each.
#' Maybe this isn't the write way of looking at things. a lot of sequences might map
#' equally well to multiple sequences
#' Maybe just plot a barplot describing the number of sequences in which each species present as a tophit is present in the db
#'
#' @param blast blast object. see [blast_parse()]
#'
#' @return ggiraph plot
#' @export
blast_top_microbes_barplot <- function(blast, interactive=TRUE){
  assert_is_blast_object(blast)
  tophits_df = blast_top_hits(blast)

  species_level_summary_df <- blast_species_level_summary(tophits_df)

  gg = species_level_summary_df |>
    #dplyr::count(scientific_name_genus, scientific_name_species) |>
    ggplot2::ggplot(ggplot2::aes(
      y= scientific_name_species,
      x = n_queries_tophit,
      fill = scientific_name_species,
      tooltip = paste0(scientific_name_species, "\nTotal Seqs Blasted = ", blast@blast_config$num_seqs_blasted, "\nPercentage of total reads: ", utilitybeltfmt::fmtpercent(n_queries_tophit/as.numeric(blast@blast_config$num_seqs_blasted) * 100))
      )) +
    ggiraph::geom_col_interactive() +
    ggplot2::ylab("Genus") +
    ggplot2::xlab("# Query Sequences Classified") +
    ggplot2::ggtitle("Top Hit Per Read") +
    ggplot2::scale_fill_brewer(palette = "Set1") +
    utilitybeltgg::theme_fivethirtyeight_two() +
    utilitybeltgg::theme_legend_right() +
    ggplot2::guides("fill" = ggplot2::guide_legend(title = "Species"))
  ggi = ggiraph::girafe(ggobj = gg)

  if(interactive)
    return(ggi)
  else
    return(gg)

}

#' Top Hits 2d-scatter
#'
#' @param blast blast object. See [blast_parse()]
#' @param annotate_strong_hits Draw ggforce annotation box around strong hits (boolean)
#'
#' @return ggiraph object
#' @export
#'
blast_top_microbes_2dscatter <- function(blast, annotate_strong_hits = TRUE, zoomed_out = FALSE, legend_hide = FALSE, size = 1, svg_width = 6, svg_height = 5){
  assert_is_blast_object(blast)
  tophits_df = blast_top_hits(blast, add_unmapped_query_seqs = FALSE)

  total_blasted_seqs = as.numeric(blast@blast_config$num_seqs_blasted)

  gg = tophits_df |>
    ggplot2::ggplot(ggplot2::aes(x = evalue, y=pident, shape = strong_hit))

  if(annotate_strong_hits){
    gg <- gg + ggforce::geom_mark_rect(
      ggplot2::aes(
        fill = strong_hit,
        filter = strong_hit == TRUE,
        label = "Strong Hits",
        description = describe_frequency(scientific_name_species, n=3)
        ),
      show.legend = FALSE,
      label.fontsize = 9,
      expand = 0.03)
  }

  gg = gg +
    ggiraph::geom_point_interactive(
      ggplot2::aes(
        color=scientific_name_species,
        tooltip = paste0(scientific_name, "\nQuery Coverage:  ", utilitybeltfmt::fmtpercent(qcovs))
        ),
      size = size
      ) +
    #ggplot2::geom_point(ggplot2::aes(color=scientific_name_species, tooltip = paste0(scientific_name, "\n", qcovs))) +
    #ggforce::facet_zoom(y==strong_hit) +
    ggplot2::scale_x_continuous(trans = "log10", oob = scales::oob_squish_infinite) +
    ggplot2::xlab("E-value") +
    ggplot2::ylab("Percent Identity") +
    ggplot2::ggtitle("Quality of Top Hits", subtitle = paste0("for each of the ", nrow(tophits_df), " input sequences \nwhich had at least one hit")) +
    utilitybeltgg::theme_fivethirtyeight_two() +
    ggplot2::guides("shape" = ggplot2::guide_legend(title = "Strong Hit")) +
    ggplot2::guides("color" = ggplot2::guide_legend(title = "Species")) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::scale_fill_brewer(palette = "Set3", direction = 1) +
    ggplot2::scale_shape_manual(values =  c("TRUE"= 20,"FALSE"=21)) +
    utilitybeltgg::theme_legend_right()

  if(zoomed_out){
    gg <- gg + ggplot2::coord_cartesian(ylim=c(0,100), xlim=c(NA, 10))
  }

  if(legend_hide){
    gg <-  gg + ggplot2::theme(legend.position = "none")
  }

  ggi <- ggiraph::girafe(ggobj=gg,  width_svg = svg_width, height_svg= svg_height)

    return(ggi)
}

describe_frequency <- function(character_vector, n = min(3, dplyr::n_distinct(character_vector))){

  if(n == 0 || length(character_vector) == 0){ return(character(0))}

  character_vector = forcats::fct_lump_n(character_vector,n = n)

  tb <- sort(table(character_vector), decreasing = TRUE)
  tb_proportion <- tb/sum(tb)*100

  #browser()
  paste0(
    utilitybeltfmt::fmtpercent(tb_proportion, decimal_places = 0)," ", names(tb_proportion), "\n", collapse = ""
    )
}

#' Find top hits
#' Takes a blast object and finds the most appopriate top hit for each.
#' Maybe this isn't the write way of looking at things. a lot of sequences might map
#' equally well to multiple sequences
#' Maybe just plot a barplot describing the number of sequences in which each species present as a tophit is present in the db


#' Find top hits
#'
#' Takes a blast object and finds the absolute top hit for each.
#' Looks for minimum hit_rank (summed ranks of evalue, -perc_identity and -qcovs).
#' If ties exist will mark as 'multiple'
#'
#' @param blast Blast object (see [blast_parse])
#'
#' @return tibble with 1 row per query sequence
#' @export
blast_top_hits <- function(blast, add_unmapped_query_seqs = TRUE){
  blast_tophits_df <- blast@blast_df |>
    dplyr::group_by(qaccver) |>
    dplyr::slice_min(hit_rank, with_ties = TRUE) |>
    dplyr::distinct(saccver, .keep_all = TRUE) |>
    dplyr::mutate(multiple_species_in_tophits = dplyr::n_distinct(scientific_name_species) > 1) |>
    dplyr::mutate(multiple_genus_in_tophits = dplyr::n_distinct(scientific_name_genus) > 1) |>
    # dplyr::summarise(
    #   dplyr::across(dplyr::everything(), .fns = function(vec){
    #     if(length(vec) > 1)
    #         return(paste(unique(vec), collapse = ","))
    #     else
    #       return(vec)
    #     }),
    #   multiple = dplyr::n() > 1
    #   ) |>
    dplyr::ungroup()

#
#     dplyr::mutate(
#       multiple_best_hits = dplyr::n_distinct(scientific_name) > 1,
#       equally_good_hits = ifelse(multiple_best_hits, yes=paste0(unique(scientific_name), collapse = ", "), no = character(0))
#     ) |>
#     #dplyr::filter(scientific_name_species == get_most_frequent_value(scientific_name_species)) |>
#     #dplyr::slice_min(hit_rank, with_ties = FALSE) |>
#     dplyr::ungroup()

  if(add_unmapped_query_seqs){

    sequences_with_no_hits = blast@names_of_query_sequences[! blast@names_of_query_sequences %in% blast_tophits_df[["qaccver"]]]

    if (length(sequences_with_no_hits) > 0){
      blast_tophits_df = dplyr::add_row(
        blast_tophits_df,
        qaccver = sequences_with_no_hits,
        #scientific_name_genus = "No Hits",
        scientific_name = "No Hits",
        #scientific_name_species = "No Hits",
        lineage = "No Hits"
      )
    }
  }
  return(blast_tophits_df)
}

#' Find top hits
#'
#' Takes a blast object and finds the absolute top hit for each.
#' Looks for minimum hit_rank (summed ranks of evalue, -perc_identity and -qcovs).
#' Breaks ties based on which species has the greatest number of distinct subjects at the given hit_rank.
#'
#' @param blast Blast object (see [blast_parse])
#'
#' @return tibble with 1 row per query sequence
#' @export
blast_top_hits_old <- function(blast, add_unmapped_query_seqs = TRUE){
  blast_tophits_df <- blast@blast_df |>
    dplyr::group_by(qaccver) |>
    dplyr::slice_min(hit_rank, with_ties = TRUE) |>
    dplyr::distinct(saccver, .keep_all = TRUE) |>
    dplyr::mutate(equally_good_hits = ifelse(dplyr::n_distinct(scientific_name) > 1, yes=paste0(unique(scientific_name), collapse = ", "), no = character(0))) |>
    dplyr::filter(scientific_name_species == get_most_frequent_value(scientific_name_species)) |>
    dplyr::slice_min(hit_rank, with_ties = FALSE) |>
    dplyr::ungroup()

  if(add_unmapped_query_seqs){

    sequences_with_no_hits = blast@names_of_query_sequences[! blast@names_of_query_sequences %in% blast_tophits_df[["qaccver"]]]

    if (length(sequences_with_no_hits) > 0){
      blast_tophits_df = dplyr::add_row(
        blast_tophits_df,
        qaccver = sequences_with_no_hits,
        #scientific_name_genus = "No Hits",
        scientific_name = "No Hits",
        #scientific_name_species = "No Hits",
        lineage = "No Hits"
     )
    }
  }
  return(blast_tophits_df)
}

blast_top_hits_simplify_for_sunburst <- function(blast_top_hits){

  add_prefix <- function(prefix,vector) {ifelse(is.na(vector), yes = NA_character_, no=paste0(prefix, vector))}

  counts_df <- blast_top_hits |>
    dplyr::count(scientific_name, scientific_name_species, scientific_name_genus, scientific_name_superkingdom, name = "n_tophits") |>
    dplyr::mutate(
      scientific_name = add_prefix("", scientific_name),
      scientific_name_species = add_prefix("(s) ", scientific_name_species),
      scientific_name_genus = add_prefix("(g) ", scientific_name_genus),
      scientific_name_superkingdom = add_prefix("(ss) ", scientific_name_superkingdom)
      )

  ultimate_parent = "Blasted\nSequences"
  #browser()
  dplyr::tibble(
    labels = vctrs::vec_c(
      counts_df$scientific_name,
      counts_df$scientific_name_species,
      counts_df$scientific_name_genus,
      counts_df$scientific_name_superkingdom,
      .ptype = character()),
    parents = vctrs::vec_c(
      counts_df$scientific_name_species,
      counts_df$scientific_name_genus,
      counts_df$scientific_name_superkingdom,
      rep(ultimate_parent, times=length(counts_df$scientific_name)),
      .ptype = character()),
    values = vctrs::vec_c(
      counts_df$n_tophits,
      rep(0, times = length(counts_df$scientific_name) *3),
      .ptype = numeric())
    ) |>
    dplyr::distinct(labels, parents, values) |>
    dplyr::filter(!is.na(labels)) |>
    tidyr::replace_na(replace = list(parents = ultimate_parent))
  #browser()
  # blast_top_hits |>
  #   dplyr::mutate(
  #     lineage = paste0(",", lineage)
  #     )
}

blast_simplify_for_sunburst_to_sunburst <- function(simplified_tophits){

  fig <- plotly::plot_ly(
    labels = simplified_tophits$labels,
    parents = simplified_tophits$parents,
    values = simplified_tophits$values,
    type = "sunburst",
    textinfo='label+percent entry',
    hoverinfo='label+percent entry+parent'
    )

  fig
}

get_most_frequent_value <- function(vec){
  names(sort(table(vec), decreasing = TRUE))[1]
}


#' Blast Tabular tophits
#'
#' @param blast Blast object (see [blast_parse])
#'
#' @return reactable
#' @export
#'
blast_tabular_tophits_reactable <- function(blast){
  assert_is_blast_object(blast)
  blast_tophits_df = blast_top_hits(blast, add_unmapped_query_seqs = FALSE)
  data = blast_tophits_df |>
    dplyr::select(
      qaccver,
      saccver,
      stitle,
      length,
      pident,
      evalue,
      bitscore,
      qcovs,
      qcovhsp,
      staxid,
      strong_hit,
      scientific_name,
      equally_good_hits
    ) |>
    dplyr::mutate(
      evalue_color_ref = dplyr::case_when(
        evalue < 1e-50 ~ "lightgreen",
        evalue < 1e-10 ~ "lightblue",
        TRUE ~ "grey"
        ),
      pident_color_ref = dplyr::case_when(
        pident > 97 ~ "lightgreen",
        pident > 90 ~ "lightblue",
        pident > 70 ~ "grey",
        TRUE ~ "lightpink",
        ),
      qcovhsp_color_ref = dplyr::case_when(
        qcovhsp > 97 ~ "lightgreen",
        qcovhsp > 90 ~ "lightblue",
        qcovhsp > 70 ~ "grey",
        TRUE ~ "lightpink",
      ),
      qcovs_color_ref = dplyr::case_when(
        qcovs > 97 ~ "lightgreen",
        qcovs > 90 ~ "lightblue",
        qcovs > 70 ~ "grey",
        TRUE ~ "lightpink",
      ),
      strong_hit_color_ref = ifelse(strong_hit, "lightgreen", "lightpink"),
      strength = ifelse(strong_hit, "confident", "putative")
      ) |>
    dplyr::select(-strong_hit)

  reactable::reactable(
    data,
    theme = reactablefmtr::fivethirtyeight(),
    filterable = TRUE,
    searchable = TRUE,
    wrap=FALSE,
    columnGroups = list(
      reactable::colGroup(name = "Quality", columns = c("length", "pident","evalue","bitscore","qcovs","qcovhsp", "equally_good_hits")),
      reactable::colGroup(name = "Subject", columns = c("saccver","stitle","scientific_name","staxid"))#,
      #reactable::colGroup(name = "Query", columns = c("qaccver")),
      ),
    columns = list(
      pident = reactable::colDef(cell = reactablefmtr::data_bars(data, fill_color_ref = "pident_color_ref")),
      qcovhsp = reactable::colDef(cell = reactablefmtr::data_bars(data, fill_color_ref = "qcovhsp_color_ref")),
      qcovs = reactable::colDef(cell = reactablefmtr::data_bars(data, fill_color_ref = "qcovs_color_ref")),
      strength = reactable::colDef(cell = reactablefmtr::pill_buttons(data = data, color_ref = 'strong_hit_color_ref', box_shadow = TRUE)),
      evalue = reactable::colDef(cell = reactablefmtr::pill_buttons(data = data, color_ref = 'evalue_color_ref', box_shadow = TRUE)),
      strong_hit_color_ref = reactable::colDef(show = FALSE),
      evalue_color_ref = reactable::colDef(show = FALSE),
      pident_color_ref = reactable::colDef(show=FALSE),
      qcovhsp_color_ref = reactable::colDef(show=FALSE),
      qcovs_color_ref = reactable::colDef(show=FALSE)
      ),
    resizable = TRUE, pagination = FALSE
    )
}



blast_sunburst <- function(blast){

  blast_tophits_df <- blast_top_hits(blast)

  browser()
  # blast_species_counts_df <- dplyr::count(x = blast_tophits_df,
  #   scientific_name_genus,
  #   scientific_name_species,
  #   scientific_name, sort = TRUE, name = "number_of_tophits"
  #   )
  #
  # blast_species_counts_df <- dplyr::mutate(.data = blast_species_counts_df,
  #     scientific_name_parent = scientific_name_species,
  #     scientific_name_species_parent = scientific_name_genus,
  #     scientific_name_genus_parent = ""
  #     )
  #
  # nspecies = nrow(blast_species_counts_df)
#
#   longform_df = dplyr::tibble(
#     labels = vctrs::vec_c(
#       blast_species_counts_df$scientific_name,
#       blast_species_counts_df$scientific_name_species,
#       blast_species_counts_df$scientific_name_genus,
#       .ptype = character()),
#     labeltype = rep(c("scientific_name", "scientific_name_species", "scientific_name_genus"), each=nspecies),
#     parent = c(blast_species_counts_df$scientific_name_parent, blast_species_counts_df$scientific_name_species_parent, blast_species_counts_df$scientific_name_genus_parent),
#     value = c(blast_species_counts_df$number_of_tophits, rep(10, times = nspecies), rep(10, times = nspecies))
#     )


  fig <- plotly::plot_ly(
    labels = as.character(longform_df$labels),
    parents = as.character(longform_df$parent),
    values = as.character(longform_df$value),
    type = 'sunburst'
  )

  # labels = c("Listeria monocytogenes", "Escherichia coli",
  #            "Escherichia coli O157:H7", "Salmonella enterica subsp. enterica serovar Typhimurium",
  #            "Salmonella enterica subsp. enterica serovar Hissar", "Listeria monocytogenes EGD-e",
  #            "Salmonella enterica subsp. enterica serovar Gallinarum", "Salmonella enterica subsp. enterica serovar Indiana",
  #            "Salmonella enterica", "Listeria monocytogenes", "Escherichia coli",
  #            "Escherichia coli", "Salmonella enterica", "Salmonella enterica",
  #            "Listeria monocytogenes", "Salmonella enterica", "Salmonella enterica",
  #            "Salmonella enterica", "Listeria", "Escherichia", "Escherichia",
  #            "Salmonella", "Salmonella", "Listeria", "Salmonella", "Salmonella",
  #            "Salmonella")
  #
  # labeltype = c("scientific_name", "scientific_name",
  #               "scientific_name", "scientific_name", "scientific_name", "scientific_name",
  #               "scientific_name", "scientific_name", "scientific_name", "scientific_name_species",
  #               "scientific_name_species", "scientific_name_species", "scientific_name_species",
  #               "scientific_name_species", "scientific_name_species", "scientific_name_species",
  #               "scientific_name_species", "scientific_name_species", "scientific_name_genus",
  #               "scientific_name_genus", "scientific_name_genus", "scientific_name_genus",
  #               "scientific_name_genus", "scientific_name_genus", "scientific_name_genus",
  #               "scientific_name_genus", "scientific_name_genus")
  # parent = c(
  #   "Listeria monocytogenes",
  #   "Escherichia coli", "Escherichia coli", "Salmonella enterica",
  #   "Salmonella enterica", "Listeria monocytogenes", "Salmonella enterica",
  #   "Salmonella enterica", "Salmonella enterica", "Listeria", "Escherichia",
  #   "Escherichia", "Salmonella", "Salmonella", "Listeria", "Salmonella",
  #   "Salmonella", "Salmonella", "", "", "", "", "", "", "", "", ""
  #   )
  # value = c(40, 29, 21, 19, 14, 10, 8, 8, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  # dput(longform_df)

  browser()
  #blast_species_counts_df
  #blast_species_counts_df
}
