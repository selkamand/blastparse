assert_is_blast_object <- function(object){
  assertthat::assert_that(c("Blast") %in% class(object), msg = utilitybeltfmt::fmterror("Expected Blast object, got [", class(object), "]"))
}

#' S4 class
#'
#' An S4 representatio of blast objects. Please create using instances using [blast_parse()].
#'
#' @slot blast_df data.frame.
#' @slot blast_config list. list of paramaters related to running of blast. This field is parsed from the blastn.config.tsv file produced by [blast_run()]
#' @slot names_of_query_sequences character. A vector containing the names of the sequences analysed using blast.
#' @slot strong_hit_evalue numeric. The maximum threshold used to filter used to classify 'strong_hits'.
#' @slot strong_hit_perc_identity numeric. The minimum threshold used to filter used to classify 'strong_hits'
#' @slot strong_hit_query_coverage numeric. The minimum threshold used to filter used to classify 'strong_hits'
#' @details Strong Hit = evalue < strong_hit_evalue & perc_identity > strong_hit_perc_identity & qcovs > strong_hit_query_coverage.
#'
#' Required config file fields
#' percentage_identity,evalue_threshold,max_target_seqs,query,outfmt,num_seqs_blasted,max_hsps,subject_besthit,qcov_hsp_perc
#'
BlastInstantiate <- setClass(
  Class = "Blast",
  slots=list(
    blast_df="data.frame",
    blast_config="list",
    names_of_query_sequences="character",
    strong_hit_evalue="numeric", # stronghit = evalue < max_evalue & perc_identity > min_perc_identity & qcovs > min_percent_query_coverage
    strong_hit_perc_identity="numeric", # stronghit = evalue < max_evalue & perc_identity > min_perc_identity & qcovs > min_percent_query_coverage
    strong_hit_query_coverage="numeric" #stronghit = evalue < max_evalue & perc_identity > min_perc_identity & qcovs > min_percent_query_coverage
    ),
  )

setMethod("show",
          "Blast",
          function(object) {
            message("-----------------\nBlast Summary\n-----------------")
            message()
            message("-----------------\nBlast Config\n-----------------")
            message(paste0(names(object@blast_config), "\t", utilitybeltfmt::fmtbold(unlist(object@blast_config)), "\n"))
          }
)

check_blast_config <- function(object){
  errors = character()

  # Check all names are valid config fields
  invalid_config_fields <- object@blast_config[!names(object@blast_config) %in% blast_valid_config_fields()]
  missing_config_fields <- blast_valid_config_fields()[! blast_valid_config_fields() %in% names(object@blast_config)]

  if(length(invalid_config_fields) > 0){
    errors <- c(errors, paste0("Invalid field/s in blast config file: \n", paste0(invalid_config_fields, collapse=",")))
  }

  if(length(missing_config_fields) > 0){
    errors <- c(errors, paste0("Missing field/s in blast config file: \n", paste0(missing_config_fields, collapse=",")))
  }
}

setValidity(Class = "Blast", method = check_blast_config)

blast_valid_config_fields <- function(){
  blast_valid_config_fields = c(
    "percentage_identity",
    "evalue_threshold",
    "max_target_seqs",
    "query",
    "outfmt",
    "num_seqs_blasted",
    "max_hsps",
    "subject_besthit",
    "qcov_hsp_perc"
  )
  return(blast_valid_config_fields)
}


# BlastSummary <- function(blast){
#   assert_is_blast_object(blast)
#   number_of_query_sequences_with_at_least_one_hit = unique(blast@blast_df[["qaccver"]])
#   total_query_sequences = blast@blast_config$num_seqs_blasted
#   message(number_of_query_sequences_with_at_least_one_hit, " / ", total_query_sequences, " seqeunces blasted had at least one hit (see Blast Config for thresholds)")
# }
