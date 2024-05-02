#' Run Commandline Blast
#'
#' Runs a commandline blast against full NCBI nucleotide database (remote).
#' Only run on a small number of sequences (no more than a few hundred) or it will take forever.
#'
#' @param query path to fasta file with sequences to blast (string)
#' @param db the database to blast against (defaults to the NCBI nucleotide database)
#' @param outfile_prefix prefix to prepend onto blast output files. Outputs will look like: <outfile_prefix>.blastn.tsv & <outfile_prefix>.blastn.config.tsv  (string)
#' @param evalue BLASTN parameter. evalue threshold to use (if unsure leave as is) (string)
#' @param max_target_seqs BLASTN parameter. maximum target seqs to use (also somewhat informs search space of heuristics) (int)
#' @param perc_identity BLASTN parameter. minimum threshold for percentage identity between subject and query sequences for a match to count as a 'hit/hsp' (int)
#' @param qcov_hsp_perc BLASTN parameter. minimum threshold for percentage of query sequence covered by a hit (hsp) for inclusion in output report (int)
#' @param outfmt BLASTN parameter. Columns to include in the blast output (string)
#' @param max_hsps BLASTN parameter. Maximum number of high scoring pairs (hsps/hits) to report per subject sequence (int)
#' @param subject_besthit BLASTN paramater. If two HSPs are from one subject and one fully encapsulates the other, return only the larger HSP (boolean)
#' @param overwrite Overwrite existing blastn output files? (boolean)
#' @return Path to main blastn output file
#' @export
#'
#' @examples
#' \dontrun{
#' query = system.file("testfiles/simulated_fasta/e_coli_1.100seqs.fasta",package = "blastparse")
#' blast_run(query, "my_outfile_prefix")
#' }
blast_run <- function(query,
                      db = "nt",
                      remote = TRUE,
                      outfile_prefix,
                      evalue = "1e-10",
                      max_target_seqs = 20,
                      perc_identity = 90,
                      qcov_hsp_perc = 0,
                      outfmt = "6 std staxid qcovs qcovhsp stitle",
                      max_hsps = 10,
                      subject_besthit = TRUE,
                      overwrite = FALSE){

  utilitybeltassertions::assert_program_exists_in_path("blastn")
  utilitybeltassertions::assert_files_exist(query)
  utilitybeltassertions::assert_filenames_have_valid_extensions(query, valid_extensions = c("fasta", "fa", "fna"))

  # Set up output paths
  blast_output_path = paste0(outfile_prefix, ".blastn.tsv")
  blast_config_path = paste0(outfile_prefix, ".blastn.config.tsv")
  blast_seqnames_path = paste0(outfile_prefix, ".blastn.seqnames.tsv")

  # Make sure user doesn't overwrite files without knowing
  if(!overwrite & any(file.exists(blast_output_path, blast_config_path))){
    utilitybeltfmt::message_warning("There already exists a blast output/config files with prefix matching:\n\t[",outfile_prefix,"]\n.\nTo overwrite these files rerun this function with `overwrite=TRUE`")
    return(NA_character_)
  }

  # Make sure user doesn't try and blast too many sequences
  utilitybeltfmt::message_info("Counting sequences in input fasta ...")

  max_seqs_allowed_in_query = 500
  query_seqs=readLines(query)
  sequence_names = grep(x=query_seqs, pattern = "^>", value = TRUE)
  sequence_names = sub(x = sequence_names, pattern = "> ?", replacement = "")
  sequence_names = sub(x = sequence_names, pattern = "^ +", replacement = "")
  sequence_names = sub(x = sequence_names, pattern = " .*$", replacement = "")
  numseqs_input=length(sequence_names)
  message(numseqs_input, " sequences")

  assertthat::assert_that(numseqs_input <= max_seqs_allowed_in_query, msg = utilitybeltfmt::fmterror(
    "Woah there, [", numseqs_input,"] is a lot of input sequences. Maybe try 100-200 short ~ 100-200bp sequences instead. Or at least use fewer than [",max_seqs_allowed_in_query,"] sequences"
  ))

  # Writing names of all query sequences to file
  utilitybeltfmt::message_info("Writing names of input sequences to file [", blast_seqnames_path ,"] ...")
  write(x=sequence_names,file = blast_seqnames_path, ncolumns = 1, append = FALSE)



  # Write config file logging the settings we used
  config_file_tibble = dplyr::tibble(
    Property = c("evalue_threshold","max_target_seqs","percentage_identity","outfmt","max_hsps","subject_besthit", "query", "num_seqs_blasted", "qcov_hsp_perc"),
    Value = c(evalue,max_target_seqs,perc_identity,outfmt,max_hsps,subject_besthit, query, numseqs_input, qcov_hsp_perc)
    )

  # Make a config file describing some of the parameters we used
  utilitybeltfmt::message_info("Writing config file to  [", blast_config_path, "] ...")
  write.table(config_file_tibble, file = blast_config_path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


  # Blast sequences
  #---build command
  command = paste0(
    "blastn -db ", db," -out ", blast_output_path," -evalue ", evalue, " -max_target_seqs ", max_target_seqs, " -outfmt '",outfmt, "' -perc_identity " ,perc_identity, " -max_hsps ", max_hsps,
    ifelse(subject_besthit, " -subject_besthit", ""),
    ifelse(remote, " -remote", ""),
    " -query ", query
    )

  #--- run comand
  utilitybeltfmt::message_info("Running: \n", command)
  utilitybeltfmt::message_info("\n\tPlease be patient, this might take a moment ... ", symbol = FALSE)
  exit_code = system(command)

  #--- check it worked
  assertthat::assert_that(exit_code==0, msg = utilitybeltfmt::fmterror("Failed to run blastn"))
  utilitybeltassertions::assert_files_exist(c(blast_output_path, blast_config_path, blast_seqnames_path))

  # Final message
  utilitybeltfmt::message_success("Successfuly wrote blast file to ", blast_output_path)
  utilitybeltfmt::message_success("Successfuly wrote config file to ", blast_config_path)
  utilitybeltfmt::message_success("Successfuly wrote seqnames file to ", blast_seqnames_path)

  # return path to blast file
  return(blast_output_path)
}
