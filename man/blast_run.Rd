% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blastrun.R
\name{blast_run}
\alias{blast_run}
\title{Run Commandline Blast}
\usage{
blast_run(
  query,
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
  overwrite = FALSE
)
}
\arguments{
\item{query}{path to fasta file with sequences to blast (string)}

\item{db}{the database to blast against (defaults to the NCBI nucleotide database)}

\item{outfile_prefix}{prefix to prepend onto blast output files. Outputs will look like: <outfile_prefix>.blastn.tsv & <outfile_prefix>.blastn.config.tsv  (string)}

\item{evalue}{BLASTN parameter. evalue threshold to use (if unsure leave as is) (string)}

\item{max_target_seqs}{BLASTN parameter. maximum target seqs to use (also somewhat informs search space of heuristics) (int)}

\item{perc_identity}{BLASTN parameter. minimum threshold for percentage identity between subject and query sequences for a match to count as a 'hit/hsp' (int)}

\item{qcov_hsp_perc}{BLASTN parameter. minimum threshold for percentage of query sequence covered by a hit (hsp) for inclusion in output report (int)}

\item{outfmt}{BLASTN parameter. Columns to include in the blast output (string)}

\item{max_hsps}{BLASTN parameter. Maximum number of high scoring pairs (hsps/hits) to report per subject sequence (int)}

\item{subject_besthit}{BLASTN paramater. If two HSPs are from one subject and one fully encapsulates the other, return only the larger HSP (boolean)}

\item{overwrite}{Overwrite existing blastn output files? (boolean)}
}
\value{
Path to main blastn output file
}
\description{
Runs a commandline blast against full NCBI nucleotide database (remote).
Only run on a small number of sequences (no more than a few hundred) or it will take forever.
}
\examples{
\dontrun{
query = system.file("testfiles/simulated_fasta/e_coli_1.100seqs.fasta",package = "blastparse")
blast_run(query, "my_outfile_prefix")
}
}
