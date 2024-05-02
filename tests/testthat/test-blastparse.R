test_that("blast_parse runs without errors", {
  blast_results = system.file("testfiles/simulated_fasta/e_coli_1.100seqs.blastn.tsv",package = "blastparse")
  expect_error(blast_parse(blast_path = blast_results), regexp = NA)
})

test_that("blast_parse produces expected results", {
  blast_results = system.file("testfiles/simulated_fasta/e_coli_1.100seqs.blastn.tsv",package = "blastparse")
  blast_parsed = blast_parse(blast_path = blast_results)

  # Produces blast object
  expect_s4_class(blast_parsed, class = "Blast")

  # Blast dataframe is nonempty
  expect_gt(nrow(blast_parsed@blast_df), expected = 0)

})
