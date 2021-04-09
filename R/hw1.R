#!/usr/bin/env Rscript

### Usage: Rscript --vanilla hw1.R <input file> <score file>
### Example: Rscript --vanilla hw1.R input.txt blosum62.txt

### Note: Smith-Waterman Algorithm

### This is one way to read in arguments in R
args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("At least two arguments must be supplied (inputFile, scoreFile).n",
       call.=FALSE) } else if (length(args)>=2) {
         # default gap penalties
         args[3] = -2
         args[4] = -1 }

## Specifying author and email
p <- c(person("Sarah", "Dudgeon", role = "aut", email =
                "sarah.dudgeon@yale.edu"))


## Load input file & make sequences
make.seqs <- function(inputFile) {
  # Load input file
  inF <- read.delim(inputFile, header = FALSE, stringsAsFactor = FALSE)
  # Assign sequences
  seq1 <<- inF[1,]
  seq2 <<- inF[2,]
}


## Make empty matrix
make.mat <- function(s1, s2) {
  # add space to top of string and tokenize components
  seq1_mat <- c(" ", strsplit(s1,"")[[1]])
  seq2_mat <- c(" ", strsplit(s2,"")[[1]])
  # make empty matrix with named row [seq1] and col [seq2]
  mat_empty <<- matrix(0,
                      nrow = length(seq1_mat),
                      ncol = length(seq2_mat),
                      dimnames = list(seq1_mat, seq2_mat))
}


## Update matrix
update.mat <- function(mat_input) {
  # make new matrix
  mat_full <<- mat_input
  for (i in rownames(mat_full)) {
    for (j in colnames(mat_full)) {
      # compute diag score
      diag_score <-((rownames(mat_full)[i] == colnames(mat_full)[j])*1) + mat_full[(i-1),(j-1)]
      # save diag score and corresponding seqs in a list of length=3
      diag_score <<- list(diag_score,
                          rownames(mat_full)[i],
                          colnames(mat_full)[j])
      # compute row scores
      # save row scores and corresponding seqs in a list of length=[[i]]

      # compute col scores
      # save col scores and corresponding seqs in a list of length=[[j]]

      # make final score list
      score_list <<- c(diag_score,
                      row_score,
                      col_score)
      ## Attach score_list as the 2nd value in the z plane of the updated matrix
      mat_full[i,j,2] <<- max(c(diag_score,))
      # make a list of values to check for max (every third)
      mat_check <- mat_full[i,j,2][seq(2, length(mat_full[i,j,2]), 3)]
      # take the max value of the list
      mat_full[i,j,3] <<- max(mat_check)
    }
  }

}


## Make output file
make.outputFile <- function(s1, s2, m, as, ar) {
  # create file
  file.create("output.txt")
  # set up a connector, where we will write
  fileConn <- file("output.txt")
  # write to the connector
  writeLines(c("-----------",
               "|Sequences|",
               "-----------",
               "sequence1",
               s1,
               "sequence2",
               s2,
               "--------------",
               "|Score Matrix|",
               "--------------",
               m,
               "----------------------",
               "|Best Local Alignment|",
               "----------------------",
               paste0("Alignment Score:",as),
               "Alignment Results:",
               ar), fileConn)
  close(fileConn)
}


## Implement your Smith-Waterman Algorithm
runSW <- function(inputFile, scoreFile, openGap = -2, extGap = -1) {
  ### load files
  make.seqs(inputFile)
  make.mat(seq1, seq2)
  ### calculation
  update.mat(mat_empty)
  ### write output
  make.outputFile()

}









## Run the main function and generate results
runSW(inputFile=args[1], scoreFile=args[2], openGap=args[3], extGap=args[4])
