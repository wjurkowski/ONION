
# TRUE helpers.
readWithoutDuplicates <- function(filePath, header = TRUE) {
    transcriptomics <- read.table(filePath, header = header, row.names = NULL);
    transcriptomicsWithoutDuplicates <- data.frame(transcriptomics[!duplicated(transcriptomics[1]), ], row.names = 1)
    transcriptomicsWithoutDuplicates
}
