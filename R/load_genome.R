#' @title load_genome
#' @description Load the genome from BSgenome given the genome name
#' @param genome_name The name of the genome, can be abbreviations for some
#'   commonly used genomes like "hg38"
#' @importFrom BSgenome available.genomes
#' @importFrom BSgenome getBSgenome
#' @return A BSgenome object
load_genome <- function(genome_name) {
  genome_list = BSgenome::available.genomes()

  if(genome_name %in% genome_list){
    genome_package = genome_name
  }
  # Handle prefix dynamically based on the genome_name and organism
  else {
    # Infer organism prefix from genome_name (a simplified heuristic)
    if (startsWith(genome_name, "hg")) {
      organism_prefix <- "BSgenome.Hsapiens.UCSC."
    } else if (startsWith(genome_name, "mm")) {
      organism_prefix <- "BSgenome.Mmusculus.UCSC."
    } else if (startsWith(genome_name, "rn")) {
      organism_prefix <- "BSgenome.Rnorvegicus.UCSC."
    } else if (startsWith(genome_name, "dm")) {
      organism_prefix <- "BSgenome.Dmelanogaster.UCSC."
    } else if (startsWith(genome_name, "ce")) {
      organism_prefix <- "BSgenome.Celegans.UCSC."
    } else if (startsWith(genome_name, "dr")) {
      organism_prefix <- "BSgenome.Drerio.UCSC."
    } else {
      stop("Unknown genome name or unsupported organism. Please use the full name for the genome_name or check if the genome exists in BSgenome::available.genomes()")
    }

    # Construct the full package name
    genome_package <- paste0(organism_prefix, genome_name)
  }

  # Fail fast if the BSgenome package is not installed — never attempt to
  # install packages at runtime, which is fragile on HPC compute nodes
  # (no internet, read-only library paths) and non-reproducible.
  if (!requireNamespace(genome_package, quietly = TRUE)) {
    stop(
      "The BSgenome package '", genome_package, "' is not installed.\n",
      "  Install it before running the pipeline with:\n",
      "    BiocManager::install(\"", genome_package, "\")"
    )
  }

  # Load the genome using getBSgenome()
  genome <- tryCatch({
    BSgenome::getBSgenome(genome_package)
  }, error = function(e) {
    stop(paste("Failed to load genome package:", genome_package, "\nError:", e$message))
  })

  return(genome)
}
