#' @title seq_config
#' @description Preset parameters for different sequencing protocol
#' @param protocol The sequence protocol
#' @param toolkit The kit to build the library, should be 5, 3, or "3lax"
#' @param adapter The sequence of the adapter which is aside the cell barcode
#' @param left_flank,right_flank The length of the left/right part aside the adapter to be extracted
#' @param drop_adapter A flag to indicate if the adapter region should be removed from the tag region.
#' @param barcode_len,UMI_len The length of the cell barcode/UMI
seq_config = function(protocol,toolkit,
                      adapter,
                      left_flank,right_flank,drop_adapter,
                      barcode_len,UMI_len){
  if(protocol == "10X"){
    if(toolkit == 5){
      adapter = "TTTCTTATATGGG"
      UMI_len = 10
    }
    else if(toolkit == "3lax"){
      adapter = "GCGTCGTG"
      UMI_len = 12
    }
    else{
      adapter = "GCGTCGTGTAG"
      UMI_len = 12
    }
    left_flank = 55
    right_flank = -10
    drop_adapter = FALSE
    barcode_len = 16

  }
  else if(protocol == "VISIUM"){
    if(toolkit == 5){
      adapter = "TTTCTTATATGGG"
    }
    else{
      adapter = "GCGTCGTGTAG"
    }
    left_flank = 55
    right_flank = -10
    drop_adapter = FALSE
    barcode_len = 16
    UMI_len = 12
  }
  else if(protocol == "Curio"){
    adapter = "TCTCGGGAACGCTGAAGA"
    left_flank = 25
    right_flank = 25
    drop_adapter = TRUE
    barcode_len = 14
    UMI_len = 9
  }
  else if(protocol == "other"){
    if(is.null(adapter)){
      stop("Please sepcify the adapter, left_flank, right_flank for your protocol to search and extract the adapter sequence!")
    }
  }
  else{
    stop("The protocol you specified is not incorporated yet, please specify your config!")
  }

  out = list(adapter,left_flank,right_flank,drop_adapter,barcode_len,UMI_len)
  names(out) = c("adapter","left_flank","right_flank","drop_adapter","barcode_len","UMI_len")
  return(out)
}


#' @title extractTagBc
#'
#' @description Extract the tag region and identify cell barcode in the tag region
#' @details Extract the tag region for each read and trim it out to generate a polished read.
#' Then the tag region is used to identify cell barcode and UMI.
#' When cores > 1 the input FASTQ is split into chunks via \code{fastqSplit}, each chunk is
#' processed in parallel via \code{future_lapply}, and the polished FASTQs are concatenated
#' with \code{cat} before the barcode-matching step.
#'
#' @param fastq_path The path of the input fastq
#' @param out_name The path for the polished fastq
#' @param barcode_path The path for the barcode whitelist
#' @param toolkit The kit to build the library, should be 5, 3, or "3lax" (10X 3-prime with
#'   shortened 8 bp adapter "GCGTCGTG", window = 7, step = 1)
#' @param adapter The sequence of the adapter which is aside the cell barcode
#' @param window The window size to search the substring of adapter in the read.
#'   Defaults to 7 for toolkit "3lax", 10 otherwise.
#' @param step The step size to search the substring of adapter in the read.
#'   Defaults to 1 for toolkit "3lax", 2 otherwise.
#' @param left_flank,right_flank The length of the left/right part aside the adapter to be extracted
#' @param drop_adapter A flag to indicate if the adapter region should be removed from the tag region.
#' @param polyA_bin The window size to search for polyA.
#' @param polyA_base_count The minimum threshold of the number of A within the polyA search window
#' @param polyA_len The maximum length of polyA to be preserved in the read
#' @param mu The initial expectation of the start positions for cell barcodes in the tag region
#' @param sigma The initial variance of the start positions for cell barcodes in the tag region
#' @param k The length of barcode substring search in the tag region
#' @param batch The number of reads to be processed to update the distribution of
#' start position of cell barcodes
#' @param top The number of candidate barcodes which are further compared with edit distance
#' @param cos_thresh The minimum threshold for the cosine similarity between the barcode and
#' the tag region, only barcode with cosine similarity larger than the thresh can be preserved as
#' candidate barcodes.
#' @param alpha The size of the confidence interval for the start position of cell barcode alignment
#' @param edit_thresh The maximum threshold for the edit distance of the barcode alignment, alignment with
#' edit distance over the threshold would be filtered out.
#' @param barcode_len,UMI_len The length of the cell barcode/UMI
#' @param flank The length of flank when extract UMI, which is used to be tolerant of indels and dels.
#' @param fastq_batch The number of reads per FASTQ chunk when splitting for parallel tag extraction.
#'   Defaults to 10000000L. Set to cores * fastq_batch >= total_reads to minimise I/O round trips.
#' @param cores The number of cores to use for parallelization
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom stringi stri_reverse
#' @importFrom spgs reverseComplement
#' @importFrom future sequential
#' @importFrom future multisession
#' @importFrom future multicore
#' @importFrom future.apply future_lapply
#' @importFrom data.table rbindlist
#' @export
#'
extractTagBc = function(fastq_path, barcode_path, out_name,
                        # parameters to extract the tag region
                        toolkit, protocol = "10X", adapter = NULL,
                        window = NULL, step = NULL,
                        left_flank = 0, right_flank = 0, drop_adapter = FALSE,
                        polyA_bin = 20, polyA_base_count = 15, polyA_len = 10,
                        # parameters for barcode match
                        barcode_len = 16, mu = 20, sigma = 10, k = 6, batch = 100,
                        top = 5, cos_thresh = 0.25, alpha = 0.05,
                        edit_thresh = 3, mean_edit_thresh = 1.5,
                        UMI_len = 10, UMI_flank = 1,
                        # parameter for parallel
                        fastq_batch = 10000000L,
                        cores = 1){

  config = seq_config(protocol, toolkit,
                      adapter,
                      left_flank, right_flank, drop_adapter,
                      barcode_len, UMI_len)
  barcode_len = config$barcode_len
  UMI_len = config$UMI_len

  # resolve window/step defaults based on toolkit
  if(is.null(window)) window = if(toolkit == "3lax") 7L else 10L
  if(is.null(step))   step   = if(toolkit == "3lax") 1L else  2L

  # normalise toolkit to integer for the C++ layer (must be 3 or 5)
  toolkit_int = if(toolkit == "3lax") 3L else as.integer(toolkit)

  # P10: parallel FASTQ tag extraction
  if(cores > 1L){
    # Split the FASTQ into chunks; each chunk becomes a separate .fq.gz file
    split_dir = tempfile("longcell_split_")
    dir.create(split_dir)
    on.exit(unlink(split_dir, recursive = TRUE), add = TRUE)

    Longcellsrc::fastqSplit(fastq_path, split_dir, fastq_batch)

    chunk_files = sort(list.files(split_dir, pattern = "\\.fq\\.gz$", full.names = TRUE))
    if(length(chunk_files) == 0){
      stop("fastqSplit produced no output files. Check that fastq_path is readable.")
    }

    # Process each chunk in parallel; each worker writes its own polished FASTQ
    chunk_out_files = file.path(split_dir, paste0("out_", seq_along(chunk_files), ".fq.gz"))

    tag_list = future_lapply(
      seq_along(chunk_files),
      function(i){
        Longcellsrc::extractTagFastq(chunk_files[i], chunk_out_files[i],
                                     config$adapter, toolkit_int, window, step,
                                     config$left_flank, config$right_flank, config$drop_adapter,
                                     polyA_bin, polyA_base_count, polyA_len)
      },
      future.globals = list(chunk_files = chunk_files,
                            chunk_out_files = chunk_out_files,
                            config = config,
                            toolkit_int = toolkit_int,
                            window = window, step = step,
                            polyA_bin = polyA_bin,
                            polyA_base_count = polyA_base_count,
                            polyA_len = polyA_len),
      future.packages = c("Longcellsrc"),
      future.seed = TRUE
    )

    # Concatenate polished FASTQ chunks into the final output file
    produced_chunks = chunk_out_files[file.exists(chunk_out_files)]
    if(length(produced_chunks) > 0){
      system(paste("cat", paste(shQuote(produced_chunks), collapse = " "),
                   ">", shQuote(out_name)))
    }

    reads = as.data.frame(data.table::rbindlist(tag_list))
  } else {
    reads = Longcellsrc::extractTagFastq(fastq_path, out_name,
                                         config$adapter, toolkit_int, window, step,
                                         config$left_flank, config$right_flank, config$drop_adapter,
                                         polyA_bin, polyA_base_count, polyA_len)
  }

  if(length(reads) == 0 || nrow(reads) == 0){
    stop("Please check if your adapter sequence is correct!")
  }

  barcode = read.table(barcode_path)[,1]
  barcode = substr(barcode,1,barcode_len)

  if(toolkit_int == 5L){
    reads$tag = stringi::stri_reverse(reads$tag)
    barcode = stringi::stri_reverse(barcode)
  }
  else if(toolkit_int == 3L){
    barcode = spgs::reverseComplement(barcode,case = "upper")
  }

  bc = BarcodeMatch(reads$tag, barcode,
                    mu = mu, sigma = sigma, sigma_start = 10,
                    k = k, batch = batch,top = top, cos_thresh = cos_thresh, alpha = alpha,
                    edit_thresh = edit_thresh,mean_edit_thresh = mean_edit_thresh,
                    UMI_len = UMI_len,UMI_flank = UMI_flank,cores = cores)

  if(toolkit_int == 5L){
    bc$barcode = stringi::stri_reverse(bc$barcode)
  }
  else if(toolkit_int == 3L){
    bc$barcode = spgs::reverseComplement(bc$barcode,case = "upper")
  }

  bc = bc %>% filter(nchar(umi) == UMI_len+2*UMI_flank)
  orig_col = setdiff(colnames(bc),"id")
  bc = cbind(bc,reads[bc$id,c("name","tag","polyA")]) %>% dplyr::select(-id)

  bc = bc[,c("name","tag",orig_col,"polyA")]
  return(bc)
}


#' @title adapter_dis
#'
#' @description Calculate the needleman score between original adapter and adapter in reads.
#' @details Calculate the edit distance between original adapter and adapter in reads to represent the data quality.
#' Only adapters upstream confidently identified barcodes will be used.
#'
#' @param data The output from the barcode match step. Should be a dataframe contains edit distance and
#' the adapter sequence
#' @param flank The length of flank when extract UMI, which is used to be tolerant of indels and dels.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr filter
#' @importFrom dplyr summarise
#' @importFrom NameNeedle needleScores
#' @return A dataframe with two columns, the first is the edit distance and the secons is
#' how many adapters have such an edit distance to the original adapter
adapter_dis = function(data){
  len = unique(nchar(data$umi))
  adapter = data %>%
            filter(edit == 0) %>%
            filter(nchar(adapter) == len)

  adapter_count = table(adapter$adapter)
  adapter_seq = names(adapter_count)[adapter_count == max(adapter_count)][1]

  adapter_dis = data.frame(
    needle = needleScores(adapter_seq, names(adapter_count)),
    count = as.integer(adapter_count)
  )
  adapter_dis = adapter_dis %>% group_by(needle) %>% summarise(count = sum(count), .groups = "drop")

  return(adapter_dis)
}