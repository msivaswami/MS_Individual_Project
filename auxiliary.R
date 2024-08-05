##################################################
# FUNCTION FOR HEAD-TAIL MS2 PLOT + DOT PROD COSINE
###################################################

ms2_plot = function(query, ref, precursor_mz = 0, ppm = 5){
  
  #run as follows: ms2_plot(query = query, ref = ref, ppm = 5)
  
  #query = data.frame(mz = c(100.0000, 101.0033, 250.0012, 251),
  #                      int = c(1000, 110, 3000, 3050))
  
  #ref = data.frame(mz = c(150.0000, 151.0033, 250.0014, 255),
  #                      int = c(1000, 110, 3000, 3050))
  
  mtx = data.frame(data_source = c('precursor', rep('query', nrow(query)), rep('ref', nrow(ref))),
                   mz = c(precursor_mz, query$mz, ref$mz),
                   int = c(precursor_mz, query$int, ref$int))
  
  mdist = stats::dist(mtx$mz)
  averageMzPair = stats::as.dist(outer(mtx$mz, mtx$mz, "+")/2)
  relativeErrors = averageMzPair * 0.000001
  m_massTolerance = mdist / relativeErrors
  clh <- fastcluster::hclust(m_massTolerance)
  # cut a desired level
  cl <- stats::cutree(clh, h = ppm)
  mtx$cl = cl
  
  mtx$key = ''
  mtx[mtx$cl==1,]$key = 'precursor'
  
  mtx = mtx %>% filter(data_source != 'precursor') %>% group_by(cl) %>% dplyr::mutate(key = ifelse(n_distinct(data_source) > 1, 'matched', 'unmatched')) %>%
    dplyr::mutate(int = ifelse(data_source != 'ref', int, int*-1))
  #mtx = mtx %>% dplyr::mutate(int = ifelse(data_source != 'ref', int, int*-1)) %>% filter(data_source != 'precursor')
  #return(p1)
  
  
  #calculate dot product
  add_ref = setdiff(dplyr::pull(mtx[mtx$data_source=='query','cl']), dplyr::pull(mtx[mtx$data_source=='ref','cl']))
  add_query = setdiff(dplyr::pull(mtx[mtx$data_source=='ref','cl']), dplyr::pull(mtx[mtx$data_source=='query','cl']))
  
  add_ref = mtx[mtx$cl %in% add_ref,]
  add_ref$data_source = 'ref'
  add_ref$int = 0
  
  add_query = mtx[mtx$cl %in% add_query,]
  add_query$data_source = 'query'
  add_query$int = 0
  
  mtx <- mtx %>% rbind(add_ref) %>% rbind(add_query)
  
  #create weighted intensity values
  mtx[mtx$data_source == 'query', 'mz']^2 * mtx[mtx$data_source == 'query', 'int']
  
  #
  mtx$w = (abs(mtx$int)^0.5) * mtx$mz^2
  mtx <- mtx %>% dplyr::arrange(data_source, cl)
  q_vec = pull(mtx[mtx$data_source=='query',])
  r_vec = pull(mtx[mtx$data_source=='ref',])
  
  q_vec * r_vec
  
  dpc = 100* (sum( mtx[mtx$data_source=='query', 'w'] * mtx[mtx$data_source=='ref', 'w']) / (sqrt(sum(mtx[mtx$data_source=='query', 'w']^2)) * sqrt(sum(mtx[mtx$data_source=='ref', 'w']^2))))
  n_matched <- mtx %>% filter(key == 'matched', data_source == 'ref') %>% nrow()
  n_total <- mtx %>% filter(w > 0, data_source == 'ref') %>% nrow()
  
  #head to tail plot for MS2 (upper = query, lower = reference)
  p1 <- mtx %>% ggplot(aes(x = mz, y = int)) + geom_segment(aes(x = mz, xend = mz, y = 0, yend=int, col = key), size = 1) +
    geom_point(size=ifelse(mtx$key=='matched', 2, 0), fill=ifelse(mtx$key=='matched', alpha("orange", 0.3), alpha("orange", 0.3)),
               alpha=ifelse(mtx$key=='matched', 0.7, 0), shape=ifelse(mtx$key=='matched', 16, 16), stroke=ifelse(mtx$key=='matched', 2, 0)) +
    geom_hline(yintercept = 0) +
    ggtitle('MS2 spectral comparison', subtitle = paste('Upper = query, Lower = reference, dpc = ', round(dpc,3),
                                                        ', n_matched = ', n_matched , ', n_total = ', n_total ,sep = '')) +
    theme_bw() +
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  return(list('plot' = p1, 'dpc' = dpc))
}




yaml_concatPath <- function(x) {

  #handler for yaml, concatenate paths based on !concat_path tag
  #make those dirs if they do not exist

    x = lapply(x, FUN=function(x){
      if(is.null(x)) x=assay
      return(x)})
    path = paste(x, collapse="/")

    print(path)

    if(file.exists(path)){
      print("The file exists")
    } else{
      print("The file does not exist")
    }
    return(path)
  }

yaml_extractYamlChunk <- function(yaml_params, level, name){
  
  names = unlist(lapply(yaml_params[[level]], FUN=function(x){
    return(x$name)}))
  
  name_idx = which(names == name)
  
  return(yaml_params[[level]][[name_idx]])
}

get_full_peak_width_xcms3 <- function(peaklist, obj){
  ###########################################
  # Get full peak width
  ###########################################
  # Args:
  #   peaklist: the peak list generated from either XCMS or CAMERA.
  #              Use the CAMERA peak list for this pipeline
  #   xsa: The CAMERA annotation object
  #
  # Returns:
  #   An updated peaklist with the full retention window ranges (and full mz ranges)
  #
  # See also:
  #   full_minmax, getpeaks, ldply (from the plyr library)
  
  message("Getting full peak widths")
  # Get 'peaks' (xcms features) from the XCMSnSet object
  # in the camera annotation object
  
  message("Get 'individual' peaks from camera-xcms object")
  
  #Modification by M Jones to allow for retrieving full chromatographic peak widths of the excluded peaks before removal from peaklist.
  #if(attributes(xsa)$class[1] != "xcmsSet"){
  #    obj = xsa@xcmsSet
  #}else{
  #    obj = xsa
  #}
  
  #rt.min = xcms::groupval(obj, method = "medret", value = "rtmin", intensity = "into")
  rt.min = xcms::featureValues(obj, method = "medret", value = "rtmin", intensity = "into")
  rt.min = apply(rt.min, 1, min, na.rm = TRUE)
  
  #rt.max = xcms::groupval(obj, method = "medret", value = "rtmax", intensity = "into")
  #rt.max = apply(rt.max, 1, min, na.rm = TRUE)
  rt.max = xcms::featureValues(obj, method = "medret", value = "rtmax", intensity = "into")
  #updated to use MAX (previously used min, which is opposite of what we're trying to find)
  rt.max = apply(rt.max, 1, max, na.rm = TRUE)
  
  #mz.min = xcms::groupval(obj, method = "medret", value = "mzmin", intensity = "into")
  mz.min = xcms::featureValues(obj, method = "medret", value = "mzmin", intensity = "into")
  mz.min = apply(mz.min, 1, min, na.rm = TRUE)
  
  #mz.max = xcms::groupval(obj, method = "medret", value = "mzmax", intensity = "into")
  #mz.max = apply(mz.max, 1, min, na.rm = TRUE)
  mz.max = xcms::featureValues(obj, method = "medret", value = "mzmax", intensity = "into")
  #updated to use MAX (previously used min, which is opposite of what we're trying to find)
  mz.max = apply(mz.max, 1, max, na.rm = TRUE)
  
  peaklist_full = cbind(peaklist, 
                        "mzmin_full" = unname(mz.min), 
                        "mzmax_full" = unname(mz.max), 
                        "rtmin_full" = unname(rt.min), 
                        "rtmax_full" = unname(rt.max))
  return(peaklist_full)
  
}


#@author G. Le Corguille
# This function convert if it is required the Retention Time in minutes
RTSecondToMinute <- function(variableMetadata, convertRTMinute) {
  if (convertRTMinute) {
    #converting the retention times (seconds) into minutes
    print("converting the retention times into minutes in the variableMetadata")
    variableMetadata[, "rt"] <- variableMetadata[, "rt"] / 60
    variableMetadata[, "rtmin"] <- variableMetadata[, "rtmin"] / 60
    variableMetadata[, "rtmax"] <- variableMetadata[, "rtmax"] / 60
  }
  return(variableMetadata)
}

#@author G. Le Corguille
# This function format ions identifiers
formatIonIdentifiers <- function(variableMetadata, numDigitsRT = 0, numDigitsMZ = 0) {
  splitDeco <- strsplit(as.character(variableMetadata$name), "_")
  idsDeco <- sapply(splitDeco,
                    function(x) {
                      deco <- unlist(x)[2]; if (is.na(deco)) return("") else return(paste0("_", deco))
                    }
  )
  namecustom <- make.unique(paste0("M", round(variableMetadata[, "mz"], numDigitsMZ), "T", round(variableMetadata[, "rt"], numDigitsRT), idsDeco))
  variableMetadata <- cbind(name = variableMetadata$name, namecustom = namecustom, variableMetadata[, !(colnames(variableMetadata) %in% c("name"))])
  return(variableMetadata)
}

#@author G. Le Corguille
# This function convert the remain NA to 0 in the dataMatrix
naTOzeroDataMatrix <- function(dataMatrix, naTOzero) {
  if (naTOzero) {
    dataMatrix[is.na(dataMatrix)] <- 0
  }
  return(dataMatrix)
}


getTables <- function(xdata, intval = "into", convertRTMinute = F, numDigitsMZ = 4, numDigitsRT = 0, naTOzero = T, sampleNamesList) {
  dataMatrix <- featureValues(xdata, method = "medret", value = intval)
  colnames(dataMatrix) <- make.names(tools::file_path_sans_ext(colnames(dataMatrix)))
  # dataMatrix <- cbind(name = groupnames(xdata), dataMatrix)
  rownames(dataMatrix) = groupnames(xdata)
  variableMetadata <- featureDefinitions(xdata)
  colnames(variableMetadata)[1] <- "mz"; colnames(variableMetadata)[4] <- "rt"
  variableMetadata <- data.frame(name = groupnames(xdata), variableMetadata)
  
  variableMetadata <- RTSecondToMinute(variableMetadata, convertRTMinute)
  variableMetadata <- formatIonIdentifiers(variableMetadata, numDigitsRT = numDigitsRT, numDigitsMZ = numDigitsMZ)
  dataMatrix <- naTOzeroDataMatrix(dataMatrix, naTOzero)
  
  # FIX: issue when the vector at peakidx is too long and is written in a new line during the export
  variableMetadata[, "peakidx"] <- vapply(variableMetadata[, "peakidx"], FUN = paste, FUN.VALUE = character(1), collapse = ",")
  
  return(list(variableMetadata, dataMatrix))
}


getPeakTable <- function(xset, intval="into") {
  variableMetadata_dataMatrix <- xcms::peakTable(xset, method="medret",
                                                 value=intval)
  variableMetadata_dataMatrix <- cbind(name=xcms::groupnames(xset),
                                       variableMetadata_dataMatrix)
  
  variableMetadata <-
    variableMetadata_dataMatrix[, !(
      colnames(variableMetadata_dataMatrix) %in% c(make.names(xcms::sampnames(xset))))]
  
  return(variableMetadata)
}


get_av_spectra <- function(x) {
  
  if (length(x$av_intra) > 0) {
    av_intra_df <- plyr::ldply(x$av_intra)
    
    if (nrow(av_intra_df) == 0) {
      av_intra_df <- NULL
    } else{
      av_intra_df$method <- "intra"
    }
    
  }else{
    av_intra_df <- NULL
  }
  
  if ((is.null(x$av_inter)) || (nrow(x$av_inter) == 0)) {
    av_inter_df <- NULL
  }else{
    av_inter_df <- x$av_inter
    av_inter_df$method <- "inter"
  }
  
  if ((is.null(x$av_all)) || (nrow(x$av_all) == 0)) {
    av_all_df <- NULL
  }else{
    av_all_df <- x$av_all
    av_all_df$method <- "all"
  }
  
  combined <- plyr::rbind.fill(av_intra_df, av_inter_df, av_all_df)
  
  return(combined)
}


tableAverageMS2spec = function(pa){
  
  if (length(pa) > 0) {
    
    av_spectra <- plyr::ldply(pa@av_spectra, get_av_spectra)
    
    if (nrow(av_spectra) == 0) {
      message("No average spectra available")
    } else {
      colnames(av_spectra)[1] <- "grpid"
      av_spectra$grpid <- names(pa@av_spectra)[av_spectra$grpid]
      
      if ((length(pa@av_intra_params) > 0) || (length(pa@av_inter_params) > 0)) {
        # Add some extra info (only required if av_intra or av_inter performed)
        colnames(av_spectra)[2] <- "fileid"
        av_spectra$avid <- seq_len(nrow(av_spectra))
        
        filenames <- sapply(av_spectra$fileid,
                            function(x) names(pa@fileList)[as.integer(x)])
        # filenames_galaxy <- sapply(
        #        av_spectra$fileid, function(x) basename(pa@fileList[as.integer(x)]))
        
        av_spectra <- as.data.frame(
          append(av_spectra, list(filename = filenames), after = 2))
      }
      return(av_spectra)
    }
  }
}


get_spectrum = function(raw_data_spec, rt){
  filt <- unlist(lapply(raw_data_spec, function(x) x@rt - as.numeric(rt)))
  name_spec = names(which.min(abs(filt)))
  return(raw_data_spec[[name_spec]])
}


get_peak_summary = function(pgs, rtr, mzr, exact_mass_st){
  rtr = rtr_st
  mzr = mzr_st
  fd = featureDefinitions(pgs, rt=rtr, mz=mzr)
  fv = featureValues(pgs, value = "into")
  cp_pgs = chromPeaks(pgs)
  l_features = vector(mode = "list", length = nrow(fd))
  
  if (nrow(fd) == 0){
    return(l_features)
  }
  
  for (i in 1:nrow(fd)){
    
    peak_idxs = unlist(fd[i, "peakidx"])
    cp = cp_pgs[peak_idxs, ]
    mz_values = as.numeric(cp[,"mz"])
    samples_idxs = cp[,"sample"]
    
    ppm_errors = (mz_values - exact_mass_st) / (0.000001 * exact_mass_st)
    da_errors = mz_values - exact_mass_st
    
    df_fv = data.frame(feature=rownames(fd)[1],
                       peak_idx=peak_idxs,
                       sample_name=colnames(fv)[samples_idxs],
                       sample_idx=samples_idxs,
                       mz = mz_values,
                       intensity = as.numeric(fv[1,][samples_idxs]),
                       ppm_error = ppm_errors,
                       da_error = da_errors)
    
    missing_samples = setdiff(1:length(fv[1,]), as.numeric(samples_idxs))
    if (length(missing_samples) > 0){
      for (j in 1:length(missing_samples)) {
        df_fv[nrow(df_fv) + 1,] = c(NA, NA, colnames(fv)[missing_samples[j]], 
                                    missing_samples[j], NA, NA, NA, NA)
      }
    }
    
    df_fv = df_fv[order(df_fv$sample_idx), ]
    df_fv$peak_idx = df_fv$peak_idx
    df_fv$mz = as.numeric(df_fv$mz)
    df_fv$intensity = as.numeric(df_fv$intensity)
    df_fv$ppm_error = as.numeric(df_fv$ppm_error)     
    df_fv$da_error = as.numeric(df_fv$da_error)
    l_features[[i]] = df_fv
  }
  
  return(l_features)
}


is.blank <- function(x, false.triggers=FALSE){
  if(is.function(x)) return(FALSE) # Some of the tests below trigger
  # warnings when used on functions
  return(
    is.null(x) ||        # Actually this line is unnecessary since
      length(x) == 0 ||    # length(NULL) = 0, but I like to be clear
      all(is.na(x)) ||
      all(x=="") ||
      (false.triggers && all(!x))
  )
}


write_tables_to_xlsx = function(pgs, fn_xlsx){
  
  tbls = getTables(pgs, naTOzero=FALSE) # variableMetadata, dataMatrix
  
  wb <- createWorkbook()
  addWorksheet(wb, "variableMetaData")
  writeData(wb, "variableMetaData", tbls[[1]], colNames=TRUE, keepNA=FALSE)
  
  addWorksheet(wb, "dataMatrix")
  writeData(wb, "dataMatrix", tbls[[2]], rowNames=TRUE, colNames=TRUE, keepNA=FALSE)
  
  addWorksheet(wb, "metaData")
  writeData(wb, "metaData", phenoData(pgs)@data, colNames=TRUE, keepNA=FALSE)
  
  saveWorkbook(wb, file = fn_xlsx, overwrite = TRUE)
  
}


update_file_paths = function(obj, path_raw_data){
  for (i in 1:length(obj@processingData@files)){
    print(paste("Set mzML file", paste(path_raw_data, basename(obj@processingData@files[i]), sep=path_sep)))
    obj@processingData@files[i] = paste(path_raw_data, basename(obj@processingData@files[i]), sep=path_sep)
    # pgs@processingData@files[i] = paste(path_raw_data, basename(pgs@processingData@files[i]), sep=path_sep)
  }
  return(obj)
}


groupNames_custom <- function(object, mzdec = 0, rtdec = 0){
  if('mz' %in% unlist(attributes(object)) == FALSE){
    stop('object does not have attribute "mz"')
  }
  if('rt' %in% unlist(attributes(object)) == FALSE){
    stop('object does not have attribute "rt')
  }
  mzfmt <- paste0("%.", mzdec, "f")
  rtfmt <- paste0("%.", rtdec, "f")
  gnames <- paste0("M", sprintf(mzfmt, object$mz),"T", sprintf(rtfmt, object$rt))
  if (any(dup <- duplicated(gnames))){
    for (dupname in unique(gnames[dup])) {
      dupidx <- which(gnames == dupname)
      gnames[dupidx] <- paste(gnames[dupidx], seq(along = dupidx), sep = "_")
    }
    object$name = gnames
  }
  return(object)
}

list_to_df <- function(listfordf){
  #https://gist.github.com/aammd/9ae2f5cce9afd799bafb?permalink_comment_id=1241571#gistcomment-1241571
  if(!is.list(listfordf)) stop("it should be a list")
  
  df <- list(keep = listfordf)
  class(df) <- c("tbl_df", "data.frame")
  attr(df, "row.names") <- .set_row_names(length(listfordf))
  
  if (!is.null(names(listfordf))) {
    df$grpid <- names(listfordf)
  }
  
  df
}


get_all_ms2_data <- function(pa, grp_idxs){
  
  #adapted from https://github.com/computational-metabolomics/msPurity/blob/ae3fcde2a7e67c6a78c0e8c21702401465b483a8/R/purityA-av-spectra.R#L386
  
  out_df = data.frame()
  
  for(grp_idx in grp_idxs){
    
    grped_info <- pa@grped_df[pa@grped_df==as.numeric(grp_idx),]
    grped_spectra <- pa@grped_ms2[as.character(grp_idx)][[1]]
    
    grped_info$index <- 1:nrow(grped_info)
    names(grped_spectra) <- 1:length(grped_spectra)
    
    ##############################################################################
    # Create a new dataframe with only the valid info and frag spectra
    ##############################################################################
    grped_spectra <- plyr::llply(grped_spectra, data.frame)
    df <- data.frame(do.call("rbind", grped_spectra))
    
    colnames(df)[1:2] <- c('mz', 'i')
    df$index <-   rep(seq_along(grped_spectra), sapply(grped_spectra, nrow))
    
    #if (!length(pa@filter_frag_params)==0){
    #  # if prior filtering performed only use those that have passed
    #  df<-df[df$pass_flag==1,]
    #}
    
    spectra_to_average <- merge(df, grped_info[, c('grpid', 'sample', 'cid', 'index', 'inPurity')], by = "index")
    if(nrow(out_df) == 0){
      out_df = spectra_to_average
    }else{
      out_df = rbind(out_df, spectra_to_average)  
    }
    
  }
  return(out_df)
}