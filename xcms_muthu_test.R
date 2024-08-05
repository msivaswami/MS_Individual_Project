library(BiocParallel)
library(msPurity)
library(xcms)
library(MSnbase)
library(ggplot2)
library(openxlsx)
library(dplyr)
library(tidyr)

path_wd = getwd()
# path_wd = "/rds/projects/2015/viantm-01/users/weberrj/run-mspurity-metfrag"

# set path to auxiliary.R
source(file.path(path_wd, "workflows", "scripts", "auxiliary.R"))

set.seed(57475)
cores = 2

# Register cores
print(sprintf("%s system", .Platform$OS.type))
if (.Platform$OS.type == "unix") {
  register(bpstart(MulticoreParam(cores)))
} else {
  register(bpstart(SnowParam(cores)))
}

##################################
##################################
#### Get mzML paths

assay = "HILIC_POS"
suffix = ""
path_mzml_files = file.path(path_wd, "data", assay, "mzml")
path_results = file.path(path_wd, "results", "xcms_mspurity")

multiple_samples <- c("MS2_IL_70_140_F_QC", "MS2_IL_115_185_F_QC", "MS2_IL_150_220_F_QC",
                      "MS2_IL_200_310_F_QC", "MS2_IL_300_510_F_QC", "MS2_IL_500_1050_F_QC")

select_samples = multiple_samples

paths_mzml = list.files(path_mzml_files, pattern=glob2rx(paste(suffix, "*.mzML", sep="")), 
                        full.names = TRUE)
bname = sub(".mzML$", "", basename(paths_mzml))

# Check and match files with select_samples
matched_paths_mzml = paths_mzml[match(select_samples, bname)]
missing_samples = setdiff(select_samples, bname)

# Check for any missing samples
if (length(missing_samples) > 0) {
  stop(paste("Missing samples in the specified directory:", paste(missing_samples, collapse=", ")))
}

fn_rdata = file.path(path_results, paste(assay, ".RData", sep = ""))
fn_msPurity_msp = file.path(path_results, paste("mspurity_MS2_", assay, ".msp", sep=""))
fn_msPurity_grpeddf = file.path(path_results, paste("mspurity_grpeddf_", assay, ".xlsx", sep=""))
matched_fn_xlsx = file.path(path_results, paste('matched_', assay, ".xlsx", sep=''))
matched_fn_txt = file.path(path_results, paste('matched_', assay, ".txt", sep=''))

print("====================================================")
print(assay)
print(matched_paths_mzml)
print(cores)
print("====================================================")

##################################
##################################
#### Load cw parameters

cw_params = CentWaveParam(ppm = 10.0,
                          mzdiff = 0.001,
                          peakwidth = c(3.0,30.0),
                          snthresh = 10.0,
                          prefilter = c(3,100),
                          fitgauss = FALSE,
                          noise = 0,
                          verboseColumns = FALSE,
                          firstBaselineCheck = TRUE)
print("")
print(cw_params)
print("")

##################################
##################################
#### Set grouping parameters

# Set sample groups to match the length of matched_paths_mzml
sampleGroups = rep("sample", length(matched_paths_mzml))

pd_params = PeakDensityParam(sampleGroups = sampleGroups,
                             bw = 0.5,
                             minFraction = 0.01,
                             binSize = 0.005,
                             maxFeatures = 50)
print("")
print(sampleGroups)
print(pd_params)
print("")

##################################
## Processing steps
##################################

# Read data
print(paste("Reading .mzML files ....................."))
print(matched_paths_mzml)
raw_data <- readMSData(files = matched_paths_mzml, mode = "onDisk", msLevel = 1) # mslevel set to 1 and skip ms2
print(raw_data)

##################################
##################################
#### Find chromatographic peaks

print("Running findChromPeaks .....................")
xdata <- findChromPeaks(raw_data, cw_params)

print("Running groupChromPeaks .....................")
xdata = groupChromPeaks(xdata, pd_params)

# print("Writing tables to .xlsx .....................")
# write_tables_to_xlsx(xdata, fn_xlsx)

######################################################################
######################################################################

# tbls = getTables(xdata, naTOzero=FALSE) # variableMetadata, dataMatrix
# 
# wb <- createWorkbook()
# addWorksheet(wb, "variableMetaData")
# writeData(wb, "variableMetaData", tbls[[1]], colNames=TRUE, keepNA=FALSE)
# 
# addWorksheet(wb, "dataMatrix")
# writeData(wb, "dataMatrix", tbls[[2]], rowNames=TRUE, colNames=TRUE, keepNA=FALSE)
# 
# addWorksheet(wb, "metaData")
# writeData(wb, "metaData", phenoData(xdata)@data, colNames=TRUE, keepNA=FALSE)
# 
# saveWorkbook(wb, file = fn_xlsx, overwrite = TRUE)

pa <- purityA(matched_paths_mzml, interpol = "linear") %>%
  frag4feature(xcmsObj=xdata, createDb=FALSE) %>%
  filterFragSpectra(ilim=0,
                    plim=0.5,
                    ra=0,
                    snr=3.0,
                    rmp=FALSE,
                    snmeth='median'
  ) %>%  
  averageIntraFragSpectra(minfrac = 0.5,
                          minnum = 1,
                          ppm = 5.0,
                          snr = 0.0,
                          ra = 0.0,
                          av = 'median',
                          sumi = TRUE,
                          rmp = FALSE,
                          cores = cores
  ) %>%
  averageInterFragSpectra(minfrac = 0.5, 
                          minnum = 1, 
                          ppm = 5.0, 
                          snr = 0.0, 
                          ra = 0.0, 
                          av = 'median',
                          sumi = TRUE, 
                          rmp = FALSE, 
                          cores = cores)

# write msp to file
msPurity::createMSP(pa, fn_msPurity_msp, method="av_inter")

fdefs = data.frame(xcms::featureDefinitions(xdata)) %>% 
  dplyr::rename(mz = mzmed, rt = rtmed) %>%  
  groupNames_custom() %>%
  relocate(name, .before = mz)

# for each xcms feature for which MS2 data was collected, return information regarding precursor used to generate MS2
grped_df = pa@grped_df %>% 
  mutate(name = fdefs$name[pa@grped_df$grpid]) %>% 
  relocate(name, .before = grpid)

##################################
# msPurity to BEAMSpy
##################################

# pull out all (i.e. not average and not filtered) MS2 data associated with a given grpid (grpid = featureDefinitions(xdata) row) 
ms2Info_raw = pa %>% 
  get_all_ms2_data(grp_idxs = names(pa@grped_ms2)) %>% 
  dplyr::rename(scan = index) %>% 
  relocate(c(grpid, sample), .before = scan)

# update sample column to include the .mzml file name
ms2Info_raw$sample = unique(pa@grped_df$filename)[ms2Info_raw$sample]

# link ms2Info_raw to featureDefinitions(xcms) based on grpid and then add 'name' column as M###T### format
ms2Info_raw = ms2Info_raw %>% 
  mutate(name = fdefs$name[as.numeric(ms2Info_raw$grpid)]) %>% 
  relocate(name, .before = grpid)

# get averaged MS2 spectra
features_with_avrg_MS2 = pa@av_spectra %>% sapply(function(x) !is.null(x$av_inter))
grpids_with_avrg_MS = names(features_with_avrg_MS2[features_with_avrg_MS2 != FALSE])

ms2Info_avrg = lapply(grpids_with_avrg_MS, function(x, pa){
  ms2 = pa@av_spectra[[x]]$av_inter
  ms2$grpid = x
  return(ms2)
}, pa=pa)

ms2Info_avrg <- do.call('rbind', ms2Info_avrg) %>% 
  dplyr::relocate(grpid, .before = 'cl') 

ms2Info_avrg <- ms2Info_avrg %>% 
  mutate(name = fdefs[as.numeric(ms2Info_avrg$grpid), 'name']) %>% 
  dplyr::relocate(name, .before = 'grpid')

# write outputs to Excel
wb <- createWorkbook() %T>% 
  addWorksheet(paste("MS2precInfo", sep = '')) %T>%
  writeData(paste("MS2precInfo", sep = ''), grped_df, colNames=TRUE, keepNA=FALSE) %T>%
  addWorksheet(paste("rawMS2data", sep = '')) %T>%
  writeData(paste("rawMS2data", sep = ''), ms2Info_raw, colNames=TRUE, keepNA=FALSE) %T>%
  addWorksheet(paste("avrgMS2spectra", sep = '')) %T>%
  writeData(paste("avrgMS2spectra", sep = ''), ms2Info_avrg, colNames=TRUE, keepNA=FALSE) %>%
  saveWorkbook(file = matched_fn_xlsx, overwrite = TRUE)

# write BEAMSpy output to .txt file
write.table(ms2Info_raw, file = matched_fn_txt, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#save(list=ls(), file=fn_rdata)
save(xdata, grped_df, fdefs, ms2Info_avrg, file=fn_rdata)
