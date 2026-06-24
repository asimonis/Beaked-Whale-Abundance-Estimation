#Generate a single acoustic study for each drift and then combine them all into 
#single acoustic study


library(beepr)
library(PAMpal)

db <- 'F:/Beaked_whale_data/Databases/CalCurCEAS_027.sqlite3'
bin <- 'F:/Beaked_whale_data/Binaries/CalCurCEAS_027'
pps <- PAMpalSettings(db=db, binaries=bin, sr_hz='auto', winLen_sec=.0025,
                      filterfrom_khz=10, filterto_khz=NULL)

data <- processPgDetections(pps, mode='db', id='RoboJStudy')
data <- setSpecies(data, method='pamguard')

data@ancillary[["warnings"]][["message"]]
beep(1)

saveRDS(data,file='F:/Beaked_whale_data/AcousticStudy/CalCurCEAS027_EventsOnly.rds')


#####################


# ── Configuration ─────────────────────────────────────────────────────────────
rds_dir      <- "F:/Beaked_whale_data/AcousticStudy"   # <-- set this
output_path  <- file.path(rds_dir, "CalCurCEAS_Combined_EventsOnly.rds")
# ──────────────────────────────────────────────────────────────────────────────

# Discover files
rds_files <- list.files(
  path       = rds_dir,
  pattern    = "^CalCurCEAS.*EventsOnly\\.rds$",
  full.names = TRUE
)

if (length(rds_files) == 0) {
  stop("No matching .rds files found in: ", rds_dir)
}

message("Found ", length(rds_files), " AcousticStudy file(s):")
message(paste0("  ", basename(rds_files), collapse = "\n"))

# Load all AcousticStudy objects
study_list <- lapply(rds_files, function(f) {
  message("Loading: ", basename(f))
  readRDS(f)
})

# Combine using PAMpal::bindStudies
message("\nCombining studies...")
combined <- bindStudies(study_list)


# Save
message("Saving combined study to: ", output_path)
saveRDS(combined, file = output_path)

message("Done. Combined study contains ", 
        length(events(combined)), " event(s).")

