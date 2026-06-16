library(DBI)
library(RSQLite)

#Rename beaked whale species labels used in Pamguard databases according to new 
#NOAA standardized naming 

# ── Configuration ─────────────────────────────────────────────────────────────
db_dir       <- "path/to/your/folder"   # <-- set this
table_name   <- "YourTable"             # <-- table containing the variable
column_name  <- "YourColumn"            # <-- column to correct
old_value    <- "wrong_string"          # <-- value to find
new_value    <- "correct_string"        # <-- replacement value
# ──────────────────────────────────────────────────────────────────────────────

db_files <- list.files(
  path       = db_dir,
  pattern    = "\\.sqlite3$",
  full.names = TRUE
)

if (length(db_files) == 0) stop("No .sqlite3 files found in: ", db_dir)

message("Found ", length(db_files), " database(s). Running corrections...")

results <- lapply(db_files, function(f) {
  con <- dbConnect(SQLite(), f)
  on.exit(dbDisconnect(con))  # ensures connection closes even if error occurs
  
  # Count matching rows before update
  n_before <- dbGetQuery(con,
                         sprintf("SELECT COUNT(*) AS n FROM \"%s\" WHERE \"%s\" = ?",
                                 table_name, column_name),
                         params = list(old_value)
  )$n
  
  if (n_before > 0) {
    dbExecute(con,
              sprintf("UPDATE \"%s\" SET \"%s\" = ? WHERE \"%s\" = ?",
                      table_name, column_name, column_name),
              params = list(new_value, old_value)
    )
    message(sprintf("  %-50s  %d row(s) updated", basename(f), n_before))
  } else {
    message(sprintf("  %-50s  no match, skipped", basename(f)))
  }
  
  data.frame(file = basename(f), rows_updated = n_before)
})

summary_df <- do.call(rbind, results)
message("\nTotal rows updated across all databases: ", sum(summary_df$rows_updated))