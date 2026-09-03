# =============================================================================
# 10_load_reference_metric.R -- Loads the reference implementation AND the eleven data sets,
#              without implementing the metric a second time and without
#              altering 12_metric_reference.R.
# =============================================================================
# The S3 brief, section 0.3: "do not rewrite it, call it."
#
# THE PROCEDURE, and why it has to be this way:
# 12_metric_reference.R defines at the top everything that is needed (the gene sets,
# kern(), mk_zieh(), mk_zieh_L(), einzel_f(), kontrast_f(), the protein lengths
# AS), then loads the eleven data sets below and calls auswerten() for each --
# which takes about 25 minutes and is not needed here.
#
# Instead of copying the loading block out (that WOULD be a second
# implementation and could quietly drift apart), the file is parsed and
# evaluated in two sections:
#   (1) every expression BEFORE the line "# Die elf Datensaetze" -- the
#       definitions. After that the computational core is in memory.
#   (2) auswerten() is replaced by a stub that collects its arguments and
#       returns NULL.
#   (3) every expression between that line and the line "TAB <- do.call" --
#       the loading block. It calls the stub and therefore computes nothing.
# The result: DATEN, a list holding Z, meta, expr, label, arm, rolle and
# klasse for all eleven points, produced by exactly the code that also
# produced §8.5.
#
# The line numbers are not hard-wired but found through the text marks -- if
# the file shifts, this does not break silently.
# =============================================================================

# --- The path parameter: the session tree ---------------------------------
# This script used to point at a directory that survived only as a Windows
# directory junction onto a backup -- a junction cannot be cloned or
# versioned, and its target was documented nowhere. It is replaced by an
# explicit, overridable parameter: SESSIONS points at the tree of session
# directories in which the reference implementation of the metric and the
# loading block live. That tree holds the raw analysis sessions and is not
# part of the public archive; set PAPER_V2_SESSIONS to point at it.
SESSIONS <- Sys.getenv("PAPER_V2_SESSIONS")
if (!nzchar(SESSIONS))
  stop("Set PAPER_V2_SESSIONS to the tree of raw analysis sessions. ",
       "That tree is not part of the public archive; see README.md.")
stopifnot(dir.exists(SESSIONS))
FP   <- file.path(SESSIONS, "03_Metrik_Elf_Punkte",
                  "reference_implementations", "12_metric_reference.R")

zeilen <- readLines(FP)
mark_a <- grep("^# Die elf Datensaetze", zeilen)
mark_b <- grep("^TAB <- do\\.call", zeilen)
stopifnot(length(mark_a) == 1, length(mark_b) == 1, mark_a < mark_b)

EX  <- parse(FP, keep.source = TRUE)
srf <- utils::getSrcref(EX)
start <- vapply(srf, function(s) as.integer(s)[1], integer(1))

cat(sprintf("10_load_reference_metric.R: %d expressions; definitions < line %d, loading block %d..%d\n",
            length(EX), mark_a, mark_a, mark_b))

# --- (1) The definitions ---------------------------------------------------
for (i in which(start < mark_a)) eval(EX[[i]], envir = globalenv())
stopifnot(is.function(kern), is.function(einzel_f), is.function(auswerten),
          length(S_NEUTRAL) > 100, length(AS) > 10000)

# --- (2) The stub ----------------------------------------------------------
auswerten_echt <- auswerten
DATEN <- list()
auswerten <- function(Z, meta, expr, label, arm, rolle, klasse) {
  DATEN[[label]] <<- list(Z = Z, meta = meta, expr = expr, label = label,
                          arm = arm, rolle = rolle, klasse = klasse)
  cat(sprintf("  loaded: %-30s %-11s class %s  (%d genes x %d samples)\n",
              label, arm, klasse, nrow(Z), ncol(Z)))
  NULL
}

# --- (3) The loading block -------------------------------------------------
for (i in which(start >= mark_a & start < mark_b)) eval(EX[[i]], envir = globalenv())

auswerten <- auswerten_echt
cat(sprintf("10_load_reference_metric.R: %d data sets loaded.\n", length(DATEN)))
stopifnot(length(DATEN) == 11)
