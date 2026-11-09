# Test runMethylTFR

library(methylTFR)

load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))
load(system.file("extdata", "example_data.rda", package = "methylTFR"))

# Write the bundled sites as bismarkCov files, one per sample
writeExampleSamples <- function(sample_dir, n = 2) {
    dir.create(sample_dir, showWarnings = FALSE)
    set.seed(1)
    files <- paste0("sample_", seq_len(n), ".cov")
    cov <- round(msites$coverage)
    keep <- is.finite(msites$score) & msites$score >= 0 &
        msites$score <= 1 & is.finite(cov) & cov >= 1
    m <- msites[keep]
    cov <- cov[keep]
    for (f in files) {
        n_meth <- stats::rbinom(length(m), cov, m$score)
        write.table(
            data.frame(
                chr = as.character(GenomicRanges::seqnames(m)),
                start = GenomicRanges::start(m),
                end = GenomicRanges::end(m),
                percent = 100 * n_meth / cov,
                meth = n_meth,
                unmeth = cov - n_meth
            ),
            file.path(sample_dir, f),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
        )
    }
    write.table(
        data.frame(sampleName = sub(".cov", "", files), bedFile = files),
        file.path(sample_dir, "samples.tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE
    )
    files
}

test_that("runMethylTFR validates its input", {
    expect_error(runMethylTFR(
        sample_ann = "sample_ann.tsv", sample_dir = ".",
        filetype = "invalid"
    ), "Please provide a valid file type")

    expect_error(runMethylTFR(
        sample_ann = "sample_ann.tsv", sample_dir = ".",
        filetype = NULL
    ), "Please provide a valid file type")

    expect_error(runMethylTFR(
        sample_ann = "sample_ann.tsv", sample_dir = ".",
        filetype = "epp", sampleColName = NULL
    ), "Please provide a valid sample column name")

    expect_error(runMethylTFR(
        sample_ann = "sample_ann.tsv", sample_dir = ".",
        filetype = "epp", full_path = "yes"
    ), "Invalid full path flag")

    expect_error(runMethylTFR(
        sample_ann = "sample_ann.tsv", sample_dir = ".",
        filetype = "epp"
    ), "Please load the annotation objects for given genome.")

    expect_error(runMethylTFR(
        sample_ann = "sample_ann.tsv", sample_dir = "no_such_dir",
        tf_bindsites = tf_bindsites, gcfreqs = gcfreqs, gc_dist = gcdist,
        filetype = "epp"
    ), "Sample directory does not exist")

    expect_error(runMethylTFR(
        sample_ann = "sample_ann.txt", sample_dir = ".",
        tf_bindsites = tf_bindsites, gcfreqs = gcfreqs, gc_dist = gcdist,
        filetype = "epp"
    ), "does not exist")
})

test_that("annotation inputs are checked against each other", {
    expect_error(runMethylTFR(
        sample_ann = "s.tsv", sample_dir = ".",
        tf_bindsites = tf_bindsites, gcfreqs = list(OTHER = gcfreqs$BATF),
        gc_dist = gcdist, filetype = "epp"
    ), "do not match")

    noBin <- gcdist
    S4Vectors::mcols(noBin) <- NULL
    expect_error(runMethylTFR(
        sample_ann = "s.tsv", sample_dir = ".",
        tf_bindsites = tf_bindsites, gcfreqs = gcfreqs,
        gc_dist = noBin, filetype = "epp"
    ), "GC_bin")

    expect_error(runMethylTFR(
        sample_ann = "s.tsv", sample_dir = ".",
        tf_bindsites = tf_bindsites, gcfreqs = unname(gcfreqs),
        gc_dist = gcdist, filetype = "epp"
    ), "gcfreqs must be named")

    expect_error(runMethylTFR(
        sample_ann = "s.tsv", sample_dir = ".",
        tf_bindsites = tf_bindsites, gcfreqs = gcfreqs,
        gc_dist = gcdist, filetype = "epp", enhancer = "chr1"
    ), "enhancer")
})

test_that("runMethylTFR runs end to end on files", {
    sample_dir <- tempfile("methylTFR_run")
    tmpDir <- tempfile("methylTFR_h5")
    on.exit(unlink(c(sample_dir, tmpDir), recursive = TRUE), add = TRUE)
    files <- writeExampleSamples(sample_dir, n = 3)

    devs <- runMethylTFR(
        sample_ann = "samples.tsv",
        sample_dir = sample_dir,
        tf_bindsites = tf_bindsites,
        gcfreqs = gcfreqs,
        gc_dist = gcdist,
        filetype = "bismarkCov",
        tempDir = tmpDir
    )
    expect_s4_class(devs, "methylTFRdeviations")
    expect_equal(dim(deviations(devs)), c(1L, 3L))
    expect_equal(colnames(deviations(devs)), files)
    expect_equal(dim(deviationZScores(devs)), c(1L, 3L))
    expect_true("sampleName" %in% colnames(SummarizedExperiment::colData(devs)))
    expect_length(list.files(tmpDir, pattern = "\\.h5$"), 0)

    # annfile takes a full path to the annotation file
    devs2 <- runMethylTFR(
        annfile = file.path(sample_dir, "samples.tsv"),
        sample_dir = sample_dir,
        tf_bindsites = tf_bindsites,
        gcfreqs = gcfreqs,
        gc_dist = gcdist,
        filetype = "bismarkCov",
        tempDir = tmpDir
    )
    expect_equal(deviations(devs), deviations(devs2))
})

test_that("invalid run options fall back to the defaults", {
    expect_warning(
        methylTFR:::check_run_options(chunkSize = "a"), "chunk size"
    )
    expect_warning(
        methylTFR:::check_run_options(threads = 0), "thread count"
    )
    expect_warning(
        methylTFR:::check_run_options(ignoreStrand = "x"), "strand option"
    )
    expect_warning(
        methylTFR:::check_run_options(cov_threshold = -1), "cov_threshold"
    )
    opts <- suppressWarnings(methylTFR:::check_run_options(
        chunkSize = "a", threads = 0, ignoreStrand = "x", cov_threshold = -1
    ))
    expect_equal(opts$chunkSize, 20)
    expect_equal(opts$threads, 1)
    expect_true(opts$ignoreStrand)
    expect_equal(opts$cov_threshold, 1)
})

test_that("bpparam_from_threads picks a serial back-end for one thread", {
    expect_s4_class(
        methylTFR:::bpparam_from_threads(1), "SerialParam"
    )
    expect_s4_class(
        methylTFR:::bpparam_from_threads(NULL), "SerialParam"
    )
})

test_that("valid_core_motifs drops motifs without sites", {
    gf <- list(BATF = gcfreqs$BATF, EMPTY = gcfreqs$BATF)
    tb <- c(as.list(tf_bindsites), list(EMPTY = GenomicRanges::GRanges()))
    expect_equal(suppressMessages(methylTFR:::valid_core_motifs(tb, gf)), "BATF")
    expect_error(
        suppressMessages(methylTFR:::valid_core_motifs(
            list(EMPTY = GenomicRanges::GRanges()),
            list(EMPTY = gcfreqs$BATF)
        )),
        "No valid motifs"
    )
})
