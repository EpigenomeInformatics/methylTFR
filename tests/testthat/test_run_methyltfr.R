# Test runMethylTFR

library(methylTFR)

load(system.file("extdata", "gcdist_subset.rda", package = "methylTFR"))
load(system.file("extdata", "BATF_gcfreqs.rda", package = "methylTFR"))
load(system.file("extdata", "BATF_tf_bindsites.rda", package = "methylTFR"))
load(system.file("extdata", "example_data.rda", package = "methylTFR"))

# Write the bundled sites as bismarkCov files, one per sample
writeExampleSamples <- function(sampleDir, n = 2) {
    dir.create(sampleDir, showWarnings = FALSE)
    set.seed(1)
    files <- paste0("sample_", seq_len(n), ".cov")
    for (f in files) {
        nMeth <- stats::rbinom(length(msites), msites$coverage, msites$score)
        write.table(
            data.frame(
                chr = as.character(GenomicRanges::seqnames(msites)),
                start = GenomicRanges::start(msites),
                end = GenomicRanges::end(msites),
                percent = 100 * nMeth / msites$coverage,
                meth = nMeth,
                unmeth = msites$coverage - nMeth
            ),
            file.path(sampleDir, f),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
        )
    }
    write.table(
        data.frame(sampleName = sub(".cov", "", files), bedFile = files),
        file.path(sampleDir, "samples.tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE
    )
    files
}

test_that("runMethylTFR validates its input", {
    expect_error(runMethylTFR(
        sampleAnn = "sample_ann.tsv", sampleDir = ".",
        filetype = "invalid"
    ), "Please provide a valid file type")

    expect_error(runMethylTFR(
        sampleAnn = "sample_ann.tsv", sampleDir = ".",
        filetype = NULL
    ), "Please provide a valid file type")

    expect_error(runMethylTFR(
        sampleAnn = "sample_ann.tsv", sampleDir = ".",
        filetype = "epp", sampleColName = NULL
    ), "Please provide a valid sample column name")

    expect_error(runMethylTFR(
        sampleAnn = "sample_ann.tsv", sampleDir = ".",
        filetype = "epp", fullPath = "yes"
    ), "Invalid full path flag")

    expect_error(runMethylTFR(
        sampleAnn = "sample_ann.tsv", sampleDir = ".",
        filetype = "epp"
    ), "Please load the annotation objects for given genome.")

    expect_error(runMethylTFR(
        sampleAnn = "sample_ann.tsv", sampleDir = "no_such_dir",
        tfBindsites = tf_bindsites, gcfreqs = gcfreqs, gcDist = gcdist,
        filetype = "epp"
    ), "Sample directory does not exist")

    expect_error(runMethylTFR(
        sampleAnn = "sample_ann.txt", sampleDir = ".",
        tfBindsites = tf_bindsites, gcfreqs = gcfreqs, gcDist = gcdist,
        filetype = "epp"
    ), "does not exist")
})

test_that("annotation inputs are checked against each other", {
    expect_error(runMethylTFR(
        sampleAnn = "s.tsv", sampleDir = ".",
        tfBindsites = tf_bindsites, gcfreqs = list(OTHER = gcfreqs$BATF),
        gcDist = gcdist, filetype = "epp"
    ), "do not match")

    noBin <- gcdist
    S4Vectors::mcols(noBin) <- NULL
    expect_error(runMethylTFR(
        sampleAnn = "s.tsv", sampleDir = ".",
        tfBindsites = tf_bindsites, gcfreqs = gcfreqs,
        gcDist = noBin, filetype = "epp"
    ), "GC_bin")

    expect_error(runMethylTFR(
        sampleAnn = "s.tsv", sampleDir = ".",
        tfBindsites = tf_bindsites, gcfreqs = unname(gcfreqs),
        gcDist = gcdist, filetype = "epp"
    ), "gcfreqs must be named")

    expect_error(runMethylTFR(
        sampleAnn = "s.tsv", sampleDir = ".",
        tfBindsites = tf_bindsites, gcfreqs = gcfreqs,
        gcDist = gcdist, filetype = "epp", enhancer = "chr1"
    ), "enhancer")
})

test_that("runMethylTFR runs end to end on files", {
    sampleDir <- tempfile("methylTFR_run")
    tmpDir <- tempfile("methylTFR_h5")
    on.exit(unlink(c(sampleDir, tmpDir), recursive = TRUE), add = TRUE)
    files <- writeExampleSamples(sampleDir, n = 3)

    devs <- runMethylTFR(
        sampleAnn = "samples.tsv",
        sampleDir = sampleDir,
        tfBindsites = tf_bindsites,
        gcfreqs = gcfreqs,
        gcDist = gcdist,
        filetype = "bismarkCov",
        tempDir = tmpDir
    )
    expect_s4_class(devs, "methylTFRdeviations")
    expect_equal(dim(deviations(devs)), c(1L, 3L))
    expect_equal(colnames(deviations(devs)), files)
    expect_equal(dim(deviationZScores(devs)), c(1L, 3L))
    expect_true("sampleName" %in% colnames(SummarizedExperiment::colData(devs)))
    expect_length(list.files(tmpDir, pattern = "\\.h5$"), 0)

    # A full path to the annotation file works as well
    devs2 <- runMethylTFR(
        sampleAnn = file.path(sampleDir, "samples.tsv"),
        sampleDir = sampleDir,
        tfBindsites = tf_bindsites,
        gcfreqs = gcfreqs,
        gcDist = gcdist,
        filetype = "bismarkCov",
        tempDir = tmpDir
    )
    expect_equal(deviations(devs), deviations(devs2))
})

test_that("invalid run options fall back to the defaults", {
    expect_warning(
        methylTFR:::checkRunOptions(chunkSize = "a"), "chunk size"
    )
    expect_warning(
        methylTFR:::checkRunOptions(threads = 0), "thread count"
    )
    expect_warning(
        methylTFR:::checkRunOptions(ignoreStrand = "x"), "strand option"
    )
    expect_warning(
        methylTFR:::checkRunOptions(covThreshold = -1), "covThreshold"
    )
    opts <- suppressWarnings(methylTFR:::checkRunOptions(
        chunkSize = "a", threads = 0, ignoreStrand = "x", covThreshold = -1
    ))
    expect_equal(opts$chunkSize, 20)
    expect_equal(opts$threads, 1)
    expect_true(opts$ignoreStrand)
    expect_equal(opts$covThreshold, 1)
})

test_that("bpparamFromThreads picks a serial back-end for one thread", {
    expect_s4_class(
        methylTFR:::bpparamFromThreads(1), "SerialParam"
    )
    expect_s4_class(
        methylTFR:::bpparamFromThreads(NULL), "SerialParam"
    )
})

test_that("validCoreMotifs drops motifs without sites", {
    gf <- list(BATF = gcfreqs$BATF, EMPTY = gcfreqs$BATF)
    tb <- c(as.list(tf_bindsites), list(EMPTY = GenomicRanges::GRanges()))
    expect_equal(suppressMessages(methylTFR:::validCoreMotifs(tb, gf)), "BATF")
    expect_error(
        suppressMessages(methylTFR:::validCoreMotifs(
            list(EMPTY = GenomicRanges::GRanges()),
            list(EMPTY = gcfreqs$BATF)
        )),
        "No valid motifs"
    )
})
