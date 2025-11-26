# Create simulated data formats for the methylTFR tests
output_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR/inst/extdata/"

# Encode data
Encode <- data.frame(
  chr = c("chr1", "chr1", "chr1", "chr1", "chr1", "chr1"),
  start = c(1000170, 1000190, 1000191, 1000198, 1000199, 1000206),
  end = c(1000171, 1000191, 1000192, 1000199, 1000200, 1000207),
  name = rep("HepG2_B1__GC_", 6),
  score = c(62, 62, 31, 62, 31, 31),
  strand = c("+", "+", "-", "+", "-", "-"),
  thickStart = c(1000170, 1000190, 1000191, 1000198, 1000199, 1000206),
  thickEnd = c(1000171, 1000191, 1000192, 1000199, 1000200, 1000207),
  itemRgb = rep("255,255,0", 6),
  blockCount = c(62, 62, 31, 62, 31, 31),
  blockSizes = c(6, 3, 0, 10, 0, 10)
)

# bismarkCytosine data
bismarkCytosine <- data.frame(
  chr = rep("chr2", 7),
  #start = rep(2, 7),
  end = c(16050097, 16050098, 16050114, 16050115, 16115591, 16117938, 16122790),
  strand = c("+", "-", "+", "-", "+", "-", "+"),
  score1 = c(0, 1, 0, 0, 1, 0, 0),
  score2 = c(3, 4, 4, 0, 1, 2, 1),
  context = rep("CG", 7),
  trinucleotide = c("CGG", "CGA", "CGG", "CGT", "CGC", "CGT", "CGC")
)

# bissnp data
bissnp <- data.frame(
  chr = c("chr1", "chr1", "chr1", "chr1"),
  start = c(10496, 10524, 864802, 864803),
  end = c(10497, 10525, 864803, 864804),
  score = c(79.69, 90.62, 58.70, 50.00),
  coverage = c(64, 64, 46, 4),
  strand = c("+", "+", "+", "-"),
  thickStart = c(10496, 10524, 864802, 864803),
  thickEnd = c(10497, 10525, 864803, 864804),
  itemRgb = c("180,60,0", "210,0,0", "120,120,0", "90,150,0"),
  blockCount = c(64, 64, 46, 4),
  blockSizes = c(0, 0, 5, 145)
)

# EPP data
epp <- data.frame(
  chr = c("chr1", "chr1", "chr1", "chr1", "chr1", "chr1"),
  start = c(3010957, 3010959, 3010971, 3010973, 3011025, 3011027),
  end = c(3010958, 3010960, 3010972, 3010974, 3011026, 3011028),
  meth = c("27/27", "2/7", "10/20", "1/20", "57/70", "81/100"),
  score = c(1000, 500, 1000, 500, 814, 500),
  strand = c("+", "-", "+", "-", "+", "-")
)

# bismarkCov
bismarkCov <- data.frame(
  chr = c("chr1", "chr1", "chr1", "chr1", "chr1", "chr1"),
  start = c(73252, 73253, 73256, 73260, 73262, 73269),
  end = c(73253, 73254, 73257, 73261, 73263, 73270),
  meth = c(100, 0, 100, 0, 100, 100),
  mcov = c(1, 0, 1, 0, 1, 1),
  ucov = c(0, 1, 0, 1, 0, 0)
)

# allc
allc <- data.frame(
  chromosome = c("chr12", "chr12", "chr12"),
  position = c(18283342, 18283343, 18283344),
  strand = c("+", "-", "+"),
  sequence_context = c("CGT", "CGA", "CGG"),
  mc = c(1, 2, 1),
  cov = c(2, 3, 2),
  methylated = c(1, 0, 1)
)

# Save as a compressed TSV files
write.table(Encode, file = gzfile(paste0(output_dir, "encode.tsv.gz")), sep = "\t", row.names = FALSE, quote = FALSE)

write.table(bismarkCytosine, file = gzfile(paste0(output_dir, "bismarkCytosine.tsv.gz")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

write.table(bissnp, file = gzfile(paste0(output_dir, "bissnp.tsv.gz")), sep = "\t", row.names = FALSE, quote = FALSE)

write.table(epp, file = gzfile(paste0(output_dir, "epp.tsv.gz")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

write.table(bismarkCov, file = gzfile(paste0(output_dir, "bismarkCov.tsv.gz")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

write.table(allc, file = gzfile(paste0(output_dir, "allc.tsv.gz")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
