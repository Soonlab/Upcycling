#!/usr/bin/env Rscript
# A7: Predict minimum doubling time for 111 MAGs via gRodon2
# Uses ribosomal protein codon usage bias
suppressMessages({
    library(gRodon)
    library(Biostrings)
})

bakta <- "/data/data/Upcycling/MAGs_FASTA_files/bakta_results"
outdir <- "/data/data/Upcycling/research/additional/A7_grodon"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

hero <- c("S13","S16","S23","C22","M1","S26")
mags <- list.dirs(bakta, recursive = FALSE, full.names = FALSE)
cat("[A7] processing", length(mags), "MAGs\n")

results <- data.frame(MAG=character(), group=character(),
                      CUBHE=numeric(), ConsistencyHE=numeric(),
                      CPB=numeric(), FilteredSequences=integer(),
                      nHE=integer(), nGenes=integer(),
                      d_hours=numeric(),
                      LowerCI=numeric(), UpperCI=numeric(),
                      stringsAsFactors=FALSE)

for (m in mags) {
    ffn <- file.path(bakta, m, paste0(m, ".ffn"))   # nucleotide CDS
    tsv <- file.path(bakta, m, paste0(m, ".tsv"))   # annotations
    if (!file.exists(ffn) | !file.exists(tsv)) { next }
    genes <- tryCatch(readDNAStringSet(ffn), error=function(e) NULL)
    if (is.null(genes)) next
    ann <- read.table(tsv, sep="\t", comment.char="#", header=FALSE,
                      col.names=c("contig","type","start","end","strand",
                                  "locus","gene","product","dbxrefs"),
                      quote="", fill=TRUE)
    # highly expressed = ribosomal protein large/small subunit
    # Match either "ribosomal ... protein" or short gene names rplA / rpsB / rpmC
    # Exclude rRNA methyltransferases, ribosome-binding, etc.
    he_mask <- (grepl("ribosomal.*protein|ribosomal subunit protein|[0-9]{2}S ribosomal",
                      ann$product, ignore.case=TRUE) &
                !grepl("methyltransferase|kinase|acetyltransferase|acetylase|rimI|rimM|rimO",
                       ann$product, ignore.case=TRUE)) |
               grepl("^(rpl|rps|rpm)[A-Z0-9]", ann$gene, ignore.case=TRUE)
    he_loci <- ann$locus[he_mask]
    gene_names <- sub(" .*", "", names(genes))
    he_logi <- gene_names %in% he_loci
    if (sum(he_logi) < 10) {
        cat("[A7]", m, "- <10 ribosomal genes found, skipping\n")
        next
    }
    # gRodon: predictGrowth(genes, highly_expressed, temperature=20, mode="partial")
    result <- tryCatch(
        predictGrowth(genes, he_logi, mode="partial"),
        error=function(e) { cat("[A7]", m, "ERROR:", conditionMessage(e), "\n"); NULL }
    )
    if (is.null(result)) next
    grp <- ifelse(m %in% hero, "MICP_complete", "rest")
    results <- rbind(results, data.frame(
        MAG=m, group=grp,
        CUBHE=result$CUBHE, ConsistencyHE=result$ConsistencyHE,
        CPB=result$CPB,
        FilteredSequences=ifelse(!is.null(result$FilteredSequences), result$FilteredSequences, NA),
        nHE=sum(he_logi), nGenes=length(genes),
        d_hours=result$d,
        LowerCI=result$LowerCI, UpperCI=result$UpperCI
    ))
}

write.csv(results, file.path(outdir, "gRodon_growth_rates_per_MAG.csv"), row.names=FALSE)

cat("\n[A7] growth rate summary:\n")
cat(sprintf("  MICP-complete (n=%d): median %.2f h (IQR %.2f-%.2f)\n",
    sum(results$group=="MICP_complete"),
    median(results$d_hours[results$group=="MICP_complete"], na.rm=TRUE),
    quantile(results$d_hours[results$group=="MICP_complete"], 0.25, na.rm=TRUE),
    quantile(results$d_hours[results$group=="MICP_complete"], 0.75, na.rm=TRUE)))
cat(sprintf("  Rest          (n=%d): median %.2f h (IQR %.2f-%.2f)\n",
    sum(results$group=="rest"),
    median(results$d_hours[results$group=="rest"], na.rm=TRUE),
    quantile(results$d_hours[results$group=="rest"], 0.25, na.rm=TRUE),
    quantile(results$d_hours[results$group=="rest"], 0.75, na.rm=TRUE)))

mwu <- wilcox.test(
    results$d_hours[results$group=="MICP_complete"],
    results$d_hours[results$group=="rest"]
)
cat(sprintf("  Mann-Whitney p = %.4f\n", mwu$p.value))

cat("\n[A7] per-MAG MICP-complete:\n")
print(results[results$group=="MICP_complete",
              c("MAG","d_hours","LowerCI","UpperCI","nHE","ConsistencyHE")])
cat("[A7] DONE\n")
