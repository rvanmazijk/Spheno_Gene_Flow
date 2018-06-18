# cytb database
dim(xx)

# spheno database
s = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/metadata/sphenomorphine_species.csv", stringsAsFactors = F)

s = s[complete.cases(s$sanger_sample), ]
s = s[s$cytB == FALSE, ]

# for (i in 1:nrow(s)) {
#   sample = s[i, "sanger_sample"]
#   sample = gsub('[A-Z]+', '', sample)
#   a = grep(sample, xx$SAMPLE_ID.x)
#   b = grep(sample, xx$ABTC_FIELD)
#   rows = c(a, b)
#   if (length(rows) > 0) {
#     cat(s[i, "OTU"], "\n")
#     cat(s[i, "sanger_sample"], "\n")
#     cat(xx[rows, "SAMPLE_ID.x"], "\n")
#     cat("***\n")
#   }
# }

d = read.csv("~/Desktop/cytb_add.csv", stringsAsFactors = F)
# spheno database
s = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/metadata/sphenomorphine_species.csv", stringsAsFactors = F)
s[s$sanger_sample %in% d$name, "cytB"] = TRUE
write.csv(s, "~/Dropbox/Sphenomorphine_Gene_Flow/data/metadata/sphenomorphine_species.csv", row.names = F)

sink("~/Desktop/cytb_to_add.fasta")
for (i in 1:nrow(d)) {
  cat('>', d[i, "name"], "\n", sep="")
  cat(xx[which(xx$SAMPLE_ID.x == d[i, "SAMPLE_ID"]), "SEQ"], "\n", sep="")
}
sink()