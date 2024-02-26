The human genome is composed of four "building blocks", A, T, C and G. By chance, each should be present 25% each. In total we have 3,100 Nillion of these.

In reality we have a bias distribution of CG being in fewer numbers: C 19%, G 19%, A 28%, T 28%

<details>
  <summary>Show code</summary>

```{.r code-line-numbers="|1-3|5-8|5-13|15-16|18-26|28-31|"}
# Import the human genome sequence, version hg38
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38

# Only use annotated chromosomes (drop those not aligned properly)
chrom <- rownames(as.data.frame(seqinfo(genome)))
annotated_chrom <- grepl("_", chrom)
annotated_chrom <- chrom[!annotated_chrom]

#> annotated_chrom
# [1] "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9" 
#[10] "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18"
#[19] "chr19" "chr20" "chr21" "chr22" "chrX"  "chrY"  "chrM" 

# Choose pattern to match
nucleotides <- c("A", "C", "G", "T", "N")

# Count the pattern per chromosome, save to a list
count <- list()
for (chrom in annotated_chrom){
	print(paste("counting chrom:", chrom))
	
	count[[chrom]] <- sapply(nucleotides, function(nucleotide){
		countPattern(nucleotide, genome[[chrom]])
	})
}

# Convert list into a nice looking table
dat <- Reduce(function(x,y) rbind(x,y), count)
dat <- data.frame(dat)
rownames(dat) <- annotated_chrom
```

```{.r code-line-numbers=""}
# Calculate the column sums
x <- colSums(dat)
x <- round( x / sum(x), 2 )

labels <- paste(x * 100, "%" , names(x))
pie(x, labels = labels)
```

</details>

![](img/piechart_base_dist.png){width="500" fig-align="left"}

The distribution is NOT 25% each, AT > CG
