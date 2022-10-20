library(DESeq2)

## good example ##
# https://lashlock.github.io/compbio/R_presentation.html

# get the counts file
args = commandArgs(trailingOnly=TRUE)
# check to see that two argument is given
if (length(args)!=2) {
  stop("Two arguments must be supplied (input file).n", call.=FALSE)
}

infilename = args[1]
outfilename = args[2]

# read count file
tabla <- read.table(infilename, sep=",", row.names=1, header=TRUE)
# group columns for the design matrix TODO: automate this by checking the number of columns
group <- factor(rep(c("A","B"),each=3))
# design matrix or treatment data construction
design_df <- data.frame(
    id=colnames(tabla), 
    ko=group
    )
# construct deseq2 object
dds <- DESeqDataSetFromMatrix(countData = tabla,
                              colData = design_df,
                              design = ~ ko)
# run DESeq pipeline, will normalize using median of ratios
dds <- DESeq(dds)
# get the results
res <- results(dds)
# save to file .. 
write.table(res, file=outfilename, sep=",", row.names=TRUE, col.names=TRUE)
