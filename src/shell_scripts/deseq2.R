library(DESeq2)

## good example ##
# https://lashlock.github.io/compbio/R_presentation.html

# get the counts file
args = commandArgs(trailingOnly=TRUE)
# check to see that two argument is given
if (length(args)!=4) {
  stop("Fours arguments must be supplied (input file).n", call.=FALSE)
}

infilename = args[1]
outfilename = args[2]
input_reps = args[3]
output_reps = args[4]

# read count file
tabla <- read.table(infilename, sep=",", row.names=1, header=TRUE)
# group columns for the design matrix TODO: automate this by checking the number of columns
group <- factor(rep(c("A","B"), times=c(as.numeric(input_reps), as.numeric(output_reps))))
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
