##############################################
##############################################
##############################################
# script for analysing qPCR data
##############################################
##############################################
##############################################

# import libraries
library("plyr")
library("ggplot2")
library("optparse")
library("gtools")

# make options list
option_list <- list(
               make_option(c("-i", "--infile"),
                           help="provide infile of Ct values per sample per gene"),
               make_option(c("--housekeeper"),
                           help="specify housekeeping gene to normalise against"),
               make_option(c("--deltact-file"),
                           help="outfile to place 2^deltaCT"),
               make_option(c("--anova-file"),
                           help="outfile to place anova results"),
               make_option(c("--plot-order"),
                           help="order of conditions in plot")
               )

##############################
# get command line options
##############################

opt <- parse_args(OptionParser(option_list=option_list))

# read in sample and Ct information
dat <- read.csv(opt$`infile`, 
                header = T, 
                stringsAsFactors = F, 
                sep = "\t")

# change any "Undetermined" to NA
dat$Ct[dat$Ct == "Undetermined"] <- NA
dat$Ct <- as.numeric(dat$Ct)

# take mean of technical replicates
dat.ave <- ddply(dat, 
                 c("sample", "gene"), 
                 summarise, 
                 meanCt = mean(Ct, na.rm = T))

# housekeeping gene expression
hk <- dat.ave[dat.ave$gene == opt$`housekeeper`,]

# normalise to hk expression
deltacts = c()
for (i in 1:nrow(dat.ave)){
     sample <- dat.ave[i,]$sample
     deltact <- hk$meanCt[hk$sample == sample] - dat.ave[i,]$meanCt 
     deltacts <- append(deltacts, deltact)
     }

# add deltaCt to data frame
dat.ave$deltaCt <- deltacts

# add 2^deltaCt column
dat.ave$power2.deltaCt <- 2^dat.ave$deltaCt

# remove hk gene
dat.ave <- dat.ave[dat.ave$gene != opt$`housekeeper`,]

# write out table
write.table(dat.ave,
            file=opt$`deltact-file`,
            sep="\t",
            row.names=F,
	    quote=F)

# remove NAs
dat.ave <- na.omit(dat.ave)

# get the conditions
conds.a <- unlist(strsplit(dat.ave$sample, ".R"))
conds <- conds.a[seq(1, length(conds.a), 2)]
reps <- conds.a[seq(2, length(conds.a), 2)]
dat.ave$cond <- conds
dat.ave$rep <- reps

# plot 2^deltaCt
genes <- unique(dat.ave$gene)

if (!(is.null(opt$`plot-order`))){
   plot.order <- unlist(strsplit(opt$`plot-order`, ","))}else{
   plot.order <- mixedsort(unique(dat.ave$cond))
}

for (gene in genes){
    outname <- paste(gene, "pdf", sep = ".")
    res <- dat.ave[dat.ave$gene == gene,]
    plot1 <- ggplot(res, aes(x = factor(cond, levels = plot.order), 
                             y = power2.deltaCt, 
                             colour = factor(cond, levels = plot.order)))
    plot2 <- plot1 + geom_boxplot() + geom_point(position = position_jitter(width = 0.2))
    plot3 <- plot2 + ggtitle(gene)
    plot4 <- plot3 + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
    plot5 <- plot4 + xlab("Condition") + ylab("2^deltaCT")
    plot5 + scale_colour_discrete(name="Condition")
    ggsave(outname)
    }

# peform statistical testing
result <- matrix(ncol=3, nrow=length(genes))

for (i in 1:length(genes)){
    gene <- genes[i]
    deltact <- dat.ave$power2.deltaCt[dat.ave$gene == gene]
    conditions <- factor(dat.ave$cond[dat.ave$gene == gene])
    aov1 <- aov(deltact ~ conditions)
    s.aov1 <- summary(aov1)
    f = s.aov1[[1]]$F[1]
    p = s.aov1[[1]]$P[1]
    result[i, 1] <-gene
    result[i,2] <- f
    result[i,3] <- p
    }

colnames(result) <- c("gene", "F-statistic", "P-value")
result <- as.data.frame(result)
write.table(result,
            file=opt$`anova-file`,
            sep="\t",
            quote=F,
	    row.names = F)









