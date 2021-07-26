
library(DESeq2)

TnSeqDataSet <- function(metadata, name, qpath, mpath, mincount=10, sample=''){
    cat("Building dataset for set ", name, "\n", sep="")
    # Slice metadata to relevant sample
    metadata <- metadata[metadata$Inoculum==name,]

    # Import quantified data from qpath
    countdata <- read_quantified_files(paste(qpath, paste(sample, rownames(metadata), ".count.txt", sep=""), sep="/"), mincount)

    # Import mapped library
    mapdata <- read_mapping_file(paste(mpath, metadata$Library[1], "barcode_map.txt", sep="/"))

    # Return combined dataset
    dataset = list(name=name, meta=metadata, counts=countdata, map=mapdata)

    return(dataset)
}

read_beat_mapping_file <- function(file){
    map <- read.table(file,stringsAsFactors=F)
    libs <- split(map,map$V4)
    for(i in 1:length(libs)){
        lib <- libs[[i]]
        rownames(lib) <- lib[,5]
        lib <- lib[,1:4]
        colnames(lib) <- c("Element","Position","Strand","Library")
        libs[[i]] <- lib
    }
    return(libs)
}

OLD_read_mapping_file <- function(file){
    map <- read.table(file, stringsAsFactors=F)
    colnames(map) <- c("Barcode", "Count", "Position", "Locus", "Strand", "Count2", "Feature", "ShortName")
    return(map)
}

merge_by_rowname <- function(df1,df2){
    new <- merge(df1,df2,by=0,all=T)
    rownames(new) <- new[,1]
    new <- new[,2:ncol(new)]
    cat(paste("Merged:",ncol(new),"x",nrow(new),"\n"))
    return(new)
}

merge_list <- function(dflist){
    new <- dflist[[1]]
    for(i in 2:length(dflist)){
        new <- merge_by_rowname(new,dflist[[i]])
    }
    colnames(new) <- names(dflist)
    return(new)
}

read_quantified_file <- function(file){
    df <- read.table(file,fill=T,stringsAsFactors=F)
    return(df)
}

read_quantified_files <- function(files,min_count){
    fulldfs <- list()
    dfs <- list()
    for(file in files){
        name <- strsplit(file,"/")[[1]]
        name <- name[length(name)]
        name <- sub(".txt", "", name)
        fulldfs[[name]] <- read_quantified_file(file)
        fulldfs[[name]] <- fulldfs[[name]][fulldfs[[name]][,2]>=min_count,]
        dfs[[name]] <- fulldfs[[name]][,2,F]
        rownames(dfs[[name]]) <- fulldfs[[name]][,1]
        cat(paste(file,"-",nrow(dfs[[name]]),"\n"))
    }
    qdf <- merge_list(dfs)
    qdf[is.na(qdf)] <- 0
    return(qdf)
}

read_metadata <- function(file){
    meta <- read.table(file,sep="\t",stringsAsFactors=F,row.names=1)
    colnames(meta) <- c("Library","TVNo","Inoculum","Mouse","Day","Locus")
    return(meta)
}

calculate_fitness <- function(dataset){
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=dataset$counts, colData=dataset$meta, design=~Day)
    dds <- DESeq2::DESeq(dds)
    fitness <- list()
    for(day in levels(dds$Day)[-1]){
        fitness[[day]] = as.data.frame(DESeq2::results(dds,contrast=c("Day",day,"d0")))
    }
    dataset$fitness <- fitness
    return(dataset)
}

read_mapping_file <- function(file){
    df <- read.table(file,stringsAsFactors=F,row.names=1)
    df <- df[,c(2,3,4,6,7)]
    colnames(df) <- c("Position","Element","Strand","Feature", "ShortName")
    return(df)
}

build_output <- function(res,map,control_tags){
    barcodes <- rownames(res$foldChanges)
    output <- map[barcodes,]
    rownames(output) <- barcodes
    output <- cbind(output,res$foldChanges,res$pvalues,res$lfcSEs)

    median_controls <- lapply(res,function(x) apply(x[control_tags,],2,median,na.rm=T))
    zscores <- sapply(1:ncol(res$foldChanges),function(x) calculate_2_dist_zscore(res$foldChanges[,x],res$lfcSEs[,x],median_controls$foldChanges[x],median_controls$lfcSEs[x]))
    output <- cbind(output,zscores)

    colnames(output) <- c(colnames(map),paste(colnames(res$foldChanges),"_log2FC",sep=""),paste(colnames(res$pvalues),"_padj",sep=""),paste(colnames(res$lfcSEs),"_lfcSE",sep=""),paste(colnames(res$foldChanges),"_zscore",sep=""))

    return(output)
}

map_gene <- function(gene,gffs){
    gff <- gffs[[gene[1]]]
    pos <- gene[2]
    hits <- gff[(gff$V4 <= as.numeric(pos)) & (as.numeric(pos) <= gff$V5),c("cds","gene")]
    cds <- paste(hits$cds,collapse=";")
    name <- paste(hits$gene,collapse=";")
    gene <- c(gene,CDS=cds,Gene=name)
    return(gene)
}

map_genes <- function(lib,gff_files){
    gffs <- list()
    for(file in gff_files){
        name <- sub(".*/","",sub(".gff","",file))
        gff <- read.table(file,sep="\t",stringsAsFactors=F)
        gff <- gff[gff$V3=="CDS",]
        gff$cds <- sub(";","",substr(regmatches(gff$V9,regexec("ID=.*?;",gff$V9)),8,999))
        gff$gene <- sub(";","",substr(regmatches(gff$V9,regexec("gene=.*?;",gff$V9)),6,999))
        gffs[[name]] <- gff
    }
    map <- as.data.frame(t(apply(lib,1,function(y) map_gene(y,gffs))))

    return(map)
}

calculate_2dist_zscore <- function(u1,s1,u2,s2){
    z <- (u1-u2)/sqrt((s1^2)+(s2^2))
    return(z)
}

read_controls <- function(file){
    controls <- read.table(file,stringsAsFactors=F)
    colnames(controls) <- c("Control","Tag","Type","Concentration")
    return(controls)
}

map_controls <- function(controls,maps){
    names(maps) <- unlist(lapply(maps,rownames))
    map <- do.call(rbind,maps)
    control_barcodes <- controls[,2]
    control_map <- map[control_barcodes,]
    return(control_map)
}

check_controls <- function(dataset, controls, plot=TRUE, sample = "", outDir="scratch/results/"){
    control_counts <- dataset$counts[controls$Tag,]
    control_counts[is.na(control_counts)] <- 0

    # Anna's renaming bit
    n <- function(s, dataset) {
        code <- (unlist(strsplit(s, "_|\\."))[2])
        return (paste(as.character(dataset$meta[code,c('Mouse', 'Day', 'Locus')]), collapse='-'))
    }

    titles <- sapply(colnames(control_counts), n, dataset=dataset)


    if(plot){
        pdf(paste(outDir, "/", sample, dataset$name, "_controls.pdf", sep=""), width=3*length(unique(controls$Type)), height=3)
        par(mfrow=c(1, length(unique(controls$Type))))
    }

    linearity <- matrix(0, nrow=ncol(control_counts), ncol=length(unique(controls$Type)))
    rownames(linearity) <- colnames(control_counts)
    colnames(linearity) <- unique(controls$Type)

    for(sample in 1:ncol(control_counts)){
        for(type in unique(controls$Type)){
            set <- which(controls$Type==type)
            subcontrols <- controls[set,]
            subcounts <- control_counts[set,]
            linearity[sample, type] <- cor(controls[set,]$Concentration, control_counts[set,sample])
            plot(controls[set,]$Concentration, control_counts[set,sample]+1, pch=20, main=paste(titles[sample], type), sub=paste("R =",round(linearity[sample, type],3)), log="xy", xlab="Nominal Concentration", ylim=c(1, 100000), ylab="Counts")
        }
    }

    dev.off()

    dataset$control_stats <- as.data.frame(linearity)
    return(dataset)
}




filter_dataset <- function(dataset, type, cutoff){
    pass <- (dataset$control_stats[,type] >= cutoff) & (!is.na(dataset$control_stats[,type]))
    dataset$control_stats$pass <- pass
    dataset$meta <- dataset$meta[pass,]
    dataset$counts <- dataset$counts[,pass]
    return(dataset)
}

calculate_comps <- function(dataset, controls, type){
    genes <- unique(dataset$map$Feature)
    genes <- genes[genes!="-"]
    other_barcodes <- rownames(dataset$map)[dataset$map$Feature=="-"]

    control_tags <- controls[controls$Type==type,]$Tag

    comps <- list()
    for(day in names(dataset$fitness)){
        control_fits <- dataset$fitness[[day]][control_tags,]
        control_fits <- control_fits[!is.na(control_fits$baseMean),]
        control_mu <- mean(control_fits$log2FoldChange)
        control_sigma <- sqrt(sum(control_fits$lfcSE^2))/nrow(control_fits)
        gene_comps <- list()
        used_barcodes <- c()
        for(gene in genes){
            gene_tags <- rownames(dataset$map[dataset$map$Feature==gene,])
            gene_fits <- dataset$fitness[[day]][gene_tags,,F]
            gene_fits <- gene_fits[!is.na(gene_fits$baseMean),,F]
            if(nrow(gene_fits)>0){
                gene_mu <- mean(gene_fits$log2FoldChange)
                gene_sigma <- sqrt(sum(gene_fits$lfcSE^2))/nrow(gene_fits)
                gene_comps[[gene]] <- calculate_2dist_zscore(gene_mu, gene_sigma, control_mu, control_sigma)
            }
        }
        for(barcode in other_barcodes){
            other_fit <- dataset$fitness[[day]][barcode,]
            if(!is.na(other_fit$baseMean)){
                gene_comps[[barcode]] <- calculate_2dist_zscore(other_fit$log2FoldChange, other_fit$lfcSE, control_mu, control_sigma)
            }
        }
        comps[[day]] <- unlist(gene_comps)
    }

    dataset$competitiveness <- comps
    return(dataset)
}

comp_stats <- function(comp, name){
    pvalues <- 2*pnorm(-abs(comp), 0, 1)
    padj <- p.adjust(pvalues, method="BH")
    stats <- cbind(comp, pvalues, padj)
    colnames(stats) <- paste(name, c("zscore", "pvalue", "padj"), sep="_")
    return(stats)
}

build_output <- function(dataset){
    # Stats
    stats <- lapply(names(dataset$competitiveness), function(x) comp_stats(dataset$competitiveness[[x]], x))
    stats <- do.call(cbind, stats)

    # Tag counts
    tag_tab <- table(dataset$map$Feature)
    tag_counts <- tag_tab[rownames(stats)]
    tag_counts[is.na(tag_counts)] <- 1

    # Short names
    short_name <- dataset$map[match(rownames(stats), dataset$map$Feature),]$ShortName

    output <- cbind(short_name, tag_counts, stats)
    rownames(output) <- rownames(stats)

    dataset$output <- as.data.frame(output)
    return(dataset)
}

save_controls <- function(dataset, controls){
    control_fit <- lapply(dataset$fitness, function(x) x[controls$Tag,])
    dataset$controls <- cbind(controls, do.call(cbind, control_fit))
    return(dataset)
}

plot_fitness <- function(fitness){
    par(mfrow=c(1,2))

    subset <- fitness[fitness$baseMean>100,]
    cor1 <- lm(subset$log2FoldChange ~ log10(subset$baseMean))
    plot(log10(subset$baseMean), subset$log2FoldChange, xlim=c(2,log10(max(fitness$baseMean))), ylim=c(-2,2), xlab="log10 Normalised Mean Count", ylab="log2 Fold Change", pch= "..")
    abline(cor1, col=2)
    cor2 <- lm(subset$lfcSE ~ log10(subset$baseMean))
    plot(log10(subset$baseMean), subset$lfcSE, xlim=c(2,log10(max(fitness$baseMean))), xlab="Normalised Mean Count", ylab="LFC SE", pch= "..")
    abline(cor2, col=2)
}


demux_code_to_meta <- function(count_name, dataset){ # Example of count_name dnaid2023_66.count
        code <- (unlist(strsplit(count_name, "_|\\."))[2])
        return (paste(as.character(dataset$meta[code,c('Mouse', 'Day', 'Locus')]), collapse='-'))
}

rename_counts <- function(dataset){
    titles <- sapply(colnames(dataset$counts), demux_code_to_meta, dataset=dataset)
    names(dataset$counts)[match(names(titles), names(dataset$counts))] <- titles
    merged_counts <- merge(dataset$counts, dataset$map, by=0, all.x=TRUE)
    rownames(merged_counts) <- merged_counts$Row.names
    merged_counts <- merged_counts[, names(merged_counts) != "Row.names"]
    dataset$merged_counts <- merged_counts
    return(dataset)

}


args <- commandArgs(trailingOnly=TRUE)

metadataFile <- args[1]
controlsFile <- args[2]
countsDir <- args[3]
mapDir <- args[4]
outDir <- args[5]
sample <- args[6]


#
#controlsFile <- "../data/metadata/controls.txt"
#mapDir <- "/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/chris/hardt/nguyenb/tnseq_mapping_2031/scratch"
#outDir <- "../scratch/test/"
#
## Sample 2023
#metadataFile <- "../data/metadata/dnaid2023_metadata.txt"
#countsDir <- "../scratch/counts/dnaid2023"
#sample <- "dnaid2023_"
#

# Sample 2029
#metadataFile <- "../data/metadata/dnaid2029_metadata.txt"
#countsDir <- "../scratch/counts/dnaid2029"
#sample <- "dnaid2029_"

resultsPrefix <- paste0(outDir, "/", sample, "results_")
controlsPrefix <- paste0(outDir, "/", sample, "controls_")
countsPrefix <- paste0(outDir, "/", sample, "counts_")
featurePrefix <- paste0(outDir, "/", sample, "features_")



# Read in metadata
meta <- read_metadata(metadataFile)

# Determine datasets
dataset_names <- unique(meta$Inoculum)

# Build datasets
datasets <- lapply(dataset_names, function(x) TnSeqDataSet(meta, x, qpath=countsDir, mpath=mapDir, sample=sample))

# Check controls and filter
controls <- read_controls(controlsFile)
datasets <- lapply(datasets, function(x) check_controls(x, controls, sample=sample, outDir=outDir))
datasets <- lapply(datasets, function(x) filter_dataset(x, "wt", 0.8))

# Calculate fitnesses
datasets <- lapply(datasets, calculate_fitness)

# Update controls
datasets <- lapply(datasets, function(x) save_controls(x, controls))

# Calculate competitiveness scores
datasets <- lapply(datasets, function(x) calculate_comps(x, controls, "wt"))

# Construct a nice output table
datasets <- lapply(datasets, build_output)
datasets <- lapply(datasets, rename_counts)
names(datasets) <- dataset_names

# Write output
lapply(dataset_names, function(x) write.table(datasets[[x]]$output, paste(resultsPrefix, x, ".txt", sep="")))
lapply(dataset_names, function(x) write.table(datasets[[x]]$controls, paste(controlsPrefix, x, ".txt", sep="")))
lapply(dataset_names, function(x) write.table(datasets[[x]]$merged_counts, paste(countsPrefix, x, ".txt", sep="")))
lapply(dataset_names, function(x) write.table(datasets[[x]]$map, paste(featurePrefix, x, ".txt", sep="")))