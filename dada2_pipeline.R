library(dada2)
library(ggplot2)
library(dplyr)

path = "./Demultiplexed/Trimmed"
Filter_out = "./Filttrim"


#---load files
# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fq.gz and SAMPLENAME_R2.fq.gz
fwd_files <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
rev_files <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fwd_files), "_"), `[`, 1)
sample.names = substr(sample.names, start = 9, stop = nchar(sample.names))

#---Filter and Trim
filtF1 <- file.path("./Filttrim", paste0(sample.names, "_F_filt.fastq.gz"))
filtR1 <- file.path("./Filttrim", paste0(sample.names, "_R_filt.fastq.gz"))

filtered_seqs <- filterAndTrim(fwd = fwd_files, #Enter the details of the filter
                               filt = filtF1, 
                               rev = rev_files, 
                               filt.rev = filtR1,
                               #truncLen = c(0,0), #TruncLen
                               trimLeft = 10, #Trim_left
                               maxN = 0, #MaxN
                               maxEE = 2, #MaxEE
                               compress = TRUE, 
                               verbose = TRUE)
                               
saveRDS(filtered_seqs, file = "./Pipeline_output/filtered_seqs.RData")

#---Dereplicate
derepF1 <- derepFastq(filtF1, verbose=TRUE)
derepR1 <- derepFastq(filtR1, verbose=TRUE)

names(derepF1) = sample.names #doesn't work for single files, also not necessary
names(derepR1) = sample.names #doesn't work for single files, also not necessary

saveRDS(derepF1, file = "./Pipeline_output/derepF1.RData") #To save for later use.
saveRDS(derepF1, file = "./Pipeline_output/derepR1.RData") #To save for later use.

#---ErrorRate
errF <- learnErrors(derepF1, multithread=TRUE)
errR <- learnErrors(derepR1, multithread=TRUE)

saveRDS(errF, file = "./Pipeline_output/errF.RData") #To save for later use.
saveRDS(errR, file = "./Pipeline_output/errR.RData") #To save for later use.

#---sample inference
dadaFs = dada(derepF1, err=errF, multithread = TRUE)
dadaRs = dada(derepR1, err=errR, multithread = TRUE)

saveRDS(dadaFs, file = "./Pipeline_output/dadaFs.RData")
saveRDS(dadaRs, file = "./Pipeline_output/dadaRs.RData")

#---merging paired reads
mergers <- mergePairs(dadaFs, derepF1, dadaRs, derepR1, verbose=TRUE, trimOverhang = TRUE)

ASVs = mergers
saveRDS(ASVs, file = "./Pipeline_output/ASVs.RData")

#---Constructing a sequence table
seqtab = makeSequenceTable(ASVs)
dim(seqtab)
table(nchar(getSequences(seqtab)))

#---Remove the chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

saveRDS(seqtab.nochim, file = "./Pipeline_output/seqtabnochim.RData")

#---Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(filtered_seqs, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

write.csv(track, file = "./Pipeline_output/Track.csv")

#---Assigning Taxonomy using dada2 and PR²
taxa <- assignTaxonomy(seqtab.nochim, "pr2_version_5.0.0_SSU_dada2.fasta.gz", taxLevels = c("Domain","Supergroup","Division","Subdivision", "Class","Order","Family","Genus","Species"), multithread=TRUE, verbose=TRUE) 

saveRDS(taxa, file =  "./Pipeline_output/taxaPR2.RData")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#---Making a read abundance table
ASV_tab = cbind(ASV_tab, taxa) #Add The taxonomic annotation
rownames(ASV_tab) = NULL
head(ASV_tab)

write.csv(ASV_tab, file = "./Pipeline_output/read_abundance_table_PR2.csv")
