#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : WES PID                                                                    |
#  Data Owner  : Newcastle University - Prof. Sophie Hambleton                              |
#  Description : Parse ExomeDepth results and produce a Karyogram with annotated CNV sites  |
#-------------------------------------------------------------------------------------------#



##'load essential packages and set directories
##'-----------------------------------------------------------------------------------------#
library(readr)
library(dplyr)
library(ggbio)
library(openxlsx)
library(GenomicRanges)
library(tools)

pack_dir <- "/Volumes/IGM_BRC-Genomics/Sophie_Exomes/Andrew_ExomeDepth/"
out_dir  <- "/Volumes/IGM_BRC-Genomics/Sophie_Exomes/Andrew_ExomeDepth/Karyograms/"
files_in <- list.files(path = pack_dir, pattern = "*.xlsx", recursive = T, full.names = T)
##'-----------------------------------------------------------------------------------------#




##'Loop through results, parse ExomeDepth output, produce labelled karyogram
##'-----------------------------------------------------------------------------------------#
for(i in files_in) {
  f_name   <- i %>% basename %>% file_path_sans_ext %>% gsub("ExomeDepth_","",.)
  ss       <- read.xlsx(i) %>%
              as.data.frame

  # Subset by a particular chomosome
  #%>%
  # filter(chromosome == 19)
  ss.range <- GRanges(seqnames = paste0("chr",ss$chromosome),
                      ranges   = IRanges(start = ss$start,
                                         end   = ss$end),
                      strand   = "*")

  data(hg19IdeogramCyto, package = "biovizBase")
  hg19     <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", c(1:22, "X", "Y")))

  gg       <- ggplot(hg19) +
              layout_karyogram(cytoband = TRUE) +
              layout_karyogram(ss.range, geom = "rect", ylim = c(11, 21), color = "red") +
              ggtitle(f_name) + xlab("")
  # print(gg)

  png(paste0(out_dir,"/",f_name,".png"),
      width  = 6.98,height = 6.98,
      units  = "in",res    = 600)
  print(gg); dev.off()
}
##'-----------------------------------------------------------------------------------------#
