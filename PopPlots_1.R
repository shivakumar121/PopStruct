#source("http://bioconductor.org/biocLite.R")
#biocLite("SNPRelate")
library ("SNPRelate")
library ("VariantAnnotation")
library ("ggplot2")
library("ggrepel")
library("ape")
##################################################################################################
PFChromosomes <- c("M76611", "PFC10_API_IRAB", "Pf3D7_01_v3", "Pf3D7_02_v3", "Pf3D7_03_v3", "Pf3D7_04_v3", "Pf3D7_05_v3", "Pf3D7_06_v3", "Pf3D7_07_v3", "Pf3D7_08_v3", "Pf3D7_09_v3", "Pf3D7_10_v3", "Pf3D7_11_v3", "Pf3D7_12_v3", "Pf3D7_13_v3", "Pf3D7_14_v3")
InputVCFFile <- "/media/rathodlab/Data_disk/NGS_home/India_Patient_NGS_home/AllSamplesPlusDd2_Oct242016/Mpileup/PopStudyOct212016_1.AllChr.mpileup.calls.annotated.vcf.gz"
OutDir <- "/media/rathodlab/Data_disk/NGS_home/India_Patient_NGS_home/GenomicsPaper/Feb072017_WG/R_OutDir"
PopInfo <- read.csv ("/media/rathodlab/Data_disk/NGS_home/India_Patient_NGS_home/AllSamplesPlusDd2_Oct242016/PopInfo.csv")
PopInfo_GMC_RMRC <- read.csv ("/media/rathodlab/Data_disk/NGS_home/India_Patient_NGS_home/AllSamplesPlusDd2_Oct242016/PopInfo_GMC_RMRC.csv")
ColorReference <- read.csv ("/media/rathodlab/Data_disk/NGS_home/India_Patient_NGS_home/AllSamplesPlusDd2_Oct242016/ColorReference.csv")
##################################################################################################
for (ChrName in PFChromosomes)
{
  if (!(file.exists (paste0 (OutDir, "/", ChrName, ".filtered_Feb072017.vcf"))))
  {
    VCFRegion_Location <- GRanges (seqnames = ChrName, ranges = IRanges (start = 0, end = 29900000, names = ChrName))
    MyVCF <- readVcf (file = InputVCFFile, param = VCFRegion_Location, genome = "Pf3D7_v3")
    MyVCF <- MyVCF[,-c(9,10)]  ## Remove samples that are bad
    MyVCF <- MyVCF[,-c(22,25,29,30,32,34,74,75,82,93,99,107,110,112)]
    ReadDepth <- geno(MyVCF)$DP
    RowIndex <- which (ReadDepth <= 10, arr.ind = T)[,1]
    RowIndex <- sort (unique (RowIndex))
    MyVCF <- MyVCF[-RowIndex,]
    HighQualIndex <- which (rowRanges(MyVCF)$QUAL >= 100)
    MyVCF <- MyVCF[HighQualIndex,]
    ### Find Index of Non Intergenic #####
    #Index_HIGH_MODERATE <- grep (as.list (info (MyVCF)$EFF), pattern = "HIGH", ignore.case = F)
    #Index_HIGH_MODERATE <- c(Index_HIGH_MODERATE, grep (as.list (info (MyVCF)$EFF), pattern = "MODERATE", ignore.case = F))
    #Index_HIGH_MODERATE <- sort(unique(Index_HIGH_MODERATE))
    #MyVCF <- MyVCF[Index_HIGH_MODERATE,]
    ### Find Index of VARs ###
    IndexRifin = grep (as.list (info (MyVCF)$EFF), pattern = "rifin")
    #IndexIntergenic = grep (as.list (info (MyVCF)$EFF), pattern = "intergenic_region")
    IndexStevor = grep (as.list (info (MyVCF)$EFF), pattern = "stevor")
    IndexEMP = grep (as.list (info (MyVCF)$EFF), pattern = "erythrocyte\\+membrane\\+protein")
    IndexVARs = c(IndexRifin, IndexStevor, IndexEMP)
    IndexVARs = sort (unique (IndexVARs))
    if (length (IndexVARs) > 0)
    {
      MyVCF = MyVCF[-(IndexVARs),]
    }
    writeVcf (obj = MyVCF, filename = paste0 (OutDir, "/", ChrName, ".filtered_Feb072017.vcf"))
    remove (MyVCF)
    gc()
  }
}
####### Concat All Filtered VCFs #######################################
#AllVCFFiles <- list.files (path = OutDir, pattern = ".filtered_NoVARs.Oct282016.vcf$", full.names = T)
if (!(file.exists(paste0 (OutDir, ".filtered_Feb072017.vcf"))))
{
  AllInputVCFFileNames = NULL
  for (i in 1:length(PFChromosomes))
  {
    AllInputVCFFileNames = c(AllInputVCFFileNames, paste0 (OutDir, "/", PFChromosomes[i], ".filtered_Feb072017.vcf"))
  }
  AllInputVCFFileNames = paste (AllInputVCFFileNames, collapse = " ")
  system (paste0 ("bcftools concat ", AllInputVCFFileNames, " -o " , OutDir, "/AllChr.filtered_Feb072017.vcf"))
}
###### Convert VCF to GDS ##################
if (!(file.exists (paste0 (OutDir, "/", "AllChr.filtered_Feb072017.gds"))))
{
  snpgdsVCF2GDS(vcf.fn = paste0 (OutDir, "/AllChr.filtered_Feb072017.vcf"), out.fn = paste0 (OutDir, "/", "AllChr.filtered_Feb072017.gds"), method="copy.num.of.ref")
}
#Read GDS file
genofile <- snpgdsOpen(paste0 (OutDir, "/", "AllChr.filtered_Feb072017.gds"))
pca <- snpgdsPCA(genofile, autosome.only = F)
snpset <- snpgdsLDpruning(genofile, ld.threshold = 0.3, autosome.only = F, method = "corr", remove.monosnp = F)
snpset.id <- unlist(snpset)
Anti.snp.id <- read.gdsn(index.gdsn(genofile, "snp.id"))
Anti.snp.id <- setdiff (Anti.snp.id, snpset.id)
pca <- snpgdsPCA(genofile, autosome.only = F, snp.id = snpset.id)
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  EV5 = pca$eigenvect[,5],
                  EV6 = pca$eigenvect[,6],
                  EV7 = pca$eigenvect[,7],
                  EV8 = pca$eigenvect[,8],
                  stringsAsFactors = FALSE,
                  pop = factor(PopInfo$Region)[match(pca$sample.id, PopInfo$SampleName)])
head(tab)
ShortPop <- as.character (tab$pop)
Index <- grep (ShortPop, pattern = "Cambodia", ignore.case = T)
ShortPop[Index] <- "Cambodia"
Index <- grep (ShortPop, pattern = "Thailand", ignore.case = T)
ShortPop[Index] <- "Thailand"
tab$pop <- as.factor (ShortPop)
## Add Number to Locations
tab$pop <- as.character (tab$pop)
for (i in 1:length(tab$pop))
{
  if (tab$pop[i] == "India")
  {
    tab$pop[i] <- "1 India"
  } else if (grepl ("Bangladesh", pattern = tab$pop[i]))
  {
    tab$pop[i] <- "2 Bangladesh"
  } else if (grepl ("CentralMyanmar", pattern = tab$pop[i]))
  {
    tab$pop[i] <- "3 Myanmar"
  } else if (grepl ("Thailand", pattern = tab$pop[i]))
  {
    tab$pop[i] <- "4 Thailand"
  } else if (grepl ("Cambodia", pattern = tab$pop[i]))
  {
    tab$pop[i] <- "5 Cambodia"
  } else if (grepl ("Laos,Vietnam", pattern = tab$pop[i]))
  {
    tab$pop[i] <- paste0 ("6 ", tab$pop[i])
  } else if (grepl ("Mozambique", pattern = tab$pop[i]))
  {
    tab$pop[i] <- "7 Mozambique"
  } else if (grepl ("D.R.Congo,Gambia,Ghana,Mali,Nigeria", pattern = tab$pop[i]))
  {
    tab$pop[i] <- paste0 ("8 ", tab$pop[i])
  } else if (grepl ("Brazil,French Guiana", pattern = tab$pop[i]))
  {
    tab$pop[i] <- paste0 ("9 ", tab$pop[i])
  } else if (grepl ("Dd2", pattern = tab$pop[i]))
  {
    tab$pop[i] <- paste0 ("Z ", tab$pop[i])
  }
  
}  
## Add Color Reference
MyColors <- NULL
MyColors <- as.character (MyColors)
for (i in 1:nrow (ColorReference))
{
  TempIndex <- grep (tab$pop, pattern = ColorReference$Region[i])
  MyColors[TempIndex] <- as.character (ColorReference$Color[i])
}
tab$MyColors <- MyColors
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")
ColorReference$Region <- as.character (ColorReference$Region)
ColorReference$Color <- as.character (ColorReference$Color)
###### Make GGplot2 Plots ####################
p1 <- ggplot (data = tab, aes(EV1,EV2))
p1 <- p1 + geom_point(size = 5, aes (color = factor (pop)))
p1 <- p1 + scale_color_manual (values = ColorReference[,2], name = "Region")
p1 <- p1 + xlim (-0.2,0.16) + ylim (-0.1, 0.2)
p1 <- p1 + theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_line(colour = "#f2f2f2"), 
                 plot.title = element_text(size = rel(1.5)), axis.title.x = element_text(size = rel(1.2)), axis.title.y = element_text(size = rel(1.2)),
                 legend.title = element_text(size = rel(1.2)), legend.text=element_text(size = rel(1.2)))
#p1 <- p1 + geom_text_repel (data = subset (tab, pop == "1 India") , aes(label = sample.id), force = 1, size = 6)
p1
###### make plots without India samples ##################
MyCols <- ColorReference[,2]
MyCols[1] <- "white"
p1 <- ggplot (data = tab, aes(EV1,EV2))
p1 <- p1 + geom_point(size = 5, aes (color = factor (pop)))
p1 <- p1 + scale_color_manual (values = MyCols, name = "Region")
#p1 <- p1 + xlim (-0.2,0.16) + ylim (-0.1, 0.2)
p1 <- p1 + theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_line(colour = "#f2f2f2"), 
                 plot.title = element_text(size = rel(1.5)), axis.title.x = element_text(size = rel(1.2)), axis.title.y = element_text(size = rel(1.2)),
                 legend.title = element_text(size = rel(1.2)), legend.text=element_text(size = rel(1.2)))
#p1 <- p1 + geom_text_repel (data = subset (tab, pop == "1 India") , aes(label = sample.id), force = 1, size = 6)
p1


if (!(file.exists(paste0 (OutDir, "/PCA_SNPRelateAllSamples_NoVars_Feb082017.pdf"))))
{
  ggsave (filename = paste0 (OutDir, "/PCA_SNPRelateAllSamples_NoVars_Feb082017.pdf"), plot = p1, units = "in", width = 11, height = 8.5)
}
##### Color NEIndia and SWIndia separately ###################
EastIndiaSamples <- c("akxy","diby","brdz","ixyr","sxbz","monj","gguj","aaam","hfff","ottv","faea","hrxz")
WestIndiaSamples <- c("2NRN","5TBX","durf","bngc","cgiy","czum","gabv","hpve","hxjy","istd","JZT6","kmml","kphn","lwct","NUSP","oego","oirv","rwht","TCFT","ulpg","XGK2","yyfl","znmf")
Pop_NE_SW <- rep("3 Non-India", times = nrow(tab))
Pop_NE_SW[match(EastIndiaSamples, tab$sample.id)] <- "1 N.E. India"
Pop_NE_SW[match(WestIndiaSamples, tab$sample.id)] <- "2 S.W. India"
tab$Pop_NE_SW <- Pop_NE_SW
ColorReference_NE_SW <- c("#e60000", "#e60000", "#cccccc")
p1 <- ggplot (data = tab, aes(EV1,EV2))
p1 <- p1 + geom_point(size = 5, aes (color = factor (Pop_NE_SW))) 
p1 <- p1 + scale_color_manual (values = ColorReference_NE_SW, name = "Region")
p1 <- p1 + geom_point (data = subset (tab, Pop_NE_SW == "2 S.W. India"), size = 2, color = "yellow")
#p1 <- p1 + geom_point (data = subset (tab, Pop_NE_SW == "2 S.W. India"), size = 3, color = "black", show.legend = TRUE)
p1 <- p1 + theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_line(colour = "#f2f2f2"), 
                 plot.title = element_text(size = rel(1.5)), axis.title.x = element_text(size = rel(1.2)), axis.title.y = element_text(size = rel(1.2)),
                 legend.title = element_text(size = rel(1.2)), legend.text=element_text(size = rel(1.2)))
p1

