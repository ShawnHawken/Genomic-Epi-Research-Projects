# MLST ANALYSIS

#MAKE A TABLE OF MLST TYPES OF ISOLATES IN OUR STUDY

#mlst types
mlst_types<-read.table("../data/2016-10-27_genomes_eclo_mlst.txt", sep = "\t", header = TRUE, quote = "")
colnames(mlst_types)[1:2]<-c("genome","ST")

#read in phenotype df to get genome names in UM study

carb_phenotype<-read.table("2016-11-29-carb_phenotype", sep=",", header = T)

UM_mlst<-rbind(mlst_types[grepl("UM-CRE.*",mlst_types$genome), ], mlst_types[grepl("Eclo_436.*" , mlst_types$genome), ],mlst_types[grepl("Eclo_437.*" , mlst_types$genome), ])
UM_mlst$genome<-sub("_genome","", UM_mlst$genome)
UM_mlst$genome<-sub("Eclo_","", UM_mlst$genome)
UM_mlst<-UM_mlst[UM_mlst$genome %in% as.vector(carb_phenotype$genome),] #subset based on genomes in UM analysis

sort(table(UM_mlst$ST))
