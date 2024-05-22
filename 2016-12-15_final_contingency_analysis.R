# These are the final funciton calls for the Enterobacter project carbR vs. CarbS analysis and blood vs other source analysis

#setwd("/Users/shawnwhitefield/Desktop/Snitkin_lab/UM_CRE-1/final_manuscript_analyses/contingency_analysis/data/")
#setwd("../data/")
library(ape)
library(phytools)
library(Biostrings)
data(BLOSUM62)
data(BLOSUM80)
data(BLOSUM100)
library(wesanderson)
library(ggplot2)
source("https://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap")
library(ComplexHeatmap)
library(colorRamps)
library(circlize)
library(heatmap3)
library(rateratio.test)

source("../scripts/contingency_analysis_functions.R")

################# READ IN THE DATA #######################
##########################################################
#READ IN THE SNP MATRIX, GENE LENGTH, TREE, DELETERIOUS MUTATION SPACE, PHENOTYPE DATA

  #READ IN FINAL INDEL MAT #81 genomes 11/29/17 FROM EVAN'S PIPELINE
    indel_mat<-read.table("../data/genomes_eclo_2016-11-25_10-04-01_mutRelRef_minDist2_RecFilt_SNPs_and_indels_with_12p5_perc_ID_of_ECNIH2.SNPmat",sep = "\t", header = FALSE, quote = "") #, row.names = 1);)

  #READ IN THE GENE LOCATION AND LENGTH INFORMATION
    gene_loc<-read.table('../data/ECNIH2_genes_location.txt',sep = "\t", header = TRUE, row.names = 1)

  #READ IN THE TREE
    tree = read.tree('../data/2016-06-22_22-23-34_SNPs_among_UM_Eclo_within_12p5_perc_ID_of_ECNIH2_minDist2_noRecFilt.tree');

  # READ IN DELETERIOUS MUTATION SPACE PER GENE TO CALCULATED EXPECTED MUTATIONS
      #this is the output of /nfs/esnitkin/Shawn/UM_CRE-1/2016-12-06_calculate_deleterious_mut_space/calculate_del_mut_space.py 
      #using ECNIH2 as the input
    del_mut_space<-read.csv("../data/2016-12-12_ECNIH2.out", header = F)
    
  # READ IN PHENOTYPE DATA THAT CORRESPONDS TO THE INDEL MATRIX
    # using carb_phenotype df from 2016-10-25_antibiogram_verification script
    #this calls genomes R if they are R by any method
    carb_phenotype<-read.table("../data/2017-01-03-carb_phenotype", sep=",", header = T)
    fulltree_carb_phenotype<-read.table("../data/2017-01-03-carb_phenotype", sep=",", header = T)
    
  # READ IN THE PATHWAY ANALYSIS COMPARISON DATABASES
    
    #from ftp://ftp.patricbrc.org/patric2/patric3/genomes/1333850.3/
    pathway_annot = read.table('../data/1333850.3.PATRIC.pathway.tab', sep = "\t", header = TRUE, quote = "");
    special_annot = read.table('../data/1333850.3.PATRIC.spgene.tab', sep = "\t", header = TRUE, quote = "");
    
    
    
##### CLEAN UP THE ROW AND COLUMN NAMES AND LABELS #######
########################################################## 
    
  #FIX THE COLUMN NAMES OF THE DELETERIOUS MUTATION SPACE DATAFRAME
    
    del_mut_space<-fix.colnames_delmutspace(deleterious_mutation_space_df = del_mut_space)
    
  #FIX THE TREE TIP LABELS SO THEY ARE EASIER TO READ
    
    tree<-fix.tree_labels(a_tree = tree)
    
  #REMOVE INTERGENIC ROWS FROM VARIANT MATRIX AND ADD ROWNAMES THAT ARE NOT DUPLICATED
    
    indel_mat_includes_intgenics<-indel_mat #full dataset that includes intergenic rows
    
    indel_mat<-remove.intergenic_rows(indel_mat)
    
  #FIX THE COLUMN NAMES OF THE VARIANT MATRIX
    
    indel_mat<-fix.colnames_var_mat(indel_mat)
  
  # MAKE THE PHENOTYPE COLUMNS INTO VECTORS
    carb_phenotype$genome<-as.vector(carb_phenotype$genome)
    carb_phenotype$carbR<-as.vector(carb_phenotype$carbR)
    # ADD PHENOTYPES FOR BODY SITE OF ISOLATION
      carb_phenotype$blood<-TRUE
      carb_phenotype$blood[grepl('^UM', carb_phenotype$genome)]<-FALSE# carb_phenotype$blood
      carb_phenotype$blood[carb_phenotype$genome == "UM-CRE-24"]<-TRUE #24,40 and 1 are the only blood isolates that arent from vic UM_CRE_45  is JP drain and17 is "fluid" may be similar
      carb_phenotype$blood[carb_phenotype$genome == "UM-CRE-40"]<-TRUE
      carb_phenotype$blood[carb_phenotype$genome == "UM-CRE-1"]<-TRUE
        #carb_phenotype$blood[carb_phenotype$genome == "UM-CRE-45"]<-TRUE # JP drain
       #carb_phenotype$blood[carb_phenotype$genome == "UM-CRE-17"]<-TRUE # "fluid" 
  
  # CLEAN UP ROWNAMES FOR GENE LOCATTION DF AND CALCULATE GENE LENGTH
      # DETERMINE LENGTH OF GENES
      rownames(gene_loc)<-sub("Eclo_ECNIH2_genome:","", rownames(gene_loc)) #clean up row names so they match the genes
      gene_loc$length<- abs(gene_loc$start - gene_loc$end)
      


##### GIVE VARIANT MATRIX A NEW REFERENCE GENOME#######
######################################################## 
    
  # GIVE THE VARIANT MATRIX A NEW REFERENCE GENOME
    
   indel_mat<- give.new_reference(indel_mat,1,"43677")
    
    
##### MAKE MATRICES THAT ARE +/- KPC CONTAINING GENOMES & GENOMES THAT ARE CLOSELY RELATED#######
#################################################################################################
   
  # VECTOR OF KPC POSITIVE GENOMES TO USE FOR SUBSETTING 
    #remove KPC #"UM-CRE-14" "UM-CRE-15" "UM-CRE-20" "UM-CRE-31" "UM-CRE-32" "UM-CRE-35" "UM-CRE-45" "UM-CRE-7"  "43693" 
    KPC<-c("UM-CRE-14","UM-CRE-15" ,"UM-CRE-20" ,"UM-CRE-31","UM-CRE-32" ,"UM-CRE-35", "UM-CRE-45", "UM-CRE-7","43693")
  
  # KPC - SUBSET VARIANT MATRIX
    KPC_minus_var_mat<-indel_mat[,c(!(colnames(indel_mat) %in% KPC))] #71 genomes included
    
  # KPC + CONTAINING VARIANT MATRIX (ALL GENOMES)
    
    
########################IDENTIFY CORE GENOME GENES############################################
############################################################################################## 
    
    #GET LIST OF ALL GENES IN MATRIX TO USE FOR SUBSETTING BASED ON CORE GENOME
      #GENES total matrix
      var_desc = as.matrix(rownames(indel_mat))
      genes = var_desc
      genes = gsub("^.+Eclo_ECNIH2_genome:","",genes)
      genes = gsub(" \\(.+\\)$" ,"", genes, perl = TRUE);
      #unique genes for subsetting on core genome
      uniq_genes<-unique(genes) #unique genes
    
    
    #IDENTIFY NUMBER OF POSITIONS IN EACH GENE THAT ARE NOT -1S (NOT MISSING) IN EACH GENOME
      #this takes forever!
      genes_no_neg_ones<-get.non_missing_positions(indel_mat) # core genome genes in all genomes in variant matrix
    
      KPC_minus_genes_no_neg_ones<-genes_no_neg_ones[ ,!(colnames(genes_no_neg_ones) %in% KPC)] #remove KPC genome columns 
      
    # GET LIST OF CORE GENOME GENES
      
      core_genome_genes<-get.core_genome_genes(genes_no_neg_ones, core_percentage = .95) # vector of core genome genes in dataset with KPC+ genomes included
        #3069 core genome genes
      non_KPC_core_genome_genes<-get.core_genome_genes(KPC_minus_genes_no_neg_ones, core_percentage = .95) # vector of core genome genes in dataset without KPC+ genomes
        #3091 core genome genes without KPC genomes included
 
      
########################SUBSET VARIANT MATRIX BASED ON CORE GENOME############################
##############################################################################################      
      
    # SUBSET INDEL MAT BASED ON CORE GENOME
      
      #FULL DATASET
      indel_mat<-indel_mat[genes %in% core_genome_genes, ]
      
      #KPC - GENOMES
      KPC_minus_var_mat<-KPC_minus_var_mat[genes %in% non_KPC_core_genome_genes, ]
      
     
    
########################SUBSET PHENOTYPES TO MATCH MATRICES WITH AND WITHOUT KPC##############
##############################################################################################  
      
    #PHENOTYPE FOR KPC- GENOMES
      KPC_minus_carb_phenotype<-carb_phenotype[as.vector(carb_phenotype$genome) %in% colnames(KPC_minus_var_mat),] #phenotypes match dataset that excludes KPC genomes
      #SORT PHENOTYPE BY MATRIX
      KPC_minus_carb_phenotype<-KPC_minus_carb_phenotype[order(KPC_minus_carb_phenotype$genome[colnames(KPC_minus_var_mat)]),]
      
    #PHENOTYPE FOR ALL GENOMES 
      carb_phenotype<-carb_phenotype[as.vector(carb_phenotype$genome) %in% colnames(indel_mat),] #phenotypes match full dataset
      #SORT PHENOTYPE BY MATRIX
      carb_phenotype<-carb_phenotype[order(carb_phenotype$genome[colnames(indel_mat)]),]
      

########################IDENTIFY VARIANT TYPES INCLUDED IN VARIANT MATRIX##################### 
##############################################################################################
    
      #IDENTIFY TYPES OF VARIANTS
      variant_types_list<-get.variant_types(indel_mat) # full dataset
      #var_desc, genes, annots, annot_by_gene, stops, start_aa, end_aa, ns_mut, syn_mut, del_mut,pos
        all_genes<-as.vector(unlist(variant_types_list["genes"]))
        all_var_desc<-as.vector(unlist(variant_types_list["var_desc"]))
        all_annots<-as.vector(unlist(variant_types_list["annots"]))
        all_annot_by_gene<-as.vector(unlist(variant_types_list["annot_by_gene"]))
        all_stops<-as.vector(unlist(variant_types_list["stops"]))
        all_start_aa<-as.vector(unlist(variant_types_list["start_aa"]))
        all_end_aa<-as.vector(unlist(variant_types_list["end_aa"]))
        all_ns_mut<-as.vector(unlist(variant_types_list["ns_mut"]))
        all_syn_mut<-as.vector(unlist(variant_types_list["syn_mut"]))
        all_del_mut<-as.vector(unlist(variant_types_list["del_mut"]))
        all_pos<-as.vector(unlist(variant_types_list["pos"]))
        all_indel_row<-as.vector(unlist(variant_types_list["indel_row"]))
      
      # unique genes vector full dataset
        all_unique_genes<-unique(all_genes)
      
      
      KPC_minus_variant_types_list<-get.variant_types(KPC_minus_var_mat) #only KPC - genomes
      #var_desc, genes, annots, annot_by_gene, stops, start_aa, end_aa, ns_mut, syn_mut, del_mut,pos
        KPC_minus_genes<-as.vector(unlist(KPC_minus_variant_types_list["genes"]))
        KPC_minus_var_desc<-as.vector(unlist(KPC_minus_variant_types_list["var_desc"]))
        KPC_minus_annots<-as.vector(unlist(KPC_minus_variant_types_list["annots"]))
        KPC_minus_annot_by_gene<-as.vector(unlist(KPC_minus_variant_types_list["annot_by_gene"]))
        KPC_minus_stops<-as.vector(unlist(KPC_minus_variant_types_list["stops"]))
        KPC_minus_start_aa<-as.vector(unlist(KPC_minus_variant_types_list["start_aa"]))
        KPC_minus_end_aa<-as.vector(unlist(KPC_minus_variant_types_list["end_aa"]))
        KPC_minus_ns_mut<-as.vector(unlist(KPC_minus_variant_types_list["ns_mut"]))
        KPC_minus_syn_mut<-as.vector(unlist(KPC_minus_variant_types_list["syn_mut"]))
        KPC_minus_del_mut<-as.vector(unlist(KPC_minus_variant_types_list["del_mut"]))
        KPC_minus_pos<-as.vector(unlist(KPC_minus_variant_types_list["pos"]))
        KPC_minus_indel_row<-as.vector(unlist(KPC_minus_variant_types_list["indel_row"]))
      
      #unique genes vector KPC- genomes
        KPC_minus_unique_genes<-unique(KPC_minus_genes)
  
  ########################DEFINE PAIRS OF MAXIMUM CLOSELY RELATED GENOMES##################### 
  ############################################################################################
    
    #DISTANCE MATRIX FOR FULL DATASET
        #distance in this matrix represents the number of shared variants between genomes 
          max_shared_var_dist_mat<-get.max_shared_variants(indel_mat)
        # pairs df containing rows with pair-mates
           neighbor_pairs_df<-make.pairs_df(max_shared_var_dist_mat)
        #find genetic distances between pairs
           neighbor_pairs_dist<-find.pair_genetic_dist(neighbor_pairs_df, var_mat = indel_mat)
          
           
    #DISTANCE MATRIX FOR ONLY KPC- GENOMES
        #distance in this matrix represents the number of shared variants between genomes
          KPC_minus_max_shared_var_dist_mat<-get.max_shared_variants(KPC_minus_var_mat)
        # pairs df containing rows with pair-mates
          KPC_minus_neighbor_pairs_df<-make.pairs_df(KPC_minus_max_shared_var_dist_mat)
        #find genetic distances between pairs
          KPC_minus_neighbor_pairs_dist<-find.pair_genetic_dist(KPC_minus_neighbor_pairs_df, KPC_minus_var_mat)


########################SUBSET PAIRS BASED ON GENETIC DISTANCES ############################ 
############################################################################################
     
    #SUBSET PAIRS BASED ON THOSE THAT ARE CLOSELY RELATED ENOUGH TO BE INFORMATIVE
      # full dataset
        neighbor_pairs_df<-neighbor_pairs_df[rownames(neighbor_pairs_df)[neighbor_pairs_dist < 1000],]
        neighbor_pairs_dist<-neighbor_pairs_dist[neighbor_pairs_dist < 1000]     
        #37 genomes left
        #subset phenotype to match this matrix
        carb_phenotype<-carb_phenotype[carb_phenotype$genome %in% rownames(neighbor_pairs_df),]
        
        
      #KPC minus genomes 
        KPC_minus_neighbor_pairs_df<-KPC_minus_neighbor_pairs_df[rownames(KPC_minus_neighbor_pairs_df)[KPC_minus_neighbor_pairs_dist < 1000],]
        KPC_minus_neighbor_pairs_dist<-KPC_minus_neighbor_pairs_dist[KPC_minus_neighbor_pairs_dist < 1000]
      
        #31 genomes left
        #subset phenotype to match this matrix
        KPC_minus_carb_phenotype<-KPC_minus_carb_phenotype[KPC_minus_carb_phenotype$genome %in% rownames(KPC_minus_neighbor_pairs_df),]
        
########################IDENTIFY UNIQUE PAIR VARIANTS ##################################### 
############################################################################################
      
    #INCLUDING KPC+ GENOMES
        
        # unique stops and del_muts & ns_muts in dataset including KPC
          unique_pair_muts_with_ns_muts<-find.unique_pair_mutants_ns(neighbor_pairs_df,indel_mat,del_mut = all_del_mut, stops = all_stops, ns_mut = all_ns_mut)
        
        unique_pair_muts<-find.unique_pair_mutants(neighbor_pairs_df,indel_mat,del_mut = all_del_mut, stops = all_stops)
        
        # unique indels in dataset including KPC
          unique_pair_indels<-find.unique_pair_indels(neighbor_pairs_df, indel_mat,indel_row = all_indel_row)
          
       # unique indels stops and del_muts (for fisher tests)
          unique_pair_indels_muts<-find.unique_pair_indels_del_muts(neighbor_pairs_df, indel_mat,indel_row = all_indel_row,stops = all_stops,del_mut = all_del_mut  )
          
          #FOR PLOTTING PATHWAY ANALYSIS WITH INDELS AND MUTATIONS
          unique_pair_indels_muts_with_ns_muts<-find.unique_pair_indels_del_muts_ns(neighbor_pairs_df, indel_mat,indel_row = all_indel_row,stops = all_stops,del_mut = all_del_mut, ns_mut = all_ns_mut  )
          
          
    #EXCLUDING KPC+ GENOMES
        
        # unique stops and del_muts in dataset excluding KPC
          KPC_minus_unique_pair_muts<-find.unique_pair_mutants(KPC_minus_neighbor_pairs_df, KPC_minus_var_mat, 
                                              del_mut = KPC_minus_del_mut,stops = KPC_minus_stops)
        
        # unique indels in dataset excluding KPC
          KPC_minus_unique_pair_indels<-find.unique_pair_indels(KPC_minus_neighbor_pairs_df, KPC_minus_var_mat, indel_row = KPC_minus_indel_row)
          
        # unique indels stops and del_muts (for fisher tests)
          KPC_minus_unique_pair_indels_muts<-find.unique_pair_indels_del_muts(KPC_minus_neighbor_pairs_df, KPC_minus_var_mat, indel_row = KPC_minus_indel_row, stops = KPC_minus_stops, del_mut = KPC_minus_del_mut  )
          
          
##########################PLOT GENETIC DISTANCES TO EXAMINE DIFFERENCES BETWEEN PHENOTYPE STRATUM###########################
############################################################################################################################          
          
  #PLOT GENETIC DISTANCE BOX PLOT BY PHENOTYPE
        #full dataset
          #total distance BLOOD
          blood_and_other_dist<-as.data.frame(neighbor_pairs_dist)
          blood_and_other_dist$blood<-carb_phenotype$blood 
          colnames(blood_and_other_dist)<-c("distance","blood")
          p <- ggplot(blood_and_other_dist, aes(factor(blood), distance))
          p1 <- p + geom_boxplot(outlier.shape = 3)
          p1 + geom_point(position = position_jitter(width = 0.2))
          p2 <- p + geom_boxplot(outlier.colour = NA)
          p3<-p2 + geom_point(position = position_jitter(width = 0.2)) 
          p3 + ggtitle("Blood vs.Other genetic distance between pairs full dataset")
          
          #total distance CARB RESISTANCE
          blood_and_other_dist<-as.data.frame(neighbor_pairs_dist)
          blood_and_other_dist$blood<-carb_phenotype$carbR 
          colnames(blood_and_other_dist)<-c("distance","carbapenem_resistant")
          p <- ggplot(blood_and_other_dist, aes(factor(carbapenem_resistant), distance))
          p1 <- p + geom_boxplot(outlier.shape = 3)
          p1 + geom_point(position = position_jitter(width = 0.2))
          p2 <- p + geom_boxplot(outlier.colour = NA)
          p3<-p2 + geom_point(position = position_jitter(width = 0.2)) 
          p3 + ggtitle("carbR vs carbS genetic distance between pairs full dataset")
        
        #KPC- dataset
          #total distance blood vs other
          blood_and_other_dist<-as.data.frame(KPC_minus_neighbor_pairs_dist)
          blood_and_other_dist$blood<-KPC_minus_carb_phenotype$blood 
          colnames(blood_and_other_dist)<-c("distance","blood")
          p <- ggplot(blood_and_other_dist, aes(factor(blood), distance))
          p1 <- p + geom_boxplot(outlier.shape = 3)
          p1 + geom_point(position = position_jitter(width = 0.2))
          p2 <- p + geom_boxplot(outlier.colour = NA)
          p3<-p2 + geom_point(position = position_jitter(width = 0.2)) 
          p3 + ggtitle("Blood vs.Other genetic distance between pairs KPC minus dataset")
          
          # total distance carbR vs carbS
          blood_and_other_dist<-as.data.frame(KPC_minus_neighbor_pairs_dist)
          blood_and_other_dist$blood<-KPC_minus_carb_phenotype$blood 
          colnames(blood_and_other_dist)<-c("distance","carbapenem_resistant")
          p <- ggplot(blood_and_other_dist, aes(factor(carbapenem_resistant), distance))
          p1 <- p + geom_boxplot(outlier.shape = 3)
          p1 + geom_point(position = position_jitter(width = 0.2))
          p2 <- p + geom_boxplot(outlier.colour = NA)
          p3<-p2 + geom_point(position = position_jitter(width = 0.2)) 
          p3 + ggtitle("carbR vs carbS genetic distance between pairs KPC minus dataset")
          
##########################SUBSET UNIQUE VARIANT MATRIX TO SPEED-UP CONTINGENCY TESTS###########################
###############################################################################################################
          
# SUBSET GENE LIST TO CHECK FOR VARIANTS, THIS SPEEDS THINGS UP QUITE A BIT.
    
    
    # full dataset del_muts and stops and del_muts
    all_genes_with_vars_ns_muts<-all_genes[rowSums(unique_pair_muts_with_ns_muts) > 0 ]
    all_genes_with_vars<-all_genes[rowSums(unique_pair_muts) > 0 ] #the genes that are relevant to neighbor analysis SNPS 
    all_genes_with_vars_indels<-all_genes[rowSums(unique_pair_indels) > 0 ] #indels
    all_genes_with_vars_and_indels<-all_genes[rowSums(unique_pair_indels_muts) > 0]# indels and del_muts and stops
    all_genes_with_vars_ns_muts_and_indels<-all_genes[rowSums(unique_pair_indels_muts_with_ns_muts) > 0]
    
    #subset matrix to match
    #full dataset
    unique_pair_muts_with_ns_muts<-unique_pair_muts_with_ns_muts[rowSums(unique_pair_muts_with_ns_muts) > 0,]
    unique_pair_muts<-unique_pair_muts[rowSums(unique_pair_muts) > 0,] # subsetting the matrix based on genes that have aSNP
    unique_pair_indels<-unique_pair_indels[rowSums(unique_pair_indels) > 0,] #matrix based on genes with an indel
    unique_pair_indels_muts<-unique_pair_indels_muts[rowSums(unique_pair_indels_muts) > 0,] # indels,and del_muts and stops
    unique_pair_indels_muts_with_ns_muts<-unique_pair_indels_muts_with_ns_muts[rowSums(unique_pair_indels_muts_with_ns_muts) >0,]
    
    
    
    
    #KPC- dataset
    #genes and indels
    KPCminus_genes_with_vars<-KPC_minus_genes[rowSums(KPC_minus_unique_pair_muts) > 0 ] #the genes that are relevant to neighbor analysis SNPS 
    KPCminus_with_vars_indels<-KPC_minus_genes[rowSums(KPC_minus_unique_pair_indels) > 0 ] #indels 
    KPCminus_with_vars_indels_muts<-KPC_minus_genes[rowSums(KPC_minus_unique_pair_indels_muts) > 0 ] #indels & del muts and stops
    
    # KPC- matrix to match
    KPC_minus_unique_pair_muts<-KPC_minus_unique_pair_muts[rowSums(KPC_minus_unique_pair_muts) > 0,] # subsetting the matrix based on genes that have aSNP
    KPC_minus_unique_pair_indels<-KPC_minus_unique_pair_indels[rowSums(KPC_minus_unique_pair_indels) > 0,] #matrix based on genes with an indel
    KPC_minus_unique_pair_indels_muts<-KPC_minus_unique_pair_indels_muts[rowSums(KPC_minus_unique_pair_indels_muts) > 0,]
    
# get matrix with counts for genes with unique variants in pairs
    #full dataset
    #del_muts and ns muts
    all_unique_mut_ns_mut_test_mat<-get.var_mat_for_test(uniq_genes = all_unique_genes, pair1_mut_genes = all_genes_with_vars_ns_muts, pair1_mut_sub = unique_pair_muts_with_ns_muts, pairs_df = neighbor_pairs_df)
    #ns_muts
    all_unique_mut_test_mat<-get.var_mat_for_test(uniq_genes = all_unique_genes, pair1_mut_genes = all_genes_with_vars, pair1_mut_sub = unique_pair_muts, pairs_df = neighbor_pairs_df)
    #indels
    all_unique_indel_test_mat<-get.var_mat_for_test(uniq_genes = all_unique_genes, pair1_mut_genes = all_genes_with_vars_indels, pair1_mut_sub = unique_pair_indels, pairs_df = neighbor_pairs_df)
    #indels and del muts
    all_unique_indel_del_mut_test_mat<-get.var_mat_for_test(uniq_genes = all_unique_genes, pair1_mut_genes = all_genes_with_vars_and_indels, pair1_mut_sub = unique_pair_indels_muts, pairs_df = neighbor_pairs_df)
    # indels del_ muts and ns_muts
    all_unique_indel_mut_ns_mut_test_mat<-get.var_mat_for_test(uniq_genes = all_unique_genes, pair1_mut_genes = all_genes_with_vars_ns_muts_and_indels, pair1_mut_sub = unique_pair_indels_muts_with_ns_muts, pairs_df = neighbor_pairs_df)
    
    #KPC minus dataset
    #genes
    KPC_minus_unique_mut_test_mat<-get.var_mat_for_test(uniq_genes = KPC_minus_unique_genes, pair1_mut_genes = KPCminus_genes_with_vars, pair1_mut_sub = KPC_minus_unique_pair_muts, pairs_df = KPC_minus_neighbor_pairs_df)
    #indels
    KPC_minus_unique_indel_test_mat<-get.var_mat_for_test(uniq_genes = KPC_minus_unique_genes, pair1_mut_genes = KPCminus_with_vars_indels, pair1_mut_sub = KPC_minus_unique_pair_indels, pairs_df = KPC_minus_neighbor_pairs_df)
    # indels and del muts
    KPC_minus_unique_indel_and_del_mut_test_mat<-get.var_mat_for_test(uniq_genes = KPC_minus_unique_genes, pair1_mut_genes = KPCminus_with_vars_indels_muts, pair1_mut_sub = KPC_minus_unique_pair_indels_muts, pairs_df = KPC_minus_neighbor_pairs_df)
    
    
    
##########################PERFORM CONTINGENCY TESTS############################################################
###############################################################################################################
    
    # FOR FISHER.TEST USE MATRIX WITH BOTH INDELS AND DEL MUTS INCLUDED
    #all mutations carbR vs S fisher exact test
    all_mut_carbR_fisher_results<-carbR_fisher.test(unique_genes = all_unique_genes, positive_smaller_sub = all_unique_indel_del_mut_test_mat, carb_phenotype = carb_phenotype, annot_by_gene = all_annot_by_gene)
    
    # KPC - dataset, this makes more sense to do since we know why KPC+ genomes are resistant.
    KPC_minus_carbR_fisher_results<-carbR_fisher.test(unique_genes = KPC_minus_unique_genes, positive_smaller_sub = KPC_minus_unique_indel_and_del_mut_test_mat, carb_phenotype = KPC_minus_carb_phenotype, annot_by_gene = KPC_minus_annot_by_gene)
    
    #Blood vs non-blood
    #full datset
    all_mut_blood_fisher_results<-blood_fisher.test(unique_genes = all_unique_genes, positive_smaller_sub = all_unique_indel_del_mut_test_mat, carb_phenotype = carb_phenotype, annot_by_gene = all_annot_by_gene)
    
    
##########################PLOT FISHER.TEST RESULTS############################################################
###############################################################################################################    
    
    #Blood:
    blood_hits<-all_unique_indel_del_mut_test_mat[,all_mut_blood_fisher_results$p_value <.2]
    gene_annots_blood<-all_mut_blood_fisher_results$annotations[all_mut_blood_fisher_results$p_value <.2] #gene annotaitons for plotting
    colnames(blood_hits)<-gene_annots_blood
    
    #carbR:
    carb_hits<-KPC_minus_unique_indel_and_del_mut_test_mat[,KPC_minus_carbR_fisher_results$p_value <.25]
    gene_annots_carb<-KPC_minus_carbR_fisher_results$annotations[KPC_minus_carbR_fisher_results$p_value <.25] #gene annotaitons for plotting
    colnames(carb_hits)<-gene_annots_carb
    
    
# MAKE 
S
    
    # prepare tree for clustering heatmaps
    #subset tree to include only genomes in the matrix
    blood_tree<-drop.tip(tree, tree$tip.label[!(tree$tip.label %in% rownames(blood_hits))]) #tree with isolates corresponding to matrix
    carb_tree<-drop.tip(tree, tree$tip.label[!(tree$tip.label %in% rownames(carb_hits))]) #tree with isolates corresponding to matrix
    
    # sort the trace by the tree
    
    blood_hits<-blood_hits[blood_tree$tip.label,]
    carb_hits<-carb_hits[carb_tree$tip.label,]
    
    # genome phenotype annotations for plotting heatmaps
    blood_cols<-carb_phenotype$blood
    names(blood_cols)<-carb_phenotype$genome
    blood_cols<-blood_cols[blood_tree$tip.label] #order by tree
    
    #blood:
    blood_cols[blood_cols == TRUE]<-"red"
    blood_cols[blood_cols == FALSE]<-"blue"
    
    #carbR:
    carb_cols<-KPC_minus_carb_phenotype$carbR
    names(carb_cols)<-KPC_minus_carb_phenotype$genome
    carb_cols<-carb_cols[carb_tree$tip.label] #order by tree
    carb_cols[carb_cols == TRUE]<-"red"
    carb_cols[carb_cols == FALSE]<-"blue"
    
    
    # BLOOD HEATMAP
    file = paste(format(Sys.time(), "%Y-%m-%d"), 'blood_hits_fishertest.pdf');
    pdf(file, height = 12, width = 12)
    
    heatmap3(t(blood_hits), 
             scale="none",
             cexRow=.7,
             cexCol=.7,
             #Rowv= NA,
             Colv= NA,
             col=c("white",blue2red(9)),
             breaks=c(-1,0:9),
             ColSideColors =  blood_cols,
             ColSideLabs = "blood phenotype",
             main= "del stops indels blood p < .2 fisher.test")
    
    dev.off();
    
    #CARBAPENEM R HEATMAP
    file = paste(format(Sys.time(), "%Y-%m-%d"), 'carb_hits_fishertest.pdf');
    pdf(file, height = 12, width = 12)
    
    heatmap3(t(carb_hits), 
             scale="none",
             cexRow=.5,
             cexCol=.5,
             #Rowv= NA,
             Colv= NA,
             col=c("white",blue2red(14)),
             breaks=c(-1,0:14),
             ColSideColors =  carb_cols,
             ColSideLabs = "resistance phenotype",
             main= "del stops indels carbR vs S p < .25") 
    
    dev.off();
    
##########################PATHWAY ANALYSIS FOR BLOOD VS OTHER VARIANT CONTINGENCY TESTS########################
###############################################################################################################
    
   # Blood_pathway_poission_test<-function(positive_smaller_sub = matrix,
   #                                       positive_smaller_sub_indel = matrix,
   #                                       carb_phenotype = data.frame, 
   #                                       gene_loc = data.frame, 
   #                                       uniq_genes = vector,
   #                                       del_mut_space = data.frame,
   #                                       pathway_annot = data.frame)
    
  #pathway results with ns_muts
    
    all_blood_pathway_results_ns_muts<-Blood_pathway_poission_test(positive_smaller_sub = all_unique_mut_ns_mut_test_mat, 
                                                           positive_smaller_sub_indel = all_unique_indel_test_mat, 
                                                           carb_phenotype = carb_phenotype,
                                                           gene_loc = gene_loc,
                                                           uniq_genes = all_unique_genes,
                                                           del_mut_space = del_mut_space,
                                                           pathway_annot = pathway_annot)
    
  #no ns_muts  

  all_blood_pathway_results<-Blood_pathway_poission_test(positive_smaller_sub = all_unique_mut_test_mat, 
                                 positive_smaller_sub_indel = all_unique_indel_test_mat, 
                                 carb_phenotype = carb_phenotype,
                                 gene_loc = gene_loc,
                                 uniq_genes = all_unique_genes,
                                 del_mut_space = del_mut_space,
                                 pathway_annot = pathway_annot)
    
  write.table(all_blood_pathway_results,"blood_pathway_results.txt", sep=",", col.names = T)
 
  
  write.table(all_blood_pathway_results_ns_muts,"with_ns_muts_blood_pathway_results.txt", sep=",", col.names = T)
  
  
  #GET MATRIX WITH PATHWAYS, SHOWING WHAT GENOMES HAVE A HIT IN THAT PATHWAY         
    
    #get.pathway_matrix<-function(pathway_annot = data.frame, positive_smaller_sub = matrix 
    
   pathway_hits_mat<-get.pathway_matrix(pathway_annot = pathway_annot, positive_smaller_sub = all_unique_indel_test_mat , pathway_results = all_blood_pathway_results)
    
    
    pathway_hits_mat_ns_muts<-get.pathway_matrix(pathway_annot = pathway_annot, positive_smaller_sub = all_unique_indel_mut_ns_mut_test_mat , pathway_results = all_blood_pathway_results_ns_muts )
    
    
  # PLOT THE PATHWAY MATRIX
    # genome phenotype annotations for plotting heatmaps
    blood_cols<-carb_phenotype$blood
    names(blood_cols)<-carb_phenotype$genome
    #blood_cols<-blood_cols[trace_tree$tip.label] #order by tree
    
    
    #All pathways blood:
    blood_cols[blood_cols == TRUE]<-"red"
    blood_cols[blood_cols == FALSE]<-"blue"
    
    file = paste0(format(Sys.time(), "%Y-%m-%d"),'pathway_blood_vs_other.pdf');
    pdf(file, height = 12, width = 12)
    
    heatmap3(t(pathway_hits_mat), 
             scale="none",
             cexRow=.6,
             cexCol=.5,
             Rowv= NA,
             Colv= NA,
             col=c("white",blue2red(10)),
             breaks=c(-1,0:10),
             ColSideColors =  blood_cols,
             ColSideLabs = "blood phenotype",
             main= "pathway_hits_matrix blood vs other")
    dev.off();
    
    file = paste0(format(Sys.time(),"%Y-%m-%d"),'pathway_blood_vs_other_ns_muts.pdf');
    pdf(file, height = 12, width = 12)

    heatmap3(t(pathway_hits_mat_ns_muts), 
             scale="none",
             cexRow=.6,
             cexCol=.5,
             Rowv= NA,
             Colv= NA,
             col=c("white",blue2red(10)),
             breaks=c(-1,0:10),
             ColSideColors =  blood_cols,
             ColSideLabs = "blood phenotype",
             main= "pathway_hits_matrix with ns_muts blood vs other")
    dev.off();
    
    
    
    
    # subset on significant pathways
    pathway_hits<-pathway_hits_mat[,all_blood_pathway_results$two_sample_p_value <.05 & !is.na(all_blood_pathway_results$two_sample_p_value)]
    gene_annots_blood<-all_blood_pathway_results$pathway_name[all_blood_pathway_results$two_sample_p_value <.05 & !is.na(all_blood_pathway_results$two_sample_p_value)] #gene annotaitons for plotting
    colnames(pathway_hits)<-gene_annots_blood
    
     
    file = paste0(format(Sys.time(),"%Y-%m-%d"),'pathway_blood_vs_other_ns_mutsp_lessp1.pdf');
    pdf(file, height = 12, width = 12)
    
    heatmap3(t(pathway_hits), 
             scale="none",
             cexRow=.6,
             cexCol=.5,
             Rowv= NA,
             Colv= NA,
             col=c("white",blue2red(10)),
             breaks=c(-1,0:10),
             ColSideColors =  blood_cols,
             ColSideLabs = "blood phenotype",
             main= "pathway_hits_matrix with ns_muts blood vs other P < 0.1")
    dev.off();
    
    
  # MAKE HEATMAP WITH ALL PATHWAYS SORTED BY SIGNIFICANCE IN 2 SAMPLE POISSION TEST
    ##sort carb_phenotype by order of genomes in variant matrix
    #carb_phenotype<-carb_phenotype[order(carb_phenotype$genome[colnames(re_ref_indel_mat)]),] 
    #full set of pathways
    colnames(pathway_hits_mat_ns_muts)<-all_blood_pathway_results_ns_muts$pathway_name
    pathway_hits_mat_ns_muts<-pathway_hits_mat_ns_muts[,order(all_blood_pathway_results_ns_muts$two_sample_p_value)]
    
    all_blood_pathway_results_ns_muts<-all_blood_pathway_results_ns_muts[order(all_blood_pathway_results_ns_muts$two_sample_p_value),]
    
    file = paste0(format(Sys.time(),"%Y-%m-%d"),'pathway_blood_vs_other_all_ordered_pvalue.pdf');
    pdf(file, height = 12, width = 12)
    
    heatmap3(t(pathway_hits_mat_ns_muts), 
             scale="none",
             cexRow=.6,
             cexCol=.5,
             Rowv= NA,
             Colv= NA,
             col=c("white",blue2red(10)),
             breaks=c(-1,0:10),
             ColSideColors =  blood_cols,
             ColSideLabs = "blood phenotype",
             main= "pathway_hits_matrix with ns_muts blood vs other ordered by pvalue")
    dev.off();
    
    
       