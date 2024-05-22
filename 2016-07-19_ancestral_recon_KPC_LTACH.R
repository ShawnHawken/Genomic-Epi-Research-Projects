#FIRST ANALYSIS OF GEOGRAPHIC (facility location) DISSEMINATION OF KPC IN LTACH North

library(ape);
library(phytools);
library(RColorBrewer)
library(outbreaker)
#library(Heatplus)
library(heatmap3)
library(ggplot2)
library(reshape2)
library(igraph)
library(heatmap.plus)
library(igraph)
library(xlsx)
library(XLConnect)
library(ggplot2)
library(wesanderson)

source('../../Penn Analysis/epi_network_functions_v4.R')
setwd("/Users/shawnwhitefield/Desktop/Snitkin_lab/KPC_LTACH/LTACH_north_trace_analysis/Rush KPC Analysis For Shawn/Rush KPC LTACH Trace Analysis")

#####READ IN DATA SIMPLE######

##READ IN SEQUENCE DATA AND MAKE TREE
evan_sample_key = read.xlsx('../data/2015-9-19_Rush_KPC_LTACH_isolates.xlsx', sheetIndex = 1);
#fix date format
evan_sample_key_csv<-read.csv("../data/2015-9-19_Rush_KPC_LTACH_isolates.csv", sep=",", header=T)
evan_sample_key$Culture.Date<-evan_sample_key_csv$Culture.Date
rm(evan_sample_key_csv)
evan_sample_key$Culture.Date<-as.Date(evan_sample_key$Culture.Date,format = "%m/%d/%y")

#dna = read.dna('../data/2015-12-13_Rush_KPCallref_filter2_noUnmapped_noFiltered.fa',
#                format = 'fasta')
load('../data/2015-12-13_Rush_KPCallref_filter2_noUnmapped_noFiltered.Rdata'); #this dnabin is in Rdata?

dna_labels = gsub("Rush_KPC_|Sample_", "", row.names(dna));
dna_labels = gsub('_R1_', '', dna_labels);
dna_labels = gsub('_R_', '', dna_labels);
dna_labels = gsub('_R_', '', dna_labels);
dna_labels = gsub('__', '', dna_labels);
dna_labels[grepl("KPNIH1", dna_labels)] = "KPNIH1";

dna_labels = gsub('43714', '1', dna_labels)
dna_labels = gsub('43715', '2', dna_labels)
dna_labels = gsub('43716', '3', dna_labels)
dna_labels = gsub('43717', '4', dna_labels)
dna_labels = gsub('43718', '8', dna_labels)
dna_labels = gsub('43719', '14', dna_labels)
dna_labels = gsub('43720', '18', dna_labels)
dna_labels = gsub('43721', '19', dna_labels)
dna_labels = gsub('43722', '29', dna_labels)
dna_labels = gsub('43723', '30', dna_labels)

row.names(dna) = dna_labels;

dna_pt_lookup = sapply(dna_labels, FUN = function(x){evan_sample_key$Unique.Patient.Identifier[evan_sample_key$Sample.ID == x]});
dna_pt_labels = unlist(sapply(dna_pt_lookup, FUN = function(x){ if(length(x) == 1){as.character(x[1])}else{'NA'}} ))


#IDENTIFY VARIABLE POSITIONS AND SUBSET SEQEQUNCE
dna_var_pos = apply(dna, 2, FUN = function(x){sum(x==x[1]) < nrow(dna)})
dna_var = dna[, dna_var_pos];


#DETERMINE PAIRWISE SNP DISTANCES
dna_dist = dist.dna(dna_var, model = "raw");
dna_dist_JC = dist.dna(dna_var, model = "JC69");

snp_dist = as.matrix(dna_dist) * ncol(dna_var);
snp_dist_1kSNPmax = snp_dist;
snp_dist_1kSNPmax[snp_dist_1kSNPmax > 1000] = 1000;  # for plotting

snp_dist = as.matrix(dna_dist) * ncol(dna_var);
tree = nj(as.matrix(dna_dist_JC))

#CREATE TRACE MATRIX
#READ IN TRACE AND CULTURE DATA AND PUT IN PROPER FORMAT
surv_file = '../data/UPDATED2_North_ADM & PPS culture results_master KPC.LTACH db.xlsx';
trace_file = '../data/UPDATED_MASTER LIST_KPC(+) bed trace_clean_w study ID.xlsx';

trace_data = read.xlsx2(trace_file, 2, header = TRUE, check.names = FALSE, endRow = 315)
surv_data_full = read.xlsx2(surv_file, 1, header = TRUE, colIndex = 1:3);
surv_data = surv_data_full[surv_data_full[,2] %in% colnames(trace_data), ];

trace_mat_raw = as.matrix(trace_data[, 3:ncol(trace_data)]);

#remove non facility rooms
#determine facility locations
unique_obs_trace_mat<-as.data.frame(sort(unique(c(trace_mat_raw)))) #unique list of locations 
out_of_facility<-c("AMA","discharged","DTLT","DTRF","EXP","HHS","HOME","HOSP","HSMF","SNF","TEMP","TMPLOA1"," ")
trace_mat_raw[trace_mat_raw %in% out_of_facility]<-"" #make these 0

pts = trace_data$StudyID;
dates = as.Date(as.numeric(colnames(trace_data)[3:ncol(trace_data)]), origin = '1899-12-30"')
trace_mat = matrix(data = as.numeric(trace_mat_raw != ""), nrow = length(dates), ncol = length(pts), dimnames = list(as.character(dates), as.character(pts)), byrow = TRUE)
trace_mat2 = matrix(data = as.numeric(trace_mat_raw != ""), nrow = length(dates), ncol = length(pts), dimnames = list((dates), as.character(pts)), byrow = TRUE)

#ADD IN POSITIVE AND NEGATIVE SURVEILLANCE
pos_surv_dates = as.Date(as.numeric(as.character(surv_data[grepl("^P", surv_data[, 3]), 2])), origin = '1899-12-30')  
pos_surv_pts = surv_data[grepl("^P", surv_data[, 3]), 1]  

neg_surv_dates = as.Date(as.numeric(as.character(surv_data[grepl("^N", surv_data[, 3]), 2])), origin = '1899-12-30')  
neg_surv_pts = surv_data[grepl("^N", surv_data[, 3]), 1]  

dna_w_trace = dna_labels[as.character(dna_pt_labels[dna_labels]) %in% colnames(trace_mat)]
#upgma_tree = midpoint.root(upgma(as.matrix(dna_dist_JC)[dna_w_trace, dna_w_trace]));

#FIND THE DATE ROOM AND FLOOR AND WARD WHERE PATIENTS CONVERTED
#MAKE DIFFERENT LEVELS OF LOCATION INFORMATION EG DIFFERENT TRACE_MATS

t_trace_mat_raw<-t(trace_mat_raw) #transpose raw locaiton matrix so it matches trace matrix

#FLOOR
floor_trace<-trace_mat
floor_trace[grepl("*^1", t_trace_mat_raw)]<-5 #floor 1 = 5
floor_trace[grepl("*^2", t_trace_mat_raw)]<-6 #floor 2 = 6
floor_trace[grepl("*^3", t_trace_mat_raw)]<-7 #floor 3 = 7
floor_trace[grepl("*^4", t_trace_mat_raw)]<-8 #floor 4 = 8
floor_trace[grepl("*^S", t_trace_mat_raw)]<-9 #SCU = 9
floor_trace[grepl("*^H", t_trace_mat_raw)]<-9 #HAU = 9 #need 5 cols
floor_loc_mat<-floor_trace

#add surveillence cultures back in
#add surveillence cultures back in
floor_trace[cbind(as.character(pos_surv_dates), as.character(pos_surv_pts))] = 1.5  #assigning 1.5 to positive surveillence culture
floor_trace[cbind(as.character(neg_surv_dates), as.character(neg_surv_pts))] = 1.25 #assigning 1.25 to negative surveillence culture, all others are 0

#Ward trace matrix
#FLOOR REGION
#read in floor region data
floor_region_key = read.xlsx('../data/LTACH_room_lookup.xlsx', sheetIndex = 1);

#floor_region_mat
floor_region_mat<-trace_mat
t_trace_mat_raw<-t(trace_mat_raw) #transpose raw locaiton matrix so it matches trace matrix

for(room in as.vector(floor_region_key$unique_rooms)){
  room_loc<-t_trace_mat_raw == room
  floor_region_mat[room_loc]<-floor_region_key$plot_number[floor_region_key[,1]==room]
}

floor_region_loc_mat<-floor_region_mat

#add the surveillence cultures back in
floor_region_mat[cbind(as.character(pos_surv_dates), as.character(pos_surv_pts))] = 1.5  #assigning 1.5 to positive surveillence culture
floor_region_mat[cbind(as.character(neg_surv_dates), as.character(neg_surv_pts))] = 1.25 #assigning 1.25 to negative surveillence culture, all others are 0


#SUBSET SAMPLE KEY LOOKUP TABLE for 1. isolates that were sequenced
#                                   2. patients in the trace (e.g. not 2008)

#SUBSET SAMPLE KEY BASED ON SEQUENCED GENOMES
sequenced_samples<-evan_sample_key$Sample.ID[evan_sample_key$Sample.ID %in% colnames(snp_dist)]
sample_key<-evan_sample_key[evan_sample_key$Sample.ID %in% sequenced_samples,] #samples that are sequenced

#SUBSET SAMPLE KEY BASED ON PATIENTS BEING IN TRACE (remove 2008 isolates and patients)
patient_in_floor_trace<- sample_key$Unique.Patient.Identifier %in% colnames(floor_trace) #patients in floor trace
sample_key<-sample_key[patient_in_floor_trace,]


#DETERMINE WHERE A PATIENT ISOLATE CAME FROM AND ADD THIS INFORMATION TO THE SAMPLE KEY

sample_key$floor<-NA #set up columns to hold floor and ward locations
sample_key$ward<-NA
sample_key$pt_convert<-NA #to hold converts vs not. logical

#change the dates of patients whose culture date is after the end of the trace to the last date in the trace
#   these cultures likely grew later than the final day of the study but
#   were likely collected during the time the patient was in their final location in the trace
sample_key$Culture.Date[sample_key$Culture.Date == "2013-06-25"] <-"2013-06-18"


# DETERMINE WHICH PATIENTS WERE CONVERTS VS INTAKE POSITIVE

trace_converts = apply(floor_trace, 2, FUN = function(x){ neg = which(x == 1.25)
first_pos = min(which(x == 1.5))
last_neg = max(neg[neg < first_pos])
if(!is.infinite(first_pos) && !is.infinite(last_neg))
{
  
  (first_pos > last_neg) && sum(x[last_neg:first_pos] >=1) ==  (first_pos - last_neg + 1);
  
  #}else if(!is.infinite(first_pos) && is.infinite(last_neg)){ #CHANGED BECAUSE THERE ARE SOME PATIENTS IN TRACE WITHOUT A POSITIVE
}else{
  
  
  FALSE;
  
}
} )

#add the convert, and location data to the sample_key df

for(rows in 1:nrow(sample_key)){
  patient<-sample_key[rows,2]
  culture_date<-as.character(sample_key[rows,3])
  if(! is.na(culture_date)){  #SKIP THE NEXT 4 LINES IF ITS NA DATE
    floor_location<-floor_loc_mat[rownames(floor_loc_mat) == culture_date,colnames(floor_loc_mat) == patient] #floor
    ward_location<-floor_region_loc_mat[rownames(floor_region_loc_mat) == culture_date, colnames(floor_region_loc_mat) == patient]
    sample_key[rows,29]<-floor_location
    sample_key[rows,30]<-ward_location
    sample_key[rows,31]<-trace_converts[names(trace_converts) == patient]
    
  }
}

# ASSIGN ALL OF THE PATIENTS WHO WERE NOT CONVERTS TO WARD/FLOOR 18
sample_key$floor[sample_key$pt_convert ==FALSE ]<-18
sample_key$ward[sample_key$pt_convert == FALSE ]<-18

# GET COLOR CODE FOR EACH LOCATION
#named vector where names are isolates and values are the colors
#floor
floors<-as.character(sample_key$floor)
names(floors)<-sample_key$Sample.ID
floors[is.na(floors)]<-20

#floor_colors<-sample_key$floor
#names(floor_colors)<-sample_key$Sample.ID #7 colors, 6 "floors", 1 NA
#floor_colors[floor_colors == 5]<-"#85D4E3"
#floor_colors[floor_colors == 6]<-"#F4B5BD"
#floor_colors[floor_colors == 7]<-"#9C964A"
#floor_colors[floor_colors == 8]<-"#CDC08C"
#floor_colors[floor_colors == 9]<-"#798E87" #SCU #grey
#floor_colors[floor_colors == 18]<-"red" #intake positive
#floor_colors[is.na(floor_colors)]<-"black"

#named color vector for plotting
uniq_floors<-sort(unique(floors))
floor_cols<-(c("#85D4E3","#F4B5BD","#9C964A","#CDC08C","#798E87","red","black"))
names(floor_cols)<-uniq_floors

#ward
ward<-sample_key$ward
names(ward)<-sample_key$Sample.ID
ward[is.na(ward)]<-20

#ward_colors<-sample_key$ward
#names(ward_colors)<-sample_key$Sample.ID #14 ward colors including a NA ward (1 patient isolate without a culture date)
#ward_colors[ward_colors == 5]<-"#85D4E3" #blues
#ward_colors[ward_colors == 6]<-"#9FCCDA"
#ward_colors[ward_colors == 7]<-"#E2B9C3" #pinks
#ward_colors[ward_colors == 8]<-"#EFB3B6"
#ward_colors[ward_colors == 9]<-"#DCAC9E"
#ward_colors[ward_colors == 10]<-"#A19A50" #greens
#ward_colors[ward_colors == 11]<-"#ABA35E"
#ward_colors[ward_colors == 12]<-"#B5AC6C"
#ward_colors[ward_colors == 13]<-"#C0B47A" #strange rooms that arent on floor plan prob psych unit?
#ward_colors[ward_colors == 14]<-"#E7CD82"
#ward_colors[ward_colors == 15]<-"#F0D27E"
#ward_colors[ward_colors == 16]<-"#FAD77B"
#ward_colors[ward_colors == 17]<-"#798E87" #SCU grey
#ward_colors[ward_colors == 18]<-"red"
#ward_colors[is.na(ward_colors)]<-"black"

#named color vector for plotting
uniq_ward<-sort(unique(ward))
ward_cols<-(c("#85D4E3","#9FCCDA","#E2B9C3","#EFB3B6","#DCAC9E",
              "#A19A50","#ABA35E","#B5AC6C","#E7CD82",
              "#F0D27E","#FAD77B","#798E87","red","black"))
names(ward_cols)<-uniq_ward


#SUBSET TREE BASED ON THE TRACE
tips_to_drop<-!(tree$tip.label %in% sample_key$Sample.ID) 
tree_sub<-drop.tip(tree,tip = tree$tip.label[tips_to_drop])


# CHANGE TIP LABEL ON TREE TO PATIENT CODE RATHER THAN ISOLATE CODE SO WE KNOW WHAT ISOLATES CAME FROM SAME PATIETN
#CHANGE TIP LABEL TO PATIENT CODE

tree_sub$tip.label = sapply(tree_sub$tip.label, FUN = function(x) { if(sum(x %in% sample_key$Sample.ID) > 0)
{
  sample_key$Unique.Patient.Identifier[sample_key$Sample.ID == x];
}else{
  
  x;
  
}#end if
});

#PLOT TREE
#PHYLOGRAM, OR RADIAL, OR CLADOGRAM
  #floor
plot(tree_sub, cex = .15, type = "phylogram", label.offset = 2, open.angle = 15, no.margin = TRUE, , root.edge = TRUE, rotate.tree = -90, use.edge.length = FALSE)

tiplabels(pie = to.matrix(floors, names(floor_cols)),
         piecol = floor_cols, cex = 0.15)

legend('topleft', legend = names(floor_cols), col = "black", pt.bg = floor_cols, pch = 21)

  #ward
plot(tree_sub, cex = .15, type = "phylogram", label.offset = 2, open.angle = 15, no.margin = TRUE, , root.edge = TRUE, rotate.tree = -90, use.edge.length = FALSE)

tiplabels(pie = to.matrix(ward, names(ward_cols)),
          piecol = ward_cols, cex = 0.04)

legend('topleft', legend = names(ward_cols), col = "black", pt.bg = ward_cols, pch = 21)



#SUBSET ONLY ST258
d_mat = cophenetic.phylo(tree_sub);
not_ST258<-c(row.names(d_mat)[which(rowMeans(d_mat) > 0.2)])
tree_root = root(tree_sub, '128') #root by an arbitrary ST258
tree_ST258 = root(drop.tip(tree_root, not_ST258), '128')


#PREPARE PIE COLORS
ST258_floors<-as.character(sample_key$floor[sample_key$Sample.ID %in% intersect(tree_ST258$tip.label, sample_key$Sample.ID)])

names(ST258_floors)<-sample_key$Sample.ID[sample_key$Sample.ID %in% intersect(tree_ST258$tip.label, sample_key$Sample.ID)]

ST258_uniq_floors<-sort(as.numeric(unique(ST258_floors)))
ST258_uniq_floors[is.na(ST258_uniq_floors)]<-20
ST258_floor_cols<-(c("#85D4E3","#F4B5BD","#9C964A","#CDC08C","#798E87","red","grey"))
names(ST258_floor_cols)<-ST258_uniq_floors



#GET LTACH FOR EACH SAMPLE IN TREE
plot(tree_ST258, cex = .15, type = "phylogram", label.offset = 2, open.angle = 15, no.margin = TRUE, , root.edge = TRUE, rotate.tree = -90, use.edge.length = FALSE)

tiplabels(pie = to.matrix(ST258_floors, names(ST258_floor_cols)),
          piecol = ST258_floor_cols, cex = 0.15)



#PLOT TREE WITH ANCESTRAL STATES
floor_fitER <- rerootingMethod(tree_ST258, ST258_floors, model = "ER")







##################################################################################################################
# EXTRA STUFF ABOUT CONVERSIONS THAT MIGHT WANT LATER 
##################################################################################################################
#floors
floor_trace_sub<-floor_trace[,colnames(floor_trace) %in% sample_key_sub$Unique.Patient.Identifier]
floor_trace_sub[sample_key_sub$Unique.Patient.Identifier %in% colnames(floor_trace_sub)

convert_floors<-sapply(1:ncol(floor_trace_sub),function(x){
  first_pos<-min(which(floor_trace_sub[ ,x]==1.5))
  first_pos<-loc_mat[first_pos,x]
}) #warnings are ok these are patients who never convert


#floor region (ward)
convert_floor_region<-sapply(1:ncol(floor_region_mat), function(x){
  first_pos<-min(which(floor_region_mat[ ,x]==1.5))
  first_pos<-floor_region_mat_loc[first_pos,x]
}) #warnings are ok these are patients who never convert


#add patient IDs back to convert floors
convert_floor_df<-cbind(convert_floors,colnames(floor_trace)) #floor
colnames(convert_floor_df)<-c("convert_floor","Unique.Patient.Identifier")

convert_floor_region_df<-cbind(convert_floor_region, colnames(floor_region_mat)) #floor region (ward)
colnames(convert_floor_region_df)<-c("convert_ward","Unique.Patient.Identifier")

convert_loc_df<-merge(convert_floor_df,convert_floor_region_df, by = "Unique.Patient.Identifier") #DF WITH PATIENTS AND FLOOR AND WARD OF CONVERSION

#MERGE CONVERSION INFO WITH SAMPLE KEY
sample_key<-merge(evan_sample_key,convert_loc_df) #this only contains patients who had a conversion or intake positive

#make histograms of ward intakes
convert_loc_df$convert_floor<-as.numeric(as.vector(convert_loc_df$convert_floor))#fix factor levels for plotting
convert_loc_df$convert_ward<-as.numeric(as.vector(convert_loc_df$convert_ward))

#make hist of conversions by ward
hist(convert_loc_df$convert_ward[!trace_converts],breaks=c(5.00, 6.00,  7.00,  8.00,  9.00,10.00, 11.00, 12.00, 13.00, 14.00, 15.00, 16.00, 17.00),
    col="steelblue",labels=T, freq = T)
par(new=TRUE)
hist(convert_loc_df$convert_ward[trace_converts], col="lightpink")
par(new=TRUE)
hist(convert_loc_df$convert_ward[!trace_converts], col="red")

#conversions by floor
hist(convert_loc_df$convert_floor,breaks=c(-1,0,1,1.25,1.5,5,6,7,8,9,18))


# ASSIGN ALL PATIENTS WHO WERE INTAKE POSITIVE TO A "WARD" 18
convert_loc_df$convert_floor[!trace_converts]<-18 #floor
convert_loc_df$convert_ward[!trace_converts]<-18 #ward

#######

