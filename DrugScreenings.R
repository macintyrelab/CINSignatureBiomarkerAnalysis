
#Run trifecta_analysis prior to this script

#################################################
##### SCREENING OF DRUG RESPONSE BIOMARKERS #####
#####         & NOVEL DRUG TARGEST          #####
#################################################

##### LOAD PACKAGES AND ENVIRONMENT ##### 
#LIBRARIES
library(dplyr)
library(ggplot2)
library(ppcor)
library(ggalluvial)
library(ggrepel)
library(GGally)
library(geomnet)
library(ggnetwork)
library(igraph)
library(this.path)

#PATHS
BASE=dirname(this.path())
OUTFIGURES=file.path(BASE,"figures")
dir.create(OUTFIGURES, showWarnings = FALSE, recursive = TRUE)


#### LOAD DRUG DATA ######
#download from https://depmap.org/repurposing/
drug_screen<-read.csv(paste0(BASE,"/data/secondary-screen-dose-response-curve-parameters.csv"),stringsAsFactors = F)
drug_screen<-drug_screen[drug_screen$screen_id=="HTS002",]

drug_lookup<-c()
unique_target<-unique(drug_screen[,c("name","moa","target")])
for(i in 1:nrow(unique_target))
{
  drug<-unique_target[i,"name"]
  moa<-unique_target[i,"moa"]
  targets<-unlist(strsplit(unique_target[i,"target"],", "))
  drug_lookup<-rbind(drug_lookup,cbind(drug,moa,targets))
}

drug_lookup<-unique(drug_lookup)
drug_lookup<-data.frame(drug_lookup)
colnames(drug_lookup)<-c("drug","moa","target")
drug_lookup<-drug_lookup[!is.na(drug_lookup$target),]


#### LOAD SCREENING RESULTS #####
#Output files from trifecta_analysis.Rmd --> saved in the main folder
setwd(paste0(BASE,"/data"))
comb_drug<-read.table("comb_drug_response_qval.txt",sep="\t",header=T) #CRISPR+RNAi+drug
KO_drug<-read.table("KO_drug_response_qval.txt",sep="\t",header=T) #CRISPR+drug
KD_drug<-read.table("KD_drug_response_qval.txt",sep="\t",header=T) #RNAi+drug


#### CRISPR + RNAi + DRUG DATA #####
#List of genes by drugs (I have merged info related to mechanism of action)
comb_drug<-merge(comb_drug, drug_lookup, by=c("target", "drug"))
comb_drug$Analysis <- "Combined"
#List of significant genes by mechanism of action
comb_moa<-comb_drug[!duplicated(comb_drug[c(1,3,7)]),]

#### CRISPR + DRUG DATA #####
#List of genes by drugs (I have merged info related to mechanism of action)
KO_drug<-merge(KO_drug, drug_lookup, by=c("target", "drug"))
KO_drug$Analysis <- "CRISPR"
#List of significant genes by mechanism of action
KO_moa<-KO_drug[!duplicated(KO_drug[c(1,3,7)]),]

#### RNAi + DRUG DATA #####
#List of genes by drugs (I have merged info related to mechanism of action)
KD_drug<-merge(KD_drug, drug_lookup, by=c("target", "drug"))
KD_drug$Analysis <- "RNAi"
#List of significant genes by mechanism of action
KD_moa<-KD_drug[!duplicated(KD_drug[c(1,3,7)]),]


#### COMBINATION OF THREE LISTS #######
sig_moa <- rbind(comb_moa, KD_moa)
sig_moa <- rbind(sig_moa, KO_moa)
sig_moa <- sig_moa[!duplicated(sig_moa[c(1,3,7)]),] #104 hits
sig_moa$signature <- paste("CX", sig_moa$signature, sep = "")

#Save
#write.table(sig_moa, file="sig_moa.txt", sep="\t", row.names=F)


#### EXPLORING MOA GROUPS ####
#MOA list is manually explored in order to perform groups --> sig_moa_groups.txt is created
#Load MOA groups
sig_moa_group<-read.table("sig_moa_groups.txt",sep="\t",header=T) #MOA groups
sig_moa_group$group_moa[sig_moa_group$group_moa=="Microtubule inhibitor "] <- "Microtubule inhibitor"

#CX7 is removed since we don't have any information about the putative aetiology
sig_moa_group_noCX7 <- sig_moa_group[sig_moa_group$signature != "CX7",]


##### ALLUVIAL PLOT DRUG RESPONSE BIOMARKERS (PANEL B) #####
#Prepare data for alluvial plot --> Group by CS + Drug + MOA group
d <- sig_moa_group_noCX7 %>% 
  group_by(signature, drug, group_moa) %>% 
  summarise(n = n())

d$signature <- factor(d$signature, levels=c("CX1", "CX2", "CX3", "CX4", "CX5", "CX9",
                                            "CX10","CX11","CX13","CX14","CX15","CX17"))

#Optimize the order of moa groups in the right column in order to minimize crossing connection lines in the Alluvial plots
ggplot(data = sig_moa_group_noCX7, aes(from_id = group_moa, to_id = signature)) +
  geom_net(aes(colour = as.character(signature)), layout.alg = "kamadakawai", 
           size = 2, labelon = TRUE, vjust = -0.6, ecolour = "grey60",
           directed =FALSE, fontsize = 3, ealpha = 0.5)

#Order levels of moa according to network groups
d$group_moa <- factor(d$group_moa, levels=c("RAF inhibitor", "Microtubule inhibitor", "Neurotransmission inhibitor", 
                                            "Guanylate cyclase stimulant", "Immunosuppressant", "Progesterone receptor antagonist",
                                            "JAK inhibitor", "Microtubule disrupting agent", "Anti-infective agent", "PARP inhibitor",
                                            "Microtubule stabilizing agent", "PKC activator",
                                            "Multikinase inhibition", "Neuroprotective agent", 
                                            "KIT inhibitor", "PDGFR tyrosine kinase receptor inhibitor", "HDAC inhibitor", 
                                            "ALK tyrosine kinase receptor inhibitor", "MET inhibitor", "PI3K inhibitor",
                                            "EGFR inhibitor","Apoptotic agent", "CDK inhibitor", "PKC inhibitor"))
d <- d[order(d$group_moa),]

d$drug <- factor(d$drug, levels=c("AZ-628", "cabazitaxel", "paclitaxel", "colchicine",
                                    "mebendazole", "vinorelbine", "thiocolchicoside", "epinephrine", "maprotiline",
                                    "methyldopa", "CGS-15943", "sertindole", "riociguat", "kynurenic-acid","gestrinone",
                                    "AZD1480", "ruxolitinib","filgotinib", "2-methoxyestradiol", "puromycin",
                                    "olaparib","ixabepilone","epothilone-d","ingenol-mebutate", "imatinib",
                                    "sunitinib" ,"ponatinib","RAF265" , "regorafenib" , "dasatinib",
                                    "vandetanib", "salidroside", "masitinib", "pazopanib", "crenolanib",
                                    "givinostat", "JNJ-26481585", "AP26113", "LY2801653" , "taselisib",
                                    "afatinib", "PETCM", "arcyriaflavin-a", "bisindolylmaleimide-ix"))


#Alluvial plot
setwd(OUTFIGURES)
col <- c(c("CX1"="red1", "CX2"="grey50", "CX3"="royalblue3", "CX4"="skyblue1", "CX5"="tan2", "CX9"="gold4",
           "CX10"="pink","CX11"="red4","CX13"="palegreen3","CX14"="plum4","CX15"="violet","CX17"="gold1"))

#Save as pdf
pdf("AlluvialPlot_DrugBiomarkers.pdf", width = 8, height = 12)
pb <- ggplot(d, aes(axis1 = drug, axis2 = signature)) +
  geom_alluvium(aes(fill = signature), width = 1/20) +
  geom_stratum(width = 1/20, fill = "white", color = "grey30") +
  geom_text(stat = "stratum", colour="black", size=3, aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Drug", "CIN Signature"),
                   expand = c(.3, .3)) +
  scale_fill_manual(values = col) +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(colour = "black", size = 14))+
  theme(axis.text.y = element_blank())+
  theme(title = element_blank())+
  theme(axis.line=element_blank())+
  theme(panel.background=element_blank())+
  theme(panel.border=element_blank())+
  theme(plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "cm"))
pb
dev.off()


##### ALLUVIAL PLOT TARGETS #####
#We manually inspected the tractability of targets that have no targeted therapies.
#Target tractability was explored using canSAR database (https://cansarblack.icr.ac.uk/)
#We also selected those genes that have clear implication in CIN --> genes_nodrug_selected.txt is created

#Load data
setwd(paste0(BASE,"/data"))
sig_target_group<-read.table("genes_nodrug_selected.txt",sep="\t",header=T) #target groups

sig_target_group$signature <- factor(sig_target_group$signature, levels=c("CX1", "CX2", "CX3", "CX4", "CX5", "CX6", "CX9",
                                            "CX10","CX11","CX13","CX14","CX15","CX17"))

#Optimize the order of target groups in the right column in order to minimize crossing connection lines
sig_target_group$GroupGenes[sig_target_group$GroupGenes=="DNA repair (NER)"] <- "DNA repair"
ggplot(data = sig_target_group, aes(from_id = GroupGenes, to_id = signature)) +
  geom_net(aes(colour = as.character(signature)), layout.alg = "kamadakawai", 
           size = 2, labelon = TRUE, vjust = -0.6, ecolour = "grey60",
           directed =FALSE, fontsize = 3, ealpha = 0.5)


#Order levels of target groups according to network groups
sig_target_group$GroupGenes <- factor(sig_target_group$GroupGenes, levels=c("SWI/SNF", "Telomere maintenance", 
                                                                            "DNA repair", "Histone modification",
                                                                            "ROS stress response", "DNA damage and replication stress",
                                                                            "Cytoeskeletal estructure", "Chromatin architecture", "Cell cycle checkpoints",
                                                                            "Rho GTPase"))

sig_target_group <- sig_target_group[order(c(sig_target_group$GroupGenes,sig_target_group$signatures)),]
ord <- unique(sig_target_group$genes)
sig_target_group$genes <- factor(sig_target_group$genes, levels=ord)


#Alluvial Plot
setwd(OUTFIGURES)
col <- c(c("CX1"="red1", "CX2"="grey50", "CX3"="royalblue3", "CX4"="skyblue1", "CX5"="tan2", "CX6" = "springgreen4",
           "CX9"="gold4","CX10"="pink","CX11"="red4","CX13"="palegreen3","CX14"="plum4","CX15"="violet","CX17"="gold1"))

#Save as pdf
pdf("AlluvialPlot_DrugTargets.pdf", width = 8, height = 12)
pt <- ggplot(sig_target_group, aes(axis1 = signature, axis2 = genes)) +
  geom_alluvium(aes(fill = signature), width = 1/20) +
  geom_stratum(width = 1/20, fill = "white", color = "grey30") +
  geom_text(stat = "stratum", colour="black", size=3, aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("CIN Signature", "Target"),
                   expand = c(.3, .3)) +
  scale_fill_manual(values = col) +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(colour = "black", size = 14))+
  theme(axis.text.y = element_blank())+
  theme(title = element_blank())+
  theme(axis.line=element_blank())+
  theme(panel.background=element_blank())+
  theme(panel.border=element_blank())+
  theme(plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "cm"))
pt
dev.off()


#### ARRANGE BOTH ALLUVIAL PLOTS #####
## Figure 3a includes a diagram of the procedure --> it has been created using affinity designer
## Figure 3b includes both alluvial plots --> both plots has been maerged and modified using affinity designer

