## Gene oncoprint script
## Author: Lolita Lecompte

## import packages
require(tidyverse)
require(ComplexHeatmap)
require(RColorBrewer)

## import data
mat <- readRDS("matrix.rds")

fq <- apply(mat, 1, function (x) sum((x != ""), na.rm = TRUE) / ncol(mat))
fqOrdered <- fq[order(fq, decreasing = TRUE)]
newmat <- mat[which(fq >= 0.03),]
new_order <- names(fqOrdered[(fqOrdered >= 0.03)])

################################################################################
## OncoKB annotation
oncogene <- read.csv("oncoKB_cancerGeneList.tsv", header = TRUE, sep = "\t")
oncogene <- mutate(oncogene, class = case_when(
  Is.Oncogene == "Yes" & Is.Tumor.Suppressor.Gene == "Yes"  ~ "Both",
  Is.Oncogene == "Yes" & Is.Tumor.Suppressor.Gene == "No"  ~ "Oncogene",
  Is.Oncogene == "No" & Is.Tumor.Suppressor.Gene == "Yes"  ~ "TSG",
  Is.Oncogene == "No" & Is.Tumor.Suppressor.Gene == "No"  ~ "Unknown",
))

oncogene <- oncogene %>% mutate(class = ifelse(Hugo.Symbol %in% c("LRP1B"), "TSG", class))
oncogene$class <- factor(oncogene$class, levels = c("Oncogene", "TSG", "Both","Unknown"))
oncogeneTSG <- oncogene %>% filter(class != "Unknown") %>% pull(Hugo.Symbol)

oncogenemat <- newmat[which(rownames(newmat) %in% oncogeneTSG),]
oncogene_order <- names(fqOrdered[which(names(fqOrdered) %in% rownames(oncogenemat))])

################################################################################
## Clinical data
sclinical <- read.csv("clinical_data.csv", row.names = 1) 

################################################################################
## Table 1. WES N = 291
sclinical %>% 
  mutate(Age.cat = factor(ifelse(Age <= 50, "<=50",">50"))) %>% 
  mutate(figo.cat = factor(ifelse(Figo.18..4.classes. == "I" | Figo.18..4.classes. == "II", "I/II","III/IV"))) %>%
  mutate(HPV.pos.neg = factor(ifelse(HPV == "NA","NA",
                                     ifelse(HPV =="NEG", "NEG","POS")))) %>% 
  mutate(statut.PFS = factor(statut.PFS)) %>% 
  mutate(Type.histo.CRF = factor(Type.histo.CRF)) %>% 
  select(Age.cat, figo.cat, Type.histo.CRF, HPV.pos.neg, necrose.cat, TILs.cat,  treatment, statut.PFS) %>% summary()

################################################################################
## Oncoprint colors
col = c("background"="#CCCCCC",
        "frameshift indel" ="orange", 
        "inframe indel"="yellowgreen", 
        "splicing"="hotpink", #"pink", 
        unknown_splice = "darkorchid4",
        "missense"="forestgreen", 
        synonymous_variant="mediumorchid2",
        "stop gain"="gold", 
        start_lost = "yellow2", 
        stop_lost = "red",
        "unknown"="darkgrey",
        DEL = 'deepskyblue',
        AMP = 'red3',
        'TERT' = 'dodgerblue4', 
        'alt' = 'grey38', 
        fusion = "black")

################################################################################
## Heatmap annotation
column_ha = HeatmapAnnotation(
  `Histological type` = sclinical[colnames(mat), "type.histo.4.classes"],
  HPV = sclinical[colnames(mat), "cat.HPV"],
  `Figo stage` = sclinical[colnames(mat), "Figo.18..4.classes."],
  Necrosis = sclinical[colnames(mat),"necrose.cat"],
  TILs = sclinical[colnames(mat),"TILs.cat"],
  Treatment = sclinical[colnames(mat), "treatment"],
  `Mutational signature` = sclinical[colnames(mat), "Signature"],
  TMB = sclinical[colnames(mat),"TMB"],
  MSI = sclinical[colnames(mat),"MSI"],
  HRD = sclinical[colnames(mat),"shallowHRD.status"],

  cbar = anno_oncoprint_barplot(height = unit (2, "cm")),
  
  gp = gpar(width = unit(14, "in")),
  
  col = list(`Histological type` =
               c("Adenocarcinoma" = rev(brewer.pal(3, "Greens"))[1] ,
                 "Adenosquamous" = rev(brewer.pal(3, "Greens"))[2],
                 "Squamous" = rev(brewer.pal(3, "Greens"))[3],
                 "Unknown/Clear cell"="grey"),
             
             HPV = 
               c("NEG"="#7d458d",
                 "16"="#ffe6f2", 
                 "18"="#b95990", 
                 "Other"="#ed9ab4",
                 "Unknown"="grey"),
             
             `Figo stage` =
               c("I" = "#d90000", 
                 "II" = "#f1701e", 
                 "III" = "#ffb366",
                 "IV" = "#fef0d9"),
             
             Necrosis = c("Presence" = "#b96b35", "Absence"="ivory1", "Unknown"="grey"),
             TILs = c("Presence" = "#90befa", "Absence" = "ivory1","Unknown"="grey"), 
             Treatment =
               c("Chemoradiation" = "#ffeded", 
                 "Neoadjuvant chemotherapy" = "#f38c75",
                 "Surgery" = "#cb0000"),
             
             `Mutational signature` = 
               c("AID"="#661100", "APOBEC"="#f5d487","Artifact"="#DDCC77", "Aza"="#999933","Chemotherapy/Treatment"="#117733","Deamin_5MC"="#44AA99","HRD"="#88CCEE",
                 "MMR"="#6699CC","NTHL1"="#332288","POL"="#882255","ROS"="#AA4499", "Tobacco"="#CC6677","Unknown"="ivory","UV"="#888888"),
             
             TMB = c("High" = "chartreuse4", "Low" = "ivory"), 
             MSI = c("High" = "orange", "Low" = "ivory"),
             HRD = c("HRD"="red2", "HRP"="ivory"),
             PFS = c("PFS >= 12 months" = "#DEEBF7", "PFS < 12 months" = "#3182BD")
  ))

heatmap_legend_param = list(title = "Alterations", 
                            at = c("AMP", "DEL", "TERT", "stop gain", "frameshift indel", "inframe indel", "missense","splicing","fusion"), 
                            labels = c("AMP", "DEL", "Upstream gene variant (only TERT)", "Stop gain", "Frameshift indel", "Inframe indel","Missense", "Splicing","Fusion"),
                            ncol = 1,
                            legend_gp = gpar(fontsize = 8))

gene_order <- c("PIK3CA","TERT","MAPK1","FGFR3","YAP1","KRAS","ERBB2","NFE2L2","PGR","KLF5","FBXW7","KMT2D","FAT1","KMT2C","PTEN","TP53","STK11","RB1","ARID1A","CASP8","EP300","ATRX","BAP1","RUNX1","SPEN","CREBBP","BRCA2","RASA1","ZFHX3","PTPRD","IFNGR1","ROBO1","TGFBR2","ANKRD11","LRP1B","AJUBA","APC","EPHA3","LATS1","GRIN2A","KEAP1","NOTCH1","NSD1")

op_oncogene <- oncoPrint(oncogenemat,
                         alter_fun = list(
                           
                           background = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                                       gp = gpar(fill = "#CCCCCC", col = NA)),
                           
                           AMP = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                                gp = gpar(fill = col["AMP"], col = NA)), 
                           
                           DEL = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                                gp = gpar(fill = col["DEL"], col = NA)),
                           
                           `TERT` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.8, 
                                                                   gp = gpar(fill = col["TERT"], col = NA)), 
                           
                           `stop gain` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.8, 
                                                                        gp = gpar(fill = col["stop gain"], col = NA)),
                           
                           `frameshift indel` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.65, 
                                                                               gp = gpar(fill = col["frameshift indel"], col = NA)),
                           
                           `inframe indel` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.6, 
                                                                            gp = gpar(fill = col["inframe indel"], col = NA)),
                           
                           missense = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.5, 
                                                                     gp = gpar(fill = col["missense"], col = NA)),
                           
                           splicing = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.3, 
                                                                     gp = gpar(fill = col["splicing"], col = NA)),
                           
                           stop_lost = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.2, 
                                                                      gp = gpar(fill = col["stop_lost"], col = NA)),
                           
                           start_lost = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.2, 
                                                                       gp = gpar(fill = col["start_lost"], col = NA)),
                           
                           fusion = function(x, y, w, h) {
                             grid.segments(x - w*0.3, y - h*0.3, x + w*0.3, y + h*0.3, gp = gpar(lwd = 1))
                             grid.segments(x + w*0.3, y - h*0.3, x - w*0.3, y + h*0.3, gp = gpar(lwd = 1))}
                           
                         ), col = col, 
                         row_order = oncogene_order,
                         heatmap_legend_param = heatmap_legend_param,
                         row_names_gp = gpar(fontsize = 12), pct_gp = gpar(fontsize = 8),
                         row_split = oncogene[match(rownames(oncogenemat), oncogene$Hugo.Symbol),"class"],
                         gap = unit(c(5),"mm"),
                         remove_empty_columns = FALSE, show_column_names = FALSE, alter_fun_is_vectorized = TRUE) 

ht_OP_oncogene = draw(op_oncogene)
patient_order <- column_order(ht_OP_oncogene)
new_patient_order <- colnames(oncogenemat)[patient_order][c(0:18, 20:50, 52:88, 90:91, 19,51,89, 92:291)] 

op2 <- oncoPrint(oncogenemat,
                 alter_fun = list(
                   
                   background = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                               gp = gpar(fill = "#CCCCCC", col = NA)),
                   
                   AMP = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                        gp = gpar(fill = col["AMP"], col = NA)), 
                   
                   DEL = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                        gp = gpar(fill = col["DEL"], col = NA)),
                   
                   `TERT` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.8, 
                                                           gp = gpar(fill = col["TERT"], col = NA)), 
                   
                   `stop gain` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.8, 
                                                                gp = gpar(fill = col["stop gain"], col = NA)),
                   
                   `frameshift indel` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.65, 
                                                                       gp = gpar(fill = col["frameshift indel"], col = NA)),
                   
                   `inframe indel` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.6, 
                                                                    gp = gpar(fill = col["inframe indel"], col = NA)),
                   
                   missense = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.5, 
                                                             gp = gpar(fill = col["missense"], col = NA)),
                   
                   splicing = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.3, 
                                                             gp = gpar(fill = col["splicing"], col = NA)),
                   
                   stop_lost = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.2, 
                                                              gp = gpar(fill = col["stop_lost"], col = NA)),
                   
                   start_lost = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.2, 
                                                               gp = gpar(fill = col["start_lost"], col = NA)),
                   
                   fusion = function(x, y, w, h) {
                     grid.segments(x - w*0.3, y - h*0.3, x + w*0.3, y + h*0.3, gp = gpar(lwd = 1))
                     grid.segments(x + w*0.3, y - h*0.3, x - w*0.3, y + h*0.3, gp = gpar(lwd = 1))}
                   
                   
                 ), col = col, 
                 top_annotation = column_ha,
                 row_order = oncogene_order, column_order = new_patient_order,
                 heatmap_legend_param = heatmap_legend_param,
                 row_names_gp = gpar(fontsize = 12), pct_gp = gpar(fontsize = 8),
                 row_split = oncogene[match(rownames(oncogenemat), oncogene$Hugo.Symbol),"class"],
                 gap = unit(c(5),"mm"),
                 remove_empty_columns = FALSE, show_column_names = FALSE, alter_fun_is_vectorized = TRUE) 

ht = draw(op2, annotation_legend_side = "bottom", heatmap_legend_side = "right", merge_legend = FALSE)

pdf("Figure_1._S.Barraud_et_al_Oncoprint_global.pdf",width = 20, height = 16)
ht
dev.off()
