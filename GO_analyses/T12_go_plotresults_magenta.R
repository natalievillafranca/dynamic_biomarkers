# this is the final script in GO analysis using MWU test, as in
# Voolstra et al PLoS ONE 2011, 6(5): e20392.
rm(list=ls()) 

####setting data up####
#Note: must run this script while results of GO_MWU are loaded into R
#must run these to initiate, then can ignore...
gg=read.table("Acer_iso2gene.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)  #the iso2gene file for your species
#dfull=read.csv("MAY_unmerged_hostGO_MM_magenta_ASV.csv") #read in proper VSD file for genes in module - generated in wgcna script 
dfull=read.csv("T12unmerged_hostGO_MM_magenta_ASV.csv") #read in proper VSD file for genes in module - generated in wgcna script 
rownames(dfull)<-dfull$X #make gene names your rownames
efulldata=c(2:186) #the columns where your expression data are (2:186 or 269) 
library(pheatmap)

###############plotting what significant GO genes are
#"BP_unmerged_hostGO_MM_purple_numeric.csv" is without ASV 
#MWU_BP_unmerged_hostGO_MM_purple_numeric.csv is without ASV
golist=read.table("BP_T12unmerged_hostGO_MM_magenta_ASV.csv",sep="	",header=T) #read in proper GO list from gomwu output - BP/MF/CC and module name
in.mwu=read.table("MWU_BP_T12unmerged_hostGO_MM_magenta_ASV.csv",sep=" ",header=T)
d=read.csv("T12host_vsd_MM_asvmagenta.csv")#read in proper VSD file for genes in module - generated in wgcna script 
rownames(d)<-d$X #make gene names your rownames
edata=c(2:187) #only columns with your expression data 

#gene=subset(golist,name=="response to organic cyclic compound") #write GO term of interest here from your sig list#
#t=rownames(gene)
#t

#is <- golist[rownames(golist) %in% t, c("seq", "term")]
#length(is) #then pulls out the isogroups with this annotation

#### OR pulling from alternate lists (skip this if using GO term above)
is<-c("isogroupDN72717_c4_g1", "isogroupDN74671_c2_g1", "isogroupDN57704_c0_g1", "isogroupDN74674_c1_g3") 
#in order: hub gene IRF, IRF T12vinter and prochlor, e. selectin, allograft inflammatory factor, tyrosine protein kinase receptor 


#####################first loop through genes matching GO term in module or in "good expression" subset
sel=c();gnms=c()
for ( i in is){
  if (i %in% d$X){
    sel=rbind(sel,d[d$X==i,])
    gnms=append(gnms,substr(paste(as.character(gg$V2[gg$V1==i]),collapse="."),1,50))
  }
}
row.names(sel)=paste(gnms,sel$X,sep=".")
nrow(sel)
rownames(sel)

exp=sel[,edata]
if (length(exp[,1])==1) { exp=rbind(exp,exp) }
nrow(exp)
#should match 'good genes' tabulation in figure
rownames(exp)

##May 
#sel <- data.frame()
#gnms <- c()

#for (i in is$seq) {
#  if (i %in% d$X) {
#    sel <- rbind(sel, d[d$X == i, ])
#    gnms <- append(gnms, substr(paste(as.character(gg$V2[gg$V1 == i]), collapse = "."), 1, 50))
#  }
#}

#row.names(sel)=paste(gnms,sel$X,sep=".")
#nrow(sel)
#rownames(sel)

#exp=sel[,edata]
#if (length(exp[,1])==1) { exp=rbind(exp,exp) }
#nrow(exp)
#should match 'good genes' tabulation in figure
#rownames(exp)
#------------------ heatmap of GOgenes stacked by symbiont density (or other trait...)
library(pheatmap)

#reading in traits 
#lnames = load(file = "MAY25_BiomarkersTraitsAndInteractionsAssigned.RData") 
lnames=load(file="MAY25_BiomarkersTraitsAndInteractionsAssigned_T12.RData")

#fixing name error, refer to merged_data in WGCNA code for correct number (2)
#datTraits$replicate[datTraits$SRR == "SRR10600259"] <- 2

#Changing 34 to 36 in genotype column 
datTraits <- datTraits %>%
  mutate(genotype = ifelse(genotype == 34, 36, genotype))

datTraitsT12 <- datTraitsT12 %>%
  mutate(genotype = ifelse(genotype == 34, 36, genotype))

#datTraits <- datTraits %>%
#  mutate(sample = paste0("G", genotype, "r", replicate)) %>%
#  select(sample, everything())  # Move 'sample' to the first column

#lnames=load(file="datTraitsT12.RData")
#datTraitsT12 <- datTraitsT12 %>%
#  mutate(sample = paste0("G", genotype, "r", replicate)) %>%
#  select(sample, everything())  # Move 'sample' to the first column


library(dplyr)

names(datTraitsT12)
# rename the traits to look nice for the heatmap, and reorder them so they are grouped by symb etc
datTraitsT12=datTraitsT12%>%relocate("50", .after = "T12_Vinter")
datTraitsT12=datTraitsT12%>%relocate("3", .after ="50")
datTraitsT12=datTraitsT12%>%relocate("1", .after ="3")
datTraitsT12=datTraitsT12%>%relocate("31", .after ="1")
datTraitsT12=datTraitsT12%>%relocate("44", .after ="31")
datTraitsT12=datTraitsT12%>%relocate("36", .after ="44")
datTraitsT12=datTraitsT12%>%relocate("7", .after ="36")
datTraitsT12=datTraitsT12%>%relocate("41", .after ="7")
datTraitsT12=datTraitsT12%>%relocate("62", .after ="41")
datTraitsT12=datTraitsT12%>%relocate("13", .after ="62")

#genotype only 
datTraits1=datTraitsT12[,c(1)]

totresval_genotype_order <- c(50, 3, 1, 31, 44, 36, 7, 41, 62, 13)
#datTraits$genotype <- factor(datTraits$genotype, levels = totresval_genotype_order)
#datTraits <- datTraits[order(datTraits$genotype),]
datTraitsT12$genotype <- factor(datTraitsT12$genotype, levels = totresval_genotype_order)
datTraitsT12 <- datTraitsT12[order(datTraitsT12$genotype),]

####
sorted=exp[,rownames(datTraitsT12)] #sort rows by values of a trait *** if you don't want to sort, just say sorted=exp

#sorted<-exp #do this if you do not want to sort by predefined trait

############# Center scale
means=apply(sorted,1,mean) # means of rows
#hist(means)
expcT0=sorted-means #rescale expression data
#save(expcT12, file= "T12expc.RData")
#save(expcT0, file="T0expc.RData")
#lnames=load(file="T12expc.RData")
#lnames=load(file="T0expc.RData")

#colnames(expcT0) <- datTraits$sample[match(colnames(expcT0), rownames(datTraits))]
#colnames(expcT12) <- datTraitsT12$sample[match(colnames(expcT12), rownames(datTraitsT12))]

#all_columns <- union(colnames(expcT0), colnames(expcT12))

# Add missing columns to expcT0
#for (col in setdiff(all_columns, colnames(expcT0))) {
#  expcT0[[col]] <- NA
#}

# Add missing columns to expcT12
#for (col in setdiff(all_columns, colnames(expcT12))) {
#  expcT12[[col]] <- NA
#}

# Reorder columns to ensure both dataframes have the same column order
#expcT0 <- expcT0[, all_columns]
#expcT12 <- expcT12[, all_columns]

#column_numbers <- as.numeric(gsub("G([0-9]+).*", "\\1", colnames(expcT0)))

# Reorder columns based on the desired genotype order
#expcT0 <- expcT0[, order(match(column_numbers, totresval_genotype_order))]
#expcT12 <- expcT12[, order(match(column_numbers, totresval_genotype_order))]

# Verify the column order
#colnames(expcT0)
#colnames(expcT12)

#names <- data.frame(matrix(ncol = 0, nrow = length(colnames(expcT0)))) 
# Assign the column names of expcT0 as row names of the names dataframe
#rownames(names) <- colnames(expcT0)
# Add the genotype column, extracting the number after 'G' from the row names
#names$genotype <- as.numeric(gsub("G([0-9]+).*", "\\1", rownames(names)))

# View the updated dataframe
#head(names)


library(RColorBrewer)
col=color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")),bias=1.2)(30) #adjust bias to center 0

#setting up annotation_col option to get legends saying what genotype and status each sample is 
#only if you are trying to see gene correlation with genotype as a trait 

Genotype <- datTraitsT12$genotype
annotation_col <- data.frame(Genotype = Genotype)
rownames(annotation_col) <- colnames(expcT0)  # Ensure rownames match the column names of expc
annotation_colors=list(Genotype=c("50"="#053061", "3"="#2166ac", "1"="#4393c3", "31"="#92c5de", 
                                    "44"="#d1e5f0", "36"="#fddbc7", "7"="#f4a582", "41"="#d6604d",
                                     "62"="#b2182b", "13"="#67001f"))
genotype_breaks <- cumsum(table(Genotype))

pheatmap(
  expcT0, 
  color = col,
  cluster_cols = FALSE, 
  clustering_distance_rows = "binary",
  show_colnames = FALSE,
  border_color = "black",  # Adds black borders to all cells
  gaps_col = genotype_breaks # Adds lines between genotype groups
)


gene=as.data.frame(t(sorted))
texp=data.frame(t(sorted))
sexp=stack(texp)
unique(sexp$ind) #unique returns a list of all the gene names in your GOterm for the module of choice

final <- bind_cols(gene,
                   genotype = datTraitsT12$genotype, 
                   Prochlorococcus=datTraitsT12$`Prochlorococcus MIT9313`,
                   Synechococcus=datTraitsT12$`Synechococcus CC9902`,
                   MD355=datTraitsT12$`MD3-55`, 
                   NS5 = datTraitsT12$`NS5 marine group`,
                   HIMB11 = datTraitsT12$HIMB11, 
                   CandidatusActinomarina = datTraitsT12$`Candidatus Actinomarina`, 
                   Shannon = datTraitsT12$Shannon)

####GE correlation with genotype ####
#running stats to see how correlated GE is with genotype and how different exp is between genotypes 
#with anova and post hoc tukey, not goint to check for normalization because we normalized counts in WGCNA with DESEQ2
names(final)
#IRF2 HUB
IRF2hubAOV <- aov(`Interferon regulatory factor 2 OS=Gallus gallus OX.isogroupDN72717_c4_g1` ~ genotype, data = final)
summary(IRF2hubAOV)
AIC(IRF2hubAOV)  
IRF2hubAOV<- TukeyHSD(IRF2hubAOV)
tukey_table <- IRF2hubAOV$genotype

#e selectin
Esel.aov <- aov(final$`E-selectin OS=Bos taurus OX=9913 GN=SELE PE=2 SV=1.isogroupDN74671_c2_g1` ~ final$genotype)
summary(Esel.aov)
AIC(Esel.aov)
Esel.Tuk <- TukeyHSD(Esel.aov)
tukey_results <- Esel.Tuk$`final$genotype`
pvals <- tukey_results[, "p adj"]
names(pvals) <- rownames(tukey_results)

cld_vec <- multcompLetters(pvals)$Letters
letter_df <- data.frame(
  genotype = names(cld_vec),
  cld = as.character(cld_vec),
  stringsAsFactors = FALSE
)

letter_df$genotype <- factor(letter_df$genotype, levels = totresval_genotype_order)
letter_df <- letter_df[order(letter_df$genotype), ]
print(letter_df)
#genotype cld
#10       50  ab
#1         3   a
#2         1  ab
#3        31 abc
#4        44   a
#5        36 bcd
#6         7   a
#7        41   e
#8        62  cd
#9        13   d

#allograft inflammatory factor 1 
Allograft.aov <-aov(final$`Allograft inflammatory factor 1 OS=Homo sapiens OX.isogroupDN57704_c0_g1` ~ final$genotype)
summary(Allograft.aov)
AIC(Allograft.aov) ## not great... 54.09729. 
Allograft.Tuk <- TukeyHSD(Allograft.aov)
tukey_results <- Allograft.Tuk$`final$genotype`
pvals <- tukey_results[, "p adj"]
names(pvals) <- rownames(tukey_results)

cld_vec <- multcompLetters(pvals)$Letters
letter_df <- data.frame(
  genotype = names(cld_vec),
  cld = as.character(cld_vec),
  stringsAsFactors = FALSE
)

letter_df$genotype <- factor(letter_df$genotype, levels = totresval_genotype_order)
letter_df <- letter_df[order(letter_df$genotype), ]
print(letter_df)
#genotype cld
#      50  ab
#      3   a
#      1  ab
#      31   c
#      44   a
#      36  cd
#      7  bc
#      41   d
#      62  cd
#      13   d

#tyrosine protein kinase receptor tie 1 
PTKaov <- aov(final$`Tyrosine-protein kinase receptor Tie-1 (Fragment) .isogroupDN74674_c1_g3` ~ final$genotype)
summary(PTKaov)
AIC(PTKaov)
PTK.Tuk <- TukeyHSD(PTKaov)
tukey_results <- PTK.Tuk$`final$genotype`
pvals <- tukey_results[, "p adj"]
names(pvals) <- rownames(tukey_results)

cld_vec <- multcompLetters(pvals)$Letters
letter_df <- data.frame(
  genotype = names(cld_vec),
  cld = as.character(cld_vec),
  stringsAsFactors = FALSE
)

letter_df$genotype <- factor(letter_df$genotype, levels = totresval_genotype_order)
letter_df <- letter_df[order(letter_df$genotype), ]
print(letter_df)

#genotype cld
#     50   a
#     3  ab
#     1   a
#     31   c
#     44   a
#     36   a
#      7  bc
#     41   c
#     62   c
#     13   c
####

####microbe ~ GE + genotype ####
library(lme4)
library(lmerTest)
#only want to use samples where there is both ASV and GE data 
final <- final[!is.na(final$Synechococcus), ]
#IRF2 and syn 

# IRF (nonhub) and syn
result <- lm(`Interferon regulatory factor 2 OS=Gallus gallus OX.isogroupDN72717_c4_g1` ~ Synechococcus + genotype, data = final) #only significant with G7
summary(result)
anova(result)
car::Anova(result)

##### IRF nonhub and pro 
#result<- lm(`Interferon regulatory factor 2 OS=Gallus gallus OX.isogroupDN72717_c4_g1` ~ Prochlorococcus + genotype, data = final) 
#summary(result)

###IRF hub and pro 
#result <- lm(`Interferon regulatory factor 2 OS=Sigmodon hispidu.isogroupDN102558_c0_g1` ~ Prochlorococcus +  genotype, data = final) 
#summary(result) 

##### IRF hub and syn 
result <- lm(`Interferon regulatory factor 2 OS=Sigmodon hispidu.isogroupDN102558_c0_g1` ~ Synechococcus +  genotype, data = final) 
summary(result) 

result<- lm(`Interferon regulatory factor 2 OS=Gallus gallus OX.isogroupDN72717_c4_g1` ~ Synechococcus + genotype, data = final) 
summary(result) 



result <- lm(`Interferon regulatory factor 2 OS=Sigmodon hispidu.isogroupDN102558_c0_g1` ~ MD355 +  genotype, data = final) 
summary(result) 

result <- lm(`Interferon regulatory factor 2 OS=Gallus gallus OX.isogroupDN72717_c4_g1` ~ MD355 +  genotype, data = final) 
summary(result) 

result <- lm(`E-selectin OS=Bos taurus OX=9913 GN=SELE PE=2 SV=1.isogroupDN74671_c2_g1` ~ MD355 +  genotype, data = final) 
summary(result) 

####



