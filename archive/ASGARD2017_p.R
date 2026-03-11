### ASGARD processing site
library(tidyverse)
library(janitor)
library(here)
library(gplots) # for heatmap.2()
library(viridis)
library(scales) # for rescale()
library(vegan)

### 全サンプル含むdf作る
seqtab_16Spropcol_mradec = sort(seqtab_16Spropcol_mra, decreasing = TRUE) # sort it with descending order (length is 3076)

seqtab_16Spropcol_3076 <- seqtab_16Sprop[,names(seqtab_16Spropcol_mradec)] # 1633*3076
# top100を抽出するmradec2ではなく、3076サンプルを指定するmradecを使う

ODV.otu.pick_full = match(colnames(seqtab_16Spropcol_3076), colnames(seqtab_filt)) 
colnames(seqtab_16Spropcol_3076) = shorternames[ODV.otu.pick_full] #オブジェクト作るみたいな構造に見えるが、colnames() = で列名を変更できる

seqtab_ODV_full = merge(meta_denovo_2, seqtab_16Spropcol_3076, by = 0, all.x = FALSE)
# 960*3122
# 3122 = 3076 + 46 (metadata_denovo_2は1163*45, Row.namesが含まれて+1)
#rownames(seqtab_ODV_full) = seqtab_ODV_full$Row.names

seqtab_ODV_Processing_full = seqtab_ODV_full[seqtab_ODV_full$station_type=="P",]
# 241*3122
# S = 689, P = 241, E = 30

ASGARD_p <- seqtab_ODV_Processing_full %>%
  filter(project == "ASGARD2017") #81*3122

# Extract data from Processing stations 
asgard_vec_p <- ASGARD_p$Row.names #length is 81
seqtab_prop_asgard_p <- seqtab_16Sprop[asgard_vec_p, ] #81*3076

seqtab_proprow_asgard_p <- rowSums(seqtab_prop_asgard_p) 
hist(seqtab_proprow_asgard_p, breaks = 100) 

# ASGARDのmeta data抽出
meta_asgard_p <- meta_denovo_2[asgard_vec_p, ] #81*45

# prop tableの種類！
# 普通のprop tableはseqtab_16Spropcol_mra100mat
# そこからasgardだけ抽出したのは、seqtab_prop_asgard (181*3076)

# Set a minimum read count threshold (e.g., ASVs with <5 total reads are removed)
mincutoff_p = (apply(seqtab_prop_asgard_p, 2, max) > 0.001)
asgard_filtered_p <- seqtab_prop_asgard_p[, mincutoff_p]
asgard_filtered_p <- asgard_filtered_p[, colSums(asgard_filtered_p > 0) > 2]
hist(rowSums(asgard_filtered_p),breaks=20)
dim(asgard_filtered_p) #81*223

ODV.otu.pick_p = match(colnames(asgard_filtered_p), colnames(seqtab_filt)) #length is 258
colnames(asgard_filtered_p) = shorternames[ODV.otu.pick_p] 

# Merge ASV table and meta data
asgard_filtered_p <- as.data.frame(asgard_filtered_p)
asgard_filtered_p$Sample <- rownames(asgard_filtered_p)

meta_asgard_p$Sample <- rownames(meta_asgard_p)

# Add " µm" behind 0.2, 3 and 20
meta_asgard_p <- meta_asgard_p %>%
  mutate(filter = paste0(filter, " µm")) # Add " µm" behind 0.2, 3 and 20

asgard_processing <- left_join(meta_asgard_p, asgard_filtered_p, by = "Sample")
#81*269


### Ternary Plot

## df with only "filter" column and ASV columns
asgard_processing_frac <- asgard_processing %>% 
  select(filter, contains("ESV")) #81*224

 df_for_ternary_ave <- asgard_processing_frac %>%
   group_by(filter) %>%
   summarise(across(everything(), mean, na.rm = TRUE)) %>%
   select(-filter)
# everything(): select all columns
# across(): for all selected cols, process something all at once
# in this case, calculate mean = mean, na.rm = TRUE: 平均値を計算（欠損値は除外）)

## make a prop table
rownames(df_for_ternary_ave) <- c("0.2 µm", "20 µm", "3 µm") 
df_for_ternary_ave <- t(df_for_ternary_ave) #transposition
ternary_prop <- prop.table(df_for_ternary_ave, margin = 1) #margin = 1 calculates proportion of each row!
 
### Using ggtern
install.packages("ggtern")
library(ggtern)

# データ例（dfは223×3で、行名はASV）
# 列名は "0.2 um", "3 um", "20 um"

## ggternで三角ダイアグラムを作る

ggtern(data = ternary_prop, 
       aes(x = `0.2 µm`, y = `3 µm`, z = `20 µm`)) +
  geom_point(alpha = 0.6, size = 2, color = "blue") +  # 点の大きさ・色・透明度
  labs(L = "0.2 µm",  # 左 = x
      R = "3 µm",    # 右 = y
      T = "20 µm") + # 上 = z
  theme_bw() +  # 白背景
  theme_showarrows()

## different colors 
ternary_prop_color=as.data.frame(ternary_prop)
ternary_prop_color$dominant <- apply(ternary_prop_color, 1, function(x) {
  colnames(ternary_prop_color)[which.max(x)]
  })

## make dominant column factor
ternary_prop_color$dominant <- factor(ternary_prop_color$dominant, levels = c("0.2 µm", "3 µm", "20 µm"))
ggtern(data = ternary_prop_color,
       aes(x = `0.2 µm`, y = `3 µm`, z = `20 µm`, color = dominant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c(`0.2 µm`="red", `3 µm`="green", `20 µm`="blue")) +
  labs(x = "0.2 µm", y = "3 µm", z = "20 µm") +
  theme_bw() +
  theme_showarrows()+
  theme_rotate(-90)

# in aes(), x = left, y = top, z = right


# ### Heatmap
# 
# ### ASGARDのhm
# library(viridis)
# 
# # delete Sample column (that is non-numeric)
# asgard_filtered_p_hm <- asgard_filtered_p %>%
#   select(-Sample)
# 
# # We need matrix for hm
# asgard_filtered_p_hm <- as.matrix(asgard_filtered_p_hm)
#   
# pdf(file="ASGARD_hm_processing.pdf") #ASV全表示したいなら, width = 20, height = 10
# 
# h1=heatmap((asgard_filtered_p_hm)^.25, 
#            distfun=function(x) vegdist(x, method="bray"),
#            hclustfun=function(x) hclust(x, method="ward.D"),
#            keep.dendro=TRUE, 
#            scale="none",
#            margins = c(20,20),
#            main="ASGARD_bray/ward.D2")
# 
# nclus_p = 4 #枝が5くらいに分けられそうだった
# oldclus_p = cutree(as.hclust(h1$Rowv),k=nclus_p)
# oldorder_p = unname(rle(oldclus_p[as.hclust(h1$Rowv)$order])$values)
# neworder_p = (1:nclus_p)
# names(neworder_p) = oldorder_p
# clusnum_p = unname(neworder_p[as.character(oldclus_p)])
# names(clusnum_p) = names(oldclus_p)
# 
# colnclus_p = 6 
# cololdclus_p = cutree(as.hclust(h1$Colv),k=colnclus_p)
# cololdorder_p = unname(rle(cololdclus_p[as.hclust(h1$Colv)$order])$values)
# colneworder_p = (1:colnclus_p)
# names(colneworder_p) = cololdorder_p
# colclusnum_p = unname(colneworder_p[as.character(cololdclus_p)])
# names(colclusnum_p) = names(cololdclus_p)
# 
# #rsc=viridis(nclus)[clusnum] #rscではなくviridisにしてもいい
# #colrsc=viridis(colnclus)[colclusnum] #rscではなくviridisにしてもいい
# 
# rsc_p=hue_pal()(nclus_p)[clusnum_p]
# colrsc_p=hue_pal()(colnclus_p)[colclusnum_p]
# colrsc_p=plasma(colnclus_p)[colclusnum_p] #列クラスターを行クラスターの色を変えたいなら
# 
# 
# h2=heatmap.2((asgard_filtered_p_hm)^.25, 
#              distfun=function(x) vegdist(x, method="bray"),
#              hclustfun=function(x) hclust(x, method="ward.D"),
#              col = viridis, 
#              RowSideColors = rsc_p,
#              ColSideColors = colrsc_p,
#              margins = c(15,15),
#              scale="none",
#              main="ASGARD_bray/ward.D2",
#              trace="none",
#              cexCol = 0.2, #列の文字サイズ
#              cexRow = 0.5 #行の文字サイズ
#              #labRow = NA,        
#              #labCol = NA
# )
# 
# ## change the column color, corresponding to rgb in Ternary Plot
# asv_rgb <- rgb(ternary_prop_color$`0.2 µm`,
#                   ternary_prop_color$`3 µm`,
#                   ternary_prop_color$`20 µm`)
# # length is 223. 189 color codes. 
# 
# ## change the row color
# filter_colors <- c("0.2 µm"="red",
#                    "3 µm"="green",
#                    "20 µm"="blue")
# 
# sample_rgb <- filter_colors[as.character(meta_asgard_p$filter)]
# # length is 81. 3 color codes (red, green and blue only).
# 
# h3=heatmap.2((asgard_filtered_p_hm)^.25, 
#              distfun=function(x) vegdist(x, method="bray"),
#              hclustfun=function(x) hclust(x, method="ward.D"),
#              col = viridis, 
#              RowSideColors = sample_rgb,
#              ColSideColors = asv_rgb,
#              margins = c(15,15),
#              scale="none",
#              main="ASGARD_bray/ward.D2",
#              trace="none",
#              cexCol = 0.2, #列の文字サイズ
#              cexRow = 0.2, #行の文字サイズ
#              labRow = meta_asgard_p$filter      
#              #labCol = NA
# )
# 
# dev.off()

## ASV extraction from Ternary_prop_color
# greater than or equal to 0.7 in 0.2 um
ext_0.2 <- ternary_prop_color %>%
  filter(`0.2 µm` >= 0.7)
asv_names_0.2 <- rownames(ext_0.2) # 72 ASVs 

# greater than or equal to 0.7 in 3 um
ext_3 <- ternary_prop_color %>%
  filter(`3 µm` >= 0.7)
asv_names_3 <- rownames(ext_3) # 9 ASVs

# greater than or equal to 0.7 in 20 um
ext_20 <- ternary_prop_color %>%
  filter(`20 µm` >= 0.7)
asv_names_20 <- rownames(ext_20) # 36 ASVs


### sequence count
rownames(ASGARD_p) <- ASGARD_p$Row.names

common_asv <- intersect(rownames(ASGARD_p),
                        rownames(seqtab_16Smat)) # length is 81

asgard_seqcount <- seqtab_16Smat[common_asv, ] # dim is 81*3076
# dim(seqtab_16Smat) is 1633*3076
#, drop=FALSE

seqtab_asgard_16Srow <- rowSums(asgard_seqcount) # length is 81

hist(seqtab_asgard_16Srow, breaks = 100) #xlim = c(0, 10000)

### 5000で切ってみる
seqtab_16S_asgard_10000 <- seqtab_matrow >=5000 #length is 2193
seqtab_16Smat_asgard_10000 <- seqtab_mat[seqtab_16S_asgard_10000,] #1238*3076
seqtab_16Srow_asgard_10000 <- rowSums(seqtab_16Smat_asgard_10000) #length is 1238
common_asv2 <- intersect(rownames(ASGARD_p),
                        rownames(seqtab_16Smat_asgard_10000)) #length is 78

asgard_seqcount2 <- seqtab_16Smat[common_asv2, ] # dim is 78*3076
seqtab_asgard_16Srow_2 <- rowSums(asgard_seqcount2) # length is 78
hist(seqtab_asgard_16Srow_2, breaks = 100)

## hmをrerunする(missingが消えるかも！？)
asgard_seqcount_16Sprop <- prop.table(asgard_seqcount2, 1) # 78*3076
asgard_vec_p2 <- rownames(asgard_seqcount2) #length is 78
seqtab_prop_asgard_p2 <- asgard_seqcount_16Sprop[asgard_vec_p2, ] #78*3076

# Set a minimum read count threshold (e.g., ASVs with <5 total reads are removed)
mincutoff_p2 = (apply(seqtab_prop_asgard_p2, 2, max) > 0.001)
asgard_filtered_p2 <- seqtab_prop_asgard_p2[, mincutoff_p2]
asgard_filtered_p2 <- asgard_filtered_p2[, colSums(asgard_filtered_p2 > 0) > 2]
hist(rowSums(asgard_filtered_p2),breaks=20)
dim(asgard_filtered_p2) #78*221

ODV.otu.pick_p2 = match(colnames(asgard_filtered_p2), colnames(seqtab_filt)) 
colnames(asgard_filtered_p2) = shorternames[ODV.otu.pick_p2] 

# Merge ASV table and meta data
asgard_filtered_p2 <- as.data.frame(asgard_filtered_p2)
asgard_filtered_p2$Sample <- rownames(asgard_filtered_p2) #78*222
vec78 <- rownames(asgard_filtered_p2)

meta_asgard_p2 <- meta_asgard_p[vec78, ] #78*46

asgard_processing2 <- left_join(meta_asgard_p2, asgard_filtered_p2, by = "Sample")
#78*267

asgard_processing_frac2 <- asgard_processing2 %>% 
  select(filter, contains("ESV")) #78*222

df_for_ternary_ave2 <- asgard_processing_frac2 %>%
  group_by(filter) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  select(-filter)

rownames(df_for_ternary_ave2) <- c("0.2 µm", "20 µm", "3 µm") 
df_for_ternary_ave2 <- t(df_for_ternary_ave2) #transposition
ternary_prop2 <- prop.table(df_for_ternary_ave2, margin = 1) #margin = 1 calculates proportion of each row!

# different colors 
ternary_prop_color2=as.data.frame(ternary_prop2)
ternary_prop_color2$dominant <- apply(ternary_prop_color2, 1, function(x) {
  colnames(ternary_prop_color2)[which.max(x)]
})

# make dominant column factor
ternary_prop_color2$dominant <- factor(ternary_prop_color2$dominant, levels = c("0.2 µm", "3 µm", "20 µm"))

ggtern(data = ternary_prop_color2,
       aes(x = `0.2 µm`, y = `3 µm`, z = `20 µm`, color = dominant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c(`0.2 µm`="red", `3 µm`="green", `20 µm`="blue")) +
  labs(x = "0.2 µm", y = "3 µm", z = "20 µm") +
  theme_bw() +
  theme_showarrows()+
  theme_rotate(-90)

### Make hm >5000 reads

# change the column color, corresponding to rgb in Ternary Plot
asv_rgb2 <- rgb(ternary_prop_color2$`0.2 µm`,
               ternary_prop_color2$`3 µm`,
               ternary_prop_color2$`20 µm`)
# length is 221. 189 color codes. 

## change the row color
filter_colors <- c("0.2 µm"="red",
                   "3 µm"="green",
                   "20 µm"="blue")

sample_rgb2 <- filter_colors[as.character(meta_asgard_p2$filter)]
# length is 78 3 color codes (red, green and blue only).

#必要なら、asgard_filtered_p_hm2 <- asgard_filtered_p2 %>% 
  #select(-Sample) #78*221

asgard_filtered_p_hm2 <- as.matrix(asgard_filtered_p_hm2)

asgard_filtered_49asvs <- asgard_filtered_p_hm2[, zero_cols] #78*49
asgard_filtered_49asvs2 <- asgard_filtered_49asvs[rowSums(asgard_filtered_49asvs) > 0, colSums(asgard_filtered_49asvs) > 0]

meta_asgard_p2_49 <- meta_asgard_p2[rownames(asgard_filtered_49asvs2),]

#asgard_fl <- asgard_filtered_p_hm2[, freeliving_no_na]
#asgard_suspended <- asgard_filtered_p_hm2[, suspended_no_na]
#asgard_sinking <- asgard_filtered_p_hm2[, particleass_no_na2]

#meta_asgard_fl <- meta_asgard_p[rownames(asgard_fl),]
#meta_asgard_suspended <- meta_asgard_p[rownames(asgard_suspended),]
#meta_asgard_sinking <- meta_asgard_p[rownames(asgard_sinking),]


pdf(file="ASGARD_hm_processing_5000over.pdf", width = 20, height = 20) #ASV全表示したいなら, width = 20, height = 10

h3=heatmap.2((asgard_filtered_p_hm2)^.25, 
             distfun=function(x) vegdist(x, method="bray"),
             hclustfun=function(x) hclust(x, method="ward.D"),
             col = viridis, 
             RowSideColors = sample_rgb2,
             ColSideColors = asv_rgb2,
             margins = c(15,15),
             scale="none",
             main="ASGARD_bray/ward.D2",
             trace="none",
             cexCol = 0.2, #列の文字サイズ
             cexRow = 0.2, #行の文字サイズ
            
             labRow = meta_asgard_p2$side     
             #labCol = NA
)

h4=heatmap.2((asgard_filtered_49asvs2)^.25, 
             distfun=function(x) vegdist(x, method="bray"),
             hclustfun=function(x) hclust(x, method="ward.D"),
             col = viridis, 
             #RowSideColors = sample_rgb2,
             #ColSideColors = asv_rgb2,
             margins = c(15,15),
             scale="none",
             main="ASGARD_bray/ward.D2",
             trace="none",
             cexCol = 0.8, #列の文字サイズ
             cexRow = 0.8, #行の文字サイズ
             labRow = meta_asgard_p2_49$side     
             #labCol = NA
)

#h5=heatmap.2((asgard_fl)^.25, 
             #distfun=function(x) vegdist(x, method="bray"),
             #hclustfun=function(x) hclust(x, method="ward.D"),
             #col = viridis, 
             #RowSideColors = sample_rgb2,
             #ColSideColors = asv_rgb2,
             #margins = c(15,15),
             #scale="none",
             #main="ASGARD_bray/ward.D2_0.2µm",
             #trace="none",
             #cexCol = 0.8, #列の文字サイズ
             #cexRow = 0.8, #行の文字サイズ
             #labRow = meta_asgard_fl$side     
             #labCol = NA
#)

#h6=heatmap.2((asgard_suspended)^.25, 
             #distfun=function(x) vegdist(x, method="bray"),
             #hclustfun=function(x) hclust(x, method="ward.D"),
             #col = viridis, 
             #RowSideColors = sample_rgb2,
             #ColSideColors = asv_rgb2,
             #margins = c(15,15),
             #scale="none",
             #main="ASGARD_bray/ward.D2_3µm",
             #trace="none",
             #cexCol = 0.8, #列の文字サイズ
             #cexRow = 0.8, #行の文字サイズ
             #labRow = meta_asgard_suspended$side     
             #labCol = NA
#)

#h7=heatmap.2((asgard_sinking)^.25, 
             #distfun=function(x) vegdist(x, method="bray"),
             #hclustfun=function(x) hclust(x, method="ward.D"),
             #col = viridis, 
             #RowSideColors = sample_rgb2,
             #ColSideColors = asv_rgb2,
             #margins = c(15,15),
             #scale="none",
             #main="ASGARD_bray/ward.D2_20µm",
             #trace="none",
             #cexCol = 0.8, #列の文字サイズ
             #cexRow = 0.8, #行の文字サイズ
             #labRow = meta_asgard_sinking$side     
             #labCol = NA
#)

# 全体の長さが181になっているので、ここ要修正
nclus_p = 4 
oldclus_p = cutree(as.hclust(h3$rowDendrogram),k=nclus_p)
oldorder_p = unname(rle(oldclus_p[as.hclust(h3$rowDendrogram)$order])$values)
neworder_p = (1:nclus_p)
names(neworder_p) = oldorder_p
clusnum_p = unname(neworder_p[as.character(oldclus_p)])
names(clusnum_p) = names(oldclus_p)

rsc_p <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
sample_rgb3 <- rsc_p[clusnum_p]

h8=heatmap.2((asgard_filtered_p_hm2)^.25, 
             distfun=function(x) vegdist(x, method="bray"),
             hclustfun=function(x) hclust(x, method="ward.D"),
             col = viridis, 
             RowSideColors = sample_rgb3,
             Rowv = h3$rowDendrogram,
             ColSideColors = asv_rgb2,
             margins = c(15,15),
             scale="none",
             main="ASGARD_bray/ward.D2",
             trace="none",
             cexCol = 0.2, #列の文字サイズ
             cexRow = 0.2, #行の文字サイズ
             labRow = meta_asgard_p2$side     
             #labCol = NA
)


h9=heatmap.2((asgard_filtered_p_hm2)^.25, 
             distfun=function(x) vegdist(x, method="bray"),
             hclustfun=function(x) hclust(x, method="ward.D"),
             col = viridis, 
             RowSideColors = sample_rgb3,
             Rowv = h3$rowDendrogram,
             ColSideColors = asv_rgb2,
             margins = c(15,15),
             scale="none",
             main="ASGARD_bray/ward.D2",
             trace="none",
             cexCol = 0.2, #列の文字サイズ
             cexRow = 0.2, #行の文字サイズ
             labRow = meta_asgard_p2$Sample    
             #labCol = NA
)

dev.off()

### どのサンプルがどこで撮られたかチェックする地図
# サンプル番号を付与
a = meta_asgard_p2
a$ID <- 1:nrow(a)

a$cluster = as.factor(clusnum_p)

a$depth_type <- factor(a$depth_type, levels = c("surf", "mid", "bottom"))

a <- a %>%
  filter(Sample != "BOX_6_26") #77*48

# 地図＋点＋サンプル番号ラベル
pdf("processing_map.pdf", width = 20, height = 20)

map_plot_3 <- ggmap(mapz) +
  geom_point(data = a, 
             aes(x = lon, y = lat, color = cluster),
             alpha = 1,
             size = 3) +
  geom_text_repel(data = a,
                  aes(x = lon, y = lat, label = ID),
                  size = 3,
                  max.overlaps = Inf) +  # 重なりを自動調整
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(depth_type ~ cluster)+
  labs(color = "cluster")

print(map_plot_3)

map_plot_4 <- ggmap(mapz) +
  geom_point(data = a, 
             aes(x = lon, y = lat, color = cluster),
             alpha = 1,
             size = 3) +
  geom_text_repel(data = a,
                  aes(x = lon, y = lat, label = ID),
                  size = 3,
                  max.overlaps = Inf) +  # 重なりを自動調整
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(. ~ cluster)+
  labs(color = "cluster")

print(map_plot_4)

dev.off()



# 番号とサンプル名の対応表を作成
legend_table <- a[, c("ID", "side")]
print(legend_table)



### FL, PAのabundance map

asgard_fl <- asgard_filtered_p_hm2[, freeliving_no_na]
asgard_suspended <- asgard_filtered_p_hm2[, suspended_no_na]
asgard_sinking <- asgard_filtered_p_hm2[, particleass_no_na2]

meta_asgard_fl <- meta_asgard_p[rownames(asgard_fl),]
meta_asgard_suspended <- meta_asgard_p[rownames(asgard_suspended),]
meta_asgard_sinking <- meta_asgard_p[rownames(asgard_sinking),] #78*46

asgard_sinking_sum <- as.data.frame(asgard_sinking)
asgard_sinking_sum$sum <- rowSums(asgard_sinking_sum) #78*59
asgard_sinking_sum$Sample <- rownames(asgard_sinking_sum) #78*60

asgard_sinking_map <- left_join(meta_asgard_sinking, asgard_sinking_sum, by = "Sample")

min(asgard_sinking_map$sum) #0.0130147
max(asgard_sinking_map$sum) #0.7256913

map_plot_2 <- ggmap(mapz) +
  geom_point(data = asgard_sinking_map, aes(x = lon, y = lat, size = sum)
             ,alpha = 1) +
  #scale_size(range = c(3, 15))+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  facet_grid(. ~ filter)

print(map_plot_2)

# 念の為、Surveyデータも見る
particleass_no_na3 <- particleass_no_na[
  particleass_no_na %in% colnames(asgard_filtered)
] #59から42になる

asgard_sinking_Survey <- asgard_filtered[, particleass_no_na3] #181*42
meta_asgard_sinking_Survey <- meta_asgard[rownames(asgard_sinking_Survey),] #181*45

asgard_sinking_sum_Survey <- as.data.frame(asgard_sinking_Survey)
asgard_sinking_sum_Survey$sum <- rowSums(asgard_sinking_sum_Survey) #181*44
asgard_sinking_sum_Survey$Sample <- rownames(asgard_sinking_sum_Survey) #181*44
meta_asgard_sinking_Survey$Sample <- rownames(asgard_sinking_sum_Survey)

asgard_sinking_map_Survey <- left_join(meta_asgard_sinking_Survey, asgard_sinking_sum_Survey, by = "Sample")

map_plot_3 <- ggmap(mapz) +
  geom_point(data = asgard_sinking_map_Survey, aes(x = lon, y = lat, size = sum^2)
             ,alpha = 1) +
  #scale_size(range = c(3, 15))+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  facet_grid(. ~ depth_type)

print(map_plot_3)

### sinkingのbox plot

#########asgard_filteredのASVを減らす、asgard_filteredを使おう！

asgard_sinking_Survey <- asgard_filtered[, particleass_no_na3] #181*42

asgard_beta_sinking = asgard_filtered^.25
asgard_beta = asgard_beta[,colSums(asgard_beta) > 0]

# euclidian, bray, jaccardの3つ作る
asgard_eucmat = vegdist(asgard_beta, method="euclidean")
asgard_braymat = vegdist(asgard_beta, method="bray")
asgard_jacmat = vegdist(asgard_beta, method="jaccard")

# hmも作れる
library(pheatmap)  # For heatmap visualization
library(RColorBrewer)  # For color gradients

pheatmap(
  asgard_braymat, 
  clustering_method = "ward",  # Ward clustering
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), 
  main = "Bray-Curtis Distance Heatmap",
  border_color = NA
)

### PCoA計算する
# dist型に変換
asgard_braydist <- as.dist(asgard_braymat)
asgard_jacdist <- as.dist(asgard_jacmat)
asgard_eucdist <- as.dist(asgard_eucmat)

# PCoAする (dimensional reduction)
asgard_pcoa_bray <- cmdscale(asgard_braydist, k = 2, eig = TRUE)
asgard_pcoa_jaccard <- cmdscale(asgard_jacdist, k = 2, eig = TRUE)
asgard_pcoa_euclidean <- cmdscale(asgard_eucdist, k = 2, eig = TRUE)

# Create a dataframe for ggplot
asgard_sinking_Survey_pcoa <- asgard_filtered[, particleass_no_na3] #181*42


suspended_no_na2 <- suspended_no_na[
  suspended_no_na %in% colnames(asgard_filtered)
]
asgard_suspended_Survey_pcoa <- asgard_filtered[, suspended_no_na2] #181*34

freeliving_no_na2 <- freeliving_no_na[
  freeliving_no_na %in% colnames(asgard_filtered)
]
asgard_freeliving_Survey_pcoa <- asgard_filtered[, freeliving_no_na2] #181*71


asgard_pcoa_df_sinking <- data.frame(
  Sample = rownames(asgard_sinking_Survey_pcoa), #行名をとるだけなので、asgard_filteredでも可
  PCoA1_Bray = asgard_pcoa_bray$points[, 1],
  PCoA2_Bray = asgard_pcoa_bray$points[, 2],
  PCoA1_Jaccard = asgard_pcoa_jaccard$points[, 1],
  PCoA2_Jaccard = asgard_pcoa_jaccard$points[, 2],
  PCoA1_Euclidean = asgard_pcoa_euclidean$points[, 1],
  PCoA2_Euclidean = asgard_pcoa_euclidean$points[, 2]
) #181*7

asgard_pcoa_df_suspended <- data.frame(
  Sample = rownames(asgard_suspended_Survey_pcoa), #行名をとるだけなので、asgard_filteredでも可
  PCoA1_Bray = asgard_pcoa_bray$points[, 1],
  PCoA2_Bray = asgard_pcoa_bray$points[, 2],
  PCoA1_Jaccard = asgard_pcoa_jaccard$points[, 1],
  PCoA2_Jaccard = asgard_pcoa_jaccard$points[, 2],
  PCoA1_Euclidean = asgard_pcoa_euclidean$points[, 1],
  PCoA2_Euclidean = asgard_pcoa_euclidean$points[, 2]
) 

asgard_pcoa_df_freeliving <- data.frame(
  Sample = rownames(asgard_freeliving_Survey_pcoa), #行名をとるだけなので、asgard_filteredでも可
  PCoA1_Bray = asgard_pcoa_bray$points[, 1],
  PCoA2_Bray = asgard_pcoa_bray$points[, 2],
  PCoA1_Jaccard = asgard_pcoa_jaccard$points[, 1],
  PCoA2_Jaccard = asgard_pcoa_jaccard$points[, 2],
  PCoA1_Euclidean = asgard_pcoa_euclidean$points[, 1],
  PCoA2_Euclidean = asgard_pcoa_euclidean$points[, 2]
) 

meta_asgard <- meta_asgard %>%
  mutate(Sample = rownames(asgard_filtered)) 

meta_asgard_sinking = meta_asgard
meta_asgard_sinking$pct.sinking <- rowSums(asgard_sinking_Survey_pcoa)

meta_asgard_suspended = meta_asgard
meta_asgard_suspended$pct.suspended <- rowSums(asgard_suspended_Survey_pcoa)

meta_asgard_freeliving = meta_asgard
meta_asgard_freeliving$pct.freeliving <- rowSums(asgard_freeliving_Survey_pcoa)


asgard_pcoa_df_sinking <- left_join(meta_asgard_sinking, asgard_pcoa_df_sinking, by = "Sample")
asgard_pcoa_df_sinking$cluster = as.factor(clusnum)

asgard_pcoa_df_suspended <- left_join(meta_asgard_suspended, asgard_pcoa_df_suspended, by = "Sample")
asgard_pcoa_df_suspended$cluster = as.factor(clusnum)

asgard_pcoa_df_freeliving <- left_join(meta_asgard_freeliving, asgard_pcoa_df_freeliving, by = "Sample")
asgard_pcoa_df_freeliving$cluster = as.factor(clusnum)


pdf("PA_boxplots.pdf")

g1 = ggplot(asgard_pcoa_df_sinking, aes(x = cluster, y = pct.sinking, fill = cluster)) +
  geom_boxplot() +
  scale_fill_manual(values = rsc) +   # Clusterごとに色を割り当て
  ylab("% Sinking") +
  theme_bw()

g2 = ggplot(asgard_pcoa_df_suspended, aes(x = cluster, y = pct.suspended, fill = cluster)) +
  geom_boxplot() +
  scale_fill_manual(values = rsc) +   # Clusterごとに色を割り当て
  ylab("% Suspended") +
  theme_bw()

g3 = ggplot(asgard_pcoa_df_freeliving, aes(x = cluster, y = pct.freeliving, fill = cluster)) +
  geom_boxplot() +
  scale_fill_manual(values = rsc) +   # Clusterごとに色を割り当て
  ylab("% Freeliving") +
  theme_bw()

print(g1)
print(g2)
print(g3)

dev.off()


#pdf(file="asgard_sinking_boxplots.pdf")
#rsc = hue_pal()(nclus)[asgard_pcoa_df_sinking$cluster]

#plot(asgard_pcoa_df_sinking$PCoA1_Bray,
     #asgard_pcoa_df_sinking$PCoA2_Bray, col = rsc, pch = 19)

#for(var in colnames(asgard_pcoa_df_sinking)) {
  #gg = ggplot(asgard_pcoa_df_sinking, aes(x = cluster, y = .data[[var]])) +
    #geom_boxplot(aes(fill=cluster),
                # outlier.shape = NA
  #  )+
   # geom_jitter(width=.4,height=0)+
   # theme(text = element_text(size = 24))
  
 # print(gg)
  
  #gp = ggplot(asgard_pcoa_df) + geom_point(aes(x=PCoA1_Bray, y= PCoA2_Bray, col=cluster, size=.data[[var]]))
  #  plot(asgard_pcoa_df$PCoA1_Bray,asgard_pcoa_df$PCoA2_Bray, col = rsc, cex=pch = 19)
 # print(gp) 
#}
#dev.off()


### map

asgard_processing_ggmap <- asgard_processing2 %>% 
  select(filter, lat, lon, contains("ESV")) #81*226

bbox <- make_bbox(lon = asgard_processing_ggmap$lon, lat = asgard_processing_ggmap$lat, f = 0.1)
mapz = get_stadiamap(bbox, maptype = "stamen_terrain", zoom =4)
asgard_processing_ggmap$filter <- factor(
  asgard_processing_ggmap$filter,
  levels = c("0.2 µm", "3 µm", "20 µm")  
)

## 単独ESVの地図を見る場合
# ESV_番号で抽出パターン
apesv91 <- asgard_processing_ggmap %>%
  select(lat, lon, filter, contains("ESV_501"))
# ESVの名前で抽出パターン
apesv91 <- asgard_processing %>%
  select(lat, lon, filter, "Flavobacteriales (100); Flavobacteriaceae (100); uncultured (100); ESV_34")

colnames(apesv91)[4] <- "ESV_501"

apesv91$filter <- factor(apesv91$filter, levels = c("0.2 µm", "3 µm", "20 µm"))

map_plot_1 <- ggmap(mapz) +
  geom_point(data = apesv91, aes(x = lon, y = lat, size = ESV_501)
             ,alpha = 1) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  facet_grid(. ~ filter)

  print(map_plot_1)
  
## loopで全抽出
esv <- asgard_processing_ggmap %>%
  select(lat, lon, filter, contains("ESV_")) #78*224

esv0.2 <- esv %>%
  filter(filter == "0.2 µm") #23*224

esv_only <- esv0.2 %>%
  select(contains("ESV_")) #23*221, ESV列だけ抽出

# 全て0の列を特定
zero_cols <- names(esv_only)[colSums(esv_only != 0) == 0]

# 元データからその列を削除
esv_zero_only <- asgard_processing_ggmap %>%
  select(lat, lon, filter, any_of(zero_cols)) #78*52 (lat+lo+filter+49ASVs)
# 元は221 ASVs. zero_colsのlengthは49. よって0.2umで0以上だった175 ASVsが除かれた.
# 3と20umでabundantだった49ASVsのみ抽出

pdf("maps_pa.pdf", width = 8, height = 6)

asv_cols <- colnames(esv_zero_only)[4:ncol(esv_zero_only)]

for (asv in asv_cols) {
  
  map_plot <- ggmap(mapz) +
    geom_point(data = esv_zero_only, aes(x = lon, y = lat, size = .data[[asv]]),
      alpha = 1
    ) +
    theme(text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right",        # legendの位置を右に
          legend.key.size = unit(0.5, "cm"),# legend key（マーカー）のサイズ
          legend.text = element_text(size = 2), # legendの文字サイズ
          legend.title = element_text(size = 2)
          ) +
    facet_grid(. ~ filter)
  
  print(map_plot)  
}

dev.off()  
















# ## filterサイズの順序が合ってるか確認
# df<-asgard_processing%>%
#   select(filter, "Flavobacteriales (100); Flavobacteriaceae (100); Flavobacterium (100); ESV_501")
# 
# df2 <- df %>%
#   group_by(filter) %>%
#   summarise(mean = mean(`Flavobacteriales (100); Flavobacteriaceae (100); Flavobacterium (100); ESV_501`), na.rm = TRUE)
# #20が一番少ない
# 
# df3 <- df2 %>%
#   select(mean)
# df3 <- as.matrix(df3)
# df3 <- prop.table(df3, 2)
# view(df3)

### presence or absense
## DBO3におけるPAのabundance調べたい
## 49 ASVsだけのhm作ってみたら？
asgard_presence <- asgard_processing2 %>%
  count(station, filter, depth_type) %>%
  pivot_wider(
    names_from = filter,
    values_from = n,
    values_fill = 0
  )
  
### db-RDA (Distance-Based Redundancy Analysis)と、boxplots
# Beta diversityの行列を作る
# 必要なら、asgard_filtered_p2 <- asgard_filtered_p2 %>%
# select(-Sample)
asgard_beta_p = asgard_filtered_p2^.25
asgard_beta_p = asgard_beta_p[,colSums(asgard_beta_p) > 0] #78*221

# euclidian, bray, jaccardの3つ作る
asgard_eucmat_p = vegdist(asgard_beta_p, method="euclidean")
asgard_braymat_p = vegdist(asgard_beta_p, method="bray")
asgard_jacmat_p = vegdist(asgard_beta_p, method="jaccard")

# hmも作れる
library(pheatmap)  # For heatmap visualization
library(RColorBrewer)  # For color gradients

pheatmap(
  asgard_braymat_p, 
  clustering_method = "ward",  # Ward clustering
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), 
  main = "Bray-Curtis Distance Heatmap",
  border_color = NA
)

### PCoA計算する
# dist型に変換
asgard_braydist_p <- as.dist(asgard_braymat_p)
asgard_jacdist_p <- as.dist(asgard_jacmat_p)
asgard_eucdist_p <- as.dist(asgard_eucmat_p)

# PCoAする (dimensional reduction)
asgard_pcoa_bray_p <- cmdscale(asgard_braydist_p, k = 2, eig = TRUE)
asgard_pcoa_jaccard_p <- cmdscale(asgard_jacdist_p, k = 2, eig = TRUE)
asgard_pcoa_euclidean_p <- cmdscale(asgard_eucdist_p, k = 2, eig = TRUE)

# Create a dataframe for ggplot
asgard_pcoa_df_p <- data.frame(
  Sample = rownames(asgard_filtered_p2), #行名をとるだけなので、asgard_filteredでも可
  PCoA1_Bray = asgard_pcoa_bray_p$points[, 1],
  PCoA2_Bray = asgard_pcoa_bray_p$points[, 2],
  PCoA1_Jaccard = asgard_pcoa_jaccard_p$points[, 1],
  PCoA2_Jaccard = asgard_pcoa_jaccard_p$points[, 2],
  PCoA1_Euclidean = asgard_pcoa_euclidean_p$points[, 1],
  PCoA2_Euclidean = asgard_pcoa_euclidean_p$points[, 2]
) #78*7

meta_asgard_p2 <- meta_asgard_p2 %>%
  mutate(Sample = rownames(asgard_filtered_p2)) #78*46

asgard_pcoa_df_p <- left_join(asgard_pcoa_df_p, meta_asgard_p2, by = "Sample") #78*52
asgard_pcoa_df_p$filter = factor(asgard_pcoa_df_p$filter, levels = c("0.2 µm", "3 µm", "20 µm")) #78*52

### boxplotを作成して、hmのクラスターを分けた要因を究明する
pdf(file="asgard_boxplots_processing.pdf")
# asgard_pcoa_df作ってから！

# 4色の場合、clusnum_p(長さ78)
rsc_p <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
sample_rgb3 <- rsc_p[clusnum_p]

clusnum_p_bp <- factor(clusnum_p, levels = c("1","2","3","4"))
# clusnum_pは1-4の数字
# そのままだとinteger扱いになるので、ファクターにする

plot(asgard_pcoa_df_p$PCoA1_Bray,asgard_pcoa_df_p$PCoA2_Bray, col = sample_rgb3, pch = 19)

for(var in colnames(asgard_pcoa_df_p)) {
  gg = ggplot(asgard_pcoa_df_p, aes(x = clusnum_p_bp, y = .data[[var]])) +
    geom_boxplot(aes(fill=clusnum_p_bp),
                 #outliers = FALSE
                 outlier.shape = NA
                 #outlier.alpha = 0,
                 #outlier.size = 0, 
                 #outlier.colour = NULL,
    )+
    geom_jitter(width=.4,height=0)+
    theme(text = element_text(size = 24))
  
  print(gg)
  
  gp = ggplot(asgard_pcoa_df_p) +
    geom_point(aes(x=PCoA1_Bray, y= PCoA2_Bray, col=clusnum_p_bp, size=.data[[var]]))
  #  plot(asgard_pcoa_df$PCoA1_Bray,asgard_pcoa_df$PCoA2_Bray, col = rsc, cex=pch = 19)
  print(gp) 
}
dev.off()




### db-RDA plot for processing
# asgard_pcoa_dfの各列内の、NA数える
for(var in  c("salinity", "temp", "depth_m", "DO", "POC (ug/L)","chl (ug/l)",
              "PON (ug/L)", "SPM (ug/L)", "FlECO-AFL(mg/m^3)", "chl depth",
              "phaeo (ug/l)", "PO4(uM)", "Sil(uM)", "NO2(uM)", "NH4(uM)",
              "N+N (umol/L)", "NO3(uM)")) { 
  asgard_na_count_p <- sum(is.na(asgard_pcoa_df_p[["NO3(uM)"]])) # forループで用いるdf[[変数]]はdf$変数と同じ意味
  print(asgard_na_count_p)
}
# 現状のNAの数: POC+PON+SPM+ = 27, NO3 = 3, DO+temp+sal+FlECO = 1, depth = 0

# 必要な列だけを備えたsubset作る
asgard_pcoa_df_sub_p <- asgard_pcoa_df_p %>% # 78*12
  select(-c(secondary_number,cellid, side, project, station, station_type
            , date, CTD_number, filter, delO18, sedchla_mg_m2, mass, cast_depth, 
            depth_original, `depth_m`,`depth_type`, `bottom_depth`, lat, lon, `log_name`, `bottle`, 
            `niskin_id_meta`, `volume_filtered_l`, `transmission_pct`, `niskin_id_nutrients`,
            `poc_pon_ratio`, `poc_pon_ratio_redfield`, `notes`, `notes2`,
            `POC (ug/L)`, `chl (ug/l)`, `PON (ug/L)`, `SPM (ug/L)`, `chl depth`,
            `phaeo (ug/l)`, `PO4(uM)`, `Sil(uM)`, `NO2(uM)`, `NH4(uM)`, `N+N (umol/L)`))

# db-RDAやる前に、モデルに扱うパラメーター(PCoAのデータ)を標準化しないといけない
# scale()で標準化する。なおPCoAの値などには手を加えない(sub1)
# 1/4 transformationは、sequencing dataの標準化
asgard_pcoa_df_sub1_p = asgard_pcoa_df_p %>% 
  select(c("Sample", "PCoA1_Bray","PCoA2_Bray","PCoA1_Jaccard",
           "PCoA2_Jaccard","PCoA1_Euclidean","PCoA2_Euclidean"))

asgard_pcoa_df_sub2_p = asgard_pcoa_df_p %>% 
  select(c("temp", "salinity","DO","NO3(uM)", "FlECO-AFL(mg/m^3)")) %>% 
  scale() 

asgard_pcoa_df_sub_p = cbind(asgard_pcoa_df_sub1_p, asgard_pcoa_df_sub2_p)
#78*12

# can only use data with no NAs
# 全てがNAの列を削除
### POC等はNAを64含む。よって行を消すと、64サンプルを失う(64/181なので、およそ1/3のサンプルを失う)
### chl-aはNAを19もち、かつtemp, sal, depth, DO, NO3を含んだモデルで*だったので削除した

### ここが大問題！！
asgard_filtered_p_frt <- asgard_filtered_p2^.25
asgard_complete_cases_p = complete.cases(asgard_pcoa_df_sub_p) #NAあればFALSE, なければTRUEでlogi vecを返す
asgard_pcoa_cc_p = asgard_pcoa_df_sub_p[asgard_complete_cases_p,] #75*12, TRUEの行だけ抽出(181行から174行にまで減る)
asgard_bray_cc_p = vegdist(asgard_filtered_p_frt[asgard_complete_cases_p,],method="bray") 

# chl (ug/l)という名前だと、スペースあるせいでエラー
# clean_names()使うか、列名を変更するか、バッククォート(escの下)で囲む
colnames(asgard_pcoa_cc_p)[colnames(asgard_pcoa_cc_p) == "chl (ug/l)"] <- "chl_ug_l"

asgard_dbrda_model_p <- capscale(asgard_bray_cc_p ~ temp + salinity 
                               + DO + `NO3(uM)` + `FlECO-AFL(mg/m^3)`,
                               data = asgard_pcoa_cc_p) # `FlECO-AFL(mg/m^3)`

# 結果のプリント
summary(asgard_dbrda_model_p)

# Test significance of the db-RDA model (999 permutations)
# anova.cca()で、CCA(Canonical Correspondence Analysis, 正準対応分析)やRDAを順列検定(follow-up test)
# CCAとRDAはともに、多変量データと環境要因の相関を見るが、CCAは非線形、RDAは線形
asgard_anova_dbrda_p <- anova.cca(asgard_dbrda_model_p, permutations = 999)
print(asgard_anova_dbrda_p) #0.001***

###大事!! Test significance of each environmental variable
# by = "margin"で各説明変数の寄与を評価するために、順列検定を実行
asgard_anova_env_p <- anova.cca(asgard_dbrda_model_p, by = "margin", permutations = 999)
print(asgard_anova_env_p) # salは0.001 ***, flECOは.

# Test significance of db-RDA axes
# by = "axis"で各軸に対して個別に順列検定を実行
# axesなので、CAP1-4が出る
# CAP = Canonical Axis Principal
asgard_anova_axes_p <- anova.cca(asgard_dbrda_model_p, by = "axis", permutations = 999)
print(asgard_anova_axes_p)

# Adjusted R² value (accounts for number of predictors)
RsquareAdj(asgard_dbrda_model_p) #r^2は0.1616333, adjは0.100882

# Plot results (PCAみたい）
# Extract db-RDA site scores (sample coordinates)
asgard_dbrda_scores_p <- as.data.frame(scores(asgard_dbrda_model_p, display = "sites"))
# scores()は、多変量解析の結果からスコア（座標）を抽出するために使用
asgard_dbrda_scores_p$MDS1 = asgard_dbrda_model_p$CA$u[,1]
asgard_dbrda_scores_p$MDS2 = asgard_dbrda_model_p$CA$u[,2]

asgard_dbrda_scores_p$Sample = rownames(asgard_dbrda_scores_p)

# Extract sample scores (for points)
asgard_dbrda_vectors_p <- as.data.frame(scores(asgard_dbrda_model_p, display = "bp"))  # "bp" = biplot arrows
asgard_dbrda_vectors_p$Variable <- rownames(asgard_dbrda_vectors_p)  # 列追加, Add variable names

asgard_dbrda_df_p = merge(asgard_dbrda_scores_p, asgard_pcoa_cc_p, by="Sample", sort=FALSE)

# add filter column

df_75_for_dbrda <- meta_asgard_p2 %>%
  filter(!Sample %in% c("BOX_6_6", "BOX_6_7", "BOX_6_26"))

filter_75 <- df_75_for_dbrda$filter #length is 75
names(filter_75) <- rownames(df_75_for_dbrda)
filter_75 <- factor(filter_75, levels = c("0.2 µm", "3 µm", "20 µm"))
#identical(rownames(df_75_for_dbrda), names(filter_75)) 行名が一致しているかチェック

asgard_dbrda_merged_p = left_join(asgard_dbrda_df_p, rownames_to_column(as.data.frame(filter_75)), by = c("Sample" = "rowname"))
asgard_dbrda_merged_p$filter = factor(asgard_dbrda_merged_p$filter, levels = c("0.2 µm", "3 µm", "20 µm"))

# 4色db-RDA用 (長さは75)
clusnum_p_db <- clusnum_p[!names(clusnum_p) %in% c("BOX_6_6", "BOX_6_7", "BOX_6_26")]
clusnum_p_db <- factor(clusnum_p_db, levels = c("1", "2", "3", "4")) # length is 75


GroupVariable = "salinity" 

# Plot db-RDA in ggplot2
ggplot() +
  # Plot points, colored by group
  geom_point(data = asgard_dbrda_merged_p, 
             aes(x = MDS1, y = MDS2, color = clusnum_p_db, size = `NO3(uM)`)
  ) + 
  #size = `NO3(uM)`みたいにdotのサイズ指定可
  #.data[[GroupVariable]]
  # xyを、CAP1と2ではなく、MDS1と2に変更可能
  
  # Add environmental arrows
  geom_segment(data = asgard_dbrda_vectors_p, 
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black", size = 1) +
  
  # Add environmental variable labels near arrow tips
  geom_text(data = asgard_dbrda_vectors_p, 
            aes(x = CAP1, y = CAP2, label = Variable), 
            vjust = -0.5, hjust = 0.5, size = 5) +
  
  # Labels and theme
  labs(x = "MDS1", y = "MDS2", color = "filter") +
  theme_minimal()
# title = "db-RDA Ordination with Environmental Vectors"

### temp vs. sal plots
a = meta_asgard_p2 
a$cluster = clusnum_p

plot(x = a$salinity, 
     y = a$temp,
     col = sample_rgb3,
     pch = 19,
     cex = 1.5
     )

plot(x = a$salinity, 
     y = a$`NO3(uM)`,
     col = sample_rgb3,
     pch = 19,
     cex = 1.5
)

### importing 18S data
euk_class = read.table("class_relabund_by_station_depth.tsv",sep="\t",head=T,row.names=1) # 851*122
euk_phylum = read.table("phylum_relabund_by_station_depth.tsv",sep="\t",head=T,row.names=1) # 851*61

# extract P
euk_class_p <- euk_class %>%
  filter(station_type == "P") #291*122
euk_phylum_p <- euk_phylum %>%
  filter(station_type == "P") #291*61

# remove "station", "depth_m", "station_type", "filter" columns from the original table
euk_class_p <- euk_class_p %>%
  select(-c("station", "depth_m", "station_type", "filter")) # 291*118

a = euk_class_p %>%
  select(-Sample)

rowSums(a)

b = euk_phylum_p %>%
  select(-Sample)

rowSums(b)

euk_phylum_p <- euk_phylum_p %>%
  select(-c("station", "depth_m", "station_type", "filter")) # 291*57


# extract ASGARD sample (same as meta_asgard_p2)
#vec = rownames(meta_asgard_p2) #78
#euk_class_p <- euk_class_p[vec, ] #78*122
#euk_phylum_p <- euk_phylum_p[vec, ] #78*61

# make Sample column to left_join()
euk_class_p$Sample <- rownames(euk_class_p)
euk_phylum_p$Sample <- rownames(euk_phylum_p)

meta_asgard_p_euk_class <- 
  left_join(meta_asgard_p2, euk_class_p, 
            by = "Sample") # 78*164

rownames(meta_asgard_p_euk_class) <- meta_asgard_p_euk_class$Sample

meta_asgard_p_euk_phylum <- 
  left_join(meta_asgard_p2, euk_phylum_p, 
            by = "Sample") # 78*103

rownames(meta_asgard_p_euk_phylum) <- meta_asgard_p_euk_phylum$Sample

## boxplots
euk_clusnum_p = factor(clusnum_p, levels = c("1", "2", "3", "4"))

# ASVだけ引っ張る
num_cols_class <- sapply(meta_asgard_p_euk_class, is.numeric)
num_df_class <- meta_asgard_p_euk_class[, num_cols_class]

num_cols_phylum <- sapply(meta_asgard_p_euk_phylum, is.numeric)
num_df_phylum <- meta_asgard_p_euk_phylum[, num_cols_phylum]

# 18Sのclassのboxplots
pdf(file="18Sclass_boxplots_processing.pdf")

for(var in colnames(num_df_class)) {
  gg = ggplot(num_df_class, aes(x = euk_clusnum_p, y = .data[[var]])) +
    geom_boxplot(aes(fill=euk_clusnum_p),
                 #outliers = FALSE
                 outlier.shape = NA
                 #outlier.alpha = 0,
                 #outlier.size = 0, 
                 #outlier.colour = NULL,
    )+
    geom_jitter(width=.4,height=0)+
    scale_y_log10()+
    theme(text = element_text(size = 24))
  
  print(gg)
}

dev.off()

# 18Sのphylumのboxplots
pdf(file="18Sphylum_boxplots_processing.pdf")

for(var in colnames(num_df_phylum)) {
  gg = ggplot(num_df_phylum, aes(x = euk_clusnum_p, y = .data[[var]])) +
    geom_boxplot(aes(fill=euk_clusnum_p),
                 #outliers = FALSE
                 outlier.shape = NA
                 #outlier.alpha = 0,
                 #outlier.size = 0, 
                 #outlier.colour = NULL,
    )+
    geom_jitter(width=.4,height=0)+
    scale_y_log10()+
    theme(text = element_text(size = 24))
  
  print(gg)
}

dev.off()


## 18Sのhm
asgard_euk_class_hm <- num_df_class[, -(1:26)]
asgard_euk_class_hm <- as.matrix(asgard_euk_class_hm) #78*118

# 各列の合計を計算
col_sums <- colSums(asgard_euk_class_hm, na.rm = TRUE)

# 0より大きい列名だけ抽出
asgard_euk_class_hm_filtered <- asgard_euk_class_hm[, col_sums > 0] #78*66

# NAの行を削除(36, 38, 41, 73行目)
asgard_euk_class_hm_filtered <- asgard_euk_class_hm_filtered[complete.cases(asgard_euk_class_hm_filtered), ] # 74*66
sample_rgb4 <- sample_rgb3[-c(36, 38, 41, 73)] #74


asgard_filtered_p_hm3 <- asgard_filtered_p_hm2[rownames(asgard_euk_class_hm_filtered), ]

pdf(file="ASGARD_hm_processing_18S.pdf", width = 20, height = 20) #ASV全表示したいなら, width = 20, height = 10

h3=heatmap.2((asgard_filtered_p_hm3)^.25, 
             distfun=function(x) vegdist(x, method="bray"),
             hclustfun=function(x) hclust(x, method="ward.D"),
             col = viridis, 
             #RowSideColors = sample_rgb2,
             #ColSideColors = asv_rgb2,
             margins = c(15,15),
             scale="none",
             main="ASGARD_bray/ward.D2",
             trace="none",
             cexCol = 0.2, #列の文字サイズ
             cexRow = 0.2, #行の文字サイズ
             labRow = meta_asgard_p2$side     
             #labCol = NA
)
nclus_p_euk = 4 
oldclus_p_euk = cutree(as.hclust(h3$rowDendrogram),k=nclus_p_euk)
oldorder_p_euk = unname(rle(oldclus_p_euk[as.hclust(h3$rowDendrogram)$order])$values)
neworder_p_euk = (1:nclus_p_euk)
names(neworder_p_euk) = oldorder_p_euk
clusnum_p_euk = unname(neworder_p_euk[as.character(oldclus_p_euk)])
names(clusnum_p_euk) = names(oldclus_p_euk)

rsc_p_euk <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
sample_rgb4 <- rsc_p_euk[clusnum_p_euk]

h10=heatmap.2((asgard_euk_class_hm_filtered)^.25, 
             distfun=function(x) vegdist(x, method="bray"),
             hclustfun=function(x) hclust(x, method="ward.D"),
             col = viridis, 
             RowSideColors = sample_rgb4,
             Rowv = h3$rowDendrogram,
             #ColSideColors = asv_rgb2,
             margins = c(15,15),
             scale="none",
             main="ASGARD_bray/ward.D2",
             trace="none",
             cexCol = 0.8, #列の文字サイズ
             cexRow = 0.8 #行の文字サイズ
             #labRow = meta_asgard_p2_49$side
             #labCol = NA
)

dev.off()


### h11: New heatmap from esv_relabund_by_station_depth.tsv (in the style of h10)

# Read ESV-level relative abundance table (16S, by station and depth)
esv_relabund <- read.table("esv_relabund_by_station_depth.tsv",
                            sep = "\t", header = TRUE, row.names = 1)
# Expected: ~858 rows x 5443 columns (5 metadata + ESVs)

# Filter for Processing stations and keep only ESV columns
esv_relabund_p <- esv_relabund %>%
  filter(station_type == "P") %>%
  select(contains("ESV"))

# Keep only ASGARD samples already used in the 18S heatmap (asgard_euk_class_hm_filtered: 74 samples)
common_esv_samples <- intersect(rownames(asgard_euk_class_hm_filtered),
                                 rownames(esv_relabund_p))
esv_asgard_p <- esv_relabund_p[common_esv_samples, ]

# Apply the same minimum-abundance filtering used in the 16S processing analysis
esv_mincutoff <- apply(esv_asgard_p, 2, max) > 0.001
esv_asgard_filt <- esv_asgard_p[, esv_mincutoff]
esv_asgard_filt <- esv_asgard_filt[, colSums(esv_asgard_filt > 0) > 2]

# Convert to matrix
esv_asgard_mat <- as.matrix(esv_asgard_filt)

# Row side colors aligned to existing cluster assignments from h10 (sample_rgb4)
sample_rgb_esv <- sample_rgb4[common_esv_samples]

pdf(file = "ASGARD_hm_processing_esv_relabund.pdf", width = 20, height = 20)

h11 <- heatmap.2(
  (esv_asgard_mat)^.25,
  distfun      = function(x) vegdist(x, method = "bray"),
  hclustfun    = function(x) hclust(x, method = "ward.D"),
  col          = viridis,
  RowSideColors = sample_rgb_esv,
  Rowv         = h3$rowDendrogram,
  margins      = c(15, 15),
  scale        = "none",
  main         = "ASGARD ESV relabund by station/depth (Bray/ward.D)",
  trace        = "none",
  cexCol       = 0.8,
  cexRow       = 0.8
)

dev.off()

