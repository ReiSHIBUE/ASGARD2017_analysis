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

seqtab_ODV_Survey_full = seqtab_ODV_full[seqtab_ODV_full$station_type=="S",]
# 689*3122
# 0.2 filterはstation SとP含むが、3と20はPしかない
# よってstation typeをSでフィルターすれば、必然とSだけになる


### 1クルーズ(ASGARD)の全微生物の解析をする
# ASGARDだけのdf作る
# seqtab_ODV_Surveyなので、フィルターサイズは0.22umに限定、かつtop100だけ
# ASVの名前はfamilyとgenus
ASGARD <- seqtab_ODV_Survey_full %>%
  filter(project == "ASGARD2017") # 181*3122

# ASGARDだけを抽出 (他のクルーズでも適用可！)
asgard_vec <- ASGARD$Row.names #length is 181
seqtab_prop_asgard <- seqtab_16Sprop[asgard_vec, ] #181*3076

# ASGARDのmeta data抽出
meta_asgard <- meta_denovo_2[asgard_vec, ] #181*45

### prop tableの種類！
# 普通のprop tableはseqtab_16Spropcol_mra100mat
# そこからasgardだけ抽出したのは、seqtab_prop_asgard (181*3076)

# Set a minimum read count threshold (e.g., ASVs with <5 total reads are removed)
mincutoff = (apply(seqtab_prop_asgard, 2, max) > 0.001)
asgard_filtered <- seqtab_prop_asgard[, mincutoff ]
asgard_filtered <- asgard_filtered[, colSums(asgard_filtered > 0) > 2]
hist(rowSums(asgard_filtered),breaks=20)
dim(asgard_filtered) #181*258

ODV.otu.pick_3 = match(colnames(asgard_filtered), colnames(seqtab_filt)) #length is 258
colnames(asgard_filtered) = shorternames[ODV.otu.pick_3] 

###### Make a fourth-root transformed proportional table
# ヒートマップに用いるmatrixは標準化されてないとダメ
# hm()の中で、df^.25としてればok!!


### シーケンス深度
# prop tableではなく、リード数が必要
# asgardの各配列の数(propではなく、実際のseqのカウント数)
asgard_seqcount <- seqtab_16Smat[asgard_vec, ]

# Compute total reads per sample
seq_depth <- rowSums(asgard_seqcount)
# シーケンス深度は何回シーケンスするか
# コイントスを3回やった、表が2回でた、よって表が出る確率は2/3だ！
# これはおかしい、100回くらいやって初めて統計学的に信頼できる
# シーケンスも同じ、数回ではエラーが起きたり、ある遺伝子だけが多く読まれたりしちゃう
# よって、シーケンス深度が高いほど、ゲノムのリード回数が多くなり、より正確で信頼性の高い情報が得られる
# 30倍なら、30回読むこと。

# Create a sequencing depth plot
ggplot(data.frame(Sample = rownames(asgard_seqcount), Reads = seq_depth), 
       aes(x = Sample, y = Reads)) +
  geom_point(color = "blue") +
  geom_hline(yintercept = min(seq_depth), linetype = "dashed", color = "red") +
  labs(title = "Sequencing Depth per Sample",
       x = "Sample",
       y = "Total Reads") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# What is an appropriate sequencing depth?
hist(seq_depth,breaks=20) # ~くらいでちょうどいい

# Compute rarefaction curves
# Rarefaction Curveは、十分なシーケンスリード数を示す
# アルファ多様性がプラトーに達した = それ以上リード数を増やしてもアルファ多様性は増えない
rarecurve(asgard_seqcount, step = 100, col = "blue", cex = 0.5, label = FALSE)
abline(v = min(seq_depth), col = "red", lty = 2)  # Show minimum sequencing depth
# sequencing depthが~くらいあれば、シーケンシングの信頼度が高いと言える

### ASGARDのhm
library(viridis)

pdf(file="ASGARD_hm.pdf")
h1=heatmap((asgard_filtered)^.25, 
           distfun=function(x) vegdist(x, method="bray"),
           hclustfun=function(x) hclust(x, method="ward.D"),
           keep.dendro=TRUE, 
           scale="none",
           margins = c(20,20),
           main="ASGARD_bray/ward.D2")

nclus = 5 #枝が5くらいに分けられそうだった
oldclus = cutree(as.hclust(h1$Rowv),k=nclus)
oldorder = unname(rle(oldclus[as.hclust(h1$Rowv)$order])$values)
neworder = (1:nclus)
names(neworder) = oldorder
clusnum = unname(neworder[as.character(oldclus)])
names(clusnum) = names(oldclus)

colnclus = 8 #枝が5くらいに分けられそうだった
cololdclus = cutree(as.hclust(h1$Colv),k=colnclus)
cololdorder = unname(rle(cololdclus[as.hclust(h1$Colv)$order])$values)
colneworder = (1:colnclus)
names(colneworder) = cololdorder
colclusnum = unname(colneworder[as.character(cololdclus)])
names(colclusnum) = names(cololdclus)



#rsc=viridis(nclus)[clusnum] #rscではなくviridisにしてもいい
#colrsc=viridis(colnclus)[colclusnum] #rscではなくviridisにしてもいい

rsc=hue_pal()(nclus)[clusnum]
colrsc=hue_pal()(colnclus)[colclusnum]
colrsc=plasma(colnclus)[colclusnum] #列クラスターを行クラスターの色を変えたいなら


h2=heatmap.2((asgard_filtered)^.25, 
             distfun=function(x) vegdist(x, method="bray"),
             hclustfun=function(x) hclust(x, method="ward.D"),
             col = viridis, 
             RowSideColors = rsc,
             ColSideColors = colrsc,
             margins = c(15,15),
             scale="none",
             main="ASGARD_bray/ward.D2",
             trace="none",
             cexCol = 0.2, #列の文字サイズ
             cexRow = 0.5 #行の文字サイズ
             #labRow = NA,        
             #labCol = NA
)


h5=heatmap.2((asgard_filtered[,colclusnum==8])^.25, 
             distfun=function(x) vegdist(x, method="euc"),
             hclustfun=function(x) hclust(x, method="ward.D"),
             col = viridis, 
             RowSideColors = rsc,
             Rowv = h2$rowDendrogram,
             labRow = meta_asgard$depth_type,
#             ColSideColors = colrsc,
             margins = c(15,15),
             scale="none",
             main="ASGARD_bray/ward.D2",
             trace="none",
             cexCol = 0.2, #列の文字サイズ
             cexRow = 0.1, #行の文字サイズ
             labCol = NULL
)


newclus=cutree(as.hclust(h5$colDendrogram),k = 3)
newrsc=hue_pal()(3)[newclus]

h5=heatmap.2((asgard_filtered[,colclusnum==8])^.25, 
             distfun=function(x) vegdist(x, method="euc"),
             hclustfun=function(x) hclust(x, method="ward.D"),
             col = viridis, 
             RowSideColors = rsc,
             Rowv = h2$rowDendrogram,
             ColSideColors = newrsc,
             labRow = meta_asgard$depth_type,
             margins = c(15,15),
             scale="none",
             #main="ASGARD_bray/ward.D2",
             trace="none",
             cexCol = 0.2, #列の文字サイズ
             cexRow = 0.5, #行の文字サイズ
             labCol = NULL
             )

missing_asvs = colnames(asgard_filtered[,colclusnum==8][,newclus==2])

dev.off() 

# PNG (背景透明色で作る)
png("transparent_hm.png", 
    width = 800, 
    height = 600, 
    bg = "transparent")  # 背景を透明に

h2=heatmap.2((asgard_filtered)^.25, 
             distfun=function(x) vegdist(x, method="bray"),
             hclustfun=function(x) hclust(x, method="ward.D"),
             col = viridis, 
             RowSideColors = rsc,
             ColSideColors = colrsc,
             margins = c(15,15),
             scale="none",
             #main="ASGARD_bray/ward.D2",
             trace="none",
             key = FALSE,
             cexCol = 0.2, #列の文字サイズ
             cexRow = 0.5, #行の文字サイズ
             labRow = NA,        
             labCol = NA
)

dev.off()


# h5 = missingの、各クラスターにどのASVがいるか把握する
# これを使うなら、asgard_filteredの列名は配列であるべき！
fullnameboot[missing_asvs,]

View(fullnameboot[colnames(asgard_filtered[,colclusnum==8][,newclus==1]),])
View(fullnameboot[colnames(asgard_filtered[,colclusnum==8][,newclus==3]),])

View(fullnameboot[colnames(asgard_filtered[,colclusnum==5][,newclus==2]),])
View(fullnameboot[colnames(asgard_filtered[,colclusnum==5][,newclus==3]),])


# nclusは5, newclusはmissing ASVのhmの列用で3しかない

### h2の縦横ASVs比較
h2_colclus <- cutree(as.hclust(h2$colDendrogram), k = 8)  # 列のクラスター
h2_rowclus <- cutree(as.hclust(h2$rowDendrogram), k = 5)  # 列のクラスター
View(fullnameboot[colnames(asgard_filtered[,h2_colclus == 1][h2_rowclus == 1,]),])
View(fullnameboot[colnames(asgard_filtered[,h2_colclus == 2][h2_rowclus == 2,]),])
View(fullnameboot[colnames(asgard_filtered[,h2_colclus == 5][h2_rowclus == 3,]),])
View(fullnameboot[colnames(asgard_filtered[,h2_colclus == 4][h2_rowclus == 4,]),])
View(fullnameboot[colnames(asgard_filtered[,h2_colclus == 1][h2_rowclus == 1,]),])
View(fullnameboot[colnames(asgard_filtered[,h2_colclus == 8][h2_rowclus == 1,]),])

# ピンク(during), 5C
View(fullnameboot[colnames(asgard_filtered[,h2_colclus == 3][h2_rowclus == 5,]),])
# 緑 (pre), 3C
View(fullnameboot[colnames(asgard_filtered[,h2_colclus == 3][h2_rowclus == 3,]),])
View(fullnameboot[colnames(asgard_filtered[,h2_colclus == 5][h2_rowclus == 3,]),])


# 青 (post), 4C
View(fullnameboot[colnames(asgard_filtered[,h2_colclus == 3][h2_rowclus == 4,]),])
# オレンジ, misiing
View(fullnameboot[colnames(asgard_filtered[,h2_colclus == 1][h2_rowclus == 1,]),])

# 3C
View(fullnameboot[colnames(asgard_filtered[,h2_colclus == 3][h2_rowclus == 3,]),])
# 4C
View(fullnameboot[colnames(asgard_filtered[,h2_colclus == 3][h2_rowclus == 4,]),])



t1 = fullnameboot[colnames(asgard_filtered[,h2_colclus == 3][h2_rowclus == 3,]),]
t2 = fullnameboot[colnames(asgard_filtered[,h2_colclus == 3][h2_rowclus == 4,]),]
t3 = fullnameboot[colnames(asgard_filtered[,h2_colclus == 3][h2_rowclus == 5,]),]
# 行結合
combined_df <- rbind(t1, t2, t3)
# 重複除去
unique_df <- combined_df[!duplicated(combined_df$ESV), ]


# ESVも1000番台なので、少しはASVが変わる！
# 行は下から1-5, 列は右から1-8





#列6
View(fullnameboot[colnames(asgard_filtered[,h2_newclus==6]),])
View(fullnameboot[colnames(asgard_filtered[,colclusnum==6]),])
#列1
View(fullnameboot[colnames(asgard_filtered[,colclusnum==8]),])


dt <- as.POSIXct(meta_asgard$date, format = "%b %d %Y %H:%M", tz = "UTC")
jday <- as.numeric(format(dt, "%j"))


### JPEGとして保存する場合
# x軸の名前を削除する
library(viridis)
library(scales) #hue_pal()のため。ggplotと同じ色。

jpeg(filename = "ASGARD_hm.jpeg",
     width = 8,                         # 幅（インチ）
     height = 6,                        # 高さ（インチ）
     units = "in",                      # 単位（インチ）
     res = 2000,                        # 解像度（dpi）
     quality = 100                      # 画質（0-100）
)

h2=heatmap.2((asgard_filtered)^.25, 
             distfun=function(x) vegdist(x, method="bray"),
             hclustfun=function(x) hclust(x, method="ward.D"),
             col = viridis, 
             RowSideColors = rsc,
             ColSideColors = colrsc,
             margins = c(15,15),
             scale="none",
             #main="ASGARD_bray/ward.D2",
             trace="none",
             #cexCol = 0.2, #列の文字サイズ
             #cexRow = 0.5, #行の文字サイズ
             labRow = NA,        
             labCol = NA,
             key = FALSE
)

dev.off() 

# jpeg(filename = "ASGARD_hm.jpeg",
#      width = 8,                         # 幅（インチ）
#      height = 6,                        # 高さ（インチ）
#      units = "in",                      # 単位（インチ）
#      res = 2000,                        # 解像度（dpi）
#      quality = 100                      # 画質（0-100）
#      )
# 
# rsc = hue_pal()(5)[clusnum]
# # rsc=rsc(nclus)[clusnum]が元のコード！
# 
# h2=heatmap.2(asgard_filtered^.25, 
#              distfun=function(x) vegdist(x, method="chord"),
#              hclustfun=function(x) hclust(x, method="ward.D"),
#              col = viridis, 
#              RowSideColors = rsc,
#              margins = c(15,15),
#              scale="none",
#              #main="Microbial Abundance in Each Samples",
#              trace="none",
#              cexCol = 0.2, #列の文字サイズ
#              cexRow = 0.5, #行の文字サイズ
#              labCol = NA 
# )
# 
# dev.off() 
# 
# 
# # 新しいhmのjpeg
# jpeg(filename = "ASGARD_hm_2.jpeg",
#      width = 8,                         # 幅（インチ）
#      height = 6,                        # 高さ（インチ）
#      units = "in",                      # 単位（インチ）
#      res = 2000,                        # 解像度（dpi）
#      quality = 100                      # 画質（0-100）
# )
# 
# h2=heatmap.2((asgard_filtered)^.25, 
#              distfun=function(x) vegdist(x, method="bray"),
#              hclustfun=function(x) hclust(x, method="ward.D"),
#              col = viridis, 
#              RowSideColors = rsc,
#              ColSideColors = colrsc,
#              margins = c(15,15),
#              scale="none",
#              #main="ASGARD_bray/ward.D2",
#              trace="none",
#              #cexCol = 0.2, #列の文字サイズ
#              #cexRow = 0.5, #行の文字サイズ
#              labRow = NA,        
#              labCol = NA,
#              key=FALSE
# )
# dev.off() 
             
             
             
### ASGARDのcorhm
# 相関係数行列をもとに、hm
#library(vegan)
#library(fastcluster)
#library(gplots)

#transform using the the fourth root (4乗根, 1/4はone fourth) transformation (値のばらつきを小さくする！)
#ASGARD_ASV_filtered_frt <- ASGARD_ASV_filtered^.25
#ASGARD_ASV_filtered_frt <- (asgard_filtered)^.25
#ASGARD_cor <- cor(ASGARD_ASV_filtered_frt,method = "spearman") 
#ASGARD_cor <- cor(ASGARD_ASV_filtered_frt,method = "pearson") 


#pdf(file="ASGARD_corhm.pdf",width=24,height=24)

#my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 99) #カラーパレット作成
#breaks <- seq(-1, 1, length.out = 100)
#h3 = heatmap.2(ASGARD_cor,
          # trace = "none",
          # key = FALSE,
          # margins = c(40,40), #余白広めにする
          # col = my_palette,
          # breaks = breaks,
          # distfun=function(x) vegdist(x,method="gower"), 
          # hclustfun=function(x) hclust(x,method="ward"),
          # scale='none',
          # cexRow = 0.5,
          # cexCol = 0.5)

# nclus = 8 #枝が5くらいに分けられそうだった
# oldclus = cutree(as.hclust(h3$rowDendrogram),k=nclus)
# oldorder = unname(rle(oldclus[as.hclust(h3$rowDendrogram)$order])$values)
# neworder = (1:nclus)
# names(neworder) = oldorder
# clusnum = unname(neworder[as.character(oldclus)])
# names(clusnum) = names(oldclus)
# 
# colrsc=hue_pal()(nclus)[clusnum]
# 
# h3 = heatmap.2(ASGARD_cor,
#                trace = "none",
#                key = FALSE,
#                margins = c(40,40), #余白広めにする
#                col = my_palette,
#                breaks = breaks,
#                RowSideColors = colrsc,
#                distfun=function(x) vegdist(x,method="gower"), 
#                hclustfun=function(x) hclust(x,method="ward"),
#                scale='none',
#                cexRow = 0.5,
#                cexCol = 0.5)
# 
# 
# 
# h2=heatmap.2(asgard_filtered^.25, 
#              distfun=function(x) vegdist(x, method="bray"),
#              hclustfun=function(x) hclust(x, method="ward.D2"),
#              col = viridis, 
# #             Colv = h3$colDendrogram,
#              RowSideColors = rsc,
#              ColSideColors = colrsc,
#              margins = c(15,15),
#              scale="none",
#              main="ASGARD_bray/ward.D2",
#              trace="none",
#              cexCol = 0.2, #列の文字サイズ
#              cexRow = 0.5 #行の文字サイズ
# )
# dev.off()
# 
# # brayではなく、gowerにした！
# # brayだとWarning: results may be meaningless because data have negative entries in method “bray”
# 
# # 相関なので、範囲を-1から1にする
# # ASGARD_cor_df <- as.data.frame(ASGARD_cor)
# 
# pdf(file="ASGARD_corhm.pdf",width=24,height=24)
# my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 99) #カラーパレット作成
# breaks <- seq(-1, 1, length.out = 100)
# heatmap.2(ASGARD_cor,
#           trace = "none",
#           key = FALSE,
#           margins = c(40,40), #余白広めにする
#           col = my_palette,
#           breaks = breaks,
#           distfun=function(x) vegdist(x,method="gower"),
#           hclustfun=function(x) hclust(x,method="ward"),
#           scale='none',
#           cexRow = 0.5,
#           cexCol = 0.5)
# 
# dev.off()
# 
# # 色が青 = 相関係数が大きい = 正の相関
# # 色が赤 = 相関係数が小さい = 負の相関
# 
# # ASGARDのstacked bar plots作る
# # 温度や塩分濃度を因子にして、各パラメーター毎のabundance調べる
# # ODVで良くね？
# # 既存のdfを使って、facet_wrapでクルーズ毎の図を作れないかな？
# 
# ### corhmのJPEG
# jpeg(filename = "ASGARD_corhm_key.jpeg",
#      width = 8,                         # 幅（インチ）
#      height = 6,                        # 高さ（インチ）
#      units = "in",                      # 単位（インチ）
#      res = 2000,                        # 解像度（dpi）
#      quality = 100                      # 画質（0-100）
# )

## 重要なASVのみ表示したい場合
# 全てのASV名
#all_labels <- colnames(ASGARD_cor)

# 7ASVsだけを行列名に表示！
# keep_labels <- c("Flavobacteriales (100); Flavobacteriaceae (100); Polaribacter (94); ESV_7",
#                  "Flavobacteriales (100); Flavobacteriaceae (100); Ulvibacter (99); ESV_9",
#                  "Flavobacteriales (100); Cryomorphaceae (100); uncultured (100); ESV_32",
#                  "Flavobacteriales (100); Flavobacteriaceae (100); Tenacibaculum (59); ESV_48",
#                  "Flavobacteriales (100); Flavobacteriaceae (100); Ulvibacter (98); ESV_139",
#                  "Flavobacteriales (100); Crocinitomicaceae (100); Fluviicola (88); ESV_201",
#                  "Sphingobacteriales (100); NS11−12_marine_group (100); NS11−12_marine_group_ge (100); ESV_68"
#                  )

# 表示用ラベル（その他は空白に）
# custom_labels <- ifelse(all_labels %in% keep_labels, all_labels, "")
# 
# my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 99) #カラーパレット作成
# breaks <- seq(-1, 1, length.out = 100)
# heatmap.2(ASGARD_cor,
#           trace = "none",
#           key = TRUE,
#           margins = c(5,5), #余白広めにする
#           col = my_palette,
#           breaks = breaks,
#           distfun=function(x) vegdist(x,method="gower"),
#           hclustfun=function(x) hclust(x,method="ward"),
#           scale='none',
#           cexRow = 0.5,
#           cexCol = 0.5,
#           labRow = NA,  # 行ラベルに反映
#           labCol = NA # 列ラベルに反映
#           )
# 
# dev.off()
# 
# ### familyとgenusのみの短いver. for corhm
# ASGARD_short <- seqtab_ODV_Survey_short %>%
#   filter(project == "ASGARD2017")
# 
# ASGARD_ASV_short <- ASGARD_short %>%
#   select(Row.names, contains("ESV_"))
# 
# asgard_vec_short <- ASGARD_short$Row.names #length is 181
# 
# seqtab_prop_asgard_short <- seqtab_16Spropcol_mra100mat2[asgard_vec_short, ] #181*100
# 
# # ASGARDのmeta data抽出
# meta_asgard_short <- meta_denovo_2[asgard_vec_short, ] #181*45
# 
# # Set a minimum read count threshold (e.g., ASVs with <5 total reads are removed)
# asgard_filtered_short <- seqtab_prop_asgard_short[, colSums(seqtab_prop_asgard_short > 0) > 5]
# rowSums(asgard_filtered_short)
# 
# asgard_frtprop_short = asgard_filtered_short^.25
# 
# ### seqではなく、ESVの名前が欲しいならこれ使う！
# ASGARD_ASV_filtered_short <- ASGARD_ASV_short[, colSums(ASGARD_ASV_short > 0) > 5]
# ASGARD_ASV_filtered_short <- ASGARD_ASV_filtered_short %>%
#   select(-Row.names)
# 
# ASGARD_cor_short <- cor(ASGARD_ASV_filtered_short)

### 環境因子毎のabundance調べる
# 微生物同士の相関行列作る、cor()

# Sum ASV abundances at the Class level
class_abundance_asgard <- aggregate(t(seqtab_prop_asgard), by=list(fullnameboot[colnames(seqtab_prop_asgard),"Order"]), FUN = sum)

# long formatにする
class_long_asgard <- pivot_longer(class_abundance_asgard, 
                                  cols = -Group.1, # Group.1以外の列をロングにする
                                  names_to = "Sample", #列名(リード配列)がSample列にまとまる
                                  values_to = "Abundance") #数値(行列の要素)がAbundance列にまとまる

# Make a fourth-root transformed proportional table
# 249行目で、hm作るとき、mat2を^.25してる
asgard_frtprop = asgard_filtered^.25

# Choose which table to go forward with
asgard_beta = asgard_frtprop
asgard_beta = asgard_beta[,colSums(asgard_beta) > 0]

# bray-curtis matrix作る
asgard_bcm = vegdist(asgard_beta, method="bray")

# Convert Bray-Curtis matrix to distance format (bray-curtisの行列をdist型に変換)
# asgard_bcm <- as.dist(asgard_bcm)

# Perform hierarchical clustering (UPGMA method)
asgard_bh <- hclust(asgard_bcm, method = "ward.D")

# Extract ordered sample names
asgard_ordered_samples <- asgard_bh$labels[asgard_bh$order]
# labelsはサンプル名、orderはindicesでその番号の位置にいるサンプルを抽出
# ordered_samplesはサンプル名の順番を変えたもの

# Convert Sample column to factor with correct order
class_long_asgard$Sample <- factor(class_long_asgard$Sample, levels = asgard_ordered_samples)
class_long_asgard$Sample <- as.factor(class_long_asgard$Sample)
# levelsはfactor()の引数で、因子のレベル(カテゴリの順序)を指定

# raを計算
class_long_asgard <- class_long_asgard %>%
  group_by(Sample) %>%
  mutate(Relative_Abundance = Abundance / sum(Abundance))

# Stacked bar plot with relative abundance
ggplot(class_long_asgard, aes(x = Sample, y = Relative_Abundance, fill = Group.1)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +  # Show as percentage
  theme_minimal() +
  labs(title = "Microbial Community Composition at Class Level from ASGARD 2017",
       x = "Sample (Ordered by Bray-Curtis Clustering)", 
       y = "Relative Abundance (%)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("ASGARD_barchart.pdf",width=96,height=96, limitsize=F)
# ピンクと青緑が多い = SAR11系統

### network analysis

# Compute Spearman correlation matrix between ASVs
# Spearman is useful because it uses ranks rather than absolute values
# ピアソン相関係数はactual valueを使う一方、スピアマンは順位使う
asgard_correlation <- cor(ASGARD_ASV_filtered, method = "spearman")
hist(asgard_correlation, breaks=100) # normally distributed

# make hm
library(pheatmap)

pdf(file="asgard_cor_hm.pdf",width=24,height=24)
pheatmap(asgard_correlation, 
         clustering_method = "average", 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cex=0.1,
         main = "ASGARD ASV Co-occurrence Heatmap")
dev.off()

# Replace weak correlations with zero (set threshold)
asgard_correlation[abs(asgard_correlation) < 0.6] <- 0  # Adjust threshold as needed

# Convert correlation matrix to adjacency matrix(隣接行列)
# 隣接行列はネットワークアナリシスの各サンプルの頂点の接続関係(ノードを結ぶ線の長さ)を示してくれる
asgard_adjacency_matrix <- asgard_correlation #dfの名前を変える
diag(asgard_adjacency_matrix) <- 0  # Remove self-connections
# diag() = diagonal = 対角成分を作る

# Create a network graph object
# グラフ理論では、ノードが点々、エッジがノードを結ぶ接続を表す
library(igraph)
asgard_network <- graph_from_adjacency_matrix(asgard_adjacency_matrix, mode = "undirected", weighted = TRUE) #グラフオブジェクトを作る
# mode = undirectedで無向グラフ(ASVから別のASVへ矢印向かない。ASV同士のノードを線で結ぶもの = ノード間の接続が双方向的)
# weighted = TRUEで重みあり(エッジに数値がついたもの)

asgard_network_layout = layout_nicely(asgard_network)
# layout_nicely()は、レイアウト(グラフの頂点を2次元または3次元空間に配置するための座標)を計算する

# Save network graphs to file
pdf(file="ASGARD_network.pdf",width=8,height=8)

# Plot the network
plot(asgard_network, 
     layout = asgard_network_layout,
     vertex.size = 2,       # ノードのサイズ
     vertex.color = "red",  # ノードの色
     vertex.label = NA,
     edge.arrow.size = 1,  # 矢印のサイズ
     edge.width = 2,         # エッジの太さ
     edge.color = "green",    # エッジの色
     main = "ASGARD_Co-occurrence Network of ASVs")

dev.off()

### Example 1: Color by dominant season
# season以外のvariable使うので、depth_typeにしてみる
# Q: これはseasonのようなカテゴリカルのみに使える？

# Aggregate ASV abundances by Season
asgard_depth_type <- aggregate(ASGARD_ASV_filtered, by=list(meta_asgard$depth_type), FUN = mean) #181*78

rownames(asgard_depth_type) = asgard_depth_type$Group.1
asgard_depth_type$Group.1 = NULL

##Find dominant season for each ASV
asgard_pick_max_dt = as.numeric(sapply(asgard_depth_type,which.max)) #ながさ78
asgard_dominant_dt = rownames(asgard_depth_type)[asgard_pick_max_dt] #ながさ78

## Define colors for seasons
asgard_dt_colors <- c("surf" = "orange", "mid" = "blue", "bottom" = "black")
#ながさ78

asgard_vertex_colors_dt <- asgard_dt_colors[asgard_dominant_dt] #ながさ78

# Ensure colors are properly assigned
asgard_vertex_colors_dt[is.na(asgard_vertex_colors_dt)] <- "gray"  # Default for unassigned ASVs

plot(asgard_network, 
     layout = asgard_network_layout,
     vertex.size = 10,  # Scale node size by abundance or degree
     vertex.color = asgard_vertex_colors_dt,  # Color nodes by temp
     vertex.frame.color = NA,
     vertex.label = NA,  # Remove labels for clarity
     main = "Network Graph with Depth Type")

# 78はASVの数, 181はサンプル数
# networkは長さ78, レイアウトは長さ156(78*2, xy軸)

### Example 2: Plot abundance as point size

# take the sqrt of the mean relative abundances of each ASV
asgard_node_size = sqrt(colMeans(ASGARD_ASV_filtered))

plot(asgard_network, 
     layout = asgard_network_layout,
     vertex.size = 100*asgard_node_size,  # Scale node size by abundance or degree
     vertex.color = vertex_colors,  # Color nodes by season
     vertex.frame.color = NA,
     vertex.label = NA,  # Remove labels for clarity
     main = "Network Graph with Dominant Season colors")

### Example 3: Assign colors to ASVs based on heirarchical clustering

# Compute Jaccard distance between ASVs
asgard_jaccard_dist <- vegdist(t(ASGARD_ASV_filtered), method = "jaccard")  # Transpose: ASVs as rows

# Perform hierarchical clustering using Ward's method
asgard_hclust <- hclust(asgard_jaccard_dist, method = "ward.D2")

# Define number of clusters
k <- 8  # Adjust as needed

# Assign cluster memberships
asgard_hclusters <- cutree(asgard_hclust, k = k)

# Generate distinct colors for clusters
hcluster_colors <- rainbow(k)

# Assign colors to ASVs
asgard_color_mapping <- hcluster_colors[asgard_hclusters]

# Ensure ASV names match network vertices
vertex_colors <- asgard_color_mapping

# Ensure colors are properly assigned
vertex_colors[is.na(vertex_colors)] <- "gray"  # Default for unassigned ASVs

plot(asgard_network, 
     layout = asgard_network_layout,
     vertex.size = 100*asgard_node_size,  # Scale node size by abundance or degree
     vertex.color = vertex_colors,  # Color nodes by season
     vertex.frame.color = NA,
     vertex.label = NA,  # Remove labels for clarity
     main = "Network Graph with Heirarchical Clustering colors")

### Example 4: Use tSNE to layout ASVs 

# Run tSNE dimensional reduction to get a nice graph layout
library(Rtsne)
asgard_tsne = Rtsne(t(ASGARD_ASV_filtered)^.25, perplexity=10, dims = 2, check_duplicates = F,
                    pca = T, max_iter = 2500, verbose = T)
asgard_network_layout2 <- as.matrix(asgard_tsne$Y)
### perplexityが長いとerror出た

# Plot network using tSNE-based layout
# t-SNE（t-Distributed Stochastic Neighbor Embedding）は、高次元データを2次元または3次元に可視化するための次元削減手法
plot(asgard_network, 
     layout = asgard_network_layout2,  # Use custom layout
     vertex.size = asgard_node_size*30,  # Scale node size by abundance or degree
     vertex.color = vertex_colors,  # Color nodes by season
     vertex.frame.color = NA,
     edge.width = E(asgard_network)$weight * 0.05, 
     vertex.label = NA,  # Remove labels for clarity
     main = "Network Graph with tSNE Layout")

### Example 5: Plot Seasonal variation in the network

# Plot network
node_size <- sqrt(asv_season["spring",])

plot(network, 
     layout = network_layout,  # Use custom layout
     vertex.size = node_size*30,  # Scale node size by abundance
     vertex.color = vertex_colors,  # Color nodes by cluster
     vertex.frame.color = NA,
     edge.width = E(network)$weight * 0.05, 
     vertex.label = NA,  # Remove labels for clarity
     main = "Spring")

node_size <- sqrt(asv_season["summer",])

plot(network, 
     layout = network_layout,  # Use custom layout
     vertex.size = node_size*30,  # Scale node size by abundance
     vertex.color = vertex_colors,  # Color nodes by cluster
     vertex.frame.color = NA,
     edge.width = E(network)$weight * 0.05, 
     vertex.label = NA,  # Remove labels for clarity
     main = "Summer")


node_size <- sqrt(asv_season["fall",])

plot(network, 
     layout = network_layout,  # Use custom layout
     vertex.size = node_size*30,  # Scale node size by abundance
     vertex.color = vertex_colors,  # Color nodes by cluster
     vertex.frame.color = NA,
     edge.width = E(network)$weight * 0.05, 
     vertex.label = NA,  # Remove labels for clarity
     main = "Fall")

dev.off()

###alpha diversity

# Compute alpha diversity indices
asgard_alpha_diversity <- data.frame(
  Sample = rownames(asgard_seqcount),
  Observed_ASVs = rowSums(asgard_seqcount > 0),  # Count of nonzero ASVs per sample
  Shannon = vegan::diversity(asgard_seqcount, index = "shannon"),  
  Simpson = vegan::diversity(asgard_seqcount, index = "simpson"),  
  Chao1 = estimateR(asgard_seqcount)[1, ]  # First row of estimateR gives Chao1
)
# estimateR()は一行目にChao1, 二行目にACE, 三行目にShannon, 四行目にSimpsonを表示する
# veganとigraphの両方にdiversity()が存在するため、vegan::と明示したほうがいい

# Use Sample column as a key and merge it to metadata
meta_asgard <- meta_asgard %>%
  mutate(Sample = rownames(ASGARD_ASV_filtered_frt))

asgard_alpha_meta <- merge(asgard_alpha_diversity, meta_asgard, by = "Sample", sort=F)

# choose a variable and index
GroupVariable = "DO"  #"NO3(uM)"
DiversityVariable = "Simpson"

ggplot(asgard_alpha_meta, aes(x = .data[[GroupVariable]], y = .data[[DiversityVariable]], fill = .data[[GroupVariable]])) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  
  labs(title = paste(DiversityVariable, "by", GroupVariable),
       x = GroupVariable,
       y = DiversityVariable) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Kruskal Wallis test
# rank-based, nonparametric
# 対応のない3つ以上のグループ間で中央値に差があるかどうかを検定するクラスカル・ウォリス検定(H検定)
# ノンパラメトリック検定であるため、データが正規分布していない場合や、等分散性が仮定できない場合に適している
kruskal.test(asgard_alpha_meta[[DiversityVariable]] ~ asgard_alpha_meta[[GroupVariable]])
# dfは168 (169-1)
# Kruskal-Wallis chi-squaredは大きければ、グループ間の中央値に有位な差がある(ただ有位水準なし)
# p値が有位水準(0.05)より小さければ、グループ間の中央値に有位な差がある
# Simpsonでは、p-value = 0.3458なので、グループ間のalpha diversityに有位な差なし
# sal & Shannonでも0.3905, Chao1も0.417

### make a transformation
# rarefied tableやlog10
asgard_asv_rarefied = rrarefy(asgard_seqcount, sample=min(rowSums(asgard_seqcount)))
# ランダム希釈(random rarefication)は、サンプル間のシーケンシング深度の違いを補正するための手法
# 例えば、サンプルAが10,000リード、サンプルBが5,000リードの場合、両方を5,000リードにランダムに減らす

# Make a rarefied proportional table
asgard_prop_rarefied = prop.table(as.matrix(asgard_asv_rarefied), margin=1)

# Make a log-transformed rarefied proportional table
asgard_logprop_rarefied = log1p(asgard_asv_rarefied) #log10(0) = -Infになるので、全てに1足す

### db-RDA (Distance-Based Redundancy Analysis)
# 環境変数（説明変数）の、生物群集の組成（応答変数）への影響を評価
# 距離行列求めて、主座標分析(PCoA)して、重回帰
# PCoAはサンプル間のユークリッド距離を用いるPCA。PCAは数値行列使う
# PERMANOVAはcategorical, db-RDAはcontinuous

# Beta diversityの行列を作る
# asgard_beta = ASGARD_ASV_filtered_frt
asgard_beta = asgard_filtered^.25
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
asgard_pcoa_df <- data.frame(
  Sample = rownames(ASGARD_ASV_filtered_frt), #行名をとるだけなので、asgard_filteredでも可
  PCoA1_Bray = asgard_pcoa_bray$points[, 1],
  PCoA2_Bray = asgard_pcoa_bray$points[, 2],
  PCoA1_Jaccard = asgard_pcoa_jaccard$points[, 1],
  PCoA2_Jaccard = asgard_pcoa_jaccard$points[, 2],
  PCoA1_Euclidean = asgard_pcoa_euclidean$points[, 1],
  PCoA2_Euclidean = asgard_pcoa_euclidean$points[, 2]
) #181*7

meta_asgard <- meta_asgard %>%
  mutate(Sample = rownames(ASGARD_ASV_filtered_frt)) #181*46

asgard_pcoa_df <- left_join(asgard_pcoa_df, meta_asgard, by = "Sample") #181*52
asgard_pcoa_df$cluster = as.factor(clusnum) #181*53

### boxplotを作成して、hmのクラスターを分けた要因を究明する
pdf(file="asgard_boxplots.pdf")
# asgard_pcoa_df作ってから！
rsc = hue_pal()(nclus)[asgard_pcoa_df$cluster]

plot(asgard_pcoa_df$PCoA1_Bray,asgard_pcoa_df$PCoA2_Bray, col = rsc, pch = 19)

for(var in colnames(asgard_pcoa_df)) {
  gg = ggplot(asgard_pcoa_df, aes(x = cluster, y = .data[[var]])) +
    geom_boxplot(aes(fill=cluster),
                 #outliers = FALSE
                 outlier.shape = NA
                 #outlier.alpha = 0,
                 #outlier.size = 0, 
                 #outlier.colour = NULL,
    )+
    geom_jitter(width=.4,height=0)+
    theme(text = element_text(size = 24))
  
  print(gg)
  
  gp = ggplot(asgard_pcoa_df) + geom_point(aes(x=PCoA1_Bray, y= PCoA2_Bray, col=cluster, size=.data[[var]]))
#  plot(asgard_pcoa_df$PCoA1_Bray,asgard_pcoa_df$PCoA2_Bray, col = rsc, cex=pch = 19)
 print(gp) 
  }
dev.off()



# temp: cluster 5 has relatively high temp (>5). cluster 1 and 2 has low temp (<4)
# salinity: cluster 4 have highest salinity. 5 has lowest average salinity. 1 and 2 are similar
# DO: 2 has highest primary productivity. 4 and 5 has the highest respiration

# 背景透明なboxplot
png("transparent_fleco.png", 
    width = 800, 
    height = 600, 
    bg = "transparent") 

rsc = hue_pal()(nclus)[asgard_pcoa_df$cluster]


  ggplot(asgard_pcoa_df, aes(x = cluster, y = `FlECO-AFL(mg/m^3)`)) +
  geom_boxplot(aes(fill=cluster),
                 #outliers = FALSE
                 outlier.shape = NA
                 #outlier.alpha = 0,
                 #outlier.size = 0, 
                 #outlier.colour = NULL,
    )+
  geom_jitter(width=.4,height=0)+
  theme(text = element_text(size = 24),
          # 背景を透明に設定
          #panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", color = NA),
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "transparent"))

dev.off()


# pcoa_dfをmeta dataとmerge
# mergeのために、Sample列を追加し、これをキー列とする
#meta_asgard <- meta_asgard %>%
  #mutate(Sample = rownames(ASGARD_ASV_filtered_frt))

#asgard_pcoa_df <- merge(asgard_pcoa_df, meta_asgard, by = "Sample", sort = FALSE)
#asgard_pcoa_df <- left_join(asgard_pcoa_df, meta_asgard, by = "Sample")

#181*53

### Parameter selection before db-RDA (db-RDAに用いる変数を決める)
## ggpairs
# ggpairs()は変数*変数のスピアマン相関係数を示す
# 同じ変数同士なら分布
# 異なる変数は相関係数と散布図で示される
# install.packages("GGally")
library(GGally)
asgard_ggpairs <- asgard_pcoa_df %>%
  select(-c("chl depth", "log_name", "bottle", "bottom_depth", "lat", "lon", "N+N (umol/L)"))
# NO3とN+Nは同じだったので削除

ggpairs(asgard_ggpairs,
        columns = 19:33,
        lower = list(continuous = "smooth") # add regression line
)

# poc*pon = 0.976, chl*FLECO = 0.867, Sil*PO4 = 0.817, NO3*PO4 = 0.903, NO2*Sil = 0.850, NO2*NO3 = 0.925
# NH4*NO3 = 0.562, NH4*NO3 = 0.693, NH4*PO4 = 0.73
# poc, chl, NO3を選択すればいいかも

## VIF (Variance Inflation Factor, 分散拡大要因)
#重回帰分析を行なった際に各変数に多重共線性が無いかを調べるための統計量
# VIF は、ある説明変数が他の説明変数とどの程度強く相関しているかを示す指標
# 一般的に、VIF が 5 または 10 以上の場合、多重共線性が強いと判断され、回帰モデルの信頼性が低下する可能性がある
# VIF統計量は一般的にに10以下であれば多重共線性がないとされる。理想値は2以下である。VIF統計量が10を超えた変数がある場合にはモデルからその変数を外してもう一度VIF統計量を計算するなど、モデルを再構成する必要がある。

# install.packages("car")
library(car) # Loading required package: carDataと出るが無視してok

# lm()でmulti linear regression計算
# lm(y ~ x1 + x2 + x3...)
# yには目的変数、xには説明変数(影響を及ぼす変数)
# yにはとりあえず、一番大事そうなsalinityを入れてみる
asgard_lm <- lm(`NO3(uM)` ~ salinity + depth_m + DO + `POC (ug/L)` + `PON (ug/L)`
                + `SPM (ug/L)` + `FlECO-AFL(mg/m^3)` + `chl (ug/l)` + `phaeo (ug/l)`
                + `PO4(uM)` + `Sil(uM)` + temp + `NO2(uM)` + `NH4(uM)`,
                data = asgard_ggpairs)

# 重回帰モデルの統計的有意性をテスト
# salinity (y)に、各説明変数(xたち)が与える影響をp値で表す
summary(asgard_lm)
# モデルは有意: p-value: < 2.2e-16

#vif()に重回帰モデルを入れる
asgard_vif <- vif(asgard_lm)
print(asgard_vif)
# POCが26.54, PONが27.92, PO4が10.86, Silが14.6, NO3が19.2, NO2が12.6
# POCとNO3を選択する。従って、変数は、temp, sal, DO, POC, NO3, chl, depth


# PCoA Bray-Curtis plot
GroupVariable = "salinity"

ggplot(asgard_pcoa_df, aes(x = PCoA1_Bray, y = PCoA2_Bray, color = .data[[GroupVariable]])) +
  geom_point(size = 3) +
  labs(title = "PCoA - Bray-Curtis", x = "PCoA 1", y = "PCoA 2") +
  theme_minimal()+
  scale_color_viridis()

# クラスター毎に色分けした
# asgard_dbrda_mergedを先に作ってから
ggplot(asgard_dbrda_merged, aes(x = PCoA1_Bray, y = PCoA2_Bray, color = clusnum)) +
  geom_point(size = 3) +
  labs(title = "PCoA - Bray-Curtis", x = "PCoA 1", y = "PCoA 2") +
  theme_minimal()
  #scale_color_viridis()


### db-RDA
# 非ユークリッド距離尺度を用いたデータに対して制約付き順序付け（constrained ordination）を実行する手法
# CCAやRDAはユークリッド距離を用いるが、community count dataには適応できない
# 距離行列を計算、行列にPCoA、PCoAで得られたeigenvaluesをRDAに適応
# NAのないデータしか使えない

# asgard_pcoa_dfの各列内の、NA数える
for(var in  c("salinity", "temp", "depth_m", "DO", "POC (ug/L)","chl (ug/l)",
              "PON (ug/L)", "SPM (ug/L)", "FlECO-AFL(mg/m^3)", "chl depth",
              "phaeo (ug/l)", "PO4(uM)", "Sil(uM)", "NO2(uM)", "NH4(uM)",
              "N+N (umol/L)", "NO3(uM)")) { 
  asgard_na_count <- sum(is.na(asgard_pcoa_df[["Sil(uM)"]])) # forループで用いるdf[[変数]]はdf$変数と同じ意味
  print(asgard_na_count)
}
# 現状のNAの数: POC+PON+SPM+ = 64, chl+phaeo = 19, NO3+PO4+NO2+NH4+N+N+Sil = 7, DO+temp+sal+depth = 0
# FlECOはNAなし！

# 必要な列だけを備えたsubset作る
asgard_pcoa_df_sub <- asgard_pcoa_df %>% # 181*13
  select(-c(secondary_number,cellid, side, project, station, station_type
            , date, CTD_number, filter, delO18, sedchla_mg_m2, mass, cast_depth, 
            depth_original, `depth_m`, `depth_type`, `bottom_depth`, lat, lon, `log_name`, `bottle`, 
            `niskin_id_meta`, `volume_filtered_l`, `transmission_pct`, `niskin_id_nutrients`,
            `poc_pon_ratio`, `poc_pon_ratio_redfield`, `notes`, `notes2`,
            `POC (ug/L)`, `chl (ug/l)`, `PON (ug/L)`, `SPM (ug/L)`, `chl depth`,
            `phaeo (ug/l)`, `PO4(uM)`, `Sil(uM)`, `NO2(uM)`, `NH4(uM)`, `N+N (umol/L)`))

# db-RDAやる前に、モデルに扱うパラメーター(PCoAのデータ)を標準化しないといけない
# scale()で標準化する。なおPCoAの値などには手を加えない(sub1)
# 1/4 transformationは、sequencing dataの標準化
asgard_pcoa_df_sub1 = asgard_pcoa_df %>% 
  select(c("Sample", "PCoA1_Bray","PCoA2_Bray","PCoA1_Jaccard",
           "PCoA2_Jaccard","PCoA1_Euclidean","PCoA2_Euclidean"))

asgard_pcoa_df_sub2 = asgard_pcoa_df %>% 
  select(c("temp", "salinity","DO","NO3(uM)", "FlECO-AFL(mg/m^3)")) %>% 
  scale() 

asgard_pcoa_df_sub = cbind(asgard_pcoa_df_sub1, asgard_pcoa_df_sub2)
#181*12

# can only use data with no NAs
# 全てがNAの列を削除
### POC等はNAを64含む。よって行を消すと、64サンプルを失う(64/181なので、およそ1/3のサンプルを失う)
### chl-aはNAを19もち、かつtemp, sal, depth, DO, NO3を含んだモデルで*だったので削除した

### ここが大問題！！
asgard_complete_cases = complete.cases(asgard_pcoa_df_sub) #NAあればFALSE, なければTRUEでlogi vecを返す
asgard_pcoa_cc = asgard_pcoa_df_sub[asgard_complete_cases,] #174*12, TRUEの行だけ抽出(181行から174行にまで減る)
asgard_bray_cc = vegdist(ASGARD_ASV_filtered_frt[asgard_complete_cases,],method="bray") #asv_frtpropは4乗根でtransformedしたprop table

# chl (ug/l)という名前だと、スペースあるせいでエラー
# clean_names()使うか、列名を変更するか、バッククォート(escの下)で囲む
colnames(asgard_pcoa_cc)[colnames(asgard_pcoa_cc) == "chl (ug/l)"] <- "chl_ug_l"

asgard_dbrda_model <- capscale(asgard_bray_cc ~ temp + salinity 
                               + DO + `NO3(uM)` + `FlECO-AFL(mg/m^3)`,
                               data = asgard_pcoa_cc) # `FlECO-AFL(mg/m^3)`
# `N+N (umol/L)`はNO3と全く同じデータ=重複してるため(共線性高い)、除外される

# 結果のプリント
summary(asgard_dbrda_model)

# Test significance of the db-RDA model (999 permutations)
# anova.cca()で、CCA(Canonical Correspondence Analysis, 正準対応分析)やRDAを順列検定(follow-up test)
# CCAとRDAはともに、多変量データと環境要因の相関を見るが、CCAは非線形、RDAは線形
asgard_anova_dbrda <- anova.cca(asgard_dbrda_model, permutations = 999)
print(asgard_anova_dbrda) #0.001***

###大事!! Test significance of each environmental variable
# by = "margin"で各説明変数の寄与を評価するために、順列検定を実行
asgard_anova_env <- anova.cca(asgard_dbrda_model, by = "margin", permutations = 999)
print(asgard_anova_env) # temp, sal, NO3のp-valueは0.001 ***, DOは**, flECOは*

# Test significance of db-RDA axes
# by = "axis"で各軸に対して個別に順列検定を実行
# axesなので、CAP1-4が出る
# CAP = Canonical Axis Principal
asgard_anova_axes <- anova.cca(asgard_dbrda_model, by = "axis", permutations = 999)
print(asgard_anova_axes)

# Adjusted R² value (accounts for number of predictors)
RsquareAdj(asgard_dbrda_model) 

# Plot results (PCAみたい）
# Extract db-RDA site scores (sample coordinates)
asgard_dbrda_scores <- as.data.frame(scores(asgard_dbrda_model, display = "sites"))
# scores()は、多変量解析の結果からスコア（座標）を抽出するために使用
asgard_dbrda_scores$MDS1 = asgard_dbrda_model$CA$u[,1]
asgard_dbrda_scores$MDS2 = asgard_dbrda_model$CA$u[,2]

asgard_dbrda_scores$Sample = rownames(asgard_dbrda_scores)

# Extract sample scores (for points)
asgard_dbrda_vectors <- as.data.frame(scores(asgard_dbrda_model, display = "bp"))  # "bp" = biplot arrows
asgard_dbrda_vectors$Variable <- rownames(asgard_dbrda_vectors)  # 列追加, Add variable names

asgard_dbrda_df = merge(asgard_dbrda_scores, asgard_pcoa_cc, by="Sample", sort=FALSE)

# add cluster column

asgard_dbrda_merged = left_join(asgard_dbrda_df, rownames_to_column(as.data.frame(clusnum)), by = c("Sample" = "rowname"))
asgard_dbrda_merged$clusnum = as.factor(asgard_dbrda_merged$clusnum)

GroupVariable = "salinity" 

# Plot db-RDA in ggplot2
ggplot() +
  # Plot points, colored by group
  geom_point(data = asgard_dbrda_merged, 
             aes(x = CAP1, y = CAP2, color = clusnum, size = `NO3(uM)`)
             ) + 
#size = `NO3(uM)`みたいにdotのサイズ指定可
#.data[[GroupVariable]]
# xyを、CAP1と2ではなく、MDS1と2に変更可能
  
   # Add environmental arrows
  geom_segment(data = asgard_dbrda_vectors, 
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black", size = 1) +
  
  # Add environmental variable labels near arrow tips
  geom_text(data = asgard_dbrda_vectors, 
            aes(x = CAP1, y = CAP2, label = Variable), 
            vjust = -0.5, hjust = 0.5, size = 5) +
  
  # Labels and theme
  labs(x = "CAP1", y = "CAP2", color = "cluster") +
  theme_minimal()
# title = "db-RDA Ordination with Environmental Vectors"

### dbRDAのJPEG
jpeg(filename = "ASGARD_dbrda_mds.jpeg",
     width = 8,                         # 幅（インチ）
     height = 6,                        # 高さ（インチ）
     units = "in",                      # 単位（インチ）
     res = 2000,                        # 解像度（dpi）
     quality = 100                      # 画質（0-100）
)

ggplot() +
  # Plot points, colored by group
  geom_point(data = asgard_dbrda_merged, 
             aes(x = MDS1, y = MDS2, color = clusnum, size = `NO3(uM)`)
  ) + 
  #size = `NO3(uM)`みたいにdotのサイズ指定可
  #.data[[GroupVariable]]
  # xyを、CAP1と2ではなく、MDS1と2に変更可能
  
  # Add environmental arrows
  geom_segment(data = asgard_dbrda_vectors, 
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black", size = 1) +
  
  # Add environmental variable labels near arrow tips
  geom_text(data = asgard_dbrda_vectors, 
            aes(x = CAP1, y = CAP2, label = Variable), 
            vjust = -0.5, hjust = 0.5, size = 5) +
  
  # Labels and theme
  labs(x = "MDS1", y = "MDS2", color = "cluster") +
  theme_minimal()

dev.off()

# db-rdaの透明png
png("transparent_dbrda.png", 
    width = 800, 
    height = 600, 
    bg = "transparent") 

ggplot() +
  # Plot points, colored by group
  geom_point(data = asgard_dbrda_merged, 
             aes(x = MDS1, y = MDS2, color = clusnum, size = `NO3(uM)`)
  ) + 
  #size = `NO3(uM)`みたいにdotのサイズ指定可
  #.data[[GroupVariable]]
  # xyを、CAP1と2ではなく、MDS1と2に変更可能
  
  # Add environmental arrows
  geom_segment(data = asgard_dbrda_vectors, 
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black", size = 1) +
  
  # Add environmental variable labels near arrow tips
  geom_text(data = asgard_dbrda_vectors, 
            aes(x = CAP1, y = CAP2, label = Variable), 
            vjust = -0.5, hjust = 0.5, size = 5) +
  
  # Labels and theme
  labs(x = "MDS1 (19.21 %)", y = "MDS2 (9.21 %)", color = "cluster") +
  theme(text = element_text(size = 30),
        # 背景を透明に設定
        #panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"))


dev.off()

### ggplot maps
# Stadia mapsからAPI得る
# API Key = 85bfb0bd-23a7-437e-9e98-1f1fd9e283f7
# new API: 6430985f-a969-4816-81f3-bafee5bd9da9
#install.packages(("ggmap"))
library(ggmap)
register_stadiamaps("6430985f-a969-4816-81f3-bafee5bd9da9")

bbox <- make_bbox(lon = meta_asgard$lon, lat = meta_asgard$lat, f = 0.1)
mapz = get_stadiamap(bbox, maptype = "stamen_terrain", zoom =4)

asgard_ggmap = left_join(asgard_pcoa_df,
                         rownames_to_column(as.data.frame(meta_asgard)),
                         by = c("Sample" = "rowname"))

# cluster列を追加する
asgard_ggmap = ASGARD_full
asgard_ggmap$cluster = as.factor(clusnum)

### cluster付きdfを別名にする！
asgard_filtered_cluster = asgard_ggmap


# surf, mid, bottomの順に定義
asgard_ggmap$depth_type <- factor(
  asgard_ggmap$depth_type,
  levels = c("surf", "mid", "bottom")  
)

# rowがdepth
map_plot_1 <- ggmap(mapz) +
  geom_point(data = asgard_ggmap, aes(x = lon, y = lat, color = cluster),
             size = 3, alpha = 1) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1)
        )+
  facet_grid(rows = . ~ depth_type)+
  labs(color = "cluster")

print(map_plot_1)

# rowがcluster, columnがdepth
map_plot_2 <- ggmap(mapz) +
  geom_point(data = asgard_ggmap, aes(x = lon, y = lat, color = cluster),
             size = 3, alpha = 1) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  facet_grid(depth_type ~ cluster)+
  labs(color = "cluster")

print(map_plot_2)

#透明ggmap
png("transparent_map.png", 
    width = 800, 
    height = 600, 
    bg = "transparent")

map_plot_2 <- ggmap(mapz) +
  geom_point(data = asgard_ggmap, aes(x = lon, y = lat, color = cluster),
             size = 3, alpha = 1) +
  facet_grid(depth_type ~ cluster)+
  labs(color = "cluster")+
  theme(text = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1),
        #panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"))
 

print(map_plot_2)

dev.off()

### うまくできない箇所！！！！！
#Visualizing Environmental Gradients on PCoA
#Error in text.default(asgard_pcoa_bray_cc$points, labels = rownames(asgard_pcoa_bray_cc$points),  : 
#plot.new has not been called yet というエラー出る
asgard_pcoa_bray_cc <- cmdscale(asgard_bray_cc, k = 2, eig = TRUE)
asgard_env_fit <- envfit(asgard_pcoa_bray_cc, asgard_pcoa_cc[, c("temp", "salinity", "depth_m","DO", "chl_ug_l", "NO3(uM)", "POC (ug/L)")], perm = 999, na.rm=T)
plot(asgard_env_fit, col = "red")

### Which ASVs are significantly correlated with environmental parameters?
# ASVと環境パラメーターの相関！

# fullnamebootは全クルーズのデータ含むので、ASGARDだけ抽出
ASGARD_taxa <- fullnameboot[colnames(seqtab_prop_asgard), ] #100*24

#pcoa_ccは、db-RDAの際に、NA含むサンプルを除外している(64から34に減った)
asgard_filtered_cc = asgard_frtprop[asgard_pcoa_cc$Sample,]
asgard_filtered_cc = asgard_filtered_cc[,colSums(asgard_filtered_cc)>0]

asgard_env_asvs = NULL
for(envvar in  c("salinity", "temp", "depth_m", "DO", "NO3(uM)")) { #envvarに環境パラメーターを代入
  print(envvar)
  asgard_cor_results <- apply(asgard_filtered_cc, 2, function(asv) cor(asv, asgard_pcoa_cc[,envvar], method = "spearman"))
  #pcoa_ccのenvvarと一致する列を抽出
  # apply()は行列やデータフレームの行や列に対して、指定した関数を適用
  # 1で行、2で列に関数を適用
  # よって、asv_filtered_ccの2 = 各列(リード配列)をfunction(asv)に適用する=各列を渡す
  asgard_significant_asvs <- names(which(abs(asgard_cor_results) > 0.3))  # Adjust threshold
  # 閾値は個別にsignificant_asvsを実行して、値を一通り見てから設定する！
  # abs()は相関係数の絶対値を計算
  # cor_resultsのabs()が閾値より大きいものの位置を、which()でindicesで返す
  # cor_resultsはnamed vectorなので、namesでindicesに相当するリード配列を抽出してくれる
  asgard_env_asvs = rbind(asgard_env_asvs,cbind(ASGARD_taxa[asgard_significant_asvs,],envvar))
  # taxaから、significant_asvsを抽出
  # 上をenvvarをcbind = 列を結合(横方向)
  # それと、env_asvsをrbind = 行を結合(下方向)
  # envvarには、ASVとの相関が、閾値を超える複数の環境パラメーターがすべて記録される

}
View(asgard_env_asvs)
# 0.3以上入れると、行数が一致しないエラー！
# うまく読み込めてない
print(asgard_cor_results) #ASVと環境パラメーターのスピアマン相関係数が表示される！
hist(asgard_cor_results) #ヒストグラムでチェック


#library(pheatmap)
pdf(file="asgard_envhm.pdf",width=12,height=12)
env_taxa_correlation <- cor(asgard_filtered_cc, asgard_pcoa_cc[, c("temp", "salinity","DO")], method = "spearman")
pheatmap(env_taxa_correlation, main = "Correlation of ASVs with Environmental Variables", cex=0.1)
dev.off()
# pcoa_ccから[, c()]で4列を抽出

### PERMANOVAとPERMDISP

# PERMANOVAは、群間差を検定するための統計手法
# ANOVAは1次元データ、こいつは多次元(多変量)データで、距離行列に基づく
# 距離行列~グループ変数を渡す
# permutations = 順列検定の回数(デフォは999回)

GroupVariable = "temp" #NO3(uM)
# asgard_pcoa_dfからNA除去する
# asgard_pcoa_df[is.na(asgard_pcoa_df)] <- 0 

asgard_permanova_bray <- adonis2(asgard_braymat ~ asgard_pcoa_df[[GroupVariable]], permutations = 999)
asgard_permanova_jaccard <- adonis2(asgard_jacmat ~ asgard_pcoa_df[[GroupVariable]], permutations = 999)
asgard_permanova_euclidean <- adonis2(asgard_eucmat ~ asgard_pcoa_df[[GroupVariable]], permutations = 999)

# temp, sal, DOはNAなし、一方NO3はNAを7つ含むので、NO3には別のobject作る
asgard_pcoa_df_no3 <- asgard_pcoa_df %>%
  filter(!is.na(`NO3(uM)`)) #174*52

# 距離行列も行と列が同じ長さでないといけない
asgard_no3_vec <- !is.na(asgard_pcoa_df$`NO3(uM)`) #長さ181, !is.naでNA含まない行をTRUEにする

asgard_beta_no3 <- asgard_beta[asgard_no3_vec, ] #174*78

asgard_eucmat_no3 = vegdist(asgard_beta_no3, method="euclidean")
asgard_braymat_no3 = vegdist(asgard_beta_no3, method="bray")
asgard_jacmat_no3 = vegdist(asgard_beta_no3, method="jaccard")

asgard_permanova_bray <- adonis2(asgard_braymat_no3 ~ asgard_pcoa_df_no3[[GroupVariable]], permutations = 999)
asgard_permanova_jaccard <- adonis2(asgard_jacmat_no3 ~ asgard_pcoa_df_no3[[GroupVariable]], permutations = 999)
asgard_permanova_euclidean <- adonis2(asgard_eucmat_no3 ~ asgard_pcoa_df_no3[[GroupVariable]], permutations = 999)

# salinity, temp, NO3, DOは全てで、0.001***
print(asgard_permanova_bray)
print(asgard_permanova_jaccard)
print(asgard_permanova_euclidean)

# PERMANOVA assumes similar dispersion (variance) = 分散 between groups.
# We test this with PERMDISP (Beta Dispersion Test).
# 多変量データのグループ内の分散が等しいかどうかを検定
# PERMANOVA がグループ間の差異を検定するのに対し、PERMDISP はグループ内の分散を比較する

# Test dispersion homogeneity (PERMDISP = Permutational Analysis of Multivariate Dispersions)
asgard_permdisp_bray <- betadisper(asgard_braymat, asgard_pcoa_df[[GroupVariable]])
# salinityで上をやるとこのエラー発生: Warning in betadisper(asgard_braymat, asgard_pcoa_df[[GroupVariable]]) :some squared distances are negative and changed to zero
# PO4だと、some squared distances are negative and changed to zero
asgard_permdisp_jaccard <- betadisper(asgard_jacmat, asgard_pcoa_df[[GroupVariable]])
asgard_permdisp_euclidean <- betadisper(asgard_eucmat, asgard_pcoa_df[[GroupVariable]])
# betadisper()はadonis2()と文法がほぼ同じ

# NO3用
GroupVariable = "temp" #NO3(uM)
asgard_permdisp_bray <- betadisper(asgard_braymat_no3, asgard_pcoa_df_no3[[GroupVariable]])
asgard_permdisp_jaccard <- betadisper(asgard_jacmat_no3, asgard_pcoa_df_no3[[GroupVariable]])
asgard_permdisp_euclidean <- betadisper(asgard_eucmat_no3, asgard_pcoa_df_no3[[GroupVariable]])

# Run statistical test
anova(asgard_permdisp_bray) # salは0.03879*, tempは6.391e-16 ***, NO3は1.344e-12 ***, DOは8.154e-16 ***
anova(asgard_permdisp_jaccard) # salは0.0148*, tempは8.667e-15 ***, NO3は2.2e-16 ***,  DOは2.2e-16 ***
anova(asgard_permdisp_euclidean) # 0.0005124 ***, tempは2.2e-16 ***、NO3は2.2e-16 ***,  DOは1.378e-15 ***
# PERMDISP の順列検定と ANOVA の結果を比較することで、結果の頑健性を確認し、結果が合致すれば信頼性高まる
# tempで実行すると、Warning in anova.lm(lm(Distances ~ Groups, data = model.dat)) :ANOVA F-tests on an essentially perfect fit are unreliableとエラー出る

### Mantel test for continuous variables (like temperature, salinity)
# マンテルテストは、2つの距離行列の間の相関を評価
# continuous variableは、euclidian matrixにする
# ASVはbrays

# Compute Euclidean distance matrix for a continuous variable (e.g., pH)
asgard_salinity_dist <- dist(asgard_pcoa_df$salinity) #デフォルトでeuclidian
asgard_temp_dist <- dist(asgard_pcoa_df$temp) #デフォルトでeuclidian
asgard_do_dist <- dist(asgard_pcoa_df$DO) #デフォルトでeuclidian
asgard_no3_dist <- dist(asgard_pcoa_df_no3$`NO3(uM)`) #デフォルトでeuclidian

# Perform Mantel test
asgard_mantel_result_sal <- mantel(asgard_braymat, asgard_salinity_dist, permutations = 999)
asgard_mantel_result_temp <- mantel(asgard_braymat, asgard_temp_dist, permutations = 999)
asgard_mantel_result_do <- mantel(asgard_braymat, asgard_do_dist, permutations = 999)
asgard_mantel_result_no3 <- mantel(asgard_braymat_no3, asgard_no3_dist, permutations = 999)

# Print result
print(asgard_mantel_result_sal) #0.3902, pvalue is 0.001***, 正の相関あり
print(asgard_mantel_result_temp) #0.3026, pvalue is 0.001***, 正の相関あり
print(asgard_mantel_result_do) #0.1123, pvalue is 0.002, 正の相関あり
print(asgard_mantel_result_no3) #0.1006, pvalue is 0.013, 正の相関あり

# Upper quantiles of permutations (null model)は、順列検定の上位分位数（90%、95%、97.5%、99%）が表示
# Mantel統計量 0.3902 は、これらの分位数よりも大きいため、統計的に有意

### ASGARDからは、ASVsは幾つ見つかった？
### サンプルの数は幾つ？

# 各クルーズがどれだけASVを含んでいるか知りたい
b = seqtab_filt
bprop = prop.table(b, 1)

dim(bprop) #2193*18535

colnames(bprop) = shorternames
view(bprop)

bprop_merged = merge(meta_denovo_2, bprop, by = 0, all.x = FALSE)
view(bprop_merged)

dim(meta_denovo_2) #1163*45

dim(bprop_merged) #1158*18581

bprop_asgard = bprop_merged %>%
  filter(project == "ASGARD2017")

view(bprop_asgard)
dim(bprop_asgard) #290*18581

#bprop_asgardのメタデータの列は45, Row.namesが1列
#よって、ASV列は18581 (bprop_mergedの列数)から46を引いた、18535

# 本当に正しいかチェックする
bprop_check = bprop_asgard %>%
  select(-c(1:46))
view(bprop_check) 
dim(bprop_check) #290*18535

bprop_counts = colSums(bprop_check) #vector
length(bprop_counts) #18535
bprop_0 = bprop_counts[bprop_counts == 0]
view(bprop_0) #14575

## サンプル数
a <- meta_denovo_2 %>%
  filter(project == "ASGARD2017") #290*45

#つまり、ASGARDは290サンプル、そこから、18535-14575 = 3960 ASVsが見つかった


### クラスター毎のscatter plot
  geom_point(stat = "identity")


ggplot(asgard_pcoa_df, aes(x = salinity, y = temp, colour = cluster, size=`NO3(uM)`)) +
  geom_point()


### クラスター毎のASVのリストアップ
names(colclusnum[colclusnum==8])


# waffle chart in terms of abundance (class)
#taxaout = readRDS("/Users/shibuerei/Desktop/taxout_edit.rds")

install.packages("waffle")
x=cbind(fullnamemat[which(shorternames %in% names(colclusnum)),"Phylum"],colclusnum,colSums(asgard_filtered))
x = as.data.frame(x)
x$V3 = as.numeric(x$V3)

#df = x %>% count(V1, colclusnum)

library(waffle)


df2 = x %>% 
  group_by(V1, colclusnum) %>%
  summarize(weight = sum(V3), .groups="drop")

df2 = df2 %>% mutate(units = round(weight*100))

total_units = sum(df2$units)

df3 = df2 %>% mutate(units = round(weight / sum(weight) * 100))

df4 = df2 %>% group_by(colclusnum) %>%
  mutate(units = round(weight / sum(weight) * 100)) %>%
  ungroup()


# 
# gg = ggplot(df2, aes(fill=V1, values=units)) +
#   geom_waffle(n_rows = 10, size = 0.00, colour = NA) +
#   facet_wrap(~colclusnum)
# 
# 
# gg2 = ggplot(df2, aes(fill=V1, values=units)) +
#   geom_waffle(n_rows = ceiling(sqrt(total_units)), size = 0.00, colour = NA) +
#   facet_wrap(~colclusnum)
# 
# gg3 = ggplot(df2, aes(fill=V1, values=units)) +
#   geom_waffle(n_rows = ceiling(sqrt(total_units)), size = 0.00, colour = NA) +
#   facet_wrap(~colclusnum) +
#   coord_equal()


# gg4 = ggplot(df3, aes(fill=V1, values=units)) +
#   geom_waffle(n_rows = 10, size = 0.00, colour = NA) +
#   facet_wrap(~colclusnum)


gg5 = ggplot(df4, aes(fill=V1, values=units)) +
  geom_waffle(n_rows = 10, size = 0.00, colour = NA) +
  facet_wrap(~colclusnum) + coord_equal() +
  scale_fill_viridis_d(option = "turbo")+
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1)
  )
#+ scale_fill_brewer(palette = "Set1")


pdf(file="waffle.pdf", width=24, height=24)
# gg
# gg2
# gg3
# gg4
gg5
dev.off()

#waffleの透明png
png("transparent_waffle.png", 
    width = 800, 
    height = 600, 
    bg = "transparent") 
gg5 = ggplot(df4, aes(fill=V1, values=units)) +
  geom_waffle(n_rows = 10, size = 0.00, colour = NA) +
  facet_wrap(~colclusnum) + coord_equal() +
  scale_fill_viridis_d(option = "turbo")+
  labs(fill="")+
  theme(text = element_text(size = 30),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(1, "cm"),
      # 背景を透明に設定
      #panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", color = NA),
      #panel.grid.major = element_blank(),
      #panel.grid.minor = element_blank(),
      legend.background = element_rect(fill = "transparent"),
      legend.box.background = element_rect(fill = "transparent"))
gg5
dev.off()

# DBO3を各クラスターが含むか調べる
# フィルタリング
dbo3check <- asgard_pcoa_df[grep("DBO3", asgard_pcoa_df$station), ]

# クラスターごとのカウント
dbo3count <- aggregate(station ~ cluster, data = dbo3check, FUN = length)
colnames(dbo3count) <- c("cluster", "count")
print(dbo3count)
# クロス集計表作る = クラスターが幾つのサンプルを持っているか調べる
cluster_xtabs <- xtabs(~ cluster, data = asgard_pcoa_df)
print(cluster_xtabs)

# 各クラスターが、surf, mid, depthをいくつ持ってるか？
# 数式形式での集計
cluster_depth_xtabs <- xtabs(~ cluster + depth_type, data = asgard_pcoa_df)
print(cluster_depth_xtabs)

#サンプル数把握
dim(seqtab_mat) #2193*3076
dim(meta_denovo_2) #1163*45

meta_asgard_new <-meta_denovo_2 %>%
  filter(project == "ASGARD2017")

dim(meta_asgard_new) #290*45
# ASGARDは290サンプル

##次にASVの数を調べる
new290 <- meta_asgard_new$cellid
length(new290)

asgard_ASVnums <- seqtab_mat[new290,]
dim(asgard_ASVnums)

# 列数は変化してない
# 290行のデータの中で、すべての行が0の列を削除すれば、ASGARDで最低でも1は見つかったASVだけ抽出
asgard_ASVnums_filtered <- asgard_ASVnums[, colSums(asgard_ASVnums) > 0]
dim(asgard_ASVnums_filtered)
# 290 534