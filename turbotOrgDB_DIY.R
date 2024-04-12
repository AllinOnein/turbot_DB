#install.packages("")
library(openxlsx)
library(dplyr)
library(stringr)
library(jsonlite)
library(AnnotationForge)
#顺手设置一下options
options(stringsAsFactors = F)

#D读取注释文件
emapper <- read.xlsx("D:/BaiduSyncdisk/RNA-seq-dataDB-DIY/turbot2021CWJ.annotations2.xlsx", sep = "\t")
#将空值替换为NA，方便后续使用na.omit()函数提出没有注释到的行
emapper[emapper=="-"]<-NA

gene_info <- emapper %>% dplyr::select(GID = query, GENENAME = Preferred_name) %>% na.omit()
gos <- emapper %>% dplyr::select(query, GOs) %>% na.omit()

gene2go = data.frame(GID = character(),
                     GO = character(),
                     EVIDENCE = character())

gos_list <- function(x){
  the_gos <- str_split(x[2], ",", simplify = FALSE)[[1]]
  df_temp <- data.frame(GID = rep(x[1], length(the_gos)),
                        GO = the_gos,
                        EVIDENCE = rep("IEA", length(the_gos)))
  return(df_temp)
}
gene2gol <- apply(as.matrix(gos),1,gos_list)
gene2gol_df <- do.call(rbind.data.frame, gene2gol)
gene2go <- gene2gol_df
gene2go$GO[gene2go$GO=="-"]<-NA
gene2go<-na.omit(gene2go)


gene2ko <- emapper %>% dplyr::select(GID = query, Ko = KEGG_ko)
gene2ko$Ko[gene2ko$Ko=="-"]<-NA
gene2ko <- na.omit(gene2ko)
gene2kol <- apply(as.matrix(gene2ko),1,gos_list)
gene2kol_df <- do.call(rbind.data.frame, gene2kol)
gene2ko <- gene2kol_df[,1:2]
colnames(gene2ko) <- c("GID","Ko")
gene2ko$Ko <- gsub("ko:","",gene2ko$Ko)


# library(jsonlite)
# 下面的json = "ko00001.json"，如果你下载到其他地方，记得加上路径
update_kegg <- function(json = "D:/Desktop/turbotDB_DIY/ko00001.json") {
  pathway2name <- tibble(Pathway = character(), Name = character())
  ko2pathway <- tibble(Ko = character(), Pathway = character())
  kegg <- fromJSON(json)
  for (a in seq_along(kegg[["children"]][["children"]])) {
    A <- kegg[["children"]][["name"]][[a]]
    for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
      B <- kegg[["children"]][["children"]][[a]][["name"]][[b]] 
      for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
        pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
        pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
        pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
        pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
        kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
        kos <- str_match(kos_info, "K[0-9]*")[,1]
        ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))}}}
  save(pathway2name, ko2pathway, file = "kegg_info.RData")}
# 调用函数后在本地创建kegg_info.RData文件，以后只需要载入 "kegg_info.RData"即可
update_kegg()
# 载入kegg_info.RData文件
load(file = "kegg_info.RData")

gene2pathway <- gene2ko %>% left_join(ko2pathway, by = "Ko") %>% dplyr::select(GID, Pathway) %>% na.omit()


#可以去重，以防万一，或者等下面报错提示去重再去跑代码也行。原理是，makeOrgPackage不允许有重复的行，因此需要删除，除了gene2pathway可能有重复的行，其他几乎不可能有重复的情况。
gene2go <- unique(gene2go)
gene2go <- gene2go[!duplicated(gene2go),]
gene2ko <- gene2ko[!duplicated(gene2ko),]
gene2pathway <- gene2pathway[!duplicated(gene2pathway),]
gene_info <- gene_info[!duplicated(gene_info),]


makeOrgPackage(gene_info=gene_info,
               go=gene2go,
               ko=gene2ko,
               pathway=gene2pathway,
               version="1.1",  #版本，使用？makeOrgPackage，拉到最下面查看
               maintainer = "sunfei <sunfei_1997@163.com>",  #修改为你的名字和邮箱
               author = "sunfei <sunfei_1997@163.com>",  #修改为你的名字和邮箱
               outputDir = "D:/Desktop/turbotDB_DIY",  #输出文件位置
               tax_id = "52904",  #你在NCBI上查并定义的id
               genus = "Scophthalmus",  #你在NCBI上查并定义的属名
               species = "maximus",  #你在NCBI上查并定义的种名
               goTable="go")


install.packages("D:/Desktop/turbotDB_DIY/org.Smaximus.eg.db", repos=NULL, type="sources")
library(org.Smaximus.eg.db)
# 查看所有列的信息
columns(org.Smaximus.eg.db)
# 查看所有基因
keys(org.Smaximus.eg.db)
# 查看特定基因的信息
# library(dplyr)
select(org.Smaximus.eg.db, keys = "ENSSMAT00000016768", columns = c("GO"))


###因为我们在后面进行kegg分析的时候需要背景文件，因此为了方便即可导出pathway2name和pathway2gene
# 做一些常规格式化
pathway2name$Name <- gsub(" \\[BR:ko[0-9]{5}\\]", "",pathway2name$Name)
pathway2name <- na.omit(pathway2name)
pathway2gene <- gene2pathway[, c("Pathway","GID")]
# 输出
write.table(pathway2name,file = "D:/Desktop/turbotDB_DIY/pathway2name", sep = "\t", quote = F, row.names = F)
write.table(pathway2gene,file = "D:/Desktop/turbotDB_DIY/pathway2gene", sep = "\t", quote = F, row.names = F)

library(ClusterProfiler)
pathway2gene <- read.table("D:/Desktop/turbotDB_DIY/pathway2gene",header = T,sep = "\t")
pathway2name <- read.table("D:/Desktop/turbotDB_DIY/pathway2name",header = T,sep = "\t")
geneIDMap <- bitr(pathway2gene$GID, fromType="GID", toType="GENENAME", OrgDb="org.Smaximus.eg.db")
pathway2gene2 <- merge(pathway2gene, geneIDMap, by = "GID")
pathway2gene2 <- pathway2gene2 [,-1]
write.table(pathway2gene2,file = "D:/Desktop/turbotDB_DIY/pathway2gene2", sep = "\t", quote = F, row.names = F)
