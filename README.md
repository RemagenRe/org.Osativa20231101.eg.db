水稻GO富集分析-R语言自建orgdb库

Rice GO enrichment analysis

数据来源：

https://rapdb.dna.affrc.go.jp/download/irgsp1.html

http://eggnog-mapper.embl.de/

已建好的水稻orgdb库(可直接下载导入)：
https://wwi.lanzoup.com/iMVhB1dg7uyh
密码:e5ia

①下载org.Osativa20231101.eg.db.7z文件并解压，
假设解压到路径"D:/workdata/test"下，
通过运行以下R语言代码导入orgdb：
```
rm(list=ls())
setwd("D:/workdata/test")
install.packages("org.Osativa20231101.eg.db",repos=NULL,type="source")
```
②使用Transcript ID进行R语言GO分析(例如：Os08t0560900-01,Os08t0560900-02,Os08t0560900-03)

附上水稻GID转换网站：https://rapdb.dna.affrc.go.jp/converter/

批量处理含Transcript ID文件的R语言代码：

```
rm(list=ls())
setwd("D:/workdata/test") #此脚本批量处理文件夹内所有*_TranscriptID文件
library(org.Osativa20231101.eg.db)
library(clusterProfiler)
library(ggplot2)

#id_txt="*_TranscriptID"

riceGO <- function(id_txt) {
  gene_ids <- readLines(id_txt)
  
  ego <- enrichGO(gene_ids,
                  OrgDb = org.Osativa20231101.eg.db,
                  keyType = 'GID',
                  ont = 'ALL',
                  qvalueCutoff = 0.05,
                  pvalueCutoff = 0.05,
                  pool = TRUE)
  
  df_ego <- as.data.frame(ego)
  write.csv(df_ego, paste0("GO of ", id_txt,".csv"), row.names = FALSE)
  
  dotp <- enrichplot::dotplot(ego,font.size =8,split = 'ONTOLOGY')+
    facet_grid(ONTOLOGY~., scale="free")+     
    theme(legend.key.size = unit(10, "pt"),#调整图例大小
          plot.margin=unit(c(1,1,1,1),'lines'))#调整四周留白大小
  ggsave(dotp,filename = paste0("riceGO_dotplot_", id_txt, ".png"),width =10,height =14)
}

my_list <- list.files(pattern = "\\_TranscriptID$")

lapply(my_list, function(file) {
  tryCatch({
    riceGO(file)
  }, error=function(err) {
    cat(paste("Error occurred while processing file:", file, "\n"))
    cat("Error message:", conditionMessage(err), "\n")
  })
})
```

该脚本运行后输出图表形式和表格形式的GO结果

示例的输入文件test.example_TranscriptID：

![image](https://github.com/RemagenRe/org.Osativa20231101.eg.db/assets/114082077/1082cfcd-e497-449c-9089-db8994f25be8)


示例输出：

![image](https://github.com/RemagenRe/org.Osativa20231101.eg.db/assets/114082077/23e363cc-feb4-4328-9b7e-12be420c41ff)

![image](https://github.com/RemagenRe/org.Osativa20231101.eg.db/assets/114082077/c9ea0b9e-b795-46e9-8b95-eae806a6aecf)


同时发布在：https://mp.weixin.qq.com/s/Pp7eUPY00of7A11SULvs1Q


DOI: 10.5281/zenodo.11483546 
