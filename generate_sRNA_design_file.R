# Generate small RNA design file 

library(stringr)

# read small RNA file
md <- read.csv("~/rimod/smallRNA/Files_Names_14092016_Tenzin_Javi.csv",  stringsAsFactors = F, header = F)
colnames(md) <- c("sample", "group", "uk", "note", "bamfile", "file")
md <- md[!md$file == "",]

# read RiMod file
rimod <- read.csv("~/rimod/files/FTD_Brain.csv", stringsAsFactors = F)
rimod <- rimod[rimod$REGION == "frontal",]
rimod$SAMPLEID <- str_pad(rimod$SAMPLEID, 5, side="left", pad="0")


 # store samples not in rimod design sheet
only.mirna <- md[!md$sample %in% rimod$SAMPLEID,]

# match both sheets
md <- md[md$sample %in% rimod$SAMPLEID,]
rimod <- rimod[rimod$SAMPLEID %in% md$sample,]
rimod <- rimod[match(md$sample,rimod$SAMPLEID),]

md$age <- rimod$AGE
md$pmd <- rimod$PMD.MIN.
md$gender <- rimod$GENDER
md$mutation <- rimod$gene

only.mirna$age <- c(NA, NA)
only.mirna$pmd <- c(NA, NA)
only.mirna$gender <- c(NA, NA)

md.all <- rbind(md, only.mirna)


write.table(md.all, "~/rimod/smallRNA/smallRNA_design_file.txt", sep="\t", row.names = F, quote=F)
