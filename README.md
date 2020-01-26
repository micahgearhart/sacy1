``` r
library("ggplot2")
library("magrittr")
library("dplyr")
library("magrittr")
library("rtracklayer")
library("GenomicFeatures")
library("GenomicAlignments")
library("RColorBrewer")
library("DESeq2")
library("goseq")
library("ComplexHeatmap")
library("Gviz")
library("BSgenome.Celegans.UCSC.ce11")
ce11<-BSgenome.Celegans.UCSC.ce11
suppressWarnings(seqlevelsStyle(ce11)<-"Ensembl")

ts<-format(Sys.time(), "%a_%b_%d_%Y_%H%M")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
```

``` r
gg_plotcounts <- function(x="WBGene00004984",d=cds,returnData=F) {
  if (substr(x,1,6)=="WBGene") {
    title<-symbols[grep(x,symbols$gene_id),"gene_name"][1]
  } else {
    x<-tolower(x)
    title<-x
    x<-symbols[grep(paste0("^",title,"$"),symbols$gene_name),"gene_id"][1]
  }
  
  if(returnData) {return(plotCounts(d,x,intgroup=c("tir1","aid"),returnData=T))}
  
  plotCounts(d,x,intgroup=c("tir1","aid"),returnData=T) %>%
    tibble::rownames_to_column() %>%
    mutate(aid=as.character(aid)) %>% 
    mutate(tir1=as.character(tir1)) %>% 
    mutate(strain=factor(paste0(aid,"_",tir1),levels=c("control_soma","experimental_soma","control_germline","experimental_germline"))) %>% 
    ggplot(aes(x=strain, y=count,colour=aid,shape=tir1)) +
    geom_point(position=position_jitter(w=0.1,h=0),size=3) + ggtitle(paste0(title," : ",x)) +
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                 geom = "crossbar", width = 0.35,size=0.4) +
    expand_limits(x=0, y = 0) + xlab("") + ylab("Normalized Counts") +
    scale_color_manual(values=c("#E69F00","#56B4E9")) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(legend.position="none")
}

listGO<-function(goid,degs) {
  #print(go_terms[goid])
  tg<-wb_go[grep(goid,wb_go$`GO ID`),"ensembl",drop=F]
  tg$go<-goid
  idx<-match(tg$ensembl,symbols$gene_id)
  tg$symbol<-symbols[idx,"gene_name"]
  tg$deg<-degs[tg$ensembl]
  as.data.frame(tg[tg$deg==1,])
}
```

Make Txdb from GTF file
=======================

``` r
gff3_file<-"data/Caenorhabditis_elegans.WBcel235.97.gtf"

  v97<-import(gff3_file,format="gtf")
  symbols<-as.data.frame(mcols(v97))
  symbols<-symbols[symbols$type=="transcript",c("gene_id","gene_name","gene_biotype","transcript_id")]
  table(symbols$gene_biotype)
```

    ## 
    ##  antisense_RNA        lincRNA          miRNA          ncRNA          piRNA 
    ##            104            184            718           7779          15363 
    ## protein_coding     pseudogene           rRNA         snoRNA          snRNA 
    ##          34214           1958             22            346            129 
    ##           tRNA 
    ##            634

``` r
#make txdb from gencode BASIC annotatios
suppressWarnings(txdb<-makeTxDbFromGFF(gff3_file,format="gtf"))
```

    ## Import genomic features from the file as a GRanges object ... OK
    ## Prepare the 'metadata' data frame ... OK
    ## Make the TxDb object ... OK

``` r
#seqlevelsStyle(txdb)<-"NCBI"
exons<-exonsBy(txdb,by="tx",use.names=T)
exons_per_transcript<-data.frame(num_exons=lengths(exons))
exons_per_transcript$gene<-symbols[match(rownames(exons_per_transcript),symbols$transcript_id),"gene_id"]
exons_per_transcript<-exons_per_transcript[with(exons_per_transcript,order(-num_exons)),]
exons_per_gene<-exons_per_transcript[!duplicated(exons_per_transcript$gene),]
rownames(exons_per_gene)<-exons_per_gene$gene
#exons<-exons[!names(exons) %in% symbols[symbols$gene_biotype=='rRNA',]$transcript_id]
introns<-intronsByTranscript(txdb,use.names=T)
introns_per_transcript<-data.frame(wid_intron=width(introns))
colnames(introns_per_transcript)<-c("group_num","transcript_id","width")
introns_per_transcript<-introns_per_transcript %>% 
  group_by(transcript_id) %>% 
  summarize(intron_sum=sum(width),intron_max=max(width)) %>% 
  ungroup()
introns_per_transcript$gene<-symbols[match(introns_per_transcript$transcript_id,symbols$transcript_id),"gene_id"]
introns_per_transcript$gene_name<-symbols[match(introns_per_transcript$transcript_id,symbols$transcript_id),"gene_name"]
introns_per_transcript<-introns_per_transcript[with(introns_per_transcript,order(-intron_sum)),]
introns_per_gene<-data.frame(introns_per_transcript[!duplicated(introns_per_transcript$gene),])
rownames(introns_per_gene)<-introns_per_gene$gene


exonsUL<-unlist(exons)
names(exonsUL)<-exonsUL$exon_name

#transcipts for Gviz
transcripts<-transcripts(txdb,use.names=TRUE)
transcripts$gene_id<-symbols[match(transcripts$tx_name,symbols$transcript_id),"gene_id"]

save(symbols,introns_per_gene,exons_per_gene,transcripts,file="symbols97.rdata")
```

Load Gene Lists Strome Soma/Germline and OMA1/LIN41 IPs
=======================================================

``` r
strome<-readxl::read_excel("data/TableS1.xlsx",sheet="H3K9me2 data") %>% dplyr::select(Gene.WB.ID,Gene.Public.Name,`soma-specific`,`germline-specific`) %>% as.data.frame
length(soma_genes<-strome[strome$`soma-specific`==1,"Gene.WB.ID"])
```

    ## [1] 1181

``` r
length(germline_genes<-strome[strome$`germline-specific`==1,"Gene.WB.ID"])
```

    ## [1] 169

``` r
#oma1/lin41 Enriched list
rnp_file<-"../lin41/enrichment_relative_to_LysateRZ_withGonadEnrichment_Fri_Apr_21_2017_1638.csv"
rnp<-read.csv(rnp_file,stringsAsFactors = F,header = T, row.names = 1)

length(lin41_list<-rownames(subset(rnp,Lin.41.IP > 2 & Lin41_padj < 0.05)))
```

    ## [1] 1115

``` r
length(oma1_list<-rownames(subset(rnp,Oma.1.IP > 2 & Oma1_padj < 0.05)))
```

    ## [1] 2259

\#DESeq2

``` r
load("featureCounts_SACY1_92_q55_Thu_Aug_15_2019_0844.rdata")
load("symbols97.rdata")

(coldata<-data.frame(counts=apply(unstranded_counts$counts,2,sum)))
```

    ##                                     counts
    ## WBcel235.CA1200.1.sort.bam.repair 26363866
    ## WBcel235.CA1200.2.sort.bam.repair 26260859
    ## WBcel235.CA1200.3.sort.bam.repair 28533206
    ## WBcel235.CA1352.1.sort.bam.repair 31004610
    ## WBcel235.CA1352.2.sort.bam.repair 33492636
    ## WBcel235.CA1352.3.sort.bam.repair 26381561
    ## WBcel235.DG4700.1.sort.bam.repair 25489209
    ## WBcel235.DG4700.2.sort.bam.repair 26399873
    ## WBcel235.DG4700.3.sort.bam.repair 30441183
    ## WBcel235.DG4703.1.sort.bam.repair 30131816
    ## WBcel235.DG4703.2.sort.bam.repair 31971368
    ## WBcel235.DG4703.3.sort.bam.repair 19177158

``` r
str(unstranded_counts)
```

    ## List of 4
    ##  $ counts    : int [1:46904, 1:12] 0 0 20 2393 19 32 1822 0 2 121 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:46904] "WBGene00197333" "WBGene00198386" "WBGene00015153" "WBGene00002061" ...
    ##   .. ..$ : chr [1:12] "WBcel235.CA1200.1.sort.bam.repair" "WBcel235.CA1200.2.sort.bam.repair" "WBcel235.CA1200.3.sort.bam.repair" "WBcel235.CA1352.1.sort.bam.repair" ...
    ##  $ annotation:'data.frame':  46904 obs. of  6 variables:
    ##   ..$ GeneID: chr [1:46904] "WBGene00197333" "WBGene00198386" "WBGene00015153" "WBGene00002061" ...
    ##   ..$ Chr   : chr [1:46904] "V" "V" "V;V;V;V;V;V" "V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V" ...
    ##   ..$ Start : chr [1:46904] "180" "180" "1663;2851;5538;5548;6024;6024" "6588;6589;6589;6589;6593;6593;6913;6913;6913;7158;7158;7158;7158;7158;7158;7433;7433;7433;7433;7433;7433;7651;7"| __truncated__ ...
    ##   ..$ End   : chr [1:46904] "329" "329" "1782;3019;5966;5966;6632;6632" "7110;6754;6754;6754;7110;7110;7110;7110;7110;7390;7384;7393;7390;7384;7393;7609;7609;7609;7609;7609;7609;7822;7"| __truncated__ ...
    ##   ..$ Strand: chr [1:46904] "+" "-" "+;+;+;+;+;+" "-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-" ...
    ##   ..$ Length: int [1:46904] 150 150 1327 1108 165 220 5597 21 106 553 ...
    ##  $ targets   : chr [1:12] "WBcel235.CA1200.1.sort.bam.repair" "WBcel235.CA1200.2.sort.bam.repair" "WBcel235.CA1200.3.sort.bam.repair" "WBcel235.CA1352.1.sort.bam.repair" ...
    ##  $ stat      :'data.frame':  12 obs. of  13 variables:
    ##   ..$ Status                           : chr [1:12] "Assigned" "Unassigned_Unmapped" "Unassigned_MappingQuality" "Unassigned_Chimera" ...
    ##   ..$ WBcel235.CA1200.1.sort.bam.repair: int [1:12] 26363866 2058088 0 0 1643280 0 1454265 0 0 240850 ...
    ##   ..$ WBcel235.CA1200.2.sort.bam.repair: int [1:12] 26260859 2355099 0 0 1351843 0 1658555 0 0 231443 ...
    ##   ..$ WBcel235.CA1200.3.sort.bam.repair: int [1:12] 28533206 2261523 0 0 1763922 0 1647373 0 0 257960 ...
    ##   ..$ WBcel235.CA1352.1.sort.bam.repair: int [1:12] 31004610 2483360 0 0 1856281 0 2082649 0 0 267507 ...
    ##   ..$ WBcel235.CA1352.2.sort.bam.repair: int [1:12] 33492636 2433381 0 0 2062518 0 2247710 0 0 321617 ...
    ##   ..$ WBcel235.CA1352.3.sort.bam.repair: int [1:12] 26381561 1939525 0 0 1681331 0 1895562 0 0 257004 ...
    ##   ..$ WBcel235.DG4700.1.sort.bam.repair: int [1:12] 25489209 1877600 0 0 1501024 0 1699790 0 0 248594 ...
    ##   ..$ WBcel235.DG4700.2.sort.bam.repair: int [1:12] 26399873 2078087 0 0 1588112 0 1813273 0 0 258322 ...
    ##   ..$ WBcel235.DG4700.3.sort.bam.repair: int [1:12] 30441183 2392910 0 0 1864281 0 2052135 0 0 309049 ...
    ##   ..$ WBcel235.DG4703.1.sort.bam.repair: int [1:12] 30131816 2314935 0 0 1763565 0 1677411 0 0 293198 ...
    ##   ..$ WBcel235.DG4703.2.sort.bam.repair: int [1:12] 31971368 2382282 0 0 1892245 0 1791443 0 0 307364 ...
    ##   ..$ WBcel235.DG4703.3.sort.bam.repair: int [1:12] 19177158 1815542 0 0 1159260 0 1109330 0 0 197397 ...

``` r
cds<-DESeqDataSetFromMatrix(unstranded_counts$counts,colData = coldata,design = ~1)
all.equal(rownames(cds),unstranded_counts$annotation$GeneID)
```

    ## [1] TRUE

``` r
mcols(cds)$basepairs<-unstranded_counts$annotation$Length
cds <- estimateSizeFactors(cds)
cds$filename<-rownames(colData(cds))
colnames(cds)<- cds$filename %>% gsub("WBcel235\\.","",.) %>% gsub("\\.sort\\.bam\\.repair","",.) %>% gsub("\\.","_rep",.)

cds$strain<-subseq(colnames(cds),1,6)
cds$rep<-subseq(colnames(cds),8,11)

cds$tir1<-factor(case_when(
  cds$strain == "CA1352" ~ "germline",
  cds$strain == "DG4700" ~ "germline",
  cds$strain == "CA1200" ~ "soma", 
  cds$strain == "DG4703" ~ "soma",
))

cds$aid<-factor(case_when(
  cds$strain == "CA1352" ~ "control",
  cds$strain == "DG4700" ~ "experimental",
  cds$strain == "CA1200" ~ "control",
  cds$strain == "DG4703" ~ "experimental",
))

cds$group<-factor(paste0(cds$tir1,"_",cds$aid))


cds<-cds[,!stringr::str_detect("CA1200_rep2",colnames(cds))]
as.data.frame(colData(cds))
```

    ##               counts sizeFactor                          filename strain  rep
    ## CA1200_rep1 26363866  1.0385536 WBcel235.CA1200.1.sort.bam.repair CA1200 rep1
    ## CA1200_rep3 28533206  1.1088140 WBcel235.CA1200.3.sort.bam.repair CA1200 rep3
    ## CA1352_rep1 31004610  1.0678817 WBcel235.CA1352.1.sort.bam.repair CA1352 rep1
    ## CA1352_rep2 33492636  1.2304098 WBcel235.CA1352.2.sort.bam.repair CA1352 rep2
    ## CA1352_rep3 26381561  0.9292652 WBcel235.CA1352.3.sort.bam.repair CA1352 rep3
    ## DG4700_rep1 25489209  0.9272360 WBcel235.DG4700.1.sort.bam.repair DG4700 rep1
    ## DG4700_rep2 26399873  0.9266389 WBcel235.DG4700.2.sort.bam.repair DG4700 rep2
    ## DG4700_rep3 30441183  1.0978403 WBcel235.DG4700.3.sort.bam.repair DG4700 rep3
    ## DG4703_rep1 30131816  1.0911330 WBcel235.DG4703.1.sort.bam.repair DG4703 rep1
    ## DG4703_rep2 31971368  1.1920995 WBcel235.DG4703.2.sort.bam.repair DG4703 rep2
    ## DG4703_rep3 19177158  0.7178937 WBcel235.DG4703.3.sort.bam.repair DG4703 rep3
    ##                 tir1          aid                 group
    ## CA1200_rep1     soma      control          soma_control
    ## CA1200_rep3     soma      control          soma_control
    ## CA1352_rep1 germline      control      germline_control
    ## CA1352_rep2 germline      control      germline_control
    ## CA1352_rep3 germline      control      germline_control
    ## DG4700_rep1 germline experimental germline_experimental
    ## DG4700_rep2 germline experimental germline_experimental
    ## DG4700_rep3 germline experimental germline_experimental
    ## DG4703_rep1     soma experimental     soma_experimental
    ## DG4703_rep2     soma experimental     soma_experimental
    ## DG4703_rep3     soma experimental     soma_experimental

``` r
  plotPCA(normTransform(cds),intgroup=c("tir1","aid"),returnData=T)
```

    ##                    PC1        PC2                 group     tir1          aid
    ## CA1200_rep1  -1.754614  15.787932          soma:control     soma      control
    ## CA1200_rep3  -5.172494  17.934475          soma:control     soma      control
    ## CA1352_rep1  19.809312  -7.660557      germline:control germline      control
    ## CA1352_rep2  20.014277  -6.457336      germline:control germline      control
    ## CA1352_rep3  24.800092 -13.973031      germline:control germline      control
    ## DG4700_rep1  13.835096   5.902847 germline:experimental germline experimental
    ## DG4700_rep2  13.748825   6.022152 germline:experimental germline experimental
    ## DG4700_rep3   9.465578   1.103153 germline:experimental germline experimental
    ## DG4703_rep1 -28.508233 -11.446225     soma:experimental     soma experimental
    ## DG4703_rep2 -31.209559  -1.545414     soma:experimental     soma experimental
    ## DG4703_rep3 -35.028280  -5.667996     soma:experimental     soma experimental
    ##                    name
    ## CA1200_rep1 CA1200_rep1
    ## CA1200_rep3 CA1200_rep3
    ## CA1352_rep1 CA1352_rep1
    ## CA1352_rep2 CA1352_rep2
    ## CA1352_rep3 CA1352_rep3
    ## DG4700_rep1 DG4700_rep1
    ## DG4700_rep2 DG4700_rep2
    ## DG4700_rep3 DG4700_rep3
    ## DG4703_rep1 DG4703_rep1
    ## DG4703_rep2 DG4703_rep2
    ## DG4703_rep3 DG4703_rep3

``` r
 plotPCA(normTransform(cds),intgroup=c("tir1","aid"),returnData=F,ntop=1000)+
    scale_color_manual(values=cbPalette[c(6,7,4,2)]) + theme_bw()
```

![](sacy1_polyA_files/figure-markdown_github/deseq2-1.png)

``` r
#svglite::svglite(file=paste0("2019 Manuscript/pca_",ts,".svg"),width=6,height=4)  
  plotPCA(normTransform(cds),intgroup=c("tir1","aid"),returnData=F,ntop=1000)+
    scale_color_manual(values=cbPalette[c(6,7,4,2)]) + theme_bw()
```

![](sacy1_polyA_files/figure-markdown_github/deseq2-2.png)

``` r
#dev.off()  
colData(cds)
```

    ## DataFrame with 11 rows and 8 columns
    ##                counts        sizeFactor                          filename
    ##             <integer>         <numeric>                       <character>
    ## CA1200_rep1  26363866  1.03855360354668 WBcel235.CA1200.1.sort.bam.repair
    ## CA1200_rep3  28533206  1.10881403155905 WBcel235.CA1200.3.sort.bam.repair
    ## CA1352_rep1  31004610   1.0678816600844 WBcel235.CA1352.1.sort.bam.repair
    ## CA1352_rep2  33492636  1.23040977851531 WBcel235.CA1352.2.sort.bam.repair
    ## CA1352_rep3  26381561 0.929265180890905 WBcel235.CA1352.3.sort.bam.repair
    ## DG4700_rep1  25489209 0.927236015333653 WBcel235.DG4700.1.sort.bam.repair
    ## DG4700_rep2  26399873 0.926638899480908 WBcel235.DG4700.2.sort.bam.repair
    ## DG4700_rep3  30441183  1.09784026274059 WBcel235.DG4700.3.sort.bam.repair
    ## DG4703_rep1  30131816  1.09113296535249 WBcel235.DG4703.1.sort.bam.repair
    ## DG4703_rep2  31971368  1.19209953957569 WBcel235.DG4703.2.sort.bam.repair
    ## DG4703_rep3  19177158 0.717893722227499 WBcel235.DG4703.3.sort.bam.repair
    ##                  strain         rep     tir1          aid                 group
    ##             <character> <character> <factor>     <factor>              <factor>
    ## CA1200_rep1      CA1200        rep1     soma      control          soma_control
    ## CA1200_rep3      CA1200        rep3     soma      control          soma_control
    ## CA1352_rep1      CA1352        rep1 germline      control      germline_control
    ## CA1352_rep2      CA1352        rep2 germline      control      germline_control
    ## CA1352_rep3      CA1352        rep3 germline      control      germline_control
    ## DG4700_rep1      DG4700        rep1 germline experimental germline_experimental
    ## DG4700_rep2      DG4700        rep2 germline experimental germline_experimental
    ## DG4700_rep3      DG4700        rep3 germline experimental germline_experimental
    ## DG4703_rep1      DG4703        rep1     soma experimental     soma_experimental
    ## DG4703_rep2      DG4703        rep2     soma experimental     soma_experimental
    ## DG4703_rep3      DG4703        rep3     soma experimental     soma_experimental

``` r
gg_plotcounts("tra-2")
```

![](sacy1_polyA_files/figure-markdown_github/deseq2-3.png)

``` r
gg_plotcounts("her-1")
```

![](sacy1_polyA_files/figure-markdown_github/deseq2-4.png)

``` r
f <- fpkm(cds,robust=TRUE)

f_mean <- f %>% as.data.frame() %>% 
  tibble::rownames_to_column(var="gene_id") %>% 
  tidyr::gather(sample,fpkm,-gene_id) %>%
  dplyr::mutate(group=colData(cds)[sample,]$group) %>% 
  dplyr::group_by(gene_id,group) %>% 
    dplyr::summarize(mean_fpkm=round(mean(fpkm),3)) %>% 
    #dplyr::summarize(mean_fpkm=round(log2(mean(fpkm)+0.01),3)) %>% 
    ungroup() %>% 
  dplyr::select(gene_id,group,mean_fpkm) %>% 
  tidyr::spread(group,mean_fpkm) %>%  as.data.frame() 

rownames(f_mean)<-f_mean$gene_id
f_mean<-f_mean[rownames(cds),-1]
head(f_mean)
```

    ##                germline_control germline_experimental soma_control
    ## WBGene00197333            0.000                 0.000        0.000
    ## WBGene00198386            0.000                 0.000        0.000
    ## WBGene00015153            0.771                 0.887        0.529
    ## WBGene00002061           95.109                82.935       76.691
    ## WBGene00255704            4.537                 5.825        4.543
    ## WBGene00235314            6.705                 5.642        4.429
    ##                soma_experimental
    ## WBGene00197333             0.000
    ## WBGene00198386             0.000
    ## WBGene00015153             0.572
    ## WBGene00002061            88.874
    ## WBGene00255704             5.669
    ## WBGene00235314             4.193

``` r
f_mean[unique(symbols[grep("her-1",symbols$gene_name),"gene_id"]),]
```

    ##                germline_control germline_experimental soma_control
    ## WBGene00001842            4.102                16.995        9.228
    ##                soma_experimental
    ## WBGene00001842             8.794

``` r
f_mean[unique(symbols[grep("gld-1",symbols$gene_name),"gene_id"]),]
```

    ##                germline_control germline_experimental soma_control
    ## WBGene00001595           66.756                55.791       52.638
    ##                soma_experimental
    ## WBGene00001595            61.056

``` r
f_mean[unique(symbols[grep("tra-2",symbols$gene_name),"gene_id"]),]
```

    ##                germline_control germline_experimental soma_control
    ## WBGene00006605            7.487                 6.962        7.107
    ##                soma_experimental
    ## WBGene00006605              6.88

Define sets of Differentially expressed genes
=============================================

``` r
colData(cds)
```

    ## DataFrame with 11 rows and 8 columns
    ##                counts        sizeFactor                          filename
    ##             <integer>         <numeric>                       <character>
    ## CA1200_rep1  26363866  1.03855360354668 WBcel235.CA1200.1.sort.bam.repair
    ## CA1200_rep3  28533206  1.10881403155905 WBcel235.CA1200.3.sort.bam.repair
    ## CA1352_rep1  31004610   1.0678816600844 WBcel235.CA1352.1.sort.bam.repair
    ## CA1352_rep2  33492636  1.23040977851531 WBcel235.CA1352.2.sort.bam.repair
    ## CA1352_rep3  26381561 0.929265180890905 WBcel235.CA1352.3.sort.bam.repair
    ## DG4700_rep1  25489209 0.927236015333653 WBcel235.DG4700.1.sort.bam.repair
    ## DG4700_rep2  26399873 0.926638899480908 WBcel235.DG4700.2.sort.bam.repair
    ## DG4700_rep3  30441183  1.09784026274059 WBcel235.DG4700.3.sort.bam.repair
    ## DG4703_rep1  30131816  1.09113296535249 WBcel235.DG4703.1.sort.bam.repair
    ## DG4703_rep2  31971368  1.19209953957569 WBcel235.DG4703.2.sort.bam.repair
    ## DG4703_rep3  19177158 0.717893722227499 WBcel235.DG4703.3.sort.bam.repair
    ##                  strain         rep     tir1          aid                 group
    ##             <character> <character> <factor>     <factor>              <factor>
    ## CA1200_rep1      CA1200        rep1     soma      control          soma_control
    ## CA1200_rep3      CA1200        rep3     soma      control          soma_control
    ## CA1352_rep1      CA1352        rep1 germline      control      germline_control
    ## CA1352_rep2      CA1352        rep2 germline      control      germline_control
    ## CA1352_rep3      CA1352        rep3 germline      control      germline_control
    ## DG4700_rep1      DG4700        rep1 germline experimental germline_experimental
    ## DG4700_rep2      DG4700        rep2 germline experimental germline_experimental
    ## DG4700_rep3      DG4700        rep3 germline experimental germline_experimental
    ## DG4703_rep1      DG4703        rep1     soma experimental     soma_experimental
    ## DG4703_rep2      DG4703        rep2     soma experimental     soma_experimental
    ## DG4703_rep3      DG4703        rep3     soma experimental     soma_experimental

``` r
design(cds)<- ~group
cds<-DESeq(cds)
```

    ## using pre-existing size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
#plotMA(cds,ylim=c(-10,10)) 
resultsNames(cds)
```

    ## [1] "Intercept"                                      
    ## [2] "group_germline_experimental_vs_germline_control"
    ## [3] "group_soma_control_vs_germline_control"         
    ## [4] "group_soma_experimental_vs_germline_control"

Define DEGs in Soma
-------------------

``` r
summary(res_soma<-results(cds,contrast=c("group","soma_experimental","soma_control"),alpha=0.05))
```

    ## 
    ## out of 27700 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 1352, 4.9%
    ## LFC < 0 (down)     : 1735, 6.3%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 10349, 37%
    ## (mean count < 8)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
all.equal(rownames(res_soma),rownames(f_mean))
```

    ## [1] TRUE

``` r
res_soma<-cbind(as.data.frame(res_soma),f_mean)
res_soma$gene_name<-symbols[match(rownames(res_soma),symbols$gene_id),"gene_name"]
res_soma$gene_biotype<-symbols[match(rownames(res_soma),symbols$gene_id),"gene_biotype"]


dim(temp_soma<-subset(res_soma, baseMean > 25 & padj<0.05 & abs(log2FoldChange) > 1 &
                    (germline_control >= 2.5 | germline_experimental >= 2.5 |
                       soma_control >= 2.5 | soma_experimental >= 2.5 )))
```

    ## [1] 484  12

``` r
table(sign(temp_soma$log2FoldChange))
```

    ## 
    ##  -1   1 
    ## 242 242

``` r
length(soma_up<-rownames(subset(temp_soma,log2FoldChange > 1)))
```

    ## [1] 242

``` r
length(soma_dn<-rownames(subset(temp_soma,log2FoldChange < -1)))
```

    ## [1] 242

``` r
all(c("col-17", "col-41", "col-46", "col-47", "col-90", "col-128", "col-149", "dpy-3", "dpy-4", 
  "dpy-5", "dpy-6", "dpy-8", "dpy-9", "dpy-13", "dpy-20", "lon-3", "mlt-7", "qua-1", "rol-6", "rol-8", "sqt-1","sqt-2") %in% res_soma[soma_dn,"gene_name"])
```

    ## [1] TRUE

Define DEGs in germline
-----------------------

``` r
summary(res_germ<-results(cds,contrast=c("group","germline_experimental","germline_control"),alpha=0.05))
```

    ## 
    ## out of 27700 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 2700, 9.7%
    ## LFC < 0 (down)     : 2292, 8.3%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 10349, 37%
    ## (mean count < 8)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
stopifnot(all.equal(rownames(res_germ),rownames(f_mean)))
res_germ<-cbind(as.data.frame(res_germ),f_mean)
res_germ$gene_name<-symbols[match(rownames(res_germ),symbols$gene_id),"gene_name"]
res_germ$gene_biotype<-symbols[match(rownames(res_germ),symbols$gene_id),"gene_biotype"]


dim(temp_germ<-subset(res_germ, baseMean > 25 & padj<0.05 & abs(log2FoldChange) > 1 &
                    (germline_control >= 2.5 | germline_experimental >= 2.5 |
                       soma_control >= 2.5 | soma_experimental >= 2.5 )))
```

    ## [1] 437  12

``` r
table(sign(temp_germ$log2FoldChange))
```

    ## 
    ##  -1   1 
    ## 126 311

``` r
length(germline_up<-rownames(subset(temp_germ,log2FoldChange > 1)))
```

    ## [1] 311

``` r
res_germ[res_germ$gene_name %in% c("her-1","xol-1","tra-2"),] #sex-determination pathway
```

    ##                  baseMean log2FoldChange     lfcSE       stat       pvalue
    ## WBGene00001842  144.51609      2.0486952 0.2529907  8.0979084 5.591215e-16
    ## WBGene00006962   41.14594      1.6918722 0.3951274  4.2818396 1.853546e-05
    ## WBGene00006605 1004.17507     -0.1035522 0.1625503 -0.6370472 5.240941e-01
    ##                        padj germline_control germline_experimental soma_control
    ## WBGene00001842 1.347405e-13            4.102                16.995        9.228
    ## WBGene00006962 2.041961e-04            0.268                 0.859        0.691
    ## WBGene00006605 6.694807e-01            7.487                 6.962        7.107
    ##                soma_experimental gene_name   gene_biotype
    ## WBGene00001842             8.794     her-1 protein_coding
    ## WBGene00006962             0.629     xol-1 protein_coding
    ## WBGene00006605             6.880     tra-2 protein_coding

``` r
length(germline_dn<-rownames(subset(temp_germ,log2FoldChange < -1)))
```

    ## [1] 126

Combine soma/germline results for supplementary table
-----------------------------------------------------

``` r
temp <-full_join(tibble::rownames_to_column(res_soma,var="WormbaseID"),tibble::rownames_to_column(res_germ,var="WormbaseID"),
                 by=c("WormbaseID","gene_name","gene_biotype","baseMean",
                 "soma_control","soma_experimental","germline_control","germline_experimental"
                 ),suffix=c(".soma",".germline")) %>% 
  dplyr::select(WormbaseID,gene_name,gene_biotype,baseMean,soma_control,soma_experimental,germline_control,germline_experimental,
                log2FoldChange.soma,padj.soma,log2FoldChange.germline,padj.germline) %>% 
  dplyr::filter(baseMean > 25 & (germline_control >= 2.5 | germline_experimental >= 2.5 | soma_control >= 2.5 | soma_experimental >= 2.5 )) %>% 
  dplyr::filter((padj.soma <0.05 & abs(log2FoldChange.soma) > 1) | (padj.germline <0.05 & abs(log2FoldChange.germline) > 1)) %>% 
  dplyr::arrange(padj.germline) 

#Check it
stopifnot(all.equal(sort(soma_up),sort(subset(temp,log2FoldChange.soma > 1 & padj.soma < 0.05)$WormbaseID)))
stopifnot(all.equal(sort(soma_dn),sort(subset(temp,log2FoldChange.soma < -1 & padj.soma < 0.05)$WormbaseID)))
stopifnot(all.equal(sort(germline_up),sort(subset(temp,log2FoldChange.germline > 1 & padj.germline < 0.05)$WormbaseID)))
stopifnot(all.equal(sort(germline_dn),sort(subset(temp,log2FoldChange.germline < -1 & padj.germline < 0.05)$WormbaseID)))

#Rename Columns
  temp2 <- temp %>% 
  dplyr::rename("Gene Name" = gene_name) %>% 
  dplyr::rename("Gene Type" = gene_biotype) %>% 
  dplyr::rename("Mean Counts Across Samples" = baseMean) %>% 
  dplyr::rename("Soma Control FPKM" = soma_control) %>% 
  dplyr::rename("Soma Deplete FPKM" = soma_experimental) %>% 
  dplyr::rename("Germline Control FPKM" = germline_control) %>% 
  dplyr::rename("Germline Deplete FPKM" = germline_experimental) %>% 
  dplyr::rename("Log2 Fold Change Soma"= log2FoldChange.soma) %>%
  dplyr::rename("BH Corrected p-value Soma"= padj.soma) %>%
  dplyr::rename("Log2 Fold Change Germline"= log2FoldChange.germline) %>%
  dplyr::rename("BH Corrected p-value Germline"= padj.germline)


write.csv(temp2,file=paste0("2019 Manuscript/sacy1_deplete_degs_",ts,".csv"),quote=F,row.names=F)
```

Four way Venn on global\_results
--------------------------------

``` r
#svglite::svglite(file=paste0("2019 Manuscript/venn_",ts,".svg"),width=5,height=5)
grid.newpage()
VennDiagram::draw.quad.venn(area1=length(germline_up),
                            area2=length(germline_dn),
                            area3=length(soma_up),
                            area4=length(soma_dn),
                            n12=length(intersect(germline_up,germline_dn)),
                            n13=length(intersect(germline_up,soma_up)),
                            n14=length(intersect(germline_up,soma_dn)),
                            n23=length(intersect(germline_dn,soma_up)),
                            n24=length(intersect(germline_dn,soma_dn)),
                            n34=length(intersect(soma_up,soma_dn)),
                            n123=length(intersect(intersect(germline_up,germline_dn),soma_up)),
                            n124=length(intersect(intersect(germline_up,germline_dn),soma_dn)),
                            n134=length(intersect(intersect(germline_up,soma_up),soma_dn)),
                            n234=length(intersect(intersect(germline_dn,soma_up),soma_dn)),
                            n1234=length(intersect(intersect(germline_up,germline_dn),intersect(soma_up,soma_dn))),
                            fill=cbPalette[c(7,6,2,4)],
                            category=c("germline_up","germline_dn","soma_up","soma_dn"))
```

![](sacy1_polyA_files/figure-markdown_github/four_way_venn-1.png)

    ## (polygon[GRID.polygon.270], polygon[GRID.polygon.271], polygon[GRID.polygon.272], polygon[GRID.polygon.273], polygon[GRID.polygon.274], polygon[GRID.polygon.275], polygon[GRID.polygon.276], polygon[GRID.polygon.277], text[GRID.text.278], text[GRID.text.279], text[GRID.text.280], text[GRID.text.281], text[GRID.text.282], text[GRID.text.283], text[GRID.text.284], text[GRID.text.285], text[GRID.text.286], text[GRID.text.287], text[GRID.text.288], text[GRID.text.289], text[GRID.text.290], text[GRID.text.291], text[GRID.text.292], text[GRID.text.293], text[GRID.text.294], text[GRID.text.295], text[GRID.text.296])

``` r
#dev.off()

temp_soma[intersect(soma_up,germline_up),] #17
```

    ##                  baseMean log2FoldChange     lfcSE      stat       pvalue
    ## WBGene00017413   45.47375       1.033115 0.3627446  2.848052 4.398773e-03
    ## WBGene00016361   50.40427       1.017239 0.3804612  2.673701 7.501925e-03
    ## WBGene00006097   62.58709       1.755602 0.3614292  4.857389 1.189436e-06
    ## WBGene00016927   44.88284       1.193389 0.3808829  3.133217 1.729013e-03
    ## WBGene00219905   32.27106       2.709835 0.4781176  5.667717 1.447128e-08
    ## WBGene00008267   81.35764       1.137052 0.3002092  3.787533 1.521503e-04
    ## WBGene00014069   37.10703       1.232075 0.4395744  2.802882 5.064821e-03
    ## WBGene00044638 4413.07895       1.010652 0.1977980  5.109516 3.229857e-07
    ## WBGene00000675 4944.68697       1.972059 0.2053820  9.601909 7.847744e-22
    ## WBGene00235263   34.02978       1.406309 0.4372210  3.216471 1.297774e-03
    ## WBGene00018643  103.92132       1.361498 0.3709340  3.670458 2.421164e-04
    ## WBGene00020905   38.45928       2.131551 0.4753538  4.484137 7.320962e-06
    ## WBGene00003448   29.15098       2.348230 0.4889165  4.802926 1.563636e-06
    ## WBGene00003432   45.46353       1.595444 0.5242192  3.043467 2.338688e-03
    ## WBGene00002103   30.40456       1.905697 0.5590904  3.408567 6.530505e-04
    ## WBGene00004965   42.78701       2.036328 0.4193070  4.856412 1.195322e-06
    ## WBGene00017270  109.24651       3.263450 0.2856571 11.424361 3.159709e-30
    ##                        padj germline_control germline_experimental soma_control
    ## WBGene00017413 2.995404e-02            1.801                 3.842        1.941
    ## WBGene00016361 4.406429e-02            2.241                 4.936        2.222
    ## WBGene00006097 3.394392e-05            0.380                 3.183        1.099
    ## WBGene00016927 1.477110e-02            0.438                 1.835        1.143
    ## WBGene00219905 6.767954e-07            3.922                 9.142        2.648
    ## WBGene00008267 2.151556e-03            1.100                 3.466        0.967
    ## WBGene00014069 3.309970e-02            1.462                 3.616        1.594
    ## WBGene00044638 1.098848e-05          431.562              1011.085      416.585
    ## WBGene00000675 2.669926e-19           95.992               213.434       60.863
    ## WBGene00235263 1.180791e-02            3.699                 7.708        4.416
    ## WBGene00018643 3.118754e-03            2.716                 6.252        1.629
    ## WBGene00020905 1.664823e-04            1.094                 2.518        0.720
    ## WBGene00003448 4.390071e-05            1.011                 2.339        0.872
    ## WBGene00003432 1.859697e-02            1.869                 4.499        1.851
    ## WBGene00002103 6.836071e-03            1.848                 4.334        0.715
    ## WBGene00004965 3.400004e-05            0.720                 1.696        0.674
    ## WBGene00017270 2.108619e-27            0.664                 1.324        1.465
    ##                soma_experimental gene_name   gene_biotype
    ## WBGene00017413             3.932   F13A2.4 protein_coding
    ## WBGene00016361             4.511    drd-10 protein_coding
    ## WBGene00006097             3.738    str-31 protein_coding
    ## WBGene00016927             2.627   nhr-172 protein_coding
    ## WBGene00219905            17.123 F19F10.13          ncRNA
    ## WBGene00008267             2.137   C53A5.9 protein_coding
    ## WBGene00014069             3.791   ZK678.3 protein_coding
    ## WBGene00044638           839.733   F23A7.8 protein_coding
    ## WBGene00000675           238.741   col-101 protein_coding
    ## WBGene00235263            11.855  Y66H1A.8 protein_coding
    ## WBGene00018643             4.142    drd-50 protein_coding
    ## WBGene00020905             3.122  T28H11.7 protein_coding
    ## WBGene00003448             4.468    msp-55 protein_coding
    ## WBGene00003432             5.586    msp-36 protein_coding
    ## WBGene00002103             2.675    ins-20 protein_coding
    ## WBGene00004965             2.758    spe-11 protein_coding
    ## WBGene00017270            14.041    numr-1 protein_coding

``` r
temp_soma[intersect(soma_dn,germline_dn),] #8
```

    ##                  baseMean log2FoldChange     lfcSE      stat       pvalue
    ## WBGene00044109  114.40863      -2.614063 0.5492786 -4.759085 1.944727e-06
    ## WBGene00000618 1154.40697      -3.579762 0.3811196 -9.392752 5.845246e-21
    ## WBGene00019662   34.56442      -2.454711 0.4524887 -5.424911 5.798352e-08
    ## WBGene00019495 5919.32239      -1.387550 0.1471454 -9.429792 4.109044e-21
    ## WBGene00021690  285.45012      -1.473027 0.3131185 -4.704375 2.546453e-06
    ## WBGene00017068  170.90377      -1.076005 0.3111763 -3.457862 5.444810e-04
    ## WBGene00013701  152.75592      -1.138598 0.2471573 -4.606775 4.089621e-06
    ## WBGene00013130   34.20326      -3.968491 0.9596154 -4.135501 3.541808e-05
    ##                        padj germline_control germline_experimental soma_control
    ## WBGene00044109 5.305496e-05            8.340                 3.565        0.664
    ## WBGene00000618 1.748635e-18           70.938                34.240        6.952
    ## WBGene00019662 2.367228e-06            2.678                 1.212        6.852
    ## WBGene00019495 1.288098e-18          388.238               172.813      248.199
    ## WBGene00021690 6.654142e-05           67.561                31.832       10.919
    ## WBGene00017068 5.949175e-03           25.265                11.556        5.173
    ## WBGene00013701 1.010812e-04            7.386                 2.447        4.326
    ## WBGene00013130 6.382609e-04            9.590                 4.085        1.858
    ##                soma_experimental  gene_name   gene_biotype
    ## WBGene00044109             0.104  K02E11.10 protein_coding
    ## WBGene00000618             0.587     col-41 protein_coding
    ## WBGene00019662             1.268   K11H12.6 protein_coding
    ## WBGene00019495            94.869     sdz-24 protein_coding
    ## WBGene00021690             3.886 Y48G8AL.12 protein_coding
    ## WBGene00017068             2.435    D2092.8 protein_coding
    ## WBGene00013701             1.966  Y106G6D.4 protein_coding
    ## WBGene00013130             0.100  Y52B11B.1 protein_coding

``` r
temp_soma[intersect(soma_up,germline_dn),] #13
```

    ##                 baseMean log2FoldChange     lfcSE      stat       pvalue
    ## WBGene00002016 2463.3726       1.527342 0.1303672 11.715699 1.059141e-31
    ## WBGene00002018 3946.0711       1.622787 0.1277816 12.699694 5.936118e-37
    ## WBGene00021601  266.0784       1.092049 0.2207976  4.945929 7.578147e-07
    ## WBGene00002015  765.5775       1.256667 0.1531833  8.203682 2.331336e-16
    ## WBGene00002019  679.6127       1.167491 0.2126217  5.490931 3.998213e-08
    ## WBGene00002020 1569.1360       1.404208 0.1440446  9.748427 1.873487e-22
    ## WBGene00002017  384.4505       1.238306 0.1952043  6.343643 2.243943e-10
    ## WBGene00219327  131.1207       2.252773 0.2483913  9.069450 1.196193e-19
    ## WBGene00009691 2433.5848       1.333108 0.1474873  9.038803 1.583971e-19
    ## WBGene00009692 1865.6592       1.411346 0.1537723  9.178151 4.385563e-20
    ## WBGene00012603  162.5798       1.204195 0.2317010  5.197192 2.023211e-07
    ## WBGene00002026 2057.9675       1.291487 0.1521479  8.488371 2.095544e-17
    ## WBGene00000667  250.4192       1.264959 0.2734734  4.625527 3.736465e-06
    ##                        padj germline_control germline_experimental soma_control
    ## WBGene00002016 7.657149e-29          308.565               102.989       38.415
    ## WBGene00002018 5.420926e-34          493.773               141.999       47.013
    ## WBGene00021601 2.319020e-05           10.110                 3.933        2.975
    ## WBGene00002015 3.816133e-14          128.666                41.316       19.675
    ## WBGene00002019 1.697636e-06          126.895                32.845       16.263
    ## WBGene00002020 7.223750e-20          268.734                90.060       36.088
    ## WBGene00002017 1.426178e-08           57.087                19.603        9.804
    ## WBGene00219327 2.965020e-17           19.201                 9.572        3.535
    ## WBGene00009691 3.817150e-17           86.064                28.184       11.835
    ## WBGene00009692 1.207840e-17           68.481                21.353        8.083
    ## WBGene00012603 7.164232e-06           28.015                11.729       19.065
    ## WBGene00002026 3.787477e-15           72.537                20.282       12.197
    ## WBGene00000667 9.368699e-05           15.837                 7.098        3.549
    ##                soma_experimental gene_name   gene_biotype
    ## WBGene00002016           110.797  hsp-16.2 protein_coding
    ## WBGene00002018           144.831 hsp-16.41 protein_coding
    ## WBGene00021601             6.334  Y46H3A.5 protein_coding
    ## WBGene00002015            47.076  hsp-16.1 protein_coding
    ## WBGene00002019            36.585 hsp-16.48 protein_coding
    ## WBGene00002020            95.468 hsp-16.49 protein_coding
    ## WBGene00002017            23.133 hsp-16.11 protein_coding
    ## WBGene00219327            16.834   ZK822.9 protein_coding
    ## WBGene00009691            29.816   F44E5.4 protein_coding
    ## WBGene00009692            21.492   F44E5.5 protein_coding
    ## WBGene00012603            43.986    nspe-6 protein_coding
    ## WBGene00002026            29.844    hsp-70 protein_coding
    ## WBGene00000667             8.486    col-92 protein_coding

``` r
temp_soma[intersect(soma_dn,germline_up),] #1
```

    ##                baseMean log2FoldChange     lfcSE      stat      pvalue
    ## WBGene00021162 25.80832      -1.889059 0.5040611 -3.747678 0.000178479
    ##                       padj germline_control germline_experimental soma_control
    ## WBGene00021162 0.002434583            0.877                 3.111        2.869
    ##                soma_experimental gene_name   gene_biotype
    ## WBGene00021162             0.776   Y5H2A.1 protein_coding

her-1 and tra-2 bigWig diagrams
===============================

``` r
options(ucscChromosomeNames=FALSE)

plotTranscript<-function(common_name="her-1",ylim=500,pad=c(100,500),rev=FALSE,logT=FALSE) {

gene<-unique(symbols[grep(common_name,symbols$gene_name),"gene_id"])

gene_start<-start(transcripts[transcripts$gene_id==gene])
if (length(gene_start) > 1) { gene_start<-min(gene_start) }
gene_end<-end(transcripts[transcripts$gene_id==gene])
if (length(gene_end) > 1) { gene_end<-max(gene_end) }
gene_chr<-unique(as.character(seqnames(transcripts[transcripts$gene_id==gene])))

myImportFun <- function(file, selection){
  (fls<-list.files("data","*.bam$",full.names=T) %>% stringr::str_subset("CA1352|DG4700"))
  s_names<-sapply(strsplit(basename(fls),"\\."),function(x) x[1]) %>% gsub("-","_rep",.)
  weights<-colData(cds[,s_names])$sizeFactor
  param_total <- ScanBamParam(what="mapq",flag = scanBamFlag(isUnmappedQuery = FALSE,isMinusStrand = NA,isMateMinusStrand = NA,isProperPair = TRUE),which=selection)
  gr<-GRanges(seqnames=seqnames(selection),IRanges(start=start(selection):end(selection),width=1),strand="*")
  mcols(gr)<-do.call(cbind,lapply(fls,function(file) as.numeric(coverage(suppressWarnings(GenomicAlignments::readGAlignmentPairs(file,use.names=F,param=param_total)))[[as.character(seqnames(selection))]])[start(selection):end(selection)]))
  mcols(gr)<-sweep(as.matrix(mcols(gr)),MARGIN=2,weights,"/")
  return(gr)
}


#myImportFun2("data/DG4700-1.Aligned.out.sort.dedup.unique.bam",GRanges(seqnames=gene_chr,IRanges(start=gene_start,end=gene_end),strand="*"))
paddedLog2<-function(x) {
  if(logT) return(log2(x+0.1))
  else(return(x))
  }
gtrack <- GenomeAxisTrack(col="darkgray")
txTr <- GeneRegionTrack(txdb, chromosome =  gene_chr, start = gene_start-1, end = gene_end+1,fill="gray",col="black",fontcolor.group="black",fill="darkgray",name="WS273 Transcripts")
dT<-DataTrack(range="data/DG4700-1.Aligned.out.sort.dedup.unique.bam", genome="ce11", type="p", name="Log2 Coverage", chromosome=gene_chr,importFunction = myImportFun,stream=T,col=cbPalette[6:7])
plotTracks(list(gtrack,dT,txTr), from=gene_start-pad[1], to = gene_end+pad[2],transformation=paddedLog2,reverseStrand=rev,
           cex=0.8,add53=TRUE,ylim=c(0,ylim),type=c("a","confint"),
           groups=rep(c("Control", "Deplete"), each=3),
           col.axis="black",background.title="transparent",fontcolor.title="black") 

}

svglite::svglite(file=paste0("2019 Manuscript/her1_log2_",ts,".svg"),width=5,height=3) 
plotTranscript("her-1",ylim=10,pad=c(200,500),logT=T,rev=F)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
svglite::svglite(file=paste0("2019 Manuscript/tra2_log2_",ts,".svg"),width=5,height=3) 
plotTranscript("tra-2",ylim=10,pad=c(25,100),logT=T,rev=T)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

Look at alternative splicing events
===================================

``` r
import_rmats_dataset<-function(dataset) {
  results<-list()
  results[["A3SS"]]<-read.table(paste0("rmats/rmats_4.0.2.",dataset,".kimble/A3SS.MATS.JCEC.txt"),stringsAsFactors = F,header=T)[,c("GeneID","geneSymbol","FDR","IncLevelDifference","chr","longExonStart_0base","longExonEnd","IncLevelDifference")]
  colnames(results[["A3SS"]])<-c("GeneID","geneSymbol","FDR","IncLevelDifference","chr","ExonStart","ExonEnd","IncLevelDifference")
  results[["A5SS"]]<-read.table(paste0("rmats/rmats_4.0.2.",dataset,".kimble/A5SS.MATS.JCEC.txt"),stringsAsFactors = F,header=T)[,c("GeneID","geneSymbol","FDR","IncLevelDifference","chr","longExonStart_0base","longExonEnd","IncLevelDifference")]
  colnames(results[["A5SS"]])<-c("GeneID","geneSymbol","FDR","IncLevelDifference","chr","ExonStart","ExonEnd","IncLevelDifference")
  results[["MXE"]]<-read.table(paste0("rmats/rmats_4.0.2.",dataset,".kimble/MXE.MATS.JCEC.txt"),stringsAsFactors = F,header=T)[,c("GeneID","geneSymbol","FDR","IncLevelDifference","chr","X1stExonStart_0base","X1stExonEnd","IncLevelDifference")]
  colnames(results[["MXE"]])<-c("GeneID","geneSymbol","FDR","IncLevelDifference","chr","ExonStart","ExonEnd","IncLevelDifference")
  results[["RI"]]<-read.table(paste0("rmats/rmats_4.0.2.",dataset,".kimble/RI.MATS.JCEC.txt"),stringsAsFactors = F,header=T)[,c("GeneID","geneSymbol","FDR","IncLevelDifference","chr","riExonStart_0base","riExonEnd","IncLevelDifference")]
  colnames(results[["RI"]])<-c("GeneID","geneSymbol","FDR","IncLevelDifference","chr","ExonStart","ExonEnd","IncLevelDifference")
  results[["SE"]]<-read.table(paste0("rmats/rmats_4.0.2.",dataset,".kimble/SE.MATS.JCEC.txt"),stringsAsFactors = F,header=T)[,c("GeneID","geneSymbol","FDR","IncLevelDifference","chr","exonStart_0base","exonEnd","IncLevelDifference")]
  colnames(results[["SE"]])<-c("GeneID","geneSymbol","FDR","IncLevelDifference","chr","ExonStart","ExonEnd","IncLevelDifference")
  return(bind_rows(results,.id="event"))
}

rmats_soma<-import_rmats_dataset("soma")
rmats_soma<-rmats_soma[with(rmats_soma,order(FDR)),] %>%  distinct(event,GeneID,chr,ExonStart,ExonEnd,.keep_all=TRUE) # removes less significant flanking exon hits, eg
rmats_germline<-import_rmats_dataset("germline")
rmats_germline<-rmats_germline[with(rmats_germline,order(FDR)),] %>%  distinct(event,GeneID,chr,ExonStart,ExonEnd,.keep_all=TRUE)
rmats_ortiz<-import_rmats_dataset("ortiz")
rmats_ortiz<-rmats_ortiz[with(rmats_ortiz,order(FDR)),] %>%  distinct(event,GeneID,chr,ExonStart,ExonEnd,.keep_all=TRUE)



temp <-full_join(rmats_soma,rmats_germline,by=c("event","GeneID","geneSymbol","chr","ExonStart","ExonEnd"),suffix=c(".soma",".germline"))
#rmats_results <- left_join(temp,rmats_ortiz,by=c("event","GeneID","geneSymbol","chr","ExonStart","ExonEnd"),suffix=c(".sacy1",".ortiz"))
rmats_results <- full_join(temp,rmats_ortiz,by=c("event","GeneID","geneSymbol","chr","ExonStart","ExonEnd"),suffix=c(".sacy1",".ortiz"))

head(rmats_results)
```

    ##   event      GeneID geneSymbol FDR.soma IncLevelDifference.soma    chr
    ## 1  A3SS  MSTRG.8077       <NA>        0                  -0.486 chrIII
    ## 2  A3SS  MSTRG.9065       <NA>        0                  -0.443 chrIII
    ## 3  A3SS  MSTRG.7825       <NA>        0                  -0.350 chrIII
    ## 4  A3SS   MSTRG.925    F27C1.2        0                  -0.328   chrI
    ## 5  A3SS  MSTRG.7066       <NA>        0                  -0.283 chrIII
    ## 6  A3SS MSTRG.17325       <NA>        0                  -0.438   chrV
    ##   ExonStart  ExonEnd FDR.germline IncLevelDifference.germline          FDR
    ## 1   7332939  7333431  0.047382532                      -0.105 0.000000e+00
    ## 2  11688125 11688336  0.494178016                      -0.030 8.542106e-09
    ## 3   6225716  6225850  0.438506510                      -0.041 2.104121e-11
    ## 4   5430760  5430992  0.002916894                      -0.083 0.000000e+00
    ## 5   2544834  2544895  1.000000000                      -0.011 0.000000e+00
    ## 6    754242   754415  1.000000000                      -0.006           NA
    ##   IncLevelDifference
    ## 1              0.529
    ## 2              0.303
    ## 3              0.466
    ## 4              0.200
    ## 5              0.065
    ## 6                 NA

``` r
table(rmats_results$sig<-factor(case_when(
  rmats_results$FDR.germline < 0.05 & rmats_results$FDR.soma < 0.05  ~ "both", 
  rmats_results$FDR.germline < 0.05 & (rmats_results$FDR.soma >= 0.05 | is.na(rmats_results$FDR.soma)) ~ "germline", 
  (rmats_results$FDR.germline >= 0.05 | is.na(rmats_results$FDR.germline) ) & rmats_results$FDR.soma < 0.05  ~ "soma", 
  (rmats_results$FDR.germline >= 0.05 | is.na(rmats_results$FDR.germline) ) & (rmats_results$FDR.soma >= 0.05 | is.na(rmats_results$FDR.soma)) & rmats_results$FDR < 0.05 ~ "Ortiz_sig",
  (rmats_results$FDR.germline >= 0.05 | is.na(rmats_results$FDR.germline) ) & (rmats_results$FDR.soma >= 0.05 | is.na(rmats_results$FDR.soma)) ~ "notsig"
)),useNA="always")
```

    ## 
    ##      both  germline    notsig Ortiz_sig      soma      <NA> 
    ##       236       560     17041       930      1370         0

``` r
nrow(rmats_results<-rmats_results[rmats_results$sig !="notsig",])
```

    ## [1] 3096

``` r
nrow(rmats_results[rmats_results$sig =="both" | rmats_results$sig =="soma",])
```

    ## [1] 1606

``` r
nrow(rmats_results[rmats_results$sig =="both" | rmats_results$sig =="germline",])
```

    ## [1] 796

``` r
table(rmats_results$ortiz_sig<- factor(case_when(
  (rmats_results$sig == "germline" | rmats_results$sig == "both") & sign(rmats_results$IncLevelDifference.germline) < 0 & sign(rmats_results$IncLevelDifference) > 0 & rmats_results$FDR < 0.05 ~ "Oocyte Enriched",
  (rmats_results$sig == "germline" | rmats_results$sig == "both") & sign(rmats_results$IncLevelDifference.germline) > 0 & sign(rmats_results$IncLevelDifference) < 0 & rmats_results$FDR < 0.05 ~ "Oocyte Enriched",
  (rmats_results$sig == "germline" | rmats_results$sig == "both") & sign(rmats_results$IncLevelDifference.germline) < 0 & sign(rmats_results$IncLevelDifference) < 0 & rmats_results$FDR < 0.05 ~ "Sperm Enriched",
  (rmats_results$sig == "germline" | rmats_results$sig == "both") & sign(rmats_results$IncLevelDifference.germline) > 0 & sign(rmats_results$IncLevelDifference) > 0 & rmats_results$FDR < 0.05 ~ "Sperm Enriched",
  rmats_results$sig == "soma" & sign(rmats_results$IncLevelDifference.soma) < 0 & sign(rmats_results$IncLevelDifference) > 0 & rmats_results$FDR < 0.05 ~ "Oocyte Enriched",
  rmats_results$sig == "soma" & sign(rmats_results$IncLevelDifference.soma) > 0 & sign(rmats_results$IncLevelDifference) < 0 & rmats_results$FDR < 0.05 ~ "Oocyte Enriched",
  rmats_results$sig == "soma" & sign(rmats_results$IncLevelDifference.soma) < 0 & sign(rmats_results$IncLevelDifference) < 0 & rmats_results$FDR < 0.05 ~ "Sperm Enriched",
  rmats_results$sig == "soma" & sign(rmats_results$IncLevelDifference.soma) > 0 & sign(rmats_results$IncLevelDifference) > 0 & rmats_results$FDR < 0.05 ~ "Sperm Enriched",
  rmats_results$sig ==  "Ortiz_sig" ~ "TRUE",
  is.na(rmats_results$FDR) ~ "not sig",
  TRUE ~ "not sig"
  )),useNA="always")
```

    ## 
    ##         not sig Oocyte Enriched  Sperm Enriched            TRUE            <NA> 
    ##            1496             614              56             930               0

``` r
#how many in Oritz
nrow(rmats_results[rmats_results$ortiz_sig =="Sperm Enriched" | rmats_results$ortiz_sig =="Oocyte Enriched" | rmats_results$ortiz_sig=="TRUE",])
```

    ## [1] 1600

``` r
sum(rmats_ortiz$FDR < 0.05)
```

    ## [1] 1600

``` r
#Identify overlapping Transcript (Stringtie Merge doesn't always pull the right name from reference annotation)
get_overlapping_transcripts<-function(gr) { paste(unique(symbols[symbols$gene_id %in% unique(transcripts[subjectHits(findOverlaps(gr,transcripts))]$gene_id),"gene_name"]),collapse=";") }
rmats_ranges<-GRanges(seqnames=gsub("chr","",rmats_results$chr),IRanges(start=rmats_results$ExonStart,end=rmats_results$ExonEnd),strand="*")
rmats_results$overlapping_transcript<-"unknown"

rmats_ranges$name<-rmats_results$sig
export(rmats_ranges,"rmats_ranges.bed")
for(i in 1:length(rmats_ranges)) {rmats_results$overlapping_transcript[i]<-get_overlapping_transcripts(rmats_ranges[i])}  #non-vectorized slow step

#svglite::svglite(file=paste0("2019 Manuscript/splicing_effects_",ts,".svg"),width=8,height=4) 
rmats_results %>% 
  dplyr::filter(sig != "Ortiz_sig") %>%  #remove ortiz data from figure
  distinct(event,GeneID,overlapping_transcript,.keep_all=TRUE) %>% 
  dplyr::mutate(event=factor(event,levels=c("A5SS","A3SS","RI","MXE","SE"))) %>% 
  dplyr::mutate(sig=factor(sig,levels=c("soma","both","germline"))) %>% 
  dplyr::mutate(ortiz_sig=factor(ortiz_sig,levels=c("Oocyte Enriched","Sperm Enriched","not sig"))) %>% 
  ggplot(aes(event)) + geom_bar(aes(fill=ortiz_sig),color="black") + facet_grid(~sig) +
  scale_fill_manual(values=c("magenta","navyblue","gray"),name="Different\nin \nSperm \nvs \nOocyte",) +
  ylab("Number of Genes Affected") + xlab("") + ggtitle("Genomewide SACY-1 Splicing Effects") + theme_bw() 
```

![](sacy1_polyA_files/figure-markdown_github/import_rmats-1.png)

``` r
#dev.off()


temp<-rmats_results %>% 
  dplyr::mutate(gr=paste0(chr,":",ExonStart,"-",ExonEnd)) %>% 
  dplyr::select(overlapping_transcript,GeneID,geneSymbol,gr,event,sig,FDR.soma,IncLevelDifference.soma,
                                      FDR.germline,IncLevelDifference.germline,ortiz_sig,FDR,IncLevelDifference) %>% 
  dplyr::rename(FDR.ortiz=FDR,IncLevelDifference.ortiz=IncLevelDifference) %>% 
  dplyr::mutate(sig=factor(sig,levels=c("germline","both","soma","Ortiz_sig"))) %>% 
  dplyr::mutate(event=factor(event,levels=c("A3SS","A5SS","RI","SE","MXE"))) %>% 
  arrange(sig,event,FDR.germline,FDR.ortiz) 

temp2<-temp %>% 
   dplyr::rename("Overlapping Transcript" = overlapping_transcript) %>%
   dplyr::rename("StringTie Gene ID" = GeneID) %>% 
   dplyr::rename("StringTie Defined Symbol" = geneSymbol) %>% 
   dplyr::rename("Affected Exon Coordinates" = gr) %>% 
   dplyr::rename("Splicing Event Type" = event) %>% 
   dplyr::rename("Significant Dataset" = sig) %>% 
   dplyr::rename("False Discovery Rate q-value Soma" = FDR.soma) %>% 
   dplyr::rename("Inclusion Level Difference Soma" = IncLevelDifference.soma) %>% 
   dplyr::rename("False Discovery Rate q-value Germline" = FDR.germline) %>% 
   dplyr::rename("Inclusion Level Difference Germline" = IncLevelDifference.germline) %>% 
   dplyr::rename("Relative Enrichment Sperm/Oocyte" = ortiz_sig) %>% 
   dplyr::rename("False Discovery Rate q-value Sperm/Oocyte" = FDR.ortiz) %>% 
   dplyr::rename("Inclusion Level Difference Sperm/Oocyte" = IncLevelDifference.ortiz)  
  
write.csv(temp2,file=paste0("sacy1_alternative_splicing_events_",ts,".csv"),quote=F,row.names=F)

#find out how many unique genes are alternatively spliced in Ortiz dataset
nrow(rmats_ortiz_subset<-subset(rmats_ortiz,FDR<0.05))
```

    ## [1] 1600

``` r
rmats_ortiz_subset_ranges<-GRanges(seqnames=gsub("chr","",rmats_ortiz_subset$chr),IRanges(start=rmats_ortiz_subset$ExonStart,end=rmats_ortiz_subset$ExonEnd),strand="*")
rmats_ortiz_subset$overlapping_transcript<-"unknown"
for(i in 1:length(rmats_ortiz_subset)) {rmats_ortiz_subset$overlapping_transcript[i]<-get_overlapping_transcripts(rmats_ortiz_subset_ranges[i])}
table(rmats_ortiz_subset$overlapping_transcript)
```

    ## 
    ##          cdk-8       F56C9.10         him-18         nkcc-1          sax-3 
    ##              1              2              1              1              1 
    ##          ula-1        unknown         usp-50 ZK1127.3;cth-2 
    ##              1           1591              1              1

``` r
#table(rmats_ortiz_subset$geneSymbol,useNA="always")


nrow(x<- rmats_ortiz_subset %>%  distinct(GeneID,.keep_all=TRUE) )
```

    ## [1] 1071

``` r
nrow(y<- rmats_results %>%  distinct(GeneID,.keep_all=TRUE) )
```

    ## [1] 2074

``` r
(x[paste(x$GeneID,x$event,x$overlapping_transcript) %in% paste(y$GeneID,y$event,y$overlapping_transcript),])
```

    ##   event      GeneID geneSymbol FDR IncLevelDifference    chr ExonStart  ExonEnd
    ## 1  A3SS MSTRG.14364       <NA>   0              0.314  chrIV  14115406 14115556
    ## 2  A3SS  MSTRG.7483       <NA>   0              0.263 chrIII   4712820  4713104
    ## 3  A3SS  MSTRG.8077       <NA>   0              0.289 chrIII   7332412  7332835
    ## 4  A3SS MSTRG.21348     usp-50   0              0.293   chrV  17929614 17929780
    ## 5  A3SS MSTRG.22917       <NA>   0              0.686   chrX   3452617  3452939
    ## 6  A3SS  MSTRG.1671       <NA>   0              0.349   chrI   8704532  8704760
    ## 7  A3SS  MSTRG.4774       <NA>   0              0.422  chrII   7055925  7056185
    ##   overlapping_transcript
    ## 1                 nkcc-1
    ## 2                 him-18
    ## 3               F56C9.10
    ## 4                 usp-50
    ## 5                  sax-3
    ## 6                  cdk-8
    ## 7         ZK1127.3;cth-2

``` r
(y[paste(y$GeneID,y$event,y$overlapping_transcript) %in% paste(x$GeneID,x$event,x$overlapping_transcript),])
```

    ##      event      GeneID geneSymbol     FDR.soma IncLevelDifference.soma    chr
    ## 1     A3SS  MSTRG.8077       <NA> 0.000000e+00                  -0.486 chrIII
    ## 245   A3SS MSTRG.22917       <NA> 5.144261e-08                  -0.343   chrX
    ## 473   A3SS MSTRG.21348     usp-50 1.031779e-04                  -0.503   chrV
    ## 602   A3SS  MSTRG.1671       <NA> 8.451629e-04                  -0.134   chrI
    ## 832   A3SS  MSTRG.4774       <NA> 8.416536e-03                  -0.080  chrII
    ## 1415  A3SS MSTRG.14364       <NA> 2.717880e-01                  -0.245  chrIV
    ## 1477  A3SS  MSTRG.7483       <NA> 3.876641e-01                  -0.059 chrIII
    ##      ExonStart  ExonEnd FDR.germline IncLevelDifference.germline        FDR
    ## 1      7332939  7333431 0.0473825316                      -0.105 0.00000000
    ## 245    3452617  3452939 0.5181523008                      -0.051 0.00000000
    ## 473   17928542 17928653 0.7816319579                       0.008 0.01217625
    ## 602    8704532  8704760 0.0001294857                      -0.164 0.00000000
    ## 832    7055925  7056185 0.0138506093                      -0.079 0.00000000
    ## 1415  14115406 14115556 0.2658660878                      -0.180 0.00000000
    ## 1477   4712820  4713104 0.3811244119                      -0.077 0.00000000
    ##      IncLevelDifference       sig       ortiz_sig overlapping_transcript
    ## 1                 0.529      both Oocyte Enriched               F56C9.10
    ## 245               0.686      soma Oocyte Enriched                  sax-3
    ## 473               0.269      soma Oocyte Enriched                 usp-50
    ## 602               0.349      both Oocyte Enriched                  cdk-8
    ## 832               0.422      both Oocyte Enriched         ZK1127.3;cth-2
    ## 1415              0.314 Ortiz_sig            TRUE                 nkcc-1
    ## 1477              0.263 Ortiz_sig            TRUE                 him-18

``` r
sum(paste(x$GeneID) %in% paste(y$GeneID))
```

    ## [1] 1071

``` r
goi<-c("sel-10", "nos-3", "cye-1", "nos-2", "gld-1","cpb-1", "cdk-8", "rde-11", "ptc-1", "ife-3", "lin-54", "egl-45", "eif-3.2", "wdr-4", "taf-10", "picc-1")

#View(rmats_results[rmats_results$overlapping_transcript %in% goi,])
```

import counts for sperm/oocyte data for normalization
=====================================================

``` r
load("featureCounts_SACY1_92_with_kimble_q55_Wed_Dec_11_2019_1120.rdata")

(coldata2<-data.frame(counts=apply(unstranded_counts$counts,2,sum)))
```

    ##                                                 counts
    ## WBcel235.CA1200.1.Aligned.out.sort.bam.repair 24909469
    ## WBcel235.CA1200.2.Aligned.out.sort.bam.repair 25398667
    ## WBcel235.CA1200.3.Aligned.out.sort.bam.repair 27024526
    ## WBcel235.CA1352.1.Aligned.out.sort.bam.repair 29373050
    ## WBcel235.CA1352.2.Aligned.out.sort.bam.repair 31319499
    ## WBcel235.CA1352.3.Aligned.out.sort.bam.repair 24699836
    ## WBcel235.DG4700.1.Aligned.out.sort.bam.repair 23890410
    ## WBcel235.DG4700.2.Aligned.out.sort.bam.repair 24956332
    ## WBcel235.DG4700.3.Aligned.out.sort.bam.repair 28743897
    ## WBcel235.DG4703.1.Aligned.out.sort.bam.repair 28348437
    ## WBcel235.DG4703.2.Aligned.out.sort.bam.repair 30067629
    ## WBcel235.DG4703.3.Aligned.out.sort.bam.repair 18319285
    ## ...kimble..SRX527951.q4.bam.repair            40557963
    ## ...kimble..SRX527952.q4.bam.repair            27394835
    ## ...kimble..SRX527953.q4.bam.repair            32281817
    ## ...kimble..SRX527954.q4.bam.repair            39458193
    ## ...kimble..SRX527955.q4.bam.repair            26345161
    ## ...kimble..SRX527956.q4.bam.repair            35923324
    ## ...kimble..SRX527957.q4.bam.repair            31981458
    ## ...kimble..SRX527958.q4.bam.repair            43778644
    ## ...kimble..SRX527959.q4.bam.repair            37788135
    ## ...kimble..SRX527960.q4.bam.repair            30248649
    ## ...kimble..SRX527961.q4.bam.repair            27265518
    ## ...kimble..SRX527962.q4.bam.repair            32140380
    ## ...kimble..SRX527963.q4.bam.repair            34391135
    ## ...kimble..SRX527964.q4.bam.repair            34904021
    ## ...kimble..SRX527965.q4.bam.repair            35960311
    ## ...kimble..SRX527966.q4.bam.repair            33720283

``` r
str(unstranded_counts)
```

    ## List of 4
    ##  $ counts    : int [1:46904, 1:28] 0 0 21 2463 19 19 1299 0 2 38 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:46904] "WBGene00197333" "WBGene00198386" "WBGene00015153" "WBGene00002061" ...
    ##   .. ..$ : chr [1:28] "WBcel235.CA1200.1.Aligned.out.sort.bam.repair" "WBcel235.CA1200.2.Aligned.out.sort.bam.repair" "WBcel235.CA1200.3.Aligned.out.sort.bam.repair" "WBcel235.CA1352.1.Aligned.out.sort.bam.repair" ...
    ##  $ annotation:'data.frame':  46904 obs. of  6 variables:
    ##   ..$ GeneID: chr [1:46904] "WBGene00197333" "WBGene00198386" "WBGene00015153" "WBGene00002061" ...
    ##   ..$ Chr   : chr [1:46904] "V" "V" "V;V;V;V;V;V" "V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V;V" ...
    ##   ..$ Start : chr [1:46904] "180" "180" "1663;2851;5538;5548;6024;6024" "6588;6589;6589;6589;6593;6593;6913;6913;6913;7158;7158;7158;7158;7158;7158;7433;7433;7433;7433;7433;7433;7651;7"| __truncated__ ...
    ##   ..$ End   : chr [1:46904] "329" "329" "1782;3019;5966;5966;6632;6632" "7110;6754;6754;6754;7110;7110;7110;7110;7110;7390;7384;7393;7390;7384;7393;7609;7609;7609;7609;7609;7609;7822;7"| __truncated__ ...
    ##   ..$ Strand: chr [1:46904] "+" "-" "+;+;+;+;+;+" "-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-;-" ...
    ##   ..$ Length: int [1:46904] 150 150 1327 1108 165 220 5597 21 106 553 ...
    ##  $ targets   : chr [1:28] "WBcel235.CA1200.1.Aligned.out.sort.bam.repair" "WBcel235.CA1200.2.Aligned.out.sort.bam.repair" "WBcel235.CA1200.3.Aligned.out.sort.bam.repair" "WBcel235.CA1352.1.Aligned.out.sort.bam.repair" ...
    ##  $ stat      :'data.frame':  12 obs. of  29 variables:
    ##   ..$ Status                                       : chr [1:12] "Assigned" "Unassigned_Unmapped" "Unassigned_MappingQuality" "Unassigned_Chimera" ...
    ##   ..$ WBcel235.CA1200.1.Aligned.out.sort.bam.repair: int [1:12] 24909469 1789 0 0 4195742 0 1046496 0 0 252118 ...
    ##   ..$ WBcel235.CA1200.2.Aligned.out.sort.bam.repair: int [1:12] 25398667 1300 0 0 3377479 0 1153263 0 0 244643 ...
    ##   ..$ WBcel235.CA1200.3.Aligned.out.sort.bam.repair: int [1:12] 27024526 1638 0 0 4478357 0 1146187 0 0 269753 ...
    ##   ..$ WBcel235.CA1352.1.Aligned.out.sort.bam.repair: int [1:12] 29373050 2330 0 0 4828065 0 1631497 0 0 279551 ...
    ##   ..$ WBcel235.CA1352.2.Aligned.out.sort.bam.repair: int [1:12] 31319499 2335 0 0 5686569 0 1646973 0 0 336568 ...
    ##   ..$ WBcel235.CA1352.3.Aligned.out.sort.bam.repair: int [1:12] 24699836 1893 0 0 4527165 0 1363479 0 0 269613 ...
    ##   ..$ WBcel235.DG4700.1.Aligned.out.sort.bam.repair: int [1:12] 23890410 1608 0 0 4236421 0 1275471 0 0 259711 ...
    ##   ..$ WBcel235.DG4700.2.Aligned.out.sort.bam.repair: int [1:12] 24956332 2146 0 0 4146523 0 1288288 0 0 269935 ...
    ##   ..$ WBcel235.DG4700.3.Aligned.out.sort.bam.repair: int [1:12] 28743897 2191 0 0 4870085 0 1490647 0 0 327264 ...
    ##   ..$ WBcel235.DG4703.1.Aligned.out.sort.bam.repair: int [1:12] 28348437 1964 0 0 4859437 0 1099761 0 0 308020 ...
    ##   ..$ WBcel235.DG4703.2.Aligned.out.sort.bam.repair: int [1:12] 30067629 2417 0 0 5244734 0 1204625 0 0 322287 ...
    ##   ..$ WBcel235.DG4703.3.Aligned.out.sort.bam.repair: int [1:12] 18319285 1205 0 0 2917716 0 722606 0 0 206976 ...
    ##   ..$ ...kimble..SRX527951.q4.bam.repair           : int [1:12] 40557963 0 41141231 0 0 0 0 0 0 167856 ...
    ##   ..$ ...kimble..SRX527952.q4.bam.repair           : int [1:12] 27394835 0 27816058 0 0 0 0 0 0 81722 ...
    ##   ..$ ...kimble..SRX527953.q4.bam.repair           : int [1:12] 32281817 0 32754313 0 0 0 0 0 0 128234 ...
    ##   ..$ ...kimble..SRX527954.q4.bam.repair           : int [1:12] 39458193 0 40026338 0 0 0 0 0 0 153684 ...
    ##   ..$ ...kimble..SRX527955.q4.bam.repair           : int [1:12] 26345161 0 26734929 0 0 0 0 0 0 108850 ...
    ##   ..$ ...kimble..SRX527956.q4.bam.repair           : int [1:12] 35923324 0 36434295 0 0 0 0 0 0 131709 ...
    ##   ..$ ...kimble..SRX527957.q4.bam.repair           : int [1:12] 31981458 0 32499798 0 0 0 0 0 0 143270 ...
    ##   ..$ ...kimble..SRX527958.q4.bam.repair           : int [1:12] 43778644 0 44379511 0 0 0 0 0 0 138264 ...
    ##   ..$ ...kimble..SRX527959.q4.bam.repair           : int [1:12] 37788135 0 38432915 0 0 0 0 0 0 193894 ...
    ##   ..$ ...kimble..SRX527960.q4.bam.repair           : int [1:12] 30248649 0 30740830 0 0 0 0 0 0 131891 ...
    ##   ..$ ...kimble..SRX527961.q4.bam.repair           : int [1:12] 27265518 0 27709002 0 0 0 0 0 0 132492 ...
    ##   ..$ ...kimble..SRX527962.q4.bam.repair           : int [1:12] 32140380 0 32602484 0 0 0 0 0 0 90132 ...
    ##   ..$ ...kimble..SRX527963.q4.bam.repair           : int [1:12] 34391135 0 34970281 0 0 0 0 0 0 173095 ...
    ##   ..$ ...kimble..SRX527964.q4.bam.repair           : int [1:12] 34904021 0 35489546 0 0 0 0 0 0 174416 ...
    ##   ..$ ...kimble..SRX527965.q4.bam.repair           : int [1:12] 35960311 0 36550438 0 0 0 0 0 0 149825 ...
    ##   ..$ ...kimble..SRX527966.q4.bam.repair           : int [1:12] 33720283 0 34292534 0 0 0 0 0 0 164388 ...

``` r
cds2<-DESeqDataSetFromMatrix(unstranded_counts$counts,colData = coldata2,design = ~1)
stopifnot(all.equal(rownames(cds2),unstranded_counts$annotation$GeneID))
mcols(cds2)$basepairs<-unstranded_counts$annotation$Length
cds2 <- estimateSizeFactors(cds2)
cds2$filename<-rownames(colData(cds2))
colnames(cds2)<- cds2$filename %>% gsub("WBcel235\\.","",.) %>% 
  gsub("\\.Aligned\\.out","",.) %>% 
  gsub("\\.\\.\\.kimble\\.\\.","",.) %>% 
  gsub("\\.q4\\.bam\\.repair","",.) %>% 
  gsub("\\.sort\\.bam\\.repair","",.) %>% gsub("\\.","_rep",.)

cds2$strain<-sapply(strsplit(colnames(cds2),"_"), function(x) x[1])
cds2$rep<-sapply(strsplit(colnames(cds2),"_"), function(x) x[2])
#cds2$strain<-subseq(colnames(cds2),1,6)
#cds2$rep<-subseq(colnames(cds2),8,11)

cds2$tir1<-factor(case_when(
  cds2$strain == "CA1352" ~ "germline", 
  cds2$strain == "DG4700" ~ "germline",
  cds2$strain == "CA1200" ~ "soma", 
  cds2$strain == "DG4703" ~ "soma",
  TRUE ~ "ortiz"
))

cds2$aid<-factor(case_when(
  cds2$strain == "CA1352" ~ "control",
  cds2$strain == "DG4700" ~ "experimental",
  cds2$strain == "CA1200" ~ "control",
  cds2$strain == "DG4703" ~ "experimental",
  TRUE ~ "ortiz"
))

cds2$group<-paste0(cds2$tir1,"_",cds2$aid)
# oocytes [fog-2(q71)]
colData(cds2)[rownames(colData(cds2)) %>% stringr::str_subset("27951|27952|27953|27954|27955|27956|27957|27958"),]$group<-"fog2"
# sperm [fem-3(q96)] 
colData(cds2)[rownames(colData(cds2)) %>% stringr::str_subset("27959|27960|27961|27962|27963|27964|27965|27966"),]$group<-"fem3"
cds2$group<-factor(cds2$group,levels=c("soma_control","soma_experimental","germline_control","germline_experimental","fem3", "fog2"))

cds2<-cds2[,!stringr::str_detect("CA1200_rep2",colnames(cds2))]
as.data.frame(colData(cds2))
```

    ##               counts sizeFactor                                      filename
    ## CA1200_rep1 24909469  1.0669938 WBcel235.CA1200.1.Aligned.out.sort.bam.repair
    ## CA1200_rep3 27024526  1.1892308 WBcel235.CA1200.3.Aligned.out.sort.bam.repair
    ## CA1352_rep1 29373050  1.0720565 WBcel235.CA1352.1.Aligned.out.sort.bam.repair
    ## CA1352_rep2 31319499  1.2002190 WBcel235.CA1352.2.Aligned.out.sort.bam.repair
    ## CA1352_rep3 24699836  0.8800891 WBcel235.CA1352.3.Aligned.out.sort.bam.repair
    ## DG4700_rep1 23890410  1.0418887 WBcel235.DG4700.1.Aligned.out.sort.bam.repair
    ## DG4700_rep2 24956332  1.0017141 WBcel235.DG4700.2.Aligned.out.sort.bam.repair
    ## DG4700_rep3 28743897  1.1609627 WBcel235.DG4700.3.Aligned.out.sort.bam.repair
    ## DG4703_rep1 28348437  1.2234244 WBcel235.DG4703.1.Aligned.out.sort.bam.repair
    ## DG4703_rep2 30067629  1.3517324 WBcel235.DG4703.2.Aligned.out.sort.bam.repair
    ## DG4703_rep3 18319285  0.8238747 WBcel235.DG4703.3.Aligned.out.sort.bam.repair
    ## SRX527951   40557963  1.2368743            ...kimble..SRX527951.q4.bam.repair
    ## SRX527952   27394835  0.8392385            ...kimble..SRX527952.q4.bam.repair
    ## SRX527953   32281817  1.0137049            ...kimble..SRX527953.q4.bam.repair
    ## SRX527954   39458193  1.2324900            ...kimble..SRX527954.q4.bam.repair
    ## SRX527955   26345161  0.7633368            ...kimble..SRX527955.q4.bam.repair
    ## SRX527956   35923324  1.0225077            ...kimble..SRX527956.q4.bam.repair
    ## SRX527957   31981458  0.7317396            ...kimble..SRX527957.q4.bam.repair
    ## SRX527958   43778644  1.4907323            ...kimble..SRX527958.q4.bam.repair
    ## SRX527959   37788135  0.9956877            ...kimble..SRX527959.q4.bam.repair
    ## SRX527960   30248649  0.7374596            ...kimble..SRX527960.q4.bam.repair
    ## SRX527961   27265518  0.6839448            ...kimble..SRX527961.q4.bam.repair
    ## SRX527962   32140380  0.7693382            ...kimble..SRX527962.q4.bam.repair
    ## SRX527963   34391135  0.8609340            ...kimble..SRX527963.q4.bam.repair
    ## SRX527964   34904021  0.8988578            ...kimble..SRX527964.q4.bam.repair
    ## SRX527965   35960311  1.0417226            ...kimble..SRX527965.q4.bam.repair
    ## SRX527966   33720283  0.8969055            ...kimble..SRX527966.q4.bam.repair
    ##                strain  rep     tir1          aid                 group
    ## CA1200_rep1    CA1200 rep1     soma      control          soma_control
    ## CA1200_rep3    CA1200 rep3     soma      control          soma_control
    ## CA1352_rep1    CA1352 rep1 germline      control      germline_control
    ## CA1352_rep2    CA1352 rep2 germline      control      germline_control
    ## CA1352_rep3    CA1352 rep3 germline      control      germline_control
    ## DG4700_rep1    DG4700 rep1 germline experimental germline_experimental
    ## DG4700_rep2    DG4700 rep2 germline experimental germline_experimental
    ## DG4700_rep3    DG4700 rep3 germline experimental germline_experimental
    ## DG4703_rep1    DG4703 rep1     soma experimental     soma_experimental
    ## DG4703_rep2    DG4703 rep2     soma experimental     soma_experimental
    ## DG4703_rep3    DG4703 rep3     soma experimental     soma_experimental
    ## SRX527951   SRX527951 <NA>    ortiz        ortiz                  fog2
    ## SRX527952   SRX527952 <NA>    ortiz        ortiz                  fog2
    ## SRX527953   SRX527953 <NA>    ortiz        ortiz                  fog2
    ## SRX527954   SRX527954 <NA>    ortiz        ortiz                  fog2
    ## SRX527955   SRX527955 <NA>    ortiz        ortiz                  fog2
    ## SRX527956   SRX527956 <NA>    ortiz        ortiz                  fog2
    ## SRX527957   SRX527957 <NA>    ortiz        ortiz                  fog2
    ## SRX527958   SRX527958 <NA>    ortiz        ortiz                  fog2
    ## SRX527959   SRX527959 <NA>    ortiz        ortiz                  fem3
    ## SRX527960   SRX527960 <NA>    ortiz        ortiz                  fem3
    ## SRX527961   SRX527961 <NA>    ortiz        ortiz                  fem3
    ## SRX527962   SRX527962 <NA>    ortiz        ortiz                  fem3
    ## SRX527963   SRX527963 <NA>    ortiz        ortiz                  fem3
    ## SRX527964   SRX527964 <NA>    ortiz        ortiz                  fem3
    ## SRX527965   SRX527965 <NA>    ortiz        ortiz                  fem3
    ## SRX527966   SRX527966 <NA>    ortiz        ortiz                  fem3

``` r
  plotPCA(normTransform(cds2),intgroup=c("group"),returnData=F)
```

![](sacy1_polyA_files/figure-markdown_github/deseq2_for_sperm_oocyte_dataset-1.png)

make transcript plots with sperm/oocyte data
============================================

``` r
options(ucscChromosomeNames=FALSE)
paddedLog2<-function(x) { if(logT) return(log2(x+0.1)) else(return(x)) }

plotTranscript<-function(x="fog-1",ylim=500,pad=c(100,500),rev=FALSE,logT=FALSE) {

if (stringr::str_detect(x,":")) {
    x<-gsub("-",":",x) %>% gsub(",","",.)
    gene_chr<-sapply(strsplit(x,":"),function(y) y[1]) %>% gsub("chr","",.)
    gene_start<-as.numeric(sapply(strsplit(x,":"),function(y) y[2]))
    gene_end<-as.numeric(sapply(strsplit(x,":"),function(y) y[3]))
  } else {
    gene<-unique(symbols[grep(common_name,symbols$gene_name),"gene_id"])
    gene_start<-start(transcripts[transcripts$gene_id==gene])
    if (length(gene_start) > 1) { gene_start<-min(gene_start) }
    gene_end<-end(transcripts[transcripts$gene_id==gene])
    if (length(gene_end) > 1) { gene_end<-max(gene_end) }
    gene_chr<-unique(as.character(seqnames(transcripts[transcripts$gene_id==gene])))
  }
  
myImportFun<- function(file, selection){
  (fls<-list.files("data","*.bam$",full.names=T) %>% stringr::str_subset("CA1200-1|CA1200-3|DG4703|CA1352|DG4700"))
  s_names<-sapply(strsplit(basename(fls),"\\."),function(x) x[1]) %>% gsub("-","_rep",.)
  weights<-colData(cds[,s_names])$sizeFactor
  param_total <- ScanBamParam(what="mapq",flag = scanBamFlag(isUnmappedQuery = FALSE,isMinusStrand = NA,isMateMinusStrand = NA,isProperPair = TRUE),which=selection)
  gr<-GRanges(seqnames=seqnames(selection),IRanges(start=start(selection):end(selection),width=1),strand="*")
  mcols(gr)<-do.call(cbind,lapply(fls,function(file) as.numeric(coverage(suppressWarnings(GenomicAlignments::readGAlignmentPairs(file,use.names=F,param=param_total)))[[as.character(seqnames(selection))]])[start(selection):end(selection)]))
  mcols(gr)<-sweep(as.matrix(mcols(gr)),MARGIN=2,weights,"/")
  return(gr)
}

paddedLog2<-function(x) { if(logT) return(log2(x+0.1)) else(return(x)) }
#myImportFun("data/DG4700-1.Aligned.out.sort.dedup.unique.bam",GRanges(seqnames=gene_chr,IRanges(start=gene_start,end=gene_end),strand="*"))
gtrack <- GenomeAxisTrack(col="darkgray")
txTr <- GeneRegionTrack(txdb, chromosome =  gene_chr, start = gene_start-1, end = gene_end+1,fill="gray",col="black",fontcolor.group="black",fill="darkgray",name="WS273 Transcripts")
dT<-DataTrack(range="data/DG4700-1.Aligned.out.sort.dedup.unique.bam", genome="ce11", type="p", name="Log2 Coverage", chromosome=gene_chr,importFunction = myImportFun,stream=T,col=cbPalette[c(4,6,7,2)])

group_factor=factor(c(rep("Soma Control",2),rep("Germline Control",3),rep("Germline Deplete",3), rep("Soma Deplete",3)),levels=c("Soma Control","Germline Control","Germline Deplete","Soma Deplete"))


plotTracks(list(gtrack,dT,txTr), from=gene_start-pad[1], to = gene_end+pad[2],transformation=paddedLog2,reverseStrand=rev,
           cex=0.8,add53=TRUE,ylim=c(0,ylim),type=c("a","confint"),
           groups=group_factor,
           col.axis="black",background.title="transparent",fontcolor.title="black") 
}


plotTranscriptSO<-function(x="fog-1",ylim=500,pad=c(100,500),rev=FALSE,logT=FALSE) {

if (stringr::str_detect(x,":")) {
    x<-gsub("-",":",x) %>% gsub(",","",.)
    gene_chr<-sapply(strsplit(x,":"),function(y) y[1]) %>% gsub("chr","",.)
    gene_start<-as.numeric(sapply(strsplit(x,":"),function(y) y[2]))
    gene_end<-as.numeric(sapply(strsplit(x,":"),function(y) y[3]))
  } else {
    gene<-unique(symbols[grep(common_name,symbols$gene_name),"gene_id"])
    gene_start<-start(transcripts[transcripts$gene_id==gene])
    if (length(gene_start) > 1) { gene_start<-min(gene_start) }
    gene_end<-end(transcripts[transcripts$gene_id==gene])
    if (length(gene_end) > 1) { gene_end<-max(gene_end) }
    gene_chr<-unique(as.character(seqnames(transcripts[transcripts$gene_id==gene])))
  }

myImportFunSO <- function(file, selection){
  (fls<-list.files("../data","*.bam$",full.names=T) %>% stringr::str_subset("ortiz"))
  s_names<-sapply(strsplit(basename(fls),"_"),function(x) x[5]) %>% gsub("\\.q4\\.bam","",.)
  weights<-colData(cds2[,s_names])$sizeFactor
  group<-as.character(colData(cds2[,s_names])$group)
  param_total <- ScanBamParam(what="mapq",flag = scanBamFlag(isUnmappedQuery = FALSE,isMinusStrand = NA,isMateMinusStrand = NA,isProperPair = FALSE),which=selection)
  gr<-GRanges(seqnames=seqnames(selection),IRanges(start=start(selection):end(selection),width=1),strand="*")
  mcols(gr)<-do.call(cbind,lapply(fls,function(file) as.numeric(coverage(suppressWarnings(GenomicAlignments::readGAlignments(file,use.names=F,param=param_total)))[[as.character(seqnames(selection))]])[start(selection):end(selection)]))
  mcols(gr)<-sweep(as.matrix(mcols(gr)),MARGIN=2,weights,"/")
  colnames(mcols(gr))<-group
  return(gr)
}

paddedLog2<-function(x) { if(logT) return(log2(x+0.1)) else(return(x)) }
#x<-myImportFunSO("data/DG4700-1.Aligned.out.sort.dedup.unique.bam",GRanges(seqnames=gene_chr,IRanges(start=gene_start,end=gene_end),strand="*"))
gtrack <- GenomeAxisTrack(col="darkgray")
txTr <- GeneRegionTrack(txdb, chromosome =  gene_chr, start = gene_start-1, end = gene_end+1,fill="gray",col="black",fontcolor.group="black",fill="darkgray",name="WS273 Transcripts")
dT<-DataTrack(range="data/DG4700-1.Aligned.out.sort.dedup.unique.bam", genome="ce11", type="p", name="Log2 Coverage", chromosome=gene_chr,importFunction = myImportFunSO,stream=T,col=cbPalette[c(3,8)])
plotTracks(list(gtrack,dT,txTr), from=gene_start-pad[1], to = gene_end+pad[2],transformation=paddedLog2,reverseStrand=rev,
           cex=0.8,add53=TRUE,ylim=c(0,ylim),type=c("a","confint"),
           groups=factor(rep(c("Sperm [fem-3(q96)]","Oocytes [fog-2(q71)]"), each=8),levels=c("Sperm [fem-3(q96)]","Oocytes [fog-2(q71)]")),
           col.axis="black",background.title="transparent",fontcolor.title="black") 

}

#plotTranscriptSO("fog-1",ylim=10,pad=c(100,10),logT=T,rev=F)
#plotTranscriptSO("chrI:3,210,365-3,212,646",ylim=10,pad=c(100,10),logT=T,rev=F)
#plotTranscript("fog-1",ylim=10,pad=c(100,10),logT=T,rev=F)
#plotTranscript("chrI:3,210,365-3,212,646",ylim=10,pad=c(100,10),logT=T,rev=F)

#svglite::svglite(file=paste0("2019 Manuscript/ptc1_log2_",ts,".svg"),width=5,height=3) 
plotTranscript("chrII:7,893,466-7,893,573",ylim=10,pad=c(0,0),logT=T,rev=F)
```

![](sacy1_polyA_files/figure-markdown_github/transcript_plots_sperm_oocyte-1.png)

``` r
#dev.off()
#svglite::svglite(file=paste0("2019 Manuscript/ptc1_log2_kimble_",ts,".svg"),width=5,height=3) 
plotTranscriptSO("chrII:7,893,466-7,893,573",ylim=10,pad=c(0,0),logT=T,rev=F)
```

![](sacy1_polyA_files/figure-markdown_github/transcript_plots_sperm_oocyte-2.png)

``` r
#dev.off()

#lin-54
#svglite::svglite(file=paste0("2019 Manuscript/lin54_log2_",ts,".svg"),width=5,height=3) 
plotTranscript("chrIV:13,241,150-13,241,500",ylim=10,pad=c(0,0),logT=T,rev=F)
```

![](sacy1_polyA_files/figure-markdown_github/transcript_plots_sperm_oocyte-3.png)

``` r
#dev.off()
#svglite::svglite(file=paste0("2019 Manuscript/lin54_log2_kimble_",ts,".svg"),width=5,height=3) 
plotTranscriptSO("chrIV:13,241,150-13,241,500",ylim=10,pad=c(0,0),logT=T,rev=F)
```

![](sacy1_polyA_files/figure-markdown_github/transcript_plots_sperm_oocyte-4.png)

``` r
#dev.off()

#wdr-4
#svglite::svglite(file=paste0("2019 Manuscript/wdr4_log2_",ts,".svg"),width=5,height=3) 
plotTranscript("chrIII:6,726,380-6,726,559",ylim=10,pad=c(0,0),logT=T,rev=F)
```

![](sacy1_polyA_files/figure-markdown_github/transcript_plots_sperm_oocyte-5.png)

``` r
#dev.off()
#svglite::svglite(file=paste0("2019 Manuscript/wdr4_log2_kimble_",ts,".svg"),width=5,height=3) 
plotTranscriptSO("chrIII:6,726,380-6,726,559",ylim=10,pad=c(0,0),logT=T,rev=F)
```

![](sacy1_polyA_files/figure-markdown_github/transcript_plots_sperm_oocyte-6.png)

``` r
#dev.off()

#C17H12.2
#svglite::svglite(file=paste0("2019 Manuscript/C17H12_log2_",ts,".svg"),width=5,height=3) 
plotTranscript("chrIV:6,788,665-6,788,810",ylim=10,pad=c(0,0),logT=T,rev=T)
```

![](sacy1_polyA_files/figure-markdown_github/transcript_plots_sperm_oocyte-7.png)

``` r
#dev.off()
#svglite::svglite(file=paste0("2019 Manuscript/C17H12_log2_kimble_",ts,".svg"),width=5,height=3) 
plotTranscriptSO("chrIV:6,788,665-6,788,810",ylim=10,pad=c(0,0),logT=T,rev=T)
```

![](sacy1_polyA_files/figure-markdown_github/transcript_plots_sperm_oocyte-8.png)

``` r
#dev.off()

#tra-2
plotTranscript("chrII:6,955,397-6,956,404",ylim=10,pad=c(0,0),logT=T,rev=T)
```

![](sacy1_polyA_files/figure-markdown_github/transcript_plots_sperm_oocyte-9.png)

``` r
plotTranscriptSO("chrII:6,955,397-6,956,404",ylim=10,pad=c(0,0),logT=T,rev=T)
```

![](sacy1_polyA_files/figure-markdown_github/transcript_plots_sperm_oocyte-10.png)

Figures - GO
============

``` r
#download.file("ftp://ftp.wormbase.org/pub/wormbase/releases/WS273/ONTOLOGY/gene_association.WS273.wb", destfile="gene_association.WS273.wb")
#download.file("ftp://ftp.wormbase.org/pub/wormbase/releases/WS273/ONTOLOGY/gene_ontology.WS273.obo",destfile= "gene_ontology.WS273.obo")

ga<-readr::read_delim("gene_association.WS273.wb",skip=3,delim="\t",col_names=F)
#Scan in OBO file to get Biological Process terms
t1<-readLines("gene_ontology.WS273.obo",n=-1L)
terms<-list()
for (i in seq_along(1:length(t1))) {
  if (t1[i]=="[Term]") {
    gt<-substr(t1[i+1],5,14)
    name<-gsub("name: ","",t1[i+2])
    if (t1[i+3]=="namespace: biological_process") { 
      #print (paste0(gt,":  ",name)) 
      terms[[gt]]<-name
    }
    }
}
go_terms<-unlist(terms)
head(ga)
wb_go<-ga[ga$'X5' %in% names(go_terms),c('X2','X5')]
colnames(wb_go)<-c("ensembl","GO ID")
dim(wb_go<-as.data.frame(wb_go))
#remove duplicated entries (possibly from multiple lines of evidence)
dim(wb_go<-wb_go[!duplicated(wb_go),])
wb_go.list<-split(wb_go$`GO ID`,wb_go$ensembl)
#List Terms for a gene
symbols[grep("sacy-1",symbols$gene_name),]
go_terms[wb_go.list[["WBGene00019245"]]]
symbols[grep("lin-41",symbols$gene_name),]
go_terms[wb_go.list[["WBGene00003026"]]]
symbols[grep("lin-29",symbols$gene_name),]
go_terms[wb_go.list[["WBGene00003015"]]]
#List Genes for a Term
symbols[symbols$gene_id %in% wb_go[grep("GO:0045138",wb_go$`GO ID`),"ensembl"],]

go_terms["GO:0008150"] #bp category
wb_go<-wb_go[wb_go$`GO ID`!="GO:0008150",]

save(wb_go,go_terms,file="wb_273_goterms.rdata")
```

GO Analysis
-----------

``` r
load("wb_273_goterms.rdata")
length(expressed_genes<-rownames(res_soma[res_soma$baseMean > 1,]))
```

    ## [1] 21248

``` r
dim(wb_go<-wb_go[wb_go$ensembl %in% expressed_genes,])
```

    ## [1] 31808     2

``` r
wb_go.list<-split(wb_go$`GO ID`,wb_go$ensembl)

#bias data
bd<-mcols(cds)$basepairs
names(bd)<-rownames(cds)
bd<-bd[names(bd) %in% expressed_genes]

table(degs_germline_up<-as.numeric(expressed_genes %in% germline_up ))
```

    ## 
    ##     0     1 
    ## 20937   311

``` r
names(degs_germline_up)<-expressed_genes
germline_up_GO<-goseq(nullp(degs_germline_up,bias.data=bd),gene2cat=wb_go.list)
```

    ## Using manually entered categories.

    ## For 12333 genes, we could not find any categories. These genes will be excluded.

    ## To force their use, please run with use_genes_without_cat=TRUE (see documentation).

    ## This was the default behavior for version 1.15.1 and earlier.

    ## Calculating the p-values...

    ## 'select()' returned 1:1 mapping between keys and columns

![](sacy1_polyA_files/figure-markdown_github/GO-1.png)

``` r
germline_up_GO[1:10,c("category","term","numDEInCat","numInCat","over_represented_pvalue")]
```

    ##        category
    ## 1543 GO:0031146
    ## 455  GO:0006511
    ## 1238 GO:0016567
    ## 2829 GO:0061158
    ## 2624 GO:0051260
    ## 778  GO:0007369
    ## 767  GO:0007275
    ## 3345 GO:1902742
    ## 2674 GO:0051574
    ## 363  GO:0006357
    ##                                                                         term
    ## 1543 SCF-dependent proteasomal ubiquitin-dependent protein catabolic process
    ## 455                            ubiquitin-dependent protein catabolic process
    ## 1238                                                  protein ubiquitination
    ## 2829                                    3'-UTR-mediated mRNA destabilization
    ## 2624                                             protein homooligomerization
    ## 778                                                             gastrulation
    ## 767                                       multicellular organism development
    ## 3345                               apoptotic process involved in development
    ## 2674                        positive regulation of histone H3-K9 methylation
    ## 363                         regulation of transcription by RNA polymerase II
    ##      numDEInCat numInCat over_represented_pvalue
    ## 1543          8       36            5.543578e-10
    ## 455           8      117            1.696841e-06
    ## 1238          8      163            1.053716e-05
    ## 2829          4       19            1.515674e-05
    ## 2624          5       52            1.308216e-04
    ## 778           3       22            1.899177e-04
    ## 767           8      429            9.750197e-04
    ## 3345          2       12            1.840957e-03
    ## 2674          1        1            3.461719e-03
    ## 363           6      210            4.682993e-03

``` r
listGO("GO:0061158",degs_germline_up)
```

    ##              ensembl         go     symbol deg
    ## 26217 WBGene00009537 GO:0061158     ccch-2   1
    ## 31249 WBGene00013794 GO:0061158     dct-13   1
    ## 31251 WBGene00013796 GO:0061158 Y116A8C.19   1
    ## 31252 WBGene00013797 GO:0061158 Y116A8C.20   1

``` r
gg_plotcounts("ccch-2")
```

![](sacy1_polyA_files/figure-markdown_github/GO-2.png)

``` r
listGO("GO:0051574",degs_germline_up)
```

    ##             ensembl         go symbol deg
    ## 2873 WBGene00000499 GO:0051574  chk-2   1

``` r
gg_plotcounts("chk-2")
```

![](sacy1_polyA_files/figure-markdown_github/GO-3.png)

``` r
listGO("GO:0007369",degs_germline_up)
```

    ##              ensembl         go symbol deg
    ## 6411  WBGene00001311 GO:0007369  end-3   1
    ## 28020 WBGene00011124 GO:0007369 sdz-27   1
    ## 37451 WBGene00020073 GO:0007369 sdz-28   1

``` r
table(degs_germline_dn<-as.numeric(expressed_genes %in% germline_dn))
```

    ## 
    ##     0     1 
    ## 21122   126

``` r
names(degs_germline_dn)<-expressed_genes
germline_dn_GO<-goseq(nullp(degs_germline_dn,bias.data=bd),gene2cat=wb_go.list)
```

    ## Using manually entered categories.

    ## For 12333 genes, we could not find any categories. These genes will be excluded.

    ## To force their use, please run with use_genes_without_cat=TRUE (see documentation).

    ## This was the default behavior for version 1.15.1 and earlier.

    ## Calculating the p-values...

    ## 'select()' returned 1:1 mapping between keys and columns

![](sacy1_polyA_files/figure-markdown_github/GO-4.png)

``` r
germline_dn_GO[1:10,c("category","term","numDEInCat","numInCat","over_represented_pvalue")]
```

    ##        category
    ## 2574 GO:0050976
    ## 948  GO:0009408
    ## 833  GO:0007638
    ## 1907 GO:0036498
    ## 966  GO:0009612
    ## 1511 GO:0030968
    ## 3186 GO:0098542
    ## 2902 GO:0070588
    ## 1785 GO:0034620
    ## 3538 GO:1905789
    ##                                                                                                 term
    ## 2574                        detection of mechanical stimulus involved in sensory perception of touch
    ## 948                                                                                 response to heat
    ## 833                                                                          mechanosensory behavior
    ## 1907                                                         IRE1-mediated unfolded protein response
    ## 966                                                                  response to mechanical stimulus
    ## 1511                                                 endoplasmic reticulum unfolded protein response
    ## 3186                                                              defense response to other organism
    ## 2902                                                             calcium ion transmembrane transport
    ## 1785                                                           cellular response to unfolded protein
    ## 3538 positive regulation of detection of mechanical stimulus involved in sensory perception of touch
    ##      numDEInCat numInCat over_represented_pvalue
    ## 2574          4       12            2.202977e-07
    ## 948           5       46            2.220603e-06
    ## 833           4       28            9.617514e-06
    ## 1907          6      117            2.527469e-05
    ## 966           3       18            7.198375e-05
    ## 1511          4       53            1.417511e-04
    ## 3186          3       20            1.718241e-04
    ## 2902          3       35            4.708037e-04
    ## 1785          2        9            7.536101e-04
    ## 3538          2       10            1.027480e-03

``` r
listGO("GO:0031146",degs_germline_dn)
```

    ## [1] ensembl go      symbol  deg    
    ## <0 rows> (or 0-length row.names)

``` r
listGO("GO:0036498",degs_germline_dn)
```

    ##             ensembl         go    symbol deg
    ## 9256 WBGene00002015 GO:0036498  hsp-16.1   1
    ## 9261 WBGene00002017 GO:0036498 hsp-16.11   1
    ## 9264 WBGene00002018 GO:0036498 hsp-16.41   1
    ## 9266 WBGene00002019 GO:0036498 hsp-16.48   1
    ## 9271 WBGene00002020 GO:0036498 hsp-16.49   1
    ## 9289 WBGene00002026 GO:0036498    hsp-70   1

``` r
table(degs_soma_up<-as.numeric(expressed_genes %in% soma_up ))
```

    ## 
    ##     0     1 
    ## 21006   242

``` r
names(degs_soma_up)<-expressed_genes
soma_up_GO<-goseq(nullp(degs_soma_up,bias.data=bd),gene2cat=wb_go.list)
```

    ## Using manually entered categories.

    ## For 12333 genes, we could not find any categories. These genes will be excluded.

    ## To force their use, please run with use_genes_without_cat=TRUE (see documentation).

    ## This was the default behavior for version 1.15.1 and earlier.

    ## Calculating the p-values...

    ## 'select()' returned 1:1 mapping between keys and columns

![](sacy1_polyA_files/figure-markdown_github/GO-5.png)

``` r
soma_up_GO[1:10,c("category","term","numDEInCat","numInCat","over_represented_pvalue")]
```

    ##        category                                            term numDEInCat
    ## 948  GO:0009408                                response to heat          6
    ## 1785 GO:0034620           cellular response to unfolded protein          3
    ## 1779 GO:0034605                       cellular response to heat          3
    ## 1907 GO:0036498         IRE1-mediated unfolded protein response          7
    ## 639  GO:0006986                    response to unfolded protein          3
    ## 2590 GO:0051085  chaperone cofactor-dependent protein refolding          3
    ## 2208 GO:0045087                          innate immune response          9
    ## 1954 GO:0042026                               protein refolding          3
    ## 1511 GO:0030968 endoplasmic reticulum unfolded protein response          4
    ## 249  GO:0006082                  organic acid metabolic process          3
    ##      numInCat over_represented_pvalue
    ## 948        46            7.674518e-07
    ## 1785        9            1.180552e-05
    ## 1779       11            2.249635e-05
    ## 1907      117            6.226842e-05
    ## 639        16            8.860262e-05
    ## 2590       17            2.552244e-04
    ## 2208      224            3.694736e-04
    ## 1954       21            8.322296e-04
    ## 1511       53            9.984178e-04
    ## 249        40            2.842254e-03

``` r
listGO("GO:0009408",degs_soma_up)
```

    ##              ensembl         go    symbol deg
    ## 9257  WBGene00002016 GO:0009408  hsp-16.2   1
    ## 9258  WBGene00002017 GO:0009408 hsp-16.11   1
    ## 9262  WBGene00002018 GO:0009408 hsp-16.41   1
    ## 9268  WBGene00002020 GO:0009408 hsp-16.49   1
    ## 9284  WBGene00002026 GO:0009408    hsp-70   1
    ## 13163 WBGene00003474 GO:0009408     mtl-2   1

``` r
listGO("GO:0036498",degs_soma_up)
```

    ##              ensembl         go    symbol deg
    ## 9256  WBGene00002015 GO:0036498  hsp-16.1   1
    ## 9261  WBGene00002017 GO:0036498 hsp-16.11   1
    ## 9264  WBGene00002018 GO:0036498 hsp-16.41   1
    ## 9266  WBGene00002019 GO:0036498 hsp-16.48   1
    ## 9271  WBGene00002020 GO:0036498 hsp-16.49   1
    ## 9289  WBGene00002026 GO:0036498    hsp-70   1
    ## 24550 WBGene00007992 GO:0036498   fipr-24   1

``` r
table(degs_soma_dn<-as.numeric(expressed_genes %in% soma_dn ))
```

    ## 
    ##     0     1 
    ## 21006   242

``` r
names(degs_soma_dn)<-expressed_genes
soma_dn_GO<-goseq(nullp(degs_soma_dn,bias.data=bd),gene2cat=wb_go.list)
```

    ## Using manually entered categories.

    ## For 12333 genes, we could not find any categories. These genes will be excluded.

    ## To force their use, please run with use_genes_without_cat=TRUE (see documentation).

    ## This was the default behavior for version 1.15.1 and earlier.

    ## Calculating the p-values...

    ## 'select()' returned 1:1 mapping between keys and columns

![](sacy1_polyA_files/figure-markdown_github/GO-6.png)

``` r
soma_dn_GO[1:10,c("category","term","numDEInCat","numInCat","over_represented_pvalue")]
```

    ##        category
    ## 1991 GO:0042338
    ## 1392 GO:0019934
    ## 290  GO:0006182
    ## 916  GO:0009190
    ## 596  GO:0006869
    ## 2208 GO:0045087
    ## 1332 GO:0019346
    ## 1851 GO:0035434
    ## 2890 GO:0070482
    ## 527  GO:0006644
    ##                                                                                    term
    ## 1991 cuticle development involved in collagen and cuticulin-based cuticle molting cycle
    ## 1392                                                            cGMP-mediated signaling
    ## 290                                                           cGMP biosynthetic process
    ## 916                                              cyclic nucleotide biosynthetic process
    ## 596                                                                     lipid transport
    ## 2208                                                             innate immune response
    ## 1332                                                                   transsulfuration
    ## 1851                                                 copper ion transmembrane transport
    ## 2890                                                          response to oxygen levels
    ## 527                                                      phospholipid metabolic process
    ##      numDEInCat numInCat over_represented_pvalue
    ## 1991          9       34            4.235427e-11
    ## 1392          4        7            1.286810e-07
    ## 290           4       32            6.955488e-05
    ## 916           4       36            1.085473e-04
    ## 596           3       34            2.504732e-03
    ## 2208          8      224            3.512709e-03
    ## 1332          2       10            4.873109e-03
    ## 1851          2        9            5.211530e-03
    ## 2890          2       14            5.947091e-03
    ## 527           3       34            6.163753e-03

``` r
listGO("GO:0042338",degs_soma_dn)
```

    ##              ensembl         go symbol deg
    ## 5141  WBGene00001065 GO:0042338  dpy-3   1
    ## 5145  WBGene00001066 GO:0042338  dpy-4   1
    ## 5149  WBGene00001067 GO:0042338  dpy-5   1
    ## 5154  WBGene00001070 GO:0042338  dpy-8   1
    ## 5167  WBGene00001074 GO:0042338 dpy-13   1
    ## 16896 WBGene00004397 GO:0042338  rol-6   1
    ## 18825 WBGene00005016 GO:0042338  sqt-1   1
    ## 18826 WBGene00005017 GO:0042338  sqt-2   1
    ## 40309 WBGene00022743 GO:0042338  mlt-7   1

``` r
listGO("GO:0070482",degs_soma_dn)
```

    ##             ensembl         go symbol deg
    ## 7609 WBGene00001555 GO:0070482 gcy-35   1
    ## 7618 WBGene00001556 GO:0070482 gcy-36   1

``` r
g5<-head(soma_dn_GO,10) %>% 
  dplyr::mutate(term=glue::glue("{term} ({category})")) %>% 
  dplyr::mutate(term=factor(term,levels=rev(term))) %>%
  ggplot(aes(x=term,y=-log10(over_represented_pvalue))) +
  geom_point(color=cbPalette[4],aes(size=numDEInCat)) +
  scale_radius(limits=c(1,12)) +
 # geom_bar(stat="identity",fill=cbPalette[8]) + ylim(c(0,10)) +
  coord_flip() + xlab("") + ggtitle("Soma DN GO Terms") + ylim(0,12) +
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),aspect.ratio=4/3)

g6<-head(soma_up_GO,10) %>% 
  dplyr::mutate(term=glue::glue("{term} ({category})")) %>% 
  dplyr::mutate(term=factor(term,levels=rev(term))) %>%
  ggplot(aes(x=term,y=-log10(over_represented_pvalue))) +
  geom_point(color=cbPalette[2],aes(size=numDEInCat)) +
   scale_radius(limits=c(1,12)) +
 # geom_bar(stat="identity",fill=cbPalette[8]) + ylim(c(0,10)) +
  coord_flip() + xlab("") + ggtitle("Soma UP GO Terms") + ylim(0,12) +
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),aspect.ratio=4/3)

g7<-head(germline_dn_GO,10) %>% 
  dplyr::mutate(term=glue::glue("{term} ({category})")) %>% 
  dplyr::mutate(term=factor(term,levels=rev(term))) %>%
  ggplot(aes(x=term,y=-log10(over_represented_pvalue))) +
  geom_point(color=cbPalette[6],aes(size=numDEInCat)) +
   scale_radius(limits=c(1,12)) +
 # geom_bar(stat="identity",fill=cbPalette[8]) + ylim(c(0,10)) +
  coord_flip() + xlab("") + ggtitle("Germline DN GO Terms") + ylim(0,12) +
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),aspect.ratio=4/3)

g8<-head(germline_up_GO,10) %>% 
  dplyr::mutate(term=glue::glue("{term} ({category})")) %>% 
  dplyr::mutate(term=factor(term,levels=rev(term))) %>%
  ggplot(aes(x=term,y=-log10(over_represented_pvalue))) +
  geom_point(color=cbPalette[7],aes(size=numDEInCat)) +
   scale_radius(limits=c(1,12)) +
 # geom_bar(stat="identity",fill=cbPalette[8]) + ylim(c(0,10)) +
  coord_flip() + xlab("") + ggtitle("Germline UP GO Terms") + ylim(0,12) +
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),aspect.ratio=4/3)

#pdf(paste0("GO_Overlap_",ts,".pdf"),width=3,height=3);
#gridExtra::grid.arrange(g5,g6,g7,g8,ncol=1)
#dev.off()
g5
```

![](sacy1_polyA_files/figure-markdown_github/GO-7.png)

``` r
g6
```

![](sacy1_polyA_files/figure-markdown_github/GO-8.png)

``` r
g7
```

![](sacy1_polyA_files/figure-markdown_github/GO-9.png)

``` r
g8
```

![](sacy1_polyA_files/figure-markdown_github/GO-10.png)

``` r
#svglite::svglite(file=paste0("2019 Manuscript/GO_soma_dn_",ts,".svg"),width=10,height=3); g5; dev.off()
#svglite::svglite(file=paste0("2019 Manuscript/GO_soma_up_",ts,".svg"),width=10,height=3); g6; dev.off()
#svglite::svglite(file=paste0("2019 Manuscript/GO_germ_dn_",ts,".svg"),width=10,height=3); g7; dev.off()
#svglite::svglite(file=paste0("2019 Manuscript/GO_germ_up_",ts,".svg"),width=10,height=3); g8; dev.off()
```

\#Remake volcano plots with some of the major GO category genes plotted

``` r
df_germ<-res_germ[res_germ$gene_name %in% c("Y54G2A.36","her-1","hsp-16.41","hsp-70","hsp-16.1","sdz-27","skr-8","xol-1"),]
(df_germ<-subset(df_germ,abs(log2FoldChange) > 1 & padj < 1e-10))
```

    ##                 baseMean log2FoldChange     lfcSE       stat        pvalue
    ## WBGene00002018 3946.0711      -1.797509 0.1107486 -16.230535  3.067750e-59
    ## WBGene00002015  765.5775      -1.637007 0.1278731 -12.801807  1.601794e-37
    ## WBGene00001842  144.5161       2.048695 0.2529907   8.097908  5.591215e-16
    ## WBGene00021898 1961.4116       2.405268 0.1054269  22.814552 3.288027e-115
    ## WBGene00011124 4641.1750       1.079828 0.1007539  10.717481  8.426667e-27
    ## WBGene00004814  763.2603       1.932052 0.2045432   9.445692  3.530630e-21
    ## WBGene00002026 2057.9675      -1.838308 0.1335676 -13.763129  4.247377e-43
    ##                         padj germline_control germline_experimental
    ## WBGene00002018  2.661426e-55          493.773               141.999
    ## WBGene00002015  4.632121e-34          128.666                41.316
    ## WBGene00001842  1.347405e-13            4.102                16.995
    ## WBGene00021898 5.705055e-111           36.246               192.161
    ## WBGene00011124  1.207203e-23          123.376               260.849
    ## WBGene00004814  2.268887e-18           13.165                50.239
    ## WBGene00002026  1.842406e-39           72.537                20.282
    ##                soma_control soma_experimental gene_name   gene_biotype
    ## WBGene00002018       47.013           144.831 hsp-16.41 protein_coding
    ## WBGene00002015       19.675            47.076  hsp-16.1 protein_coding
    ## WBGene00001842        9.228             8.794     her-1 protein_coding
    ## WBGene00021898       43.044            52.228 Y54G2A.36 protein_coding
    ## WBGene00011124      252.792           294.932    sdz-27 protein_coding
    ## WBGene00004814       57.337            49.516     skr-8 protein_coding
    ## WBGene00002026       12.197            29.844    hsp-70 protein_coding

``` r
#svglite::svglite(file=paste0("2019 Manuscript/res_germ_volcano_",ts,".svg"),width=5,height=4)
res_germ %>% 
  tibble::rownames_to_column(var = "gene_id") %>%
  dplyr::filter(!(is.na(padj) | is.na(log2FoldChange)) & baseMean > 25) %>% 
  mutate(neglog10=-1*log10(padj)) %>% 
  mutate(cat=case_when(
    gene_id %in% germline_up ~ "UP",
    gene_id %in% germline_dn ~ "DN",
    TRUE ~ "not_sig"  
  )) %>%
  ggplot(aes(x=log2FoldChange,y=neglog10,color=cat)) + 
  scale_color_manual(values=cbPalette[c(6,1,7)]) +
  geom_point() + ylim(c(0,150)) + xlim(c(-5,5)) +  theme_bw() + 
  annotate("text",x=df_germ$log2FoldChange,y=(-1*log10(df_germ$padj)),label=df_germ$gene_name) +
  ggtitle("Gene Expression Changes in Germline Deplete") + ylab("- Log10 Adjusted P-value") + xlab("Log2 Fold Change")
```

![](sacy1_polyA_files/figure-markdown_github/volcano_plots_with_GO_highlights-1.png)

``` r
#dev.off()

df_soma<-res_soma[res_soma$gene_name %in% c("hsp-16.41","hsp-70","hsp-16.1","col-179","pgp-8","gst-24","dpy-9","dpy-8","dpy-20","sqt-1","rol-6"),]
(df_soma<-subset(df_soma,abs(log2FoldChange) > 1 & padj < 1e-10))
```

    ##                 baseMean log2FoldChange     lfcSE       stat        pvalue
    ## WBGene00002018 3946.0711       1.622787 0.1277816  12.699694  5.936118e-37
    ## WBGene00002015  765.5775       1.256667 0.1531833   8.203682  2.331336e-16
    ## WBGene00001070  746.4521      -4.150061 0.2973960 -13.954666  2.946793e-44
    ## WBGene00004002  263.7100       4.470347 0.2414085  18.517765  1.484740e-76
    ## WBGene00000752 5657.6893       1.576601 0.1563737  10.082268  6.618053e-24
    ## WBGene00001071 1191.4255      -4.115705 0.2193440 -18.763701  1.496111e-78
    ## WBGene00001079  165.1514      -2.904942 0.3459904  -8.396019  4.618716e-17
    ## WBGene00004397 1505.5258      -3.719881 0.2396206 -15.524047  2.385198e-54
    ## WBGene00005016 1297.0857      -3.660847 0.3834606  -9.546868  1.336770e-21
    ## WBGene00001772 1174.8487       3.720493 0.1560138  23.847200 1.082701e-125
    ## WBGene00002026 2057.9675       1.291487 0.1521479   8.488371  2.095544e-17
    ##                         padj germline_control germline_experimental
    ## WBGene00002018  5.420926e-34          493.773               141.999
    ## WBGene00002015  3.816133e-14          128.666                41.316
    ## WBGene00001070  3.652129e-41           38.638                21.771
    ## WBGene00004002  5.152345e-73            0.432                 0.529
    ## WBGene00000752  2.944355e-21          133.754               225.402
    ## WBGene00001071  6.489753e-75           75.842                56.348
    ## WBGene00001079  8.261788e-15            5.856                 5.206
    ## WBGene00004397  5.173196e-51           99.762                57.269
    ## WBGene00005016  4.295239e-19           97.894                53.441
    ## WBGene00001772 1.878594e-121           11.302                13.236
    ## WBGene00002026  3.787477e-15           72.537                20.282
    ##                soma_control soma_experimental gene_name   gene_biotype
    ## WBGene00002018       47.013           144.831 hsp-16.41 protein_coding
    ## WBGene00002015       19.675            47.076  hsp-16.1 protein_coding
    ## WBGene00001070        6.312             0.352     dpy-8 protein_coding
    ## WBGene00004002        0.348             7.691     pgp-8 protein_coding
    ## WBGene00000752      111.636           332.924   col-179 protein_coding
    ## WBGene00001071       20.331             1.161     dpy-9 protein_coding
    ## WBGene00001079        1.861             0.247    dpy-20 protein_coding
    ## WBGene00004397       17.903             1.358     rol-6 protein_coding
    ## WBGene00005016       16.478             1.289     sqt-1 protein_coding
    ## WBGene00001772       13.695           180.615    gst-24 protein_coding
    ## WBGene00002026       12.197            29.844    hsp-70 protein_coding

``` r
#svglite::svglite(file=paste0("2019 Manuscript/res_soma_volcano_",ts,".svg"),width=5,height=4)
res_soma %>% 
  tibble::rownames_to_column(var = "gene_id") %>% 
  dplyr::filter(!(is.na(padj) | is.na(log2FoldChange)) & baseMean > 25) %>% 
  mutate(neglog10=-1*log10(padj)) %>% 
  mutate(cat=case_when(
    gene_id %in% soma_up ~ "UP",
    gene_id %in% soma_dn ~ "DN",
    TRUE ~ "not_sig"  
  )) %>% 
  ggplot(aes(x=log2FoldChange,y=neglog10,color=cat)) + 
  scale_color_manual(values=cbPalette[c(4,1,2)]) +
  geom_point() + ylim(c(0,150)) + xlim(c(-5,5)) +  theme_bw() + 
  annotate("text",x=df_soma$log2FoldChange,y=(-1*log10(df_soma$padj)),label=df_soma$gene_name) +
  ggtitle("Gene Expression Changes in somaline Deplete") + ylab("- Log10 Adjusted P-value") + xlab("Log2 Fold Change")
```

![](sacy1_polyA_files/figure-markdown_github/volcano_plots_with_GO_highlights-2.png)

``` r
#dev.off()
```

\#Heatmap from top 5 GO terms

``` r
nterms<-5

length(go_gois<- unique(c(do.call(rbind,lapply(soma_up_GO[1:nterms,c("category")],listGO,degs_soma_up))$ensembl,
                   do.call(rbind,lapply(soma_dn_GO[1:nterms,c("category")],listGO,degs_soma_dn))$ensembl,
                   do.call(rbind,lapply(germline_up_GO[1:nterms,c("category")],listGO,degs_germline_up))$ensembl,
                   do.call(rbind,lapply(germline_dn_GO[1:nterms,c("category")],listGO,degs_germline_dn))$ensembl)))
```

    ## [1] 49

``` r
temp<-as.matrix(t(scale(t(log2(f[go_gois,]+0.1)),center=T,scale=T)))
rownames(temp)<-symbols[match(rownames(temp),symbols$gene_id),"gene_name"]
col_fun = circlize::colorRamp2(c(min(temp),0, max(temp)), c("blue","white", "red"))

#svglite::svglite(file=paste0("2019 Manuscript/heatmap_",ts,".svg"),width=3.5,height=10)  

Heatmap(temp[,c(1,2,9,10,11,3:8)], column_title = "GO GOIs", rect_gp = gpar(col = "black", lwd = 0.5),
                  cluster_rows = TRUE, cluster_columns=FALSE,
                  row_names_side = "left",col = col_fun, show_heatmap_legend = TRUE)
```

![](sacy1_polyA_files/figure-markdown_github/supplemental_heatmap-1.png)

``` r
#dev.off()  
```

Figure 7B etr-1
===============

chrII:167,000-168,600
---------------------

``` r
options(ucscChromosomeNames=FALSE)
common_name<-"etr-1"
 gene_chr<-"II"
 gene_start<-167000
 gene_end<-168600
 logT<-TRUE
 
 
gene<-unique(symbols[grep(common_name,symbols$gene_name),"gene_id"])

myImportFun<- function(file, selection){
  (fls<-list.files("data","*.bam$",full.names=T) %>% stringr::str_subset("CA1200-1|CA1200-3|DG4703|CA1352|DG4700"))
  s_names<-sapply(strsplit(basename(fls),"\\."),function(x) x[1]) %>% gsub("-","_rep",.)
  weights<-colData(cds[,s_names])$sizeFactor
  param_total <- ScanBamParam(what="mapq",flag = scanBamFlag(isUnmappedQuery = FALSE,isMinusStrand = NA,isMateMinusStrand = NA,isProperPair = TRUE),which=selection)
  gr<-GRanges(seqnames=seqnames(selection),IRanges(start=start(selection):end(selection),width=1),strand="*")
  mcols(gr)<-do.call(cbind,lapply(fls,function(file) as.numeric(coverage(suppressWarnings(GenomicAlignments::readGAlignmentPairs(file,use.names=F,param=param_total)))[[as.character(seqnames(selection))]])[start(selection):end(selection)]))
  mcols(gr)<-sweep(as.matrix(mcols(gr)),MARGIN=2,weights,"/")
  return(gr)
}

#myImportFun("data/DG4700-1.Aligned.out.sort.dedup.unique.bam",GRanges(seqnames=gene_chr,IRanges(start=gene_start,end=gene_end),strand="*"))
paddedLog2<-function(x) { if(logT) return(log2(x+0.1)) else(return(x)) }
gtrack <- GenomeAxisTrack(col="darkgray")
txTr <- GeneRegionTrack(txdb, chromosome =  gene_chr, start = gene_start, end = gene_end,fill="gray",col="black",fontcolor.group="black",fill="darkgray",name="",transcriptAnnotation="transcript")
length(txTr@range <- txTr@range[stringr::str_detect("T01D1.2m.3|T01D1.2a.5|T01D1.2n.3|T01D1.2l.5",txTr@range$transcript)])
```

    ## [1] 56

``` r
#symbol(txTr)<-c("T01D1.2a.1","T01D1.2e.5")

dT<-DataTrack(range="data/DG4700-1.Aligned.out.sort.dedup.unique.bam", genome="ce11", type="p", name="Log2 Coverage", chromosome=gene_chr,importFunction = myImportFun,stream=T,col=cbPalette[c(4,6,7,2)])

group_factor=factor(c(rep("Soma Control",2),rep("Germline Control",3),rep("Germline Deplete",3), rep("Soma Deplete",3)),levels=c("Soma Control","Germline Control","Germline Deplete","Soma Deplete"))


#svglite::svglite(file=paste0("2019 Manuscript/etr1_log2_",ts,".svg"),width=8,height=3) 
plotTracks(list(gtrack,dT,txTr), from=gene_start, to=gene_end,transformation=paddedLog2,reverseStrand=FALSE,
           cex=0.8,add53=TRUE,ylim=c(0,10),type=c("a","confint"),
           groups=group_factor,
           col.axis="black",background.title="transparent",fontcolor.title="black") 
```

![](sacy1_polyA_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
#dev.off()
```

\#etr-1 zoom for insets

``` r
options(ucscChromosomeNames=FALSE)

common_name<-"etr-1"
 gene_start<-167880
 gene_end<-167920
 logT<-TRUE


gene<-unique(symbols[grep(common_name,symbols$gene_name),"gene_id"])
gene_chr<-unique(as.character(seqnames(transcripts[transcripts$gene_id==gene])))
```

    ## Warning in transcripts$gene_id == gene: longer object length is not a multiple
    ## of shorter object length

``` r
myImportFun <- function(file, selection){
  (fls<-list.files("data","*.bam$",full.names=T) %>% stringr::str_subset("CA1200-1|CA1200-3|DG4703|CA1352|DG4700"))
  s_names<-sapply(strsplit(basename(fls),"\\."),function(x) x[1]) %>% gsub("-","_rep",.)
  weights<-colData(cds[,s_names])$sizeFactor
  param_total <- ScanBamParam(what="mapq",flag = scanBamFlag(isUnmappedQuery = FALSE,isMinusStrand = NA,isMateMinusStrand = NA,isProperPair = TRUE),which=selection)
  gr<-GRanges(seqnames=seqnames(selection),IRanges(start=start(selection):end(selection),width=1),strand="*")
  mcols(gr)<-do.call(cbind,lapply(fls,function(file) as.numeric(coverage(suppressWarnings(GenomicAlignments::readGAlignmentPairs(file,use.names=F,param=param_total)))[[as.character(seqnames(selection))]])[start(selection):end(selection)]))
  mcols(gr)<-sweep(as.matrix(mcols(gr)),MARGIN=2,weights,"/")
  return(gr)
}


#myImportFun("data/DG4700-1.Aligned.out.sort.dedup.unique.bam",GRanges(seqnames=gene_chr,IRanges(start=gene_start,end=gene_end),strand="*"))
paddedLog2<-function(x) { if(logT) return(log2(x+0.1)) else(return(x)) }

sTrack <- SequenceTrack(ce11,complement=F)
dT<-DataTrack(range="data/DG4700-1.Aligned.out.sort.dedup.unique.bam", genome="ce11", type="p", name="Log2 Coverage", chromosome=gene_chr,importFunction = myImportFun,stream=T,col=cbPalette[c(4,6,7,2)])

group_factor=factor(c(rep("Soma Control",2),rep("Germline Control",3),rep("Germline Deplete",3), rep("Soma Deplete",3)),levels=c("Soma Control","Germline Control","Germline Deplete","Soma Deplete"))


#svglite::svglite(file=paste0("2019 Manuscript/etr1_zoom1_",ts,".svg"),width=3.5,height=3) 
plotTracks(list(dT,sTrack), from=gene_start, to = gene_end,transformation=paddedLog2,reverseStrand=FALSE,
           cex=0.4,add53=TRUE,ylim=c(0,10),type=c("a","confint"),
           groups=group_factor,
           col.axis="black",background.title="transparent",fontcolor.title="black") 
```

    ## Warning in plotTracks(list(dT, sTrack), from = gene_start, to = gene_end, : The
    ## track chromosomes in 'trackList' differ. Setting all tracks to chromosome 'II'

![](sacy1_polyA_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
#dev.off()

#svglite::svglite(file=paste0("2019 Manuscript/etr1_zoom2_",ts,".svg"),width=3.5,height=3) 
 gene_start<-168523
 gene_end<-168563
 plotTracks(list(dT,sTrack), from=gene_start, to = gene_end,transformation=paddedLog2,reverseStrand=FALSE,
           cex=0.4,add53=TRUE,ylim=c(0,10),type=c("a","confint"),
           groups=group_factor,
           col.axis="black",background.title="transparent",fontcolor.title="black") 
```

    ## Warning in plotTracks(list(dT, sTrack), from = gene_start, to = gene_end, : The
    ## track chromosomes in 'trackList' differ. Setting all tracks to chromosome 'II'

![](sacy1_polyA_files/figure-markdown_github/unnamed-chunk-3-2.png)

``` r
 #dev.off()
```

Figure 7C prdx-6 (soma RI)
==========================

Four condition plot for prdx-6
------------------------------

chrIV:219,273-221,018
---------------------

``` r
options(ucscChromosomeNames=FALSE)
common_name="prdx-6"
rev=FALSE
logT=TRUE
gene<-unique(symbols[grep(common_name,symbols$gene_name),"gene_id"])
gene_start<-219273
gene_end<-220990
pad<-c(0,0)
gene_chr<-unique(as.character(seqnames(transcripts[transcripts$gene_id==gene])))

myImportFun_prdx6 <- function(file, selection){
  (fls<-list.files("data","*.bam$",full.names=T) %>% stringr::str_subset("CA1200-1|CA1200-3|DG4703|CA1352|DG4700"))
  s_names<-sapply(strsplit(basename(fls),"\\."),function(x) x[1]) %>% gsub("-","_rep",.)
  weights<-colData(cds[,s_names])$sizeFactor
  param_total <- ScanBamParam(what="mapq",flag = scanBamFlag(isUnmappedQuery = FALSE,isMinusStrand = NA,isMateMinusStrand = NA,isProperPair = TRUE),which=selection)
  gr<-GRanges(seqnames=seqnames(selection),IRanges(start=start(selection):end(selection),width=1),strand="*")
  mcols(gr)<-do.call(cbind,lapply(fls,function(file) as.numeric(coverage(suppressWarnings(GenomicAlignments::readGAlignmentPairs(file,use.names=F,param=param_total)))[[as.character(seqnames(selection))]])[start(selection):end(selection)]))
  mcols(gr)<-sweep(as.matrix(mcols(gr)),MARGIN=2,weights,"/")
  return(gr)
}

#myImportFun_prdx6("data/DG4700-1.Aligned.out.sort.dedup.unique.bam",GRanges(seqnames=gene_chr,IRanges(start=gene_start,end=gene_end),strand="*"))
paddedLog2<-function(x) { if(logT) return(log2(x+0.1)) else(return(x)) }
gtrack <- GenomeAxisTrack(col="darkgray")
txTr <- GeneRegionTrack(txdb, chromosome =  gene_chr, start = gene_start-1, end = gene_end+1,fill="gray",col="black",fontcolor.group="black",fill="darkgray",name="WS273 Transcripts")
dT<-DataTrack(range="data/DG4700-1.Aligned.out.sort.dedup.unique.bam", genome="ce11", type="p", name="Log2 Coverage", chromosome=gene_chr,importFunction = myImportFun_prdx6,stream=T,col=cbPalette[c(4,6,7,2)])

group_factor=factor(c(rep("Soma Control",2),rep("Germline Control",3),rep("Germline Deplete",3), rep("Soma Deplete",3)),levels=c("Soma Control","Germline Control","Germline Deplete","Soma Deplete"))


#svglite::svglite(file=paste0("2019 Manuscript/prdx6_log2_",ts,".svg"),width=5,height=3) 
plotTracks(list(gtrack,dT,txTr), from=gene_start, to = gene_end,transformation=paddedLog2,reverseStrand=FALSE,
           cex=0.8,add53=TRUE,ylim=c(0,10),type=c("a","confint"),
           groups=group_factor,
           col.axis="black",background.title="transparent",fontcolor.title="black") 
```

![](sacy1_polyA_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
#dev.off()
```

Figure 7D lin-28 all reads
==========================

chrI:8,408,476-8,408,526

``` r
options(ucscChromosomeNames=FALSE)
common_name<-"lin-28"
 gene_start<-8408480
 gene_end<-8408520
 ylim<-1700
 logT<-FALSE
 
gene<-unique(symbols[grep(common_name,symbols$gene_name),"gene_id"])
gene_chr<-unique(as.character(seqnames(transcripts[transcripts$gene_id==gene])))

myImportFun_lin28b <- function(file, selection){
  (fls<-list.files("data","*.bam$",full.names=T) %>% stringr::str_subset("CA1200-1|CA1200-3|DG4703|CA1352|DG4700"))
  s_names<-sapply(strsplit(basename(fls),"\\."),function(x) x[1]) %>% gsub("-","_rep",.)
  weights<-colData(cds[,s_names])$sizeFactor
  param_total <- ScanBamParam(what="mapq",flag = scanBamFlag(isUnmappedQuery = FALSE,isMinusStrand = NA,isMateMinusStrand = NA,isProperPair = TRUE),which=selection)
  gr<-GRanges(seqnames=seqnames(selection),IRanges(start=start(selection):end(selection),width=1),strand="*")
  mcols(gr)<-do.call(cbind,lapply(fls,function(file) as.numeric(coverage(suppressWarnings(GenomicAlignments::readGAlignmentPairs(file,use.names=F,param=param_total)))[[as.character(seqnames(selection))]])[start(selection):end(selection)]))
  mcols(gr)<-sweep(as.matrix(mcols(gr)),MARGIN=2,weights,"/")
  return(gr)
}

#myImportFun_lin28b("data/DG4700-1.Aligned.out.sort.dedup.unique.bam",GRanges(seqnames=gene_chr,IRanges(start=gene_start,end=gene_end),strand="*"))
paddedLog2<-function(x) {
  if(logT) return(log2(x+0.1))
  else(return(x))
  }
gtrack <- GenomeAxisTrack(col="darkgray")
sTrack <- SequenceTrack(ce11,complement=T)
txTr <- GeneRegionTrack(txdb, chromosome =  gene_chr, start = gene_start-1, end = gene_end+1,fill="gray",col="black",fontcolor.group="black",fill="darkgray",name="",transcriptAnnotation="transcript")
#length(txTr@range <- txTr@range[stringr::str_detect("T01D1.2m.3|T01D1.2a.5|T01D1.2n.3|T01D1.2l.5",txTr@range$transcript)])


dT<-DataTrack(range="data/DG4700-1.Aligned.out.sort.dedup.unique.bam", genome="ce11", type="p", name="Log2 Coverage", chromosome=gene_chr,importFunction = myImportFun_lin28b,stream=T,col=cbPalette[c(4,6,7,2)])

group_factor=factor(c(rep("Soma Control",2),rep("Germline Control",3),rep("Germline Deplete",3), rep("Soma Deplete",3)),levels=c("Soma Control","Germline Control","Germline Deplete","Soma Deplete"))


#svglite::svglite(file=paste0("2019 Manuscript/lin28_full_log2_",ts,".svg"),width=3.5,height=3) 
plotTracks(list(gtrack,dT,sTrack,txTr), from=gene_start, to = gene_end,transformation=paddedLog2,reverseStrand=TRUE,
           cex=0.4,add53=TRUE,ylim=c(0,1700),type=c("a","confint"),
           groups=group_factor,
           col.axis="black",background.title="transparent",fontcolor.title="black") 
```

![](sacy1_polyA_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
#dev.off()
```
