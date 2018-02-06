#!/usr/bin/env littler
# 三组或更多组的画法，只需添加数据即可
# https://mp.weixin.qq.com/s/FZWinr14RTs6YSUE_juaug
# Galaxy实测，两组也行。

library(tidyverse)

df <- data.frame(
  Phylum=c("Ruminococcaceae","Bacteroidaceae","Eubacteriaceae","Lachnospiraceae","Porphyromonadaceae"),
  GroupA=c(37.7397,31.34317,222.08827,5.08956,3.7393),
  GroupB=c(113.2191,94.02951,66.26481,15.26868,11.2179),
  GroupC=c(123.2191,94.02951,46.26481,35.26868,1.2179),
  GroupD=c(37.7397,31.34317,222.08827,5.08956,3.7393)
)

df.long <- df %>% gather(group, abundance, -Phylum)

## 组间连线数据：
## 假设第一列是Phylum
link_dat <- df %>% 
  arrange(by=desc(Phylum)) %>% 
  mutate_if(is.numeric, cumsum) 
bar.width <- 0.7
link_dat <- link_dat[, c(1,2,rep(3:(ncol(link_dat)-1),each=2), ncol(link_dat))]
link_dat <- data.frame(y=t(matrix(t(link_dat[,-1]), nrow=2)))
link_dat$x.1 <- 1:(ncol(df)-2)+bar.width/2
link_dat$x.2 <- 1:(ncol(df)-2)+(1-bar.width/2)

p <- ggplot(df.long, aes(x=group, y=abundance, fill=Phylum)) + 
     geom_bar(stat = "identity", width=bar.width, col='black') + 
     geom_segment(data=link_dat, aes(x=x.1, xend=x.2, y=y.1, yend=y.2), inherit.aes = F)

pdf(file="linkedSegments.pdf")
print(p)
dev.off()
