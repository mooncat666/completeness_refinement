require(ggplot2)
require(ggalluvial)
require(scales)
require(munsell)
require(dplyr)

df <- read.csv('results.tsv',sep='\t')

# arrow plot

arrow_plot <- ggplot(df)+
  geom_segment(aes(x=Completeness,y=Contamination,xend=NewCompl,yend=NewCont), arrow=arrow(length = unit(0.1,"inches")),col='gray85')+
  geom_point(aes(x=Completeness,y=Contamination,shape=init_GUNC),col=mnsl('2.5P 4/12'))+
  geom_point(aes(x=NewCompl,y=NewCont,shape=new_GUNC),col=mnsl('5YR 7/12'))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base=10),breaks = c(0,1,10,50,100),labels = c(0,1,10,50,100))+
  coord_flip()+
  theme_minimal()+
  scale_shape_manual(name="Pass GUNC", labels=c('False','True'),values=c('False'='circle','True'='triangle'))  
ggsave("arrow_plot.pdf", plot=arrow_plot)


# arrow plot 3: add 1.1 gunc refined step to a plot

df3 <- read.csv('results3.tsv',sep='\t')

arrow_plot3 <- ggplot(df3)+
  geom_segment(aes(x=Completeness,y=Contamination,xend=GuncRefCompl,yend=GuncRefCont), arrow=arrow(length = unit(0.07,"inches")),col='gray85',size=0.4)+
  geom_segment(aes(x=GuncRefCompl,y=GuncRefCont,xend=NewCompl,yend=NewCont), arrow=arrow(length = unit(0.07,"inches")),col='gray85',size=0.4)+
  geom_point(aes(x=Completeness,y=Contamination,shape=init_GUNC),col=mnsl('2.5P 4/12'),size=1)+
  geom_point(aes(x=GuncRefCompl,y=GuncRefCont,shape=guncref_GUNC),col=mnsl('2.5P 7/12'),size=1)+
  geom_point(aes(x=NewCompl,y=NewCont,shape=new_GUNC),col=mnsl('5YR 7/12'),size=1)+
  scale_y_continuous(trans=scales::pseudo_log_trans(base=10),breaks = c(0,1,10,50,100),labels = c(0,1,10,50,100))+
  coord_flip()+
  theme_minimal()+
  scale_shape_manual(name="Pass GUNC", labels=c('False','True'),values=c('False'='circle','True'='triangle'))
ggsave("arrow_plot_3.pdf", plot=arrow_plot3)


# box plots

n <- length(df$Completeness)
box_df <- data.frame(Completeness=c(df$Completeness,df$NewCompl), Contamination=c(df$Contamination, df$NewCont),State=c(rep('Initial', n),rep('Refined', n)))

completeness_box_plot <- ggplot(box_df,aes(x=State, y=Completeness))+
  geom_boxplot()+
  theme_minimal()+
  xlab('')
contamination_box_plot <- ggplot(box_df,aes(x=State, y=Contamination))+
  geom_boxplot()+
  theme_minimal()
ggsave("completeness_box_plot.pdf", plot=completeness_box_plot)
ggsave("contamination_box_plot.pdf", plot=contamination_box_plot)
 
# alluvial plot

df <- mutate(df, init_Quality = case_when(Completeness < 50 | Contamination >= 10 ~ "low", Completeness>90 & Contamination < 5 ~ "high", Completeness >= 50 & Contamination < 10 ~ "medium"), 
            ref_Quality = case_when(NewCompl < 50 | NewCont >= 10 ~ "low", NewCompl>90 & NewCont < 5 ~ "high", NewCompl >= 50 & NewCont < 10 ~ "medium"))
df$init_Quality <- factor(df$init_Quality, levels=c('high','medium','low'))

high=mnsl('2.5G 9/8')
medium=mnsl('10B 7/8')
low=mnsl('2.5YR 7/8')

# funny colours
#high='#A8DCD9'
#medium='#FFD670'
#low='#E8C7DE'

alluvial_plot <- ggplot(df, aes(axis1=init_Quality,axis2=ref_Quality))+
  geom_alluvium(aes(fill=init_Quality))+
  geom_stratum(width=0.12, col='white', fill='grey55')+
  scale_x_discrete(limits = c("Initial", "Refined"), expand = c(0, 0))+
  theme_minimal()+
  scale_fill_manual(name='Quality', values=c('high'=high, 'medium'=medium, 'low'=low))+
  scale_color_manual(values=c('high'=high, 'medium'=medium, 'low'=low),guide='none')+
  geom_text(stat='stratum',aes(label = after_stat(stratum),col=after_stat(stratum)), size=3, angle=90)+
  theme(axis.text.y=element_blank(),panel.grid=element_blank())
ggsave("alluvial_plot.pdf", plot=alluvial_plot)
