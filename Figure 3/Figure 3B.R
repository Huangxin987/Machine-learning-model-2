library(ggplot2)
#Figure3B_bacteria
data<-read.csv("ACE.CSV",row.names=1)
data$Day<-as.numeric(as.character(data$Day))
ggplot(data, aes(x = Day, y = Shannon, color = Bodysite, group = Bodysite)) +
  geom_smooth(size = 1) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_color_manual(values = c('#E64A34', '#4DBBD5', '#2C69B0', '#FFCA55', '#FF9FF3')) +
  theme_bw() +
  ylab('ACE index') +
  xlab("Postmortem interval (day)") +
  scale_x_continuous(breaks = c(2, 3, 5, 7, 10)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
