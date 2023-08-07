#TMEJ_predictor
library(ggplot2)
library(ggrepel)
library(ggiraph)

data<-read.csv("TMEJ_predictor_R26A.csv", header = T)
ldel<-c(data$ldel)
rdel<-c(data$rdel)
MH1<-c(data$inv.TINS.location)
MH1size<-c(data$seed.MH.size)
MH2size<-c(data$MH.size)
resolve_MH_seq<-c(data$resolve.MH.seq)
flank<-c(data$template.flank)
seed_seq<-c(data$seed.MH.seq)
seed_location<-c(data$seed.MH.location)
hairpin<-c(data$hairpin.distance)
TINS_seq<-c(data$TINS.seq)
ID<-c(data$TINS.class)
family<-c(data$junction.family)
x.index<-c(data$x.position)
y.index<-c(data$y.position)

TMEJ<-ggplot() +
  aes(ldel,rdel,color=ID,data_id=seed_seq,text=paste
      ("seed MH:",seed_seq, "\nseed pos:",seed_location,"\nhairpin:",hairpin,
        "\ninsertion:",TINS_seq,"\nresolve MH:",resolve_MH_seq,"\nflank:",flank)) +
  annotate("text",x=median(x.index),y=max(y.index)+3,label="Flank.Priming MH.Location") +
  annotate("text",x=20,y=35,label="TMEJ_predictor_R26A") +
  geom_point_interactive(data=data,x=ldel,y=rdel,size=MH2size) +
  geom_path_interactive(group=family,aes(color=family),size=0.1) +
  geom_label_interactive(label = family,aes(color=family),x=x.index,y=y.index) +
  scale_x_continuous(limit = c(0,40),breaks = c(0,10,20,30,40), minor_breaks = c(5,15,25,35)) +
  scale_y_continuous(limit = c(0,40),breaks = c(0,10,20,30,40), minor_breaks = c(5,15,25,35)) +
  scale_color_manual(name="Product Class",breaks=ID,values=c('MH.del'='#0000FF','inv'='#FF00FF','ss'='#FF9100'))

girafe(ggobj = TMEJ,
       width_svg = 8, height_svg = 6, #sizes the output plot
       options = list(
         opts_tooltip(
           opacity = 0.8, 
           css = "background-color:#4c6061; color:white; padding:10px; border-radius:5px;"),
         opts_hover_inv(css = "opacity:0.02;"),
         opts_hover(css = "stroke-width: 2; opacity: 1;")))