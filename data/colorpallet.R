
setwd("~/Documents/forecast_vassarlabs/data")
# #brk=c(seq(0,40,1),seq(45,100,5),seq(110,250,10))/3.6
# brk=c(seq(0,30,1),seq(32,60,2),seq(63,100,3))/3.6
# # brk=c(seq(0,40,1),seq(,40,1),seq(0,40,1))
# # col=colorRampPalette(c("#FFFFFF", "#4795B1","#F5E75B","#D82928","#921519","darkorchid","beige"))(length(brk)-1) 
# #col=colorRampPalette(c("cornflowerblue","#4795B1", "seagreen","#F5E75B","#D82928","#921519","darkorchid","beige"))(length(brk)-1) 
# 
# #windy
# col=colorRampPalette(c("#6271B7","#6271B7","#6271B7","#6271B7","#3D6EA2","#4A93A9","#4A9193","#4D8D7C","#4CA34C","#67A336","#A18640","#A16D5C","#8C3F5C","#964B90", "#5F649F", "#5B87A0","#5B87A0"))(length(brk)-1) 
# 
# #col=colorRampPalette(c("#4a47e8","#2F2CE0","#007BFF","#1DBEC6","#068556","#039603","#53D304","#E0C900","#C4032D","#dd8444","#795548"))(length(brk)-1) 
# 
# #col=colorRampPalette(c("#5E5777","#5D6095","#4184BB","#665E7F","#59C049","#90CB4D","#CDD83C","#DAC03A","#D99D44","#D7794D","#C94A72","#A43C5C","#903250","#631D1B"))(length(brk)-1) 
# 
# 
# plot(brk*3.6,brk/brk,col=col,type = "p",lwd = 10)
# #"limegreen",
# #bakup col=colorRampPalette(c("#FFFFFF", "#4795B1","lightgreen","#F5E75B","#D82928","#921519","red4","darkorchid","beige"))(length(brk)-1) 
# 
# write.csv(col,file = "col.csv")
# write.csv(brk,file = "brk.csv")
  

#RAIN
#brk=c(0,10,20,35,55,70,100)/3.6 #windy
brk=c(0,15,30,45,60,75,100)/3.6 #windy

brk_hr=c(brk[1])
N=20
for (i in 1:(length(brk)-1)) {
  print(i)
  brk_hr=c(brk_hr,seq(brk[i],brk[(i+1)],(brk[i+1]-brk[i])/N)[-1])
}

#windy
#col=colorRampPalette(c("#6271B7","#3D6EA2","#4A93A9","#4A9193","#4D8D7C","#4CA34C","#67A336","#A18640","#A16D5C","#8C3F5C","#964B90", "#5F649F", "#5B87A0"))(length(brk_hr)-1) 

col=colorRampPalette(c("purple4","#3D6EA2","#4A93A9","#4A9193","#4D8D7C","#4CA34C","#67A336","#A18640","#A16D5C","#8C3F5C","#964B90", "#5F649F", "#5B87A0"))(length(brk_hr)-1) 

plot(brk_hr[-1],brk_hr[-1]/brk_hr[-1],col=col,type = "p",lwd = 10,ylim=c(0.99,1.1),xlab = "Wind speed(m/s)")

write.csv(col,file = "col.csv")
write.csv(brk_hr,file = "brk.csv")

#temperature
brk=c(-20,-10,0,10,20,30,40,50) #windy

brk_hr=c(brk[1])
N=20
for (i in 1:(length(brk)-1)) {
  print(i)
  brk_hr=c(brk_hr,seq(brk[i],brk[(i+1)],(brk[i+1]-brk[i])/N)[-1])
}


col=colorRampPalette(c(rgb(149, 137, 211,maxColorValue = 256), 
rgb(150, 209, 216,maxColorValue = 256) ,
rgb(129, 204, 197,maxColorValue = 256) ,
rgb(103, 180, 186,maxColorValue = 256) ,
rgb(95, 143, 197,maxColorValue = 256) ,
rgb(80, 140, 62,maxColorValue = 256),
rgb(121, 146, 28,maxColorValue = 256), 
rgb(171, 161, 14,maxColorValue = 256) ,
rgb(223, 177, 6,maxColorValue = 256) ,
rgb(243, 150, 6,maxColorValue = 256) ,
rgb(236, 95, 21,maxColorValue = 256) ,
rgb(190, 65, 18,maxColorValue = 256) ,
rgb(138, 43, 10,maxColorValue = 256) ,
rgb(255, 51, 249,maxColorValue = 256)))(length(brk_hr)-1) 

#rgb(138, 43, 10,maxColorValue = 256) last backup

plot(brk_hr[-1],brk_hr[-1]/brk_hr[-1],col=col,type = "p",lwd = 10,ylim=c(0.99,1.1),xlab = ("Temperature (deg C)"))

write.csv(col,file = "col_tem.csv")
write.csv(brk_hr,file = "brk_tem.csv")

#RH
brk=c(0,30,50,80,90,100) #windy

brk_hr=c(brk[1])
N=500
for (i in 1:(length(brk)-1)) {
  print(i)
  brk_hr=c(brk_hr,seq(brk[i],brk[(i+1)],(brk[i+1]-brk[i])/N)[-1])
}



col=colorRampPalette(c(rgb(174, 110, 56,maxColorValue = 256),
rgb(174, 146, 56,maxColorValue = 256), 
rgb(105, 174, 56,maxColorValue = 256), 
rgb(54, 175, 148,maxColorValue = 256), 
rgb(56, 158, 174,maxColorValue = 256), 
rgb(56, 142, 174,maxColorValue = 256), 
rgb(56, 133, 174,maxColorValue = 256), 
rgb(56, 111, 166,maxColorValue = 256), 
rgb(56, 71, 116,maxColorValue = 256)))(length(brk_hr)-1) 

plot(brk_hr[-1],brk_hr[-1]/brk_hr[-1],col=col,type = "p",lwd = 10,ylim=c(0.99,1.1),xlab = ("RH (%)"))

write.csv(col,file = "col_rh.csv")
write.csv(brk_hr,file = "brk_rh.csv")


