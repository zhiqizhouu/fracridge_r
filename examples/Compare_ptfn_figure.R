library(reshape2)
library(tidyverse)

# create legend
data_p = read.delim("p_param.txt",sep = " ")
x = c(10,20,50,100,200,300,500,750,1000,1500,2000,3000)
name = c("Standard Ridge Regression", "SVD Ridge Regression", "Fractional Ridge Regression", "Glmnet Ridge Regression", "x")
data_p = t(data_p) %>% cbind(x)
colnames(data_p) = name
rownames(data_p) = NULL
dat1 = as.data.frame(data_p) 
dat1 = melt(dat1,id.vars = c("x"))
dat1

ggplot(dat1)+
  geom_line(aes(x = dat1[,1],y = dat1[,3],color = dat1[,2]))+
  geom_point(aes(x = dat1[,1],y = dat1[,3],color = dat1[,2],shape = dat1[,2],size = dat1[,2]))+
  scale_size_manual("",values=c(2.5,4,2.5,3.5))+
  scale_color_manual("", values = c("red","orange","green","blue"))+
  labs(x = "p", y = "Execution time (s)", title = "")+
  scale_shape_manual("", values = c(15, 16, 17, 18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_y_continuous(breaks = seq(0, 600, len = 7))+
  geom_vline(xintercept=100,color = "gray") +
  annotate("text", x=500, y=500, label="p = 100", angle=0, size=4.5, color="black")

##### p 
data_p = read.delim("p_param.txt",sep = " ")
a = t(data_p) %>% cbind(x) 
colnames(a) = name
rownames(a) = NULL
a = as.data.frame(a)
p1 = ggplot(a)+
  geom_line(aes(x = a[,5],y = a[,1]),color = "red")+
  geom_line(aes(x = a[,5],y = a[,2]),color = "orange")+
  geom_line(aes(x = a[,5],y = a[,3]),color = "green")+
  geom_line(aes(x = a[,5],y = a[,4]),color = "blue")+
  geom_point(aes(x = a[,5],y = a[,1]),color = "red", shape = 15, size = 2.5)+
  geom_point(aes(x = a[,5],y = a[,2]),color = "orange", shape = 16, size = 4.)+
  geom_point(aes(x = a[,5],y = a[,3]),color = "green", shape = 17, size = 2.5)+
  geom_point(aes(x = a[,5],y = a[,4]),color = "blue", shape = 18, size = 3.5)+
  labs(x = "p", y = "Execution time (s)", title = "")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_y_continuous(breaks = seq(0, 600, len = 7))+
  geom_vline(xintercept=100,color = "gray") +
  annotate("text", x=300, y=500, label="p = 100", angle=0, size=4.5, color="black")
p1






######## f
data_f = read.delim("num_params.txt",sep = " ")
x = c(1,5,10,20,50,75,100)
name = c("Standard Ridge Regression", "SVD Ridge Regression", "Fractional Ridge Regression", "Glmnet Ridge Regression","x")
data_f = t(data_f) %>% cbind(x) 
data_f = as.data.frame(data_f)
colnames(data_f) = name
rownames(data_f) = NULL
data_f


p2 = ggplot(data_f)+
  geom_line(aes(x = data_f[,5],y = data_f[,1]),color = "red")+
  geom_line(aes(x = data_f[,5],y = data_f[,2]),color = "orange")+
  geom_line(aes(x = data_f[,5],y = data_f[,3]),color = "green")+
  geom_line(aes(x = data_f[,5],y = data_f[,4]),color = "blue")+
  geom_point(aes(x = data_f[,5],y = data_f[,1]),color = "red", shape = 15, size = 2.5)+
  geom_point(aes(x = data_f[,5],y = data_f[,2]),color = "orange", shape = 16, size = 4.)+
  geom_point(aes(x = data_f[,5],y = data_f[,3]),color = "green", shape = 17, size = 2.5)+
  geom_point(aes(x = data_f[,5],y = data_f[,4]),color = "blue", shape = 18, size = 3.5)+
  labs(x = "f", y = "", title = "")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #scale_y_continuous(breaks = seq(0, 8, len = 20))+
  geom_vline(xintercept=10,color = "gray") +
  annotate("text", x=15, y=7.37, label="f = 10", angle=0, size=4.5, color="black")
p2







#####  n
data_n = read.delim("num_size_n.txt",sep = " ")
x = c(25,50,100,500,1000,5000,10000,15000,20000)
name = c("Standard Ridge Regression", "SVD Ridge Regression", "Fractional Ridge Regression", "Glmnet Ridge Regression","x")
c = t(data_n) %>% cbind(x) 
c = as.data.frame(c)
colnames(c) = name
rownames(c) = NULL
p3 = ggplot(c)+
  geom_line(aes(x = c[,5],y = c[,1]),color = "red")+
  geom_line(aes(x = c[,5],y = c[,2]),color = "orange")+
  geom_line(aes(x = c[,5],y = c[,3]),color = "green")+
  geom_line(aes(x = c[,5],y = c[,4]),color = "blue")+
  geom_point(aes(x = c[,5],y = c[,1]),color = "red", shape = 15, size = 2.5)+
  geom_point(aes(x = c[,5],y = c[,2]),color = "orange", shape = 16, size = 4.)+
  geom_point(aes(x = c[,5],y = c[,3]),color = "green", shape = 17, size = 2.5)+
  geom_point(aes(x = c[,5],y = c[,4]),color = "blue", shape = 18, size = 3.5)+
  labs(x = "n", y = "", title = "")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #scale_y_continuous(breaks = seq(0, 8, len = 20))+
  geom_vline(xintercept=5000,color = "gray") +
  annotate("text", x=6500, y=2.958, label="n = 5,000", angle=0, size=4.5, color="black")
p3



#######   t
data_t = read.delim("num_targets.txt",sep = " ")
x = c(5,25,50,100,500,1000,2500,5000)
name = c("Standard Ridge Regression", "SVD Ridge Regression", "Fractional Ridge Regression", "Glmnet Ridge Regression","x")
d = t(data_t) %>% cbind(x) 
d = as.data.frame(d)
colnames(d) = name
rownames(d) = NULL
p4 = ggplot(d)+
  geom_line(aes(x = d[,5],y = d[,1]),color = "red")+
  geom_line(aes(x = d[,5],y = d[,2]),color = "orange")+
  geom_line(aes(x = d[,5],y = d[,3]),color = "green")+
  geom_line(aes(x = d[,5],y = d[,4]),color = "blue")+
  geom_point(aes(x = d[,5],y = d[,1]),color = "red", shape = 15, size = 2.5)+
  geom_point(aes(x = d[,5],y = d[,2]),color = "orange", shape = 16, size = 4.)+
  geom_point(aes(x = d[,5],y = d[,3]),color = "green", shape = 17, size = 2.5)+
  geom_point(aes(x = d[,5],y = d[,4]),color = "blue", shape = 18, size = 3.5)+
  labs(x = "t", y = "", title = "")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #scale_y_continuous(breaks = seq(0, 8, len = 20))+
  geom_vline(xintercept=1,color = "gray") +
  annotate("text", x=250, y=200, label="t = 1", angle=0, size=4.5, color="black")
p4


p1
p2
p3
p4
