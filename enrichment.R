data1 <- data.frame(read_excel("miR-target.xlsx",sheet=1))
data2 <- data.frame(read_excel("miR-target.xlsx",sheet=2))
data3 <- data.frame(read_excel("miR-target.xlsx",sheet=3))
data4 <- data.frame(read_excel("miR-target.xlsx",sheet=4))
data5 <- data.frame(read_excel("miR-target.xlsx",sheet=5))
data6 <- data.frame(read_excel("miR-target.xlsx",sheet=6))

temp2 <- data.frame(read_excel("Table-S5.xlsx",sheet=1,col_names=F))
dim(temp2)
temp2<- read.table("InnateDB_genes_list.txt",head=T)
dim(temp2)
temp2<- read.table("FDA_approved_drugtarget_list.txt",head=T)
dim(temp2)

x1<-data1[which(data1[,5] %in% temp2[,1]),5]
x2<-data2[which(data2[,5] %in% temp2[,1]),5]
x3<-data3[which(data3[,5] %in% temp2[,1]),5]
x4<-data4[which(data4[,5] %in% temp2[,1]),5]
x5<-data5[which(data5[,5] %in% temp2[,1]),5]
x6<-data6[which(data6[,5] %in% temp2[,1]),5]

x<-unique(c(x1,x2,x3,x4,x5,x6))
length(x)
nrow(temp2)

q<-1-phyper(792,4723,37875,4427, lower.tail=TRUE)
q<-1-phyper(23,123,37875,4427, lower.tail=TRUE)
q<-1-phyper(12,123,37875,2001, lower.tail=TRUE)
q<-1-phyper(143,672,37875,4427, lower.tail=TRUE)

(792/4723)/(4427/37875)
(143/672)/(4427/37875)
(23/123)/(4427/37875)

input<-unique(as.character(read.table("ENSG.ENST.ENSP.Symbol.hg19.bed.txt",head=F)[,4]))
P<-c()
for(j in 1:50000){
x<-input[sample(1:length(input),4427)]
p<-sum(x %in% temp2[,1])
P<-c(P,p)
}
sum(P>792)/length(P)
hist(P,breaks = 200,col="blue",xlim=c(min(P),800),main="IRG: immnue related genes")
abline(v=782,col="red",lty=3,lwd=2)

P<-c()
for(j in 1:50000){
  x<-input[sample(1:length(input),4427)]
  p<-sum(x %in% temp2[,1])
  P<-c(P,p)
}
sum(P>23)/length(P)
hist(P,breaks = 200,col="blue",xlim=c(min(P),max(P)),main="GWAS")
abline(v=23,col="red",lty=3,lwd=2)

P<-c()
for(j in 1:50000){
  x<-input[sample(1:length(input),4427)]
  p<-sum(x %in% temp2[,1])
  P<-c(P,p)
}
sum(P>143)/length(P)
hist(P,breaks = 200,col="blue",xlim=c(min(P),max(P)),main="FDA drug targets")
abline(v=143,col="red",lty=3,lwd=2)
