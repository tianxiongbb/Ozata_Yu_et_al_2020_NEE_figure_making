setwd("/Users/bigbear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
# annotations----
human_pi_list=read.table("/Users/bigbear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/piRNA_cluster_annotation/piRNA.list.txt",header=F,row.names=1)
human_gene_list=read.table("/Users/bigbear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/hg19.gene.category.bed",header=F,row.names=NULL)
pd=read.table("/Users/bigbear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/human.amybPeakMerged.distance.txt",header=F,row.names=1)
human_pi=row.names(human_pi_list)
human_pachy=human_pi[which(human_pi_list[,1]=="Pachytene")]
human_pachy=setdiff(human_pachy,"17-q25-201")
human_pachy_amyb=human_pachy[pd[human_pachy,1]<=500]
human_pachy_noamyb=human_pachy[pd[human_pachy,1]>500]
human_prepachy=human_pi[which(human_pi_list[,1]=="Prepachytene")]
human_prepachy=setdiff(human_prepachy,"pi-ACADM")
human_prepachy_genic1=as.vector(read.table("/Users/bigbear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/genic_prepachy.list",header=F,row.names=NULL)[,1])
human_prepachy_genic=as.vector(read.table("/Users/bigbear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/genic_prepachy.list1",header=F,row.names=NULL)[,1])
human_prepachy_genic=intersect(human_prepachy_genic,human_prepachy)
human_prepachy_inter=setdiff(human_prepachy,human_prepachy_genic)
human_pi=setdiff(human_pi,c("pi-ACADM","17-q25-201"))
human_pi=setdiff(human_pi,c("pi-ACADM"))
human_hybrid=human_pi[which(human_pi_list[,1]=="Hybrid")]
human_bid=human_pi[which(human_pi_list[,2]=="bidirectional")]
human_pc=as.vector(human_gene_list[which(human_gene_list[,2]=="protein_coding"),1])
human_linc=as.vector(human_gene_list[which(human_gene_list[,2]=="lincRNA"),1])
pl=c("DDX4","DDX39A","FKBP6","GPAT2","HENMT1","MAEL","MOV10L1","MYBL1","NFYA",
     "PIWIL1","PIWIL2","PIWIL4","PLD6","TDRD1","TDRD3","TDRD5","TDRD6","TDRD9",
     "TDRD12","TDRKH")
srna_rn=c("N1643-5YO-Rep1-Ox","N1643-5YO-Rep2-Ox","TTC-Ox","LCH-008-Ox-ID5",
         "N7377-Ox","S8+","s904-Ox-ID3","Ambion-Ox",
         "S6+","P75-Ox-ID1","8289NT_Ox","6615NT_Ox",
         "N7358-Ox","S7+","N8150-Ox","9268NT_Ox","8357NT_Ox","8377NT_Ox")
rna_rn=c("D1643_5YO_Rep1","D1643_5YO_Rep2","TTC-Juv-ID4_S4","LCH-008-Juv-ID3_S3",
         "D7377NT","S8_NT","s904-Ad-ID2_S1","AmbionNT",
         "S6_NT","P75-Ad-ID1_S2","8289NT-ID4_S6","6615NT-ID3_S4",
         "D7358NT","S7_NT","D8150NT","9268NT-ID8_S5","8357NT-ID7_S1","8377NT-ID5_S3")

# colors
csl_aj=c("#c51b7d","#4d9221")
tcs=c(rep("#8dd3c7",4),rep("#fccde5",3),rep("#b3de69",5),rep("#fdb462",6))

# data reading----
df_srna=read.table("merged.rpm",header=T,row.names=1)
df_rna=read.table("merged.rnaseq.gene+pi.rpkm",header=T,row.names=1, check.names=F)
df_rna_rebuttal=read.table("merged.rnaseq.gene+pi.rpkm",header=T,row.names=1, check.names=F)

# Jan.2018----
setwd("/Users/BigBear/My_Publish/Human_piRNA/CUT_RUN/Jan.2018/")

mat_amyb_4N=read.table("4N-AMYB.mat",header=F,row.names=1)
mat_amyb_Prepachy=read.table("PrePACH-AMYB.mat",header=F,row.names=1)
mat_amyb_rST=read.table("Rs-AMYB.mat",header=F,row.names=1)
mat_pol2_4N=read.table("4N-RNAP2.mat",header=F,row.names=1)
mat_pol2_Prepachy=read.table("PrePACH-RNAP2.mat",header=F,row.names=1)
mat_pol2_rST=read.table("Rs-RNAP2.mat",header=F,row.names=1)
mat_nfya_4N=read.table("4N-NFYA.mat",header=F,row.names=1)
mat_nfya_Prepachy=read.table("PrePACH-NFYA.mat",header=F,row.names=1)
mat_nfya_rST=read.table("Rs-NFYA.mat",header=F,row.names=1)
mat_tcfl5_4N=read.table("4N-TCFL5.mat",header=F,row.names=1)
mat_tcfl5_Prepachy=read.table("PrePACH-TCFL5.mat",header=F,row.names=1)
mat_tcfl5_rST=read.table("Rs-TCFL5.mat",header=F,row.names=1)
mat_neg_4N=read.table("4N-Neg-CNT.mat",header=F,row.names=1)
mat_neg_Prepachy=read.table("PrePACH-Neg-CNT.mat",header=F,row.names=1)
mat_neg_rST=read.table("Rs-Neg-CNT.mat",header=F,row.names=1)

al=row.names(mat_amyb_4N)
al=setdiff(al,pi)
list1=al[which(Exp_a[al,1]<0.1 & Exp_a[al,2]>10)]
list2=al[which(Exp_a[al,2]<0.1 & Exp_a[al,1]>10)]

fun_plot=function(data_in,fn,m){
  annoR=data.frame(Location=factor(rep(c("Intergenic","Genic"),c(114,100))),Type=factor(rep(c("Pachytene","Prepachytene","Hybrid"),c(100,84,30))))
  row.names(annoR)=c(Pachy,Prepachy,Hybrid)
  annoR[grep("pi-",c(Pachy,Prepachy,Hybrid)),1]="Genic"
  annoR[grep("pi-",c(Pachy,Prepachy,Hybrid),invert=T),1]="Intergenic"
  ann_colors=list(Type=c(Pachytene="#E41A1C",Hybrid="#984ea3",Prepachytene="#4daf4a"),
                  position=c(beforeTSS="grey",afterTSS="navy"),Location=c(Intergenic="black",Genic="grey"))
  x=as.matrix(data_in[c(Pachy,Prepachy,Hybrid),])
  colnames(x)=c("-5kb",rep("",99),"TSS",rep("",98),"5kb")
  pheatmap(x,
           col=colorRampPalette(c(brewer.pal(9,"Blues")[9:1],"white"))(100),
           breaks=0:100/100*m,
           show_rownames=FALSE,
           cluster_col=FALSE,cluster_row=FALSE,cellwidth=1,cellheight=3,
           annotation_row=annoR,annotation_colors=ann_colors,cex=1,
           filename=paste("../../figures/Jan23.2018.CUTRUN/",fn,sep=""),
           legend_labels=c("low","high"),legend_breaks=c(0,m))
}
fun_plot1=function(data_in,fn,m){
  annoR=data.frame(Location=factor(rep(c("Intergenic","Genic"),c(114,100))),Type=factor(rep(c("Pachytene","Prepachytene","Hybrid"),c(100,84,30))))
  row.names(annoR)=c(Pachy,Prepachy,Hybrid)
  annoR[grep("pi-",c(Pachy,Prepachy,Hybrid)),1]="Genic"
  annoR[grep("pi-",c(Pachy,Prepachy,Hybrid),invert=T),1]="Intergenic"
  ann_colors=list(Type=c(Pachytene="#E41A1C",Hybrid="#984ea3",Prepachytene="#4daf4a"),
                  position=c(beforeTSS="grey",afterTSS="navy"),Location=c(Intergenic="black",Genic="grey"))
  x=as.matrix(data_in[c(Pachy,Prepachy,Hybrid),])
  for(i in 1:214){
    rn=max(x[i,])-min(x[i,])+1e-6
    x[i,]=(x[i,]-min(x[i,]))/rn
  }
  colnames(x)=c("-5kb",rep("",99),"TSS",rep("",98),"5kb")
  pheatmap(x,
           col=colorRampPalette(c(brewer.pal(9,"Blues")[9:1],"white"))(100),
           breaks=0:100/100,
           show_rownames=FALSE,
           cluster_col=FALSE,cluster_row=FALSE,cellwidth=1,cellheight=3,
           annotation_row=annoR,annotation_colors=ann_colors,cex=1,
           filename=paste("../../figures/Jan23.2018.CUTRUN/Rowscaled_",fn,sep=""),
           legend_labels=c("low","high"))
}
fun_plot2=function(data_in,fn,m){
  annoR=data.frame(Stage=factor(rep(c("Prepachy_Gene","Pachy_Gene"),c(159,703))))
  row.names(annoR)=c(list2,list1)
  ann_colors=list(Stage=c(Prepachy_Gene="black",Pachy_Gene="grey"))
  x=as.matrix(data_in[c(list2,list1),])
  colnames(x)=c("-5kb",rep("",99),"TSS",rep("",98),"5kb")
  pheatmap(x,
           col=c(brewer.pal(9,"Blues")[9:1],"white"),
           breaks=0:10/10*m,
           show_rownames=FALSE,
           cluster_col=FALSE,cluster_row=FALSE,cellwidth=1,cellheight=1,
           annotation_row=annoR,annotation_colors=ann_colors,cex=1,
           filename=paste("../../figures/genes_",fn,sep=""),
           legend_labels=c("low","high"),legend_breaks=c(0,m))
}

fun_plot1(mat_amyb_4N,"amyb_4N_allreads.pdf",5)
fun_plot1(mat_amyb_Prepachy,"amyb_prePachy_allreads.pdf",5)
fun_plot1(mat_amyb_rST,"amyb_rST_allreads.pdf",5)
fun_plot1(mat_pol2_4N,"pol2_4N_allreads.pdf",5)
fun_plot1(mat_pol2_Prepachy,"pol2_prePachy_allreads.pdf",5)
fun_plot1(mat_pol2_rST,"pol2_rST_allreads.pdf",5)
fun_plot1(mat_nfya_4N,"nfya_4N_allreads.pdf",5)
fun_plot1(mat_nfya_Prepachy,"nfya_prePachy_allreads.pdf",5)
fun_plot1(mat_nfya_rST,"nfya_rST_allreads.pdf",5)
fun_plot1(mat_tcfl5_4N,"tcfl5_4N_allreads.pdf",5)
fun_plot1(mat_tcfl5_Prepachy,"tcfl5_prePachy_allreads.pdf",5)
fun_plot1(mat_tcfl5_rST,"tcfl5_rST_allreads.pdf",5)
fun_plot1(mat_neg_4N,"neg_4N_allreads.pdf",5)
fun_plot1(mat_neg_Prepachy,"neg_prePachy_allreads.pdf",5)
fun_plot1(mat_neg_rST,"neg_rST_allreads.pdf",5)

fun_plot(mat_amyb_4N,"amyb_4N_allreads.pdf",5)
fun_plot(mat_amyb_Prepachy,"amyb_prePachy_allreads.pdf",5)
fun_plot(mat_amyb_rST,"amyb_rST_allreads.pdf",5)
fun_plot(mat_pol2_4N,"pol2_4N_allreads.pdf",5)
fun_plot(mat_pol2_Prepachy,"pol2_prePachy_allreads.pdf",5)
fun_plot(mat_pol2_rST,"pol2_rST_allreads.pdf",5)
fun_plot(mat_nfya_4N,"nfya_4N_allreads.pdf",5)
fun_plot(mat_nfya_Prepachy,"nfya_prePachy_allreads.pdf",5)
fun_plot(mat_nfya_rST,"nfya_rST_allreads.pdf",5)
fun_plot(mat_tcfl5_4N,"tcfl5_4N_allreads.pdf",5)
fun_plot(mat_tcfl5_Prepachy,"tcfl5_prePachy_allreads.pdf",5)
fun_plot(mat_tcfl5_rST,"tcfl5_rST_allreads.pdf",5)
fun_plot(mat_neg_4N,"neg_4N_allreads.pdf",5)
fun_plot(mat_neg_Prepachy,"neg_prePachy_allreads.pdf",5)
fun_plot(mat_neg_rST,"neg_rST_allreads.pdf",5)

# Oct.2017----
setwd("/Users/BigBear/My_Publish/Human_piRNA/CUT_RUN/Oct.2017/")

mat_amyb_pachy=read.table("Pach-Amyb.mat",header=F,row.names=1)
mat_nfya_pachy=read.table("Pach-NFYA.mat",header=F,row.names=1)
mat_tcfl5_pachy=read.table("Pach-TCFL5.mat",header=F,row.names=1)
mat_neg_pachy=read.table("Pach-Neg-CNT.mat",header=F,row.names=1)
mat_nfya_prepachy=read.table("PrePach-NFYA.mat",header=F,row.names=1)

al=row.names(mat_amyb_pachy)
al=setdiff(al,pi)
list1=al[which(Exp_a[al,1]<0.1 & Exp_a[al,2]>10)]
list2=al[which(Exp_a[al,2]<0.1 & Exp_a[al,1]>10)]

fun_plot=function(data_in,fn,m){
  annoR=data.frame(Location=factor(rep(c("Intergenic","Genic"),c(114,100))),Type=factor(rep(c("Pachytene","Prepachytene","Hybrid"),c(100,84,30))))
  row.names(annoR)=c(Pachy,Prepachy,Hybrid)
  annoR[grep("pi-",c(Pachy,Prepachy,Hybrid)),1]="Genic"
  annoR[grep("pi-",c(Pachy,Prepachy,Hybrid),invert=T),1]="Intergenic"
  ann_colors=list(Type=c(Pachytene="#E41A1C",Hybrid="#984ea3",Prepachytene="#4daf4a"),
                  position=c(beforeTSS="grey",afterTSS="navy"),Location=c(Intergenic="black",Genic="grey"))
  x=as.matrix(data_in[c(Pachy,Prepachy,Hybrid),])
  colnames(x)=c("-5kb",rep("",99),"TSS",rep("",98),"5kb")
  pheatmap(x,
           col=colorRampPalette(c(brewer.pal(9,"Blues")[9:1],"white"))(100),
           breaks=0:100/100*m,
           show_rownames=FALSE,
           cluster_col=FALSE,cluster_row=FALSE,cellwidth=1,cellheight=3,
           annotation_row=annoR,annotation_colors=ann_colors,cex=1,
           filename=paste("../../figures/Oct24.2017.cutrun2/",fn,sep=""),
           legend_labels=c("low","high"),legend_breaks=c(0,m))
}
fun_plot1=function(data_in,fn,m){
  annoR=data.frame(Location=factor(rep(c("Intergenic","Genic"),c(114,100))),Type=factor(rep(c("Pachytene","Prepachytene","Hybrid"),c(100,84,30))))
  row.names(annoR)=c(Pachy,Prepachy,Hybrid)
  annoR[grep("pi-",c(Pachy,Prepachy,Hybrid)),1]="Genic"
  annoR[grep("pi-",c(Pachy,Prepachy,Hybrid),invert=T),1]="Intergenic"
  ann_colors=list(Type=c(Pachytene="#E41A1C",Hybrid="#984ea3",Prepachytene="#4daf4a"),
                  position=c(beforeTSS="grey",afterTSS="navy"),Location=c(Intergenic="black",Genic="grey"))
  x=as.matrix(data_in[c(Pachy,Prepachy,Hybrid),])
  for(i in 1:214){
    rn=max(x[i,])-min(x[i,])+1e-6
    x[i,]=(x[i,]-min(x[i,]))/rn
  }
  colnames(x)=c("-5kb",rep("",99),"TSS",rep("",98),"5kb")
  pheatmap(x,
           col=colorRampPalette(c(brewer.pal(9,"Blues")[9:1],"white"))(100),
           breaks=0:100/100,
           show_rownames=FALSE,
           cluster_col=FALSE,cluster_row=FALSE,cellwidth=1,cellheight=3,
           annotation_row=annoR,annotation_colors=ann_colors,cex=1,
           filename=paste("../../figures/Oct24.2017.cutrun2/Rowscaled_",fn,sep=""),
           legend_labels=c("low","high"))
}
fun_plot2=function(data_in,fn,m){
  annoR=data.frame(Stage=factor(rep(c("Prepachy_Gene","Pachy_Gene"),c(159,703))))
  row.names(annoR)=c(list2,list1)
  ann_colors=list(Stage=c(Prepachy_Gene="black",Pachy_Gene="grey"))
  x=as.matrix(data_in[c(list2,list1),])
  colnames(x)=c("-5kb",rep("",99),"TSS",rep("",98),"5kb")
  pheatmap(x,
           col=c(brewer.pal(9,"Blues")[9:1],"white"),
           breaks=0:10/10*m,
           show_rownames=FALSE,
           cluster_col=FALSE,cluster_row=FALSE,cellwidth=1,cellheight=1,
           annotation_row=annoR,annotation_colors=ann_colors,cex=1,
           filename=paste("../../figures/genes_",fn,sep=""),
           legend_labels=c("low","high"),legend_breaks=c(0,m))
}
fun_plot_aggre=function(data_in,fn,m){
  plot(apply(data_in[Pachy,],2,mean,trim=0.05),type="l",lwd=2,
       xaxt="n",ylab="RPM",xlab="",main=m)
  axis(1,c(0,100,200),label=c("-5kb","TSS","5kb"))
}
fun_plot_aggre1=function(data_in,fn,m){
  plot(apply(data_in[Prepachy,],2,mean,trim=0.05),type="l",lwd=2,
       xaxt="n",ylab="RPM",xlab="",main=m)
  axis(1,c(0,100,200),label=c("-5kb","TSS","5kb"))
}
fun_plot_aggre2=function(data_in,fn,m){
  plot(apply(data_in[Hybrid,],2,mean,trim=0.05),type="l",lwd=2,
       xaxt="n",ylab="RPM",xlab="",main=m)
  axis(1,c(0,100,200),label=c("-5kb","TSS","5kb"))
}
fun_plot_aggre3=function(data_in,fn,m){
  plot(apply(data_in[list1,],2,mean,trim=0.05),type="l",lwd=2,
       xaxt="n",ylab="RPM",xlab="",main=m)
  axis(1,c(0,100,200),label=c("-5kb","TSS","5kb"))
}
fun_plot_aggre4=function(data_in,fn,m){
  plot(apply(data_in[list2,],2,mean,trim=0.05),type="l",lwd=2,
       xaxt="n",ylab="RPM",xlab="",main=m)
  axis(1,c(0,100,200),label=c("-5kb","TSS","5kb"))
}

pdf("../figures/Oct24.2017.cutrun/Pachytene_cluster.pdf",width=3,height=2.5)
par(mar=c(3,4,4,1),cex=0.6)
fun_plot_aggre(mat_amyb_pachy,"amyb_pachy_allreads.pdf","AMYB Pachytene")
fun_plot_aggre(mat_nfya_pachy,"nfya_pachy_allreads.pdf","NFYA Pachytene")
fun_plot_aggre(mat_tcfl5_pachy,"tcfl5_pachy_allreads.pdf","TCLF5 Pachytene")
fun_plot_aggre(mat_neg_pachy,"neg_pachy_allreads.pdf","Control Pachytene")
fun_plot_aggre(mat_nfya_prepachy,"nfya_prepachy_allreads.pdf","NFYA Prepachytene")
dev.off()
pdf("../figures/Oct24.2017.cutrun/Prepachytene_cluster.pdf",width=3,height=2.5)
par(mar=c(3,4,4,1),cex=0.6)
fun_plot_aggre1(mat_amyb_pachy,"amyb_pachy_allreads.pdf","AMYB Pachytene")
fun_plot_aggre1(mat_nfya_pachy,"nfya_pachy_allreads.pdf","NFYA Pachytene")
fun_plot_aggre1(mat_tcfl5_pachy,"tcfl5_pachy_allreads.pdf","TCLF5 Pachytene")
fun_plot_aggre1(mat_neg_pachy,"neg_pachy_allreads.pdf","Control Pachytene")
fun_plot_aggre1(mat_nfya_prepachy,"nfya_prepachy_allreads.pdf","NFYA Prepachytene")
dev.off()
pdf("../figures/Oct24.2017.cutrun/Hybrid_cluster.pdf",width=3,height=2.5)
par(mar=c(3,4,4,1),cex=0.6)
fun_plot_aggre2(mat_amyb_pachy,"amyb_pachy_allreads.pdf","AMYB Pachytene")
fun_plot_aggre2(mat_nfya_pachy,"nfya_pachy_allreads.pdf","NFYA Pachytene")
fun_plot_aggre2(mat_tcfl5_pachy,"tcfl5_pachy_allreads.pdf","TCLF5 Pachytene")
fun_plot_aggre2(mat_neg_pachy,"neg_pachy_allreads.pdf","Control Pachytene")
fun_plot_aggre2(mat_nfya_prepachy,"nfya_prepachy_allreads.pdf","NFYA Prepachytene")
dev.off()
pdf("../figures/Oct24.2017.cutrun/Pachy_genes.pdf",width=3,height=2.5)
par(mar=c(3,4,4,1),cex=0.6)
fun_plot_aggre3(mat_amyb_pachy,"amyb_pachy_allreads.pdf","AMYB Pachytene")
fun_plot_aggre3(mat_nfya_pachy,"nfya_pachy_allreads.pdf","NFYA Pachytene")
fun_plot_aggre3(mat_tcfl5_pachy,"tcfl5_pachy_allreads.pdf","TCLF5 Pachytene")
fun_plot_aggre3(mat_neg_pachy,"neg_pachy_allreads.pdf","Control Pachytene")
fun_plot_aggre3(mat_nfya_prepachy,"nfya_prepachy_allreads.pdf","NFYA Prepachytene")
dev.off()
pdf("../figures/Oct24.2017.cutrun/Prepachy_genes.pdf",width=3,height=2.5)
par(mar=c(3,4,4,1),cex=0.6)
fun_plot_aggre4(mat_amyb_pachy,"amyb_pachy_allreads.pdf","AMYB Pachytene")
fun_plot_aggre4(mat_nfya_pachy,"nfya_pachy_allreads.pdf","NFYA Pachytene")
fun_plot_aggre4(mat_tcfl5_pachy,"tcfl5_pachy_allreads.pdf","TCLF5 Pachytene")
fun_plot_aggre4(mat_neg_pachy,"neg_pachy_allreads.pdf","Control Pachytene")
fun_plot_aggre4(mat_nfya_prepachy,"nfya_prepachy_allreads.pdf","NFYA Prepachytene")
dev.off()

fun_plot(mat_amyb_pachy,"amyb_pachy_allreads.pdf",5)
fun_plot(mat_nfya_pachy,"nfya_pachy_allreads.pdf",5)
fun_plot(mat_nfya_pachy,"nfya_pachy_allreads1.pdf",5)
fun_plot(mat_tcfl5_pachy,"tcfl5_pachy_allreads.pdf",10)
fun_plot(mat_neg_pachy,"neg_pachy_allreads.pdf",5)
fun_plot(mat_nfya_prepachy,"nfya_prepachy_allreads.pdf",5)

fun_plot1(mat_amyb_pachy,"amyb_pachy_allreads.pdf",5)
fun_plot1(mat_nfya_pachy,"nfya_pachy_allreads.pdf",5)
fun_plot1(mat_tcfl5_pachy,"tcfl5_pachy_allreads.pdf",10)
fun_plot1(mat_neg_pachy,"neg_pachy_allreads.pdf",5)
fun_plot1(mat_nfya_prepachy,"nfya_prepachy_allreads.pdf",5)

fun_plot2(mat_amyb_pachy,"amyb_pachy_allreads.pdf",5)
fun_plot2(mat_nfya_pachy,"nfya_pachy_allreads.pdf",5)
fun_plot2(mat_tcfl5_pachy,"tcfl5_pachy_allreads.pdf",10)
fun_plot2(mat_neg_pachy,"neg_pachy_allreads.pdf",5)
fun_plot2(mat_nfya_prepachy,"nfya_prepachy_allreads.pdf",5)

# conserved peaks Jan22----
# library
library(VennDiagram)
setwd("/Users/BigBear/My_Publish/Human_piRNA/conserved_peaks/")
# vennplot
three_mo_1=as.vector(read.table("3MO-1-AMYB_peaks.list",header=F,row.names=NULL)[,1])
three_mo_2=as.vector(read.table("3MO-2-AMYB_peaks.list",header=F,row.names=NULL)[,1])
four_mo=as.vector(read.table("4MO-AMYB_peaks.list",header=F,row.names=NULL)[,1])
cutrun=as.vector(read.table("Pach-Amyb-Cut-Run-ID11_rep1_peaks.list",header=F,row.names=NULL)[,1])
xin=as.vector(read.table("mouse_adult_amyb_peaks.list",header=F,row.names=NULL)[,1])
human_amyb=as.vector(read.table("HUMAN_AMYB_peaks.hg38.mm10.list",header=F,row.names=1))
human_nfya=as.vector(read.table("Human_NFYA_peaks.hg38.list",header=F,row.names=NULL)[,1])
k=as.vector(human_amyb[which(human_amyb[,2]!="none"),2])

pdf("../figures/Jan19.2018.conservedPeaks/pek_overlap_amyb.pdf",width=4,height=4)
par(mfrow=c(1,1),mar=c(1,1,1,1))
plot.new()
grid.draw(venn.diagram(list("NFYA"=human_nfya,"AMYB"=as.vector(row.names(human_amyb))),
                       filename=NULL,fill=csl[c(2,5)],alpha=c(0.5,0.5),margin=0.1))
text(0.5,1,label="human_AMYB\nvs human_NFYA",pos=1)
g=c()
eiTogn=read.table("human_EID_GN.txt",header=F,row.names=1)
for(i in intersect(human_nfya,as.vector(row.names(human_amyb)))){
  g=c(g,as.character(eiTogn[i,1]))
}
write.table(g,"../figures/Jan19.2018.conservedPeaks/human_AMYB_vs_human_NFYA.txt",
            row.names=F,col.names=F,sep="\t",quote=F)
plot.new()
text(0.5,0.5,label="AMYB in\nmouse and human",font=2,cex=1.5)
plot.new()
grid.draw(venn.diagram(list("mouse"=three_mo_2,"human"=k),
                       filename=NULL,fill=csl[c(2,5)],alpha=c(0.5,0.5),margin=0.1))
text(0.5,1,label="mouse_deniz_3MO_2\nvs human_deniz",pos=1)
g=c()
for(i in intersect(k,three_mo_2)){
  g=c(g,as.character(human_amyb[which(human_amyb[,2]==i),1]))
}
write.table(g,"../figures/Jan19.2018.conservedPeaks/human_AMYB_vs_mouse_AMYB_3MO_2.txt",
            row.names=F,col.names=F,sep="\t",quote=F)
plot.new()
grid.draw(venn.diagram(list("mouse"=three_mo_1,"human"=k),
                       filename=NULL,fill=csl[c(2,5)],alpha=c(0.5,0.5),margin=0.1))
text(0.5,1,label="mouse_deniz_3MO_1\nvs human_deniz",pos=1)
g=c()
for(i in intersect(k,three_mo_1)){
  g=c(g,as.character(human_amyb[which(human_amyb[,2]==i),1]))
}
write.table(g,"../figures/Jan19.2018.conservedPeaks/human_AMYB_vs_mouse_AMYB_3MO_1.txt",
            row.names=F,col.names=F,sep="\t",quote=F)
plot.new()
grid.draw(venn.diagram(list("mouse"=four_mo,"human"=k),
                       filename=NULL,fill=csl[c(2,5)],alpha=c(0.5,0.5),margin=0.1))
text(0.5,1,label="mouse_deniz_4MO\nvs human_deniz",pos=1)
g=c()
for(i in intersect(k,four_mo)){
  g=c(g,as.character(human_amyb[which(human_amyb[,2]==i),1]))
}
write.table(g,"../figures/Jan19.2018.conservedPeaks/human_AMYB_vs_mouse_AMYB_4MO.txt",
            row.names=F,col.names=F,sep="\t",quote=F)
plot.new()
grid.draw(venn.diagram(list("mouse"=xin,"human"=k),
                       filename=NULL,fill=csl[c(2,5)],alpha=c(0.5,0.5),margin=0.1))
text(0.5,1,label="mouse_xin_42week\nvs human_deniz",pos=1)
g=c()
for(i in intersect(k,xin)){
  g=c(g,as.character(human_amyb[which(human_amyb[,2]==i),1]))
}
write.table(g,"../figures/Jan19.2018.conservedPeaks/human_AMYB_vs_mouse_AMYB_Xin_6Week.txt",
            row.names=F,col.names=F,sep="\t",quote=F)
plot.new()
grid.draw(venn.diagram(list("mouse"=cutrun,"human"=k),
                       filename=NULL,fill=csl[c(2,5)],alpha=c(0.5,0.5),margin=0.1))
text(0.5,1,label="mouse_deniz_PachyCell_CUTRUN\nvs human_deniz",pos=1)
g=c()
for(i in intersect(k,cutrun)){
  g=c(g,as.character(human_amyb[which(human_amyb[,2]==i),1]))
}
write.table(g,"../figures/Jan19.2018.conservedPeaks/human_AMYB_vs_mouse_AMYB_PachyCell_CUTRUN.txt",
            row.names=F,col.names=F,sep="\t",quote=F)
dev.off()

# peak distance
setwd("/Users/BigBear/My_Publish/Human_piRNA/conserved_peaks/")
# for human
pd=read.table("Human_NFYA_peaks_distance.txt",header=F,row.names=1)
pd1=read.table("Human_AMYB_peaks_distance.txt",header=F,row.names=1)
pl=c("DDX4","DDX39A","FKBP6","GPAT2","HENMT1","MAEL","MOV10L1","MYBL1","NFYA",
     "PIWIL1","PIWIL2","PIWIL4","PLD6","TDRD1","TDRD3","TDRD5","TDRD6","TDRD9",
     "TDRD12","TDRKH")
# NFYA
l1=pd[which(pd[,2]=="protein_coding"),3]
l2=c()
for(i in pl){l2=c(l2,pd[which(pd[,1]==i),3])}
l3=pd[which(pd[,2]=="pachytene"),3]
l4=pd[which(pd[,2]=="prepachytene"),3]
l5=pd[which(pd[,2]=="hybrid"),3]
pdf("peak_distance.nfya.pdf",width=4,height=4)
par(mar=c(4,4,1,1),cex=5/6,bty="n",xpd=T)
plot(NA,xlim=c(0,9),ylim=c(0,7),ylab="distance (kbp)",xlab="",xaxt="n",yaxt="n")
ph=20
boxplot(log10(l1+1),log10(l2+1),log10(l3+1),log10(l4+1),log10(l5+1),outline=F,
        at=c(1,3,5,7,9)-0.5,boxwex=1,staplewex=0,add=T,yaxt="n",xaxt="n")
#points(runif(length(l1),0+0.05,1-0.05),log10(l1+1),col=csl[1],pch=ph,lwd=1,cex=0.2)
points(runif(length(l2),2+0.05,3-0.05),log10(l2+1),col=csl[2],pch=ph,lwd=1,cex=0.5)
points(runif(length(l3),4+0.05,5-0.05),log10(l3+1),col=csl[3],pch=ph,lwd=1,cex=0.5)
points(runif(length(l4),6+0.05,7-0.05),log10(l4+1),col=csl[4],pch=ph,lwd=1,cex=0.5)
points(runif(length(l5),8+0.05,9-0.05),log10(l5+1),col=csl[5],pch=ph,lwd=1,cex=0.5)
text(c(1,3,5,7,9)-0.5,-1,label=c("mRNA genes","piRNA\nbiogenesis","pachytene\nclusters",
                                 "prepachytene\nclusters","hybrid\nclusters"),xpd=T,
     srt=45,col=csl[1:5])
axis(2,0:7,c(expression(10^0,10^1,10^2,10^3,10^4,10^5,10^6,10^7)))
dev.off()
# AMYB
l1=pd1[which(pd1[,2]=="protein_coding"),3]
l2=c()
for(i in pl){l2=c(l2,pd1[which(pd1[,1]==i),3])}
l3=pd1[which(pd1[,2]=="pachytene"),3]
l4=pd1[which(pd1[,2]=="prepachytene"),3]
l5=pd1[which(pd1[,2]=="hybrid"),3]
pdf("peak_distance.amyb.pdf",width=4,height=4)
par(mar=c(4,4,1,1),cex=5/6,bty="n",xpd=T)
plot(NA,xlim=c(0,9),ylim=c(0,7),ylab="distance (kbp)",xlab="",xaxt="n",yaxt="n")
ph=20
boxplot(log10(l1+1),log10(l2+1),log10(l3+1),log10(l4+1),log10(l5+1),outline=F,
        at=c(1,3,5,7,9)-0.5,boxwex=1,staplewex=0,add=T,yaxt="n",xaxt="n")
#points(runif(length(l1),0+0.05,1-0.05),log10(l1+1),col=csl[1],pch=ph,lwd=1,cex=0.2)
points(runif(length(l2),2+0.05,3-0.05),log10(l2+1),col=csl[2],pch=ph,lwd=1,cex=0.5)
points(runif(length(l3),4+0.05,5-0.05),log10(l3+1),col=csl[3],pch=ph,lwd=1,cex=0.5)
points(runif(length(l4),6+0.05,7-0.05),log10(l4+1),col=csl[4],pch=ph,lwd=1,cex=0.5)
points(runif(length(l5),8+0.05,9-0.05),log10(l5+1),col=csl[5],pch=ph,lwd=1,cex=0.5)
text(c(1,3,5,7,9)-0.5,-1,label=c("mRNA genes","piRNA\nbiogenesis","pachytene\nclusters",
                                 "prepachytene\nclusters","hybrid\nclusters"),xpd=T,
     srt=45,col=csl[1:5])
axis(2,0:7,c(expression(10^0,10^1,10^2,10^3,10^4,10^5,10^6,10^7)))
dev.off()
# for mouse
pl=c("Ddx4","Ddx39a","Fkbp6","Gpat2","Henmt1","Mael","Mov10l1","Mybl1","Nfya",
     "Piwil1","Piwil2","Piwil4","Pld6","Tdrd1","Tdrd3","Tdrd5","Tdrd6","Tdrd9",
     "Tdrd12","Tdrkh")
fun_plot_dis=function(pd, fn){
  par(mar=c(4,4,1,1),cex=5/6,bty="n",xpd=T,tcl=0.3)
  plot(NA,xlim=c(0,9),ylim=c(0,7),ylab="distance (kbp)",xlab="",xaxt="n",yaxt="n")
  ph=20
  x1=log10(pd[proteincoding,]+1); x2=log10(pd[pl,]+1);
  x3=log10(pd[Pachy,]+1); x4=log10(pd[Prepachy,]+1); x5=log10(pd[Hybrid,]+1)
  boxplot(x1,x2,x3,x4,x5,outline=F,
          at=c(1,3,5,7,9)-0.5,boxwex=1,staplewex=0,add=T,yaxt="n",xaxt="n")
  #points(runif(length(l1),0+0.05,1-0.05),log10(l1+1),col=csl[1],pch=ph,lwd=1,cex=0.2)
  points(runif(length(x2),2+0.05,3-0.05),x2,col=csl[2],pch=ph,lwd=1,cex=0.5)
  points(runif(length(x3),4+0.05,5-0.05),x3,col=csl[3],pch=ph,lwd=1,cex=0.5)
  points(runif(length(x4),6+0.05,7-0.05),x4,col=csl[4],pch=ph,lwd=1,cex=0.5)
  points(runif(length(x5),8+0.05,9-0.05),x5,col=csl[5],pch=ph,lwd=1,cex=0.5)
  text(c(1,3,5,7,9)-0.5,-1,label=c("mRNA genes","piRNA\nbiogenesis","pachytene\nclusters",
                                   "prepachytene\nclusters","hybrid\nclusters"),xpd=T,
       srt=45,col=csl[1:5])
  text(4.5,7,label=fn,font=2,cex=1.5)
  axis(2,0:7,c(expression(10^0,10^1,10^2,10^3,10^4,10^5,10^6,10^7)))
}
pdf("mouse_peaks_distance.pdf",width=4,height=4)
pd=read.table("3MO-2-AMYB_peaks.distance.txt",header=F,row.names=1)
fun_plot_dis(pd, "AMYB-3MO")
pd=read.table("NFYA_3MO-2.distance.txt",header=F,row.names=1)
fun_plot_dis(pd, "NFYA-3MO-2")
pd=read.table("NFYA_3MO-3.distance.txt",header=F,row.names=1)
fun_plot_dis(pd, "NFYA-3MO-3")
pd=read.table("NFYA_3MO-4.distance.txt",header=F,row.names=1)
fun_plot_dis(pd, "NFYA-3MO-4")
pd=read.table("NFYA_4MO.distance.txt",header=F,row.names=1)
fun_plot_dis(pd, "NFYA-4MO")
pd=read.table("NFYA_kidney.distance.txt",header=F,row.names=1)
fun_plot_dis(pd, "NFYA-kidney")
pd=read.table("NFYA_spleen.distance.txt",header=F,row.names=1)
fun_plot_dis(pd, "NFYA-spleen")
pd=read.table("TCFL5_2MO.distance.txt",header=F,row.names=1)
fun_plot_dis(pd, "TCFL5-2MO")
pd=read.table("TCFL5_3MO.distance.txt",header=F,row.names=1)
fun_plot_dis(pd, "TCFL5-3MO")
pd=read.table("TCFL5_4MO.distance.txt",header=F,row.names=1)
fun_plot_dis(pd, "TCFL5-4MO")
dev.off()

# TCFL5, AMYB, NFYA for marmoset and rat----
setwd("/Users/BigBear/My_Publish/Human_piRNA/aggreMatrix/")
pl=c("DDX4","DDX39A","FKBP6","GPAT2","HENMT1","MAEL","MOV10L1","MYBL1","NFYA",
     "PIWIL1","PIWIL2","PIWIL4","PLD6","TDRD1","TDRD3","TDRD5","TDRD6","TDRD9",
     "TDRKH")
pl1=c("Ddx4","Ddx39a","Fkbp6","Gpat2","Henmt1","Mael","Mov10l1","Mybl1","Nfya",
      "Piwil1","Piwil2","Piwil4","Pld6","Tdrd1","Tdrd3","Tdrd5","Tdrd6","Tdrd9",
      "Tdrkh")
rat_tcfl5=read.table("Rat-TCFL5-ID6_S1.rpm.uniq.mat.t5.u5000.d5000.b200.mat",header=F,row.names=1)
rat_amyb=read.table("Rat-AMYB-ID5_S7.rpm.uniq.mat.t5.u5000.d5000.b200.mat",header=F,row.names=1)
rat_nfya=read.table("Rat-NFYA-ID7_S3.rpm.uniq.mat.t5.u5000.d5000.b200.mat",header=F,row.names=1)
rat_input=read.table("Rat-Input-ID8_S5.rpm.uniq.mat.t5.u5000.d5000.b200.mat",header=F,row.names=1)
marmoset_tcfl5=read.table("Marmoset-TCFL5-ID2_S4.rpm.uniq.mat.t5.u5000.d5000.b200.mat",header=F,row.names=1)
marmoset_amyb=read.table("Marmoset-AMYB-ID1_S8.rpm.uniq.mat.t5.u5000.d5000.b200.mat",header=F,row.names=1)
marmoset_nfya=read.table("Marmoset-NFYA-ID3_S2.rpm.uniq.mat.t5.u5000.d5000.b200.mat",header=F,row.names=1)
marmoset_input=read.table("Marmoset-Input-ID4_S6.rpm.uniq.mat.t5.u5000.d5000.b200.mat",header=F,row.names=1)

rn=c("-5kb",rep("",99),"TSS",rep("",98),"+5kb")
pdf("../chipseq_heatmap/rat.chip.pdf",width=4,height=3)
pheatmap(rat_tcfl5[pl1,],cluster_rows=F,cluster_col=F,main="tcfl5",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(rat_tcfl5[pl1,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(rat_amyb[pl1,],cluster_rows=F,cluster_col=F,main="amyb",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(rat_amyb[pl1,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(rat_nfya[pl1,],cluster_rows=F,cluster_col=F,main="nfya",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(rat_nfya[pl1,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(rat_input[pl1,],cluster_rows=F,cluster_col=F,main="input",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(rat_input[pl1,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
dev.off()

pdf("../chipseq_heatmap/marmoset.chip.pdf",width=4,height=3)
pheatmap(marmoset_tcfl5[pl1,],cluster_rows=F,cluster_col=F,main="tcfl5",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(marmoset_tcfl5[pl1,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(marmoset_amyb[pl1,],cluster_rows=F,cluster_col=F,main="amyb",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(marmoset_amyb[pl1,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(marmoset_nfya[pl1,],cluster_rows=F,cluster_col=F,main="nfya",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(marmoset_nfya[pl1,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(marmoset_input[pl1,],cluster_rows=F,cluster_col=F,main="input",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(marmoset_input[pl1,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
dev.off()

# Dec CUTRUN correlation----
setwd("/Users/BigBear/My_Publish/Human_piRNA/CUT_RUN/Oct.2017/")
rpm=read.table("RPM.txt",header=T,row.names=1)*1000
fun_plot=function(i,j){
  par(bty="n",tcl=0.3)
  plot(-100,xlim=c(0,4),ylim=c(0,4),xlab=colnames(rpm)[i],ylab=colnames(rpm)[j],
       main="log10 ( RPKM around TSS )")
  abline(0,1); abline(-log10(2),1,lty=2); abline(log10(2),1,lty=2)
  points(log10(rpm[Hybrid,i]+1),log10(rpm[Hybrid,j]+1),lwd=1.5,col=csl[4])
  points(log10(rpm[Prepachy,i]+1),log10(rpm[Prepachy,j]+1),lwd=1.5,col=csl[3])
  points(log10(rpm[Pachy,i]+1),log10(rpm[Pachy,j]+1),lwd=1.5,col=csl[1])
  cor1=round(cor(rpm[Pachy,i],rpm[Pachy,j],method="spearman"),2)
  cor2=round(cor(rpm[Prepachy,i],rpm[Prepachy,j],method="spearman"),2)
  cor3=round(cor(rpm[Hybrid,i],rpm[Hybrid,j],method="spearman"),2)
  text(c(0,0,0),c(4,3.7,3.4),
       label=paste(c("pachytene","prepachytene","hybrid"),c(cor1,cor2,cor3),sep=" "),
       col=csl[c(1,3,4)],pos=4)
}
pdf("../../figures/Oct24.2017.cutrun2/coorelation.pdf",width=8,height=8)
par(mar=c(4,4,4,2),mfrow=c(2,2))
fun_plot(1,2)
fun_plot(1,3)
fun_plot(2,3)
fun_plot(2,4)
dev.off()
Pachy[which(rpm[Pachy,1]>1000 & rpm[Pachy,3]<1000)]
Pachy[which(rpm[Pachy,3]>1000 & rpm[Pachy,1]<1000)]

# Feb15.2018, motif conservation ------
setwd("/Users/BigBear/My_Publish/Human_piRNA/motifs/motif_score/")
amyb=read.table("myba.mm10.in_gene+pi.result",header=F,row.names=1)
tcfl5=read.table("tcfl5.mm10.in_gene+pi.result",header=F,row.names=1)
nfy=read.table("nfy.mm10.in_gene+pi.result",header=F,row.names=1)
sig_amyb=read.table("3MO-2-AMYB.FE",header=F,row.names=1)
sig_tcfl5=read.table("TCFL5_3MO_rep1.FE",header=F,row.names=1)
sig_nfya=read.table("NFYA_3MO-3_rep1.FE",header=F,row.names=1)
sig_mat=cbind(sig_amyb,sig_tcfl5,sig_nfya)
colnames(sig_mat)=c("AMYB 3MO","TCFL5 3MO","NFYA 3MO")
amyb_pachy=intersect(row.names(amyb),Pachy)
tcfl5_pachy=intersect(row.names(tcfl5),Pachy)
nfy_pachy=intersect(row.names(nfy),Pachy)
motif_score=matrix(0,46676,3);row.names(motif_score)=row.names(sig_mat);colnames(motif_score)=c("amyb","tcfl5","nfya")
motif_score[row.names(amyb),1]=amyb[,1]
motif_score[row.names(tcfl5),2]=tcfl5[,1]
motif_score[row.names(nfy),3]=nfy[,1]

# scatterplot
pdf("correlation_chipseq_3month.pdf",width=5,height=5,useDingbats=F)
par(mar=c(4,4,1,1),bty="n")
fun_scat=function(i,j){
  lb=expression(2^0,2^1,2^2,2^3,2^4,2^5,2^6,2^7,2^8)
  x=log2(sig_mat+1)
  plot(x[Pachy,i],x[Pachy,j],lwd=2,col=csl[1],
       xlim=c(0,max(x[Pi,i])),ylim=c(0,max(x[Pi,j])),
       xlab=colnames(sig_mat)[i],ylab=colnames(sig_mat)[j],xaxt="n",yaxt="n")
  axis(1,0:8,label=lb);axis(2,0:8,label=lb)
  points(x[Prepachy,i],x[Prepachy,j],lwd=2,col=csl[2])
  points(x[Hybrid,i],x[Hybrid,j],lwd=2,col=csl[3])
  cr=c(cor(x[Pachy,i],x[Pachy,j]),cor(x[Prepachy,i],x[Prepachy,j]),
       cor(x[Pachy,i],x[Pachy,j]))
  cr=round(cr,2)
  text(max(x[Pi,i]),max(x[Pi,j]),pos=2,label=paste("cor = ",cr[1],sep=""),col=csl[1],font=2)
  text(max(x[Pi,i]),max(x[Pi,j])*9.3/10,pos=2,label=paste("cor = ",cr[2],sep=""),col=csl[2],font=2)
  text(max(x[Pi,i]),max(x[Pi,j])*8.6/10,pos=2,label=paste("cor = ",cr[3],sep=""),col=csl[3],font=2)
}
fun_scat(1,2)
fun_scat(2,3)
fun_scat(1,3)
dev.off()

# barplot
pn=c(setdiff(Pachy,c(amyb_pachy,tcfl5_pachy,nfy_pachy)),tcfl5_pachy,
     setdiff(nfy_pachy,tcfl5_pachy),setdiff(amyb_pachy,c(tcfl5_pachy,nfy_pachy)))[100:1]
pdf("barplot.pdf",width=8,height=1.3)
par(mar=c(0,0,0,0))
mat=as.matrix(log2(sig_mat[Pachy,]+1))
for(i in 1:3){mat[Pachy,i]=mat[Pachy,i]/max(mat[Pachy,i])}
barplot(t(mat[pn,]),
        beside=T,col=csl[1:3],border=csl[1:3],yaxt="n",xaxt="n")
dev.off()
# heatmap2
mat=as.matrix(log2(sig_mat[Pachy,]+1))
for(i in 1:3){mat[Pachy,i]=mat[Pachy,i]/max(mat[Pachy,i])}
pheatmap(t(mat[pn,]),cluster_rows=F,cluster_cols=F,labels_row=c("AMYB","TCFL5","NFYA"),
        filename="heatmap2.pdf",cellwidth=8,cellheight=24)
# heatmap
motif_score=matrix(0,100,3)
row.names(motif_score)=pn
motif_score[tcfl5_pachy,2]=tcfl5[tcfl5_pachy,1]
motif_score[amyb_pachy,1]=amyb[amyb_pachy,1]
motif_score[nfy_pachy,3]=nfy[nfy_pachy,1]
pheatmap(t(motif_score),
         cluster_rows=F,cluster_cols=F,
         color=c("grey",colorRampPalette(brewer.pal(9,"Reds"))(100)),
         breaks=c(0,1:100/15+4.5),
         labels_row=c("AMYB","TCFL5","NFYA"),fontsize_col=7,
         filename="heatmap.pdf",cellwidth=8,cellheight=24)

pdf("AMYB_score.pdf",width=3.5,height=4)
par(tcl=0.3,bty="n",cex=5/6)
x1=intersect(amyb_pachy,tcfl5_pachy)
x2=intersect(amyb_pachy,nfy_pachy)
x3=setdiff(amyb_pachy,c(tcfl5_pachy,nfy_pachy))
boxplot(motif_score[x1,1],motif_score[x2,1],motif_score[x3,1],ylab="motif score",
        xaxt="n",main="89 AMYB promoters\np-value = 0.03")
points(runif(length(x1),0.6,1.4),motif_score[x1,1],pch=19,col=csl[1])
points(runif(length(x2),1.6,2.4),motif_score[x2,1],pch=19,col=csl[2])
points(runif(length(x3),2.6,3.4),motif_score[x3,1],pch=19,col=csl[3])
axis(1,1:3,label=c("with TCFL5\npromoters","with NFYA\npromoters",
                   "AMYB-specific\npromoters"),lwd=0,cex.axis=0.7)
dev.off()
pdf("TCFL5_score.pdf",width=3,height=4)
par(tcl=0.3,bty="n",cex=5/6)
x1=tcfl5_pachy[which(amyb[tcfl5_pachy,]>=8)]
x2=setdiff(tcfl5_pachy,x1)
boxplot(motif_score[x1,2],motif_score[x2,2],ylab="motif score",
        xaxt="n",main="37 TCFL5 promoters\np-value = 0.25")
points(runif(length(x1),0.6,1.4),motif_score[x1,2],pch=19,col=csl[1])
points(runif(length(x2),1.6,2.4),motif_score[x2,2],pch=19,col=csl[3])
axis(1,1:2,label=c("with high\nAMYB score\n>=8","with low\nAMYB score\n<8"),lwd=0,cex.axis=0.7)
dev.off()

# Feb22.2018, ChIPseq Heatmap and Aggregation plot----
setwd("/Users/BigBear/My_Publish/Human_piRNA/chipseq_heatmap/")
am_amyb=read.table("AMYB_3MO.aggre.mat",header=F,row.names=1)
am_tcfl5=read.table("TCFL5_3MO.aggre.mat",header=F,row.names=1)
am_nfya=read.table("NFYA_3MO.aggre.mat",header=F,row.names=1)
am=list(AMYB=am_amyb,TCFL5=am_tcfl5,NFYA=am_nfya)

fun_plot=function(i,fn,mn){
  pdf(fn,width=4,height=8)
  m=as.matrix(log2(i+1))
  row.names(m)=rep("",dim(m)[1])
  colnames(m)=c("-3kb",rep("",49),"TSS",rep("",98),"TES",rep("",49),"3kb")
  par(mar=c(2,2,2,2),bty="n",mfrow=c(3,1))
  plot(apply(i,2,mean),type="l",col="#313695",lwd=2,ylab="",xlab="",xaxt="n")
  mtext(mn,side=3,font=2,cex=1.2)
  axis(1,c(0,50,150,200),label=c("-3kb","TSS","TES","3kb"))
  pheatmap(m[o,],cluster_rows=F,cluster_cols=F,
           colorRampPalette(brewer.pal(11,"RdYlBu"))(100))
  dev.off()
}
o=order(apply(am_amyb[Pachy,],1,sum),decreasing=T)
fun_plot(am_amyb[Pachy,],"amyb.pachy.pdf","mouse_AMYB")
fun_plot(am_tcfl5[Pachy,],"tcfl5.pachy.pdf","mouse_TCFL5")
fun_plot(am_nfya[Pachy,],"nfya.pachy.pdf","mouse_NFYA")
o=order(apply(am_tcfl5[Prepachy,],1,sum),decreasing=T)
fun_plot(am_amyb[Prepachy,],"amyb.prepachy.pdf","mouse_AMYB")
fun_plot(am_tcfl5[Prepachy,],"tcfl5.prepachy.pdf","mouse_TCFL5")
fun_plot(am_nfya[Prepachy,],"nfya.prepachy.pdf","mouse_NFYA")
o=order(apply(am_tcfl5[Hybrid,],1,sum),decreasing=T)
fun_plot(am_amyb[Hybrid,],"amyb.hybrid.pdf","mouse_AMYB")
fun_plot(am_tcfl5[Hybrid,],"tcfl5.hybrid.pdf","mouse_TCFL5")
fun_plot(am_nfya[Hybrid,],"nfya.hybrid.pdf","mouse_NFYA")

# Mar08.2018, CAGE and PAS Heatmap and Aggregation plot, also for H3K4me3----
# read table
setwd("/Users/BigBear/My_Publish/Human_piRNA/aggreMatrix/")
mat_cage=list()
mat_cage$CAGE_1643_Rep1_Id2=read.table("CAGE_1643_Rep1_Id2.t53.u5000.d5000.b200.mat",header=F,row.names=1)
mat_cage$CAGE_1643_Rep2_Id3=read.table("CAGE_1643_Rep2_Id3.t53.u5000.d5000.b200.mat",header=F,row.names=1)
mat_cage$CAGE_1643_Rep3_Id1=read.table("CAGE_1643_Rep3_Id1.t53.u5000.d5000.b200.mat",header=F,row.names=1)
mat_cage$CAGE_5431_Id7=read.table("CAGE_5431_Id7.t53.u5000.d5000.b200.mat",header=F,row.names=1)
mat_cage$CAGE_7358_Id6=read.table("CAGE_7358_Id6.t53.u5000.d5000.b200.mat",header=F,row.names=1)
mat_cage$CAGE_7377_Id5=read.table("CAGE_7377_Id5.t53.u5000.d5000.b200.mat",header=F,row.names=1)
mat_cage$CAGE_8150_Id4=read.table("CAGE_8150_Id4.t53.u5000.d5000.b200.mat",header=F,row.names=1)
mat_pas=list()
mat_pas$ACTGT.Dedup_Atrimmed=read.table("ACTGT.Dedup_Atrimmed.t53.u5000.d5000.b200.mat",header=F,row.names=1)
mat_pas$ATCTG.Dedup_Atrimmed=read.table("ATCTG.Dedup_Atrimmed.t53.u5000.d5000.b200.mat",header=F,row.names=1)
mat_pas$CGGGA.Dedup_Atrimmed=read.table("CGGGA.Dedup_Atrimmed.t53.u5000.d5000.b200.mat",header=F,row.names=1)
mat_pas$TGACT.Dedup_Atrimmed=read.table("TGACT.Dedup_Atrimmed.t53.u5000.d5000.b200.mat",header=F,row.names=1)
# piRNA gene names
human_pi=row.names(mat_cage$CAGE_5431_Id7)
human_pachy=human_pi[grep("^pachy",human_pi)]
#human_pachy=setdiff(human_pachy,"pachytene_78")
human_prepachy=human_pi[grep("^prepa",human_pi)]
human_hybrid=human_pi[grep("^hybr",human_pi)]
human_pi=c(human_pachy,human_hybrid,human_prepachy)
annoRow=data.frame(Type=factor(rep(c("pachytene","hybrid","prepachytene"),c(155,17,78))))
row.names(annoRow)=human_pi
annoColor=list(Type=c(pachytene=csl[1],hybrid=csl[3],prepachytene=csl[2]))
# function
fun_rownor=function(x){
  if(min(x)==max(x)){
    return(rep(0,length(x)))
    }
  else{
    return((x-min(x))/(max(x)-min(x)))
  }
}
fun_plot=function(i,fn,mn,ma){
  # pdf(fn,width=4,height=8)
  m=as.matrix(log10(i+1))
  for(i in 1:dim(m)[1]){
    m[i,which(m[i,]>ma)]=ma
  }
  # m=as.matrix(i)
  # row.names(m)=rep("",dim(m)[1])
  colnames(m)=c("-5kb",rep("",49),"TSS",rep("",98),"TES",rep("",49),"5kb")
  # par(mar=c(2,2,2,2),bty="n",mfrow=c(3,1))
  # plot(apply(i,2,mean),type="l",col="#313695",lwd=2,ylab="",xlab="",xaxt="n")
  # mtext(mn,side=3,font=2,cex=1.2)
  # axis(1,c(0,50,150,200),label=c("-3kb","TSS","TES","3kb"))
  pheatmap(m,cluster_rows=F,cluster_cols=F,
           colorRampPalette(brewer.pal(11,"RdYlBu"))(100),
           cellwidth=1.5,cellheight=4,filename=fn,main=mn,
           fontsize_row=4,annotation_row=annoRow,annotation_colors=annoColor)
  dev.off()
}
# plot cage
i="CAGE_7358_Id6"
fun_plot(mat_cage[[i]][human_pi,],paste("human_pi",i,"pdf",sep="."),i,1.5)
i="CAGE_8150_Id4"
fun_plot(mat_cage[[i]][human_pi,],paste("human_pi",i,"pdf",sep="."),i,1.5)
i="TGACT.Dedup_Atrimmed"
fun_plot(mat_pas[[i]][human_pi,],paste("human_pi",i,"pdf",sep="."),i,1.5)

mat_h3k4me3=read.table("human_H3K4me3_testis.t53.u5000.d5000.b200.mat",header=F,row.names=1)
fun_plot(mat_h3k4me3[human_pi,],"human_pi_h3k4me3.pdf","human_H3K4me3_testis",0.7)

# Mar12.2018, srna in 5utr, 3utr, cds, intron Heatmap----
setwd("/Users/BigBear/My_Publish/Human_piRNA/srna_across_prepachy/")
srna_list=list()
srna_list$juvenile_1=read.table("hg19.prepachytene.1643-5YO-Rep1-Ox.mat",header=T,row.names=1)
srna_list$juvenile_2=read.table("hg19.prepachytene.1643-5YO-Rep2-Ox.mat",header=T,row.names=1)
srna_list$N2692=read.table("hg19.prepachytene.2692_NT.mat",header=T,row.names=1)
srna_list$N5125=read.table("hg19.prepachytene.5125-Ox.mat",header=T,row.names=1)
srna_list$N7358=read.table("hg19.prepachytene.7358_NT.mat",header=T,row.names=1)
srna_list$N7377=read.table("hg19.prepachytene.7377-Ox.mat",header=T,row.names=1)
srna_list$N8150=read.table("hg19.prepachytene.8150_NT-Ox.mat",header=T,row.names=1)
pdf("summary.reads.pdf",width=2.5,height=5)
for(i in names(srna_list)){
  pheatmap(log10(srna_list[[i]][,c(3,5,7,1)]+1),main=i,
           cluster_rows=F,cluster_cols=F,fontsize_row=6)
}
dev.off()
pdf("summary.readsPer1k.pdf",width=2.5,height=5)
for(i in names(srna_list)){
  pheatmap(log10(srna_list[[i]][,c(4,6,8,2)]+1),main=i,
           cluster_rows=F,cluster_cols=F,fontsize_row=6)
}
dev.off()

# Mar23.2018, picluster comparison----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/Bear_Mar23.2018_piclusterComparison/")
sm=read.table("summary.txt",header=T,row.names=1)
pdf("summary.pdf",width=8,height=5)
layout(matrix(c(1,2,2,2,2,1,3,3,3,3),5,2))
par(mar=c(0,0,0,0))
plot.new()
legend("center",legend=c("shared","specific","pachytene","prepachytene","hybrid"),
       density=c(-1,20,-1,-1,-1),ncol=5,fill=c("black","black",csl[1:3]))
par(mar=c(2,2,2,2))
pie(c(sm[1:3,1],sm[1:3,2]),border=F,col=csl[1:3],
    density=c(20,20,20,-1,-1,-1),label=c(sm[1:3,1],sm[1:3,2]),
    main="Deniz's Annotation")
pie(as.integer(sm[4,]),border=F,col="black",
    density=c(20,-1),label=sm[4,],
    main="Previous Annotation")
dev.off()

# Mar23.2018, tpm of different piclusters----
setwd("/Users/BigBear/My_Publish/Human_piRNA/tpm/")
st=read.table("genes_expression_tpm.tsv",header=T,row.names=1)
i30nt=c("D8150NT","D7358NT","S7_NT","S12_NT")
i26nt=c("S6_NT","D7377NT","AmbionNT","S8_NT")
human_piRNA_bed=read.table("../piRNA_list/human_piRNA.hg19.bed9",header=F,row.names=4)
piN_STRG_map=read.table("../piRNA_list/piName_STRG.map",header=F,row.names=1)
strg.pi.uni=as.vector(piN_STRG_map[row.names(human_piRNA_bed)[which(human_piRNA_bed[,6]=="unidirectional")],1])
strg.pi.bi=as.vector(piN_STRG_map[row.names(human_piRNA_bed)[which(human_piRNA_bed[,6]=="bidirectional")],1])
strg.pi.genic=as.vector(piN_STRG_map[row.names(human_piRNA_bed)[which(human_piRNA_bed[,7]=="genic")],1])
strg.pi.inter=as.vector(piN_STRG_map[row.names(human_piRNA_bed)[which(human_piRNA_bed[,7]=="intergenic")],1])
strg.pachy=as.vector(piN_STRG_map[row.names(human_piRNA_bed)[grep("^pachy",row.names(human_piRNA_bed))],1])
strg.pi.uni2=intersect(strg.pi.uni,strg.pachy)
strg.pi.bi2=intersect(strg.pi.bi,strg.pachy)
strg.pi.genic2=intersect(strg.pi.genic,strg.pachy)
strg.pi.inter2=intersect(strg.pi.inter,strg.pachy)
# function
fun_plot=function(n1,l1,n2,l2,yl){
  layout(matrix(c(1,2,2,2,2,2,2),7,1))
  par(mar=c(0,0,0,0))
  plot.new()
  legend("center",
         legend=c(paste(l1,",n=",length(n1),sep=""),paste(l2,",n=",length(n2),sep="")),
         fill=csl[1:2],ncol=2,bty="n")
  par(mar=c(4,4,0.5,0.5),tcl=0.3,bty="n")
  boxplot(log10(st[n1,c(i30nt,i26nt)]+1),xlim=c(0.5,8.5),col=csl[1],
          at=1:8-0.2,boxwex=0.4,xaxt="n",yaxt="n",ylab=yl)
  boxplot(log10(st[n2,c(i30nt,i26nt)]+1),xlim=c(0,8),col=csl[2],
          at=1:8+0.2,boxwex=0.4,xaxt="n",yaxt="n",add=T)
  axis(1,1:8,label=c(i30nt,i26nt),lwd=0,cex.axis=0.8)
  axis(2,0:4,label=expression(10^0,10^1,10^2,10^3,10^4))
  pv=c()
  for(i in c(i30nt,i26nt)){
    pv=c(pv,signif(wilcox.test(st[n1,i],st[n2,i])$p.value,2))
  }
  df_pv=data.frame(pv=pv,cs="black",stringsAsFactors=F);df_pv[df_pv[,1]<0.05,2]=csl[1]
  text(1:8,2.8,
       pos=3,label=df_pv[,1],col=df_pv[,2],font=2)
}
# plot
pdf("picluster.uni.bi.genic.inter.pdf",width=5,height=4)
fun_plot(strg.pi.uni2,"unidirection",strg.pi.bi2,"bidirection","TPM of pachytene piRNA genes")
fun_plot(strg.pi.genic2,"genic",strg.pi.inter2,"intergenic","TPM of pachytene piRNA genes")
fun_plot(strg.pi.uni,"unidirection",strg.pi.bi,"bidirection","TPM of piRNA genes")
fun_plot(strg.pi.genic,"genic",strg.pi.inter,"intergenic","TPM of piRNA genes")
dev.off()

# May06.2018, TCFL5, AMYB, NFYA for marmoset, rat, mouse and human----
# read data
setwd("/Users/BigBear/My_Publish/Human_piRNA/aggreMatrix/")
pl=c("DDX4","DDX39A","FKBP6","GPAT2","HENMT1","MAEL","MOV10L1","MYBL1","NFYA",
     "PIWIL1","PIWIL2","PIWIL4","PLD6","TDRD1","TDRD3","TDRD5","TDRD6","TDRD9",
     "TDRKH")
pl1=c("Ddx4","Ddx39a","Fkbp6","Gpat2","Henmt1","Mael","Mov10l1","Mybl1","Nfya",
      "Piwil1","Piwil2","Piwil4","Pld6","Tdrd1","Tdrd3","Tdrd5","Tdrd6","Tdrd9",
      "Tdrkh")
pl2=c("DDX4","DDX39A","FKBP6","GPAT2","HENMT1","MAEL","MOV10L1","MYBL1","NFYA",
     "PIWIL1","PIWIL2","PIWIL4","PLD6","TDRD1","TDRD3","TDRD5","TDRD6","TDRD9",
     "TDRKH")
rat_tcfl5=read.table("Rat-TCFL5-ID6_S1.rpm.uniq.mat.t5.u5000.d5000.b200.mat",header=F,row.names=1)
rat_amyb=read.table("Rat-AMYB-ID5_S7.rpm.uniq.mat.t5.u5000.d5000.b200.mat",header=F,row.names=1)
rat_nfya=read.table("Rat-NFYA-ID7_S3.rpm.uniq.mat.t5.u5000.d5000.b200.mat",header=F,row.names=1)
rat_input=read.table("Rat-Input-ID8_S5.rpm.uniq.mat.t5.u5000.d5000.b200.mat",header=F,row.names=1)
marmoset_tcfl5=read.table("Marmoset-TCFL5-ID2_S4.rpm.uniq.mat.t5.u5000.d5000.b200.mat",header=F,row.names=1)
marmoset_amyb=read.table("Marmoset-AMYB-ID1_S8.rpm.uniq.mat.t5.u5000.d5000.b200.mat",header=F,row.names=1)
marmoset_nfya=read.table("Marmoset-NFYA-ID3_S2.rpm.uniq.mat.t5.u5000.d5000.b200.mat",header=F,row.names=1)
marmoset_input=read.table("Marmoset-Input-ID4_S6.rpm.uniq.mat.t5.u5000.d5000.b200.mat",header=F,row.names=1)
mouse_tcfl5=read.table("mouse_3MO-TCFL5-IP-ID02_S2.rpm.uniq.t5.u5000.d5000.b200.mat",header=F,row.names=1)
mouse_tcfl5_input=read.table("mouse_3MO-TCFL5-input-ID04_S6.rpm.uniq.t5.u5000.d5000.b200.mat",header=F,row.names=1)
mouse_nfya=read.table("mouse_3MO-3-NFYA-ID1_S3.rpm.uniq.t5.u5000.d5000.b200.mat",header=F,row.names=1)
mouse_nfya_input=read.table("mouse_3MO-3-NFYA-input-ID1_S3.rpm.uniq.t5.u5000.d5000.b200.mat",header=F,row.names=1)
mouse_amyb=read.table("mouse_3MO-2-AMYB-Id8_S7.rpm.uniq.mat.t5.u5000.d5000.b200.mat",header=F,row.names=1)
human_amyb=read.table("human_AMYB.t5.u5000.d5000.b200.mat",header=F,row.names=1)
human_amyb_input=read.table("human_AMYB_input.t5.u5000.d5000.b200.mat",header=F,row.names=1)
human_nfya=read.table("human_NFYA.t5.u5000.d5000.b200.mat",header=F,row.names=1)
human_nfya_input=read.table("human_NFYA_input.t5.u5000.d5000.b200.mat",header=F,row.names=1)
human_tcfl5_ID2=read.table("5745NT-TCFL5-ID2_S6.rpm.uniq.t5.u5000.d5000.b200.mat",header=F,row.names=1)
human_tcfl5_input_ID4=read.table("5745NT-input-ID4_S5.rpm.uniq.t5.u5000.d5000.b200.mat",header=F,row.names=1)
human_tcfl5_ID1=read.table("5952NT-TCFL5-ID1_S1.rpm.uniq.t5.u5000.d5000.b200.mat",header=F,row.names=1)
human_tcfl5_input_ID3=read.table("5952NT-input-ID3_S3.rpm.uniq.t5.u5000.d5000.b200.mat",header=F,row.names=1)
macaque_tcfl5=read.table("Macaque-TCFL5-ID6_S4.rpm.uniq.t5.u5000.d5000.b200.mat",header=F,row.names=1)
macaque_amyb=read.table("Macaque-AMYB-ID5_S8.rpm.uniq.t5.u5000.d5000.b200.mat",header=F,row.names=1)
macaque_nfya=read.table("Macaque-NFYA-ID7_S2.rpm.uniq.t5.u5000.d5000.b200.mat",header=F,row.names=1)
macaque_input=read.table("Macaque-Input-ID8_S7.rpm.uniq.t5.u5000.d5000.b200.mat",header=F,row.names=1)

# rat
rn=c("-5kb",rep("",99),"TSS",rep("",98),"+5kb")
pdf("../chipseq_heatmap/rat.chip.pdf",width=4,height=3)
pheatmap(rat_tcfl5[pl1,],cluster_rows=F,cluster_col=F,main="tcfl5",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(rat_tcfl5[pl1,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(rat_amyb[pl1,],cluster_rows=F,cluster_col=F,main="amyb",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(rat_amyb[pl1,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(rat_nfya[pl1,],cluster_rows=F,cluster_col=F,main="nfya",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(rat_nfya[pl1,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(rat_input[pl1,],cluster_rows=F,cluster_col=F,main="input",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(rat_input[pl1,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
dev.off()
# marmoset
pdf("../chipseq_heatmap/marmoset.chip.pdf",width=4,height=3)
pheatmap(marmoset_tcfl5[pl,],cluster_rows=F,cluster_col=F,main="tcfl5",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(marmoset_tcfl5[pl,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(marmoset_amyb[pl,],cluster_rows=F,cluster_col=F,main="amyb",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(marmoset_amyb[pl,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(marmoset_nfya[pl,],cluster_rows=F,cluster_col=F,main="nfya",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(marmoset_nfya[pl,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(marmoset_input[pl,],cluster_rows=F,cluster_col=F,main="input",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(marmoset_input[pl,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
dev.off()
# macaque
pdf("../chipseq_heatmap/macaque.chip.pdf",width=4,height=3)
pheatmap(macaque_tcfl5[pl,],cluster_rows=F,cluster_col=F,main="tcfl5",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(macaque_tcfl5[pl,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(macaque_amyb[pl,],cluster_rows=F,cluster_col=F,main="amyb",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(macaque_amyb[pl,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(macaque_nfya[pl,],cluster_rows=F,cluster_col=F,main="nfya",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(macaque_nfya[pl,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(macaque_input[pl,],cluster_rows=F,cluster_col=F,main="input",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(macaque_input[pl,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
dev.off()
# mouse
pdf("../chipseq_heatmap/mouse.chip.pdf",width=4,height=3)
pl1[2]="Ddx39"
pheatmap(mouse_tcfl5[pl1,],cluster_rows=F,cluster_col=F,main="tcfl5",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(mouse_tcfl5[pl1,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(mouse_tcfl5_input[pl1,],cluster_rows=F,cluster_col=F,main="tcfl5_input",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(mouse_tcfl5_input[pl1,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(mouse_amyb[pl1,],cluster_rows=F,cluster_col=F,main="amyb",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(mouse_amyb[pl1,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(mouse_nfya[pl1,],cluster_rows=F,cluster_col=F,main="nfya",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(mouse_nfya[pl1,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(mouse_nfya_input[pl1,],cluster_rows=F,cluster_col=F,main="nfya_input",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(mouse_nfya_input[pl1,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
dev.off()
# human
pdf("../chipseq_heatmap/human.chip.pdf",width=4,height=3)
pheatmap(human_tcfl5_ID1[pl,],cluster_rows=F,cluster_col=F,main="tcfl5_ID1",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(human_tcfl5_ID1[pl,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(human_tcfl5_ID2[pl,],cluster_rows=F,cluster_col=F,main="tcfl5_ID2",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(human_tcfl5_ID2[pl,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(human_tcfl5_input_ID3[pl,],cluster_rows=F,cluster_col=F,main="tcfl5_input_ID3",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(human_tcfl5_input_ID3[pl,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(human_tcfl5_input_ID4[pl,],cluster_rows=F,cluster_col=F,main="tcfl5_input_ID4",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(human_tcfl5_input_ID4[pl,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(human_amyb[pl,],cluster_rows=F,cluster_col=F,main="amyb",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(human_amyb[pl,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(human_amyb_input[pl,],cluster_rows=F,cluster_col=F,main="amyb_input",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(human_amyb_input[pl,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(human_nfya[pl,],cluster_rows=F,cluster_col=F,main="nfya",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(human_nfya[pl,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
pheatmap(human_nfya_input[pl,],cluster_rows=F,cluster_col=F,main="nfya_input",labels_col=rn)
par(mar=c(3,3,1,1),tcl=0.3)
plot(colMeans(human_nfya_input[pl,]),type="l",ylab="",xlab="",xaxt="n")
axis(1,c(1,100,200),label=c("-5kb","TSS","+5kb"))
dev.off()


# May06.2018, TCFL5 and AMYB in proteincoding genes----
setwd("/Users/BigBear/My_Publish/Human_piRNA/motifs/motif_score/")
amyb=read.table("myba.mm10.in_gene+pi.result",header=F,row.names=1)
tcfl5=read.table("tcfl5.mm10.in_gene+pi.result",header=F,row.names=1)
nfy=read.table("nfy.mm10.in_gene+pi.result",header=F,row.names=1)
sig_amyb=read.table("3MO-2-AMYB.FE",header=F,row.names=1)
sig_tcfl5=read.table("TCFL5_3MO_rep1.FE",header=F,row.names=1)
sig_nfya=read.table("NFYA_3MO-3_rep1.FE",header=F,row.names=1)
sig_mat=cbind(sig_amyb,sig_tcfl5,sig_nfya)
peak_distance=read.table("../../conserved_peaks/gene+pi.distance.txt",header=T,row.names=1)
tpd=peak_distance[proteincoding,]
tpd=tpd[which(sig_mat[row.names(tpd),1]>0),]
amyb_list=row.names(tpd[tpd[,1]<500,])
tcfl5_list=row.names(tpd[tpd[,2]<500,])
nfya_list=row.names(tpd[tpd[,3]<500,])
plt=venn.diagram(list(AMYB=amyb_list,TCFL5=tcfl5_list,NFYA=nfya_list),
             filename=NULL,margin=0.1,fill=csl[1:3],alpha=c(0.5,0.5,0.5))
pdf("../../figures/May06.2018.TCFL5vsAMYB/venn.tcfl5_amyb_nfya.pdf",width=4,height=4)
par(mar=c(0,0,0,0))
plot.new()
text(0.5,0.5,label="protein-coding genes with AMYB,\nTCFL5 or NFYA peaks within \nTSS+-500bp")
plot.new()
grid.draw(plt)
dev.off()
plt=venn.diagram(list(AMYB=amyb_list,TCFL5=tcfl5_list),
                 filename=NULL,margin=0.1,fill=csl[1:2],alpha=c(0.5,0.5))
pdf("../../figures/May06.2018.TCFL5vsAMYB/venn.tcfl5_amyb.pdf",width=4,height=4)
par(mar=c(0,0,0,0))
plot.new()
text(0.5,0.5,label="protein-coding genes with AMYB,\nor TCFL5 peaks within TSS+-500bp")
plot.new()
grid.draw(plt)
dev.off()

# heatmap for tcfl5, amyb and nfya
t=intersect(amyb_list,intersect(tcfl5_list,nfya_list));t=t[hclust(dist(sig_mat[t,]))$order];t1=t
t=intersect(amyb_list,setdiff(tcfl5_list,nfya_list));t=t[hclust(dist(sig_mat[t,]))$order];t2=t
t=intersect(amyb_list,setdiff(nfya_list,tcfl5_list));t=t[hclust(dist(sig_mat[t,]))$order];t3=t
t=intersect(tcfl5_list,setdiff(nfya_list,amyb_list));t=t[hclust(dist(sig_mat[t,]))$order];t4=t
t=setdiff(amyb_list,c(tcfl5_list,nfya_list));t=t[hclust(dist(sig_mat[t,]))$order];t5=t
t=setdiff(tcfl5_list,c(amyb_list,nfya_list));t=t[hclust(dist(sig_mat[t,]))$order];t6=t
t=setdiff(nfya_list,c(amyb_list,tcfl5_list));t=t[hclust(dist(sig_mat[t,]))$order];t7=t

annoR=data.frame(Type=factor(rep(c("AMYB_TCFL5_NFYA","AMYB_TCFL5","AMYB_NFYA",
                                   "TCFL5_NFYA","AMYB","TCFL5","NFYA"),
                                 c(length(t1),length(t2),length(t3),length(t4),
                                   length(t5),length(t6),length(t7)))))
ann_colors=list(Type=c(AMYB_TCFL5_NFYA="black",AMYB_TCFL5=csl[4],AMYB_NFYA=csl[5],
                       TCFL5_NFYA=csl[7],AMYB=csl[1],TCFL5=csl[2],NFYA=csl[3]))
row.names(annoR)=c(t1,t2,t3,t4,t5,t6,t7)
sig_mat_nor=normalize.quantiles(as.matrix(sig_mat))
row.names(sig_mat_nor)=row.names(sig_mat);colnames(sig_mat_nor)=colnames(sig_mat)
pheatmap(log10(t(sig_mat_nor[c(t1,t2,t3,t4,t5,t6,t7),])+1),cluster_rows=F,cluster_cols=F,
         show_colnames=F,labels_col=c("AMYB","TCFL5","NFYA"),annotation_col=annoR,
         filename="pc.signal.pdf",cellwidth=0.1,cellheight=32,annotation_colors=ann_colors,
         labels_row=c("amyb","tcfl5","nfya"))
pheatmap(t(motif_score[c(t1,t2,t3,t4,t5,t6,t7),]),cluster_rows=F,cluster_cols=F,
         show_colnames=F,labels_col=c("AMYB","TCFL5","NFYA"),annotation_col=annoR,
         filename="pc.motif_score.pdf",cellwidth=0.1,cellheight=32,
         annotation_colors=ann_colors,color=c("grey",colorRampPalette(brewer.pal(9,"Reds"))(100)),
         breaks=c(0,1:100/15+4.5))

# heatmap for only tcfl5 and amyb
t=intersect(amyb_list,tcfl5_list);t=t[hclust(dist(sig_mat[t,]))$order];t1=t
t=setdiff(amyb_list,tcfl5_list);t=t[hclust(dist(sig_mat[t,]))$order];t2=t
t=setdiff(tcfl5_list,amyb_list);t=t[hclust(dist(sig_mat[t,]))$order];t3=t

annoR=data.frame(Type=factor(rep(c("AMYB_TCFL5","AMYB","TCFL5"),
                                 c(length(t1),length(t2),length(t3)))))
ann_colors=list(Type=c(AMYB_TCFL5=csl[4],AMYB=csl[1],TCFL5=csl[2]))
row.names(annoR)=c(t1,t2,t3)
sig_mat_nor=normalize.quantiles(as.matrix(sig_mat))
row.names(sig_mat_nor)=row.names(sig_mat);colnames(sig_mat_nor)=colnames(sig_mat)
pheatmap(log10(t(sig_mat_nor[c(t1,t2,t3),1:2])+1),cluster_rows=F,cluster_cols=F,
         show_colnames=F,annotation_col=annoR,
         filename="pc.signal.pdf",cellwidth=0.1,cellheight=32,annotation_colors=ann_colors,
         labels_row=c("amyb","tcfl5"))
pheatmap(t(motif_score[c(t1,t2,t3),1:2]),cluster_rows=F,cluster_cols=F,
         show_colnames=F,annotation_col=annoR,
         filename="pc.motif_score.pdf",cellwidth=0.1,cellheight=32,
         annotation_colors=ann_colors,color=c("grey",colorRampPalette(brewer.pal(9,"Reds"))(100)),
         breaks=c(0,1:100/15+4.5))

# boxplot
pdf("pc.AMYB_score.pdf",width=3,height=4)
par(tcl=0.3,bty="n",cex=5/6)
boxplot(motif_score[t1,1],motif_score[t2,1],ylab="motif score",
        xaxt="n",main="5806 AMYB protein-coding\npromoters p-value<2.2e-16",
        border=csl[c(4,1)],lwd=2)
dev.off()
#points(runif(length(t1),0.6,1.4),motif_score[t1,1],pch=19,col=csl[1])
#points(runif(length(t2),1.6,2.4),motif_score[t2,1],pch=19,col=csl[3])
axis(1,1:2,label=c("with TCFL5\npromoters","AMYB specific\npromoters"),lwd=0,cex.axis=0.7)
dev.off()
pdf("pc.TCFL5_score.pdf",width=3,height=4)
par(tcl=0.3,bty="n",cex=5/6)
boxplot(motif_score[t1,2],motif_score[t3,2],ylab="motif score",
        xaxt="n",main="11010 TCFL5 protein-coding\npromoters p-value<1.5e-8",
        border=csl[c(4,2)],lwd=2)
#points(runif(length(t1),0.6,1.4),motif_score[t1,1],pch=19,col=csl[1])
#points(runif(length(t2),1.6,2.4),motif_score[t2,1],pch=19,col=csl[3])
axis(1,1:2,label=c("with AMYB\npromoters","TCFL5 specific\npromoters"),lwd=0,cex.axis=0.7)
dev.off()





write.table(amyb_list,"amyb_proteincoding.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(tcfl5_list,"tcfl5_proteincoding.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(nfya_list,"nfya_proteincoding.txt",row.names=F,col.names=F,sep="\t",quote=F)




# May10.2018, TCFL5 and AMYB in proteincoding genes; for fc5 q0.00001----
setwd("/Users/BigBear/My_Publish/Human_piRNA/motifs/motif_score/")
amyb=read.table("myba.mm10.in_gene+pi.result",header=F,row.names=1)
tcfl5=read.table("tcfl5.mm10.in_gene+pi.result",header=F,row.names=1)
nfy=read.table("nfy.mm10.in_gene+pi.result",header=F,row.names=1)
sig_amyb=read.table("3MO-2-AMYB.FE",header=F,row.names=1)
sig_tcfl5=read.table("TCFL5_3MO_rep1.FE",header=F,row.names=1)
sig_nfya=read.table("NFYA_3MO-3_rep1.FE",header=F,row.names=1)
sig_mat=cbind(sig_amyb,sig_tcfl5,sig_nfya)
peak_distance=read.table("../../conserved_peaks/gene+pi.distance.fc10.q0.00001.txt",header=T,row.names=1)
tpd=peak_distance[proteincoding,]
tpd=tpd[which(sig_mat[row.names(tpd),1]>0),]
amyb_list=row.names(tpd[tpd[,1]<500,])
tcfl5_list=row.names(tpd[tpd[,2]<500,])
nfya_list=row.names(tpd[tpd[,3]<500,])
plt=venn.diagram(list(AMYB=amyb_list,TCFL5=tcfl5_list,NFYA=nfya_list),
                 filename=NULL,margin=0.1,fill=csl[1:3],alpha=c(0.5,0.5,0.5))
pdf("../../figures/May06.2018.TCFL5vsAMYB/venn.tcfl5_amyb_nfya.fc10.pdf",width=4,height=4)
par(mar=c(0,0,0,0))
plot.new()
text(0.5,0.5,label="protein-coding genes with AMYB,\nTCFL5 or NFYA peaks within \nTSS+-500bp")
plot.new()
grid.draw(plt)
dev.off()
plt=venn.diagram(list(AMYB=amyb_list,TCFL5=tcfl5_list),
                 filename=NULL,margin=0.1,fill=csl[1:2],alpha=c(0.5,0.5))
pdf("../../figures/May06.2018.TCFL5vsAMYB/venn.tcfl5_amyb.fc10.pdf",width=4,height=4)
par(mar=c(0,0,0,0))
plot.new()
text(0.5,0.5,label="protein-coding genes with AMYB,\nor TCFL5 peaks within TSS+-500bp")
plot.new()
grid.draw(plt)
dev.off()

# heatmap for tcfl5, amyb and nfya
t=intersect(amyb_list,intersect(tcfl5_list,nfya_list));t=t[hclust(dist(sig_mat[t,]))$order];t1=t
t=intersect(amyb_list,setdiff(tcfl5_list,nfya_list));t=t[hclust(dist(sig_mat[t,]))$order];t2=t
t=intersect(amyb_list,setdiff(nfya_list,tcfl5_list));t=t[hclust(dist(sig_mat[t,]))$order];t3=t
t=intersect(tcfl5_list,setdiff(nfya_list,amyb_list));t=t[hclust(dist(sig_mat[t,]))$order];t4=t
t=setdiff(amyb_list,c(tcfl5_list,nfya_list));t=t[hclust(dist(sig_mat[t,]))$order];t5=t
t=setdiff(tcfl5_list,c(amyb_list,nfya_list));t=t[hclust(dist(sig_mat[t,]))$order];t6=t
t=setdiff(nfya_list,c(amyb_list,tcfl5_list));t=t[hclust(dist(sig_mat[t,]))$order];t7=t

annoR=data.frame(Type=factor(rep(c("AMYB_TCFL5_NFYA","AMYB_TCFL5","AMYB_NFYA",
                                   "TCFL5_NFYA","AMYB","TCFL5","NFYA"),
                                 c(length(t1),length(t2),length(t3),length(t4),
                                   length(t5),length(t6),length(t7)))))
ann_colors=list(Type=c(AMYB_TCFL5_NFYA="black",AMYB_TCFL5=csl[4],AMYB_NFYA=csl[5],
                       TCFL5_NFYA=csl[7],AMYB=csl[1],TCFL5=csl[2],NFYA=csl[3]))
row.names(annoR)=c(t1,t2,t3,t4,t5,t6,t7)
sig_mat_nor=normalize.quantiles(as.matrix(sig_mat))
row.names(sig_mat_nor)=row.names(sig_mat);colnames(sig_mat_nor)=colnames(sig_mat)
pheatmap(log10(t(sig_mat_nor[c(t1,t2,t3,t4,t5,t6,t7),])+1),cluster_rows=F,cluster_cols=F,
         show_colnames=F,labels_col=c("AMYB","TCFL5","NFYA"),annotation_col=annoR,
         filename="pc.signal.fc10.pdf",cellwidth=0.1,cellheight=32,annotation_colors=ann_colors,
         labels_row=c("amyb","tcfl5","nfya"))
pheatmap(t(motif_score[c(t1,t2,t3,t4,t5,t6,t7),]),cluster_rows=F,cluster_cols=F,
         show_colnames=F,labels_col=c("AMYB","TCFL5","NFYA"),annotation_col=annoR,
         filename="pc.motif_score.fc10.pdf",cellwidth=0.1,cellheight=32,
         annotation_colors=ann_colors,color=c("grey",colorRampPalette(brewer.pal(9,"Reds"))(100)),
         breaks=c(0,1:100/15+4.5))

# heatmap for only tcfl5 and amyb
t=intersect(amyb_list,tcfl5_list);t=t[hclust(dist(sig_mat[t,]))$order];t1=t
t=setdiff(amyb_list,tcfl5_list);t=t[hclust(dist(sig_mat[t,]))$order];t2=t
t=setdiff(tcfl5_list,amyb_list);t=t[hclust(dist(sig_mat[t,]))$order];t3=t

annoR=data.frame(Type=factor(rep(c("AMYB_TCFL5","AMYB","TCFL5"),
                                 c(length(t1),length(t2),length(t3)))))
ann_colors=list(Type=c(AMYB_TCFL5=csl[4],AMYB=csl[1],TCFL5=csl[2]))
row.names(annoR)=c(t1,t2,t3)
sig_mat_nor=normalize.quantiles(as.matrix(sig_mat))
row.names(sig_mat_nor)=row.names(sig_mat);colnames(sig_mat_nor)=colnames(sig_mat)
pheatmap(log10(t(sig_mat_nor[c(t1,t2,t3),1:2])+1),cluster_rows=F,cluster_cols=F,
         show_colnames=F,annotation_col=annoR,
         filename="pc.signal.fc10.pdf",cellwidth=0.1,cellheight=32,annotation_colors=ann_colors,
         labels_row=c("amyb","tcfl5"))
pheatmap(t(motif_score[c(t1,t2,t3),1:2]),cluster_rows=F,cluster_cols=F,
         show_colnames=F,annotation_col=annoR,
         filename="pc.motif_score.fc10.pdf",cellwidth=0.1,cellheight=32,
         annotation_colors=ann_colors,color=c("grey",colorRampPalette(brewer.pal(9,"Reds"))(100)),
         breaks=c(0,1:100/15+4.5))

# boxplot
pdf("pc.AMYB_score.fc10.pdf",width=3,height=4)
par(tcl=0.3,bty="n",cex=5/6)
boxplot(motif_score[t1,1],motif_score[t2,1],ylab="motif score",
        xaxt="n",main="5806 AMYB protein-coding\npromoters p-value<2.2e-16",
        border=csl[c(4,1)],lwd=2)
#points(runif(length(t1),0.6,1.4),motif_score[t1,1],pch=19,col=csl[1])
#points(runif(length(t2),1.6,2.4),motif_score[t2,1],pch=19,col=csl[3])
axis(1,1:2,label=c("with TCFL5\npromoters","AMYB specific\npromoters"),lwd=0,cex.axis=0.7)
dev.off()
pdf("pc.TCFL5_score.fc10.pdf",width=3,height=4)
par(tcl=0.3,bty="n",cex=5/6)
boxplot(motif_score[t1,2],motif_score[t3,2],ylab="motif score",
        xaxt="n",main="11010 TCFL5 protein-coding\npromoters p-value<1.5e-8",
        border=csl[c(4,2)],lwd=2)
#points(runif(length(t1),0.6,1.4),motif_score[t1,1],pch=19,col=csl[1])
#points(runif(length(t2),1.6,2.4),motif_score[t2,1],pch=19,col=csl[3])
axis(1,1:2,label=c("with AMYB\npromoters","TCFL5 specific\npromoters"),lwd=0,cex.axis=0.7)
dev.off()
write.table(amyb_list,"amyb_proteincoding.fc10.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(tcfl5_list,"tcfl5_proteincoding.fc10.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(nfya_list,"nfya_proteincoding.fc10.txt",row.names=F,col.names=F,sep="\t",quote=F)

# May10.2018, merged NFYA heatmap----
setwd("/Users/BigBear/My_Publish/Human_piRNA/aggreMatrix/")
merged_nfya=read.table("prepachy_merged_nfya.rpm.uniq.t5.u5000.d5000.b200.mat",header=F,row.names=1)
annoR=data.frame(Type=factor(rep(c("Pachytene","Prepachytene","Hybrid"),c(100,84,30))))
row.names(annoR)=c(Pachy,Prepachy,Hybrid)
ann_colors=list(Type=c(Pachytene="#E41A1C",Hybrid="#984ea3",Prepachytene="#4daf4a"),
                position=c(beforeTSS="grey",afterTSS="navy"),Location=c(Intergenic="black",Genic="white"))
pheatmap(merged_nfya[c(Pachy,Prepachy,Hybrid),],cluster_rows=F,cluster_col=F,
         main="nfya_in_prepachyteneCell",labels_col=rn,fontsize_row=8,
         annotation_row=annoR,annotation_colors=ann_colors,
         filename="../chipseq_heatmap/mouse_merged_nfya_prepachy.pdf",
         cellwidth=1.5,height=18,
         col=colorRampPalette(c(brewer.pal(9,"Blues")[9:1],"white"))(100))
pheatmap(merged_nfya[c(Pachy),],cluster_rows=F,cluster_col=F,
         main="nfya_in_prepachyteneCell",labels_col=rn,fontsize_row=8,
         annotation_row=annoR,annotation_colors=ann_colors,
         filename="../chipseq_heatmap/mouse_merged_nfya_prepachy1.pdf",
         cellwidth=1.5,height=18,
         col=colorRampPalette(c(brewer.pal(9,"Blues")[9:1],"white"))(100))





# May11.2018, correlation_chipseq_FE_3month----
setwd("/Users/BigBear/My_Publish/Human_piRNA/motifs/motif_score/")
amyb=read.table("myba.mm10.in_gene+pi.result",header=F,row.names=1)
tcfl5=read.table("tcfl5.mm10.in_gene+pi.result",header=F,row.names=1)
nfy=read.table("nfy.mm10.in_gene+pi.result",header=F,row.names=1)
sig_amyb=read.table("3MO-2-AMYB.FE",header=F,row.names=1)
sig_tcfl5=read.table("TCFL5_3MO_rep1.FE",header=F,row.names=1)
sig_nfya=read.table("NFYA_3MO-3_rep1.FE",header=F,row.names=1)
x=log2(sig_amyb[Pachy,1]+1);y=log2(sig_tcfl5[Pachy,1]+1)

pdf("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/Bear_Feb21.2018.peakConservation/correlation_chipseq_FE_3month.pdf",
    width=4,height=4,useDingbats=F)
par(mar=c(4,4,2,2))
plot(x,y,col=csl[1],xlim=c(min(x),max(x)),ylim=c(min(y),max(y)),lwd=2,
     xlab="AMYB 3MO",ylab="TCFL5 3MO")
text(max(x),max(y),label=paste("cor = ",round(cor(x,y,method="spearman"),2),sep=""),
     col=csl[1],pos=2)
dev.off()






# May23.2018, Figure1----
# amyb heatmaps
setwd("/Users/BigBear/My_Publish/Human_piRNA/aggreMatrix/")
human_amyb=read.table("human_AMYB_piRNA.t53.u3000.d3000.b200.mat",head=F,row.names=1)
fun_plot=function(i,fn,mn){
  pdf(fn,width=4,height=8)
  par(tcl=0.3)
  m=as.matrix(i)
  row.names(m)=rep("",dim(m)[1])
  colnames(m)=c("-3kb",rep("",49),"TSS",rep("",98),"TES",rep("",49),"3kb")
  par(mar=c(2,2,2,2),bty="n",mfrow=c(3,1))
  plot(apply(i,2,mean,trim=0.05),type="l",col="#313695",lwd=2,ylab="",xlab="",xaxt="n",ylim=c(0,30))
  mtext(mn,side=3,font=2,cex=1.2)
  axis(1,c(0,50,150,200),label=c("-3kb","TSS","TES","3kb"))
  for(i in 1:dim(m)[2]){m[m[,i]>40,i]=40}
  pheatmap(m,cluster_rows=F,cluster_cols=F,
           col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
           breaks=0:100*40/100)
  dev.off()
}
hg_pachy=as.vector(read.table("../piRNA_list/pachytene.order.list",header=F,row.names=NULL)[,1])
hg_prepachy=as.vector(read.table("../piRNA_list/prepachytene.order.list",header=F,row.names=NULL)[,1])
hg_hybrid=as.vector(read.table("../piRNA_list/hybrid.order.list",header=F,row.names=NULL)[,1])
human_amyb_copy1=human_amyb[hg_pachy,]
human_amyb_copy2=human_amyb[hg_prepachy,]
o1=order(rowMeans(human_amyb_copy1),decreasing=T)
o2=order(rowMeans(human_amyb_copy2),decreasing=T)
human_amyb_copy1=human_amyb_copy1[o1,]
human_amyb_copy2=human_amyb_copy2[o2,]
fun_plot(human_amyb_copy1,"../figures/May23.2018.FigureModification/Figure1a.amyb.pdf",
         "human AMYB; pachytene")
fun_plot(human_amyb_copy2,"../figures/May23.2018.FigureModification/Figure1b.amyb.pdf",
         "human AMYB; prepachytene")
# srnaseq heatmap
setwd("/Users/BigBear/My_Publish/Human_piRNA/ppm/")
cn=c(1,2,5,7,12,10,6,8,11,13)
mat=as.matrix(read.table("allSample.25_31nt.allmap.rpm",header=T,row.names=1))
mat=mat[,cn]
mat_ppm=mat
annoC=data.frame(Sample=factor(rep(c("Juvenile","Adult_30nt","Adult_26nt"),c(2,4,4))))
ann_colors=list(Type=c(Pachytene="#E41A1C",Hybrid="#984ea3",
                       Prepachytene="#4daf4a"),
                Sample=c(Juvenile="#66c2a5",Adult_30nt="#fc8d62",
                         Adult_26nt="#8da0cb"))
row.names(annoC)=colnames(mat)
pdf("../figures/May23.2018.FigureModification/figure1a.srna.1.pdf",width=4,height=8)
pheatmap(mat[hg_pachy[o1],],cluster_rows=F,cluster_col=F,show_rownames=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         annotation_col=annoC,annotation_colors=ann_colors,
         scale="row",
         border="white",lwd=0.3)
pheatmap(log10(mat[hg_pachy[o1],]+1),cluster_rows=F,cluster_col=F,show_rownames=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         annotation_col=annoC,annotation_colors=ann_colors,
         #scale="row",
         border="white",lwd=0.3)
dev.off()
pdf("../figures/May23.2018.FigureModification/figure1b.srna.pdf",width=4,height=8)
pheatmap(mat[hg_prepachy[o2],],cluster_rows=F,cluster_col=F,show_rownames=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         annotation_col=annoC,annotation_colors=ann_colors,scale="row",
         border="white",lwd=0.3)
dev.off()
# rnaseq heatmap
tpm_deniz=read.table("../tpm/genes_expression_tpm.tsv",header=T,row.names=1)
piN_STRG_map=read.table("../piRNA_list/piName_STRG.map",header=F,row.names=1)
mat=tpm_deniz[as.vector(piN_STRG_map[hg_pachy[o1],1]),
              c("D1643_5YO_Rep1","D1643_5YO_Rep2","D7358NT",
                "D8150NT","S7_NT","S12_NT","D7377NT","AmbionNT",
                "S6_NT","S8_NT")]
row.names(annoC)=colnames(mat)
pdf("../figures/May23.2018.FigureModification/figure1a.rna.pdf",width=4,height=8)
pheatmap(mat,cluster_rows=F,cluster_col=F,show_rownames=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         annotation_colors=ann_colors,
         annotation_col=annoC,scale="row",
         border="white",lwd=0.3)
dev.off()
mat=tpm_deniz[as.vector(piN_STRG_map[hg_prepachy[o2],1]),
              c("D1643_5YO_Rep1","D1643_5YO_Rep2","D7358NT",
                "D8150NT","S7_NT","S12_NT","D7377NT","AmbionNT",
                "S6_NT","S8_NT")]
pdf("../figures/May23.2018.FigureModification/figure1b.rna.pdf",width=4,height=8)
pheatmap(mat,cluster_rows=F,cluster_col=F,show_rownames=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         annotation_colors=ann_colors,annotation_col=annoC,scale="row",
         border="white",lwd=0.3)
dev.off()
# srna heatmap for all piRNA genes: Figure1C
oo1=hclust(dist(log10(mat_ppm[c(hg_prepachy),]+10)))$order
oo2=hclust(dist(log10(mat_ppm[c(hg_hybrid),]+10)))$order
oo3=hclust(dist(log10(mat_ppm[c(hg_pachy),]+10)))$order
row.names(annoC)=colnames(mat_ppm)
annoR=data.frame(Type=factor(rep(c("Prepachytene","Hybrid","Pachytene"),c(80,17,153))))
row.names(annoR)=c(hg_prepachy[oo1],hg_hybrid[oo2],hg_pachy[oo3])
pheatmap(log10(t(mat_ppm[c(hg_prepachy[oo1],hg_hybrid[oo2],hg_pachy[oo3]),c(1:2,7:10,3:6)])+10),
         labels_col=rep("",250),labels_row=rep("",10),cluster_rows=F,cluster_cols=F,
         annotation_colors=ann_colors,annotation_row=annoC,annotation_col=annoR,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         filename="../figures/May23.2018.FigureModification/figure1c.pdf",
         cellwidth=4,cellheight=16,border="white",lwd=0.3)


pheatmap(t(mat_ppm[c(hg_prepachy[oo1],hg_hybrid[oo2],hg_pachy[oo3]),c(1:2,7:10,3:6)]),
         labels_col=rep("",250),labels_row=rep("",10),cluster_rows=F,cluster_cols=F,
         annotation_colors=ann_colors,annotation_row=annoC,annotation_col=annoR,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         filename="../figures/May23.2018.FigureModification/figure1c.1.pdf",
         cellwidth=4,cellheight=16,border="white",lwd=0.3,scale="column")


plot.new()
plot(log10(mat_ppm[hg_pachy,]+1))
for(i in 1:dim(mat_ppm[hg_pachy,])[1]){if(length(which(mat_ppm[hg_pachy[i],]>10))==0){print(hg_pachy[i])}}

l=c()
for(i in 1:dim(mat_ppm[hg_pachy,])[1]){
  if(length(which(mat_ppm[hg_pachy[i],]>1))==0){l=c(l,hg_pachy[i])}}
pheatmap(log10(t(mat_ppm[l,c(1:2,7:10,3:6)])+1),
         labels_col=rep("",250),labels_row=rep("",10),cluster_rows=F,cluster_cols=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         filename="../figures/May23.2018.FigureModification/figure1c.pdf",
         cellwidth=4,cellheight=16,border="white",lwd=0.3)



# May25.2018, Annotate piRNA genes, version1.0----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
mat_ppm=read.table("human_ox.TPM0.1.rpm",header=T,row.names=1)
mat_ppm=mat_ppm[,c(1,2,4,5,7,12,3,10,6,8,9,11,13,14)]
for(ct in c(100,200,500,1000)){
print(paste(ct, length(which(rowMeans(mat_ppm[,3:4])>=ct | rowMeans(mat_ppm[,5:12])>=ct)), sep=" "))
}
mat_ppm_pi1=mat_ppm[which(rowMeans(mat_ppm[,3:4])>=100 | rowMeans(mat_ppm[,5:12])>=100),]
ll=c()
for(i in 1:dim(mat_ppm_pi1)[1]){
  s=0
  if(rowMeans(mat_ppm_pi1[i,3:4])<100){
    for(j in 5:12){
      if(mat_ppm_pi1[i,j]>=100){
        s=s+1
      }
    }
    if(s<4){
      ll=c(ll,i)
    }
  }
}
pheatmap(t(mat_ppm_pi1[ll,3:12]),
         cluster_rows=F,cluster_cols=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),filename="../piRNA_cluster_annotation/leaved_piclusters.pdf",
         cellwidth=16,cellheight=16,border="white",lwd=0.3)
mat_ppm_pi2=mat_ppm_pi1[setdiff(1:1141,ll),]
pheatmap(log10(t(mat_ppm_pi2[,3:12])+1),
         labels_col=rep("",250),labels_row=rep("",10),cluster_rows=F,cluster_cols=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         cellwidth=4,cellheight=16,border="white",lwd=0.3,
         filename="../piRNA_cluster_annotation/allpiRNAGenes.heatmap.pdf")
human_pachy=row.names(mat_ppm_pi2)[which(rowMeans(mat_ppm_pi2[,5:12])>4*rowMeans(mat_ppm_pi2[,3:4]))]
pheatmap(log10(t(mat_ppm_pi2[human_pachy,3:12])+1),
         labels_col=rep("",250),labels_row=rep("",10),cluster_rows=F,cluster_cols=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         cellwidth=4,cellheight=16,border="white",lwd=0.3,
         filename="../piRNA_cluster_annotation/Pachytene.heatmap.pdf")
human_hybrid=row.names(mat_ppm_pi2)[which(rowMeans(mat_ppm_pi2[,5:12])<=4*rowMeans(mat_ppm_pi2[,3:4]) &
                                            rowMeans(mat_ppm_pi2[,5:12])>2*rowMeans(mat_ppm_pi2[,3:4]))]
pheatmap(log10(t(mat_ppm_pi2[human_hybrid,3:12])+1),
         labels_col=rep("",250),labels_row=rep("",10),cluster_rows=F,cluster_cols=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         cellwidth=4,cellheight=16,border="white",lwd=0.3,
         filename="../piRNA_cluster_annotation/Hybrid.heatmap.pdf")
human_prepachy=row.names(mat_ppm_pi2)[which(rowMeans(mat_ppm_pi2[,5:12])<=2*rowMeans(mat_ppm_pi2[,3:4]))]
pheatmap(log10(t(mat_ppm_pi2[human_prepachy[1:100],3:12])+1),
         labels_col=rep("",250),labels_row=rep("",10),cluster_rows=F,cluster_cols=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         cellwidth=4,cellheight=16,border="white",lwd=0.3,
         filename="../piRNA_cluster_annotation/Prepachytene.heatmap.top100.pdf")

annoC=data.frame(Type=factor(rep(c("Prepachytene","Hybrid","Pachytene"),c(100,19,116))))
annoR=data.frame(Sample=factor(rep(c("Juvenile","Adult_26nt","Adult_30nt"),c(2,4,4))))
ann_colors=list(Type=c(Prepachytene="#4daf4a",Hybrid="#984ea3",Pachytene="#E41A1C"),
                Sample=c(Juvenile="#66c2a5",Adult_26nt="#fc8d62",
                         Adult_30nt="#8da0cb"))
mt=t(mat_ppm_pi2[c(human_prepachy[1:100],human_hybrid,human_pachy),3:12])
row.names(annoR)=row.names(mt)
row.names(annoC)=colnames(mt)
pheatmap(log10(mt+1),
         labels_col=rep("",235),cluster_rows=F,cluster_cols=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         cellwidth=4,cellheight=16,border="white",lwd=0.3,
         filename="../piRNA_cluster_annotation/Pachy_Prepachy_Hybrid.heatmap.pdf",
         annotation_row=annoR,annotation_col=annoC,annotation_colors=ann_colors)
mt=t(mat_ppm_pi2[c(human_prepachy[1:100],human_hybrid,human_pachy),])
mt=cbind(t(mt),annoC)
write.table(mt,"../piRNA_cluster_annotation/piRNA.list.txt",quote=F,sep="\t")




# May31.2018, Annotate piRNA genes, version2.0----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/piRNA_cluster_annotation/Version_2.0")
mat_ppm=read.table("stringtie.June1.rpm",header=T,row.names=1)
#mat_ppm=transform(mat_ppm,mean=rowMeans(mat_ppm),max=apply(mat_ppm,1,max))
#write.table(mat_ppm,"tmp",quote=F,sep="\t")
mat_ppm=mat_ppm[,c(1,3,4,2,6,9,11,5,7,8,10,12,13)]
for(ct in c(100,200,500,1000)){
  print(paste(ct, length(which(rowMeans(mat_ppm[,2:3])>=ct | rowMeans(mat_ppm[,4:11])>=ct)), sep=" "))
}
mat_ppm_pi1=mat_ppm[which(rowMeans(mat_ppm[,2:3])>=100 | rowMeans(mat_ppm[,4:11])>=100),]
ll=c()
for(i in 1:dim(mat_ppm_pi1)[1]){
  s=0
  if(rowMeans(mat_ppm_pi1[i,2:3])<100){
    for(j in 4:11){
      if(mat_ppm_pi1[i,j]>=100){
        s=s+1
      }
    }
    if(s<4){
      ll=c(ll,i)
    }
  }
}
pdf("leaved_piclusters.pdf",height=4,width=6)
pheatmap(t(mat_ppm_pi1[ll,2:11]),
         cluster_rows=F,cluster_cols=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         cellwidth=16,cellheight=16,border="white",lwd=0.3)
pheatmap(t(mat_ppm_pi1[ll,2:11]),
         cluster_rows=F,cluster_cols=F,
         col=c("grey",csl[1]),
         breaks=c(0,100,1500),
         cellwidth=16,cellheight=16,border="white",lwd=0.3)
dev.off()
mat_ppm_pi2=mat_ppm_pi1[setdiff(1:dim(mat_ppm_pi1)[1],ll),]

pdf("ss.pdf",width=45,height=4.5)
out=pheatmap(log10(t(mat_ppm_pi2[,2:11])+1),
         clustering_distance_rows="euclidean",clustering_method="complete",
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         cellwidth=4,cellheight=16,border="white",lwd=0.3)
par(mar=c(1,4,2,10))
barplot(log2(apply(mat_ppm_pi2[out$tree_col[["order"]],4:11],1,mean)/
               apply(mat_ppm_pi2[out$tree_col[["order"]],2:3],1,mean)),
        border="white",col="black",lwd=0.02,xlim=c(10,850))
abline(h=c(1,2),lwd=1,lty=2,col=csl[c(1,2)])
dev.off()

par(op)
t1=(rowMeans(df_srna[human_pachy,c(1,2,3,4,5,13,19)]))/(rowMeans(df_srna[human_pachy,c(7,18)]))
t1[which(t1>8)]=8
t2=(rowMeans(df_srna[human_prepachy,c(1,2,3,4,5,13,19)]))/(rowMeans(df_srna[human_prepachy,c(7,18)]))
t2[which(t2>8)]=8
t3=(rowMeans(df_srna[human_hybrid,c(1,2,3,4,5,13,19)]))/(rowMeans(df_srna[human_hybrid,c(7,18)]))
t3[which(t3>8)]=8
hist(t1,breaks=100,xlim=c(0,8),col=csl[1])
hist(t2,breaks=50,xlim=c(0,8),col=csl[2],add=T)
hist(t3,breaks=200,xlim=c(0,8),col=csl[3],add=T)

pheatmap(log10(t(mat_ppm_pi2[,2:11])+1),
         labels_col=rep("",250),labels_row=rep("",10),cluster_rows=F,cluster_cols=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         cellwidth=4,cellheight=16,border="white",lwd=0.3,
         filename="allpiRNAGenes.heatmap.pdf")
human_pachy=row.names(mat_ppm_pi2)[which(rowMeans(mat_ppm_pi2[,4:11])>4*rowMeans(mat_ppm_pi2[,2:3]))]
human_pachy=human_pachy[hclust(dist(mat_ppm_pi2[human_pachy,2:11]))$order]
pheatmap(log10(t(mat_ppm_pi2[human_pachy,2:11])+1),
         labels_col=rep("",250),labels_row=rep("",10),cluster_rows=F,cluster_cols=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         cellwidth=4,cellheight=16,border="white",lwd=0.3,
         filename="Pachytene.heatmap.pdf")
human_hybrid=row.names(mat_ppm_pi2)[which(rowMeans(mat_ppm_pi2[,4:11])<=4*rowMeans(mat_ppm_pi2[,2:3]) &
                                            rowMeans(mat_ppm_pi2[,4:11])>2*rowMeans(mat_ppm_pi2[,2:3]))]
human_hybrid=human_hybrid[hclust(dist(mat_ppm_pi2[human_hybrid,2:11]))$order]
pheatmap(log10(t(mat_ppm_pi2[human_hybrid,2:11])+1),
         labels_col=rep("",250),labels_row=rep("",10),cluster_rows=F,cluster_cols=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         cellwidth=4,cellheight=16,border="white",lwd=0.3,
         filename="Hybrid.heatmap.pdf")
human_prepachy_all=row.names(mat_ppm_pi2)[which(rowMeans(mat_ppm_pi2[,4:11])<=2*rowMeans(mat_ppm_pi2[,2:3]))]
human_prepachy_all=human_prepachy_all[order(mat_ppm_pi2[human_prepachy_all,13],decreasing=T)]
human_prepachy=human_prepachy_all[1:100]
human_prepachy=human_prepachy[hclust(dist(mat_ppm_pi2[human_prepachy[1:100],2:11]))$order]
write.table(human_prepachy,"prepachy_100.txt",row.names=F,col.names=F,sep="\t",quote=F)
pheatmap(log10(t(mat_ppm_pi2[human_prepachy,2:11])+1),
         labels_col=rep("",250),labels_row=rep("",10),cluster_rows=F,cluster_cols=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         cellwidth=4,cellheight=16,border="white",lwd=0.3,
         filename="Prepachytene.heatmap.top100.pdf")

annoC=data.frame(Type=factor(rep(c("Prepachytene","Hybrid","Pachytene"),c(100,9,96))))
annoR=data.frame(Sample=factor(rep(c("Juvenile","Adult_26nt","Adult_30nt"),c(2,4,4))))
ann_colors=list(Type=c(Prepachytene="#4daf4a",Hybrid="#984ea3",Pachytene="#E41A1C"),
                Sample=c(Juvenile="#66c2a5",Adult_26nt="#fc8d62",
                         Adult_30nt="#8da0cb"))
mt=t(mat_ppm_pi2[c(human_prepachy,human_hybrid,human_pachy),2:11])
row.names(annoR)=row.names(mt)
row.names(annoC)=colnames(mt)
pheatmap(log10(mt+1),
         labels_col=rep("",235),cluster_rows=F,cluster_cols=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         cellwidth=4,cellheight=16,border="white",lwd=0.3,
         filename="Pachy_Prepachy_Hybrid.heatmap.pdf",
         annotation_row=annoR,annotation_col=annoC,annotation_colors=ann_colors)
mt=t(mat_ppm_pi2[c(human_prepachy,human_hybrid,human_pachy),])
mt=cbind(t(mt),annoC)
write.table(mt,"piRNA.list.txt",quote=F,sep="\t")

t=mat_ppm_pi2[human_prepachy_all,13]
t[t>1000]=1000
plot(hist(t,breaks=(10:100)*10))
abline(v=460)

# June12.2018, heatmaps for curated piG----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
mat_ppm=read.table("merged.rpm",header=T,row.names=1,check.names=F)
#mat_ppm=transform(mat_ppm,mean=rowMeans(mat_ppm),max=apply(mat_ppm,1,max))
#write.table(mat_ppm,"tmp",quote=F,sep="\t")
# mat_ppm=mat_ppm[,c(3,4,13,2,14,1,6,10,12,5,7,8,9,11)]
# pheatmap(log10(t(mat_ppm)+1),
#          labels_col=rep("",250),labels_row=rep("",10),cluster_rows=F,cluster_cols=F,
#          col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
#          cellwidth=4,cellheight=16,border="white",lwd=0.3,
#          filename="allpiRNAGenes.heatmap.pdf")
# #human_pachy=row.names(mat_ppm)[which(rowMeans(mat_ppm[,3:10])>4*rowMeans(mat_ppm[,1:2]))]
# human_pachy=human_pachy[hclust(dist(mat_ppm[human_pachy,1:14]))$order]
# pheatmap(log10(t(mat_ppm[human_pachy,1:14])+1),
#          labels_col=rep("",250),labels_row=rep("",10),cluster_rows=F,cluster_cols=F,
#          col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
#          cellwidth=4,cellheight=16,border="white",lwd=0.3,
#          filename="Pachytene.heatmap.pdf")
# #human_hybrid=row.names(mat_ppm)[which(rowMeans(mat_ppm[,3:10])<=4*rowMeans(mat_ppm[,1:2]) &
# #                                            rowMeans(mat_ppm[,3:10])>2*rowMeans(mat_ppm[,1:2]))]
# human_hybrid=human_hybrid[hclust(dist(mat_ppm[human_hybrid,1:14]))$order]
# pheatmap(log10(t(mat_ppm[human_hybrid,1:14])+1),
#          labels_col=rep("",250),labels_row=rep("",10),cluster_rows=F,cluster_cols=F,
#          col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
#          cellwidth=4,cellheight=16,border="white",lwd=0.3,
#          filename="Hybrid.heatmap.pdf")
#human_prepachy=row.names(mat_ppm)[which(rowMeans(mat_ppm[,3:10])<=2*rowMeans(mat_ppm[,1:2]))]
#human_prepachy_all=human_prepachy_all[order(mat_ppm[human_prepachy_all,2],decreasing=T)]
#human_prepachy=human_prepachy_all[1:100]
# human_prepachy=human_prepachy[hclust(dist(mat_ppm[human_prepachy,1:14]))$order]
# pheatmap(log10(t(mat_ppm[human_prepachy,1:14])+1),
#          labels_col=rep("",250),labels_row=rep("",10),cluster_rows=F,cluster_cols=F,
#          col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
#          cellwidth=4,cellheight=16,border="white",lwd=0.3,
#          filename="Prepachytene.heatmap.pdf")

annoC=data.frame(Type=factor(rep(c("Pachytene","Hybrid","Prepachytene"),c(89,10,83))))
annoR=data.frame(Sample=factor(rep(c("Juvenile","Adult_26nt","Adult_30nt"),c(3,5,5))))
ann_colors=list(Type=c(Prepachytene=csl[2],Hybrid="#984ea3",Pachytene="#E41A1C"),
                Sample=c(Juvenile="#66c2a5",Adult_26nt="#fc8d62",
                         Adult_30nt="#8da0cb"))
mt=t(cbind(apply(mat_ppm[c(human_pachy,human_hybrid,human_prepachy),1:2],1,mean),
     mat_ppm[c(human_pachy,human_hybrid,human_prepachy),3:14]))
row.names(mt)[1]="Juv"
row.names(annoR)=row.names(mt)
row.names(annoC)=colnames(mt)
pheatmap(log10(mt+1),
         labels_col=rep("",235),cluster_rows=F,cluster_cols=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         cellwidth=4,cellheight=16,border="white",lwd=0.3,
         filename="Pachy_Prepachy_Hybrid.heatmap.pdf",
         annotation_row=annoR,annotation_col=annoC,annotation_colors=ann_colors)
mt=t(mat_ppm[c(human_pachy,human_hybrid,human_prepachy),])
mt=cbind(t(mt),annoC)
write.table(mt,"piRNA.list.txt",quote=F,sep="\t")

#rnaseq
mat_rpkm=read.table("human.piG.rpkm",header=T,row.names=1,check.names=F)
mat_rpkm=mat_rpkm[,c(1,2,6,8,12,14,5,7,11,13)]
annoC=data.frame(Type=factor(rep(c("Pachytene","Hybrid","Prepachytene"),c(90,10,84))))
annoR=data.frame(Sample=factor(rep(c("Juvenile","Adult_26nt","Adult_30nt"),c(2,4,4))))
ann_colors=list(Type=c(Prepachytene="#4daf4a",Hybrid="#984ea3",Pachytene="#E41A1C"),
                Sample=c(Juvenile="#66c2a5",Adult_26nt="#fc8d62",
                         Adult_30nt="#8da0cb"))
mt=t(mat_rpkm[c(human_pachy,human_hybrid,human_prepachy),1:10])
row.names(annoR)=row.names(mt)
row.names(annoC)=colnames(mt)
pheatmap(log10(mt+1),
         labels_col=rep("",235),cluster_rows=F,cluster_cols=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         cellwidth=4,cellheight=16,border="white",lwd=0.3,
         filename="rnaseq.Pachy_Prepachy_Hybrid.heatmap.pdf",
         annotation_row=annoR,annotation_col=annoC,annotation_colors=ann_colors)

# June15.2018: heatmap for all chipseq----
pl=c("DDX4","DDX39A","FKBP6","GPAT2","HENMT1","MAEL","MOV10L1","MYBL1","NFYA",
     "PIWIL1","PIWIL2","PIWIL4","PLD6","TDRD1","TDRD3","TDRD5","TDRD6","TDRD9",
     "TDRKH")
pl1=c("Ddx4","Ddx39a","Fkbp6","Gpat2","Henmt1","Mael","Mov10l1","Mybl1","Nfya",
      "Piwil1","Piwil2","Piwil4","Pld6","Tdrd1","Tdrd3","Tdrd5","Tdrd6","Tdrd9",
      "Tdrkh")
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/aggreMatrix/matrix_FE/")
list_hs=list()
l_hs=c("5958NT-AMYB-Id2_S2","6249NT-AMYB-Id4_S5","6615NT-AMYB-Id1_S10","7358-AMYB-Id7_S1",
       "8246NT-AMYB-Id3_S8","P75-Hs-AMYB-ID1_S5","s904-Hs-AMYB-ID4_S4",
       "s904-Hs-NFYA-ID5_S3","P75-Hs-NFYA-ID3_S2","P75-Hs-TCFL5-ID2_S4","hs-kidney-AMYB-ID8_S3","hs-liver-AMYB-ID3_S4",
       "hs-spleen-AMYB-ID5_S11")
for(i in l_hs){
  list_hs[[i]]=read.table(paste(i,"FE.t5.u5000.d5000.b200.mat",sep="."),
                          header=F,row.names=1)
}

list_mac=list()
l_mac=c("Macaque-Input-ID8_S7","Macaque-AMYB-ID5_S8")
for(i in l_mac){
  list_mac[[i]]=read.table(paste("../matrix_rpm/",i,".rpm.uniq.t5.u5000.d5000.b200.mat",sep=""),
                          header=F,row.names=1)
}
list_rn=list()
l_rn=c("Rat-AMYB-ID5_S7","Rat-AMYB-ID7_S8","Rat-K27Ac-ID9_S3_K27ac",
        "Rat-NFYA-ID7_S3","Rat-NFYA-ID8_S7","Rat-TCFL5-ID6_S1")
for(i in l_rn){
  list_rn[[i]]=read.table(paste(i,"FE.t5.u5000.d5000.b200.mat",sep="."),
                           header=F,row.names=1)
}
list_marmo=list()
l_marmo=c("Marmoset-AMYB-ID1_S8","Marmoset-AMYB-ID4_S1","Marmoset-NFYA-ID3_S2",
       "Marmoset-NFYA-ID6_S3","Marmoset-TCFL5-ID2_S4","Marmoset-TCFL5-ID5_S6")
for(i in l_marmo){
  list_marmo[[i]]=read.table(paste(i,"FE.t5.u5000.d5000.b200.mat",sep="."),
                          header=F,row.names=1)
}
# human
pdf("../figures_FE/human_all.pdf",width=3.5,height=7)
for(i in names(list_hs)){
  annoC=data.frame(Type=factor(rep(c("Pachytene","Hybrid","Prepachytene"),c(90,10,84))))
  row.names(annoC)=c(human_pachy,human_hybrid,human_prepachy)
  ann_colors=list(Type=c(Prepachytene="#4daf4a",Hybrid="#984ea3",Pachytene="#E41A1C"),
                  Sample=c(Juvenile="#66c2a5",Adult_26nt="#fc8d62",
                           Adult_30nt="#8da0cb"))
  pheatmap(log10(list_hs[[i]][c(human_pachy,human_hybrid,human_prepachy),]+1),
           cluster_rows=F,cluster_cols=F,labels_row="",main=i,
           #cellwidth=0.7,cellheight=2,
           labels_col=c("-5kb",rep("",99),"TSS",rep("",98),"+5kb"),
           annotation_row=annoC,annotation_colors=ann_colors,
           col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100))
  #filename=paste("heatmap.amyb.piRNA",i,"pdf",sep=".")
}
dev.off()
# human rpm
for(i in l_hs){
  list_hs[[i]]=read.table(paste(i,"rpm.t5.u5000.d5000.b200.mat",sep="."),
                          header=F,row.names=1)
}
pdf("../figures_rpm/human_all.pdf",width=3.5,height=7)
for(i in names(list_hs)){
  annoC=data.frame(Type=factor(rep(c("Pachytene","Hybrid","Prepachytene"),c(90,10,84))))
  row.names(annoC)=c(human_pachy,human_hybrid,human_prepachy)
  ann_colors=list(Type=c(Prepachytene="#4daf4a",Hybrid="#984ea3",Pachytene="#E41A1C"),
                  Sample=c(Juvenile="#66c2a5",Adult_26nt="#fc8d62",
                           Adult_30nt="#8da0cb"))
  pheatmap(log10(list_hs[[i]][c(human_pachy,human_hybrid,human_prepachy),]+1),
           cluster_rows=F,cluster_cols=F,labels_row="",main=i,
           #cellwidth=0.7,cellheight=2,
           labels_col=c("-5kb",rep("",99),"TSS",rep("",98),"+5kb"),
           annotation_row=annoC,annotation_colors=ann_colors,
           breaks=c((0:99/100)*0.4,0.5),
           col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100))
  #filename=paste("heatmap.amyb.piRNA",i,"pdf",sep=".")
}
dev.off()
# macaque
pdf("../figures_rpm/macaque_all.pdf",width=3.5,height=7,useDingbats=F)
tn=row.names(list_mac[[1]])
pi_pc=tn[grep("_PC_",tn)]
pi_ig=tn[grep("_IG_",tn)]
tn=c(pi_ig,pi_pc,pl)
for(i in names(list_mac)){
  annoC=data.frame(Type=factor(rep(c("intergenic_piG","genic_piG","pathway"),
                                   c(length(pi_ig),length(pi_pc),length(pl1)))))
  row.names(annoC)=tn
  ann_colors=list(Type=c(genic_piG="#4daf4a",pathway="#984ea3",intergenic_piG="#E41A1C"),
                  Sample=c(Juvenile="#66c2a5",Adult_26nt="#fc8d62",
                           Adult_30nt="#8da0cb"))
  pheatmap(log10(list_mac[[i]][tn,]/43/43+0.1),
           cluster_rows=F,cluster_cols=F,labels_row="",main=i,
           #cellwidth=0.7,cellheight=2,
           labels_col=c("-5kb",rep("",99),"TSS",rep("",98),"+5kb"),
           annotation_row=annoC,annotation_colors=ann_colors,
           col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
           breaks=(0:100)/100*2.5-1)
  #filename=paste("heatmap.amyb.piRNA",i,"pdf",sep=".")
}
dev.off()
# rat
pdf("../figures_FE/rat_all.pdf",width=3.5,height=7)
tn=row.names(list_rn[[1]])
pi_pc=tn[grep("_PC_",tn)]
pi_ig=tn[grep("_IG_",tn)]
tn=c(pi_ig,pi_pc,pl1)
for(i in names(list_rn)){
  annoC=data.frame(Type=factor(rep(c("intergenic_piG","genic_piG","pathway"),
                                   c(length(pi_ig),length(pi_pc),length(pl1)))))
  row.names(annoC)=tn
  ann_colors=list(Type=c(genic_piG="#4daf4a",pathway="#984ea3",intergenic_piG="#E41A1C"),
                  Sample=c(Juvenile="#66c2a5",Adult_26nt="#fc8d62",
                           Adult_30nt="#8da0cb"))
  pheatmap(log10(list_rn[[i]][tn,]+1),
           cluster_rows=F,cluster_cols=F,labels_row="",main=i,
           #cellwidth=0.7,cellheight=2,
           labels_col=c("-5kb",rep("",99),"TSS",rep("",98),"+5kb"),
           annotation_row=annoC,annotation_colors=ann_colors,
           col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100))
  #filename=paste("heatmap.amyb.piRNA",i,"pdf",sep=".")
}
dev.off()
# marmoset
pdf("../figures_FE/marmoset_all.pdf",width=3.5,height=7)
tn=row.names(list_marmo[[1]])
pi_pc=tn[grep("_PC_",tn)]
pi_ig=tn[grep("_IG_",tn)]
tn=c(pi_ig,pi_pc,pl1)
for(i in names(list_marmo)){
  annoC=data.frame(Type=factor(rep(c("intergenic_piG","genic_piG","pathway"),
                                   c(length(pi_ig),length(pi_pc),length(pl1)))))
  row.names(annoC)=tn
  ann_colors=list(Type=c(genic_piG="#4daf4a",pathway="#984ea3",intergenic_piG="#E41A1C"),
                  Sample=c(Juvenile="#66c2a5",Adult_26nt="#fc8d62",
                           Adult_30nt="#8da0cb"))
  pheatmap(log10(list_marmo[[i]][tn,]+1),
           cluster_rows=F,cluster_cols=F,labels_row="",main=i,
           #cellwidth=0.7,cellheight=2,
           labels_col=c("-5kb",rep("",99),"TSS",rep("",98),"+5kb"),
           annotation_row=annoC,annotation_colors=ann_colors,
           col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100))
  #filename=paste("heatmap.amyb.piRNA",i,"pdf",sep=".")
}
dev.off()

# July12.2018: metaplot and heatmap for CAGE and PAS in piRNA genes----
setwd("/Users/BigBear/My_Publish/Human_piRNA/aggreMatrix/")
list_cage=list()
list_cage[["CAGE_1643_Rep1"]]=read.table("CAGE_1643_Rep1_Id2.t53.u5000.d5000.b200.mat",header=F,row.names=1)
list_cage[["CAGE_1643_Rep2"]]=read.table("CAGE_1643_Rep2_Id3.t53.u5000.d5000.b200.mat",header=F,row.names=1)
list_cage[["CAGE_1643_Rep3"]]=read.table("CAGE_1643_Rep3_Id1.t53.u5000.d5000.b200.mat",header=F,row.names=1)
list_cage[["CAGE_5431"]]=read.table("CAGE_5431_Id7.t53.u5000.d5000.b200.mat",header=F,row.names=1)
list_cage[["CAGE_7358"]]=read.table("CAGE_7358_Id6.t53.u5000.d5000.b200.mat",header=F,row.names=1)
list_cage[["CAGE_7377"]]=read.table("CAGE_7377_Id5.t53.u5000.d5000.b200.mat",header=F,row.names=1)
list_cage[["CAGE_8150"]]=read.table("CAGE_8150_Id4.t53.u5000.d5000.b200.mat",header=F,row.names=1)
list_pas=list()
list_pas[["PAS_ACTGT"]]=read.table("ACTGT.Dedup_Atrimmed.t53.u5000.d5000.b200.mat",header=F,row.names=1)
list_pas[["PAS_ATCTG"]]=read.table("ATCTG.Dedup_Atrimmed.t53.u5000.d5000.b200.mat",header=F,row.names=1)
list_pas[["PAS_TGACT"]]=read.table("TGACT.Dedup_Atrimmed.t53.u5000.d5000.b200.mat",header=F,row.names=1)
list_pas[["PAS_CGGGA"]]=read.table("CGGGA.Dedup_Atrimmed.t53.u5000.d5000.b200.mat",header=F,row.names=1)
# plot metaplot
fun_meta=function(m1,m2,mn,lp,cl){
  par(mar=c(2,4,4,1),bty="n",tcl=0.3,las=1)
  x1=apply(m1,2,mean,rm=0.05);x2=apply(m2,2,mean,rm=0.05)
  ym=max(c(x1,x2))
  plot(0:199,apply(m1,2,mean,rm=0.05),type="l",lwd=2,xaxt="n",xlab="",ylab="RPM",
       col=cl,ylim=c(0,ym),main=mn)
  lines(0:199,apply(m2,2,mean,rm=0.05),lwd=2,col=cl,lty=2)
  legend(lp,legend=c("adult","juvenil"),col="black",lwd=3,bty="n",lty=c(1,2))
  axis(1,c(1,51,151,200),label=c("-5kb","TSS","TES","+5kb"))
}
pdf("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/Bear_Jul12.2018/metaplot.cage_pas.pdf",width=6,height=4)
fun_meta(list_cage$CAGE_8150[human_pachy,],list_cage$CAGE_1643_Rep1[human_pachy,],
         "CAGE in Pachytene piRNA Genes","topright",csl[1])
fun_meta(list_cage$CAGE_8150[human_prepachy,],list_cage$CAGE_1643_Rep2[human_prepachy,],
         "CAGE in Prepachytene piRNA Genes","topright",csl[3])
fun_meta(list_cage$CAGE_8150[human_hybrid,],list_cage$CAGE_1643_Rep2[human_hybrid,],
         "CAGE in Hybrid piRNA Genes","topright",csl[4])
fun_meta(list_pas$PAS_ATCTG[human_pachy,],list_pas$PAS_ACTGT[human_pachy,],
         "PAS in Pachytene piRNA Genes","topright",csl[1])
fun_meta(list_pas$PAS_ATCTG[human_prepachy,],list_pas$PAS_ACTGT[human_prepachy,],
         "PAS in Prepachytene piRNA Genes","topright",csl[3])
fun_meta(list_pas$PAS_ATCTG[human_hybrid,],list_pas$PAS_ACTGT[human_hybrid,],
         "PAS in Hybrid piRNA Genes","topright",csl[4])
dev.off()
# plot heatmap
pdf("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/Bear_Jul12.2018/heatmap.cage_pas.pdf",width=3.5,height=7)
for(i in names(list_cage)){
  annoC=data.frame(Type=factor(rep(c("Pachytene","Hybrid","Prepachytene"),c(90,10,84))))
  row.names(annoC)=c(human_pachy,human_hybrid,human_prepachy)
  ann_colors=list(Type=c(Prepachytene="#4daf4a",Hybrid="#984ea3",Pachytene="#E41A1C"),
                  Sample=c(Juvenile="#66c2a5",Adult_26nt="#fc8d62",
                           Adult_30nt="#8da0cb"))
  pheatmap(log10(list_cage[[i]][c(human_pachy,human_hybrid,human_prepachy),]+1),
           cluster_rows=F,cluster_cols=F,labels_row="",main=i,
           #cellwidth=0.7,cellheight=2,
           labels_col=c("-5kb",rep("",49),"TSS",rep("",99),"TES",rep("",48),"+5kb"),
           annotation_row=annoC,annotation_colors=ann_colors,
           breaks=c((0:99/100)*1,2),
           col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100))
  #filename=paste("heatmap.amyb.piRNA",i,"pdf",sep=".")
}
for(i in names(list_pas)){
  annoC=data.frame(Type=factor(rep(c("Pachytene","Hybrid","Prepachytene"),c(90,10,84))))
  row.names(annoC)=c(human_pachy,human_hybrid,human_prepachy)
  ann_colors=list(Type=c(Prepachytene="#4daf4a",Hybrid="#984ea3",Pachytene="#E41A1C"),
                  Sample=c(Juvenile="#66c2a5",Adult_26nt="#fc8d62",
                           Adult_30nt="#8da0cb"))
  pheatmap(log10(list_pas[[i]][c(human_pachy,human_hybrid,human_prepachy),]+1),
           cluster_rows=F,cluster_cols=F,labels_row="",main=i,
           #cellwidth=0.7,cellheight=2,
           labels_col=c("-5kb",rep("",49),"TSS",rep("",99),"TES",rep("",48),"+5kb"),
           annotation_row=annoC,annotation_colors=ann_colors,
           #breaks=c((0:99/100)*0.3,2),
           col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100))
  #filename=paste("heatmap.amyb.piRNA",i,"pdf",sep=".")
}
dev.off()

# July12.2018: use quantsmooth to plot piRNA genes on each chromosome----
library(quantsmooth)
human_pi_bed=read.table("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/piRNA_cluster_annotation/final_cluster/piG.bed.nochr",header=F,row.names=4)
df1=data.frame(human_pi_bed[human_pachy,1],rowMeans(human_pi_bed[human_pachy,2:3]))
df2=data.frame(human_pi_bed[human_prepachy,1],rowMeans(human_pi_bed[human_prepachy,2:3]))
df3=data.frame(human_pi_bed[human_hybrid,1],rowMeans(human_pi_bed[human_hybrid,2:3]))
colnames(df1)=c("CHR","MapInfo")
colnames(df2)=c("CHR","MapInfo")
colnames(df3)=c("CHR","MapInfo")
chrompos1=prepareGenomePlot(df1,paintCytobands=TRUE,organism="hsa",sexChromosomes=T,
                           units="hg19",cols="grey50")
chrompos2=prepareGenomePlot(df2,paintCytobands=TRUE,organism="hsa",sexChromosomes=T,
                            units="hg19",cols="grey50")
pdf("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/Bear_Jul12.2018/karyoplot.pdf",width=6,height=5)
chrompos3=prepareGenomePlot(df3,paintCytobands=TRUE,organism="hsa",sexChromosomes=T,
                            units="hg19",cols="grey50",bleach=0,topspace=2)
points(chrompos1[,2],chrompos1[,1],col=csl[1],pch="l",font=1)
points(chrompos2[,2],chrompos2[,1]+0.15,col=csl[3],pch="l",font=1)
points(chrompos3[,2],chrompos3[,1]-0.15,col=csl[4],pch="l",font=1)
legend("top",legend=c("Pachytene","Prepachytene","Hybrid"),lwd=2,col=csl[c(1,3,4)],ncol=3,cex=5/6)
dev.off()

kp <-plotKaryotype()
ss=human_pi_bed;ss[,3]=ss[,2]
kpPlotMarkers(kp,data=toGRanges(ss[1:10,]),labels=row.names(human_pi_bed)[1:10])



# July12.2018: merged AMYB peak distance to all genes----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/")
pd=read.table("matrix/human.amybPeakMerged.distance.txt",header=F,row.names=1)
#pd=read.table("matrix/t",header=F,row.names=1)
#pd[which(pd[,1]==(-1)),1]=500000
# AMYB
pdf("Bear_Jul12.2018/peak_distance.amyb.pdf",width=4,height=4,useDingbats=F)
par(mar=c(4,4,1,1),cex=5/6,bty="n",xpd=T)
plot(NA,xlim=c(0,13),ylim=c(0,7),ylab="distance (bp)",xlab="",xaxt="n",yaxt="n")
ph=21
boxplot(log10(pd[human_pachy_amyb,1]+1),log10(pd[human_pachy_noamyb,1]+1),
        log10(pd[pl,1]+1),log10(pd[human_prepachy,1]+1),
        log10(pd[human_hybrid,1]+1),log10(pd[human_pc,1]+1),log10(pd[human_linc,1]+1),outline=F,
        at=c(1,3,5,7,9,11,13)-0.5,boxwex=1,staplewex=0,add=T,yaxt="n",xaxt="n")
#points(runif(length(l1),0+0.05,1-0.05),rnalog10(l1+1),col=csl[1],pch=ph,lwd=1,cex=0.2)
points(runif(length(pl),4+0.05,5-0.05),log10(pd[pl,1]+1),bg=csl[2],pch=ph,lwd=0.6,cex=0.7,col="white")
points(runif(length(human_pachy_amyb),0+0.05,1-0.05),log10(pd[human_pachy_amyb,1]+1),bg=csl[1],pch=ph,lwd=0.6,cex=0.7,col="white")
points(runif(length(human_pachy_noamyb),2+0.05,3-0.05),log10(pd[human_pachy_noamyb,1]+1),bg=csl[1],pch=ph,lwd=0.6,cex=0.7,col="white")
points(runif(length(human_prepachy),6+0.05,7-0.05),log10(pd[human_prepachy,1]+1),bg=csl[3],pch=ph,lwd=0.6,cex=0.7,col="white")
points(runif(length(human_hybrid),8+0.05,9-0.05),log10(pd[human_hybrid,1]+1),bg=csl[4],pch=ph,lwd=0.6,cex=0.7,col="white")
points(runif(length(human_pc),10+0.05,11-0.05),log10(pd[human_pc,1]+1),bg="black",pch=ph,lwd=0.6,cex=0.7,col="white")
points(runif(length(human_linc),12+0.05,13-0.05),log10(pd[human_linc,1]+1),bg="black",pch=ph,lwd=0.6,cex=0.7,col="white")
text(c(1,3,5,7,9,11,13)-0.5,-1,label=c("amyb-positive\npachytene\nclusters",
                                       "amyb-negative\npachytene\nclusters","piRNA\nbiogenesis",
                                 "prepachytene\nclusters","hybrid\nclusters","mRNA genes",
                                 "lincRNA genes"),
     xpd=T,srt=45,col=c(csl[c(1,1,2,3,4)],"black","black"))
axis(2,0:7,c(expression(10^0,10^1,10^2,10^3,10^4,10^5,10^6,10^7)))
dev.off()

# get gene list within amyb peaks
rn=row.names(pd)[which(pd[,1]<=500)]
write.table(intersect(rn,human_pc),"Bear_Jul12.2018/amyb_gene_list/protein_coding.txt",quote=F,row.names=F,col.names=F,sep="\t")
write.table(intersect(rn,human_linc),"Bear_Jul12.2018/amyb_gene_list/lincRNA.txt",quote=F,row.names=F,col.names=F,sep="\t")
write.table(intersect(rn,human_pachy),"Bear_Jul12.2018/amyb_gene_list/pachytene.txt",quote=F,row.names=F,col.names=F,sep="\t")
write.table(intersect(rn,human_prepachy),"Bear_Jul12.2018/amyb_gene_list/prepachytene.txt",quote=F,row.names=F,col.names=F,sep="\t")
write.table(intersect(rn,human_hybrid),"Bear_Jul12.2018/amyb_gene_list/hybrid.txt",quote=F,row.names=F,col.names=F,sep="\t")
write.table(intersect(rn,pl),"Bear_Jul12.2018/amyb_gene_list/pathway.txt",quote=F,row.names=F,col.names=F,sep="\t")

# July12.2018: buckets for amyb in pathway genes----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/")
bucket_testis=read.table("aggreMatrix/matrix_rpm/merged.AMYB.t5.u50000.d100000.b10000.mat",header=F,row.names=1)
bucket_kidney=read.table("aggreMatrix/matrix_rpm/kidney.AMYB.t5.u50000.d100000.b10000.mat",header=F,row.names=1)
pdf("Bear_Jul12.2018/coveragePlot.pdf",width=5,height=3)
barplot(as.matrix(bucket_testis["PIWIL1",3134:5650]),space=0,col="blue",
        border="blue",xaxt="n",main="HIWI",ylab="RPM",ylim=c(0,0.5))
axis(1,c(1,200,2516),label=c("-3000bp","TSS","34750bp"))
barplot(as.matrix(bucket_testis["MYBL1",3134:6939]),space=0,col="blue",
        border="blue",xaxt="n",main="AMYB",ylab="RPM",ylim=c(0,0.5))
axis(1,c(1,200,3805),label=c("-3000bp","TSS","52072bp"))
barplot(as.matrix(bucket_testis["15-q22-56093",1477:7219]),space=0,col="blue",
        border="blue",xaxt="n",main="15-q22",ylab="RPM",ylim=c(0,0.25))
axis(1,c(1,1878,5743),label=c("-28174bp","TSS","58286bp"))
barplot(as.matrix(bucket_testis["pi-TMEM181",3134:9933]),space=0,col="blue",
        border="blue",xaxt="n",main="pi-TMEM181",ylab="RPM",ylim=c(0,0.5))
axis(1,c(1,200,6799),label=c("-3000bp","TSS","98993bp"))
dev.off()

# setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/")
# bucket_amyb=read.table("aggreMatrix/matrix_rpm/hg19.pathway.amyb.t53.u5000.d5000.b2000.mat",header=F,row.names=1)
# bucket_input=read.table("aggreMatrix/matrix_rpm/hg19.pathway.amyb.kidney.t53.u5000.d5000.b2000.mat",header=F,row.names=1)
# pdf("Bear_Jul12.2018/AMYB.buckets.pathway.pdf",width=5,height=4)
# par(mfrow=c(2,1))
# for(i in row.names(bucket_amyb)){
#   par(mar=c(1,4,3,1))
#   ym=max(c(max(bucket_amyb[i,]),max(bucket_input[i,]/100)))
#   barplot(as.matrix(bucket_amyb[i,]),xaxt="n",ylab="RPM",col=csl[2],border=csl[2],
#           ylim=c(0,ym),main=i,space=0)
#   #axis(1,c(1,501,1501,2000),label=c("-5kb","TSS","TES","+5kb"))
#   par(mar=c(4,4,0,1))
#   barplot(as.matrix(bucket_input[i,]/100),xaxt="n",ylab="RPM",col="grey",
#           border="grey",ylim=c(0,ym),space=0)
#   axis(1,c(1,501,1501,2000),label=c("-5kb","TSS","TES","+5kb"))
# }
# dev.off()

# July12.2018: buckets for amyb in pathway genes for rhesus----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/")
bucket_mouse1=read.table("aggreMatrix/matrix_rpm/mouse_3MO-2-AMYB-Id8_S7.rpm.uniq.AMYB_HIWI.t5.u50000.d100000.b10000.mat",header=F,row.names=1)
bucket_mouse2=read.table("aggreMatrix/matrix_rpm/mouse_3MO-2-AMYB-input-Id10_S3.rpm.uniq.AMYB_HIWI.t5.u50000.d100000.b10000.mat",header=F,row.names=1)
bucket_rhesus1=read.table("aggreMatrix/matrix_rpm/Macaque-AMYB-ID5_S8.rpm.uniq.AMYB_HIWI.t5.u50000.d100000.b10000.mat",header=F,row.names=1)
bucket_rhesus2=read.table("aggreMatrix/matrix_rpm/Macaque-Input-ID8_S7.rpm.uniq.AMYB_HIWI.t5.u50000.d100000.b10000.mat",header=F,row.names=1)
pdf("Bear_Jul12.2018/AMYB_PIWI1.coveragePlot.rhesus.pdf",width=5,height=6)
par(mar=c(2,4,3,1),mfrow=c(2,1),tcl=0.3,cex=5/6)
# barplot(as.matrix(bucket_mouse1["Piwil1",3134:5650])/500000,space=0,col="blue",
#         border="blue",xaxt="n",main="mouse PIWIL1",ylab="RPM",ylim=c(0,0.5))
# axis(1,c(1,200,2516),label=c("-3000bp","TSS","34750bp"))
# barplot(as.matrix(bucket_mouse2["Piwil1",3134:5650])/500000,space=0,col="blue",
#         border="blue",xaxt="n",main="mouse input",ylab="RPM",ylim=c(0,0.5))
# axis(1,c(1,200,2516),label=c("-3000bp","TSS","34750bp"))
# barplot(as.matrix(bucket_mouse1["Mybl1",3134:6939])/43/43,space=0,col="blue",
#         border="blue",xaxt="n",main="mouse AMYB",ylab="RPM",ylim=c(0,10))
# axis(1,c(1,200,3805),label=c("-3000bp","TSS","52072bp"))
# barplot(as.matrix(bucket_mouse2["Mybl1",3134:6939])/43/43,space=0,col="blue",
#         border="blue",xaxt="n",main="mouse input",ylab="RPM",ylim=c(0,10))
# axis(1,c(1,200,3805),label=c("-3000bp","TSS","52072bp"))
barplot(as.matrix(bucket_rhesus1["ENSMMUG00000013974",3134:6590])/43/43,space=0,col="blue",
        border="blue",xaxt="n",main="rhesus PIWIL1",ylab="RPM",ylim=c(0,1))
axis(1,c(1,200,2516),label=c("-3000bp","TSS","48851bp"))
barplot(as.matrix(bucket_rhesus2["ENSMMUG00000013974",3134:6590])/43/43,space=0,col="blue",
        border="blue",xaxt="n",main="rhesus input",ylab="RPM",ylim=c(0,1))
axis(1,c(1,200,2516),label=c("-3000bp","TSS","48851p"))
barplot(as.matrix(bucket_rhesus1["ENSMMUG00000005629",3134:5294])/43/43,space=0,col="blue",
        border="blue",xaxt="n",main="rhesus AMYB",ylab="RPM",ylim=c(0,1))
axis(1,c(1,200,3805),label=c("-3000bp","TSS","29421bp"))
barplot(as.matrix(bucket_rhesus2["ENSMMUG00000005629",3134:5294])/43/43,space=0,col="blue",
        border="blue",xaxt="n",main="rhesus input",ylab="RPM",ylim=c(0,1))
axis(1,c(1,200,3805),label=c("-3000bp","TSS","29421bp"))
dev.off()

# July18.2018: H3K4me3 heatmap for piRNA genes----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/aggreMatrix/matrix_rpm/")
l_hs=c("5431-H3K4me3-Id2_S9","5431-input-H3K4me3-Id5_S5","7358-H3K4me3-Id1_S2",
       "7358-input-H3K4me3-Id4_S8","8246NT-H3K4me3-Id10_S6","8246NT-H3K4me3-input-Id12_S7",
       "Hs-Kidney-H3K4me3-ID1_S9","Hs-Liver-H3K4me3-ID3_S6",
       "Hs-Liver-H3K4me3-input-ID4_S4","Hs-Spleen-H3K4me3-ID2_S5")
list_hs=list()
for(i in l_hs){
  list_hs[[i]]=read.table(paste(i,"rpm.uniq.t5.u5000.d5000.b200.mat",sep="."),
                          header=F,row.names=1)
}
pdf("../figures_rpm/human_all.H3K4me3.pdf",width=3.5,height=7)
for(i in l_hs){
  annoC=data.frame(Type=factor(rep(c("Pachytene","Hybrid","Prepachytene"),c(90,10,84))))
  row.names(annoC)=c(human_pachy,human_hybrid,human_prepachy)
  ann_colors=list(Type=c(Prepachytene="#4daf4a",Hybrid="#984ea3",Pachytene="#E41A1C"),
                  Sample=c(Juvenile="#66c2a5",Adult_26nt="#fc8d62",
                           Adult_30nt="#8da0cb"))
  pheatmap(log10(list_hs[[i]][c(human_pachy,human_hybrid,human_prepachy),]+1),
           cluster_rows=F,cluster_cols=F,labels_row="",main=i,
           #cellwidth=0.7,cellheight=2,
           labels_col=c("-5kb",rep("",99),"TSS",rep("",98),"+5kb"),
           annotation_row=annoC,annotation_colors=ann_colors,
           breaks=c((0:99/100)*0.5,0.8),
           col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100))
  #filename=paste("heatmap.amyb.piRNA",i,"pdf",sep=".")
}
dev.off()

# July19.2018: boxplot of rnaseq, srnaseq for piRNA genes----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
df_srna=read.table("merged.rpm",header=T,row.names=1,check.names=F)
df_rna=read.table("merged.rnaseq.gene+pi.rpkm",header=T,row.names=1,check.names=F)
df_rna=df_rna[,rna_rn]
pdf("../Bear_Jul12.2018/boxplot.rna.srna.ForpiRNA.pdf",width=5,height=5.5,useDingbats=F)
par(mfrow=c(2,1),mar=c(4,4,1,1),tcl=0.3,bty="n",cex=5/6)
boxplot(log10(df_srna[human_pachy,srna_rn[1]]+1),log10(df_srna[human_prepachy,srna_rn[1]]+1),
        log10(df_srna[human_pachy,srna_rn[3]]+1),log10(df_srna[human_prepachy,srna_rn[3]]+1),
        log10(df_srna[human_pachy,srna_rn[4]]+1),log10(df_srna[human_prepachy,srna_rn[4]]+1),
        log10(df_srna[human_pachy,srna_rn[7]]+1),log10(df_srna[human_prepachy,srna_rn[7]]+1),
        log10(df_srna[human_pachy,srna_rn[10]]+1),log10(df_srna[human_prepachy,srna_rn[10]]+1),
        log10(df_srna[human_pachy,srna_rn[15]]+1),log10(df_srna[human_prepachy,srna_rn[15]]+1),
        col=csl[c(1,3)],at=c(1,2,4,5,7,8,10,11,13,14,16,17),pch=20,cex=0.3,xaxt="n",
        ylab="piRNA RPM",yaxt="n")
axis(2,0:4,label=c(0,10,100,1000,"10,000"))
boxplot(log10(df_rna[human_pachy,rna_rn[1]]+1),log10(df_rna[human_prepachy,rna_rn[1]]+1),
        log10(df_rna[human_pachy,rna_rn[3]]+1),log10(df_rna[human_prepachy,rna_rn[3]]+1),
        log10(df_rna[human_pachy,rna_rn[4]]+1),log10(df_rna[human_prepachy,rna_rn[4]]+1),
        log10(df_rna[human_pachy,rna_rn[5]]+1),log10(df_rna[human_prepachy,rna_rn[5]]+1),
        log10(df_rna[human_pachy,rna_rn[10]]+1),log10(df_rna[human_prepachy,rna_rn[10]]+1),
        log10(df_rna[human_pachy,rna_rn[15]]+1),log10(df_rna[human_prepachy,rna_rn[15]]+1),
        col=csl[c(1,3)],at=c(1,2,4,5,7,8,10,11,13,14,16,17),pch=20,cex=0.3,xaxt="n",
        ylab="Transcript RPKM",yaxt="n")
axis(1,c(1.5,4.5,7.5,10.5,13.5,16.5),label=c("5YO","12YO","13YO","26nt","inter","30nt"),lwd=0)
axis(2,0:4,label=c(0,10,100,1000,"10,000"))
dev.off()
pdf("../Bear_Jul12.2018/boxplot.rna.srna.ForpiRNA.linear.pdf",width=5,height=5.5,useDingbats=F)
par(mfrow=c(2,1),mar=c(4,4,1,1),tcl=0.3,bty="n",cex=5/6)
boxplot(df_srna[human_pachy,srna_rn[1]],df_srna[human_prepachy,srna_rn[1]],
        df_srna[human_pachy,srna_rn[3]],df_srna[human_prepachy,srna_rn[3]],
        df_srna[human_pachy,srna_rn[4]],df_srna[human_prepachy,srna_rn[4]],
        df_srna[human_pachy,srna_rn[7]],df_srna[human_prepachy,srna_rn[7]],
        df_srna[human_pachy,srna_rn[10]],df_srna[human_prepachy,srna_rn[10]],
        df_srna[human_pachy,srna_rn[15]],df_srna[human_prepachy,srna_rn[15]],
        col=csl[c(1,3)],at=c(1,2,4,5,7,8,10,11,13,14,16,17),pch=20,cex=0.3,xaxt="n",
        ylab="piRNA RPM",ylim=c(0,20000))
boxplot(df_rna[human_pachy,rna_rn[1]],df_rna[human_prepachy,rna_rn[1]],
        df_rna[human_pachy,rna_rn[3]],df_rna[human_prepachy,rna_rn[3]],
        df_rna[human_pachy,rna_rn[4]],df_rna[human_prepachy,rna_rn[4]],
        df_rna[human_pachy,rna_rn[5]],df_rna[human_prepachy,rna_rn[5]],
        df_rna[human_pachy,rna_rn[10]],df_rna[human_prepachy,rna_rn[10]],
        df_rna[human_pachy,rna_rn[15]],df_rna[human_prepachy,rna_rn[15]],
        col=csl[c(1,3)],at=c(1,2,4,5,7,8,10,11,13,14,16,17),pch=20,cex=0.3,xaxt="n",
        ylab="Transcript RPKM",ylim=c(0,30))
axis(1,c(1.5,4.5,7.5,10.5,13.5,16.5),label=c("5YO","12YO","13YO","26nt","inter","30nt"),lwd=0)
dev.off()
# heatmap for correlation of rnaseqCluster
annoR=data.frame(Length_Type=factor(rep(c("Juv","l26nt","inter","l30nt"),c(3,3,5,6))))
rna_rn1=c("D1643_5YO",rna_rn[3:18])
row.names(annoR)=rna_rn1
df_rnat=cbind(apply(df_rna[,1:2],1,mean),df_rna[,-c(1,2)])
colnames(df_rnat)[1]="D1643_5YO"
ann_colors=list(Length_Type=c(Juv="#8dd3c7",l26nt="#fccde5",inter="#b3de69",l30nt="#fdb462"))
pheatmap(cor(df_rnat[setdiff(row.names(df_rnat),human_pi),rna_rn1],
             method="spearman"),annotation_row=annoR,annotation_col=annoR,
         annotation_colors=ann_colors,
         filename="../Bear_Jul12.2018/heatmap.Cor.rnaseqCluster.pdf",
         cellwidth=24,cellheight=24,main="gene expression correlation between samples")
pheatmap(cor(df_rnat[human_pi,rna_rn1],
             method="spearman"),annotation_row=annoR,annotation_col=annoR,
         annotation_colors=ann_colors,
         filename="../Bear_Jul12.2018/heatmap.Cor.rnaseqCluster.basedOnpiRNA.pdf",
         cellwidth=24,cellheight=24,main="gene expression correlation between samples")

# boxplot of pathway gene, pachy, prepachy expression for all samples
pdf("../Bear_Jul12.2018/boxplot.rna.srna.ForPathwayAndPi.pdf",width=8,height=4)
par(mar=c(8,4,2,2),bty="n")
boxplot(log10(df_rna[pl,c(9,10,20,14,17,19,12,21,8,11,13,15,18,16,1,2,3,4,5,6,7)]+1),
        xaxt="n",pch=20,cex=0.2,col=csl[2],at=(1:21),xlim=c(1,22),boxwex=0.2)
boxplot(log10(df_rna[human_pachy,c(9,10,20,14,17,19,12,21,8,11,13,15,18,16,1,2,3,4,5,6,7)]+1),
        xaxt="n",pch=20,cex=0.2,col=csl[1],at=(1:21)+0.2,boxwex=0.2,add=T)
boxplot(log10(df_rna[human_prepachy,c(9,10,20,14,17,19,12,21,8,11,13,15,18,16,1,2,3,4,5,6,7)]+1),
        xaxt="n",pch=20,cex=0.2,col=csl[3],at=(1:21)+0.4,boxwex=0.2,add=T)
axis(1,(1:21)+0.2,label=colnames(df_rna)[c(9,10,20,14,17,19,12,21,8,11,13,15,18,16,1,2,3,4,5,6,7)],
     las=3,lwd=0)
dev.off()

# boxplot for AMYB and HIWI
ma1=c(max(df_rna["MYBL1",1:4]),max(df_rna["MYBL1",c(5:7,9)]),max(df_rna["MYBL1",c(8,10,11,12)]),max(df_rna["MYBL1",13:18]))
ma2=c(max(df_rna["PIWIL1",1:4]),max(df_rna["PIWIL1",c(5:7,9)]),max(df_rna["PIWIL1",c(8,10,11,12)]),max(df_rna["PIWIL1",13:18]))
mi1=c(min(df_rna["MYBL1",1:4]),min(df_rna["MYBL1",c(5:7,9)]),min(df_rna["MYBL1",c(8,10,11,12)]),min(df_rna["MYBL1",13:18]))
mi2=c(min(df_rna["PIWIL1",1:4]),min(df_rna["PIWIL1",c(5:7,9)]),min(df_rna["PIWIL1",c(8,10,11,12)]),min(df_rna["PIWIL1",13:18]))
ma3=c(max(df_rna["TCFL5",1:4]),max(df_rna["TCFL5",c(5:7,9)]),max(df_rna["TCFL5",c(8,10,11,12)]),max(df_rna["TCFL5",13:18]))
mi3=c(min(df_rna["TCFL5",1:4]),min(df_rna["TCFL5",c(5:7,9)]),min(df_rna["TCFL5",c(8,10,11,12)]),min(df_rna["TCFL5",13:18]))
ma4=c(max(df_rna["PIWIL2",1:4]),max(df_rna["PIWIL2",c(5:7,9)]),max(df_rna["PIWIL2",c(8,10,11,12)]),max(df_rna["PIWIL2",13:18]))
mi4=c(min(df_rna["PIWIL2",1:4]),min(df_rna["PIWIL2",c(5:7,9)]),min(df_rna["PIWIL2",c(8,10,11,12)]),min(df_rna["PIWIL2",13:18]))
pdf("../../../Ozata et al. 1st Paper/Rebutal Cell/Figure3C.pdf",width=5,height=5)
par(mfrow=c(2,1),mar=c(1,4,1,1),bty="n",tcl=0.3)
plot(NA,xlim=c(0.5,4.5),ylim=c(0,max(df_rna[c("MYBL1"),rna_rn])),
     xlab="",ylab="AMYB",xaxt="n")
for(i in 1:4){arrows(i,ma1[i],i,mi1[i],length=0.3,angle=90,code=3)}
for(i in 1:4){a=runif(1,0.8,1.2);points(a,df_rna["MYBL1",i],pch=21,col="white",bg="black");text(a,df_rna["MYBL1",i],label=rna_rn[i],cex=0.6)}
for(i in c(5:7,9)){a=runif(1,1.8,2.2);points(a,df_rna["MYBL1",i],pch=21,col="white",bg="black");text(a,df_rna["MYBL1",i],label=rna_rn[i],cex=0.6)}
for(i in c(8,10,11,12)){a=runif(1,2.8,3.2);points(a,df_rna["MYBL1",i],pch=21,col="white",bg="black");text(a,df_rna["MYBL1",i],label=rna_rn[i],cex=0.6)}
for(i in 13:18){a=runif(1,3.8,4.2);points(a,df_rna["MYBL1",i],pch=21,col="white",bg="black");text(a,df_rna["MYBL1",i],label=rna_rn[i],cex=0.6)}
plot(NA,xlim=c(0.5,4.5),ylim=c(0,max(df_rna[c("PIWIL1"),rna_rn])),
     xlab="",ylab="HIWI",xaxt="n")
for(i in 1:4){arrows(i,ma2[i],i,mi2[i],length=0.3,angle=90,code=3)}
for(i in 1:4){a=runif(1,0.8,1.2);points(a,df_rna["PIWIL1",i],pch=21,col="white",bg="black");text(a,df_rna["PIWIL1",i],label=rna_rn[i],cex=0.6)}
for(i in c(5:7,9)){a=runif(1,1.8,2.2);points(a,df_rna["PIWIL1",i],pch=21,col="white",bg="black");text(a,df_rna["PIWIL1",i],label=rna_rn[i],cex=0.6)}
for(i in c(8,10,11,12)){a=runif(1,2.8,3.2);points(a,df_rna["PIWIL1",i],pch=21,col="white",bg="black");text(a,df_rna["PIWIL1",i],label=rna_rn[i],cex=0.6)}
for(i in 13:18){a=runif(1,3.8,4.2);points(a,df_rna["PIWIL1",i],pch=21,col="white",bg="black");text(a,df_rna["PIWIL1",i],label=rna_rn[i],cex=0.6)}
axis(1,1:4,label=c("Juv","26nt","inter","30nt"),lwd=0)
plot(NA,xlim=c(0.5,4.5),ylim=c(0,max(df_rna[c("TCFL5"),rna_rn])),
     xlab="",ylab="TCFL5",xaxt="n")
for(i in 1:4){arrows(i,ma3[i],i,mi3[i],length=0.3,angle=90,code=3)}
for(i in 1:4){a=runif(1,0.8,1.2);points(a,df_rna["TCFL5",i],pch=21,col="white",bg="black");text(a,df_rna["TCFL5",i],label=rna_rn[i],cex=0.6)}
for(i in c(5:7,9)){a=runif(1,1.8,2.2);points(a,df_rna["TCFL5",i],pch=21,col="white",bg="black");text(a,df_rna["TCFL5",i],label=rna_rn[i],cex=0.6)}
for(i in c(8,10,11,12)){a=runif(1,2.8,3.2);points(a,df_rna["TCFL5",i],pch=21,col="white",bg="black");text(a,df_rna["TCFL5",i],label=rna_rn[i],cex=0.6)}
for(i in 13:18){a=runif(1,3.8,4.2);points(a,df_rna["TCFL5",i],pch=21,col="white",bg="black");text(a,df_rna["TCFL5",i],label=rna_rn[i],cex=0.6)}
axis(1,1:4,label=c("Juv","26nt","inter","30nt"),lwd=0)
plot(NA,xlim=c(0.5,4.5),ylim=c(0,max(df_rna[c("PIWIL2"),rna_rn])),
     xlab="",ylab="HILI",xaxt="n")
for(i in 1:4){arrows(i,ma4[i],i,mi4[i],length=0.3,angle=90,code=3)}
for(i in 1:4){a=runif(1,0.8,1.2);points(a,df_rna["PIWIL2",i],pch=21,col="white",bg="black");text(a,df_rna["PIWIL1",i],label=rna_rn[i],cex=0.6)}
for(i in c(5:7,9)){a=runif(1,1.8,2.2);points(a,df_rna["PIWIL2",i],pch=21,col="white",bg="black");text(a,df_rna["PIWIL1",i],label=rna_rn[i],cex=0.6)}
for(i in c(8,10,11,12)){a=runif(1,2.8,3.2);points(a,df_rna["PIWIL2",i],pch=21,col="white",bg="black");text(a,df_rna["PIWIL1",i],label=rna_rn[i],cex=0.6)}
for(i in 13:18){a=runif(1,3.8,4.2);points(a,df_rna["PIWIL2",i],pch=21,col="white",bg="black");text(a,df_rna["PIWIL1",i],label=rna_rn[i],cex=0.6)}
axis(1,1:4,label=c("Juv","26nt","inter","30nt"),lwd=0)
dev.off()

# boxplot for all pathway genes
pdf("../../../Ozata et al. 1st Paper/Rebutal Cell/Figure3B.pdf",width=4,height=4,useDingbats=F)
par(mfrow=c(1,1),tcl=0.3,mar=c(6,4,1,1),bty="n")
boxplot(as.vector(unlist(df_rna[pl,1:4])),as.vector(unlist(df_rna[pl,13:18])),
        as.vector(unlist(df_rna[pl,c(8,10,11,12)])),as.vector(unlist(df_rna[pl,c(5:7,9)])),
        col=tcs[c(1,5,8,13)],cex=0.3,pch=20,ylab="RPKM",names=c("Juv","30nt","inter","26nt"))
#boxplot(df_rna[pl,1:18],at=c(1:4,7:9,12:16,19:24),col=tcs,pch=20,cex=0.3,las=3,ylab="RPKM")
dev.off()

# July19.2018: syntenic analysis----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
syn=read.table("syntenic/hg19.sytenicAndPi.mat",header=T,row.names=1)
o1=order(apply(syn[human_pachy,],1,mean),decreasing=T)
o2=order(apply(syn[human_hybrid,],1,mean),decreasing=T)
o3=order(apply(syn[human_prepachy,],1,mean),decreasing=T)
annoC=data.frame(Type=factor(rep(c("Pachytene","Hybrid","Prepachytene"),c(89,10,83))))
row.names(annoC)=c(human_pachy,human_hybrid,human_prepachy)
ann_colors=list(Type=c(Prepachytene="#4daf4a",Hybrid="#984ea3",Pachytene="#E41A1C"),
                Sample=c(Juvenile="#66c2a5",Adult_26nt="#fc8d62",
                         Adult_30nt="#8da0cb"))
# heatmap
pheatmap(as.matrix(syn[c(human_pachy[o1],human_hybrid[o2],human_prepachy[o3]),
                        c(6,1,2,5,3,4,7)]),cluster_cols=F,cluster_rows=F,
         annotation_row=annoC,annotation_colors=ann_colors,col=c("grey","#9ecae1","#3182bd"),
         #filename="../Bear_Jul12.2018/syntenic.heatmap.pdf",
         cellwidth=16,cellheight=4,legend_labels=c("not_syntenic","syntenic","conserved"),
         legend_breaks=c(0,1,4),fontsize_row=4,breaks=c(0,0.9,1.9,4))
# define eutherian-conservation...
for(i in 1:ncol(syn)){
  syn[which(syn[,i]==1),i]=0
}
human_euthe=human_pachy[which(apply(syn[human_pachy,c(1,6)],1,sum)>0 & apply(syn[human_pachy,c(2,5)],1,sum)>0)]
human_primate=human_pachy[which(apply(syn[human_pachy,c(1,6)],1,sum)>0)]
human_primate=setdiff(human_primate,human_euthe)
human_non=setdiff(human_pachy,c(human_euthe,human_primate))
write.table(human_euthe,"../../Metagene_Plots/matrix/pachytene_eutherian.txt",col.names=F,row.names=F,sep="\t",quote=F)
write.table(human_primate,"../../Metagene_Plots/matrix/pachytene_primates.txt",col.names=F,row.names=F,sep="\t",quote=F)
write.table(human_non,"../../Metagene_Plots/matrix/pachytene_noncon.txt",col.names=F,row.names=F,sep="\t",quote=F)

# chromosome plot
cs_hs=read.table("syntenic/hg19.chrom.size",header=F,row.names=1,stringsAsFactors=F)/100
cs_mm=read.table("syntenic/mm10.chrom.size",header=F,row.names=1,stringsAsFactors=F)/100
cs_rm=read.table("syntenic/rheMac8.chrom.size",header=F,row.names=1,stringsAsFactors=F)/100
bed_hs=read.table("syntenic/piG.bed",header=F,row.names=4,stringsAsFactors=F)
bed_mm=read.table("syntenic/hg19.mm10.piRNA.overlap",header=F,row.names=4,stringsAsFactors=F)
bed_rm=read.table("syntenic/hg19.rheMac8.piRNA.overlap",header=F,row.names=4,stringsAsFactors=F)
pdf("../Bear_Jul12.2018/karyoplot.syntenic.pdf",width=5,height=15)
par(mar=c(2,0,4,0),bty="n",lwd=2)
plot(NA,xlim=c(0.5,3.5),ylim=c(0,max(c(cs_hs["chrY",2],cs_mm["chrY",2],cs_rm["chrY",2]))),
     xlab="",ylab="",xaxt="n",yaxt="n",main="conservation of piRNA genes")
cs=rep(c("#d9d9d9","#525252"),20);cs[23:24]=c("#fcbba1","#c6dbef")
for(i in 1:dim(cs_hs)[1]){
  polygon(c(1.98,1.98,2.02,2.02),c(cs_hs[i,1],cs_hs[i,2],cs_hs[i,2],cs_hs[i,1]),
          col=cs[i],border=cs[i])
}
cs=rep(c("#d9d9d9","#525252"),20);cs[20:21]=c("#fcbba1","#c6dbef")
for(i in 1:dim(cs_mm)[1]){
  if(i==1){p=0}else{p=sum(cs_mm[1:(i-1),1])}
  polygon(c(2.98,2.98,3.02,3.02),c(cs_mm[i,1],cs_mm[i,2],cs_mm[i,2],cs_mm[i,1]),
          col=cs[i],border=cs[i])
}
cs=rep(c("#d9d9d9","#525252"),20);cs[21:22]=c("#fcbba1","#c6dbef")
for(i in 1:dim(cs_rm)[1]){
  if(i==1){p=0}else{p=sum(cs_rm[1:(i-1),1])}
  polygon(c(0.98,0.98,1.02,1.02),c(cs_rm[i,1],cs_rm[i,2],cs_rm[i,2],cs_rm[i,1]),
          col=cs[i],border=cs[i])
}
mtext("hg19",1);mtext("rheMac8",1,at=1);mtext("mm10",1,at=3)
for(i in human_hybrid){
  lines(c(1.97,2.03),c(cs_hs[bed_hs[i,1],1]+bed_hs[i,2]/100,cs_hs[bed_hs[i,1],1]+bed_hs[i,2]/100),
        col=csl[4],lwd=1)
}
for(i in human_prepachy){
  lines(c(1.97,2.03),c(cs_hs[bed_hs[i,1],1]+bed_hs[i,2]/100,cs_hs[bed_hs[i,1],1]+bed_hs[i,2]/100),
        col=csl[3],lwd=1)
}
for(i in human_pachy){
  lines(c(1.97,2.03),c(cs_hs[bed_hs[i,1],1]+bed_hs[i,2]/100,cs_hs[bed_hs[i,1],1]+bed_hs[i,2]/100),
        col=csl[1],lwd=1)
}
for(i in human_hybrid){
  lines(c(1.97,2.03),c(cs_hs[bed_hs[i,1],1]+bed_hs[i,2]/100,cs_hs[bed_hs[i,1],1]+bed_hs[i,2]/100),
        col=csl[4],lwd=1)
  lines(c(2.08,2.92),c(cs_hs[bed_hs[i,1],1]+bed_hs[i,2]/100,cs_mm[bed_mm[i,1],1]+bed_mm[i,2]/100),
        col=csl[4])
  lines(c(2.97,3.03),c(cs_mm[bed_mm[i,1],1]+bed_mm[i,2]/100,cs_mm[bed_mm[i,1],1]+bed_mm[i,2]/100),
        col=csl[4],lwd=1)
  lines(c(1.92,1.08),c(cs_hs[bed_hs[i,1],1]+bed_hs[i,2]/100,cs_rm[bed_rm[i,1],1]+bed_rm[i,2]/100),
        col=csl[4])
  lines(c(0.97,1.03),c(cs_rm[bed_rm[i,1],1]+bed_rm[i,2]/100,cs_rm[bed_rm[i,1],1]+bed_rm[i,2]/100),
        col=csl[4],lwd=1)
}
for(i in human_prepachy){
  lines(c(1.97,2.03),c(cs_hs[bed_hs[i,1],1]+bed_hs[i,2]/100,cs_hs[bed_hs[i,1],1]+bed_hs[i,2]/100),
        col=csl[3],lwd=1)
  lines(c(2.08,2.92),c(cs_hs[bed_hs[i,1],1]+bed_hs[i,2]/100,cs_mm[bed_mm[i,1],1]+bed_mm[i,2]/100),
        col=csl[3])
  lines(c(2.97,3.03),c(cs_mm[bed_mm[i,1],1]+bed_mm[i,2]/100,cs_mm[bed_mm[i,1],1]+bed_mm[i,2]/100),
        col=csl[3],lwd=1)
  lines(c(1.92,1.08),c(cs_hs[bed_hs[i,1],1]+bed_hs[i,2]/100,cs_rm[bed_rm[i,1],1]+bed_rm[i,2]/100),
        col=csl[3])
  lines(c(0.97,1.03),c(cs_rm[bed_rm[i,1],1]+bed_rm[i,2]/100,cs_rm[bed_rm[i,1],1]+bed_rm[i,2]/100),
        col=csl[3],lwd=1)
}
for(i in human_pachy){
  lines(c(1.97,2.03),c(cs_hs[bed_hs[i,1],1]+bed_hs[i,2]/100,cs_hs[bed_hs[i,1],1]+bed_hs[i,2]/100),
        col=csl[1],lwd=1)
  lines(c(2.08,2.92),c(cs_hs[bed_hs[i,1],1]+bed_hs[i,2]/100,cs_mm[bed_mm[i,1],1]+bed_mm[i,2]/100),
        col=csl[1])
  lines(c(2.97,3.03),c(cs_mm[bed_mm[i,1],1]+bed_mm[i,2]/100,cs_mm[bed_mm[i,1],1]+bed_mm[i,2]/100),
        col=csl[1],lwd=1)
  lines(c(1.92,1.08),c(cs_hs[bed_hs[i,1],1]+bed_hs[i,2]/100,cs_rm[bed_rm[i,1],1]+bed_rm[i,2]/100),
        col=csl[1])
  lines(c(0.97,1.03),c(cs_rm[bed_rm[i,1],1]+bed_rm[i,2]/100,cs_rm[bed_rm[i,1],1]+bed_rm[i,2]/100),
        col=csl[1],lwd=1)
}
dev.off()


# Figure 1b,c----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
mat_ppm=read.table("merged.rpm",header=T,row.names=1,check.names=F)
df_rna=read.table("merged.rnaseq.gene+pi.rpkm",header=T,row.names=1,check.names=F)
mat_ppm=mat_ppm[,srna_rn]
mat_ppm1=mat_ppm[,c(1,2,5,6,9,8,14,13,15)]
mat_ppm2=mat_ppm[,setdiff(1:18,c(1,2,5,6,9,8,14,13,15))]
# Figure 1b
annoC=data.frame(Type=factor(rep(c("Pachytene","Hybrid","Prepachytene"),c(89,10,83))))
annoR=data.frame(Sample=factor(rep(c("Juv","l26nt","inter","l30nt"),c(2,2,2,3))))
ann_colors=list(Type=c(Prepachytene="#4daf4a",Hybrid="#984ea3",Pachytene="#E41A1C"),
                Sample=c(Juv="#8dd3c7",l26nt="#fccde5",inter="#b3de69",l30nt="#fdb462"))
o1=order(apply(mat_ppm1[human_pachy,],1,mean),decreasing=T)
o2=order(apply(mat_ppm1[human_hybrid,],1,mean),decreasing=T)
o3=order(apply(mat_ppm1[human_prepachy,],1,mean),decreasing=T)
mt=mat_ppm1[c(human_pachy[o1],human_hybrid[o2],human_prepachy[o3]),]
row.names(annoR)=colnames(mt)
row.names(annoC)=row.names(mt)
pheatmap(log10(mt+1),
         cluster_rows=F,cluster_cols=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         cellwidth=16,cellheight=4,border="white",lwd=0.3,fontsize_row=4,fontsize_col=16,
         filename="../Bear_Jul12.2018/Figure1.srna.heatmap.pdf",
         annotation_row=annoC,annotation_col=annoR,annotation_colors=ann_colors)
mt=df_rna[c(human_pachy,human_hybrid,human_prepachy),rna_rn[c(1,2,5,6,9,8,14,13,15)]]
row.names(annoR)=colnames(mt)
pheatmap(log10(mt+1),
         cluster_rows=F,cluster_cols=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         cellwidth=16,cellheight=4,border="white",lwd=0.3,fontsize_row=4,fontsize_col=16,
         filename="../Bear_Jul12.2018/Figure1.rna.heatmap.pdf",
         annotation_row=annoC,annotation_col=annoR,annotation_colors=ann_colors)
# Figure 1b additional
annoC=data.frame(Type=factor(rep(c("Pachytene","Hybrid","Prepachytene"),c(89,10,83))))
annoR=data.frame(Sample=factor(rep(c("Juv","l26nt","inter","l30nt"),c(2,1,3,3))))
ann_colors=list(Type=c(Prepachytene="#4daf4a",Hybrid="#984ea3",Pachytene="#E41A1C"),
                Sample=c(Juv="#8dd3c7",l26nt="#fccde5",inter="#b3de69",l30nt="#fdb462"))
o1=order(apply(mat_ppm1[human_pachy,],1,mean),decreasing=T)
o2=order(apply(mat_ppm1[human_hybrid,],1,mean),decreasing=T)
o3=order(apply(mat_ppm1[human_prepachy,],1,mean),decreasing=T)
mt=mat_ppm2[c(human_pachy[o1],human_hybrid[o2],human_prepachy[o3]),]
row.names(annoR)=colnames(mt)
row.names(annoC)=row.names(mt)
pheatmap(log10(mt+1),
         cluster_rows=F,cluster_cols=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         cellwidth=16,cellheight=4,border="white",lwd=0.3,fontsize_row=4,fontsize_col=16,
         filename="../Bear_Jul12.2018/Figure1.srna.heatmap.newSample.pdf",
         annotation_row=annoC,annotation_col=annoR,annotation_colors=ann_colors)
mt=df_rna[c(human_pachy,human_hybrid,human_prepachy),setdiff(rna_rn,rna_rn[c(1,2,5,6,9,8,14,13,15)])]
row.names(annoR)=colnames(mt)
pheatmap(log10(mt+1),
         cluster_rows=F,cluster_cols=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         cellwidth=16,cellheight=4,border="white",lwd=0.3,fontsize_row=4,fontsize_col=16,
         filename="../Bear_Jul12.2018/Figure1.rna.heatmap.newSample.pdf",
         annotation_row=annoC,annotation_col=annoR,annotation_colors=ann_colors)

# Figure 1c
annoC=data.frame(Type=factor(rep(c("Pachytene","Hybrid","Prepachytene"),c(89,10,83))))
annoR=data.frame(Sample=factor(rep(c("Juv","l26nt","inter","l30nt"),c(3,3,5,6))))
ann_colors=list(Type=c(Prepachytene=csl[2],Hybrid="#984ea3",Pachytene="#E41A1C"),
                Sample=c(Juv="#8dd3c7",l26nt="#fccde5",inter="#b3de69",l30nt="#fdb462"))
o1=hclust(dist(mat_ppm[human_pachy,]))$order
o2=hclust(dist(mat_ppm[human_hybrid,]))$order
o3=hclust(dist(mat_ppm[human_prepachy,]))$order
mt=cbind(apply(mat_ppm[c(human_pachy[o1],human_hybrid[o2],human_prepachy[o3]),1:2],1,mean),
         mat_ppm[c(human_pachy[o1],human_hybrid[o2],human_prepachy[o3]),-c(1,2)])
colnames(mt)[1]="N1643-5YO-ox"
row.names(annoR)=colnames(mt)
row.names(annoC)=row.names(mt)
pheatmap(log10(t(mt)+1),
         cluster_rows=F,cluster_cols=F,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),
         cellwidth=5,cellheight=12,border="white",lwd=0.3,fontsize_row=10,fontsize_col=4,
         filename="../Bear_Jul12.2018/Figure1.srna.horizontal.heatmap.pdf",
         annotation_row=annoR,annotation_col=annoC,annotation_colors=ann_colors)

# Aug01.2018: [average] lendis of all samples; rebuttal in May05, 2019----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
# average lendis for pachytene and prepachytene
aver.pachy.lendis=read.table("lendis/merged.pachytene.lendis",header=T,row.names=1,check.names=F)
aver.prepachy.lendis=read.table("lendis/merged.prepachy.lendis",header=T,row.names=1,check.names=F)
fun_t=function(l){
  return(sum(l*c(15:35))/sum(l))
}

pdf("../Bear_Jul12.2018/avergae_lendis_pachyAndPrepachy.pdf",width=6,height=3,useDingbats = F)
par(mfrow=c(1,2),mar=c(2,2,3,1),tcl=0.3,cex=5/6,bty="n",xpd=F)
plot(c(1:3,5:7,9:13,15:20),apply(aver.pachy.lendis[,srna_rn[-2]],2,fun_t),
     lwd=1.5,xaxt="n",xlab="",main="pachytene piRNAs",ylim=c(26,30))
axis(1,c(2,6,11,17.5),label=c("juv","groupIII","groupII","groupI"),lwd=0)
abline(v=c(4,8,14),lty=2)
plot(c(1:3,5:7,9:13,15:20),apply(aver.prepachy.lendis[,srna_rn[-2]],2,fun_t),
     lwd=1.5,xaxt="n",xlab="",main="prepachytene piRNAs",ylim=c(26,30))
axis(1,c(2,6,11,17.5),label=c("juv","groupIII","groupII","groupI"),lwd=0)
abline(v=c(4,8,14),lty=2)
dev.off()



pdf("../../../Ozata et al. 1st Paper/Rebutal Cell/Figure3A.pdf",width=12,height=6,useDingbats = F)
par(mfrow=c(2,4),mar=c(2,2,3,1),tcl=0.3,cex=5/6,bty="n",xpd=F)
fun_t=function(df,n,m,mm){
  barplot(apply(df[4:21,n],1,mean),space=0,border="white",col="black",names=18:35,
          ylab="RPM",main=m,ylim=c(0,mm))
}
fun_t(aver.pachy.lendis,srna_rn[1:4],"pachy juv",150000)
fun_t(aver.pachy.lendis,srna_rn[13:18],"pachy 30nt",150000)
fun_t(aver.pachy.lendis,srna_rn[8:12],"pachy inter",150000)
fun_t(aver.pachy.lendis,srna_rn[5:7],"pachy 26nt",150000)
fun_t(aver.prepachy.lendis,srna_rn[1:4],"prepachy juv",15000)
fun_t(aver.prepachy.lendis,srna_rn[13:18],"prepachy 30nt",15000)
fun_t(aver.prepachy.lendis,srna_rn[8:12],"prepachy inter",15000)
fun_t(aver.prepachy.lendis,srna_rn[5:7],"prepachy 26nt",15000)
dev.off()

# avergae lendis
aver.lendis=read.table("merged.average.lendis",header=F,row.names=1)
pdf("../Bear_Jul12.2018/average.lendis.pdf",width=5,height=4)
par(bty="n",cex=5/6,mar=c(9,4,1,1))
plot(NA,xlim=c(1,18),ylim=c(26.5,29),xlab="",ylab="average length of piRNAs",
     xaxt="n")
polygon(c(0,0,4.5,4.5),c(0,30,30,0),border="#8dd3c7",col="#8dd3c7")
polygon(c(4.5,4.5,7.5,7.5),c(0,30,30,0),border="#fccde5",col="#fccde5")
polygon(c(7.5,7.5,12.5,12.5),c(0,30,30,0),border="#b3de69",col="#b3de69")
polygon(c(12.5,12.5,20,20),c(0,30,30,0),border="#fdb462",col="#fdb462")
points(aver.lendis[srna_rn,1],pch=21,bg="white",col="black",lwd=2)
abline(h=c(27.35,28),col="black",lty=2)
text(c(2.5,6,10,15.5),c(29,29,29,29)-0.1,label=c("Juv","26nt","inter","30nt"))
axis(1,1:18,label=srna_rn,las=3,lwd=0)
dev.off()

# lendis with normalization
human_lendis=list()
m=0
for(i in srna_rn){
  human_lendis[[i]]=read.table(paste("lendis/",i,".lendis",sep=""),header=F,row.names=1)
  m=max(m,human_lendis[[i]][,1])
}
pdf("../Bear_Jul12.2018/allSample.lendis.pdf",width=4,height=3)
par(mar=c(3,4,4,1),cex=5/6)
n=1;tcs=c(rep("#8dd3c7",4),rep("#fccde5",3),rep("#b3de69",5),rep("#fdb462",6))
for(i in srna_rn){
  barplot(as.vector(human_lendis[[i]][18:35,1]),ylab="RPM",names=18:35,main=i,
          ylim=c(0,m),col=tcs[n])
  n=n+1
}
dev.off()

# prepachytene piRNA lendis with normalization
human_lendis=list()
m=0
for(i in srna_rn){
  human_lendis[[i]]=read.table(paste("lendis/",i,".prepachy.lendis",sep=""),header=F,row.names=1)
  m=max(m,human_lendis[[i]][,1])
}
pdf("../Bear_Jul12.2018/allSample.prepachy.lendis.pdf",width=4,height=3)
par(mar=c(3,4,4,1),cex=5/6)
n=1;tcs=c(rep("#8dd3c7",4),rep("#fccde5",3),rep("#b3de69",5),rep("#fdb462",6))
for(i in srna_rn){
  barplot(as.vector(human_lendis[[i]][4:21,1]),ylab="RPM",names=18:35,main=i,
          ylim=c(0,m),col=tcs[n])
  n=n+1
}
dev.off()


# Aug01.2018: examples for syntenic piRNA genes----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
human_c=read.table("syntenic/examples/human.c",header=F,row.names=NULL)
human_w=read.table("syntenic/examples/human.w",header=F,row.names=NULL)
mouse_c=read.table("syntenic/examples/mouse.c",header=F,row.names=NULL)
mouse_w=read.table("syntenic/examples/mouse.w",header=F,row.names=NULL)
rat_c=read.table("syntenic/examples/rat.c",header=F,row.names=NULL)
rat_w=read.table("syntenic/examples/rat.w",header=F,row.names=NULL)
rhesus_c=read.table("syntenic/examples/rhesus.c",header=F,row.names=NULL)
rhesus_w=read.table("syntenic/examples/rhesus.w",header=F,row.names=NULL)
marmoset_c=read.table("syntenic/examples/marmoset.c",header=F,row.names=NULL)
marmoset_w=read.table("syntenic/examples/marmoset.w",header=F,row.names=NULL)
human_bed=read.table("syntenic/examples/human.bed",header=F,row.names=NULL)
mouse_bed=read.table("syntenic/examples/mouse.bed",header=F,row.names=NULL)
rat_bed=read.table("syntenic/examples/rat.bed",header=F,row.names=NULL)
rhesus_bed=read.table("syntenic/examples/rhesus.bed",header=F,row.names=NULL)
marmoset_bed=read.table("syntenic/examples/marmoset.bed",header=F,row.names=NULL)
human_phast=read.table("syntenic/examples/human.phastCon.examples.smoothed.bdg",header=F,row.names=NULL)
rhesus_amyb=read.table("syntenic/examples/rhesus.amyb.examples.bdg",header=F,row.names=NULL)
rhesus_input=read.table("syntenic/examples/rhesus.input.examples.bdg",header=F,row.names=NULL)
mouse_amyb=read.table("syntenic/examples/mouse.amyb.examples.bdg",header=F,row.names=NULL)
mouse_input=read.table("syntenic/examples/mouse.input.examples.bdg",header=F,row.names=NULL)

# function
fun_plot=function(df_w,df_c,chrom,sp,df_bed){
  tc=df_c[which(df_c[,1]==chrom),]
  tw=df_w[which(df_w[,1]==chrom),]
  tb=df_bed[which(df_bed[,1]==chrom),]
  plot(NA,xlim=c(min(tb[,2]-10000),max(tb[,3])+10000),xlab="",ylab=paste(sp,chrom,sep=" "),ylim=c(min(tc[,4]),max(tw[,4])))
  for(i in 1:dim(tc)[1]){polygon(c(tc[i,2],tc[i,2],tc[i,3],tc[i,3]),c(0,tc[i,4],tc[i,4],0),col=csl[2],border=csl[2])}
  for(i in 1:dim(tw)[1]){polygon(c(tw[i,2],tw[i,2],tw[i,3],tw[i,3]),c(0,tw[i,4],tw[i,4],0),col=csl[1],border=csl[1])}
}
fun_plot1=function(df,chrom,sp,df_bed){
  tb=df_bed[which(df_bed[,1]==chrom),]
  tc=df[which(df[,1]==chrom & df[,2]>=min(tb[,2]-10000) & df[,3]<max(tb[,3])+10000),]
  plot(NA,xlim=c(min(tb[,2]-10000),max(tb[,3])+10000),xlab="",ylab=paste(sp,chrom,sep=" "),ylim=c(0,max(tc[,4])))
  for(i in 1:dim(tc)[1]){polygon(c(tc[i,2],tc[i,2],tc[i,3],tc[i,3]),c(0,tc[i,4],tc[i,4],0),col=csl[2],border=csl[2])}
}
# for human piRNA gene 15-q22-56093 and 15-q22-8218
par(mfrow=c(5,1),mar=c(2,4,1,1))
fun_plot(human_w,human_c,"chr15","human",human_bed)
fun_plot(rhesus_w,rhesus_c,"chr7","rhesus",rhesus_bed)
fun_plot(marmoset_w,marmoset_c,"chr10","marmoset",marmoset_bed)
fun_plot(mouse_w,mouse_c,"chr9","mouse",mouse_bed)
fun_plot(rat_w,rat_c,"chr8","rat",rat_bed)

### example1
pdf("../Bear_Jul12.2018/syntenic.example1.pdf",width=8,height=8*6/5)
par(tcl=0.3,bty="n",xpd=T,mfrow=c(6,1),mar=c(1,4,0.5,1))
tl=70000
# human
par(xpd=F)
fun_plot1(human_phast,"chr15","human",human_bed)
par(xpd=T)
x1=62488208;x2=62516382;x3=62516075;x4=62574361;chrom="chr15"
tc=human_c[which(human_c[,1]==chrom),]
tw=human_w[which(human_w[,1]==chrom),]
plot(NA,xlim=c(min(x3-tl),x3+tl),xlab="",ylab="human",ylim=c(min(tc[,4]),max(tw[,4])),xaxt="n")
for(i in 1:dim(tc)[1]){polygon(c(tc[i,2],tc[i,2],tc[i,3],tc[i,3]),c(0,tc[i,4],tc[i,4],0),col=csl[2],border=csl[2])}
for(i in 1:dim(tw)[1]){polygon(c(tw[i,2],tw[i,2],tw[i,3],tw[i,3]),c(0,tw[i,4],tw[i,4],0),col=csl[1],border=csl[1])}
polygon(c(x1,x1,x2,x2),
        c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
        col=csl[2])
polygon(c(x3,x3,x4,x4),
        c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
        col=csl[1])
text((x1+x2)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="15-q22-8218",col="white")
text((x3+x4)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="15-q22-56093",col="white")
text(x3-tl,max(tw[,4]),pos=4,
     label=paste(chrom,":",format(x1,big.mark=","),"-",format(x4,big.mark=","),sep=""))
# rhesus
x1=38758450;x2=38797823;x3=38798201;x4=38847100;chrom="chr7"
tc=rhesus_c[which(rhesus_c[,1]==chrom),]
tw=rhesus_w[which(rhesus_w[,1]==chrom),]
plot(NA,xlim=c(min(x3-tl),x3+tl),xlab="",ylab="rhesus",ylim=c(min(tc[,4]),max(tw[,4])),xaxt="n")
for(i in 1:dim(tc)[1]){polygon(c(tc[i,2],tc[i,2],tc[i,3],tc[i,3]),c(0,tc[i,4],tc[i,4],0),col=csl[2],border=csl[2])}
for(i in 1:dim(tw)[1]){polygon(c(tw[i,2],tw[i,2],tw[i,3],tw[i,3]),c(0,tw[i,4],tw[i,4],0),col=csl[1],border=csl[1])}
# polygon(c(x1,x1,x2,x2),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[2])
# polygon(c(x3,x3,x4,x4),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[1])
#text((x1+x2)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="15-q22-8218",col="white")
#text((x3+x4)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="15-q22-56093",col="white")
text(x3-tl,max(tw[,4]),pos=4,
     label=paste(chrom,":",format(x1,big.mark=","),"-",format(x4,big.mark=","),sep=""))
# marmoset
x1=4512401;x2=4570300;x3=4570601;x4=4639500;chrom="chr10"
tc=marmoset_c[which(marmoset_c[,1]==chrom),]
tw=marmoset_w[which(marmoset_w[,1]==chrom),]
plot(NA,xlim=c(min(x3-tl),x3+tl),xlab="",ylab="marmoset",ylim=c(min(tc[,4]),max(tw[,4])),xaxt="n")
for(i in 1:dim(tc)[1]){polygon(c(tc[i,2],tc[i,2],tc[i,3],tc[i,3]),c(0,tc[i,4],tc[i,4],0),col=csl[2],border=csl[2])}
for(i in 1:dim(tw)[1]){polygon(c(tw[i,2],tw[i,2],tw[i,3],tw[i,3]),c(0,tw[i,4],tw[i,4],0),col=csl[1],border=csl[1])}
# polygon(c(x1,x1,x2,x2),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[2])
# polygon(c(x3,x3,x4,x4),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[1])
#text((x1+x2)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="15-q22-8218",col="white")
#text((x3+x4)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="15-q22-56093",col="white")
text(x3-tl,max(tw[,4]),pos=4,
     label=paste(chrom,":",format(x1,big.mark=","),"-",format(x4,big.mark=","),sep=""))
# mouse
x1=67691251;x2=67733786;x3=67733944;x4=67760929;chrom="chr9"
tc=mouse_c[which(mouse_c[,1]==chrom),]
tw=mouse_w[which(mouse_w[,1]==chrom),]
plot(NA,xlim=c(min(x3+tl),x3-tl),xlab="",ylab="mouse",ylim=c(min(tc[,4]),max(tw[,4])),xaxt="n")
for(i in 1:dim(tc)[1]){polygon(c(tc[i,2],tc[i,2],tc[i,3],tc[i,3]),c(0,tc[i,4],tc[i,4],0),col=csl[2],border=csl[2])}
for(i in 1:dim(tw)[1]){polygon(c(tw[i,2],tw[i,2],tw[i,3],tw[i,3]),c(0,tw[i,4],tw[i,4],0),col=csl[1],border=csl[1])}
# polygon(c(x1,x1,x2,x2),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[2])
# polygon(c(x3,x3,x4,x4),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[1])
#text((x1+x2)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="15-q22-56093",col="white")
#text((x3+x4)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="15-q22-8218",col="white")
text(x3+tl,max(tw[,4]),pos=4,
     label=paste(chrom,":",format(x4,big.mark=","),"-",format(x1,big.mark=","),sep=""))
# rat
x1=73533101;x2=73572100;x3=73572401;x4=73599100;chrom="chr8"
tc=rat_c[which(rat_c[,1]==chrom),]
tw=rat_w[which(rat_w[,1]==chrom),]
plot(NA,xlim=c(min(x3+tl),x3-tl),xlab="",ylab="rat",ylim=c(min(tc[,4]),max(tw[,4])),xaxt="n")
for(i in 1:dim(tc)[1]){polygon(c(tc[i,2],tc[i,2],tc[i,3],tc[i,3]),c(0,tc[i,4],tc[i,4],0),col=csl[2],border=csl[2])}
for(i in 1:dim(tw)[1]){polygon(c(tw[i,2],tw[i,2],tw[i,3],tw[i,3]),c(0,tw[i,4],tw[i,4],0),col=csl[1],border=csl[1])}
# polygon(c(x1,x1,x2,x2),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[2])
# polygon(c(x3,x3,x4,x4),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[1])
#text((x1+x2)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="15-q22-56093",col="white")
#text((x3+x4)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="15-q22-8218",col="white")
text(x3+tl,max(tw[,4]),pos=4,
     label=paste(chrom,":",format(x4,big.mark=","),"-",format(x1,big.mark=","),sep=""))
dev.off()

### example2
pdf("../Bear_Jul12.2018/syntenic.example2.pdf",width=8,height=8*6/5)
par(tcl=0.3,bty="n",xpd=T,mfrow=c(5,1),mar=c(1,4,0.5,1))
tl=70000
# human
par(xpd=F)
fun_plot1(human_phast,"chr6","human",human_bed)
par(xpd=T)
x1=33835283;x2=33860229;x3=33860941;x4=33885986;chrom="chr6"
tc=human_c[which(human_c[,1]==chrom),]
tw=human_w[which(human_w[,1]==chrom),]
plot(NA,xlim=c(min(x3-tl),x3+tl),xlab="",ylab="human",ylim=c(min(tc[,4]),max(tw[,4])),xaxt="n")
for(i in 1:dim(tc)[1]){polygon(c(tc[i,2],tc[i,2],tc[i,3],tc[i,3]),c(0,tc[i,4],tc[i,4],0),col=csl[2],border=csl[2])}
for(i in 1:dim(tw)[1]){polygon(c(tw[i,2],tw[i,2],tw[i,3],tw[i,3]),c(0,tw[i,4],tw[i,4],0),col=csl[1],border=csl[1])}
polygon(c(x1,x1,x2,x2),
        c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
        col=csl[2])
polygon(c(x3,x3,x4,x4),
        c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
        col=csl[1])
text((x1+x2)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="6-p21-16923",col="white")
text((x3+x4)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="6-p21-43244",col="white")
text(x3-tl,max(tw[,4]),pos=4,
     label=paste(chrom,":",format(x1,big.mark=","),"-",format(x4,big.mark=","),sep=""))
# rhesus
x1=34726601;x2=34752100;x3=34752101;x4=34810600;chrom="chr4"
tc=rhesus_c[which(rhesus_c[,1]==chrom),]
tw=rhesus_w[which(rhesus_w[,1]==chrom),]
plot(NA,xlim=c(min(x3-tl),x3+tl),xlab="",ylab="rhesus",ylim=c(min(tc[,4]),max(tw[,4])),xaxt="n")
for(i in 1:dim(tc)[1]){polygon(c(tc[i,2],tc[i,2],tc[i,3],tc[i,3]),c(0,tc[i,4],tc[i,4],0),col=csl[2],border=csl[2])}
for(i in 1:dim(tw)[1]){polygon(c(tw[i,2],tw[i,2],tw[i,3],tw[i,3]),c(0,tw[i,4],tw[i,4],0),col=csl[1],border=csl[1])}
# polygon(c(x1,x1,x2,x2),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[2])
# polygon(c(x3,x3,x4,x4),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[1])
#text((x1+x2)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="6-p21-16923",col="white")
#text((x3+x4)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="6-p21-43244",col="white")
text(x3-tl,max(tw[,4]),pos=4,
     label=paste(chrom,":",format(x1,big.mark=","),"-",format(x4,big.mark=","),sep=""))
# marmoset
x1=34626374;x2=34719301;x3=34719501;x4=34882900;chrom="chr4"
tc=marmoset_c[which(marmoset_c[,1]==chrom),]
tw=marmoset_w[which(marmoset_w[,1]==chrom),]
plot(NA,xlim=c(min(x3-tl),x3+tl),xlab="",ylab="marmoset",ylim=c(min(tc[,4]),max(tw[,4])),xaxt="n")
for(i in 1:dim(tc)[1]){polygon(c(tc[i,2],tc[i,2],tc[i,3],tc[i,3]),c(0,tc[i,4],tc[i,4],0),col=csl[2],border=csl[2])}
for(i in 1:dim(tw)[1]){polygon(c(tw[i,2],tw[i,2],tw[i,3],tw[i,3]),c(0,tw[i,4],tw[i,4],0),col=csl[1],border=csl[1])}
# polygon(c(x1,x1,x2,x2),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[2])
# polygon(c(x3,x3,x4,x4),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[1])
#text((x1+x2)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="6-p21-16923",col="white")
#text((x3+x4)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="6-p21-43244",col="white")
text(x3-tl,max(tw[,4]),pos=4,
     label=paste(chrom,":",format(x1,big.mark=","),"-",format(x4,big.mark=","),sep=""))
# mouse
x1=27288275;x2=27325028;x3=27325196;x4=27367483;chrom="chr17"
tc=mouse_c[which(mouse_c[,1]==chrom),]
tw=mouse_w[which(mouse_w[,1]==chrom),]
plot(NA,xlim=c(min(x3-tl),x3+tl),xlab="",ylab="mouse",ylim=c(min(tc[,4]),max(tw[,4])),xaxt="n")
for(i in 1:dim(tc)[1]){polygon(c(tc[i,2],tc[i,2],tc[i,3],tc[i,3]),c(0,tc[i,4],tc[i,4],0),col=csl[2],border=csl[2])}
for(i in 1:dim(tw)[1]){polygon(c(tw[i,2],tw[i,2],tw[i,3],tw[i,3]),c(0,tw[i,4],tw[i,4],0),col=csl[1],border=csl[1])}
# polygon(c(x1,x1,x2,x2),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[2])
# polygon(c(x3,x3,x4,x4),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[1])
#text((x1+x2)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="6-p21-16923",col="white")
#text((x3+x4)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="6-p21-43244",col="white")
text(x3-tl,max(tw[,4]),pos=4,
     label=paste(chrom,":",format(x1,big.mark=","),"-",format(x4,big.mark=","),sep=""))
# rat
x1=6526001;x2=6652300;x3=6652501;x4=6692800;chrom="chr20"
tc=rat_c[which(rat_c[,1]==chrom),]
tw=rat_w[which(rat_w[,1]==chrom),]
plot(NA,xlim=c(min(x3-tl),x3+tl),xlab="",ylab="rat",ylim=c(min(tc[,4]),max(tw[,4])),xaxt="n")
for(i in 1:dim(tc)[1]){polygon(c(tc[i,2],tc[i,2],tc[i,3],tc[i,3]),c(0,tc[i,4],tc[i,4],0),col=csl[2],border=csl[2])}
for(i in 1:dim(tw)[1]){polygon(c(tw[i,2],tw[i,2],tw[i,3],tw[i,3]),c(0,tw[i,4],tw[i,4],0),col=csl[1],border=csl[1])}
# polygon(c(x1,x1,x2,x2),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[2])
# polygon(c(x3,x3,x4,x4),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[1])
#text((x1+x2)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="6-p21-16923",col="white")
#text((x3+x4)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="6-p21-43244",col="white")
text(x3-tl,max(tw[,4]),pos=4,
     label=paste(chrom,":",format(x1,big.mark=","),"-",format(x4,big.mark=","),sep=""))
dev.off()

### example3
pdf("../Bear_Jul12.2018/syntenic.example3.pdf",width=8,height=8*6/5)
par(tcl=0.3,bty="n",xpd=T,mfrow=c(6,1),mar=c(1,4,0.5,1))
tl=70000
# human
par(xpd=F)
#fun_plot1(human_phast,"chr12","human",human_bed)
par(xpd=T)
x1=3546566;x2=3570665;x3=3570736;x4=3600637;chrom="chr12"
tc=human_c[which(human_c[,1]==chrom),]
tw=human_w[which(human_w[,1]==chrom),]
plot(NA,xlim=c(min(x3-tl),x3+tl),xlab="",ylab="human",ylim=c(min(tc[,4]),max(tw[,4])),xaxt="n")
for(i in 1:dim(tc)[1]){polygon(c(tc[i,2],tc[i,2],tc[i,3],tc[i,3]),c(0,tc[i,4],tc[i,4],0),col=csl[2],border=csl[2])}
for(i in 1:dim(tw)[1]){polygon(c(tw[i,2],tw[i,2],tw[i,3],tw[i,3]),c(0,tw[i,4],tw[i,4],0),col=csl[1],border=csl[1])}
polygon(c(x1,x1,x2,x2),
        c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
        col=csl[2])
polygon(c(x3,x3,x4,x4),
        c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
        col=csl[1])
text((x1+x2)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="15-q22-8218",col="white")
text((x3+x4)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="15-q22-56093",col="white")
text(x3-tl,max(tw[,4]),pos=4,
     label=paste(chrom,":",format(x1,big.mark=","),"-",format(x4,big.mark=","),sep=""))
# rhesus
x1=3562601;x2=3605801;x3=3606101;x4=3676600;chrom="chr11"
tc=rhesus_c[which(rhesus_c[,1]==chrom),]
tw=rhesus_w[which(rhesus_w[,1]==chrom),]
plot(NA,xlim=c(min(x3-tl),x3+tl),xlab="",ylab="rhesus",ylim=c(min(tc[,4]),max(tw[,4])),xaxt="n")
for(i in 1:dim(tc)[1]){polygon(c(tc[i,2],tc[i,2],tc[i,3],tc[i,3]),c(0,tc[i,4],tc[i,4],0),col=csl[2],border=csl[2])}
for(i in 1:dim(tw)[1]){polygon(c(tw[i,2],tw[i,2],tw[i,3],tw[i,3]),c(0,tw[i,4],tw[i,4],0),col=csl[1],border=csl[1])}
# polygon(c(x1,x1,x2,x2),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[2])
# polygon(c(x3,x3,x4,x4),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[1])
#text((x1+x2)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="15-q22-8218",col="white")
#text((x3+x4)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="15-q22-56093",col="white")
text(x3-tl,max(tw[,4]),pos=4,
     label=paste(chrom,":",format(x1,big.mark=","),"-",format(x4,big.mark=","),sep=""))
# marmoset
x1=14996301;x2=15096300;x3=15096601;x4=15192600;chrom="chr9"
tc=marmoset_c[which(marmoset_c[,1]==chrom),]
tw=marmoset_w[which(marmoset_w[,1]==chrom),]
plot(NA,xlim=c(min(x3-tl),x3+tl),xlab="",ylab="marmoset",ylim=c(min(tc[,4]),max(tw[,4])),xaxt="n")
for(i in 1:dim(tc)[1]){polygon(c(tc[i,2],tc[i,2],tc[i,3],tc[i,3]),c(0,tc[i,4],tc[i,4],0),col=csl[2],border=csl[2])}
for(i in 1:dim(tw)[1]){polygon(c(tw[i,2],tw[i,2],tw[i,3],tw[i,3]),c(0,tw[i,4],tw[i,4],0),col=csl[1],border=csl[1])}
# polygon(c(x1,x1,x2,x2),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[2])
# polygon(c(x3,x3,x4,x4),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[1])
#text((x1+x2)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="15-q22-8218",col="white")
#text((x3+x4)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="15-q22-56093",col="white")
text(x3-tl,max(tw[,4]),pos=4,
     label=paste(chrom,":",format(x1,big.mark=","),"-",format(x4,big.mark=","),sep=""))
# mouse
x1=127776075;x2=127796372;x3=127796430;x4=127841890;chrom="chr6"
tc=mouse_c[which(mouse_c[,1]==chrom),]
tw=mouse_w[which(mouse_w[,1]==chrom),]
plot(NA,xlim=c(min(x3-tl),x3+tl),xlab="",ylab="mouse",ylim=c(min(tc[,4]),max(tw[,4])),xaxt="n")
for(i in 1:dim(tc)[1]){polygon(c(tc[i,2],tc[i,2],tc[i,3],tc[i,3]),c(0,tc[i,4],tc[i,4],0),col=csl[2],border=csl[2])}
for(i in 1:dim(tw)[1]){polygon(c(tw[i,2],tw[i,2],tw[i,3],tw[i,3]),c(0,tw[i,4],tw[i,4],0),col=csl[1],border=csl[1])}
# polygon(c(x1,x1,x2,x2),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[2])
# polygon(c(x3,x3,x4,x4),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[1])
#text((x1+x2)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="6-p21-16923",col="white")
#text((x3+x4)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="6-p21-43244",col="white")
text(x3-tl,max(tw[,4]),pos=4,
     label=paste(chrom,":",format(x1,big.mark=","),"-",format(x4,big.mark=","),sep=""))
# rat
x1=160341301;x2=160362300;x3=160362301;x4=160417227;chrom="chr4"
tc=rat_c[which(rat_c[,1]==chrom),]
tw=rat_w[which(rat_w[,1]==chrom),]
plot(NA,xlim=c(min(x3-tl),x3+tl),xlab="",ylab="rat",ylim=c(min(tc[,4]),max(tw[,4])),xaxt="n")
for(i in 1:dim(tc)[1]){polygon(c(tc[i,2],tc[i,2],tc[i,3],tc[i,3]),c(0,tc[i,4],tc[i,4],0),col=csl[2],border=csl[2])}
for(i in 1:dim(tw)[1]){polygon(c(tw[i,2],tw[i,2],tw[i,3],tw[i,3]),c(0,tw[i,4],tw[i,4],0),col=csl[1],border=csl[1])}
# polygon(c(x1,x1,x2,x2),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[2])
# polygon(c(x3,x3,x4,x4),
#         c(min(tc[,4]),min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.05,min(tc[,4])),
#         col=csl[1])
#text((x1+x2)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="6-p21-16923",col="white")
#text((x3+x4)/2,min(tc[,4])-(max(tw[,4])-min(tc[,4]))*0.025,label="6-p21-43244",col="white")
text(x3-tl,max(tw[,4]),pos=4,
     label=paste(chrom,":",format(x1,big.mark=","),"-",format(x4,big.mark=","),sep=""))
dev.off()

# Aug02.2018: piechart of rmsk content----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
rmsk_mat=read.table("piG.rmsk.mat",header=T,row.names=1)
fun_pie=function(df){
  p1=sum(df[,2]);p2=sum(df[,3]);p3=sum(df[,4]);p4=sum(df[,5]);
  p5=sum(df[,6]);p6=sum(df[,7]);p7=sum(df[,8]);p8=sum(df[,9])
  to=sum(df[,1])-p1-p2-p3-p4-p5-p6-p7-p8;tt=sum(df[,1])
  pie(c(p1,p2,p3,p4,p5,p6,p7,p8,to),
      density=c(-1,20,-1,20,-1,20,-1,20,-1),
      col=c(csl[c(1,1,2,2,3,3,4,4)],"grey"),
      label=paste(round(c(p1,p2,p3,p4,p5,p6,p7,p8,to)/tt*100,1),"%"))
}
pdf("../Bear_Jul12.2018/human_rmsk_pie.pdf",width=6,height=8)
par(mfrow=c(3,2),mar=c(3,3,3,3),xpd=T)
plot.new()
legend("left",fill=csl[1:4],legend=c("LINE","SINE","LTR","DNA"),cex=1.3)
legend("right",density=c(20,-1),legend=c("anti-sense","sense"),cex=1.3)
fun_pie(rmsk_mat[human_pachy,])
text(-1.4,-0.3,label="pachytene",srt=90,pos=4,font=2,cex=1.5)
fun_pie(rmsk_mat[human_pi,])
text(-1.4,-0.3,label="piRNA genes",srt=90,pos=4,font=2,cex=1.5)
fun_pie(rmsk_mat[human_prepachy,])
text(-1.4,-0.3,label="prepachytene",srt=90,pos=4,font=2,cex=1.5)
plot.new()
fun_pie(rmsk_mat[human_hybrid,])
text(-1.4,-0.3,label="hybrid",srt=90,pos=4,font=2,cex=1.5)
dev.off()


# Aug05.2018: differential expression analysis of 26nt and 30nt RNAseq----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
library(DESeq2)
df_rna1=read.table("merged.rnaseq.gene+pi.readCount",header=T,row.names=1,check.names=F)
df_rna1=df_rna1[,rna_rn]
fun_deseq=function(cn1,cn2){
  cd=data.frame(condition=factor(rep(c("untreated","treated"),c(length(cn1),length(cn2)))),
                row.names=colnames(df_rna1)[c(cn1,cn2)])
  dds=DESeqDataSetFromMatrix(countData=df_rna1[,c(cn1,cn2)],cd,design=~condition)
  keep=rowSums(counts(dds))>=10
  dds=dds[keep,]
  dds=DESeq(dds)
  res=results(dds)
  res=res[order(res$padj),]
  print(sizeFactors(dds))
  return(res)
}
der=fun_deseq(c(5,6,7,9),c(13,14,15,16,17,18))
t1=der[which(der$padj<0.001 & der$log2FoldChange<(-1)),]
t2=der[which(der$padj<0.001 & der$log2FoldChange>1),]
write.table(t1,"../../../Ozata et al. 1st Paper/Rebutal Cell/30vs26.upregulated.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(t2,"../../../Ozata et al. 1st Paper/Rebutal Cell/30vs26.downregulated.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(t1[intersect(row.names(t1),human_pc),],"../../../Ozata et al. 1st Paper/Rebutal Cell/30vs26.upregulated,mRNA.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(t2[intersect(row.names(t2),human_pc),],"../../../Ozata et al. 1st Paper/Rebutal Cell/30vs26.downregulated.mRNA.txt",row.names=T,col.names=T,sep="\t",quote=F)

# plot scatterplot
df_rna1=read.table("merged.rnaseq.gene+pi.rpkm",header=T,row.names=1,check.names=F)
df_rna1=df_rna1[,rna_rn]
pdf("../../../Ozata et al. 1st Paper/Rebutal Cell/Figure4B.pdf",width=5,height=5,useDingbats=F)
par(mar=c(4,4,2,2),tcl=0.3,cex=5/6,bty="n")
plot(log10(apply(df_rna1[human_pc,c(5,6,7)],1,mean)+1),
     log10(apply(df_rna1[human_pc,c(13,14,15,16,17,18)],1,mean)+1),
     bg="grey",col="white",pch=21,cex=0.6,
     xlab="mean RPKM in 26nt samples",ylab="mean RPKM in 30nt samples",
     xaxt="n",yaxt="n")
axis(1,0:6,label=expression(10^0,10^1,10^2,10^3,10^4,10^5,10^6))
axis(2,0:6,label=expression(10^0,10^1,10^2,10^3,10^4,10^5,10^6))
points(log10(apply(df_rna1[human_prepachy,c(5,6,7,9)],1,mean)+1),
       log10(apply(df_rna1[human_prepachy,c(13,14,15,16,17,18)],1,mean)+1),
       bg=csl[3],col="white",pch=21)
points(log10(apply(df_rna1[human_hybrid,c(5,6,7,9)],1,mean)+1),
       log10(apply(df_rna1[human_hybrid,c(13,14,15,16,17,18)],1,mean)+1),
       bg=csl[4],col="white",pch=21)
points(log10(apply(df_rna1[human_pachy,c(5,6,7,9)],1,mean)+1),
       log10(apply(df_rna1[human_pachy,c(13,14,15,16,17,18)],1,mean)+1),
       bg=csl[1],col="white",pch=21)
points(log10(apply(df_rna1[pl,c(5,6,7,9)],1,mean)+1),
       log10(apply(df_rna1[pl,c(13,14,15,16,17,18)],1,mean)+1),
       bg=csl[2],col="white",pch=21)
text(log10(apply(df_rna1[c("MYBL1","PIWIL1"),c(5,6,7,9)],1,mean)+1),
     log10(apply(df_rna1[c("MYBL1","PIWIL1"),c(13,14,15,16,17,18)],1,mean)+1),
     col="black",font=2,pos=3,label=c("AMYB","PIWIL1"))
abline(0,1)
plot(log10(apply(df_rna1[human_pc,c(8,10,11,12)],1,mean)+1),
     log10(apply(df_rna1[human_pc,c(13,14,15,16,17,18)],1,mean)+1),
     bg="grey",col="white",pch=21,cex=0.6,
     xlab="mean RPKM in intermediate samples",ylab="mean RPKM in 30nt samples",
     xaxt="n",yaxt="n")
axis(1,0:6,label=expression(10^0,10^1,10^2,10^3,10^4,10^5,10^6))
axis(2,0:6,label=expression(10^0,10^1,10^2,10^3,10^4,10^5,10^6))
points(log10(apply(df_rna1[human_prepachy,c(8,10,11,12)],1,mean)+1),
       log10(apply(df_rna1[human_prepachy,c(13,14,15,16,17,18)],1,mean)+1),
       bg=csl[3],col="white",pch=21)
points(log10(apply(df_rna1[human_hybrid,c(8,10,11,12)],1,mean)+1),
       log10(apply(df_rna1[human_hybrid,c(13,14,15,16,17,18)],1,mean)+1),
       bg=csl[4],col="white",pch=21)
points(log10(apply(df_rna1[human_pachy,c(8,10,11,12)],1,mean)+1),
       log10(apply(df_rna1[human_pachy,c(13,14,15,16,17,18)],1,mean)+1),
       bg=csl[1],col="white",pch=21)
points(log10(apply(df_rna1[pl,c(8,10,11,12)],1,mean)+1),
       log10(apply(df_rna1[pl,c(13,14,15,16,17,18)],1,mean)+1),
       bg=csl[2],col="white",pch=21)
text(log10(apply(df_rna1[c("MYBL1","PIWIL1"),c(8,10,11,12)],1,mean)+1),
     log10(apply(df_rna1[c("MYBL1","PIWIL1"),c(13,14,15,16,17,18)],1,mean)+1),
     label=c("AMYB","PIWIL1"),
     col="black",font=2,pos=3)
abline(0,1)
dev.off()

# Aug15.2018: rnaseq vs srnaseq scatterplot----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
df_rna=read.table("merged.rnaseq.gene+pi.rpkm",header=T,row.names=1,check.names=F)
df_rna=df_rna[,rna_rn]
df_srna=read.table("merged.rpm",header=T,row.names=1,check.names=F)
df_srna=df_srna[,srna_rn]
human_gene_len=read.table("hg19.gene+pi.exon.length",header=F,row.names=1)
df_srna=df_srna/human_gene_len[row.names(df_srna),]*1000

fun_temp=function(i1,i2){
  x1=log10(df_rna[c(human_pachy,human_prepachy,human_hybrid),i1]+1)
  x2=log10(df_srna[c(human_pachy,human_prepachy,human_hybrid),i2]+1)
  xm=2.6;ym=5
  plot(x1,x2,bg=c(rep(csl[1],89),rep(csl[3],83),rep(csl[4],10)),col="white",pch=21,
       xlab=paste(i1," RPKM",sep=""),ylab=paste(i2," RPKM",sep=""),
       xaxt="n",yaxt="n",xlim=c(0,xm),ylim=c(0,ym))
  axis(1,0:6,label=expression(10^0,10^1,10^2,10^3,10^4,10^5,10^6))
  axis(2,0:6,label=expression(10^0,10^1,10^2,10^3,10^4,10^5,10^6))
  c1=round(cor(df_rna[human_pachy,i1],df_srna[human_pachy,i2],method="spearman"),2)
  c2=round(cor(df_rna[human_prepachy,i1],df_srna[human_prepachy,i2],method="spearman"),2)
  text(xm*19/20,0,pos=3,label=c1,col=csl[1])
  text(xm*19/20,ym/20,pos=3,label=c2,col=csl[3])
}
pdf("../Bear_Jul12.2018/scatter.rna_vs_srna.pdf",width=4,height=4,useDingbats=F)
par(mar=c(4,4,1,1),tcl=0.3,bty="n",cex=5/6)
for(i in 1:18){fun_temp(rna_rn[i],srna_rn[i])}
dev.off()

# Aug15.2018: GO analysis
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
go_up=read.table("../../../Ozata et al. 1st Paper/Rebutal Cell/GO.upregulated.txt",header=T,row.names=1,sep="\t",stringsAsFactors=F)
go_down=read.table("../../../Ozata et al. 1st Paper/Rebutal Cell/GO.downregulated.txt",header=T,row.names=1,sep="\t",stringsAsFactors=F)
ti=which(go_up[,4]=="+");ti=ti[length(ti):1]
pdf("../../../Ozata et al. 1st Paper/Rebutal Cell/Figure4C.pdf",width=8,height=4)
par(mar=c(4,1,4,1),tcl=0.3,cex=5/6,bty="n",mfrow=c(1,2))
barplot(go_up[ti,5],horiz=T,xlim=c(0,8),space=0.2,border="white",xaxt="n",
        xlab="fold enrichment of GO term",main="up-regulated genes")
axis(1,0:3,label=0:3)
text(go_up[ti,5],c(1:length(ti))*1.2-0.6,
     label=paste(row.names(go_up)[ti]," q=",go_up[ti,7],sep=""),pos=4,cex=0.7)
ti=which(go_down[,4]=="+");ti=ti[length(ti):1]
par(mar=c(4,1,4,1),tcl=0.3,cex=5/6,bty="n")
barplot(as.numeric(go_down[ti,5]),horiz=T,xlim=c(0,8),space=0.2,border="white",xaxt="n",
        xlab="fold enrichment of GO term",main="down-regulated genes")
axis(1,0:3,label=0:3)
text(as.numeric(go_down[ti,5]),c(1:length(ti))*1.2-0.6,
     label=paste(row.names(go_down)[ti]," q=",go_down[ti,7],sep=""),pos=4,cex=0.7)
dev.off()



pdf("../Bear_Jul12.2018/GO.barplot.pdf",width=8,height=4)
par(mar=c(4,1,4,1),tcl=0.3,cex=5/6,bty="n",mfrow=c(1,2))
barplot(go_up[ti,5],horiz=T,xlim=c(0,8),space=0.2,border="white",xaxt="n",
        xlab="fold enrichment of GO term",main="up-regulated genes")
axis(1,0:3,label=0:3)
text(go_up[ti,5],c(1:length(ti))*1.2-0.6,
     label=paste(row.names(go_up)[ti]," q=",go_up[ti,7],sep=""),pos=4,cex=0.7)
ti=which(go_down[,4]=="+");ti=ti[length(ti):1]
par(mar=c(4,1,4,1),tcl=0.3,cex=5/6,bty="n")
barplot(as.numeric(go_down[ti,5]),horiz=T,xlim=c(0,8),space=0.2,border="white",xaxt="n",
        xlab="fold enrichment of GO term",main="down-regulated genes")
axis(1,0:3,label=0:3)
text(as.numeric(go_down[ti,5]),c(1:length(ti))*1.2-0.6,
     label=paste(row.names(go_down)[ti]," q=",go_down[ti,7],sep=""),pos=4,cex=0.7)
dev.off()



# Aug20.2018: boxplot for rnaseq: gene vs prepachytene genes----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
df_rna=read.table("merged.rnaseq.gene+pi.rpkm",header=T,row.names=1,check.names=F)
df_rna=df_rna[,rna_rn]
pdf("../Bear_Jul12.2018/boxplot.gene_vs_prepachy.rpkm.pdf",width=8,height=6,useDingbats=F)
boxplot(log10(df_rna[as.vector(human_pc),]+0.1),at=1:18,xlab="",ylab="log10(RPKM+0.1)",
        boxwex=0.3,outline=F,las=3)
boxplot(log10(df_rna[human_prepachy,]+0.1),at=1:18+0.3,add=T,xlab="",ylab="",xaxt="n",
        yaxt="n",boxwex=0.3,col=csl[3],pch=20,cex=0.3)
dev.off()

pdf("../Bear_Jul12.2018/boxplot2.gene_vs_prepachy.rpkm.pdf",width=6,height=5,useDingbats=F)
par(mar=c(4,4,2,2),tcl=0.3,bty="n")
o=order(apply(log10(df_rna[as.vector(human_pc),]+0.1),1,mean),decreasing=T)
plot(apply(log10(df_rna[as.vector(human_pc[o]),]+0.1),1,mean),pch=20,xlab="rank",ylab="log10(RPKM+0.1)")
boxplot(apply(log10(df_rna[human_prepachy,]+0.1),1,mean),at=10000,boxwex=1000,add=T,
        xlab="",ylab="",xaxt="n",yaxt="n",pch=20,cex=0.3,col=csl[3])
abline(h=c(0,1),
       col=csl[3],lty=2)
dev.off()

t1=human_prepachy[which(apply(df_rna[human_prepachy,],1,mean)<10)]
t2=human_prepachy[which(apply(df_rna[human_prepachy,],1,mean)>=10)]
cor(df_rna[t1,"D1643_5YO_Rep1"],df_srna[t1,"N1643-5YO-Rep1-Ox"],method="spearman")
cor(df_rna[t2,"D1643_5YO_Rep1"],df_srna[t2,"N1643-5YO-Rep1-Ox"],method="spearman")
cor(df_rna[t1,"D8150NT"],df_srna[t1,"N8150-Ox"],method="spearman")
cor(df_rna[t2,"D8150NT"],df_srna[t2,"N8150-Ox"],method="spearman")

# abline(h=c(summary(apply(log10(df_rna[human_prepachy,]+0.1),1,mean))[[2]],
#            summary(apply(log10(df_rna[human_prepachy,]+0.1),1,mean))[[5]]),
#        col=csl[3],lty=2)



o=order(apply(log10(df_rna[as.vector(human_pc),]+0.1),1,mean),decreasing=T)

# Aug21.2018: Xin's MC data analysis----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
# lendis
amyb_het=read.table("Xin_MC/mouse_amybHet_17.5dpp.lendis",header=F,row.names=1)
amyb_mut=read.table("Xin_MC/mouse_amybMut_17.5dpp.lendis",header=F,row.names=1)
miwi_het=read.table("Xin_MC/mouse_miwiHet_17.5dpp.lendis",header=F,row.names=1)
miwi_mut=read.table("Xin_MC/mouse_miwiMut_17.5dpp.lendis",header=F,row.names=1)
wt=read.table("Xin_MC/mouse_wt_17.5dpp.lendis",header=F,row.names=1)
tdf=cbind(wt/6.42166,amyb_het/62.9994,amyb_mut/2.5799,miwi_het/14.0531,miwi_mut/14.8815)
colnames(tdf)=c("wt","amyb_het","amyb_mut","miwi_het","miwi_mut")
tm=max(apply(tdf,1,max))
f=c()
pdf("../Bear_Jul12.2018/Xin.lendis.pdf",width=9,height=5,useDingbats=F)
par(mar=c(2,4,3,1),bty="n",cex=5/6,tcl=0.3,mfrow=c(2,3))
for(i in 1:5){
  barplot(tdf[,i],names=15:35,ylim=c(0,tm),col="black",main=colnames(tdf)[i])
}
dev.off()
# heatmap
df_xin_srna=read.table("Xin_MC/Xin.pirna.rpm",header=T,row.names=1)
df_xin_rna=read.table("Xin_MC/Xin_MC.rpkm",header=T,row.names=1)
o1=order(apply(df_xin_srna[Pachy,],1,sum),decreasing=T)
o2=order(apply(df_xin_srna[Hybrid,],1,sum),decreasing=T)
o3=order(apply(df_xin_srna[Prepachy,],1,sum),decreasing=T)
pdf("../Bear_Jul12.2018/Xin.heatmap.pdf",width=5,height=10,useDingbats=F)
annoR=data.frame(Type=factor(rep(c("Pachytene","Prepachytene","Hybrid"),c(100,84,30))))
row.names(annoR)=c(Pachy,Prepachy,Hybrid)
ann_colors=list(Type=c(Pachytene="#E41A1C",Hybrid="#984ea3",Prepachytene="#4daf4a"),
                position=c(beforeTSS="grey",afterTSS="navy"),Location=c(Intergenic="black",Genic="white"))
pheatmap(log10(df_xin_srna[c(Pachy[o1],Hybrid[o2],Prepachy[o3]),c(5,1,2,3,4)]+1),
         cluster_cols=F,cluster_rows=F,cellwidth=16,cellheight=2,
         fontsize_row=2,main="piRNA RPM",annotation_row=annoR,annotation_colors=ann_colors,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),lwd=0.3)
pheatmap(log10(df_xin_rna[c(Pachy[o1],Hybrid[o2],Prepachy[o3]),c(5,1,2,3,4)]+1),
         cluster_cols=F,cluster_rows=F,cellwidth=16,cellheight=2,
         fontsize_row=2,main="piRNA RPM",annotation_row=annoR,annotation_colors=ann_colors,
         col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100),lwd=0.3)
dev.off()



# Aug22.2018: 26-27nt and 30nt piRNA signal analysis----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
df_srna_26_27=read.table("merged.26-27nt.rpm",header=T,row.names=1,check.names=F)
df_srna_30=read.table("merged.30nt.rpm",header=T,row.names=1,check.names=F)
df_srna_26_27=df_srna_26_27[,srna_rn]
df_srna_30=df_srna_30[,srna_rn]
tcs=c(rep("#8dd3c7",4),rep("#fccde5",3),rep("#b3de69",5),rep("#fdb462",6))
pdf("../Bear_Jul12.2018/boxplot.26_27.vs.30.pdf",width=6,height=4,useDingbats=F)
par(mfrow=c(1,2),mar=c(4,4,4,1),tcl=0.3,bty="n",cex=5/6)
boxplot(log10(apply(df_srna_26_27[human_pachy,1:4],1,mean)+1),
        log10(apply(df_srna_26_27[human_pachy,5:7],1,mean)+1),
        log10(apply(df_srna_26_27[human_pachy,8:12],1,mean)+1),
        log10(apply(df_srna_26_27[human_pachy,13:18],1,mean)+1),col=tcs[c(1,5,8,13)],
        pch=20,cex=0.3,names=c("Juv","26nt","inter","30nt"),ylab="RPM",
        main="26-27nt piRNAs",yaxt="n")
axis(2,0:5,label=expression(10^0,10^1,10^2,10^3,10^4,10^5))
boxplot(log10(apply(df_srna_30[human_pachy,1:4],1,mean)+1),
        log10(apply(df_srna_30[human_pachy,5:7],1,mean)+1),
        log10(apply(df_srna_30[human_pachy,8:12],1,mean)+1),
        log10(apply(df_srna_30[human_pachy,13:18],1,mean)+1),col=tcs[c(1,5,8,13)],
        pch=20,cex=0.3,names=c("Juv","26nt","inter","30nt"),ylab="RPM",
        main="30nt piRNAs",yaxt="n")
axis(2,0:5,label=expression(10^0,10^1,10^2,10^3,10^4,10^5))
boxplot(log10(apply(df_srna_26_27[human_prepachy,1:4],1,mean)+1),
        log10(apply(df_srna_26_27[human_prepachy,5:7],1,mean)+1),
        log10(apply(df_srna_26_27[human_prepachy,8:12],1,mean)+1),
        log10(apply(df_srna_26_27[human_prepachy,13:18],1,mean)+1),col=tcs[c(1,5,8,13)],
        pch=20,cex=0.3,names=c("Juv","26nt","inter","30nt"),ylab="RPM",
        main="26-27nt piRNAs",yaxt="n")
axis(2,0:5,label=expression(10^0,10^1,10^2,10^3,10^4,10^5))
boxplot(log10(apply(df_srna_30[human_prepachy,1:4],1,mean)+1),
        log10(apply(df_srna_30[human_prepachy,5:7],1,mean)+1),
        log10(apply(df_srna_30[human_prepachy,8:12],1,mean)+1),
        log10(apply(df_srna_30[human_prepachy,13:18],1,mean)+1),col=tcs[c(1,5,8,13)],
        pch=20,cex=0.3,names=c("Juv","26nt","inter","30nt"),ylab="RPM",
        main="30nt piRNAs",yaxt="n")
axis(2,0:5,label=expression(10^0,10^1,10^2,10^3,10^4,10^5))
dev.off()

# Aug30.2018: ildar's srnaseq data----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
ildar.srna=read.table("mouse.ildar.merged.rpm",header=T,row.names=1,check.names=F)
ildar.srna1=read.table("mouse.ildar.merged.26-27nt.rpm",header=T,row.names=1,check.names=F)
ildar.srna2=read.table("mouse.ildar.merged.30nt.rpm",header=T,row.names=1,check.names=F)
ildar.lendis=list()
for(i in colnames(ildar.srna)){
  ildar.lendis[[i]]=read.table(paste("lendis/",i,".lendis",sep=""),header=T,row.names=1)
}
pdf("../Bear_Jul12.2018/ildar.lendis.pdf",width=7,height=5,useDingbats=F)
par(mfrow=c(2,3),mar=c(2,2,3,1),tcl=0.3,cex=5/6,bty="n",cex.main=0.6)
for(i in colnames(ildar.srna)){
  barplot(ildar.lendis[[i]][4:21,1],main=i,col="black",names=18:35)
}
dev.off()

fun_t=function(td,fn){
  annoR=data.frame(Type=factor(rep(c("Pachytene","Prepachytene","Hybrid"),c(100,84,30))))
  row.names(annoR)=c(Pachy,Prepachy,Hybrid)
  ann_colors=list(Type=c(Pachytene="#E41A1C",Hybrid="#984ea3",Prepachytene="#4daf4a"))
  o1=order(apply(td[Pachy,],1,mean),decreasing=T)
  o2=order(apply(td[Hybrid,],1,mean),decreasing=T)
  o3=order(apply(td[Prepachy,],1,mean),decreasing=T)
  t=td[c(Pachy[o1],Hybrid[o2],Prepachy[o3]),]
  t1=cbind(rowMeans(t[,c(4,14)]),rowMeans(t[,c(6,8)]),t[,2],t[,1],
           rowMeans(t[,c(19,20)]),rowMeans(t[,c(21,22)]))
  colnames(t1)=c("WT.Prepachy","WT.4N","AMYB--.Prepachy","AMYB--.4N","MIWI.HET.4N","MIWI-KO.4N")
  pheatmap(log10(t1+1),
           annotation_row=annoR,,annotation_colors=ann_colors,cluster_cols=F,cluster_rows=F,
           cellwidth=16,cellheight=2,fontsize_row=2,
           filename=fn,
           col=colorRampPalette(c("white",brewer.pal(9,"Reds")))(100))
  return(t1)
}
t1=fun_t(ildar.srna,"../Bear_Jul12.2018/ildar.heatmap.pdf")
t1=fun_t(ildar.srna1,"../Bear_Jul12.2018/ildar.heatmap.26-27nt.pdf")
t1=fun_t(ildar.srna2,"../Bear_Jul12.2018/ildar.heatmap.30nt.pdf")

# Sep04.2018: AMYB peak distance for rhesus----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
amyb_distance_rhesus=read.table("../Bear_Jul12.2018/amyb_gene_list/A-MYB_ChIP_positive_genes/Rhesus/rhesus.amyb.distance.txt",header=F,row.names=1)
rhesus_igp=row.names(amyb_distance_rhesus)[grep("_IG_",row.names(amyb_distance_rhesus))]
rhesus_pcp=row.names(amyb_distance_rhesus)[grep("_PC_",row.names(amyb_distance_rhesus))]
rhesus_pc=setdiff(row.names(amyb_distance_rhesus),c(rhesus_igp,rhesus_pcp))
t=c(rhesus_igp,rhesus_pcp)
rhesus_pi_amyb=t[amyb_distance_rhesus[t,1]<=500]
rhesus_pi_noamyb=t[amyb_distance_rhesus[t,1]>500]
pdf("../Bear_Jul12.2018/peak_distance_rhesus.pdf",width=4,height=4,useDingbats=F)
x1=log10(amyb_distance_rhesus[rhesus_pi_amyb,1]+1)
x2=log10(amyb_distance_rhesus[rhesus_pi_noamyb,1]+1)
x3=log10(amyb_distance_rhesus[pl,1]+1)
x4=log10(amyb_distance_rhesus[rhesus_pc,1]+1)
boxplot(x1,x2,x3,x4,yaxt="n",
        names=c("amyb\npositive\npachy\npiCluster","amyb\nnegative\npachy\npiCluster","pathway","protein\ncoding"),
        border=csl[c(1,1,2,7)],outline=F)
points(runif(length(x1),0.7,1.3),x1,pch=21,col="white",bg=csl[1],lwd=0.7,cex=0.7)
points(runif(length(x2),1.7,2.3),x2,pch=21,col="white",bg=csl[1],lwd=0.7,cex=0.7)
points(runif(length(x3),2.7,3.3),x3,pch=21,col="white",bg=csl[2],lwd=0.7,cex=0.7)
points(runif(length(x4),3.7,4.3),x4,pch=21,col="white",bg=csl[7],lwd=0.7,cex=0.7)
add_axis(2)
dev.off()



# Sep05.2018: piRNA annotation comparison----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/piRNA_cluster_annotation/comparison/")
anno_bear=read.table("piG.exon.sum",header=F,row.names=1)
anno_xxx=read.table("previous_annotation.XXX.sum",header=F,row.names=1)
anno_jinchuan=read.table("previous_annotation.Jinchuan.sum",header=F,row.names=1)
anno_william=read.table("previous_annotation.willian.sum",header=F,row.names=1)
anno_bear1=read.table("piG.exon.uniq.sum",header=F,row.names=1)
anno_xxx1=read.table("previous_annotation.XXX.uniq.sum",header=F,row.names=1)
anno_jinchuan1=read.table("previous_annotation.Jinchuan.uniq.sum",header=F,row.names=1)
anno_william1=read.table("previous_annotation.willian.uniq.sum",header=F,row.names=1)


fun_t=function(d,ct){
  x1=c();x2=c();tl=0;ta=0
  for(i in 1:dim(d)[1]){
    if(tl<ct){tl=tl+d[i,3];ta=ta+d[i,1]/10000;x1=c(x1,tl);x2=c(x2,ta)}
  }
  om=cbind(x1,x2);return(om)
}

pdf("../../Bear_Jul12.2018/piRNA_comparison.pdf",width=8,height=4)
par(mar=c(4,4,3,1),tcl=0.3,bty="n",cex=5/6,mfrow=c(1,2))
ct=3000000
ss=fun_t(anno_bear,ct=ct)
plot(ss[,1],ss[,2],type="l",xlab="Total Mbp in piclusters",ylab="piRNA percentage(%)",
     xlim=c(0,3000000),xaxt="n",lwd=3,col="black",main="AllMappers",ylim=c(0,100))
axis(1,c(0,10,20,30)*100000,label=c(0,1,2,3))
text(c(3,3,3,3)*10^6,c(2,7,12,17),
     label=c("This study","Previous: Girard","Previous: Jinchuan","Previous: Williams"),
     col=c("black",cs_gg[2:4]),font=c(2,1,1,1),pos=2)
ss=fun_t(anno_xxx,ct=ct)
lines(ss[,1],ss[,2],lwd=2,col=cs_gg[2])
ss=fun_t(anno_jinchuan,ct=ct)
lines(ss[,1],ss[,2],lwd=2,col=cs_gg[3])
ss=fun_t(anno_william,ct=ct)
lines(ss[,1],ss[,2],lwd=2,col=cs_gg[4])
ss=fun_t(anno_bear1,ct=ct)
plot(ss[,1],ss[,2],type="l",xlab="Total Mbp in piclusters",ylab="piRNA percentage(%)",
     xlim=c(0,3000000),xaxt="n",lwd=3,col="black",main="UniqueMappers",ylim=c(0,100))
axis(1,c(0,10,20,30)*100000,label=c(0,1,2,3))
ss=fun_t(anno_xxx1,ct=ct)
lines(ss[,1],ss[,2],lwd=2,col=cs_gg[2])
ss=fun_t(anno_jinchuan1,ct=ct)
lines(ss[,1],ss[,2],lwd=2,col=cs_gg[3])
ss=fun_t(anno_william1,ct=ct)
lines(ss[,1],ss[,2],lwd=2,col=cs_gg[4])
text(c(3,3,3,3)*10^6,c(2,7,12,17),
     label=c("This study","Previous: Girard","Previous: Jinchuan","Previous: Williams"),
     col=c("black",cs_gg[2:4]),font=c(2,1,1,1),pos=2)
dev.off()




# Sep06.2018: sorted cells RNAseq----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
tdf=read.table("Sorted_Cells.rpkm",header=T,row.names=1)
tdf=tdf[,c(5,1,4,2,3,6)]
pdf("../Bear_Jul12.2018/boxplot.celltypes.pdf",width=6,height=4,useDingbats=F)
boxplot(log10(tdf[human_pachy,]+0.1),at=c(1:6)*5-4,col=csl[1],
        xaxt="n",xlab="",ylab="RPKM",xlim=c(0,30),pch=20,cex=0.3)
boxplot(log10(tdf[human_prepachy,]+0.1),at=c(1:6)*5-3,col=csl[3],add=T,
        xaxt="n",xlab="",ylab="",yaxt="n",pch=20,cex=0.3)
boxplot(log10(tdf[human_hybrid,]+0.1),at=c(1:6)*5-2,col=csl[4],add=T,
        xaxt="n",xlab="",ylab="",yaxt="n",pch=20,cex=0.3)
# boxplot(log10(tdf[pl,]+0.1),at=c(1:6)*5-1,col=csl[2],add=T,
#         xaxt="n",xlab="",ylab="",yaxt="")
axis(1,(1:6)*5-2,label=colnames(tdf),lwd=0,las=3)
dev.off()
pdf("../Bear_Jul12.2018/boxplot.celltypes.linear.pdf",width=6,height=4,useDingbats=F)
boxplot(tdf[human_pachy,],at=c(1:6)*5-4,col=csl[1],
        xaxt="n",xlab="",ylab="RPKM",xlim=c(0,30),pch=20,cex=0.3,ylim=c(0,100))
boxplot(tdf[human_prepachy,],at=c(1:6)*5-3,col=csl[3],add=T,
        xaxt="n",xlab="",ylab="",yaxt="n",pch=20,cex=0.3)
boxplot(tdf[human_hybrid,],at=c(1:6)*5-2,col=csl[4],add=T,
        xaxt="n",xlab="",ylab="",yaxt="n",pch=20,cex=0.3)
# boxplot(log10(tdf[pl,]+0.1),at=c(1:6)*5-1,col=csl[2],add=T,
#         xaxt="n",xlab="",ylab="",yaxt="")
axis(1,(1:6)*5-2,label=colnames(tdf),lwd=0,las=3)
dev.off()

pdf("../Bear_Jul12.2018/X_inactive.pdf",width=4,height=4,useDingbats=F)
tn=c("X-p22-3627","pi-KANTR")
plot(unlist(tdf[tn[1],]),col=cs_gg[3],pch=19,ylim=c(0,3),xaxt="n",xlab="")
points(unlist(tdf[tn[2],]),col=cs_gg[4],pch=19)
text(c(1,1),c(2.3,2.7),label=tn,pos=4,col=cs_gg[3:4])
axis(1,1:6,label=colnames(tdf),lwd=0,las=3)
dev.off()

# Sep06.2018: phastCon metaplot----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
meta_phast=read.table("phastCon.t5.u5000.d10000.b150.mat",header=F,row.names=1)
meta_phast1=read.table("phastCon.random.intergenic.t5.u5000.d10000.b150.mat",header=F,row.names=1)
meta_phylop=read.table("phylop.t5.u5000.d10000.b150.mat",header=F,row.names=1)
meta_phylop1=read.table("phylop.random.intergenic.t5.u5000.d10000.b150.mat",header=F,row.names=1)
tn=human_pachy[order(df_srna[human_pachy,5],decreasing=T)]
fun_t=function(tm,x,tcs,mn){
  par(mar=c(2,2,3,1),tcl=0.3,cex=5/6,bty="n")
  t=apply(tm[x,],2,mean,trim=0.05)
  ti=c(49:51)[order(t[49:51])[3]]
  zs=round((t[ti]-mean(t[-ti]))/sd(t[-ti]),2)
  ymi=floor(min(t)*100)/100
  ym=ceiling(max(t)*100)/100
  plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(ymi,ym),xlim=c(0,150),main=mn)
  lines(t,col=tcs)
  axis(2,c(ymi,ym),label=c(ymi,ym))
  text(50,ym*19/20,label=paste("zscore = ",zs,sep=""),pos=4)
}
pdf("../Bear_Jul12.2018/metaplot.phastCon.pdf",width=3,height=10*7/5,useDingbats=F)
par(mfrow=c(7,1))
fun_t(meta_phast,human_pc,csl[5],"protein-coding")
fun_t(meta_phast,human_linc,csl[7],"lincrna")
fun_t(meta_phast,human_prepachy_genic,csl[2],"genicPre")
fun_t(meta_phast,human_prepachy_inter,csl[2],"interPre")
fun_t(meta_phast,human_pachy_amyb,csl[1],"amybPachy")
fun_t(meta_phast,human_pachy_noamyb,csl[1],"noamybPachy")
fun_t(meta_phast1,row.names(meta_phast1),csl[9],"random intergenic")
dev.off()
pdf("../Bear_Jul12.2018/metaplot.phylop.pdf",width=3,height=10*7/5,useDingbats=F)
par(mfrow=c(7,1))
fun_t(meta_phylop,human_pc,csl[5],"protein-coding")
fun_t(meta_phylop,human_linc,csl[7],"lincrna")
fun_t(meta_phylop,human_prepachy_genic,csl[2],"genicPre")
fun_t(meta_phylop,human_prepachy_inter,csl[2],"interPre")
fun_t(meta_phylop,human_pachy_amyb,csl[1],"amybPachy")
fun_t(meta_phylop,human_pachy_noamyb,csl[1],"noamybPachy")
fun_t(meta_phylop1,row.names(meta_phast1),csl[9],"random intergenic")
dev.off()


# Sep16.2018: mouse SNP analysis----
### mouse
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/snp_analysis/")
mrna.cds=read.table("mm10.mRNA.CDS.result",header=F,row.names=NULL)
mrna.5utr=read.table("mm10.mRNA.5UTR.result",header=F,row.names=NULL)
mrna.3utr=read.table("mm10.mRNA.3UTR.result",header=F,row.names=NULL)
mrna.intron=read.table("mm10.mRNA.intron.result",header=F,row.names=NULL)
lincrna.exon=read.table("mm10.lincRNA.exon.result",header=F,row.names=NULL)
lincrna.intron=read.table("mm10.lincRNA.intron.result",header=F,row.names=NULL)
pachy1.exon=read.table("mm10.top20pachy.exon.result",header=F,row.names=NULL)
pachy2.exon=read.table("mm10.bottom80pachy.exon.result",header=F,row.names=NULL)
pachy2.intron=read.table("mm10.bottom80pachy.intron.result",header=F,row.names=NULL)
prepachy.exon=read.table("mm10.prepachy.exon.result",header=F,row.names=NULL)
prepachy.intron=read.table("mm10.prepachy.intron.result",header=F,row.names=NULL)
random.inter=read.table("mm10.random.intergenic.result",header=F,row.names=NULL)

pdf("../../Bear_Jul12.2018/SNP.density.pdf",width=9,height=6,useDingbats=F)
par(mar=c(8,4,1,1),tcl=0.3,cex=5/6,bty="n")
o=order(mrna.cds[,1])
plot(NA,xlim=c(0.5,21),ylim=c(-1,1.5),
     xlab="",ylab="SNP/kb",xaxt="n",yaxt="n")
points(rep(1,36),log10(mrna.cds[o,1]+0.1),bg=colorRampPalette(brewer.pal(9,"RdYlBu"))(36),
       pch=21,lwd=0.3,col="grey")
points(rep(2,36),log10(mrna.5utr[o,1]+0.1),bg=colorRampPalette(brewer.pal(9,"RdYlBu"))(36),
       pch=21,lwd=0.3,col="grey")
points(rep(3,36),log10(mrna.3utr[o,1]+0.1),bg=colorRampPalette(brewer.pal(9,"RdYlBu"))(36),
       pch=21,lwd=0.3,col="grey")
points(rep(4,36),log10(mrna.intron[o,1]+0.1),bg=colorRampPalette(brewer.pal(9,"RdYlBu"))(36),
       pch=21,lwd=0.3,col="grey")
points(rep(5,36),log10(lincrna.exon[o,1]+0.1),bg=colorRampPalette(brewer.pal(9,"RdYlBu"))(36),
       pch=21,lwd=0.3,col="grey")
points(rep(6,36),log10(lincrna.intron[o,1]+0.1),bg=colorRampPalette(brewer.pal(9,"RdYlBu"))(36),
       pch=21,lwd=0.3,col="grey")
points(rep(7,36),log10(prepachy.exon[o,1]+0.1),bg=colorRampPalette(brewer.pal(9,"RdYlBu"))(36),
       pch=21,lwd=0.3,col="grey")
points(rep(8,36),log10(prepachy.intron[o,1]+0.1),bg=colorRampPalette(brewer.pal(9,"RdYlBu"))(36),
       pch=21,lwd=0.3,col="grey")
points(rep(9,36),log10(pachy1.exon[o,1]+0.1),bg=colorRampPalette(brewer.pal(9,"RdYlBu"))(36),
       pch=21,lwd=0.3,col="grey")
points(rep(10,36),log10(pachy2.exon[o,1]+0.1),bg=colorRampPalette(brewer.pal(9,"RdYlBu"))(36),
       pch=21,lwd=0.3,col="grey")
points(rep(11,36),log10(pachy2.intron[o,1]+0.1),bg=colorRampPalette(brewer.pal(9,"RdYlBu"))(36),
       pch=21,lwd=0.3,col="grey")
points(rep(12,36),log10(random.inter[o,1]+0.1),bg=colorRampPalette(brewer.pal(9,"RdYlBu"))(36),
       pch=21,lwd=0.3,col="grey")
axis(1,1:12,label=c("mrna.cds","mrna.5utr","mrna.3utr","mrna.intron",
                    "lincrna.exon","lincrna.intron","prepachy.exon","prepachy.intron",
                    "top20pachy.exon","otherpachy.exon","otherpachy.intron","random.inter"),lwd=0,las=3)
names=c("129P2_OlaHsd","129S1_SvImJ","129S5SvEvBrd","AKR_J",
        "A_J","BALB_cJ","BTBR_T+_Itpr3tf_J","BUB_BnJ","C3H_HeH",
        "C3H_HeJ","C57BL_10J","C57BL_6NJ","C57BR_cdJ","C57L_J",
        "C58_J","CAST_EiJ","CBA_J","DBA_1J","DBA_2J",
        "FVB_NJ","I_LnJ","KK_HiJ","LEWES_EiJ","LP_J",
        "MOLF_EiJ","NOD_ShiLtJ","NZB_B1NJ","NZO_HlLtJ","NZW_LacJ",
        "PWK_PhJ","RF_J","SEA_GnJ","SPRET_EiJ","ST_bJ",
        "WSB_EiJ","ZALENDE_EiJ")
legend("right",legend=names[o],col=colorRampPalette(brewer.pal(9,"RdYlBu"))(36),bty="n",
       pch=19,ncol=2)
axis(2,c(log10(c(0,1,10,100)+0.1)),label=c(0,1,10,100))
axis(2,log10(c((1:9)/10,2:9,(2:9)*10)),lty=2,label=rep("",25))
dev.off()

# Aug22.2018: conservation analysis of piRNA genes----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
pha_mrna_cds=read.table("conservation_analysis/hg19.mRNA.CDS.phastCon",header=F,row.names=1,check.names=F)
pha_mrna_5utr=read.table("conservation_analysis/hg19.mRNA.5UTR.phastCon",header=F,row.names=1,check.names=F)
pha_mrna_3utr=read.table("conservation_analysis/hg19.mRNA.3UTR.phastCon",header=F,row.names=1,check.names=F)
pha_mrna_intron=read.table("conservation_analysis/hg19.mRNA.intron.phastCon",header=F,row.names=1,check.names=F)
pha_orna_exon=read.table("conservation_analysis/hg19.otherRNA.exon.phastCon",header=F,row.names=1,check.names=F)
pha_orna_intron=read.table("conservation_analysis/hg19.otherRNA.intron.phastCon",header=F,row.names=1,check.names=F)
pha_pirna_exon=read.table("conservation_analysis/piG.exon.phastCon",header=F,row.names=1,check.names=F)
pha_pirna_intron=read.table("conservation_analysis/piG.intron.phastCon",header=F,row.names=1,check.names=F)
pha_random=read.table("conservation_analysis/hg19.random.intergenic.phastCon",header=F,row.names=1,check.names=F)
phy_mrna_cds=read.table("conservation_analysis/hg19.mRNA.CDS.phylop",header=F,row.names=1,check.names=F)
phy_mrna_5utr=read.table("conservation_analysis/hg19.mRNA.5UTR.phylop",header=F,row.names=1,check.names=F)
phy_mrna_3utr=read.table("conservation_analysis/hg19.mRNA.3UTR.phylop",header=F,row.names=1,check.names=F)
phy_mrna_intron=read.table("conservation_analysis/hg19.mRNA.intron.phylop",header=F,row.names=1,check.names=F)
phy_orna_exon=read.table("conservation_analysis/hg19.otherRNA.exon.phylop",header=F,row.names=1,check.names=F)
phy_orna_intron=read.table("conservation_analysis/hg19.otherRNA.intron.phylop",header=F,row.names=1,check.names=F)
phy_pirna_exon=read.table("conservation_analysis/hg19.piRNA.exon.phylop",header=F,row.names=1,check.names=F)
phy_pirna_intron=read.table("conservation_analysis/hg19.piRNA.intron.phylop",header=F,row.names=1,check.names=F)
phy_random=read.table("conservation_analysis/hg19.random.intergenic.phylop",header=F,row.names=1,check.names=F)

setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/snp_analysis/")
human.af=read.table("merged.DAF.stats.af",header=T,row.names=NULL)
human.af_afr=read.table("merged.DAF.stats.af_afr",header=T,row.names=NULL)
human.af_amr=read.table("merged.DAF.stats.af_amr",header=T,row.names=NULL)
human.af_eas=read.table("merged.DAF.stats.af_eas",header=T,row.names=NULL)
human.af_eur=read.table("merged.DAF.stats.af_eur",header=T,row.names=NULL)
human.af_sas=read.table("merged.DAF.stats.af_sas",header=T,row.names=NULL)
human.af1=read.table("../population/merged.piG.af",header=T,row.names=NULL,check.names = F)
to=c(16,15,14,17,12,13,7,6,5,9,1,2,18,19,28)
to1=c(16,15,14,17,12,13,7,6,5,9,3,4,26,27,20,21,28)

pdf("../../Bear_Jul12.2018/conservation.analysis.pdf",width=5,height=4,useDingbats=F)
par(mar=c(8,4,2,2),tcl=0.3,bty="n",cex=5/6)
boxplot(pha_mrna_cds[,1],pha_mrna_5utr[,1],pha_mrna_3utr[,1],pha_mrna_intron[,1],
        #pha_mrna_cds[pl,1],pha_mrna_5utr[pl,1],pha_mrna_3utr[pl,1],pha_mrna_intron[pl,1],
        pha_orna_exon[human_linc,1],pha_orna_intron[human_linc,1],
        pha_mrna_cds[human_prepachy_genic1,1],pha_mrna_5utr[human_prepachy_genic1,1],
        pha_mrna_3utr[human_prepachy_genic1,1],pha_mrna_intron[human_prepachy_genic1,1],
        pha_pirna_exon[human_prepachy_inter,1],pha_pirna_intron[human_prepachy_inter,1],
        pha_pirna_exon[human_pachy_amyb,1],pha_pirna_intron[human_pachy_amyb,1],
        pha_pirna_exon[human_pachy_noamyb,1],pha_pirna_intron[human_pachy_noamyb,1],
        pha_random[,1],
        col=c(rep(csl[5],4),rep(csl[7],2),rep(csl[3],6),rep(csl[1],4),csl[9]),
        outline=F,xaxt="n",ylab="phastCon score",
        at=c(1:4,6:7,9:12,14:15,17:18,20:21,23))
abline(v=c(5,8,13,16,19,22),lty=2)
axis(1,at=c(1:4,6:7,9:12,14:15,17:18,20:21,23),
     label=c("mRNA.CDS","mRNA.5UTR","mRNA.3UTR","mRNA.intron",
             "lincRNA.exon","lincRNA.intron","genicPre.CDS","genicPre.5UTR",
             "genicPre.3UTR","genicPre.intron","interPre.exon","interPre.intron",
             "amybPachy.exon","amybPachy.intron","noamybPachy.exon","noamybPachy.intron",
             "random.inter"),las=3,lwd=0)
abline(h=(0:5)/5,lty=2,lwd=0.8,col="grey")

boxplot(phy_mrna_cds[,1],phy_mrna_5utr[,1],phy_mrna_3utr[,1],phy_mrna_intron[,1],
        #phy_mrna_cds[pl,1],phy_mrna_5utr[pl,1],phy_mrna_3utr[pl,1],phy_mrna_intron[pl,1],
        phy_orna_exon[human_linc,1],phy_orna_intron[human_linc,1],
        phy_mrna_cds[human_prepachy_genic1,1],phy_mrna_5utr[human_prepachy_genic1,1],
        phy_mrna_3utr[human_prepachy_genic1,1],phy_mrna_intron[human_prepachy_genic1,1],
        phy_pirna_exon[human_prepachy_inter,1],phy_pirna_intron[human_prepachy_inter,1],
        phy_pirna_exon[human_pachy_amyb,1],phy_pirna_intron[human_pachy_amyb,1],
        phy_pirna_exon[human_pachy_noamyb,1],phy_pirna_intron[human_pachy_noamyb,1],
        phy_random[,1],
        col=c(rep(csl[5],4),rep(csl[7],2),rep(csl[3],6),rep(csl[1],4),csl[9]),
        outline=F,xaxt="n",ylab="phylop score",
        at=c(1:4,6:7,9:12,14:15,17:18,20:21,23))
abline(v=c(5,8,13,16,19,22),lty=2)
axis(1,at=c(1:4,6:7,9:12,14:15,17:18,20:21,23),
     label=c("mRNA.CDS","mRNA.5UTR","mRNA.3UTR","mRNA.intron",
             "lincRNA.exon","lincRNA.intron","genicPre.CDS","genicPre.5UTR",
             "genicPre.3UTR","genicPre.intron","interPre.exon","interPre.intron",
             "amybPachy.exon","amybPachy.intron","noamybPachy.exon","noamybPachy.intron",
             "random.inter"),las=3,lwd=0)
abline(h=(0:5)/5,lty=2,lwd=0.8,col="grey")

par(mar=c(6,4,3,1),tcl=0.3,bty="n",cex=5/6)
plot(NA,xlim=c(0,21),ylim=c(0.02,0.08),
     xlab="",ylab="DAF",xaxt="n",main="derived allele frequency")
arrows(c(1:4,6:7,9:12,14:15,17:18,20),unlist(human.af[2,to]),
       c(1:4,6:7,9:12,14:15,17:18,20),unlist(human.af[3,to]),length=0,angle=90,
       col=c(rep(csl[5],4),rep(csl[7],2),rep(csl[3],4),rep(csl[1],4),csl[9]))
points(c(1:4,6:7,9:12,14:15,17:18,20),
       unlist(human.af[1,to]),pch=21,col="white",
       bg=c(rep(csl[5],4),rep(csl[7],2),rep(csl[3],4),rep(csl[1],4),csl[9]))
abline(v=c(5,8,13,16,19,22),lty=2)
axis(1,c(1:4,6:7,9:12,14:15,17:18,20),label=colnames(human.af)[to],lwd=0,las=3)

par(mar=c(6,4,3,1),tcl=0.3,bty="n",cex=5/6)
plot(NA,xlim=c(0,24),ylim=c(0.02,0.10),
     xlab="",ylab="DAF",xaxt="n",main="derived allele frequency")
arrows(c(1:4,6:7,9:12,14:15,17:18,20:21,23),unlist(human.af[2,to1]),
       c(1:4,6:7,9:12,14:15,17:18,20:21,23),unlist(human.af[3,to1]),length=0,angle=90,
       col=c(rep(csl[5],4),rep(csl[7],2),rep(csl[3],4),rep(csl[1],6),csl[9]))
points(c(1:4,6:7,9:12,14:15,17:18,20:21,23),
     unlist(human.af[1,to1]),pch=21,col="white",
     bg=c(rep(csl[5],4),rep(csl[7],2),rep(csl[3],4),rep(csl[1],6),csl[9]))
abline(v=c(5,8,13,16,19,22),lty=2)
axis(1,c(1:4,6:7,9:12,14:15,17:18,20:21,23),label=colnames(human.af)[to1],lwd=0,las=3)
dev.off()

pdf("../../Bear_Jul12.2018/conservation.analysis.pachytene.pdf",width=3.5,height=4,useDingbats = F)
par(mar=c(6,4,3,1),tcl=0.3,bty="n",cex=5/6)
plot(NA,xlim=c(0.5,7.5),ylim=c(0.05,0.10),
     xlab="",ylab="DAF",xaxt="n",main="derived allele frequency")
arrows(c(1:3,5:7),unlist(human.af1[2,]),
       c(1:3,5:7),unlist(human.af1[3,]),length=0,angle=90)
points(c(1:3,5:7),
       unlist(human.af1[1,]),pch=21,col="white",bg="black")
abline(v=c(4),lty=2)
axis(1,c(1:3,5:7),label=colnames(human.af1),lwd=0,las=3)
dev.off()

# Oct05.2018: exon-intron metaplot for conservation----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/exon_intron_conservation/")
ei_daf=read.table("around_pachy_intron.DAF.t53.u2000.d2000.b100.mat",header=F,row.names=1)
ei_phastCon=read.table("around_pachy_intron.phasCon.t53.u2000.d2000.b100.mat",header=F,row.names=1)
ei_rpm=read.table("around_pachy_intron.rpm.t53.u2000.d2000.b100.mat",header=F,row.names=1)
#ei_daf=ei_daf[!is.na(apply(ei_daf,1,mean)),]
#ei_rpm=ei_rpm[!is.na(apply(ei_rpm,1,mean)),]
#ei_phastCon=ei_phastCon[!is.na(apply(ei_phastCon,1,mean)),]
fun_t=function(d,yl,tm){
  par(bty="n",cex=5/6,mar=c(2,4,1,1))
  td=apply(d,2,mean,trim=tm)
  ym=max(td)
  plot(td,xaxt="n",ylab=yl,ylim=c(0,ym),type="l")
  lines(c(35,66),c(0,0))
  lines(c(1,35),c(0,0),lwd=3)
  lines(c(66,100),c(0,0),lwd=3)
  axis(1,c(35/2,50,65/2+50),label=c("exon","intron","exon"),lwd=0)
  abline(v=c(35,66),lwd=1,lty=2)
}
pdf("../../Bear_Jul12.2018/exon_intron_conservation.pdf",width=5,height=9)
par(mfrow=c(3,1))
fun_t(ei_daf,"DAF",tm=0)
fun_t(ei_rpm,"ppm",tm=0.05)
fun_t(ei_phastCon,"phastCon",tm=0.05)
dev.off()




# Oct09.2018: population conservation1----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
chop.pi=read.table("merged.picluster.chopped.rpm",header=T,row.names=1,check.names=F)
t=unlist(strsplit(row.names(chop.pi),"_"))
t=t[(1:(length(t)/2))*2-1]
n1=c();n2=c()
for(i in human_pachy){n1=c(n1,which(t==i))}
for(i in human_prepachy){n2=c(n2,which(t==i))}
fun_t=function(i1,i2){
  x1=log10(chop.pi[n1,i1]+1)
  x2=log10(chop.pi[n1,i2]+1)
  m=max(c(x1,x2))
  plot(x1,x2,xlim=c(0,m),ylim=c(0,m),xaxt="n",yaxt="n",xlab=colnames(chop.pi)[i1],
       ylab=colnames(chop.pi)[i2])
  c=cor(chop.pi[n1,i1],chop.pi[n1,i2])
  add_axis(1);add_axis(2)
  text(0,m*9/10,label=paste("pearson cor = ",round(c,3)),pos=4)
}
pdf("../Bear_Jul12.2018/chopped.picluster.ppm.pdf",width=5,height=5)
par(mar=c(4,4,1,1),cex=5/6,tcl=0.3,bty="n")
for(i in 1:6){
  for(j in 1:6){
    if(i<j){fun_t(i,j)}
  }
}
dev.off()



# Oct11.2018: population conservation2-16nt----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/population/")
list_highA=list()
list_lowA=list()
elements=c("mRNA.CDS","mRNA.3UTR","mRNA.5UTR","mRNA.intron","lincRNA.exon",
           "lincRNA.intron","prepachy.exon","prepachy.intron","pachy.exon",
           "pachy.intron","others")
for(i in elements){
  list_highA[[i]]=read.table(paste("final.highAbun.hg19",i,"result",sep="."),header=F,row.names=NULL)
  list_lowA[[i]]=read.table(paste("final.lowAbun.hg19",i,"result",sep="."),header=F,row.names=NULL)
}
cdf=as.vector(read.table("cdf.tab",header=F,row.names=NULL)[,1])
abun_highA=read.table("final.highAbun.abundance.result",header=F,row.names=NULL)
abun_lowA=read.table("final.lowAbun.abundance.result",header=F,row.names=NULL)

### CDF
pdf("../../Bear_Jul12.2018/CDF.piRNA.perloci.pdf",width=4,height=4,useDingbats=F)
par(mar=c(4,4,3,1),tcl=0.3,cex=5/6,bty="n")
plot(cdf*100/cdf[length(cdf)],type="l",xlab="seeds producing loci number",
     ylab="percentage of produced piRNAs",main="CDF of piRNAs in seed producing loci")
abline(h=90.59,lty=2)
abline(v=200000,lty=2)
text(200000,50,label="RPM cutoff for low and high\nabundance loci: 0.3284",pos=4)
text(300000,90.59,label="90.59%",pos=1)
dev.off()
### seed conservation
pdf("../../Bear_Jul12.2018/seedConservation.pdf",width=7,height=8,useDingbats = F)
number=c(nrow(list_highA$mRNA.CDS),nrow(list_lowA$mRNA.CDS),
         nrow(list_highA$mRNA.3UTR),nrow(list_lowA$mRNA.3UTR),
         nrow(list_highA$mRNA.5UTR),nrow(list_lowA$mRNA.5UTR),
         nrow(list_highA$mRNA.intron),nrow(list_lowA$mRNA.intron),
         nrow(list_highA$lincRNA.exon),nrow(list_lowA$lincRNA.exon),
         nrow(list_highA$lincRNA.intron),nrow(list_lowA$lincRNA.intron),
         nrow(list_highA$prepachy.exon),nrow(list_lowA$prepachy.exon),
         nrow(list_highA$prepachy.intron),nrow(list_lowA$prepachy.intron),
         nrow(list_highA$pachy.exon),nrow(list_lowA$pachy.exon),
         nrow(list_highA$pachy.intron),nrow(list_lowA$pachy.intron),
         nrow(list_highA$others),nrow(list_lowA$others))
number=prettyNum(number,big.mark=",")
par(mfrow=c(2,1),mar=c(1,4,1,1),tcl=0.3,bty="n",cex=5/6)
boxplot(list_highA$mRNA.CDS[,5],list_lowA$mRNA.CDS[,5],
        list_highA$mRNA.3UTR[,5],list_lowA$mRNA.3UTR[,5],
        list_highA$mRNA.5UTR[,5],list_lowA$mRNA.5UTR[,5],
        list_highA$mRNA.intron[,5],list_lowA$mRNA.intron[,5],
        list_highA$lincRNA.exon[,5],list_lowA$lincRNA.exon[,5],
        list_highA$lincRNA.intron[,5],list_lowA$lincRNA.intron[,5],
        list_highA$prepachy.exon[,5],list_lowA$prepachy.exon[,5],
        list_highA$prepachy.intron[,5],list_lowA$prepachy.intron[,5],
        list_highA$pachy.exon[,5],list_lowA$pachy.exon[,5],
        list_highA$pachy.intron[,5],list_lowA$pachy.intron[,5],
        list_highA$others[,5],list_lowA$others[,5],
        outline=F,staplewex=0,lty=1,
        at=c(1:2,4:5,7:8,10:11,13:14,16:17,19:20,22:23,25:26,28:29,31:32),
        col=csl[1:2],xaxt="n",ylab="cosine similarity")
legend("bottomleft",fill=csl[1:2],legend=c("high abundance","low abundance"))

mean_vector=c(mean(list_highA$mRNA.CDS[,5]),mean(list_lowA$mRNA.CDS[,5]),
              mean(list_highA$mRNA.3UTR[,5]),mean(list_lowA$mRNA.3UTR[,5]),
              mean(list_highA$mRNA.5UTR[,5]),mean(list_lowA$mRNA.5UTR[,5]),
              mean(list_highA$mRNA.intron[,5]),mean(list_lowA$mRNA.intron[,5]),
              mean(list_highA$lincRNA.exon[,5]),mean(list_lowA$lincRNA.exon[,5]),
              mean(list_highA$lincRNA.intron[,5]),mean(list_lowA$lincRNA.intron[,5]),
              mean(list_highA$prepachy.exon[,5]),mean(list_lowA$prepachy.exon[,5]),
              mean(list_highA$prepachy.intron[,5]),mean(list_lowA$prepachy.intron[,5]),
              mean(list_highA$pachy.exon[,5]),mean(list_lowA$pachy.exon[,5]),
              mean(list_highA$pachy.intron[,5]),mean(list_lowA$pachy.intron[,5]),
              mean(list_highA$others[,5]),mean(list_lowA$others[,5]))
# sd_vector=c(sd(list_highA$mRNA.CDS[,5]),sd(list_lowA$mRNA.CDS[,5]),
#             sd(list_highA$mRNA.3UTR[,5]),sd(list_lowA$mRNA.3UTR[,5]),
#             sd(list_highA$mRNA.5UTR[,5]),sd(list_lowA$mRNA.5UTR[,5]),
#             sd(list_highA$mRNA.intron[,5]),sd(list_lowA$mRNA.intron[,5]),
#             sd(list_highA$lincRNA.exon[,5]),sd(list_lowA$lincRNA.exon[,5]),
#             sd(list_highA$lincRNA.intron[,5]),sd(list_lowA$lincRNA.intron[,5]),
#             sd(list_highA$prepachy.exon[,5]),sd(list_lowA$prepachy.exon[,5]),
#             sd(list_highA$prepachy.intron[,5]),sd(list_lowA$prepachy.intron[,5]),
#             sd(list_highA$pachy.exon[,5]),sd(list_lowA$pachy.exon[,5]),
#             sd(list_highA$pachy.intron[,5]),sd(list_lowA$pachy.intron[,5]),
#             sd(list_highA$others[,5]),sd(list_lowA$others[,5]))
plot(c(1:2,4:5,7:8,10:11,13:14,16:17,19:20,22:23,25:26,28:29,31:32),
     mean_vector,col=csl[1:2],lwd=2,
     ylab="average cosine similarity",xaxt="n",ylim=c(0.98,1),yaxt="n")
axis(2,c(0.9,0.92,0.94,0.96,0.98,1),label=c(0.9,0.92,0.94,0.96,0.98,1))
text(c(1:2,4:5,7:8,10:11,13:14,16:17,19:20,22:23,25:26,28:29,31:32),rep(0.88,22),
     label=paste("n=",number,sep=""),srt=45,cex=0.6)
text((1:11)*3-1.5,rep(0.86,11),label=elements,srt=45,cex=5/6)
dev.off()



# Oct11.2018: population conservation2-11nt----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/population/2-11nt/")
list_highA=list()
list_lowA=list()
elements=c("mRNA.CDS","mRNA.3UTR","mRNA.5UTR","mRNA.intron","lincRNA.exon",
           "lincRNA.intron","prepachy.exon","prepachy.intron","pachy.exon",
           "pachy.intron","others")
for(i in elements){
  list_highA[[i]]=read.table(paste("final.highAbun.hg19",i,"result",sep="."),header=F,row.names=NULL)
  list_lowA[[i]]=read.table(paste("final.lowAbun.hg19",i,"result",sep="."),header=F,row.names=NULL)
}
#abun_highA=read.table("final.highAbun.abundance.result",header=F,row.names=NULL)
#abun_lowA=read.table("final.lowAbun.abundance.result",header=F,row.names=NULL)
### seed conservation
pdf("../../../Bear_Jul12.2018/seedConservation.2-11nt.pdf",width=7,height=8,useDingbats = F)
number=c(nrow(list_highA$mRNA.CDS),nrow(list_lowA$mRNA.CDS),
         nrow(list_highA$mRNA.3UTR),nrow(list_lowA$mRNA.3UTR),
         nrow(list_highA$mRNA.5UTR),nrow(list_lowA$mRNA.5UTR),
         nrow(list_highA$mRNA.intron),nrow(list_lowA$mRNA.intron),
         nrow(list_highA$lincRNA.exon),nrow(list_lowA$lincRNA.exon),
         nrow(list_highA$lincRNA.intron),nrow(list_lowA$lincRNA.intron),
         nrow(list_highA$prepachy.exon),nrow(list_lowA$prepachy.exon),
         nrow(list_highA$prepachy.intron),nrow(list_lowA$prepachy.intron),
         nrow(list_highA$pachy.exon),nrow(list_lowA$pachy.exon),
         nrow(list_highA$pachy.intron),nrow(list_lowA$pachy.intron),
         nrow(list_highA$others),nrow(list_lowA$others))
number=prettyNum(number,big.mark=",")
par(mfrow=c(2,1),mar=c(1,4,1,1),tcl=0.3,bty="n",cex=5/6)
boxplot(list_highA$mRNA.CDS[,5],list_lowA$mRNA.CDS[,5],
        list_highA$mRNA.3UTR[,5],list_lowA$mRNA.3UTR[,5],
        list_highA$mRNA.5UTR[,5],list_lowA$mRNA.5UTR[,5],
        list_highA$mRNA.intron[,5],list_lowA$mRNA.intron[,5],
        list_highA$lincRNA.exon[,5],list_lowA$lincRNA.exon[,5],
        list_highA$lincRNA.intron[,5],list_lowA$lincRNA.intron[,5],
        list_highA$prepachy.exon[,5],list_lowA$prepachy.exon[,5],
        list_highA$prepachy.intron[,5],list_lowA$prepachy.intron[,5],
        list_highA$pachy.exon[,5],list_lowA$pachy.exon[,5],
        list_highA$pachy.intron[,5],list_lowA$pachy.intron[,5],
        list_highA$others[,5],list_lowA$others[,5],
        outline=F,staplewex=0,lty=1,
        at=c(1:2,4:5,7:8,10:11,13:14,16:17,19:20,22:23,25:26,28:29,31:32),
        col=csl[1:2],xaxt="n",ylab="cosine similarity")
legend("bottomleft",fill=csl[1:2],legend=c("high abundance","low abundance"))

mean_vector=c(mean(list_highA$mRNA.CDS[,5]),mean(list_lowA$mRNA.CDS[,5]),
              mean(list_highA$mRNA.3UTR[,5]),mean(list_lowA$mRNA.3UTR[,5]),
              mean(list_highA$mRNA.5UTR[,5]),mean(list_lowA$mRNA.5UTR[,5]),
              mean(list_highA$mRNA.intron[,5]),mean(list_lowA$mRNA.intron[,5]),
              mean(list_highA$lincRNA.exon[,5]),mean(list_lowA$lincRNA.exon[,5]),
              mean(list_highA$lincRNA.intron[,5]),mean(list_lowA$lincRNA.intron[,5]),
              mean(list_highA$prepachy.exon[,5]),mean(list_lowA$prepachy.exon[,5]),
              mean(list_highA$prepachy.intron[,5]),mean(list_lowA$prepachy.intron[,5]),
              mean(list_highA$pachy.exon[,5]),mean(list_lowA$pachy.exon[,5]),
              mean(list_highA$pachy.intron[,5]),mean(list_lowA$pachy.intron[,5]),
              mean(list_highA$others[,5]),mean(list_lowA$others[,5]))
# sd_vector=c(sd(list_highA$mRNA.CDS[,5]),sd(list_lowA$mRNA.CDS[,5]),
#             sd(list_highA$mRNA.3UTR[,5]),sd(list_lowA$mRNA.3UTR[,5]),
#             sd(list_highA$mRNA.5UTR[,5]),sd(list_lowA$mRNA.5UTR[,5]),
#             sd(list_highA$mRNA.intron[,5]),sd(list_lowA$mRNA.intron[,5]),
#             sd(list_highA$lincRNA.exon[,5]),sd(list_lowA$lincRNA.exon[,5]),
#             sd(list_highA$lincRNA.intron[,5]),sd(list_lowA$lincRNA.intron[,5]),
#             sd(list_highA$prepachy.exon[,5]),sd(list_lowA$prepachy.exon[,5]),
#             sd(list_highA$prepachy.intron[,5]),sd(list_lowA$prepachy.intron[,5]),
#             sd(list_highA$pachy.exon[,5]),sd(list_lowA$pachy.exon[,5]),
#             sd(list_highA$pachy.intron[,5]),sd(list_lowA$pachy.intron[,5]),
#             sd(list_highA$others[,5]),sd(list_lowA$others[,5]))
plot(c(1:2,4:5,7:8,10:11,13:14,16:17,19:20,22:23,25:26,28:29,31:32),
     mean_vector,col=csl[1:2],lwd=2,
     ylab="average cosine similarity",xaxt="n",ylim=c(0.96,1),yaxt="n")
axis(2,c(0.9,0.92,0.94,0.96,0.98,1),label=c(0.9,0.92,0.94,0.96,0.98,1))
text(c(1:2,4:5,7:8,10:11,13:14,16:17,19:20,22:23,25:26,28:29,31:32),rep(0.98,22),
     label=paste("n=",number,sep=""),srt=45,cex=0.6)
text((1:11)*3-1.5,rep(0.97,11),label=elements,srt=45,cex=5/6)
dev.off()



# Oct11.2018: population conservation17-25nt----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/population/12-25nt/")
list_highA=list()
list_lowA=list()
elements=c("mRNA.CDS","mRNA.3UTR","mRNA.5UTR","mRNA.intron","lincRNA.exon",
           "lincRNA.intron","prepachy.exon","prepachy.intron","pachy.exon",
           "pachy.intron","others")
for(i in elements){
  list_highA[[i]]=read.table(paste("final.highAbun.hg19",i,"result",sep="."),header=F,row.names=NULL)
  list_lowA[[i]]=read.table(paste("final.lowAbun.hg19",i,"result",sep="."),header=F,row.names=NULL)
}
abun_highA=read.table("final.highAbun.abundance.result",header=F,row.names=NULL)
abun_lowA=read.table("final.lowAbun.abundance.result",header=F,row.names=NULL)
### seed conservation
pdf("../../../Bear_Jul12.2018/seedConservation.12-25nt.pdf",width=7,height=8,useDingbats = F)
number=c(nrow(list_highA$mRNA.CDS),nrow(list_lowA$mRNA.CDS),
         nrow(list_highA$mRNA.3UTR),nrow(list_lowA$mRNA.3UTR),
         nrow(list_highA$mRNA.5UTR),nrow(list_lowA$mRNA.5UTR),
         nrow(list_highA$mRNA.intron),nrow(list_lowA$mRNA.intron),
         nrow(list_highA$lincRNA.exon),nrow(list_lowA$lincRNA.exon),
         nrow(list_highA$lincRNA.intron),nrow(list_lowA$lincRNA.intron),
         nrow(list_highA$prepachy.exon),nrow(list_lowA$prepachy.exon),
         nrow(list_highA$prepachy.intron),nrow(list_lowA$prepachy.intron),
         nrow(list_highA$pachy.exon),nrow(list_lowA$pachy.exon),
         nrow(list_highA$pachy.intron),nrow(list_lowA$pachy.intron),
         nrow(list_highA$others),nrow(list_lowA$others))
number=prettyNum(number,big.mark=",")
par(mfrow=c(2,1),mar=c(1,4,1,1),tcl=0.3,bty="n",cex=5/6)
boxplot(list_highA$mRNA.CDS[,5],list_lowA$mRNA.CDS[,5],
        list_highA$mRNA.3UTR[,5],list_lowA$mRNA.3UTR[,5],
        list_highA$mRNA.5UTR[,5],list_lowA$mRNA.5UTR[,5],
        list_highA$mRNA.intron[,5],list_lowA$mRNA.intron[,5],
        list_highA$lincRNA.exon[,5],list_lowA$lincRNA.exon[,5],
        list_highA$lincRNA.intron[,5],list_lowA$lincRNA.intron[,5],
        list_highA$prepachy.exon[,5],list_lowA$prepachy.exon[,5],
        list_highA$prepachy.intron[,5],list_lowA$prepachy.intron[,5],
        list_highA$pachy.exon[,5],list_lowA$pachy.exon[,5],
        list_highA$pachy.intron[,5],list_lowA$pachy.intron[,5],
        list_highA$others[,5],list_lowA$others[,5],
        outline=F,staplewex=0,lty=1,
        at=c(1:2,4:5,7:8,10:11,13:14,16:17,19:20,22:23,25:26,28:29,31:32),
        col=csl[1:2],xaxt="n",ylab="cosine similarity")
legend("bottomleft",fill=csl[1:2],legend=c("high abundance","low abundance"))

mean_vector=c(mean(list_highA$mRNA.CDS[,5]),mean(list_lowA$mRNA.CDS[,5]),
              mean(list_highA$mRNA.3UTR[,5]),mean(list_lowA$mRNA.3UTR[,5]),
              mean(list_highA$mRNA.5UTR[,5]),mean(list_lowA$mRNA.5UTR[,5]),
              mean(list_highA$mRNA.intron[,5]),mean(list_lowA$mRNA.intron[,5]),
              mean(list_highA$lincRNA.exon[,5]),mean(list_lowA$lincRNA.exon[,5]),
              mean(list_highA$lincRNA.intron[,5]),mean(list_lowA$lincRNA.intron[,5]),
              mean(list_highA$prepachy.exon[,5]),mean(list_lowA$prepachy.exon[,5]),
              mean(list_highA$prepachy.intron[,5]),mean(list_lowA$prepachy.intron[,5]),
              mean(list_highA$pachy.exon[,5]),mean(list_lowA$pachy.exon[,5]),
              mean(list_highA$pachy.intron[,5]),mean(list_lowA$pachy.intron[,5]),
              mean(list_highA$others[,5]),mean(list_lowA$others[,5]))
# sd_vector=c(sd(list_highA$mRNA.CDS[,5]),sd(list_lowA$mRNA.CDS[,5]),
#             sd(list_highA$mRNA.3UTR[,5]),sd(list_lowA$mRNA.3UTR[,5]),
#             sd(list_highA$mRNA.5UTR[,5]),sd(list_lowA$mRNA.5UTR[,5]),
#             sd(list_highA$mRNA.intron[,5]),sd(list_lowA$mRNA.intron[,5]),
#             sd(list_highA$lincRNA.exon[,5]),sd(list_lowA$lincRNA.exon[,5]),
#             sd(list_highA$lincRNA.intron[,5]),sd(list_lowA$lincRNA.intron[,5]),
#             sd(list_highA$prepachy.exon[,5]),sd(list_lowA$prepachy.exon[,5]),
#             sd(list_highA$prepachy.intron[,5]),sd(list_lowA$prepachy.intron[,5]),
#             sd(list_highA$pachy.exon[,5]),sd(list_lowA$pachy.exon[,5]),
#             sd(list_highA$pachy.intron[,5]),sd(list_lowA$pachy.intron[,5]),
#             sd(list_highA$others[,5]),sd(list_lowA$others[,5]))
plot(c(1:2,4:5,7:8,10:11,13:14,16:17,19:20,22:23,25:26,28:29,31:32),
     mean_vector,col=csl[1:2],lwd=2,
     ylab="average cosine similarity",xaxt="n",ylim=c(0.96,1),yaxt="n")
axis(2,c(0.9,0.92,0.94,0.96,0.98,1),label=c(0.9,0.92,0.94,0.96,0.98,1))
text(c(1:2,4:5,7:8,10:11,13:14,16:17,19:20,22:23,25:26,28:29,31:32),rep(0.98,22),
     label=paste("n=",number,sep=""),srt=45,cex=0.6)
text((1:11)*3-1.5,rep(0.97,11),label=elements,srt=45,cex=5/6)
dev.off()



# Oct22.2018: rmsk content----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
exon.rmsk=read.table("hg19.gene+pi.exon.rmsk.content.txt",header=F,row.names=1)
intron.rmsk=read.table("hg19.gene+pi.intron.rmsk.content.txt",header=F,row.names=1)
fun_t=function(d,n){
  t=d[intersect(row.names(d),n),]
  return(((t[,2]+t[,3])/t[1])[,1]*100)
}
pdf("../Bear_Jul12.2018/rmsk.content.boxplot.pdf",width=4,height=4)
par(mar=c(6,4,3,1),bty="n",tcl=0.3,cex=5/6)
boxplot(fun_t(exon.rmsk,human_pachy),fun_t(exon.rmsk,human_prepachy),
        fun_t(exon.rmsk,human_hybrid),fun_t(exon.rmsk,human_pc),
        fun_t(exon.rmsk,human_linc),col=csl[c(1,2,3,7,9)],xaxt="n",
        ylab="percentage (%)",main="rmsk content in exons",outline=F,lty=1)
axis(1,1:5,label=c("pachytene","prepachy","hybrid","mRNA","lincRNA"),las=3,lwd=0)
abline(h=0.4472413*100,lty=2)
boxplot(fun_t(intron.rmsk,human_pachy),fun_t(intron.rmsk,human_prepachy),
        fun_t(intron.rmsk,human_hybrid),fun_t(intron.rmsk,human_pc),
        fun_t(intron.rmsk,human_linc),col=csl[c(1,2,3,7,9)],xaxt="n",
        ylab="percentage (%)",main="rmsk content in introns",outline=F,lty=1)
axis(1,1:5,label=c("pachytene","prepachy","hybrid","mRNA","lincRNA"),las=3,lwd=0)
abline(h=0.4472413*100,lty=2)
dev.off()

pdf("../Bear_Jul12.2018/rmsk.content.pdf",width=8,height=3.5)
par(mar=c(8,4,1,1),tcl=0.3,cex=1/2,bty="n",las=3,lwd=0.5)
n=c(human_pachy,human_hybrid,human_prepachy)[182:1]
ci=as.integer(unlist(piG.rmsk[n,2]/rowSums(piG.rmsk[n,2:3]+0.001)*100))+1
barplot(rowSums(piG.rmsk[n,2:3])/piG.rmsk[n,1]*100,
        col=colorRampPalette(brewer.pal(9,"RdYlBu"))(101)[ci],ylab="rmsk percentage",
        border="grey")
abline(h=0.4472413*100,lty=2)
legend("topleft",fill=brewer.pal(9,"RdYlBu"),legend=c("100%",rep("",7),"0%"),
       title="anti-sense",y.intersp = 0.5)
dev.off()

# Oct24.2018: DAF in seeds----
setwd("/Users/Bigbear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/population/new")
cdf=read.table("cdf.txt",header=F,row.names=NULL)[,1]
daf.h1=read.table("merged.2-16.high.stats.af",header=T,row.names=NULL)
daf.l1=read.table("merged.2-16.low.stats.af",header=T,row.names=NULL)
daf.o1=read.table("merged.2-16.other.stats.af",header=T,row.names=NULL)
daf.h2=read.table("merged.17-25.high.stats.af",header=T,row.names=NULL)
daf.l2=read.table("merged.17-25.low.stats.af",header=T,row.names=NULL)
daf.o2=read.table("merged.17-25.other.stats.af",header=T,row.names=NULL)

cdf1=cdf[c(1:1000)*as.integer(length(cdf)/1000)]/10000
pdf("../../../Bear_Jul12.2018/CDF.piRNA.perloci.pdf",width=4,height=4,useDingbats=F)
par(mar=c(4,4,3,1),tcl=0.3,cex=5/6,bty="n")
# plot(cdf*100/cdf[length(cdf)],type="l",xlab="seeds producing loci number",
#      ylab="percentage of produced piRNAs",main="CDF of piRNAs in seed producing loci")
# abline(h=66.666,lty=2)
# abline(v=481321,lty=2)
# text(481321,50,label="RPM cutoff for low and high\nabundance loci: 0.132014",pos=4)
# text(681321,66.666,label="66.66%",pos=1)
t=length(cdf)/1000
plot(c(0,cdf1),type="l",xlab="",ylab="",xaxt="n",ylim=c(0,100))
axis(1,c(0,5,10,15,20,25)*10^6/t,label=c("0","5M","10M","15M","20M","25M"))
abline(h=66.666,lty=2)
abline(v=481321/t,lty=2)
dev.off()

to=c(7,6,5,8,3,4,19,20,11,12,13,14,17,18,15,16,21)
to1=c(7,6,5,8,3,4,19,20,11,12,1,2,9,10,21)
pdf("../../../Bear_Jul12.2018/seedConservation.pdf",width=6,height=4,useDingbats = F)
par(mar=c(6,4,3,1),tcl=0.3,bty="n",cex=5/6)
plot(NA,xlim=c(0,23),ylim=c(0.02,0.12),
     xlab="",ylab="DAF",xaxt="n",main="derived allele frequency for 2-16nt seed")
arrows(c(1:4,6:7,9:10,12:13,15:20,22)-0.2,unlist(daf.h1[2,to]),
       c(1:4,6:7,9:10,12:13,15:20,22)-0.2,unlist(daf.h1[3,to]),length=0,angle=90,
       col=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,1,1,9)])
points(c(1:4,6:7,9:10,12:13,15:20,22)-0.2,
       unlist(daf.h1[1,to]),pch=21,col="white",bg=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,1,1,9)])
arrows(c(1:4,6:7,9:10,12:13,15:20,22),unlist(daf.l1[2,to]),
       c(1:4,6:7,9:10,12:13,15:20,22),unlist(daf.l1[3,to]),length=0,angle=90,
       col=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,1,1,9)],lty=2)
points(c(1:4,6:7,9:10,12:13,15:20,22),
       unlist(daf.l1[1,to]),pch=21,col="white",bg=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,1,1,9)])
arrows(c(1:4,6:7,9:10,12:13,15:20,22)+0.2,unlist(daf.o1[2,to]),
       c(1:4,6:7,9:10,12:13,15:20,22)+0.2,unlist(daf.o1[3,to]),length=0,angle=90,
       col=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,1,1,9)])
points(c(1:4,6:7,9:10,12:13,15:20,22)+0.2,pch=21,
       unlist(daf.o1[1,to]),col=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,1,1,9)])
abline(v=c(5,8,11,14),lty=2)
axis(1,c(1:4,6:7,9:10,12:13,15:20,22),label=colnames(daf.h1)[to],lwd=0,las=3)
legend("topleft",lty=c(1,2,1),pch=c(19,19,1),legend=c("highAbun","lowABun","others"))
plot(NA,xlim=c(0,23),ylim=c(0.02,0.12),
     xlab="",ylab="DAF",xaxt="n",main="derived allele frequency for 17-25nt seed")
arrows(c(1:4,6:7,9:10,12:13,15:20,22)-0.2,unlist(daf.h2[2,to]),
       c(1:4,6:7,9:10,12:13,15:20,22)-0.2,unlist(daf.h2[3,to]),length=0,angle=90,
       col=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,1,1,9)])
points(c(1:4,6:7,9:10,12:13,15:20,22)-0.2,
       unlist(daf.h2[1,to]),pch=21,col="white",bg=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,1,1,9)])
arrows(c(1:4,6:7,9:10,12:13,15:20,22),unlist(daf.l2[2,to]),
       c(1:4,6:7,9:10,12:13,15:20,22),unlist(daf.l2[3,to]),length=0,angle=90,
       col=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,1,1,9)],lty=2)
points(c(1:4,6:7,9:10,12:13,15:20,22),
       unlist(daf.l2[1,to]),pch=21,col="white",bg=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,1,1,9)])
arrows(c(1:4,6:7,9:10,12:13,15:20,22)+0.2,unlist(daf.o2[2,to]),
       c(1:4,6:7,9:10,12:13,15:20,22)+0.2,unlist(daf.o2[3,to]),length=0,angle=90,
       col=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,1,1,9)])
points(c(1:4,6:7,9:10,12:13,15:20,22)+0.2,pch=21,
       unlist(daf.o2[1,to]),col=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,1,1,9)])
abline(v=c(5,8,11,14),lty=2)
axis(1,c(1:4,6:7,9:10,12:13,15:20,22),label=colnames(daf.h1)[to],lwd=0,las=3)
legend("topleft",lty=c(1,2,1),pch=c(19,19,1),legend=c("highAbun","lowABun","others"))

par(mar=c(6,4,3,1),tcl=0.3,bty="n",cex=5/6)
plot(NA,xlim=c(0,23),ylim=c(0.02,0.12),
     xlab="",ylab="DAF",xaxt="n",main="derived allele frequency for 2-16nt seed")
arrows(c(1:4,6:7,9:10,12:13,15:18,20)-0.2,unlist(daf.h1[2,to1]),
       c(1:4,6:7,9:10,12:13,15:18,20)-0.2,unlist(daf.h1[3,to1]),length=0,angle=90,
       col=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,9)])
points(c(1:4,6:7,9:10,12:13,15:18,20)-0.2,
       unlist(daf.h1[1,to1]),pch=21,col="white",bg=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,,1,9)])
arrows(c(1:4,6:7,9:10,12:13,15:18,20),unlist(daf.l1[2,to1]),
       c(1:4,6:7,9:10,12:13,15:18,20),unlist(daf.l1[3,to1]),length=0,angle=90,
       col=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,9)],lty=2)
points(c(1:4,6:7,9:10,12:13,15:18,20),
       unlist(daf.l1[1,to1]),pch=21,col="white",bg=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,9)])
arrows(c(1:4,6:7,9:10,12:13,15:18,20)+0.2,unlist(daf.o1[2,to1]),
       c(1:4,6:7,9:10,12:13,15:18,20)+0.2,unlist(daf.o1[3,to1]),length=0,angle=90,
       col=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,9)])
points(c(1:4,6:7,9:10,12:13,15:18,20)+0.2,pch=21,
       unlist(daf.o1[1,to1]),col=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,9)])
abline(v=c(5,8,11,14),lty=2)
axis(1,c(1:4,6:7,9:10,12:13,15:18,20),label=colnames(daf.h1)[to1],lwd=0,las=3)
legend("topleft",lty=c(1,2,1),pch=c(19,19,1),legend=c("highAbun","lowABun","others"))
plot(NA,xlim=c(0,23),ylim=c(0.02,0.12),
     xlab="",ylab="DAF",xaxt="n",main="derived allele frequency for 17-25nt seed")
arrows(c(1:4,6:7,9:10,12:13,15:18,20)-0.2,unlist(daf.h2[2,to1]),
       c(1:4,6:7,9:10,12:13,15:18,20)-0.2,unlist(daf.h2[3,to1]),length=0,angle=90,
       col=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,9)])
points(c(1:4,6:7,9:10,12:13,15:18,20)-0.2,
       unlist(daf.h2[1,to1]),pch=21,col="white",bg=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,9)])
arrows(c(1:4,6:7,9:10,12:13,15:18,20),unlist(daf.l2[2,to1]),
       c(1:4,6:7,9:10,12:13,15:18,20),unlist(daf.l2[3,to1]),length=0,angle=90,
       col=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,9)],lty=2)
points(c(1:4,6:7,9:10,12:13,15:18,20),
       unlist(daf.l2[1,to1]),pch=21,col="white",bg=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,9)])
arrows(c(1:4,6:7,9:10,12:13,15:18,20)+0.2,unlist(daf.o2[2,to1]),
       c(1:4,6:7,9:10,12:13,15:18,20)+0.2,unlist(daf.o2[3,to1]),length=0,angle=90,
       col=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,9)])
points(c(1:4,6:7,9:10,12:13,15:18,20)+0.2,pch=21,
       unlist(daf.o2[1,to1]),col=csl[c(5,5,5,5,7,7,3,3,1,1,1,1,1,1,9)])
abline(v=c(5,8,11,14),lty=2)
axis(1,c(1:4,6:7,9:10,12:13,15:18,20),label=colnames(daf.h1)[to1],lwd=0,las=3)
legend("topleft",lty=c(1,2,1),pch=c(19,19,1),legend=c("highAbun","lowABun","others"))
dev.off()


# Oct24.2018: GTEx----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
df_rna=read.table("merged.rnaseq.gene+pi.rpkm",header=T,row.names=1,check.names=F)
df_rna=df_rna[,rna_rn]
df_gtex=read.table("GTEx.testis.tpm.txt",header=T,row.names=1,check.names=F)
gn=intersect(row.names(df_rna),row.names(df_gtex))
gn=intersect(gn,human_pc)
tdf=cbind(df_rna[gn,],df_gtex[gn,])
tdf1=as.matrix(tdf)
tdf1=normalize.quantiles(tdf1)
row.names(tdf1)=row.names(tdf);colnames(tdf1)=colnames(tdf)
pc=prcomp(tdf1)
pdf("../Bear_Jul12.2018/GTEx.PCA.pdf",width=5,height=5,useDingbats = F)
par(mar=c(4,4,1,1),bty="n",tcl=0.3,cex=5/6)
plot(pc[[2]][,1:2],
     bg=c(rep(csl[5],4),rep(csl[3],5),rep(csl[2],3),rep(csl[1],6),rep("black",184)),
     pch=21,col="white",lwd=0.2)
legend("top",bty="n",pch=19,col=c(csl[c(5,3,2,1)],"black"),
       legend=c("Juv","26nt","inter","30nt","GTEx samples"))
dev.off()

c=cor(tdf,method="spearman")
pheatmap(c,cluster_rows=T,cluster_cols=T,cellwidth=4,cellheight=4,
         fontsize=4,filename="../Bear_Jul12.2018/GTEx.heatmap.cor.pdf")

plot(log10(tdf1[,c(11,19)]+0.1),lwd=0.2)
points(log10(tdf1[pl,c(11,19)]+0.1),pch=21,col="white",bg=csl[1])
abline(0,1)
plot(log10(tdf1[,c(11,1)]+0.1),lwd=0.2)
points(log10(tdf1[pl,c(11,5)]+0.1),pch=21,col="white",bg=csl[1])
abline(0,1)
plot(hclust(dist(t(tdf1[pl,]))),cex=3/6)


plot(unlist(tdf["MYBL1",]),pch=21,col="white",
     bg=c(rep(csl[5],4),rep(csl[3],5),rep(csl[2],3),rep(csl[1],6),rep("black",184)))
plot(unlist(tdf["TDRD12",]),pch=21,col="white",
     bg=c(rep(csl[5],4),rep(csl[3],5),rep(csl[2],3),rep(csl[1],6),rep("black",184)))



# Nov28.2018: AMYB across human and mouse----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/Bear_Jul12.2018/amyb_gene_list/")
ths=as.vector(read.table("compare_with_mouse/human.to.mouse.gene.amyb.list",header=F,row.names=NULL)[,1])
tms=as.vector(read.table("compare_with_mouse/mouse.amyb.DEG.txt",header=F,row.names=NULL)[,1])
pdf("compare_with_mouse/venn.pdf",width=5,height=5,useDingbats = F)
t=venn.diagram(list(human=ths,mouse=tms),
               filename=NULL,fill=csl[1:2],main=paste("AMYB overlap with previous paper"),margin=0.1)
par(mar=c(0,0,0,0))
plot.new();grid.draw(t)
t=venn.diagram(list(human=ths,mouse=c(amyb_gene,amyb_linc,amyb_other)),
               filename=NULL,fill=csl[1:2],main=paste("AMYB overlap with Xin's paper"),margin=0.1)
par(mar=c(0,0,0,0))
plot.new();grid.draw(t)
dev.off()

write.table(intersect(ths,tms),"compare_with_mouse/human_schemanti.common.txt",row.names=F,
            col.names=F,sep="\t",quote=F)
write.table(intersect(ths,c(amyb_gene,amyb_linc,amyb_other)),
            "compare_with_mouse/human_Xin.common.txt",row.names=F,
            col.names=F,sep="\t",quote=F)



# Dec02.2018: how many piRNAs after filtering in picluster----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
t=as.vector(read.table("N8150-Ox.trimmed.x_rRNA.x_hairpin.hg19v2.all.25_31nt.uniq.cdf",
             header=F,row.names=NULL)[,1])
t[which(t>10)]=10
plot(hist(t,breaks=0:10),main="Histogram of piRNAs",
     xlab="piRNA abundance",ylab="piRNA species",xaxt="n")
axis(1,0:9+0.5,label=c(1:9,">=10"),lwd=0)





# Dec18.2018: exon intron piRNA density----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
tdf=read.table("N8150.piC.exon_intron.rpkm",header=F,row.names=1)
tn=row.names(tdf)[which(tdf[,2]==(-1))]
pdf("../Bear_Jul12.2018/exon_vs_intron.pirna_density.pdf",width=3.5,height=3,useDingbats = F)
par(mar=c(5,4,1,1),tcl=0.3,bty="n",cex=5/6)
boxplot(log10(tdf[setdiff(human_prepachy,tn),1]+1),log10(tdf[setdiff(human_prepachy,tn),2]+1),
        log10(tdf[setdiff(human_pachy,tn),1]+1),log10(tdf[setdiff(human_pachy,tn),2]+1),
        log10(tdf[intersect(human_prepachy,tn),1]+1),log10(tdf[intersect(human_pachy,tn),1]+1),
        pch=20,cex=0.3,lty=1,xaxt="n",taxt="",at=c(1,2,4,5,7,9),col=csl[c(1,2,1,2,1,1)],
        ylab="RPKM",outline=T,range=10000,staplewex=0)
legend("topleft",bty="n",legend=c("exon","intron"),fill=csl[1:2])
axis(1,at=c(1.5,4.5,7,9),label=c("prepachy\nwith intron","pachy\nwith intron",
     "prepachy\nwithout intron","pachy\nwithout intron"),las=3,lwd=0)
x=log10(tdf[setdiff(human_prepachy,tn),1]+1);ti=1
points(runif(length(x),ti-0.2,ti+0.2),x,pch=21,col="white",bg="black")
x=log10(tdf[setdiff(human_prepachy,tn),2]+1);ti=2
points(runif(length(x),ti-0.2,ti+0.2),x,pch=21,col="white",bg="black")
x=log10(tdf[setdiff(human_pachy,tn),1]+1);ti=4
points(runif(length(x),ti-0.2,ti+0.2),x,pch=21,col="white",bg="black")
x=log10(tdf[setdiff(human_pachy,tn),2]+1);ti=5
points(runif(length(x),ti-0.2,ti+0.2),x,pch=21,col="white",bg="black")
x=log10(tdf[intersect(human_prepachy,tn),1]+1);ti=7
points(runif(length(x),ti-0.2,ti+0.2),x,pch=21,col="white",bg="black")
x=log10(tdf[intersect(human_pachy,tn),1]+1);ti=9
points(runif(length(x),ti-0.2,ti+0.2),x,pch=21,col="white",bg="black")
dev.off()

# Dec18.2018: pirna and longrna correlation in prepachy 3utr----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
tdf=read.table("N8150.pi3UTR.rpkm",header=F,row.names=1)
pdf("../Bear_Jul12.2018/utr.correlation.pdf",width=4,height=4,useDingbats = F)
par(mar=c(4,4,3,1),tcl=0.3,bty="n",cex=5/6)
plot(log10(tdf[,1:2]+1),pch=21,col="white",bg="black",lwd=0.3,xlim=c(0,2.5),ylim=c(0,2.5),
     xaxt="n",yaxt="n",xlab="longRNA RPKM",ylab="piRNA RPKM",main="prepachy UTR")
add_axis(1);add_axis(2)
text(0,2.3,pos=4,label="spearman cor=0.345")
dev.off()

# Dec18.2018: pirna abundance in euthe, primate, noncon piclusters; based on July19 syntenic analysis----
pdf("../Bear_Jul12.2018/pirnaAbun.eutherian_primate_noncon_clusters.pdf",width=2.5,height=2.2,useDingbats = F)
par(mar=c(4,4,1,1),tcl=0.3,cex=5/6,bty="n")
boxplot(log10(df_srna[human_euthe,12]+1),log10(df_srna[human_primate,12]+1),
        log10(df_srna[human_non,12]+1),col=cs_gg[1:3],lty=1,ylab="RPM",xaxt="n")
axis(1,1:3,label=c("eutherian","primates","non-conserved"),las=3,lwd=0)
dev.off()
# Dec18.2018: amyb motif distance in other mammals----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
tl=list()
tl[["marmoset"]]=read.table("calJac3.myba.motif.distance",header=F,row.names=1)
tl[["rat"]]=read.table("rn6.myba.motif.distance",header=F,row.names=1)
tl[["opossum"]]=read.table("monDom5.myba.motif.distance",header=F,row.names=1)
tl[["platypus"]]=read.table("ornAna1.myba.motif.distance",header=F,row.names=1)
pdf("../Bear_Jul12.2018/amyb.motif.distance.pdf",width=2.2,height=2.5,useDingbats = F)
par(mar=c(5,4,3,1),tcl=0.3,cex=5/6,bty="n")
for(i in names(tl)){
  r1=row.names(tl[[i]])[grep("_IG_",row.names(tl[[i]]))]
  r2=row.names(tl[[i]])[which(tl[[i]][,1]=="genicPi")];r2=setdiff(r2,r1)
  r3=row.names(tl[[i]])[which(tl[[i]][,1]=="protein_coding")]
  boxplot(log10(tl[[i]][r1,2]+1),log10(tl[[i]][r2,2]+1),log10(tl[[i]][r3,2]+1),
          xaxt="n",ylab="AMYB motif distance",yaxt="n",main=i,col=cs_gg[1:3],
          pch=20,cex=0.3,lty=1)
  add_axis(2)
  axis(1,1:3,label=c("intergenic\npiC","genic\npiC","mRNA"),lwd=0,las=3)
}
dev.off()


# Dec19.2018: rmsk content for utr, cds and intron----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
t3utr=read.table("hg19.mRNA.3UTR.rmsk.content.txt",header=F,row.names=1)
t5utr=read.table("hg19.mRNA.5UTR.rmsk.content.txt",header=F,row.names=1)
tcds=read.table("hg19.mRNA.CDS.rmsk.content.txt",header=F,row.names=1)
tintron=read.table("hg19.mRNA.intron.rmsk.content.txt",header=F,row.names=1)
fun_t=function(d,n){
  t=d[intersect(row.names(d),n),]
  return(((t[,2]+t[,3])/t[1])[,1]*100)
}
pdf("../Bear_Jul12.2018/rmsk.content.genicPre.utrCdsIntron.boxplot.pdf",width=4,height=4)
par(mar=c(6,4,3,1),bty="n",tcl=0.3,cex=5/6)
boxplot(fun_t(t3utr,row.names(t3utr)),fun_t(t3utr,intersect(row.names(t3utr),human_prepachy_genic1)),
        fun_t(t5utr,row.names(t5utr)),fun_t(t5utr,intersect(row.names(t5utr),human_prepachy_genic1)),
        fun_t(tcds,row.names(tcds)),fun_t(tcds,intersect(row.names(tcds),human_prepachy_genic1)),
        fun_t(tintron,row.names(tintron)),fun_t(tintron,intersect(row.names(tintron),human_prepachy_genic1)),
        col=csl[c(7,2)],xaxt="n",at=c(1,2,4,5,7,8,10,11),
        ylab="percentage (%)",main="rmsk content in exons",outline=F,lty=1)
axis(1,c(1,4,7,10)+0.5,label=c("3'UTR","5'UTR","CDS","intron"),las=3,lwd=0)
abline(h=0.4472413*100,lty=2)
legend("topleft",legend=c("mRNA","genicPre"),fill=csl[c(7,2)],bty="n")
dev.off()

pdf("../Bear_Jul12.2018/rmsk.content.pdf",width=8,height=3.5)
par(mar=c(8,4,1,1),tcl=0.3,cex=1/2,bty="n",las=3,lwd=0.5)
n=c(human_pachy,human_hybrid,human_prepachy)[182:1]
ci=as.integer(unlist(piG.rmsk[n,2]/rowSums(piG.rmsk[n,2:3]+0.001)*100))+1
barplot(rowSums(piG.rmsk[n,2:3])/piG.rmsk[n,1]*100,
        col=colorRampPalette(brewer.pal(9,"RdYlBu"))(101)[ci],ylab="rmsk percentage",
        border="grey")
abline(h=0.4472413*100,lty=2)
legend("topleft",fill=brewer.pal(9,"RdYlBu"),legend=c("100%",rep("",7),"0%"),
       title="anti-sense",y.intersp = 0.5)
dev.off()


# Updated Figure4B with amyb genes----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/")
pd=read.table("matrix/human.amybPeakMerged.distance.txt",header=F,row.names=1)
human_amyb_pc=human_pc[which(pd[human_pc,1]<500)]
human_amyb_linc=human_pc[which(pd[human_linc,1]<500)]
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
df_srna=read.table("merged.rpm",header=T,row.names=1,check.names=F)
df_rna=read.table("merged.rnaseq.gene+pi.rpkm",header=T,row.names=1,check.names=F)
df_rna=df_rna[,rna_rn]
pdf("../Bear_Jul12.2018/boxplot.amybGene.pdf",width=4,height=4,useDingbats=F)
par(tcl=0.3,mar=c(6,4,1,1),bty="n")
boxplot(as.vector(unlist(df_rna[human_amyb_pc,1:4])),as.vector(unlist(df_rna[human_amyb_pc,5:7])),
        as.vector(unlist(df_rna[human_amyb_pc,8:12])),as.vector(unlist(df_rna[human_amyb_pc,13:18])),
        col=tcs[c(1,5,8,13)],cex=0.3,pch=20,ylab="RPKM",names=c("Juv","26nt","inter","30nt"),
        ylim=c(0,60))
dev.off()

# Jan11.2019: TCFL5 and NFYA peak distance in human and mouse----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/")
pd_mm_nfya=read.table("matrix/mouse_3MO-3-NFYA-ID1_S3.peak.distance.txt",header=F,row.names=1)
pd_mm_tcfl5=read.table("matrix/mouse_3MO-TCFL5-IP-ID02_S2.peak.distance.txt",header=F,row.names=1)
pd_hg_nfya1=read.table("matrix/P75-Hs-NFYA-ID3_S2.peak.distance.txt",header=F,row.names=1)
pd_hg_nfya2=read.table("matrix/s904-Hs-NFYA-ID5_S3.peak.distance.txt",header=F,row.names=1)
# function
fun_plot_mm=function(pd,m){
  plot(NA,xlim=c(0,11),ylim=c(0,7),ylab="distance (bp)",xlab="",xaxt="n",yaxt="n")
  ph=21
  boxplot(log10(pd[Pachy,1]+1),log10(pd[pathway_list1,1]+1),
          log10(pd[Prepachy,1]+1),log10(pd[Hybrid,1]+1),
          log10(pd[proteincoding,1]+1),log10(pd[lincRNA,1]+1),outline=F,
          at=c(1,3,5,7,9,11)-0.5,boxwex=1,staplewex=0,add=T,yaxt="n",xaxt="n",main=m)
  #points(runif(length(l1),0+0.05,1-0.05),rnalog10(l1+1),col=csl[1],pch=ph,lwd=1,cex=0.2)
  points(runif(length(Prepachy),4+0.05,5-0.05),log10(pd[Prepachy,1]+1),bg=csl[2],pch=ph,lwd=0.6,cex=0.7,col="white")
  points(runif(length(Pachy),0+0.05,1-0.05),log10(pd[Pachy,1]+1),bg=csl[1],pch=ph,lwd=0.6,cex=0.7,col="white")
  points(runif(length(pathway_list1),2+0.05,3-0.05),log10(pd[pathway_list1,1]+1),bg=csl[1],pch=ph,lwd=0.6,cex=0.7,col="white")
  points(runif(length(Hybrid),6+0.05,7-0.05),log10(pd[Hybrid,1]+1),bg=csl[3],pch=ph,lwd=0.6,cex=0.7,col="white")
  points(runif(length(proteincoding),8+0.05,9-0.05),log10(pd[proteincoding,1]+1),bg=csl[4],pch=ph,lwd=0.6,cex=0.7,col="white")
  points(runif(length(lincRNA),10+0.05,11-0.05),log10(pd[lincRNA,1]+1),bg="black",pch=ph,lwd=0.6,cex=0.7,col="white")
  text(c(1,3,5,7,9,11,13)-0.5,-1,label=c("pachytene\nclusters",
                                         "piRNA\nbiogenesis",
                                         "prepachytene\nclusters","hybrid\nclusters","mRNA genes",
                                         "lincRNA genes"),
       xpd=T,srt=45,col=c(csl[c(1,1,2,3,4)],"black","black"))
  axis(2,0:7,c(expression(10^0,10^1,10^2,10^3,10^4,10^5,10^6,10^7)))
}
fun_plot_hg=function(pd,m){
  plot(NA,xlim=c(0,13),ylim=c(0,7),ylab="distance (bp)",xlab="",xaxt="n",yaxt="n")
  ph=21
  boxplot(log10(pd[human_pachy_amyb,1]+1),log10(pd[human_pachy_noamyb,1]+1),
          log10(pd[pl,1]+1),log10(pd[human_prepachy,1]+1),
          log10(pd[human_hybrid,1]+1),log10(pd[human_pc,1]+1),log10(pd[human_linc,1]+1),outline=F,
          at=c(1,3,5,7,9,11,13)-0.5,boxwex=1,staplewex=0,add=T,yaxt="n",xaxt="n",main=m)
  #points(runif(length(l1),0+0.05,1-0.05),rnalog10(l1+1),col=csl[1],pch=ph,lwd=1,cex=0.2)
  points(runif(length(pl),4+0.05,5-0.05),log10(pd[pl,1]+1),bg=csl[2],pch=ph,lwd=0.6,cex=0.7,col="white")
  points(runif(length(human_pachy_amyb),0+0.05,1-0.05),log10(pd[human_pachy_amyb,1]+1),bg=csl[1],pch=ph,lwd=0.6,cex=0.7,col="white")
  points(runif(length(human_pachy_noamyb),2+0.05,3-0.05),log10(pd[human_pachy_noamyb,1]+1),bg=csl[1],pch=ph,lwd=0.6,cex=0.7,col="white")
  points(runif(length(human_prepachy),6+0.05,7-0.05),log10(pd[human_prepachy,1]+1),bg=csl[3],pch=ph,lwd=0.6,cex=0.7,col="white")
  points(runif(length(human_hybrid),8+0.05,9-0.05),log10(pd[human_hybrid,1]+1),bg=csl[4],pch=ph,lwd=0.6,cex=0.7,col="white")
  points(runif(length(human_pc),10+0.05,11-0.05),log10(pd[human_pc,1]+1),bg="black",pch=ph,lwd=0.6,cex=0.7,col="white")
  points(runif(length(human_linc),12+0.05,13-0.05),log10(pd[human_linc,1]+1),bg="black",pch=ph,lwd=0.6,cex=0.7,col="white")
  text(c(1,3,5,7,9,11,13)-0.5,-1,label=c("amyb-positive\npachytene\nclusters",
                                         "amyb-negative\npachytene\nclusters","piRNA\nbiogenesis",
                                         "prepachytene\nclusters","hybrid\nclusters","mRNA genes",
                                         "lincRNA genes"),
       xpd=T,srt=45,col=c(csl[c(1,1,2,3,4)],"black","black"))
  axis(2,0:7,c(expression(10^0,10^1,10^2,10^3,10^4,10^5,10^6,10^7)))
}
fun_plot_hg1=function(pd,m){
  plot(NA,xlim=c(0,11),ylim=c(0,7),ylab="distance (bp)",xlab="",xaxt="n",yaxt="n")
  ph=21
  boxplot(log10(pd[human_pachy,1]+1),
          log10(pd[pl,1]+1),log10(pd[human_prepachy,1]+1),
          log10(pd[human_hybrid,1]+1),log10(pd[human_pc,1]+1),log10(pd[human_linc,1]+1),outline=F,
          at=c(1,3,5,7,9,11)-0.5,boxwex=1,staplewex=0,add=T,yaxt="n",xaxt="n",main=m)
  #points(runif(length(l1),0+0.05,1-0.05),rnalog10(l1+1),col=csl[1],pch=ph,lwd=1,cex=0.2)
  points(runif(length(pl),2+0.05,3-0.05),log10(pd[pl,1]+1),bg=csl[2],pch=ph,lwd=0.6,cex=0.7,col="white")
  points(runif(length(human_pachy_amyb),0+0.05,1-0.05),log10(pd[human_pachy_amyb,1]+1),bg=csl[1],pch=ph,lwd=0.6,cex=0.7,col="white")
  points(runif(length(human_prepachy),4+0.05,5-0.05),log10(pd[human_prepachy,1]+1),bg=csl[3],pch=ph,lwd=0.6,cex=0.7,col="white")
  points(runif(length(human_hybrid),6+0.05,7-0.05),log10(pd[human_hybrid,1]+1),bg=csl[4],pch=ph,lwd=0.6,cex=0.7,col="white")
  points(runif(length(human_pc),8+0.05,9-0.05),log10(pd[human_pc,1]+1),bg="black",pch=ph,lwd=0.6,cex=0.7,col="white")
  points(runif(length(human_linc),10+0.05,11-0.05),log10(pd[human_linc,1]+1),bg="black",pch=ph,lwd=0.6,cex=0.7,col="white")
  text(c(1,3,5,7,9,11)-0.5,-1,label=c("pachytene\nclusters","piRNA\nbiogenesis",
                                         "prepachytene\nclusters","hybrid\nclusters","mRNA genes",
                                         "lincRNA genes"),
       xpd=T,srt=45,col=c(csl[c(1,1,2,3,4)],"black","black"))
  axis(2,0:7,c(expression(10^0,10^1,10^2,10^3,10^4,10^5,10^6,10^7)))
}

# AMYB
pdf("Bear_Jan11.2019.TCFL5/peak_distance.Ttcfl5_Nfya.pdf",width=4,height=4,useDingbats=F)
par(mar=c(4,4,3,1),cex=5/6,bty="n",xpd=T)
fun_plot_mm(pd_mm_nfya,"mouse NFYA")
fun_plot_mm(pd_mm_tcfl5,"mouse TCFL5")
fun_plot_hg(pd_hg_nfya1,"human P75 NFYA")
fun_plot_hg1(pd_hg_nfya1,"human P75 NFYA")
fun_plot_hg(pd_hg_nfya2,"human s904 NFYA")
fun_plot_hg1(pd_hg_nfya2,"human s904 NFYA")
dev.off()

# compare mouse and human AMYB signal----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/AMYB_peak_human_and_mouse/")
mouse_amyb_ortho=as.vector(read.table("mouse.amyb.list",header=F,row.names=NULL)[,1])
amyb_mouse=read.table("mouse.3MO.AMYB.uniq.rpm",header=F,row.names=1)
amyb_human=list()
for(i in list.files(pattern="^h.*m$")){
  amyb_human[[strsplit(i,".uniq")[[1]][1]]]=read.table(i,header=F,row.names=1)
}
for(i in names(amyb_human)){amyb_human[[i]][which(amyb_human[[i]][,1]>0.5),1]=0.5}
ortho=intersect(row.names(amyb_mouse),row.names(amyb_human[[1]]))
cons=read.table("human.mouse.con",header=F,row.names=NULL)
cons_mouse=read.table("human.mouse.con",header=F,row.names=2)
cons_human=read.table("human.mouse.con",header=F,row.names=1)
transname=read.table("human_mouse.genes.list",header=F,row.names=1)
mouse.amyb.motif=read.table("mouse.gene+pi.amybMotif.distance.txt",header=F,row.names=1)
mouse_amybM_ortho=row.names(mouse.amyb.motif)[which(mouse.amyb.motif[,1]<100)]
mouse_amybM_ortho1=intersect(mouse_amybM_ortho,ortho)
human.amyb.motif=read.table("human.gene+pi.amybMotif.distance.txt",header=F,row.names=1)
human_amybM_ortho=row.names(human.amyb.motif)[which(human.amyb.motif[,1]<100)]
human_amybM_ortho1=intersect(as.vector(transname[human_amybM_ortho,1]),ortho)
human_amyb_ortho=read.table("human.AMYB.union.distance.txt",header=F,row.names=1)
human_amyb_ortho=row.names(human_amyb_ortho)[which(human_amyb_ortho[,1]<500)]
human_amyb_ortho1=intersect(as.vector(transname[human_amyb_ortho,1]),ortho)
# scatterplot
pdf("../../Bear_Jul12.2018/amyb_mouse_vs_human.rpm.scatter.pdf",width=12,height=9,useDingbats=F)
par(mfrow=c(3,4),mar=c(4,4,1,1),tcl=0.3,bty="n",cex=5/6)
for(i in names(amyb_human)){
  sc=0.01
  tn=setdiff(ortho,c(human_amyb_ortho,mouse_amyb_ortho))
  xma=log10(max(amyb_mouse[ortho,1])+sc);xmi=log10(min(amyb_mouse[ortho,1])+sc)
  yma=log10(max(amyb_human[[i]][ortho,1])+sc);ymi=log10(min(amyb_human[[i]][ortho,1])+sc)
  x=log10(amyb_mouse[tn,1]+sc);y=log10(amyb_human[[i]][tn,1]+sc)
  plot(x,y,col=rgb(0,0,0,maxColorValue=255,alpha=50),xlab="log2FE; mouse 3MO testis",
       ylab=paste("log2FE; ",i,sep=""),pch=20,lwd=0.05,cex=0.7,
       xlim=c(xmi,xma),ylim=c(ymi,yma))
  tn=intersect(human_amyb_ortho,mouse_amyb_ortho);x1=log10(amyb_mouse[tn,1]+sc);y1=log10(amyb_human[[i]][tn,1]+sc)
  tn=setdiff(human_amyb_ortho,mouse_amyb_ortho);x2=log10(amyb_mouse[tn,1]+sc);y2=log10(amyb_human[[i]][tn,1]+sc)
  tn=setdiff(mouse_amyb_ortho,human_amyb_ortho);x3=log10(amyb_mouse[tn,1]+sc);y3=log10(amyb_human[[i]][tn,1]+sc)
  points(x3,y3,col="white",bg=rgb(0,0,255,maxColorValue=255,alpha=100),pch=21,lwd=0.1)
  points(x1,y1,col="white",bg=rgb(255,0,0,maxColorValue=255,alpha=100),pch=21,lwd=0.1)
  points(x2,y2,col="white",bg=rgb(0,255,0,maxColorValue=255,alpha=100),pch=21,lwd=0.1)
  c=round(cor(x,y,method="spearman"),2);c1=round(cor(x1,y1,method="spearman"),2)
  c2=round(cor(x2,y2,method="spearman"),2);c3=round(cor(x3,y3,method="spearman"),2)
  text(xmi,ymi+(yma-ymi)*9/10,pos=4,label=paste("no amyb",c),col="grey")
  text(xmi,ymi+(yma-ymi)*8/10,pos=4,label=paste("mouse+human",c1),col="red")
  text(xmi,ymi+(yma-ymi)*7/10,pos=4,label=paste("human",c2),col="green")
  text(xmi,ymi+(yma-ymi)*6/10,pos=4,label=paste("mouse",c3),col="blue")
}
dev.off()
# boxplot for mouse----
conspachy_human=intersect(as.vector(cons[,1]),human_pachy)
conspachy_human_amyb=intersect(conspachy_human,human_amyb_ortho)
conspachy_human_noamyb=setdiff(conspachy_human,human_amyb_ortho)
conspachy_mouse=intersect(as.vector(cons[,2]),Pachy)
pathway_ortho=intersect(pathway_list1,ortho)
pdf("../../Bear_Jul12.2018/amyb.boxplot.detailedGroup.pdf",width=8,height=6,useDingbats=F)
par(mar=c(15,4,3,1),tcl=0.3,cex=5/6,bty="n")
pathway_ortho=c("Piwil1","Piwil2","Piwil4","Mybl1","Tdrd1","Tdrd3","Tdrd5",
                "Tdrd6","Tdrd9","Tdrd12","Ddx4","Ddx39","Fkbp6","Gpat2",
                "Henmt1","Mael","Mov10l1","Nfya","Pld6")
human_amyb_ortho1=c(human_amyb_ortho1,c("Mybl1"))
mouse_amyb_ortho=c(mouse_amyb_ortho,c("Piwil1","Pld6"))
# boxplot1 mouse amyb peak----
bxl=list()
bxl[["conserved pachy with\namyb in human and mouse"]]=log10(amyb_mouse[as.vector(cons_human[conspachy_human_amyb,1]),1]+sc)
bxl[["conserved pachy with\namyb only in human"]]=log10(amyb_mouse[as.vector(cons_human[conspachy_human_noamyb,1]),1]+sc)
bxl[["nonconserved pachy"]]=log10(amyb_mouse[setdiff(Pachy,conspachy_mouse),1]+sc)
bxl[["hybrid"]]=log10(amyb_mouse[Hybrid,1]+sc)
bxl[["prepachy"]]=log10(amyb_mouse[Prepachy,1]+sc)
bxl[["conserved gene with\n no amyb in mouse and human"]]=log10(amyb_mouse[setdiff(ortho,c(human_amyb_ortho1,mouse_amyb_ortho)),1]+sc)
bxl[["conserved gene with\namyb only in human"]]=log10(amyb_mouse[setdiff(human_amyb_ortho1,mouse_amyb_ortho),1]+sc)
bxl[["conserved gene with\namyb only in mouse"]]=log10(amyb_mouse[setdiff(mouse_amyb_ortho,human_amyb_ortho1),1]+sc)
bxl[["conserved gene with\namyb in human and mouse"]]=log10(amyb_mouse[intersect(human_amyb_ortho1,mouse_amyb_ortho),1]+sc)
bxl[["conserved pathway with\n no amyb in mouse and human"]]=log10(amyb_mouse[setdiff(pathway_ortho,c(human_amyb_ortho1,mouse_amyb_ortho)),1]+sc)
bxl[["conserved pathway with\namyb only in human"]]=log10(amyb_mouse[setdiff(intersect(human_amyb_ortho1,pathway_ortho),mouse_amyb_ortho),1]+sc)
bxl[["conserved pathway with\namyb only in mouse"]]=log10(amyb_mouse[setdiff(intersect(mouse_amyb_ortho,pathway_ortho),human_amyb_ortho1),1]+sc)
bxl[["conserved pathway with\namyb in human and mouse"]]=log10(amyb_mouse[intersect(intersect(human_amyb_ortho1,pathway_ortho),mouse_amyb_ortho),1]+sc)
tcs=csl[c(1,1,1,9,9,2,2,2,2,3,3,3,3)]
yma=-100;ymi=100000
tl=c()
for(i in names(bxl)){yma=max(yma,max(bxl[[i]]));ymi=min(ymi,min(bxl[[i]]));tl=c(tl,length(bxl[[i]]))}
boxplot(bxl,names=names(bxl),las=3,pch=20,cex=0.3,staplewex=0,lty=1,boxwex=0.5,
        ylim=c(ymi,yma+(yma-ymi)/8),col=tcs,ylab="mouse AMYB RPM",
        main="amyb stands for amyb peak")
text(1:length(bxl),rep(yma,length(bxl)),srt=45,format(tl,big.mark=",",trim=T),pos=3)
# boxplot2 mouse amyb motif----
human_amybM_ortho1=c(human_amybM_ortho1,"Pld6")
mouse_amybM_ortho1=c(mouse_amybM_ortho1,"Pld6")
bxl=list()
bxl[["conserved pachy with\nno amyb in mouse and human"]]=log10(amyb_mouse[setdiff(conspachy_mouse,c(mouse_amybM_ortho,as.vector(cons_human[human_amybM_ortho,1]))),1]+sc)
bxl[["conserved pachy with\namyb only in human"]]=log10(amyb_mouse[setdiff(intersect(conspachy_mouse,as.vector(cons_human[human_amybM_ortho,1])),mouse_amybM_ortho),1]+sc)
bxl[["conserved pachy with\namyb only in mouse"]]=log10(amyb_mouse[setdiff(intersect(conspachy_mouse,mouse_amybM_ortho),as.vector(cons_human[human_amybM_ortho,1])),1]+sc)
bxl[["conserved pachy with\namyb in human and mouse"]]=log10(amyb_mouse[intersect(c(conspachy_mouse,mouse_amybM_ortho),as.vector(cons_human[human_amybM_ortho,1])),1]+sc)
bxl[["nonconserved pachy with\nno amyb in mouse"]]=log10(amyb_mouse[setdiff(setdiff(Pachy,conspachy_mouse),mouse_amybM_ortho),1]+sc)
bxl[["nonconserved pachy with\namyb in mouse"]]=log10(amyb_mouse[intersect(setdiff(Pachy,conspachy_mouse),mouse_amybM_ortho),1]+sc)
bxl[["hybrid"]]=log10(amyb_mouse[Hybrid,1]+sc)
bxl[["prepachy"]]=log10(amyb_mouse[Prepachy,1]+sc)
bxl[["conserved gene with\n no amyb in mouse and human"]]=log10(amyb_mouse[setdiff(ortho,c(human_amybM_ortho1,mouse_amybM_ortho1)),1]+sc)
bxl[["conserved gene with\namyb only in human"]]=log10(amyb_mouse[setdiff(human_amybM_ortho1,mouse_amybM_ortho1),1]+sc)
bxl[["conserved gene with\namyb only in mouse"]]=log10(amyb_mouse[setdiff(mouse_amybM_ortho1,human_amybM_ortho1),1]+sc)
bxl[["conserved gene with\namyb in human and mouse"]]=log10(amyb_mouse[intersect(human_amybM_ortho1,mouse_amybM_ortho1),1]+sc)
bxl[["conserved pathway with\n no amyb in mouse and human"]]=log10(amyb_mouse[setdiff(pathway_ortho,c(human_amybM_ortho1,mouse_amybM_ortho1)),1]+sc)
bxl[["conserved pathway with\namyb only in human"]]=log10(amyb_mouse[setdiff(intersect(human_amybM_ortho1,pathway_ortho),mouse_amybM_ortho1),1]+sc)
bxl[["conserved pathway with\namyb only in mouse"]]=log10(amyb_mouse[setdiff(intersect(mouse_amybM_ortho1,pathway_ortho),human_amybM_ortho1),1]+sc)
bxl[["conserved pathway with\namyb in human and mouse"]]=log10(amyb_mouse[intersect(intersect(human_amybM_ortho1,pathway_ortho),mouse_amybM_ortho1),1]+sc)
tcs=csl[c(1,1,1,1,1,1,9,9,2,2,2,2,3,3,3,3)]
yma=-100;ymi=100000
tl=c()
for(i in names(bxl)){yma=max(yma,max(bxl[[i]]));ymi=min(ymi,min(bxl[[i]]));tl=c(tl,length(bxl[[i]]))}
boxplot(bxl,names=names(bxl),las=3,pch=20,cex=0.3,staplewex=0,lty=1,boxwex=0.5,
        ylim=c(ymi,yma+(yma-ymi)/8),col=tcs,ylab="mouse AMYB RPM",
        main="amyb stands for amyb motif with score higher than 7")
text(1:length(bxl),rep(yma,length(bxl)),srt=45,format(tl,big.mark=",",trim=T),pos=3)
# boxplot3 human amyb peak----
bxl=list()
tdf=amyb_human[["human.7358-AMYB-Id7_S1"]]
bxl[["conserved pachy with\namyb in human and mouse"]]=log10(tdf[conspachy_human_amyb,1]+sc)
bxl[["conserved pachy with\namyb in only in mouse"]]=log10(tdf[conspachy_human_noamyb,1]+sc)
bxl[["nonconserved pachy with\nno amyb in human"]]=log10(tdf[setdiff(setdiff(human_pachy,conspachy_human),human_amyb_ortho),1]+sc)
bxl[["nonconserved pachy with\namyb in human"]]=log10(tdf[intersect(setdiff(human_pachy,conspachy_human),human_amyb_ortho),1]+sc)
bxl[["hybrid"]]=log10(tdf[human_hybrid,1]+sc)
bxl[["prepachy"]]=log10(tdf[human_prepachy,1]+sc)
bxl[["conserved gene with\n no amyb in mouse and human"]]=log10(tdf[setdiff(ortho,c(human_amyb_ortho1,mouse_amyb_ortho)),1]+sc)
bxl[["conserved gene with\namyb only in mouse"]]=log10(tdf[setdiff(mouse_amyb_ortho,human_amyb_ortho1),1]+sc)
bxl[["conserved gene with\namyb only in human"]]=log10(tdf[setdiff(human_amyb_ortho1,mouse_amyb_ortho),1]+sc)
bxl[["conserved gene with\namyb in human and mouse"]]=log10(tdf[intersect(human_amyb_ortho1,mouse_amyb_ortho),1]+sc)
bxl[["conserved pathway with\n no amyb in mouse and human"]]=log10(tdf[setdiff(pathway_ortho,c(human_amyb_ortho1,mouse_amyb_ortho)),1]+sc)
bxl[["conserved pathway with\namyb only in mouse"]]=log10(tdf[setdiff(intersect(mouse_amyb_ortho,pathway_ortho),human_amyb_ortho1),1]+sc)
bxl[["conserved pathway with\namyb only in human"]]=log10(tdf[setdiff(intersect(human_amyb_ortho1,pathway_ortho),mouse_amyb_ortho),1]+sc)
bxl[["conserved pathway with\namyb in human and mouse"]]=log10(tdf[intersect(intersect(human_amyb_ortho1,pathway_ortho),mouse_amyb_ortho),1]+sc)
tcs=csl[c(1,1,1,1,9,9,2,2,2,2,3,3,3,3)]
yma=-100;ymi=100000
tl=c()
for(i in names(bxl)){yma=max(yma,max(bxl[[i]]));ymi=min(ymi,min(bxl[[i]]));tl=c(tl,length(bxl[[i]]))}
boxplot(bxl,names=names(bxl),las=3,pch=20,cex=0.3,staplewex=0,lty=1,boxwex=0.5,
        ylim=c(ymi,yma+(yma-ymi)/8),col=tcs,ylab="human AMYB RPM",
        main="amyb stands for amyb peak")
text(1:length(bxl),rep(yma,length(bxl)),srt=45,format(tl,big.mark=",",trim=T),pos=3)
# boxplot4 human amyb motif----
bxl=list()
tdf=amyb_human[["human.7358-AMYB-Id7_S1"]]
bxl[["conserved pachy with\nno amyb in mouse and human"]]=log10(tdf[setdiff(conspachy_human,c(human_amybM_ortho,as.vector(cons_mouse[mouse_amybM_ortho,1]))),1]+sc)
bxl[["conserved pachy with\namyb only in human"]]=log10(tdf[setdiff(intersect(conspachy_human,human_amybM_ortho),as.vector(cons_mouse[mouse_amybM_ortho,1])),1]+sc)
bxl[["conserved pachy with\namyb only in mouse"]]=log10(tdf[setdiff(intersect(conspachy_human,as.vector(cons_mouse[mouse_amybM_ortho,1])),human_amybM_ortho),1]+sc)
bxl[["conserved pachy with\namyb in human and mouse"]]=log10(tdf[intersect(c(conspachy_human,as.vector(cons_mouse[mouse_amybM_ortho,1])),human_amybM_ortho),1]+sc)
bxl[["nonconserved pachy with\nno amyb in human"]]=log10(tdf[setdiff(setdiff(human_pachy,conspachy_human),human_amybM_ortho),1]+sc)
bxl[["nonconserved pachy with\namyb in human"]]=log10(tdf[intersect(setdiff(human_pachy,conspachy_human),human_amybM_ortho),1]+sc)
bxl[["hybrid"]]=log10(tdf[human_pachy,1]+sc)
bxl[["prepachy"]]=log10(tdf[human_prepachy,1]+sc)
bxl[["conserved gene with\n no amyb in mouse and human"]]=log10(tdf[setdiff(ortho,c(human_amybM_ortho1,mouse_amybM_ortho1)),1]+sc)
bxl[["conserved gene with\namyb only in mouse"]]=log10(tdf[setdiff(mouse_amybM_ortho1,human_amybM_ortho1),1]+sc)
bxl[["conserved gene with\namyb only in human"]]=log10(tdf[setdiff(human_amybM_ortho1,mouse_amybM_ortho1),1]+sc)
bxl[["conserved gene with\namyb in human and mouse"]]=log10(tdf[intersect(human_amybM_ortho1,mouse_amybM_ortho1),1]+sc)
bxl[["conserved pathway with\n no amyb in mouse and human"]]=log10(tdf[setdiff(pathway_ortho,c(human_amybM_ortho1,mouse_amybM_ortho1)),1]+sc)
bxl[["conserved pathway with\namyb only in mouse"]]=log10(tdf[setdiff(intersect(mouse_amybM_ortho1,pathway_ortho),human_amybM_ortho1),1]+sc)
bxl[["conserved pathway with\namyb only in human"]]=log10(tdf[setdiff(intersect(human_amybM_ortho1,pathway_ortho),mouse_amybM_ortho1),1]+sc)
bxl[["conserved pathway with\namyb in human and mouse"]]=log10(tdf[intersect(intersect(human_amybM_ortho1,pathway_ortho),mouse_amybM_ortho1),1]+sc)
tcs=csl[c(1,1,1,1,1,1,9,9,2,2,2,2,3,3,3,3)]
yma=-100;ymi=100000
tl=c()
for(i in names(bxl)){yma=max(yma,max(bxl[[i]]));ymi=min(ymi,min(bxl[[i]]));tl=c(tl,length(bxl[[i]]))}
boxplot(bxl,names=names(bxl),las=3,pch=20,cex=0.3,staplewex=0,lty=1,boxwex=0.5,
        ylim=c(ymi,yma+(yma-ymi)/8),col=tcs,ylab="human AMYB RPM",
        main="amyb stands for amyb motif with score higher than 7")
text(1:length(bxl),rep(yma,length(bxl)),srt=45,format(tl,big.mark=",",trim=T),pos=3)
dev.off()

# compare mouse and human AMYB motif----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/AMYB_peak_human_and_mouse/")
pdf("../../Bear_Jul12.2018/amyb_mouse_vs_human.motif.boxplot.pdf",width=4,height=4,useDingbats = F)
par(mar=c(7,4,3,1),tcl=0.3,bty="n",cex=5/6,las=3)
mouse.amyb.motif=read.table("mouse.gene+pi.amybMotif.distance.txt",header=F,row.names=1)
human.amyb.motif=read.table("human.gene+pi.amybMotif.distance.txt",header=F,row.names=1)
boxplot(log10(human.amyb.motif[human_bid,1]+1),log10(mouse.amyb.motif[Bid,1]+1),
        log10(human.amyb.motif[setdiff(human_pachy,human_bid),1]+1),log10(mouse.amyb.motif[setdiff(Pachy,Bid),1]+1),
        log10(human.amyb.motif[human_prepachy,1]+1),log10(mouse.amyb.motif[Prepachy,1]+1),
        log10(human.amyb.motif[pl,1]+1),log10(mouse.amyb.motif[pathway_list1,1]+1),
        log10(human.amyb.motif[human_pc,1]+1),log10(mouse.amyb.motif[proteincoding,1]+1),
        col=csl[1:2],pch=20,yaxt="n",xaxt="n",at=c(1,2,4,5,7,8,10,11,13,14),ylim=c(0,5),
        main="AMYB motifs with motif score higher than 7",ylab="distance (bp)",cex=0.3)
text(c(1,4,7,10,13)+0.5,rep(4.4,5),pos=3,label=c("0.015","0.223","0.176","0.915","2.2e-16"),
     col=c("red","black","black","red"))
axis(1,c(1,4,7,10,13)+0.5,label=c("bidirectional\nPachy","unidirectional\nPachy","Prepachy","pathway","mRNA"))
add_axis(2)
dev.off()

plot(log10(human.amyb.motif[pl,1]+1),log10(mouse.amyb.motif[pathway_list_name,1]+1))



# new population variation analysis----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/population_new/")
fun_plot=function(df,mn){
  par(mar=c(12,4,3,1),bty="n",tcl=0.3,cex=5/6)
  ymi=min(df[1:3,]);yma=max(df[1:3,])
  if(nrow(df)==4){tyma=yma+(yma-ymi)/5;yl="DAF"}else{tyma=yma+(yma-ymi)/3;yl="π"}
  plot(NA,xlim=c(0,ncol(df)*10/9),ylim=c(ymi,tyma),xaxt="n",ylab=yl,xlab="",
       main=mn)
  points(unlist(df[1,]))
  for(i in 1:ncol(df)){lines(c(i,i),c(df[2,i],df[3,i]))}
  axis(1,1:ncol(df),labels=colnames(df),las=3,lwd=0)
  if(nrow(df)==4){
    text(1:ncol(df),yma,label=format(as.vector(df[4,]),big.mark=","),srt=45,pos=3,cex=0.7)
    text(ncol(df)*20/19,yma,pos=3,label="# of SNP",srt=45,cex=0.7)
  }else{
    text(1:ncol(df),yma+(yma-ymi)/10,srt=45,pos=3,
         label=paste(format(as.vector(df[4,]),big.mark=","),format(as.vector(df[5,]),big.mark=","),sep="/"),cex=0.7)
    text(ncol(df)*20/19,yma+(yma-ymi)/10,pos=3,label="# of SNP/# of nucleotide",srt=45,cex=0.7)
  }
}
# all DAF and pi including RMSK----
daf.all=list()
pi.all=list()
pop=c("All","AFR","AMR","EUR","SAS","EAS")
for(i in pop){
  daf.all[[i]]=read.table(paste("final_matrix/",i,".all.DAF.mat",sep=""),header=T,row.names=1)
  pi.all[[i]]=read.table(paste("final_matrix/",i,".all.pi.mat",sep=""),header=T,row.names=1)
}
pdf("all.pdf",width=7,height=5,useDingbats = F)
tcn=c(21,22,3,4,25,26,19,20,7,6,5,8,15,14,13,16,11,12,27,1,2,17,18)
fun_plot(daf.all$All[,tcn],"DAF in all population")
fun_plot(daf.all$AFR[,tcn],"DAF in Africa")
fun_plot(daf.all$AMR[,tcn],"DAF in America")
fun_plot(daf.all$EUR[,tcn],"DAF in European")
fun_plot(daf.all$SAS[,tcn],"DAF in South Asia")
fun_plot(daf.all$EAS[,tcn],"DAF in East Asia")
tcn=c(21,22,3,4,25,26,19,20,7,6,5,8,15,14,13,16,11,12,27,1,2,17,18)
fun_plot(pi.all$All[,tcn],"pi in all population")
fun_plot(pi.all$AFR[,tcn],"pi in Africa")
fun_plot(pi.all$AMR[,tcn],"pi in America")
fun_plot(pi.all$EUR[,tcn],"pi in European")
fun_plot(pi.all$SAS[,tcn],"pi in South Asia")
fun_plot(pi.all$EAS[,tcn],"pi in East Asia")
dev.off()

# all DAF and pi excluding RMSK----
daf.all=list()
pi.all=list()
pop=c("All","AFR","AMR","EUR","SAS","EAS")
for(i in pop){
  daf.all[[i]]=read.table(paste("final_matrix_noRMSK/",i,".all.DAF.mat",sep=""),header=T,row.names=1)
  pi.all[[i]]=read.table(paste("final_matrix_noRMSK/",i,".all.pi.mat",sep=""),header=T,row.names=1)
}
pdf("all.noRMSK.pdf",width=7,height=5,useDingbats = F)
tcn=c(21,22,3,4,25,26,19,20,7,6,5,8,15,14,13,16,11,12,27,1,2,17,18)
fun_plot(daf.all$All[,tcn],"DAF in all population")
fun_plot(daf.all$AFR[,tcn],"DAF in Africa")
fun_plot(daf.all$AMR[,tcn],"DAF in America")
fun_plot(daf.all$EUR[,tcn],"DAF in European")
fun_plot(daf.all$SAS[,tcn],"DAF in South Asia")
fun_plot(daf.all$EAS[,tcn],"DAF in East Asia")
tcn=c(21,22,3,4,25,26,19,20,7,6,5,8,15,14,13,16,11,12,27,1,2,17,18)
fun_plot(pi.all$All[,tcn],"pi in all population")
fun_plot(pi.all$AFR[,tcn],"pi in Africa")
fun_plot(pi.all$AMR[,tcn],"pi in America")
fun_plot(pi.all$EUR[,tcn],"pi in European")
fun_plot(pi.all$SAS[,tcn],"pi in South Asia")
fun_plot(pi.all$EAS[,tcn],"pi in East Asia")
dev.off()


# DAF in piRNA seed----
fun_plot=function(df1,df2,df3,mn){
  par(mar=c(12,4,3,1),bty="n",tcl=0.3,cex=5/6)
  ymi=min(min(df1[1:3,]),min(df2[1:3,]),min(df3[1:3,]))
  yma=max(max(df1[1:3,]),max(df2[1:3,]),max(df3[1:3,]))
  if(nrow(df1)==4){tyma=yma+(yma-ymi)/5;yl="DAF"}else{tyma=yma+(yma-ymi)/3;yl="π"}
  plot(NA,xlim=c(0,(3*ncol(df1))*10/9),ylim=c(ymi,tyma),xaxt="n",ylab=yl,xlab="",
       main=mn)
  points(1:ncol(df1)*3-2,unlist(df1[1,]),col=cs_gg[1])
  points(1:ncol(df1)*3-1,unlist(df2[1,]),col=cs_gg[2])
  points(1:ncol(df1)*3,unlist(df3[1,]),col=cs_gg[3])
  for(i in 1:ncol(df1)){lines(c(i*3-2,i*3-2),c(df1[2,i],df1[3,i]),col=cs_gg[1])}
  for(i in 1:ncol(df1)){lines(c(i*3-1,i*3-1),c(df2[2,i],df2[3,i]),col=cs_gg[2])}
  for(i in 1:ncol(df1)){lines(c(i*3,i*3),c(df3[2,i],df3[3,i]),col=cs_gg[3])}
  axis(1,1:ncol(df1)*3-1,labels=colnames(df1),las=3,lwd=0)
  if(nrow(df1)==4){
    text(1:ncol(df1)*3-2,yma,label=format(as.vector(df1[4,]),big.mark=","),srt=45,pos=3,cex=0.7,col=cs_gg[1])
    text(1:ncol(df1)*3-1,yma,label=format(as.vector(df2[4,]),big.mark=","),srt=45,pos=3,cex=0.7,col=cs_gg[2])
    text(1:ncol(df1)*3,yma,label=format(as.vector(df3[4,]),big.mark=","),srt=45,pos=3,cex=0.7,col=cs_gg[3])
    text(ncol(df1)*3*20/19,yma,pos=3,label="# of SNP",srt=45,cex=0.7)
  }else{
    text(1:ncol(df1)*3-2,yma+(yma-ymi)/10,srt=45,pos=3,col=cs_gg[1],
         label=paste(format(as.vector(df1[4,]),big.mark=","),format(as.vector(df1[5,]),big.mark=","),sep="/"),cex=0.7)
    text(1:ncol(df1)*3-1,yma+(yma-ymi)/10,srt=45,pos=3,col=cs_gg[2],
         label=paste(format(as.vector(df2[4,]),big.mark=","),format(as.vector(df2[5,]),big.mark=","),sep="/"),cex=0.7)
    text(1:ncol(df1)*3,yma+(yma-ymi)/10,srt=45,pos=3,col=cs_gg[3],
         label=paste(format(as.vector(df3[4,]),big.mark=","),format(as.vector(df3[5,]),big.mark=","),sep="/"),cex=0.7)
    text(ncol(df1)*3*20/19,yma+(yma-ymi)/10,pos=3,label="# of SNP/# of nucleotide",srt=45,cex=0.7)
  }
  legend("bottomright",pch=1,lty=1,col=cs_gg[1:3],
         legend=c("high piRNA bundance","low piRNA abundance","no piRNA"),bty="n")
}
# DAF and pi in high piA, low piA and no piA seed1(2-16nt) including RMSK----
daf.hip=list();daf.lop=list();daf.nop=list()
pi.hip=list();pi.lop=list();pi.nop=list()
pop=c("All","AFR","AMR","EUR","SAS","EAS")
for(i in pop){
  daf.hip[[i]]=read.table(paste("final_matrix/",i,".2-16.high.DAF.mat",sep=""),header=T,row.names=1)
  daf.lop[[i]]=read.table(paste("final_matrix/",i,".2-16.low.DAF.mat",sep=""),header=T,row.names=1)
  daf.nop[[i]]=read.table(paste("final_matrix/",i,".2-16.nopi.DAF.mat",sep=""),header=T,row.names=1)
  pi.hip[[i]]=read.table(paste("final_matrix/",i,".2-16.high.pi.mat",sep=""),header=T,row.names=1)
  pi.lop[[i]]=read.table(paste("final_matrix/",i,".2-16.low.pi.mat",sep=""),header=T,row.names=1)
  pi.nop[[i]]=read.table(paste("final_matrix/",i,".2-16.nopi.pi.mat",sep=""),header=T,row.names=1)
}
pdf("seed2-16.pdf",width=14,height=5,useDingbats = F)
tcn=c(21,22,3,4,25,26,19,20,7,6,5,8,15,14,13,16,11,12,27,1,2,17,18)
fun_plot(daf.hip$All[,tcn],daf.lop$All[,tcn],daf.nop$All[,tcn],"DAF in all population for 2-16nt seed")
fun_plot(daf.hip$AFR[,tcn],daf.lop$AFR[,tcn],daf.nop$AFR[,tcn],"DAF in Africa population for 2-16nt seed")
fun_plot(daf.hip$AMR[,tcn],daf.lop$AMR[,tcn],daf.nop$AMR[,tcn],"DAF in America population for 2-16nt seed")
fun_plot(daf.hip$EUR[,tcn],daf.lop$EUR[,tcn],daf.nop$EUR[,tcn],"DAF in European population for 2-16nt seed")
fun_plot(daf.hip$SAS[,tcn],daf.lop$SAS[,tcn],daf.nop$SAS[,tcn],"DAF in South Asia population for 2-16nt seed")
fun_plot(daf.hip$EAS[,tcn],daf.lop$EAS[,tcn],daf.nop$EAS[,tcn],"DAF in East Asia population for 2-16nt seed")
tcn=c(21,22,3,4,25,26,19,20,7,6,5,8,15,14,13,16,11,12,27,1,2,17,18)
fun_plot(pi.hip$All[,tcn],pi.lop$All[,tcn],pi.nop$All[,tcn],"pi in all population for 2-16nt seed")
fun_plot(pi.hip$AFR[,tcn],pi.lop$AFR[,tcn],pi.nop$AFR[,tcn],"pi in Africa population for 2-16nt seed")
fun_plot(pi.hip$AMR[,tcn],pi.lop$AMR[,tcn],pi.nop$AMR[,tcn],"pi in America population for 2-16nt seed")
fun_plot(pi.hip$EUR[,tcn],pi.lop$EUR[,tcn],pi.nop$EUR[,tcn],"pi in European population for 2-16nt seed")
fun_plot(pi.hip$SAS[,tcn],pi.lop$SAS[,tcn],pi.nop$SAS[,tcn],"pi in South Asia population for 2-16nt seed")
fun_plot(pi.hip$EAS[,tcn],pi.lop$EAS[,tcn],pi.nop$EAS[,tcn],"pi in East Asia population for 2-16nt seed")
dev.off()

# DAF and pi in high piA, low piA and no piA seed1(2-16nt) excluding RMSK----
daf.hip=list();daf.lop=list();daf.nop=list()
pi.hip=list();pi.lop=list();pi.nop=list()
pop=c("All","AFR","AMR","EUR","SAS","EAS")
for(i in pop){
  daf.hip[[i]]=read.table(paste("final_matrix_noRMSK/",i,".2-16.high.DAF.mat",sep=""),header=T,row.names=1)
  daf.lop[[i]]=read.table(paste("final_matrix_noRMSK/",i,".2-16.low.DAF.mat",sep=""),header=T,row.names=1)
  daf.nop[[i]]=read.table(paste("final_matrix_noRMSK/",i,".2-16.nopi.DAF.mat",sep=""),header=T,row.names=1)
  pi.hip[[i]]=read.table(paste("final_matrix_noRMSK/",i,".2-16.high.pi.mat",sep=""),header=T,row.names=1)
  pi.lop[[i]]=read.table(paste("final_matrix_noRMSK/",i,".2-16.low.pi.mat",sep=""),header=T,row.names=1)
  pi.nop[[i]]=read.table(paste("final_matrix_noRMSK/",i,".2-16.nopi.pi.mat",sep=""),header=T,row.names=1)
}
pdf("seed2-16.noRMSK.pdf",width=14,height=5,useDingbats = F)
tcn=c(21,22,3,4,25,26,19,20,7,6,5,8,15,14,13,16,11,12,27,1,2,17,18)
fun_plot(daf.hip$All[,tcn],daf.lop$All[,tcn],daf.nop$All[,tcn],"DAF in all population for 2-16nt seed")
fun_plot(daf.hip$AFR[,tcn],daf.lop$AFR[,tcn],daf.nop$AFR[,tcn],"DAF in Africa population for 2-16nt seed")
fun_plot(daf.hip$AMR[,tcn],daf.lop$AMR[,tcn],daf.nop$AMR[,tcn],"DAF in America population for 2-16nt seed")
fun_plot(daf.hip$EUR[,tcn],daf.lop$EUR[,tcn],daf.nop$EUR[,tcn],"DAF in European population for 2-16nt seed")
fun_plot(daf.hip$SAS[,tcn],daf.lop$SAS[,tcn],daf.nop$SAS[,tcn],"DAF in South Asia population for 2-16nt seed")
fun_plot(daf.hip$EAS[,tcn],daf.lop$EAS[,tcn],daf.nop$EAS[,tcn],"DAF in East Asia population for 2-16nt seed")
tcn=c(21,22,3,4,25,26,19,20,7,6,5,8,15,14,13,16,11,12,27,1,2,17,18)
fun_plot(pi.hip$All[,tcn],pi.lop$All[,tcn],pi.nop$All[,tcn],"pi in all population for 2-16nt seed")
fun_plot(pi.hip$AFR[,tcn],pi.lop$AFR[,tcn],pi.nop$AFR[,tcn],"pi in Africa population for 2-16nt seed")
fun_plot(pi.hip$AMR[,tcn],pi.lop$AMR[,tcn],pi.nop$AMR[,tcn],"pi in America population for 2-16nt seed")
fun_plot(pi.hip$EUR[,tcn],pi.lop$EUR[,tcn],pi.nop$EUR[,tcn],"pi in European population for 2-16nt seed")
fun_plot(pi.hip$SAS[,tcn],pi.lop$SAS[,tcn],pi.nop$SAS[,tcn],"pi in South Asia population for 2-16nt seed")
fun_plot(pi.hip$EAS[,tcn],pi.lop$EAS[,tcn],pi.nop$EAS[,tcn],"pi in East Asia population for 2-16nt seed")
dev.off()

# DAF and pi in high piA, low piA and no piA seed2(17-25nt) including RMSK----
daf.hip=list();daf.lop=list();daf.nop=list()
pi.hip=list();pi.lop=list();pi.nop=list()
pop=c("All","AFR","AMR","EUR","SAS","EAS")
for(i in pop){
  daf.hip[[i]]=read.table(paste("final_matrix/",i,".17-25.high.DAF.mat",sep=""),header=T,row.names=1)
  daf.lop[[i]]=read.table(paste("final_matrix/",i,".17-25.low.DAF.mat",sep=""),header=T,row.names=1)
  daf.nop[[i]]=read.table(paste("final_matrix/",i,".17-25.nopi.DAF.mat",sep=""),header=T,row.names=1)
  pi.hip[[i]]=read.table(paste("final_matrix/",i,".17-25.high.pi.mat",sep=""),header=T,row.names=1)
  pi.lop[[i]]=read.table(paste("final_matrix/",i,".17-25.low.pi.mat",sep=""),header=T,row.names=1)
  pi.nop[[i]]=read.table(paste("final_matrix/",i,".17-25.nopi.pi.mat",sep=""),header=T,row.names=1)
}
pdf("seed17-25.pdf",width=14,height=5,useDingbats = F)
tcn=c(21,22,3,4,25,26,19,20,7,6,5,8,15,14,13,16,11,12,27,1,2,17,18)
fun_plot(daf.hip$All[,tcn],daf.lop$All[,tcn],daf.nop$All[,tcn],"DAF in all population for 17-25nt seed")
fun_plot(daf.hip$AFR[,tcn],daf.lop$AFR[,tcn],daf.nop$AFR[,tcn],"DAF in Africa population for 17-25nt seed")
fun_plot(daf.hip$AMR[,tcn],daf.lop$AMR[,tcn],daf.nop$AMR[,tcn],"DAF in America population for 17-25nt seed")
fun_plot(daf.hip$EUR[,tcn],daf.lop$EUR[,tcn],daf.nop$EUR[,tcn],"DAF in European population for 17-25nt seed")
fun_plot(daf.hip$SAS[,tcn],daf.lop$SAS[,tcn],daf.nop$SAS[,tcn],"DAF in South Asia population for 17-25nt seed")
fun_plot(daf.hip$EAS[,tcn],daf.lop$EAS[,tcn],daf.nop$EAS[,tcn],"DAF in East Asia population for 17-25nt seed")
tcn=c(21,22,3,4,25,26,19,20,7,6,5,8,15,14,13,16,11,12,27,1,2,17,18)
fun_plot(pi.hip$All[,tcn],pi.lop$All[,tcn],pi.nop$All[,tcn],"pi in all population for 17-25nt seed")
fun_plot(pi.hip$AFR[,tcn],pi.lop$AFR[,tcn],pi.nop$AFR[,tcn],"pi in Africa population for 17-25nt seed")
fun_plot(pi.hip$AMR[,tcn],pi.lop$AMR[,tcn],pi.nop$AMR[,tcn],"pi in America population for 17-25nt seed")
fun_plot(pi.hip$EUR[,tcn],pi.lop$EUR[,tcn],pi.nop$EUR[,tcn],"pi in European population for 17-25nt seed")
fun_plot(pi.hip$SAS[,tcn],pi.lop$SAS[,tcn],pi.nop$SAS[,tcn],"pi in South Asia population for 17-25nt seed")
fun_plot(pi.hip$EAS[,tcn],pi.lop$EAS[,tcn],pi.nop$EAS[,tcn],"pi in East Asia population for 17-25nt seed")
dev.off()

# DAF and pi in high piA, low piA and no piA seed2(17-25nt) excluding RMSK----
daf.hip=list();daf.lop=list();daf.nop=list()
pi.hip=list();pi.lop=list();pi.nop=list()
pop=c("All","AFR","AMR","EUR","SAS","EAS")
for(i in pop){
  daf.hip[[i]]=read.table(paste("final_matrix_noRMSK/",i,".17-25.high.DAF.mat",sep=""),header=T,row.names=1)
  daf.lop[[i]]=read.table(paste("final_matrix_noRMSK/",i,".17-25.low.DAF.mat",sep=""),header=T,row.names=1)
  daf.nop[[i]]=read.table(paste("final_matrix_noRMSK/",i,".17-25.nopi.DAF.mat",sep=""),header=T,row.names=1)
  pi.hip[[i]]=read.table(paste("final_matrix_noRMSK/",i,".17-25.high.pi.mat",sep=""),header=T,row.names=1)
  pi.lop[[i]]=read.table(paste("final_matrix_noRMSK/",i,".17-25.low.pi.mat",sep=""),header=T,row.names=1)
  pi.nop[[i]]=read.table(paste("final_matrix_noRMSK/",i,".17-25.nopi.pi.mat",sep=""),header=T,row.names=1)
}
pdf("seed17-25.noRMSK.pdf",width=14,height=5,useDingbats = F)
tcn=c(21,22,3,4,25,26,19,20,7,6,5,8,15,14,13,16,11,12,27,1,2,17,18)
fun_plot(daf.hip$All[,tcn],daf.lop$All[,tcn],daf.nop$All[,tcn],"DAF in all population for 17-25nt seed")
fun_plot(daf.hip$AFR[,tcn],daf.lop$AFR[,tcn],daf.nop$AFR[,tcn],"DAF in Africa population for 17-25nt seed")
fun_plot(daf.hip$AMR[,tcn],daf.lop$AMR[,tcn],daf.nop$AMR[,tcn],"DAF in America population for 17-25nt seed")
fun_plot(daf.hip$EUR[,tcn],daf.lop$EUR[,tcn],daf.nop$EUR[,tcn],"DAF in European population for 17-25nt seed")
fun_plot(daf.hip$SAS[,tcn],daf.lop$SAS[,tcn],daf.nop$SAS[,tcn],"DAF in South Asia population for 17-25nt seed")
fun_plot(daf.hip$EAS[,tcn],daf.lop$EAS[,tcn],daf.nop$EAS[,tcn],"DAF in East Asia population for 17-25nt seed")
tcn=c(21,22,3,4,25,26,19,20,7,6,5,8,15,14,13,16,11,12,27,1,2,17,18)
fun_plot(pi.hip$All[,tcn],pi.lop$All[,tcn],pi.nop$All[,tcn],"pi in all population for 17-25nt seed")
fun_plot(pi.hip$AFR[,tcn],pi.lop$AFR[,tcn],pi.nop$AFR[,tcn],"pi in Africa population for 17-25nt seed")
fun_plot(pi.hip$AMR[,tcn],pi.lop$AMR[,tcn],pi.nop$AMR[,tcn],"pi in America population for 17-25nt seed")
fun_plot(pi.hip$EUR[,tcn],pi.lop$EUR[,tcn],pi.nop$EUR[,tcn],"pi in European population for 17-25nt seed")
fun_plot(pi.hip$SAS[,tcn],pi.lop$SAS[,tcn],pi.nop$SAS[,tcn],"pi in South Asia population for 17-25nt seed")
fun_plot(pi.hip$EAS[,tcn],pi.lop$EAS[,tcn],pi.nop$EAS[,tcn],"pi in East Asia population for 17-25nt seed")
dev.off()

# rebutal; cell; review1.2----
pdf("../../rebutal_cell/review1.2.pdf",width=5,height=4,useDingbats=F)
par(mar=c(4,4,3,1),tcl=0.3,bty="n",bty="n",cex=4/6)
t1=(rowMeans(df_srna[human_pachy,c(1,2,3,4,5,6,10,12,13,14,15,16)]))/(rowMeans(df_srna[human_pachy,c(7,8,9,18)]))
t1[which(t1>8)]=8
t2=(rowMeans(df_srna[human_prepachy,c(1,2,3,4,5,6,10,12,13,14,15,16)]))/(rowMeans(df_srna[human_prepachy,c(7,8,9,18)]))
t2[which(t2>8)]=8
t3=(rowMeans(df_srna[human_hybrid,c(1,2,3,4,5,6,10,12,13,14,15,16)]))/(rowMeans(df_srna[human_hybrid,c(7,8,9,18)]))
t3[which(t3>8)]=8
hist(t1,breaks=100,xlim=c(0,8),col=csl[1],main="gourpI+groupII / juv")
hist(t2,breaks=50,xlim=c(0,8),col=csl[2],add=T)
hist(t3,breaks=200,xlim=c(0,8),col=csl[3],add=T)
dev.off()

# rebutal; cell; review2.4----
pdf("../../rebutal_cell/review2.4.pdf",width=8,height=4,useDingbats=F)
par(mfrow=c(1,2),mar=c(4,4,3,3),tcl=0.3,cex=4/6,bty="n")
# piRNA
x=df_srna[human_pi,"N1643.5YO.Rep1.Ox"];y=df_srna[human_pi,"N1643.5YO.Rep2.Ox"]
plot(x,y,main="piRNA in Juv1 for 182 piRNA genes",xlab="RPM",ylab="RPM");abline(0,1)
c1=round(cor(x,y,method="spearman"),3)
c2=round(cor(x,y,method="pearson"),3)
text(0,2600,label=paste("soearman cor = ",c1,sep=""),pos=4)
text(0,2400,label=paste("pearson cor = ",c2,sep=""),pos=4)
# longRNA
x=df_rna[setdiff(row.names(df_rna),human_pi),"D1643_5YO_Rep1"]
y=df_rna[setdiff(row.names(df_rna),human_pi),"D1643_5YO_Rep2"]
x1=df_rna[human_pi,"D1643_5YO_Rep1"];y1=df_rna[human_pi,"D1643_5YO_Rep2"]
plot(x,y,main="longRNA in Juv1 for all genes and piRNA genes",xlab="RPKM",ylab="RPKM");abline(0,1)
points(x1,y1,col=csl[1],pch=19)
c1=round(cor(x,y,method="spearman"),3)
c2=round(cor(x,y,method="pearson"),3)
text(0,2500,label=paste("soearman cor = ",c1,sep=""),pos=4)
text(0,2300,label=paste("pearson cor = ",c2,sep=""),pos=4)
c1=round(cor(x1,y1,method="spearman"),3)
c2=round(cor(x1,y1,method="pearson"),3)
text(0,2100,label=paste("soearman cor = ",c1,sep=""),pos=4,col=csl[1])
text(0,1900,label=paste("pearson cor = ",c2,sep=""),pos=4,col=csl[1])
text(1500,500,label="all genes",col="black")
text(1500,300,label="piRNA genes",col=csl[1])
dev.off()

# rebutal; cell; review2.8----
random.pi=read.table("N8150-Ox.25_31nt.randomRegion.rpm",header=F,row.names=1)
random.rna=read.table("8150NT.randomRegion.rpkm",header=F,row.names=1)
linc.pi=read.table("N8150-Ox.lincRNA.rpm",header=F,row.names=1)
mRNA.pi=read.table("N8150-Ox.mRNA.rpm",header=F,row.names=1)
pdf("../rebutal_cell/review2.8.pdf",width=5,height=4,useDingbats=F)
par(mfrow=c(1,2),bty="n",mar=c(7,4,1,1),cex=4/6,tcl=0.3)
boxplot(log10(random.pi[,1]+1),log10(linc.pi[,1]+1),log10(mRNA.pi[,1]+1),
        log10(df_srna[human_prepachy,12]+1),
        log10(df_srna[human_hybrid,12]+1),log10(df_srna[human_pachy,12]+1),
        col=csl[c(9,7,6,3,2,1)],yaxt="n",xaxt="n",ylab="piRNA; RPM",
        outline=F)
add_axis(2)
axis(1,1:6,label=c("randomRegion","lincRNA","mRNA","Prepachy","Hybrid","Pachy"),lwd=0,las=3)
boxplot(log10(random.rna[,1]+1),log10(df_rna[human_linc,11]+1),log10(df_rna[human_pc,11]+1),
        log10(df_rna[human_prepachy,11]+1),
        log10(df_rna[human_hybrid,11]+1),log10(df_rna[human_pachy,11]+1),
        col=csl[c(9,7,6,3,2,1)],ylab="longRNA; RPKM",
        outline=F,yaxt="n")
add_axis(2)
axis(1,1:6,label=c("randomRegion","lincRNA","mRNA","Prepachy","Hybrid","Pachy"),lwd=0,las=3)
dev.off()

# rebutal; cell; review3.1----
library(affy)
tn=c(human_prepachy,human_hybrid,human_pachy)
tcs=c(rep(csl[3],83),rep(csl[2],10),rep(csl[1],89))
pdf("../rebutal_cell/review3.1.pdf",width=4,height=3,useDingbats=F)
par(mar=c(4,4,3,1),tcl=0.3,bty="n",cex=5/6)
M=(rowMeans(df_srna[tn,c(1,2,3,4,5,6,10,12,13,14,16)]))/(rowMeans(df_srna[tn,c(7,8,9,18)]))
A=rowMeans(df_srna[tn,c(1,2,3,4,5,6,10,12,13,14,16,7,8,9,18)])
plot(log10(A+1),log2(M),xlim=c(1,5),ylim=c(-10,10),col=tcs,
     xlab="meam; log10RPM",ylab="fold enrichment; log2(adult/juv)",
     main="groupI+groupII / juv")
abline(0,0)
dev.off()

# rebutal; cell; review A6 issue----
pdf("../rebutal_cell/A6_issue.pdf",width=4,height=3,useDingbats=F)
par(mar=c(4,4,3,1),tcl=0.3,bty="n",cex=5/6)
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
df_rna=read.table("rebuttal.rnaseq.gene+pi.rpkm",header=T,row.names=1,check.names=F)
trna_rn=rna_rn;trna_rn[9]="S6-A7-ID9_S6"
df_rna=df_rna[,trna_rn]
annoR=data.frame(Length_Type=factor(rep(c("Juv","l26nt","inter","l30nt"),c(3,3,5,6))))
rna_rn1=c("D1643_5YO",rna_rn[3:18]);rna_rn1[8]="S6-A7-ID9_S6"
row.names(annoR)=rna_rn1
df_rnat=cbind(apply(df_rna[,1:2],1,mean),df_rna[,-c(1,2)])
colnames(df_rnat)[1]="D1643_5YO"
ann_colors=list(Length_Type=c(Juv="#8dd3c7",l26nt="#fccde5",inter="#b3de69",l30nt="#fdb462"))
pheatmap(cor(df_rnat[setdiff(row.names(df_rnat),human_pi),rna_rn1],
             method="spearman"),annotation_row=annoR,annotation_col=annoR,
         annotation_colors=ann_colors,
         filename="../rebutal_cell/A6_issue.pdf",
         cellwidth=24,cellheight=24,main="gene expression correlation between samples")
dev.off()

# rebutal; cell; correlation between piRNA abundance; long RNA abundance and DAF----
tdaf=read.table("human.pachy.exon.individulas.DAF",header=F,row.names=1)
tdaf[which(tdaf==(-1)),1]=NA
tn=human_pachy[which(!is.na(tdaf[human_pachy,1]))]
c1=round(cor(log10(df_srna[tn,11]+1),tdaf[tn,1],method="spearman"),3)
c2=round(cor(log10(df_rna[tn,15]+1),tdaf[tn,1],method="spearman"),3)
pdf("../../../Ozata et al. 1st Paper/Rebutal Cell/cor.DAFvspiRNAandLongRNAForPachy.pdf",width=8,height=4,useDingbats=F)
par(mfrow=c(1,2),mar=c(4,4,1,1),cex=5/6,bty="n",tcl=0.3)
plot(log10(df_srna[human_pachy,11]+1),tdaf[human_pachy,1],
     pch=21,col="white",bg="black",lwd=0.5,xlab="log10 piRNA RPM",ylab="average DAF",
     xlim=c(0,4),ylim=c(0,0.1))
text(0,0.09,label=paste("spearman R = ",c1,sep=""),pos=4)
plot(log10(df_rna[human_pachy,15]+1),tdaf[human_pachy,1],
     pch=21,col="white",bg="black",lwd=0.5,xlab="log10 longRNA RPKM",ylab="average DAF",
     xlim=c(0,2.2),ylim=c(0,0.1))
text(0,0.09,label=paste("spearman R = ",c2,sep=""),pos=4)
dev.off()

# rebutal; NEE; review1.2----
treads=read.table("../rebutal_NEE/tables/summary.reads.txt",header=T,row.names=1)
tspecies=read.table("../rebutal_NEE/tables/summary.species.txt",header=T,row.names=1)
treads=treads[c(1:30,55:102,109:114),]
tspecies=tspecies[c(1:30,55:102,109:114),]
r1=treads[(1:14)*6-2,2]/treads[(1:14)*6-2,1]
r2=treads[(1:14)*6-3,2]/treads[(1:14)*6-3,1]
r3=treads[(1:14)*6,2]/treads[(1:14)*6,1]
r4=treads[(1:14)*6-1,2]/treads[(1:14)*6-1,1]
r5=treads[(1:14)*6-4,2]/treads[(1:14)*6-4,1]
r6=treads[(1:14)*6-5,2]/treads[(1:14)*6-5,1]
r7=c(1707/5544,1508/5476)
r8=c(1972/9319,1613/7018)
pdf("../rebutal_NEE/figures/reviewer1.2.pdf",width=4.5,height=3.5,useDingbats=F)
par(mar=c(4,4,1,1),tcl=0.3,bty="n")
tdf=100*c(mean(r1),mean(r2),mean(r3),mean(r4),mean(r5),mean(r6),mean(r7),mean(r8))
tdf1=100*c(sd(r1),sd(r2),sd(r3),sd(r4),sd(r5),sd(r6),sd(r7),sd(r8))
barplot(tdf,col=cs_gg[1:2],space=c(0,0,1,0,1,0,1,0),border="white",
        ylab="percentage of TE piRNAs reads(%)",ylim=c(0,35))
axis(1,c(1,4,7,10),label=c("pachy","prepachy","hybrid","fetal"),lwd=0)
ps=c(1,2,4,5,7,8,10,11)-0.5
text(ps,tdf/2,label=paste(round(tdf,1),"%",sep=""),srt=90)
legend("topright",legend=c("uniqMappers","multiMappers"),col=cs_gg[1:2],
       pch=15,bty="n")
for(i in 1:8){
  lines(c(ps[i],ps[i]),c(tdf[i]-tdf1[i],tdf[i]+tdf1[i]))
}

r1=tspecies[(1:14)*6-2,2]/tspecies[(1:14)*6-2,1]
r2=tspecies[(1:14)*6-3,2]/tspecies[(1:14)*6-3,1]
r3=tspecies[(1:14)*6,2]/tspecies[(1:14)*6,1]
r4=tspecies[(1:14)*6-1,2]/tspecies[(1:14)*6-1,1]
r5=tspecies[(1:14)*6-4,2]/tspecies[(1:14)*6-4,1]
r6=tspecies[(1:14)*6-5,2]/tspecies[(1:14)*6-5,1]
r7=c(1656/5014,1468/5068)
r8=c(1852/7167,1523/4927)
tdf=100*c(mean(r1),mean(r2),mean(r3),mean(r4),mean(r5),mean(r6),mean(r7),mean(r8))
tdf1=100*c(sd(r1),sd(r2),sd(r3),sd(r4),sd(r5),sd(r6),sd(r7),sd(r8))
barplot(tdf,col=cs_gg[1:2],space=c(0,0,1,0,1,0,1,0),border="white",
        ylab="percentage of TE piRNAs species(%)",ylim=c(0,50))
axis(1,c(1,4,7,10),label=c("pachy","prepachy","hybrid","fetal"),lwd=0)
ps=c(1,2,4,5,7,8,10,11)-0.5
text(ps,tdf/2,label=paste(round(tdf,1),"%",sep=""),srt=90)
legend("topright",legend=c("uniqMappers","multiMappers"),col=cs_gg[1:2],
       pch=15,bty="n")
for(i in 1:8){
  lines(c(ps[i],ps[i]),c(tdf[i]-tdf1[i],tdf[i]+tdf1[i]))
}
dev.off()

# rebutal; NEE; review1.3----
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/TCFL5_project/")
mouse_pachy=as.vector(read.table("../../My_Projects/Epigenetic_Regulation_of_Mouse_piRNA_Genes_during_Spermatogenesis/Data/piRNA_list/pachytene.txt",header=FALSE,row.names=NULL)[,1])
amybScore=read.table("matrix/motifScore/mm10.gene+pi.myba.motifScore",header=F,row.names=1)
setwd("/Users/BigBear/Dropbox (UMass Medical School)/IGV_piRNA_Cluster_12April17/Metagene_Plots/matrix/")
amybSignal=read.table("../rebutal_NEE/tables/mm10.gene+pi.AMYB_rpm.txt",header=F,row.names=1)
tdf=cbind(log2(amybSignal[mouse_pachy,1]+1),amybScore[mouse_pachy,2])
row.names(tdf)=mouse_pachy
pheatmap(t(tdf),cellwidth=4,cellheight=16,cluster_rows=F,
         filename="../rebutal_NEE/figures/reviewer1.3.pdf",
         show_rownames=T,show_colnames=T,
         labels_row=c("AMYB signal (log2(RPM+1))","AMYB motif score"),
         fontsize_row=16,fontsize_col=4)
# rebutal; NEE; review1.5----
l140=as.vector(read.table("../rebutal_NEE/tables/list_140genes.txt",header=F,row.names=NULL)[,1])
lrepro=as.vector(read.table("../rebutal_NEE/tables/reproduction.gene.list",header=F,row.names=NULL)[,1])
limmun=as.vector(read.table("../rebutal_NEE/tables/immune_system_process.gene.list",header=F,row.names=NULL)[,1])

pdf("../rebutal_NEE/figures/reviewer1.5.pdf",width=4,height=4,useDingbats=F)
par(mar=c(4,4,3,1),tcl=0.3,bty="n")
gjuv=c(7,8,12,17);g1=c(15,9,5,3,11,4);g2=c(1,6,13,2);g3=c(10,14,16,18)
boxplot(apply(df_rna[l140,gjuv],1,mean),
        apply(df_rna[l140,g1],1,mean),
        apply(df_rna[l140,g2],1,mean),
        apply(df_rna[l140,g3],1,mean),outline=F,
        staplewex=0,lty=1,ylim=c(0,60),ylab="RPKM",main="the 140 genes (DEG)")
boxplot(apply(df_rna[lrepro,gjuv],1,mean),
        apply(df_rna[lrepro,g1],1,mean),
        apply(df_rna[lrepro,g2],1,mean),
        apply(df_rna[lrepro,g3],1,mean),outline=F,
        staplewex=0,lty=1,ylim=c(0,60),ylab="RPKM",main="reproduction genes (n=1,416)")
boxplot(apply(df_rna[limmun,gjuv],1,mean),
        apply(df_rna[limmun,g1],1,mean),
        apply(df_rna[limmun,g2],1,mean),
        apply(df_rna[limmun,g3],1,mean),outline=F,
        staplewex=0,lty=1,ylim=c(0,60),ylab="RPKM",main="immune system process gene (n=3,237)")
dev.off()

# rebutal NEE; review2.2----
ta=c(99,1)
tb=c(49,40)
tc=c(104,85)
pdf("../rebutal_NEE/figures/reviewer2.2.pdf",width=4,height=4,useDingbats=F)
pie(ta,col=cs_gg[1:2],border="white",labels=c("AMYB 99%","non-AMYB 1%"),main="mouse")
pie(tb,col=cs_gg[1:2],border="white",labels=c("AMYB 55%","non-AMYB 45%"),main="human")
pie(tc,col=cs_gg[1:2],border="white",labels=c("AMYB 55%","non-AMYB 45%"),main="rhesus")
dev.off()
pdf("../rebutal_NEE/figures/reviewer2.2.barplot.pdf",width=2,height=4,useDingbats=F)
par(mar=c(10,4,1,1))
barplot(ta,col=cs_gg[1:2],border="white",names=c("AMYB 99%","non-AMYB 1%"),las=3,ylab="mouse")
barplot(tb,col=cs_gg[1:2],border="white",names=c("AMYB 55%","non-AMYB 45%"),las=3,ylab="human")
barplot(tc,col=cs_gg[1:2],border="white",names=c("AMYB 55%","non-AMYB 45%"),las=3,ylab="rhesus")
dev.off()

# rebutal NEE; review2.3----
# mouse MAF and pi analysis
# 1. set cutoff to devied high and low expressed regions
tdf=read.table("../rebutal_NEE/tables/mouse.25_31nt.uniq.accum",header=F,row.names=NULL)
pdf("../rebutal_NEE/figures/reviewer2.3_cutoff.pdf",width=4,height=4,useDingbats=F)
par(mar=c(4,4,1,1),tcl=0.3,cex=5/6,bty="n")
plot(tdf[,2],ylim=c(0,1000000),type="l",ylab="percentage of total piRNA (%)",
     xlab="number of piRNA-production loci",yaxt="n")
axis(2,c(0,200000,400000,600000,800000,1000000),label=c(0,20,40,60,80,100))
abline(h=666666,v=149985,lty=2)
text(500000,666666,pos=1,label="cutoff = 0.891 RPM")
dev.off()
# 2. MAF and pi
pdf("../rebutal_NEE/figures/reviewer2.3.MAFandpi.pdf",width=6,height=6,useDingbats=F)
all.daf=read.table("../rebutal_NEE/tables/All.DAF.mat",header=T,row.names=1)
all.pi=read.table("../rebutal_NEE/tables/All.pi.mat",header=T,row.names=1)
par(mar=c(10,4,1,1),tcl=0.3,bty="n",cex=5/6)
plot(unlist(all.daf[1,]),xaxt="n",ylim=c(0.05,0.1),xlab="",ylab="MAF")
for(i in 1:15){lines(c(i,i),c(all.daf[2,i],all.daf[3,i]))}
axis(1,1:15,labels=colnames(all.daf),las=3)
plot(unlist(all.pi[1,]),xaxt="n",ylim=c(0.08,0.15),xlab="",ylab="pi")
for(i in 1:15){lines(c(i,i),c(all.pi[2,i],all.pi[3,i]))}
axis(1,1:15,labels=colnames(all.daf),las=3)
dev.off()

# rebutal NEE; review2.7----
treads=read.table("../rebutal_NEE/tables/summary.reads.txt",header=T,row.names=1)
tspecies=read.table("../rebutal_NEE/tables/summary.species.txt",header=T,row.names=1)
pdf("../rebutal_NEE/figures/reviewer2.5.pdf",width=12,height=5,useDingbats=F)
par(mar=c(8,4,1,1),tcl=0.3,bty="n")
barplot(as.vector(treads[,2]/treads[,1]),space=rep(c(1,0,0,0,0,0),19),
        col=csl[c(3,3,1,1,2,2)],density=c(-1,30),
        ylab="percentage of TE piRNAs reads(%)")
text(1,c(0.6,0.55,0.5),pos=4,
     label=c("hybrid","pachytene","prepachytene"),col=csl[c(3,1,2)])
lbs=unlist(strsplit(row.names(treads),".25_31nt."))[(1:19)*12-1]
axis(1,(1:19)*7-2.5,label=lbs,lwd=0,las=3)
legend("topright",legend=c("unique mappers","multiple mappers"),
       col="black",density=c(-1,30),bty="n")

barplot(as.vector(tspecies[,2]/tspecies[,1]),space=rep(c(1,0,0,0,0,0),19),
        col=csl[c(3,3,1,1,2,2)],density=c(-1,30),
        ylab="percentage of TE piRNAs species(%)")
text(1,c(0.6,0.55,0.5),pos=4,
     label=c("hybrid","pachytene","prepachytene"),col=csl[c(3,1,2)])
axis(1,(1:19)*7-2.5,label=lbs,lwd=0,las=3)
legend("topright",legend=c("unique mappers","multiple mappers"),
       col="black",density=c(-1,30),bty="n")
dev.off()

# rebutal NEE; review2.4----
# >=2 reads
anno_bear=read.table("../piRNA_cluster_annotation/comparison/piG.exon.bed.MT2R.sum",header=F,row.names=1)
anno_xxx=read.table("../piRNA_cluster_annotation/comparison/previous_annotation.XXX.bed.MT2R.sum",header=F,row.names=1)
anno_jinchuan=read.table("../piRNA_cluster_annotation/comparison/previous_annotation.Jinchuan.bed.MT2R.sum",header=F,row.names=1)
anno_william=read.table("../piRNA_cluster_annotation/comparison/previous_annotation.willian.bed.MT2R.sum",header=F,row.names=1)
anno_bear1=read.table("../piRNA_cluster_annotation/comparison/piG.exon.bed.MT2R.uniq.sum",header=F,row.names=1)
anno_xxx1=read.table("../piRNA_cluster_annotation/comparison/previous_annotation.XXX.bed.MT2R.uniq.sum",header=F,row.names=1)
anno_jinchuan1=read.table("../piRNA_cluster_annotation/comparison/previous_annotation.Jinchuan.bed.MT2R.uniq.sum",header=F,row.names=1)
anno_william1=read.table("../piRNA_cluster_annotation/comparison/previous_annotation.willian.bed.MT2R.uniq.sum",header=F,row.names=1)
fun_t=function(d,ct){
  x1=c();x2=c();tl=0;ta=0
  for(i in 1:dim(d)[1]){
    if(tl<ct){tl=tl+d[i,3];ta=ta+d[i,1]/10000;x1=c(x1,tl);x2=c(x2,ta)}
  }
  om=cbind(x1,x2);return(om)
}
pdf("../rebutal_NEE/figures/reviewer2.4.>=2reads.pdf",width=8,height=4)
par(mar=c(4,4,3,1),tcl=0.3,bty="n",cex=5/6,mfrow=c(1,2))
ct=3000000
ss=fun_t(anno_bear,ct=ct)
plot(ss[,1],ss[,2],type="l",xlab="Total Mbp in piclusters",ylab="piRNA percentage(%)",
     xlim=c(0,3000000),xaxt="n",lwd=3,col="black",main="AllMappers",ylim=c(0,100))
axis(1,c(0,10,20,30)*100000,label=c(0,1,2,3))
text(c(3,3,3,3)*10^6,c(2,7,12,17),
     label=c("This study","Previous: Girard","Previous: Jinchuan","Previous: Williams"),
     col=c("black",cs_gg[2:4]),font=c(2,1,1,1),pos=2)
ss=fun_t(anno_xxx,ct=ct)
lines(ss[,1],ss[,2],lwd=2,col=cs_gg[2])
ss=fun_t(anno_jinchuan,ct=ct)
lines(ss[,1],ss[,2],lwd=2,col=cs_gg[3])
ss=fun_t(anno_william,ct=ct)
lines(ss[,1],ss[,2],lwd=2,col=cs_gg[4])
ss=fun_t(anno_bear1,ct=ct)
plot(ss[,1],ss[,2],type="l",xlab="Total Mbp in piclusters",ylab="piRNA percentage(%)",
     xlim=c(0,3000000),xaxt="n",lwd=3,col="black",main="UniqueMappers",ylim=c(0,100))
axis(1,c(0,10,20,30)*100000,label=c(0,1,2,3))
ss=fun_t(anno_xxx1,ct=ct)
lines(ss[,1],ss[,2],lwd=2,col=cs_gg[2])
ss=fun_t(anno_jinchuan1,ct=ct)
lines(ss[,1],ss[,2],lwd=2,col=cs_gg[3])
ss=fun_t(anno_william1,ct=ct)
lines(ss[,1],ss[,2],lwd=2,col=cs_gg[4])
text(c(3,3,3,3)*10^6,c(2,7,12,17),
     label=c("This study","Previous: Girard","Previous: Jinchuan","Previous: Williams"),
     col=c("black",cs_gg[2:4]),font=c(2,1,1,1),pos=2)
dev.off()
# all
anno_bear=read.table("../piRNA_cluster_annotation/comparison/piG.exon.bed.ALL.sum",header=F,row.names=1)
anno_xxx=read.table("../piRNA_cluster_annotation/comparison/previous_annotation.XXX.bed.ALL.sum",header=F,row.names=1)
anno_jinchuan=read.table("../piRNA_cluster_annotation/comparison/previous_annotation.Jinchuan.bed.ALL.sum",header=F,row.names=1)
anno_william=read.table("../piRNA_cluster_annotation/comparison/previous_annotation.willian.bed.ALL.sum",header=F,row.names=1)
anno_bear1=read.table("../piRNA_cluster_annotation/comparison/piG.exon.bed.ALL.uniq.sum",header=F,row.names=1)
anno_xxx1=read.table("../piRNA_cluster_annotation/comparison/previous_annotation.XXX.bed.ALL.uniq.sum",header=F,row.names=1)
anno_jinchuan1=read.table("../piRNA_cluster_annotation/comparison/previous_annotation.Jinchuan.bed.ALL.uniq.sum",header=F,row.names=1)
anno_william1=read.table("../piRNA_cluster_annotation/comparison/previous_annotation.willian.bed.ALL.uniq.sum",header=F,row.names=1)
fun_t=function(d,ct){
  x1=c();x2=c();tl=0;ta=0
  for(i in 1:dim(d)[1]){
    if(tl<ct){tl=tl+d[i,3];ta=ta+d[i,1]/10000;x1=c(x1,tl);x2=c(x2,ta)}
  }
  om=cbind(x1,x2);return(om)
}
pdf("../rebutal_NEE/figures/reviewer2.4.ALL.pdf",width=8,height=4)
par(mar=c(4,4,3,1),tcl=0.3,bty="n",cex=5/6,mfrow=c(1,2))
ct=3000000
ss=fun_t(anno_bear,ct=ct)
plot(ss[,1],ss[,2],type="l",xlab="Total Mbp in piclusters",ylab="piRNA percentage(%)",
     xlim=c(0,3000000),xaxt="n",lwd=3,col="black",main="AllMappers",ylim=c(0,100))
axis(1,c(0,10,20,30)*100000,label=c(0,1,2,3))
text(c(3,3,3,3)*10^6,c(2,7,12,17),
     label=c("This study","Previous: Girard","Previous: Jinchuan","Previous: Williams"),
     col=c("black",cs_gg[2:4]),font=c(2,1,1,1),pos=2)
ss=fun_t(anno_xxx,ct=ct)
lines(ss[,1],ss[,2],lwd=2,col=cs_gg[2])
ss=fun_t(anno_jinchuan,ct=ct)
lines(ss[,1],ss[,2],lwd=2,col=cs_gg[3])
ss=fun_t(anno_william,ct=ct)
lines(ss[,1],ss[,2],lwd=2,col=cs_gg[4])
ss=fun_t(anno_bear1,ct=ct)
plot(ss[,1],ss[,2],type="l",xlab="Total Mbp in piclusters",ylab="piRNA percentage(%)",
     xlim=c(0,3000000),xaxt="n",lwd=3,col="black",main="UniqueMappers",ylim=c(0,100))
axis(1,c(0,10,20,30)*100000,label=c(0,1,2,3))
ss=fun_t(anno_xxx1,ct=ct)
lines(ss[,1],ss[,2],lwd=2,col=cs_gg[2])
ss=fun_t(anno_jinchuan1,ct=ct)
lines(ss[,1],ss[,2],lwd=2,col=cs_gg[3])
ss=fun_t(anno_william1,ct=ct)
lines(ss[,1],ss[,2],lwd=2,col=cs_gg[4])
text(c(3,3,3,3)*10^6,c(2,7,12,17),
     label=c("This study","Previous: Girard","Previous: Jinchuan","Previous: Williams"),
     col=c("black",cs_gg[2:4]),font=c(2,1,1,1),pos=2)
dev.off()

# rebutal NEE; review2.5----
# rebutal NEE; review2.6----
# rebutal NEE; review3.6----
sps=c("human","rhesus","marmoset");genes=c("CBL","PDPK1","USP49","SYNM")
yms=list(CBL=20,PDPK1=10,USP49=2,SYNM=2)
regions=list();exons=list();crick=list();watson=list()
for(sp in sps){
  regions[[sp]]=read.table(paste("../rebutal_NEE/tables/review3.6/",sp,".region.bed",sep=""),header=F,row.names=4)
  exons[[sp]]=read.table(paste("../rebutal_NEE/tables/review3.6/",sp,".exon.bed",sep=""),header=F,row.names=NULL)
  crick[[sp]]=read.table(paste("../rebutal_NEE/tables/review3.6/",sp,".crick.bdg",sep=""),header=F,row.names=NULL)
  watson[[sp]]=read.table(paste("../rebutal_NEE/tables/review3.6/",sp,".watson.bdg",sep=""),header=F,row.names=NULL)
}
pdf("../rebutal_NEE/figures/reviewer3.6.pdf",width=6,height=6,useDingbats=F)
par(mfrow=c(3,1),bty="n",cex=5/6,tcl=0.3,mar=c(2,4,1,1))
for(g in genes){
  for(sp in sps){
    tr=regions[[sp]];te=exons[[sp]];tc=crick[[sp]];tw=watson[[sp]]
    chrom=as.vector(tr[g,1])
    start=as.vector(tr[g,2])
    end=as.vector(tr[g,3])
    te=te[which(te[,1]==chrom & te[,2]>=start & te[,3]<=end),]
    tc=tc[which(tc[,1]==chrom & tc[,2]>=start & tc[,3]<=end),]
    tw=tw[which(tw[,1]==chrom & tw[,2]>=start & tw[,3]<=end),]
    ym=yms[[g]]
    plot(NA,xlim=c(start,end),ylim=c(-ym,ym),xaxt="n",ylab="RPM",xlab="")
    text(start,ym*18/20,label=paste(sp,"; ",chrom,":",start,"-",end,sep=""),pos=4)
    if(te[1,6]=="-"){tcs=cs_gg[2]}else{tcs=cs_gg[1]}
    lines(c(te[1,2],te[nrow(te),3]),rep(-ym,2),col=tcs)
    axis(1,(start+end)/2,label=g,lwd=0)
    for(i in 1:nrow(te)){
      polygon(c(te[i,2],te[i,3],te[i,3],te[i,2]),
              c(-ym*20/19,-ym*20/19,-ym*18/19,-ym*18/19),
              col=tcs,border=tcs)
    }
    for(i in 1:nrow(tw)){
      polygon(c(tw[i,2],tw[i,3],tw[i,3],tw[i,2]),c(0,0,tw[i,4],tw[i,4]),
              col=cs_gg[1],border=cs_gg[1])
    }
    for(i in 1:nrow(tc)){
      polygon(c(tc[i,2],tc[i,3],tc[i,3],tc[i,2]),c(0,0,tc[i,4],tc[i,4]),
              col=cs_gg[2],border=cs_gg[2])
    }
  }
}
dev.off()

# rebutal NEE; review3.7----
tdf=read.table("../rebutal_NEE/tables/hg19.gene+pi.len",header=T,row.names=1)
pdf("../rebutal_NEE/figures/reviewer3.7.pdf",width=3,height=4,useDingbats=F)
par(mar=c(7,4,2,1),tcl=0.3,bty="n",las=3)
boxplot(tdf[human_pachy,1],tdf[human_pachy,2],
        tdf[human_prepachy,1],tdf[human_prepachy,2],
        tdf[human_hybrid,1],tdf[human_hybrid,2],
        tdf[human_pc,1],tdf[human_pc,2],
        tdf[human_linc,1],tdf[human_linc,2],
        at=c(1,2,4,5,7,8,10,11,13,14),outline=F,lty=1,
        staplewex=0,col=cs_gg[1:2],xaxt="n",ylab="length (bp)")
axis(1,c(1,4,7,10,13)+0.5,label=c("pachy","prepachy","hybrid","mRNA",
                                  "lincRNA"),lwd=0)
text(c(0.5,0.5),c(200000,220000),label=c("intron","exon"),col=cs_gg[2:1],pos=4)
tc=c(median((tdf[human_pachy,2]+1)/(tdf[human_pachy,1]+1)),
     median((tdf[human_prepachy,2]+1)/(tdf[human_prepachy,1]+1)),
     median((tdf[human_hybrid,2]+1)/(tdf[human_hybrid,1]+1)),
     median((tdf[human_pc,2]+1)/(tdf[human_pc,1]+1)),
     median((tdf[human_linc,2]+1)/(tdf[human_linc,1]+1)))
par(mar=c(6,4,1,1),tcl=0.3,bty="n",las=3)
barplot(tc,space=0,border="white",
        col=csl[1:5],ylab="median intron_length/exon_length",
        names=c("pachy","prepachy","hybrid","mRNA","lincRNA"))
text(1:5-0.5,c(1,tc[2:5]/2),label=round(tc,2),srt=90)
dev.off()

# rebutal NEE; pre-pachytene piRNA definition question----
t2=(rowMeans(df_srna[human_prepachy,c(6,11,15,17,10,12,14,16)]))/(rowMeans(df_srna[human_prepachy,c(8,9)]))
n1=names(t2)[which(t2<=1)] # FC≤1
n2=names(t2)[which(t2>1)] # FC between 1 and 2
write.table(n1,"../matrix/Prepachy.FC0-1.list",row.names=F,col.names=F,sep="\t",
            quote=F)
write.table(n2,"../matrix/Prepachy.FC1-2.list",row.names=F,col.names=F,sep="\t",
            quote=F)
  # lendis for pachy, prepachy, hybrid, FC0-1 prepachy and FC1-2 prepachy----
aver.pachy.lendis=read.table("lendis/merged.pachytene.lendis",header=T,row.names=1,check.names=F)
aver.prepachy.lendis=read.table("lendis/merged.prepachy.lendis",header=T,row.names=1,check.names=F)
aver.prepachy.1.lendis=read.table("lendis/merged.prepachy.FC0-1.lendis",header=T,row.names=1,check.names=F)
aver.prepachy.2.lendis=read.table("lendis/merged.prepachy.FC1-2.lendis",header=T,row.names=1,check.names=F)
aver.hybrid.lendis=read.table("lendis/merged.hybrid.lendis",header=T,row.names=1,check.names=F)
pdf("../Bear_Jul12.2018/Lendis.AllGroup.Allpi.pdf",width=12,height=15,useDingbats = F)
par(mfrow=c(5,4),mar=c(2,4,3,1),tcl=0.3,cex=5/6,bty="n",xpd=F)
fun_t=function(df,n,m,mm,yl){
  barplot(apply(df[4:21,n],1,mean),space=0,border="white",col="black",names=18:35,
          ylab=yl,main=m,ylim=c(0,mm))
}
fun_t(aver.pachy.lendis,srna_rn[1:4],"juv",150000,"pachytene")
fun_t(aver.pachy.lendis,srna_rn[13:18],"groupI (30nt)",150000,"")
fun_t(aver.pachy.lendis,srna_rn[8:12],"groupI (intermediate)",150000,"")
fun_t(aver.pachy.lendis,srna_rn[5:7],"groupI (26nt)",150000,"")
fun_t(aver.prepachy.lendis,srna_rn[1:4],"",15000,"pre-pachytene")
fun_t(aver.prepachy.lendis,srna_rn[13:18],"",15000,"")
fun_t(aver.prepachy.lendis,srna_rn[8:12],"",15000,"")
fun_t(aver.prepachy.lendis,srna_rn[5:7],"",15000,"")
fun_t(aver.prepachy.1.lendis,srna_rn[1:4],"",15000,"pre-pachytene (FC0-1)")
fun_t(aver.prepachy.1.lendis,srna_rn[13:18],"",15000,"")
fun_t(aver.prepachy.1.lendis,srna_rn[8:12],"",15000,"")
fun_t(aver.prepachy.1.lendis,srna_rn[5:7],"",15000,"")
fun_t(aver.prepachy.2.lendis,srna_rn[1:4],"",1500,"pre-pachytene (FC1-2)")
fun_t(aver.prepachy.2.lendis,srna_rn[13:18],"",1500,"")
fun_t(aver.prepachy.2.lendis,srna_rn[8:12],"",1500,"")
fun_t(aver.prepachy.2.lendis,srna_rn[5:7],"",1500,"")
fun_t(aver.hybrid.lendis,srna_rn[1:4],"",500,"hybrid")
fun_t(aver.hybrid.lendis,srna_rn[13:18],"",500,"")
fun_t(aver.hybrid.lendis,srna_rn[8:12],"",500,"")
fun_t(aver.hybrid.lendis,srna_rn[5:7],"",500,"")
dev.off()
  # individual lendis ordered by FC for piRNA genes----
fun_t=function(pn,fc){
  tdf=read.table(paste("lendis/detailed/",pn,".merged.lendis",sep=""),check.names=F,header=T,row.names=1)
  mm=max(apply(tdf[4:21,srna_rn[1:4]],1,mean),apply(tdf[4:21,srna_rn[13:18]],1,mean),
         apply(tdf[4:21,srna_rn[8:12]],1,mean),apply(tdf[4:21,srna_rn[5:7]],1,mean))
  barplot(apply(tdf[4:21,srna_rn[1:4]],1,mean),space=0,border="white",col=cs_gg[4],names=18:35,
          ylab=paste(pn," FC=",fc,sep=""),main="",ylim=c(0,mm))
  barplot(apply(tdf[4:21,srna_rn[13:18]],1,mean),space=0,border="white",col=cs_gg[1],names=18:35,
          ylab="",main="",ylim=c(0,mm))
  barplot(apply(tdf[4:21,srna_rn[8:12]],1,mean),space=0,border="white",col=cs_gg[2],names=18:35,
          ylab="",main="",ylim=c(0,mm))
  barplot(apply(tdf[4:21,srna_rn[5:7]],1,mean),space=0,border="white",col=cs_gg[3],names=18:35,
          ylab="",main="",ylim=c(0,mm))
}
pdf("../Updated_Figures_NEE/lendis.all_prepachy.pdf",width=12,height=3,useDingbats = F)
par(mfrow=c(1,4),mar=c(2,4,3,1),tcl=0.3,cex=5/6,bty="n",xpd=F)
t2=(rowMeans(df_srna[human_prepachy,c(6,11,15,17,10,12,14,16)]))/(rowMeans(df_srna[human_prepachy,c(8,9)]))
o1=order(t2,decreasing=T)
for(i in human_prepachy[o1]){
  fc=round(t2[i],2)
  fun_t(i,fc)
}
dev.off()
pdf("../Updated_Figures_NEE/lendis.all_hybrid.pdf",width=12,height=3,useDingbats = F)
par(mfrow=c(1,4),mar=c(2,4,3,1),tcl=0.3,cex=5/6,bty="n",xpd=F)
t2=(rowMeans(df_srna[human_hybrid,c(6,11,15,17,10,12,14,16)]))/(rowMeans(df_srna[human_hybrid,c(8,9)]))
o1=order(t2,decreasing=T)
for(i in human_hybrid[o1]){
  fc=round(t2[i],2)
  fun_t(i,fc)
}
dev.off()
pdf("../Updated_Figures_NEE/lendis.all_pachy.pdf",width=12,height=3,useDingbats = F)
par(mfrow=c(1,4),mar=c(2,4,3,1),tcl=0.3,cex=5/6,bty="n",xpd=F)
t2=(rowMeans(df_srna[human_pachy,c(6,11,15,17,10,12,14,16)]))/(rowMeans(df_srna[human_pachy,c(8,9)]))
o1=order(t2,decreasing=T)
for(i in human_pachy[o1]){
  fc=round(t2[i],2)
  fun_t(i,fc)
}
dev.off()
  # individual lendis ordered by FC for piRNA genes in all humans----
fun_t=function(pn,fc){
  tdf=read.table(paste("lendis/detailed/",pn,".merged.lendis",sep=""),check.names=F,header=T,row.names=1)
  mm=max(tdf)
  tcs=rep(cs_gg[c(4,1,2,3)],c(4,6,5,3))
  barplot(tdf[4:21,srna_rn[1]],space=0,border="white",col=cs_gg[4],names=18:35,
          ylab=paste(pn," FC=",fc,sep=""),main="",ylim=c(0,mm))
  for(i in c(2:4,13:18,8:12,5:7)){
    barplot(tdf[4:21,srna_rn[i]],space=0,border="white",col=tcs[i],names=18:35,
            ylab="",main="",ylim=c(0,mm))
  }
}
pdf("../Updated_Figures_NEE/lendis.all_prepachy.inAllSamples.pdf",width=48,height=3,useDingbats = F)
par(mfrow=c(1,18),mar=c(2,4,3,1),tcl=0.3,cex=5/6,bty="n",xpd=F)
t2=(rowMeans(df_srna[human_prepachy,c(6,11,15,17,10,12,14,16)]))/(rowMeans(df_srna[human_prepachy,c(8,9)]))
o1=order(t2,decreasing=T)
for(i in human_prepachy[o1]){
  fc=round(t2[i],2)
  fun_t(i,fc)
}
dev.off()
pdf("../Updated_Figures_NEE/lendis.all_hybrid.inAllSamples.pdf",width=48,height=3,useDingbats = F)
par(mfrow=c(1,18),mar=c(2,4,3,1),tcl=0.3,cex=5/6,bty="n",xpd=F)
t2=(rowMeans(df_srna[human_hybrid,c(6,11,15,17,10,12,14,16)]))/(rowMeans(df_srna[human_hybrid,c(8,9)]))
o1=order(t2,decreasing=T)
for(i in human_hybrid[o1]){
  fc=round(t2[i],2)
  fun_t(i,fc)
}
dev.off()
pdf("../Updated_Figures_NEE/lendis.all_pachy.inAllSamples.pdf",width=48,height=3,useDingbats = F)
par(mfrow=c(1,18),mar=c(2,4,3,1),tcl=0.3,cex=5/6,bty="n",xpd=F)
t2=(rowMeans(df_srna[human_pachy,c(6,11,15,17,10,12,14,16)]))/(rowMeans(df_srna[human_pachy,c(8,9)]))
o1=order(t2,decreasing=T)
for(i in human_pachy[o1]){
  fc=round(t2[i],2)
  fun_t(i,fc)
}
dev.off()

