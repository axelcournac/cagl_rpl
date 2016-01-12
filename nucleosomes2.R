# script to display nuclesomes pattern around a set of positions in Glabrata 
# 2: double graphs 

require(fields)
require('matrixStats')

# nucleosome data set:
d=read.table("GSM552914_Cgla_CM_Dec0908_nucCount_100.wig5")

# set of positions of interest:
#e= read.table("early_ori_GF.dat3")
#e= read.table("ARS_GF.dat4")
#e= read.table("late_ori_GF.dat3")
#e= read.table("mega_sat.dat3")

#e= read.table("genes.dat5")
#e= read.table("Cagl.txt5")
#e= read.table("Cagl.txt35")

#e= read.table("manip2_acs.pos3");
#e=subset(e,e$V4>6)

e= read.table("/run/media/axel/RSG3/BACK_UP/data/yeasts_species_project/glabrata/bio_data/GF_march2015/bonnass_acs.pos_sorted3");

bin_nuc = 20;   # bin of the nucleosome signal
area = 800 ;    # lenght of the area -+ around a ori

# Put in memory with lists of matrices :
d2=list()
for(i in 1:13) {print(i);d2[[i]]=subset(d,d$V1==i);}

M=matrix();
for(i in 1:dim(e)[1]  )
#for(i in 1:100  )
{
print(i);  
c=e[i,1];
#e[i,2] = floor(runif(1) * dim(d2[[c]])[1] ); 
ei = subset(d2[[c]]$V3,e[i,2]-area <= d2[[c]]$V2 & d2[[c]]$V2 <= e[i,2]+area);

if(i==1) {M=ei;}
else {M=rbind(M,ei);}
}


image.plot(t(M));

# reorderning the columns : 
h=heatmap(M,Colv=NA)


par(mfrow=c(2,1) )
# layout(matrix(c(1,1,2,3), 2, 1, byrow = TRUE), widths=c(1,1.5), heights=c(4,1))

#image.plot(t(M[h$rowInd,]**0.5));
h=heatmap(M,Colv=NA)

m=colMeans(M, na.rm = FALSE, dims = 1)

plot(m,type="l",lwd=3,xlab="Position (in bp)",ylab="Nucleosome Density");
points(800,30,type="h",col="green");

