setwd("~/github/CNEA-AQ/megan/src/db")

#Original Files
EF=read.csv("./orig/EFv210806.csv", sep=",")  

SpCrop=read.csv("./orig/SpeciationCrop210806.csv" , sep=",")  
SpHerb=read.csv("./orig/SpeciationHerb210806.csv" , sep=",")
SpShru=read.csv("./orig/SpeciationShrub210806.csv", sep=",")  
SpTree=read.csv("./orig/SpeciationTree210725.csv" , sep=","); SpTree=SpTree[,c(1:3)]  
#----------------------------------------------------------
SpCrop$Gtyp="Crop"
SpHerb$Gtyp="Herb"
SpShru$Gtyp="Shrub"
SpTree$Gtyp="Tree"
SpColumnNames=c("EcotypeID","VegID","SpFrac","Gtyp")

colnames(SpCrop)=SpColumnNames
colnames(SpHerb)=SpColumnNames
colnames(SpShru)=SpColumnNames
colnames(SpTree)=SpColumnNames

#Create new Speciation{Ntr,Btr}.csv files
df=merge(EF[c("VegID","GrowthForm")],SpTree, by="VegID",all.y = T)              #Add growthForm column

   #Needleleaf trees
   SpNtr=subset(df,GrowthForm=="Ntr")                                           #Get only Ntr species.
   NtEtFracs=aggregate(x = SpNtr$SpFrac, by=list(et=SpNtr$EcotypeID),FUN=sum)   #Compute partial-fraction: Ntr/tree
   for (i in c(1:nrow(NtEtFracs))){
     k=max(NtEtFracs$x[i],1e-23)
     et=NtEtFracs$et[i]
     SpNtr$SpFrac[which(SpNtr$EcotypeID==et)]=SpNtr$SpFrac[which(SpNtr$EcotypeID==et)]/k
   }
   #Broadleaf trees
   SpBtr=subset(df,GrowthForm=="Btr")                                           #Get only Btr species.
   BtEtFracs=aggregate(SpBtr$SpFrac,by=list(et=SpBtr$EcotypeID),FUN=sum)        #Compute partial fraction: Btr/Tree
   for (i in c(1:nrow(NtEtFracs))){
     k =max(1e-23,BtEtFracs$x[i])   
     et=BtEtFracs$et[i]   
     SpBtr$SpFrac[which(SpBtr$EcotypeID==et)]=SpBtr$SpFrac[which(SpBtr$EcotypeID==et)]/k   
   }
   ##I can write the new files:
   #write.csv(SpNtr[c("EcotypeID","VegID","SpFrac")],"./SpeciationNtree.csv")
   #write.csv(SpBtr[c("EcotypeID","VegID","SpFrac")],"./SpeciationBtree.csv")
   
SpNtr=SpNtr[SpColumnNames]; SpNtr$Gtyp="Ntr"
SpBtr=SpBtr[SpColumnNames]; SpBtr$Gtyp="Btr"
# All in one File:
   
df_list=list(crop=SpCrop,herb=SpHerb,shrub=SpShru,neddle=SpNtr,broad=SpBtr)
SpDF=do.call(rbind,df_list)
rownames(SpDF)=NULL
write.csv(SpDF,"SpeciationAll.csv")
rm(list=ls())
#===============================================================================
#
#Calculo el EF de cada "growth-type" (GT) agrupado por Ecotype

#Original Files
EF=read.csv("./orig/EFv210806.csv", sep=",");
SP=read.csv("./SpeciationAll.csv" , sep=",");

#Merge by VegID
M=merge(EF,SP,by = "VegID")
#Select only used columns
M[c("References","Comment","X","Family","GenusGroup","CommonName","Type","GrowthForm","VegID")]=NULL
#Multiply EF by Fraction and then Sum
efs=c(paste("EF",c(1:19),sep=""),paste("LDF",c(3:6),sep=""))
M[efs]=M[efs]*M$SpFrac

DF=aggregate(M[efs],by=list(gt=M$Gtyp,et=M$EcotypeID),FUN=sum)
DF=DF[order(DF$gt),]
DF=format(DF,digits=3,nsmall=3)
write.table(DF,"GtEFbyEcotype_new.csv",sep=",",dec=".",quote = F,row.names = F,col.names = F)
  
rm(list=ls())