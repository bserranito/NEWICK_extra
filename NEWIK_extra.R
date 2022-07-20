# extraction de Données NEWIK
library(ape)
library(data.table)
library(tcltk2)
library(dplyr)
library(geiger)
library(gplots)
library(metaSEM)
library(plyr)
library(reshape)
library(xlsx)


rm(list=ls())

choix=choose.files()


d = read.tree(file=choix)
d<- di2multi(d, 0.00009)


Edge=d$edge



# On extrait les premiers caracteres des noms
Names=d$tip.label
Names2=gsub("_.*", "", Names,perl=TRUE)

#Names2=substr(Names, start = 1, stop = 3)


# On identifie les niveaux terminaux
n <- Ntip(d)
Nods=d$edge[which(d$edge[, 2] <= n), 1]
Nods.ter=d$edge[which(d$edge[, 2] <= n), 2]

# ... et les autres niveaux
o_Nods=d$edge[which(d$edge[, 2] > n),]


#longueur des branches N2 inférieures au seuil
#L_Nods=d$edge.length[which(Nods)]


# On recherche les n replicats ...
df=data.frame(Label=Names2,Nods.ter,Nods)
dt = data.table(df)

#Replicats du premier degré
D1=df[duplicated(df$Nods) | duplicated(df$Nods, fromLast=TRUE),] 
tN1=as.data.frame(table(D1$Nods))
tN1_un=unique(tN1$Freq)


# Extraction des faux degré 1 ---------------------------------------------

# Evaluer si valeur de noeud est reliée a un niveau terminal

for (i in 1:dim(tN1)[1]){
  aa=as.numeric(levels(tN1$Var1))[i]
  ifelse(aa %in% o_Nods[,1],tN1$ev[i]<-c('N'),tN1$ev[i]<-c('O'))
}

tN1_r=(tN1[which(tN1$ev %in% c("O")),])
tN1_r[] <- lapply(tN1_r, function(x) if(is.factor(x)) factor(x) else x)


# Boucle de selection Degré 1-----------------------------------------------------

mat_NL=matrix(nrow=dim(tN1_r)[1], ncol=2)
N1=list()
T_N1=list()

for (j in 1:dim(tN1_r)[1]){  #Nombre de Dichotomie et Trichotomie ...
  a=as.numeric(levels(tN1_r$Var1))[j]
  grp=as.character(sort(factor(df[which(df$Nods==a),1])))
  
  mat_NL[j,2]=a
  mat_NL[j,1]=paste(grp, collapse='*')

  hh=as.matrix(combn(grp,2, FUN=paste, collapse='*'))
  N1<-rbind(N1,hh) #Incrémentation dans temp
}

# Boucle selection degré 2 ------------------------------------------------

# Récupération des numéros de Noeuds répliqués du premier degré
mat_NL=data.frame(mat_NL)
colnames(mat_NL)<- c("Label","Nods")
mat_NL$Nods.2=Edge[which(Edge[,2] %in% mat_NL$Nods),1]


# On recherche les niveaux terminaux non repliqué au premier degré
occ <- data.frame(table(df$Nods))
n_occ<-occ[occ$Freq ==1,1]
mat_n1nr=df[which(df$Nods %in% n_occ),]

on1nr=df[which(df$Nods %in% tN1$Var1[which(tN1$ev==c("N"))]),]
N1NR=rbind(mat_n1nr,on1nr)

#Recherche des valeurs communes chez les deux matrices
u_no=unique(N1NR$Nods)
mmm=intersect(mat_NL$Nods.2, u_no)

N2_1=list()
for (k in 1:length(mmm)){
  x=mmm[k]
 # m1=paste0("(",mat_NL$Label[which(mat_NL$Nods.2 %in% x)],")")
  m1.2=as.matrix(mat_NL$Label[which(mat_NL$Nods.2 %in% x)])
  m2=as.matrix(N1NR$Label[which(N1NR$Nods %in% x)])
  M=rbind(m1.2,m2)
  N2_temp=as.matrix(combn(M,2, FUN=paste, collapse='*'))
  N2_1=rbind(N2_1,N2_temp)
}

# Connection entre deux groupes de niveau 1
#Replicats du second degré 
D2=mat_NL[duplicated(mat_NL$Nods.2) | duplicated(mat_NL$Nods.2, fromLast=TRUE),] 

if (dim(D2)[1]==0){
  N2=N2_1
}else{
  tN2=as.data.frame(table(D2$Nods.2))
  
  N2_2=list()
  for (z in 1:dim(tN2)[1]){  #Nombre de Dichotomie et Trichotomie ...
    a2=as.numeric(levels(tN2$Var1))[z]
    grp2=as.character(sort(factor(mat_NL[which(mat_NL$Nods.2==a2),1])))
  #####################  grp2=paste0("(",grp2,")")
    # mat_NL[j,2]=a
    #  mat_NL[j,1]=paste(grp, collapse='*')
    
    hh2=as.matrix(combn(grp2,2, FUN=paste, collapse='*'))
    N2_2<-rbind(N2_2,hh2) #Incrémentation dans temp
  }
  N2=rbind(N2_1,N2_2)
}


# Extraction des données --------------------------------------------------

# Calcule des fréquences
Niveau1 <- matrix(unlist(N1), byrow = TRUE)
Niveau2 <- matrix(unlist(N2), byrow = TRUE)


N1_table=as.data.frame(sort(table(Niveau1)))
N2_table=as.data.frame(sort(table(Niveau2)))


#print(N1_table)
#print(N2_table)

#tk_choose.dir(default = "", caption = "Selection du dossier")

#Export en csv
#sink('output.csv')
##cat('Niveau 1')
#cat('\n')
#write.csv(N1_table)
#cat('____________________________')
#cat('\n')
#cat('Niveau 2')
#cat('\n')
#write.csv(N2_table)

#sink()

# Table de contingence----------------------------------------------------------
# Niveau 1 ----------------------------------------------------------------


temp=strsplit(as.character(N1), "*", fixed=TRUE)
mat_temp  <- matrix(unlist(temp), ncol=2, byrow=TRUE)

#rr=as.data.frame(cbind(mat_temp,N1_table$Freq))
rr=as.data.frame(mat_temp)

N1c_sort=rr[order(rr$V1,rr$V2 ),]

N1.df=count(N1c_sort, c('V1', 'V2')) 

un1.1 <- unique(unlist(N1.df[1:2]))
N1.df[1:2] <- lapply(N1.df[1:2], factor, levels=un1.1)



smat1 <- cast(N1.df, V1 ~ V2,
              fill=0, add.missing=TRUE,
              drop=TRUE)

Smat_fin1=subset(smat1,select=-c(V1))
rownames(Smat_fin1)= smat1$V1
Smat_fin1=as.matrix(Smat_fin1)

N1_Smat=as.matrix(Smat_fin1+t(Smat_fin1))
#diag(N1_Smat)<-0   #On remplace la diagonale par 0

colnames(N1_Smat)<- smat1$V1
rownames(N1_Smat)<- smat1$V1

# On remplace la diagonale par des zéros
N1_Smat[ row(N1_Smat) == col(N1_Smat) ] <- 0


Prop.N1=sweep(N1_Smat, 2, colSums(N1_Smat), FUN="/")

#my_palette <- colorRampPalette(c("black","orange", "red"))(n = 10)
#heatmap.2(N1_Smat, col = my_palette,
#          Rowv = NA,Colv = NA,  trace="none", keysize=0.75, 
#          key.par = list(cex=0.5))


# Niveau 2 ----------------------------------------------------------------

temp=strsplit(as.character(N2), "*", fixed=TRUE)
N2_couple=list()
for (z in 1:length(temp)){
  temp22=t(as.data.frame(unlist(temp[z])))
           
  d1 <- combn(temp22,2,simplify=FALSE)
  d1.1=as.data.frame(do.call(rbind, d1))
#  rownames(d1)<-NULL
  
  N2_couple=rbind(N2_couple,d1.1)
}

N2c_sort=N2_couple[order(N2_couple$V1,N2_couple$V2 ),]

N2.df=count(N2c_sort, c('V1', 'V2')) 

un1 <- unique(unlist(N2.df[1:2]))
N2.df[1:2] <- lapply(N2.df[1:2], factor, levels=un1)


#Pour obtenir une matrice symmetrique
smat2 <- cast(N2.df, V1 ~ V2,
                       fill=0, add.missing=TRUE,
                       drop=TRUE)

Smat_fin2=subset(smat2,select=-c(V1))
rownames(Smat_fin2)= smat2$V1
Smat_fin2=as.matrix(Smat_fin2)

N2_Smat=as.matrix(Smat_fin2+t(Smat_fin2))
#diag(N2_Smat)<-0   #On remplace la diagonale par 0

colnames(N2_Smat)<- smat2$V1
rownames(N2_Smat)<- smat2$V1

# On remplace la diagonale par des zéros
N2_Smat[ row(N2_Smat) == col(N2_Smat) ] <- 0

Prop.N2=sweep(N2_Smat, 2, colSums(N2_Smat), FUN="/")

#my_palette <- colorRampPalette(c("black","orange", "red"))(n = 40)
#col_breaks = c(seq(0,0,length=1),  # for red
#               seq(0,max(N2_Smat)/2,length=3),           # for yellow
#               seq(max(N2_Smat)/2,max(N2_Smat),length=3)) 

#heatmap.2(N2_Smat, col = my_palette,
#          Rowv = NA,Colv = NA,  trace="none", keysize=0.75, 
#          key.par = list(cex=0.5))


#  Export dans Excel ------------------------------------------------------

# Table
write.xlsx(N1_table, file="Output.xlsx", sheetName="N1_table", row.names=TRUE)
write.xlsx(N2_table, file="Output.xlsx", sheetName="N2_table", append=TRUE, row.names=TRUE)

# Matrix
write.xlsx(N1_Smat, file="Output.xlsx", sheetName="N1_matrix", append=TRUE, row.names=TRUE)
write.xlsx(N2_Smat, file="Output.xlsx", sheetName="N2_matrix", append=TRUE, row.names=TRUE)

#Prop Matrix
write.xlsx(Prop.N1, file="Output.xlsx", sheetName="N1_prop", append=TRUE, row.names=TRUE)
write.xlsx(Prop.N2, file="Output.xlsx", sheetName="N2_prop", append=TRUE, row.names=TRUE)


#Export en csv
#sink('output2.csv')
#cat('Niveau 1')
#cat('\n')
#write.csv(N1_Smat)
#cat('____________________________')
#cat('\n')
#cat('Niveau 2')
#cat('\n')
#write.csv(N2_Smat)
#sink()

