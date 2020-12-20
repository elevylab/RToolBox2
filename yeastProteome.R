###
### General functions to load function yeast data as well as to benchmark data
### 
library(org.Sc.sgd.db)
library(Biostrings)

##
## Return a matrix where yeast proteins are annotated according to a particular physical position
##
get.REF = function(){

  allGenes = read.csv("/data/elevy/70_R_Data/ALLGENES_R2.csv", stringsAsFactors=FALSE)
   tmp = read.csv("/data/elevy/70_R_Data/DHFR-LibraryOligosPosList.csv")
  allGenes$PCR = tmp$Fragment.size
  
  ## fakes the 16th plate
  fake.16 = cbind( paste( rep(NA, 384),1:384, sep=""), rep(NA, 384), rep(16,384), allGenes[ 1:384, c(4:6)], rep(0,384), rep(0,384), rep(0,384), rep(0,384),rep(0,384), rep(0,384), rep(0,384), rep(0,384) , rep(0,384) )
  colnames(fake.16) = colnames(allGenes)
  allGenes = rbind(allGenes, fake.16)
  imageJ.order = read.table("/data/elevy/62_ADE4/bin/imageJ.order", header=TRUE)
  new.order = order(imageJ.order$plate, imageJ.order$cadran, imageJ.order$pos)
  all.Genes.order = allGenes
  all.Genes.order[new.order,] = all.Genes.order
  REF = cbind(all.Genes.order,imageJ.order)
  REF[ which(REF[,1]==""),1]= paste( rep("Empty",4), 1:4, sep="")
  colnames(REF)[1]= c("ORF")
  REF$Desc = as.vector(unlist(sapply(as.list(REF[,1]), mget, env=org.Sc.sgdDESCRIPTION, ifnotfound=NA)))
  REF$row2 = rep(NA,6144)
  REF$col2 = rep(NA,6144)

  REF$row2[ REF$plate %in% c(1,3,9,11) ] =     REF$row[REF$plate %in% c(1,3,9,11)]*2
  REF$row2[ REF$plate %in% c(2,4,10,12) ] = 1+ (REF$row[REF$plate %in% c(2,4,10,12)]*2)

  REF$row2[ REF$plate %in% c(5,7,13,15) ] =     REF$row[REF$plate %in% c(5,7,13,15)]*2
  REF$row2[ REF$plate %in% c(6,8,14,16) ] = 1+ (REF$row[REF$plate %in% c(6,8,14,16)]*2)

  REF$col2[ REF$plate %in% c(1,3,9,11) ] =     REF$col[REF$plate %in% c(1,3,9,11)]*2
  REF$col2[ REF$plate %in% c(2,4,10,12) ] = 1+ (REF$col[REF$plate %in% c(2,4,10,12)]*2)

  REF$col2[ REF$plate %in% c(5,7,13,15) ] =     REF$col[REF$plate %in% c(5,7,13,15)]*2
  REF$col2[ REF$plate %in% c(6,8,14,16) ] = 1+ (REF$col[REF$plate %in% c(6,8,14,16)]*2)

  REF$P1536 = 0
  REF$P1536[ REF$plate %in% c(1,3,9,11) ] = 1
  REF$P1536[ REF$plate %in% c(2,4,10,12) ] = 2
  REF$P1536[ REF$plate %in% c(5,7,13,15) ] = 3
  REF$P1536[ REF$plate %in% c(6,8,14,16) ] = 4
  
  REF$ord=1:6144
  REF$Prot.Name = as.vector(unlist(sapply(as.list(REF[,1]), mget, env=org.Sc.sgdGENENAME, ifnotfound=NA)))
  REF$Prot.Name[is.na(REF$Prot.Name)]=REF$ORF[is.na(REF$Prot.Name)]
  REF$SGD = as.vector(unlist(sapply(as.list(REF[,1]), mget, env=org.Sc.sgdSGD, ifnotfound=NA)))
  
  return(REF)
}

##
## Return a matrix where yeast proteins are annotated according to a particular physical position
##
get.REF.full = function(){

  allGenes = read.csv("/data/elevy/70_R_Data/MichnickLablist_DHFRPCAcollectionList.csv", stringsAsFactors=FALSE)
  #tmp = read.csv("/data/elevy/70_R_Data/MichnickLablist_DHFRPCAcollectionList.csv")
  allGenes$PCR = tmp$Fragment.size
  
  ## fakes the 16th plate
  fake.16 = cbind( paste( rep(NA, 384),1:384, sep=""), rep(NA, 384), rep(16,384), allGenes[ 1:384, c(4:6)], rep(0,384), rep(0,384), rep(0,384), rep(0,384),rep(0,384), rep(0,384), rep(0,384), rep(0,384) , rep(0,384)  , rep(0,384)  , rep(0,384)  , rep(0,384)  , rep(0,384)  , rep(0,384)  , rep(0,384)  , rep(0,384)  , rep(0,384) )
  colnames(fake.16) = colnames(allGenes)
  allGenes = rbind(allGenes, fake.16)
  imageJ.order = read.table("/data/elevy/62_ADE4/bin/imageJ.order", header=TRUE)
  new.order = order(imageJ.order$plate, imageJ.order$cadran, imageJ.order$pos)
  all.Genes.order = allGenes
  all.Genes.order[new.order,] = all.Genes.order
  REF = cbind(all.Genes.order,imageJ.order)
  REF[ which(REF[,1]==""),1]= paste( rep("Empty",4), 1:4, sep="")
  colnames(REF)[1]= c("ORF")
  REF$Desc = as.vector(unlist(sapply(as.list(REF[,1]), mget, env=org.Sc.sgdDESCRIPTION, ifnotfound=NA)))
  REF$row2 = rep(NA,6144)
  REF$col2 = rep(NA,6144)

  REF$row2[ REF$plate %in% c(1,3,9,11) ] =     REF$row[REF$plate %in% c(1,3,9,11)]*2
  REF$row2[ REF$plate %in% c(2,4,10,12) ] = 1+ (REF$row[REF$plate %in% c(2,4,10,12)]*2)

  REF$row2[ REF$plate %in% c(5,7,13,15) ] =     REF$row[REF$plate %in% c(5,7,13,15)]*2
  REF$row2[ REF$plate %in% c(6,8,14,16) ] = 1+ (REF$row[REF$plate %in% c(6,8,14,16)]*2)

  REF$col2[ REF$plate %in% c(1,3,9,11) ] =     REF$col[REF$plate %in% c(1,3,9,11)]*2
  REF$col2[ REF$plate %in% c(2,4,10,12) ] = 1+ (REF$col[REF$plate %in% c(2,4,10,12)]*2)

  REF$col2[ REF$plate %in% c(5,7,13,15) ] =     REF$col[REF$plate %in% c(5,7,13,15)]*2
  REF$col2[ REF$plate %in% c(6,8,14,16) ] = 1+ (REF$col[REF$plate %in% c(6,8,14,16)]*2)

  REF$ord=1:6144
  REF$Prot.Name = as.vector(unlist(sapply(as.list(REF[,1]), mget, env=org.Sc.sgdGENENAME, ifnotfound=NA)))
  REF$SGD = as.vector(unlist(sapply(as.list(REF[,1]), mget, env=org.Sc.sgdSGD, ifnotfound=NA)))
  
  return(REF)
}


###
### Creates a matrix containing many descriptors of the yeast proteome.
###
get.proteome.table = function(REF=get.REF()){

  ORD = data.frame( ORF = as.vector(REF[,1]), ord=1:6144)
  
  SC.abund = read.csv("/data/elevy/70_R_Data/TABLES/sc_abund_christine.csv")[,1:16]
  colnames(SC.abund) = c("ORF","ab.apex.YPD","ab.apex.YMD","ab.GFP.YPD","ab.GFP.YMD","ab.western","ab.tap","mrna.sage","mrna.HDA","mrna.wang","mrna.av","cod.cai","cod.bias","cod.fop","pest.1","pest.2")
  SC.abund = SC.abund[,-c(11)]

  SC.tAI = read.table("/data/elevy/70_R_Data/TABLES/sc.tAI.values.txt", sep="\t", head=FALSE)
  colnames(SC.tAI) = c("ORF","cod.tAI")
  SC.tAI[,2]= round(SC.tAI[,2],3)
  
  SC.proteom = read.table("/data/elevy/48_subho_sc/data/Sc_expression_matrix.dat", sep="\t", header=TRUE)
  colnames(SC.proteom) = c("ORF","len","cod.tAI2","mrna.sde","ab.sde","mrna12","prot12","ribo.occu")
  SC.proteom = SC.proteom[,-c(3,4,5)]

  SC.addition = read.table("/data/elevy/70_R_Data/TABLES/sc.ABUND.paxDB.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)  
  colnames(SC.addition) = c("ORF","ab.pax")
  SC.addition[,2] = round(as.numeric(SC.addition[,2]),2)
  ## Subho has been using the rna HDA and the westerm for prots.
           
  SC.data = merge(SC.abund, SC.tAI, all.x=TRUE)
  SC = merge(SC.data,SC.proteom)
  SC = merge(SC,SC.addition, all.x=TRUE)

  SC.2.uniprot = read.csv("/data/elevy/70_R_Data/TABLES/yeast_2_UNIPROT.csv", sep="\t", stringsAsFactors=FALSE, header=FALSE)[,c(2,5)]
  colnames(SC.2.uniprot) = c("ORF","UPROT")  
  SC = merge(SC, SC.2.uniprot,all.x=TRUE)
    
  SC.maya = read.csv("/data/elevy/70_R_Data/TABLES/sc.ab.GFP.maya.csv")[,c(1,3,5,6,9,10,13,14,17)];
  colnames(SC.maya) = c("ORF","ab.GFP.maya.ymd","loc.ymd","ab.GFP.maya.dtt","loc.dtt","ab.GFP.maya.h2o2","loc.h2o2","ab.GFP.maya.starv","loc.starv")
  #index.to.replace = which(SC.maya[, grep("ab", colnames(SC.maya)) ] < 0)
  #SC.maya[index.to.replace, grep("ab", colnames(SC.maya)) ] = NA
  SC = merge(SC,SC.maya, all.x=TRUE)

  
  SC.dosage = read.csv("/data/elevy/70_R_Data/TABLES/sc.DOSAGE.sopko.txt", header=TRUE)
  SC.dosage = data.frame(ORF=  SC.dosage[,1], over.tox= 1 )
  SC = merge(SC,SC.dosage, all.x=TRUE)
  SC$over.tox[is.na(SC$over.tox)]=0

  ## intensity.H = Haploid??
  SC.mann = read.csv("/data/elevy/70_R_Data/TABLES/sc.MannAbund.csv")
  SC.mann.ins = SC.mann[,c(1,34,38,39)]
  colnames(SC.mann.ins)=c("ORF","ab.ms.ratio","ab.ms.haploid","ab.ms.diploid")
  SC = merge(SC,SC.mann.ins, all.x=TRUE)
  #SC$over.tox[is.na(SC$over.tox)]=0
  
  SC.diso = read.table("/data/elevy/70_R_Data/TABLES/sc.diso.txt", header=FALSE)
  colnames(SC.diso) = c("ORF", "len2", "diso05","diso1","diso2","sticky05","sticky1","sticky2","KRord","KRdes")
  SC = merge(SC,SC.diso, all.x=TRUE)

  SC.genom = read.csv("/data/elevy/70_R_Data/TABLES/sc.GENOMFEATURES.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
  colnames(SC.genom) = c("SGDID","ORF","chr","ch-start","ch-end","strand","chr2")
  SC.genom = SC.genom[,-c(1,7)]
  SC.genom$chr = gsub("chr", "", SC.genom$chr)
  SC = merge(SC, SC.genom, all.x=TRUE, by=c(1))

  #SC.phenom = read.csv("/data/elevy/70_R_Data/TABLES/sc.multidrug.table.csv")[,-c(2,6,8,9,10:88)];
  #colnames(SC.phenom) = c("ORF","essential","slow.hete","slow.homo","num.cond")
  #SC = merge(SC, SC.phenom, all.x=TRUE, by=c(1))  

  ### There are a few lines with B" and the " need to be deleted
  SC.phenom = read.csv("/data/elevy/70_R_Data/TABLES/sc.genotypes.SGD.131008.tsv", sep="\t");
  colnames(SC.phenom) = c("SGD.num", "ORF","PROT","PROT.synonym","ORF.qualif","exp.type", "disturbance.type", "readout.type", "readout.outcome", "notsure", "strain", "drug", "media","remark","remark2","pmid","publi")
  SC.phenom = SC.phenom[, -c(1,3,4,5,10,14,15)]

  ### Viable / inviable 
  SC.phenom1 = data.frame(ORF = SC.phenom$ORF[ SC.phenom[,9]==12140549 & SC.phenom[,8]==""& SC.phenom[,7]==""],
      viable = as.vector(SC.phenom$readout.type[ SC.phenom[,9]==12140549 & SC.phenom[,8]==""& SC.phenom[,7]==""]),
      stringsAsFactors=FALSE)
  SC.phenom1$viable[which(SC.phenom1$viable=="Viable")]=1
  SC.phenom1$viable[which(SC.phenom1$viable=="Inviable")]=0
  SC.phenom1$viable[SC.phenom1$viable != 1 & SC.phenom1$viable != 0]= NA

  SC = merge(SC, SC.phenom1, all.x=TRUE, by=c(1))
  
  ### slow growth Giaver 2002
  SC.phenom2 = table(SC.phenom$ORF[SC.phenom[,9]==12140549 & SC.phenom[,8]=="" & SC.phenom[,5]=="decreased"])
  SC.phenom2 = SC.phenom2[SC.phenom2 > 0]
  SC.phenom2 = data.frame( ORF = names(SC.phenom2), pheno.giav = SC.phenom2, stringsAsFactors=FALSE)

  SC = merge(SC, SC.phenom2, all.x=TRUE, by=c(1))
  
  ### slow growth Dudley / carbon source
  SC.phenom3 = table(SC.phenom$ORF[SC.phenom$pmid==16729036 & SC.phenom[,4]== "Utilization of carbon source"])
  SC.phenom3 = SC.phenom3[SC.phenom3 > 0]
  SC.phenom3 = data.frame( ORF = names(SC.phenom3), pheno.growth.dudley = SC.phenom3, stringsAsFactors=FALSE)

  SC = merge(SC, SC.phenom3, all.x=TRUE, by=c(1))

  ### slow growth Dudley / drug resistance
  SC.phenom4 = table(SC.phenom$ORF[SC.phenom$pmid==16729036 & SC.phenom[,4]== "Resistance to chemicals"])
  SC.phenom4 = SC.phenom4[SC.phenom4 > 0]
  SC.phenom4 = data.frame( ORF = names(SC.phenom4), pheno.resist.dudley = SC.phenom4, stringsAsFactors=FALSE)

  SC = merge(SC, SC.phenom4, all.x=TRUE, by=c(1))

  ### slow growth all --> pheno.growth.less / pheno.growth.more / pheno.growth.same
  SC.pheno.growth.increase = table(SC.phenom$ORF[SC.phenom[,4]=="Competitive fitness" & SC.phenom[,5]=="increased" & grepl("S288", SC.phenom[,6])] )
  SC.pheno.growth.decrease = table(SC.phenom$ORF[SC.phenom[,4]=="Competitive fitness" & SC.phenom[,5]=="decreased" & grepl("S288", SC.phenom[,6])] )
  SC.pheno.growth.normal   = table(SC.phenom$ORF[SC.phenom[,4]=="Competitive fitness" & SC.phenom[,5]=="normal"    & grepl("S288", SC.phenom[,6])] )

  SC.pheno.growth.increase = SC.pheno.growth.increase[SC.pheno.growth.increase>0]
  SC.pheno.growth.decrease = SC.pheno.growth.decrease[SC.pheno.growth.decrease>0]
  SC.pheno.growth.normal   = SC.pheno.growth.normal  [SC.pheno.growth.normal  >0]

  SC.pheno.growth.increase = data.frame(ORF = names(SC.pheno.growth.increase),  pheno.growth.up   = SC.pheno.growth.increase, stringsAsFactors=FALSE )
  SC.pheno.growth.decrease = data.frame(ORF = names(SC.pheno.growth.decrease),  pheno.growth.down = SC.pheno.growth.decrease, stringsAsFactors=FALSE )
  SC.pheno.growth.normal   = data.frame(ORF = names(SC.pheno.growth.normal  ),  pheno.growth.stab = SC.pheno.growth.normal  , stringsAsFactors=FALSE )

  SC = merge(SC, SC.pheno.growth.increase, all.x=TRUE, by=c(1))
  SC = merge(SC, SC.pheno.growth.decrease, all.x=TRUE, by=c(1))
  SC = merge(SC, SC.pheno.growth.normal, all.x=TRUE, by=c(1))
  
  ### drug resistance --> pheno.resist.less / pheno.resist.more / pheno.resist.same
  SC.pheno.resist.increase = table(SC.phenom$ORF[SC.phenom[,4]=="Resistance to chemicals" & SC.phenom[,5]=="increased" & grepl("S288", SC.phenom[,6])] )
  SC.pheno.resist.decrease = table(SC.phenom$ORF[SC.phenom[,4]=="Resistance to chemicals" & SC.phenom[,5]=="decreased" & grepl("S288", SC.phenom[,6])] )
  SC.pheno.resist.normal   = table(SC.phenom$ORF[SC.phenom[,4]=="Resistance to chemicals" & SC.phenom[,5]=="normal"    & grepl("S288", SC.phenom[,6])] )

  SC.pheno.resist.increase = SC.pheno.resist.increase[SC.pheno.resist.increase >0]
  SC.pheno.resist.decrease = SC.pheno.resist.decrease[SC.pheno.resist.decrease >0]
  SC.pheno.resist.normal   = SC.pheno.resist.normal  [SC.pheno.resist.normal   >0]
  
  SC.pheno.resist.increase = data.frame(ORF = names(SC.pheno.resist.increase),  pheno.resist.up   = SC.pheno.resist.increase , stringsAsFactors=FALSE)
  SC.pheno.resist.decrease = data.frame(ORF = names(SC.pheno.resist.decrease),  pheno.resist.down = SC.pheno.resist.decrease , stringsAsFactors=FALSE)
  SC.pheno.resist.normal   = data.frame(ORF = names(SC.pheno.resist.normal  ),  pheno.resist.stab = SC.pheno.resist.normal   , stringsAsFactors=FALSE)

  SC = merge(SC, SC.pheno.resist.increase, all.x=TRUE, by=c(1))
  SC = merge(SC, SC.pheno.resist.decrease, all.x=TRUE, by=c(1))
  SC = merge(SC, SC.pheno.resist.normal, all.x=TRUE, by=c(1))

  ### homozygote hap
  ##SC = merge(SC, SC.phenom, all.x=TRUE, by=c(1))  
  
  SC.ord = merge(ORD, SC, all.x=TRUE, by=c(1))  
  SC = SC.ord[ order(SC.ord$ord),]  
  return(SC)
}


## Loads all the yeast GO annotations
##
load.sc.GO = function(){
  go.table = read.csv("/data/elevy/70_R_Data/TABLES/sc.GO.table")[,c(2,3,5,7,10,12)]
  colnames(go.table)=c("ORF","Prot.name","name","type","evidence","desc")
  TYPE=rep(0,length(go.table[,1]))
  TYPE[go.table$type=="cellular_component"]="CC"
  TYPE[go.table$type=="biological_process"]="BP"
  TYPE[go.table$type=="molecular_function"]="MF"
  go.table$type=TYPE

  TYPE=rep(0,length(go.table[,1]))
  TYPE[go.table$evidence=="manually curated"]=100
  TYPE[go.table$evidence=="high-throughput"]=10
  TYPE[go.table$evidence=="computational"]=1
  go.table$evidence=TYPE  
  return(go.table)
}

## Format the object given by load.sc.GO, and create CC/BP/MF specific objects
## Min gives the minimum number of entries required per annotation to consider the annotation as valid.
##
format.GO = function(SC.GO, REF=get.REF(), type="CC", min=50){

  SC.GO.sub = SC.GO[SC.GO$type == type,]

  ### Retrieves all the names that have over X components
  GO.names = -sort(-table(SC.GO.sub$name))
  GO.names = names(GO.names[GO.names>min])

  GO.names = GO.names[which(GO.names != "cellular_component")]
  
  SC.GO.sub = SC.GO.sub[SC.GO.sub$name %in% GO.names,]

  SC.GO.sub = SC.GO.sub[SC.GO.sub$ORF %in% REF$ORF,]

  ### Final table:
  ### row=ORF, col=CCname, give=evidence code (0=nothing, 1=pred, 2=HT, 4=exp, 3=pred+HT, 5=pred+exp, 6=HT+exp, 7=all).
  ###
  SC.res = matrix(0, ncol=length(GO.names), nrow=NROW(REF))
  rownames(SC.res) = REF$ORF
  colnames(SC.res) = GO.names

  SC.GO.sub$name = as.vector(SC.GO.sub$name)
  SC.GO.sub$ORF = as.vector(SC.GO.sub$ORF)  
  
  for (i in 1:length(SC.GO.sub[,1])){
    if(i %% 1000 == 0){ 
      print(paste(i, "entries processed"))
    }
    SC.res[ SC.GO.sub$ORF[i], SC.GO.sub$name[i] ] = SC.res[ SC.GO.sub$ORF[i], SC.GO.sub$name[i] ] + SC.GO.sub$evidence[i]
  }
  SC.res.counts = colSums(SC.res>0)
  SC.res = SC.res[ , which(SC.res.counts>min)]
  return(SC.res)
}


###
### Loads localization data from Maya
###
get.loc.maya = function(REF){
    maya.loc = read.csv(file="/data/elevy/70_R_Data/MATRICES/sc.maya.localization.csv")
    colnames(maya.loc)[6] = "Control.Localization"
    maya.loc.aligned = merge(x=REF[,c("ORF","ord")], y=maya.loc, by="ORF", all.x=TRUE)
    maya.loc.aligned = maya.loc.aligned[ order(maya.loc.aligned$ord),]
    return(maya.loc.aligned)
}


## Loads PPI data.
## PCA PPI => 18467557
##
load.PPIs = function(REF = get.REF(), recalculate=FALSE, remove.PMID=c(), only=FALSE, MYFILE="/data/elevy/70_R_Data/TABLES/sc.BIOGRID.3.2.103.tab"){

  if( file.exists("/data/elevy/70_R_Data/PPI.mat.aligned") & !recalculate ){

    print(paste(date(), "Reading the file ..."))
    matrix.square = as.matrix(read.table(file="/data/elevy/70_R_Data/PPI.mat.aligned", header=TRUE, row.names=c(1), sep=" "))
    print(paste(date(), "Done"))
    
  } else {
    
    ORD = data.frame( ORF = as.vector(REF[,1]), ord=1:6144)
    matrix.square = matrix(ncol=6144,nrow=6144, 0)
    colnames(matrix.square) = REF[,1]
    rownames(matrix.square) = REF[,1]

    PPI = read.delim(file=MYFILE, header = TRUE, sep = "\t", quote="\"")
    PPI = PPI[ PPI$Experimental.System.Type == "physical",]
    PPI = PPI[ PPI$Experimental.System != "Affinity Capture-RNA",]

    ### 80802 interactions

    ###
    ### Here I want a high confidence dataset - first, order all interactions.
    ### 
    node.from = as.vector(PPI$Systematic.Name.Interactor.A)
    node.to   = as.vector(PPI$Systematic.Name.Interactor.B)
    node.PMID = as.vector(PPI$Pubmed.ID)

    #### 1 - sort by alphabetical order each interaction
    ft <- cbind(node.from, node.to)
    node.direction = apply( ft, 1, function(x) return(x[1]>x[2]))

    #### 2 - creates a list of directed interactions, i.e., each pair has always the same direction.
    node.from.ord = node.from
    node.to.ord   = node.to
    node.from.ord[node.direction] = node.to[node.direction]    # When TO is smaller, replace FROM with TO
    node.to.ord[node.direction]    = node.from[node.direction] #     and             replace TO with FROM --> HERE all INTs are ordered

    ft.ord <- cbind(node.from.ord, node.to.ord, node.PMID)

    ft.tmp = c()
    
    if( length(remove.PMID) > 0){

      for( PMID in remove.PMID){

        if(only){
          ft.tmp = rbind(ft.tmp, ft.ord[ which(ft.ord[,3]==PMID),])
        } else {
          ft.ord = ft.ord[ which(ft.ord[,3]!=PMID),]
        }
      }
    }

    if(only){
      ft.ord=ft.tmp
    }
    
    all.ints = apply( ft.ord, 1, function(x) paste(x, collapse="%") )
    all.ints = unique(all.ints)
    ft.ord.uniq = matrix(unlist(strsplit(all.ints,"%")), ncol=3, byrow=TRUE)

     # length(ft.ord.uniq) 73106

    ft.ppi = ft.ord.uniq[,1:2]
    all.ints.pairs = apply( ft.ppi, 1, function(x) paste(x, collapse="%") )
    ppi.table = table(all.ints.pairs)
    all.pairs = names(ppi.table)
    all.pairs.ppi = matrix(unlist(strsplit(all.pairs,"%")), ncol=2, byrow=TRUE)
    all.pairs.ppi = cbind(all.pairs.ppi, as.vector(ppi.table))
    
    table = all.pairs.ppi
    present.1 = table[,1] %in% REF[,1]
    present.2 = table[,2] %in% REF[,1]

    table = table[ present.1 & present.2,]

    col.i = data.frame( ORF=table[,1], pos.i = 1:(length(table[,1])) , sign= table[,3] )
    col.j = data.frame( ORF=table[,2], pos.j = 1:(length(table[,1])) , sign= table[,3] )

    col.i = merge(col.i, ORD);  col.i = col.i[ order(col.i$pos.i),]
    col.j = merge(col.j, ORD);  col.j = col.j[ order(col.j$pos.j),]

    matrix.square[ (col.i$ord-1)*length(ORD[,1]) + col.j$ord ]= as.numeric(as.vector(col.i$sign))
    matrix.square[ (col.j$ord-1)*length(ORD[,1]) + col.i$ord ]= as.numeric(as.vector(col.i$sign))

    if( length(remove.PMID)==0){
      write.table(matrix.square, file="/data/elevy/70_R_Data/PPI.mat.aligned", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=" ")
    }
  }
  return(matrix.square)
}


##
## Takes a SQUARE binary matrix with rows carrying labels --> returns a two column table
## e.g.,
## prot1 prot4
## prot1 prot10
## prot2 prot34
## ... etc
##
MAT2LIST = function(MAT){

  if(NROW(MAT) != NCOL(MAT)){
    print("Error, the matrix must be square")
    return()   
  } else {
    LABS = rownames(MAT)
    MAT = MAT*lower.tri(MAT)
    all.ints = c()    
    tmp=sapply(1:NCOL(MAT), function(x){ if( sum(MAT[x,]==1) > 0) { all.ints <<- rbind(all.ints, cbind(LABS[x], LABS[ MAT[x,]==1]) )} ;1 } )
  }
  return(all.ints)
}

##
## Takes a SQUARE binary matrix with rows carrying labels --> returns a list for all pairs > 0 and also gives associated value (number of papers)
##
MAT2LIST2 = function(MAT){

  if(NROW(MAT) != NCOL(MAT)){
    print("Error, the matrix must be square")
    return()   
  } else {
    LABS = rownames(MAT)
    MAT = MAT*lower.tri(MAT)
    all.ints = c()    
    tmp=sapply(1:NCOL(MAT), function(x){ if( sum(MAT[x,]>0) > 0) { all.ints <<- rbind(all.ints, cbind(LABS[x], LABS[ MAT[x,]>=1], MAT[x,MAT[x,]>=1] ) )} ;1 } )
  }
  return(all.ints)
}

##
## Takes a list carrying labels --> returns a SQUARE binary matrix
##
LIST2MAT = function(REF = get.REF(), LIST){

  vals=c()
  if(NCOL(LIST) == 3){
    vals=LIST[,3]
  } else {
    vals = rep(1,NROW(LIST))
  }
  
  LIST = LIST[ LIST[,1] %in% REF[,1] & LIST[,2] %in% REF[,1] ,]
  my.matrix = matrix(ncol=NROW(REF), nrow=NROW(REF),0)
  colnames(my.matrix) = REF[,1]
  rownames(my.matrix) = REF[,1]
  tmp = sapply(1:NROW(LIST), function(x){
    if(x %% 1000 == 0){print(x)};
    my.matrix[LIST[x,1], LIST[x,2]] <<- vals[x] ;
    my.matrix[LIST[x,2], LIST[x,1]] <<- vals[x] } )

  ## SETS all non defined lines/columns to NA
  #na.cols = rowSums(my.matrix)==0
  #my.matrix[na.cols,na.cols]=NA
    my.matrix = as.numeric(my.matrix)
  print(sum(my.matrix, na.rm=TRUE))
  return(my.matrix)
}



get.exp.cors = function(met=""){

  if( file.exists(paste("/data/elevy/70_R_Data/MATRICES/COR",met,".exp",sep="")) ){

    print(paste(date(), "Reading the exp.cor file ..."))
    rosetta.cor = as.matrix(read.table(file=paste("/data/elevy/70_R_Data/MATRICES/COR",met,".exp",sep=""), header=TRUE, row.names=c(1), sep=" "))
    print(paste(date(), "Done"))

  } else {

    ORD = data.frame( ORF = as.vector(get.REF()[,1]), ord=1:NROW(REF))
    
    rosetta.exp = read.table("/data/elevy/08_chipAnalysis/rosetta.mat")
    rosetta.names = read.table("/data/elevy/08_chipAnalysis/rosetta.names")
    rosetta.names = as.vector(rosetta.names[,1])
    rosetta.exp$ORF = rosetta.names
    rosetta.REF  = merge(ORD, rosetta.exp, all.x=TRUE)
    rosetta.REF2 = rosetta.REF[ order(rosetta.REF$ord), 3:302]
  
    print(paste(date(), "calculating correlations for expressions ..."))
    if(met==""){
      rosetta.cor = cor(t(rosetta.REF2), use="pairwise.complete.obs")
    } else {
      rosetta.cor = cor(t(rosetta.REF2), use="pairwise.complete.obs", method="spearman")
    }
    print(paste("Done", date()))
    diag(rosetta.cor)=0  
    colnames(rosetta.cor)=ORD$ORF #rosetta.names
    rownames(rosetta.cor)=ORD$ORF #rosetta.names

    rosetta.cor = round(rosetta.cor,4)
    
    write.table(rosetta.cor, file=paste("/data/elevy/70_R_Data/MATRICES/COR",met,".exp",sep=""), quote=FALSE, col.names=TRUE, row.names=TRUE, sep=" ")
  }
  return(rosetta.cor)
}


## Loads co-expression data
##
##
get.exp.cors.gasch = function(REF=get.REF(), met=""){

  if( file.exists(paste("/data/elevy/70_R_Data/MATRICES/COR",met,"gasch.exp",sep="")) ){

    print(paste(date(), "Reading the exp.cor file ..."))
    gasch.cor = as.matrix(read.table(file=paste("/data/elevy/70_R_Data/MATRICES/COR",met,"gasch.exp",sep=""), header=TRUE, row.names=c(1), sep=" "))
    print(paste(date(), "Done"))

  } else {
    
    ORD = data.frame( ORF = as.vector(REF[,1]), ord=1:6144)
    
    gasch.exp = read.csv("/data/elevy/66_YPCAN_ana/data/GASCH_data.csv")
    gasch.names = gasch.exp[,1]
    gasch.exp   = gasch.exp[, -c(1,2,3)]

    gasch.names = as.vector(gasch.names)
    gasch.exp$ORF = gasch.names
    gasch.REF  = merge(ORD, gasch.exp, all.x=TRUE)
    gasch.REF2 = gasch.REF[ order(gasch.REF$ord), 3:175]

    print(paste(date(), "calculating correlations for expressions ..."))
    if(met==""){
      gasch.cor = cor(t(gasch.REF2), use="pairwise.complete.obs")
    } else {
      gasch.cor = cor(t(gasch.REF2), use="pairwise.complete.obs", method="spearman")
    }
    print(paste("Done", date()))
    diag(gasch.cor)=0  
    colnames(gasch.cor)=ORD$ORF #gasch.names
    rownames(gasch.cor)=ORD$ORF #gasch.names

    gasch.cor = round(gasch.cor,4)
    
    write.table(gasch.cor, file=paste("/data/elevy/70_R_Data/MATRICES/COR",met,"gasch.exp",sep=""), quote=FALSE, col.names=TRUE, row.names=TRUE, sep=" ")
  }
  return(gasch.cor)
}


PPI.pmid.2.names = function(PMIDS){

  PPI.desc = get.PPI.studies.desc()
  my.names = c()

  for (each.pmid in PMIDS){

    if(each.pmid %in% PPI.desc[[2]]){

      pos = which(PPI.desc[[2]]==each.pmid)
      my.names = c(my.names, PPI.desc[[1]][pos])
    } else {
      my.names = c(my.names, each.pmid)
    }        
  }
  return(my.names)
}

get.PPI.studies.desc = function(){
  
  a1=c("Gilmore et al.","Starita et al.","Malsburg et al.", "Sahasranaman et al.", "Garcia Gomez et al.","Lee et al. 2011", "Schwer et al. 2011","Fasolo et al. 2011")
  p1=c("22199229","22106047","21944719","21926967","21825077", "21734642", "21558325","21460040")
  d1=c("Histone binding prots","ubiquitin attachment","mitofilin","pre-rRNA processing", "RNA helicase Spb4","SAGA and ADA complexes","snRNPs and snoRNPs", "kinase PPIs by prot microarray")
  
  a2=c("Ziv et al. 2011"," Ren et al. 2011", "Panasenko et al. 2011","Sammons et al. 2011", "Akiyoshi et al. 2010")
  p2=c("21427232","21386897","21321079", "21277287", "21107429")
  d2=c("ubiquitin","pre-mRNA splicing factors connected to Prp19","proteasome assembly and functional integrity", "Gis2 and translation machinery","kinetochore-microtubule attachments")
  
  a3=c("Ranjitkar et al. 2010", "Muller et al 2010", "Ossareh-Nazari et al. 2010","Breitkreutz et al. 2010","Kaake et al. 2010","Tonikian et al. 2009")
  p3=c("21070971", "20595233", "20508643","20489023","20170199","19841731")
  d3=c("ubiquitin and histone localization", "Orc1 and DNA replication origins", "new partners of the ubiquitin protease Ubp3","Signaling network","Cell cycle specific 26S interactions","SH3 domain interactome")
  
  a4=c("Gong et al. 2009","Lambert et al. 2009", "Guerrero et al. 2008", "Yu et al. 2008","Colomina et al.", "Tarassov et al. 2008", "Alber et al 2007","Oeffinger et al. 2007")
  p4=c("19536198", "19106085", "18757749", "18719252", "18667435", "18467557", "18046405", "17922018")
  d4=c("chaperone-protein interactions", "chromatin-associated protein network", "proteasome interaction network by QTAX", "HQ yeast network", "Whi3 with mRNAs in ER", "PCA Network", "Nuclear pore cplx", "diverse ribonucleoprotein complexes")
  
  a5=c("Wong et al. 2007", "Collins et al. 2007", "Hesselberth et al. 2006","Krogan et al. 2006","Gavin et al. 2006","Ptacek et al. 2005","Mo et al. 2005", "Miller et al. 2005")
  p5=c("17634282","17200106","16606443", "16554755","16429126","16319894","16300994", "16093310")
  d5=c("interaction map of the mitotic spindle", "Large HQ network by MS", "WW domains and their interacting proteins", "MS Network Krog","MS Network Gav","Kinase in vitro","sterol biosynthetic protein-protein interactions","membrane interactions")
  
  a6=c("Zhao et al. 2005", "Wohlschlegel et al. 2004", "Krogan et al. 2004","Graumann et al. 2004","Hitchcock et al. 2003","Peng et al. 2003", "Grandi et al. 2003","Sanders et al. 2002","Ho et al. 2002","Gavin et al. 2002", "Tong et al. 2002","Drees et al. 2001","Ito et al. 2001","Fromont-Racine et al. 2000","Uetz et al. 2000")
  p6=c("15766533", "15326169","14759368","14660704","14557538","12872131","12150911","12052880","11805837","11805826","11743162","11489916","11283351","10900456","10688190")
  d6=c("chaperone network", "protein sumoylation network", "yeast RNA-processing complexes", "TAP MudPIT for pathways","membrane, ubiquitin, ER", "protein ubiquitination","ribosome related","transcription machinery","Old MS Network 2","Old MS Network 1", "LS yeast","protein interaction map for cell polarity","Y2H interactome","Sm-like proteins network","Y2H first")

  a7=c("Wang et al., 2012", "Babu et al, 2012", "Kolawa et al, 2013")
  p7=c("22875988", "22940862", "23793018")
  d7=c("Coiled-coil interactions",  "Yeast membrane proteins", "Ubiquitin network")
  
  PPI.desc = list()
  PPI.desc[[1]] = c(a1,a2,a3,a4,a5,a6,a7)
  PPI.desc[[2]] = c(p1,p2,p3,p4,p5,p6,p7)
  PPI.desc[[3]] = c(d1,d2,d3,d4,d5,d6,d7)

  return(PPI.desc)
}



## Takes scores obtained for 2 sets of interactions,
## Return a ROC curve of how well they agree with each other.
##
##
get.FP.TP.NEW = function(bench.data, values, datapoints = c(seq(0.1,2.4,by=0.1),seq(2.5,10,by=0.5),11:20,seq(25,100,by=5) ) , rand=FALSE){
  
  TP = c()
  FP = c()
  ACC = c()
  PPV = c()
  REC = c()
  ROC1 =c()
  ROC2 =c()

  bench.data = bench.data[ !is.na(values),]
  values = values[!is.na(values)]

  L  = length(bench.data$i)
  L.pos = sum(bench.data$i==1)
  L.neg = sum(bench.data$i==0)
  
  reorder = sample(1:L, replace=FALSE)
  bench.data = bench.data[reorder ,]

  if(rand==FALSE){              # If random, values are not realigned to get a "randomized dataset".
    values = values[reorder]
  }

  bench.data = bench.data[order(-values),]

  for ( CUT in datapoints){
    
    X = as.integer(0.01 * CUT * L )
    
    TP = sum( bench.data$i[1:X] >  0)
    FP = sum( bench.data$i[1:X] == 0)

    FN = sum( bench.data$i[X:L] >  0)
    TN = sum( bench.data$i[X:L] == 0)

    ACC = c(ACC, (TP/L.pos) / ( (TP/L.pos)  +  (FP/L.neg) ))
    PPV = c(PPV, (TP / (TP  + FP) ) )
    REC = c(REC,CUT)
    #ROC1 = c(ROC1, TP/L1)
    #ROC2 = c(ROC2, FP/L2)
  }

  return( data.frame(REC, ACC, PPV))
}


## 
##
get.CO.info = function(CO.mat, ORF.1, ORF.2){

  ### Keeps only those ORFs that are defined in CO.mat.
  def.row = which(rowSums(is.na(CO.mat)) != NCOL(CO.mat))
  def.col = which(colSums(is.na(CO.mat)) != NROW(CO.mat))
  #print(paste("DEF COL = ",def.row))

  present.1 = ORF.1 %in% colnames(CO.mat)[def.row]
   present.2 = ORF.2 %in% colnames(CO.mat)[def.col]
  
  #present.1 = ORF.1 %in% colnames(CO.mat)[!is.na(CO.mat[,def.col]) & !is.na(CO.mat[def.col,]) ]
  #present.2 = ORF.2 %in% colnames(CO.mat)[!is.na(CO.mat[,def.col]) & !is.na(CO.mat[def.col,]) ]

  print( paste( sum( !present.1 | !present.2)," interactions out of",length(present.1)," are not in the CO.mat and were removed"))
  
  ORF.1b = ORF.1[ present.1 & present.2]
  ORF.2b = ORF.2[ present.1 & present.2]

  LOOP = 1:(length(ORF.1b))

  CORS.b = sapply(LOOP, function(x){ return (CO.mat[ ORF.1b[x], ORF.2b[x] ])})

  CORS = rep(NA,length(ORF.1))
  CORS[ present.1 & present.2 ] = CORS.b
  
  return(CORS)
}

get.CO.info.rank = function(CO.mat, ORF.1, ORF.2){

  ### Keeps only those ORFs that are defined in CO.mat.
  present.1 = ORF.1 %in% colnames(CO.mat) & !is.na(ORF.1) #[which(!is.na(CO.mat[,2449]) & !is.na(CO.mat[2449,]) )]
  present.2 = ORF.2 %in% colnames(CO.mat) & !is.na(ORF.2) #[which(!is.na(CO.mat[,2449]) & !is.na(CO.mat[2449,]) )]

  #print( paste( sum( !present.1 | !present.2)," interactions out of",length(present.1)," are not in the CO.mat and were removed"))
  
  ORF.1b = ORF.1[ present.1 & present.2]
  ORF.2b = ORF.2[ present.1 & present.2]

  LOOP = 1:(length(ORF.1b))

  CORS.b = sapply(LOOP, function(x){ if(ORF.1b[x] != ORF.2b[x]){ return (  min( CO.mat[ ORF.1b[x], ORF.2b[x] ], CO.mat[ ORF.2b[x], ORF.1b[x] ])  ) } else { return(NA) }})

  CORS = rep(NA,length(ORF.1))
  CORS[ present.1 & present.2 ] = CORS.b
  
  return(CORS)
}

## Takes a matrix of similarity and return the N partners corresponding to the highest scores in the matrix
##
## /!\ if matrix is symmetric, each number is returned twice.
get.COR.highest = function(CO.mat, N){

  diag(CO.mat)=NA
  
  n.def = sum(!is.na(CO.mat))
  Q = quantile(CO.mat, na.rm=TRUE, prob=c(0,1-(N/n.def),1))
  cut.off = Q[[2]]

  goods = which(CO.mat > cut.off)
  I = ceiling(goods / NROW(REF))
  J = goods %% NROW(REF)

  names = colnames(CO.mat)
  
  pairs = data.frame( ORF.1 = as.vector(names[I]), ORF.2 = as.vector(names[J]), val= as.vector(CO.mat[goods]) )

  return(pairs)
}


## Takes a matrix of similarity and return the closest partner for each entry in the form of: rbind(ENTRIES, CLOSESTS)
##
get.COR.closest = function(CO.mat){

  diag(CO.mat)=0
  closest = apply(CO.mat, 1, function(x){ return( which(x == max(x, na.rm=TRUE))[1])} )
  closest.val = apply(CO.mat, 1, function(x){ return( max(x, na.rm=TRUE))} )

  pairs = data.frame( ORF.1 = colnames(CO.mat), ORF.2 = colnames(CO.mat)[closest], val=closest.val )

  return(pairs)
}

## Takes a matrix of similarity and return the N closest partners for each entry in the form of: list`ENTRY` = c(CLOSESTS)
##
get.COR.closest.N = function(CO.mat, N, ORF.to.exclude=c()){
    
  diag(CO.mat)=NA
  ## Removes all lines where there are only NAs
  undefined = apply(CO.mat, 1, function(x){ return( sum(is.na(x)))})
  undefined = undefined == length(undefined)
  CO.mat = CO.mat[ !undefined, !undefined]
  ORFs.REF = colnames(CO.mat)

  closest = apply(CO.mat, 1, function(x){ return( ORFs.REF[order(-x)[1:(N)]] )})   #which(x == max(x, na.rm=TRUE))[1])} )
  if( length(closest) %% N != 0){
    print("Problem, with get.COR.closest.N")
    return()
  }

  closest.val = apply(CO.mat, 1, function(x){ return( x[order(-x)[1:(N)]] )})   #which(x == max(x, na.rm=TRUE))[1])} )
  if( length(closest) %% N != 0){
    print("Problem, with get.COR.closest.N")
    return()
  }

  ORFs.rep = rep( ORFs.REF, rep(N, length(closest)/N ))
  closest.2 = cbind(ORFs.rep, as.vector(closest), as.vector(closest.val))
  closest.2 = closest.2[order(- as.numeric(closest.2[,3])),]

  if(length(ORF.to.exclude)>0 ){
          
      closest.2 = closest.2[ which(!closest.2[,1] %in%  ORF.to.exclude & !closest.2[,2] %in%  ORF.to.exclude),]
          
  }
  
  #closest = rbind(ORFs.REF, closest)
  #closest = split(t(closest), rep(1:(length(ORFs.REF) ),N ) )
  #names(closest) = ORFs.REF
  
  return(closest.2)
}

## Measures the predictive power of MAT1 for MAT2.
##
benchmark.data = function(MAT1, MAT2, CUT=NA){
  
  ## First, creates a positive set
  ##
  if(is.na(CUT)){
    positive.set  = get.COR.closest.N(CO.mat=MAT2, N=6)
  } else {

    n.pos = sum(MAT2 >= CUT, na.rm=TRUE)
    positive.set = get.COR.highest(CO.mat=MAT2, N=n.pos)

  }
  benchmark.set = get.benchmark.from.ORFs.NEW(REF, MAT=MAT2, ORF.1 = positive.set[,1], ORF.2 = positive.set[,2], NEG.factor=1)

  ROC.data = get.FP.TP.NEW( bench.data=benchmark.set, values=get.mat.indexes(MAT= 1*MAT1, indexes = benchmark.set[,1:2]) )

  return(ROC.data)
}


## Takes a matrix of similarity and return the SAME matrix where each line/column is transformed into a rank
##
## byROW --> rank is by ROW /// bySIM = similarity matrix (as opposed to distance matrix).
##
get.COR.closest.rank = function(CO.mat, byROW=TRUE, bySIM=TRUE){

  diag(CO.mat)=NA
  
  if(byROW==TRUE){

    closest = t(apply(-CO.mat,1,rank))
    
  } else {

    closest = apply(-CO.mat,2,rank)
  }

  closest[is.na(CO.mat)]=NA

  colnames(closest)=colnames(CO.mat)
  rownames(closest)=rownames(CO.mat)
  
  return(closest)
}


## Takes a matrix of similarity and return the N closest partners for each entry in the form of: list`ENTRY` = c(CLOSESTS)
##
get.COR.closest.N.rand = function(CO.mat, N){

  diag(CO.mat)=0
  ## Removes all lines where there are only NAs
  undefined = apply(CO.mat, 1, function(x){ return( sum(is.na(x)))})
  undefined = undefined == length(undefined)
  CO.mat = CO.mat[ !undefined, !undefined]
  ORFs.REF = colnames(CO.mat)

  closest = apply(CO.mat, 1, function(x){ return( ORFs.REF[sample(1:(length(ORFs.REF)), N-1) ] )})   #which(x == max(x, na.rm=TRUE))[1])} )
  closest = rbind(ORFs.REF, closest)
  closest = split(t(closest), rep(1:(length(ORFs.REF) ),N ) )
  names(closest) = ORFs.REF
  
  return(closest)
}

##
## Returns a matrix CODON USAGE correlations
##
get.codon.info = function(remove.non.info = TRUE){

  ORD = data.frame( ORF = as.vector(get.REF()[,1]), ord=1:NROW(REF))

  cod.table = read.table("/data/elevy/51_stickyness_lastAna/results/sc_tAI_table.txt", sep="\t", header=TRUE)
  colnames(cod.table)[1] = "ORF"
  cod.table.tmp = matrix( ncol=64, nrow=length(cod.table[,1]),  as.numeric(as.matrix(cod.table[, c(2:65)]))  )
  cod.table.ord = merge(ORD,cod.table,all.x=TRUE,sort=FALSE)
  cod.table.ord = cod.table.ord[order(cod.table.ord$ord),]
  
  cod.table.f = matrix( ncol=64, nrow=NROW(REF),  as.numeric(as.matrix(cod.table.ord[, c(3:66)])) )
  colnames(cod.table.f) = colnames(cod.table.ord)[3:66]
    
  ## Codon frequency table
  cod.table.f.norm = round(100*cod.table.f / rowSums(cod.table.f),2)
  rownames(cod.table.f.norm) = as.vector(ORD$ORF)

  if( remove.non.info){
    indexes.to.remove = c(
      grep("S3",colnames(cod.table.f.norm)), grep("F1",colnames(cod.table.f.norm)),
      grep("Y1",colnames(cod.table.f.norm)),
      grep("Z",colnames(cod.table.f.norm)),
      grep("C",colnames(cod.table.f.norm)),
      grep("W",colnames(cod.table.f.norm)),
      grep("L[4-6]",colnames(cod.table.f.norm)),
      grep("P[2-3]",colnames(cod.table.f.norm)),
      grep("H",colnames(cod.table.f.norm)),
      grep("R[1-4]",colnames(cod.table.f.norm)),
      grep("I[2-3]",colnames(cod.table.f.norm)),
      grep("M1",colnames(cod.table.f.norm)),
      grep("T[2-3]",colnames(cod.table.f.norm)),
      grep("A[2-3]",colnames(cod.table.f.norm)),
      grep("G[1-3]",colnames(cod.table.f.norm))
      )
  } else {
    indexes.to.remove = c(
      grep("Z",colnames(cod.table.f.norm))
      )
  }

  cod.table.f.norm = cod.table.f.norm[ , - indexes.to.remove ]

  return(cod.table.f.norm)
}



###
### Returns sets of positively/negatively correlated genes with indexes aligned with the REF.
###
get.benchmark.data.EXP = function(CO.mat= CO.exp, N.closest=10, MIN.info = 5, tAI="all"){

  rosetta.exp = read.table("/data/elevy/08_chipAnalysis/rosetta.mat")
  rosetta.names = read.table("/data/elevy/08_chipAnalysis/rosetta.names")
  rosetta.names = as.vector(rosetta.names[,1])
  rosetta.exp$ORF = rosetta.names

  rosetta.REF  = merge(REF, rosetta.exp, all.x=TRUE)
  rosetta.REF2 = rosetta.REF[ order(rosetta.REF$ord), 23:322]
  rownames(rosetta.REF2) = as.vector(rosetta.REF[ order(rosetta.REF$ord), 1])

  ##############################  co-expression
  
  sum.under.1 = apply( rosetta.REF2,1, function(x){ return(sum(x < 0.15 & x > -0.15,na.rm=TRUE))}) ## remove all proteins with no info
  sum.over.1  = apply( rosetta.REF2,1, function(x){ return(sum(x > 0.2 | x < -0.2,na.rm=TRUE))}) ## to keep all proteins with info
  little.info = which(sum.under.1 >= 295 | sum.under.1 == 0) ## do not change in over 295/300 conditions
  
  if(tAI == "all"){
    more.info   = which(sum.over.1 >= MIN.info)    ## do change in at least   5/300 conditions
  } else {
    
    if(tAI == "low"){
      more.info   = which(sum.over.1 >= MIN.info & SC$cod.tAI < 0.42)    ## do change in at least   5/300 conditions      
    } else {
      more.info   = which(sum.over.1 >= MIN.info & SC$cod.tAI > 0.42)    ## do change in at least   5/300 conditions      
    }
  }
  
  L = length(sum.over.1)
  
  ############################## Positive
  indexes.j = apply( CO.mat, 1, function(x){ return( (order(-x))[1:N.closest]    )} )
  indexes.i = rep(1:L, rep(N.closest,L ) )
    
  indexes.j2 = indexes.j[indexes.i %in% more.info & indexes.j %in% more.info & indexes.j > indexes.i]
  indexes.i2 = indexes.i[indexes.i %in% more.info & indexes.j %in% more.info & indexes.j > indexes.i]

  print(paste("There are ",length(indexes.i2)," protein pairs in the positive dataset"))
  #plot(density(CO.mat[(indexes.j2-1)*L+indexes.i2],na.rm=TRUE, from=-1, to=1,bw=0.05), lty=1)

  ############################## Negative
  indexes.j.neg = sample(indexes.j2,50000, replace=TRUE)
  indexes.i.neg = sample(indexes.i2,50000, replace=TRUE)

  indexes.j.neg2 = indexes.j.neg[ !is.na(CO.mat[(indexes.j.neg-1)*L+indexes.i.neg]) & CO.mat[(indexes.j.neg-1)*L+indexes.i.neg] > -0.15 & CO.mat[(indexes.j.neg-1)*L+indexes.i.neg] < 0.15]
  indexes.i.neg2 = indexes.i.neg[ !is.na(CO.mat[(indexes.j.neg-1)*L+indexes.i.neg]) & CO.mat[(indexes.j.neg-1)*L+indexes.i.neg] > -0.15 & CO.mat[(indexes.j.neg-1)*L+indexes.i.neg] < 0.15]

  print(paste("There are ",length(indexes.i.neg2)," protein pairs in the negative dataset"))
  
  #lines(density(CO.mat[(indexes.j.neg2-1)*L+indexes.i.neg2],na.rm=TRUE, from=-1, to=1,bw=0.05), lty=2)
  
  result=list()
  result[[1]]= data.frame(POS.i = indexes.i2, POS.j = indexes.j2)
  result[[2]]= data.frame(NEG.i = indexes.i.neg2, NEG.j = indexes.j.neg2)
  return(result)
}


get.CODON.interact.table = function(benchmark.data = benchmark.exp, cod.table.f.norm = get.codon.info(),
  my.breaks.mean = c(0,0.1,0.25,0.5,1,2,5,10,50),
  my.breaks.diff = c(0,0.1,0.25,0.5,1,2,5,10,50)
  ){

  #### 
  pos.diff = apply( benchmark.data[[1]] ,1, function(x){return(abs(cod.table.f.norm[x[1],] - cod.table.f.norm[x[2],]) )})
  pos.mean = apply( benchmark.data[[1]] ,1, function(x){return( apply( rbind( cod.table.f.norm[x[1],], cod.table.f.norm[x[2],]),2,max) )})
  
  neg.diff = apply( benchmark.data[[2]] ,1, function(x){return(abs(cod.table.f.norm[x[1],] - cod.table.f.norm[x[2],]) )})
  neg.mean = apply( benchmark.data[[2]] ,1, function(x){return( apply( rbind( cod.table.f.norm[x[1],], cod.table.f.norm[x[2],]),2,max) )})

  my.matrices = list()
  
  for(i in (1: (length(cod.table.f.norm[1,]) ) ) ){      ## For each codon

    print( colnames(cod.table.f.norm)[i])

    table.pos = table( cut( pos.mean[i,], breaks=my.breaks.mean),  cut( pos.diff[i,], breaks=my.breaks.diff))
    table.neg = table( cut( neg.mean[i,], breaks=my.breaks.mean),  cut( neg.diff[i,], breaks=my.breaks.diff))

    table.prob = round( ( table.pos/sum(table.pos) ) /  ( table.pos/sum(table.pos) + table.neg/sum(table.neg) ), 4)

    table.prob2 = table.prob
    table.prob2[ is.infinite(table.prob2) | is.na(table.prob2) | table.prob2==0] = 0.5
    
    table.prob[ is.infinite(table.prob)] = max(table.prob2, na.rm=TRUE)
    table.prob[ table.prob == 0] = min(table.prob2,na.rm=TRUE)
    table.prob[ is.na(table.prob)] = 0.5
          
    my.matrices[[i]] = table.prob
  }
  names(my.matrices) = colnames(cod.table.f.norm)
  return(my.matrices)
}

get.CODON.interact.bench.fast = function(my.matrices = my.mats.1, benchmark.data = benchmark.exp[[1]], cod.table.f.norm = get.codon.info(),

  my.breaks.mean = c(0,0.1,0.25,0.5,1,2,5,10,50),
  my.breaks.diff = c(0,0.1,0.25,0.5,1,2,5,10,50),
  LOG= FALSE
  ){

  print(date())
  diff = apply( benchmark.data ,1, function(x){return(abs(cod.table.f.norm[x[1],] - cod.table.f.norm[x[2],]) )})
  mean = apply( benchmark.data ,1, function(x){return( apply( rbind( cod.table.f.norm[x[1],], cod.table.f.norm[x[2],]),2,max) )})

  print(date())

  diff.cut = matrix( as.numeric(factor(cut(diff, my.breaks.diff))), dim(diff), byrow=FALSE)
  mean.cut = matrix( as.numeric(factor(cut(mean, my.breaks.mean))), dim(diff), byrow=FALSE)
  
  all.probs = c()

  all.matrices = array(dim=c( dim(my.matrices[[1]]), length(my.matrices) ) ) # ROW * COL * DEEP

  COD.NUMBER = length(cod.table.f.norm[1,])
  
  for ( COD.NUM in 1:COD.NUMBER  ){      ## For each codon
    all.matrices[,,COD.NUM] = my.matrices[[COD.NUM]]
    #all.matrices[0,0,COD.NUM]=0.5
  }

  LOOP = 1:(length(diff[1,]))
  print(date())

  COD.INDEXES = matrix( ncol=length(diff[1,]) , rep( seq(0, 64*(COD.NUMBER-1),by=64), length(diff[1,])), byrow=FALSE)
  COL.INDEXES = (diff.cut-1)*( length(my.breaks.mean) - 1)
  ROW.INDEXES = mean.cut
  ALL.INDEXES = COD.INDEXES + COL.INDEXES + ROW.INDEXES

  probs = list()
  
  probs[[1]] = sapply( LOOP, function(x){ val=all.matrices[ ALL.INDEXES[,x]]; val[is.na(val)]=0.5;  return( sum(log(val))) } )
  probs[[2]] = sapply( LOOP, function(x){ return( mean(    all.matrices[ ALL.INDEXES[,x] ]  , na.rm=TRUE)  ) })

  print(date())
  return(probs)
}

get.CODON.interact.one.gene = function(my.matrices = my.mats.1, gene.num=1, cod.table.f.norm = COD.freqs,
  my.breaks.mean = c(-0.1,0.1,0.25,0.5,1,2,5,10,50),
  my.breaks.diff = c(-0.1,0.1,0.25,0.5,1,2,5,10,50)
  ){

  #my.matrices = my.mats.COD.TABLE
  
  gene.table = matrix( ncol= length(cod.table.f.norm[1,]), rep( cod.table.f.norm[gene.num,], length(cod.table.f.norm[,1])), byrow=TRUE)
  diff.tmp = cod.table.f.norm-gene.table

  diff = abs(diff.tmp)
  mean = 0*diff.tmp  

  mean[ diff.tmp > 0 & !is.na(mean)] = cod.table.f.norm[ diff.tmp > 0  & !is.na(mean)]
  mean[ diff.tmp < 0 & !is.na(mean)] = gene.table[ diff.tmp < 0  & !is.na(mean)]
    
  diff.cut = matrix( as.numeric(factor(cut(diff, my.breaks.diff))), dim(diff), byrow=FALSE)
  mean.cut = matrix( as.numeric(factor(cut(mean, my.breaks.mean))), dim(diff), byrow=FALSE)
  
  all.probs = c()

  all.matrices = array(dim=c( dim(my.matrices[[1]]), length(my.matrices) ) ) # ROW * COL * DEEP

  COD.NUMBER = length(cod.table.f.norm[1,])
  
  for ( COD.NUM in 1:COD.NUMBER  ){      ## For each codon
    all.matrices[,,COD.NUM] = my.matrices[[COD.NUM]]
  }

  LOOP = 1:(length(diff[,1]))

  #LOOP = 1:20 #(length(diff[,1]))

  COD.INDEXES = matrix( nrow=length(diff[,1]) , rep( seq(0, 64*(COD.NUMBER-1),by=64), length(diff[,1])), byrow=TRUE)
  
  COL.INDEXES = (diff.cut-1)*( length(my.breaks.mean) - 1)
  ROW.INDEXES = mean.cut
  ALL.INDEXES = COD.INDEXES + COL.INDEXES + ROW.INDEXES

  probs = list()
  
  probs[[1]] = sapply( LOOP, function(x){ return( sum( log(all.matrices[ ALL.INDEXES[x,] ]) , na.rm=TRUE)  ) })
  probs[[2]] = sapply( LOOP, function(x){ return( mean(    all.matrices[ ALL.INDEXES[x,] ]  , na.rm=TRUE)  ) })
  return(probs)
}

get.CODON.all.genes = function(my.MATs = my.mats.COD.TABLE, cod.freqs = COD.freqs, LOG=TRUE, recalculate=FALSE, addName=""){

  
  if( ! file.exists(paste("/data/elevy/70_R_Data/MATRICES/COR.codon",addName,".log",sep="") ) | recalculate){

    L = length(COD.freqs[,1])
    COR.RES = list()
    
    COR.RES[[1]] = matrix(ncol=L, nrow=L, NA) ## LOG
    COR.RES[[2]] = matrix(ncol=L, nrow=L, NA)
       
    for(i in c(1:L) ){

      if(i %% 500 == 0){
        print(paste(i, date()))
      }
      probs = get.CODON.interact.one.gene(gene.num=i, cod.table.f.norm=COD.freqs, my.matrices=my.MATs)
      COR.RES[[1]][i,] = probs[[1]]
      COR.RES[[2]][i,] = probs[[2]]
    }
    colnames(COR.RES[[1]]) = as.vector(REF$ORF)
    colnames(COR.RES[[2]]) = as.vector(REF$ORF)
    rownames(COR.RES[[1]]) = as.vector(REF$ORF)
    rownames(COR.RES[[2]]) = as.vector(REF$ORF)

    NA.pos = which( is.na(cod.table.f.norm[,1])  )
    COR.RES[[1]][ NA.pos,       ]=NA
    COR.RES[[1]][       ,NA.pos ]=NA

    COR.RES[[2]][ NA.pos,       ]=NA
    COR.RES[[2]][       ,NA.pos ]=NA

    print(paste(length(NA.pos)," NA positions have been added"))
    
    write.table(round(COR.RES[[1]],2), file=paste("/data/elevy/70_R_Data/MATRICES/COR.codon",addName,".log",sep="") , quote=FALSE, col.names=TRUE, row.names=TRUE, sep=" ")
    write.table(round(COR.RES[[2]],4), file=paste("/data/elevy/70_R_Data/MATRICES/COR.codon",addName,".lin",sep="") , quote=FALSE, col.names=TRUE, row.names=TRUE, sep=" ")    
  } else {

  
    if(LOG){

      if( file.exists(paste("/data/elevy/70_R_Data/MATRICES/COR.codon",addName,".log",sep="")) ){
        print(paste(date(), "Reading the file ..."))
        COD.cor = read.table(file=paste("/data/elevy/70_R_Data/MATRICES/COR.codon",addName,".log",sep=""), header=TRUE, row.names=c(1), sep=" ")
        print(paste(date(), "Done"))
      } 
      
    } else {
      
      if( file.exists(paste("/data/elevy/70_R_Data/MATRICES/COR.codon",addName,".lin",sep="")) ){
        print(paste(date(), "Reading the file ..."))
        COD.cor = read.table(file=paste("/data/elevy/70_R_Data/MATRICES/COR.codon",addName,".lin",sep=""), header=TRUE, row.names=c(1), sep=" ")
      print(paste(date(), "Done"))
      }
    }
    return(COD.cor)    
  }
}


##
## Returns a matrix CODON USAGE correlations
##
get.codon.cors = function(met=""){

  if( file.exists(paste("/data/elevy/70_R_Data/MATRICES/COR",met,".codon",sep="") ) ){

    print(paste(date(), "Reading the codon.cor file ..."))
    cod.cor = as.matrix(read.table(file=paste("/data/elevy/70_R_Data/MATRICES/COR",met,".codon",sep=""), header=TRUE, row.names=c(1), sep=" "))
    print(paste(date(), "Done"))

  } else {

    ORD = data.frame( ORF = as.vector(get.REF()[,1]), ord=1:6144)

    cod.table = read.table("/data/elevy/51_stickyness_lastAna/results/sc_tAI_table.txt", sep="\t", header=TRUE)
    colnames(cod.table)[1] = "ORF"
    cod.table.tmp = matrix( ncol=64, nrow=length(cod.table[,1]),  as.numeric(as.matrix(cod.table[, c(2:65)]))  )

    cod.table.ord = merge(ORD,cod.table,all.x=TRUE,sort=FALSE)
    cod.table.ord = cod.table.ord[order(cod.table.ord$ord),]
  
    cod.table.f = matrix( ncol=64, nrow=6144,  as.numeric(as.matrix(cod.table.ord[, c(3:66)])) )
    colnames(cod.table.f) = colnames(cod.table.ord)[3:66]
    
    ## Codon frequency table
    cod.table.f.norm = round(100*cod.table.f / rowSums(cod.table.f),2)
    rownames(cod.table.f.norm) = as.vector(ORD$ORF)

    indexes.to.remove = c(  grep("S3",colnames(cod.table.f.norm)), grep("F1",colnames(cod.table.f.norm)),
      grep("Y1",colnames(cod.table.f.norm)),
      grep("Z",colnames(cod.table.f.norm)),
      grep("C",colnames(cod.table.f.norm)),
      grep("W",colnames(cod.table.f.norm)),
      grep("L[4-6]",colnames(cod.table.f.norm)),
      grep("P[2-3]",colnames(cod.table.f.norm)),
      grep("H",colnames(cod.table.f.norm)),
      grep("R[1-4]",colnames(cod.table.f.norm)),
      grep("I[2-3]",colnames(cod.table.f.norm)),
      grep("M1",colnames(cod.table.f.norm)),
      grep("T[2-3]",colnames(cod.table.f.norm)),
      grep("A[2-3]",colnames(cod.table.f.norm)),
      grep("G[1-3]",colnames(cod.table.f.norm))
      )

    cod.table.f.norm = cod.table.f.norm[ , - indexes.to.remove ]

    print(paste(date(), "calculating correlations for codons ..."))
    if(met==""){
      cod.cor = cor(t(cod.table.f.norm), method="pearson", use="pairwise.complete")
    } else {
      cod.cor = cor(t(cod.table.f.norm), method="spearman", use="pairwise.complete")
    }    
    print(paste("Done", date()))
    diag(cod.cor)=0  
    cod.cor = round(cod.cor,3)    
    write.table(cod.cor, file=paste("/data/elevy/70_R_Data/MATRICES/COR",met,".codon",sep=""), quote=FALSE, col.names=TRUE, row.names=TRUE, sep=" ")
  }
  return(cod.cor)
}



## Loads SGA data
##
load.SGA = function(REF=get.REF(), recalculate=FALSE){

  if( file.exists("/data/elevy/70_R_Data/SGA.mat-len1-interm2-string4.aligned") & !recalculate ){

    print(paste(date(), "Reading the file ..."))
    matrix.square = as.matrix(read.table(file="/data/elevy/70_R_Data/SGA.mat-len1-interm2-string4.aligned", header=TRUE, row.names=c(1), sep=" "))
    print(paste(date(), "Done"))
    
  } else {

    ORD = data.frame( ORF = as.vector(REF[,1]), ord=1:6144)
    matrix.square = matrix(ncol=6144,nrow=6144, 0)
    colnames(matrix.square) = REF[,1]
    rownames(matrix.square) = REF[,1]

    types= c("lenient","intermediate","stringent")
    multipl = c(1,2,4)
    
    for ( type in 1:3){

      print(paste("loading type",types[type]," ..."))
      table = read.table(file=paste("/data/elevy/66_YPCAN_ana/data/sgadata_costanzo2009_",types[type],"Cutoff_101120.txt",sep=""),sep="\t", header=FALSE)
      table = table[, c(1,3,5)]
      present.1 = table[,1] %in% REF[,1]
      present.2 = table[,2] %in% REF[,1]

      table = table[ present.1 & present.2,]

      col.i = data.frame( ORF=table[,1], pos.i = 1:(length(table[,1])) , sign= -1*(table[,3]<0) + (1*table[,3]>0) )
      col.j = data.frame( ORF=table[,2], pos.j = 1:(length(table[,1])) , sign= -1*(table[,3]<0) + (1*table[,3]>0) )

      col.i = merge(col.i, ORD);  col.i = col.i[ order(col.i$pos.i),]
      col.j = merge(col.j, ORD);  col.j = col.j[ order(col.j$pos.j),]

      matrix.square[ (col.i$ord-1)*length(ORD[,1]) + col.j$ord ]= multipl[type] * col.j$sign
      matrix.square[ (col.j$ord-1)*length(ORD[,1]) + col.i$ord ]= multipl[type] * col.j$sign
    }
    
    write.table(matrix.square, file="/data/elevy/70_R_Data/SGA.mat-len1-interm2-string4.aligned", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=" ")
  }
     
  return(matrix.square)
}


## Loads SGA CORRELATION data
##
get.SGA.cors = function(met="", recalculate=FALSE){

  if( file.exists(paste("/data/elevy/70_R_Data/MATRICES/COR",met,".SGA",sep="")) & !recalculate ){

    print(paste(date(), "Reading the file ..."))
    SGA.cor = read.table(file=paste("/data/elevy/70_R_Data/MATRICES/COR",met,".SGA",sep=""), header=TRUE, row.names=c(1), sep=" ")
    print(paste(date(), "Done"))

  } else {

    SGA.MAT = load.SGA()
    SGA.MAT[SGA.MAT>0]=0
    
    print(paste(date(), "calculating correlations for SGA ..."))

    if(met==""){
      SGA.cor = cor(SGA.MAT)
    } else {
      SGA.cor = cor(SGA.MAT, method="spearman")
    }
    
    print(paste("Done", date()))
    diag(SGA.cor)=0  
    colnames(SGA.cor)=ORD$ORF #SGA.names
    rownames(SGA.cor)=ORD$ORF #SGA.names
    
    SGA.cor = round(SGA.cor,4)
    
    write.table(SGA.cor, file=paste("/data/elevy/70_R_Data/MATRICES/COR",met,".SGA",sep=""), quote=FALSE, col.names=TRUE, row.names=TRUE, sep=" ")
  }
  return(SGA.cor)
}


## Loads YPCAN network
##
##
## Returns the matrix of raw intensities generated during the large-scale screen.
##
get.YPCAN.matrix = function(REF=get.REF()){

  if( file.exists("/data/elevy/70_R_Data/MATRICES/YPCAN.mat.aligned") ){

    print(paste(date(), "Reading the file ..."))
    matrix.square = read.table(file="/data/elevy/70_R_Data/MATRICES/YPCAN.mat.aligned", header=TRUE, row.names=c(1), sep=" ")
    colnames(matrix.square) = gsub("\\.","-",colnames(matrix.square))
    print(paste(date(), "Done"))
    
  } else {
    matrix.raw = as.matrix(read.table("/data/elevy/66_YPCAN_ana/YPCAN/YPCAN_raw_intensities.txt", sep="\t"))

    row.a     = cbind(read.table("/data/elevy/66_YPCAN_ana/YPCAN/YPCAN_a_array_row_names.txt", stringsAsFactors=FALSE),1:6144)
    row.alpha = cbind(read.table("/data/elevy/66_YPCAN_ana/YPCAN/YPCAN_alpha_array_row_names.txt", stringsAsFactors=FALSE),1:6144)
    
    colnames(row.a) = c("ORF","l")
    matrix.raw[as.vector(row.a$ORF)=="PC",]=0;
    matrix.raw[as.vector(row.a$ORF)=="NC",]=0;
    matrix.raw[as.vector(row.a$ORF)=="-",]=0;

    preys     = read.table("/data/elevy/66_YPCAN_ana/YPCAN/YPCAN_column_names.txt", header=TRUE, stringsAsFactors=FALSE)
    colnames(preys)=c("col","ORF","array");
    
    #preys2 = preys[ which(as.vector(preys$Prey) %in% as.vector(row.a[,1])),]
    preys2 = preys

    matrix.raw2 = matrix.raw[ , which(as.vector(preys$Prey) %in% as.vector(row.a[,1]))]
    preys2[,1] = 1:(length(preys2[,1]))
    colnames(preys2)=c("col","ORF","array");

    ORD = data.frame( ORF = as.vector(REF[,1]), ord=1:6144)
  
    reorder.a = merge(ORD, row.a, all.x=TRUE)
    reorder.a = reorder.a[ order(reorder.a$ord),]

    reorder.alpha = merge(ORD, preys2[,c(1,2)], all.x=TRUE)
    reorder.alpha = reorder.alpha[ order(reorder.alpha$ord),]
  
    matrix.square = matrix(ncol=6144,nrow=6144, 1)

    row.indexes = which(!is.na(reorder.a$l))
    col.indexes = which(!is.na(reorder.alpha$col))

    matrix.square[ row.indexes , col.indexes] = matrix.raw2[ reorder.a$l[!is.na(reorder.a$l)] , reorder.alpha$col[!is.na(reorder.alpha$col)] ]
  
    colnames(matrix.square) = as.vector(REF[,1])
    rownames(matrix.square) = as.vector(REF[,1])

    write.table(matrix.square, file="/data/elevy/70_R_Data/MATRICES/YPCAN.mat.aligned", quote=FALSE, col.names=TRUE, row.names=TRUE, sep=" ")
  }
  return(matrix.square)
}

## Loads all the yeast GO annotations
##
load.sc.GO = function(){
  go.table = read.csv("/data/elevy/70_R_Data/TABLES/sc.GO.table")[,c(2,3,5,7,10,12)]
  colnames(go.table)=c("ORF","Prot.name","name","type","evidence","desc")
  TYPE=rep(0,length(go.table[,1]))
  TYPE[go.table$type=="cellular_component"]="CC"
  TYPE[go.table$type=="biological_process"]="BP"
  TYPE[go.table$type=="molecular_function"]="MF"
  go.table$type=TYPE

  TYPE=rep(0,length(go.table[,1]))
  TYPE[go.table$evidence=="manually curated"]=100
  TYPE[go.table$evidence=="high-throughput"]=10
  TYPE[go.table$evidence=="computational"]=1
  go.table$evidence=TYPE  
  return(go.table)
}

## Format the object given by load.sc.GO, and create CC/BP/MF specific objects
## Min gives the minimum number of entries required per annotation to consider the annotation as valid.
##
format.GO = function(SC.GO, REF=get.REF(), type="CC", min=50){

  SC.GO.sub = SC.GO[SC.GO$type == type,]

  ### Retrieves all the names that have over X components
  GO.names = -sort(-table(SC.GO.sub$name))
  GO.names = names(GO.names[GO.names>min])

  GO.names = GO.names[which(GO.names != "cellular_component")]
  
  SC.GO.sub = SC.GO.sub[SC.GO.sub$name %in% GO.names,]

  SC.GO.sub = SC.GO.sub[SC.GO.sub$ORF %in% REF$ORF,]

  ### Final table:
  ### row=ORF, col=CCname, give=evidence code (0=nothing, 1=pred, 2=HT, 4=exp, 3=pred+HT, 5=pred+exp, 6=HT+exp, 7=all).
  ###
  SC.res = matrix(0, ncol=length(GO.names), nrow=NROW(REF))
  rownames(SC.res) = REF$ORF
  colnames(SC.res) = GO.names

  SC.GO.sub$name = as.vector(SC.GO.sub$name)
  SC.GO.sub$ORF = as.vector(SC.GO.sub$ORF)  
  
  for (i in 1:length(SC.GO.sub[,1])){
    if(i %% 1000 == 0){ 
      print(paste(i, "entries processed"))
    }
    SC.res[ SC.GO.sub$ORF[i], SC.GO.sub$name[i] ] = SC.res[ SC.GO.sub$ORF[i], SC.GO.sub$name[i] ] + SC.GO.sub$evidence[i]
  }
  SC.res.counts = colSums(SC.res>0)
  SC.res = SC.res[ , which(SC.res.counts>min)]
  return(SC.res)
}

##
## Calculates GO enrichment given:
## - a GO matrix containing CC/BP/MF information
## - a list of protein IDs of interest
## - a list of proteins IDs making the universe
## - a number of iterations to calculate empirical Z-scores.
get.GO.enrichment = function(GO.mat, list, univ, N, REF){

    if(length(list)>1){
    
        ## First checks that all proteins of interest are in the universe.
        list.in.univ = list %in% univ
        if(sum(list.in.univ) != length(list.in.univ)){
            print("Warning: some proteins are not in the universe")
            print(list[ !list %in% univ ])
        }
        list = list[list.in.univ]
        list.N = length(unique(list))
        
        GO.mat.univ = 1* (GO.mat[ rownames(GO.mat) %in% univ,]>0)
        univ.size = NROW(GO.mat.univ)
        
        index = rownames(GO.mat.univ) %in% list
        
        GO.obs = apply( GO.mat.univ[index,],2,sum )
        
        GO.rand.obs = matrix(NA, ncol=NCOL(GO.mat.univ), nrow=N)
        for (i in 1:N){
            list.rand = sample(1:univ.size, sum(index), replace=FALSE)
            GO.rand.obs[i,] = apply(GO.mat.univ[list.rand,],2,sum)
        }
        GO.exp  = apply(GO.rand.obs, 2, mean)
        GO.size = apply(GO.mat.univ, 2, sum)
        GO.sd  = apply(GO.rand.obs, 2, sd)
        GO.zscore1 = (GO.obs-GO.exp)/GO.sd
        GO.zscore2 = (GO.obs-GO.exp-1)/GO.sd
        new.order = order( -abs(GO.zscore1))
        GO.enrich = data.frame(
            genes = "",
            size = GO.size[new.order],
            zscore1 = round(GO.zscore1[  new.order ],2) ,
            zscore2 = round(GO.zscore2[  new.order ],2),
            obs=GO.obs[new.order],
            exp=round(GO.exp[new.order],3),
            sd=round(GO.sd[new.order],3), 
            cat = names(GO.zscore1[new.order]))
        GO.enrich.good = GO.enrich$zscore1 > 3 & GO.enrich$zscore2 > 2
        GO.enrich$good = 1*GO.enrich.good
        GO.enrich.final = GO.enrich[c(which(GO.enrich.good)),] #, which(!GO.enrich.good)
        GO.enrich.final = GO.enrich[c(which(GO.enrich.good), which(!GO.enrich.good)),] #, which(!GO.enrich.good)
        GO.enrich.final$genes=""
        for (each.good in 1:sum(GO.enrich.final$good)){
            GO.term = rownames(GO.enrich.final)[each.good]
            GO.col = which(colnames(GO.mat)==GO.term)
            gene.list = REF$Prot.Name[ which(GO.mat[ ,GO.col]>0 & REF$ORF %in% list)]
            if(length(gene.list) < 100){            
                GO.enrich.final$genes[each.good] = paste(gene.list, collapse=", ")
            }
        }
    } else {

        GO.enrich.final =  data.frame(
            genes = c(0),
            size = c(0),
            zscore1 = c(0),
            zscore2 = c(0),
            obs= c(0),
            exp= c(0),
            sd= c(0),
            cat = c(0),
            good=c(0)
            )

    }
    return(GO.enrich.final)
}


## Loads GO-sharing information
##
## 1. localization (a.) using all classes
get.co.GO.cc.l = function(){  
  
  if( length(grep("SC.GO.cc",ls()))>0){
    n.GO = colSums(1* (SC.GO.cc>0))
    GO.cc.l  = SC.GO.cc %*% t(SC.GO.cc)
  } else {

    SC.GO = load.sc.GO()
    SC.GO.cc = format.GO(SC.GO, REF=get.REF(), type="CC", min=15)
    n.GO = colSums(1* (SC.GO.cc>0))
    GO.cc.l  = SC.GO.cc %*% t(SC.GO.cc)    
  }
  return(GO.cc.l)
}

## 2 localization (b.) using only classes which size < 500
get.co.GO.cc.s = function(){

  if( length(grep("SC.GO.cc",ls()))>0){
    n.GO = colSums(1* (SC.GO.cc>0))
    tmp  = SC.GO.cc[,which(n.GO<500)]
    GO.cc.s  = tmp %*% t(tmp)
  } else {
    SC.GO = load.sc.GO()
    SC.GO.cc = format.GO(SC.GO, REF=get.REF(), type="CC", min=15)
    n.GO = colSums(1* (SC.GO.cc>0))
    tmp  = SC.GO.cc[,which(n.GO<500)]
    GO.cc.s  = tmp %*% t(tmp)
  }  
  return(GO.cc.s)
}

## 3 Biological process
get.co.GO.bp = function(){

  if( length(grep("SC.GO.bp",ls()))>0){
    n.GO = colSums(1* (SC.GO.bp>0))
    tmp  = SC.GO.bp[,which(n.GO<500)]
    GO.bp.s  = tmp %*% t(tmp)
  } else {
    SC.GO = load.sc.GO()
    SC.GO.bp = format.GO(SC.GO, REF=get.REF(),type="BP", min=15)
    n.GO = colSums(1* (SC.GO.bp>0))
    tmp  = SC.GO.bp[,which(n.GO<500)]
    GO.bp.s  = tmp %*% t(tmp)
  }  
  return(GO.bp.s)
}

## Get.CYC.data
## ribosomal --> whether ribosomal proteins are loaded or not
## bigMachines --> whether ribosome/DNA polymerase are loaded or not
## bigComplexes --> maximal size of a loaded complex (larger complexes are not considered).
##
get.CYC.data = function(REF = REF, ribosomal=TRUE, bigMachines=TRUE, bigComplexes=1000){

  CYC.data = read.csv("/data/elevy/66_YPCAN_ana/data/CYC2008_complex.csv")
  colnames(CYC.data)[2] = "Prot.Name"
  CYC.data = CYC.data[ order(CYC.data$Complex),]
  CYC.complexes.table = table(as.numeric(factor(CYC.data$Complex)))

  CYC.ord = data.frame( CYC.data$Complex, num = rep( CYC.complexes.table, CYC.complexes.table) )

  CYC.data = cbind(CYC.data, size=CYC.ord[,2])
  
  CYC.ord.nums = order(CYC.ord[,2])
  CYC.data = CYC.data[ order(CYC.ord[,2]),]

  CYC.data$cyc.ord= 1:(length(CYC.data[,1]))

  ORD = data.frame( ORF = as.vector(REF[,1]), ord=1:NROW(REF))
  CYC.tmp = merge(CYC.data, ORD, all.x=TRUE)
  CYC.data = CYC.tmp[ order(-CYC.tmp$cyc.ord),]
  
  ### First sorts the complexes by size
  
  if(bigComplexes == 1000){
  
    if(bigMachines){
      if(ribosomal){
        return(CYC.data)
      } else {
        return(CYC.data[ -grep("ibosom", CYC.data$Complex),])
      }
    } else {
      return(CYC.data[ - c( grep("ibosom", CYC.data$Complex), grep("olymeras", CYC.data$Complex) ),])
    }
  } else {

    comps = table(CYC.data$Complex)
    comps = comps[order(-comps)]
    comps.to.keep = names(comps[comps <= bigComplexes])

    return( CYC.data [ CYC.data$Complex %in% comps.to.keep ,])
  }
}

get.CYC.pairs = function(CYC.data=CYC.data){

  all.ORFs = split(CYC.data$ORF, CYC.data$Complex)

  bench.data=list()

  all.ints = c()
  
  ## First generate all possible pairs among the complexes
  LOOP1 = 1:(length(all.ORFs))
  
  for (complex in LOOP1){

    ORFs.1 = as.vector(all.ORFs[[complex]] )  
    all.ints = rbind( all.ints,
      matrix(
             unlist(sapply( 1:(length(ORFs.1)-1), function(x1) {
               sapply( (x1+1):(length(ORFs.1)) , function(x2){ return(c(ORFs.1[x1],ORFs.1[x2]))}) } )), ncol=2, byrow=TRUE))
  }

  return( all.ints )
}


## Loads PIC data, which is the data from the original study predicting PPIs from codon usage.
## If PIC = "T", then it loads PICT data.
##
get.PIC.mat = function(REF = get.REF(), PIC=""){

  co.PIC = matrix( nrow=6719, ncol=6719, readBin(paste("/data/elevy/66_YPCAN_ana/data/Yeast_PIC",PIC,"/Scerev_PIC",PIC,".bin", sep=""), double(), n = 6719^2)  )

  PIC.names = data.frame ( ORF = as.vector(read.table("/data/elevy/66_YPCAN_ana/data/Yeast_PIC/ordered_genes.txt", header=FALSE)[,1]), pic.ord=1:6719)
   
  ORD = data.frame( ORF = as.vector(REF[,1]), ord=1:6144)
 
  ord.1 = merge( ORD, PIC.names, all.x=TRUE)
  ord.1 = ord.1[order(ord.1$ord),]
  
  co.PIC = co.PIC[ ord.1$pic.ord, ord.1$pic.ord ]

  colnames(co.PIC) = ord.1$ORF
  rownames(co.PIC) = ord.1$ORF

  return(co.PIC)
}

##

##
## Load Marcotte's YeastNet network V3.
##
load.yeastNet.V3 = function(REF, recalculate=FALSE){

    if( file.exists(paste("/data/elevy/70_R_Data/MATRICES/sc.yeastNet.V3.mat",sep="")) & recalculate==FALSE){

        print(paste(date(), "Reading the yeastNet file ..."))
        matrix.square = as.matrix(read.table(file=paste("/data/elevy/70_R_Data/MATRICES/sc.yeastNet.V3.mat",sep=""), header=TRUE, row.names=c(1), sep=" "))
        print(paste(date(), "Done"))

    } else {
    
        net = read.table(file="/data/elevy/70_R_Data/MATRICES/sc.YeastNet.v3.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)
        colnames(net) = c("ORF1","ORF2","all")

        net = net[ net$ORF1 %in% REF[,1] & net$ORF2 %in% REF[,1],]
        
        matrix.square = matrix(ncol=6144,nrow=6144, 0)
        colnames(matrix.square) = REF[,1]
        rownames(matrix.square) = REF[,1]
        tmp = sapply(1:NROW(net), function(x){ matrix.square[net$ORF1[x], net$ORF2[x]] <<- net$all[x] ; matrix.square[net$ORF2[x], net$ORF1[x]] <<- 1+net[x,3]  })

        #net = read.table(file="/data/elevy/70_R_Data/MATRICES/sc.YeastNet.v3.benchmark.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)
        #colnames(net) = c("ORF1","ORF2")

        #net = net[ net$ORF1 %in% REF[,1] & net$ORF2 %in% REF[,1],]

        #tmp = sapply(1:NROW(net), function(x){ matrix.square[net$ORF1[x], net$ORF2[x]] <<- 100 ; matrix.square[net$ORF2[x], net$ORF1[x]] <<- 100 })
        
        write.table(matrix.square, file=paste("/data/elevy/70_R_Data/MATRICES/sc.yeastNet.V3.mat",sep=""), quote=FALSE, col.names=TRUE, row.names=TRUE, sep=" ")
    }
  
    return(matrix.square)
}


