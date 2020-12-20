
source("/data/elevy/70_R_Data/bin/RToolBox/MicroscopeToolBox.R")


HSG.design = microscope.get.design(
  F=c("/media/elmicro/or/11HSG/161103-image_HSG01_06"), # PIC FOLDER
  D=c("HSG01_t0"),                                      ## Imaged immediately after HS
  PDEF=c("/media/elmicro/or/11HSG/platedef_PAB1.csv"),  # Plate def
  FORMAT=384,
  OUT=c("_output"),
  CHANELS=c("GFP","RFP"),                               # Channels
  MEASURE.COL = "_int_b7",                             
  DIR.res = "/media/elusers/users/benjamin/Microscope/Working/GitHub/Microscope2016/HSG_Results/"
)

data.1 = microscope.load.data(HSG.design) 

lsos() ## check how much memory the object takes
 
HSG.design = uscope.process.estimate.background(data.1, HSG.design)
data.1    = uscope.process.reorder(data.1, design=HSG.design)
data.1    = uscope.process.remove.first.pic(data.1)
data.1    = uscope.process.remove.background(data.1, HSG.design)

uscope.count.cells(data.1)

data.1    = uscope.process.remove.small(data.1, MIN.size=800,MAX.size=2000)

data.1   = uscope.process.BF(data.1)
data.12   = uscope.process.remove.BF.outliers(data.1, cutoff=0.8)

data.12   = uscope.process.add.ncells(data = data.12)

uscope.count.cells(data.12)

init.result.folders(HSG.design)

##
## First we make a couple of diagnostic plots to make sure things are OK 
##
## BACKGROUND INTENSITY ~ PLATE 
pdf(file= paste(HSG.design$DIR.res,"plate_diag_HSG_int_b9_new.pdf",sep=""),width=20, height=15)
diagnostic.intensity(data=data.12, design=HSG.design, col.of.interest="GFP_int_b9")
dev.off()
## CELL HIGHEST INTENSITY ~ PLATE 
pdf(file= paste(HSG.design$DIR.res,"plate_diag_GFP_int_b1_new.pdf",sep=""),width=20, height=15)
diagnostic.intensity(data=data.12, design=HSG.design, col.of.interest="GFP_int_b1$")
dev.off()
## NUMBER OF CELLS ~ PLATE 
pdf(file= paste(HSG.design$DIR.res,"plate_diag_ncells_new.pdf",sep=""),width=20, height=15)
diagnostic.intensity(data=data.12, design=HSG.design, col.of.interest="GFP_int_b1$", fun2use=length)
dev.off()
