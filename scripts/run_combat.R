args <- commandArgs()
file<-args[6]
source(file)
EIF='EIF.dat'
SIF='SIF.dat'
ComBat(EIF,SIF,write=T)
