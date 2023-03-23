#################
#inputfile
rm(list = ls())
library(xcms)
library(gtools)
path_data<-'E:/NTA workflow/negative'
polarity<--1
ppmshift<-0*10^(-6) 
path_results<-'E:/NTA workflow/Result'
setwd(path_data)
msfiles<-list.files()
msfiles<-mixedsort(msfiles)###sort the file name according to number 
electron<-0.000549
mzwin<-3###2.5ppm for mz cutoff
timewin<-0.5###30 sec for rt cutoff
timewin2<-2####60 sec for rt cutoff, since library was established for a long time
xset.raw<-xcmsSet(msfiles,method='centWave',ppm=3,peakwidth=c(10,30),snthresh=10,nSlaves=1)##peak width, the min and max range of chromatographic peaks in seconds
xtest<-xcmsRaw(msfiles[1])
mzrange<-xtest@mzrange
####delete those peaks with extreme m/z values, if this step is not done, some issues may happen during getEIC step#######
xset<-xset.raw
index<-which(xset.raw@peaks[,1]<mzrange[1]+1)
if (length(index)>0){
  xset@peaks<-xset.raw@peaks[-index,]}
index<-which(xset@peaks[,1]>mzrange[2]-1)
if (length(index)>0){
  xset@peaks<-xset@peaks[-index,]}
peaklist<-xset@peaks#into is the peak area, maxo is the max intensity                                                                      
len<-length(xset@peaks[,1])
setwd(path_results)

##################group peaks
xset1<-group(xset,bw=60,minsamp=1,minfrac=1/length(msfiles),mzwid=0.002)

######################retention time correction
#xset2<-retcor(xset1,family="symmetric")
#xset2<-group(xset2,bw=30,minsamp=1,minfrac=1/length(msfiles),mzwid=0.001)###re-group using
xset2<-xset1 


###################function to fill peaks######################
fillpeak<-function(xset.input,ppm,btw){
  ppm<-ppm/10^6
  idsave<-matrix(rep(0,length(xset.input@groupidx)*length(msfiles)*4),nrow=length(xset.input@groupidx)*length(msfiles),ncol=4)
  k<-1
  for (i in 1:length(xset.input@groupidx)){
    index<-unlist(xset.input@groupidx[i])
    sampleid<-xset.input@peaks[index,11]
    mz.value<-mean(xset.input@peaks[index,1])
    rt.value<-mean(xset.input@peaks[index,4])
    for (j in 1:length(unlist(phenoData(xset.input)))){
      index<-which(sampleid==j)
      if (length(index)<1){
        idsave[k,]<-c(i,j,mz.value,rt.value)##save groupidx,sampleid
        k<-k+1
      }}}
  index<-which(idsave[,1]==0)
  idsave<-idsave[-index,]
  newpeak<-matrix(rep(0,nrow(idsave)*5),nrow=nrow(idsave),ncol=5)
  kk<-1
  msfile<-filepaths(xset.input)
  minmz<-min(xset.input@peaks[,1])
  maxmz<-max(xset.input@peaks[,1])
  minrt<-min(xset.input@peaks[,4])
  maxrt<-max(xset.input@peaks[,4])
  if (length(idsave)>0){
    for (k in 1:length(unlist(phenoData(xset.input)))){
      print(c('filling peakID...',k,'of...',length(unlist(phenoData(xset.input)))))
      index<-which(idsave[,2]==k)
      if (length(index)==0){next}
      xraw<-xcmsRaw(msfiles[k])
      for (n in 1:length(index)){
        mz.value<-idsave[index[n],3]
        mzmin<-max(minmz,mz.value-mz.value*ppm)
        mzmax<-min(maxmz,mz.value+mz.value*ppm)
        rt.value<-idsave[index[n],4]
        rtmin<-max(minrt,rt.value-btw)
        rtmax<-min(maxrt,rt.value+btw)
        peak<-rawEIC(xraw,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))
        intensity<-max(peak$intensity)
        newpeak[kk,]<-c(mz.value,rt.value,intensity,idsave[index[n],1],k)##mz, rt, intensity,groupidx,sampleid
        kk<-kk+1
      }}}
  
  
  peakid<-length(xset.input@peaks[,1])
  for (i in 1:nrow(newpeak)){
    groupid<-newpeak[i,4]
    tempvalue<-unlist(xset.input@groupidx[groupid])
    peakid<-peakid+1
    xset.input@groupidx[groupid]<-list(c(tempvalue,peakid))}
  peak.combine<-matrix(0,ncol=11,nrow=nrow(newpeak))
  peak.combine[,1]<-newpeak[,1]
  peak.combine[,4]<-newpeak[,2]
  peak.combine[,9]<-newpeak[,3]
  peak.combine[,11]<-newpeak[,5]
  xset.input@peaks<-rbind(xset.input@peaks,peak.combine)
  return(xset.input)
}

######################
setwd(path_data)
xset3<-fillpeak(xset2,10,20)
test<-unlist(xset3@groupidx)##peak ID for each group
#peak_target<-which(test>=min(range_target)&&test<=max(range_target2))
len<-length(xset3@groupidx)#group number
len2<-length(msfiles)#data files
final<-array(rep(0,len*(len2+2)),dim=c(len,(len2+2)))##columns are m/z, rt, averaged intensity, and then 8 ratios
for (i in 1:len){
  print(c('PeakID...',i,'of...',len))
  temp<-unlist(xset3@groupidx[i])
  len3<-length(temp)
  for (j in 1:len3){
    index1<-xset3@peaks[temp[j],11]
    final[i,1]<-xset3@peaks[temp[j],1]##mz
    final[i,2]<-xset3@peaks[temp[j],4]/60##rt 
    final[i,index1+2]<-max(xset3@peaks[temp[j],9],final[i,index1+2])##intensity,select the maximum one if there are multiple peaks for one sample
  }}
setwd(path_results)
final[which(final==0)]<-100##replace 0 values
write.table(final, file="48NRs_2_Neg.csv", sep = ',',row.names=FALSE,col.names=cbind("mz","rt",t(msfiles)))

