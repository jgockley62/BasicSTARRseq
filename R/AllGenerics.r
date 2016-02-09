setGeneric(name="STARRseqData",
           def=function(sample, control, pairedEnd){
               standardGeneric("STARRseqData")
           }
)

setGeneric(name="sample",
           def=function(object){
               standardGeneric("sample")
           }
)

setGeneric(name="sample<-",
           def=function(object, value){
               standardGeneric("sample<-")
           }
)

setGeneric(name="control",
           def=function(object){
               standardGeneric("control")
           }
)

setGeneric(name="control<-",
           def=function(object, value){
               standardGeneric("control<-")
           }
)

setGeneric(name="getPeaks",
           def=function(object, minQuantile=0.9, peakWidth=500,
                        maxPval=0.001, deduplicate=TRUE, model=1){
               standardGeneric("getPeaks")
           }
)
