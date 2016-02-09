# STARRseqData: class definition for the STARRseqData class
setClass("STARRseqData",
         representation=list(
             sample="GRanges",
             control="GRanges"),
         prototype=prototype(
             sample=GRanges(),
             control=GRanges()),
         validity=function(object){
             if (length(seqinfo(sample(object))) == 0 ||
                 length(seqinfo(control(object))) == 0){
                 return("Sample and/or control does not contain any seqinfo.")
             }
             if (sum(seqnames(seqinfo(sample(object))) !=
                     seqnames(seqinfo(control(object)))) != 0){
                 return("Sample and control does not have the same seqinfo.")
             }
             return(TRUE)
         }
)
