################################################################################
################ Constructors of STARRseqData class ############################

setMethod(f="STARRseqData",
          signature(sample="character", control="character",
                    pairedEnd="logical"),
          function(sample, control, pairedEnd=TRUE){
              if (!file.exists(sample)){
                  stop("file not found: ", sample)
              }
              if (!file.exists(control)){
                  stop("file not found: ", control)
              }
              if (pairedEnd){
                  message("read sample data")
                  sampleData <- readGAlignmentPairs(sample)
                  message("read control data")
                  controlData <- readGAlignmentPairs(control)
              } else {
                  message("read sample data")
                  sampleData <- readGAlignments(sample)
                  message("read control data")
                  controlData <- readGAlignments(control)
              }
              message("create GRanges from sample data")
              sampleRanges <- granges(sampleData)
              message("create GRanges from control data")
              controlRanges <- granges(controlData)
              new("STARRseqData", sample=sampleRanges, control=controlRanges)
          }
)

setMethod(f="STARRseqData",
          signature(sample="GRanges", control="GRanges"),
          function(sample, control){
              new("STARRseqData", sample=sample, control=control)
          }
)

################################################################################
################# Show method of STARRseqData class ############################

setMethod(f="show",
          signature="STARRseqData",
          definition=function(object){
              cat("STARRseqData object with", 
                  length(sample(object)), "STARR-seq fragments and", 
                  length(control(object)), "input fragments")
          }
)

################################################################################
################ Getter of STARRseqData class slots ############################

setMethod(f="sample",
          signature="STARRseqData",
          definition=function(object){
              object@sample
          }
)

setMethod(f="control",
          signature="STARRseqData",
          definition=function(object){
              object@control
          }
)

################################################################################
################ Setter of STARRseqData class slots ############################

setReplaceMethod(f="sample",
                 signature(object="STARRseqData", value="GRanges"),
                 function(object, value){
                     object@sample <- value
                     validObject(object)
                     object
                 }
)

setReplaceMethod(f="control",
                 signature(object="STARRseqData", value="GRanges"),
                 function(object, value){
                     object@control <- value
                     validObject(object)
                     object
                 }
)
