
# This file is a generated template, your changes will not be overwritten

modmedClass <- if (requireNamespace('jmvcore')) R6::R6Class(
    "modmedClass",
    inherit = modmedBase,
    private = list(
        .run = function() {

           # get the data

          data <- self$data

          data <- na.omit(data)

          if (self$results$isFilled())
            return()

          ready <- !is.null(self$options$dep) &&  !is.null(self$options$pred) &&  !is.null(self$options$meds)

          if (ready) {

            suppressWarnings({

            results <- moderatedMediationSem (
                       dat = data,
                       yvar  = self$options$dep,
                       mvars = self$options$meds,
                       xvar = self$options$pred,
                       xmmod = self$options$xmmod,
                       mymod = self$options$mymod,
                       cmvars = self$options$covsm,
                       cyvars = self$options$covsy,
                       nboot = self$options$bootstrap)
            
            

            })   # suppressWarnings
            
            
            private$.populateRsqTable(results)

            

        }

        } , # end run
        
        .populateRsqTable = function(results) {
          
          table <- self$results$rsq <- results$output$Rsq
          
          for (rowKey in table$rowKeys) {
            
            row <- table[rowKey]
            
            table$setRow(rowKey=rowKey, values=list(
              value=row)
            )
          }
          
          
        }


        )
)
