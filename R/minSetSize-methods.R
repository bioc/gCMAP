setMethod(
          "minSetSize",
          signature( sets = "CMAPCollection" ),
          function(sets, universe=NULL, min.members = 5) {
            
            if (is.null(universe)) {
              universe = featureNames(sets)
            }
            sets[, Matrix::colSums( abs( members( sets ) ) ) >= min.members ]
          }
          )

