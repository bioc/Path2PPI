setMethod("show", 
          signature=".ReferenceSpecies",
          definition=function(object)
          {
            titel <- paste(
              object@taxName," (TaxId: ",object@taxId,")\n",collapse="",sep="")
            cat(titel)
            cat(paste(rep("-",nchar(titel)),collapse=""),"\n")
            
            cat(nrow(object@proteins)," proteins (",
                nrow(object@proteins[!object@proteins$used,]),
                " not used)\n",sep="")
            cat(nrow(object@interactions)," interactions:\n",sep="")
            
            both <- nrow(object@interactions[
              object@interactions$A.in.prot.list & 
                object@interactions$B.in.prot.list,])
            one <- nrow(object@interactions[
              object@interactions$A.in.prot.list | 
                object@interactions$B.in.prot.list,])
            complex.interactions <- nrow(object@interactions[
              object@interactions$A.db=="complex",])
            complexes <- length(unique(as.character(object@interactions[
              object@interactions$A.db=="complex","ref"])))
            cat("- ",both,
                " interactions have both interactors in protein list.\n",
                sep="")
            cat("- ",
                one," interactions have at least one ", 
                "interactor in protein list.\n",sep="")
            cat("- ",complex.interactions," interactions in ",
                complexes," protein complexes.\n",sep="")
            #}
            
            
          }
)

setMethod("show",
          signature=".ReferenceContainer",
          definition=function(object)
          {
            
            if (length(object@species) > 0) {
              cat("This reference container contains ",length(object@species),
                  " species:\n\n")
              for (i in 1:length(object@species)) {
                show(object@species[[i]])
                cat("\n")
              }
            }
            else
              cat("No reference-species in container yet.\n")
            
          }
)

setMethod("show",
          signature=".TargetSpecies",
          definition=function(object)
          {
            cat("Target species: ",object@taxName,"\n")
            cat("----------------------------------------------------\n")
            if (length(object@homologs) > 0) {
              for (i in 1:length(object@homologs)) {
                cat(nrow(object@homologs[[i]]), "homologous relationships to 
                    species with taxonomy id ", names(object@homologs)[i],"\n")
              }
            }
            else
              cat("No homologous relationships in container yet.\n")
            if (length(object@proteins) > 0) 
              cat("Relevant proteins: ",length(object@proteins),"\n")
            else
              cat("No relevant proteins yet.\n")
          }
)

setMethod("show",
          signature="Path2PPI",
          definition=function(object)
          {
            titel <- paste(object@pathway," in ", 
                           .getTaxName(object@targetSpecies)," (",
                           .getTaxId(object@targetSpecies),")\n",collapse="",
                           sep="")
            cat(titel)
            cat(paste(rep("-",nchar(titel)),collapse=""),"\n")
            r.species.n <- .countSpecies(object@referenceContainer)
            if (r.species.n > 0) 
              cat(.countSpecies(object@referenceContainer),
                  " reference species: ",
                  paste(.getAllTaxIds(object@referenceContainer),
                        collapse=", "),"\n",sep="")
            else
              cat("No reference species yet.\n")
            cat(paste(rep("-",nchar(titel)),collapse=""),"\n")
            if (nrow(object@ppi)>0) {
              #cat(paste(rep("-",nchar(titel)),collapse=""),"\n")
              hybrid <- getHybridNetwork(object)
              n.predicted <- 
                which(hybrid$t.species.id==.getTaxId(object@targetSpecies))
              hybrid <- hybrid[-n.predicted,]
              n.predicted <- length(n.predicted)
              r.species <- unique(hybrid$t.species.id)
              cat("Number of predicted interactions: ",n.predicted,"\n",sep="")
              cat("Predicted PPI based on ",
                  .countSpecies(object@referenceContainer),
                  " reference species:\n",sep="")
              for (i in 1:length(r.species))
                cat(r.species[i]," (",
                    length(which(hybrid$t.species.id==
                                   r.species[i]&hybrid$type=="interaction")),
                    " interactions and ",
                    length(which(hybrid$t.species.id==
                                   r.species[i]&hybrid$type=="homology")),
                    " homologous relations)\n",sep="")
              cat(paste(rep("-",nchar(titel)),collapse=""),"\n")
              cat("Settings:\n",sep="")
              cat("Homology threshold: ",object@h.thresh,"\n",sep="")
              cat("Homology range: [",object@h.range[1],",",
                  object@h.range[2],"]\n",sep="")
              cat("Interactions threshold: ",object@i.thresh,"\n",sep="")
              cat("Consider complexes: ",object@consider.complexes,"\n",sep="")
              if (object@consider.complexes)
                cat("Max. complex size: ",object@max.complex.size,"\n",sep="")
              
            } else
              cat("No predicted PPI yet.\n")
            
          }
)