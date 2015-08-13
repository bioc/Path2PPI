setMethod(".getProteins", ".ReferenceSpecies",
          function(object) {
            
            return (object@proteins)
          }
)


setMethod(".getTaxId", ".ReferenceSpecies",
          function(object) {
            return(object@taxId)
          }
)


setMethod(".getTaxName", ".ReferenceSpecies",
          function(object) {
            return(object@taxName)
          }
)



setMethod(".getReferenceInteractions", ".ReferenceSpecies",
          function(object, ref, detailed=FALSE) {
            are.all.in <- ref %in% object@interactions$ref
            
            if (!all(are.all.in))
              message("Could not find some of the references: ",
                  paste(ref[!are.all.in], collapse=", "), sep="")
            if (detailed)
              return (object@irefindex[object@irefindex$irigid %in% ref,])
            else
              return (object@interactions[object@interactions$ref %in% ref,])
          }
)

setMethod(".getInteractions", ".ReferenceSpecies",
          function(object, mode="both", consider.complexes=FALSE,
                   max.complex.size=10) {
            warning.number <- 1000
            
            who.is.complex <- object@interactions$A.db == "complex"
            complexes <- object@interactions[who.is.complex,]
            interactions <- object@interactions[!who.is.complex,]
            interactions$complex <- FALSE
            
            if (mode == "both") {
              interactions <- 
                interactions[interactions$A.in.prot.list &
                               interactions$B.in.prot.list,]
            }
            
            if (nrow(interactions)>=warning.number) {
              message("Warning: Large network, may take a long time to 
                      compute!\nYour settings will lead to approximately ",
                      nrow(interactions), " interactions in this reference 
                      species.\nIt is recommended to reduce the size of the 
                      protein list or by using the \"both\"-mode.")
              line <- readline("Proceed anyway (y/n)?")
              
              if (line == "n")
                stop("Canceled by user.")
              else
                warning.number <- 100000000 #Ignore second warning
            }            
            
            if (consider.complexes) {
              #Get each complex size
              all.complexes <- table(as.character(complexes$A.accession))
              #Remove complexes which are to large
              all.complexes <- 
                names(all.complexes[all.complexes <= max.complex.size]) 
              
              if (mode == "both") {
                #Keep only those where all participants are in prot list
                which.complexes.have.all <- sapply(all.complexes, function(x) 
                  all(complexes[complexes$A.accession == x,"B.in.prot.list"]))
                all.complexes <- all.complexes[which.complexes.have.all]
              }
              
              complexes <- complexes[as.character(complexes$A.accession) 
                                     %in% all.complexes,]
              
              if (nrow(complexes)>0) {
                approx.number.of.interactions <- 
                  sum(choose(table(as.character(complexes$A.accession)), 2)) +
                  nrow(interactions)
                
                if (approx.number.of.interactions>=warning.number) {
                  message("Warning: Large network, may take a long time to 
                          compute!\n",
                          paste("Your settings will lead to approximately ",
                                approx.number.of.interactions,
                                " interactions in this reference species.\nIt 
                                is recommended to reduce complex size or use 
                                the \"both\"-mode.", collapse=""))
                  line <- readline("Proceed anyway (y/n)?")
                  
                  if (line == "n")
                    stop("Canceled by user.")
                }
                
                #Now translate to interaction list and add to interactions
                all.complexes <- unique(as.character(complexes$A.accession))
                for (a in 1:length(all.complexes)) {
                  temp.complexes <- 
                    complexes[complexes$A.accession == all.complexes[a],]
                  combinations <- 
                    t(combn(as.character(temp.complexes$B.accession),2))
                  temp.complex.df <- data.frame(ref=temp.complexes$ref[1],
                                                A.db="",
                                                A.accession=combinations[,1],
                                                A.in.prot.list=FALSE,B.db="",
                                                B.accession=combinations[,2],
                                                B.in.prot.list=NA,
                                                stringsAsFactors=FALSE)
                  for (b in 1:nrow(temp.complexes)) {
                    fill <- temp.complexes[b, c("B.db", "B.in.prot.list")]
                    fill[, 1] <- as.character(fill[, 1])
                    temp.complex.df[which(temp.complex.df$A.accession ==
                                            temp.complexes$B.accession[b]),
                                    c("A.db", "A.in.prot.list")] <- fill
                    temp.complex.df[which(temp.complex.df$B.accession ==
                                            temp.complexes$B.accession[b]),
                                    c("B.db", "B.in.prot.list")] <- fill
                  }
                  
                  temp.complex.df$complex <- TRUE
                  interactions <- rbind(interactions, temp.complex.df)
                }
              }
              
            }
            
            return (interactions)
          }           
)

setMethod(".addReferenceSpecies", ".ReferenceContainer",
          function(object, taxName, taxId, proteins, irefindex.file) {
            if (!is.character(proteins)) #Not a character vector
              stop ("Proteins list is not a character vector!")
            
            if (is.null(names(proteins))) { #Not a named vector
              if (any(duplicated(proteins))) 
                stop ("Ambiguous entries in your proteins list detected!", 
                      " Please remove duplicated IDs first.")
              proteins <- data.frame(proteins=proteins, row.names=proteins,
                                     stringsAsFactors=FALSE)
            }
            else { #Is a named vector
              proteins <- data.frame(proteins=proteins,
                                     row.names=names(proteins),
                                     stringsAsFactors=FALSE)
            }
            
            if (taxId %in% names(object@species))
              stop ("TaxID ",taxId," already in container!")
            
            species <- .ReferenceSpecies(
              taxId, taxName, proteins, irefindex.file)
            object@species[[taxId]] <- species
            
            return (object)
          }
)


setMethod(".removeReferenceSpecies", ".ReferenceContainer",
          function(object, species) {
            if ((class(species) == "numeric") && 
                  ((species < 0) || (species > length(object@species))))
              stop ("Invalid index!")
            
            if ((class(species) == "character") && 
                  !(species %in% names(object@species)))
              stop ("Species with taxId ", species, " not in container!")
            
            object@species[[species]] <- NULL
            
            return (object)
          }
)


setMethod(".getInteractions", ".ReferenceContainer",
          function(object, taxId, mode="both", consider.complexes=FALSE,
                   max.complex.size=10) {
            return (.getInteractions(object@species[[taxId]],mode,
                                     consider.complexes,max.complex.size))
          }
)

setMethod(".getReferenceInteractions", ".ReferenceContainer",
          function(object,taxId,ref,detailed=FALSE) {
            return (.getReferenceInteractions(object@species[[taxId]],
                                              ref,detailed))
          }
)


setMethod(".getAllTaxIds", ".ReferenceContainer",
          function(object) {
            return (names(object@species))
          }
)

setMethod(".getTaxName", ".ReferenceContainer",
          function(object,taxId) {
            return(.getTaxName(object@species[[taxId]]))
          }
)


setMethod(".getReferenceSpecies", ".ReferenceContainer",
          function(object,species,info,detailed) {
            if (is.numeric(species) && 
                  ((species > length(object@species)) || (species < 0)))
              stop("Wrong species number!")
            if (is.character(species) && 
                  !(species %in% names(object@species)))
              stop("Wrong TaxId!")
            
            if (info)
              show(object@species[[species]])
            else
              return (object@species[[species]])
          }
)


setMethod(".getSourceDBs", ".ReferenceContainer",
          function(object, irefindex) {
            temp <- table(irefindex[,"sourcedb"])
            dbs <- data.frame(db.id=gsub(".*:(.*)\\(.*", "\\1", 
                                         names(temp)),
                              db.name=gsub(".*\\((.*)\\).*", "\\1", 
                                           names(temp)),count=as.numeric(temp))
            
            return (dbs)
          }
)


setMethod(".countSpecies", ".ReferenceContainer",
          function(object) {
            return (length(object@species))
          }
)


setMethod(".getHomologs", ".TargetSpecies",
          function(object,taxId,proteins=NA,is.source.or.target="t") {
            if (is.na(proteins))
              return (object@homologs[[taxId]])
            else {
              if (is.source.or.target == "t")
                return (object@homologs[[taxId]]
                        [object@homologs[[taxId]][,2] %in% proteins,c(1,11)])
              else
                return (object@homologs[[taxId]]
                        [object@homologs[[taxId]][,1] %in% proteins,c(2,11)])
            }
          }
)


setMethod(".addHomologs", ".TargetSpecies",
          function(object,taxId,homologs) {
            object@homologs[[taxId]] <- homologs
            
            return (object)
          }
)



setMethod(".removeHomologs", ".TargetSpecies",
          function(object,species) {
            if ((class(species) == "numeric") 
                && ((species < 0) || (species > length(object@homologs))))
              stop ("Invalid index!")
            
            if ((class(species) == "character") 
                && !(species %in% names(object@homologs)))
              stop ("Species with taxId ",species," not in container!")
            
            object@homologs[[species]] <- NULL
            
            return (object)
          }
)


setMethod(".getTaxId", ".TargetSpecies",
          function(object) {
            return (object@taxId)
          }  
)

setMethod(".getTaxName", ".TargetSpecies",
          function(object) {
            return (object@taxName)
          }  
)