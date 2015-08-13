setMethod ("initialize", ".ReferenceSpecies", 
           function(.Object, taxId, taxName, proteins, irefindex) {
             
             #Load file
             if (!is.data.frame(irefindex)) {
               message("Load interaction file and search for all relevant 
                   interactions:")
               message("0%")
               irefindex <- read.table(irefindex, sep="\t", comment.char = "",
                                       header=TRUE, quote="",
                                       stringsAsFactors=FALSE)
             } else {
               message("Search for all relevant interactions:")
               message("0%", appendLF = FALSE)
             }

             #Only interspecies interactions
             taxA <- grepl(taxId, irefindex$taxa)
             taxB <- grepl(taxId, irefindex$taxb)
             irefindex <- 
               irefindex[(taxA & taxB) | (grepl("-", irefindex$taxa) & taxB),]
             #Temporary remove duplicates (due different DBs)
             
             if (colnames(irefindex)[1]=="X.uidA") 
               interactions <- 
               irefindex[!duplicated(irefindex[, c("X.uidA", "uidB")]),]
             else
               interactions <- 
               irefindex[!duplicated(irefindex[, c("uidA", "uidB")]),]
             
             rownames(interactions) <- 1:nrow(interactions)
             message("--25%", appendLF = FALSE)             
             
             #Interactor A
             
             if (colnames(irefindex)[1]=="X.uidA") 
               idA <- apply(interactions[, 
                                         c("X.uidA", "altA", 
                                           "OriginalReferenceA", 
                                           "FinalReferenceA", 
                                           "aliasA")], 1, paste, collapse="|")
             else
               idA <- apply(interactions[, 
                                         c("uidA", "altA", 
                                           "OriginalReferenceA", 
                                           "FinalReferenceA", 
                                           "aliasA")], 1, paste, collapse="|")
             
             idA <- strsplit(idA, "|", fixed=TRUE)
             idA <- data.frame(ref=rep(rownames(interactions), 
                                       lapply(idA,length)), temp=unlist(idA),
                               stringsAsFactors=FALSE)  
             idA <- cbind(idA, do.call(rbind, 
                                       strsplit(idA$temp, ":", fixed=TRUE)))
             idA <- idA[, -2]
             colnames(idA) <- c("ref", "database", "accession")
             idA <- idA[!duplicated(idA[, c("ref", "accession")]),]
             message("--50%", appendLF = FALSE)
             #Interactor B
             idB <- apply(interactions[, c("uidB", "altB", 
                                           "OriginalReferenceB", 
                                           "FinalReferenceB","aliasB")],
                          1, paste,collapse="|")
             idB <- strsplit(idB, "|", fixed=TRUE)
             idB <- data.frame(
               ref=rep(rownames(interactions), lapply(idB,length)), 
               temp=unlist(idB), stringsAsFactors=FALSE)  
             idB <- cbind(idB, do.call(rbind, 
                                       strsplit(idB$temp, ":", fixed=TRUE)))
             idB <- idB[,-2]
             colnames(idB) <- c("ref", "database", "accession") 
             idB <- idB[!duplicated(idB[, c("ref", "accession")]),]  
             message("--75%", appendLF = FALSE)
             
             #Find relevant complex ids using B interactors
             
             inB <- idB[idB$accession %in% 
                          c(rownames(proteins), proteins[,1]), "ref"]
             inB <- unique(inB)
             complexes <- as.character(
               idA[(idA$ref %in% inB) &
                     (idA$database == "complex"), "accession"])
             
             #All relevant lines in A
             inA <- idA[idA$accession %in% 
                          c(rownames(proteins),proteins[,1], complexes), "ref"]
             
             #Combine all relavant lines
             cmbn <- unique(c(inA, inB))
             
             #Reduce A and B to the relevant lines...
             idA <- idA[idA$ref %in% cmbn,]
             idB <- idB[idB$ref %in% cmbn,]
             
             idA[, 3] <- as.character(idA[, 3])
             idA[, 2] <- as.character(idA[, 2])
             idB[, 3] <- as.character(idB[, 3]) 
             idB[, 2] <- as.character(idB[, 2])
             
             #Additionally search for the aliases and replace with the 
             #original ids
             for (i in 1:nrow(proteins)) {
               temp <- which(idA[, 3] == proteins[i, 1])
               idA[temp, 3] <- rownames(proteins)[i]
               idA[temp, 2] <- "replaced"
               temp <- which(idB[, 3] == proteins[i, 1])
               idB[temp, 3] <- rownames(proteins)[i]
               idB[temp, 2] <- "replaced"              
             }
             
             #...and prepare for the filtering one id per entry
             idA$inProtList <- idA$accession %in% rownames(proteins)
             idB$inProtList <- idB$accession %in% rownames(proteins)
             
             #First only keep non-duplicates, complexes and those in orgin list
             idA <- idA[!duplicated(idA$ref) | (idA$database == "complex") | 
                          idA$inProtList,]
             idB <- idB[!duplicated(idB$ref) | idB$inProtList,]
             
             #Then remove duplicates from last, this will keep those which are 
             #in origin list but not the major id (assign by irefindex)
             idA <- idA[!duplicated(idA$ref, fromLast=TRUE),]
             idB <- idB[!duplicated(idB$ref, fromLast=TRUE),]
             
             #Combine
             both <- cbind(idA, idB)
             both[, 1] <- interactions[both[, 1], "irigid"]
             both <- both[, -5]
             colnames(both) <- c("ref", "A.db", "A.accession", 
                                 "A.in.prot.list", "B.db", "B.accession", 
                                 "B.in.prot.list")
             
             irefindex <- irefindex[irefindex$irigid %in% both$ref,]
             
             proteins$used <- rownames(proteins) %in% 
               c(as.character(both$A.accession), 
                 as.character(both$B.accession))
             message("--100%")
             
             not.used <- rownames(proteins[!proteins$used,])
             if (length(not.used) > 0) {
              message("Could not find interactions for",
                      " some of the proteins:\n",
                      paste(not.used,collapse=", "))
             }
             
             .Object@taxId <- taxId
             .Object@taxName <- taxName
             .Object@proteins <- proteins
             .Object@interactions <- both
             .Object@irefindex <- irefindex
             
             .Object
           }
)

setMethod ("initialize", ".ReferenceContainer", 
           function(.Object) {
             .Object
           }
)

setMethod ("initialize", ".TargetSpecies", 
           function(.Object, taxName, taxId) {
             if (is.na(taxName))
               stop("Please enter a title for the target species.")
             .Object@taxName <- taxName
             .Object@taxId <- taxId
             
             .Object
           }
)

setMethod ("initialize", "Path2PPI", 
           function(.Object, pathway, targetName, targetId) {
             
             if(missing(pathway)) 
               stop("Please specify a name for the pathway.")
             if(missing(targetName)) 
               stop("Please specify a name for the target species.")
             if(missing(targetId)) 
               stop("Please specify an identifier for the target species.")
             
             .Object@pathway=pathway
             .Object@referenceContainer <- .ReferenceContainer()
             .Object@targetSpecies <- .TargetSpecies(targetName, targetId)
             
             .Object
           }
)