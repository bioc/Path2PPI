setMethod("addReference",
          signature="Path2PPI",
          definition=function(path2ppi,taxName,taxId,
                              proteins,irefindex,homologs)
          {
            
            path2ppi@referenceContainer <- .addReferenceSpecies(
              path2ppi@referenceContainer,taxName,taxId,proteins,irefindex)
            
            #Load homology file and remove irrelevant homologs
            if (!is.data.frame(homologs)) {
              message("Load homologs file and remove irrelevant homologs.") 
              homologs <- read.table(file=homologs, sep="\t",quote="",
                                     header=FALSE,stringsAsFactors=FALSE)
            } else
              message("Remove irrelevant homologs.")
            
            proteins <- showReferences(path2ppi,taxId,"proteins")
            
            path2ppi@targetSpecies <- .addHomologs(
              path2ppi@targetSpecies,taxId,
              homologs[homologs[,2]%in%rownames(proteins),])
            
            path2ppi
          }
)

setMethod("showReferences",
          signature="Path2PPI",
          definition=function(path2ppi,species,returnValue)
          {
            if (.countSpecies(path2ppi@referenceContainer)>0) {
              if (!is.na(returnValue) && (returnValue!="proteins") && 
                    (returnValue!="interactions") && 
                    (returnValue!="irefindex"))
                stop("Wrong value for 'returnValue' parameter!")
              
              if (is.na(species)) {
                for (i in 1:.countSpecies(path2ppi@referenceContainer)) {
                  .getReferenceSpecies(path2ppi@referenceContainer, 
                                       i, info=TRUE)
                  cat("\n\n")
                }
              } else {
                if (is.na(returnValue))
                  .getReferenceSpecies(path2ppi@referenceContainer,species,
                                       info=TRUE)
                else {
                  species <- .getReferenceSpecies(path2ppi@referenceContainer,
                                                  species,info=FALSE)
                  if (returnValue=="proteins") 
                    return (species@proteins)
                  if (returnValue=="interactions") 
                    return (species@interactions)
                  if (returnValue=="irefindex") 
                    return (species@irefindex)
                }
                
              }
              
            } else {
              message("No reference species yet.")
            }
          }
)

setMethod("plot.Path2PPI",
          signature="Path2PPI",
          definition=function(x,type="ppi",multiple.edges=FALSE,scores=FALSE,
                              species.colors=c(),vertices.opacity=0.8,
                              use.identifiers=FALSE,protein.labels=NA,
                              show.legend=TRUE,vertices.coordinates=NA,
                              return.coordinates=FALSE,tkplot=FALSE,...)
          {
            path2ppi <- x
            
            if (nrow(getPPI(path2ppi))==0) stop("No predicted PPI yet. Please 
                                                call predictPPI-method first!")
            if (type!="hybrid" & type!="ppi") stop("Wrong network type!")
            
            #Some genereal plot setting
            if (!exists("vertex.size")) vertex.size = 10
            if (!exists("vertex.label.dist")) vertex.label.dist = 0
            if (!exists("vertex.label.degree")) vertex.label.degree = -pi/4
            if (!exists("vertex.label.cex")) vertex.label.cex = 0.8
            if (!exists("vertex.label.color")) vertex.label.color=rgb(
              0.1,0.1,0.1)
            if (!exists("edge.label.cex")) edge.label.cex = 0.8
            if (!exists("edge.label.color")) edge.label.color=rgb(0.2,0.2,0.2)
            if (!exists("vertex.frame.color")) vertex.frame.color=NA
            if (!exists("vertex.label.font")) vertex.label.font=2
            if (!exists("vertex.label.family")) vertex.label.family="sans"
            if (!exists("edge.label.family")) edge.label.family="sans"
            
            
            
            edges <- c("interaction"="solid","homology"="dotted")
            
            if (type=="hybrid") { #Hybrid
              titel <- paste(path2ppi@pathway," in ", .getTaxName(
                path2ppi@targetSpecies)," (",.getTaxId(path2ppi@targetSpecies),
                ") - Hybrid",collapse="",sep="")
              
              g <- getHybridNetwork(path2ppi,igraph=TRUE)
              species <- unique(get.vertex.attribute(g,"Species"))
              
              #Check if user defined COLORS
              if (length(species.colors)>0) {
                species <- setdiff(species,names(species.colors))
                if (length(species) > 0)
                  warning("Species colors vector is not complete. 
                          Default colors were added.")
                temp <- species
                species <- rainbow(length(species))
                names(species) <- temp
                species.colors <- c(species.colors,species)
              }
              else {
                species.colors <- rainbow(length(species))
                names(species.colors) <- species
              }
              
              vertex.color <- as.character(
                species.colors[get.vertex.attribute(g,"Species")])
              edges.color <- as.character(species.colors[get.edge.attribute(
                g,"t.species.id")])
              edges.type <- as.character(edges[get.edge.attribute(g,"type")])
              
              
              
              if (scores) {
                edges.labels <- round(get.edge.attribute(g,"score"),digits=2)
                
                if (type=="hybrid")
                  edges.labels[get.edge.attribute(g,"t.species.id")!=
                                 .getTaxId(path2ppi@targetSpecies)&
                                 get.edge.attribute(g,"type")==
                                 "interaction"]<-NA
                
              } else
                edges.labels <- NA
              
              par(cex.main = 1)
              
              if (is.matrix(vertices.coordinates))
                coordinates <- vertices.coordinates
              else
                coordinates <- layout.auto(g,weight=E(g)$score)
              
            } else { #PPI
              
              if (multiple.edges) {
                titel <- paste(path2ppi@pathway," in ", 
                               .getTaxName(path2ppi@targetSpecies)," (",
                               .getTaxId(path2ppi@targetSpecies),
                               ") - PPI detailed",collapse="",sep="")
                
                g <- getPPI(path2ppi,raw=TRUE,igraph=TRUE)
                
                species <- unique(get.edge.attribute(g,"r.species"))
                species <- c(.getTaxId(path2ppi@targetSpecies),c(species))
                
                #Check if user defined COLORS
                if (length(species.colors)>0) {
                  species <- setdiff(species,names(species.colors))
                  if (length(species) > 0)
                    warning("Species colors vector is not complete. Default 
                            colors were added.")
                  temp <- species
                  species <- rainbow(length(species))
                  names(species) <- temp
                  species.colors <- c(species.colors,species)
                }
                else {
                  species.colors <- rainbow(length(species))
                  names(species.colors) <- species
                }
                vertex.color <- species.colors[1]
                edges.color <- as.character(
                  species.colors[get.edge.attribute(g,"r.species")])
                
              }
              else {
                titel <- paste(path2ppi@pathway," in ", .getTaxName(
                  path2ppi@targetSpecies),
                  " (",.getTaxId(path2ppi@targetSpecies)
                  ,") - PPI",collapse="",sep="")
                
                g <- getPPI(path2ppi,raw=FALSE,igraph=TRUE)
                
                if ((length(species.colors)==0) || !((.getTaxId(
                  path2ppi@targetSpecies) %in% names(species.colors)))) {
                  if (length(species.colors)>0) 
                    warning("Color for target species is not defined in color 
                            vector. Default color is used.")
                  vertex.color <- edges.color <- rainbow(1)
                } else
                  vertex.color <- edges.color <- 
                  species.colors[.getTaxId(path2ppi@targetSpecies)]
              }
              
              edges.type <- "solid"
              
              par(cex.main = 1)
              if (is.matrix(vertices.coordinates))
                coordinates <- vertices.coordinates
              else
                coordinates <- layout.auto(g,weight=E(g)$score)
              
              if (scores)
                edges.labels <- round(get.edge.attribute(g,"score"),digits=2)
              else
                edges.labels <- NA
              
            }
            
            vertex.color <- adjustcolor(vertex.color, vertices.opacity) 
            
            #Get labels
            labels <- get.vertex.attribute(g,"name")
            names(labels) <- labels
            if (!use.identifiers) {
              #Replace with available reference species labels
              for (i in 1:length(path2ppi@referenceContainer@species)) {
                proteins <- showReferences(path2ppi,i,"proteins")
                labels[names(labels) %in% rownames(proteins)] <- 
                  proteins[
                    names(labels[names(labels) %in% rownames(proteins)]), 1]
              }
              
              if (!is.null(names(protein.labels))) {
                labels[names(labels) %in% names(protein.labels)] <- 
                  protein.labels[names(labels[names(labels) %in% 
                                                names(protein.labels)])]
              }
            }
            
            if (tkplot)
              tkplot(g,vertex.label=labels,vertex.label.font=vertex.label.font,
                     layout=coordinates,vertex.color=vertex.color,
                     vertex.label.color=vertex.label.color,
                     vertex.label.dist=vertex.label.dist,
                     vertex.label.degree=vertex.label.degree,
                     vertex.size=vertex.size,vertex.label.cex=vertex.label.cex,
                     vertex.label.family=vertex.label.family,
                     edge.label.family=edge.label.family,
                     vertex.frame.color=vertex.frame.color, 
                     edge.color=edges.color,edge.lty=edges.type,
                     edge.width=ceiling(E(g)$score),edge.label=edges.labels,
                     edge.label.cex=edge.label.cex,
                     edge.label.color=edge.label.color,main=titel,asp=0,
                     edge.loop.angle=pi,margin=0,...)
            else {
              plot(g,vertex.label=labels,vertex.label.font=vertex.label.font,
                   layout=coordinates,vertex.color=vertex.color,
                   vertex.label.color=vertex.label.color,
                   vertex.label.dist=vertex.label.dist,
                   vertex.label.degree=vertex.label.degree,
                   vertex.size=vertex.size,vertex.label.cex=vertex.label.cex,
                   vertex.label.family=vertex.label.family,
                   edge.label.family=edge.label.family,
                   vertex.frame.color=vertex.frame.color,
                   edge.color=edges.color,
                   edge.lty=edges.type,edge.width=ceiling(E(g)$score),
                   edge.label=edges.labels,edge.label.cex=edge.label.cex,
                   edge.label.color=edge.label.color,main=titel,asp=0,
                   edge.loop.angle=pi,margin=0,...)
              if ((type=="ppi") && !multiple.edges) show.legend=FALSE
              if (show.legend) legend("bottomright",names(species.colors),
                                      col=species.colors, pch=19, cex=0.8,
                                      bty="o", text.col=rgb(0.2,0.2,0.2))
            }
            
            if (return.coordinates) {
              return (coordinates)
            }
          }
)


setMethod("getHybridNetwork", "Path2PPI",
          function(path2ppi, igraph=FALSE) {
            used.homologues <- data.frame(
              source.id=
                c(path2ppi@raw.ppi$source.id,path2ppi@raw.ppi$target.id), 
              target.id=c(path2ppi@raw.ppi$r.species.s, 
                          path2ppi@raw.ppi$r.species.t),
              t.species.id=rep(path2ppi@raw.ppi$r.species,2),
              score=c(path2ppi@raw.ppi$h.scoreA,path2ppi@raw.ppi$h.scoreB),
              type="homology", stringsAsFactors=FALSE)
            used.homologues <- used.homologues[!duplicated(
              used.homologues[,c(1,2)]),]
            
            #All used interactions
            used.interactions <- data.frame(
              source.id=path2ppi@raw.ppi$r.species.s,
              target.id=path2ppi@raw.ppi$r.species.t,
              t.species.id=path2ppi@raw.ppi$r.species,
              score=1,type="interaction", stringsAsFactors=FALSE)
            used.interactions <- used.interactions[!duplicated(
              used.interactions[,c(1,2)]),]
            
            hybrid <- rbind(used.homologues, used.interactions, 
                            data.frame(source.id=path2ppi@ppi$source.id,
                                target.id=path2ppi@ppi$target.id,
                                t.species.id=.getTaxId(path2ppi@targetSpecies),
                                score=path2ppi@ppi$score,type="interaction", 
                                stringsAsFactors=FALSE))
            
            if (igraph) {
              targetId <- .getTaxId(path2ppi@targetSpecies)
              taxIds <- setdiff(unique(hybrid$t.species.id),targetId)
              
              vertices <- list()
              #First, get all target species nodes
              temp <- hybrid[which(hybrid$t.species.id==targetId),
                             c("source.id","target.id")]
              vertices[[targetId]] <- 
                c(hybrid[which(hybrid$type=="homology"),"source.id"],
                  temp[,"source.id"],temp[,"target.id"])
              vertices[[targetId]] <- 
                vertices[[targetId]][!duplicated(vertices[[targetId]])]
              
              #Then each reference species node
              for (s in 1:length(taxIds)) {
                temp <- hybrid[
                  (hybrid$type=="interaction") & 
                    (hybrid$t.species.id==taxIds[s]),
                  c("source.id","target.id")]
                vertices[[taxIds[s]]] <- 
                  c(hybrid[(hybrid$type=="homology")&
                             (hybrid$t.species.id==taxIds[s]),"target.id"],
                    temp[,"source.id"],temp[,"target.id"])
                vertices[[taxIds[s]]] <- 
                  vertices[[taxIds[s]]][!duplicated(vertices[[taxIds[s]]])]
              }
              
              vertices <- data.frame(ID=do.call(c,vertices),
                                     Species=rep(names(vertices),
                                                 lapply(vertices,length)),
                                     stringsAsFactors=FALSE)
              
              return (graph.data.frame(hybrid, directed=FALSE, 
                                       vertices=vertices))
            } else
              return (hybrid)
          }
)


setMethod("getPPI", "Path2PPI",
          function(path2ppi, raw=FALSE, igraph=FALSE) {
            #, igraph=TRUE
            
            if (igraph) {
              all.vertices <- c(path2ppi@ppi[,1],path2ppi@ppi[,2])
              all.vertices <- unique(all.vertices)
              
              vertices <- data.frame(ID=all.vertices,stringsAsFactors=FALSE)
              if (raw)
                return (graph.data.frame(path2ppi@raw.ppi[,-11], 
                                         directed=FALSE,
                                         vertices=vertices))
              else
                return (graph.data.frame(path2ppi@ppi[,-5], directed=FALSE, 
                                         vertices=vertices))
            } else {
              if (raw)
                return (path2ppi@raw.ppi[,-11])
              else
                return (path2ppi@ppi[,-5])
            }
          }
          
)


setMethod("showInteraction", "Path2PPI",
          function(path2ppi, interaction, mode="default", verbose=TRUE) {
            index <- which((path2ppi@ppi$source.id==interaction[1]&
                            path2ppi@ppi$target.id==interaction[2])|
                           (path2ppi@ppi$source.id==interaction[2]&
                            path2ppi@ppi$target.id==interaction[1]))
            if (length(index) == 0)
              message("Interaction not found!")
            else {
              #default, detailed, references, references.detailed              
              
              temp.ppi <- path2ppi@ppi[index,]
              raw.indices <- unlist(strsplit(as.character(temp.ppi[,5]),",",
                                             fixed=TRUE))
              temp.raw <- path2ppi@raw.ppi[raw.indices,]
              
              counts <- table(temp.raw$r.species)
              
              if (verbose) {
                cat("Score: ",temp.ppi[,"score"],"\n",sep="")
                cat(nrow(temp.raw), " reference interactions: ",sep="")
                
                for (i in 1:nrow(counts))
                  cat(names(counts)[i]," (",counts[i],") ",sep="")
                cat("\n")
              }
              
              
              if ((mode=="default") && !verbose)
                return (temp.ppi[,-5])
              if (mode=="detailed")
                return (temp.raw)
              if (mode=="references") {
                for (i in 1:nrow(counts)) {
                  if (i==1)
                    df <- cbind(species=names(counts)[i],
                          .getReferenceInteractions(
                            path2ppi@referenceContainer, names(counts)[i],
                            temp.raw[which(temp.raw$r.species==
                                             names(counts)[i]),"ref"],
                  FALSE))
                  else
                    df <- rbind(df,cbind(species=names(counts)[i],
                    .getReferenceInteractions(path2ppi@referenceContainer,
                    names(counts)[i],
                    temp.raw[which(temp.raw$r.species==names(counts)[i]),
                             "ref"],
                    FALSE)))
                }
                
                return (df)
              }
              if (mode=="references.detailed") {
                for (i in 1:nrow(counts)) {
                  if (i==1)
                    df <- .getReferenceInteractions(
                      path2ppi@referenceContainer, names(counts)[i], 
                      temp.raw[which(temp.raw$r.species==names(counts)[i]), 
                               "ref"],
                    TRUE)
                  else
                    df <- rbind(df,
                          .getReferenceInteractions(
                            path2ppi@referenceContainer, names(counts)[i], 
                            temp.raw[which(temp.raw$r.species==
                                             names(counts)[i]),"ref"],
                    TRUE))
                }
                return (df)
              }                
              
            }
          }
          
)


setMethod("predictPPI", "Path2PPI",
          function(path2ppi,mode,h.thresh,h.range,i.thresh,consider.complexes,
                   max.complex.size,decline.self.interaction.ref,
                   decline.self.interaction.tar,verbose) {
            
            .combineScores <- function(x,k=length(x)) {
              #Recursive function, e.g.:
              #x <- c(0.6,0.7,0.9)
              #0.9 + (1-0.9)*0.7 + (1-(0.9+(1-0.9)*0.7))*0.6  
              if (k == 1) return (x[k])
              else {
                if (k == length(x)) sort(x)
                return (.combineScores(x,k-1)+(1-.combineScores(x,k-1))*x[k])
              }
            }
            
            if (max.complex.size==0) consider.complexes=FALSE
            
            taxIds <- .getAllTaxIds(path2ppi@referenceContainer)
            id.counter <- 0 
            ppi <- data.frame(id=numeric(),source.id=character(),
                              target.id=character(),score=numeric(),
                              h.scoreA=numeric(), h.scoreB=numeric(), 
                              r.species=character(),r.species.s=character(),
                              r.species.t=character(), pos.edges=numeric(), 
                              used.edges=numeric(), ref=character(), 
                              stringsAsFactors=FALSE)
            
            #For each species
            for (s in 1:length(taxIds)) { 
              taxId <- taxIds[s]
              taxName <- .getTaxName(path2ppi@referenceContainer,taxId)
              
              if (verbose) message("Begin with ",taxName,sep="")
              
              #Get relevant interactions
              interactions <- .getInteractions(path2ppi@referenceContainer,
                                               taxId,mode,consider.complexes,
                                               max.complex.size)
              
              if (decline.self.interaction.ref) 
                interactions <- interaction[
                  !(interactions$A.accession==interactions$B.accession),
                  ]
              
              n <- 0 #Counter for adopted interactions
              for (i in 1:nrow(interactions)) { 
                
                pA_r <- as.character(interactions$A.accession[i])
                pB_r <- as.character(interactions$B.accession[i])
                
                pA_t <- .getHomologs(path2ppi@targetSpecies,taxId,pA_r,"t")
                pB_t <- .getHomologs(path2ppi@targetSpecies,taxId,pB_r,"t")
                
                if (nrow(pA_t) == 1) hpA.one2one <- TRUE 
                else hpA.one2one <- FALSE
                if (nrow(pB_t) == 1) hpB.one2one <- TRUE 
                else hpB.one2one <- FALSE
                
                #At least one homologous for each protein
                if ((nrow(pA_t) > 0) && (nrow(pB_t) > 0)) { 
                  pos.edges <- 0
                  used.edges <- 0
                  for (pA_t_one in 1:nrow(pA_t)) {
                    for (pB_t_one in 1:nrow(pB_t)) {
                      if (!decline.self.interaction.tar || 
                            (pA_t[pA_t_one,1]!=pB_t[pB_t_one,1])) {
                        
                        #Compute both HomPairScores (hp.score) for each 
                        #Homologue pair if only one homologue was found than 
                        #the score is 1 otherwise it is computed by a linear 
                        #function            
                        if ((hpA.one2one) && (pA_t[pA_t_one,2] <= h.thresh)) 
                          hp.scoreA <- 1 
                        else 
                          hp.scoreA <- homologyScore(pA_t[pA_t_one,2],h.range)
                        if ((hpB.one2one) && (pB_t[pB_t_one,2] <= h.thresh)) 
                          hp.scoreB <- 1 
                        else 
                          hp.scoreB <- homologyScore(pB_t[pB_t_one,2],h.range)
                        
                        #Each score has to be above 0
                        if ((hp.scoreA > 0) && (hp.scoreB > 0)) {
                          pos.edges<-pos.edges+1
                          score <- (hp.scoreA+hp.scoreB)/2
                        }
                        else
                          score <- 0
                        
                        if (score >= i.thresh) {
                          n<-n+1
                          used.edges<-used.edges+1
                          id.counter<-id.counter+1
                          ppi <- rbind(ppi, 
                                       data.frame(id=id.counter, 
                                                  source.id=pA_t[pA_t_one,1],
                                                  target.id=pB_t[pB_t_one,1],
                                                  score=score,
                                                  h.scoreA=hp.scoreA, 
                                                  h.scoreB=hp.scoreB, 
                                                  r.species=taxId,
                                                  r.species.s=pA_r,
                                                  r.species.t=pB_r,
                                                  pos.edges=-1,
                                                  used.edges=-1, 
                                                  ref=interactions[i,"ref"], 
                                                  stringsAsFactors=FALSE))
                          
                        } 
                      }  
                    }
                  }
                  
                  if (nrow(ppi) > 0) {
                    #How many edges were possible from this refrence species 
                    #protein pair (regarding the degree of homology fits)
                    ppi[ppi$pos.edges==-1,"pos.edges"] <- pos.edges 
                    #How many edges were left after new compute score which 
                    #have to be >= i.thresh
                    ppi[ppi$used.edges==-1,"used.edges"] <- used.edges 
                  }        
                }
              }
              
              if (verbose) message(paste(i, 
                                " interactions processed. These lead",
                                " to ", n," interactions in target species.",
                                collapse="",sep=""))
              if (verbose) message("-------------------------------")
              
            } 
            
            #Combine, remove redundancies and increment counter
            save.ppi <- ppi  
            combine.equal.interactions.order <- c()
            
            if (verbose) message("Combine results to one single PPI.")
            
            #New container with combined results
            c.ppi <- data.frame(source.id=character(),target.id=character(),
                                score=numeric(),r.species=character(),
                                reference=character(), stringsAsFactors=FALSE) 
            
            repeat { #For each founded interaction in target species
              
              pA <- ppi[1,"source.id"] #Start (id)
              pB <- ppi[1,"target.id"] #Target (id)
              
              #Find each similar edge
              temp.indices <- 
                which(ppi[,"source.id"]==pA & ppi[,"target.id"]==pB, 
                      arr.ind = TRUE) 
              #vice versa (ie start = target and target = start)
              temp.indices.r <- 
                which(ppi[,"source.id"]==pB & ppi[,"target.id"]==pA, 
                      arr.ind = TRUE) 
              
              #Combine all edges
              if (length(temp.indices.r) > 0) 
                temp.indices <- c(temp.indices, temp.indices.r) 

              #In case of self-interaction (For wach self-interaction two 
              #edges a generated, hence, remove the each second:  A-A would 
              #be found in temp.indices and temp.indices.r
              temp.indices <- unique(temp.indices)
              
              temp.interactions <- ppi[temp.indices,] #Temp. save similar edges

              #In which reference species?
              in.which.species <- unique(temp.interactions$r.species) 
              
              references <- c()
              scores <- c()
              for (i in 1:length(in.which.species)) 
                scores <- c(scores,.combineScores(
                  temp.interactions[temp.interactions[,"r.species"]==
                                      in.which.species[i],"score"]))
              
              references <- 
                c(references, 
                  ((length(combine.equal.interactions.order)+1):
                     (length(combine.equal.interactions.order)+
                        nrow(temp.interactions))))
              combine.equal.interactions.order <- 
                c(combine.equal.interactions.order,temp.interactions[,"id"])
              
              avg.score <- mean(scores, na.rm = TRUE)
              inter.score <- avg.score + length(scores)-1
              
              cur.row <- data.frame(source.id=pA, target.id=pB,
                                    score=inter.score,
                                    r.species=paste(sort(in.which.species),
                                                    collapse=","),
                                    reference=paste(references,collapse=","),
                                    stringsAsFactors=FALSE)                
              
              c.ppi <- rbind(c.ppi, cur.row)

              #Remove current interaction (and all similar) from container and 
              #proceed
              ppi <- ppi[-temp.indices,] 
              
              if (nrow(ppi)==0) {
                if (verbose) 
                  message(paste("A total of ",nrow(c.ppi), 
                                " putative interactions were predicted in ",
                                "target species.",collapse="",sep=""))
                break
              }
              
            }
            
            save.ppi <- save.ppi[combine.equal.interactions.order,-1]
            c.ppi <- c.ppi[order(c.ppi$score,decreasing = TRUE),]
            rownames(save.ppi) <- 1:nrow(save.ppi)
            rownames(c.ppi) <- 1:nrow(c.ppi)
            
            path2ppi@raw.ppi <- save.ppi
            path2ppi@ppi <- c.ppi
            
            path2ppi@h.thresh <- h.thresh
            path2ppi@h.range <- h.range
            path2ppi@i.thresh <- i.thresh
            path2ppi@consider.complexes <- consider.complexes
            path2ppi@max.complex.size <- max.complex.size
            
            
            return (path2ppi)
            
          }  
)

setMethod("removeReference",
          signature="Path2PPI",
          definition=function(path2ppi,species)
          {
            path2ppi@referenceContainer <- 
              .removeReferenceSpecies(path2ppi@referenceContainer,species)
            path2ppi@targetSpecies <- 
              .removeHomologs(path2ppi@targetSpecies,species)
            
            path2ppi
          }
)

homologyScore <- function(e.value,h.range) {
  dist <- log10(h.range[1])-log10(h.range[2])
  m <- 1/(dist)
  b <- -(log10(h.range[2])/dist)
  return (sapply(e.value, function(e) if (e >= h.range[2]) 0 
                 else if (e <= h.range[1]) 1 else round(m*log10(e)+b,2)))
}