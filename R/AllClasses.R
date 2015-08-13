#PRIVATE CLASSES

# CLASS REFERENCE SPECIES

.ReferenceSpecies <- setClass(".ReferenceSpecies",
    representation(
        taxId = "character",
        taxName = "character",
        proteins = "data.frame",
        irefindex = "data.frame",
        interactions = "data.frame"
    )
)

# CLASS REFERENCE CONTAINER
.ReferenceContainer <- setClass(".ReferenceContainer",
    representation(
        species = "list"
    ),
    prototype=list(
    )
)

# CLASS TARGET SPECIES
.TargetSpecies <- setClass(".TargetSpecies",
    representation(
        taxId = "character",
        taxName = "character",
        proteins = "character",
        homologs = "list"
    )
)

#PUBLIC CLASSES

# CLASS Path2PPI

Path2PPI <- setClass("Path2PPI",
    representation(
        pathway = "character",
        targetSpecies = ".TargetSpecies",
        referenceContainer = ".ReferenceContainer",
        h.thresh = "numeric",
        h.range = "numeric",
        i.thresh = "numeric",
        consider.complexes = "logical",
        max.complex.size = "numeric",
        raw.ppi = "data.frame",
        ppi = "data.frame"
    ),

    prototype=list(
      targetSpecies = NULL,
      referenceContainer = NULL,
      h.thresh = 1e-5,
      h.range = c(1e-100,1e-20),
      i.thresh = 0.7,
      consider.complexes = FALSE,
      max.complex.size = 5
    )
)
