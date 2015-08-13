test_homologyScore <- function() {
  l <- 1e-100
  u <- 1e-20
  h.range <- c(l,u)
  e.values <- c(1e-20,1e-40,1e-60,1e-80,1e-100)  
  
  values <- homologyScore(e.values,h.range)
  
  checkTrue(all(is.numeric(values)),msg="All a number?")
  
  checkTrue(all(values>=0),msg="All values greater or equal than 0?")
  
  checkTrue(all(values<=1),msg="All values less or equal than 1?")
  
  checkEquals(length(values),length(e.values),
              msg="Same number of values as input?")
  
}

test_Path2PPI_pipeline <- function() {
  data(ai)
  
  ppi <- Path2PPI("Autophagy induction", "Podospora anserina", "5145")
  
  checkTrue(is(ppi,"Path2PPI"),
            msg="Class Path2PPI: Value an Path2PPI-object?")
  
  ppi <- addReference(ppi, "Homo sapiens", "9606", human.ai.proteins, 
                      human.ai.irefindex, pa2human.ai.homologs)
  
  ppi <- addReference(ppi, "Saccharomyces cerevisiae (S288c)", "559292", 
                      yeast.ai.proteins, yeast.ai.irefindex, 
                      pa2yeast.ai.homologs)  
  
  
  checkEquals(length(ppi@referenceContainer@species),2,
              "After addReference(): Are all ref. species stored?")
  
  checkException(addReference(ppi, "Homo sapiens", "9606", human.ai.proteins, 
                              human.ai.irefindex, pa2human.ai.homologs),
                 "Throws Exception if one reference species was added twice?")
  
  ppi <- predictPPI(ppi)
  
  checkEquals(nrow(getPPI(ppi)),8,
              "Default prediction lead to 8 interactions?")
  
  
}