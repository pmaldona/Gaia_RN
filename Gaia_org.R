library("limSolve")
library("XML")


# returns reaction network to a list object from an .sbml file
rn.xmlextract <- function(f, verbose=F) {
  p <- xmlRoot(xmlParse(f))[["model"]]
  sp.id <- NULL
  sp.name <- NULL
  if ("listOfSpecies" %in% names(p)) {
    sp.l <- xmlApply(p[["listOfSpecies"]],xmlAttrs)
    sp.id <- as.vector(sapply(sp.l,function(e) e["id"]))
    sp.name <- as.vector(sapply(sp.l,function(e) e["name"]))
  }
  sp.idn <- sp.id
  i <- which(grepl("^S[0-9]+$",sp.idn))
  sp.idn[i] <- sp.name[i]
  p <- p[["listOfReactions"]]
  n <- which(xmlSApply(p, xmlName) == "reaction")
  a <- xmlApply(p, xmlAttrs)
  l <- list()
  n <- which(xmlApply(p, xmlName) == "reaction")
  R <- list()
  for (i in n) {
    reversible <- a[[i]]["reversible"]
    reversible <- if (!is.na(reversible) && reversible=="false") F else T
    reactants <- if("listOfReactants" %in% names(p[[i]])) xmlApply(p[[i]][["listOfReactants"]], xmlAttrs) else NULL
    products <- if("listOfProducts" %in% names(p[[i]])) xmlApply(p[[i]][["listOfProducts"]], xmlAttrs) else NULL
    r.stoichiometry <- sapply(reactants,function(e) { s <- e["stoichiometry"]; if (is.na(s)) 1 else s } )
    p.stoichiometry <- sapply(products,function(e) { s <- e["stoichiometry"]; if (is.na(s)) 1 else s } )
    names(r.stoichiometry) <- NULL; names(p.stoichiometry) <- NULL
    reactants <- if(is.null(reactants)) NULL else sapply(reactants,function(e) e["species"])
    products <- if(is.null(products)) NULL else sapply(products,function(e) e["species"])
    names(reactants) <- NULL; names(products) <- NULL
    R <- c(R,list(list( reactants=reactants,
                        products=products,
                        r.stoichiometry=as.numeric(r.stoichiometry),
                        p.stoichiometry=as.numeric(p.stoichiometry) )))
    if (reversible)
      R <- c(R,list(list( reactants=products,
                          products=reactants,
                          r.stoichiometry=as.numeric(p.stoichiometry),
                          p.stoichiometry=as.numeric(r.stoichiometry) )))
    if (verbose) {
      print("----------------------------")
      print(reversible)
      print(c(reactants,"-->",products))
      print(c(r.stoichiometry,"-->",p.stoichiometry))
    }
  }
  list(sp.id=sp.id,sp.name=sp.name,sp.idn=sp.idn,reac=R)
}


# returns the reaction network (as a list) including the stoichiometric matrix. 
# Within the execution it computes whether the network is an organization or not. 
rn.proc <- function(f=rn.test,tot_org=T) {
  if (is.numeric(f)) {
    f <- paste0("./ReacNet/",dir("ReacNet","*.xml")[f])
  } 
  cat("f =",f,"\n")
  rn <- rn.xmlextract(f)
  scod <- 1:length(rn$sp.id); names(scod) <- rn$sp.id
  mr <- mp <- matrix(0,length(scod),length(rn$reac))
  rownames(mr) <- rn$sp.id; rownames(mp) <- rn$sp.id

  
  i <- 1
  for (r in rn$reac) {
    l <- length(r$reactants); lp <- length(r$products)
    j <- 1; for (s in r$reactants) { mr[s,i] <- r$r.stoichiometry[j]; j <- j + 1 } 
    j <- 1; for (s in r$products) { mp[s,i] <- r$p.stoichiometry[j]; j <- j + 1 } 
    i <- i + 1
  }
  cat(nrow(mr),"especies,",ncol(mr), "reacciones\n")
  
  if(!tot_org){
    rn <- c(rn,list(spc=scod,mr=mr,mp=mp))
    nsp <- rn.linp_org(rn,rn$sp.id,F)
    rn$spc=scod
    rn$mr=mr
    rn$mp=mp
    rn$nsp=nsp
    rn <- rn
    return(rn)
  }
  
  m=mp-mr
  
  rsp <- which(sapply(1:length(m[,1]),function(i) any(m[i,]!=0))) 
  csp <- which(sapply(1:length(m[,1]),function(i) all(m[i,]<=0))) 
  rcsp <- intersect(rsp,csp) 
  if(length(rcsp!=0)){
    cm <- matrix(0,length(scod),length(rcsp))
    mr <- cbind(mr,cm)
    
    for(i in 1:length(rcsp)){
      cm[rcsp[i],i]=1
      rn$reac <- list.append(rn$reac,list(reactants=NULL,products=rcsp[i],r.stoichiometry=numeric(0),p.stoichiometry=1))
    } 
    mp <- cbind(mp,cm)
    
    
  }
  
  
  rn$spc=scod
  rn$mr=mr
  rn$mp=mp
  nsp <- rn.linp_org(rn,rn$sp.id,F)
  
  if(length(nsp)==0) {
    rn <- rn
    return(rn)
  }
  
  cm <- matrix(0,length(scod),length(nsp))
  rn$mr <- cbind(rn$mr,cm)
  for(i in 1:length(nsp)){ 
    cm[nsp[i],i]=1
    rn$reac <- list.append(rn$reac,list(reactants=NULL,products=nsp[i],r.stoichiometry=numeric(0),p.stoichiometry=1))
  }
  rn$mp <- cbind(rn$mp,cm)
  
  rn <- rn
  return(rn)
}

# returns the required components so the proposed reaction network to be an organization. (function used un rn.proc)
rn.linp_org <- function(crn,id=crn$sp.id,verbose=T) {
  Ns <- nrow(crn$mr); Nr <- ncol(crn$mr)
  E <- cbind(diag(Ns),-diag(Ns),crn$mp-crn$mr) # virtual reactions of introduction and destruction of components are added
  f <- rep(0,Ns) # production equals 0 (since we add virtual destruction reactions)
  G <- diag(ncol(E))
  h <- rep(1.0,ncol(E)); h[1:(2*Ns)] <- 0 # virtual reactions can have rate 0, the others greater than 0 (arbitrarily 1.0)
  Cost <- rep(0,ncol(E)); Cost[1:Ns] <- 1 # ideally cost 0 (no virtual input reaction operating)
  rn.linp.r <- linp(E,f,G,h,Cost)
  X <- rn.linp.r$X
  k <- which(X[(1:Ns)]>0)
  out <- k
  if (verbose) {
    cat("required components: ",id[k],"\n")
  }
  return(out)
}

