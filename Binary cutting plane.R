#cutting plane
cuttingPlane <- function(f.obj, f.conProb, f.dirProb, f.rhsProb)
{  
numCons <- length(f.conProb[1,])
f.binCon <- rep(0,(numCons)^2) #add binary constraints (0<=x<=1)
for(i in 1:(numCons)^2) 
{
 if(i%%(numCons+1)==1) 
 {
   f.binCon[i] <- 1
 }
}
f.binCon <- matrix (f.binCon, nrow = numCons, byrow = TRUE) #binary constraints into matrix
f.con <- rbind(f.conProb,f.binCon,f.binCon) #combine with constraint matrix
zeroDir <-  rep(">=",numCons) 
oneDir <- rep("<=",numCons)
binDir <- c(oneDir,zeroDir)
f.dir <- c(f.dirProb,binDir)
numProb <- length(f.rhsProb)
zeroBound <- rep(0,numCons)
oneBound <- rep(1,numCons)
f.rhsBin <- c(oneBound,zeroBound) #binary constraint right hand sides
f.rhs <- c(f.rhsProb, f.rhsBin) #combine right hand sides
original <- lp("max", f.obj, f.con, f.dir, f.rhs)$solution # relaxed LP solution
vars <- length(original)
cutP.obj <- rep(0,vars) #initialize cutting plane LP objective vector
augmented <- rep(0,vars) #initialize augmented LP solution vector for after adding cutting plane
doneCutting <- FALSE #check if cutting planes can be added
check <- FALSE # checks if augmented is the same as original to decide if more cutting planes can be added
while(doneCutting == FALSE)
{
 checkSol <- all(augmented==0) #only store original as old augmented after 1st round
  if(checkSol == FALSE)
  {
  original <- augmented
  }
  for(j in 1:vars) # set objective values to find cutting plane
  {
  cutP.obj[j] <- (1-original[j]) 
  if (cutP.obj[j]<0) # avoid small negative values
   {
    cutP.obj[j] <- 0
   }
  }
original <- augmented #store original as augmented for 1st time through
cutP.obj <- cutP.obj*10^7 # avoid non-integer coefficients
  for(i in 1:numProb) # find cutting planes by solving a 0-1 problem for each individual constraint
  {
  f.rhs1 <-  c(f.rhsProb[i]+1)
  cutP.dir <- c(">=")
  f.con1 <- matrix(c(f.con[i,]),nrow=1,byrow=TRUE)
  cutP <- lp ("min", cutP.obj, f.con1, cutP.dir, f.rhs1, all.bin = TRUE)
  cutPsol <- lp ("min", cutP.obj, f.con1, cutP.dir, f.rhs1, all.bin = TRUE, num.bin.solns = 5)$solution #get 5 solutions since lpsolve does not always pick the best 
  cutPsol <- cutPsol[-length(cutPsol)]
  solutions <- matrix(cutPsol, nrow = 5, byrow = TRUE)
  max<-10000000
  for(m in 1:5)
  {
  posBin <- all(solutions[m,]>=0) # make sure no negative variable values (problem in lpsolve)
  allZeros <- all(solutions[m,]==0) # do not select all 0 solution when infeasible
  if (posBin == TRUE && allZeros == FALSE)
  {  
    solObj <- solutions[m,]*cutP.obj #pick best solution
    totObj <- sum(solObj)
    if(totObj <= max)
    {
      cutPsol<-solutions[m,]
      max<-totObj
    }
  }  
  }
  cutPVal <- max #save solution's lp value
  extCov <- cutPsol*c(f.con[i,]) #vectors of coefficients with used variables
  minCov <- sum(extCov) 
  if(cutPVal<10^7)
  {
    cutRhs <- sum(cutPsol)-1
    for(k in 1:numCons)
    {
      if((minCov - extCov[k]) > f.rhs1) #find minimal cover
      {
        cutPsol[k] <- 0
      }
    }
    for(k in 1:numCons)
    {
      if(f.con[i,k]>=max(extCov)) #find extended cover
      {
        cutPsol[k] <- 1
      }
    } 
    addCut <- FALSE
    addCutRhs <- FALSE
    n=1
    while((addCut == FALSE || addCutRhs == FALSE) && n < nrow(f.con)) #check to see if the constraint already exists
      {
      addCut <- identical(cutPsol, f.con[n,]) 
      addCutRhs <- identical(cutRhs, f.rhs[n])
      n <- n+1
      }
    if(addCut == FALSE)
     {
    f.con <- rbind(f.con,cutPsol) #add covers to constraints
    f.rhs <- c(f.rhs,cutRhs)
     }
   }
  }
augmented <- lp("max", f.obj, f.con, f.dir, f.rhs)$solution # new solution
check <- identical(original, augmented)
doneCutting <- check
}
print(augmented)
print(f.con)
print(f.rhs)
}
