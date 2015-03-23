steps.merge = function(stepData,thetaThreshold) {
  ### Arguments:
    # stepData is a list and is the output from steps.define function.
    # thetaThreshold sets the value in degrees used to decide on merging consecutive step lengths.
  ### Values:
    # Returns a list containing the merged step lengths, original step lengths, and the beginning of each pause in the x-y trajectory.
    # merged step lengths ('mergedSteps') is an n-by-3 matrix where column 1 is the distance of each step (after mergning), the beginning, and end positions of each merged step.
  
  x = stepData[[1]]
  y = stepData[[2]]
  steps = stepData[[3]]
  startPause = stepData[[4]]
  
  mergedSteps = NULL
  firststep = 1
  nextstep = 2
  Vec = cbind(x[steps[firststep:nextstep,3]] - x[steps[firststep:nextstep,2]], y[steps[firststep:nextstep,3]] - y[steps[firststep:nextstep,2]])
  
  # Finds angle between vectors created by consecutive steps.
  while (nextstep < nrow(steps)) {
    dotprod = Vec[1,]%*%Vec[2,]
    xnorm = norm(Vec[1,],type="2")
    ynorm = norm(Vec[2,],type="2")
    Theta = (acos(pmin(pmax(as.numeric(dotprod/(xnorm*ynorm)),-1.0),1.0))) * (180/pi)
    
    # Checks if angle is less than (or equal to) thetaThreshold.
    # If so, re-establishes a vec from the beginning of the first step to the end of the second step (merges steps). 
    # This new vector is used to compare with the vector defining the next step.
    # Otherwise, computes the distance and adds it to 'mergedSteps'.
    if ((Theta <= thetaThreshold) & (!any(startPause == steps[firststep,3]))) {
      Vec = cbind(x[steps[nextstep,3]] - x[steps[firststep,2]], y[steps[nextstep,3]] - y[steps[firststep,2]])
      nextstep = nextstep + 1
      Vec = rbind(Vec, cbind(x[steps[nextstep,3]] - x[steps[nextstep,2]], y[steps[nextstep,3]] - y[steps[nextstep,2]]))
    }  else {
      distance = sqrt((Vec[1,1])^2 + (Vec[1,2])^2)
      newStep = cbind(distance,steps[firststep,2],steps[nextstep-1,3])
      mergedSteps = rbind(mergedSteps,newStep)
      firststep = nextstep
      nextstep = nextstep + 1
      Vec = cbind(x[steps[firststep:nextstep,3]] - x[steps[firststep:nextstep,2]], y[steps[firststep:nextstep,3]] - y[steps[firststep:nextstep,2]])
    }
  }
  
  return(list(mergedSteps,steps,startPause))
}
