steps.define = function(x,y,time,r,w) {
  ### Arguments:
    # x and y are coordinates defining the trajectory.
    # time is the time of each x-y sample in the trajectory. This function assumes TIME IS CONSTANT!
    # r defines the threshold distance between consecutive samples that must be exceeded, otherwise the model assumes the person paused.
    # w defines the threshold distance between a line connecting a beginning and end point and a sample point in between the beginning and end point.
  ### Values:
    # Returns a list containing x, y, stepData, and the starting position for each pause, in this order.
    # stepData consists of a n-by-5 matrix where each row contains the length of a step, starting position of the step (in x-y trajectory), ending position of the step, r value used, and w value used.
    
  consDist = sqrt((diff(x))^2 + (diff(y))^2)
  
  # 'startPause' and 'endPause' are the starting and ending row, respectively, in x-y coords where person moved less than r. 'pauseTimes' provides the time of each pause (for later use).
  pauseIndex = 0*(1:length(consDist)-1)
  pauseIndex[which(consDist <= r )] = 1
  startPause = c((which(diff(pauseIndex)==1))+1,length(x))
  endPause = c((which(diff(pauseIndex)==-1))+1,length(x))
  pauseTimes = time[endPause] - time[startPause]
  
  anchor = 1
  stepData = NULL
  
  # Loops through segments of the trajectory. Segments are defined as portion of trajectory between paused positions.
  for (i in 1:length(startPause)) {
    lead = startPause[i]
    
    # Establishes a vector from the anchor to the lead.
    while (anchor < startPause[i]) {
      alVec = cbind((x[lead] - x[anchor]),(y[lead] - y[anchor]))
      
      # Catches any instance where person revisited the exact same location (vector (0,0)), and moves anchor-lead up one position.
      if (alVec[1]==0 && alVec[2]==0) {
        anchor = anchor + 1
        lead = lead + 1
        alVec = cbind((x[lead] - x[anchor]),(y[lead] - y[anchor]))
      }
      
      # Finds the projection lines created between the anchor-lead vector and each test vector. Test vectors are vectors between anchor and each point between the anchor and lead position. 
      # Then finds the distances between the end of the projection lines and each test vector, called 'side_dist'.
      testVecs = cbind((x[(anchor+1):(lead-1)] - x[anchor]),(y[(anchor+1):(lead-1)] - y[anchor]))
      projVecs = (apply(sweep(testVecs,2,alVec,'*'),1,sum) / ((sqrt(alVec[1]^2 + alVec[2]^2))^2)) %*% alVec 
      projVecsEnd = sweep(projVecs, 2, cbind(x[anchor],y[anchor]), '+')
      side_dist = sqrt((projVecsEnd[1:nrow(projVecsEnd),1]-x[(anchor+1):(lead-1)])^2 + (projVecsEnd[1:nrow(projVecsEnd),2]-y[(anchor+1):(lead-1)])^2)
      
      # Checks any distance greater than w threshold. 
      # If so, moves 'lead' back to this position and repeats process.
      # Otherwise, establishes a step length as the distance between anchor and lead.
      if (any(side_dist >= w)) {
        lead = anchor + min( which(side_dist >= w) )
      } else {
        stepDist = sqrt( (x[lead]-x[anchor])^2 + (y[lead]-y[anchor])^2 )
        stepData = rbind(stepData,cbind(stepDist,anchor,lead))
        anchor = lead
        lead = startPause[i]
      }
    
    }
    
    anchor = endPause[i]
  
  }
  
  stepData = cbind(stepData,rep(r,nrow(stepData)),rep(w,nrow(stepData)))
  return(list(x,y,stepData,startPause))
  
}
