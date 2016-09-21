#############################
#
# Cython Functions
#
# Functions written in cython
# are grouped here.
# 
# Written by Daniele Michilli
#
#############################

#Installation: sudo python setup.py build_ext --inplace

cimport cython

#---------------------------------
# Gives a pulse-code to each event
#---------------------------------
@cython.boundscheck(False)
@cython.wraparound(False)
def Get_Group(float[::1] DM not None,
          float[::1] Sigma not None,
          float[::1] Time not None,
          int[::1] Pulse not None,
          unsigned int n_steps,
          float durat,
          float step):
  
  cdef:
    unsigned int i, j, k, j_min, j_max, empty, SNR_max
    unsigned int code = 0
    unsigned int dim = len(DM)
    float step_min, step_max, dDM, DM_min
    float DM_new = -1.
    float float_err = 0.001

  # Assign a code to each event.
  # Events must have been ordered in DM and then in Time
  for i in range(0,dim):
    
    #Remove close events at the same DM
    j = i+1
    if DM[i] == DM[j]:
    
      if abs(Time[i]-Time[j]) < durat:
        
        if j < dim : 
        
          if Sigma[i] < Sigma[j] : Pulse[i] = -1
          else : Pulse[j] = -1
  
    if Pulse[i]==-1: continue
    
    # Set a code to the events that aren't in a pulse
    if Pulse[i]==0: 
      code += 1
      Pulse[i] = code
      
    # Defines for a certain DM a range of events that can be grouped
    if DM[i] != DM_new:
      
      j_min = 0
      j_max = dim
      
      DM_new = DM[i]
        
      step_min = step - float_err
      
      step_max = step * (n_steps + 1) + float_err
      
      
      #find the minimum and maximum event in the range
      for j in range(i+1,dim):
        
        dDM = DM[j] - DM_new

        if dDM > step_max:
          
          j_max = j
          
          break

        if dDM > step_min: 
          
          if j_min == 0: j_min = j
          
    empty = 0
    
    if j_min > 0:

      # Gives a code to the next event in the pulse
      for j in range(j_min,j_max):

        if abs(Time[i]-Time[j]) < durat:   #MAYBE add a condition on SNR (attention: dSNR depends on SNR!) 
          
          if Pulse[j] == -1: continue
          
          if Pulse[j] > 0: 
            
            Pulse[j] = -1
            continue
          
          if empty == 0:
            
            Pulse[j] = Pulse[i]
            SNR_max = j
            empty = 1
            DM_min = DM[j]
            
          else:
            
            if DM[j] > DM_min: break
                        
            if Sigma[j] > Sigma[SNR_max]:
              
              Pulse[SNR_max] = -1
              SNR_max = j
              Pulse[j] = Pulse[i]
              
            else:
              
              Pulse[j] = -1
         
  return
  
  