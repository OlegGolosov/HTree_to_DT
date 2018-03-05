## Scripts and macros to make and analyze trees for HADES experiment
1. Directory «treeMaker»
  * Directory «lists» contains lists of trees in Hades format (HTree)
  
  * Directory «scripts»:
  
    env.sh runs the install of root and HYDRA environment
    
    make.sh runs «make» for the source code makeTree.cc
    
    submitJobs.sh runs jobs which convert trees.
    
    jobScript.sh is used to display what happens when we run submitJobs.sh
    
  * Directory «src»:
  
    makeTree.cc is used to convert trees 
    
    MHConstants.h contains different constants and enumerators
    
    mhqvector.h - class for Qn vectors
    
    mhwalldivider.h - divides FW in rings
    
2. Directory «macros» contains scripts for a drawing distribution of different variables 
3. Directory «results» contains figures that were obtained with the scripts from the directory «macros» 
