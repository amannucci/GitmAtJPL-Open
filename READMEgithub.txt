This README is for github usage only. 

I've created the repository using the mac. I am amannucci on github (public, non-JPL identity).
(If there is a public repository I want to be part of, I fork it). 

On the remote machine, in projects area, the git commands will work in the github subdirectory. 
I got things started with 
>>git clone git://github.com/amannucci/GitmAtJPL-Open.git github

Copy all code changes there. First create a branch there called "clusterbranch"?
>>git branch ClusterBranch -- creates the branch
Now I need to set this as the current branch:
>>git checkout ClusterBranch
>>git status -- verifies the branch.

How to make this branch visible to the main repository online?
This fails:
>>git push https://github.com/amannucci/GitmAtJPL-Open.git ClusterBranch
fatal: Unable to find remote helper for 'https'

I also tried to created a "remote" called cluster but it fails with the same error. I think this is just
a short-hand.
>>git remote add cluster https://github.com/amannucci/GitmAtJPL-Open.git
>>git push cluster ClusterBranch
(Fails also)
