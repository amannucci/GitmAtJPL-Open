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

I've uploaded new code (tested on Mac) to the hub. Is there an easy 
way to get it on the cluster (w/o wiping out what I have on hub?).

Added a new public key to the hub. Account settings.
See https://help.github.com/articles/generating-ssh-keys

Now this works.
git push cluster ClusterBranch
git remote -v
git remote add cluster git@github.com:amannucci/GitmAtJPL-Open
cluster: git@github.com:amannucci/GitmAtJPL-Open

To update GAIM with code, do this:
git fetch
get merge origin/master

Perhaps the sequence goes like this:
master branch is "stasis"
Need to make a change, create a branch. 
Test within that branch. 
If OK, then merge branch back into master and delete the branch.
Version control is accomplished via "commits", not via branches. 

Right now, I have one branch, master. 

NEXT: how do I commit a change on cluster, and synchronize?

Here is how to commit, but no sync yet (gitpro.pdf online book):
Stage the changes (done before commit)
git add READMEgithub.txt tseriesbatchLoc350.py
git status # For information
git diff --cached # Shows the changes that are staged and ready to be committed
git commit # Does the commit. 
You will see the repository is not yet modified. github.com, login. 

Here is how to sync: push command
git push git@github.com:amannucci/GitmAtJPL-Open master
(There should be a way to define my remote with git@ as I did with cluster above, but leave for now)
(The following would probably work also:
git push cluster master
Use 
git remote -v
to see the listed remotes.
origin remote (https://) seems to work with merges
cluster remote (git@) seems to be needed for push changes from the cluster
)

NOTE: I've defined a key with my public key, this is the github PW for now.

Adding a new file to repository. Let's review steps. New file copied into 
the github directory.
git status 
git add new_file
new_file is already staged.
git commit
git push cluster master

To bring stuff in: (reminder)
git fetch
git merge origin/master
(git merge cluster/master will probably work also)
