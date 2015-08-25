# WZRun2Analysis
Instructions for check out:

git clone git@github.com:vuko-brigljevic/WZRun2Analysis.git

Instructions for commit:

git add nameOfFile

git commit -m "comment"

git push

Instructions to get newest remote repository:

git pull git@github.com:vuko-brigljevic/WZRun2Analysis.git master

UPDATE:
In some cases 'git pull' fails giving the message:
"git: 'Merge branch 'master' of github.com:vuko-brigljevic/WZRun2Analysis' is not a git command."

Roundabout way:
> git fetch origin master
> git merge FETCH_HEAD

