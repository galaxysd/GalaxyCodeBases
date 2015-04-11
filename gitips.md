Detach subdirectory into separate Git repository
======
http://stackoverflow.com/questions/359424/detach-subdirectory-into-separate-git-repository

### The Easy Way&trade;

It turns out that this is such a common and useful practice that the overlords of git made it really easy, but you have to have a newer version of git (>= 1.7.11 May 2012). See the **appendix** for how to install the latest git. Also, there's a **real-world example** in the **walkthrough** below.

0. Prepare the old repo

        pushd <big-repo>
        git subtree split -P <name-of-folder> -b <name-of-new-branch>
        popd

  **Note:** `<name-of-folder>` must NOT contain leading or trailing characters.  For instance, the folder named `subproject` MUST be passed as `subproject`, NOT `./subproject/`

  **Note for windows users:** when your folder depth is > 1, `<name-of-folder>` must have *nix style folder separator (/). For instance, the folder named `path1\path2\subproject` MUST be passed as `path1/path2/subproject`

0. Create the new repo

        mkdir <new-repo>
        pushd <new-repo>

        git init
        git pull </path/to/big-repo> <name-of-new-branch>
    
0. Link the new repo to Github or wherever

        git remote add origin <git@github.com:my-user/new-repo.git>
        git push origin -u master

0. Cleanup, *if desired*

        popd # get out of <new-repo>
        pushd <big-repo>

        git rm -rf <name-of-folder>

  **Note**: This leaves all the historical references in the repository.See the **Appendix** below if you're actually concerned about having committed a password or you need to decreasing the file size of your `.git` folder.

...

### Walkthrough

These are the **same steps as above**, but following my exact steps for my repository instead of using `<meta-named-things>`.

Here's a project I have for implementing JavaScript browser modules in node:

    tree ~/Code/node-browser-compat
    
    node-browser-compat
    ├── ArrayBuffer
    ├── Audio
    ├── Blob
    ├── FormData
    ├── atob
    ├── btoa
    ├── location
    └── navigator

I want to split out a single folder, `btoa`, into a separate git repository

    pushd ~/Code/node-browser-compat/
    git subtree split -P btoa -b btoa-only
    popd

I now have a new branch, `btoa-only`, that only has commits for `btoa` and I want to create a new repository.

    mkdir ~/Code/btoa/
    pushd ~/Code/btoa/
    git init
    git pull ~/Code/node-browser-compat btoa-only

Next I create a new repo on Github or bitbucket, or whatever and add it is the `origin` (btw, "origin" is just a convention, not part of the command - you could call it "remote-server" or whatever you like)

    git remote add origin git@github.com:node-browser-compat/btoa.git
    git push origin -u master

Happy day!

**Note:** If you created a repo with a `README.md`, `.gitignore` and `LICENSE`, you will need to pull first:

    git pull origin -u master
    git push origin -u master

Lastly, I'll want to remove the folder from the bigger repo

    git rm -rf btoa

...

### Appendix

#### Latest git on OS X

To get the latest version of git:

    brew install git

To get brew for OS X:

<http://brew.sh>

#### Latest git on Ubuntu

    sudo apt-get update
    sudo apt-get install git
    git --version

If that doesn't work (you have a very old version of ubuntu), try

    sudo add-apt-repository ppa:git-core/ppa
    sudo apt-get update
    sudo apt-get install git

If that still doesn't work, try

    sudo chmod +x /usr/share/doc/git/contrib/subtree/git-subtree.sh
    sudo ln -s \
    /usr/share/doc/git/contrib/subtree/git-subtree.sh \
    /usr/lib/git-core/git-subtree

Thanks to rui.araujo from the comments.

#### clearing your history

By default removing files from git doesn't actually remove them from git, it just commits that they aren't there anymore. If you want to actually remove the historical references (i.e. you have a committed a password), you need to do this:

    git filter-branch --tree-filter 'rm -rf <name-of-folder>' HEAD

After that you can check that your file or folder no longer shows up in the git history at all

    git log -S<name-of-folder> # should show nothing

However, you **can't "push" deletes to github** and the like. If you try you'll get an error and you'll have to `git pull` before you can `git push` - and then you're back to having everything in your history.

So if you want to delete history from the "origin" - meaning to delete it from github, bitbucket, etc - you'll need to delete the repo and re-push a pruned copy of the repo. But wait - **there's more**! - If you're really concerned about getting rid of a password or something like that you'll need to prune the backup (see below).

#### making `.git` smaller

The aforementioned delete history command still leaves behind a bunch of backup files - because git is all too kind in helping you to not ruin your repo by accident. It will eventually deleted orphaned files over the days and months, but it leaves them there for a while in case you realize that you accidentally deleted something you didn't want to.

So if you really want to *empty the trash* to **reduce the clone size** of a repo immediately you have to do all of this really weird stuff:

    rm -rf .git/refs/original/ && \
    git reflog expire --all && \
    git gc --aggressive --prune=now

    git reflog expire --all --expire-unreachable=0
    git repack -A -d
    git prune

That said, I'd recommend not performing these steps unless you know that you need to - just in case you did prune the wrong subdirectory, y'know? The backup files shouldn't get cloned when you push the repo, they'll just be in your local copy.

### Credit

  * http://psionides.eu/2010/02/04/sharing-code-between-projects-with-git-subtree/
  * http://stackoverflow.com/questions/1216733/remove-a-directory-permanently-from-git
  * http://blogs.atlassian.com/2013/05/alternatives-to-git-submodule-git-subtree/
  * http://stackoverflow.com/questions/1904860/how-to-remove-unreferenced-blobs-from-my-git-repo



How do you merge two git repositories ?
======
http://stackoverflow.com/questions/1425892/how-do-you-merge-two-git-repositories

A single branch of another repository can be easily placed under a subdirectory retaining its history. For example:

    git subtree add --prefix=rails git://github.com/rails/rails.git master

This will appear as a single commit where all files of Rails master branch are added into "rails" directory.
However the commit's title contains a reference to the old history tree.

    Add 'rails/' from commit <rev>

Where `<rev>` is a SHA-1 commit hash. You can still see the history, blame some changes.

    git log <rev>
    git blame <rev> -- README.md

Note that you can't see the directory prefix from here since this is an actual old branch left intact.
You should treat this like a usual file move commit: you will need an extra jump when reaching it.

    # finishes with all files added at once commit
    git log rails/README.md
    
    # then continue from original tree
    git log <rev> -- README.md

There are more complex solutions like doing this manually or rewriting the history as described in other answers.

The git-subtree command is a part of official git-contrib, some packet managers install it by default (OS X Homebrew).
But you might have to install it by yourself in addition to git.



Changing remote repository for a git submodule
======
http://stackoverflow.com/questions/913701/changing-remote-repository-for-a-git-submodule

What worked for me (on Windows, using git version 1.8.3.msysgit.0):

 - Update .gitmodules with the path to the new repository
 - Remove the corresponding line from the ".git/config" file
 - Delete the corresponding directory in the ".git/modules/external" directory
 - Delete the checked out submodule directory itself (unsure if this is necessary)
 - Run `git submodule init` and `git submodule update`
 - Make sure the checked out submodule is at the correct commit, and commit that, since it's likely that the hash will be different

After doing all that, everything is in the state I would expect. I imagine other users of the repository will have similar pain when they come to update though - it would be wise to explain these steps in your commit message!



How do I remove a Git submodule?
======
http://stackoverflow.com/questions/1260748/how-do-i-remove-a-git-submodule/16162000#16162000

Since [git1.8.3 (April 22d, 2013)][1]:

> There was no Porcelain way to say "I no longer am interested in this submodule", once you express your interest in a submodule with "`submodule init`".  
"**`submodule deinit`**" is the way to do so.

The deletion process also uses `git rm` (since git1.8.5 October 2013).  

The all removal process would then be:

    git submodule deinit asubmodule    
    git rm asubmodule
    # Note: asubmodule (no trailing slash)
    # or, if you want to leave it in your working tree
    git rm --cached asubmodule

But you seem to still need a:

    rm -rf .git/modules/asubmodule

This is mentioned in [Daniel Schroeder][2]'s [answer][3], and summarized by [Eonil][4] in [the comments][5]:

> This leaves `.git/modules/<path-to-submodule>/` unchanged.  
So if you once delete a submodule with this method and re-add them again, it will not be possible because repository already been corrupted.

----

`git rm`: See [commit 95c16418][6]:

> Currently using "`git rm`" on a submodule removes the submodule's work tree from that of the superproject and the gitlink from the index.  
But the submodule's section in `.gitmodules` is left untouched, which is a leftover of the now removed submodule and might irritate users (as opposed to the setting in `.git/config`, this must stay as a reminder that the user showed interest in this submodule so it will be repopulated later when an older commit is checked out).

> Let "`git rm`" help the user by not only removing the submodule from the work tree but by also removing the "`submodule.<submodule name>`" section from the .gitmodules file and stage both.

-----

`git submodule deinit`: It stems from [this patch][7]:

> With "`git submodule init`" the user is able to tell git he cares about one or more submodules and wants to have it populated on the next call to "`git submodule update`".  
But currently there is no easy way he could tell git he does not care about a submodule anymore and wants to get rid of his local work tree (except he knows a lot about submodule internals and removes the "`submodule.$name.url`" setting from `.git/config` together with the work tree himself).

> Help those users by providing a '**`deinit`**' command.  
This **removes the whole `submodule.<name>` section from `.git/config` either for the given
submodule(s)** (or for all those which have been initialized if '`.`' is given).  
Fail if the current work tree contains modifications unless forced.  
Complain when for a submodule given on the command line the url setting can't be found in `.git/config`, but nonetheless don't fail. 

This takes care if the (de)initialization steps (`.git/config` and `.git/modules/xxx`)

Since git1.8.5, the `git rm` takes *also* care of the:

- '`add`' step which records the url of a submodule in the `.gitmodules` file: it is need to removed for you.
- the submodule **[special entry][8]** (as illustrated by [this question][9]): the git rm removes it from the index:  
`git rm --cached path_to_submodule` (no trailing slash)  
That will remove that directory stored in the index with a special mode "160000", marking it as a submodule root directory.

If you forget that last step, and try to add what was a submodule as a regular directory, you would get error message like:

    git add mysubmodule/file.txt 
    Path 'mysubmodule/file.txt' is in submodule 'mysubmodule'


  [1]: https://github.com/git/git/blob/v1.8.3-rc0/Documentation/RelNotes/1.8.3.txt#L135-L137
  [2]: http://stackoverflow.com/users/2753241/daniel-schroeder
  [3]: http://stackoverflow.com/a/26505847/6309
  [4]: http://stackoverflow.com/users/246776/eonil
  [5]: http://stackoverflow.com/questions/1260748/how-do-i-remove-a-git-submodule/16162000?noredirect=1#comment41729982_16162000
  [6]: https://github.com/git/git/commit/95c16418f0375e2fc325f32c3d7578fba9cfd7ef
  [7]: http://git.661346.n2.nabble.com/PATCH-v3-submodule-add-deinit-command-td7576946.html
  [8]: http://stackoverflow.com/questions/1992018/git-submodule-update-needed-only-initially/2227598#2227598
  [9]: http://stackoverflow.com/q/16574625/6309



How to undo the last commit ?
======
http://stackoverflow.com/questions/927358/how-to-undo-the-last-commit

From the docs for [`git-reset`][91]:

> ### Undo a commit and redo
>
>     $ git commit ...              (1)
>     $ git reset --soft HEAD~1     (2)
>     << edit files as necessary >> (3)
>     $ git add ....                (4)
>     $ git commit -c ORIG_HEAD     (5)
>
> 1. This is what you want to undo
>
> 2. This is most often done when you remembered what you just committed is incomplete, or you misspelled your commit message<sup>1</sup>, or both. Leaves working tree as it was before "commit".
>
> 3. Make corrections to working tree files.
>
> 4. Stage changes for commit.
>
> 5. Commit the changes, reusing the old commit message. `reset` copied the old head to `.git/ORIG_HEAD`; `commit` with `-c ORIG_HEAD` will open an editor, which initially contains the log message from the old commit and allows you to edit it. If you do not need to edit the message, you could use the `-C` option instead.

-----

<sup>Editor's note 1</sup>: You don't need to reset to the an earlier commit if "you misspelled your commit message". If you `reset`, git will not link new activity to the previous commit in any way, giving you a blank slate for a new commit message. The easier option is [`git commit --amend`][92], which will open your default commit message editor pre-populated with the last commit message. 

Beware however that if you have added any new changes to the index, using `commit --amend` will add them to your previous commit.


  [91]: http://git-scm.com/docs/git-reset
  [92]: http://stackoverflow.com/q/179123/1146608



