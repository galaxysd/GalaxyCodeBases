Detach subdirectory into separate Git repository
======
http://stackoverflow.com/questions/359424/detach-subdirectory-into-separate-git-repository



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



