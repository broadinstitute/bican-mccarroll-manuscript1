# Creating an R package
How to create an R package for the bican-mccarroll-manuscript1 repository.

## Creating a package using RStudio

1. Get started with github and create a branch using the [basic github workflow](https://docs.google.com/document/d/1A5T8ZVhNxP0GCQmmo1R8qaPPXah9cvSLI96U3TZdYos/edit?tab=t.0) 
1. In RStudio, go to `File` -> `New Project...`
1. Pick `New Directory`
1. Pick `R Package`
1. Give your package a name that starts with `bican.mccarroll.` (e.g., `bican.mccarroll.frob`)
1. Create the project as a subdirectory `<your-git-sandbox>/bican-mccarroll-manuscript1/R`
1. Set up github actions for your package by copying two workflow files (replace `frob` with your package name):
   - `cp .github/workflow-templates/R-package.yml .github/workflows/R-frob.yml`
   - `cp .github/workflow-templates/roxygen.yml .github/workflows/roxygen-frob.yml`
1. Edit these two files and replace `<package_name>` with your package name
1. Do your package development as usual, including adding documentation using roxygen2.
1. When you are ready to commit your changes, use the [basic github workflow](https://docs.google.com/document/d/1A5T8ZVhNxP0GCQmmo1R8qaPPXah9cvSLI96U3TZdYos/edit?tab=t.0) to commit and push your changes.

Note that for now we are requiring review before PRs.  We can change if this is cumbersome.

Note also that the github actions are not active until they are merged into the main branch, so they won't run 
until they have been merged.  Subsequent PRs will trigger the R cmd check action.

## Generating documentation via github actions
You can cause roxygen2 to run via github action by creating a PR, and then commenting on the PR
with `/document.<package_name>`, e.g. `/document.frob`.  This will trigger the roxygen action, which will generate the 
documentation.  It's kinda weird, because it's not obvious that the action is running, but
you can see the action in the Actions tab of the repository, and eventually you'll see a new commit added to your branch.
Note also that if you want to add more commits to the branch after the action has run, you should `git pull` the
changes that were added to your branch by the action, otherwise you will get a merge conflict
when you try to push your changes.

This doesn't work all that well, so it may be easier just to generate documentation and NAMESPACE
locally, and then commit the changes to the branch.  