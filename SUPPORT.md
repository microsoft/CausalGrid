# Building the project
Notes on building:
- install (renv)[https://rstudio.github.io/renv/articles/renv.html] package and then use `renv::init()` and choose "1: Restore the project from the lockfile."
- when building, restart the R session before "Install and restart", otherwise it can't copy over the DLL
- When R CMD CHECK, set the test `do_load_all=F`
- Building copies everything over to temp dir and then deletes, so might want to move the large files (`project/sim.RData`) out to save time.


# Support

## How to file issues and get help  

This project uses GitHub Issues to track bugs and feature requests. Please search the existing 
issues before filing new issues to avoid duplicates.  For new issues, file your bug or 
feature request as a new Issue.

For help and questions about using this project, please use the Discussion section.

## Microsoft Support Policy  

Support for this project is limited to the resources listed above.
