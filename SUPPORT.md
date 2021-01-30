# Building the project
Requirements
- R (tested on v3.5)

Notes on building:
- You will need RTools (probably at least v3.5)
- install (renv)[https://rstudio.github.io/renv/articles/renv.html] package. Then after opening the project you should be able to use `renv::restore()`. Some packages (such as `brio, cpp11, knitr, ragg, systemfonts, textshaping`) aren't mentioned directly, but are used in building vignettes.
- Given the cpp you should use "Install and restart" (and not use "Load All") to get the new library (though you might be able to get away w/o it if you don't change the DLL). On Windows, when building, you should restart the R session before this otherwise it can't copy over the DLL (it stays in memory).
- If you want updated vignettes to show up when using "Load All", you can use `devtools::build_vignettes()` (possibly with `install=FALSE` to speed things up). They will get placed in `doc/` (not `docs`).
- To build the html help in `docs/` use `pkgdown::build_site()`.
- Building copies everything over to temp dir and then deletes, so might want to move the large files (`project/sim.RData`) out to save time.


# Support

## How to file issues and get help  

This project uses GitHub Issues to track bugs and feature requests. Please search the existing 
issues before filing new issues to avoid duplicates.  For new issues, file your bug or 
feature request as a new Issue.

For help and questions about using this project, please use the Discussion section.

## Microsoft Support Policy  

Support for this project is limited to the resources listed above.
