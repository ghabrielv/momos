# Requirements
- R Package in your O.S (Windows - Mac OS - GNU/Linux)
- For only development: libgit2 (Mac OS) - libgit2-dev (GNU/Linux)
# Requirements for development
```
install.packages('devtools')
install.packages(c('deSolve','ggplot2','gridExtra','readxl','gdata','FME'))
install.packages(c('covr','DT'))
install.packages(c("rmarkdown", "knitr"))
```
# Load
```{r}
devtools::load_all()
```
# Run tests
```{r}
devtools::test()
```
# Check coverage of the tests
```{r}
devtools::test_coverage()
```
# Generate documentation
```{r}
devtools::document()
```
# Generate vignettes
```{r}
devtools::build_vignettes()
```
# Check if the package can upload in CRAN
```{r}
devtools::check()
```
# Build
```{r}
devtools::build()
```
# Install package
```{r}
install.packages("./momos_0.2.0.tar.gz", repos = NULL, type = "source")
```
# Send to CRAN
`https://cran.r-project.org/submit.html`
