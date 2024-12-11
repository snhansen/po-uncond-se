library(tidyverse)
source("utils.R")
options(pillar.sigfig = 6)

mtcars <- mtcars |> 
 mutate(cyl = factor(cyl))

model <- glm(mpg ~ cyl*disp, data = mtcars)
var <- "cyl"
vals  <- c(6, 8)
po_se(model, var, vals)
marginaleffects::avg_predictions(model, variables = setNames(list(vals), var), vcov = "HC0") |>
 as_tibble() |>
 select(all_of(var), estimate, std.error)

model <- glm(mpg ~ cyl*disp, data = mtcars, family = gaussian(link = "log"))
var <- "cyl"
vals  <- c(6, 8)
po_se(model, var, vals, )
marginaleffects::avg_predictions(model, variables = setNames(list(vals), var), vcov = "HC0") |>
 as_tibble() |>
 select(all_of(var), estimate, std.error)

model <- glm(am ~ cyl + disp, data = mtcars, family = binomial())
var <- "cyl"
vals  <- c(6, 8)
po_se(model, var, vals, )
marginaleffects::avg_predictions(model, variables = setNames(list(vals), var), vcov = "HC0") |>
 as_tibble() |>
 select(all_of(var), estimate, std.error)