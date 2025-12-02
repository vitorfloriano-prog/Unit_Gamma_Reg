# UnitGammaReg
This package implements a Unit Gamma regression model with bootstrap bias correction. 
It includes a summary() method that behaves similarly to the summary.lm() and 
summary.glm() functions. The main function returns parameter MLEs and 
bias-corrected estimates. Users can also perform hypothesis tests for a single 
parameter ($beta_1$) or multiple-parameter hypothesis tests, retaining the intercept.

## Function Parameters

### `ugamma.fit()` function

```r
ugamma.fit <- function(
  formula = NULL, 
  data = NULL,
  Y = NULL,
  X = NULL,
  intercept = TRUE,
  q = 1, 
  B = 1000
)
```
- `formula`: Specifies the model using the standard R formula interface, similar to 
  `lm()` and `glm()`.  
  Example: `Response ~ predictor1 + predictor2 + ...`.

- `data`: The dataset containing the variables used in the formula, helping to avoid 
  conflicts between objects in the workspace.

- `Y`: If the user chooses not to use the `formula` argument, they can directly 
  provide a numeric vector for the response variable `Y`.

- `X`: If not using the `formula` argument, the user may provide a numeric vector or 
  matrix for the predictors.  
  The input must be numeric: `X = vector` or `X = matrix`.

- `intercept`: If the dataset already includes an intercept column, 
  the user may disable the automatic intercept insertion by setting 
  `intercept = FALSE`. Otherwise, the intercept is included by default.

- `q`: Number of parameters included in the hypothesis test.  
  For `q >= 2`, the tests become multiple-parameter hypothesis tests.

- `B`: Number of bootstrap replications performed during the bias-correction step.

### `summary.Unit.gamma` function

This function requires no additional arguments besides the fitted model, 
for the reasons described previously.
