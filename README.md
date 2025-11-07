# Unit_Gamma_Reg
This project offers an R script implementing a unit gamma model with bootstrap bias correction. It includes a function for the summary() command to work like in lm()/glm(). The main function returns parameter MLEs and bias-corrected estimates. Users can also perform hypothesis tests for beta1 or multiparametric tests, retaining the intercept.

## Parametros das funções 
### Função Ugamma.fit()
#### ugamma.fit <- function(formula = NULL, data = NULL,= NULL, X = NULL,intercepto = TRUE,q = 1, B = 1000) 

- formula: This is the method used to define the response variable and the predictor variables in a clearer way, similar to how it is implemented in the lm() and glm() commands. The input format is: Response $\sim$ predictor1 + predictor2 + ...;
- data: Here, the desired dataset must be specified, which helps avoid database conflicts;Y: If the user prefers not to use the formula argument, they can directly input a vector for Y, in the format Y = vector;
- X: If the user prefers not to use the formula argument, they can directly input a vector or a matrix with multiple columns for X, in the format X = vector or X = matrix. The input data for X must be numeric.
- intercepto (intercept): If the user has a dataset that already includes the intercept, they can choose to remove it by setting intercepto = FALSE. Otherwise, the intercept is inserted into the dataset by default in our code;
- q: This parameter defines how many variables will be tested in the hypothesis tests. With $q \geq 2$, the tests become multiparametric;
- B: This parameter defines the number of repetitions that will be performed during the bootstrap step.

The summary.Unit.gamma function does not require any parameters, only the fitted model, for the reasons already mentioned in the preceding paragraphs.
