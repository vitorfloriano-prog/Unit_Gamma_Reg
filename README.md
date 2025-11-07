# Unit_Gamma_Reg
This project offers an R script implementing a unit gamma model with bootstrap bias correction. It includes a function for the summary() command to work like in lm()/glm(). The main function returns parameter MLEs and bias-corrected estimates. Users can also perform hypothesis tests for beta1 or multiparametric tests, retaining the intercept.

## Parametros das funções 
### Função Ugamma.fit()
#### ugamma.fit <- function(formula = NULL, data = NULL,= NULL, X = NULL,intercepto = TRUE,q = 1, B = 1000) 

- formula:  é a maneira de definir a variavel resposta e as preditoras de maneira mais clara, similiar a maneira que é implementado nos comando lm() e glm(), a forma de entrada é: Resposta $\sim$ preditora1 + preditora 2 + \ldots;
- data: aqui deve-se especificar o dataset desejado, evita conflito de banco de dados;
- Y: caso o usuario não queira a ulizar a formula, ele pode entrar diretamente com um vetor para Y, da forma Y=vetor.;
- X: caso o usuario não queira a ulizar a formula, ele pode entrar diretamente com um vetor ou com uma matriz de varias colunas para X, da forma X=vetor ou X=matriz.;
- intercepto: Caso o usuario penha um dataset que ja possua o intercepto ele pode obtar por retira-lo, basta apenas definir, intercepto=FALSE, caso contrario, em nosso codigo o intercepto é inserido no conjunto de dados por padrão;
- q: este parametro defini quantas variaveis iremos testar nos teste de hipoteses, com q$\geq$2 os teste passam a ser multiparametricos;
- B: este parametros define o numero de repetições que seram realizados na etapa de bootstrap;


### Função summary.Unit.gamma() 
