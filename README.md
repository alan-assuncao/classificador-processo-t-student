# Primeiro Markdown
***
Arquivo será escrito **em negrito**, e essa palavra em *itálico*

Vamos começar a digitar, fazendo misturas de __*marcações*__ (negrito e itálico)

Listas numeradas

1. teste1
1. teste 2
   1. subteste

lista demarcada:

1. teste 1
2. teste 2
   1. teste 3 (coloque tres espaços a frente)

lista de tarefas:


- [x] fazer o teste

Criando link:

[Acesse meu Github](https://github.com/alan-assuncao)

Num | Nome | Nota
--- | --- | ---
1| Alan | 10.0
2| Kleunice | 10.0

Olha meu Comando `Rcpp`

Olha essa função

```
W_sparsa_phi = function(phi, W_sparse, D_sparse, W_m, rho, m, quant_t){
  	   
	  phit= as.vector(phi)
         phit_D = phit*D_sparse
         phit_W = matrix(0,1,m)
         for (i in 1:W_m) {
             phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]]
             phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]]
         }
     if(quant_t=="TRUE")
     {
        return (-0.5*(t(phit_D)%*%phi - rho*(phit_W%*%phi)))
      }    
       else{
           return (phit_D - rho*t(phit_W)) 
       }
  }
```

Colocando emoji 🐒 

Marcando pessoas @ricardo
















