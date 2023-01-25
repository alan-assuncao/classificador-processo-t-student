### IMPLEMENTAÇÃO STAN PARA O CASO MULTICLASSE

rm(list=ls())# limpar a memória do PC

# PACOTES QUE SERÃO UTILIZADOS NO DECORRER DO SCRIPT

library('rstan')
library('MASS')
library('pracma')
library('reshape2')
library('ggplot2')
#library('R.matlab')	# para ler arquivos matlab
#library('fields')	# para calcular distancias

####################################################################################
# MODELO DE CLASSIFICAÇÃO MULTICLASSE COM PROCESSOS GAUSSIANO

# MODELO DE CLASSIFICAÇÃO MULTICLASSE COM BASE NO FORUM STAN SOBRE CLASSIFICAÇÃO MULTICLASSE E NO ARTIGO DE WILLIAMS E BARBER (1998)

model_multiclasse.GP="

// função para implementar a f.c. exponencial quadrática (ARD)
functions{
	matrix L_cov_exp_quad_ARD(vector[] x, real alpha, vector rho, real delta){
		
		int N=size(x);
		matrix[N, N] K;
		real neg_half=-0.5;
		real sq_alpha=square(alpha);

		for(i in 1:(N-1)){
			K[i,i]=sq_alpha + delta;
			for(j in (i+1):N){
				K[i,j]=sq_alpha*exp(neg_half*dot_self((x[i]-x[j])./rho));
				K[j,i]=K[i,j];				
			}
		}
	
	K[N,N]=sq_alpha+delta;
	return cholesky_decompose(K);
	//return K;
	}

}

data{
	// se o banco de dados tiver entradas multivariadas, então utilizar a configuração abaixo
	
	int<lower=1> N; // quantidade de pontos de treinamento
	int<lower=1> D; // dimensão do vetor de entrada
	int<lower=1> C; // quantidade de classes
	vector[D] x[N]; // dados de entrada para treino
	int<lower=1> CN; // dimensão do vetor extendito f dos C processos latentes
	int<lower=0,upper=1> y[N,C]; // matrix de rótulos de classe (codificacao 1 de C)
}

transformed data{
	real delta=1e-9; // aqui seria o jitter (não colocar mais aspas em comentarios)do modelo de classificação similar a sigma2 na regressão
	real w=1;
	real a=1;

}

parameters{
	// se eu utilizar entradas multivariadas, utilizar essa configuração abaixo
	
	vector<lower=0>[D] rho[C]; // comprimento escala
	real<lower=0> alpha[C];	// sinal de variancia
	vector[CN] eta;	
}

model{
		
	vector[CN] f;
	matrix[N,C] fextent;

	{
		matrix[CN,CN] L_K;  // matrix em bloco diagonal formada pelas matrizes k1,...,kc
	
		for(i in 1:CN)
		for(j in 1:CN)
		 L_K[i,j]<-0;
		
		for(i in 1:C)
			L_K[((i-1)*N+1):(i*N),((i-1)*N+1):(i*N)]<-L_cov_exp_quad_ARD(x,alpha[i],rho[i],delta);
		
		f=L_K*eta;			

		for(i in 1:N)
			for(j in 1:C)
	 		 fextent[i,j]<-f[(((j-1)*N+1)+(i-1))];
	}

	// distribuições a priori

	for(i in 1:C)
		rho[i] ~ gamma((a/2),(a/(2*w)));			//rho~inv_gamma(5,5);
	alpha ~ normal(0,1);				//alpha ~ gamma((a/2),(a/(2*w)));
	//a ~ normal(0,1);
	eta ~ normal(0,1);

	for(i in 1:N)
	 y[i]~multinomial(softmax(to_vector(fextent[i]))); // tinha um a seguindo distribuição normal dentro da função softmax

}
"

#####################################################################################
# MODELO DE CLASSIFICAÇÃO DE TP MULTICLASSE COM BASE NO FORUM STAN SOBRE CLASSIFICAÇÃO MULTICLASSE E NO ARTIGO DE WILLIAMS E BARBER (1998)

model_multiclasse.TP="

// função para implementar a f.c. exponencial quadrática (ARD)
functions{
	matrix L_cov_exp_quad_ARD(vector[] x, real alpha, vector rho, real delta){
		
		int N=size(x);
		matrix[N, N] K;
		real neg_half=-0.5;
		real sq_alpha=square(alpha);

		for(i in 1:(N-1)){
			K[i,i]=sq_alpha + delta;
			for(j in (i+1):N){
				K[i,j]=sq_alpha*exp(neg_half*dot_self((x[i]-x[j])./rho));
				K[j,i]=K[i,j];				
			}
		}
	
	K[N,N]=sq_alpha+delta;
	return cholesky_decompose(K);
	}

}

data{
	// se o banco de dados tiver entradas multivariadas, então utilizar a configuração abaixo
	
	int<lower=1> N; // quantidade de pontos de treinamento
	int<lower=1> D; // dimensão do vetor de entrada
	int<lower=1> C; // quantidade de classes
	vector[D] x[N]; // dados de entrada para treino
	int<lower=1> CN; // dimensão do vetor extendito f dos C processos latentes
	int<lower=0,upper=1> y[N,C]; // matrix de rótulos de classe (codificacao 1 de C)
	real<lower=0> v;  // graus de liberdade da distribuição t-student
}

transformed data{
	real delta=1e-9; // aqui seria o jitter (não colocar mais aspas em comentarios)do modelo de classificação similar a sigma2 na regressão
	real w=1;
	real a=1;

	vector[CN] mu;
	for(i in 1:CN) mu[i]<-0; // vetor de médias nulo

}

parameters{
	
	vector<lower=0>[D] rho[C]; // comprimento escala
	real<lower=0> alpha[C];	// sinal de variancia
	vector[CN] eta;
	real<lower=0> u;	
}

model{
		
	vector[CN] f;
	matrix[N,C] fextent;
	{
	 matrix[CN,CN] L_K;  // matrix em bloco diagonal formada pelas matrizes k1,...,kc
		
	 for(i in 1:CN)
	 for(j in 1:CN)
	  L_K[i,j]<-0;

	for(i in 1:C)
		L_K[((i-1)*N+1):(i*N),((i-1)*N+1):(i*N)]<-L_cov_exp_quad_ARD(x,alpha[i],rho[i],delta);
	
	f=mu+sqrt(v/u)*(L_K*eta);	// amostrando de uma t-student mutivariada

	for(i in 1:N)
	   for(j in 1:C)
	 	 fextent[i,j]<-f[(((j-1)*N+1)+(i-1))];
	}

	// distribuições a priori
	
	for(i in 1:C)
		rho[i] ~ gamma((a/2),(a/(2*w)));			//rho~inv_gamma(5,5);
	alpha ~ normal(0,1);				//alpha ~ gamma((a/2),(a/(2*w)));
	//v ~ gamma(2,0.1);
	eta ~ normal(0,1);
	u ~ chi_square(v);

	for(i in 1:N)
	 y[i]~multinomial(softmax(to_vector(fextent[i]))); // tinha um a seguindo distribuição normal dentro da função softmax
}
"

#####################################################################################
# MODELO DE CLASSIFICAÇÃO DE TP MULTICLASSE COM BASE NO FORUM STAN SOBRE CLASSIFICAÇÃO MULTICLASSE E NO ARTIGO DE WILLIAMS E BARBER (1998)

model_multiclasse.TP.modified="

// função para implementar a f.c. exponencial quadrática (ARD)
functions{
	matrix L_cov_exp_quad_ARD(vector[] x, real alpha, vector rho, real delta){
		
		int N=size(x);
		matrix[N, N] K;
		real neg_half=-0.5;
		real sq_alpha=square(alpha);

		for(i in 1:(N-1)){
			K[i,i]=sq_alpha + delta;
			for(j in (i+1):N){
				K[i,j]=sq_alpha*exp(neg_half*dot_self((x[i]-x[j])./rho));
				K[j,i]=K[i,j];				
			}
		}
	
	K[N,N]=sq_alpha+delta;
	return cholesky_decompose(K);
	}

}

data{
	// se o banco de dados tiver entradas multivariadas, então utilizar a configuração abaixo
	
	int<lower=1> N; // quantidade de pontos de treinamento
	int<lower=1> D; // dimensão do vetor de entrada
	int<lower=1> C; // quantidade de classes
	vector[D] x[N]; // dados de entrada para treino
	int<lower=1> CN; // dimensão do vetor extendito f dos C processos latentes
	int<lower=0,upper=1> y[N,C]; // matrix de rótulos de classe (codificacao 1 de C)
}

transformed data{
	real delta=1e-9; // aqui seria o jitter (não colocar mais aspas em comentarios)do modelo de classificação similar a sigma2 na regressão
	real w=1;
	real a=1;

	vector[CN] mu;
	for(i in 1:CN) mu[i]<-0; // vetor de médias nulo
}

parameters{
	
	vector<lower=0>[D] rho[C]; // comprimento escala
	real<lower=0> alpha[C];	// sinal de variancia
	vector[CN] eta;
	real<lower=0> u;	
	real<lower=0> v;  // graus de liberdade da distribuição t-student
}

model{
		
	vector[CN] f;
	matrix[N,C] fextent;
	{
	 matrix[CN,CN] L_K;  // matrix em bloco diagonal formada pelas matrizes k1,...,kc
		
	 for(i in 1:CN)
	 for(j in 1:CN)
	  L_K[i,j]<-0;

	for(i in 1:C)
		L_K[((i-1)*N+1):(i*N),((i-1)*N+1):(i*N)]<-L_cov_exp_quad_ARD(x,alpha[i],rho[i],delta);
	
	f=mu+sqrt(v/u)*(L_K*eta);	// amostrando de uma t-student mutivariada

	for(i in 1:N)
	   for(j in 1:C)
	 	 fextent[i,j]<-f[(((j-1)*N+1)+(i-1))];
	}

	// distribuições a priori
	
	for(i in 1:C)
		rho[i] ~ gamma((a/2),(a/(2*w)));			//rho~inv_gamma(5,5);
	alpha ~ normal(0,1);				//alpha ~ gamma((a/2),(a/(2*w)));
	eta ~ normal(0,1);
	v ~ gamma(2,0.1);
	u ~ chi_square(v);

	for(i in 1:N)
	 y[i]~multinomial(softmax(to_vector(fextent[i]))); // tinha um a seguindo distribuição normal dentro da função softmax
}
"

# MODELO DE CLASSIFICAÇÃO PROBABILISTICA DE TESTE COM BASE NO FORUM STAN SOBRE CLASSIFICAÇÃO MULTICLASSE E NO ARTIGO DE WILLIAMS E BARBER (1998)

model_multiclasse_pred.GP="

// função para implementar a f.c. exponencial quadrática (ARD)
functions{
	matrix L_cov_exp_quad_ARD(vector[] x, real alpha, vector rho, real delta){
		
		int N=size(x);
		matrix[N, N] K;
		real neg_half=-0.5;
		real sq_alpha=square(alpha);

		for(i in 1:(N-1)){
			K[i,i]=sq_alpha + delta;
			for(j in (i+1):N){
				K[i,j]=sq_alpha*exp(neg_half*dot_self((x[i]-x[j])./rho));
				K[j,i]=K[i,j];				
			}
		}
	
	K[N,N]=sq_alpha+delta;
	return cholesky_decompose(K);
	//return K;
	}

}

data{
	// se o banco de dados tiver entradas multivariadas, então utilizar a configuração abaixo
	
	int<lower=1> N1; // quantidade de pontos de treinamento
	int<lower=1> D; // dimensão do vetor de entrada
	int<lower=1> C; // quantidade de classes
	vector[D] x1[N1]; // dados de treinamento
	int<lower=1> CN1; // dimensão do vetor latente extendido dos dados de treino
	int<lower=1> CN2; // dimensão do vetor latente extendido dos dados de teste
	int<lower=0,upper=1> y1[N1,C]; // matrix de rótulos 
	int<lower=1> N2; // quantidade de pontos de teste
	vector[D] x2[N2]; // dados de entrada dos pontos de teste
	vector<lower=0>[D] rho[C];
	real<lower=0> alpha[C];
}

transformed data{
	real delta=1e-9; // aqui seria o jitter (não colocar mais aspas em comentarios)do modelo de classificação similar a sigma2 na regressão
	int<lower=1> N = N1 + N2;
	int<lower=1> CN= CN1+CN2;
  	vector[D] x[N]; 
	matrix[CN,CN] L_K;
	//matrix[CN,CN] L;

  	for (n1 in 1:N1) x[n1] = x1[n1];
  	for (n2 in 1:N2) x[N1 + n2] = x2[n2];
	
	for(i in 1:CN)
	for(j in 1:CN)
	 L_K[i,j]<-0;
		
	for(i in 1:C)
	  L_K[((i-1)*N+1):(i*N),((i-1)*N+1):(i*N)]<-L_cov_exp_quad_ARD(x,alpha[i],rho[i],delta);
		
	  //L<-cholesky_decompose(L_K);

}

parameters{

	vector[CN]eta;
}

transformed parameters{
	
	vector[CN] f;
	matrix[N,C] fextent;

	f=L_K*eta;			
	for(i in 1:N)
	 for(j in 1:C)
	  fextent[i,j]<-f[(((j-1)*N+1)+(i-1))];
	
}
model{

	// priori para o modelo	
	eta ~ normal(0, 1);

	for(i in 1:N1)
	 y1[i]~multinomial(softmax(to_vector(fextent[i]))); // tinha um a seguindo distribuição normal dentro da função softmax

}
"

# MODELO DE CLASSIFICAÇÃO PROBABILISTICA DE TESTE COM BASE NO FORUM STAN SOBRE CLASSIFICAÇÃO MULTICLASSE E NO ARTIGO DE WILLIAMS E BARBER (1998)

model_multiclasse_pred.TP="

// função para implementar a f.c. exponencial quadrática (ARD)
functions{
	matrix L_cov_exp_quad_ARD(vector[] x, real alpha, vector rho, real delta){
		
		int N=size(x);
		matrix[N, N] K;
		real neg_half=-0.5;
		real sq_alpha=square(alpha);

		for(i in 1:(N-1)){
			K[i,i]=sq_alpha + delta;
			for(j in (i+1):N){
				K[i,j]=sq_alpha*exp(neg_half*dot_self((x[i]-x[j])./rho));
				K[j,i]=K[i,j];				
			}
		}
	
	K[N,N]=sq_alpha+delta;
	return cholesky_decompose(K);
	
	}

}

data{
	// se o banco de dados tiver entradas multivariadas, então utilizar a configuração abaixo
	
	int<lower=1> N1; // quantidade de pontos de treinamento
	int<lower=1> D; // dimensão do vetor de entrada
	int<lower=1> C; // quantidade de classes
	vector[D] x1[N1]; // dados de treinamento
	int<lower=1> CN1; // dimensão do vetor latente extendido dos dados de treino
	int<lower=1> CN2; // dimensão do vetor latente extendido dos dados de teste
	int<lower=0,upper=1> y1[N1,C]; // matrix de rótulos 
	int<lower=1> N2; // quantidade de pontos de teste
	vector[D] x2[N2]; // dados de entrada dos pontos de teste
	vector<lower=0>[D] rho[C];
	real<lower=0> alpha[C];
	real<lower=0> v;
}

transformed data{
	real delta=1e-9; // aqui seria o jitter (não colocar mais aspas em comentarios)do modelo de classificação similar a sigma2 na regressão
	int<lower=1> N = N1 + N2;
	int<lower=1> CN= CN1+CN2;
  	vector[D] x[N]; 
	matrix[CN,CN] L_K;
	vector[CN] mu;

	for(i in 1:CN) mu[i]<-0; // vetor de médias nulo

  	for (n1 in 1:N1) x[n1] = x1[n1];
  	for (n2 in 1:N2) x[N1 + n2] = x2[n2];
	
	for(i in 1:CN)
	for(j in 1:CN)
	 L_K[i,j]<-0;
		
	for(i in 1:C)
	  L_K[((i-1)*N+1):(i*N),((i-1)*N+1):(i*N)]=L_cov_exp_quad_ARD(x,alpha[i],rho[i],delta);
		
}

parameters{
	vector[CN] eta;
  	real<lower=0> u;
}

transformed parameters{
	
	matrix[N,C] fextent;
	vector[CN] f;
	{
	f=mu+sqrt(v/u)*(L_K*eta);	// amostrando de uma t-student mutivariada
		
	for(i in 1:N)
	 for(j in 1:C)
	  fextent[i,j]<-f[(((j-1)*N+1)+(i-1))];
	}
}
model{

	eta ~ normal(0,1);
	u ~ chi_square(v);

	for(i in 1:N1)
	 y1[i]~multinomial(softmax(to_vector(fextent[i]))); // tinha um a seguindo distribuição normal dentro da função softmax

}
"


# banco de dados da coluna vertebral 3 classes

dados=read.table("C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/coluna_vertebral_3_classe.csv",header=T,sep=";")
names(dados)
attach(dados)

y=dados[,8]		# classes dos pacientes

n_out=length(y)
C=3 	#classe

X=cbind(dados[,1:6])	# dados com as variáveis de entrada
X.mc=as.matrix(X)

n_norm=dim(X.mc)

# rotulando as classes

T=array(0,c(n_out,C))

for(i in 1:n_out)
{
	if(y[i]==1)
	{
	T[i,1]=1
	}

	if(y[i]==2)
	{
	 T[i,2]=1

	}
	if(y[i]==3)
	{
	 T[i,3]=1
	}
}

# GRAFICO DE DENSIDADE

dado=as.data.frame(X.mc)

for(i in 1:6)
{
x11()
print(ggplot(data = dado, aes(x = X.mc[,i])) +
  geom_density() + ggtitle('Densidade das variaveis de entrada'))
}

## função para avaliar a taxa de acerto por classe do classificador

ava_classifier2=function(Ytst,Ypred,Ntst)
{
	n=dim(Ytst)
	Nexemplos=array(NA,c(n[2]))

	for (i in 1:n[2]) # contando o numero de exemplos de cada classe no conjunto de dados de teste
	{
		Nexemplos[i]=sum(Ytst[,i])
	}
	Nacertos=array(0,c(n[2])) # um vetor para contar o acerto por classe
		
	for(j in 1:Ntst)
	{
		if(which.max(Ypred[j,])==which.max(Ytst[j,]))
		{
			ident=which.max(Ytst[j,])
			Nacertos[ident]=Nacertos[ident]+1
		}
	}

	return(cbind(Nacertos,Nexemplos))

}


### Rodando o Algoritmo

R=5        # R é o numero de rodadas que o algoritmo terá que fazer
alpha_GP=array(NA,c(C,R))
rho_GP=array(NA,c(C*n_norm[2],R))

v=array(c(2,4,8),c(3))            # graus de liberdade para os quais rodarei o codigo
alpha_TP=array(NA,c(C,4*R))
rho_TP=array(NA,c(C*n_norm[2],4*R))
v_TP=array(NA,c(R))
cont2=0

Pglobal_GP=array(NA,c(R))
Pglobal_TP=array(NA,c(4,R))
cont=0

# guardar os índices a cada rodada

Ind=array(NA,c(R,n_out))

# media das distribuições para gerar os outliers 

media=array(c(150,80,150,150,200,400),c(1,6))

system.time( # para calcular o tempo em que o código permanecer rodando
for(i in 1:R)
{
# dividindo os dados em conjunto de treinamento e de teste

 Ntrn=floor(0.8*n_out)	# 80% de dados para treinamento, N_trn é a quantidade de dados no conjunto de treino
 Ntst=n_out-Ntrn		# 20% dos dados restantes é para teste, N_tst é a quantidade de dados no conjunto de teste

 n=dim(X)
 #	set.seed(203040) # Guardar a semente
	I=randperm(n[1])
	Ind[i,]=I
	Xcomp=X.mc[I,]; ycomp=T[I,] # embaralhamento dos dados de treinamento
#	Xcomp=X.mc[Ind[i,],]; ycomp=y[Ind[i,]] # embaralhamento dos dados de treinamento

#########################################################################################################
########## ESCOLHENDO A P% DE OUTLIERS ###########################################

#  Nout=floor(0.05*Ntrn)	# escolhendo 5% de outliers
# serão perturbados apenas as 6 variavesi de entrada de uma vez para 'i' diferente

#Isamp=1:Ntrn

#for(m in 1:n_norm[2])
#{

# nsamp=sample(Isamp,size=Nout)
#
# for(j in 1:Nout)
# {		 
#	Xcomp[nsamp[j],m]=Xcomp[nsamp[j],m]+media[m]
# }
# x11()
# plot(Xcomp[,m]);points(nsamp,Xcomp[nsamp,m],pch=16,col=3)

#}

# reescalando as variáveis para que tenham media zero e variancia 1

X.norm=array(NA,c(n[1],n[2]))

for(j in 1:n[2])
{
	X.norm[,j]=(Xcomp[,j]-mean(Xcomp[,j]))/sd(Xcomp[,j])
}

	Xtrn=X.norm[1:Ntrn,]
	Ytrn=ycomp[1:Ntrn,]
	Xtst=X.norm[(Ntrn+1):n[1],]
	Ytst=ycomp[(Ntrn+1):n[1],]
 
     #set.seed(202020) # Guardar a semente
	Iout=randperm(Ntrn)
	Xtrn_out=Xtrn[Iout,]; Ytrn_out=Ytrn[Iout,] # embaralhamento dos dados de treinamento

# preparando os dados para entrar no Stan

  N_trn=dim(Xtrn)
  N_tst=dim(Xtst)
  CN=C*N_trn[1]
  CN1=C*N_trn[1]
  CN2=C*N_tst[1]

  dadostan<-list(N=N_trn[1],D=N_trn[2],C=C,CN=CN,x=Xtrn_out,y=Ytrn_out) # dados para classificação

  system.time(# para calcular o tempo que o modelo stan passa rodando
  ajuste.tp<-stan(model_code=model_multiclasse.GP,data=dadostan,iter=220,warmup=100,chains=1)
  )

  print(ajuste.tp,pars=c('alpha','rho'))
  print(c('Rodada',i,"GPC"))
# extraindo os valores dos hiperparâmetros

 amostra_post=extract(ajuste.tp,pars=c('rho','alpha'))
 rho_amostra=amostra_post$rho
 alpha_amostra=amostra_post$alpha

 rhoamost=array(0,c(C,N_trn[2]))
 alphaamost=array(0,c(C))
 Namost=dim(rho_amostra)
 for(j in 1:Namost[1])
  {
	for(l in 1:C)
	{
		rhoamost[l,]=rhoamost[l,]+rho_amostra[j,l,]
	}
	alphaamost=alphaamost+alpha_amostra[j,]
  }
  rho=rhoamost/Namost[1]
  alpha=alphaamost/Namost[1]  

  alpha_GP[,i]=alpha
  
  for(j in 1:C)
  {
  	rho_GP[((j-1)*N_trn[2]+1):(j*N_trn[2]),i]=rho[j,]
  }


# Salvando os gráficos diagnósticos

x11()
print(plot(ajuste.tp,plotfun='trace',pars=c('rho','alpha'),inc_warmup=TRUE)+ ggtitle("GPC"))
#x11()
#print(plot(ajuste.tp,plotfun='ac',pars=c('rho','alpha'),inc_warmup=TRUE)+ ggtitle("GPC "))

# rodando o modelo de TP para v fixo

 for(k in 1:3)
 {
  cont2=cont2+1
	dadostan<-list(N=N_trn[1],D=N_trn[2],C=C,CN=CN,v=v[k],x=Xtrn_out,y=Ytrn_out) # dados para classificação

	system.time(# para calcular o tempo que o modelo stan passa rodando
	ajuste.tp<-stan(model_code=model_multiclasse.TP,data=dadostan,iter=220,warmup=100,chains=1)
	)

	print(ajuste.tp,pars=c('alpha','rho'))
	# extraindo os valores dos hiperparâmetros

 amostra_post=extract(ajuste.tp,pars=c('rho','alpha'))
 rho_amostra=amostra_post$rho
 alpha_amostra=amostra_post$alpha

 rhoamost=array(0,c(C,N_trn[2]))
 alphaamost=array(0,c(C))
 Namost=dim(rho_amostra)
 for(j in 1:Namost[1])
  {
	for(l in 1:C)
	{
		rhoamost[l,]=rhoamost[l,]+rho_amostra[j,l,]
	}
	alphaamost=alphaamost+alpha_amostra[j,]
  }
  rho=rhoamost/Namost[1]
  alpha=alphaamost/Namost[1]  

	alpha_TP[,cont2]=alpha
  
  	for(j in 1:C)
  	{
  		rho_TP[((j-1)*N_trn[2]+1):(j*N_trn[2]),cont2]=rho[j,]
  	}

	print(c('Rodada',i,"v=",v[k]))

	# plotando os gráficos diagnósticos

	x11()
	print(plot(ajuste.tp,plotfun='trace',pars=c('rho','alpha'),inc_warmup=TRUE)+ ggtitle("TPC v fixo"))
#	x11()
#	print(plot(ajuste.tp,plotfun='ac',pars=c('rho','alpha'),inc_warmup=TRUE)+ ggtitle("TPC v fixo"))
}


# rodando o modelo de TP para v sendo estimado a partir dos dados

dadostan<-list(N=N_trn[1],D=N_trn[2],C=C,CN=CN,x=Xtrn_out,y=Ytrn_out) # dados para classificação

	system.time(# para calcular o tempo que o modelo stan passa rodando
	ajuste.tp<-stan(model_code=model_multiclasse.TP.modified,data=dadostan,iter=220,warmup=100,chains=1)
	)
  cont2=cont2+1
	print(ajuste.tp,pars=c('alpha','rho','v'))
	# extraindo os valores dos hiperparâmetros

	amostra_post=extract(ajuste.tp,pars=c('rho','alpha','v'))
	rho_amostra=amostra_post$rho
 	alpha_amostra=amostra_post$alpha
	v_amostra=amostra_post$v

 rhoamost=array(0,c(C,N_trn[2]))
 alphaamost=array(0,c(C))
 Namost=dim(rho_amostra)
 for(j in 1:Namost[1])
  {
	for(l in 1:C)
	{
		rhoamost[l,]=rhoamost[l,]+rho_amostra[j,l,]
	}
	alphaamost=alphaamost+alpha_amostra[j,]
  }
  rho=rhoamost/Namost[1]
  alpha=alphaamost/Namost[1]
  v_TP[i]=mean(v_amostra)  

	alpha_TP[,cont2]=alpha
  
  	for(j in 1:C)
  	{
  		rho_TP[((j-1)*N_trn[2]+1):(j*N_trn[2]),cont2]=rho[j,]
  	}


	print(c('Rodada',i))

# plotando os gráficos diagnósticos

	x11()
	print(plot(ajuste.tp,plotfun='trace',pars=c('rho','alpha'),inc_warmup=TRUE)+ ggtitle("TPC v estimado"))
#	x11()
#	print(plot(ajuste.tp,plotfun='ac',pars=c('rho','alpha'),inc_warmup=TRUE)+ ggtitle("TPC v estimado"))

# Fazendo as predições

 rhopred=matrix(rho_GP[,i],nrow=C,ncol=N_trn[2],byrow=T)
 dadostanpred<-list(N1=N_trn[1],D=N_trn[2],N2=N_tst[1],C=C,CN1=CN1,CN2=CN2,x1=Xtrn_out,rho=rhopred,alpha=alpha_GP[,i],x2=Xtst,y1=Ytrn_out)

 system.time(# para calcular o tempo que o modelo stan passa rodando
 ajustepred<-stan(model_code=model_multiclasse_pred.GP,data=dadostanpred,iter=400,warmup=200,chains=1)
 )

# AVALIANDO O DESEMPENHO DO MODELO

  amostra_pred=extract(ajustepred)

  f=amostra_pred$f

  dim_amostra=dim(f)
  S=dim_amostra[1]

  pi_pred=array(0,c(Ntst[1],C))
  Nacertos=0
  
  for(j in 1:S)
  {
	fpred=array(NA,c(Ntst[1],C))

	for(k in 1:C)
	{
		fpred[,k]=f[j,((k*n[1]+1)-Ntst[1]):(k*n[1])]
	}

	for(m in 1:Ntst[1])
		{	
			pi_pred[m,]=pi_pred[m,]+exp(fpred[m,])/sum(exp(fpred[m,]))	# PROBAB. ACUMULADA	
		}

  }

  pi_pred=pi_pred/S

  for(j in 1:Ntst[1])
	{
		if(which.max(pi_pred[j,])==which.max(Ytst[j,]))
		{
			Nacertos=Nacertos+1
		}	

	}


 Pglobal=100*(Nacertos/Ntst[1])
 print(Pglobal)
 Pglobal_GP[i]=Pglobal
 print(c('Rodada',i, "GPC"))
 print(ava_classifier2(Ytst,pi_pred,Ntst)) # conta o n de acerto por classe

# rodando o modelo de TP para v fixo

for(k in 1:3)
{
  cont=cont+1
  rhopred=matrix(rho_TP[,cont],nrow=C,ncol=N_trn[2],byrow=T)
  dadostanpred<-list(N1=N_trn[1],D=N_trn[2],C=C,CN1=CN1,CN2=CN2,v=v[k],N2=N_tst[1],x1=Xtrn_out,rho=rhopred,alpha=alpha_TP[,cont],x2=Xtst,y1=Ytrn_out)

	system.time(# para calcular o tempo que o modelo stan passa rodando
	ajustepred<-stan(model_code=model_multiclasse_pred.TP,data=dadostanpred,iter=400,warmup=200,chains=1)
	)

# AVALIANDO O DESEMPENHO DO MODELO

  amostra_pred=extract(ajustepred)

  f=amostra_pred$f

  dim_amostra=dim(f)
  S=dim_amostra[1]

  pi_pred=array(0,c(Ntst[1],C))
  Nacertos=0
  
  for(j in 1:S)
  {
	fpred=array(NA,c(Ntst[1],C))

	for(l in 1:C)
	{
		fpred[,l]=f[j,((l*n[1]+1)-Ntst[1]):(l*n[1])]
	}

	for(m in 1:Ntst[1])
		{	
			pi_pred[m,]=pi_pred[m,]+exp(fpred[m,])/sum(exp(fpred[m,]))	# PROBAB. ACUMULADA	
		}

  }

  pi_pred=pi_pred/S

  for(j in 1:Ntst[1])
	{
		if(which.max(pi_pred[j,])==which.max(Ytst[j,]))
		{
			Nacertos=Nacertos+1
		}	

	}

	Pglobal=100*(Nacertos/Ntst[1])

	print(Pglobal)
	Pglobal_TP[k,i]=Pglobal
	print(c('Rodada',i, "v=",v[k]))
	print(ava_classifier2(Ytst,pi_pred,Ntst)) # conta o n de acerto por classe
	}


# rodando o modelo de TP para v sendo estimado a partir dos dados

cont=cont+1
rhopred=matrix(rho_TP[,cont],nrow=C,ncol=N_trn[2],byrow=T)
dadostanpred<-list(N1=N_trn[1],D=N_trn[2],C=C,CN1=CN1,CN2=CN2,v=v_TP[i],N2=N_tst[1],x1=Xtrn_out,rho=rhopred,alpha=alpha_TP[,cont],x2=Xtst,y1=Ytrn_out)

	system.time(# para calcular o tempo que o modelo stan passa rodando
	ajustepred<-stan(model_code=model_multiclasse_pred.TP,data=dadostanpred,iter=400,warmup=200,chains=1)
	)

# AVALIANDO O DESEMPENHO DO MODELO

  amostra_pred=extract(ajustepred)

  f=amostra_pred$f

  dim_amostra=dim(f)
  S=dim_amostra[1]

  pi_pred=array(0,c(Ntst[1],C))
  Nacertos=0
  
  for(j in 1:S)
  {
	fpred=array(NA,c(Ntst[1],C))

	for(k in 1:C)
	{
		fpred[,k]=f[j,((k*n[1]+1)-Ntst[1]):(k*n[1])]
	}

	for(m in 1:Ntst[1])
		{	
			pi_pred[m,]=pi_pred[m,]+exp(fpred[m,])/sum(exp(fpred[m,]))	# PROBAB. ACUMULADA	
		}

  }

  pi_pred=pi_pred/S

  for(j in 1:Ntst[1])
	{
		if(which.max(pi_pred[j,])==which.max(Ytst[j,]))
		{
			Nacertos=Nacertos+1
		}	

	}

	Pglobal=100*(Nacertos/Ntst[1])

	print(Pglobal)
	Pglobal_TP[4,i]=Pglobal
	print(c('Rodada',i, "v estimado"))
	print(ava_classifier2(Ytst,pi_pred,Ntst)) # conta o n de acerto por classe	

	write.table(Pglobal_TP,"C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/scripts/simulacoes_multiclasse/Pglobal_TP.txt")
	write.table(Pglobal_GP,"C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/scripts/simulacoes_multiclasse/Pglobal_GP.txt")

	write.table(alpha_GP,"C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/scripts/simulacoes_multiclasse/alpha_mult_GP.txt")
	write.table(rho_GP,"C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/scripts/simulacoes_multiclasse/rho_mult_GP.txt")

	write.table(alpha_TP,"C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/scripts/simulacoes_multiclasse/alpha_mult_TP.txt")
	write.table(rho_TP,"C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/scripts/simulacoes_multiclasse/rho_mult_TP.txt")
	write.table(v_TP,"C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/scripts/simulacoes_multiclasse/v.txt")

}
)

resumo=data.frame(Media=rep(NA,5),Desv.pad=rep(NA,5),Mediana=rep(NA,5),Maximo=rep(NA,5),Minimo=rep(NA,5))

STAT=list(Media=mean(Pglobal_GP),Desv.pad=sd(Pglobal_GP),Mediana=median(Pglobal_GP),Maximo=max(Pglobal_GP),Minimo=min(Pglobal_GP))
resumo[1,]=STAT

for(i in 1:4)
{
STAT_TP=list(Media=mean(Pglobal_TP[i,]),Desv.pad=sd(Pglobal_TP[i,]),Mediana=median(Pglobal_TP[i,]),Maximo=max(Pglobal_TP[i,]),Minimo=min(Pglobal_TP[i,]))
resumo[(i+1),]=STAT_TP
}

v_degree=c("GPC",2,4,8,"vchap")
resut=cbind(resumo,v_degree)

# Salvando os dados


Pglobal_TP=read.table("C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/scripts/simulacoes_multiclasse/Pglobal_TP5p.txt")
Pglobal_GP=read.table("C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/scripts/simulacoes_multiclasse/Pglobal_GP5p.txt")
Pglobal_GP=Pglobal_GP[,1]
Pglobal_TP=as.matrix(Pglobal_TP)
v_TP=read.table("C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/scripts/simulacoes_multiclasse/v5p.txt")
v_TP=v_TP[,1]

