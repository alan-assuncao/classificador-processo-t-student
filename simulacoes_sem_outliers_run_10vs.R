### IMPLEMENTAÇÃO STAN PARA O CASO MULTICLASSE

rm(list=ls())# limpar a memória do PC

# PACOTES QUE SERÃO UTILIZADOS NO DECORRER DO SCRIPT

library('rstan')
library('MASS')
library('pracma')
library('reshape2')
library('ggplot2')
library('R.matlab')	# para ler arquivos matlab
library('fields')	# para calcular distancias

####################################################################################
# MODELO DE CLASSIFICAÇÃO BINARIA COM PROCESSOS T-STUDENT

model_logistic.TP="

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
	
	int<lower=1> N;
	int<lower=1> D;
	vector[D] x[N];
	int<lower=0,upper=1> z[N];
	real<lower=0> v; // graus de liberdade da distribuição t-student
}

transformed data{
	real delta=1e-9; // aqui seria o jitter (não colocar mais aspas em comentarios)do modelo de classificação similar a sigma2 na regressão
	real w=1;
	real a=1;

	vector[N] mu;
	for(i in 1:N) mu[i]<-0; // vetor de médias nulo
}

parameters{
	// se eu utilizar entradas multivariadas, utilizar essa configuração abaixo
	
	vector<lower=0>[D] rho;
	real<lower=0> alpha;
	vector[N] eta;
	real<lower=0> u;
}

model{
	
	vector[N] f;
	{
		matrix[N,N] L_K=L_cov_exp_quad_ARD(x,alpha,rho,delta);
		f=mu+sqrt(v/u)*(L_K*eta);	// amostrando de uma t-student mutivariada
	}


	// distribuições a priori
	
	rho ~ gamma((a/2),(a/(2*w)));			//rho~inv_gamma(5,5);
	alpha ~ normal(0,1);
	//v ~ gamma(2,0.1);
	eta ~ normal(0,1);
	u ~ chi_square(v);

	z ~ bernoulli_logit(f);
}
"

####################################################################################
# MODELO DE CLASSIFICAÇÃO BINARIA COM PROCESSOS T-STUDENT

model_logistic.TP.modified="

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
	
	int<lower=1> N;
	int<lower=1> D;
	vector[D] x[N];
	int<lower=0,upper=1> z[N];
	
}

transformed data{
	real delta=1e-9; // aqui seria o jitter (não colocar mais aspas em comentarios)do modelo de classificação similar a sigma2 na regressão
	real w=1;
	real a=1;	

	vector[N] mu;
	for(i in 1:N) mu[i]<-0; // vetor de médias nulo
}

parameters{
	// se eu utilizar entradas multivariadas, utilizar essa configuração abaixo
	
	vector<lower=0>[D] rho;
	real<lower=0> alpha;
	vector[N] eta;
	real<lower=0> u;
	real<lower=0> v; // graus de liberdade da distribuição t-student
}

model{

	
	vector[N] f;
	{
		matrix[N,N] L_K=L_cov_exp_quad_ARD(x,alpha,rho,delta);
		f=mu+sqrt(v/u)*(L_K*eta);	// amostrando de uma t-student mutivariada
	}

	// distribuições a priori
	
	rho ~ gamma((a/2),(a/(2*w)));			//rho~inv_gamma(5,5);
	alpha ~ normal(0,1);

	v ~ gamma(2,0.1);
	eta ~ normal(0,1);
	u ~ chi_square(v);

	z ~ bernoulli_logit(f);
}
"

####################################################################################
# MODELO DE CLASSIFICAÇÃO BINARIA COM PROCESSOS GAUSSIANO

model_logistic.GP="

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
	
	int<lower=1> N;
	int<lower=1> D;
	vector[D] x[N];
	int<lower=0,upper=1> z[N];

}

transformed data{
	real delta=1e-9; // aqui seria o jitter (não colocar mais aspas em comentarios)do modelo de classificação similar a sigma2 na regressão
	real w=1;
	real a=1;
}

parameters{
	// se eu utilizar entradas multivariadas, utilizar essa configuração abaixo
	
	vector<lower=0>[D] rho;
	real<lower=0> alpha;
	//real<lower=0> a;
	vector[N] eta;
}

model{

		// se as entradas x forem multivariadas, então utilizar a configuração abaio
	vector[N] f;
	{
		matrix[N,N] L_K=L_cov_exp_quad_ARD(x,alpha,rho,delta);
		f=L_K*eta;
	}

	// distribuições a priori

	rho ~ gamma((a/2),(a/(2*w)));			//rho~inv_gamma(5,5);
	alpha ~ normal(0,1);
	//a ~ normal(0,1);
	eta ~ normal(0,1);

	z ~ bernoulli_logit(f);
}
"

#############################################################################################
# Implementando a distribuição preditiva do modelo de classificação binário em stan #########

model_logistic_pred.TP="

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

data {
  int<lower=1> D; //dimensão do vetor de entrada
  int<lower=1> N1;
  vector[D] x1[N1]; // vetor de dimensão N1 x D
  int<lower=0, upper=1> z1[N1];
  int<lower=1> N2;
  vector[D] x2[N2]; // vetor de dimensão N2 x D
  real<lower=0> v;
  vector<lower=0>[D] rho;
  real<lower=0> alpha;
}
transformed data {
  real delta = 1e-9;
  int<lower=1> N = N1 + N2;
  vector[D] x[N]; // original era isso: real x[N]
  matrix[N,N] L_K;
  vector[N] mu;

  for (n1 in 1:N1) x[n1] = x1[n1];
  for (n2 in 1:N2) x[N1 + n2] = x2[n2];

  for(i in 1:N) mu[i]<-0; // vetor de médias nulo

  L_K=L_cov_exp_quad_ARD(x,alpha,rho,delta);
}

parameters {
  vector[N] eta;
  real<lower=0> u;
}

transformed parameters {

	vector[N] f;
	{
	f=mu+sqrt(v/u)*(L_K*eta);	// amostrando de uma t-student mutivariada
	}
}

model {

	eta ~ normal(0,1);
	u ~ chi_square(v);

  z1 ~ bernoulli_logit(f[1:N1]); // tinha um a seguindo uma N(0,1) dentro do logit
}
"

#############################################################################################
# Implementando a distribuição preditiva do modelo de classificação binário em stan #########

model_logistic_pred.GP="

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

data {
  int<lower=1> D; //dimensão do vetor de entrada
  int<lower=1> N1;
  vector[D] x1[N1]; // vetor de dimensão N1 x D
  int<lower=0, upper=1> z1[N1];
  int<lower=1> N2;
  vector[D] x2[N2]; // vetor de dimensão N2 x D
  vector<lower=0>[D] rho;
  real<lower=0> alpha;
}
transformed data {
  real delta = 1e-9;
  int<lower=1> N = N1 + N2;
  vector[D] x[N]; // original era isso: real x[N]
  matrix[N,N] L_K;

  for (n1 in 1:N1) x[n1] = x1[n1];
  for (n2 in 1:N2) x[N1 + n2] = x2[n2];

  L_K=L_cov_exp_quad_ARD(x,alpha,rho,delta);
}
parameters {

 vector[N] eta;
}
transformed parameters {
  vector[N] f;

// se as entradas x forem multivariadas, então utilizar a configuração abaio

	{
		f=L_K*eta;
	}
}

model {

  // priori para o modelo
	eta ~ normal(0, 1);

  z1 ~ bernoulli_logit(f[1:N1]); // tinha um a seguindo uma N(0,1) dentro do logit
}
"

# banco de dados da coluna vertebral 2 classes

dados=read.table("C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/coluna_vertebral_2_classe.csv",header=T,sep=";")
names(dados)
attach(dados)

y=dados[,8]		# classes dos pacientes

n_out=length(y)

X=cbind(dados[,1:6])	# dados com as variáveis de entrada
X.mc=as.matrix(X)


# banco de dados da Wisconsin Diagnostic Breast Cancer (WDBC) 2 classes

#dados=read.table("C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/Dados_WDBC_diagnostico_medico/WDBC.csv",header=T,sep=";")
#names(dados)
#attach(dados)

#y=dados[,3]		# classes dos pacientes

#n_out=length(y)

#X=cbind(dados[,4:33])	# dados com as variáveis de entrada
#X.mc=as.matrix(X)

#####################################################################################################################################


# reescalando as variáveis para que tenham media zero e variancia 1

n_norm=dim(X.mc)
X.norm=array(NA,c(n_norm[1],n_norm[2]))

for(i in 1:n_norm[2])
{
	X.norm[,i]=(X.mc[,i]-mean(X.mc[,i]))/sd(X.mc[,i])
}

X.mc=X.norm


### Rodando o Algoritmo

R=15        # R é o numero de rodadas que o algoritmo terá que fazer
alpha_GP=array(NA,c(R))
rho_GP=array(NA,c(n_norm[2],R))

v=array(c(2,4,8),c(3))            # graus de liberdade para os quais rodarei o codigo
alpha_TP=array(NA,c(4,R))
rho_TP=array(NA,c(n_norm[2],4*R))
v_TP=array(NA,c(R))
cont2=0

Pglobal_GP=array(NA,c(R))
Pglobal_TP=array(NA,c(4,R))
cont=0

system.time( # para calcular o tempo em que o código permanecer rodando
for(i in 1:R)
{
# dividindo os dados em conjunto de treinamento e de teste

 Ntrn=floor(0.8*n_out)	# 80% de dados para treinamento, N_trn é a quantidade de dados no conjunto de treino
 Ntst=n_out-Ntrn		# 20% dos dados restantes é para teste, N_tst é a quantidade de dados no conjunto de teste

 n=dim(X)
 #	set.seed(203040) # Guardar a semente
	I=randperm(n[1])
	Xcomp=X.mc[I,]; ycomp=y[I] # embaralhamento dos dados de treinamento

	Xtrn=Xcomp[1:Ntrn,]
	Ytrn=ycomp[1:Ntrn]
	Xtst=Xcomp[(Ntrn+1):n[1],]
	Ytst=ycomp[(Ntrn+1):n[1]]

# preparando os dados para entrar no Stan

  N_trn=dim(Xtrn)
  N_tst=dim(Xtst)

  dadostan<-list(N=N_trn[1],D=N_trn[2],x=Xtrn,z=Ytrn) # dados para classificação

  system.time(# para calcular o tempo que o modelo stan passa rodando
  ajuste.tp<-stan(model_code=model_logistic.GP,data=dadostan,iter=600,warmup=200,chains=1)
  )

  print(ajuste.tp,pars=c('alpha','rho'))
  print(c('Rodada',i,"GPC"))
# extraindo os valores dos hiperparâmetros

 amostra_post=extract(ajuste.tp,pars=c('rho','alpha'))
 rho_amostra=amostra_post$rho
 alpha_amostra=amostra_post$alpha

 rho=array(NA,c(N_trn[2]))
 for(j in 1:N_trn[2])
  {
		rho[j]=mean(rho_amostra[,j])
	}

  alpha_GP[i]=mean(alpha_amostra)
  rho_GP[,i]=rho

# Salvando os gráficos diagnósticos

x11()
print(plot(ajuste.tp,plotfun='trace',pars=c('rho','alpha'),inc_warmup=TRUE)+ ggtitle("GPC"))
#x11()
#print(plot(ajuste.tp,plotfun='ac',pars=c('rho','alpha'),inc_warmup=TRUE)+ ggtitle("GPC "))

# rodando o modelo de TP para v fixo

 for(k in 1:3)
 {
  cont2=cont2+1
	dadostan<-list(N=N_trn[1],D=N_trn[2],v=v[k],x=Xtrn,z=Ytrn) # dados para classificação

	system.time(# para calcular o tempo que o modelo stan passa rodando
	ajuste.tp<-stan(model_code=model_logistic.TP,data=dadostan,iter=600,warmup=200,chains=1)
	)

	print(ajuste.tp,pars=c('alpha','rho'))
	# extraindo os valores dos hiperparâmetros

	amostra_post=extract(ajuste.tp,pars=c('rho','alpha'))

	rho_amostra=amostra_post$rho
	rho=array(NA,c(N_trn[2]))
	for(j in 1:N_trn[2])
  	{
		rho[j]=mean(rho_amostra[,j])
	}


	alpha_amostra=amostra_post$alpha
	alpha=mean(alpha_amostra)

	alpha_TP[k,i]=alpha
	rho_TP[,cont2]=rho
	print(c('Rodada',i,"v=",v[k]))

	# plotando os gráficos diagnósticos

	x11()
	print(plot(ajuste.tp,plotfun='trace',pars=c('rho','alpha'),inc_warmup=TRUE)+ ggtitle("TPC v fixo"))
#	x11()
#	print(plot(ajuste.tp,plotfun='ac',pars=c('rho','alpha'),inc_warmup=TRUE)+ ggtitle("TPC v fixo"))
}


# rodando o modelo de TP para v sendo estimado a partir dos dados

dadostan<-list(N=N_trn[1],D=N_trn[2],x=Xtrn,z=Ytrn) # dados para classificação

	system.time(# para calcular o tempo que o modelo stan passa rodando
	ajuste.tp<-stan(model_code=model_logistic.TP.modified,data=dadostan,iter=600,warmup=200,chains=1)
	)
  cont2=cont2+1
	print(ajuste.tp,pars=c('alpha','rho','v'))
	# extraindo os valores dos hiperparâmetros

	amostra_post=extract(ajuste.tp,pars=c('rho','alpha','v'))

	rho_amostra=amostra_post$rho
	rho=array(NA,c(N_trn[2]))
	for(j in 1:N_trn[2])
  	{
		rho[j]=mean(rho_amostra[,j])
	}

	alpha_amostra=amostra_post$alpha
	v_amostra=amostra_post$v
	alpha=mean(alpha_amostra)

	alpha_TP[4,i]=alpha
	rho_TP[,cont2]=rho
  	v_TP[i]=mean(v_amostra)
	print(c('Rodada',i))

# plotando os gráficos diagnósticos

	x11()
	print(plot(ajuste.tp,plotfun='trace',pars=c('rho','alpha'),inc_warmup=TRUE)+ ggtitle("TPC v estimado"))
#	x11()
#	print(plot(ajuste.tp,plotfun='ac',pars=c('rho','alpha'),inc_warmup=TRUE)+ ggtitle("TPC v estimado"))

# Fazendo as predições

 dadostanpred<-list(N1=N_trn[1],D=N_trn[2],N2=N_tst[1],x1=Xtrn,rho=rho_GP[,i],alpha=alpha_GP[i],x2=Xtst,z1=Ytrn)

 system.time(# para calcular o tempo que o modelo stan passa rodando
 ajustepred<-stan(model_code=model_logistic_pred.GP,data=dadostanpred,iter=1000,warmup=200,chains=1)
 )

# AVALIANDO O DESEMPENHO DO MODELO

  amostra_pred=extract(ajustepred)

  f=amostra_pred$f
  fpred=f

  dim_amostra=dim(f)
  S=dim_amostra[1]

  pi_pred=array(0,c(Ntst[1]))
  zpred=array(0,c(Ntst[1]))
  Nacertos=0
  f=f[,(Ntrn[1]+1):n[1]]

  for(j in 1:S)
  {

	 pi_pred=pi_pred+sigmoid(f[j,])
  }

  pi_pred=pi_pred/S

  for(j in 1:Ntst[1])
	{
		if(pi_pred[j]>0.5)
		{
			zpred[j]=1
		}else{
			zpred[j]=0
		}

		if(zpred[j]==Ytst[j])
		{
			Nacertos=Nacertos+1
		}
	}


 Pglobal=100*(Nacertos/Ntst[1])
 print(Pglobal)
 Pglobal_GP[i]=Pglobal
 print(c('Rodada',i, "GPC"))

# rodando o modelo de TP para v fixo

for(k in 1:3)
{
  cont=cont+1
  dadostanpred<-list(N1=N_trn[1],D=N_trn[2],v=v[k],N2=N_tst[1],x1=Xtrn,rho=rho_TP[,cont],alpha=alpha_TP[k,i],x2=Xtst,z1=Ytrn)

	system.time(# para calcular o tempo que o modelo stan passa rodando
	ajustepred<-stan(model_code=model_logistic_pred.TP,data=dadostanpred,iter=1000,warmup=200,chains=1)
	)

	# AVALIANDO O DESEMPENHO DO MODELO

	amostra_pred=extract(ajustepred)

	f=amostra_pred$f
	dim_amostra=dim(f)
	S=dim_amostra[1]

	pi_pred=array(0,c(Ntst[1]))
	zpred=array(0,c(Ntst[1]))
	Nacertos=0
	f=f[,(Ntrn[1]+1):n[1]]

	for(j in 1:S)
	{
		pi_pred=pi_pred+sigmoid(f[j,])
	}

	pi_pred=pi_pred/S

	for(j in 1:Ntst[1])
	{
		if(pi_pred[j]>0.5)
		{
			zpred[j]=1
		}else{
			zpred[j]=0
		}

		if(zpred[j]==Ytst[j])
		{
			Nacertos=Nacertos+1
		}
	}

	Pglobal=100*(Nacertos/Ntst[1])

	print(Pglobal)
	Pglobal_TP[k,i]=Pglobal
	print(c('Rodada',i, "v=",v[k]))
	}


# rodando o modelo de TP para v sendo estimado a partir dos dados

cont=cont+1
dadostanpred<-list(N1=N_trn[1],D=N_trn[2],v=v_TP[i],N2=N_tst[1],x1=Xtrn,rho=rho_TP[,cont],alpha=alpha_TP[4,i],x2=Xtst,z1=Ytrn)

	system.time(# para calcular o tempo que o modelo stan passa rodando
	ajustepred<-stan(model_code=model_logistic_pred.TP,data=dadostanpred,iter=1000,warmup=200,chains=1)
	)

	# AVALIANDO O DESEMPENHO DO MODELO

	amostra_pred=extract(ajustepred)

	f=amostra_pred$f

	dim_amostra=dim(f)
	S=dim_amostra[1]

	pi_pred=array(0,c(Ntst[1]))
	zpred=array(0,c(Ntst[1]))
	Nacertos=0
	f=f[,(Ntrn[1]+1):n[1]]

	for(j in 1:S)
	{
		pi_pred=pi_pred+sigmoid(f[j,])
	}

	pi_pred=pi_pred/S

	for(j in 1:Ntst[1])
	{
		if(pi_pred[j]>0.5)
		{
			zpred[j]=1
		}else{
			zpred[j]=0
		}

		if(zpred[j]==Ytst[j])
		{
			Nacertos=Nacertos+1
		}
	}

	Pglobal=100*(Nacertos/Ntst[1])

	print(Pglobal)
	Pglobal_TP[4,i]=Pglobal
	print(c('Rodada',i, "v estimado"))

	# Salvando os dados	

	write.table(Pglobal_TP,"C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/scripts/simulacoes_binario/Pglobal_TP.txt")
	write.table(Pglobal_GP,"C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/scripts/simulacoes_binario/Pglobal_GP.txt")
	
	write.table(alpha_GP,"C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/scripts/simulacoes_binario/alpha_bin_GP.txt")
	write.table(rho_GP,"C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/scripts/simulacoes_binario/rho_bin_GP.txt")

	write.table(alpha_TP,"C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/scripts/simulacoes_binario/alpha_bin_TP.txt")
	write.table(rho_TP,"C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/scripts/simulacoes_binario/rho_bin_TP.txt")
	write.table(v_TP,"C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/scripts/simulacoes_binario/v.txt")

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


Pglobal_TP=read.table("C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/scripts/simulacoes/Pglobal_TP.txt")
Pglobal_GP=read.table("C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/scripts/simulacoes/Pglobal_GP.txt")
Pglobal_GP=Pglobal_GP[,1]
Pglobal_TP=as.matrix(Pglobal_TP)
v_TP=read.table("C:/Users/ALAN/Google Drive/MESTRADO/UFC/DISSERTAÇÃO/dados_e_scripts/scripts/simulacoes/v.txt")
v_TP=v_TP[,1]

# Salvando os dados
