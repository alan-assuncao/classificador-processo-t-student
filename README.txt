# classificador-processo-t-student
Códigos do classificador de processo t-Student, com função de ligação logit, produzidos no trabalho de dissertação"Processos t-Student em classificação"

O primeiro arquivo "simulacoes_sem_outliers_run_10vs" apresenta as versões em stan do classificadores binário de processo Gaussiano e processo t-Student aplicado aos dado de coluna 
vertebral encontrado no repositório UCI irvine. O script monta um loop de 10 replicações para o uso dos classificadores, em que cada replicação, o banco de dados é partici
onado aleatóriamente em dois conjuntos, um de treino e um de teste, para avaliar a taxa de acerto de cada classificador. Ao final, as taxas de acerto são guardadas e depois
se tira a média e o desvio padrão delas para avaliar a precisão dos classificadores.

O segundo arquivo "simulacoes_multiclass_sem_input_outliers_run_5vs.R" funciona da mesma forma que no primeiro, a diferença é que os classificadores agora são multiclass.
