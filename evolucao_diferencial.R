#PROFESSOR CINIRO NAMETALA
#IFMG - CAMPUS BAMBUI
#EVOLUCAO DIFERENCIAL
#-------------------------------------------------------

#limpar workspace
rm(list=ls())

#limpar tela
cat('\014')

#bibliotecas
library("plot3D")

#FUNCOES
#ESFERA
#otimo: 0..
f.sphere <- function(x) 
{
  return(sum(x^2))
}

#RASTRIGIN
#otimo: 0..
f.rastrigin <- function(x)
{
  d <- length(x)
  sum <- sum(x^2 - 10*cos(2*pi*x))
  y <- 10*d + sum
  return(y)
}

#ROSENBROCK
#otimo: 1..
f.rosenbrock <- function(x)
{
  d <- length(x)
  xi <- x[1:(d-1)]
  xnext <- x[2:d]
  
  sum <- sum(100*(xnext-xi^2)^2 + (xi-1)^2)
  
  y <- sum
  return(y)
}

#PEAKS
#otimo: 1..
f.peaks <- function(x)
{
  return(x[1] * exp(-(x[1]^2 + x[2]^2)))
}

#-----------------------------------------------------------------------------
#PARAMETRIZACAO DO ALGORITMO
f <- 1                    #fator de escala F
fmult <- 5
chancefmult <- 0.4
txcruz <- 0.7               #taxa de cruzamento

#PARAMETRIZACAO DO EXPERIMENTO
dimindividuo <- 3           #quantidade de dimensoes do problema
#funcaofit <- f.rastrigin  #funcao objetivo
funcaofit <- f.sphere
maxgeracoes <- 300             #quantidade maxima de geracoes
tampop <- 100       #tamanho da populacao
linf <- -5.12                 #limite inferior do espaco de busca
lsup <- 5.12
res <- 0.6

#EXIBICAO DE GRAFICOS
plotaGraficos <- TRUE

#GERANDO A POPULACAO INICIAL
#gera a matriz de posicoes das particulas
populacao <- matrix(
  runif(tampop*dimindividuo, linf, lsup), 
  tampop, 
  dimindividuo)

#AVALIANDO OS INDIVIDUOS DA POPULACAO INICIAL
#avalia a qualidade de cada particula com a funcao objetivo
fitness <- rep(NA,tampop)
for (i in 1:tampop)
{
  fitness[i] <- funcaofit(populacao[i,])
}

#PLOTA A PRIMEIRA DISTRIBUICAO DE INDIVIDUOS NO ESPACO
if (plotaGraficos == TRUE)
{
  #plotando a funcoes
  x <- seq(linf,lsup,res)
  y <- x
  z <- matrix(NA, length(x), length(y))
  for (i in 1:length(x))
  {
    for (j in 1:length(y))
    {
      z[i,j] <- funcaofit(c(x[i],y[j]))
    }
  }
  
  persp3D(x,y,z,theta=45,phi=30)
  
  contour(x,y,z)
  par(new=T)
  
  #plotando as particulas
  plot(populacao, 
       xlim=c(linf,lsup), 
       ylim=c(linf,lsup), 
       pch=20, 
       xlab="", ylab="", col="blue", xaxt='n', yaxt='n')
}

#vetor para armazenar a fitness media e a fitness do melhor a cada geracao
fitnessmedia <- rep(NA,maxgeracoes)
fitnessmelhor <- rep(NA,maxgeracoes)

#INICIALIZANDO O ALGORITMO
for (igeracao in 1:maxgeracoes)
{
  print(igeracao)
  #movimenta as particulas no espaco
  for (iindividuo in 1:tampop)
  {
    
    #SELECAO----
    #SORTEIA 3 INDIVIDUOS (x0, x1 e x2) NA POPULACAO
    indsorteados <- sample(1:tampop, 4, replace=FALSE)
    ind_x0 <- populacao[indsorteados[1],]
    ind_x1 <- populacao[indsorteados[2],]
    ind_x2 <- populacao[indsorteados[3],]
    ind_target <- populacao[indsorteados[4],]
    
    #MUTACAO----
    #CALCULA O VETOR DE DIFERENCAS
    #APLICA MUTACAO COM FATOR DE CONTRICAO f e soma ao X0
    if (runif(1,0,1)<=chancefmult)
    {
      ind_mutado <- ind_x0 + ((ind_x1 - ind_x2) * (f*(fmult*runif(1,0,1))))
    }
    else
    {
      ind_mutado <- ind_x0 + ((ind_x1 - ind_x2) * f)
    }
    
    #CRUZAMENTO----
    #cruzamento ponto a ponto
    ind_trial <- rep(NA, dimindividuo)
    for (igene in 1:dimindividuo)
    {
      #se valor sorteado menor que a taxa de cruzamento
      #entao trial recebe material do target senao do mutado
      #alta taxa de cruzamento priveligia o vetor target
      if (runif(1,0,1)<=txcruz)
      {
        ind_trial[igene] <- ind_target[igene]
      }
      else
      {
        ind_trial[igene] <- ind_mutado[igene]
      }
    }
    
    #COMPETICAO----
    #trial e target competem entre si
    #o melhor sera incluido na novapopulacao
    fit_trial <- funcaofit(ind_trial)
    fit_target <- funcaofit(ind_target)
    
    #fitness menor e melhor pois estamos minimizando
    if (fit_trial < fit_target)
    {
      populacao[iindividuo,] <- ind_trial
      fitness[iindividuo] <- fit_trial
    }
    else
    {
      populacao[iindividuo,] <- ind_target
      fitness[iindividuo] <- fit_target
    }
    
  } #individuo
  
  #guarda a fitness media da iteracao iter
  fitnessmedia[igeracao] <- mean(fitness)
  fitnessmelhor[igeracao] <- min(fitness)
  
  #plotando a nova populacao no espaco de busca
  if (plotaGraficos == TRUE)
  {
    contour(x,y,z)
    par(new=T)
    plot(populacao, 
         xlim=c(linf,lsup), 
         ylim=c(linf,lsup), 
         pch=20, 
         xlab="", ylab="", col="blue")
  }
  
} #geracao

indmelhor <- which(fitness == min(fitness))[1]
melhor <- populacao[indmelhor,]
fitmelhor <- min(fitness)

#exibi informacoes da melhor particula (ultimo gbest)
print(round(melhor, 4))
print(round(fitmelhor,4))
plot(fitnessmedia, type="l", xlab="Geracoes", ylab="Fitness media", main="Fitness media a cada geracao")
plot(fitnessmelhor, type="l", xlab="Geracoes", ylab="Fitness", main="Fitness do melhor a cada geracao")
