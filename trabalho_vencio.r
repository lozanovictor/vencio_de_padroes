# Plataform

# *Affymetrix Human Gene 1.0 ST Array*


# Análise dos Dados

## Download dos dados
```{r}

# carregar bibliotecas

library(Biobase)
library(GEOquery)
library(limma)

# carregar dados do GEO
 
gset <- getGEO("GSE40392", GSEMatrix =TRUE) # faz o download do estudo ( fui utilizado um código diferente para baixar o mesmo arquivo,no caso estudado : GDS4974 )
gset <- gset[[1]] #extrai o primeiro elemento da lista
dat <- exprs(gset) #extrai os dados de expressão e atribui a uma matriz
dat.anno <- pData(gset) #extrai os dados de anotação e salva em um dataframe



## Conferir se os dados estão da forma necessaria
```{r}

dim(dat) #dimensões dos dados
View(dat.anno) #ve os metadados
colnames(dat) == rownames(dat.anno) #checa ordenacao, pois em dat.anno as 10 primeiras amostras são de oene e as 8 ultimas de expostos ao Ni
table(is.na(dat)) #checa se tem NA, pois pode atrapalhar nas analises

```

## Busca outliers
```{r}
# Busca Outliers
dat.mean <- apply(dat,1,mean) # Calcula a media
dat.sd <- sqrt(apply(dat,1,var)) # Calcula o desvio padrao
dat.cv <- dat.sd/dat.mean # Calcula o coeficiente de variacao
plot(dat.mean,dat.cv,main="Sample CV vs. Mean",xlab="Mean",ylab="CV",cex=1) #Plot do grafico

# Separa em grupos
dat.anno[1:10, "color"] <- "red" #atribui uma cor ao primeiro grupo de amostras - environmental exposure (EE)
dat.anno[11:18, "color"] <- "green" #atribui cor ao segundo grupo - occupational exposure nickel (OEN)
dat.anno[colnames(dat), "color"] #checa o atribuicao

```

## Normalizacao dos dados
```{r}

limma::plotMA(dat) # distribuicao original
abline(h=0, col="dark grey")

# Normalizacao dos dados pelo método loess
dat.norm<-normalizeCyclicLoess(dat)

limma::plotMA(dat.norm) #Plot apos normalizacao
abline(h=0, col="dark grey")

```
# Podemos ver o conjunto de pontos bem mais homogeneo e alinhado a linha central do grafico, mostrando a importancia da normalizacao


## Dendograma
```{r}

colnames(dat.norm) = dat.anno$geo_accession # altera o nome das amostras
plot(hclust(dist(t(dat.norm), method = "euclidean")),main = "Hclust Dendogram");

```
# Fica claro que 7 das 8 amostras de pessoas expostas ao Ni são agrupadas na mesma categoria e as demais amostras agrupadas separadamente (amostra GSM992710, exposta a Ni no trabalho, é agrupada com as amostras de referencia).

## Correlacao
```{r}

# busca por correlacao
dat.cor <- cor(dat.norm) # calcula a correlacao
image (dat.cor, axes=F, main = "Correlation matrix") # gera a imagem da matriz

```
# Existem poucos genes correlacionados no plot, a maioria relacionando o grupo ee vs ee (controle), mas tambem alguns relacionando os dois grupos. A maioria dos dados tem baixa correlacao.

## Plot Profile
```{r}

selected = c(1:18)
plot(c(1,ncol(subdata)), range(subdata[selected,]), type='n',
main="Profile Plot - Genes", xlab="Samples", ylab="Expression")
for (i in 1:length(selected)) {
subdata.y = as.numeric(subdata[selected[i],])
lines(c(1:ncol(subdata)), subdata.y, col=i)
}
text(subdata.y, y=NULL, cex = 1)

```
# Podemos avaliar a variacao de expressao nesse grafico, ficando claro que existem diferencas na expressao entre as amostras e que isso pode estar relacionado com a diferenca dos grupos e sera mais explorada neste estudo.

### Resultados

## Expressao Diferencial
```{r}

# separacao dos dados pelo indice
data.ee = data[,11:18]  #ee
data.oen = data[,1:10]  #oen
# calculo do pvalue pela funcao abaixo
t.test.all.genes = function(data, s1, s2){
    oen_sample = data[colnames(s1)]
    ee_sample = data[colnames(s2)]
    oen_sample = as.numeric(oen_sample)
    ee_sample = as.numeric(ee_sample)
    t.out = t.test(oen_sample, ee_sample, var.equal=T)
    out = as.numeric(t.out$p.value)
    return(out)
}
# execucao do teste t
t.test.run = apply(data ,1, t.test.all.genes, s1=data.oen, s2=data.ee)
# sumario dos dados, para observar a media, quartis, etc
summary(t.test.run)
# histograma
hist(t.test.run, col="red")

```

## Volcano Plot
```{r}

oen.m = apply(data.oen, 1, mean, na.rm=T)
ee.m = apply(data.ee, 1, mean, na.rm=T)
fold = log2(ee.m/oen.m)
p.valor = -1*log10(t.test.run)
plot(range(p.valor), range(fold), type="n", xlab="-1*log(p-value)", ylab="fold change",ylim=c(-0.6,0.6),xlim=c(0,10), main="Volcano Plot")
points(p.valor, fold, col="black")
points(p.valor[(p.valor>0.5 & fold>0.4)], fold[(p.valor>0.5 & fold>0.4)], col="red")
points(p.valor[(p.valor>0.5 & fold< -0.4)], fold[(p.valor>0.5 & fold< -0.4)], col="blue", pch=15)
abline(v=1)
abline(h=-0.4)
abline(h=0.4)

```
# Podemos ver que varios genes apresentam aumento na expressao e alguns apresentam diminuicao da expressao. Esses genes poderao ser utilizados para caracterizar o grupo estudado.

## PCA
```{r}

data.pca.n <- prcomp(t(dat.norm))
data.pca.n <- data.pca.n$x[,1:2]
plot(data.pca.n[,1],data.pca.n[,2], col=as.factor(dat.anno$color), main="PCA",xlab="PC1", ylab="PC2")
legend(1,1,c("O. exposure - black","E. exposure - red"),lty=c(1,1), col=as.factor(dat.anno$color),cex=1, bty="n")

```

# Podemos ver claramente a tendencia de agrupamento entre as amostras similares, separando o grupo oen, acumulado na parte superior direita do plot (em preto), e o grupo ee, acumulado na parte inferior esquerda (vermelho). Isso nos indica diferencas na expressao dos dois grupos e que provavelmente a exposicao ao Ni no trabalho tem influencia nesse fenomeno.

## K-means
```{r}

kmeans = kmeans(t(dat.norm), centers=2,iter.max=200)
dat.anno$tipo = c(rep(1, 10),rep(2, 8))

## Check if clustering is equal to real groups
kmeans$cluster == dat.anno$tipo
