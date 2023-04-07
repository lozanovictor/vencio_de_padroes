# Plataform

# *Affymetrix Human Gene 1.0 ST Array*

## Disciplina 5955016 - Reconhecimento de Padrões
    
## Introdução
    
#*A exposição ocupacional ao níquel (Ni) esta associada a um aumento do risco de cânceres pulmonares e nasais. 
#Os compostos de Ni exibem fraca atividade mutagênica, alteram a homeostase epigenética da célula e ativam as vias de sinalização. 
#No entanto, alterações na expressão gênica associadas à exposição ao Ni foram investigadas apenas in vitro. 
#Este estudo foi conduzido em uma população chinesa para determinar se a exposição ocupacional ao Ni foi associada com perfis diferenciais de 
#expressão gênica nas células mononucleares do sangue periférico (PBMCs) dos trabalhadores do refino Ni, 
#quando comparados à referência.*
    
    
# Métodos
    
#*Oito trabalhadores da refinaria e dez "normais" foram selecionados. 
#O RNA de PBMC foi extra?do e o perfil de expressão gênica foi realizado usando matrizes de exon Affymetrix. 
#Genes diferencialmente expressos entre os dois grupos foram identificados em uma análise global. 
#A significância das mudanças de expressão gênica entre as populações exposta e normal foi 
#avaliada usando uma abordagem de modelo linear com LIMMA, que utiliza uma abordagem empirica de Bayes para gerar estatísticas 
#t moderadas, levando em consideração os erros padrão e estimativas dobrar as alterações. 
#Os genes que atingiram os limiares de significância foram incluídos no conjunto de estudos para análise posterior e 
#são referidos como genes diferencialmente expressos.*
    
    
# Plataform
    
#*Affymetrix Human Gene 1.0 ST Array*
    
    

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


dim(dat) #dimensões dos dados
View(dat.anno) #ve os metadados
colnames(dat) == rownames(dat.anno) #checa ordenacao, pois em dat.anno as 10 primeiras amostras são de oene e as 8 ultimas de expostos ao Ni
table(is.na(dat)) #checa se tem NA, pois pode atrapalhar nas analises
head(dat)


# Busca Outliers
dat.mean <- apply(dat,1,mean) # Calcula a media
dat.sd <- sqrt(apply(dat,1,var)) # Calcula o desvio padrao
dat.cv <- dat.sd/dat.mean # Calcula o coeficiente de variacao
plot(dat.mean,dat.cv,main="Sample CV vs. Mean",xlab="Mean",ylab="CV",cex=1) #Plot do grafico

# Separa em grupos
dat.anno[1:10, "color"] <- "red" #atribui uma cor ao primeiro grupo de amostras - environmental exposure (EE)
dat.anno[11:18, "color"] <- "green" #atribui cor ao segundo grupo - occupational exposure nickel (OEN)
dat.anno[colnames(dat), "color"] #checa o atribuicao



## Normalizacao dos dados

limma::plotMA(dat) # distribuicao original
abline(h=0, col="dark grey")

# Normalizacao dos dados pelo método loess
dat.norm<-normalizeCyclicLoess(dat)

limma::plotMA(dat.norm) #Plot apos normalizacao
abline(h=0, col="dark grey")


# Podemos ver o conjunto de pontos bem mais homogeneo e alinhado a linha central do grafico, 
#mostrando a importancia da normalizacao


## Dendograma

colnames(dat.norm) = dat.anno$geo_accession # altera o nome das amostras
plot(hclust(dist(t(dat.norm), method = "euclidean")),main = "Hclust Dendogram"); #plota o cluster no metodo euclidiano

# Fica claro que 7 das 8 amostras de pessoas expostas ao Ni são agrupadas
#na mesma categoria e as demais amostras agrupadas separadamente 
#(amostra GSM992710, exposta a Ni no trabalho,
#é agrupada com as amostras de referencia).

## Correlacao

# busca por correlacao
dat.cor <- cor(dat.norm) # calcula a correlacao
image (dat.cor, axes=F, main = "Correlation matrix") # gera a imagem da matriz

# Existem poucos genes correlacionados no plot, a maioria relacionando o grupo ee vs ee (controle),
#mas tambem alguns relacionando os dois grupos. A maioria dos dados tem baixa correlacao.

## Plot Profile

## PCA
```{r}

data.pca.n <- prcomp(t(dat.norm))
data.pca.n <- data.pca.n$x[,1:2]
plot(data.pca.n[,1],data.pca.n[,2], col=as.factor(dat.anno$color), main="PCA",xlab="PC1", ylab="PC2")
legend(1,1,c("O. exposure - black","E. exposure - red"),lty=c(1,1), col=as.factor(dat.anno$color),cex=1, bty="n")

```

# Podemos ver claramente a tendencia de agrupamento entre as amostras similares, 
#separando o grupo oen, acumulado na parte superior direita do plot (em preto), 
#e o grupo ee, acumulado na parte inferior esquerda (vermelho). 
#Isso nos indica diferencas na expressao dos dois grupos e que provavelmente a exposicao ao Ni no trabalho tem influencia nesse fenomeno.

## K-means
```{r}

kmeans = kmeans(t(dat.norm), centers=2,iter.max=200)
dat.anno$tipo = c(rep(1, 10),rep(2, 8))

## Check if clustering is equal to real groups
kmeans$cluster == dat.anno$tipo
