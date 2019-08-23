# rvencio ago 2019

	D = function(A,B){  return(as.numeric(sqrt(sum((A-B)^2)))) }

	KNN = function(x, X, L, k=3){
		n = nrow(X)
		d = 0*(1:n) # inicializa vetor de zeros em n posições
		for(i in 1:n){ d[i] = D(x, X[i,]) }
		o = order(d)
		vote = sort(table(L[o[1:k]]), decreasing = TRUE)
		return(as.numeric(names(vote[1])))
	}

	x = read.delim("SkillCraft1_Dataset.csv", sep=",")
	classe = x[,2]
	IDs = x[,1]
	n = ncol(x)
	y = x[, 6:n] # queremos só de APM em diante, sem idade nem nada
	y = scale(y) 
	N = nrow(y)

	m = length(unique(classe))
	confu = matrix(0, m, m)

date()
	for(u in 1:N){
		predito = KNN(y[u,] , y[-u,] , classe[-u], k=1)
		for(i in 1:m){
		  for(j in 1:m){
		   confu[i,j] = confu[i,j] + sum((predito == i) & (classe[u] == j))
		  }
		}
	}
date()




























