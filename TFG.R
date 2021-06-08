
# Grafo corto

activities <- data.frame(
  id =       c( 0, 1, 2, 3, 4, 5, 6 ),
  duration = c( 0, 2, 7, 3, 4, 5, 0 ),
  resource = c( 0, 1, 2, 3, 3, 2, 0 )
)
relations <- data.frame(
  from = c(0, 0, 0, 1, 2, 3, 4, 5),
  to =   c(1, 2, 3, 4, 4, 5, 6, 6)
  )


# Grafo 22 actividades

activities <- data.frame(
  id =       c( 0:23 ),
  duration = c( 0, 4, 2, 4, 6, 1, 2, 3, 2, 4, 10, 3, 1, 2, 3, 2, 1, 1, 2, 3, 1, 2, 5, 0 ),
  resource = c( 0, 3.373, 4.842, 2.072, 2.300, 2.480, 1.868, 4.336, 3.618, 1.168, 2.351, 4.528, 3.662, 1.238, 1.825, 2.006, 3.674, 1.504, 3.401, 3.279, 4.286, 2.863, 3.841, 0 )
)
relations <- data.frame(
  from = c(0, 1, 2, 2, 2, 17, 21, 22, 3, 3, 8, 5, 6, 9, 5, 7, 10, 11, 11, 11, 14, 18, 12, 13, 19, 20, 3, 4, 15, 16 ),
  to =   c(1, 2, 3, 5, 17, 21, 22, 23, 8, 9, 10, 6, 9, 10, 7, 10, 11, 14, 12, 13, 18, 23, 19, 19, 20, 23, 4, 15, 16, 21 )
)


# Caso TFG

activities <- data.frame(
  id =       c( 0:14 ),
  duration = c( 0,3,3,1,3,3,5,5,1,2,2,2,1,2,0 ),
  resource = c( 0,1,3,1,3,3,5,4,1,3,3,2,1,2,0 )
)
relations <- data.frame(
  from = c(0, 0, 1, 1, 1, 2, 5, 6, 3, 4, 8, 8, 7, 7, 10, 11, 9, 12, 13),
  to =   c(1, 8, 3, 2, 8, 5, 6, 7, 4, 7, 9, 12, 10, 9, 11, 13, 13, 13, 14)
)


# Las simulaciones están al final del código


solution = function(FO, graph){
  
library(criticalpath)

schedule <- Schedule$new(activities, relations)


ES = schedule$activities$ES
LS = schedule$activities$LS
TF = schedule$activities$total_float
activities = cbind(activities, ES, LS, TF)



# Matriz auxiliar X

X <- matrix(c(rep(0, length(activities$id)*(schedule$duration + 1)*4)), ncol=4, byrow=T)

colnames(X) = c("Activity i", "Time t", "Resource", "isX")

X[,1] = c(rep(activities$id, each=schedule$duration+1))
X[,2] = c(rep(0:(schedule$duration),length(activities$id)))


for(i in 1:length(activities$id)){
  for(x in 1:nrow(X)){
   if(activities$id[i]==X[x,1]){
     X[x,3]= activities$resource[i] 
    }
  }
}


library(dplyr)

for(i in 1:length(activities$id)){
  
  for(act in 1:nrow(X)){
    
    if(activities$id[i]==X[act,1]){
    if(X[act,2] %in% activities$ES[i]:activities$LS[i]){
      
    X[act,4]=1
    
      }
    } 
  }
}


X = filter(as.data.frame(X), isX == 1)
X = X[,c(1,2,3)]



n_var_X = sum(schedule$activities$total_float + 1)
n_var_Z = schedule$duration
n_var = n_var_X + n_var_Z



library(gurobi)

model <- list()


# Objetivo según cada criterio de nivelación. 
## La numeración de las FO está según Damci y Polat: 
## (Criterio - FO): 3.2.1 - 7, 3.2.2 - 8, 3.2.3 - 1, 3.2.4 - 2, 3.2.5 - 9, 3.2.6 - 3, 3.2.7 - 4, 3.2.8 - 5, 3.2.9 - 6.  

# Matriz Q y obj para función 7
# Minimization of the sum of the square of daily resource usage

if(FO=='7'){
  
model$obj = c(rep(0,n_var))
model$Q = matrix(c(rep(0,n_var*n_var)), ncol = n_var, nrow = n_var, byrow = T)

for(v in 1:schedule$duration){

  model$Q[v+n_var_X,v+n_var_X] = 1

  }
}


# Matriz Q y obj para función 8
# Minimization of the sum of the square of the deviations in daily resource usage 

if(FO=='8'){
  
model$obj = c(rep(0,n_var))
model$Q = matrix(c(rep(0,n_var*n_var)), ncol = n_var, nrow = n_var, byrow = T)

for(v in 1:schedule$duration){
  if(v==1){
    model$Q[v+n_var_X,v+n_var_X] = 1 
  } 
  else{
    if(v==schedule$duration){
      model$Q[v+n_var_X,v+n_var_X] = 1
    } else{
      model$Q[v+n_var_X,v+n_var_X] = 2
    }
  }
}


for(v in 1:(schedule$duration-1)){
  
  model$Q[v+n_var_X,v+n_var_X+1] = -1
  
}

for(v in 2:schedule$duration){
  
  model$Q[v+n_var_X,v+n_var_X-1] = -1
  
}
}



# Matriz Q y obj para función 9
# Minimization of the sum of the square of the deviations between daily resource usage 
# and the average resource usage  

if(FO=='9'){
  
model$obj = c(rep(0,n_var_X), rep(-2*(sum(activities$resource*activities$duration)/schedule$duration),n_var_Z))
model$Q = matrix(c(rep(0,n_var*n_var)), ncol = n_var, nrow = n_var, byrow = T)
for(v in 1:schedule$duration){
  
  model$Q[v+n_var_X,v+n_var_X] = 1
  
}
}


# obj para función 1
# Minimization of the sum of the absolute deviations in daily resource usage  

if(FO=='1'){
  
model$obj = c(rep(0,n_var), rep(1,n_var_Z-1))

R00 = matrix(0, nrow = (n_var_Z-1), ncol= n_var + n_var_Z - 1)

for(i in 1:(schedule$duration-1)){
  for(j in (1:schedule$duration-1)){
    if(i==j){
      R00[i,n_var+j]=1
      
    }
  }
}



R01 = matrix(0, nrow = (n_var_Z-1), ncol= n_var + n_var_Z -1)

for(i in 1:(schedule$duration-1)){
  for(j in 1:(schedule$duration-1)){
    if(i==j){
      R01[i,n_var+j]=1
      
    }
  }
}

for(i in 1:(schedule$duration-1)){
  for(j in 1:(schedule$duration-1)){
    if(i==j){
      R01[i,n_var_X+j]=-1
      
    }
  }
}
for(i in 1:(schedule$duration-1)){
  for(j in 1:(schedule$duration-1)){
    if(i==j){
      R01[i,n_var_X+j+1]=1
      
    }
  }
}



R02 = matrix(0, nrow = (n_var_Z -1), ncol= n_var + n_var_Z -1)

for(i in 1:(schedule$duration-1)){
  for(j in 1:(schedule$duration-1)){
    if(i==j){
      R02[i,n_var+j]=1
      
    }
  }
}

for(i in 1:(schedule$duration-1)){
  for(j in 1:(schedule$duration-1)){
    if(i==j){
      R02[i,n_var_X+j]=1
      
    }
  }
}
for(i in 1:(schedule$duration-1)){
  for(j in 1:(schedule$duration-1)){
    if(i==j){
      R02[i,n_var_X+j+1]=-1
      
    }
  }
}

R0 = rbind(R00, R01, R02)
}




# obj para función 2
# Minimization of the sum of only the increases in daily resource usage  


if(FO=='2'){
  
  model$obj = c(rep(0,n_var), rep(1,n_var_Z-1))

  R00 = matrix(0, nrow = (n_var_Z-1), ncol= n_var + n_var_Z - 1)
  
  for(i in 1:(schedule$duration-1)){
    for(j in (1:schedule$duration-1)){
      if(i==j){
        R00[i,n_var+j]=1
        
      }
    }
  }
  
  R02 = matrix(0, nrow = (n_var_Z -1), ncol= n_var + n_var_Z -1)
  
  for(i in 1:(schedule$duration-1)){
    for(j in 1:(schedule$duration-1)){
      if(i==j){
        R02[i,n_var+j]=1
        
      }
    }
  }
  
  for(i in 1:(schedule$duration-1)){
    for(j in 1:(schedule$duration-1)){
      if(i==j){
        R02[i,n_var_X+j]=1
        
      }
    }
  }
  for(i in 1:(schedule$duration-1)){
    for(j in 1:(schedule$duration-1)){
      if(i==j){
        R02[i,n_var_X+j+1]=-1
        
      }
    }
  }
  
  R0 = rbind(R00,R02)
}
  
  
# obj para función 3
# Minimization of the sum of the the absolute deviations between daily resource usage
# and the average resource usage  

if(FO=='3'){
  
model$obj = c(rep(0,n_var), rep(1,n_var_Z)) 


R00 = matrix(0, nrow = (n_var_Z), ncol= n_var + n_var_Z)

for(i in 1:schedule$duration){
  for(j in 1:schedule$duration){
    if(i==j){
      R00[i,n_var+j]=1
      
    }
  }
}

R01 = matrix(0, nrow = (n_var_Z), ncol= n_var + n_var_Z)

for(i in 1:schedule$duration){
  for(j in 1:schedule$duration){
    if(i==j){
      R01[i,n_var+j]=1
      
    }
  }
}

for(i in 1:schedule$duration){
  for(j in 1:schedule$duration){
    if(i==j){
      R01[i,n_var_X+j]=-1
      
    }
  }
}


R02 = matrix(0, nrow = (n_var_Z), ncol= n_var + n_var_Z)

for(i in 1:schedule$duration){
  for(j in 1:schedule$duration){
    if(i==j){
      R02[i,n_var+j]=1
      
    }
  }
}

for(i in 1:schedule$duration){
  for(j in 1:schedule$duration){
    if(i==j){
      R02[i,n_var_X+j]=1
      
    }
  }
}

R0 = rbind(R00, R01, R02)

}



# obj para función 4
# Minimization of the maximum daily resource usage

if(FO=='4'){
  
model$obj = c(rep(0,n_var), 1)

R0 = matrix(0, ncol=n_var+1, nrow = n_var_Z)

for(i in 1:schedule$duration){
  R0[i,n_var+1]=1
  R0[i, n_var_X+i]=-1
}
}


# obj para función 5
# Minimization of the maximum deviation in daily resource usage

if(FO=='5'){
  
model$obj = c(rep(0,n_var), 1)

R00 = matrix(0, ncol=n_var+1, nrow = n_var_Z-1)

for(i in 1:(schedule$duration-1)){
    R00[i,n_var+1]=1
    R00[i, n_var_X+i]=-1
    R00[i, n_var_X+i+1]=1
}

R01 = matrix(0, ncol=n_var+1, nrow = n_var_Z-1)

for(i in 1:(schedule$duration-1)){
  R01[i,n_var+1]=1
  R01[i, n_var_X+i]=1
  R01[i, n_var_X+i+1]=-1
}


R02=c(rep(0,n_var),1)

R0 = rbind(R00, R01, R02)

}



# obj para función 6
# Minimization of the maximum  absolute deviation between daily resource usage
# and the average resource usage

if(FO=='6'){
  
model$obj = c(rep(0,n_var), 1)

R00 = matrix(0, ncol=n_var+1, nrow = n_var_Z)

for(i in 1:schedule$duration){
  R00[i,n_var+1]=1
  R00[i, n_var_X+i]=-1
}

R01 = matrix(0, ncol=n_var+1, nrow = n_var_Z)

for(i in 1:schedule$duration){
  R01[i,n_var+1]=1
  R01[i, n_var_X+i]=1
}


R02=c(rep(0,n_var),1)

R0 = rbind(R00, R01, R02)
}






# Restricciones 1 y 2

R12 = matrix(0, nrow=length(activities$id), ncol = n_var)
R12[1,1]=1
for(i in 2:length(activities$id)){
  R12[i,1+(sum(activities$TF[1:i-1]+1)):((sum(activities$TF[1:i-1]+1))+activities$TF[i])]=1
}


LD12 = c(rep(1,length(activities$id)))

sense12 = c(rep("=",length(activities$id)))



# Restricciones 3 

R3 = matrix(0, ncol=n_var, nrow = length(relations$from))

for(i in 1:length(relations$from)){
    for(x in 1:nrow(X)){
      
        if(X[x,1]==relations$from[i]){
          R3[i,x]=-X[x,2]
          
        }
       else{
        if(X[x,1]==relations$to[i]){
          R3[i,x]=X[x,2]}
        else{R3[i,x]=0}
    }  
  }
}

LD3 = c(rep(0,length(relations$from)))
for(i in 1:length(relations$from)){
  for(j in 1:length(relations$from)){
    if(relations$from[i]==j-1){
      LD3[i]=activities$duration[j]
    }
  }
  
}

sense3 = c(rep(">=",length(relations$from)))




# Restricciones 4

activities$ini[1] = 1
activities$fin[1] = 1

for(i in 2:nrow(activities)){
  activities$ini[i] = activities$fin[i-1] + 1
  activities$fin[i] = activities$fin[i-1] + 1 + activities$TF[i]
}


time.ind = data.frame(ind = rep(0,max(activities$LS)))
time.ind$ind[1] = activities$fin[nrow(activities)] + 1
for(i in 2:nrow(time.ind)){
  time.ind$ind[i] = time.ind$ind[i-1] + 1
}


R4 = matrix(0, ncol=n_var, nrow = nrow(time.ind))

diag(R4[,time.ind$ind[1]:time.ind$ind[nrow(time.ind)]]) = rep(1, nrow(time.ind)) 


for(t in 0:(nrow(time.ind)-1)){
  for(i in 1:nrow(activities)){
    tauini = max(activities$ES[i],t-activities$duration[i]+1) 
    taufin = min(t, activities$LS[i]) 
    if(tauini<=taufin){ 
      for(j in tauini:taufin){
        R4[t+1,j-activities$ES[i]+activities$ini[i]] = -activities$resource[i] 
      }
    }
  }
}

sense4 = c(rep("=",schedule$duration))
LD4 = c(rep(0,schedule$duration))



LD = c(LD12, LD3, LD4)
sense = c(sense12, sense3, sense4)

model$rhs        <- LD
model$sense      <- sense

model$vtype <- c(c(rep("B",n_var_X)), c(rep("C", n_var_Z)))

model$A <- rbind(R12, R3, R4)

model$modelsense <- "min"




if(FO=="1"){
  A2 = matrix(0, nrow = nrow(model$A), ncol = n_var_Z-1)
  model$A = cbind(model$A, A2)
  model$A = rbind(model$A, R0)
  
  model$rhs = c(model$rhs, c(rep(0, (n_var_Z-1)*3)))
  model$sense = c(model$sense, c(rep(">=", (n_var_Z-1)*3)))
  model$vtype = c(model$vtype, c(rep("C", n_var_Z-1)))
}

if(FO=="2"){
  A2 = matrix(0, nrow = nrow(model$A), ncol = n_var_Z-1)
  model$A = cbind(model$A, A2)
  model$A = rbind(model$A, R0)
  
  model$rhs = c(model$rhs, c(rep(0, (n_var_Z-1)*2)))
  model$sense = c(model$sense, c(rep(">=", (n_var_Z-1)*2)))
  model$vtype = c(model$vtype, c(rep("C", n_var_Z-1)))
}

if(FO=="3"){
  A2 = matrix(0, nrow = nrow(model$A), ncol = n_var_Z)
  model$A = cbind(model$A, A2)
  model$A = rbind(model$A, R0)
  
  model$rhs = c(model$rhs, c(rep(0, n_var_Z)), c(rep(-(sum(activities$resource*activities$duration)/schedule$duration), n_var_Z)), c(rep((sum(activities$resource*activities$duration)/schedule$duration), n_var_Z)))
  model$sense = c(model$sense, c(rep(">=", n_var_Z*3)))
  model$vtype = c(model$vtype, c(rep("C", n_var_Z)))
}

if(FO=="4"){
  model$A = cbind(model$A, c(rep(0,nrow(model$A))))
  model$A = rbind(model$A, R0)
  
  model$rhs = c(model$rhs, c(rep(0, n_var_Z)))
  model$sense = c(model$sense, c(rep(">=", n_var_Z)))
  model$vtype = c(model$vtype, "C")
}

if(FO=="5"){
  model$A = cbind(model$A, c(rep(0,nrow(model$A))))
  model$A = rbind(model$A, R0)
  
  model$rhs = c(model$rhs, c(rep(0, ((n_var_Z-1)*2)+1)))
  model$sense = c(model$sense, c(rep(">=", ((n_var_Z-1)*2)+1)))
  model$vtype = c(model$vtype, "C")
}

if(FO=="6"){
  model$A = cbind(model$A, c(rep(0,nrow(model$A))))
  model$A = rbind(model$A, R0)
  
  model$rhs = c(model$rhs, c(rep(-(sum(activities$resource*activities$duration)/schedule$duration),n_var_Z)), c(rep((sum(activities$resource*activities$duration)/schedule$duration),n_var_Z)), 0)
  model$sense = c(model$sense, c(rep(">=", ((n_var_Z)*2)+1)))
  model$vtype = c(model$vtype, "C")
}



params <- list()
params$TimeLimit <- 60*60  #segundos
result <- gurobi(model, params)

if(FO=="9"){
  result$objval = result$objval + n_var_Z*((sum(activities$resource*activities$duration)/schedule$duration))^2
}

# Solución
print(result$objval) # valor objetivo
print(result$x) # resultado variables
print(result$status)   # indica si se ha encontrado una solución óptima

Plan = cbind(X[,-3],isS=result$x[1:n_var_X])
Plan = filter(as.data.frame(Plan), isS == 1)
Plan = cbind(Plan, S = Plan[,2]*Plan[,3])
Plan = Plan[,c(1,4)]
#View(Plan)



## Visualizaciones

library(ggplot2)

recursosES = c(rep(0, schedule$duration))
for(i in 1:(length(activities$ES)-1)){
  for(j in (activities$ES[i]):(activities$ES[i]+activities$duration[i]-1)){
    recursosES[j+1] <- recursosES[j+1] + activities$resource[i] 
  }
}

recursosLS = c(rep(0, schedule$duration))

for(i in 1:(length(activities$LS)-1)){
  for(j in (activities$LS[i]):(activities$LS[i]+activities$duration[i]-1)){
    recursosLS[j+1] <- recursosLS[j+1] + activities$resource[i] 
  }
}



# Histograma consumo a lo largo del tiempo

if(graph=="hist"){
  
x = seq(0,schedule$duration-1) + 0.5  
y = result$x[(n_var_X+1): n_var]
df = data.frame(cbind(x,y))

print(ggplot(df) + 
  geom_col(aes(x = x, y = y),  width = 1) + #el width = 1 es para que pinte las barras juntas
  scale_x_continuous(breaks=seq(0,schedule$duration)) +
  labs(#title = "Consumo de recursos en cada instante" , 
       #subtitle = "¿Cómo avanza el requerimiento de recursos a lo largo del tiempo? ",       x = "t", 
       y = "Zt" ) +
  ylim(c(0,max(result$x)+3)) +
  geom_hline(yintercept=(sum(activities$resource*activities$duration)/schedule$duration), color = "red")) 
}



# Comparativa de consumo con planificaciones ES y LS

if(graph=="comp"){
  
datosOpt= data.frame(Instante=seq(0,schedule$duration-1) + 0.5, Recursos=result$x[(n_var_X+1): n_var], Planificación =c(rep("S(Opt 8)",schedule$duration)))
  
datosES = data.frame(Instante=seq(0,schedule$duration-1) + 0.5 , Recursos=recursosES, Planificación =c(rep("S(ES)",schedule$duration)))

datosLS = data.frame(Instante=seq(0,schedule$duration-1) + 0.5, Recursos=recursosLS, Planificación =c(rep("S(LS)",schedule$duration)))
  
datos=rbind(datosOpt, datosES, datosLS)
  
print(ggplot(data = datos, aes(x = Instante, y = Recursos, colour = Planificación)) + 
    labs(#title = "Consumo de recursos", 
         #subtitle = "¿Cómo avanza el requerimiento de recursos a lo largo del tiempo? ",
         x = "t", 
         y = "Zt" ) +
    ylim(c(0, 9
           #max(result$x)+3  #modificar en función de las dimensiones del problema (9 es para el Caso TFG)
           )) +
    geom_line() +
    geom_point() +
    scale_color_manual(values=c('#8BC540', '#DC5D42', "#25AAE2")) +
    scale_x_continuous(breaks=seq(0,schedule$duration))) 

  }   
#}




# Tabla valores objetivo 

if(FO=='ES'){
  z = recursosES
} else{
  if(FO=='LS'){
    z = recursosLS
  }else{
    z = result$x[(n_var_X+1):(n_var_X+n_var_Z)]
  }
}


# 3.2.1: 7
z7 = c(rep(0,schedule$duration))
for(i in 1:length(z)){
  z7[i] = z[i]^2
}
sum(z7)

# 3.2.2: 8
z8 = c(rep(0,schedule$duration-1))
for(i in 1:length(z8)){
  z8[i] = (z[i]-z[i+1])^2
}
sum(z8)

# 3.2.3: 1
z1 = c(rep(0,schedule$duration-1))
for(i in 1:length(z1)){
  z1[i] = abs(z[i]-z[i+1])
}
sum(z1)

# 3.2.4: 2
z2 = c(rep(0,schedule$duration-1))
for(i in 1:length(z2)){
  if(z[i+1]>z[i]){
    z2[i] = abs(z[i]-z[i+1])
  }else{
    z2[i] = 0
  }
}
sum(z2)

# 3.2.5: 9
z9 = c(rep(0,schedule$duration))
for(i in 1:length(z9)){
  z9[i] = (z[i] - sum(activities$resource*activities$duration)/schedule$duration)^2
}
sum(z9)

# 3.2.6: 3
z3 = c(rep(0,schedule$duration))
for(i in 1:length(z3)){
  z3[i] = abs(z[i] - sum(activities$resource*activities$duration)/schedule$duration)
}
sum(z3)

# 3.2.7: 4
z4 = max(z)

# 3.2.8: 5
z5 = c(rep(0,schedule$duration-1))
for(i in 1:length(z5)){
  z5[i] = abs(z[i]-z[i+1])
}
max(z5)


# 3.2.9: 6
z6 = c(rep(0,schedule$duration))
for(i in 1:length(z6)){
  z6[i] = abs(z[i] - sum(activities$resource*activities$duration)/schedule$duration)
}
max(z6)

obj = t(t(round(c(sum(z7), sum(z8), sum(z1), sum(z2), sum(z9), sum(z3), z4, max(z5), max(z6)),2)))
#View(obj)

}




solution("ES","hist")
solution("LS","hist")
solution(7,"hist")
solution(8,"hist")
solution(1,"hist")
solution(2,"hist")
solution(9,"hist")
solution(3,"hist")
solution(4,"hist")
solution(5,"hist")
solution(6,"hist")

solution(5,"comp")





# Generar grafos


grafo <- function(n,prob){
library(pcalg)
set.seed(0)
mydag <- randomDAG(n, prob, lB = 0.1, uB = 1, V = as.character(1:n))
plot(mydag)
mydag
} 
  
grafo(50,0.1)



## Grafo n=15 y prob = 0.2

set.seed(15)
n=15

activities <- data.frame(
  id =       c( 0:(n+1) ),
  duration = c( 0, sample(1:5, n , replace = TRUE), 0 ),
  resource = c( 0, runif(n, 0, 5), 0 )
)

# A los que no le llegan ninguno: from 0
# De los que no sale ninguno: to n+1
relations <- data.frame(
  from = c(0,1,1,1,1,1,2,3,3,0,4,4,5,0,6,0,7,8,0,9,0,10,10,11,12,13,14,15),
  to =   c(1,2,8,3,10,5,5,8,13,4,11,13,14,6,16,7,12,14,9,16,10,14,11,13,14,14,15,16)
)



## Grafo n=30, prob = 0.2, duración y recursos 1:5

set.seed(3002)
n=30

activities <- data.frame(
  id =       c( 0:(n+1) ),
  duration = c( 0, sample(1:5, n , replace = TRUE), 0 ),
  resource = c( 0, runif(n, 0, 5), 0 )
)

relations <- data.frame(
  from = c(0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,0,4,4,4,4,4,4,5,5,5,5,5,6,6,6,6,7,7,8,8,8,8,8,8,
           9,9,10,10,10,10,10,10,10,10,11,11,11,11,11,12,12,15,15,13,13,13,13,13,14,14,14,14,16,16,16,16,17,18,19,20,20,20,20,21,22,23,24,25,26,27,28,29,30),
  to =   c(1,26,19,2,5,3,15,12,8,24,7,23,30,16,17,27,11,28,13,6,27,9,23,4,27,12,11,10,16,29,23,19,18,27,11,13,18,12,14,23,21,14,21,9,20,30,24,
           20,25,27,11,12,25,17,21,20,23,19,12,14,22,29,17,21,19,24,28,21,30,19,16,16,15,20,22,19,18,17,21,30,20,21,23,26,30,29,23,30,26,29,30,28,31,31,31,31)
)



## Grafo n=30, prob = 0.2, duración y recursos 1:10

set.seed(3002)
n=30

activities <- data.frame(
  id =       c( 0:(n+1) ),
  duration = c( 0, sample(1:10, n , replace = TRUE), 0 ),
  resource = c( 0, runif(n, 0, 10), 0 )
)

relations <- data.frame(
  from = c(0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,0,4,4,4,4,4,4,5,5,5,5,5,6,6,6,6,7,7,8,8,8,8,8,8,
           9,9,10,10,10,10,10,10,10,10,11,11,11,11,11,12,12,15,15,13,13,13,13,13,14,14,14,14,16,16,16,16,17,18,19,20,20,20,20,21,22,23,24,25,26,27,28,29,30),
  to =   c(1,26,19,2,5,3,15,12,8,24,7,23,30,16,17,27,11,28,13,6,27,9,23,4,27,12,11,10,16,29,23,19,18,27,11,13,18,12,14,23,21,14,21,9,20,30,24,
           20,25,27,11,12,25,17,21,20,23,19,12,14,22,29,17,21,19,24,28,21,30,19,16,16,15,20,22,19,18,17,21,30,20,21,23,26,30,29,23,30,26,29,30,28,31,31,31,31)
)



## Grafo n=30, prob = 0.3, duración y recursos 1:5

set.seed(3003)
n=30
activities <- data.frame(
  id =       c( 0:(n+1) ),
  duration = c( 0, sample(1:5, n , replace = TRUE), 0 ),
  resource = c( 0, runif(n, 0, 5), 0 )
)

relations <- data.frame(
  from = c(0,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,5,5,5,5,5,0,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,10,10,10,10,11,11,11,12,12,12,12,
           14,14,14,14,15,15,15,15,15,16,16,16,16,16,16,16,17,17,17,17,18,18,18,19,19,20,21,21,21,21,21,22,23,24,25,26,26,27,27,28,29,30),
  to =   c(1,19,27,12,8,20,5,11,24,3,26,2,27,28,21,14,5,4,10,28,3,22,9,13,10,27,15,26,22,28,12,30,5,23,16,10,28,19,21,30,25,28,6,12,19,15,23,13,7,18,22,8,14,24,30,17,9,26,18,29,17,29,15,21,11,16,9,12,23,21,16,25,26,28,29,17,30,25,26,20,23,18,15,30,14,13,
           19,30,21,23,29,19,27,21,17,20,21,24,29,25,17,28,27,20,25,26,20,26,22,30,27,21,27,29,30,24,22,25,24,30,29,28,30,29,30,31,31,31)
)



## Grafo n=30, prob = 0.3, duración y recursos 1:10

set.seed(3003)
n=30
activities <- data.frame(
  id =       c( 0:(n+1) ),
  duration = c( 0, sample(1:10, n , replace = TRUE), 0 ),
  resource = c( 0, runif(n, 0, 10), 0 )
)

relations <- data.frame(
  from = c(0,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,5,5,5,5,5,0,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,10,10,10,10,11,11,11,12,12,12,12,
           14,14,14,14,15,15,15,15,15,16,16,16,16,16,16,16,17,17,17,17,18,18,18,19,19,20,21,21,21,21,21,22,23,24,25,26,26,27,27,28,29,30),
  to =   c(1,19,27,12,8,20,5,11,24,3,26,2,27,28,21,14,5,4,10,28,3,22,9,13,10,27,15,26,22,28,12,30,5,23,16,10,28,19,21,30,25,28,6,12,19,15,23,13,7,18,22,8,14,24,30,17,9,26,18,29,17,29,15,21,11,16,9,12,23,21,16,25,26,28,29,17,30,25,26,20,23,18,15,30,14,13,
           19,30,21,23,29,19,27,21,17,20,21,24,29,25,17,28,27,20,25,26,20,26,22,30,27,21,27,29,30,24,22,25,24,30,29,28,30,29,30,31,31,31)
)



## Grafo n=40, prob = 0.1, duración y recursos 1:5

set.seed(40)
n=40

activities <- data.frame(
  id =       c( 0:(n+1) ),
  duration = c( 0, sample(1:5, n , replace = TRUE), 0 ),
  resource = c( 0, runif(n, 0, 5), 0 )
)

relations <- data.frame(
  from = c(0,1,1,1,1,1,1,2,2,2,2,2,0,3,3,3,3,3,3,0,4,4,4,4,9,5,5,5,6,0,7,7,7,7,7,8,8,8,10,10,10,0,14,15,15,0,16,17,18,19,19,0,20,20,21,21,22,23,24,25,26,26,27,28,29,29,32,0,33,34,35,36,37,38,39,40),
  to =   c(1,24,2,5,19,15,35,17,9,23,11,12,3,37,23,28,6,18,36,4,24,10,29,27,34,11,29,37,28,7,8,26,15,39,13,12,25,19,17,23,30,14,24,22,34,16,17,41,41,29,35,20,27,28,24,32,28,32,38,31,29,40,36,36,31,36,41,33,34,41,39,41,41,41,41,41)
)



## Grafo n=40, prob = 0.1, duración y recursos 1:10

set.seed(40)
n=40

activities <- data.frame(
  id =       c( 0:(n+1) ),
  duration = c( 0, sample(1:10, n , replace = TRUE), 0 ),
  resource = c( 0, runif(n, 0, 10), 0 )
)

relations <- data.frame(
  from = c(0,1,1,1,1,1,1,2,2,2,2,2,0,3,3,3,3,3,3,0,4,4,4,4,9,5,5,5,6,0,7,7,7,7,7,8,8,8,10,10,10,0,14,15,15,0,16,17,18,19,19,0,20,20,21,21,22,23,24,25,26,26,27,28,29,29,32,0,33,34,35,36,37,38,39,40),
  to =   c(1,24,2,5,19,15,35,17,9,23,11,12,3,37,23,28,6,18,36,4,24,10,29,27,34,11,29,37,28,7,8,26,15,39,13,12,25,19,17,23,30,14,24,22,34,16,17,41,41,29,35,20,27,28,24,32,28,32,38,31,29,40,36,36,31,36,41,33,34,41,39,41,41,41,41,41)
)


## Grafo n=40, prob = 0.145, duración y recursos 1:5

set.seed(402)
n=40

activities <- data.frame(
  id =       c( 0:(n+1) ),
  duration = c( 0, sample(1:5, n , replace = TRUE), 0 ),
  resource = c( 0, runif(n, 0, 5), 0 )
)

relations <- data.frame(
  from = c(0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,0,3,3,3,3,3,3,3,0,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,6,6,6,6,6,0,7,7,7,7,7,7,7,7,0,8,8,8,8,9,9,9,9,9,9,9,10,10,10,10,10,11,11,11,11,
           0,12,12,12,12,12,12,13,13,14,14,14,14,14,15,15,16,16,17,17,17,17,18,18,19,19,19,0,20,20,20,21,21,22,22,23,23,24,24,24,25,26,27,28,29,29,30,31,32,33,34,35,35,35,36,37,38,39,40),
  to =   c(1,15,35,2,22,19,5,24,18,17,35,36,39,40,27,3,26,28,40,31,23,9,11,4,28,36,17,6,22,18,26,10,22,18,16,30,11,9,13,28,34,35,27,22,7,26,23,16,14,30,33,29,25,8,9,10,18,37,19,16,30,33,25,21,38,11,13,21,32,24,14,28,33,25,
           12,36,21,30,24,29,34,38,33,40,16,18,33,30,28,17,27,18,31,28,36,30,33,28,26,40,24,20,35,36,26,27,30,31,27,31,27,40,25,29,34,31,29,41,31,37,41,41,41,41,41,36,39,40,41,41,41,41,41)
)



## Grafo n=40, prob = 0.145, duración y recursos 1:10

set.seed(402)
n=40

activities <- data.frame(
  id =       c( 0:(n+1) ),
  duration = c( 0, sample(1:10, n , replace = TRUE), 0 ),
  resource = c( 0, runif(n, 0, 10), 0 )
)

relations <- data.frame(
  from = c(0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,0,3,3,3,3,3,3,3,0,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,6,6,6,6,6,0,7,7,7,7,7,7,7,7,0,8,8,8,8,9,9,9,9,9,9,9,10,10,10,10,10,11,11,11,11,
           0,12,12,12,12,12,12,13,13,14,14,14,14,14,15,15,16,16,17,17,17,17,18,18,19,19,19,0,20,20,20,21,21,22,22,23,23,24,24,24,25,26,27,28,29,29,30,31,32,33,34,35,35,35,36,37,38,39,40),
  to =   c(1,15,35,2,22,19,5,24,18,17,35,36,39,40,27,3,26,28,40,31,23,9,11,4,28,36,17,6,22,18,26,10,22,18,16,30,11,9,13,28,34,35,27,22,7,26,23,16,14,30,33,29,25,8,9,10,18,37,19,16,30,33,25,21,38,11,13,21,32,24,14,28,33,25,
           12,36,21,30,24,29,34,38,33,40,16,18,33,30,28,17,27,18,31,28,36,30,33,28,26,40,24,20,35,36,26,27,30,31,27,31,27,40,25,29,34,31,29,41,31,37,41,41,41,41,41,36,39,40,41,41,41,41,41)
)



## Grafo n=50, prob = 0.06, duración y recursos 1:5

set.seed(501)
n=50

activities <- data.frame(
  id =       c( 0:(n+1) ),
  duration = c( 0, sample(1:5, n , replace = TRUE), 0 ),
  resource = c( 0, runif(n, 0, 5), 0 )
)

relations <- data.frame(
  from = c(0,1,1,1,1,1,2,2,2,2,0,3,0,4,5,5,0,6,6,0,7,0,8,8,8,8,0,9,9,0,10,10,0,11,12,13,13,13,13,14,14,0,15,16,0,17,17,17,17,0,18,18,18,18,18,18,19,0,20,21,21,0,22,
           0,33,33,34,35,36,36,37,23,23,24,24,25,25,0,26,27,27,28,29,30,31,32,0,38,38,39,39,40,0,41,42,43,44,45,46,47,48,49,50),
  to =   c(1,2,5,35,24,40,50,12,23,44,3,44,4,51,39,42,6,50,21,7,13,8,46,36,28,49,9,49,34,10,16,49,11,13,51,36,19,27,14,25,31,15,40,39,17,50,45,30,39,
           18,39,35,24,27,25,31,51,20,50,28,49,22,32,33,34,47,51,49,37,48,45,50,49,48,34,49,32,26,51,28,46,29,51,36,51,51,38,40,45,48,43,43,41,51,51,41,45,51,51,51,51,51,51)
)


## Grafo n=50, prob = 0.06, duración y recursos 1:10

set.seed(501)
n=50

activities <- data.frame(
  id =       c( 0:(n+1) ),
  duration = c( 0, sample(1:10, n , replace = TRUE), 0 ),
  resource = c( 0, runif(n, 0, 10), 0 )
)

relations <- data.frame(
  from = c(0,1,1,1,1,1,2,2,2,2,0,3,0,4,5,5,0,6,6,0,7,0,8,8,8,8,0,9,9,0,10,10,0,11,12,13,13,13,13,14,14,0,15,16,0,17,17,17,17,0,18,18,18,18,18,18,19,0,20,21,21,0,22,
           0,33,33,34,35,36,36,37,23,23,24,24,25,25,0,26,27,27,28,29,30,31,32,0,38,38,39,39,40,0,41,42,43,44,45,46,47,48,49,50),
  to =   c(1,2,5,35,24,40,50,12,23,44,3,44,4,51,39,42,6,50,21,7,13,8,46,36,28,49,9,49,34,10,16,49,11,13,51,36,19,27,14,25,31,15,40,39,17,50,45,30,39,
           18,39,35,24,27,25,31,51,20,50,28,49,22,32,33,34,47,51,49,37,48,45,50,49,48,34,49,32,26,51,28,46,29,51,36,51,51,38,40,45,48,43,43,41,51,51,41,45,51,51,51,51,51,51)
)



## Grafo n=50, prob = 0.1, duración y recursos 1:5

set.seed(50)
n=50

activities <- data.frame(
  id =       c( 0:(n+1) ),
  duration = c( 0, sample(1:5, n , replace = TRUE), 0 ),
  resource = c( 0, runif(n, 0, 5), 0 )
)

relations <- data.frame(
  from = c(0,1,1,1,1,1,1,1,1,2,2,2,2,0,3,3,3,3,0,4,4,5,5,5,5,5,5,0,6,6,6,6,0,7,7,7,7,7,7,7,8,8,8,9,9,9,9,10,10,10,10,10,10,10,
           11,11,11,11,12,12,12,12,13,13,13,14,14,15,15,15,0,16,16,16,16,16,16,17,17,18,18,18,18,19,19,19,19,19,19,20,21,21,21,21,22,23,23,24,24,24,25,25,26,27,27,
           28,28,29,30,31,32,33,34,35,35,35,36,37,37,38,38,39,39,39,40,40,40,41,42,43,44,45,46,47,48,49,50),
  to =   c(1,40,15,5,2,44,19,24,35,11,23,9,17,3,40,37,45,28,4,10,14,25,50,11,45,28,33,6,38,48,12,28,7,45,8,30,46,47,13,21,31,28,33,28,31,37,42,11,45,19,29,36,42,33,
           26,40,35,37,34,40,45,36,42,35,27,17,36,50,22,42,16,20,17,35,32,18,50,19,36,21,32,49,50,37,45,43,28,29,47,25,45,35,27,49,25,40,29,45,33,36,26,46,40,50,48,
           30,44,33,35,38,51,51,51,44,46,50,51,41,43,47,50,43,44,46,41,48,50,43,51,51,51,51,47,48,50,51,51)
)



## Grafo n=50, prob = 0.1, duración y recursos 1:10

set.seed(50)
n=50

activities <- data.frame(
  id =       c( 0:(n+1) ),
  duration = c( 0, sample(1:10, n , replace = TRUE), 0 ),
  resource = c( 0, runif(n, 0, 10), 0 )
)

relations <- data.frame(
  from = c(0,1,1,1,1,1,1,1,1,2,2,2,2,0,3,3,3,3,0,4,4,5,5,5,5,5,5,0,6,6,6,6,0,7,7,7,7,7,7,7,8,8,8,9,9,9,9,10,10,10,10,10,10,10,
           11,11,11,11,12,12,12,12,13,13,13,14,14,15,15,15,0,16,16,16,16,16,16,17,17,18,18,18,18,19,19,19,19,19,19,20,21,21,21,21,22,23,23,24,24,24,25,25,26,27,27,
           28,28,29,30,31,32,33,34,35,35,35,36,37,37,38,38,39,39,39,40,40,40,41,42,43,44,45,46,47,48,49,50),
  to =   c(1,40,15,5,2,44,19,24,35,11,23,9,17,3,40,37,45,28,4,10,14,25,50,11,45,28,33,6,38,48,12,28,7,45,8,30,46,47,13,21,31,28,33,28,31,37,42,11,45,19,29,36,42,33,
           26,40,35,37,34,40,45,36,42,35,27,17,36,50,22,42,16,20,17,35,32,18,50,19,36,21,32,49,50,37,45,43,28,29,47,25,45,35,27,49,25,40,29,45,33,36,26,46,40,50,48,
           30,44,33,35,38,51,51,51,44,46,50,51,41,43,47,50,43,44,46,41,48,50,43,51,51,51,51,47,48,50,51,51)
)



# Especificaciones duraciones y consumos (sin actividades ficticias)

mean(activities$duration[-c(1,nrow(activities))])
quantile(activities$duration[-c(1,nrow(activities))])

mean(activities$resource[-c(1,nrow(activities))])
quantile(activities$resource[-c(1,nrow(activities))])

