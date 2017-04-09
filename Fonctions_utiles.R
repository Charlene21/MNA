###### Fonctions communes aux deux parties #####

#Fonction h
h <- function(x){
  res = S0*exp(x) - K;
  if (res < 0){
    res = 0;
  }
  
  return(res);
}

#Fonction f
f <- function(x, t){
  return(h(x+r*t))
}


#Calcul des points de discrétisation du modèle
calcul_points_discretisation <- function(extremite_droite, extremite_gauche, delta_x){
  nb_points = (extremite_droite - extremite_gauche)/delta_x;
  vect_discr = c(extremite_gauche);
  x = extremite_gauche;
  
  for (i in 1 : nb_points){
    x = x + delta_x;
    vect_discr = c(vect_discr, x);
  }
  
  return(vect_discr);
}

#Matrice tridiagonale
create_mat_tridiag <- function(delta_x, points_discr){
  
  delta_t = delta_x * delta_x;
  j = 1 + (sigma*sigma*delta_t)/(delta_x*delta_x);
  j_1 = -(sigma*sigma*delta_t)/(2*delta_x*delta_x);
  matrix <- diag(j,length(points_discr), length(points_discr));
  for (i in 1:length(points_discr)-1){
    matrix[i,i+1] = j_1;
  }
  
  for (i in 2:length(points_discr)){
    matrix[i,i-1] = j_1;
  }
  
  return (matrix);
}



#Pour resoudre le systeme matriciel AX = B
resol_syst <- function(X_1, A, B){
  X = solve(A, B %*% X_1);
  return(X);
  
}

#Condition aux limites pour t = 0
calculate_uo <- function(points_discr){
  uo = cbind();
  for (i in 1 : length(points_discr)){
    uo = c(uo, h(points_discr[i]));
  }
  
  return (uo);
}

