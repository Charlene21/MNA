S0 = 1;
K = 1;
T=1;
r = 0.04;
sigma = 0.15;
mu = -0.2;
delta = 0.1;
lambda = 0.5;
alpha = 0;

#Creation de la matrice bidiagonale B 
create_mat_bidiag <- function(delta_x, points_discr){
  
  delta_t = delta_x * delta_x;

  #Remplissage de la diagonale et de la sous diagonale
  j = 1 + ((sigma*sigma)/2 - r + alpha)*(delta_t/delta_x);
  j_1 = -((sigma*sigma)/2 - r + alpha) * (delta_t/delta_x);
  
  matrix <- diag(j,length(points_discr), length(points_discr));
  for (i in 2:length(points_discr)){
    matrix[i,i-1] = j_1;
  }
  
  return(matrix)

}

#Calcul de u(x,T)
calcul_u <- function(Ar, Al, delta_x, points){
  
  uo = calculate_uo(points);
  delta_t = delta_x * delta_x;
  A = create_mat_tridiag(delta_x, points);
  B = create_mat_bidiag(delta_x, points);

  uprec = uo;
  t = 0;
  while (t < T){
    ucour = resol_syst(uprec, A, B)
    
    #Conditions aux limites
    ucour[1] = h(Al + r*t);
    ucour[length(ucour)] = h(Ar + r*t);
    
    t = t + delta_t;
    uprec = ucour;
  }

  plot(points, ucour, ylab = "u(T)", type="l", main = "Tracé de u(x,T)" );
  lines(points,uo,col='red')
  
  return(ucour);
}

#Calcul de V(S,0)
calcul_v <- function(Ar, Al, delta_x){
  
  #Pour calculer le temps d'exécution
  time1<-Sys.time()
  points = calcul_points_discretisation(Ar, Al, delta_x);
  uT = calcul_u(Ar, Al, delta_x, points)
  v0 = exp(-r*T) * uT;
  
  plot (points, v0, type="l", main="Visualisation de V(S,0)");
  
  time2<-Sys.time()
  Tdiff= difftime(time2, time1)
  
  #A decommenter pour afficher le temps de calcul
  #cat("Temps de calcul : ", Tdiff, "\n")


}
