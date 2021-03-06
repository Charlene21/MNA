S0 = 1;
K = 1;
T=1;
r = 0.04;
sigma = 0.15;
mu = -0.2;
delta = 0.1;
lambda = 0.5;

Bl = -0.5
Br = 0.2

L = 50;
R = 20;

#fonction g
g <- function(y){
  return (lambda * exp(-(y-mu)^2/(2*delta^2))/(delta*sqrt(2*pi)));
}

#Calcul de la constante alpha
calcul_alpha <- function(delta_x){
  
  sum = 0;
  for (j in -L:(R-1)){
    quantityToAdd = (exp(j*delta_x)+exp((1+j)*delta_x) - 2 )*g((j+0.5)*delta_x)
    sum = sum + quantityToAdd;
  }
  
  sum = sum * (delta_x/2);
  
  return (sum);
}

#Calcul de la matrice B multidiagonale
create_mat_sauts <- function(delta_x, alpha, points_discr){
  
  delta_t = delta_x * delta_x;
  
  #Remplissage de la diagonale et la diagonale inferieure comme dans le cas sans sauts
  j = 1 + ((sigma*sigma)/2 - r + alpha) * (delta_t/delta_x); 
  j_1 = -((sigma*sigma)/2 - r + alpha) * (delta_t/delta_x);
  
  matrix <- diag(j,length(points_discr), length(points_discr));
  for (i in 2:length(points_discr)){
    matrix[i,i-1] = j_1;
  }
  
  compteur = 0;
  
  #Remplissage des diagonales supérieures
  for (k in -L:0){
    
    compteur = abs(k);
    
    for(decompte in 1 : (length(points_discr)-abs(k))) {
      
      compteur = compteur + 1;
      
      if(k == -L){
        quantityToAdd = g((k + 0.5)*delta_x)*(delta_x * delta_t/2)
      }
      else {
        quantityToAdd =  (g((k - 1 + 0.5)*delta_x) + g((k + 0.5)*delta_x))*(delta_x * delta_t/2)
      }
      matrix[decompte,compteur] = matrix[decompte,compteur] + quantityToAdd ;
    }
    
  }
  
  #Remplissage des diagonales inférieures
  for (k in 1 : R){
    
    compteur = k;
    for(decompte in 1 : (length(points_discr)-abs(k))){
      compteur = compteur + 1;
      if ( k == R) {
        quantityToAdd =  g((k + 0.5)*delta_x)*(delta_x* delta_t/2);
      }
      else{
        quantityToAdd =  (g((k- 1 + 0.5)*delta_x) + g((k + 0.5)*delta_x))*(delta_x * delta_t/2);
      }
      
      matrix[compteur,decompte] = matrix[compteur,decompte] +  quantityToAdd;
    }
    
  }
  
  return(matrix)
}  


#Calcul de u pour le cas avec sauts
calcul_u_sauts <- function(Ar, Al, delta_x, points){
  
  uo = calculate_uo(points);
  delta_t = delta_x * delta_x;
  alpha = calcul_alpha(delta_x);
  uprec = uo;
  t = 0;
  A = create_mat_tridiag(delta_x, points);
  B = create_mat_sauts(delta_x, alpha, points);
  
  while (t < T){
  
    ucour = resol_syst(uprec, A, B)
    
    #Conditions aux limites
    for (i in 1:length(points)){
      if ((points[i] <= Al )|| (points[i] >= Ar)){
        ucour[i] = h(points[i] + r*t);
      }
    }
    
    t = t + delta_t;
    uprec = ucour;
  }
  
  return(ucour);
}

#Calcul de v pour le cas avec sauts
calcul_v_sauts <- function(Ar, Al, delta_x){
  
  #pour calculer le temps d'exécution
  time1<-Sys.time()
  
  points = calcul_points_discretisation(Ar+Br, Al+Bl, delta_x);
  uT = calcul_u_sauts(Ar, Al, delta_x, points)
  plot(points, uT, type="l");
  v0 = exp(-r*T) * uT;
  
  plot (points, v0, type="l");
  
  time2<-Sys.time()
  Tdiff= difftime(time2, time1) 
  
  #A de commenter pour afficher le temps de calcul
  #cat("Temps de calcul : ", Tdiff, "\n")
  
  
}
