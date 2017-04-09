S0 = 1;
K = 1;
T=1;
r = 0.04;
sigma = 0.15;
mu = -0.2;
delta = 0.1;
lambda = 0.5;
#alpha = 0;

Bl = -0.4
Br = 0.2

L = 40;
R = 20;

g <- function(y){
  return (lambda * exp(-(y-mu)^2/(2*delta^2))/(delta*sqrt(2*pi)));
}

calcul_points_discretisation_sauts <- function(Ar, Al, delta_x){
  #cat(Ar+Br, "\n");
  #cat(Al+Bl, "\n");
  nb_points = (Ar+Br - (Al+Bl))/delta_x;
  #cat("nb poitns : ", nb_points);
  vect_discr = c(Al+Bl);
  x = Al+Bl;
  
  for (i in 1 : nb_points){
    x = x + delta_x;
    vect_discr = c(vect_discr, x);
  }
  
  return(vect_discr);
}

calcul_alpha <- function(delta_x){

  sum = 0;
  for (j in -L:(R-1)){
    quantityToAdd = (exp(j*delta_x)+exp((1+j)*delta_x) - 2 )*g((j+0.5)*delta_x)
    sum = sum + quantityToAdd;
  }
  
  sum = sum * (delta_x/2);
  
  return (sum);
}

create_mat_sauts <- function(Ar, Al, delta_x){
  
  delta_t = delta_x * delta_x;
  alpha = calcul_alpha(delta_x);
  #cat("alpha : ", alpha, "\n")
  points_discr = calcul_points_discretisation_sauts(Ar, Al, delta_x);
  j = 1 + ((sigma*sigma)/2 - r + alpha) * delta_x;
  #cat("j = ", j, "\n")
  j_1 = -((sigma*sigma)/2 - r + alpha) * delta_x;
  
  matrix <- diag(j,length(points_discr), length(points_discr));
  for (i in 2:length(points_discr)){
    matrix[i,i-1] = j_1;
  }
  
  # compteur = -L;
  # for (j in 0:(length(points_discr)-1)){
  #   incr = 0;
  #   compteur = compteur + 1;
  #   for (k in (length(points_discr)-j):length(points_discr)){
  #     incr = incr + 1;
  #     
  #     if(j==0 || j==length(points_discr)-1){
  #       matrix[incr,k] = matrix[incr,k] + g((compteur + 0.5)*delta_x)*(delta_x/2);
  #     } else {
  #       matrix[incr,k] = matrix[incr,k] + (g((compteur + 0.5)*delta_x) + g((compteur + 1 + 0.5)*delta_x))*(delta_x/2) ;
  #     }
  #   }
  # }g()
  
  compteur = 0;
  
  for (k in -L:0){
  
    compteur = abs(k);
   # cat("k = ", k, "\n")
   # cat(">>>>>", length(points_discr)-abs(k), "\n");
   # cat(">>>>> length : ", length(points_discr), "\n");
    for(decompte in 1 : (length(points_discr)-abs(k))) {
     # cat(" ----------------------------decompte = ", decompte, "\n")
      compteur = compteur + 1;
     # cat("decompte,compteur = ", decompte, " , ", compteur, "\n")
      
      if(k == -L){
        
        quantityToAdd = g((k + 0.5)*delta_x)*(delta_x * delta_t/2)
        #cat("Cas de -L, k = ", k, "\n Quantity to add = ", quantityToAdd, "\n")
      }
      else {
       # cat("matrix before : ", matrix[decompte,compteur], "\n");
        quantityToAdd =  (g((k - 1 + 0.5)*delta_x) + g((k + 0.5)*delta_x))*(delta_x * delta_t/2)
        #cat("quantity = ", quantityToAdd, "\n")
      }

      #cat("qtityToAdd / mat before : ", quantityToAdd, "/ ",  matrix[decompte,compteur], "\n")
      matrix[decompte,compteur] = matrix[decompte,compteur] + quantityToAdd ;
     # cat("matrix after = ", matrix[decompte,compteur], "\n");
          }
 
  }
  

  #cat("Deuxième partie \n")
  for (k in 1 : (R)){
   # cat("2eme boucle : ", (R-1), "\n")
    compteur = k;
    for(decompte in 1 : (length(points_discr)-abs(k))){
      compteur = compteur + 1;
      # cat("compteur,decompte = ", compteur, " , ", decompte, "\n")
      if ( k == R) {
        quantityToAdd =  g((k + 0.5)*delta_x)*(delta_x* delta_t/2);
      }
      else{
        quantityToAdd =  (g((k- 1 + 0.5)*delta_x) + g((k + 0.5)*delta_x))*(delta_x * delta_t/2);
      }
      #cat("quantity to add : ", quantityToAdd, "\n")
      #cat("matrix before : ", matrix[compteur,decompte], "\n")
      matrix[compteur,decompte] = matrix[compteur,decompte] +  quantityToAdd;
     # cat("matrix after : ", matrix[compteur,decompte], "\n")
    }
    
  }

  
  
  return(matrix)
}  

create_mat_tridiag_sauts <- function(Ar, Al, delta_x){
  
  delta_t = delta_x * delta_x;
  points_discr = calcul_points_discretisation_sauts(Ar, Al, delta_x);
  #cat("length : ", length(points_discr), "\n")
  j = 1 + (sigma*sigma*delta_t)/(delta_x*delta_x);
  j_1 = -(sigma*sigma*delta_t)/(2*delta_x*delta_x);
  #cat("j = ", j , " / ", "j_1 = ", j_1, "\n")
  matrix <- diag(j,length(points_discr), length(points_discr));
  for (i in 1:length(points_discr)-1){
    matrix[i,i+1] = j_1;
  }
  
  for (i in 2:length(points_discr)){
    matrix[i,i-1] = j_1;
  }
  
  return (matrix);
}

resol_syst_sauts <- function(Ar, Al, delta_x, X_1){
 # cat("Al : ", Al, "Ar : ", Ar)
  A = create_mat_tridiag_sauts(Ar, Al, delta_x);
 # print(A)
  #cat("taille de A : ", nrow(A), "/", ncol(A), "\n");
  B = create_mat_sauts(Ar, Al, delta_x);
  
  #cat("taille de B : ", nrow(B), "/", ncol(B), "\n");
  
  #cat("taille de X_1 : ", length(X_1), "\n");
  X = solve(A, B %*% X_1);
  return(X);
  
}

calcul_u_sauts <- function(Ar, Al, delta_x){
  points = calcul_points_discretisation_sauts(Ar, Al, delta_x);
 # cat("points_u :",  points, "/", length(points), "\n");
  uo = calculate_uo(points);
 # cat("uo :",  uo, "/", length(uo), "\n");
  #uo = rev(uo);

  plot(points, uo, type='l');
  delta_t = delta_x * delta_x;
  
  uprec = uo;
  t = 0;
  while (t < T){
  
    cat("t = ", t, "\n")
    ucour = resol_syst_sauts(Ar, Al, delta_x, uprec)
    #cat("ucour :",  ucour, "/", length(uo), "\n");
    for (i in 1:length(points)){
    # cat("point : ", points[i], "\n", "Al : ", Al , "/", "Ar : ", Ar, "\n");
      if ((points[i] <= Al )|| (points[i] >= Ar)){
       # cat("yo : ", points[i], "\n")
        ucour[i] = h(points[i] + r*t);
      }
    }
 
    t = t + delta_t;
    uprec = ucour;
  }
  
  #ucour = rev(ucour);
  plot(points, ucour, type="l");
  lines(points,uo,col='red')
  
  return(ucour);
}

calcul_v <- function(Ar, Al, delta_x){
  points = calcul_points_discretisation_sauts(Ar, Al, delta_x);
  #cat("points_v :",  points, "/", length(points), "\n");
  uT = calcul_u_sauts(Ar, Al, delta_x)
  plot(points, uT, type="l");
  v0 = exp(-r*T) * uT;
  
  vect = c();
  for (i in 1:length(points)){
    vect=c(vect, h(points[i]));
  }
  
  plot (points, v0, type="l");
  lines(vect, col="red")
  
}
