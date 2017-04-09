S0 = 1;
K = 1;
T=1;
r = 0.04;
sigma = 0.15;
mu = -0.2;
delta = 0.1;
lambda = 0.5;
#alpha = 0;

Bl = -1
Br = 1

L = 100;
R = 100;

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
  
  #cat("alpha : ", alpha, "\n")
  
  j = 1 + ((sigma*sigma)/2 - r + alpha) * delta_x;
  #cat("j = ", j, "\n")
  j_1 = -((sigma*sigma)/2 - r + alpha) * delta_x;
  
  matrix <- diag(j,length(points_discr), length(points_discr));
  for (i in 2:length(points_discr)){
    matrix[i,i-1] = j_1;
  }

  
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





calcul_u_sauts <- function(Ar, Al, delta_x, points){

 # cat("points_u :",  points, "/", length(points), "\n");
  uo = calculate_uo(points);
 # cat("uo :",  uo, "/", length(uo), "\n");
  #uo = rev(uo);

  #plot(points, uo, type='l');
  delta_t = delta_x * delta_x;
  alpha = calcul_alpha(delta_x);
  uprec = uo;
  t = 0;
  A = create_mat_tridiag(delta_x, points);
  B = create_mat_sauts(delta_x, alpha, points);
  while (t < T){
  
    cat("t = ", t, "\n")
    ucour = resol_syst(uprec, A, B)
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

calcul_v_sauts <- function(Ar, Al, delta_x){
  points = calcul_points_discretisation(Ar+Br, Al+Bl, delta_x);

  #cat("points_v :",  points, "/", length(points), "\n");
  uT = calcul_u_sauts(Ar, Al, delta_x, points)
  plot(points, uT, type="l");
  v0 = exp(-r*T) * uT;
  
  # vect = c();
  # for (i in 1:length(points)){
  #   vect=c(vect, h(points[i]));
  # }
  
  plot (points, v0, type="l");
  #lines(vect, col="red")
  
}
