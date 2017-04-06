S0 = 1;
K = 1;
T=1;
r = 0.04;
sigma = 0.15;
mu = -0.2;
delta = 0.1;
lambda = 0.5;
alpha = 0;

g <- function(y){
  return (lambda * exp(-(y-mu)^2/(2*delta^2))/(delta*sqrt(2*pi)));
}

h <- function(x){
  res = S0*exp(x) - K;
  if (res < 0){
    res = 0;
  }
  
  return(res);
}

f <- function(x, t){
  return(h(x+r*t))
}

calcul_points_discretisation <- function(Ar, Al, delta_x){
  
  nb_points = (Ar - Al)/delta_x;
  vect_discr = c(Al);
  x = Al;
  
  for (i in 1 : nb_points){
    x = x + delta_x;
    vect_discr = c(vect_discr, x);
  }
  
  return(vect_discr);
}

create_mat_tridiag <- function(Ar, Al, delta_x){
  
  delta_t = delta_x * delta_x;
  points_discr = calcul_points_discretisation(Ar, Al, delta_x);
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

create_mat_bidiag <- function(Ar, Al, delta_x){
  
  delta_t = delta_x * delta_x;
  points_discr = calcul_points_discretisation(Ar, Al, delta_x);
  j = 1 + ((sigma*sigma)/2 - r + alpha)*delta_x;
  j_1 = -((sigma*sigma)/2 - r + alpha) * delta_x;
  
  matrix <- diag(j,length(points_discr), length(points_discr));
  for (i in 2:length(points_discr)){
    matrix[i,i-1] = j_1;
  }
  
  return(matrix)

}

resol_syst <- function(Ar, Al, delta_x, X_1){
  A = create_mat_tridiag(Ar, Al, delta_x);
  B = create_mat_bidiag(Ar, Al, delta_x);
  
  X = solve(A, B %*% X_1);
  return(X);
  
}

calculate_uo <- function(points_discr){
  uo = cbind();
  for (i in 1 : length(points_discr)){
    uo = c(uo, h(points_discr[i]));
  }
  
  return (uo);
}


calcul_u <- function(Ar, Al, delta_x){
  points = calcul_points_discretisation(Ar, Al, delta_x);
  uo = calculate_uo(points);
  #uo = rev(uo);
  plot(points, uo, type='l');
  delta_t = delta_x * delta_x;

  uprec = uo;
  t = 0;
  while (t < T){
    ucour = resol_syst(Ar, Al, delta_x, uprec)
    ucour[1] = h(Al + r*t);
    ucour[length(ucour)] = h(Ar + r*t);
    #cat("t : ", t, "\n")
    t = t + delta_t;
    uprec = ucour;
  }

  #ucour = rev(ucour);
  plot(points, ucour, type="l");
  lines(points,uo,col='red')
  
  return(ucour);
}

calcul_v <- function(Ar, Al, delta_x){
  points = calcul_points_discretisation(Ar, Al, delta_x);
  uT = calcul_u(Ar, Al, delta_x)
  v0 = exp(-r*T) * uT;
  
  vect = c();
  for (i in 1:length(points)){
    vect=c(vect, h(points[i]));
  }
  
  plot (points, v0, type="l");
  lines(vect, col="red")

}
