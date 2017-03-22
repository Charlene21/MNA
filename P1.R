g <- function(r, T, S0 ,K, sigma, lambda, delta, mu, y){
  return (lambda * exp(-(y-mu)^2/2*delta)/(delta*sqrt(2*pi)));
}

h <- function(x, S0, K){
  res = S0*exp(x) - K;
  if (res < 0){
    res = 0;
  }
  
  return(res);
}

f <- function(x, r, t, S0, K){
  return(h(x+r*t, S0, K))
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

create_mat_tridiag <- function(Ar, Al, delta_x, sigma, r , alpha){
  
  delta_t = delta_x * delta_x;
  points_discr = calcul_points_discretisation(Ar, Al, delta_x);
  j = 1 + (sigma*sigma*delta_t)/(delta_x*delta_x);
  j_1 = -(sigma*sigma*delta_t)/(2*delta_x*delta_x);
  
  #cat("length : ", length(points_discr));
  matrix <- diag(j,length(points_discr), length(points_discr));
  for (i in 1:length(points_discr)-1){
    matrix[i,i+1] = j_1;
  }
  
  for (i in 2:length(points_discr)){
    matrix[i,i-1] = j_1;
  }
  
  return (matrix);
}

create_mat_bidiag <- function(Ar, Al, delta_x, sigma, r , alpha){
  
  delta_t = delta_x * delta_x;
  points_discr = calcul_points_discretisation(Ar, Al, delta_x);
  j = 1 - (((sigma*sigma)/2 - r + alpha) * delta_t)/delta_x;
  j_1 = -(((sigma*sigma)/2 - r + alpha) * delta_t)/delta_x;
  
  matrix <- diag(j,length(points_discr), length(points_discr));
  for (i in 2:length(points_discr)){
    matrix[i,i-1] = j_1;
  }
  
  return(matrix)

}

resol_syst <- function(Ar, Al, delta_x, sigma, r , alpha, X_1){
  A = create_mat_tridiag(Ar, Al, delta_x, sigma, r , alpha);
  B = create_mat_bidiag(Ar, Al, delta_x, sigma, r , alpha);
  #X = cbind(rep(0,nx));
  
  X = solve(A, B %*% X_1);
  return(X);
  #Solve.tridiag(,diag(A),)
}

calculate_uo <- function(points_discr, S0, K){
  uo = cbind();
  for (i in 1 : length(points_discr)){
    #cat("h : ", h(points_discr[i], S0, K), "points_discr : ", points_discr[i], "\n");
    uo = c(uo, h(points_discr[i], S0, K));
  }
  
  return (uo);
}


calcul_u <- function(Ar, Al, delta_x, S0, K, sigma, r, alpha, T){
  points = calcul_points_discretisation(Ar, Al, delta_x);
  uo = calculate_uo(points, S0, K);
  uo = rev(uo);
  plot(points, uo, type='l');
  delta_t = delta_x * delta_x;

  uprec = uo;
  t = 0;
  while (t < T){
    ucour = resol_syst(Ar, Al, delta_x, sigma, r , alpha, uprec)
    t = t + delta_t;
    uprec = ucour;
  }

  ucour = rev(ucour);
  plot(points, ucour, type="l");
  
  return(ucour);
}

calcul_v <- function(Ar, Al, delta_x, S0, K, sigma, r, alpha, T){
  points = calcul_points_discretisation(Ar, Al, delta_x);
  uT = calcul_u(Ar, Al, delta_x, S0, K, sigma, r, alpha, T)
  v0 = exp(-r*T) * uT;
  
  vect = c();
  for (i in 1:length(points)){
    vect=c(vect, h(points[i],S0,K));
  }
  
  plot (points, v0, type="l");
  lines(vect, col="red")

}
