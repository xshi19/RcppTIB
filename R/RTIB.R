# Triangular Input Balance

aak2 <- function(h, ord) {
  n = length(h)
  res = hankelsv(h, ord + 1)
  v = res$v[, ord + 1]
  pad = rep(0, 2 * n)
  Tv = ftm(v[n + 1 - 1:n], c(0, h[1:(n - 1)]))[n + 1 - 1:n]
  hm = fft(fft(c(Tv, pad)) / fft(c(v, pad)), inverse = TRUE) / (3 * n)
  hm = hm[3 * n + 1 - 1:n]
    
  if (all(abs(Im(hm)) < .Machine$double.eps * 1000)) {
      return (Re(hm))
  }
  else {
    return (hm)
  }
  return(hm) 
}

allpastib <- function(lambda, x, real = TRUE) {
  # All pass transfer function: H(z)H*(1 / z*) = I. For example H(z) = 1 / z.
  # All pass TIB satisfies: 
  # AA' + BB' = I; A lower tirangular;
  # A'B + C'D = 0; CA' + DB' = 0;
  # A'A + C'C = I; D'D + B'B = I.
  n = length(lambda)
  res = uhesslr(lambda)
  z = rep(0, n) 
  y = x
  
  for(k in 1:length(x)) {
    w = forwardsolve(res$M, res$N %*% as.matrix(c(x[k], z)))
    z = w[1:n]
    y[k] = w[n + 1]
  }
  if (real) {
    return(list(y = Re(y), z = Re(z)))
  }
  else {
    return(list(y = y, z = z))
  }
}

ams515 <- function() {
  h = c(1, 0.03 * maxflatresponse(0.9, 47, 2499) - 0.1 * maxflatresponse(0.8, 12, 2499));
  
  # scale to ensure it is minimum phase
  minphasescale = 0.1;
  a1 = Re(toeplog(h))
  a2 = Re(toeplog(c(h[1], minphasescale * h[2:length(h)])))
  # plot(a1, type = 'l', log = 'y')
  plot(a2, type = 'l', log = 'y')
  h[2:length(h)] = minphasescale * h[2:length(h)];
  
  # make realizations of this process
  m = 100000
  n = 10
  x = matrix(rnorm(m * n), m, n)
  y = matrix(0, n, m)
  sc = 0.5  # innovation scale
  
  for (k in 1:n) {
    y[k, ] = sc * ftm(x[, k], h)
  }
  
  # to estimate the process using TIBRLS and prescribed poles,
  # first choose some poles. Here we use the shifted Chebyshev nodes
  lambda = 0.5 + 0.5 * cos(pi * seq(1, 99, by = 2) / 100);
  # lambda = fftreduce(h, 40)
  
  Cl = ir2realtibc(lambda, h);
  hl = realtibir(lambda, Cl, length(h));
  
  # plot(abs(hl - h), type = "l")
  
  # TIB RLS
  Crls = matrix(0, n, length(lambda))
  arls = matrix(0, n, 3 * length(h))
  hrls = matrix(0, n, length(h))
  
  for (k in 1:n) {
    print(k)
    t = Sys.time()
    rls = tibrls(lambda, y[k, ], 1)
    print(Sys.time() - t)
    
    Crls[k, ] = rls$C
    
    t = Sys.time()
    arls[k, ] = realtibir(lambda, Crls[k, ], 3 * length(h))
    hrls[k, ] = Re(ar2ir(arls[k, ], length(h)))
    print(Sys.time() - t)
  }
  
  hrlsm = colSums(hrls) / n
  
  plot(hrls[1, 2:length(h)], type = "l", col = "blue", lwd = 1, ylim = c(min(hrls), max(hrls[, 2:length(h)])))
  for (k in 2:n) {
    lines(hrls[k, 2:length(h)], type = "l", col = "blue", lwd = 1)
  }
  lines(h[2:length(h)], type = "l", col = "red", lwd = 4)
  lines(hl[2:length(h)], type = "p", col = "yellow", lwd = 3)
  lines(hrlsm[2:length(h)], type = "l", col = "green", lwd = 3)
  title("TIB RLS")
  
#   # TIB LMS
#   Clms = matrix(0, n, length(lambda))
#   hlms = matrix(0, n, length(h))
#   
#   for (k in 1:n) {
#     print(k)
#     t = Sys.time()
#     lms = tiblmsreal(lambda, y[k, ])
#     print(Sys.time() - t)
#     
#     Clms[k, ] = lms$C
#     
#     t = Sys.time()
#     hlms[k, ] = realtibir(lambda, Clms[k, ], length(h))
#     print(Sys.time() - t)
#   }
#   
#   hlmsm = colSums(hlms) / n
#   
#   plot(hlms[1, 2:length(h)], type = "l", col = "blue", lwd = 1, ylim = c(min(hlms), max(hlms[, 2:length(h)])))
#   for (k in 2:n) {
#     lines(hlms[k, 2:length(h)], type = "l", col = "blue", lwd = 1)
#   }
#   lines(h[2:length(h)], type = "l", col = "red", lwd = 4)
#   lines(hl[2:length(h)], type = "p", col = "yellow", lwd = 3)
#   lines(hlmsm[2:length(h)], type = "l", col = "green", lwd = 3)
#   title("TIB LMS")
  
  # Burg
  ahat = matrix(0, n, 600)
  bhhat = matrix(0, n, length(h))
  for (k in 1:n) {
    print(k)
    t = Sys.time()
    res = burg(y[k, ], 600)
    print(Sys.time() - t)
    
    ahat[k, ] = reflect2ar(res$reflect)
    bhhat[k, ] = Re(ar2ir(ahat[k, ], length(h)))
  }
  
  bhm = colSums(bhhat) / n
  
  plot(bhhat[1, 2:length(h)], type = "l", col = "blue", lwd = 1, ylim = c(min(bhhat), max(bhhat[, 2:length(h)])))
  for (k in 2:n) {
    lines(bhhat[k, 2:length(h)], type = "l", col = "blue", lwd = 1)
  }
  lines(h[2:length(h)], type = "l", col = "red", lwd = 4)
  lines(bhm[2:length(h)], type = "l", col = "green", lwd = 3)
  title("Burg")
}

ams515new <- function(cl) {
  registerDoSNOW(cl)
  
  # generate impluse response
  h = c(1, 0.03 * maxflatresponse(0.9, 47, 2499) - 0.1 * maxflatresponse(0.8, 12, 2499));
  
  # scale to ensure it is minimum phase
  minphasescale = 0.1;
  a = Re(toeplog(c(h[1], minphasescale * h[2:length(h)])))
  # plot(a1, type = 'l', log = 'y')
  plot(a, type = 'l', log = 'y')
  h[2:length(h)] = minphasescale * h[2:length(h)];
  
  # make realizations of this process
  m = 100000
  n = 10
  sc = 0.5
  x = matrix(rnorm(m * n), m, n)
  y = foreach(sc_ = rep(sc, n), 
              x_ = x, 
              h_ = matrix(rep(h, n), length(h), n), 
              .packages = 'RcppTIB') %dopar% (sc_ * ftm(x_, h_))
  
  # lambda = 0.5 + 0.5 * cos(pi * seq(1, 99, by = 2) / 100);
  lambda = fftreduce(h, 50)
  Cl = ir2realtibc(lambda, h);
  hl = realtibir(lambda, Cl, length(h));
  
  # TIB RLS
  rls = foreach(l_ = matrix(rep(lambda, n), length(lambda), n),
                y_ = y,
                na = rep(3 * length(h), n),
                nh = rep(length(h), n),
                .packages = 'RcppTIB') %dopar% (tibrlsall(l_, y_, na, nh))
  
  inforls = foreach(r_ = rls,
                    h_ = matrix(rep(h, n), length(h), n), 
                    .combine = 'c',
                    .packages = 'RcppTIB') %dopar% sqrt(sum(abs(toeplog(r_$ir) - a)^2))
  
  hrlsm = rls[[1]]$ir
  plot(rls[[1]]$ir[2:length(h)], type = "l", col = "blue", lwd = 1, ylim = c(-0.02, 0.01))
  for (k in 2:n) {
    lines(rls[[k]]$ir[2:length(h)], type = "l", col = "blue", lwd = 1)
    hrlsm = hrlsm + rls[[k]]$ir
  }
  hrlsm = hrlsm / n
  lines(h[2:length(h)], type = "l", col = "red", lwd = 4)
  lines(hl[2:length(h)], type = "p", col = "yellow", lwd = 3)
  lines(hrlsm[2:length(h)], type = "l", col = "green", lwd = 3)
  title("TIB RLS")
  
#   # TIB LMS  
#   lms = foreach(l_ = matrix(rep(lambda, n), length(lambda), n),
#                 y_ = y,
#                 nh = rep(length(h), n),
#                 .packages = 'RcppTIB') %dopar% (tiblmsall(l_, y_, nh))
#   
#   infolms = foreach(r_ = lms,
#                     h_ = matrix(rep(h, n), length(h), n), 
#                     .combine = 'c',
#                     .packages = 'RcppTIB') %dopar% sqrt(sum(abs(toeplog(r_$ir) - a)^2))
#   
#   hlmsm = lms[[1]]$ir
#   plot(lms[[1]]$ir[2:length(h)], type = "l", col = "blue", lwd = 1, ylim = c(-0.02, 0.01))
#   for (k in 2:n) {
#     lines(lms[[k]]$ir[2:length(h)], type = "l", col = "blue", lwd = 1)
#     hlmsm = hlmsm + lms[[k]]$ir
#   }
#   hlmsm = hlmsm / n
#   lines(h[2:length(h)], type = "l", col = "red", lwd = 4)
#   lines(hl[2:length(h)], type = "p", col = "yellow", lwd = 3)
#   lines(hlmsm[2:length(h)], type = "l", col = "green", lwd = 3)
#   title("TIB LMS")
  
  # Burgs method
  burgs = foreach(y_ = y,
                  nr = rep(600, n),
                  nh = rep(length(h), n),
                  .packages = 'RcppTIB') %dopar% (burgall(y_, nr, nh))
  
  infoburg = foreach(r_ = burgs,
                     h_ = matrix(rep(h, n), length(h), n), 
                     .combine = 'c',
                     .packages = 'RcppTIB') %dopar% sqrt(sum(abs(toeplog(r_$ir) - a)^2))
  
  hburgm = burgs[[1]]$ir
  plot(burgs[[1]]$ir[2:length(h)], type = "l", col = "blue", lwd = 1, ylim = c(-0.02, 0.015))
  for (k in 2:n) {
    lines(burgs[[k]]$ir[2:length(h)], type = "l", col = "blue", lwd = 1)
    hburgm = hburgm + burgs[[k]]$ir
  }
  hburgm = hburgm / n
  lines(h[2:length(h)], type = "l", col = "red", lwd = 4)
  lines(hl[2:length(h)], type = "p", col = "yellow", lwd = 3)
  lines(hburgm[2:length(h)], type = "l", col = "green", lwd = 3)
  title("Burg")
  
  plot(inforls, type = "l", col = "blue", lwd = 1, 
       ylim = c(min(c(inforls, infoburg)), max(c(inforls, infoburg))))
  # lines(infolms, type = "l", col = "red", lwd = 1)
  lines(infoburg, type = "l", col = "green", lwd = 1)
}

annularsample <- function(r1, r2, n) {
  # Generate random sample (complex) points on an annular with radius r1 and r2
  a = min(r1, r2)^2
  b = abs(r1^2 - r2^2)
  z = sqrt(runif(n, a, b)) * exp(complex(1, 0, runif(n, 0, 2 * pi)))
  return(z)
}

approxcomp <- function(ir = NULL, n = 50) {
  
  if (is.null(ir)) {
    ir = 1 / sqrt(1:1e4)
  }
  
  # exponentially increase index from 1 to length(ir)
  j = floor(2^((0:10000) / 10000 * log(length(ir)) / log(2)))
  si = cumsum(ir)[j]
  ci = toeplog(ir)
    
  # Chebyshev poles in [-1, 1]
  lambda = cos(pi * (2 * (1:n) - 1) / (2 * n))
  Cl = ir2realtibc(lambda, ir)
  hl = realtibir(lambda, Cl, length(ir))
  sl = cumsum(hl)[j]
  cl = toeplog(hl)
  
  # moves the points in the disc to 1
  mu = (lambda + 0.999) / (1 + 0.999 * lambda)
  Cm = ir2realtibc(mu, ir)
  hm = realtibir(mu, Cm, length(ir))
  sm = cumsum(hm)[j]
  cm = toeplog(hm)
  
  # equally distributed along the circle
  nu = 0.99 * exp(complex(1, 0, 2 * pi * (0:(n - 1)) / (n - 1)))
  Cn = ir2realtibc(nu, ir)
  hn = realtibir(nu, Cn, length(ir))
  sn = cumsum(hn)[j]
  cn = toeplog(hn)
  
  xi = (nu + 0.999) / (1 + 0.999 * nu)
  Cx = ir2realtibc(xi, ir)
  hx = realtibir(xi, Cx, length(ir))
  sx = cumsum(hx)[j]
  cx = toeplog(hx)
  
  plot(complex(1, lambda, 0), col = 'blue')
  lines(complex(1, mu, 0), type = 'p', col = 'red')
  lines(nu, type = 'p', col = 'green')
  lines(xi, type = 'p', col = 'yellow')
  title('different poles')
  legend(0.7, 0.99, c('lambda','mu','nu','xi'), 
         lty = c(1, 1, 1, 1), 
         col = c('blue', 'red', 'green', 'yellow'))
  
  ymax = max(c(ir[2:length(ir)], 
               hl[2:length(ir)],
               hm[2:length(ir)],
               hn[2:length(ir)],
               hx[2:length(ir)]))
  ymin = min(c(ir[2:length(ir)], 
               hl[2:length(ir)],
               hm[2:length(ir)],
               hn[2:length(ir)],
               hx[2:length(ir)]))
  plot(ir[2:length(ir)], type = 'l', col = 'black', ylim = c(ymin, ymax))
  lines(hl[2:length(ir)], col = 'blue')
  lines(hm[2:length(ir)], col = 'red')
  lines(hn[2:length(ir)], col = 'green')
  lines(hx[2:length(ir)], col = 'yellow')
  title('impluse response')
  legend(0.7 * length(ir), 0.99 * ymax, 
         c('ir','lambda','mu','nu','xi'), lty = c(1, 1, 1, 1, 1), 
         col = c('black', 'blue', 'red', 'green', 'yellow'))
  
  ymax = max(c(si, sl, sm , sn, sx))
  ymin = min(c(si, sl, sm , sn, sx))
  plot(j, si, type = 'l', col = 'black', ylim = c(ymin, ymax))
  lines(j, sl, col = 'blue')
  lines(j, sm, col = 'red')
  lines(j, sn, col = 'green')
  lines(j, sx, col = 'yellow')
  legend(0.7 * length(si), 0.99 * ymax, 
         c('ir','lambda','mu','nu','xi'), lty = c(1, 1, 1, 1, 1), 
         col = c('black', 'blue', 'red', 'green', 'yellow'))
  title('cumulative impluse response')
  
  ymax = max(Re(c(ci, cl, cm , cn)))
  ymin = min(Re(c(ci, cl, cm , cn)))
  plot(Re(ci), type = 'l', col = 'black', ylim = c(ymin, ymax))
  lines(Re(cl), col = 'blue')
  lines(Re(cm), col = 'red')
  lines(Re(cn), col = 'green')
  legend(0.7 * length(ci), 0.99 * ymax, 
         c('ir','lambda','mu','nu'), lty = c(1, 1, 1, 1), 
         col = c('black', 'blue', 'red', 'green'))
  title('cepstrum')
  # lines(Re(cx), col = 'yellow')
  
  c(sum(abs(cl - ci)^2),
    sum(abs(cm - ci)^2),
    sum(abs(cn - ci)^2),
    sum(abs(cx - ci)^2))
}

approxcomp2 <- function() {
  d = 16
  m = 24
  n = 1e4
  h = 1 / sqrt(1:n)
  sv = hankelsv(h[2:n], m)
  
  # get the AAK approximation with d poles
  hd = aak2(h[2:n], d)
  hd = c(1, Re(hd))
  
  # sinc intepolation points
  lambda = sincpoints(d)
  
  # get poles from balanced truncation of h
  mu = fftreduce(h[2:n], d)
  
  # get poles from balanced truncation of the AAK approximation
  # we use an extra pole to accommodate numerical error in the AAK
  nu = fftreduce(hd[2:n], d + 1)
  
  # get poles from balanced truncation of log h
  rho = exp(.Machine$double.eps / n - pi / sqrt(n))
  format(rho, digits = 2^4)
  hr = h * (rho^(1:length(h) - 1))
  a = Re(toeplog(hr))
  xi = fftreduce(a[2:length(a)], d) / rho
  
  # get poles from AAK approximation of dilated logarithm
  rho = exp(log(.Machine$double.eps) / n)
  format(rho, digits = 2^4)
  ar = a * (rho^(1:length(a)))
  ak = aak2(ar[2:length(ar)], d)
  chi = fftreduce(ak, d + 1)
  
  # coefficients for TIB realizations:
  Cl = ir2realtibc(lambda, h[2:length(h)])
  Cm = ir2realtibc(mu, h[2:length(h)])
  Cn = ir2realtibc(nu, h[2:length(h)])
  Cx = ir2realtibc(xi, h[2:length(h)])
  Cc = ir2realtibc(chi, h[2:length(h)])
  
  # and the impulse responses that correspond:
  hl = c(1, realtibir(lambda, Cl, length(h) - 1))
  hm = c(1, realtibir(mu, Cm, length(h) - 1))
  hn = c(1, realtibir(nu, Cn, length(h) - 1))
  hx = c(1, realtibir(xi, Cx, length(h) - 1))
  hc = c(1, realtibir(chi, Cc, length(h) - 1))
  
  # compute the power series coefficients for the log transfer functions
  a = Re(toeplog(h))
  ad = Re(toeplog(hd))
  al = Re(toeplog(hl))
  am = Re(toeplog(hm))
  an = Re(toeplog(hn))
  ax = Re(toeplog(hx))
  ac = Re(toeplog(hc))
  
  dist = matrix(0, 6, 3)
  rownames(dist) = c('AAK', 'Sinc', 'fftreduce', 'AAK poles', 'fftreduce log', 'AAK log')
  colnames(dist) = c('Info', 'H2', 'H inf')
  
  dist[, 1] = c(sqrt(sum((a - ad)^2)), sqrt(sum((a - al)^2)),
                sqrt(sum((a - am)^2)), sqrt(sum((a - an)^2)),
                sqrt(sum((a - ax)^2)), sqrt(sum((a - ac)^2)))
  dist[, 2] = c(sqrt(sum((h - hd)^2)), sqrt(sum((h - hl)^2)),
                sqrt(sum((h - hm)^2)), sqrt(sum((h - hn)^2)),
                sqrt(sum((h - hx)^2)), sqrt(sum((h - hc)^2)))
  dist[, 3] = c(max(abs(h - hd)), max(abs(h - hl)),
                max(abs(h - hm)), max(abs(h - hn)),
                max(abs(h - hx)), max(abs(h - hc)))
  dist
  
  # Compute the first m singular values of the Hankel matrices of differences:
  cl = makeCluster(2)
  H = cbind(h - hd, h - hl, h - hm, h - hn, h - hx, h - hc)
  qq = parhankelsv(H[2:nrow(H), ], m, cl)
  
  ymax = max(c(sv$d, qq[[1]]$d, qq[[2]]$d, qq[[3]]$d, qq[[4]]$d, qq[[5]]$d, qq[[6]]$d))
  ymin = min(c(sv$d, qq[[1]]$d, qq[[2]]$d, qq[[3]]$d, qq[[4]]$d, qq[[5]]$d, qq[[6]]$d))
  plot(sv$d, log = 'y', type = 'o', col = 'red', ylim = c(ymin, ymax))
  lines(qq[[1]]$d, log = 'y', type = 'o', col = 'blue')
  lines(qq[[2]]$d, log = 'y', type = 'o', col = 'green')
  lines(qq[[3]]$d, log = 'y', type = 'o', col = 'yellow')
  lines(qq[[4]]$d, log = 'y', type = 'o', col = 'black')
  lines(qq[[5]]$d, log = 'y', type = 'o', col = 'purple')
  lines(qq[[6]]$d, log = 'y', type = 'o', col = 'orange')
  legend(0.7 * length(sv$d), 0.99 * ymax, 
         c('True', rownames(dist)), lty = c(1, 1, 1, 1, 1, 1, 1), 
         col = c('red', 'blue', 'green', 'yellow', 'black', 'purple', 'orange'))
  title('Hankel singular values')
  
  
  if (n > 1e4) {
    j = round(2^((1:1e4) / 1e4 * log(n - 1) / log(2))) + 1
  }
  else {
    j = 1:1e4 + 1
  }
  
  # ymax = max(c(h[j], hd[j], hl[j], hm[j], hn[j], hx[j], hc[j]))
  # ymin = min(c(h[j], hd[j], hl[j], hm[j], hn[j], hx[j], hc[j]))
  ymax = 10 * max(abs(h[2:length(h)]))
  ymin = 0.01 * min(abs(h[2:length(h)]))
  xmax = 1.1 * (n - 1)
  xmin = 1
  plot(h[j], log = 'xy', type = 'l', col = 'red', ylim = c(ymin, ymax))
  lines(hd[j], type = 'l', col = 'blue')
  lines(hl[j], type = 'l', col = 'green')
  lines(hm[j], type = 'l', col = 'yellow')
  lines(hn[j], type = 'l', col = 'black')
  lines(hx[j], type = 'l', col = 'purple')
  lines(hc[j], type = 'l', col = 'orange')
  legend(1, 0.02 * ymax, 
         c('True', rownames(dist)), lty = c(1, 1, 1, 1, 1, 1, 1), 
         col = c('red', 'blue', 'green', 'yellow', 'black', 'purple', 'orange'))
  title('Impluse response')
}

ar2ir <- function(a, m) {
  # Estimate the impulse reponse given the ar coefficients 
  # OR
  # Estimate the ar coefficients given the impulse reponse
  
  p = length(a)
  a = as.vector(a)
  
  if ((m - 1) > p) {
    a = c(a, rep(0, (m - 1 -p)))  
  }
  else{
    a = a[1:(m - 1)]
  }
  
  h = toepsolve(c(1, -a), c(1, rep(0, (m - 1))))
}

ar2reflect <- function(a) {
  n = length(a)
  reflect = rep(0, n)
  
  for (k in (n + 2 - 2:n)) {
    reflect[k] = -a[n]
    a[1:(n - 1)] = (a[1:(n - 1)] - reflect[k] * a[n - 1:(n - 1)]) / (1 - reflect[k]^2)
    n = n - 1
  }
  
  reflect[1] = -a[n]
  return(reflect)
}

blkhankel <- function(v, m) {
  # Generate block Hankel matrix with first block column v
  # m is the row number of the blocks
  n = ncol(v)
  p = ceiling(nrow(v) / m)
  mh = p * m
  nh = p * n
  
  if(mh > nrow(v)) {
    v = rbind(v, matrix(0, mh - nrow(v), n))
  }
  
  H = matrix(0, mh, nh)
  H[, 1:n] = v
  
  sT2 = n
  j = m + 1:mh
  while(2 * sT2 < nh) {
    H[j <= mh, sT2 + 1:sT2] = H[j[j <= mh], 1:sT2]
    j[j <= mh] = j[j[j <= mh]]
    sT2 = sT2 * 2
  }
  H[j <= mh, (sT2 + 1):nh] = H[j[j <= mh], 1:(nh - sT2)]
  return(H)
}

burg <- function(y, ord) {
  n = length(y)
  
  if (ord >= n) {
    warning("ord should not be greater than sample length")
    return(NULL)
  }
  
  reflect = rep(0, ord + 1)    # reflection coefficient
  sigma = 2 * sum(y^2)
  ferror = matrix(0, n, 2)    # forward error
  berror = matrix(0, n, 2)    # backward error
  ferror[, 1] = y
  berror[, 1] = y
  flip = 1
  
  for (k in 1:ord) {
    sigma = (1 - reflect[k]^2) * sigma - (ferror[k, flip]^2 + berror[n, flip]^2)
    tau = 2 * sum(ferror[(k + 1):n, flip] * berror[k:(n - 1), flip])
    reflect[k] = -tau / sigma
    ferror[(k + 1):n, 3 - flip] = ferror[(k + 1):n, flip] + reflect[k] * berror[k:(n - 1), flip]
    berror[(k + 1):n, 3 - flip] = berror[k:(n - 1), flip] + reflect[k] * ferror[(k + 1):n, flip]
    flip = 3 - flip
  }
  
  return(list(reflect = reflect[2:(ord + 1)], alpha = sigma / (2 * (n - ord)))) 
}

burgall <- function(y, nr, nh) {
  # wrapped function for paralell
  res = burg(y, nr)
  ar = reflect2ar(res$reflect)
  ir = Re(ar2ir(ar, nh))
  return(list(ar = ar, ir = ir))
}

cmax <- function(x, y = NULL) {
  x = c(x, y)
  if (all(Im(x) == 0)) {
    return(max(Re(x)))
  }
  else {
    return(x[order(abs(x), decreasing = TRUE)[1]])
  }
}

cmin <- function(x, y = NULL) {
  x = c(x, y)
  if (all(Im(x) == 0)) {
    return(min(Re(x)))
  }
  else {
    return(x[order(abs(x), decreasing = FALSE)[1]])
  }
}

fftreduce <- function(h, ord) {
  n = length(h)
  h = c(h, rep(0, n))
  hhat = fft(h[2 * n + 1 - 1:(2 * n)])
  
  # arpack for eigen values, same as eigs in Matlab
  options = list(n = n, nev = ord, ncv = min(n, max(2 * ord + 1, 20)),
                 which = 'LM', tol = 10 * .Machine$double.eps)
  res = arpack(ffthm, extra = hhat, sym = TRUE, options = options)
  
  q = qr.Q(qr(res$vectors))
  res = eigen(t(q[1:(n - 1), ]) %*% q[2:n, ])
  lambda = res$values
  return(lambda)
}

ffthm <- function(x, hhat) {
  #   Multiply x by the Hankel matrix with first column (flipped) h
  #   hhat is the Fourier transform of (flipped) h
  #   Example:
#       n = 1000
#       h = complex(1, runif(n), runif(n))
#       x = complex(1, runif(n), runif(n))
#       H = trihankel(h)
#     
#       s0 = system.time( {y0 = H %*% x} )
#       s1 = system.time( {hhat = fft(c(rep(0, n), h[n + 1 - 1:n]))
#                         y1 = ffthm(x, hhat)})
#       mean(abs(y0 - y1)^2)
  
  m = length(hhat)
  n = length(x)
  x = c(x, rep(0, m - n))
  xhat = mvfft(as.matrix(x))
  y = mvfft(xhat * as.matrix(hhat), inverse = TRUE) / m
  y = y[m + 1 - 1:n]
  if (all(abs(Im(y)) < .Machine$double.eps * 10000)) {
    return (Re(y))
  }
  else {
    return (y)
  }
}

ftm <- function(x, t, cutoff = 128) {
  #   Multiply x by the lower triangular Toeplitz matrix with first column t
  #   Example:
  #     n = 10000
  #     t = complex(1, runif(n), runif(n))
  #     x = complex(1, runif(n), runif(n))
  #     Toe = tritoep(t)  
  #     
  #     s0 = system.time( {y0 = Toe %*% x} )
  #     s1 = system.time( {y1 = convolve(t, rev(Conj(x)), type = "open")} )
  #     s2 = system.time( {y2 = ftm(x, t)} )
  
  n = length(x)
  m = length(t)
  if (n < cutoff) {
    y = convolve(t, rev(Conj(x)), type = "open")
  }
  else {
    t = c(t, rep(0, n))
    x = c(x, rep(0, m))
    y = fft(fft(t) * fft(x), inverse = TRUE) / (m + n)
  }
  y = y[1:n]
  if (all(abs(Im(y)) < .Machine$double.eps * 10000)) {
    return (Re(y))
  }
  else {
    return (y)
  }
}

hankel <- function(x, y = NULL) {
  # Generate Hankel matrix using outer product
  # Same as Matlab function hankel
  # http://stat.ethz.ch/R-manual/R-devel/library/base/html/outer.html
  m = length(x)
  if (is.null(y)) {
    y = rep(0, m)
  }
  n = length(y)
  h = c(x, y[2:n])
  H = outer(1:m, 1:n, function(u, v) h[u + v - 1])
  return(H)
}

hankelsv <- function(h, ord) {
#   SVD of a Hankel matrix with first column h
#   m = 5000
#   n = 100
#   # h = rnorm(m)
#   h = 1 / sqrt(1:m)
#   unix.time({r1 = hankelsv(h, n)})
#   unix.time({H = trihankel(h)
#              r2 = svds(H, n)})
#   norm(as.matrix(r1$d - r2$d))
  
  n = length(h)
  h = c(h, rep(0, n))
  hhat = fft(h[2 * n + 1 - 1:(2 * n)])
  
  # arpack for eigen values, same as eigs in Matlab
  options = list(n = n, nev = ord, ncv = min(n, max(2 * ord + 1, 20)),
                 which = 'LM', tol = 10 * .Machine$double.eps)
  res = arpack(ffthm, extra = hhat, sym = TRUE, options = options)
  
  d = abs(res$values)
  j = order(d, decreasing = TRUE)
  v = qr.Q(qr(res$vectors[, j]))
  return(list(v = v, d = d[j]))
}

ir2realtibc <- function(lambda, ir) {
  # Computes the least squares approximation to a given impulse response
  # Inputs:
  #   lambda: poles of a SISO TIB system
  #   ir: given impulse response at lags 1,2,3,...
  # Outputs:
  #   C: least squares coefficient
  #
  # NB the approximation may be limited if the impulse response
  # has not decayed sufficiently in the lags in the input
  
  n = length(lambda)
  m = length(ir)
  p = ceiling(sqrt(m))
  
  tib = realtib(lambda)
  A = solvetibM(tib$M, timestibN(tib$N, tib$U))
  B = solvetibM(tib$M, as.matrix(c(1, rep(0, n - 1))))
  if (all(Im(A) == 0)) {A = Re(A)}
  if (all(Im(B) == 0)) {B = Re(B)}
  
  # Let K = [B, AB, A^2B,...]
  # min ||h - C * K||^2
  # least squares solution: C = (K' * K) \ K * h
  # K' * K = 1 since input balance, 
  # so C = K * h = B * h0 + A * B * h1 + A^2 * B * h2 + ....
  
  K = B
  while (ncol(K) < p) {
    K = cbind(K, A %*% K)
    A = A %*% A
  }
  
  sk2 = ncol(K)
  su2 = ceiling(m / sk2)
  h = c(Conj(ir), rep(0, su2 * sk2 - m))
  h = matrix(h, sk2, su2) 
  u = K %*% h
  
  while (su2 > 1 + ceiling(n / 4)) { 
    if (mod(su2, 2) == 0) {
      ind = 2 * (1:(su2 / 2))
      u = u[, ind - 1] + A %*% u[, ind]
    }
    else {
      ind = 2 * (1:((su2 - 1) / 2))
      u = u[, c(1, ind + 1)] + cbind(A %*% u[, ind], 0)
    }
    su2 = ncol(u)
    A = A %*% A
  }
  
  C = u[, su2]
  j = su2 - 1
  while (j > 0) {
    C = u[, j] + A %*% C
    j = j - 1
  }
  return(t(C))
}

ir2tibc <- function(lambda, ir) {
  n = length(lambda)
  tib = sptib_cpp(lambda)
  v = solvetibM(tib$M, matrix(c(1, rep(0, n - 1))))
  C = rep(0, n)
  for (i in 1:length(ir)) {
    C = C + ir[i] * v
    v = solvetibM(tib$M, timestibN(tib$N, v))
  }
  return(C)
}

maxflatresponse <- function(lambda, d, m) {
  mu = -log(lambda)
  y = pgamma(mu * seq(0, m - 1), d, lower = FALSE)
  return(y)
}

parhankelsv <- function(h, ord, cl) {
  registerDoSNOW(cl)
  n = ncol(h)
  res = foreach(h_ = h, ord_ = rep(ord, n),
                .packages = 'RcppTIB') %dopar% (hankelsv(h_, ord_))
  return(res)
}

realtib <- function(lambda) {
  # Reference: Andrew P. Mullhaupt and Kurt S. Riedel, 
  #            Band Matrix Representation of Triangular Input Balanced Form, 
  #            IEEE Transactions on Automatic Control, 1998
  # A * A' + B * B' = I
  # where A = M^(-1) * N, M, N band bidiagonal matrix
  
  re = abs(Im(lambda)) <= 10 * .Machine$double.eps * abs(lambda)
  rho = Re(lambda[re])   # real eigenvalues/poles
  cc = lambda[!re]     
  cc = cc[Im(cc) > 0]    # complex eigenvalues/poles
  m = length(rho)
  
  acc = abs(cc)
  psi = acos(2 * acc / (1 + acc^2)) - acos(2 * Re(cc) / (1 + acc^2))
  cc = matrix(c(Conj(cc), cc), length(cc), 2)
  cc = c(t(cc))
  mu = c(rho, abs(cc))
  tib = sptib_cpp(mu)
  
  tib$U = diag(length(lambda))
  if (length(psi) > 0) {
    for (k in 1:length(psi)) {
      j = m + 2 * (k - 1) + 1:2
      G = matrix(c(cos(psi[k]), sin(psi[k]), -sin(psi[k]), cos(psi[k])), 2, 2)
      tib$U[j, j] = G
    }
  }
  return(tib)
}

realtibir <- function(lambda, C, m) {
  # C is the coefficient for a given impulse response  ir
  # lambda: poles of a SISO TIB system
  # ir the impulse response at lags 1 through m
  # ir = [C*B C*A*B C*A^2*B ... C*A^m*B]
  #   
  # NB: the lag 0 response is not computed by this function
  # nor does the data for this function specify the lag 0 response
  # although it is 1 for the innovations response
  
  n = length(lambda)
  tib = realtib(lambda)
  A = solvetibM(tib$M, timestibN(tib$N, tib$U))
  B = solvetibM(tib$M, as.matrix(c(1, rep(0, n - 1))))
  
  if (all(Im(A) == 0)) {A = Re(A)}
  if (all(Im(B) == 0)) {B = Re(B)}
  
  # calculate sizes of J and K
  p = ceiling(sqrt(m))
  x = 2 ^ ceiling(log(p) / log(2))      # inf {2^k >= sqrt(m)}
  y = 2 ^ ceiling(log(m / x) / log(2))  # inf {2^k >= m / x}
  J = matrix(0, x, n)
  K = matrix(0, n, y)
  
  # compute [C C*A C*A^2 ... C*A^(x-1)]', 
  J[1, ] = C
  nJ = 1
  while (nJ < p) {
    J[nJ + 1:nJ, ] = J[1:nJ, ] %*% A
    nJ = 2 * nJ
    A = A %*% A
  }
  
  # compute [B (A^x)*B (A^x)^2*B ... (A^x)^(y-1)*B]
  K[, 1] = B
  nK = 1
  while (x * nK < m) {
    K[, nK + 1:nK] = A %*% K[, 1:nK]
    nK = 2 * nK
    A = A %*% A
  }
  
  ir = as.vector(J %*% K)
  return(Re(ir[1:m]))
}

reflect2ar <- function(reflect) {
  n = length(reflect)
  a = rep(0, n + 1)
  a[1] = 1
  for (k in 1:n) {
    a[1:(k + 1)] = a[1:(k + 1)] + reflect[k] * a[k + 2 - 1:(k + 1)]
  }
  return(-a[2:(n + 1)])
}

simple_reduction <- function(ir, ord) {
  #   Computes the spectrum for the balanced truncation:
  #   ir = given impulse response
  #   ord = desired order of the truncation
  ir = as.vector(ir)
  m = length(ir)
  res = svd(trihankel(ir), nu = 0, nv = ord)
  A = t(res$v[1:(m - 1), ]) %*% res$v[2:m, ]
  lambda = eig(A)
  return(lambda)
}

example_reduction <- function(n = 20) {
  # Sample some poles
  lambda = annularsample(0, 1, n)
  
  # Pick random coefficients
  C = sqrt(1 - abs(lambda)^2)
  C = C * runif(n)
  
  # Compute the impulse response long enough to decay
  m = 2 * n
  h = Re(tibir(lambda, C, m))
  
  while (abs(h[length(h)]) > .Machine$double.eps * 100) {
    m = ceil(1.5 * m)
    h = Re(tibir(lambda, C, m))
    if (m > 3000) {
      stop("impluse response too long")  
    }
  }
  
  # compute coefficients for conjugate pole tib
  lc = c(lambda, Conj(lambda))
  lc = lc[order(abs(lc))]
  Cl = ir2tibc(lc, h)
  hl = tibir(lc, Cl, m)
  
  # find reduced poles and compute coefficients and impulse response
  mu = sort(svd_reduction(h, 2 * n));
  Cm = ir2tibc(mu, h)
  hm = tibir(mu, Cm, m)
  
  # find reduced poles and compute coefficients and impulse response
  nu = sort(simple_reduction(h, 2 * n));
  Cn = ir2tibc(nu, h)
  hn = tibir(nu, Cn, m)
  
  plot(Re(h), type = "l", col = "blue")
  lines(Re(hm), type = "l", col = "red")
  lines(Re(hn), type = "l", col = "green")
  
  plot(abs(h - hm), type = "l", col = "blue")
  lines(abs(h - hn), type = "l", col = "red")
  
  plot(lc, type = "p", col = "blue", lwd = 10)
  lines(mu, type = "p", col = "red", lwd = 5)
  lines(nu, type = "p", col = "green", lwd = 2)
  
  plot(abs(Cl - Cm), type = "l", col = "blue")
  lines(abs(Cl - Cn), type = "l", col = "red")
}

sptib <- function(lambda) {
  # Compute bidiagonal matrix representation of scalar TIB system with eigenvalues lambda. 
  # Repeated eigenvalues are accepted.
  # Reference: Andrew P. Mullhaupt and Kurt S. Riedel, 
  #            Band Matrix Representation of Triangular Input Balanced Form, 
  #            IEEE Transactions on Automatic Control, 1998
  
  lambda = as.vector(lambda)
  n = length(lambda)
  
  # equation (3.3)
  c_ =  1 / sqrt(1 - abs(lambda)^2)
  s = lambda * c_
  # if (class(s) == "complex") {
  b = diag(c_)
  diag(b[2:n, 1:(n - 1)]) = Conj(s[1:(n - 1)])
  M = b / c_ / sqrt(1 - abs(lambda[1])^2)
  b = diag(s)
  diag(b[2:n, 1:(n - 1)]) = Conj(c_[1:(n - 1)])
  N = b / c_ / sqrt(1 - abs(lambda[1])^2)
  M = M %*% diag((-1)^(0:(n - 1)))
  N = N %*% diag((-1)^(0:(n - 1)))
  #   }
  #   else  {
  #     C = .sparseDiagonal(n, c_)
  #     M = solve(C, bandSparse(n, k = c(-1, 0), diag = list(Conj(s), c_)), sparse = TRUE)
  #     N = solve(C, bandSparse(n, k = c(-1, 0), diag = list(c_, s)), sparse = TRUE)
  #     M = M / sqrt(1 - abs(lambda[1])^2)
  #     N = N / sqrt(1 - abs(lambda[1])^2)
  #     
  #     # Apply signature to conform to Meixner function definition
  #     S = .sparseDiagonal(n, (-1)^(0:(n - 1)))
  #     M = M %*% S
  #     N = N %*% S
  #   }
  #   
  return(list(M = M, N = N))
}

sincpoints <- function(n, h = NULL) {
  if (n == 1) {
    return(1 / 2)
  }
  else {
    if (is.null(h)) {
      h = pi / sqrt((n - 1) / 2)
    }
    x = 0:(n - 1) - (n - 1) / 2
    return(1 / (1 + exp(h * x)))
  }
}

svd_reduction <- function(ir, ord) {
  #   Computes the spectrum for the balanced truncation:
  #   ir = given impulse response
  #   ord = desired order of the truncation
  #   using qr to ensure orthogonality
  ir = as.vector(ir)
  m = length(ir)
  res = svd(trihankel(ir), nu = 0, nv = ord)
  res = qr(res$v)
  v = qr.Q(res)
  A = t(v[1:(m - 1), ]) %*% v[2:m, ]
  lambda = eig(A)
  return(lambda)
}

tibir <- function(lambda, C ,m) {
  n = length(lambda)
  ir = rep(0, m)
  tib = sptib_cpp(lambda)
  v = solvetibM(tib$M, matrix(c(1, rep(0, n - 1))))
  for (i in 1:m) {
    ir[i] = sum(C * Conj(v))
    v = solvetibM(tib$M, timestibN(tib$N, v))
  }
  return(ir)
}

tiblmsall <- function(lambda, y, nh) {
  # wrapped function for paralell
  res = tiblmsreal(lambda, y)
  ir = c(1, realtibir(lambda, res$C, nh - 1))
  return(list(C = res$C, ir = ir))
}

tiblmsreal <- function(lambda, y, post = 0) {
  # Least mean squares
  # min 1 / 2 * (y(k) - C(k) * z(k))^2
  # gradient = C(k) * (z(k) * z(k)') - y(k)' * z(k) 
  # update:
  #   C(k + 1) = C(k) - coeff * (C(k) * (z(k) * z(k)') - y(k)' * z(k))
  # here we assume z(k) * z(k)' = I and coeff = 1 / (k + 1)
  
  m = length(y)
  n = length(lambda)
  
  tib = realtib(lambda)
  # M = as.matrix(tib$M)
  # N = as.matrix(tib$N)
  # U = as.matrix(tib$U)
  
  C = rep(0, n)
  z = rep(0, n)
  x = rep(0, m)
  xi = rep(0, m)
  yhat = rep(0, m)
  st = rep(1, m)
  
  r2 = sqrt(qchisq(0.9999, n))
  sigma = sqrt(var(y[1:50]))
  cwgt = 1
  
  for (k in 1:length(y)) {
    yhat[k] = sum(C * z)
    w = Conj(Conj(y[k]) * z)
    C = C + (w - C) / (k + 1)
    
    # compute residuals or posterior residuals
    if (post == 1) {
      x[k] = y[k] - yhat[k]
    }
    else {
      x[k] = y[k] - sum(C * z)
    } 
    
    z1 = timestibN(tib$N, tib$U %*% z)
    # z1 = N %*% (U %*% z)
    z1[1] = z1[1] + cmin(2.5, cmax(-2.5, x[k] / sigma))
    z1 = solvetibM(tib$M, as.matrix(z1))
    # z1 = forwardsolve.cpp(M, as.matrix(z1))[[1]]
    
    if (sqrt(sum(abs(z1)^2)) > r2) {
      print("Stabilizing...")
      st[k] = 1
      
      if (sqrt(sum(abs(z1)^2)) > 0) {
        Cu = C / sqrt(sum(abs(C)^2))
        z = z - sum(Cu * z) * Conj(Cu)  # remove C' component
        if (n > sum(abs(z)^2)) {
          dz = sqrt(n - sum(abs(z)^2)) * Conj(Cu)  # add C' component to length = n
          z = z + dz
        }
        else {
          z = sqrt(n) * z / sqrt(sum(abs(z)^2)) 
        }
      }
    }
    else {
      if (k > 1) {
        st[k] = st[k - 1] + 1
      }
      z = z1
    }
    
    wgt = log(k + 1) * (1 - 1 / (2 * st[k] - 1))
    sigma = sqrt((cwgt * sigma^2 + wgt * min(20, abs(x[k])^2)) / (cwgt + wgt))
    cwgt = cwgt + wgt
  }
  return(list(C = C / sigma, yhat = yhat, z = z, sigma = sigma))
}

tibrls <- function(lambda, y, m, draw = FALSE) {
  # Recursive least squares
  # estimate C from Y = Z * C + X
  # Cholesky decomposition: [Z(k), Y(k)]' * [Z(k), Y(k)] = L(k) * L(k)'
  # update Z(k+1) = [Z(k), z(k+1)], Y(k+1) = [Y(k), y(k+1)],
  # z(k+1) = A*z(k) + B*x(k)
  # and R(k+1) is obtained by rank 1 update for Cholesky
  # here we treat y as x, so the result C is of the invese system
  
  y = as.vector(y)
  p = length(lambda)
  L = diag(p + 1) * sqrt(var(y[1:60]))
  tib = realtib(lambda)
  z = matrix(0, p, 1)
  
  #   st1 = 0
  #   st2 = 0
  for (k in (p + 1):length(y)) {
    x = c(z, y[k])
    #     t = Sys.time()
    L = cholupdate(L, x)
    #     st1 = st1 + Sys.time() - t
    
    #     t = Sys.time()
    C = backsolve(t(L[1:p, 1:p]), L[p + 1, 1:p])
    
    z = timestibN(tib$N, tib$U %*% z)
    z[1] = z[1] + y[k]
    z = solvetibM(tib$M, z)
    #     st2 = st2 + Sys.time() - t
    if (draw && (k %% 100 == 0)) {
      h = realtibir(lambda, C, m)
      plot(h)
    }
  }
  #   print(c(st1, st2))
  v = L[p + 1, p + 1]
  h = realtibir(lambda, C, m)
  return(list(C = C, h = h, v = v))
}

tibrls_bkt <- function(lambda, y, m, info = FALSE) {
  # Backtest recursive least squares
  # The predictions are Cz, CAz, CA^2z....
  # Num of 'correct' predictions = Num of (sign(perdictions) == sign(y))
  # Num of 'wrong' predictions = Num of (sign(perdictions) == -sign(y))
  # Prediction rate = Num of 'correct' / (Num of 'correct' + Num of 'wrong')
  # m is the number of multi-steps predicitons
  y = as.vector(y)
  p = length(lambda)
  L = diag(p + 1) * sqrt(var(y[1:60]))
  tib = realtib(lambda)
  z = matrix(0, p, 1)
  
  mu = matrix(0, m, 1)
  correct = matrix(0, m, 1)
  wrong = matrix(0, m, 1)
  
  for (k in (p + 1):(length(y) - m)) {
    x = c(z, y[k])
    L = cholupdate(L, x)
    C = backsolve(t(L[1:p, 1:p]), L[p + 1, 1:p])
    
    z = timestibN(tib$N, tib$U %*% z)
    z[1] = z[1] + y[k]
    z = solvetibM(tib$M, z)

    # Ainv = A - BC = M \ N * U
    # Binv = B = M \ e1
    # A = Ainv + BC = M \ (N * U + e1 * C)
    A = timestibN(tib$N, tib$U)
    A[1, ] = A[1, ] + C
    A = solvetibM(tib$M, A)
    Ap = diag(p)
    for (j in 1:m) {
      mu[j] = Re(t(C) %*% Ap %*% z)
      Ap = Ap %*% A 
    }
    g = sign(cumsum(mu)) * sign(cumsum(y[k + (1:m)]))
    correct[g == 1] = correct[g == 1] + 1
    wrong[g == -1] = wrong[g == -1] + 1
    
    if (info && (k %% 100 == 0)) {
      print(k / 100)
      print(correct / (correct + wrong))
    }
  }
  rate = correct / (correct + wrong)
  return(list(rate = rate, correct = correct, wrong = wrong, C = C))
}

tibrlsall <- function(lambda, y, na, nh) {
  # wrapped function for paralell
  res = tibrls(lambda, y, 1)
  ar = realtibir(lambda, res$C, na)
  ir = Re(ar2ir(ar, nh))
  return(list(C = res$C, ar = ar, ir = ir))
}

toeplog <- function(f) {
  # Compute power series coefficients for log(f)
  n = length(f)
  d = 1:(n - 1)
  a = c(log(f[1]), toepsolve(f[1:(n - 1)], d * f[2:n]) / d)
  return(a)
}

toepsolve <- function(t, x, cutoff = 512, Toe = NULL) {
  #   Simple fast solver for the (lower triangluar) Toeplitz system
  #   Input:
  #     t: first column of the Toeplitz matrix
  #     x: corresponding vector
  #     cutoff and Toe: auxillary variable for iteration, do not touch
  #   Output:
  #     y: solution
  #   
  #   Example:
  #   n = 5000
  #   t = complex(1, 1:n, 1)
  #   x = complex(1, runif(n), runif(n))
  #   Toe = tritoep(t)
  #   
  #   s1 = system.time( {y1 = toepsolve(t, x)} )
  #   s2 = system.time( {y2 = solve(Toe, as.matrix(x))} )
  #   s3 = system.time( {y3 = forwardsolve.cpp(Toe, as.matrix(x))[[1]]} )
  #   
  #   mean(abs(x - Toe %*% y1))
  #   mean(abs(x - Toe %*% y2))
  #   mean(abs(x - Toe %*% y3))
  
  n = length(x)
  
  # check if n = 2^k
  n2 = 2^ceiling(log(n) / log(2))
  if (n < n2) {
    # if not, then force dimension to be the power of 2 by adding n2 - n zeros
    pad = rep(0, n2 - n)
    y = toepsolve(c(t, pad), c(pad, x), cutoff)
    y = y[(1 + n2 - n):n2]
    return(y)
  }
  else {
    # here we can assume that n = 2^k
    # if dimension is no larger than cutoff, then solve it directly
    if (n < cutoff) {
      if (is.null(Toe)) {
        y = forwardsolve(tritoep(t), as.matrix(x))
      }
      else {
        y = forwardsolve(Toe, as.matrix(x))
      }
    }
    else {
      # we can divide the toeplitz matrix into two parts and solve them seperately
      # T11 * y1 = x1
      # T21 * y1 + T22 * y2 = x2
      # y1 = T11 \ x1
      # y2 = T22 \ (x2 - T21 * y1)
      # for toeplitz matrix T11 = T22 is triangular toeplitz
      # T21 is toeplitz and T2 * y1 can be computed via FFT
      top = 1:(n / 2)
      bottom = top + n / 2
      y = rep(0, n)
      if (is.null(Toe)) {
        Toe = tritoep(t[1:(cutoff / 2)])
      }
      y[top] = toepsolve(t[top], x[top], Toe = Toe)
      
      # fft for T21 * y1
      w = mvfft(mvfft(as.matrix(t)) * mvfft(as.matrix(y)), inverse = TRUE) / n
      y[bottom] = toepsolve(t[top], x[bottom] - w[bottom], cutoff, Toe = Toe)
    }
  }
  return (y)
}

trihankel <- function(h) {
  # Generate Hankel matrix with first column h
  n = length(h)
  H = matrix(0, n, n)
  H[, 1] = h
  
  sT2 = 1
  j = 2:(n + 1)
  while (sT2 * 2 < n) {
    H[j <= n, sT2 + 1:sT2] = H[j[j <= n], 1:sT2]
    j[j <= n] = j[j[j <= n]]
    sT2 = sT2 * 2
  }
  
  sn = 1:(n - sT2)
  j = j[sn]
  H[sn[j <= n], sT2 + sn] = H[j[j <= n], sn]
  return (H)
}

tritoep <- function(v) {
  # Construct the lower triangular Toeplitz matrix with first column v
  
  m = length(v)
  Toe = matrix(0, m, m)
  Toe[, 1] = v
  sT2 = 1
  j = c(m, 1:(m - 1))
  
  while (2 * sT2 < m) {
    # j = m, 1, 2, 3, ..., m-1
    # 1st loop constructs the 2nd column
    # j = m, m, 1, 2, ..., m-2
    # 2nd loop constructs the 3rd-4th columns
    # j = m, m, m, 1, ..., m-3
    # 3rd loop constructs the 5th-8th columns
    # ...
    Toe[j < m, 1:sT2 + sT2] = Toe[j[j < m], 1:sT2]
    j[j < m] = j[j[j < m]]
    sT2 = sT2 * 2
  }
  
  # construct the rest (sT2 + 1):m columns
  sm = (sT2 + 1):m
  j = j[sm]
  Toe[sm[j < m], sm] = Toe[j[j < m], 1:(m - sT2)]
  return(Toe)
}

uhesslr <- function(lambda) {
  # Generate M and N such that [B, A; D, C] = M \ N is Hessenberg
  n = length(lambda)
  c_ = 1 / sqrt(1 - abs(lambda)^2)
  s = lambda * c_
  M = diag(c(c_, 1))
  diag(M[2:(n + 1), 1:n]) = Conj(s)
  N = diag(c(1, c_))
  diag(N[1:n, 2:(n + 1)]) = s
  M = M * c_[1]
  N = N * c_[1]
  return(list(M = M, N = N))
}

ziplog <- function(f) {
  # Compute power series coefficients for log(f)
  n = length(f)
  d = 1:(n - 1)
  a = c(log(f[1]), zipsolve(f[1:(n - 1)], d * f[2:n]) / d)
  return(a)
}

zipsolve <- function(t, y = NULL, cutoff = 64) {
  #   Solve lower triangular Toeplitz system Tx = y
  #   when called with one argument, y = e1 = [1; 0; 0; ... 0]
  #   Let A, B be the toeplitz matrices with first column a, b,
  #   then A*b = B*a = conv(a, b). So x = Y*T^(-1)*e1 = T^(-1)*y.
  #   
  #   Example:
  #     n = 1000
  #     t = complex(1, 1:n, 1)
  #     x = complex(1, runif(n), runif(n))
  #     Toe = tritoep(t)
  #     
  #     s1 = system.time( {y1 = zipsolve(t, x)} )
  #     s2 = system.time( {y2 = toesolve(t, x)} )
  #     s3 = system.time( {y3 = forwardsolve.cpp(Toe, as.matrix(x))[[1]]} )
  #     
  #     mean(abs(x - Toe %*% y1))
  #     mean(abs(x - Toe %*% y2))
  #     mean(abs(x - Toe %*% y3))
  
  n = length(t)
  
  # check if n = 2^k
  n2 = 2^ceiling(log(n) / log(2))
  
  if (n < n2) {
    # if not, then force dimension to be the power of 2 by adding n2 - n zeros
    pad = rep(0, n2 - n)
    if (is.null(y)) {
      x = zipsolve(c(t, pad))
    }
    else {
      x = zipsolve(c(t, pad), c(y, pad))
    }
    return(x[1:n])
  }
  else {
    # here we can assume that n = 2^k
    # if dimension is no larger than cutoff, then solve it directly
    if (n < cutoff) {
      x = forwardsolve(tritoep(t), as.matrix(c(1, rep(0, n - 1))))
    }
    else {
      k = cutoff
      x = rep(0, n)
      x[1:k] = forwardsolve(tritoep(t[1:k]), as.matrix(c(1, rep(0, k - 1))))
      # x[1:k] = solve(tritoep(t[1:k]), as.matrix(c(1, rep(0, k - 1))))
      if (all(Im(t) == 0)) {
        while (k < n) {
          xhat = fft(x[1:(2 * k)])
          s = Re(fft(fft(t[1:(2 * k)]) * xhat, inverse = TRUE)) / (2 * k)
          s = fft(c(s[k + 1:k], rep(0, k)))
          s = -Re(fft(xhat * s, inverse = TRUE)) / (2 * k)
          x[k + 1:k] = s[1:k]
          k = 2 * k
        }
      }
      else {
        while (k < n) {
          xhat = fft(x[1:(2 * k)])
          s = fft(fft(t[1:(2 * k)]) * xhat, inverse = TRUE) / (2 * k)
          s = fft(c(s[k + 1:k], rep(0, k)))
          s = -fft(xhat * s, inverse = TRUE) / (2 * k)
          x[k + 1:k] = s[1:k]
          k = 2 * k
        }
      }
      
    }
    
    if (is.null(y)) {
      return(x)
    }
    else {
      return(ftm(x, y))
    }
  }
}
