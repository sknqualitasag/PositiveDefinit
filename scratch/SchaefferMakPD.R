### #
### # Positive Definite Matrices (Schaeffer, Chapter 7.4)
### # 2019-09-20 (skn)
### # ----------------------------------------------

### # G is not positive definite
G <- cbind(c(100,80,20,6),
           c(80,50,10,2),
           c(20,10,6,1),
           c(6,2,1,1))

matrixcalc::is.positive.definite(G)
eigen(G)

### # R-script to force a symmetric matrix to be positive definite
makPD = function(A){
  D = eigen(A)
  sr = 0
  nneg = 0
  V = D$values
  U = D$vectors
  N = nrow(A)
  for(k in 1:N){
    if(V[k] < 0){
      nneg = nneg + 1
      sr = sr + V[k] + V[k]
    }
  }
  wr = (sr*sr*100+1)
  p = V[N - nneg]
  for(m in 1:N){
    if(V[m] < 0){
      c = V[m]
      V[m] = p*(sr-c)*(sr-c)/wr
    }
  }
  A = U %*% diag(V) %*% t(U)
  return(A)
}

### # Apply to G
makPD(G)



### # ----------------------------------------------
### # Example 5.7.1
### # ----------------------------------------------
### # M
M <- cbind(c(100,95,80,40,40),
           c(95,100,95,80,40),
           c(80,95,100,95,80),
           c(40,80,95,100,95),
           c(40,40,80,95,100))

D <- c(399.48,98.52,23.65,-3.12,-18.52)

### # Matrix M is not positive definite so invalid as a covariance matrix
# the function is.positive.definite() give FALSE
#https://www.rdocumentation.org/packages/matrixcalc/versions/1.0-3/topics/is.positive.definite
install.packages("matrixcalc")
matrixcalc::is.positive.definite(M)

# For a positive definite matrix, the eigenvalues should be positive. This is not the case
eigen(M)

### # The procedure is to modify the negative eigenvalues such that they become positive
### # and have values between the smallest positive eigenvalue and zero.
### # The lowest eigenvalue in this example is 23.65

### # STEP 1
### # Add together the negative eigenvalues and multiply times 2
s <- (-3.12-18.52)*2
### # Square s, multiply by 100 and add 1
t <- (s*s)*100+1

### # STEP 2
n4 <- 23.65* (s+3.12) * (s+3.12)/t
n5 <- 23.65* (s+18.52) * (s+18.52)/t

### # STEP 3
### # Replace the negative eigenvalues with the new positive values, and reconstruct M
D1 <- c(399.48,98.52,23.65,0.20363,0.0774)

makPD(M)

