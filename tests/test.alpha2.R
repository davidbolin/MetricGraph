library(GPGraph)
loc <- c(0,0.5,1,1.5)
#loc <- seq(from=0,to=1,length.out=10)
kappa=10
tau = 1
sigma=1

Q.list <- Q.alpha2.line(loc,kappa,tau)
Q <- solve(Q.list$Sigma)
Q.partial <- Q.list$Q[c(Q.list$ind.proc,Q.list$ind.der),c(Q.list$ind.proc,Q.list$ind.der)]

print(Q)
print(Q.partial)
1.0*abs(Q.partial-Q)>1e-8
pattern <- Q
pattern[abs(pattern)<1e-10]=0
pattern[abs(pattern)>1e-10]=1


l = 0.5
x <- loc[1:2]
C <- GPGraph:::build.C.beta1(0.5, kappa=kappa, sigma=sigma, nu=3/2)
S1 <- GPGraph:::matern.neumann.free2(x, x, C, kappa=kappa, sigma=sigma, nu=3/2, L = l)
S2 <- GPGraph:::matern.neumann.free2(x, x, C, kappa=kappa, sigma=sigma, nu=3/2, L = l, deriv = c(0,1))
S3 <- GPGraph:::matern.neumann.free2(x, x, C, kappa=kappa, sigma=sigma, nu=3/2, L = l, deriv = c(1,1))
Sigma  <- rbind(cbind(S1, S2), cbind(t(S2),S3))
Q.true <- solve(Sigma)

T = kappa*l
const <- 2*kappa*tau^2/((1-2*T^2)^2 - 2*exp(2*T)*(2*T^2+1)+exp(4*T))
q1 <- exp(4*T) - (1-2*T^2)^2 + 4*T*exp(2*T)
q2 <- 4*T*exp(2*T)
q3 <- 2*exp(T)*(2*T^2*(T-1)-T-exp(2*T)*(T+1)+1)
q4 <- 2*T*exp(T)*(2*T^2-exp(2*T)-1)
q5 <- -kappa^2*(1-2*T^2)^2 + exp(2*T)*(2*kappa^2+4*(kappa^2-1)*T^2 - 4*T - 2) - (kappa^2-2)*exp(4*T)
q6 <- 2*exp(T)*(-2*T^3 - 2*T^2 + T + exp(2*T)*(T-1)+1)
Q <- -const*matrix(c(kappa^2*q1,kappa*q2,kappa^2*q3,kappa*q4,
                kappa*q2,q5,kappa*q4,q6,
                kappa^2*q3, kappa*q4, kappa^2*q1, kappa*q2,
                kappa*q4, q6, kappa*q2, q5),4,4)
Q <- Q[c(1,3,2,4),c(1,3,2,4)]
print(Q.true)
print(Q)
