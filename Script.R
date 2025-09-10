packages<-c("quantmod", "dplyr", "ggplot2", "rugarch", "stats", "moments",
            "gnorm", "tseries")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }else{
    library(pkg, character.only = TRUE)
  }
}  

#Analysis on FX rate----------------

#extract financial data for EUR/USD
eur_usd_xts <- getSymbols("EURUSD=X",
                          src = "yahoo",
                          from = "2020-09-10",
                          auto.assign = FALSE)

# Convert xts to data.frame
eur_usd_df <- data.frame(
  Date = index(eur_usd_xts),
  FX   = eur_usd_xts[, "EURUSD=X.Adjusted"]
)
colnames(eur_usd_df) <- c("Date", "FX")
#remove NA values
eur_usd_df<-na.omit(eur_usd_df)


#compute the daily returns

eur_usd_df$return <- log(eur_usd_df$FX / lag(eur_usd_df$FX, 1))
#remove the first row with NA value
eur_usd_df <- eur_usd_df[-1, ]

#squared return
eur_usd_df$r2<- eur_usd_df$return^2

#check wethere tsuqred returns are skewed
print(kurtosis(eur_usd_df$r2)) #kurtosis= 6.34>3, so we have fat tails

#from the graphs, we can infer that there is vo
#apply a simple GARCH model
p<-seq(1,4)
q<-seq(1,4)

model_selection <- data.frame(
  orders = character(16),  # empty character vector of length 16
  BIC    = numeric(16),    # empty numeric vector of length 16
  stringsAsFactors = FALSE
)
for(i in 1:length(p)){
  for (j in 1:length(q)){
    spec<-ugarchspec(
      variance.model=list(model="sGARCH", garchOrder=c(p[i], q[j])),
      mean.model = list(armaOrder=c(0,0)),
      distribution.model = "norm"
    )
    fit<-ugarchfit(
      spec=spec,
      data= eur_usd_df$r2,
  
      solver="hybrid",
    )
    bic<-infocriteria(fit)[2]
    model_selection$orders[i + (j - 1) * length(p)] <- paste0("(", p[i], ",", q[j], ")")
    model_selection$BIC[i + (j - 1) * length(p)] <- bic
  }
} #best order: (1,4)

spec<-ugarchspec(
  variance.model=list(model="sGARCH", garchOrder=c(3, 4)),
  mean.model = list(armaOrder=c(0,0)),
  distribution.model = "norm"
)
fit<-ugarchfit(
  spec=spec,
  data= eur_usd_df$r2,
  solver="hybrid",
)
#extract the standardized residuals
z <- residuals(fit, standardize = TRUE)
z<-as.numeric(na.omit(z))
#plot the model
set <- par(mfrow = c(2, 2))
plot(eur_usd_df$return, type = "l", main = "return time series", ylab = "Returns",
     col= "#60b0dc")
acf(z^2, main = "ACF of squared residuals", lag.max=20)
qqnorm(z, main = "QQ-plot of residuals vs Normal",cex= 0.5)
qqline(z, col = "red")
hist(z, breaks=40, probability=TRUE, main="Histogram of Std. Residuals", col="#60b0dc")
curve(dnorm(x), col="red", lwd=2, add=TRUE)
par(set) #reset the plot layout  

#check on the kurtosis
forth_moment <- kurtosis(z)
print(forth_moment)
#test on normality 
jb_test<- 

#plot values 
a<-par(mfrow=c(2,2))
plot(eur_usd_df$Date, eur_usd_df$return, main="EUR/USD Daily Returns",
     xlab="Date", ylab="Returns", type="l", col="blue")
plot(eur_usd_df$Date, eur_usd_df$r2, main="EUR/USD Squared Daily Returns",
     xlab="Date", ylab="Squared Returns", type="l", col="red")
plot(density(z), main="Density of Std Residuals", 
     xlab="z_t", lwd=2, col="darkgreen")
curve(dnorm(x, 0, 1), add=TRUE, col="red", lty=2)
plot(density(eur_usd_df$return), main="Scatter Plot of Returns",
     xlab="Returns", col="purple", pch=19, cex=0.5)
par(a)

#we use for proposal likelihood the generalised gaussian model
logdens_gnorm <- function(x) gnorm::dgnorm(x, mu = 0, alpha = 1, 
                                           beta = 1.5, log = TRUE)

estimate_eta <- function(residuals, eta_grid, logdens_fun) {
  # returns a data.frame with eta and the average objective Q_T(eta)
  Qvals <- sapply(eta_grid, function(e) {
    if (e <= 0) return(-Inf)
    # Q_T(eta) = mean( -log(eta) + log f(eps/eta) )
    mean( -log(e) + logdens_gnorm(residuals / e) ) # from equation 6 
  })
  out <- data.frame(eta = eta_grid, Q = Qvals)
  out$rank <- rank(-out$Q, ties.method = "first")  # 1 = best
  out
}
#compute candidate eta
eta_grid <- seq(0.005, 3, by = 0.005)


z <- as.numeric(na.omit(residuals(fit, standardize = TRUE)))

# --- run the grid search ---
eta_table <- estimate_eta(z, eta_grid, logdens)
eta_hat <- as.numeric(eta_table$eta[which.max(eta_table$Q)])

#step 3 we rescale the observation by eta_hat

y_star <- eur_usd_df$return / eta_hat

# refit GARCH(1,4) with GED likelihood
spec_ged <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1,4)),
  mean.model         = list(armaOrder = c(0,0)),
  distribution.model = "ged"
)
fit_ged <- ugarchfit(spec_ged, y_star, solver = "hybrid")

# standardized residuals from the GED fit
z_ged <- residuals(fit_ged, standardize = TRUE)

op <- par(mfrow = c(2,2))
qqnorm(z_ged, main = "QQ-plot: GED residuals") 
qqline(z_ged, col = "red")
plot(density(z_ged), main = "Density: GED residuals", xlab = "z_t", lwd = 2)
curve(dnorm(x, 0, 1), add = TRUE, lty = 2)  # N(0,1) overlay
acf(z_ged,  main = "ACF(z_t)", lag.max=20)
acf(z_ged^2, main = "ACF(z_t^2)", lag.max=20)
par(op)

# quick moment checks
if (!requireNamespace("moments", quietly = TRUE)) install.packages("moments")
sk_ged  <- moments::skewness(z_ged)
kur_ged <- moments::kurtosis(z_ged)
c(skewness = sk_ged, kurtosis = kur_ged)

# optional normality test
if (!requireNamespace("tseries", quietly = TRUE)) install.packages("tseries")
tseries::jarque.bera.test(z_ged)

#Extension


# 6.1 Choice of likelihood function---------------------------------

# we replicate the formula 23 from page 184
#for t-student distribution
A_t <- function(x, nu, eta) {
  stopifnot(nu > 2, eta > 0)
  x <- x[is.finite(x)]
  u  <- x / eta
  h1 <- (-1/eta) + ((nu + 1) * u^2) / (eta * (nu + u^2))
  h2 <- ( nu * (nu - (1 + 3*nu)*u^2 - u^4) ) / ( eta^2 * (nu + u^2)^2 )
  mean(h1^2, na.rm = TRUE) / ( eta^2 * (mean(h2, na.rm = TRUE))^2 )
}

#ratio for GED distribution
lambda_beta <- function(beta) sqrt(gamma(1/beta) / gamma(3/beta))   # Var=1 scale
c_beta      <- function(beta) beta / (lambda_beta(beta)^beta)

A_ged <- function(x, beta, eta) {
  stopifnot(beta > 0, eta > 0)
  x <- x[is.finite(x)]
  u  <- x / eta
  cb <- c_beta(beta)
  h1 <- (-1/eta) + (cb/eta) * (abs(u)^beta)
  h2 <- (1 - cb*(beta + 1) * (abs(u)^beta)) / (eta^2)
  mean(h1^2) / ( eta^2 * (mean(h2))^2 )
}
#now we compute numerically the more efficient parameter for both T and GED
#distribution
eta_t_distribution <- function(residuals, eta_grid, v_vec) {
  out <- lapply(v_vec, function(df) {
    s<-sqrt(v_vec/(v_vec-2)) #variance of t distribution
    if (df <= 2) return(NULL)                 # skip invalid dfs
    Qvals <- sapply(eta_grid, function(e) {
      if (e <= 0) return(-Inf)
      res_std<- residuals/s #standardised error
      u <- res_std / e                      # <-- your u = eps / eta
      mean(-log(e) + dt(u, df, log = TRUE))   # Eq.(6) without extra scaling
    })
    data.frame(df = df,
               eta = eta_grid,
               Q   = Qvals,
               rank = rank(-Qvals, ties.method = "first"))
  })
  do.call(rbind, out)
}

# application
eta_grid_t <- seq(0.5, 3, by = 0.05)
v <- seq(2.5, 6, by = 0.5)   # use df > 2
nu_table <- eta_t_distribution(z, eta_grid_t, v)
head(nu_table) # best choice is df=2.5, with eta=0.55

#compute parameter for general distribution

logdens_gnorm <- function(x, beta) {
  gnorm::dgnorm(x, mu = 0, alpha = 1, beta = beta, log = TRUE)
}

estimate_eta <- function(residuals, eta_grid, beta_vec) {
  out <- lapply(beta_vec, function(b){
    Qvals <- sapply(eta_grid, function(e) {
      if (e <= 0) return(-Inf)
      mean(-log(e) + logdens_gnorm(residuals / e, b))
    })
    data.frame(beta = b,
               eta  = eta_grid,
               Q    = Qvals)
  })
  do.call(rbind, out)
}

eta_grid <- seq(0.005, 3, by = 0.005)
beta_vec <- seq(0.2, 2, 0.1)

eta_table <- estimate_eta(z, eta_grid, beta_vec)
beta_table_max<- which.max(eta_table$Q)
print(eta_table[beta_table_max,]) # best beta=0.9, eta=0.81

# Now we compare with the numerical efficiency A(eta)
A_t_val <- A_t(z, nu = 2.5, eta = 0.55) # 0.3497631
a_ged_val <- A_ged(z, beta = 0.9, eta = 0.81) #1.074227


#6.2 Aggregating NGQMLE and GQMLE

#we first replicate the optimal aggragation w 
#compute expected value of h2
## inputs you picked
nu  <- 2.5
eta <- 0.55
res <- z[is.finite(z)]

## h1, h2 for Student-t (df = nu)
u  <- res / eta
h1 <- (-1/eta) + ((nu + 1) * u^2) / (eta * (nu + u^2))
h2 <- ( nu * (nu - (1 + 3*nu)*u^2 - u^4) ) / ( eta^2 * (nu + u^2)^2 )

## Eq (27): w* = E[Δ] / E[Δ^2],  Δ = (1 - ε^2)/2 - h1 / (eta * E[h2])
Eh2   <- mean(h2, na.rm = TRUE)
Delta <- (1 - res^2)/2 - h1 / (eta * Eh2)
w_star <- mean(Delta, na.rm = TRUE) / mean(Delta^2, na.rm = TRUE)
w_star

#compute the estimators for the unfeasable case when nu is known
#extract from the frist GARCH fitted model 'fit'^
theta_tilde<-coef(fit)
#adjusted estimators
#rescale dataset by our esimtated eta
y_star <- eur_usd_df$return / 0.55
#compute GARCH (1,4) model
spec_hat<-ugarchspec(
  variance.model=list(model="sGARCH", garchOrder=c(1, 4)),
  distribution.model = "std"
)
fit_hat<-ugarchfit(
  spec=spec_hat,
  data= y_star,
  solver="hybrid",
)
coef_hat<- coef(fit_hat)

#weighted coefficients
theta_w <- w_star * theta_tilde[names(coef_hat)] +
  (1 - w_star) * coef_hat[names(coef_hat)]

