

setwd( "C:/Users/James.Thorson/Desktop/Project_git/TMB_experiments/SPDE vs 2D_AR1" )

library(TMB)
library(RandomFields)
library(INLA) # FROM: http://www.r-inla.org/download

###################
# Simulate data
###################

Dim = c("n_x"=10, "n_y"=10)
loc_xy = expand.grid("x"=1:Dim['n_x'], "y"=1:Dim['n_y'])
Scale = 4
Sigma2 = (0.5) ^ 2
beta0 = 3
prob_missing = 0.2

# Simulate spatial process
RMmodel = RMgauss(var=Sigma2, scale=Scale)
epsilon_xy = array(RFsimulate(model=RMmodel, x=loc_xy[,'x'], y=loc_xy[,'y'])@data[,1], dim=Dim)
image( z=epsilon_xy )

# SImulate counts
c_xy = array(NA, dim=dim(epsilon_xy))
for(x in 1:nrow(c_xy)){
for(y in 1:ncol(c_xy)){
  c_xy[x,y] = rpois(1, exp(beta0 + epsilon_xy[x,y]) )
  if( rbinom(n=1, size=1, prob=prob_missing)==1) c_xy[x,y] = NA
}}

# Helper functions
rarray = function(dim) array(rnorm(prod(dim)), dim=dim)

# Range: 0.1=rho^t -> t = log(0.1)/log(rho)

# create mesh
mesh = inla.mesh.create( loc_xy )
# Create matrices in INLA
spde <- inla.spde2.matern(mesh)

# Create 2D AR1 elements
dist_grid = dist(loc_xy, diag=TRUE, upper=TRUE)
grid_size_km = sort(unique(dist_grid))[1]
M0 = as( ifelse(as.matrix(dist_grid)==0, 1, 0), "dgTMatrix" )
M1 = as( ifelse(as.matrix(dist_grid)==grid_size_km, 1, 0), "dgTMatrix" )
M2 = as( ifelse(as.matrix(dist_grid)==sqrt(2)*grid_size_km, 1, 0), "dgTMatrix" )

# Compile
Params = list( "beta0"=0, "ln_sigma2"=0, "logit_rho"=0, "epsilon_xy"=array(rnorm(prod(dim(loc_xy))),dim=dim(epsilon_xy)) )
compile( "comparison.cpp" )
dyn.load( dynlib("comparison") )
Data = list( "Options_z"=0, "c_i"=as.vector(c_xy), "j_i"=mesh$idx$loc-1, "spde"=spde$param.inla[c('M0','M1','M2')], "M0"=M0, "M1"=M1, "M2"=M2 )

######## SPDE
Data[["Options_z"]] = 0
Params = list( "beta0"=0, "ln_tau"=0, "ln_kappa"=0, "epsilon_j"=rarray(mesh$n) )
# Build object
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_j", DLL="comparison" )
# Optimize
Opt_spde = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
Report_spde = Obj$report()

######### Unequal distance 2D autoregressive
Data[["Options_z"]] = 1
Params = list( "beta0"=0, "ln_tau"=0, "ln_kappa"=2, "epsilon_j"=rarray(nrow(loc_xy)) )
# Build object
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_j", DLL="comparison" )
# Optimize
Opt_ar = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, upper=ifelse(names(Obj$par)=="ln_kappa",0,Inf) )
Report_ar = Obj$report()

# Comparison
unlist(Report_spde[c("Range","MargSD")])
unlist(Report_ar[c("Range","MargSD")])

