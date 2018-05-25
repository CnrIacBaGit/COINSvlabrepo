# COINS.R COntrol of INvasive Species
# Authors: A. Martiradonna, F. Diele, C. Marangi,  2018
# email: a.martiradonna@ba.iac.cnr.it
# In case of use of the model, the Authors should be cited.
#
# Details can be found in
# Baker C.M., F. Diele, C. Marangi, A. Martiradonna, and S. Ragni (2018).
# Optimal Control Governed by a Diusion PDE with Holling Type II Reaction Term
# and Budget Constraint. Natural Resource Modeling, in press.

rm(list=ls()); graphics.off()

# Load packages
library(rgeos)
library(rgdal)
library(raster)
library(rts)
library(maptools)   # for shp files
library(Matrix)     # for sparse matrices
library(pracma)     # for repmat
library(ggplot2)
library(plyr)

# Settings
presen_path <- "data/initial_density.tif"
domain_path <- "data/boundaryPA.shp"
param_path  <- "data/parameters.csv"
uds_path    <- "data/land_cover.shp"


# Read inputs files
param    <- read.csv(param_path,sep=';',header=FALSE)
pres     <- raster(presen_path)
print('Tiff file read')
suppressWarnings( domain   <- readShapeSpatial(domain_path,proj4string=CRS("+proj=longlat")) )

# Parameters of the model
par_values <- as.character(param$V2)
D      <- as.numeric(par_values[1])
r      <- as.numeric(par_values[2])
c      <- as.numeric(par_values[3])
mu     <- as.numeric(par_values[4])
h      <- as.numeric(par_values[5])
nu     <- as.numeric(par_values[6])
omega  <- as.numeric(par_values[7])
delta  <- as.numeric(par_values[8])
B      <- as.numeric(par_values[9])
T      <- as.numeric(par_values[10])
k <-1; q<- 2;
m <- 2*q-1; alpha <- c*m/(B^m)

# Numerical grid
dx   <- 100
dxkm <- dx*1e-3

print("Creating numerical grid")
rasterOptions(chunksize=1e+4, maxmemory=1e+4, progress='text', overwrite=TRUE)
pres_num <- raster::aggregate(pres, fact=dx/xres(pres), fun=max, extend=FALSE)

# Suitability

if (par_values[11]=="N") {HSI <- 1} else          {
   print("Creating suitability map")
   source("HabSuit.R")
   HSI_rast <- rho
   HSI <- getValues(HSI_rast,format="matrix")
   HSI <- t(HSI)
   HSI <- HSI[,ncol(HSI):1]
   HSI <- c(HSI)
}

# Tolerance settings
tolFB    <- 1e-5
tolsolve <- 1e-15
itmax    <- 500

# Select grid points where contol will be applied
control <- par_values[12]
temp_rast <- rasterize(domain,pres_num)
temp_mat <- getValues(temp_rast,format="matrix")
temp_mat <- t(temp_mat)
temp_mat <- temp_mat[,ncol(temp_mat):1]
temp_vec <- c(temp_mat)
if (control == "Y") {
          ind_in  <- !is.na(temp_vec)
          name_flag <- "with_control"
} else {
          ind_in  <- rep(FALSE,length(temp_vec))
          name_flag <- "without_control"
          itmax=1
}

# Time settings

date0 <- as.Date("2012/7/1")

dt <- 0.05
int <- 1/dt

# Calculate some values
Nx <- dim(pres_num)[2]
Ny <- dim(pres_num)[1]
n <- Nx*Ny
Nt <- round(T/dt) +1
t <- ((0:(Nt-1)))*dt
elo <- exp(-delta*dt)


### Assemby matrices
#print("Assemby matrices")
# Matrix L (without 1/dxkm^2 factor)
L = bandSparse(n,n, k=c(0,Nx),diagonals=list(4+numeric(n),-1+numeric(Nx*(Ny-1))), symmetric = T)
# S
L[1,1]=3; L[1,2]=-3/2; L[Nx,Nx]=6; L[Nx,Nx-1]=-3
L[cbind(2:(Nx-1),1:(Nx-2))]=-1;
L[cbind(2:(Nx-1),3:Nx)]=-1;
#Z
L[(Ny-1)*Nx+1,(Ny-1)*Nx+1]=6; L[n,n]=3
L[(Ny-1)*Nx+1,(Ny-1)*Nx+2]=-3; L[n,n-1]=-3/2
L[cbind(((Ny-1)*Nx+2):(Nx*Ny-1),((Ny-1)*Nx+1):(n-2))]=-1
L[cbind(((Ny-1)*Nx+2):(Nx*Ny-1),((Ny-1)*Nx+3):(n))]=-1
#T
L[1,Nx+1]=-3/2; L[Nx,2*Nx]=-3
L[cbind(2:(Nx-1),(Nx+2):(2*Nx-1))]=-2
#Y
L[(Ny-1)*Nx+1,(Ny-2)*Nx+1]=-3; L[n,(Ny-1)*Nx]=-3/2
L[cbind( ((Ny-1)*Nx+2):(n-1), ((Ny-2)*Nx+2):((Ny-1)*Nx-1))]=-2
#X
ind1 <- ind2 <- matrix(nrow=0,ncol=2)
for (j in 1:(Ny-2)) {
    ind1 <- rbind(ind1, cbind((j*Nx+2):((j+1)*Nx-1),(j*Nx+1):((j+1)*Nx-2)),cbind((j*Nx+2):((j+1)*Nx-1),(j*Nx+3):((j+1)*Nx)) )
    ind2 <- rbind(ind2, c(j*Nx+1,j*Nx+2), c((j+1)*Nx,(j+1)*Nx-1))
}
L[ind1]=-1; L[ind2]=-2

# Matrix B
Bmat <- bandSparse(n,n, k=0,diagonals=list(1+numeric(n))) + dt/(dxkm^2)*D*L

# Initial conditions
U0 <- getValues(pres_num,format="matrix")
U0 <- t(U0)
U0 <- U0[,ncol(U0):1]
u0 <- c(U0)

vT <- nu*exp(delta*T) #+numeric(n)

# Time stepping procedure

print("Forward-backward")

u <- v <- Evec <- numeric(n)
uout <- uold <- repmat(matrix(u0,nrow=n),1,Nt)
vout <- vold <- repmat(matrix(vT,nrow=n),1,Nt)
vout[,Nt] <- vT; vold[,Nt] <- vT
Eout <- matrix(0,nrow=n,ncol=Nt)
vl1  <- repmat(matrix(vT,nrow=n),1,Nt-1)

STOP=0; it=0;

while (STOP==0)   {
    it=it+1;

    # Forward
    u = u0

    for (nt in 2:Nt) {

        v=vl1[,nt-1];
        Evec[ind_in]=(q/2/alpha*( sqrt( 1+4*alpha*mu*u[ind_in]*v[ind_in]/q^2./(1+h*mu*u[ind_in]) ) -1))^(1/(q-1))

        Fun =  r*HSI*u-r*u^2/k - mu*u*Evec/(1+h*mu*u)
        u = u + dt*Fun

        u <- solve(Bmat,u,sparse=TRUE,tol=tolsolve)

        if (min(u)<0) {
           disp(paste('u has negative elements, nt=',nt,'it=',it))
        }

        uout[,nt] <- as.numeric(u)
    }

    # Backward
    v = vout[,Nt]
    for (nt in seq(Nt-1,1,by=-1)) {

        u=uout[,nt]
        v <- solve(Bmat,v,sparse=TRUE,tol=tolsolve)

        if (min(v)<0) {
           disp('v has negative elements')
        }

        vl1[,nt]=as.numeric(v)
        Evec[ind_in]=(q/2/alpha*(sqrt(1+4*alpha*mu*u[ind_in]*v[ind_in]/q^2./(1+h*mu*u[ind_in]))-1))^(1/(q-1))

        G= omega+r*HSI*v-2*r/k*u*v-mu*v*Evec/(1+h*mu*u)^2
        v =elo*(v+dt*G)
        vout[,nt]=as.numeric(v)
    }

    Eout[ind_in,]=( q/2/alpha*(sqrt(1+4*alpha*mu*uout[ind_in,]*vout[ind_in,]/q^2/(1+h*mu*uout[ind_in,]))-1) )^(1/(q-1))

    err_u=norm(uold-uout,type="1"); norm_u=norm(uout,type="1")
    err_v=norm(vold-vout,type="1"); norm_v=norm(vout,type="1")
    STOP= min(tolFB*norm_u-err_u,tolFB*norm_v-err_v) >0 | it==itmax
    print( paste('it=', it, '  rel.err=', max(err_u/norm_u,err_v/norm_v)) )

    uold=uout
    vold=vout

}

### Output

dates <- seq(date0, by = "year", length.out = T+1)

ab_rast <- pres_num
for (s in seq((int+1),Nt,by=int)) {
    us <- uout[,s]
    Us <- matrix(us,ncol=Ny)
    Us <- Us[,ncol(Us):1]
    Us <- t(Us)
    dens <- setValues(pres_num, Us, format='matrix')
    ab_rast <- stack(ab_rast, dens)
}
ab_rts <- rts(ab_rast, dates)
write.rts(ab_rts,paste("data/AILANTHUS_density_",T,"y",sep=""),overwrite=TRUE)

e0 <- Eout[,1]
E0 <- matrix(e0,ncol=Ny)
E0 <- E0[,ncol(E0):1]
E0 <- t(E0)
eff_rast <- setValues(pres_num, E0, format='matrix')
for (s in seq((int+1),Nt,by=int)) {
    e0 <- Eout[,s]
    E0 <- matrix(e0,ncol=Ny)
    E0 <- E0[,ncol(E0):1]
    E0 <- t(E0)
    eff <- setValues(pres_num, E0, format='matrix')
    eff_rast <- stack(eff_rast, eff)
}
eff_rts <- rts(eff_rast,dates)
write.rts(eff_rts,paste("data/AILANTHUS_effort_",T,"y",sep=""),overwrite=TRUE)
