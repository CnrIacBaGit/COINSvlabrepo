# HabSuit.R (Habitat Suitability)
# Authors: A. Martiradonna, F. Diele, C. Marangi,  2018
# email: a.martiradonna@ba.iac.cnr.it
# In case of use of the model, the Authors should be cited.
#
# Details can be found in
# Baker C.M., F. Diele, C. Marangi, A. Martiradonna, and S. Ragni (2018).
# Optimal Control Governed by a Diusion PDE with Holling Type II Reaction Term
# and Budget Constraint. Natural Resource Modeling, in press.

presHSI <- pres_num
suppressWarnings(uds   <- readShapeSpatial(uds_path,proj4string=CRS("+proj=longlat"))  )

uds_rast <- rasterize(uds,presHSI,'DESC_ENG',fun="first")

pres_mat <- getValues(presHSI, format="matrix")
uds_mat  <- getValues(uds_rast, format="matrix")
Nrow=dim(pres_mat)[1]
Ncol=dim(pres_mat)[2]

ind_pres <- which(pres_mat>0,arr.in=TRUE)
ind_pres <- ind_pres[ind_pres[,1]>1 & ind_pres[,2]>1 & ind_pres[,1]<Nrow & ind_pres[,2]<Ncol,]
pres_mat[ind_pres]

lev <- levels(uds@data$DESC_ENG)
score=numeric(length(lev))
names(score)=lev

neigh <- matrix(c(0,-1,0,1,-1,0,1,0),ncol=2,byrow=T)

for (i in 1:nrow(ind_pres))  {

    ri <- ind_pres[i,1]
    ci <- ind_pres[i,2]

    score[uds_mat[ri,ci-1]] <- score[uds_mat[ri,ci-1]]+ ( pres_mat[ri,ci-1]==0 & !is.na(uds_mat[ri,ci-1]) )
    score[uds_mat[ri,ci+1]] <- score[uds_mat[ri,ci+1]]+ ( pres_mat[ri,ci+1]==0 & !is.na(uds_mat[ri,ci+1]) )
    score[uds_mat[ri-1,ci]] <- score[uds_mat[ri-1,ci]]+ ( pres_mat[ri-1,ci]==0 & !is.na(uds_mat[ri-1,ci]) )
    score[uds_mat[ri+1,ci]] <- score[uds_mat[ri+1,ci]]+ ( pres_mat[ri+1,ci]==0 & !is.na(uds_mat[ri+1,ci]) )

}

score_norm <- score/max(score)
sort(unique(score_norm) )
length(unique(score_norm))

rho_mat <- uds_mat

for (l in 1:length(lev))   rho_mat[rho_mat==l]<- score_norm[l]

rho_mat[is.na(rho_mat)]<- 1

rho <- uds_rast
rho[] <- rho_mat

writeRaster(rho,"data/suitability100",format="GTiff", overwrite=TRUE)