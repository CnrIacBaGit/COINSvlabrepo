# COINSvlabrepo
The routine implements a modelling approach for the optimal spatiotemporal control of invasive species in natural protected areas of high conservation value. 
The model  is based on diffusion equations, is spatially explicit, and includes a functional response (Holling type II) which models the control rate as a function of the state variable, i.e. the invasive species density. The control variable is represented by the effort needed to eradicate the invasive species. Furthermore a budget constraint is imposed to the amount of effort made available. 
The growth of the species is modulated by a habitat suitability function internally  computed  by using the land cover map of the study area and the map of the initial density of the invasive species. 
The habitat suitability is evaluated through  the frequency of occurrence of each land cover class in a neighborhood of  the species points of presence.   
The routine solves a constrained otpimal control problem by searching for the optimal allocation of effort which  minimizes the density of the invasive species in both time and space.
In the current version the model is applied to the Ailanthus altissima plant species populating a National Park in the South of Italy (Murgia Alta), where there is an ongoing eradication program run within a EU LIFE project


The input files, which have been assigned conventional filenames,  are:
1) a  file containing all the parameters of the model, named "parameters.csv". 
	Content an description of the paramenters in the csv file:
	D: 	 diffusivity
	r: 	 growth rate
	c: 	 penalty weight
	mu: 	 harvesting rate
	h: 	 average time spent for the abatement/eradication
	nu: 	 weight for the final population density
	omega: 	cost due to the enviromental damage
	delta: 	discount rate 
	B:     	budget 
	T:   	project duration	
	Do you want to generate the suitability map? Y (yes) or N (no)
	Do you want to run the model with control?   Y (yes) or N (no)
2) the map of initial density of the species ( TIF file), filename: "initial_density.tif"
3) the boundary of the area (SHP file), filename: "boundaryPA.shp"
4) the land cover map (SHP file), filename: "land_cover.shp".

The output files are:
5) the habitat suitability map (TIF file);
6) the raster time series of the effort allocation strategy for the control of the species along the simulation period (rts file);
7) the raster time series of the population density under control along the simulation period (rts file).

The code is implemented in R (version 3.5.0)

The development and the implementation of the model and the routine have been made possible thanks to the H2020 project ECOPOTENTIAL: Improving Future Ecosystem Benefits Through Earth Observations' (http://www.ecopotential-project.eu) which has received funding from the European Union's Horizon 2020 research and innovation programme under grant agreement No 641762. 
All the ground data regarding Ailanthus altissima (Mill.) Swingle presence and distribution are from the EU LIFE Alta Murgia Project (LIFE12 BIO/IT/000213) titled
Eradication of the invasive exotic plant species Ailanthus altissima from the Alta Murgia National Park funded by the LIFE+ Financial instrument of the European Commission.

The routine has been implemented and developed by Angela Martiradonna, Fasma Diele and Carmela Marangi, who are the owners of the IPR.
It can be used under the conditions of CC-BY-NC 2.0

The model and the specific application implemented in the routine is the result of a joint work with:
a) Baker C.M., School of Biological Sciences, The University of Queensland, Australia and CSIRO Ecosystem Sciences (for the mathematical model); 
b) Blonda P. and Tarantino C., Institute of Atmospheric Pollution Research, CNR, Bari, Italy (for land cover map and map of initial presence);
c) Casella F., Institute of Sciences of Food Production, CNR, Bari, Italy (for the in-field data and expert knowledge on the invasive species);
d) Diele F., Marangi, C., Martiradonna A., Istituto per Applicazione del Calcolo "M.Picone", CNR, Bari, Italy (for the mathematical model, data analysis, routine); 
e) Ragni S., Department of Economics and Management, University of Ferrara, Ferrara, Italy (for the mathematical model).

The description of the model is available in: 

Baker CM, Diele F, Marangi C, Martiradonna A, Ragni S. (2018) Optimal spatiotemporal effort allocation for invasive species removal incorporating a removal handling time and budget. Natural Resource Modeling; 31:e12190 https://doi.org/10.1111/nrm.12190
