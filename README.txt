Made for interim report:
circulatory_system_components.jl		Contains funtions for the individual components of the circulatory system that can be connected to form the entire system
										Updated with docstrings for final report.
model_parameters.jl						Parameter file used for interim report. Generally the same as healthy_human_parameters.jl, but with some extra parameters for 
										windkessel investigations.
interim_report_windkessels.jl			Investigation of basic windkessel circuits.
interim_report_systemic_pulmonary.jl	Investigation of more complicated representations of the systemic and pulmonary circulations, that will be used in the final
										model of the entire cardiovascular system.
interim_report_circulatory_system.jl	Model of the entire cardiovascular system.


Model Parameters:
healthy_human_parameters.jl				Model parameters for the reference case of a healthy human
biatrial_parameters.jl					Model Parameters altered to represent human circulation with biatrial porcine transplant
bicaval_parameters.jl					Model Parameters altered to represent human circulation with bicaval porcine transplant


Investigation of biatrial porcine-to-human cardiac xenotransplantation:
total_biatrial.ipynb					Jupyter notebook comparing the total effect of biatrial cardiac xenotransplantation on some key performance indicators.
total_biatrial.ipynb					Julia script comparing the total effect of biatrial cardiac xenotransplantation on some key performance indicators. (same as above notebook)

Investigation of bicaval porcine-to-human cardiac xenotransplantation:
total_bicaval.ipynb						Jupyter notebook looking at the total effect of bicaval cardiac xenotransplantation on some key performance indicators.
total_bicaval.ipynb						Julia script looking at the total effect of bicaval cardiac xenotransplantation on some key performance indicators. (same as above notebook)


Investigation of the impact of each individual physical difference:
changes_in_isolation.ipynb				Jupyter notebook looking at the effect of each individual physical difference between pig and human hearts modelled (9 for biatrial 
										transplantation and an additional 3 for bicaval transplantation) on some key performance indicators.
changes_in_isolation.jl					Julia script looking at the effect of each individual physical difference between pig and human hearts modelled (9 for biatrial 
										transplantation and an additional 3 for bicaval transplantation) on some key performance indicators. (same as above notebook)
		
		
Note that some of these notebooks/scripts will try to save images of graphs to a folder called "plots" if a folder of this name is not present in the working folder, an
error will be thrown