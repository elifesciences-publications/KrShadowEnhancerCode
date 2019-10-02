# KrShadowEnhancerCode
Code for analyzing transcriptional traces and generating figures in Waymack, et al., 
Before running any of this code you first need to run the code for tracking the individual spots of transcription
For each movie:
Run the scripts of the LivemRNA package, available on Github at https://github.com/GarciaLab/mRNADynamics with instructions for these scripts at Image processing Nikon instructions (https://docs.google.com/document/d/1t3_gyt9EjffZVtf5o6vI93UdrGLoCPojmTlul4FEL0g/edit?usp=sharing). 
Fill out the DataStatus excel sheet for each movie that you produce and want analyzed. Create a separate sheet for each construct, naming that sheet the name of the construct. Each column from B on should be filled in for a separate movie. The Prefix row and CompileParticles rows must be filled in for that movie’s data to be used in any of the following scripts. Fill in Prefix row with ‘Prefix=’ followed by the name you have given that movie. Ex: for the movie named ‘2018-03-5-Kr1_Kr1’ ,the Prefix row would be filled out as: Prefix=’2018-03-05-Kr1_Kr1’    To have an embryo’s data used in analysis, the CompiledParticles row must be filled out “ApproveAll’ 

Running "RW" scripts for compiling and ploting data as in Waymack, et al:
1.	Run “RunComparingSpotCorrAdjRW” to create SpotDiff.m files for all approved movies – The most important thing this function does is use the ElapsedTime of your movie, the nuclear tracking information (‘XXX_lin.mat’ and APDetection.mat) , and the spot tracking information (CompiledParticles) to determine for each time frame whether a given nucleus exists in the movie or not and if so, whether or not it has one or two transcriptional spots (and the corresponding fluorescence values for those spots at that time point). 
a.	You will first need to edit the ConstructList at the beginning of the code to list the names of your constructs as you have them in DataStatus.xls
b.	You can also perform the underlying function ComparingSpotsCorrRWAdj for each individual movie- to do so you need to load the CompiledParticles.m, APDetection.m, and ‘XXX_lin.mat’ files from that movie’s folder, and manually run the save function at the end of the code. 

2.	Run “RunBurstPropertiesSlope.m” to create BurstPropertiesSlope.mat files for all approved movies. This is the code that goes through all of the fluorescence traces to determine where transcriptional bursts are happening and records the different properties of each burst (i.e. burst size, frequency, duration, etc.) 
a.	Here you must also edit the ConstructList at the beginning 
b.	If you are using a different microscope, you will also want to edit the ON and OFF thresholds in SlopeBurstCalling.m (lines 26-27) to correspond to the FRNAP (or negative FRNAP for OFF threshold) calculated for your microscope as done in Lammers, et al., 2018

Depending on what you are interested in, you may want to run all or only a subset of the following scripts. In each, the calculations section is the first section and the code to make graphs of the data follows to allow for easier manipulation of different graphs/colors/etc. that you may want to produce without having to change the calculations performed each time. For all of them, you will need to first edit the ConstructList at the beginning of the code to list your constructs as you have them in DataStatus.xlsx

Calculating correlation of allele activity in dual reporter embryos – AlleleCorrelationCalcRW
1.	This script takes the allele correlations calculated for each individual nucleus by the SpotCorrelationAdj.m script (which is performed on each approved embryo when RunComparingSpotCorrAdjRW is run) and groups by AP bin, embryo, and enhancer construct. The output is a structure, AvgCorrAllAP, that contains all of the individual correlation values, as well as the average, standard deviation, standard error, and 95% confidence intervals, organized by AP bin. 
2.	The very end of the calculations section removes data points from the final structure where there was data from only 1 allele at that AP bin for the entire construct. This was done to avoid the issue of trying to calculate error bars with only one data point. If you’re not planning on imaging many embryos (>= 3) you might want to remove this section. 
3.	The remaining sections are used for graphing the allele correlation data and therefore will likely needed to be modified for individual graphing needs. 

Comparing total mRNA produced by different constructs - AvgmRNAProdbyEmbryoRW
1.	This groups all of the integrated fluorescence values associated with each allele (recorded in CompiledParticles) by AP bin, embryo, and construct. The output is a structure, AvgProdAllAP, that contains a few different things, the most useful being an array of all integrated fluorescence values of a given construct (ie across all the movies you have of that construct) organized by AP bin (each column is a different AP bin). It also uses this data to calculate the 95% CI (as well as underlying SD and SE) for each construct at each AP bin. 
2.	The very end of the calculations section removes data points from the final structure where there was only 1 allele at that AP bin for the entire construct. This was done to avoid the issue of trying to calculate error bars with only one data point. If you’re not planning on imaging many embryos (>= 3) you might want to remove this section. 
3.	The first subsection of the plotting section sets the conversion factor (MS2Conversion) for listing fluorescence values as #’s of mRNA. If you are not imaging on the Wunderlich Nikon lab with the listed laser settings, you will need to change this value to reflect your own calculated F1 value.

Comparing burst size between constructs – AvgBurstSizeRW
1.	Groups all of the burst sizes of each transcriptional locus (ie there will often be more than one per allele if that allele had more than 1 transcriptional burst during the movie) by AP bin, embryo, and construct. The output, AvgAmpAllAP, is a structure very similar to that for the total mRNA production, but providing the burst size data. 
2.	This script also removes data points from the final structure that contain only 1 data point for an entire AP bin.    
3.	The first subsection of the plotting section sets the conversion factor (MS2Conversion) for listing fluorescence values as #’s of mRNA. If you are not imaging on the Wunderlich Nikon lab with the listed laser settings, you will need to change this value to reflect your own calculated F1 value.

Comparing burst frequency between constructs – AvgBurstFrequencyCalcRW
1.	Like above, groups all the burst frequencies from each transcriptional spot (only one value per allele) by AP bin, embryo, and construct. Creates structure, AvgFreqAllAP, that has fields containing all of the Frequency values grouped by AP bin, as well as average for a whole construct by AP bin, SD, SE, 95% CI, etc. 
2.	Note this uses the definition of frequency as the number of bursts happening between the first transcriptional burst and the end of nc14 divided by that length of time 

Comparing burst duration between constructs – AvgBurstDurationCalcRW
1.	Groups burst durations from each transcriptional spot (can be more than one per allele – ie as many bursts as that allele had) by AP bin, embryo, and construct. Creates structure, AvgDurAllAP, containing all duration values grouped by AP bin, as well as average for construct by AP bin, the SD, SE, 95% CI, etc. 

Calculating total noise, covariance, and inter-allele noise – TotalNoiseCalcRW
-	This goes through each embryo for all of your listed constructs and in each nucleus that has 2 active alleles (ie both have at least three active time points) calculates the total noise, covariance, and inter-allele noise for that nucleus. Nuclei are stored according to AP bin and combined for all of the embryos of a given construct. 
-	The WholeNoise structure that is produced also contains a field for the integrated fluorescence values of each transcriptional spot to enable determination of the AP bin of highest expression for each construct. Graphs then use the noise values from the AP bin of maximum average expression for each construct. 

Calculating coefficient of variation for transcription across time of nc14 -  CVnc14TimeRW
-	Calculates the mean and std of the fluorescence of each transcriptional spot across the time of nc14 (incorporating the SpotDiff data so we know if the nucleus exists and the allele is silent vs the nucleus just does not exist at that time point). Does this for all transcriptional spots in all of the embryos of the constructs you list in ConstructList at the beginning. Organizes this data in a structure so you can access all of the CV values grouped by AP bin for each different construct. 

The code used for producing and graphing the modeling data are anneal.zip, Simulator_main.m, main.m, readerEnergies.m, plotter.m, boxplot2.m, and multiple_boxplot.zip files  
