  __  __       _ _   _                                   _               ____                _ __  __      
 |  \/  |_   _| | |_(_)_ __  _ __ ___   ___ ___  ___ ___(_)_ __   __ _  |  _ \ ___  __ _  __| |  \/  | ___ 
 | |\/| | | | | | __| | '_ \| '__/ _ \ / __/ _ \/ __/ __| | '_ \ / _` | | |_) / _ \/ _` |/ _` | |\/| |/ _ \
 | |  | | |_| | | |_| | |_) | | | (_) | (_|  __/\__ \__ \ | | | | (_| | |  _ <  __/ (_| | (_| | |  | |  __/
 |_|  |_|\__,_|_|\__|_| .__/|_|  \___/ \___\___||___/___/_|_| |_|\__, | |_| \_\___|\__,_|\__,_|_|  |_|\___|
                      |_|                                        |___/                                     

Currently this code runs well on my personal laptop:
	16 GB RAM, Intel Core i7-10875H CPU @ 2.30GHz, 2304 Mhz, 8 Cores, 16 Logical Processors

To start (see Release Notes for more information):

Install multiprocessing and shutil packages for python

In InputParameters.txt in /Phonon-Explorer-master/Input_Files/
	Set BkgMode and QMode
	Set path for sqw file
	Set Energy steps and q binning
	Set MinPointsInDataFile. 
		I choose 40 to have a good number of points.
	Set maxY and maxX for output pdf files so they are easy to scroll through.
		Make a few initial cuts to see what a good y range is.
	Define range in h,k,l
	Set InitialAmplitude
	Set Lattice Constants
	Set maxFiles. 
		Too few and background will be bad. 
		Too many and it takes a long time to generate and fit later on.
	Set Resolution. 
		I find 4 works well for ARCS data, but may be different for others. 
		Too small and the background subtraction is too much.
	Set MinPointsInDataBackgroundFile
	Set WidthLowerBound
	Set fileWithGuesses (just file name, not path)
	
In MP_GenerateConstQCuts_New.py in /Phonon-Explorer-master/Python_Code
	Set root_path. 
		This is the path up to the output directory.
	Set input_path. 
		This is where the previously mentioned InputParameters.txt is located.
	Set run. 
		This is an experiment indicator such as a temperature or pressure.
	Example:
		.../Documents/TiO2/327K/GM/0.2/
		.../.....root...../run./GM/0.2/
		
	In CopyInput() function, edit direction names and reduced q's in if statements.
	
	In "if __name__ == '__main__':" section
		Set AllDirections and direction1, 2, 3, etc. 
			I keep them in groups of 3 to not use too much memory. 
			You may need fewer if you want to run more q's at a time.
		Set q values in for loops.

Run MP_GenerateConstQCuts_New.py
	This generates cuts and performs the background subtractions.
	It's best to run overnight if on a laptop or desktop.
	Time heavily depends on how many q points there are in the sample and how many background files are generated.
	
In MP_MultiFit.py,
	Change root_run to approproate folder (up to run folder). 
		.../Documents/TiO2/327K/GM/0.2/
		.../.....root_run....../GM/0.2/
	Change directions and q's. I can run all of them at once without issue.
	Generate 'Initial Guess' files based on pdf output and place each text file in the appropriate /{q}/ folder.

Run MP_MultiFit.py
	See fitting output in {q}/constant_q/subtr_background/xconstant_q_allPlots.pdf and _FittingParam.txt
	Adjust and rerun as needed.
	
Finished!
	
I find that some fits fail or do not attempt to fit well. 
Sometimes all of the "fitted" peaks are of very similar widths, which is inaccurate.
	In this case, try changing some initial fitting parameters such as the starting amplitude in InputParameters.txt.

Other changes:
	Set upper bound on peak width in FittingData.py by changing last number in "UB[3*ii+2]=self.WidthLowerBound*4"
		This is based on the WidthLowerBound value in InputParameters.txt
	Allow position fitting more or less leniancy by changing value in PositionDelta=0.1*self.positionGuesses[ii]. 
		Default is 0.1.