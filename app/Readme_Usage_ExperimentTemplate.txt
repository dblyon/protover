##############################################################################################################
##############################################################################################################
######################### TEMPLATE FOR EXPERIMENT_FILE.TXT ###################################################
##############################################################################################################
##### Everything after a "#" will be regarded as a comment and thus disregarded by the program. You may add as
##### many comment lines as you wish.
##### Place Experiment_template.txt in the same location as the Protein turnover scripts.
##### Start the program by 
##### A.) Navigating go the location of the Protein turnover scripts and simple double klicking the 
#####     "run_experimentfile.py" file
##### B.) Go to the Command Line, change directory to the location of the protein turnover scripts.
#####     Start the program by typing "python run_experimentfile.py" (and then pressing enter).
##### Please don't change the order of the settings in the experiment_file.txt. 
##### The settings should be in the same order as shown in this very document.
##############################################################################################################
##############################################################################################################
##############################################################################################################

##### LC/MS settings
##### Retention Time Range (positive number). Enter positive integer or float. e.g. 3
##### The rt_range serves as rt-window in which to look for the spectrum of the peptides given in the SelPEx-file. 
##### e.g. a peptide with 25min Rt (in SelPEx-file) and 3min rt_range (Experiment_file) 
##### -> program looks within 22min to 28min to find the proper spectrum.
rt_range	3

##### Mass Accuracy 
##### Enter positive integer or float. e.g. 10
##### e.g. a peptide (given in the SelPEx-file) with a theoretically calculated 400mz will be searched 
##### within 399.996mz to 400.004mz.
##### mass_accuracy_first relates only to the monoisotopic precursor.
##### mass_accuracy_rest relates to all other isotopic peaks.
mass_accuracy_first	10
mass_accuracy_rest	10

##### Input file
##### TAB delimited text file without a header, containing the following four columns in the given order: 
##### AminoAcid_Sequence (tab) Charge_State (tab) Retention_Time (tab) AccessionNumber 
##### e.g. of SelPEx-file-content: "ELVISLIVESINVIENNAK	2	23.45	ABC123"
##### e.g. of setting the path to the SelPEx-file below
selpexfilepath	C:\Users\Einstein\Documents\Selpexfile_1000Peptides.txt

##### Outputdirectory 
##### location where results-txt-file and pdf-plots will be placed as a result of the execution of the program.
outputdirectory	C:\Users\Einstein\Documents\Outputdirectory\

##### Filename/TimePoint/Experiment/Replicate settings
##### Filename: absolute path to mzML file e.g. "C:\Users\Einstein\Documents\Masterplan\ABC123.mzML"
##### TimePoint: provide a positive number (including 0) indicating the chronological order of the files (corresponding to the labeling status). 
##### e.g. 1 (for very first TimePoint, minimally labeled) or 5 (for the very last TimePoint, maximally labeled, given that 5 TimePoints used).
##### Experiment: Indicates which files belong together in order to be evaluated as a time series. 
##### e.g. "A1" or "Worldpeace_experiment1" (String without spaces or special characters, any arbitrary name can be used). 
##### Replicate: provide a positive number or string (without special characters or whitespaces) indicating technical replicates
##### e.g. "1"
Filename	TimePoint	Experiment	Replicate
C:\Users\Einstein\Documents\abc123.mzML	0	ControlA	1
C:\Users\Einstein\Documents\abc234.mzML	24	ControlA	1
C:\Users\Einstein\Documents\abc345.mzML	48	ControlA	1
C:\Users\Einstein\Documents\abc456.mzML	72	ControlA	1
C:\Users\Einstein\Documents\abc567.mzML	96	ControlA	1
C:\Users\Einstein\Documents\def123.mzML	0	ControlA	2
C:\Users\Einstein\Documents\def234.mzML	24	ControlA	2
C:\Users\Einstein\Documents\def345.mzML	48	ControlA	2
C:\Users\Einstein\Documents\def456.mzML	72	ControlA	2
C:\Users\Einstein\Documents\def567.mzML	96	ControlA	2
