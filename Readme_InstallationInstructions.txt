####################################################################################################################
####################################################################################################################
#######################################   INSTALLATION INSTRUCTIONS ################################################
####################################################################################################################
####################################################################################################################
In order to run the protein turnover scripts you need to install Python 
and all necessary packages, which are listed below. 

Please follow the installation instructions below. Preferably always chose the 64bit version if possible.

##################################################################################################################
##################################################################################################################
###### WINDOWS
##################################################################################################################
##################################################################################################################
1.) Install Python.
Use Python 2.7.x version (NOT 3.x) suitable for your operating system (preferably 64 bit, e.g. "Python 2.7.5 Windows X86-64 Installer").
http://www.python.org/download/

2.) Install Numpy for Python 2.7. At the given URL download the source (e.g. "numpy-1.8.0.zip"), unzip, and follow the installation instructions in "INSTALL.txt".
https://pypi.python.org/pypi/numpy

3.) Install Lxml for Python 2.7. At the given URL download the source (e.g. "lxml-2.3-py2.7-win-amd64.egg"), and use easy_install to install lxml (e.g. "easy_install lxml-2.3-py2.7-win-amd64.egg").
https://pypi.python.org/pypi/lxml/2.3

4.) Install Pyteomics for Python 2.7. At the given URL download and launch the windows-installer (e.g. "pyteomics-2.2.0.win-amd64_py2.7.exe").
https://pypi.python.org/pypi/pyteomics

5.) Install Pandas for Python 2.7. At the given URL download and launch the windows-installer (e.g. "pandas-0.12.0.win-amd64-py2.7.exe").
https://pypi.python.org/pypi/pandas/0.12.0

6.) Download the Protein turnover scripts ("proteinturnover_sampledata_scripts.zip") from
http://promex.pph.univie.ac.at/proteinturnoverscripts/
and place them in any location on your computer (they do not need to be installed).

7.) Install R (only necessary if you want to plot/get images of the results). At the given URL go to "Download R for Windows", chose the "base" version, download, and launch the installer (e.g. "Download R 3.0.2 for Windows").
http://cran.r-project.org/
7.1) Windows: add R to the environmental variables (this means you should be able to run R from the command line from any directory in the file system).
http://www.computerhope.com/issues/ch000549.htm

##################################################################################################################
##################################################################################################################
###### Mac OSX
##################################################################################################################
##################################################################################################################
1.) Install Python.
Use Python 2.7.x version (NOT 3.x) suitable for your operating system (preferably 64 bit, e.g. "Python 2.7.5 Mac OS X 64-bit/32-bit x86-64/i386 Installer").
http://www.python.org/download/

2.) Install Numpy for Python 2.7. At the given URL download the source (e.g. "numpy-1.8.0.zip"), unzip, and follow the installation instructions in "INSTALL.txt".
https://pypi.python.org/pypi/numpy

3.) Install Lxml for Python 2.7. At the given URL download the source (e.g. "lxml-2.3-py2.7-macosx-10.6-intel.egg"), and use easy_install to install lxml (e.g. "sudo easy_install lxml-2.3-py2.7-macosx-10.6-intel.egg").
https://pypi.python.org/pypi/lxml/2.3

4.) Install Pyteomics for Python 2.7. At the given URL download the source (e.g. "pyteomics-2.2.0.tar.gz"), unpack the archive, change to the directory and execute "sudo easy_install pip" followed by "sudo pip install pyteomics".
https://pypi.python.org/pypi/pyteomics

5.) Install Pandas for Python 2.7. At the given URL download the source (e.g. "pandas-0.12.0.tar.gz"), unpack the archive, change to the directory and execute "python setup.py build install"
https://pypi.python.org/pypi/pandas/0.12.0

6.) Download the Protein turnover scripts ("proteinturnover_sampledata_scripts.zip") from
http://promex.pph.univie.ac.at/proteinturnoverscripts/
and place them in any location on your computer (they do not need to be installed).

7.) Install R (only necessary if you want to plot/get images of the results). At the given URL go to "Download R for (Mac) OS X", download the latest version (e.g. "R-3.0.2.pkg"), and launch the installer.
http://cran.r-project.org/

##################################################################################################################
##################################################################################################################
###### Linux
##################################################################################################################
##################################################################################################################
1.) Install Python
Use any of the Python 2.7 versions (NOT 3.x) suitable for your operating system (preferably 64 bit).
http://www.python.org/

2.) Install Numpy for Python 2.7, also install "python-numpy-devel".
https://pypi.python.org/pypi/numpy

3.) Install Lxml for Python 2.7
https://pypi.python.org/pypi/lxml/2.3

4.) Install Pyteomics for Python 2.7
https://pypi.python.org/pypi/pyteomics

5.) Install Pandas for Python 2.7
https://pypi.python.org/pypi/pandas/0.12.0

6.) Download the Protein turnover scripts ("proteinturnover_sampledata_scripts.zip") from
http://promex.pph.univie.ac.at/proteinturnoverscripts/
and place them in any location on your computer (they do not need to be installed).

7.) Install R (only necessary if you want to plot/get images of the results).
http://cran.r-project.org/
####################################################################################################################
####################################################################################################################

####################################################################################################################
####################################################################################################################


The program was tested with:
Windows7 Enterprise 64bit, Intel Core i7-2600 CPU @ 3.4GHz, 16GB RAM
Mac OS X 10.9, Intel Core i5-4250U, 8GB RAM

Minimum hardware requirements:
Intel Core 2 Duo, 4GB RAM, 10MB Hard Disk Space.
(The RAM usage will vary dramatically dependant on the amount 
and size of the mzML files as well as the length of the input Selpex-list 
and the user defined rt_range value.)

EXPERIMENTFILE: 
This is the ONLY file you have to configure. Within the file you can adapt parameters e.g. the retention time 
and the mass accuracy, and set the file-paths for your mzML data files.
If you want to test the program, you will need to open the "experiment_file.txt" located in the "Scripts" folder, 
and adapt the filepath to match the current location on YOUR computer.

General Information about the Experimentfile:
Filename: absolute path to a mzML file. e.g. C:\Users\Superman\Documents\Masterplan\ABC123.mzML
TimePoint: number indicating the chronological order of the file (corresponds to the labeling). 
e.g. 1 (for very first TimePoint, minimally labeled) or 5 (for the very last TimePoint, maximally labeled, given that 5 TimePoints used).
Experiment: Indicates which files belong together in order to be evaluated as a time series. 
 any OS compatible arbitrary name can be used e.g. A1
Please do not change the order of the settings in the experiment_file.txt, meaning that "rt_range" should be 
the first, mass_accuracy the second setting etc. the last settings should be the "Filename/TimePoint/Experiment/Replicate" setting.

STARTING THE PROGRAM:
Place Experiment_template.txt in the same location as the Protein turnover scripts.
Start the program by 
A.) Navigate to the location of the protein turnover scripts (simply go to the folder you've placed the downloaded scripts).
    Double klick the "run_experimentfile.py" file to start the program.
B.) Go to the Command Line, change directory to the location of the protein turnover scripts.
    Start the program by typing "python run_experimentfile.py" (and then pressing enter).

####################################################################################################################
####################################################################################################################
