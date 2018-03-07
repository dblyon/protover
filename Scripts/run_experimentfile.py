from __future__ import print_function
import cameltoolkit, camelherder, py2r, os, peptide_results, filterresults, sys
from locale import localeconv
from time import clock

class Run(object):

    def __init__(self):
        pass      
        
    def set_platform(self, platform):
        self.platform = platform
    
    def get_platform(self):
        return self.platform
    
    def set_decimalpoint(self, decimalpoint):
        self.decimalpoint = decimalpoint
    
    def get_decimalpoint(self):
        return self.decimalpoint
        
    def run_proteinturnover(self, exp_dict_oneexperiment, plottitle, 
                            selpexfilepath, rt, mass_accuracy_first, mass_accuracy_rest, 
                            outputdir, dbl=False, escapee=True):
        fn_results_settings = os.path.basename(selpexfilepath)
        fn_results_settings = fn_results_settings.replace(".txt", "")
        filename_ancres = ("%s"+"%s"+"%s"+"%s"+"%s"+"%s") % (outputdir, "temp/", fn_results_settings, "_", plottitle, "_RESULTS.txt")
        ch.herdthem_mzml(exp_dict_oneexperiment, selpexfilepath, sr,
                            rt_range = rt,
                            mz_ppm_picking_distance_mip0 = mass_accuracy_first,
                            mz_ppm_picking_distance = mass_accuracy_rest,
                            escapee=escapee, dbl=dbl, filename_ancres=filename_ancres)
    
if __name__ == "__main__":
    start = clock()
    sr = peptide_results.Selpex_results()
    protr = peptide_results.Protein_results()
    ch = camelherder.Camelherder()
    run = Run()
    decimalpoint = localeconv()["decimal_point"]
    run.set_decimalpoint(decimalpoint)
    platform = sys.platform
    if platform == "darwin":
        run.set_platform("unix")
    elif "win" in platform:
        run.set_platform("win")
    else:
        run.set_platform("unix")

################################################################################
################################################################################
################################################################################
    dbl = False
    check_input = True
    escapee = True
    cleanTP0 = True
    RIAincreasing = True
################################################################################
################################################################################
################################################################################
    
    def _check_userinput_absorrelpath(fn):
        """
        FileName -> absolute path FileName
        """
        if os.path.isabs(fn):
            return fn
        else:
            cwd = os.getcwd()
            pdir =  os.path.abspath(os.path.join(cwd, os.pardir))
            return os.path.join(pdir, fn)
    
    def _check_userinput_experimentfile(fn):
        try:
            fh = open(fn, "r")
            fh.close()
            return fn
        except:
            print(fn)
            print("'experiment_file.txt' could not be opened. Compare your file to the provided template, it should be tab-delimited and located in the same folder as the protein turnover scripts.")
            sys.exit(0)
    
    def _check_userinput_number(number, setting_name):
        try:
            number = float(number.strip())            
        except:
            print("You provided:")
            print(number)
            str2print = ("%s"+"%s"+"%s") % ("Please provide a positive number for ", setting_name, " setting in the 'experiment_file.txt' (compare your file to the provided template, it should be tab-delimited).")
            print(str2print)
            sys.exit(0)
        if number < 0.0:
            print("You provided:")
            print(number)
            str2print = ("%s"+"%s"+"%s") % ("Please provide a positive number for ", setting_name, " setting in the 'experiment_file.txt' (compare your file to the provided template, it should be tab-delimited).")
            print(str2print)
            sys.exit(0)
        return number

    def _check_userinput_file(fn, setting_name, check_input):
        if not check_input:
            return fn
        try:
            fh = open(fn, "r")
            fh.close()
            return fn
        except:
            str2print = ("%s"+"%s") % (setting_name, " in 'experiment_file.txt' expected a file.")
            print(str2print)
            print(fn)
            print(" could not be opened. Compare your file to the provided template, it should be tab-delimited.")
            sys.exit(0)
    
    def _check_userinput_createdir(directory, setting_name):
        try:
            if not os.path.exists(directory):
                os.makedirs(directory)
            return directory
        except:
            str2print = ("%s"+"%s") % (setting_name, " in 'experiment_file.txt' expected a directory.")
            print(str2print)
            print(directory)
            print("could neither be opened nor created. (Compare your file to the provided template, it should be tab-delimited.)")
            sys.exit(0)

    def experimentfile2runturnover(fn, sep="\t"):
        """
        #exp_dict: key = Experiment, val = {key = TimePoint, val = Filename}                
        #exp_dict[Experiment] = {key = TimePoint, val = fn }
        #ancres_dict: key = Experiment, val = {key = Replicate, val = Ancres_fn}
        #ancres_dict: key = Experiment, val = {key = Replicate, val = {key = 0/1/2 (=nofilter/filter12/filter1234), val = Ancres_fn/Ancres_fn/...}}
        """
        fh = open(fn, "r")
        exp_dict = {}
        ancres_dict = {}
        settings_dict = {}
        skip = False
        for line in fh:
            line = line.strip()
            line = line.replace(",", ".")
            if line == "":
                continue
            if line[0] == "#":
                continue
            else:
                line = line.split("#")[0].strip()
                if len(line.split(" ")) > 1:
                    print("Please use tabs not spaces as delimiters in the 'experiment_file.txt'")
                    sys.exit(0)
                linesplit = line.split(sep)
                if skip == False:
                    if linesplit != ["Filename", "TimePoint", "Experiment", "Replicate"]:
                        skip = False
                        try:
                            x = linesplit[1].strip()
                        except IndexError:
                            print("Error: Expected character(s) after ", x, "in the 'experiment_file.txt', please compare your file to the provided template. The file should be tab-delimited (not space).")
                            sys.exit(0)                            
                        if linesplit[0].strip() == "rt_range":
                            rt_range = _check_userinput_number(linesplit[1], "rt_range")
                            settings_dict["rt"] = rt_range
                        elif linesplit[0].strip() == "mass_accuracy_first":
                            mass_accuracy = _check_userinput_number(linesplit[1], "mass_accuracy_first")
                            settings_dict["mass_accuracy_first"] = mass_accuracy
                        elif linesplit[0].strip() == "mass_accuracy_rest":
                            mass_accuracy = _check_userinput_number(linesplit[1], "mass_accuracy_rest")
                            settings_dict["mass_accuracy_rest"] = mass_accuracy
                        elif linesplit[0].strip() == "inputfilepath":
                            fn_selpex = _check_userinput_absorrelpath(linesplit[1].strip())
                            selpexfilepath = _check_userinput_file(fn_selpex, "inputfilepath", check_input)
                            settings_dict["selpexfilepath"] = selpexfilepath
                        elif linesplit[0].strip() == "outputdirectory":
                            outputdir = linesplit[1].strip()
                            outputdir = _check_userinput_absorrelpath(outputdir)                            
                            outputdir = _check_userinput_createdir(outputdir, "outputdirectory")                        
                            settings_dict["outputdir"] = outputdir
                    else:
                        skip = True
                else:
                    if len(linesplit) != 4:
                        print("You've provided:")
                        print(line)
                        print("Please provide 4 tab delimited entries for 'File/Experiment settings' in the 'experiment_file.txt'")
                        sys.exit(0)
                    filename = linesplit[0].strip()
                    filename = _check_userinput_absorrelpath(filename)                            
                    filename = _check_userinput_file(filename, "Filename", check_input)
                    timepoint = linesplit[1].strip()
                    timepoint = _check_userinput_number(timepoint, "TimePoint")
                    experiment = linesplit[2].strip()
                    if not experiment:
                        print("Please provide an 'Experiment' name in the 'experiment_file.txt' (tab-delimited). You've provided: ")
                        print(experiment)
                        sys.exit(0)
                    replicate = linesplit[3].strip()
                    if not replicate:
                        print("Please provide a 'Replicate' name/number in the 'experiment_file.txt' (tab-delimited). You've provided: ")
                        print(replicate)
                        sys.exit(0)
                    experiment_replicate = experiment + "_" + replicate
                    if not exp_dict.has_key(experiment_replicate):
                        exp_dict[experiment_replicate] = {}
                        exp_dict[experiment_replicate][timepoint] = filename
                    else:
                        exp_dict[experiment_replicate][timepoint] = filename
                    try:
                        selpexfilepath
                    except:
                        print("Please provide a 'selpexfilepath' in the 'experiment_file.txt'")
                        sys.exit(0)
                    try:
                        outputdir
                    except:
                        print("Please provide an 'outputdirectory' in the 'experiment_file.txt'")
                        sys.exit(0)
                    fn_results_settings = os.path.basename(selpexfilepath)
                    fn_results_settings = fn_results_settings.replace(".txt", "")        
                    plottitle = experiment_replicate
                    filename_ancres = ("%s"+"%s"+"%s"+"%s") % (fn_results_settings, "_", plottitle, "_RESULTS.txt")
                    ancres_fn_long = ("%s"+"%s") % ("temp/", filename_ancres)
                    if not ancres_dict.has_key(experiment):
                        ancres_dict[experiment] = {}
                        ancres_dict[experiment][replicate] = {}
                        ancres_dict[experiment][replicate][0] = ancres_fn_long
                    else:
                        ancres_dict[experiment][replicate] = {}
                        ancres_dict[experiment][replicate][0] = ancres_fn_long
        return(settings_dict, exp_dict, ancres_dict)
    fn = os.getcwd()
    fn = fn + "/experiment_file.txt" 
    settings_dict, exp_dict, ancres_dict = experimentfile2runturnover(_check_userinput_experimentfile(fn))
    
    if len(settings_dict.keys()) != 5:        
        num_missing_settings = 5 - len(settings_dict.keys())
        str2print = ("%s"+"%s"+"%s") % ("There are ", num_missing_settings, " missing or badly set settings in the 'experiment_file.txt', please compare your file to the provided template. The file should be tab-delimited.")
        print(str2print)
        sys.exit(0)
    print("\n\n", "All settings accepted and files accessible. Starting to run protein turnover.")
    
    selpexfilepath = settings_dict["selpexfilepath"]
    rt = settings_dict["rt"]
    mass_accuracy_first = settings_dict["mass_accuracy_first"]
    mass_accuracy_rest = settings_dict["mass_accuracy_rest"]   
    
    exp_dict_keys_sorted = sorted(exp_dict.keys())
    for experiment_replicate in exp_dict_keys_sorted:
        outputdir = settings_dict["outputdir"]
        plottitle = experiment_replicate
        dir_list = []
        outputdir = outputdir + plottitle + "/"
        dir_list.append(outputdir)
        if dbl:
            histcoveragetemplate = ("%s"+"%s") % (outputdir, "Histos/")
            dir_list.append(histcoveragetemplate)
        temp_dir = ("%s"+"%s") % (outputdir, "temp/")
        dir_list.append(temp_dir)
        twoDstacks = ("%s"+"%s") % (outputdir, "2Dstacks/")
        dir_list.append(twoDstacks)
        ria = ("%s"+"%s") % (outputdir, "RIA/")
        dir_list.append(ria)
        for directory in dir_list:
                if not os.path.exists(directory):
                    os.makedirs(directory)        
        run.run_proteinturnover(exp_dict[experiment_replicate], plottitle, selpexfilepath, rt, mass_accuracy_first, mass_accuracy_rest, outputdir, dbl=dbl, escapee=escapee)
    
    ################## Post processing
    outputdir = settings_dict["outputdir"].replace("\\", "/")
    fr = filterresults.Filter(ancres_dict, outputdir)
    fr.create_rplotbatchtxt()
    plot = py2r.Talkr()
    
    ################## clean TP0
    if cleanTP0:
        fr.cleanTP0()
    ################## calc RIA
    fr.calcRIA_relint()
    
    print("\n", "Plotting 2Dstacks", "\n")
    ################## 2Dstacks PLOT
    fr.plot_2Dstacks(plot)
        
    ################# PLOT TP vs RIA (before Filter 3 = post-processing-filter)
    if dbl:
        fr.replot_tpvsria(plot, dbl=False)

    ################## FILTERS
    ################## Filter #3 --> Peptide Level # use first since fast and easy
    if RIAincreasing:
        fr.ancres_filterriaincreasing()
        
    ################# PLOT TP vs RIA (after Filter 3)
    print("\n", "Plotting RIAs", "\n")
    fr.replot_tpvsria(plot, dbl=dbl)
    
    ################# PLOT Histos
    if dbl:
        bins = (0, 1.05, 0.05)
        bins = fr._make_bins(bins)
        fr.replot_histos(plot, bins, plotfrequency=True, plotdensity=True, dbl=dbl)
    
    ################## BATCH R PLOTS
    r_plot_batch = fr.run_rplotbatch(plot, run.get_platform())
    
    ################# move 2Dstacks to positive and negative folder
    fr.move_2dstacks()
    
    ################## Cleanup files
    if r_plot_batch:
        if not dbl:
            fr.cleanup_files(run.get_decimalpoint())
    
    stop = clock()
    print("runtime[min]: ", (stop - start)/60.0 )
    print("\n\n", "DONE", "\n")

        
    
    