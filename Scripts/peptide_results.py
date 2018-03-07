from __future__ import print_function
import cameltoolkit, mfbt, py2r #camelherder
#from time import clock

class input_results(object):
    def __init__(self):
        """
        expects nothing.
        creates a dictionary: key = aaseq, vals = peptide result
        """
        self.input_results = {}   
    
    def set_r_plottemptxt(self, outputdir):
        fn = ("%s"+"%s"+"%s") % (outputdir, "PlotTXTForR", ".txt")
        fh = open(fn, "w")
        fh.close()
        self.r_plottemptxt = fn
    
    def get_r_plottemptxt(self):
        return self.r_plottemptxt        
        
    def plot_subprocess_file2r(self, fn):
        """
        File -> None
        send content of file line by line to R via subprocess to plot
        """
        plot = py2r.Talkr()
        plot.init_subprocess()
        fh = open(fn, "r")
        for line in fh:
            line = line.lstrip()
            plot.ph.stdin.write(line)
        fh.close()
        plot.ph.stdin.write("""q()""")
        plot.kill_subprocess()          
        
    def plotit(self, expmips2D = False, cruderange2D = False, expmips3D = False, cruderange3D = False, expmipsVsTP = False, plottitle=None, plotfilename=None, plotfilelocation_2D=None, plotfilelocation_3D=None, plotfilelocation_expmipsvstps=None):
        """
        """
        plot = py2r.Talkr()
        plot.set_r_plottemptxt(self.get_r_plottemptxt())
        if plotfilename:
            plotfilename = plotfilename
        else:
            plotfilename = ""
        home = cameltoolkit.get_home_forR()
        if not plotfilelocation_2D:
           plotfilelocation_2D = ("%s"+"%s") % (home, """/Downloads/plots/2Dstacks/""")
        else:   
            plotfilelocation_2D = plotfilelocation_2D
        if not plotfilelocation_3D:
            plotfilelocation_3D = ("%s"+"%s") % (home, """/Downloads/plots/3Dplots/""")
        else:
            plotfilelocation_3D = plotfilelocation_3D
        if not plotfilelocation_expmipsvstps:
            plotfilelocation_expmipsvstps = ("%s"+"%s") % (home, """/Downloads/plots/expmipsvstps/""")
        else:
           plotfilelocation_expmipsvstps = plotfilelocation_expmipsvstps
        for aaseq in self.input_results:
            pep = self.input_results[aaseq].pep
            if expmips2D:
                plot.plot_pep2camelstacks(self, pep, plotfilelocation_2D, cruderange=False, plottitle=plottitle, plotfilename=plotfilename)
            if cruderange2D:
                plot.plot_pep2camelstacks(self, pep, plotfilelocation_2D, cruderange=True, plottitle=plottitle, plotfilename=plotfilename)
            if expmips3D:
                plot.plot_pep3Dplot(self, pep, plotfilelocation_3D, cruderange=False, plottitle=plottitle, plotfilename=plotfilename)
            if cruderange3D:
                plot.plot_pep3Dplot(self, pep, plotfilelocation_3D, cruderange=True, plottitle=plottitle, plotfilename=plotfilename)
            if expmipsVsTP:
                plot.plot_pep2numexpmipsvstp(self, aaseq, plotfilelocation_expmipsvstps, plottitle=plottitle, plotfilename=plotfilename)
        
    def calc_rias_forallpeps(self):
        """
        Self -> None
        calculate the RIA (q-ratio) for all peptides within self.input_results
        """
        for aaseq in self.input_results:
            self.input_results[aaseq].calc_rias()
                
    def get_coveragetemplate_allpeps(self, tp):
        """
        Self -> List
        produce list with numbers (sticks found / sticks template) for all peptides of given TP.
        """
        coveragetemplate_allpeps_list = []
        for aaseq in self.input_results:
            coveragetemplate_allpeps_list.append(self.input_results[aaseq].get_coverage_template(tp))
        return coveragetemplate_allpeps_list
        
    def get_coveragetheomips_allpeps(self, tp):
        coveragetheomips_allpeps_list = []
        for aaseq in self.input_results:
            coveragetheomips_allpeps_list.append(self.input_results[aaseq].get_coverage(tp))
        return coveragetheomips_allpeps_list

    def get_riaincreasing_bool_list(self):
        """
        Self -> List
        produce list of booleans, True if RIA increases for all TPs for a peptide, False otherwise.
        """
        riaincreasing_bool_list = []
        for aaseq in self.input_results:
            pep = self.input_results[aaseq]
            riaincreasing_bool_list.append(pep.ria_increasing())
        return riaincreasing_bool_list

    def plot_riaincreasing_allpeps_frequency(self, plottitle=None, plotfilename=None, plotfilelocation=None):
        """
        Self -> None
        produce plot containing information of all results
        if RIA increases for all TPs (from TPmin to TPmax) -> True, otherwise False.
        """
        plot = py2r.Talkr()
        plot.set_r_plottemptxt(self.get_r_plottemptxt())
        home = cameltoolkit.get_home_forR()
        if not plottitle:
            plottitle = "Frequency of RIA increasing for all TPs of Peptide"
        if not plotfilename:
            plotfilename = "RIA_increasing_frequency"   
        if not plotfilelocation:
            plotfilelocation = ("%s"+"%s") % (home, """/Downloads/plots/histcoveragetemplate/""")
        riaincreasing_bool_list = self.get_riaincreasing_bool_list()
        num_true  = riaincreasing_bool_list.count(True)
        num_false = riaincreasing_bool_list.count(False)
        xlab1 = "NO"
        xlab2 = "YES"
        ylab = "Frequency of RIA increasing for all TPs of Peptide"
        plot.plot_riaincreasing_allpeps_frequency(plotfilelocation, plottitle, plotfilename, num_true, num_false, xlab1, xlab2, ylab)
           
    def plot_riaincreasing_allpeps_density(self, plottitle=None, plotfilename=None, plotfilelocation=None):
        """
        Self -> None
        produce plot containing information of all results
        if RIA increases for all TPs (from TPmin to TPmax) -> True, otherwise False.
        plot the 
        """
        plot = py2r.Talkr()
        plot.set_r_plottemptxt(self.get_r_plottemptxt())
        home = cameltoolkit.get_home_forR()
        if not plottitle:
            plottitle = "Density of RIA increasing for all TPs of Peptide"
        else:
            plottitle = plottitle
        if not plotfilename:
            plotfilename = "RIA_increasing_density"
        else:
            plotfilename = plotfilename
        if not plotfilelocation:
            plotfilelocation = ("%s"+"%s") % (home, """/Downloads/plots/histcoveragetemplate/""")
        else:
            plotfilelocation=plotfilelocation
        riaincreasing_bool_list = self.get_riaincreasing_bool_list()
        num_true_abs  = riaincreasing_bool_list.count(True)
        num_false_abs = riaincreasing_bool_list.count(False)
        num_true_perc  = (float(num_true_abs)  / len(riaincreasing_bool_list)) * 100
        num_false_perc = (float(num_false_abs) / len(riaincreasing_bool_list)) * 100        
        ylab = "Density of RIA increasing for all TPs of Peptide"
        plot.plot_riaincreasing_allpeps_density(plotfilelocation, plottitle, plotfilename, num_true_abs, num_false_abs, num_true_perc, num_false_perc, ylab)
         
    def plot_hist_coveragetemplate(self, bins=None, plotfileloc_histo_freq=None, plotfileloc_histo_dens=None, plotfrequency=False, plotdensity=False, coverage_template=None):
        """
        Self, String -> None 
        produce histogram plots for each given TP
        x-axis: 0 to 100% (bin-size = 1), y-axis: number of occurrences within given results
        How many peptides 
        """
        plot = py2r.Talkr()
        plot.set_r_plottemptxt(self.get_r_plottemptxt())
        if not bins:
            bins = range(0, 101, 1)
            bins = [x/100.0 for x in bins]
            changename = False
        else:
            (start, stop, interval) = bins
            start = int(start*1000)
            stop = int(stop*1000)
            interval= int(interval*1000)
            bins = range(start, stop, interval)
            bins = [x/1000.0 for x in bins]
            changename = True
        home = cameltoolkit.get_home_forR()
        if not plotfileloc_histo_freq:
            plotfileloc_histo_freq = ("%s"+"%s") % (home, """/Downloads/plots/histcoveragetemplate/""")
        else:   
            plotfileloc_histo_freq = plotfileloc_histo_freq
        if not plotfileloc_histo_dens:
            plotfileloc_histo_dens = ("%s"+"%s") % (home, """/Downloads/plots/histcoveragetemplate/""")
        else:   
            plotfileloc_histo_dens = plotfileloc_histo_dens
        x_axis_label = "Coverage [%]"
        y_axis_label_freq = "Frequency [#] of Peptides"
        y_axis_label_dens = "Density [%] of Peptides"
        tp_fn_list = self.input_results[self.input_results.keys()[0]].get_anchor_results_keys_sorted()
        for tp in tp_fn_list:
            plottitle = ("%s%s") % ("Coverage of template   ", tp)
            if not changename:
                plotfilename = tp.replace(".mzML", "")
            else:
                plotfilename_base = tp.replace(".mzML", "_")
                plotfilename_ext = "bins"
                plotfilename = ("%s%s%s") % (plotfilename_base, "_",  plotfilename_ext)
            if coverage_template:
                fractionfound_list = self.get_coveragetemplate_allpeps(tp)
                plotfilename += "_CovTemp_"
            else:
                fractionfound_list = self.get_coveragetheomips_allpeps(tp)
                plotfilename += "_CovTheo_"
            if plotfrequency:
                plot.setplotfilelocation(plotfileloc_histo_freq)
                plotfilename = ("%s%s") % ("HistoFreq_", plotfilename)
                plot.plot_hist_coveragetemplate(bins, plotfilename, plottitle, x_axis_label, y_axis_label_freq, y_axis_label_dens, fractionfound_list, plotfrequency=plotfrequency, plotdensity=plotdensity)
            if plotdensity:
                plot.setplotfilelocation(plotfileloc_histo_dens)
                plotfilename = ("%s%s") % ("HistoDens_", plotfilename)
                plot.plot_hist_coveragetemplate(bins, plotfilename, plottitle, x_axis_label, y_axis_label_freq, y_axis_label_dens, fractionfound_list, plotfrequency=plotfrequency, plotdensity=plotdensity)
      
class Protein_results(object):    
    def __init__(self):
        self.protein_results = {}    
    
    def inputres2protres(self, input_results_object):
        """
        from input_results_object.input_results create entries in protein_results_object.prot_results,
        key = accession number, value = dictionary -->
        key = AAseq, value = a list, each element being a list 
        with TP[h] and the corresponding q value. in chronological order,
        from min to max labeled.
        """
        for aaseq in input_results_object.input_results:
            for an in input_results_object.input_results[aaseq].get_an():
                if self.protein_results.has_key(an):
                    self.protein_results[an][aaseq] = input_results_object.input_results[aaseq].get_tp_q_pairs()
                else:
                    self.protein_results[an] = {aaseq: input_results_object.input_results[aaseq].get_tp_q_pairs()}
                    
    def plot_rias(self, input_results_object, plotfilelocation = None, plottitle=None):
        """
        """
        if plotfilelocation:
            plotfilelocation_rias = plotfilelocation
        else:
            home = cameltoolkit.get_home_forR()
            plotfilelocation_rias = ("%s"+"%s") % (home, """/Downloads/plots/RIA/""")
        self.prepare_plot_rias(input_results_object, plotfilelocation = plotfilelocation_rias, plottitle=plottitle)
        
    def prepare_plot_rias(self, input_results_object, plotfilelocation=None, plottitle=None):
        """
        """
        plot = py2r.Talkr()
        sr = input_results_object
        plot.set_r_plottemptxt(sr.get_r_plottemptxt())
        for an in self.protein_results:
            if plottitle:
                plottitle_2 = ("%s"+"%s"+"%s") % (an, "_", plottitle)
            else:
                plottitle_2 = an
            aaseq_list = []
            tp_q_list = []
            for aaseq in self.protein_results[an]:
                aaseq_list.append(aaseq)
                tp_q_list.append(self.protein_results[an][aaseq])
            tp_q_list_clean = []
            if tp_q_list:
                for ele in tp_q_list:
                    if isinstance(ele[0][1], float):
                        tp_q_list_clean.append(ele)
            if tp_q_list_clean:
                if plotfilelocation: 
                    plot.setplotfilelocation(plotfilelocation)
                plot.plot_prot2rias(aaseq_list, tp_q_list_clean, plottitle_2)
        
class IterPeptide_results(object):    
    def __init__(self, keys, sorted_anchor_results):
        self.sorted_anchor_results = sorted_anchor_results
        self.keys = keys
        self.filebasename = 0
        self.length = len(sorted_anchor_results)
        
    def __next__(self):
        if self.filebasename < self.length:
            self.filebasename += 1
            return self.sorted_anchor_results[self.keys[self.filebasename - 1]]
        else:
            raise StopIteration
        
    def next(self):
        return self.__next__()
    
class Peptide_results(object):
    def __init__(self, aaseq, charge, rt, pep):
        """
        aaseq, charge and rt correspond to input input,
        chonoslist short (not absolute paths)
        """
        self.aaseq = aaseq
        self.charge = charge
        self.rt = rt 
        self.pep = pep 
        self.anchor_results_dict = {}
        self.an_dict = {}
        self.temp_anchor_results_dict = {}

    def __iter__(self):
        return IterPeptide_results(self.get_anchor_results_keys_sorted(), self.anchor_results_dict)    

    def set_tpchronosreverse_dict(self, tpchronosreverse_dict):
        """
        {1: "0_10D15NLa1_15.mzML',
         2: "24_11D15NLa1_11.mzML',
         3: "48_12D15NLa1_12.mzML',
         4: "72_13D15NLa1_21.mzML',
         5: "14D15NLa1_20.mzML'}
        # key = tp (TimePoint of Experiment File), val = filename (from 1 for TPmax to n for TPmin)
        """
        self.tpchronosreverse_dict = tpchronosreverse_dict
    
    def get_tpchronosreverse_dict(self):
        return self.tpchronosreverse_dict
    
    def get_tps_chronological(self):
        l = self.get_tpchronosreverse_dict().keys()
        return sorted(l)
    
    def get_tpfilebasename_list_chronological(self):
        key_chronos = self.get_tps_chronological()
        tpfilebasename_list_chronological = []
        d = self.get_tpchronosreverse_dict()
        for tp in key_chronos:
            filebasename = d[tp]
            tpfilebasename_list_chronological.append((tp, filebasename))
        return tpfilebasename_list_chronological
    
    def add_anchor_result(self, anchor_result, filebasename):
        """
        adds an entry to anchor_results_dict
        key = filebasename
        value = anchor_result
        index = position of fn in reverse chronological order,
            starting with 0 (for TPmax) to e.g. 4 (for TP0 with 5TPs)
        """
        if not self.anchor_results_dict.has_key(filebasename):
            self.anchor_results_dict[filebasename] = anchor_result
        else:
            if anchor_result.get_total_score() > self.anchor_results_dict[filebasename].get_total_score():
                self.anchor_results_dict[filebasename] = anchor_result
    
    def set_remcamelscoresbyrt(self, filebasename, remcamelscoresbyrt):
        self.anchor_results_dict[filebasename].set_remcamelscoresbyrt(remcamelscoresbyrt)
    
    def get_anchor_results_keys_sorted(self):
        """
        None -> List
        gets keys from self.anchor_results_keys (= FileBaseName) and sorts them chronologically,
		from min to max labeled TimePoints,
        returns sorted keys of ancres dict
        TPmin to TPmax
        """
        dict = self.get_tpchronosreverse_dict()
        keys_sorted = sorted(dict.keys())
        anchor_results_keys_sorted = []
        for key in keys_sorted:
            anchor_results_keys_sorted.append(dict[key])
        return anchor_results_keys_sorted

    def get_later_tp(self, filebasename):
        """
        Self, String -> String
        produce TP chronologically after the given TP (more labeled), if given TP last in list produce False
        """
        tp_fn_list = self.get_anchor_results_keys_sorted()
        for index, tp in enumerate(tp_fn_list):
            if tp == filebasename:
                if (index+1) == len(tp_fn_list):
                    return False
                else:
                    return tp_fn_list[index+1]

    def get_earlier_tp(self, filebasename):
        """
        Self, String -> String
        produce TP chronologically before the given TP (less labeled), if given TP first in list produce False
        """
        tp_fn_list = self.get_anchor_results_keys_sorted()
        for index, tp in enumerate(tp_fn_list):
            if tp == filebasename:
                if index == 0: 
                    return False
                else:
                    return tp_fn_list[index-1]      
            
    def get_anchor_result_ofTP(self, filebasename):
        """
        expects the basename of a rawfile, returns a dictionary
        with Rts as keys and anchor results as vals
        """
        return self.anchor_results_dict[filebasename]  
    
    def get_tp_escapee_list(self):
        """
        Self -> List
        produce nested list [filebasename, escapee_list] = [String, List]
        """
        tp_escapee_list = []
        tp_chronos_list = self.get_anchor_results_keys_sorted()
        for filebasename in tp_chronos_list:
            escapee_list = self.anchor_results_dict[filebasename].get_escapee_list()
            tp_escapee_list.append([filebasename, escapee_list])
        return tp_escapee_list
    
    def set_rawdata_frame(self, rawdata_frame):
        """
        (rt_low, rt_high, mz_low, mz_high, centroid)
        add attribute to object
        """
        self.rawdata_frame = rawdata_frame
    
    def get_rawdata_frame(self):
        """
        returns self.rawdata_frame
        """
        return self.rawdata_frame
    
    def get_rtboundaries(self):
        """
        returns a tuple with rt_low and rt_high
        """
        return (self.rawdata_frame[0], self.rawdata_frame[1])
   
    def get_aaseq(self):
        return self.aaseq
    
    def get_charge(self):
        return self.charge
    
    def get_rt(self):
        """
        This Retention Time originates from the input input,
        NOT from the spectrum of the expmips and camelscore are from.
        """
        return float(self.rt)
    
    def get_leadcamel(self):
        """
        returns anchor result of TPmax.
        """
        return self.anchor_results_dict[self.get_anchor_results_keys_sorted()[-1]]
        
    def get_theomips(self):
        return self.pep.gettheomips()
	
    def get_expmips(self, filebasename):
        return self.anchor_results_dict[filebasename].get_expmips()
    
    def get_coverage(self, filebasename):
        return self.anchor_results_dict[filebasename].get_coverage()
        
    def get_coverage_template(self, filebasename):
        return self.anchor_results_dict[filebasename].get_coverage_template()
    
    def get_coverage_theomips(self, filebasename):
        return self.anchor_results_dict[filebasename].get_coverage()
    
    def get_total_score(self, filebasename):
        return self.anchor_results_dict[filebasename].get_total_score()
    
    def get_ria(self, filebasename):
        return self.anchor_results_dict[filebasename].get_ria()
  
    def add_temp_anchor_result(self, anchor_result, filebasename):
        """
        expects an Anchor_result instance and a filename of a time point,
        adds it to self.temp_anchor_results_dict, filebasename is the key and
        another dictionary the value. The second dict has Rt as key 
        and Anchor_result as value.
        returns None
        """
        if not self.temp_anchor_results_dict.has_key(filebasename):
            self.temp_anchor_results_dict[filebasename] = {anchor_result.get_rt(): anchor_result}
        else:
            self.temp_anchor_results_dict[filebasename][anchor_result.get_rt()] = anchor_result
    
    def get_temp_anchor_results_camelscores_byrt(self, filebasename):
        camelscores_byrt = {}
        for rt in self.temp_anchor_results_dict[filebasename]:
            camelscores_byrt[rt] = self.temp_anchor_results_dict[filebasename][rt].get_camel_score()
        return camelscores_byrt
    
    def get_bestanchorresult(self, filebasename):
        """
        return the anchor result object with the hightest total score.
        that lives within self.temp_anchor_results_dict
        """
        # rt_bestancres = self.get_rttotalscorelist_sorttempancresbytotalscore(filebasename)[0][0]
        rt_bestancres = self.get_rttotalscorelist_sorttempancresbytotalscore(filebasename)
        if rt_bestancres:
            return self.temp_anchor_results_dict[filebasename][rt_bestancres[0][0]]
        else:
            return self._get_negative_default_ancres()
    
    def _get_negative_default_ancres(self):
        """
        produce Anchor_result-object with negative default values.
        """
        expmips = []
        camel_score = {"coverage": 0.0, "coverage_template": 0.0, "ppm_rms": 0.0, "logmip0intensity": 0.0, "weighted_sumppm": 0.0, "rt": 0.0, "total_score": -1}
        escapee_list = []
        raw_data_cruderange = []
        # return Anchor_result(expmips, camel_score, escapee_list, raw_data_cruderange)
        return Anchor_result(expmips, camel_score, False)
    
    def get_remcamelscoresbyrt(self, filebasename):
        """
        returns a dictionary. key = rt. value = camelscore.
        of all camelscores analyzed (except camelscore of maxscore ancres).
        """
        remcamelscoresbyrt = {}
        # rt_bestancres = self.get_rttotalscorelist_sorttempancresbytotalscore(filebasename)[0][0]
        rt_bestancres = self.get_rttotalscorelist_sorttempancresbytotalscore(filebasename)
        if rt_bestancres:
            del self.temp_anchor_results_dict[filebasename][rt_bestancres[0][0]]
            for rt in self.temp_anchor_results_dict[filebasename]:
                remcamelscoresbyrt[rt] = self.temp_anchor_results_dict[filebasename][rt].get_camel_score()
            return remcamelscoresbyrt
        else:
            return remcamelscoresbyrt
        
    def set_bestanchorresult(self, filebasename):
        """
        replaces dict containing all anchor results with one single anchor result of best score.
        set best result as the ONLY result --> value of
        self.anchor_results_dict[filebasename]
        """
        rt_bestancres = self.get_rttotalscorelist_sorttempancresbytotalscore(filebasename)[0][0]
        bestancres = self.temp_anchor_results_dict[filebasename][rt_bestancres]
        self.add_anchor_result(bestancres, filebasename)
    
    def clear_temp_anchor_results_dict(self):
        self.temp_anchor_results_dict.clear()

    def clear_anchor_results_dict(self):
        self.anchor_results_dict.clear()
        
    def get_rttotalscorelist_sorttempancresbytotalscore(self, filebasename):
        """
        for each TP of self.temp_anchor_results_dict
        get the total scores of self.temp_anchor_results_dict[filebasename]
        sort them by score
        return list with Rts and totalscores sorted by scores (max to min).
        """
        rt_tscore_list = []
        if self.temp_anchor_results_dict.has_key(filebasename):
            pass
        else:
            return None
        for rt in self.temp_anchor_results_dict[filebasename]:
            rt_tscore = [rt, self.temp_anchor_results_dict[filebasename][rt].get_total_score()]
            rt_tscore_list.append(rt_tscore)
        return sorted(rt_tscore_list, key=lambda x: x[1], reverse=True)
    
    def calc_rias(self):
        """
        calculates q ratio according to "Gustavsson et al., Proteomics,  2005"
        q = 15N / (14N + 15N)
        define which x15N are 14N by looking at TP min labeled,
        rest are assumed to be 15N.
        """
        expmips_TPmin = self.anchor_results_dict[self.get_anchor_results_keys_sorted()[0]].get_expmips()
        x14N_positions = [x[-1] for x in expmips_TPmin]
        for tp_fn in self.anchor_results_dict:
            expmips = self.anchor_results_dict[tp_fn].get_expmips()	
            sum_x14N = 0.0
            sum_x15N = 0.0
            for mzintppmx15N in expmips:
                if mzintppmx15N[-1] in x14N_positions:
                    sum_x14N+= mzintppmx15N[1]
                else: sum_x15N+= mzintppmx15N[1]
                q = (sum_x15N / (sum_x14N + sum_x15N)) # !!! shouldn't be indented?
                self.anchor_results_dict[tp_fn].set_ria(q) # !!! shouldn't be indented?
                
    def set_an(self, an_list):
        """
        expects a list with an as elements,
        for each an adds entry to self.an_dict.
        return None.
        """
        for an in an_list:
            if not self.an_dict.has_key(an):
                self.an_dict[an] = ""
    
    def get_an(self):
        """
        returns a list, each element being an accession number.
        """
        return sorted(self.an_dict.keys())
     
    def get_tp_q_pairs(self):
        """
        Self -> List
        produces a nested list of [TP, RIA] pairs in chronologically ascending order
        returns a list, each element being a list 
        with TP[h] and the corresponding q value.
        """
        tp_ria_list = []
        tp_chronos_list = self.get_tpchronosreverse_dict().keys()
        tp_chronos_list = sorted(tp_chronos_list)
        for tp in tp_chronos_list:
            fn_short = self.get_tpchronosreverse_dict()[tp]
            ria = self.anchor_results_dict[fn_short].get_ria()
            tp_ria_list.append([tp, ria])
        return tp_ria_list
                          
    def get_ppmdeviations(self, theomips = True):
        """
        calculates the ppm deviation of each MIP to template,
        if theomips is True --> theoretical mz vals used as template
        if False --> leadcamel used as template
        returns a list with tuples, 
        each entry consisting of TP (basename of file) and list with ppm_vals
        """
        ppm_results_list = []
        if theomips == True:
            template = self.get_theomips()
        else:
            leadcamel_key = self.get_anchor_results_keys_sorted()[-1]
            template = self.anchor_results_dict[leadcamel_key].get_expmips()
        keyssorted = self.get_anchor_results_keys_sorted()
        for key in keyssorted:
            ppm_list = []
            expmips = self.anchor_results_dict[key].get_expmips()
            for ele in expmips:
                expmz = ele[0]
                x15N = ele[-1]
                for x in template:
                    if x[-1] == x15N:
                        ppm = mfbt.Ppmprecision(x[0], expmz).getppmprecision()
                        ppm_list.append(ppm)
            ppm_results_list.append((key, ppm_list))
        return ppm_results_list
    
    def get_fractionfound(self, filebasename):
        """
        Self, filebasename -> Number
        produce the number of sticks found divided by the number of sticks in the template
        TPmax:   sticks_found / theo_mips
        TPmax-1: sticks_found / sticks_found_TPmax
        """
        filebasename_template = self.get_later_tp(filebasename)
        if not filebasename_template:
            num_sticks_template = len(self.get_theomips())
        else:
            num_sticks_template = len(self.get_expmips(filebasename_template)) 
        num_sticks_current =  len(self.get_expmips(filebasename))
        try:
            return float(num_sticks_current) / num_sticks_template
        except ZeroDivisionError:
            num_sticks_template =  self.get_theomips()
            return float(num_sticks_current) / num_sticks_template
        
    def ria_increasing(self):
        """
        Self -> Boolean
        produce True if RIA (q-ratio) of previousTP <= followingTP for all TPs, False otherwise
        """
        tpria_list = self.get_tp_q_pairs()
        ria_chronological_list = [x[1] for x in tpria_list]
        for index, ele in enumerate(ria_chronological_list):
            if index+1 == len(ria_chronological_list):
                break
            if ele <= ria_chronological_list[index+1]:
                continue
            else:
                return False
        return True
    
class Anchor_result(object):
    
    def __init__(self, expmips, camel_score, escapee_utilized):
        self.expmips = expmips
        self.camel_score = camel_score
        self.ria = -1
        self.escapee_utilized = escapee_utilized
    
    def get_escapee_utilized(self):
        return self.escapee_utilized
        
    def set_remcamelscoresbyrt(self, remcamelscoresbyrt):
        """
        expects a dictionary. key = rt. value = camel_score (a dict).
        remaining camelscores (excluding the best/maximum score) of other spectra.
        """
        self.remcamelscoresbyrt = remcamelscoresbyrt
    
    def get_remcamelscoresbyrt(self):
        return self.remcamelscoresbyrt
    
    def get_camelscores_byrt(self):
        """
        dict of all other/remaining spectra
        key = rt
        value = camel score
        """
        return self.camelscores_byrt
        
    def get_camelscores_byrt_rt(self, rt):
        return self.camelscores_byrt[rt]        
            
    def set_remaining_camel_scores(self, camel_scores_sorted):
        """
        add anc_scores_sorted (list with camel_scores, sorted descending)
        """
        self.remaining_camel_scores = camel_scores_sorted
    
    def get_remaining_camel_scores(self):
        """
        [ppm_rms, intensity_measure, total_score]
        """
        return self.remaining_camel_scores
    
    def set_expmips(self, expmips):
        self.expmips = expmips
    
    def get_expmips(self):
        return self.expmips
    
    def get_norm_expmips(self):
        return cameltoolkit.normtomax(self.expmips)    

    def get_camel_score(self):
        """
        calculation of total_score below
        LeadCamel_total_score = (coverage - (weighted_sumppm/100.0)) - (1 / logmip0intensity)
        CamelTrail_total_score = (logmip0intensity - weighted_sumppm + coverage - penalty)
        total_score is a dictionary with the following keys:
        rt, total_score, coverage, ppm_rms, weighted_sumppm, weighted_sumintensities_log, logmip0intensity
        """
        return self.camel_score

    def get_total_score(self):
        return self.camel_score["total_score"]

    def get_coverage(self):
        """
        number of picked sticks / number of theo sticks
        """
        return self.camel_score["coverage"]

    def get_coverage_template(self):
        """
        number of picked sticks / number of sticks in template
        """
        return self.camel_score["coverage_template"]        
        
    def get_ppm_rms(self):
        return self.camel_score["ppm_rms"]
    
    def get_weighted_sumppm(self):
        return self.camel_score["weighted_sumppm"]       
        
    def get_logmip0intensity(self):
        return self.camel_score["logmip0intensity"]
    
    def get_rt(self):
        """
        this rt corresponds to the spectrum from which the expmips
        were picked and the camel_score calculated.
        returns retention time
        """
        return self.camel_score["rt"]
        
    def set_ria(self, ria):
        self.ria = ria

    def get_ria(self):
        return self.ria	

