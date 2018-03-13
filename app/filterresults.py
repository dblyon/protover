from __future__ import print_function
import os, shutil, cameltoolkit, mfbt, itertools #, multiprocessing, peptide_results, py2r, re, 
import pandas as pd
import numpy as np
from scipy import stats


class Filter(object):
    
    def __init__(self, ancres_dict, outputdir):#!!!
        """        
        outputdir: absolute path to directory
        """
        # pass
        self.set_ancres_dict(ancres_dict)
        self.set_outputdir(outputdir)
        self.rplot_dict = {}

#old Start
    def _get_slopeinterceptmeansd(self, x_vals, y_vals): #alldatapointswithinsd(x_vals, y_vals):
        """
        List, List -> Tuple(Float, Float, Float, Float)
        """
        x_vals = np.array(x_vals)
        y_vals = np.array(y_vals)
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_vals, y_vals)
        y_mean = np.mean(y_vals)
        y_sd = np.std(y_vals, axis=0)
        return slope, intercept, y_mean, y_sd

    def get_aaseqsbydflist(self, df_list, an):
        """
        List(DataFrames), AN -> List(AAseqs)
        produce list of all AAseqs corresponding to given AccessionNumber
        """
        aaseq_list = []
        for df in df_list:
            [aaseq_list.append(aaseq) for aaseq in list(df[df["ANs"] == an]["AAseq"].drop_duplicates())]
        return list(set(aaseq_list))

    def filter_repriapeptide(self, df_list, aaseq_list, threshold):
        """
        List(DataFrame), List(AAseqs) -> Tuple(List-good, List-bad)
        produce Tuple with two Lists (good/pass and bad/no pass) of AAseqs,
        for every AAseq check if RIA_vals of technical replicates are within +/-SD
        """
        aaseq_good = []
        aaseq_bad  = []
        for aaseq in aaseq_list:
            x_vals, y_vals = self.get_tpriasofdflist(df_list, aaseq, withoutTP0=True)
            #print(x_vals, y_vals)
            if self._datapointswithinsd(x_vals, y_vals, threshold):
                aaseq_good.append(aaseq)
            else:
                aaseq_bad.append(aaseq)
        return(aaseq_good, aaseq_bad)  
        
    def filter_coveragetheo(self, df_list, aaseq_list, threshold):
        """
        List(DataFrame), List(String) -> Tuple(List-good, List-bad)
        produce Tuple with two Lists (good/pass and bad/no pass) of AAseqs,
        check if number of sticks deviates +/-10% between replicates
        """
        #aaseq_good = []
        #aaseq_bad  = []
        dict_aaseq2coverage = {}
        df = df_list[0]
        for aaseq in aaseq_list:
            dict_aaseq2coverage[aaseq] = {} # key=AAseq, val={ key=TP, val=[coverage] }
            for df in df_list:
                x15N_colnameslist = self.get_x15Ncolnamesforaaseq(df, aaseq)
                len_theomips = float(len(x15N_colnameslist))
                df_pep = df[df["AAseq"] == aaseq]
                tp_list = list(df_pep["TimePoint"])
                for tp in tp_list:
                    x15N_boollist = []
                    for mz_x15N in x15N_colnameslist:
                        x15N_boollist.append(list(df_pep[df_pep["TimePoint"] == tp][mz_x15N].isnull()))
                    len_expmips = x15N_boollist.count(False)
                    coverage_theomips = len_expmips/len_theomips
                    if not dict_aaseq2coverage[aaseq].has_key(tp):
                        dict_aaseq2coverage[aaseq][tp] = [coverage_theomips]
                    else:
                        dict_aaseq2coverage[aaseq][tp].append(coverage_theomips)
        #print(dict_aaseq2coverage)
        return self._eval_coveragetheodict(dict_aaseq2coverage, threshold)

    def _eval_coveragetheodict(self, d, threshold):
        """
        Dict -> Tuple(List-good, List-bad)
        produce Tuple with two Lists (good/pass and bad/no pass) of AAseqs,
        check if number of sticks deviates +/-10% between replicates.
        key=AAseq, val={ key=TP, val=[coverage] }
        """
        aaseq_good = []
        aaseq_bad  = []
        for aaseq in d.keys():
            deviations_list = []
            for tp in d[aaseq].keys():
                for pair in itertools.combinations(d[aaseq][tp], 2):
                    a,b = pair
                    deviations_list.append(abs(a - b))
            if deviations_list == []:
                deviations_mean = np.NAN
            else:
                deviations_mean = np.mean(deviations_list)
            if np.isnan(deviations_mean):
                aaseq_good.append(aaseq)
            elif deviations_mean < threshold:
                aaseq_good.append(aaseq)
            else:
                aaseq_bad.append(aaseq)
        return(aaseq_good, aaseq_bad)

    def ancres_replicatereproducibility(self, thresh_f1, thresh_f2):
        """
        Dictionary, Directory, Integer, Float -> None
        thresh_f1 = reproducibility of RIA for peptide level (tech./biol. replicates), threshold in StandardDeviations
        thresh_f2 = reproducibility of coverage_theomips for peptide level (tech./biol. replicates), threshold in %
        compareall = switch/flag. if False compare technical replicates (see experiment file) if True compare all replicates.
        write "good" and "bad" results tables and copy analogous 2Dstacks plots to folders
        """
        ancres_list = []
        ancres_dict = self.get_ancres_dict()
        for experiment in ancres_dict: # C_A, C_B, ...
            for replicate in ancres_dict[experiment]: # 1, 2
                ancres_list.append(self.get_ancresfnmostfilteredbyexprep(experiment, replicate))
        #   #indentation was here previously --> Technical replicates
        #indentation here now --> Biological and technical replicates
        df_list = self.ancreslist2dflist(ancres_list) # Each anchor results file made into DataFrame
        aaseq_list = self.get_aaseqlistnoduplicates(df_list, equivalent=False) # Check that all DataFrames "AAseqs" are equivalent and get aaseq_list
        aaseq_good, aaseq_bad = self.filter_repriapeptide(df_list, aaseq_list, threshold=thresh_f1) # FILTER #1: technical reproducibility of RIA
        print("Filter: reproducibility of RIA of replicates")
        print("all, pos, neg", len(aaseq_list), len(aaseq_good), len(aaseq_bad))
        aaseq_good, aaseq_bad_f2 = self.filter_coveragetheo(df_list, aaseq_good, threshold=thresh_f2) # FILTER #2: technical reproducibility of number of sticks
        aaseq_bad += aaseq_bad_f2
        print("\n", "Filter: reproducibility of coverage_theomips of replicates")
        print("all, pos, neg", len(aaseq_list), len(aaseq_good), len(aaseq_bad))
        for experiment in ancres_dict: # C_A, C_B, ... # new due to novel indentation
            ancres_list = []
            for replicate in ancres_dict[experiment]: # 1, 2
                ancres_list.append(self.get_ancresfnmostfilteredbyexprep(experiment, replicate))
            df_list = self.ancreslist2dflist(ancres_list) # Each anchor results file made into DataFrame
            for index, replicate in enumerate(ancres_dict[experiment]):
                df_aaseq_good = self.get_dfbyaaseqlist(df_list[index], aaseq_good)
                #df_aaseq_bad = self.get_dfbyaaseqlist(df_list[index], aaseq_bad)
                fn = ancres_list[index]
                fn = fn.replace(".txt", "_RepRep.txt")
                fn_basename = os.path.basename(fn)
                fn_basename = "temp/" + fn_basename
                self.add_ancresfn2dict(experiment, replicate, fn_basename)
                self.write_df2file(df_aaseq_good, fn)

    def ancres_filterrepriaprotein(self, threshold=2):
        """
        Dictionary, Directory, Integer -> None
        for each AccessionNumber get TPvsRIA of all peptides of only one replicates/measurement/DataFrame
        only calculate if at least 3 peptides per protein
        #comparereplicates: for each AccessionNumber get TPvsRIA of all peptides of all tech. replicates (=Experiments in exp.-file)
        #compareall: for each AccessionNumber get TPvsRIA of all peptides of all biol./tech. replicates
        #if both of the above False: for each AccessionNumber get TPvsRIA of all peptides of only one replicates/measurement/DataFrame
        """
        print("Filter 4: Reproducibility of TPvsRIA for all Peptides of Protein")
        ancres_dict = self.get_ancres_dict()
        for experiment in ancres_dict:
            for replicate in ancres_dict[experiment]:
                ancres_list = []
                ancres_list.append(self.get_ancresfnmostfilteredbyexprep(experiment, replicate)) 
                df_list = self.ancreslist2dflist(ancres_list) 
                df = df_list[0]
                aaseq_good, aaseq_bad = self._get_aaseqgoodbadfromfilterrepriaprotein(df, threshold)
                print(experiment, replicate)
                print("aaseq_good aaseq_bad: ", len(aaseq_good), len(aaseq_bad))
                df_aaseq_good = self.get_dfbyaaseqlist(df, aaseq_good)
                fn = ancres_list[0]
                fn = fn.replace(".txt", "_Filter4.txt")
                fn_basename = os.path.basename(fn)
                fn_basename = "temp/" + fn_basename
                self.add_ancresfn2dict(experiment, replicate, fn_basename)
                self.write_df2file(df_aaseq_good, fn)

    def _get_aaseqgoodbadfromfilterrepriaprotein(self, df, threshold):
        """
        """
        aaseq_good = []
        aaseq_bad = []
        ans_list = self.get_anslistfromdf(df)
        for an in ans_list:
            aaseq_good_temp, aaseq_bad_temp = self._filter_repriaprotein([df], an, threshold) # aaseq_good is empty list if no peps within SD
            [aaseq_good.append(aaseq) for aaseq in aaseq_good_temp]
            [aaseq_bad.append(aaseq) for aaseq in aaseq_bad_temp]
        return(aaseq_good, aaseq_bad)

    def _filter_repriaprotein(self, df_list, an, threshold):
        """
        List(DataFrame), ANs -> Tuple(List-good, List-bad)
        produce Tuple with two Lists (good/pass and bad/no pass) of AAseqs (belonging to given AccessionNumber),
        for every AAseq check if RIA_vals are within +/-SD*threshold, if multiple DataFrames
        given calculate for all DFs
        only calculate if at least 3 peptides per protein
        """
        dict_aaseq2tpria = {} # key = aaseq, val = { key = TP/RIA, val = list} # if multiple DFs append TP/RIA vals to list
        #pair_list_pos = []
        #pair_list_neg = []
        aaseq_list = self.get_aaseqsbydflist(df_list, an)
        # if len(aaseq_list) < 3: # only calculate if at least 3 peptides per protein #!!!
            # return((aaseq_list, []))
        for aaseq in aaseq_list:
            for df in df_list:
                x_vals_temp, y_vals_temp = self.get_tpriasofdflist(df_list, aaseq, withoutTP0=True)
                if not dict_aaseq2tpria.has_key(aaseq):
                    dict_aaseq2tpria[aaseq] = {}
                    dict_aaseq2tpria[aaseq]["tp"] = x_vals_temp
                    dict_aaseq2tpria[aaseq]["ria"] = y_vals_temp
                else:
                    [dict_aaseq2tpria[aaseq]["tp"].append(x) for x in x_vals_temp]
                    [dict_aaseq2tpria[aaseq]["ria"].append(y) for y in y_vals_temp]
        len_combis_reverse = range(1, len(aaseq_list)+1)[::-1]
        for combi in len_combis_reverse: # start with longest combi, then longest-1, ...
            aaseq_pair_list_samelength = []
            for pair in itertools.combinations(aaseq_list, combi): # pair = tuple of AAseqs
                pair = list(pair)
                aaseq_pair_list_samelength.append(pair) # list contains all possible combinations of same length
            aaseq_pair_good = self._find_longestaaseqcombiwithlowestsd(aaseq_pair_list_samelength, dict_aaseq2tpria, threshold)
            if aaseq_pair_good: # no more tests needed, but compare results with same length of pair
                aaseq_complement = [val for val in aaseq_list if not val in aaseq_pair_good]
                return(aaseq_pair_good, aaseq_complement)
        return([], aaseq_list)       

    def _find_longestaaseqcombiwithlowestsd(self, aaseq_pair_list_samelength, dict_aaseq2tpria, threshold):
        """
        List(List(AAseqs)), Dict(key=AAseq, val={key = TP/RIA, val = list}
        start with longest possible aaseq list, 
        compare results -> if True keep result don't test rest
        , if False, test with one AAseq less and compare results, if more than 1 True take the one with lower SD
        if True return list of AAseqs, if False return empty list
        """
        pair_positive = []
        for pair in aaseq_pair_list_samelength:
            x_vals, y_vals = self._get_xyvalsfrom_dict_aaseq2tpria(dict_aaseq2tpria, pair)
            bool, sd = self._datapointswithinsd(x_vals, y_vals, threshold, return_sd=True)
            if bool:
                pair_positive.append([pair, sd])
        if not pair_positive:
            return [] 
        return min(pair_positive, key=lambda x: x[1])[0]

    def _datapointswithinsd(self, x_vals, y_vals, threshold, return_sd=False):
        """
        List, List -> Boolean
        produce True if all y_values are within +/- SD of fit function else False.    
        """
        slope, intercept, y_mean, y_sd = self._get_slopeinterceptmeansd(x_vals, y_vals) #!!! mean and SD should be for each TP individually
        vals = zip(x_vals, y_vals)
        for xy_pair in vals:
            if self._is_datapointwithinsd(xy_pair, y_sd, slope, intercept, threshold):
                continue
            else:
                if return_sd:
                    return(False, 0.0)
                else:
                    return False
        if return_sd:
            return(True, y_sd)
        else:
            return True 

    def _is_datapointwithinsd(self, xy_pair, sd, slope, intercept, threshold):
        """
        Tuple(Float, Float), Float, Float, Float -> True or False
        produce True if y_value between lower and upper boundary else False
        boundaries are +/- SD of fit function y_value at given x_value position
        """
        x_val, y_val = xy_pair
        y_valfit = slope*x_val + intercept
        lower_bound = y_valfit - sd*threshold
        upper_bound = y_valfit + sd*threshold
        if y_val > lower_bound and y_val < upper_bound:
            return True
        else:
            return False 

#old Stop
    def set_rplottemptxt(self, experiment, replicate, fn):
        """
        rplottemp_dict: key = Experiment, val = {key = Replicate, val = R_plot_batch_fn}
        """
        if not self.rplot_dict.has_key(experiment):
            self.rplot_dict[experiment] = {}
            self.rplot_dict[experiment][replicate] = fn
        else:
            self.rplot_dict[experiment][replicate] = fn
                
    def get_rplottemptxt(self, experiment, replicate):
        return self.rplot_dict[experiment][replicate]
        
    def set_ancres_dict(self, ancres_dict):
        """
        ancres_dict: key = Experiment, val = {key = Replicate, val = Ancres_fn}
        ancres_dict: key = Experiment, val = {key = Replicate, val = {key = 0/1/2 (=nofilter/filter12/filter1234), val = Ancres_fn/Ancres_fn/...}}
        """
        self.ancres_dict = ancres_dict
        
    def get_ancres_dict(self):
        return self.ancres_dict

    def get_ancresfnbyexprepfil(self, experiment, replicate, filt):
        ancresfn_base = self.get_ancres_dict()[experiment][replicate][filt]
        outputdir = self.get_outputdir()
        ancresfn = ("%s"+"%s"+"%s"+"%s"+"%s"+"%s") % (outputdir, experiment, "_", replicate, "/", ancresfn_base)
        return ancresfn
    
    def get_ancresfnorigbyexprep(self, experiment, replicate):
        key = sorted(self.get_ancres_dict()[experiment][replicate].keys())[0]
        return self.get_ancresfnbyexprepfil(experiment, replicate, key)

    def get_ancresfnmostfilteredbyexprep(self, experiment, replicate):
        key = sorted(self.get_ancres_dict()[experiment][replicate].keys())[-1]
        return self.get_ancresfnbyexprepfil(experiment, replicate, key)
        
    def add_ancresfn2dict(self, experiment, replicate, ancresfn):
        key = sorted(self.get_ancres_dict()[experiment][replicate].keys())[-1] + 1
        self.get_ancres_dict()[experiment][replicate][key] = ancresfn
        
    def set_outputdir(self, outputdir):
        self.outputdir = outputdir

    def get_outputdir(self):
        return self.outputdir
               
    def ancreslist2dflist(self, ancres_list):
        """
        List(FileName(s)) -> List(DataFrame(s))
        produce a list with pandas.dataframe(s)
        """
        df_list = []
        for ancres in ancres_list:
            df_list.append(pd.read_csv(ancres, sep="\t"))
        return df_list

    def get_aaseqlistnoduplicates(self, df_list, equivalent=True):
        """
        List(DataFrame(s)) -> List or None
        produce list of AminoAcidSequences (non duplicate entries of column "AAseq") 
        if the AAseq_list of all DataFrames is equivalent else None
        """
        aaseq_metalist = []
        for df in df_list:
            aaseq_list = sorted(list(df["AAseq"].drop_duplicates()))
            aaseq_metalist.append(aaseq_list)
        if equivalent: # all lists have to be equivalent
            for index, aaseq_list in enumerate(aaseq_metalist):
                if (index+1) == len(aaseq_metalist):
                    break
                if aaseq_metalist[index] == aaseq_metalist[index+1]:
                    continue
                else:
                    print("Filterresults.get_aaseqlistnoduplicates: AAseq list of ancres not equal.")
                    return None
        else: # take intersection of all lists
            return list(set.intersection(*map(set, aaseq_metalist)))
        return aaseq_list
            
    def get_tpria(self, df, aaseq, withoutTP0=False):
        """
        DataFrame, String -> List(List, List)
        produce list with TimePoint_vals_list and RIA_vals_list for given DataFrame and AAseq
        """
        df_pep = df[df["AAseq"] == aaseq]
        x_vals = []
        y_vals = []
        if withoutTP0:
            [x_vals.append(float(x)) for x in list(df_pep["TimePoint"])[1:]]
            [y_vals.append(y) for y in list(df_pep["RIA"])[1:]]
        else:
            [x_vals.append(float(x)) for x in list(df_pep["TimePoint"])]
            [y_vals.append(y) for y in list(df_pep["RIA"])]
        return[x_vals, y_vals]

    def get_tpriasofdflist(self, df_list, aaseq, withoutTP0):
        """
        List(DataFrame(s)), String -> List(List, List)
        produce get_tpria for all DataFrames
        """
        x_vals_all = []
        y_vals_all = []
        for df in df_list:
            x_vals, y_vals = self.get_tpria(df, aaseq, withoutTP0)
            [x_vals_all.append(x) for x in x_vals]
            [y_vals_all.append(y) for y in y_vals]
        return[x_vals_all, y_vals_all]

    def get_x15Ncolnamesforaaseq(self, df, aaseq):
        """
        DataFrame, String -> List
        produce list of x15N column names for given AAseq and DF
        """
        df_pep = df[df["AAseq"] == aaseq]
        numN = int(df_pep["#N"].iloc[0])
        col_list = []
        numN_list = range(0, numN)
        for ele in numN_list:
            colname = "mz_" + str(ele) + "x15N"
            col_list.append(colname)
        return col_list
        
    def get_aaseqsofanbydflist(self, df_list, an):
        """
        List(DataFrames), AN -> List(AAseqs)
        produce list of all AAseqs corresponding to given AccessionNumber
        """
        aaseq_list = []
        for df in df_list:
            [aaseq_list.append(aaseq) for aaseq in list(df[df["ANs"] == an]["AAseq"].drop_duplicates())]
        return list(set(aaseq_list))
        
    def _get_xyvalsfrom_dict_aaseq2tpria(self, dict, aaseq_list):
        """
        # key = aaseq, val = { key = TP/RIA, val = list} # if multiple DFs append TP/RIA vals to list
        """
        x_vals = []
        y_vals = []
        for aaseq in aaseq_list:
            [x_vals.append(x) for x in dict[aaseq]["tp"]]
            [y_vals.append(y) for y in dict[aaseq]["ria"]]           
        return(x_vals, y_vals)
        
    def copy_files2folders(self, filename_list, folder_2copyinto, directory_parent):
        """
        List(FileNames), FolderName, FolderName -> None
        copy all files of filename_list from directory_parent to folder_2copyinto    
        """
        cwd_orig = os.getcwd()
        os.chdir(directory_parent)
        if not os.path.exists(folder_2copyinto):
            os.makedirs(folder_2copyinto)
        for filename in filename_list:
            try:
                shutil.copy(filename, folder_2copyinto)
            except IOError:
                pass
        os.chdir(cwd_orig)
        
    def get_filenamefrompartialfilename(self, filenamepartial, filename_list):
        """
        FileName, List(FileNames) -> List(FileName(s))
        produce List of FileNames for all files containing filenamepartial.
        """
        list_ofcurrentfileswithfilename = []
        for fileo in filename_list:
            fileo_list = fileo.split("_")
            for ele in fileo_list:
                if filenamepartial == ele:
                    list_ofcurrentfileswithfilename.append(fileo)
        return list_ofcurrentfileswithfilename

    def exchangepartialwithrealfilenames(self, filenamepartial_list, directory):
        """
        List(FileNames), Directory -> List(FileNames)
        produce List of FileNames of given directory if filenamepartial contained in filename.
        """
        cwd_orig = os.getcwd()
        os.chdir(directory)
        filename_list_orig = os.listdir(directory)
        filename_list_new = []
        for filenamepartial in filenamepartial_list:
            filename_new = self.get_filenamefrompartialfilename(filenamepartial, filename_list_orig)
            if len(filename_new) == 1:
                filename_list_new.append(filename_new[0])
                continue
            else:
                [filename_list_new.append(x) for x in filename_new]
        os.chdir(cwd_orig)
        return filename_list_new
        
    def get_dfbyaaseqlist(self, df, aaseq_list):
        """
        DataFrame, List(AAseq) -> DataFrame
        produce DataFrame of all rows containing AAseq of aaseq_list in column "AAseq"
        """
        df_list = []
        if not aaseq_list:
            return pd.DataFrame(columns=df.columns)
        for aaseq in aaseq_list:
            df_list.append(df[df["AAseq"] == aaseq])
        return pd.concat(df_list, ignore_index=True)

    def _sort_x15Ncolnames_natkeys(self, df):
        """
        """
        mz_x15N_list = []
        int_x15N_list = []
        colnames = list(df.columns)
        for ele in colnames:
            if "mz_" in ele:
                mz_x15N_list.append(ele)
            if "int_" in ele:
                int_x15N_list.append(ele)
        mz_x15N_list = sorted(mz_x15N_list, key=cameltoolkit.natural_keys)
        int_x15N_list = sorted(int_x15N_list, key=cameltoolkit.natural_keys)
        return(mz_x15N_list, int_x15N_list)
        
    def write_df2file(self, df, filename, sep="\t", decimalpoint=".", index=False, header=True):
        df["Rt_spec"] = df["Rt_spec"].apply(lambda x: round(x, 3))
        df = df.sort_values(["AAseq", "TimePoint"], ascending=[True, True])
        col_list_sorted = ["AAseq", "charge", "Rt_input", "ANs", "TimePoint", "Filename", "Rt_spec", "RIA", "#N", "sum14N", "sum15N"]        
        mz_x15N_list, int_x15N_list = self._sort_x15Ncolnames_natkeys(df)
        col_list_sorted += mz_x15N_list + int_x15N_list
        df = df.reindex(col_list_sorted, axis=1)
        if decimalpoint == ".":
            df.to_csv(filename, sep="\t", header=header, index=index)
        else:
            df = df.applymap(lambda x: str(x).replace(".", ","))
            df = df.applymap(lambda x: str(x).replace(",mzML", ".mzML"))
            df.to_csv(filename, sep=sep, header=header, index=index)

    def get_anslistfromdf(self, df):
        """
        DataFrame -> List(ANs)
        produce list of AccessionNumbers of given test_data frame, without duplicates.
        """
        ansstrings_list = list(df["ANs"].drop_duplicates())
        an_list = []
        for ele in ansstrings_list:
            if type(ele) == str:
                [an_list.append(x.strip()) for x in ele.split(",")]
        return list(set(an_list))
        
    def get_anscolentriesofan(self, df, an):
        """
        DataFrame, String -> List(String(s)/List)
        produce List of row entries if given an within row entry (looking at column "ANs")
        e.g. given an = I3TAJ6, col_entry = 'I3TAJ6, TC174430, contig_57177_1.1'
        return ['I3TAJ6', 'TC174430', 'contig_57177_1.1']
        """
        anscolentries_list = []
        ans_df = df["ANs"]
        for row in ans_df:
            temp_ans_list = []
            if type(row) == str:
                [temp_ans_list.append(x.strip()) for x in row.split(",")]
            if an in temp_ans_list:
                anscolentries_list.append(row)
        return list(set(anscolentries_list))

    def get_dfbyan(self, df, an):
        """
        DataFrame, String -> DataFrame
        produce reindexed dataframe of all rows if given an within "ANs" column
        """
        anscolentries_list = self.get_anscolentriesofan(df, an)
        df_list = []
        for ancolentry in anscolentries_list:
            df_list.append(df[df["ANs"] == ancolentry])
        return pd.concat(df_list, ignore_index=True)
             
    def plot_txt_tpvsria(self, df, py2r_plot_instance, plottitle):
        """
        DataFrame -> String
        produce TimePoint vs. RIA plot text for R    
        """
        aaseq_list = sorted(list(df["AAseq"].drop_duplicates()))
        tp_ria_list = []
        for aaseq in aaseq_list:
            tp, ria = self.get_tpria(df, aaseq)
            temp = zip(tp, ria)
            temp2 = []
            [temp2.append(list(ele)) for ele in temp]
            tp_ria_list.append(temp2)  
        py2r_plot_instance.plot_prot2rias(aaseq_list, tp_ria_list, plottitle)
        
    def make_directories(self, dir_list):
        """
        List -> None
        create directories in list if not already existing
        """
        for directory in dir_list:
            if not os.path.exists(directory):
                os.makedirs(directory)
                       
    def copy_goodbad2dstacks2folders(self, aaseq_list, directory_parent, subdir):
        """
        List, Directory -> None
        compare directory_parent to aaseq_list and derive true filenames 
        by splitting at "_" and thereby create list of files,
        copy list of files to subdirectory of directory_parent
        """    
        filename_list = aaseq_list
        directory_2copyinto = directory_parent+"/"+subdir   
        filename_list_correct = self.exchangepartialwithrealfilenames(filename_list, directory_parent) #!!!
        self.copy_files2folders(filename_list_correct, directory_2copyinto, directory_parent)              
        
    def move_2dstacks(self):
        """
        None -> None
        move aaseq_good to "positive" subdirectory of "2Dstacks" and aaseq_bad to "negative" subdirectory
        """
        ancres_dict = self.get_ancres_dict()
        for experiment in ancres_dict:
            for replicate in ancres_dict[experiment]:
                outputdir_2dstacks = self.get_outputdir() + experiment + "_" + replicate + "/" + "2Dstacks"
                aaseq_good, aaseq_bad = self.get_goodbadaaseqs_compdfoldnew(experiment, replicate)
                self._move_goodbad2dstacks2folders(aaseq_good, outputdir_2dstacks, "positive")
                self._move_goodbad2dstacks2folders(aaseq_bad, outputdir_2dstacks, "negative")        

    def _move_goodbad2dstacks2folders(self, aaseq_list, directory_parent, subdir):
        """
        List, Directory -> None
        compare directory_parent to aaseq_list and derive true filenames 
        by splitting at "_" and thereby create list of files,
        copy list of files to subdirectory of directory_parent
        """    
        filename_list = aaseq_list
        directory_2copyinto = directory_parent+"/"+subdir   
        filename_list_correct = self.exchangepartialwithrealfilenames(filename_list, directory_parent)
        self._move_files2folders(filename_list_correct, directory_2copyinto, directory_parent)                 

    def _move_files2folders(self, filename_list, folder_2copyinto, directory_parent):
        """
        List(FileNames), FolderName, FolderName -> None
        copy all files of filename_list from directory_parent to folder_2copyinto    
        """
        cwd_orig = os.getcwd()
        os.chdir(directory_parent)
        if not os.path.exists(folder_2copyinto):
            os.makedirs(folder_2copyinto)
        for filename in filename_list:
            fn2rm = folder_2copyinto + "/" + os.path.basename(filename)
            if os.path.isfile(fn2rm):
                os.remove(fn2rm)
            try:
                shutil.move(filename, folder_2copyinto)
            except:
                pass        
        os.chdir(cwd_orig)
        
    def ancres_filterriaincreasing(self, threshold=None):
        """
        Dictionary, Directory -> None
        check TPvsRIA increasing per peptide per replicate separately,
        """
        ancres_dict = self.get_ancres_dict()
        for experiment in ancres_dict:
            for replicate in ancres_dict[experiment]:
                print(experiment, replicate)
                ancres_list = []       
                aaseq_good = []
                aaseq_bad = []
                ancres_list.append(self.get_ancresfnmostfilteredbyexprep(experiment, replicate))
                df_list = self.ancreslist2dflist(ancres_list)
                df = df_list[0]
                aaseq_list = self.get_aaseqlistnoduplicates([df])
                for aaseq in aaseq_list:
                    tp_list, ria_list = self.get_tpriasofdflist([df], aaseq, withoutTP0=False)
                    if self._is_riaincreasing(tp_list, ria_list, threshold=threshold): 
                        aaseq_good.append(aaseq)
                    else:
                        aaseq_bad.append(aaseq)
                print("Filter RIA increasing: Pos and neg ", len(aaseq_good), len(aaseq_bad))
                df_aaseq_good = self.get_dfbyaaseqlist(df, aaseq_good)              
                fn = ancres_list[0]
                fn = fn.replace(".txt", "_RIAinc.txt")
                fn_basename = os.path.basename(fn)
                fn_basename = "temp/" + fn_basename
                self.add_ancresfn2dict(experiment, replicate, fn_basename)
                self.write_df2file(df_aaseq_good, fn)

    def plot_2Dstacks(self, plot):
        """
        """
        ancres_dict = self.get_ancres_dict()
        for experiment in ancres_dict:
            for replicate in ancres_dict[experiment]:
                df = self.ancreslist2dflist([self.get_ancresfnmostfilteredbyexprep(experiment, replicate)])[0]
                aaseq_list = self.get_aaseqlistnoduplicates([df])
                for aaseq in aaseq_list:
                    tprtmzint_list_norm = self._get_mzintpairs_alltps_norm2max(aaseq, df)
                    (xvals_list, yvals_list, yvals_orig_list, yaxis_tps_h, yaxis_positions, rt_list) = self._2Dstacks_preparedataforplot(tprtmzint_list_norm, aaseq)
                    if len(xvals_list) == 0: #!!! should I check for -1 vals since they don't appear in plots and thus no info, -->  enter in log file.
                        continue
                    plottitle = aaseq + "_" + experiment + "_" + replicate
                    plotfilename = plottitle
                    plotfilelocation = self.get_outputdir() + experiment + "_" + replicate + "/" + "2Dstacks"
                    self.make_directories([plotfilelocation])
                    plotfilelocation = plotfilelocation.replace("\\", "/") + "/"
                    plot.set_r_plottemptxt(self.get_rplottemptxt(experiment, replicate))
                    chargestate = df[df["AAseq"] == aaseq]["charge"].iloc[0]
                    theomips = mfbt.Peptide(aaseq, chargestate=chargestate).gettheomips()
                    xlimbeg = theomips[0][0] - 0.7
                    xlimend = theomips[-1][0] + 0.7
                    plot.plotcamelstacks(xvals_list, yvals_list, yvals_orig_list, plottitle, plotfilename, plotfilelocation, xlimbeg, xlimend, yaxis_tps_h, yaxis_positions, rt_list)
                    
    def _2Dstacks_preparedataforplot(self, tprtmzint_list_norm, aaseq):
        """
        """
        xvals_list = []
        yvals_list = []
        yvals_orig_list = []
        yaxis_tps_h = []
        yaxis_positions = []
        rt_list = []
        for index, tprtmzint in enumerate(tprtmzint_list_norm):
            tp, rt, mzint = tprtmzint
            yaxis_positions.append(tp)
            yaxis_tps_h.append(index*120+50)
            rt_list.append(rt)
            for ele in mzint:
                expmz = ele[0]
                expint = ele[1] + (index*120)
                yval_orig = (index*120)
                xvals_list.append(expmz)
                yvals_list.append(expint)
                yvals_orig_list.append(yval_orig)
        return(xvals_list, yvals_list, yvals_orig_list, yaxis_tps_h, yaxis_positions, rt_list)

    def _is_riamissingtoomanyvalues(self, tp_list, ria_list, threshold=None):
        """
        List, List [, Float] -> Boolean
        return False if more than threshold*possible vals missing
        (default is 50%)
        else True
        """
        if not threshold:
            threshold = 0.5
        missingvalscounter = 0
        for ria in ria_list:
            if ria < 0.0:
                missingvalscounter += 1
        if missingvalscounter >= int(len(tp_list)*threshold):
            return True
        else:
            return False
        
    def _is_riaincreasing(self, tp_list, ria_list, threshold):
        """
        List, List -> Boolean
        sort by increasing tp vals,
        return True if RIA same or higher over time, else False.
        """
        tp_ria = zip(tp_list, ria_list)
        if not tp_ria: #!!!
            return False
        if self._is_riamissingtoomanyvalues(tp_list, ria_list, threshold=0.5) == True:
            return False
        tp_ria_sorted = sorted(tp_ria, key=lambda ele: ele[0])
        for index, ria in enumerate(tp_ria_sorted):
            if index+1 == len(tp_ria_sorted): # stop when reaching the end of the iterable
                break
            if not threshold:
                if abs(ria[1] - tp_ria_sorted[index+1][1]) < 0.0000001 or ria[1] < tp_ria_sorted[index+1][1]:
                    continue
                else:
                    return False
            else:
                if abs(ria[1] - tp_ria_sorted[index+1][1]) < (threshold*tp_ria_sorted[index+1][1]) or ria[1] < tp_ria_sorted[index+1][1]:
                    continue
                else:
                    return False
        return True
        
    def writenewfile(self, fn):
        fh = open(fn, "w")
        fh.close()

    def get_goodbadaaseqs_compdfoldnew(self, experiment, replicate):
        """
        None -> List, List
        produce AAseq_good (from most filtered DataFrame) and 
        AAseq_bad list (complement of DataFrame_orig - AAseq_good).
        """
        ancresfn_orig = self.get_ancresfnorigbyexprep(experiment, replicate)
        df_orig = self.ancreslist2dflist([ancresfn_orig])[0]
        ancresfn_filt = self.get_ancresfnmostfilteredbyexprep(experiment, replicate)
        df_filt = self.ancreslist2dflist([ancresfn_filt])[0]
        aaseq_orig = self.get_aaseqlistnoduplicates([df_orig])
        aaseq_filt = self.get_aaseqlistnoduplicates([df_filt])
        aaseq_bad = [x for x in aaseq_orig if x not in aaseq_filt]
        return(aaseq_filt, aaseq_bad)
        
    def replot_tpvsria(self, plot, label="positive", dbl=False):
        """
        for every replicate of every experiment plot TP vs RIA
        #ancres_dict: key = Experiment, val = {key = Replicate, val = Ancres_fn}
        replot TPvsRIA only with good AAseqs to subdirectory of given folder
        """
        ancres_dict = self.get_ancres_dict()
        for experiment in ancres_dict:
            for replicate in ancres_dict[experiment]:
                ancres_fn = self.get_ancresfnmostfilteredbyexprep(experiment, replicate)
                df = pd.read_csv(ancres_fn, sep="\t")
                an_list = self.get_anslistfromdf(df)
                plot.set_r_plottemptxt(self.get_rplottemptxt(experiment, replicate))
                if dbl:
                    plotfilelocation = self.get_outputdir() + experiment + "_" + replicate+"/" + "RIA" + "/" + label
                else:
                    plotfilelocation = self.get_outputdir() + experiment + "_" + replicate + "/" + "RIA"
                self.make_directories([plotfilelocation])
                plotfilelocation = plotfilelocation.replace("\\", "/")
                plotfilelocation += "/"
                plot.setplotfilelocation(plotfilelocation)
                for an in an_list:
                    df_an = self.get_dfbyan(df, an)
                    plottitle = an + "_" + experiment + replicate
                    self.plot_txt_tpvsria(df_an, plot, plottitle)

    def create_rplotbatchtxt(self, label="positive"):
        """
        """
        ancres_dict = self.get_ancres_dict()
        for experiment in ancres_dict:
            for replicate in ancres_dict[experiment]:
                rplottemptxt_fn = self.get_outputdir() + experiment + "_" + replicate+"/" + "temp/" + "PlotTXTForR"+label+".txt"
                self.writenewfile(rplottemptxt_fn)
                self.set_rplottemptxt(experiment, replicate, rplottemptxt_fn)
    
    def run_rplotbatch(self, plot, winorunix):
        """
        """
        x = 1
        ancres_dict = self.get_ancres_dict()
        for experiment in ancres_dict:
            for replicate in ancres_dict[experiment]:
                fn = self.get_rplottemptxt(experiment, replicate)
                r_plot_batch = plot.run_rbatchtxtfile(fn, winorunix)
                if r_plot_batch == 0:
                    pass
                else:
                    print("An error occurred when attempting to plot using R.")
                    print("Comment out 'fr.cleanup_files()' in 'run_experimentfile.py' to save the R-batch-plot-file (which would otherwise be deleted after execution). Or set 'dbl=True' which will prevent the 'temp' folder to be deleted.")
        return x                        

    def _make_bins(self, bins):
        if not bins:
            bins = range(0, 101, 1)
            bins = [x/100.0 for x in bins]
        else:
            (start, stop, interval) = bins
            start = int(start*1000)
            stop = int(stop*1000)
            interval= int(interval*1000)
            bins = range(start, stop, interval)
            bins = [x/1000.0 for x in bins]
        return bins
        
    def replot_histos(self, plot, bins, label="positive", plotfrequency=None, plotdensity=None, dbl=False):
        """
        """
        ancres_dict = self.get_ancres_dict()
        for experiment in ancres_dict:
            for replicate in ancres_dict[experiment]:
                ancres_fn = self.get_ancresfnmostfilteredbyexprep(experiment, replicate)
                df = pd.read_csv(ancres_fn, sep="\t")                
                plot.set_r_plottemptxt(self.get_rplottemptxt(experiment, replicate))
                if dbl:
                    plotfilelocation = self.get_outputdir()+experiment + "_" + replicate+"/"+"Histos"+"/"+label
                else:
                    plotfilelocation = self.get_outputdir()+experiment + "_" + replicate+"/"+"Histos"+"/"
                self.make_directories([plotfilelocation])
                plotfilelocation = plotfilelocation.replace("\\", "/")
                plotfilelocation += "/"
                plot.setplotfilelocation(plotfilelocation)                                
                x_axis_label = "Coverage [%]" 
                y_axis_label_freq = "Frequency [#] of Peptides"
                y_axis_label_dens = "Density [%] of Peptides"
                tp_list = sorted(list(df["TimePoint"].drop_duplicates()))
                for tp in tp_list:
                    plotfilename = experiment+"_"+replicate
                    plottitle = ("%s%s%s") % ("Coverage of template ", plotfilename, tp)
                    fractionfound_list = self.get_coveragetheomips_allpeps(experiment, replicate, tp)
                    if plotdensity:
                        plotfilename_d = ("%s%s%s%s") % ("HistoDens_", plotfilename, "_", tp)
                        self._plot_hist_coveragetemplate(plot, bins, plotfilename_d, plottitle, x_axis_label, y_axis_label_freq, y_axis_label_dens, fractionfound_list, plotfrequency=False, plotdensity=True)
                    if plotfrequency:
                        plotfilename_f = ("%s%s%s%s") % ("HistoFreq_", plotfilename, "_", tp)
                        self._plot_hist_coveragetemplate(plot, bins, plotfilename_f, plottitle, x_axis_label, y_axis_label_freq, y_axis_label_dens, fractionfound_list, plotfrequency=True, plotdensity=False)
                        
    def _plot_hist_coveragetemplate(self, plot, bins, plotfilename, plottitle, x_axis_label, y_axis_label_freq, y_axis_label_dens, fractionfound_list, plotfrequency, plotdensity):
        plot.plot_hist_coveragetemplate(bins, plotfilename, plottitle, x_axis_label, y_axis_label_freq, y_axis_label_dens, fractionfound_list, plotfrequency, plotdensity)
        
    def cleanup_files(self, decimalpoint):
        """
        delete all DataFrames except for most filtered one,
        rename it to orig and set decimalpoint as given.
        Delete all R-batch-plot-files including Rout
        """
        ancres_dict = self.get_ancres_dict()
        for experiment in ancres_dict:
            for replicate in ancres_dict[experiment]:
                ancres_mostfiltered_fn = self.get_ancresfnmostfilteredbyexprep(experiment, replicate)
                ancres_orig_fn = self.get_ancresfnorigbyexprep(experiment, replicate)
                df = pd.read_csv(ancres_mostfiltered_fn, sep="\t")
                fn = self.get_outputdir() + experiment + "_" + replicate + "/" +  os.path.basename(ancres_orig_fn)
                # df["Rt_spec"] = df["Rt_spec"].apply(lambda x: round(x, 3))
                self.write_df2file(df, fn, decimalpoint=decimalpoint)
                dirname, filename = os.path.split(os.path.abspath(ancres_mostfiltered_fn))
                shutil.rmtree(dirname, ignore_errors=True)

    def get_coveragetheomips_allpeps(self, experiment, replicate, tp):
        """
        produce list of coverage_theomips of given TimePoint of given Experiment and Replicate.
        """
        coverage_theomips_list = []
        ancres_fn = self.get_ancresfnmostfilteredbyexprep(experiment, replicate)
        df = self.ancreslist2dflist([ancres_fn])[0]
        aaseq_list = list(df["AAseq"].drop_duplicates())
        for aaseq in aaseq_list:
            x15N_colnameslist = self.get_x15Ncolnamesforaaseq(df, aaseq)
            len_theomips = float(len(x15N_colnameslist))
            df_pep = df[df["AAseq"] == aaseq]
            x15N_boollist = []
            for mz_x15N in x15N_colnameslist:
                x15N_boollist.append(list(df_pep[df_pep["TimePoint"] == tp][mz_x15N].isnull()))
            len_expmips = x15N_boollist.count(False)
            coverage_theomips = len_expmips/len_theomips
            coverage_theomips_list.append(coverage_theomips)
        return coverage_theomips_list        
    
    def cleanTP0(self): #!!! slow function profile this function and find faster variant
        """
        addendum: 
        add 2 new columns to DataFrame: sum14N and sum15N
        """
        ancres_dict = self.get_ancres_dict()
        for experiment in ancres_dict:
            for replicate in ancres_dict[experiment]:
                ancres_list = []
                ancres_list.append(self.get_ancresfnorigbyexprep(experiment, replicate))
                df_list = self.ancreslist2dflist(ancres_list)
                df = df_list[0]
                df["sum14N"] = -1.0
                df["sum15N"] = -1.0
                aaseq_list = self.get_aaseqlistnoduplicates([df])            
                for aaseq in aaseq_list:
                    df_temp = df[df["AAseq"] == aaseq]
                    tp_min = sorted(list(df_temp["TimePoint"]))[0]
                    df_row_temp = df_temp[df_temp["TimePoint"] == tp_min]
                    df_row_clean = self._cleanexpmipsTP0(df_row_temp)
                    df = df.drop(df_row_temp.index)                    
                    df = pd.concat([df, df_row_clean], ignore_index=True) 
                fn = ancres_list[0]
                fn_cleanTP0 = fn.replace(".txt", "_cleanTP0.txt")
                fn_cleanTP0_basename = "temp/" + os.path.basename(fn_cleanTP0)
                self.add_ancresfn2dict(experiment, replicate, fn_cleanTP0_basename)
                self.write_df2file(df, fn_cleanTP0)                                       

    def cleanTP0_v2(self): #!!! slow function, profile this function and find faster variant
        """
        addendum: 
        add 2 new columns to DataFrame: sum14N and sum15N
        """
        ancres_dict = self.get_ancres_dict()
        for experiment in ancres_dict:
            for replicate in ancres_dict[experiment]:
                ancres_list = []
                ancres_list.append(self.get_ancresfnorigbyexprep(experiment, replicate))
                df_list = self.ancreslist2dflist(ancres_list)
                df = df_list[0]
                df["sum14N"] = -1.0 #!!! helper function        addendum: add 2 new columns to DataFrame: sum14N and sum15N
                df["sum15N"] = -1.0                                
                aaseq_list = self.get_aaseqlistnoduplicates([df])
                tp_min = sorted(list(df["TimePoint"]))[0]                
                for aaseq in aaseq_list:                    
                    dft = df[df["AAseq"] == aaseq]
                    index_aaseq_tpmin = dft[dft['TimePoint'] == tp_min].index.values                    
                    cols_list = self._cleanexpmipsTP0_v2(df, index_aaseq_tpmin)                    
                    for col in cols_list:
                        df.loc[index_aaseq_tpmin, col] = np.NaN
                fn = ancres_list[0]
                fn_cleanTP0 = fn.replace(".txt", "_cleanTP0.txt")
                fn_cleanTP0_basename = "temp/" + os.path.basename(fn_cleanTP0)
                self.add_ancresfn2dict(experiment, replicate, fn_cleanTP0_basename)
                self.write_df2file(df, fn_cleanTP0)       
            
    def _get_mzintpairs_alltps_norm2max(self, aaseq, df):
        """
        AAseq, DataFrame -> nestedList
        [[tp0, [(mz1, int1), (mz2, int2), ...]],
        [tp1, [(mz1, int1), (mz2, int2), ...]], ...]
        """
        tpmzint = []
        df = df[df["AAseq"] == aaseq]
        tp_list_ascending = sorted(list(df["TimePoint"]))
        for tp in tp_list_ascending:
            mz_list  = []
            int_list = []
            df_tp = df[df["TimePoint"] == tp].dropna(axis=1, how="all")
            rt = round(df_tp["Rt_spec"].values[0], 2)
            df_col_list = sorted(df_tp.columns.values, key=cameltoolkit.natural_keys)
            for colname in df_col_list:
                if "mz_" in colname:
                    mz_list.append(df_tp[colname].values[0])
                if "int_" in colname:
                    int_list.append(df_tp[colname].values[0])
            mzint_list_norm = self.norm2max(zip(mz_list, int_list))
            tpmzint.append([tp, rt, mzint_list_norm])
        return tpmzint
    
    def norm2max(self, mzint_list):
        """
        List -> List
        finds maximum intensity, returns copy of spectrum normalized to basepeak
        """
        returnlist = []
        try:
            max_intensity = sorted(mzint_list, key=lambda ele: int(ele[1]), reverse=True)[0][1]    # TODO: remove if next line is equivalent
        except IndexError:
            return mzint_list
        for mzint in mzint_list:
            returnlist.append([mzint[0]] + [(mzint[1]/max_intensity)*100.0])
        return returnlist
    
    def _get_expmips(self, df):
        """
        DataFrame(of one TP = one row) -> List
        expmips = [mz-float, inte-float, x15N-int]
        """
        expmips = []
        df = df.dropna(axis=1, how="all")
        for colname in df.columns.values:
            if "mz_" in colname:
                mz = df[colname].values[0]
                x15N = colname.split("_")[1].split("x")[0]
                int_colname = "int_" + x15N + "x15N"
                intensity = df[int_colname].values[0]
                expmips.append([mz, intensity, int(x15N)])
        return sorted(expmips, key=lambda ele: int(ele[-1]))
                
    def _cleanexpmipsTP0(self, df):
        """
        DataFrame -> DataFrame
        for TP0 (minimum labeling) look for first empty x15N position,
        remove all expmips following that position.
        """        
        expmips = self._get_expmips(df)
        expmips_cleanTP0 = self._findgap(expmips)
        colnames = list(df.columns)
        for ele in colnames: # set all mz and inte values to NaN
            if "mz_" in ele or "int_" in ele:
                df.loc[:, ele] = np.nan
        for mz_int_x15N in expmips_cleanTP0:
            mz, inte, x15N = mz_int_x15N
            mz_x15N  = "mz_"  + str(x15N) + "x15N"
            int_x15N = "int_" + str(x15N) + "x15N"
            df.loc[:, mz_x15N] = mz
            df.loc[:, int_x15N] = inte
        return df

    def _cleanexpmipsTP0_v2(self, df, index_aaseq_tpmin):
        """
        DataFrame -> ListOfColNames
        for TP0 (minimum labeling) look for first empty x15N position,
        return the column_names as List of all expmips following that position,
        which should be deleted.
        """
        cols_list = []
        df_row = df.iloc[index_aaseq_tpmin]                
        expmips = self._get_expmips(df_row)      
        x15N_list = self._findgap_v2(expmips)
        for x15N in x15N_list:
            mz_x15N  = "mz_"  + str(x15N) + "x15N"
            int_x15N = "int_" + str(x15N) + "x15N"
            cols_list.append(mz_x15N)
            cols_list.append(int_x15N)
        return cols_list

    def _findgap_v2(self, expmips):
        '''
        List (of expmips: [mz, int, x15N]) --> List(of x15N positions to delete)
        '''
        x15N_list = []
        for index, mips in enumerate(expmips):
            x15N = mips[-1]
            if x15N != index:                
                x15N_all = [ele[-1] for ele in expmips]
                return sorted(list(set(x15N_all) - set(x15N_list)))
            else:
                x15N_list.append(x15N)
        return []

    def _findgap(self, expmips):
        """
        [mz, x15N]
        look for first empty x15N position. remove all expmips following that position.
        """
        expmips_after = []
        for index, mips in enumerate(expmips):
            x15N = mips[-1]
            if x15N != index:
                break
            else:
                expmips_after.append(mips)
        return expmips_after        

    def calcRIA(self):
        """
        """
        ancres_dict = self.get_ancres_dict()
        for experiment in ancres_dict:
            for replicate in ancres_dict[experiment]:
                ancres_list = []
                ancres_list.append(self.get_ancresfnmostfilteredbyexprep(experiment, replicate))
                df = self.ancreslist2dflist(ancres_list)[0]
                aaseq_list = self.get_aaseqlistnoduplicates([df])
                for aaseq in aaseq_list:
                    df_aaseq = df[df["AAseq"] == aaseq]
                    df_aaseq_ria = self._calc_RIA_forAAseq(df_aaseq)
                    df = df.drop(df_aaseq.index)
                    df = pd.concat([df, df_aaseq_ria], ignore_index=True)
                fn = ancres_list[0]
                fn_RIA = fn.replace(".txt", "_RIA.txt")
                fn_RIA_basename = "temp/" + os.path.basename(fn_RIA)
                self.add_ancresfn2dict(experiment, replicate, fn_RIA_basename)
                self.write_df2file(df, fn_RIA)        
                    
    def _calc_RIA_forAAseq(self, df):
        """
        DataFrame -> DataFrame
        calculates q ratio according to "Gustavsson et al., Proteomics,  2005"
        q = 15N / (14N + 15N)
        define which x15N are 14N by looking at TP min labeled,
        rest are assumed to be 15N.
        """
        df = df.dropna(axis=1, how="all")
        tp_list = sorted(list(df["TimePoint"]))
        tpmin = tp_list[0]
        df_tpmin_row = df[df["TimePoint"] == tpmin]
        expmips_tpmin = self._get_expmips(df_tpmin_row)
        x14N_positions = [x[-1] for x in expmips_tpmin]
        for tp in tp_list:
            df_tp = df[df["TimePoint"] == tp]
            expmips = self._get_expmips(df_tp)	
            sum_x14N = 0.0
            sum_x15N = 0.0
            for mzintppmx15N in expmips:
                if mzintppmx15N[-1] in x14N_positions:
                    sum_x14N+= mzintppmx15N[1]
                else: sum_x15N+= mzintppmx15N[1]
            if sum_x14N == 0.0:
                q = -1.0
            else:
                q = (sum_x15N / (sum_x14N + sum_x15N))
            df.loc[df_tp.index, "RIA"] = q
        return df 
        
    def calcRIA_relint(self):
        """
        taking relative intensities at TP0 
        into account for Light and Heavy calculations
        """
        ancres_dict = self.get_ancres_dict()
        for experiment in ancres_dict:
            for replicate in ancres_dict[experiment]:
                ancres_list = []
                ancres_list.append(self.get_ancresfnmostfilteredbyexprep(experiment, replicate))
                df = self.ancreslist2dflist(ancres_list)[0]
                aaseq_list = self.get_aaseqlistnoduplicates([df])
                for aaseq in aaseq_list:
                    df_aaseq = df[df["AAseq"] == aaseq]
                    df_aaseq_ria = self._calc_RIA_forAAseq_relint(df_aaseq)
                    df = df.drop(df_aaseq.index)
                    df = pd.concat([df, df_aaseq_ria], ignore_index=True)
                fn = ancres_list[0]
                fn_RIA = fn.replace(".txt", "_RIA.txt")
                fn_RIA_basename = "temp/" + os.path.basename(fn_RIA)
                self.add_ancresfn2dict(experiment, replicate, fn_RIA_basename)
                self.write_df2file(df, fn_RIA)

    def calcRIA_relint_v2(self):
        """
        taking relative intensities at TP0 
        into account for Light and Heavy calculations
        """
        ancres_dict = self.get_ancres_dict()
        for experiment in ancres_dict:
            for replicate in ancres_dict[experiment]:
                ancres_list = []
                ancres_list.append(self.get_ancresfnmostfilteredbyexprep(experiment, replicate))
                df = self.ancreslist2dflist(ancres_list)[0]
                aaseq_list = self.get_aaseqlistnoduplicates([df])
                for aaseq in aaseq_list:
                    index_aaseq = df[df["AAseq"] == aaseq].index.values             
                    return_list = self._calc_RIA_forAAseq_relint_v2(df, index_aaseq)
                    for index_name_val in return_list:
                        index_tp, name, val = index_name_val
                        df.loc[index_tp, name] = val
                fn = ancres_list[0]
                fn_RIA = fn.replace(".txt", "_RIA.txt")
                fn_RIA_basename = "temp/" + os.path.basename(fn_RIA)
                self.add_ancresfn2dict(experiment, replicate, fn_RIA_basename)
                self.write_df2file(df, fn_RIA)

    def _calc_RIA_forAAseq_relint_v2(self, df, index_aaseq):
        """
        DataFrame, ArrayOfIndices -> DataFrame
        calculates q ratio according to "Gustavsson et al., Proteomics,  2005"
        q = 15N / (14N + 15N)
        define which x15N are 14N by looking at TP min labeled,
        rest are assumed to be 15N.
        take relative intensities at TPmin into account.
        if MIP0_TP0 AND expMIPS_TPn exist --> rel. RIA
        else: RIA = -1.0
        returns DataFrame with relRIA or negative default
        """
        return_list = []
        dft = df.iloc[index_aaseq]   
        dft = dft.dropna(axis=1, how="all")
        tp_list = sorted(list(dft["TimePoint"]))
        tpmin = tp_list[0]
        expmips_tpmin = self._get_expmips(dft[dft["TimePoint"] == tpmin])    
        x14N_posandrelint_TP0 = self._calc_x14N_posandrelint(expmips_tpmin)
        for tp in tp_list:
            index_tp = dft[dft["TimePoint"] == tp].index.values
            expmips_ofTP = self._get_expmips(df.iloc[index_tp])
            if x14N_posandrelint_TP0 and expmips_ofTP:
                sum_x14N, sum_x15N = self._calc_relIntSum_HeavyLight(expmips_ofTP, x14N_posandrelint_TP0)
            else:
                sum_x14N = 0.0
                sum_x15N = 0.0
            if sum_x14N == 0.0:
                q = -1.0
            else:
                q = (sum_x15N / (sum_x14N + sum_x15N))
            return_list.append([index_tp, "RIA", q])
            return_list.append([index_tp, "sum14N", sum_x14N]) 
            return_list.append([index_tp, "sum15N", sum_x15N])         
        return return_list

    def _calc_RIA_forAAseq_relint(self, df):
        """
        DataFrame -> DataFrame
        DF(of AAseq)
        calculates q ratio according to "Gustavsson et al., Proteomics,  2005"
        q = 15N / (14N + 15N)
        define which x15N are 14N by looking at TP min labeled,
        rest are assumed to be 15N.
        take relative intensities at TPmin into account.
        if MIP0_TP0 AND expMIPS_TPn exist --> rel. RIA
        else: RIA = -1.0
        """
        df = df.dropna(axis=1, how="all")
        tp_list = sorted(list(df["TimePoint"]))
        tpmin = tp_list[0]
        df_tpmin_row = df[df["TimePoint"] == tpmin]
        expmips_tpmin = self._get_expmips(df_tpmin_row)
        x14N_posandrelint_TP0 = self._calc_x14N_posandrelint(expmips_tpmin)
        for tp in tp_list:
            df_tp = df[df["TimePoint"] == tp]
            expmips_ofTP = self._get_expmips(df_tp)
            if x14N_posandrelint_TP0 and expmips_ofTP:
                sum_x14N, sum_x15N = self._calc_relIntSum_HeavyLight(expmips_ofTP, x14N_posandrelint_TP0) # if MIP0 of TP_current not found mark RIA as -1
            else: # if MIP0 of TP0 not found mark all RIAs as -1 #!!!
                sum_x14N = 0.0
                sum_x15N = 0.0
            if sum_x14N == 0.0:
                q = -1.0
            else:
                q = (sum_x15N / (sum_x14N + sum_x15N))
            df.loc[df_tp.index, "RIA"] = q
            df.loc[df_tp.index, "sum14N"] = sum_x14N
            df.loc[df_tp.index, "sum15N"] = sum_x15N
        return df
        
    def _calc_x14N_posandrelint(self, expmips):
        """
        ExpMIPs -> Dict(key=x14N-position, val=relative Intensity to MIP0)
        x15N position and ratio of all peaks to MIP0
        """
        x14N_posandrelint = {}
        try:
            mip0_int = float(expmips[0][1])#!!! catch exception???
        except:
            return None
        for mzintx15n in expmips:
            ratio2MIP0 = mzintx15n[1]/mip0_int
            x14N_pos = mzintx15n[-1]
            x14N_posandrelint[x14N_pos] = ratio2MIP0
        return x14N_posandrelint

    def _calc_relIntSum_HeavyLight(self, expmips_ofTP, x14N_posandrelint_TP0):
        """
        produce sum_x14N and sum_x15N, taking rel. int. of TP0 into account
        """
        sum_x14N = 0.0
        sum_x15N = 0.0
        if expmips_ofTP[0][-1] == 0: # x15N-position of MIP0 has to be 0
            mip0 = expmips_ofTP[0][1]
            if mip0 <= 0.0:
                return(sum_x14N, sum_x15N)
        else:
            return(sum_x14N, sum_x15N)
        for mzintx15N in expmips_ofTP: # for all mz-positions 
            if x14N_posandrelint_TP0.has_key(mzintx15N[-1]): # if in 14N spectrum of TP0 AND heavy > 0 --> rel. RIA
                relint = x14N_posandrelint_TP0[mzintx15N[-1]]
                intensity = mzintx15N[1]
                light = mip0*relint
                heavy = intensity - light
                if heavy > 0:
                    sum_x14N += light
                    sum_x15N += heavy
                else:
                    sum_x14N += intensity #!!!
            else:                                            # else non rel. RIA
                sum_x15N += mzintx15N[1]
        return(sum_x14N, sum_x15N)
