from __future__ import print_function
import anchor_15N, mzML2py_pyteomics, cameltoolkit
from time import clock
import pandas as pd
from os import path as os_path
#pd.options.mode.chained_assignment = None  # default='warn'
  
class Camelherder():
    
    def __init__(self):
        pass           
    
    def _get_expmips(self, aaseq):
        return self.dict_temp_expmips[aaseq]    
    
    def _set_df(self, df):
        self.df = df
    
    def _get_df(self):
        return self.df
    
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
    
    def _write_df2file(self, filename_ancres):
        """
        """
        df = self._get_df().sort_values(["AAseq", "TimePoint"], ascending=[True, True])
        col_list_sorted = ["AAseq", "charge", "Rt_input", "ANs", "TimePoint", "Filename", "Rt_spec", "RIA", "#N"]
        mz_x15N_list, int_x15N_list = self._sort_x15Ncolnames_natkeys(df)
        col_list_sorted += mz_x15N_list + int_x15N_list
        # df = df.reindex_axis(col_list_sorted, axis=1)
        df = df.reindex(col_list_sorted, axis=1)
        df["TimePoint"] = df["TimePoint"].astype('int32')
        df["#N"] = df["#N"].astype('int32')
        df.to_csv(filename_ancres, sep="\t", header=True, index=False)        
        
    def _create_df(self, selpex_results_object):
        maxNofdataset = self._maxNofdataset(selpex_results_object)
        dicti = {"AAseq": "", "charge": 23, "Rt_input": 12.34, "ANs": "ABC123", "TimePoint": 0, "Filename": "batman.mzML", "Rt_spec": 12.67, "RIA": 0.99, "#N": 123}
        for i in range(0, maxNofdataset):
            mz_x15N  = "mz_"  + str(i) + "x15N"
            int_x15N = "int_" + str(i) + "x15N"
            dicti[mz_x15N] = 0.0
            dicti[int_x15N] = 0.0
        df = pd.DataFrame(columns=dicti)
        self._set_df(df)
    
    def _maxNofdataset(self, selpex_results_object):
        """
        selpex_results_object -> Integer
        returns the maximum number of N of all given peptide sequences.
        """
        maxNofdataset = 0
        line2write_list = []
        for aaseq in selpex_results_object.selpex_results: 
            pepres = selpex_results_object.selpex_results[aaseq]
            numberofNs = len(pepres.get_theomips())
            if numberofNs > maxNofdataset:
                maxNofdataset = numberofNs
        return maxNofdataset

    def _ancres2dictrow(self, pepres, filebasename):
        """
        Anchor_results_object -> Dict(of one row of DataFrame output)
        """
        dict_row = {}
        dict_row["AAseq"] = pepres.get_aaseq()
        dict_row["charge"] = pepres.get_charge()
        dict_row["Rt_input"] = pepres.get_rt()
        an_list = pepres.get_an()
        ANs = ""
        for an in an_list:
            ANs += an + ", "
        ANs = ANs[:-2]
        dict_row["ANs"] = ANs
        for tp, fn in pepres.get_tpchronosreverse_dict().iteritems():
            if fn == filebasename:
                dict_row["TimePoint"] = tp
        dict_row["Filename"] = filebasename
        ancres = pepres.get_anchor_result_ofTP(filebasename)
        dict_row["Rt_spec"] = round(ancres.get_rt(), 3)
        dict_row["RIA"] = -1.0 #negative default
        dict_row["#N"] = len(pepres.get_theomips())
        expmips = ancres.get_expmips()
        for mips in expmips: #[mz, int, ppm, x15N]
            mz, int, ppm, x15N = mips
            mz_x15N  = "mz_"+str(x15N)+"x15N"
            int_x15N = "int_"+str(x15N)+"x15N"
            dict_row[mz_x15N] = float(mz)
            dict_row[int_x15N] = float(int)
        return dict_row
        
    def _add_anchor_result2df(self, pepres, filebasename):
        """
        """        
        dict_row = self._ancres2dictrow(pepres, filebasename)
        df = self._get_df().append(dict_row, ignore_index=True)
        self._set_df(df)
    
    def pickcamel(self, Rt_spec, mz_arr, intensity_arr, templatemips, theomips, leadcamel, threshold_mip0int=None, threshold_maxmip=None, mz_ppm_picking_distance_mip0 = 10.0, mz_ppm_picking_distance = 10.0, escapee=True):
        """
        expects Rt of spectrum and raw_data, as well as all other possible settings to pick expmips within anc object.
        returns anchor results object
        """
        anc = anchor_15N.Anchor_15N(templatemips=templatemips, theomips=theomips, mz_arr = mz_arr, intensity_arr = intensity_arr, Rt_spec = Rt_spec, mz_ppm_picking_distance = mz_ppm_picking_distance, mz_ppm_picking_distance_mip0 = mz_ppm_picking_distance_mip0, threshold_mip0int=threshold_mip0int)
        anc.pick_expmips()
        if escapee:
            anc.iter_getmipsexp_escapee()
        if threshold_maxmip:
            anc.remexpmipsbelowthreshold(threshold_maxmip=threshold_maxmip)
        if leadcamel == True:
            anc.calc_cameltrail_score()
        else:
            anc.calc_cameltrail_score()
        return anc.get_anchor_result()
        
    def selpexfile2resultsdict(self, selpex_fn_path, selpex_results_object, rt_range=5.0, dalton_range = 0.6, centroid=False):
        """
        creates entries in selpex_results for each AAseq (key = AAseq, value = peptide_results_object)
        also creates a rawdata_frame consisiting of (rt_low, rt_high, mz_low, mz_high, centroid)
        e.g. rawdata_frame = (35.89, 40.89, 500.12345, 520.12345, centroid=False)
        """
        cameltoolkit.selpexfile2dict(selpex_fn_path, selpex_results_object)
        for aaseq in selpex_results_object.selpex_results:
            pepres_object = selpex_results_object.selpex_results[aaseq]
            theomips = pepres_object.get_theomips()
            mz_low = theomips[0][0] - dalton_range
            mz_high = theomips[-1][0] + dalton_range
            Rt_input = pepres_object.get_rt()
            rt_low = Rt_input - rt_range
            rt_high = Rt_input + rt_range
            rawdata_frame = (rt_low, rt_high, mz_low, mz_high, centroid)
            pepres_object.set_rawdata_frame(rawdata_frame)

    def pickbestanchorresult(self, selpex_results_object, filebasename, dbl):
        """
        for every peptide in selpex_results
        sort temp_anchor_results by total score.
        set ancres with max score and clear temp_anchor_results_dict.
        """
        for aaseq in selpex_results_object.selpex_results:
            pepres = selpex_results_object.selpex_results[aaseq]
            self._add_anchor_result2df(pepres, filebasename)
            
    def findgap(self, expmips):
        """
        [mz, int, ppm, x15N]
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
        
    def herdthem_mzml(self, exp_dict_oneexperiment,
        selpex_fn_path, selpex_results_object, 
        rt_range=5, dalton_range=0.6, 
        threshold_mip0int=None, 
        threshold_maxmip=None, rt_radius_avg=15, 
        mz_ppm_picking_distance_mip0 = 10.0, mz_ppm_picking_distance = 10.0, 
        decimal = 3, escapee = True, dbl=False, filename_ancres=None):
        tps_reverse_chronos = exp_dict_oneexperiment.keys()
        tps_reverse_chronos.sort(reverse=True)
        fn_list_reverse = []
        for tp in tps_reverse_chronos:
            fn_list_reverse.append(exp_dict_oneexperiment[tp]) 
        self.fn_chronoslist_l = fn_list_reverse
        fn_chronoslist_l = self.fn_chronoslist_l
        leadcamel_fn = self.fn_chronoslist_l[0] 
        centroid = False
        self.selpexfile2resultsdict(selpex_fn_path, selpex_results_object, dalton_range = dalton_range, centroid = centroid, rt_range=rt_range)
        for key in exp_dict_oneexperiment.keys():
            exp_dict_oneexperiment[key] = os_path.basename(exp_dict_oneexperiment[key])
        rtlowhighaaseq_list = []
        for aaseq in selpex_results_object.selpex_results:
            selpex_results_object.selpex_results[aaseq].set_tpchronosreverse_dict(exp_dict_oneexperiment)
            pepres = selpex_results_object.selpex_results[aaseq]
            (rt_low, rt_high) = pepres.get_rtboundaries()
            rtlowhighaaseq_list.append([rt_low, rt_high, aaseq])
        rtlowhighaaseq_list = sorted(rtlowhighaaseq_list, key=lambda x: x[0])
        do_leadcamel = True
        self._create_df(selpex_results_object)        
        for index, filename_raw_data in enumerate(self.fn_chronoslist_l):
            start = clock()
            filebasename = os_path.basename(filename_raw_data)
            scannum2aaseq_dict = self.get_scannum2aaseq_dict(filename_raw_data, selpex_results_object, rtlowhighaaseq_list)
            mzml = mzML2py_pyteomics.Mzml(filename_raw_data)
            for rawspec in mzml: 
                spec = mzml.parserawspec(rawspec)
                if not spec:
                    continue
                (Rt_spec, scannumber, mz_arr, intensity_arr) = spec
                if scannum2aaseq_dict.has_key(scannumber):
                    aaseq_list = scannum2aaseq_dict[scannumber]
                    for aaseq in aaseq_list:
                        pepres = selpex_results_object.selpex_results[aaseq]
                        if do_leadcamel: 
                            tmips = pepres.get_theomips()
                        else:
                            filebasename_last = os_path.basename(self.fn_chronoslist_l[index-1])
                            tmips = pepres.get_expmips(filebasename_last)
                            if tmips == []: 
                                tmips = pepres.get_theomips()
                        ancres = self.pickcamel(Rt_spec, mz_arr = mz_arr, intensity_arr = intensity_arr, templatemips = tmips, theomips = pepres.get_theomips(), leadcamel = do_leadcamel,
                                                threshold_mip0int = threshold_mip0int, threshold_maxmip = threshold_maxmip, 
                                                mz_ppm_picking_distance_mip0 = mz_ppm_picking_distance_mip0, mz_ppm_picking_distance = mz_ppm_picking_distance, escapee=escapee)    
                        pepres.add_anchor_result(ancres, filebasename)
            self.pickbestanchorresult(selpex_results_object, filebasename, dbl)
            do_leadcamel = False
            stop = clock()
            runtime = float((stop - start)/60)
            print("\n", filebasename, "\n", "runtime[min]: %.3f"  % runtime)
        self._write_df2file(filename_ancres)        
        
    def get_scannum2aaseq_dict(self, filename_raw_data, selpex_results_object, rtlowhighaaseq_list):
        """
        FileName -> Dict(key = Rt_mzML, vals = [AAseq1, AAseq2, ...])
        """
        scannum2aaseq_dict = {}
        mzml = mzML2py_pyteomics.Mzml(filename_raw_data)
        scannum2rt_dict = mzml.get_rtscannumdict() #Dict(key = scannumber, val = Rt)
        temp_aaseq_list_pos = []
        for scannum in sorted(scannum2rt_dict.keys()):
            rt = scannum2rt_dict[scannum]
            aaseq_list = [rtlowhighaaseq[2] for rtlowhighaaseq in rtlowhighaaseq_list if rt > rtlowhighaaseq[0] and rt < rtlowhighaaseq[1]]
            if aaseq_list:
                scannum2aaseq_dict[scannum] = aaseq_list
                temp_aaseq_list_pos += aaseq_list
        temp_aaseq_list_pos = list(set(temp_aaseq_list_pos))
        for aaseq in selpex_results_object.selpex_results:
            pepres = selpex_results_object.selpex_results[aaseq]
            if aaseq not in temp_aaseq_list_pos:
                pepres.add_anchor_result(pepres._get_negative_default_ancres(), os_path.basename(filename_raw_data))
        return scannum2aaseq_dict
