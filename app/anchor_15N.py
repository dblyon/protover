import mfbt, math, numpy, peptide_results

class Anchor_15N():
    
    def __init__(self, templatemips, theomips, mz_arr, intensity_arr, Rt_spec, mz_ppm_picking_distance = 10.0, mz_ppm_picking_distance_mip0 = 10.0, mz_equality_threshold = 0.0000001, threshold_mip0int=None): #mz_equality_threshold = 0.0000001
        """
        anchor_15N.Anchor_15N(templatemips=templatemips, theomips=theomips, 
        raw_spectrum=spec, Rt_spec = Rt_spec, mz_ppm_picking_distance = mz_ppm_picking_distance, mz_ppm_picking_distance_mip0 = mz_ppm_picking_distance_mip0, threshold_mip0int=threshold_mip0int)

        mips_entry = [mz, int, ppm, x15N]
        templatemips: search for all mips in this list
        raw_spectrum: [[mz1, int1], [mz2, int2], ... ] # mgf_MS1 header + spectrum, iterable line by line
        rt: retention time in minutes
        mz_ppm_picking_distance: ppm deviation # means mz+/- picking_distance NOT mz+/- 1/2 picking_distance.
        mz_equality_threshold: mz_a == mz_b if abs(mz_a - mz_b) < mz_equality_threshold
        """
        self.escapee_list = []
        self.threshold_mip0int = threshold_mip0int
        self.mip0_ppm = mz_ppm_picking_distance_mip0
        self.theomips = theomips
        self.set_templatemipstheomz(templatemips)
        self.mz_equality_threshold = mz_equality_threshold
        self.mz_ppm_picking_distance = mz_ppm_picking_distance
        self.expmips = []
        self.rt = Rt_spec
        self.set_camel_score_negativedefault()     
        self.set_raw_data_filtered(self.reducearray_fromlowtohighmz(mz_arr, intensity_arr))
        self.set_raw_data_filtered(self.reducearray_intensitythreshold(self.raw_data_filtered[0], self.raw_data_filtered[1], intensity_threshold=0.0))
        self.set_raw_data_cruderange(numpy.copy(self.raw_data_filtered))
        self.set_escapee_utilized(False)

    def set_raw_data_cruderange(self, raw_data_cruderange):
        self.raw_data_cruderange = raw_data_cruderange
        
    def get_raw_data_cruderange(self):
        return self.raw_data_cruderange
        
    def set_raw_data_filtered(self, raw_data_filtered):
        self.raw_data_filtered = raw_data_filtered
        
    def get_raw_data_filtered(self):
        return self.raw_data_filtered
    
    def set_templatemipstheomz(self, templatemips):
        """
        if template_frompreviousTP True: 
        should only look for sticks/MIPs/x15N positions found in previous TP, 
        but use theoretical mz values not experimental mz values from previous TP.
        [mz_theo, int=100, ppm=0, x15N_ifinpreviousTP]
        set templatemips as atribute.
        """
        self.templatemips = []
        x15N_prevTP = [x[-1] for x in templatemips]
        for mips in self.theomips:
            if mips[-1] in x15N_prevTP:
                self.templatemips.append(mips)
                
    def reducearray_removemzintpair(self, mz_array, intensity_array, mz_rem):
        """
        Float, Array, Array -> Array, Array
        remove corresponding mz intensity pair from mz_array and intensity_array and return arrays.
        Method: produce bolean-conditional-array from mz_rem and mz_array, use on mz_array and intensity_array
        """
        cond = [ abs(mz_array - mz_rem) > self.mz_equality_threshold]
        mz_reduced = mz_array[cond]
        intensity_reduced = intensity_array[cond]
        return(mz_reduced, intensity_reduced)
        
    def reducearray_fromlowtohighmz(self, mz_array, intensity_array, mz_low=None, mz_high=None):
        """
        Array, Array, Float, Float -> Tuple(Array, Array)        
        Reduce both arrays to values within lower and upper boundary (Float, Float) of first array,
         if no mz_low/mz_high given defaults to MIP0-1.1 and MIPn+0.6
        """
        if not mz_low:
            mz_low = self.theomips[0][0] - 1.1 #!!!
            mz_high = self.theomips[-1][0] + 0.6
        cond = [ (mz_array > mz_low) & (mz_array < mz_high) ]
        mz_reduced = mz_array[cond]
        intensity_reduced = intensity_array[cond]
        return (mz_reduced, intensity_reduced) 
        
    def reducearray_tomaxintensity(self, mz_array, intensity_array):
        """
        Array, Array -> Float, Float
        Find max intensity, return corresponding mz and intensity vals as floats, else None, if len(array)==0, return None.
        ??? what if two intensity values are equal and max ???
        """
        if len(intensity_array) == 0:
            return None
        max_intensity = intensity_array.max()    
        cond = [( abs(max_intensity - intensity_array) < self.mz_equality_threshold) ]
        mz_array_red = mz_array[cond]
        intensity_array_red = intensity_array[cond]
        return(mz_array_red[0], intensity_array_red[0])

    def reducearray_intensitythreshold(self, mz_array, intensity_array, intensity_threshold=None):
        """
        Array, Array, Float -> Tuple(Array, Array)
        Reduce mz-and-int-arrays to values above intensity_threshold.
        """
        if intensity_threshold:
            intensity_threshold = intensity_threshold
        else:
            intensity_threshold = 0.0
        cond = [ (intensity_array > intensity_threshold) ]
        intensity_array_red = intensity_array[cond]
        mz_array_red = mz_array[cond]
        return(mz_array_red, intensity_array_red)        

    def find_maxintensity_withinmzrange(self, mz, mzint_arr_tuple, mz_ppm_picking_distance):
        """
        Float, Tuple(Array, Array), Float -> Float, Float
        Produce most intense mz value within mz_ppm_range, else None.
        """
        mz_lower_boundary, mz_upper_boundary = mfbt.Ppmrange(mz, mz_ppm_picking_distance).getppmrange()
        mzint_arr_tuple_red = self.reducearray_fromlowtohighmz(mzint_arr_tuple[0], mzint_arr_tuple[1], mz_lower_boundary, mz_upper_boundary)
        result = self.reducearray_tomaxintensity(mzint_arr_tuple_red[0], mzint_arr_tuple_red[1])
        if result:
            mz_ofintmax, int_max = result
        else:
            return None
        return(mz_ofintmax, int_max)
    
    def getexpmips(self):
        """
        returns all experimental MIPs found in self.raw_data_filtered or self.raw_data
        --> [[m/z, absolute intensity, delta ppm], ...]
        """
        return self.expmips

    def pick_expmips(self):
        """
        searches for all entries of self.templatemips in self.raw_data_filtered
        puts all found entries in self.expmips, ppmdeviation from template is calculated
        if self.picksticks(self.templatemips[0]) yields no result, self.expmips is empty
        """
        return_list = []
        return_list.append(self.pickstick(self.templatemips[0], ppm = self.mip0_ppm))
        if return_list[0]:
            if self.threshold_mip0int:
                mip0_intensity_threshold = (self.threshold_mip0int*return_list[0][1]) / 100.0 
                mz_arr[0], inte_arr[1] = self.get_raw_data_filtered()
                self.set_raw_data_filtered(self.reducearray_intensitythreshold(mz_arr, inte_arr, mip0_intensity_threshold))
            for mips in self.templatemips[1:]:
                mz_int_ppm = self.pickstick(mips)
                if mz_int_ppm:
                    return_list.append(mz_int_ppm)
            self.expmips = return_list[:]   

    def pickstick(self, mips, ppm=None): 
        """
        mips: a standard mips entry: (m/z, intensity, ppm_deviation, x15N).
        searches for an m/z entry in self.raw_data_filtered that's within +/- ppm_deviation of mz.
        if mip0_ppm given: searches for MIP0 with +/-mip0_ppm.
        return found m/z entry as a standard mips entry or None if none is found.
        """
        if ppm:
            mz_int_pair = self.find_maxintensity_withinmzrange(mips[0], self.get_raw_data_filtered(), ppm)
        else:
            mz_int_pair = self.find_maxintensity_withinmzrange(mips[0], self.get_raw_data_filtered(), self.mz_ppm_picking_distance)
        if mz_int_pair:
            delta_ppm = mfbt.Ppmprecision(mips[0], mz_int_pair[0]).getppmprecision()
            return [mz_int_pair[0], mz_int_pair[1], delta_ppm] + mips[3:]
        else:  
            return None
   
    def array_normtomax(self, mz_array, intensity_array):
        """
        Array, Array -> Array, Array
        spectrum: standard mips format
        finds maximum intensity, returns normalized copy of spectrum
        """
        mz_ofmaxint, intensity_maxint = self.reducearray_tomaxintensity(mz_array, intensity_array)
        if intensity_maxint == 0.0:
            return(mz_array, intensity_array)
        intensity_array_norm = (intensity_array*100)/intensity_maxint
        return(mz_array, intensity_array) 
        
    def normtomax(self, spectrum):
        """
        spectrum: standard mips format
        finds maximum intensity, returns normalized copy of spectrum
        """    
        returnlist = []
        try:
            max_mz=sorted(spectrum, key=lambda ele: int(ele[1]), reverse=True)[0][1]
        except IndexError:
            return spectrum
        for vector in spectrum:
            returnlist.append(vector[:1] + [(vector[1]/max_mz)*100.0] + vector[2:])
        return returnlist
        
    def getnormexpmips(self):
        """
        returns normalized self.expmips as standard mips list
        """
        return self.normtomax(self.expmips[:])
        
    def getnormcruderange(self):
        """
        returns raw_data_cruderange normalized to maximal intensity
        """
        mz_arr, inte_arr = self.get_raw_data_filtered()
        return self.array_normtomax(mr_arr, inte_arr)

    def is_n_times_more_intense(self, mipscurrent, mipsprev, threshold):
        """
        mipscurrent: standard mips entry 
        mipsprev: standard mips entry
        threshold: sane values are substantially greater than 1.0
        returns true if intensity of mipscurrent is more than threshold times bigger than mipsprev
        """
        return mipscurrent[1] / mipsprev[1] > threshold

    def iter_getmipsexp_escapee(self):
        """
        Repeatedly remove the first escapee in self.expmips from self.raw_data_filtered and repick self.expmips
        """
        while self.escapee_peak(self.expmips):
            self.pick_expmips()        
        
    def escapee_peak(self, mipslist, threshold=3.0):
        """
        mipslist: standard mips format
        threshold: escapee threshold - sane values are substantially greater than 1.0
        Compare adjacent entries in mipslist. If the value "to the right" is more than 
        threshold times bigger than its left neighbour, remove it from self.raw_data_filtered and return True.
        If no escapee is found in mipslist and deleted in self.raw_data_filtered, return False
        """
        for i in range(len(mipslist)-1):
            if self.is_n_times_more_intense(mipslist[i+1], mipslist[i], threshold=threshold): # #!!! SpeedUp: make lambda function
                self.add_escapee(mipslist[i+1])
                mz_arr, inte_arr = self.get_raw_data_filtered()
                self.set_raw_data_filtered(self.reducearray_removemzintpair(mz_arr, inte_arr, mipslist[i+1][0])) 
                self.set_escapee_utilized(True)
                return True
        return False
      
    def add_escapee(self, mips):
        """
        Self, Mips -> None
        add every mips entry corresponding to an escapee that was removed from self.raw_data_filtered
        """
        self.escapee_list.append(mips)
    
    def get_escapee_list(self):
        return self.escapee_list
    
    def remexpmipsbelowthreshold(self, threshold_maxmip):
        """
        finds maximum intensity of MIPs, 
        sets intensity threshold derived from maximum intensity and given threshold_maxmip
        removes all MIPs whose intensity is below threshold.
        """
        if self.expmips:
            maxint = max(self.expmips, key=lambda x: x[1])[1]
            maxint_threshold = (maxint*threshold_maxmip) / 100
            self.expmips = filter(lambda x: x[1] >= maxint_threshold, self.expmips)
            
    def find_local_min(self, mipslist):
        """
        Searches and returns the first local minimum intensity value including boundary values.
        Returns a standard mips entry if a local minimum is found, otherwise None.
        """
        for i in range(len(mipslist)-1):
            if mipslist[i][1] < mipslist[i+1][1]:
                return mipslist[i]
        if mipslist[-1][1] < mipslist[-2][1]:
            return mipslist[-1]
        else:
            return None
            
    def find_local_max(self, mipslist):
        """
        Searches and returns the first local maximum intensity value including boundary values.
        Returns a standard mips entry if a local maximum is found, otherwise None.
        """
        for i in range(len(mipslist)-1):
            if mipslist[i][1] > mipslist[i+1][1]:
                return mipslist[i]
        if mipslist[-1][1] > mipslist[-2][1]:
            return mipslist[-1]
        else:
            return None
                          
    def calc_leadcamel_score(self):
        """
        calculate score for camel sticks.
        1. Root Mean Square of ppm deviation = ppm_rms
        2. (Number of sticks found) / (number of theoretical sticks) = coverage
        3. total_score = coverage - sqrt(ppm_rms)/50 #total_score = coverage - ppm_rms  #(coverage / sqrt(ppm_rms)) * 100
        save score in self.camel_score as list --> NO should be dict!
        [ppm_rms, coverage, total_score]
        [ppm_rms, weighted_sumppm, logmip0intensity, rt, total_score]
        [coverage, ppm_rms, weighted_sumppm, weighted_sumintensities_log, self.rt, total_score]
        """
        mips_list=self.getexpmips()[:]
        sticks_found = float(len(mips_list))
        sticks_theo  = float(len(self.theomips))
        ppm_list=[]
        sum_pow_ppm=0
        for ele in mips_list:
            sum_pow_ppm += math.pow(float(ele[2]), 2)
        try:
            ppm_rms = math.sqrt((sum_pow_ppm/sticks_found))
            weighted_sumppm = self.weighted_sumbyintensity(2)
            weighted_sumintensities = self.weighted_sumbyintensity(1)
        except ZeroDivisionError:
            return None
        coverage =          sticks_found/sticks_theo
        coverage_template = sticks_found/float(len(self.templatemips))
        logmip0intensity = math.log(self.expmips[0][1])
        total_score = ( (coverage - (weighted_sumppm/100.0)) - (1 / logmip0intensity) )
        self.set_camel_score({"coverage": coverage, "coverage_template": coverage_template, "ppm_rms": ppm_rms, "logmip0intensity": logmip0intensity, "weighted_sumppm": weighted_sumppm, "rt": self.rt, "total_score": total_score})

    def calc_cameltrail_score(self):
        """
        calculate score for camel sticks.
        1. Root Mean Square of ppm deviation = ppm_rms
        2. (Number of sticks found) / (number of theoretical sticks) = coverage
        3. total_score = coverage - sqrt(ppm_rms)/50 #total_score = coverage - ppm_rms  #(coverage / sqrt(ppm_rms)) * 100
        save score in self.camel_score as list --> NO should be dict!
        [ppm_rms, coverage, total_score]
        [ppm_rms, logmip0intensity, total_score]
        [ppm_rms, weighted_sumppm, logmip0intensity, rt, total_score]
        [ppm_rms, weighted_sumppm, weighted_sumintensities_log, rt, total_score]
        [ppm_rms, weighted_sumppm, weighted_sumintensities_log, rt, total_score]
        [coverage, ppm_rms, weighted_sumppm, weighted_sumintensities_log, self.rt, total_score]
        """
        penalty = 0.0
        if not self.expmips:
            return None
        if self.is_mip0_firststick(intensity_threshold = 200) == False:
            penalty = 3
        x15N_required_to_score = [0, 1, 2]
        x15N_expmips = [x[-1] for x in self.expmips]
        for x15N in x15N_required_to_score:
            if x15N not in x15N_expmips:
                self.set_camel_score_negativedefault()
                return None
        mips_list=self.getexpmips()[:]
        sticks_found = float(len(mips_list))
        sticks_theo  = float(len(self.theomips))
        ppm_list=[]
        sum_pow_ppm=0
        for ele in mips_list:
            sum_pow_ppm += math.pow(float(ele[2]), 2)
            intensities = [x[1] for x in self.expmips]
        try:
            ppm_rms = math.sqrt((sum_pow_ppm/sticks_found))
            weighted_sumppm = self.weighted_sumbyintensity(2)
            weighted_sumintensities = self.weighted_sumbyintensity(1)            
        except ZeroDivisionError:
            return None
        logmip0intensity = math.log(self.expmips[0][1])
        coverage =          sticks_found/sticks_theo
        coverage_template = sticks_found/float(len(self.templatemips))
        total_score = (logmip0intensity - weighted_sumppm + coverage - penalty) 
        self.set_camel_score({"coverage": coverage, "coverage_template": coverage_template, "ppm_rms": ppm_rms, "logmip0intensity": logmip0intensity, "weighted_sumppm": weighted_sumppm, "rt": self.rt, "total_score": total_score})
        
    def is_mip0_firststick(self, intensity_threshold = 90):
        """
        calculates mz value by subtracting delta_mz (= MIP1 - MIP0 ) from MIP0 
        --> mz of MIP minus 1 then pickstick in raw_data,
        compares height of MIP minus 1 to MIP0 --> if below threshold (meaning probably just noise)
        then returns True
        if above threshold returns False
        """
        delta_theomip0to1 = (self.theomips[1][0] - self.theomips[0][0])
        mipminus1_mz_theo = (self.theomips[0][0] - delta_theomip0to1)
        lookup_mipminus1 = [mipminus1_mz_theo, 100.0, 0.0, -1]
        mipminus1 = self.pickstick(lookup_mipminus1, ppm=10)
        if mipminus1:
            intensity_threshold_absval = (self.expmips[0][1] / 100) * intensity_threshold 
            intensity_lower = self.expmips[0][1] - intensity_threshold_absval
            intensity_upper = self.expmips[0][1] + intensity_threshold_absval
            if mipminus1[1] > intensity_lower and mipminus1[1] < intensity_upper: 
                return False
            else:
                return True
        else:
            return True
            
    def set_camel_score_negativedefault(self):
        """
        set self.camel_score to negative default values.
        self.camel_score = {"ppm_rms": 0.0, "coverage": 0.0, "total_score": -1}
        returns None
        """
        self.camel_score = {"coverage": 0.0, "coverage_template": 0.0, "ppm_rms": 0.0, "logmip0intensity": 0.0, "weighted_sumppm": 0.0, "rt": self.rt, "total_score": -1.0}

    def weighted_sumbyintensity(self, mzintorppm):
        """
        mips_entry = [mz, int, ppm, x15N]
        calculates (intensity of MIP / summed intensities of MIPs in self.expmips) 
        and multiplies the outcome with ppm_deviation of MIP
        creating weighted_ppm_deviation of mip
        averages weighted_ppm_deviations
        returns weighted_ppm_deviation
        """
        intensities_list = [x[1] for x in self.expmips]
        intensities_sum = float(sum(intensities_list))
        weighted_ppm_deviations_list = []
        for mip in self.expmips:
            weighted_ppm_deviation_of_mip = (mip[1] / intensities_sum) * abs(mip[mzintorppm])
            weighted_ppm_deviations_list.append(weighted_ppm_deviation_of_mip)
        return sum(weighted_ppm_deviations_list)
    
    def set_camel_score(self, camel_score):
        self.camel_score = camel_score
    
    def get_camel_score(self):
        """
        return camel score
        --> [ppm_rms, coverage, total_score=(coverage/ppm_rms)*100]
        """
        return self.camel_score
    
    def set_escapee_utilized(self, bool):
        self.escapee_utilized = bool
    
    def get_escapee_utilized(self):
        return self.escapee_utilized
    
    def get_anchor_result(self):
        """
        returns an anchor_result object
        """
        return peptide_results.Anchor_result(self.getexpmips(), self.get_camel_score(), self.get_escapee_utilized())
        
    
