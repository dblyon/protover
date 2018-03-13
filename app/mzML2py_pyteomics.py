from pyteomics import mzml as mzmlpomics

class Mzml():
    def __init__(self, fn):
        """
        """
        self.set_fn(fn)
        self.rawspectra = mzmlpomics.read(self.get_fn())
        
    def __iter__(self):
        return self.rawspectra
        
    def set_fn(self, fn):
        self.fn = fn
        
    def get_fn(self):
        return self.fn  

    def rewind_fh(self):
        self.rawspectra = mzmlpomics.read(self.get_fn())
    
    def parserawspec(self, rawspec):
        """
        Self -> Tuple(rt_spec, mz-array, intensity-array)
        produces Rt, mz_arr, intensity_arr of current spectrum if self.mslevel True, otherwise returns None
        """
        if rawspec["ms level"] == 1:          
           return(rawspec["scanList"]["scan"][0]["scan start time"], rawspec["index"], rawspec["m/z array"], rawspec["intensity array"])

    def get_rtscannumdict(self):
        """
        Self -> Dict(key = scannumber, val = Rt)
        List(of Float and Int)
        """
        rtscannum_dict = {}
        self.rewind_fh()
        for rawspec in self.rawspectra:
            if rawspec["ms level"] == 1:
                rtscannum_dict[rawspec["index"]] = rawspec["scanList"]["scan"][0]["scan start time"]
        return rtscannum_dict
        
    def get_rawspec_at_rt(self, rt):
        """
        (Rt_spec, scannumber, mz_arr, intensity_arr) = spec
        """
        self.rewind_fh()
        for rawspec in self.rawspectra:
            if rawspec["ms level"] == 1:
                spec = self.parserawspec(rawspec)
                (Rt_spec, scannumber, mz_arr, intensity_arr) = spec
                if Rt_spec >= rt:
                    return spec        
    
    def get_rawspec_at_scannum(self, scannum):
        """
        (Rt_spec, scannumber, mz_arr, intensity_arr) = spec
        """
        self.rewind_fh()
        for rawspec in self.rawspectra:
            if rawspec["ms level"] == 1:
                spec = self.parserawspec(rawspec)
                (Rt_spec, scannumber, mz_arr, intensity_arr) = spec
                if scannumber >= scannum:
                    return spec     
        
        
        
        
        
        