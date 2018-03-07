from __future__ import print_function
import re, mfbt, os, peptide_results
import numpy as np
from scipy import stats
            
def get_home_forR():
    """
    Self -> String
    produce full path to home directory exchanging "\\" with "/"
    """
    home = os.path.expanduser("~")
    return home.replace("\\", "/") 

def return_int(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    """
    in order to sort by "human standard" integers in text from
    fn are returned as integers --> sorting works as expected
    """
    return [return_int(c) for c in re.split("(\d+)", text)]

def normtomax(spectrum):
    """
    spectrum: standard mips format
    finds maximum intensity, returns normalized copy of spectrum
    """
    returnlist = []
    try:
        max_intensity = sorted(spectrum, key=lambda ele: int(ele[1]), reverse=True)[0][1]    # TODO: remove if next line is equivalent
    except IndexError:
        #print("CAMELTOOLKIT: IndexError")
        return spectrum
    for vector in spectrum:
        returnlist.append(vector[:1] + [(vector[1]/max_intensity)*100.0] + vector[2:])
    return returnlist
    
def addppm(value, ppm):
    """
    calculates ppmrange of given value and ppm,
    if negative ppm given returns lower boundary of ppmrange.
    if positive ppm given returns upper boundary of ppmrange.
    """
    if ppm < 0:
        return mfbt.Ppmrange(value, abs(ppm)).getppmrange()[0]
    if ppm > 0:
	    return mfbt.Ppmrange(value, abs(ppm)).getppmrange()[1]
        
def selpexfile2dict(selpex_fn_path, selpex_results_object, sep="\t"):
    """
    expects a selpex_fn (containing AAseq; charge; rt) and 
    a selpex_results_object (containing a selpex_results dict)
    creates a Peptide instance including rt as attribute
    returns None.
    """
    fh = open(selpex_fn_path, "r")
    for line in fh:
        line = line.replace(",", ".")
        aaseq, charge, rt, an = line.strip().split(sep)
        charge = int(charge)
        rt = float(rt)
        an = an.split(";")
        peptide_instance = mfbt.Peptide(aaseq, chargestate = charge)
        peptide_instance.rt = rt
        selpex_results_object.selpex_results[aaseq] = peptide_results.Peptide_results(aaseq, charge, rt, peptide_instance) # create instance of Peptide_results  # containing solely selpex infos, to be filled with anchor_results later
        selpex_results_object.selpex_results[aaseq].set_an(an)
        
def files2chronoslist(source_path):
    """
    make list of all filenames in source_path directory,
    sort in "human" way, and oldest=max.labeled is first then less labeled files,
    --> saves sorted list of filenames as attribute.
    --> saves long names (including absolute path) as attribute.
    TPmax to TPmin
    """
    fn_list = os.listdir(source_path)
    fn_list.sort(key = natural_keys, reverse=True)
    fn_chronoslist = fn_list[:]
    fn_chronoslist_l = []
    for fn in fn_chronoslist:
        if fn.endswith("mzML"):
            x = source_path + "\\" + fn
            fn_chronoslist_l.append(x)
    return fn_chronoslist_l
            
# def fn2fh(source_path): # not used, I think
    # """
    # expects a list of filenames
    # returns a list of sorted filehandles.
    # """
    # fh_list = []
    # fn_list = files2chronoslist(source_path)
    # for fn in fn_list:
        # fn_path = source_path + "\\" + fn
        # fh = open(fn_path, "r")
        # fh_list.append(fh)
    # return fh_list

# def fn2spec_dict(fn): # not used, I think
    # """
    # expects filename, opens file and parses spectra.
    # entire spectrum = header, mz-int, end ions
    # returns dictionary[rt] = [entire spectrum]
    # """
    # spec_dict = {}
    # fh = open(fn, "r")
    # first = True
    # spec = []
    # for line in fh:
        # if first == True:
            # if "RTINSECONDS=" in line:
                # rt = float(line.strip().replace("RTINSECONDS=", ""))
                # first = False
        # spec.append(line)
        # if line.strip() == "END IONS":
            # spec_dict[rt] = spec[:]
            # spec = []
            # first=True
    # fh.close()
    # return spec_dict

def camelcroppedspec2spec_dict(camelcroppedspec, spec_rt, spec_dict): # not used, I think
    """
    """
    spec_dict[spec_rt] = camelcroppedspec #self.spec_dict_col[spec_dict_name][spec_rt] = camelcroppedspec

def peptideresults2csv(filename, selpex_results_object):
    """
    e.g. for one peptide:
    AAseq, TP0, Score#1, Score#2, ..., Score#n
    AAseq, TP1, Score#1, Score#2, ..., Score#n
    ...
    AAseq, TPmax, Score#1, Score#2, ..., Score#n
    tab delimited. comma delimited within one score triplet.
    """
    fh = open(filename, "w")        
    for key in selpex_results_object.selpex_results:
        TP_sorted = selpex_results_object.selpex_results[key].get_anchor_results_keys_sorted()
        for tp in TP_sorted:
            csv_line = ""
            csv_line += key
            csv_line += "\t"+tp
            camel_score = selpex_results_object.selpex_results[key].anchor_results_dict[tp].get_camel_score()
            camel_score_string = str(camel_score).replace("[", "").replace("]", "")
            csv_line+="\t"+camel_score_string
            remaining_camel_scores = selpex_results_object.selpex_results[key].anchor_results_dict[tp].get_remaining_camel_scores()
            for ele in remaining_camel_scores:
                remscorestring = str(ele).replace("[", "").replace("]", "")
                csv_line+="\t"+remscorestring
            csv_line+="\n"
            fh.write(csv_line)
    fh.close()

def selpexevaluation2csv(filename, selpex_results_object):
    """
    itta shouldda looka likea thisa:
    AAseq	charge	Rt	mz_MIP0	AN(s)	MIP0_found_Yes/No	forotherTPs	Coverage[%]	forotherTPs	TotalScore_Yes/No_(No if -1)	forotherTPs	ria	forotherTPs
    tab delimited, comma separated, dot as decimal
    """
    sr = selpex_results_object.selpex_results
    fh = open(filename, "w")
    header = "AAseq"+"\t"+"charge"+"\t"+"Rt"+"\t"+"mz_MIP0"+"\t"+"AN(s)"
    first = True
    for aaseq in sr:
        pepres = sr[aaseq]
        tp_fn_short_maxtomin = pepres.get_anchor_results_keys_sorted()[::-1]
        line2write = ""
        line2write += aaseq + "\t" + str(pepres.get_charge()) + "\t" + str(pepres.get_rt()) + "\t" + str(pepres.get_theomips()[0][0])
        an_list = pepres.get_an()
        an_string = ""
        for an in an_list:
            an_string += an + ", "
        an_string = an_string[:-2] # remove last ", "
        line2write += "\t" + an_string
        MIP0_found_YesNo = ""
        Coverage = ""
        totalscoreYesNo = ""
        rias = ""
        MIP0_found_tp = ""
        Coverage_tp = ""
        TotalScore_tp = ""
        rias_tp = ""
        for tp in tp_fn_short_maxtomin:
            if first:
                MIP0_found_tp += str(tp) + "_" + "MIP0_found" + "\t"
                Coverage_tp += str(tp) + "_" + "Coverage" + "\t"
                TotalScore_tp += str(tp) + "_" + "TotalScore" + "\t"
                rias_tp += str(tp) + "_" + "rias" + "\t"
            if pepres.get_expmips(tp): # is True if expmips NOT an empty list (which means that MIP0 has to be there).
                MIP0_found_YesNo += "yes" + "\t"
            else:
                MIP0_found_YesNo += "no" + "\t"
            Coverage += str(pepres.get_coverage(tp)) + "\t"
            if pepres.get_total_score(tp) == -1:
                totalscoreYesNo += "no" + "\t"
            else:
                totalscoreYesNo += "yes" + "\t"
            rias += str(pepres.get_ria(tp)) + "\t"                                
        if first:
            header += "\t" + MIP0_found_tp + Coverage_tp + TotalScore_tp + rias_tp + "\n"
            fh.write(header)
        first = False
        line2write += "\t" + MIP0_found_YesNo + Coverage + totalscoreYesNo + rias + "\n"
        fh.write(line2write)
    fh.close()

def get_sumintensities(expmips):
    """
    Expmips -> Float
    produce sum of all intensities
    """
    mz_list = [x[1] for x in expmips]
    return sum(mz_list)

def get_RIA0(expmips):
    """
    Expmips -> Float[0,1]
    produce ratio of MIP0/rest-MIPs
    if rest-MIPs == 0, return MIP0
    """
    MIP0_mz = [x[0] for x in expmips if x[-1] == 0][0] # get mz-val of x15N == 0
    rest_mz = [x[0] for x in expmips if x[-1] != 0]
    sum_rest_mz = sum(rest_mz)
    if sum_rest_mz == 0:
        return MIP0_mz
    else:
        return MIP0_mz / sum_rest_mz
        
def ancres2file(filename, selpex_results_object):
    """
    AAseq	charge	Rt_SelPeX   ANs	TimePoint    Filename	
        rt_spec	RIA	#N mz_0x15N	mz_1x15N	int_0x15N	int_1x15N    
    if nothing found "n.a." written to cell (NOT a Zero)
    Selpex[AAseq, charge, rt_selpex], number of Ns, AN(s), one time per TP (separate files for each TP)
    [tp_fn_short, rt_spec, ria, total_score, coverage, ppm_rms, weighted_sumppm, weighted_sumintensities_log, logmip0intensity, mz_0, mz_1, ..., mz_n, int_1, int_2, ..., int_n, ]
    """
    fh = open(filename, "w")
    sr = selpex_results_object.selpex_results
    header = "AAseq"+"\t"+"charge"+"\t"+"Rt_SelPeX"+"\t"+"ANs"+"\t"+ "TimePoint"+"\t"+"Filename" +"\t"+ "rt_spec"+"\t"+"RIA"+"\t"+"#N"
    first = True
    maxNofdataset = 0
    line2write_list = []
    for aaseq in sr: # find max number of N
        pepres = sr[aaseq]
        numberofNs = len(pepres.get_theomips())
        if numberofNs > maxNofdataset:
            maxNofdataset = numberofNs
    mz_x15N = ""
    int_x15N= ""
    for i in range(0, maxNofdataset):
        mz_x15N += "mz_"+str(i)+"x15N"+"\t"
        int_x15N+= "int_"+str(i)+"x15N"+"\t"
    header += "\t" + mz_x15N + int_x15N + "\n" #"\t" + "MZ_0x15N_to_nx15N" + ("\t"*maxNofdataset) + "\t" + "Intensity_0x15N_to_nx15N" + ("\t"*maxNofdataset) +"\n"  
    for aaseq in sr:
        pepres = sr[aaseq]
        numberofNs = len(pepres.get_theomips())
        an_list = pepres.get_an()
        ANs = ""
        for an in an_list:
            ANs += an + ", "
        ANs = ANs[:-2] # remove last ", "
        tpfilebasename_list_chronological = pepres.get_tpfilebasename_list_chronological()
        line2write_begin = ""
        line2write_begin += aaseq + "\t" + str(pepres.get_charge()) + "\t" + str(pepres.get_rt()) + "\t" +  ANs ###
        for tpfilebasename in tpfilebasename_list_chronological:
            line2write = ""
            tp, filebasename = tpfilebasename
            tp = str(tp)
            line2write += line2write_begin + "\t" + tp + "\t" + filebasename ###
            expmips = pepres.get_expmips(filebasename)
            # if expmips: # is True if expmips NOT an empty list (which means that MIP0 has to be there).
            ancres = pepres.anchor_results_dict[filebasename]
            x15N = 0
            MZ_0x15N_to_nx15N = ""
            Intensity_0x15N_to_nx15N = ""
            expmips_dict = {} # key = x15N value = [mz, int]
            for mips in pepres.get_expmips(filebasename):
                expmips_dict[mips[-1]] = mips[0:2]
            for mips in pepres.get_theomips():
                if expmips_dict.has_key(mips[-1]):
                    MZ_0x15N_to_nx15N += str(expmips_dict[mips[-1]][0]) + "\t"
                    Intensity_0x15N_to_nx15N += str(expmips_dict[mips[-1]][1]) + "\t"
                else:
                    MZ_0x15N_to_nx15N += "NaN" + "\t" #no value --> "n.a." NOT Zero
                    Intensity_0x15N_to_nx15N += "NaN" + "\t"
            MZ_0x15N_to_nx15N = MZ_0x15N_to_nx15N.rstrip() # remove last "\t"
            Intensity_0x15N_to_nx15N = Intensity_0x15N_to_nx15N.rstrip()
            rt_spec = str(ancres.get_rt())
            ria = str(ancres.get_ria())
            logmip0intensity = str(ancres.get_logmip0intensity())
            line2write += "\t" +  rt_spec  + "\t" +  ria #+ "\t" +  total_score +  "\t" +  coverage +  "\t" +  ppm_rms +  "\t" +  weighted_sumppm + "\t" + logmip0intensity
            line2write += "\t" + str(numberofNs) + "\t" +  MZ_0x15N_to_nx15N + "\t" + ("\t"*(maxNofdataset - numberofNs)) + Intensity_0x15N_to_nx15N
            line2write = line2write.rstrip()
            # else:
                # rt_spec = 0.0
                # ria = -1
                # line2write += "\t" +  rt_spec  + "\t" +  ria
                # MZ_0x15N_to_nx15N
                # Intensity_0x15N_to_nx15N
                # line2write += "\t" + str(numberofNs) + "\t" +  MZ_0x15N_to_nx15N + "\t" + ("\t"*(maxNofdataset - numberofNs)) + Intensity_0x15N_to_nx15N
            line2write += "\n"
            line2write_list.append(line2write)  
    line2write_list.insert(0, header)
    for line in line2write_list:
        fh.write(line)
    fh.close()

