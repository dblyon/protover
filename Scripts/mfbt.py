from __future__ import print_function
from aa_reference import aa_reference
from isotope_reference import isotope_reference
import subprocess, shlex, sys, os
#from pprint import pprint
#import anchor_15N, py2r##
from time import clock
import numpy as np
import pandas as pd

class Proteincollection(object):

    def __init__(self):
        self.prot_col_dict = {} # key = AN, val = Protein-Object
        self.pep_dict = {}
    
    def readfile(self, filename, fastatype):
        """
        creates filehandle for fasta file, sets fastatype
        """
        self.set_fastatype(fastatype)
        self.fh_fasta=open(filename, "r")
        self.fastafilename = filename
    
    def get_fastafilename(self):
        return self.fastafilename
        
    def get_pcd(self):
        """
        return protein collection dictionary
        """
        return self.prot_col_dict

    def readentry(self):
        """
        reads self.fh_fasta line by line,
        removes "*", strips leading and trailing whitespaces, converts to upper case
        creates list with elements ["AA_seq", "header"]
        closes self.fh_fasta
        """
        liste = []
        seq = ""
        header = ""
        didfirst = False              
        for line in self.fh_fasta:       
            if line[0] == ">":
                if didfirst == True:
                    liste.append([seq, header])
                else:
                    didfirst = True
                header = line.rstrip()
                seq = ""
            else:
                seq = seq + line.replace("*", "").strip().upper()
        liste.append([seq, header])
        self.fh_fasta.close()
        return liste
    
    def fasta2simpledict(self, filename, fastatype, pepseq_list):
        """
        FileName, FastaType, List[AAseqs] -> Dict
        Produce a dictionary with AN as key and List[aaseq, header] as value, 
        only proteins containing pepseqs are included in the dict.
        parsing proteins from readentry_list
        creating temporary instances of the class protein
        adding only those entries to fasta_simple_dict that contain elements of pepseq_list
        fasta_simple_dict[AN] = [aaseq, header]
        returns fasta_simple_dict
        """
        fasta_simple_dict = {}
        self.readfile(filename, fastatype)
        fastaliste = []
        seq = ""
        header = ""
        didfirst = False              
        for line in self.fh_fasta:       
            if line[0] == ">":
                if didfirst == True:
                    for ele in pepseq_list:
                        if ele in seq:
                            fastaliste.append([seq, header])
                else:
                    didfirst = True
                header = line.rstrip()
                seq = ""
            else:
                seq = seq + line.replace("*", "").strip().upper()
        for ele in pepseq_list:
            if ele in seq:
                fastaliste.append([seq, header])
        self.fh_fasta.close()
        for ele in fastaliste:
            prot = Protein(ele[1], ele[0], self.get_fastatype())
            an = prot.get_an() # prot.getan(ele[1])
            aaseq = prot.get_sequence()
            header = ele[1]
            fasta_simple_dict[an] = [aaseq, header]
        return fasta_simple_dict 
        
    def pepseq2an(self, pepseq_list, fastafilename, fastatype):
        """
        List, FileName, FastaType -> Dict(key=pepseq, val=ANs)
        produces a dictionary with Peptide sequences as keys and AccessionNumbers as values in a list.
        expects a list of peptide sequences, a fastafile and fastatype.
        creates fasta_simple_dict[AN] = [aaseq, header]
        and searches entries in pepseq_list within aaseq.
        returns a dictionary with pepseq as key and accession number(s) as values
        """
        fasta_simple_dict = self.fasta2simpledict(fastafilename, fastatype, pepseq_list)
        pepseq2an_dict = {}
        for pepseq in pepseq_list:
            pepseq2an_dict[pepseq] = []
            for key in fasta_simple_dict:
                if pepseq in fasta_simple_dict[key][0]:
                    pepseq2an_dict[pepseq].append(key)
        return pepseq2an_dict
    
    def add_protein(self, protein):
        """
        expects a protein object
        adds it to self.pc 
        """
        an = protein.get_an()
        self.prot_col_dict[an] = protein # dict: key=an, value=protein_object
    
    def makeproteins(self, filename, fastatype): #remember in protein collection which proteins generated here?
        """
        parsing proteins from readentry_list
        creating instances of the class protein
        self.prot_col_dict[AN] = proteinObject
        returns None
        """
        self.readfile(filename, fastatype)
        fastalist = self.readentry()
        for ele in fastalist:
            prot = Protein(ele[1], ele[0], self.get_fastatype()) # header, AAseq, fastatype
            an = prot.get_an()
            self.prot_col_dict[an] = prot #dict: key=an, value=protein_object           
        
    def makepep_dict(self):
        """
        for each protein/entry in self.prot_col_dict
        for each cleaved peptide sequence
        create an entry in self.pep_dict using pepseq as key
        and accession number as value.
        return None
        """
        for an in self.prot_col_dict: #key=an value=protein object
            for pep in self.prot_col_dict[an].pep_list: # for every pep_seq in protein
                if not self.pep_dict.has_key(pep): # if dict doesn't have pepseq as key
                    self.pep_dict[pep] = [[an]]
                else:
                    self.pep_dict[pep].append(an)
    
    def set_fastatype(self, fastatype="uniprot"):
        self.fastatype = fastatype # legt fest welche parserfunktione verwendet wird
    
    def change_fastatypeforallproteins(self, fastatype):
        """
        sets fastatype for PC and for each protein in PC.
        """
        self.set_fastatype(fastatype)
        an_list = self.get_ans()
        for an in an_list:
            self.get_protein(an).set_fastatype(fastatype)
    
    def get_fastatype(self):
        return self.fastatype
    
    def get_ans(self):
        """
        returns a list with all AccessionNumbers
        """
        return self.prot_col_dict.keys()

    def get_protseq(self, an):
        """
        return protein sequence of an.
        """
        prot = self.get_protein(an)
        if prot:
            return prot.get_sequence()
        else:
            return None
    
    def get_header(self, an):
        """
        return header of an.
        """
        return self.prot_col_dict[an].get_header()
    
    def get_protein(self, an):
        """
        expects an AccessionNumber,
        returns a protein_instance of the corresponding an.
        """
        try:
            return self.prot_col_dict[an]
        except KeyError:
            return None
        
    def write2file(self, filename, linelength=60):
        pass

    def get_aaseqlistofallpeptides(self, missedcleavages=0, enzyme="trypsin"):
        """
        """
        aaseqs_list = []
        for an in self.get_pcd():
            prot = self.get_pcd()[an]
            aaseqs_list += prot.cutprotein(enzyme=enzyme, missedcleavages=missedcleavages)
        return aaseqs_list
     
    def get_mzlistofallpeptides(self, chargestate=2):
        """
        """
        mz_list = []
        aaseqs_list = self.get_aaseqlistofallpeptides()
        for aaseq in aaseqs_list:
            print(Peptide(aaseq, chargestate=chargestate).gettheomip(asfloat=True))
            break
            # mz_list += Peptide(aaseq, chargestate=chargestate).gettheomip(asfloat=True)
        return mz_list

    def write2file_subsetfasta_fromanslist(self, ans_list, fn_out, linelength=60):
        fh = open(fn_out, "w")
        for an in ans_list:
            header = self.get_header(an)  + "\n"
            aaseq  = self.get_protseq(an)
            fh.write(header)
            if linelength:
                while aaseq:
                    fh.write(aaseq[:linelength])
                    fh.write("\n")
                    aaseq = aaseq[linelength:]
            else:
                aaseq += "\n"
                fh.write(aaseq)
        fh.close()
        
class Protein(object):
    def __init__(self, header, aaseq, fastatype):
        self.pep_inst_list = []
        self.set_fastatype(fastatype)
        #self.header_dict={}
        self.set_sequence(aaseq)
        self.set_header(header)
        #self.addheader(header, fastatype)
    def get_pep_inst_list(self):
        return self.pep_inst_list
    
    def set_fastatype(self, fastatype):
        self.fastatype = fastatype
        
    def get_fastatype(self):
        return self.fastatype

    #def addheader(self, header, fastatype):
    #    """
    #    dict = {"accessionnumber": "header"}
    #    """
    #    an = self.getan(header)
    #    self.header_dict[an] = header
        
    def set_header(self, header):
        self.header = header
    
    def get_header(self):#, accessionnumber=None):
        return self.header
        #if accessionnumber==None:
        #    return self.header_dict.items()
        #else:
        #    return self.header_dict[accessionnumber]

    def set_sequence(self, aaseq):
        self.sequence = aaseq

    def get_sequence(self, linelength=0):
        return self.sequence
        
    def get_an(self):
        """
        glyme: >Glyma11g24934.1|PACid:26297255|PF05970|AT3G51690.1|PIF1 helicase
        -->Glyma11g24934.1
        """
        ft = self.get_fastatype()
        header = self.get_header()
        if ft.upper() == "UNIPROT": #>tr|Q5F4N8|Q5F4N8_TOBAC Putative ...
            an = header.split("|")[1]
        elif ft.upper() == "6FT_INCLFRAME":
            an = header[1:11]
        elif ft.upper() == "EVERYTHING":
            an = header.split(">")[1]    
        elif ft.upper() == "UREF100": #>UniRef100_A0EJF0 Beta-glucan-binding...
            an = header[11:17]
        elif ft.upper() == "MTGI": #>TC395598 homologue to UniRef100_A7P1E8 Cluster: Chromosome ...
            an = header.split(" ")[0]
            an = an[1:]
        elif ft.upper() == "MTGI_6FT":
            an = header[1:9]
            index = header.rfind("rframe")
            rframe = header[index:index+8]
            an = an + "_" + rframe
        elif ft.upper() == "IMGA": #>IMGA|contig_65682_1.1 Lipase contig_65682 902-62 ...
            an = header.split("|")[1].split(" ")[0]
        elif ft.upper() == "IMGANEW": #>Medtr1g004930.1 | hypothetical protein | LC | chr1:7332-689 | 20130731
            an = header.split("|")[0].split(">")[1].strip()
        elif ft.upper() == "MTGI_6FT_ANU": # AccessionNumberUnique --> see class SLH: def cat_6ftan # >BE204420_rframe-1 homologue to UniRef100_A9TVA3 Cluster: Predicted protein; n=5; Physcomitrella patens subsp. patens Rep: Predicted protein - Physcomitrella patens subsp. patens, partial (81%)_rframe-1
            an = header.split(" ")[0].replace(">", "")
        elif ft.upper() == "GLYMA": 
            an = header.split("|")[0][1:]
        elif ft.upper() == "MTGI_LSAN":
            an = header[1:9]
        elif ft.upper() == "COMBI": #>G7KQP3 | UniRef100_G7KQP3 Putative uncharact... __***__ Medtr7g089500.1| IMGA|Medtr7g089500.1 hy
            an = header.split("|")[0][1:].strip()
        elif ft.upper() == "UNIGENE_CONTIG":
            an = header.split("_")[-2]
            #an += "_" + header.split("_")[-1]
        elif ft.upper() == "UNIGENE_NCBI":
            an = header.split()[-1].split("_")[0]
        elif ft.upper() == "UNIGENE_MORE":
            temp = header.split()[1].split("_")[:-1]
            an = ""
            for ele in temp:
                an+= ele + "_"
            an = an[:-1]
        elif ft.upper() == "NTGI":
            an = header.split(" ")[0]
            an = an.split("_")[0]
            an = an[1:]
        elif ft.upper() == "TAIR": #>AT1G77850.1 | Symbols: ARF17 | auxin response f
            an = header.split("|")[0].split(".")[0][1:]
        elif ft.upper() == "UBIQUITIN":
            an = header.split("|")[0][1:]
        elif ft.upper() == "PHYSCO": # >Pp1s1_2V6.1 PEP --> Pp1s1_2V6.1
            an = header.split(" ")[0][1:]
        # elif ft.upper() == "EMBOSS6FT": 
            # # >TC1_1 UniRef100_A9TJ92 Cluster: Cellulose s --> "TC1_1 UniRef100_A9TJ92"
            # # >NP13121314_4 GB|XM_001752902.1|XP_001752954.1 hypothetical protein --> "NP13121314_4"
            # an = header.split
        elif ft.upper() == "GI": #>gi|53828153|gb|AAU94356.1| iron reductase [Pisum sativum] --> "AAU94356.1"
            an = header.split("|")[3]
        else:
            an=None
        return an
        
    #def getans(self):
    #    """
    #    returns all accession numbers of protein as list
    #    """
    #    return self.header_dict.keys()
    
    def cutprotein(self, enzyme="trypsin", missedcleavages=2):
        """
        Signature: protein_sequence, enzyme, missedcleavages --> List
        Purpose: digests AAseq into Peptides
        """
        protein_sequence = self.get_sequence()
        protease = Protease(protein_sequence, enzyme=enzyme, missedcleavages=missedcleavages)
        return protease.get_peplist()
        
    def make_peptides(self, enzyme="trypsin", missedcleavages=2):
        """
        Signature:
        Purpose: 
        """
        self.pep_list = self.cutprotein(enzyme=enzyme, missedcleavages=missedcleavages)
        for pepseq in self.pep_list:
            peptide = Peptide(pepseq)
            self.pep_inst_list.append(peptide)

class Peptide(object):
    def __init__(self, aaseq, protein=None, chargestate=2):
        self.aaseq = aaseq.upper().rstrip()
        self.labellingmanager = Labellingmanager()
        self.setchargestate(chargestate) #mod
        self.camelresults = []

    def getsequence(self):
        """
        returns aminoacid sequence of peptide
        in upper case one letter code
        """
        return self.aaseq
                    
    def setchargestate(self, chargestate):
        self.chargestate = chargestate
        
    def getchargestate(self):
        """
        returns charge state integer value
        """
        return self.chargestate

    def gettheomips(self):
        """
        returns a list, each element a list with
        [mz, relative Intensity = 100, ppm = 0, "XYZ"(x 15N)]
        """
        masslist = self.labellingmanager.gettheomips(self.aaseq)
        charge = int(self.chargestate)
        for element in masslist:
            element[0] = (element[0] + charge*1.00727)/charge
            element.insert(2, 0.0)
        return masslist

    def gettheomip(self, x15N=0, asfloat=False):
        """
        returns theoretically calculated
        MonoIsotopicPrecursor 0 x 15N
        -> String
        """
        if asfloat:
            return self.gettheomips()[x15N][0]
        else:
            return "%.5f" % self.gettheomips()[x15N][0]
    
    def spectra2filetable(self, filename):
        """
        creates tab delimited file of
        monoisotopic precursor - relative intensity -
        number of 15N atoms
        """
        fh = open(filename, "w")
        for element in self.gettheomips():
            fh.write(("""%.5f\t%.2f\t%s\n""") %
                     (float(element[0]), float(element[1]), str(element[2])))
        fh.close()

    def spectra2targetlist(self, filename):
        """
        creates file with monoisotopic precursors
        of 0 x 15N to sumN x 15N
        """
        fh=open(filename, "w")
        for element in self.gettheomips():
            fh.write(("""%.5f\n""") % (float(element[0])))
        fh.close()
        
# class Isotopomermasses():
    # """
    # calculate masses according to Martin etal JPR 2012
    # """
    # def __init__(self):
        # pass
    # def gettheomips(self, aaseq):
        # composition = self._getcomposition(aaseq)
        # return self._calculatespectra(composition)
    # def _getcomposition(self, aaseq):
        # atomicdict={"C":0, "H": -(len(aaseq)-1)*2, "N":0, "O": -(len(aaseq)-1), "S":0}
        # (len(aaseq)-1)
        # for c in aaseq:
            # for key in atomicdict:
                # atomicdict[key] += aa_reference[c]["comp"][key]
        # return atomicdict
    # def _calculatespectra(self, atomicdict):
        # """
        # """
        # MIP0 = comp_dict["C"]*12.0 + comp_dict["H"]*1.007825 + comp_dict["N"]*14.003074 + comp_dict["O"]*15.994915 + comp_dict["S"]*31.972072
        
class Labellingmanager(object):
    def __init__(self):
        pass

    def gettheomips(self, aaseq):
        composition = self._getcomposition(aaseq)
        return self._calculatespectra(composition)

    def _getcomposition(self, aaseq):
        """
        writes elemental composition from aaseq to atomicdict,
        returns atomicdict
        """
        atomicdict={"C":0, "H": -(len(aaseq)-1)*2, "N":0, "O": -(len(aaseq)-1), "S":0}
        (len(aaseq)-1)
        for c in aaseq:
            for key in atomicdict:
                atomicdict[key] += aa_reference[c]["comp"][key]
        return atomicdict

    def _calculatespectra(self, atomicdict):
        mass = 0.0
        resultslist=[]
        deltam = isotope_reference["N"]["15N"][1]-isotope_reference["N"]["14N"][1]
        probability = 100.0
        for c in atomicdict:
            mass += atomicdict[c]*isotope_reference[c][sorted(isotope_reference[c])[0]][1]            
        for i in range(0, atomicdict["N"]+1, 1):
            resultslist.append([mass+deltam*i, probability, i])#"%s" % (int(i))])
        return resultslist

class Protease(object): # exceptIons NOT included yet !!! (e.g. PK/PR)
    def __init__(self, aaseq, enzyme="trypsin", missedcleavages=2):
        self.set_enzyme(enzyme)
        self.set_missedcleavages(missedcleavages)
        self.aaseq = aaseq
        self.enzyme_properties = {
         "trypsin": {"cuttingsites": ["K", "R"], "exceptions_cleaves_not": ["KP", "RP", "CKD", "DKD", "CKH", "CKY", "CRK", "RRH", "RRR"], "exceptions_cleaves": ["WKP", "MRP"]},
         "pepsin": {"cuttingsites": ["F", "Y", "W", "L"], "exceptions_cleaves_not": [], "exceptions_cleaves": []}
         }
        self.cleaver()
        
    def set_missedcleavages(self, missedcleavages):
        self.missedcleavages = missedcleavages
        
    def get_missedcleavages(self):
        return self.missedcleavages
        
    def set_enzyme(self, enzyme):
        self.enzyme = enzyme
        
    def get_enzyme(self):
        return self.enzyme
        
    def cleaver(self):
        for cs in  self.enzyme_properties[self.get_enzyme()]["cuttingsites"]: # cs is K, R
            self.aaseq = self.aaseq.replace(cs, cs+"*")
        pep_list_temp = self.aaseq.split("*")
        pep_list_temp = filter(None, pep_list_temp) # remove " " entries from list, which can occur when two cutting sites next to each other
        self.pep_list=[]
        for ele in pep_list_temp:
            index = pep_list_temp.index(ele)
            mcpep = ""
            for i in range(0, self.get_missedcleavages() + 1):
                try:
                    mcpep = mcpep + pep_list_temp[index+i]
                    self.pep_list.append(mcpep)
                except IndexError:
                    pass
                    
    def get_peplist(self):
        return self.pep_list

class Ppmrange(object):
    def __init__(self, mz, ppm=5):
        """
        calculates +/- ppm range of given mz value
        default = 5 ppm
        returns list with MZlow and MZhigh
        """
        self.mz = float(mz)
        self.ppm = float(ppm)
        self.ppmrange = []
        self.ppmrange = self.calcppmrange()
    def calcppmrange(self):
        deltam = ((self.mz/(10.0**6))*self.ppm)
        self.ppmrange.append(self.mz - deltam)
        self.ppmrange.append(self.mz + deltam)
        return self.ppmrange
    def getppmrange(self):
        return self.ppmrange
    
class Ppmprecision(object):  
    def __init__(self, mz_theo, mz_exp):
        """
        calculates ppm-deviation of experimental
        to theoretical mz values,
        returns float with delta ppm
        """
        self.mz_theo = float(mz_theo)
        self.mz_exp = float(mz_exp)
        self.ppmprecision = self.calcppmprecision()
    def calcppmprecision(self):
        try:
            dppm = (((self.mz_exp - self.mz_theo)/self.mz_theo)*(10**6))
        except TypeError:
            return "ND"
        dppm_2decimals="""%.2f""" % float(dppm)
        return float(dppm_2decimals)
    def getppmprecision(self):
        return self.ppmprecision

class SLH(object): # Santas Little Helpers
    def __init__(self):
        pass
    
    def cat_6ftan(self, proteincollection, fastafnout, linelength=None):
        """
        expects a proteincollection and an output filename,
        concatenates the AccessionNumber with the rframe
        producing a "unique AN".
        e.g.>TC203013 similar to UniRef100_Q9M5L0 Cluster: 60S ribosomal protein L35; n=1; Euphorbia esula Rep: 60S ribosomal protein L35 - Euphorbia esula (Leafy spurge), complete_rframe-3
        >TC203013_rframe-3 similar to UniRef100_Q9M5L0 Cluster: 60S ribosomal protein L35; n=1; Euphorbia esula Rep: 60S ribosomal protein L35 - Euphorbia esula (Leafy spurge), complete_rframe-3
        writes fasta to file.
        returns None.
        """
        fh = open(fastafnout, "w")
        pc = proteincollection
        an_list = pc.get_ans()
        for an in an_list:
            header = pc.get_header(an)
            header = header.replace(an.split("_")[0], an) # maybe too specific?
            header += "\n"
            fh.write(header)
            aaseq = pc.get_protseq(an)
            if linelength:
                while aaseq:
                    fh.write(aaseq[:linelength])
                    fh.write("\n")
                    aaseq = aaseq[linelength:]
            else:
                aaseq += "\n"
                fh.write(aaseq)
        fh.close()
    
    def add_proteotypic2file(self, aaseqfn_in, aaseqfn_out, fastafn, fastatype):
        """
        Signature: Fasta-file, fasta-type, txt-file-first-column-AAseq -> txt-file-second-column-proteotypic
        Purpose: add a column to txt file, if proteotypic "yes" if not "AN1, AN2, ..."
        Stub: 
        Unit-tests: test with big csv containing headers and other columns.
        Template:
        """            
        fh_in = open(aaseqfn_in, "r")
        pepseq_list = [line.split("\t")[0].strip() for line in fh_in]
        pc = Proteincollection()
        pepseq2an_dict = pc.pepseq2an(pepseq_list, fastafn, fastatype) # returns a dictionary with pepseq as key and accession number(s) as values           
        fh_out = open(aaseqfn_out, "w")
        fh_in.seek(0)
        for line in fh_in:
            line2write = ""
            pepseq = line.split("\t")[0].strip()
            line2write += pepseq+"\t" #1. write pepseq from original first column
            val = pepseq2an_dict[pepseq]
            #2. write "yes"/nothing/"AN1, AN2, ..." 
            if pepseq == "": # empty
                line2write = line.strip()+"\t"
            elif val == []: # means that no AN found, --> ?header? or not a peptide? or empty?
                line2write = line.strip()+"\t"
            elif len(val) == 1: # is proteotypic
                line2write += "yes"+"\t"
            else: # associated with many ANs
                ans_string = ""
                for ele in val:
                    ans_string += ele+", "
                line2write += ans_string[:-2]+"\t"
            #3. write rest of original remaining columns
            #line2append = line.strip().replace(pepseq, "")
            #line2write += line2append+"\n"
            line2write += "\n"
            fh_out.write(line2write)
        fh_in.close()
        fh_out.close()
        return pepseq2an_dict # ToDo comment out
            
    #def clean_6ftfasta(self, proteincollection, fastafnout, linelenth=None):
    #    fh = open(fastafnout, "w")
    #    pc = proteincollection
    #    an_list = pc.get_ans()
    #    anshort_dict = {}
    #    for an in an_list:
    #        an_noframe = an.split("_")[0]
    #        anshort_dict[an_noframe] = []
    #    for anshort in anshort_dict:
    #        for an in an_list:
    #            if anshort in an:
    #                #### blabla ToDo
    #        
    #        header = pc.get_header(an)
    #        header = header.replace(an.split("_")[0], an) # maybe too specific?
    #        header += "\n"
    #        fh.write(header)
    #        aaseq = pc.get_protseq(an)
    #        if linelength:
    #            while aaseq:
    #                fh.write(aaseq[:linelength])
    #                fh.write("\n")
    #                aaseq = aaseq[linelength:]
    #        else:
    #            aaseq += "\n"
    #            fh.write(aaseq)
    #    fh.close()
    
    def get_mzabundancelist(self, fasta_fn, fasta_type, min_aaseq_length=6, missedcleavages=2, charges=[2, 3]):
        """
        Signature: Fasta_fn, fasta_type and restrictions -> Tuple(List, List)
        Purpose: produce list of mz values (monoisotopic precursors within restrictions)
        Stub:
        Template:
        for every protein in proteincollection
        digest the protein into peptides with Trypsin using max. 2 missed cleavages
        filter AAsequences < 6 
        calculate MIP0 for charge 2 and 3
        add to list
        """
        mz_list = []
        ans_larger10k = []
        pc = Proteincollection()
        pc.makeproteins(fasta_fn, fasta_type)
        ans_list = pc.get_ans()
        for an in ans_list: ## for every protein in proteincollection
            prot = pc.get_protein(an)
            #prot.make_peptides() ## digest the protein into peptides with Trypsin using max. 2 missed cleavages
            pep_list = prot.cutprotein()
            #for pep in prot.get_pep_inst_list(): ## filter AAsequences < 6
            for pepseq in pep_list:
                #if len(pep.getsequence()) >= 6: ## filter AAsequences < 6 
                if len(pepseq) >= 6:
                    pep = Peptide(pepseq)
                    mz_list.append(pep.gettheomip()) ## calculate MIP0 for charge 2 and 3
                    if float(pep.gettheomip()) > 10000:
                        ans_larger10k.append(an)
                    pep.setchargestate(3)
                    mz_list.append(pep.gettheomip())
                    if float(pep.gettheomip()) > 10000:
                        ans_larger10k.append(an)
        return(mz_list, ans_larger10k)
    
    def mzlist2file(self, mz_list, fn):
        fh = open(fn, "w")
        for mz in mz_list:
            line2write = mz+"\n"
            fh.write(line2write)
        fh.close()

    def _istrueallexpnames(self, df, expnames):
        """
        return True if all expnames present in Experiment column
        """
        for ele in expnames:
            if True == df["Experiment"].isin([ele]).any(): # boolean vector checking if one expname is in any row of Expname column of given pepseq, then if a single True is there
                continue
            else:
                return False
        return True

    def _getdataframeforgivenaaseq(self, df, aaseq):
        return df[df["Sequence"] == aaseq]
        
    def maxquantdf2dfallreps(self, df, reps_list):
        """
        DataFrame, List(replicate names list) -> DataFrame
        produce DF if AAseq present in all replicates of reps_list
        """
        aaseqs_noduplicates = df["Sequence"].drop_duplicates()
        aaseq_list_true = []
        for aaseq in aaseqs_noduplicates:
            df_foraaseq = self._getdataframeforgivenaaseq(df, aaseq)
            if self._istrueallexpnames(df_foraaseq, reps_list):
                aaseq_list_true.append(aaseq)
        return df[df["Sequence"].isin(aaseq_list_true)]

    def df2selpexdf(self, df, an=False):
        """
        DataFrame -> DataFrame
        produce DataFrame with AAseq, Charge, Rt, AN (Selpex/Input for CamelCropper)
        """
        aaseqs_noduplicates = df["Sequence"].drop_duplicates()
        dict_selpex = {}
        for index, aaseq in enumerate(aaseqs_noduplicates):
            df_foraaseq = self._getdataframeforgivenaaseq(df, aaseq)
            sequence = aaseq
            rt = np.mean(df_foraaseq["Calibrated retention time"])
            charge = int(round(np.mean(df_foraaseq["Charge"])))
            if an:
                an_lrp = df_foraaseq["Leading Razor Protein"].iloc[0]
               # an_lrp = x.ix[x.index[0]].split("|")[0][1:]
                dict_selpex[index] = {"Sequence": sequence, "Charge": charge, "Rt": rt, "AN": an_lrp}
            else:
                dict_selpex[index] = {"Sequence": sequence, "Charge": charge, "Rt": rt}
        df_selpex = pd.DataFrame(dict_selpex)
        df_selpex = df_selpex.transpose()
        return df_selpex
    
    def remove_NAN_intensity_score_PEP(self, df):
        """
        DataFrame -> DataFrame
        removed all rows containing NANs in Intensity, Score, and PEP column
        """
        print(len(df))
        df = df[np.isfinite(df["Intensity"])] # remove NaN
        print(len(df))
        df = df[np.isfinite(df["Score"])] # remove NaN
        print(len(df))
        df = df[np.isfinite(df["PEP"])] # remove NaN
        print(len(df))
        return df
    
    def remove_NAN_reverse_contaminant_modification(self, df):
        """
        DataFrame -> DataFrame
        removed all rows containing NANs in Intensity, Score, and PEP column
        """
        print(len(df))
        df = df[df["Reverse"] != "+"]
        print(len(df))
        df = df[df["Contaminant"] != "+"]
        print(len(df))
        df = df[df["Modifications"] == "Unmodified"]
        print(len(df))
        return df
        
if __name__ == "__main__":
    slh = SLH()
##################### MaxQuant -> Selpex list    
    # fn=r"C:\Users\English\Documents\Uni\Medicago\15N\LCMS\20121105_AscEx_15N\MaxQuant\leaf14Ncontrol_techrepssameexp_matchbetruns_2PepsPerProt_CombiFASTA20131127\combined\txt\evidence.txt"
    fn=r"C:\Users\English\Documents\Uni\Medicago\15N\LCMS\20121105_AscEx_15N\MaxQuant\leaf14Ncontrol_techrepssameexp_matchbetruns_2PepsPerProt_CombiFASTA20131127_2\combined\txt\evidence.txt"
    df = pd.read_csv(fn, sep="\t")
    print("len DataFrame: ", len(df))
    df = slh.remove_NAN_reverse_contaminant_modification(df)
    print("len DataFrame: ", len(df))
    # df = slh.remove_NAN_intensity_score_PEP(df)
    # print("len DataFrame: ", len(df))
    
############## AAseq has to occur in all TPs of one replicate, sum of all AASeqs in Selpex-list
    # reps_list_a = ['0_10C14NLa_14', '24_11C14NLa_5', '48_12C14NLa_18', '72_13C14NLa_17', '96_14C14NLa_4']
    # reps_list_b = ['0_10C14NLb_29', '48_12C14NLb_27', '24_11C14NLb_30', '72_13C14NLb_39a', '96_14C14NLb_36']
    # reps_list_c = ['0_10C14NLc_43', '24_11C14NLc_57', '48_12C14NLc_63a', '72_13C14NLc_48_121111232214', '96_14C14NLc_53']
    # df_allreps_a = slh.maxquantdf2dfallreps(df, reps_list_a)
    # df_allreps_selpex_a = slh.df2selpexdf(df_allreps_a, an=True)
    # print("DataFrame after filtering AAseq in all replicates: ", len(df_allreps_selpex_a))
    # df_allreps_b = slh.maxquantdf2dfallreps(df, reps_list_b)
    # df_allreps_selpex_b = slh.df2selpexdf(df_allreps_b, an=True)
    # print("DataFrame after filtering AAseq in all replicates: ", len(df_allreps_selpex_a))
    # df_allreps_c = slh.maxquantdf2dfallreps(df, reps_list_c)
    # df_allreps_selpex_c = slh.df2selpexdf(df_allreps_c, an=True)
    # print("DataFrame after filtering AAseq in all replicates: ", len(df_allreps_selpex_a))    
    # list_a =  list(df_allreps_selpex_a["Sequence"])
    # list_b =  list(df_allreps_selpex_b["Sequence"])
    # list_c =  list(df_allreps_selpex_c["Sequence"])
    # list_abc = list_a + list_b + list_c
    # print("Len listABC redundant: ", len(list_abc))
    # list_abc_nonred = list(set(list_abc))
    # print("Len listABC non-redundant: ", len(list_abc_nonred))
    # df_list = [df_allreps_selpex_a, df_allreps_selpex_b, df_allreps_selpex_c]
    # concat = pd.concat(df_list, axis=0)
    # df_selpex_pepinalltpssumrep = concat.drop_duplicates(cols="Sequence")
    # num_aaseqs = len(df_selpex_pepinalltpssumrep)
    # print("Number of AAseqs in Selpex: ", num_aaseqs)
    # num_ans = len(list(df_selpex_pepinalltpssumrep["AN"].drop_duplicates()))
    # print("Number of AccessionNumbers: ", num_ans)
    # fn_out_selpex = fn.replace("evidence.txt", "MQselpex_PepinallTPsofRep.txt") # Write to file
    # df_selpex_pepinalltpssumrep = pd.DataFrame(df_selpex_pepinalltpssumrep, columns=["Sequence", "Charge", "Rt", "AN"])
    # df_selpex_pepinalltpssumrep.to_csv(fn_out_selpex, sep="\t", header=False, index=False)
### NOT restricted to occurrence of AAseq in all TPs
    df = slh.df2selpexdf(df, an=True)    
    num_aaseqs = len(df)
    print("Number of AAseqs in Selpex: ", num_aaseqs)
    num_ans = len(list(df["AN"].drop_duplicates()))
    print("Number of AccessionNumbers: ", num_ans)
    fn_out_selpex = fn.replace("evidence.txt", "MQselpex_PepinallTPsofRep.txt") # Write to file
    df = pd.DataFrame(df, columns=["Sequence", "Charge", "Rt", "AN"])
    df.to_csv(fn_out_selpex, sep="\t", header=False, index=False)
##############



    
# ############################### MTGI 6 frame translated fasta --> unique AN fasta
# ############################### cat_6ftan: from a sixframetranslated fasta --> create unique AN fasta.
    # home = os.path.expanduser("~")
    # fasta_fn = ("%s"+"%s") % (home, "\\Documents\\Uni\\fasta_mapping\\Medicago\\origin\\MTGI_032511_6frames.fasta")
    # fastatype = "mtgi_6ft"
    # pc_mt_gi6ft = Proteincollection()
    # pc_mt_gi6ft.makeproteins(fasta_fn, fastatype)
    # slh = SLH()
    # fastafnout = fasta_fn.replace(".fasta", "_ANunique.fasta")
    # slh.cat_6ftan(pc_mt_gi6ft, fastafnout, linelength=80)
    
    
################################ for Christiana Staudinger. mz abundance, fasta and restrictions given -> calculate mz values and produce list 
################################ used for plots: histogram mz vs abundance.
    # start = clock()
    # slh = SLH()
    # #fasta_fn = r"C:\Users\English\Documents\Magic Briefcase\ProteotypicPeptides_ChristianaStaudinger\UniR-taxonomy3880+AND+identity1.0.fasta"
    # #fasta_type = "uref100"
    # fasta_fn= r"C:\Users\English\Documents\Uni\fasta_mapping\Medicago\output\mt_uniprot_imga_mtgi_20130513.fasta"
    # #fasta_fn= r"C:\Users\English\Documents\Uni\fasta_mapping\Medicago\output\test.fasta"
    # fasta_type = "combi"
    # mz_list, ans_larger10k = slh.get_mzabundancelist(fasta_fn, fasta_type)
    # mz_list_fn = r"C:\Users\English\Documents\Magic Briefcase\ProteotypicPeptides_ChristianaStaudinger\mz_list_mt_uniprot_imga_mtgi_20130513.txt"
    # slh.mzlist2file(mz_list, mz_list_fn)
    # mz_list_arr = np.array(mz_list, dtype=float32)
    # ans_larger10k_arr = np.array(ans_larger10k, dtype=float32)
    # np.save("mz_list_arr", mz_list_arr)
    # np.save("ans_larger10k_arr", ans_larger10k_arr)
    # stop = clock()
    # print("runtime[min]: ", (stop-start)/60)
################################
################################ for Christiana Staudinger. add information to AAseq txt file if peptide proteotypic or not.
################################ if proteotypic add "yes" if not "AN1, AN2, ..."
    #start = clock()
    #slh = SLH()
    #aaseq_fn_in = r"C:\Users\English\Documents\Magic Briefcase\ProteotypicPeptides_ChristianaStaudinger\put-proteotypic-pep-Archetype130619_C10D11-14.txt"
    ##aaseq_fn_out = r"C:\Users\English\Documents\Magic Briefcase\ProteotypicPeptides_ChristianaStaudinger\put-proteotypic-pep-Archetype130619_C10D11-14_proteotypic.txt"
    #aaseq_fn_out = r"C:\Users\English\Documents\Magic Briefcase\ProteotypicPeptides_ChristianaStaudinger\put-proteotypic-pep-Archetype130619_C10D11-14_proteotypic_UP_IMGA_MTGI.txt"
    ##fasta_fn = r"C:\Users\English\Documents\Magic Briefcase\ProteotypicPeptides_ChristianaStaudinger\UniR-taxonomy3880+AND+identity1.0.fasta"    
    ##fasta_type = "uref100"    
    #fasta_fn= r"C:\Users\English\Documents\Uni\fasta_mapping\Medicago\output\mt_uniprot_imga_mtgi_20130513.fasta"
    #fasta_type = "combi"
    #x=slh.add_proteotypic2file(aaseq_fn_in, aaseq_fn_out, fasta_fn, fasta_type)
    #stop = clock()
    #print("runtime[min]: ", (stop-start)/60)
############################################    
    
    
##    #set fastafile and make protein collection
##    #fasta_fn = r"C:\Users\English\Documents\Uni\Medicago\fasta\uniprot\uniprot-organism%3AMedicago.fasta"
##    fasta_fn = r"C:\Users\English\Documents\Uni\Medicago\fasta\MTGI_11_LORFperFrame.fasta"
##    #fasta_fn = r"C:\Users\English\Documents\Uni\Medicago\fasta\mtgi_6ft_test.txt"
##    
####    fn="uniprot-organism%3A%22Sinorhizobium+medicae+%28strain+WSM419%29+%5B366394%5D%22.fasta"
####    fn_directory="C:\\bin\\fasta\\original\\"
####    fn_fasta="%s%s" % (fn_directory, fn)
##    start = clock()
##    pc=Proteincollection()
##    pc.readfile(fasta_fn, "mtgi_6ft")
##    pc.makeproteins()
##    stop = clock()
##    print("runtime: ", stop-start)
####    for an in pc.prot_col_dict: #key=an value=protein object
####        for pep in pc.prot_col_dict[an].pep_list: # for every pep_seq in protein
####            if not pc.pep_dict.has_key(pep): # if dict doesn't have key 
####                pc.pep_dict[pep] = [[an]]
####            else:
####                pc.pep_dict[pep].append(an)

    # # cat_6ftan: from a sixframetranslated fasta --> create unique AN fasta.
    # home = os.path.expanduser("~")
    # fasta_fn = ("%s"+"%s") % (home, "\\Documents\\Magic Briefcase\\fasta\\MTGI_11_LORFperFrame.fasta")
    # fastatype = "mtgi_6ft"
    # pc_mt_gi6ft = Proteincollection()
    # pc_mt_gi6ft.makeproteins(fasta_fn, fastatype)
    # slh = SLH()
    # fastafnout = ("%s"+"%s") % (home, "\\Documents\\Magic Briefcase\\fasta\\MTGI_11_LORFperFrame_ANunique.fasta")
    # slh.cat_6ftan(pc_mt_gi6ft, fastafnout, linelength=80)
    
    #fasta_fn = r"C:\Users\English\Documents\Uni\Medicago\fasta\uniprot\uniprot_20130327_MedicagoORSinorhizobiumMedicaeWSM419ORRhizobiumMeliloti1021.fasta"
    #fastatype = "uniprot"
    ##fasta_fn = r"C:\Users\english\Documents\Magic Briefcase\MTGI_11_LORFperFrame.fasta"
    ##fastatype = "mtgi_6ft"
    #pep_list = ['IFVGNLPFDVDSEK', 'STLTDSLVAAAGIIAQEVAGDVR', 'VVDLLAPYQR', 'LPLFGATDSSQVLK', 'VLYEVFPMSFLMEQAGGQAFTGK', 'VALVYGQMNEPPGAR', 'LTSVFGGAAEPPK', 'YDSTLGIFDADVKPVGTDGISVDGK', 'DSPLDVIAINDTGGVK', 'IIGVSVDSSGKPALR', 'QYADAVIEVLPTQLIPDDNEGK', 'VVDLADIVANNWK', 'WSPELAAACEVWK', 'VPTPNVSVVDLVVQVSK', 'RLVYTNDAGEVVK', 'LPLFGATDASQVLK', 'IIGFDNVR', 'VVGAFLEGGTPDENK', 'IEISGQNSWVGK', 'ALPTYTPETPADATR', 'IGLFGGAGVGK', 'IKPSDPSELESLLGAK', 'LLEATGISTVPGSGFGQK', 'GLAYDISDDQQDITR', 'VPVFLDGGVR', 'SGDDMTSLKDYVTR', 'FAQVTNPAIDPLR', 'TISKPVEHPSELPK', 'KFETLSYLPPLTR', 'GLLSFPFDGAYSTDFIEEWVK', 'YLEGAALGDANQDAIK', 'APITQQLPGESDTDFADFSSK', 'GILAADESTGTIGKR', 'IVNTGTVLQVGDGIAR', 'VGGPPAPSGGLPGTLNSDEAR', 'TIAMDATEGVVR', 'GLVPLAGSNDESWCQGLDGLASR', 'IGGIGTVPVGR', 'GGLDFTKDDENVNSQPFMR', 'LWGENFFDPATK', 'FYGEVTQQMLK', 'IPTDTQPGLLR', 'VVNALAKPIDGR', 'NNAGFPHNVIFDEDEIPSGVDAAK', 'LGNDISDLTLSQAK', 'IISGPAEKPLIGVNYK', 'ALDSQIAALSQDIVNK', 'FAPDANSQIVPASAIPDGWMGLDIGPDSIK', 'GHYLNATAGTCEDMMK', 'SYGVLIPDQGIALR', 'DGFEYITLR', 'MGINPIMMSAGELESGNAGEPAK', 'RLTFDEIQSK', 'FETLSYLPPLTEDQLAK', 'IVDVLIEQNIVPGIK', 'FETLSYLPPLTR', 'VIGQDEAVEAISR', 'ATGIIAQIPVSEGYLGR', 'KFETLSYLPPLTEDQLAK', 'GASTGYDNAVALPAGGR']
#    pc = Proteincollection()
#    #pc.readfile()
#    pc.makeproteins(fasta_fn, fastatype) # --> pc.prot_col_dict
#


##    #pick an accessionnumber
##    print(pc.prot_col_dict.keys()[:9])
##    k=pc.prot_col_dict.keys()[4] 
##    print("\n"+k+"\n")
##    
##    #show fastatype
##    print(pc.prot_col_dict[k].getfastatype(), "\n")
##
##    #get all and then specific header
##    print(pc.prot_col_dict[k].getheader(), "\n")
##    print(pc.prot_col_dict[k].getheader(k), "\n")
##
##    #print aaseq
##    print(pc.prot_col_dict[k].getsequence(), "\n")
##    aaseq = pc.prot_col_dict[k].getsequence()
##
##    #show proteasome/tryptic peptides of accessionnumber
##    print(pc.prot_col_dict[k].pep_list)
##
##    x=pc.prot_col_dict[k].pep_inst_list

