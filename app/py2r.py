import subprocess, os

class Talkr():
    """
    communicates with R via subprocess / piping
    predefined plots can easily be "filled" with experimentally derived test_data
    """
    def __init__(self, arial=False):
        self.type = "h"
        self.plotfilelocation = """./"""
        self.plotfilename = "plotfilename_temp.pdf"
        if arial:
            with open(self.plotfilename, "a") as fh:
                fh.write("library(extrafont)\nloadfonts()")
            fh.close()
        self.plottitle= "plottitle_temp"
        self.xlab = "m/z"
        self.ylab = "Intensity"
          
    def replace_orig_short(self, fn_orig, fn_short):
        fn_orig_newname = fn_orig.replace(".txt", "_orig.txt")
        if os.path.isfile(fn_orig_newname):
            os.remove(fn_orig_newname)
        os.rename(fn_orig, fn_orig_newname)
        os.rename(fn_short, fn_orig)
    
    def concat100(self, liste):
        """
        List -> (List, String)
        produce Tuple of String with <=100 characters, and List with elements removed that are within String
        """
        line2write = ""
        for index, ele in enumerate(liste):
            if len(line2write) <= 100:
                line2write +=  ele + ","
            else:
                liste_red = liste[index:]
                line2write += "\n"
                return(liste_red, line2write)
        liste_red = []
        if line2write.endswith(","):
            line2write = line2write[:-1]
        return(liste_red, line2write)
    
    def splitlineandwrite(self, line, fh_short):
        """
        String -> None
        Cut line into <100 character bits and write them to file
        """
        line_list = line.split(",")
        while line_list:
            line_list, line2write = self.concat100(line_list)
            fh_short.write(line2write.lstrip()) #!!!
    
    def shortline(self, fn):
        """
        FileName -> None
        rename input file to "_orig.txt" and
        create new file with same name as fn and same content
        but split lines >100 characters long to new line.
        """
        fh_orig = open(fn, "r")
        fn_short = fn.replace(".txt", "_short.txt")
        fh_short = open(fn_short, "w")
        for line in fh_orig:
            if len(line) <= 80:
                fh_short.write(line)
            else:
                self.splitlineandwrite(line, fh_short)
        fh_orig.close()
        fh_short.close()
        self.replace_orig_short(fn, fn_short)
        
    def run_rbatchtxtfile(self, fn, winorunix):
        text = """R CMD BATCH %s""" % fn
        if winorunix == "win":
            try:
                ph=subprocess.Popen(text, stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
            except:
                print("R is not installed, thus pdf plots will not be plotted.")
                return 1
        else:
            try:
                ph=subprocess.Popen(text, stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, shell=True)
            except:
                print("R is not installed, thus pdf plots will not be plotted.")
                return 1
        for line in ph.stdout:
            print line,
        for line in ph.stderr:
            print line,
        return ph.wait()
            
    def set_r_plottemptxt(self, fn):
        self.r_plottemptxt = fn
    
    def get_r_plottemptxt(self):
        return self.r_plottemptxt
  
    def plot_riaincreasing_allpeps_density(self, plotfilelocation, plottitle, plotfilename, num_true_abs, num_false_abs, num_true_perc, num_false_perc, ylab):
        """
        """
        self.setplotfilelocation(plotfilelocation)
        string2r = ("""
dataraw <- test_data.frame(row.names = c("%s"),  num_true_abs=c(%s), num_false_abs=c(%s))
dataperc <- test_data.frame(row.names = c("%s"), true=c(%s),  false=c(%s))
pdf("%s%s%s", paper="a4r")
x <- barplot(t(as.matrix(dataperc)), col=c("darkgreen", "darkred"), legend=TRUE, border=NA, xlim=c(0,8), args.legend=list(bty="n", border=NA, x=3), ylab="%s")
text(x, dataperc$true-3, labels=round(dataraw$num_true_abs), col="black")
text(x, dataperc$true+3, labels=round(dataraw$num_false_abs))
dev.off()

""" % (plottitle, num_true_abs, num_false_abs,
                plottitle, num_true_perc, num_false_perc,
                self.plotfilelocation, plotfilename, ".pdf",
                ylab))
        fn = self.get_r_plottemptxt()
        with open(fn, "a") as fh:
            fh.write(string2r)
  
    def plot_riaincreasing_allpeps_frequency(self, plotfilelocation, plottitle, plotfilename, num_true, num_false, xlab1, xlab2, ylab):
        """
        """
        self.setplotfilelocation(plotfilelocation)
        string2r = ("""
main<-"%s"
num_true<-c(%s)
num_false <-c(%s)
height<-c(num_false, num_true)
ylab<-"%s"
height.names<-c("%s", "%s")
pdf("%s%s%s", paper="a4r")
barplot(height, names.arg=height.names, main=main, ylab=ylab)
dev.off()

""" % (plottitle,
               num_true,
               num_false,
               ylab,
               xlab1,
               xlab2,
               self.plotfilelocation, plotfilename, ".pdf"))
        fn = self.get_r_plottemptxt()
        with open(fn, "a") as fh:
            fh.write(string2r)
               
    def plot_hist_coveragetemplate(self, bins, plotfilename, plottitle, x_axis_label, y_axis_label_freq, y_axis_label_dens, fractionfound_list, plotfrequency=False, plotdensity=False):
        """
        """
        fractionfound_list = self.setvals(fractionfound_list)
        bins = self.setvals(bins)
        if plotfrequency:
            self.plot_hist_coveragetemplate_frequency(fractionfound_list, bins, x_axis_label, y_axis_label_freq, plottitle, plotfilename)
        if plotdensity:
            self.plot_hist_coveragetemplate_density(fractionfound_list, bins, x_axis_label, y_axis_label_dens, plottitle, plotfilename)
  
    def plot_hist_coveragetemplate_frequency(self, fractionfound_list, bins, x_axis_label, y_axis_label_freq, plottitle, plotfilename):
        string2r = ("""
fractionfound_list<-c(%s)
bins<-c(%s)
xlab<-"%s" 
ylab<-"%s"
plottitle<-"%s" 
h = hist(fractionfound_list, breaks=bins)
h$density = h$counts/sum(h$counts)
h$density = h$density*100
pdf("%s%s%s", paper="a4r")
plot(h, freq=T, main=plottitle, xlab=xlab, ylab=ylab, col="gray", labels=TRUE)
dev.off()

""" % (fractionfound_list,
                   bins,
                   x_axis_label, 
                   y_axis_label_freq,
                   plottitle,
                   self.plotfilelocation, plotfilename, ".pdf"))
        fn = self.get_r_plottemptxt()
        with open(fn, "a") as fh:
            fh.write(string2r)      
            
    def plot_hist_coveragetemplate_density(self, fractionfound_list, bins, x_axis_label, y_axis_label_dens, plottitle, plotfilename):
        """
        #pdf("%s%s%s", paper="a4r", family="Arial")
        """
        string2r = ("""
plottitle<-"%s"
xlab<-"%s" 
ylab<-"%s"
fractionfound_list<-c(%s)
bins<-c(%s)
h = hist(fractionfound_list, breaks=bins)
h$density = h$counts/sum(h$counts)
h$density = h$density*100
pdf("%s%s%s", paper="a4r")
plot(h, freq=F, main=plottitle, xlab=xlab, ylab=ylab, col="gray", labels=TRUE)
dev.off()

""" % (plottitle,
                   x_axis_label, 
                   y_axis_label_dens,
                   fractionfound_list,
                   bins,
                   self.plotfilelocation, plotfilename, ".pdf"))
        fn = self.get_r_plottemptxt()
        with open(fn, "a") as fh:
            fh.write(string2r)
  
    def plotcamelstacks(self, xvals, yvals, yvals_orig, plottitle, plotfilename, plotfilelocation,  xlimbeg, xlimend, yaxis_tps_h, yaxis_positions, rt_list):
        """
        plot points, then connect with segments
        if crude = True --> plot raw_data_cruderange
        """
        self.xvals = self.setvals(xvals)
        self.yvals = self.setvals(yvals)
        self.yvals_orig = self.setvals(yvals_orig)
        self.yaxis_tps_h = self.setvals(yaxis_tps_h)
        self.yaxis_positions = self.setvals(yaxis_positions)
        self.rt_labels = self.setvals(rt_list)
        string2r = ("""
xvals<-c(%s)
yvals<-c(%s)
yvals_orig<-c(%s)
ylim_max<-max(yvals)
ylim_dbl<-c(0, ylim_max)
xlim_dbl<-c(%s, %s)
xlab<-"m/z"
ylab<-"Timepoints [h] / Intensity"
plottitle<-"%s" 
pdf("%s%s%s", paper="a4r")
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 0, 0, 2))
yaxis_positions<-c(%s)
yaxis_tps_h<-c(%s)
rt_labels<-c(%s)
plot(xvals, yvals, xlab=xlab, ylab=ylab, main = plottitle, ylim=ylim_dbl, xlim=xlim_dbl, yaxt="n")
s<-seq(length(xvals))
segments(xvals[s], yvals_orig[s], xvals[s], yvals[s])
axis(2, at=yaxis_positions, labels=yaxis_tps_h, las=2)
axis(4, at=yaxis_positions, labels=rt_labels, las=0)
mtext("Rt [min]", side=4, line=3)
dev.off()

""" % (self.xvals, self.yvals, self.yvals_orig,
                   xlimbeg, xlimend,
                   plottitle,
                   plotfilelocation, plotfilename, ".pdf", 
                   self.yaxis_tps_h, self.yaxis_positions,
                   self.rt_labels))
        fn = self.get_r_plottemptxt()
        with open(fn, "a") as fh:
            fh.write(string2r)
            
    def plot_prot2rias(self, aaseq_list, tp_q_list, plottitle):
        """
        expects:
        tp_q_list = TimePoint in hours and RIA, concatenated list
        aaseq_list = AAsequence, concatenated list,
        plottitle = string
        plots a simple scatterplot, each tp-q pair of a series gets a different symbol,
        legend with sequences and their corresponding symbol
        """ 
        self.pepseqs=""
        self.pchs=""
        pch_dbl=[19, 17, 8, 6, 1, 2, 3, 4, 5, 7, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 17, 8, 6, 1, 2, 3, 4, 5, 7, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 17, 8, 6, 1, 2, 3, 4, 5, 7, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 17, 8, 6, 1, 2, 3, 4, 5, 7, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 17, 8, 6, 1, 2, 3, 4, 5, 7, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 17, 8, 6, 1, 2, 3, 4, 5, 7, 9, 10, 11, 12, 13, 14, 15, 16, 18]
        tp_main=[]
        q_main=[]
        for vals in tp_q_list[0]:
            tp_main.append(vals[0])
            q_main.append(vals[1])
        tp_main = self.setvals(tp_main)
        q_main = self.setvals(q_main)
        pepseq_main = aaseq_list[0]
        ylim_max = 1.0
        plottitle_main = plottitle
        plotfilelocation = self.getplotfilelocation()
        plotfilename = plottitle_main 
        pepseqs = aaseq_list.__repr__()[1:-1]
        pchs = pch_dbl[0:len(tp_q_list)].__repr__()[1:-1]
        points_r = ""
        for i in range(1, len(tp_q_list)):
            tp=[]
            q=[]
            for vals in tp_q_list[i]:
                tp.append(int(vals[0]))
                q.append(vals[1])
            pepseq = aaseq_list[i]
            tp = self.setvals(tp)
            tp = "c(%s)" % tp
            q = self.setvals(q)
            q = "c(%s)" % q
            pch = pch_dbl[i]
            points = "points(%s, %s, pch=%s)" % (tp, q, pch)
            if points_r == "":
                points_r += points
            else:
                points_r += "; " + points
        plotfilelocation_name = ("%s%s%s" % (plotfilelocation, plotfilename, ".pdf"))
        self.plotrias(tp_main, q_main, plottitle_main, plotfilelocation_name, ylim_max, points_r, pepseqs, pchs)
   
    def plotrias(self, tp_main, q_main, plottitle_main, plotfilelocation_name, ylim_max, points_r, pepseqs, pchs):
        """
        ###tiff(file="%s%s%s", width = 20, height = 20, units = "cm", res = 600) old
        ##tiff(file="%s%s", width = 20, height = 20, units = "cm", res = 600)  --> plotfilelocation_name, ".tiff",
        ##pdf("%s", paper="a4r") --> plotfilelocation_name,
        """
        string2r = ("""
tp_main<-c(%s)
q_main<-c(%s)
plottitle_main<-"%s"
pdf("%s", paper="a4r")
plot(tp_main, q_main, pch=19, xlab="TimePoints[h]", ylab="RIA", ylim=c(0, %s), main = plottitle_main)
%s 
legend(bty=0, box.lty=NULL, "topleft", c(%s), pch = c(%s))
dev.off()

""" % (tp_main, 
               q_main, 
               plottitle_main,
               plotfilelocation_name,
               ylim_max,
               points_r,
               pepseqs, pchs))
        fn = self.get_r_plottemptxt()
        with open(fn, "a") as fh:
            fh.write(string2r)
       
    def settype(self, t):
        self.type = t
        
    def setplotfilelocation(self, plotfilelocation):
        self.plotfilelocation = plotfilelocation
        
    def setplotfilename(self, plotfilename):
        self.plotfilename = plotfilename
        
    def setplottitle(self, plottitle):
        self.plottitle = plottitle
        
    def setxlab(self, xlab):
        self.xlab = xlab
        
    def setylab(self, ylab):
        self.ylab = ylab
        
    def pair2seplists(self, mz_int_list):
        """
        expects list containing list with m/z and intensity values as pairs
        e.g.[[1085.500928, 16383930.895546], [1086.001046, 19037900.789467], ['nd', 'nd'], ...]
        creates separate lists for m/z and intensity
        returns one list with first element m/z-list and second intensity-list
        """
        xvals=[]
        yvals=[]
        for ele in mz_int_list:
            x, y = ele[0], ele[1]
            xvals.append(x)
            yvals.append(y)
        return [xvals, yvals]
        
    def setvals(self, vals):
        """
        expects xvals as list,
        returns comma separated string, "nd" replaced with
        nothing just skipped
        e.g. input [1, 2, nd, 4, 5, 6, 7, 8, 9, 10]
        returns '1, 2, 4, 5, 6, 7, 8, 9, 10'
        """
        rvals = []
        for i in vals:
            if i == "nd":
                continue
            else:
                rvals.append(i)
        r0vals = ", ".join(str(x) for x in rvals)
        return r0vals  
        
    def getplotfilelocation(self):
        return self.plotfilelocation
        
    def getplotfilename(self):
        return self.plotfilename
        
    def getplottitle(self):
        return self.plottitle                           
