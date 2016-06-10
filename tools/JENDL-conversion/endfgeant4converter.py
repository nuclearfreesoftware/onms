#numpy
import numpy as np

#pyne
from pyne.xs.data_source import ENDFDataSource
from pyne.endf import Evaluation as eva 
from pyne import nucname
from pyne.bins import stair_step

#periodictable
from periodictable import *

#debug
from pprint import pprint

#file handling
import os

# other generated files
import newmf
import geantdata

# Converter Class for files
# implemented to read all the entries occuring in the JENDL AN 2005 evaluation
class ENDFGeant4Converter:
    def __init__(self,
                 basedir = "",
                 inputdir="input/",
                 outputdir="output/",
                 datasource="JENDL",
                 librarytype="NMSAlphaLE",
                 ending=".dat",
                 maincrosssection=201):
        self.basedir = basedir
        if(basedir == ""):
            self.outputdir = outputdir
            self.csdir = inputdir
        else:
            self.outputdir = os.path.join(basedir, outputdir)
            self.csdir = os.path.join(basedir, inputdir)
        self.datasource = datasource
        self.librarytype = librarytype
        self.ending = ending
        self.maincrosssection = maincrosssection
        self.mttogeant = {4: {'dirname': "AN", 'infoType': 0},
                          16: {'dirname': "A2N", 'infoType': 3},
                          22: {'dirname': "AAN", 'infoType': 5},
                          28: {'dirname': "APN", 'infoType': 9},
                          50: {'dirname': "AN_0", 'infoType': 0},
                          51: {'dirname': "AN_1", 'infoType': 0},
                          52: {'dirname': "AN_2", 'infoType': 0},
                          53: {'dirname': "AN_3", 'infoType': 0},
                          54: {'dirname': "AN_4", 'infoType': 0},
                          55: {'dirname': "AN_5", 'infoType': 0},
                          56: {'dirname': "AN_6", 'infoType': 0},
                          57: {'dirname': "AN_7", 'infoType': 0},
                          58: {'dirname': "AN_8", 'infoType': 0},
                          59: {'dirname': "AN_9", 'infoType': 0},
                          60: {'dirname': "AN_10", 'infoType': 0},
                          61: {'dirname': "AN_11", 'infoType': 0},
                          62: {'dirname': "AN_12", 'infoType': 0},
                          63: {'dirname': "AN_13", 'infoType': 0},
                          64: {'dirname': "AN_14", 'infoType': 0},
                          65: {'dirname': "AN_15", 'infoType': 0},
                          66: {'dirname': "AN_16", 'infoType': 0},
                          67: {'dirname': "AN_17", 'infoType': 0},
                          68: {'dirname': "AN_18", 'infoType': 0},
                          69: {'dirname': "AN_19", 'infoType': 0},
                          70: {'dirname': "AN_20", 'infoType': 0},
                          71: {'dirname': "AN_21", 'infoType': 0},
                          72: {'dirname': "AN_22", 'infoType': 0},
                          73: {'dirname': "AN_23", 'infoType': 0},
                          74: {'dirname': "AN_24", 'infoType': 0},
                          75: {'dirname': "AN_25", 'infoType': 0},
                          76: {'dirname': "AN_26", 'infoType': 0},
                          77: {'dirname': "AN_27", 'infoType': 0},
                          78: {'dirname': "AN_28", 'infoType': 0},
                          79: {'dirname': "AN_29", 'infoType': 0},
                          80: {'dirname': "AN_30", 'infoType': 0},
                          81: {'dirname': "AN_31", 'infoType': 0},
                          82: {'dirname': "AN_32", 'infoType': 0},
                          83: {'dirname': "AN_33", 'infoType': 0},
                          84: {'dirname': "AN_34", 'infoType': 0},
                          85: {'dirname': "AN_35", 'infoType': 0},
                          86: {'dirname': "AN_36", 'infoType': 0},
                          87: {'dirname': "AN_37", 'infoType': 0},
                          88: {'dirname': "AN_38", 'infoType': 0},
                          89: {'dirname': "AN_39", 'infoType': 0},
                          90: {'dirname': "AN_40", 'infoType': 0},
                          91: {'dirname': "AN_C", 'infoType': 0},
                          201: {'dirname': "ANprod", 'infoType': 201}
             }

    def convert(self):
        onlyfiles = [ f for f in os.listdir(self.csdir) if (os.path.isfile(os.path.join(self.csdir,f)) & f.endswith(self.ending)) ]
        for datafile in onlyfiles:
            print "====== Loading JENDL (alpha,n) File: " + datafile
            nuclidefilename = self.readendf(datafile)
            self.writegeant(nuclidefilename)
            print "====== Finished JENDL (alpha,n) File: " + datafile


    def readendf(self, datafile):
        endfds = ENDFDataSource(os.path.join(self.csdir, datafile))
        endfds.load()

        #    print endfds.library.intdict
        print "++++ The File contains data for the following nuclides:"
        self.lines = {}
        for nuc, value in endfds.library.mat_dict.iteritems() :
            cname = nucname.name(nuc)
            print "+++ Data for nuclide " + cname
            elname = elements[nucname.znum(nuc)].name
            elname = elname.capitalize()
            nuclidefilename = str(nucname.znum(nuc)) + "_" + str(nucname.anum(nuc)) + "_" + elname
            if (nucname.snum(nuc) != 0):
                print "Found nuclide in exited state, will add _<state> to filename"
                nuclidefilename = nuclidefilename + "_" + str(nucname.snum(nuc))
            for keyb, valueb in value['mfs'].iteritems() :
                print "+++ The nuclide has data for MF=", keyb[0], " and MT=", keyb[1]
                # MF=3 (cs)            
                if(keyb[0] == 3) :
                    print "Convert to geant compatible lines for MF=3"
                    reactiondata = endfds.reaction(nuc, keyb[1])
                    if (len(reactiondata['intschemes']) == 1):
                        if (reactiondata['intschemes'] == 2) :
                            writeline = ""
                            writeline += "  {:12d}".format(len(reactiondata['xs']))
                            writeline += "\n"
                            zippedxs = zip(reactiondata['xs'], reactiondata['e_int'])
                            count = 0
                            for sxs, se in zippedxs:
                                count += 1
                                writeline += "  {:8.6e}  {:8.6e}".format(se, sxs)
                                if (count == 3):
                                    writeline += "\n"
                                    count = 0
                            if(count != 0):
                                writeline += "\n"
                            if(not keyb[1] in self.lines):
                                self.lines.update({keyb[1]: {}})
                            self.lines[keyb[1]].update({keyb[0]: writeline})

                        else:
                            print "Non-linear interpolation scheme (Type {}) is currently not supported.".format(reactiondata['intschemes'])
                    else:
                        print "Multiple interpolation schemes are currently not supprted."
                        #print reactiondata
                        #exit()
                if(keyb[0] == 1) :
                    #                xsdata =  endfds.library.get_rx(key, 1, keyb[1]).reshape(-1,6)
                    #                print xsdata
                    print "Conversion for MF=1 not yet possible"
                if(keyb[0] == 6) :
                    print "Convert to geant compatible lines for MF=6"
                    mf6data = newmf.mf6()
                    mf6data.setdatasource(endfds)
                    mf6data.loadmtenangdata(nuc, keyb[1])
                    enangdist = mf6data.getdata()

                    gdata = geantdata.geantdata()
                    gdata.importenangdist(enangdist)
                    gdata.writesection(6, keyb[1])

                    if(not gdata.producterror):
                        if(not keyb[1] in self.lines):
                            self.lines.update({keyb[1]: {}})
                        self.lines[keyb[1]].update({keyb[0]: gdata.getline()})
                    else:
                        print("No output for MT/MF combination because of errors")
        return nuclidefilename

    def writegeant(self, nuclidefilename):
        for mt in self.lines:
            if(mt == self.maincrosssection):
                fullfilename = "Inelastic/CrossSection/" + nuclidefilename
                fullfilename = self.outputdir + fullfilename
                if not os.path.exists(os.path.dirname(fullfilename)):
                    os.makedirs(os.path.dirname(fullfilename))
                    print os.path.dirname(fullfilename) + " has been created (did not exist before)."
                outf = open(fullfilename, 'w')
                outf.write(self.librarytype + '\n')
                outf.write(self.datasource + '\n')
                outf.write(self.lines[mt][3])
                outf.close()
                print "Main cross section has been written to " + fullfilename
            if(mt in self.mttogeant):
                if(6 in self.lines[mt]):
                    fullfilename = "Inelastic/" + self.mttogeant[mt]['dirname'] + "/" + nuclidefilename
                    fullfilename = self.outputdir + fullfilename
                    if not os.path.exists(os.path.dirname(fullfilename)):
                        os.makedirs(os.path.dirname(fullfilename))
                        print os.path.dirname(fullfilename) + " has been created (did not exist before)."
                    outf = open(fullfilename, 'w')
                    outf.write("  {:12d}".format(self.mttogeant[mt]['infoType']))
                    outf.write("  {:12d}".format(3))
                    outf.write("\n")
                    outf.write("  {:12d}".format(mt)) 
                    outf.write("\n")
                    outf.write("  {:12d}".format(0))
                    outf.write("\n")
                    outf.write(self.lines[mt][3])
                    outf.write(self.lines[mt][6])
                    outf.close()
                    print "Reaction cross section and angular distribution data have been written to " + fullfilename
                else:
                    print "There is no angular distribution data for MT=%d, hence currently no output is generated." % mt
            else:
                print "MT=%d is not part of conversion dictionary" % mt

            
    def convertF0(self):
        onlyfiles = [ f for f in os.listdir(self.csdir) if (os.path.isfile(os.path.join(self.csdir,f)) & f.endswith(self.ending)) ]
        for datafile in onlyfiles:
            print "====== Loading JENDL (alpha,n) File: " + datafile
            endfeval = eva(os.path.join(self.csdir, datafile), verbose = False)
            endfeval.read()
            nuc = nucname.zzzaaa_to_id(endfeval.target['ZA'])
            elname = elements[nucname.znum(nuc)].name
            elname = elname.capitalize()
            nuclidefilename = str(nucname.znum(nuc)) + "_" + str(nucname.anum(nuc)) + "_" + elname
            if (endfeval.target['isomeric_state'] != 0):
                print "Found nuclide in exited state, will add _<state> to filename"
                nuclidefilename = nuclidefilename + "_" + str(endfeval.target['isomeric_state'])

            mf3data = {}
            fullfilename = "Inelastic/F01/" + nuclidefilename
            fullfilename = self.outputdir + fullfilename
            if not os.path.exists(os.path.dirname(fullfilename)):
                os.makedirs(os.path.dirname(fullfilename))
                print os.path.dirname(fullfilename) + " has been created (did not exist before)."
            outf = open(fullfilename, 'w')
            for mt in [4] + range(50, 92):
                if(mt in endfeval.reactions):
                    mf3data[mt] = endfeval.reactions[mt]
                    if(mf3data[mt].xs.nbt[0] != len(mf3data[mt].xs.x)):
                       print("Problem with interpolation")
                    #print mf3data[mt].xs.interp
                    outf.write("  {:12d}".format(1))
                    outf.write("\n")
                    outf.write("  {:12d}".format(3))
                    outf.write("\n")
                    outf.write("  {:12d}".format(mt)) 
                    outf.write("  {:12d}".format(0))
                    outf.write("\n")
                    outf.write("  {:12d}".format(int(mf3data[mt].Q_reaction)))
                    outf.write("  {:12d}".format(0))
                    outf.write("  {:12d}".format(len(mf3data[mt].xs.x)))
                    outf.write("\n")
                    writeline = ""
                    zippedxs = zip(mf3data[mt].xs.x, mf3data[mt].xs.y)
                    count = 0
                    for se, sxs in zippedxs:
                        count += 1
                        writeline += "  {:8.6e}  {:8.6e}".format(se, sxs)
                        if (count == 3):
                            writeline += "\n"
                            count = 0
                    if(count != 0):
                        writeline += "\n"
                    outf.write(writeline)

            outf.close()
