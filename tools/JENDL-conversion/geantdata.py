#!/usr/bin/python

def outputtab(length, x, y, xformat=" {: 8.6e}", yformat=" {: 8.6e}"):
    outputline = ""
    for i in range(length):
        outputline += xformat.format(x[i])
        outputline += yformat.format(y[i])
        lastnewline = False
        if(i % 3 == 2):
            lastnewline = True
            outputline += "\n"
    if(not lastnewline):
        outputline += "\n"
    return outputline

class geantdata:
    def importenangdist(self, enangdist):
        self.enangdist = enangdist

    def getline(self):
        return self.writeline
        
    def writesection(self, mf, mt):

        self.writeline = ""
        self.writeline += "  {:12d}\n".format(mt)
        self.writeline += "  {:12d}\n".format(mf)

        if(mf == 6):
            self.writeenangdist()

    def writeenangdist(self):
        #  	<AWR> 	<LCT> 	<NK>
        self.writeline += " {: 8.6e}".format(self.enangdist['head_flags']['AWR'])
        self.writeline += "  {:12d}".format(int(self.enangdist['head_flags']['LCT']))
        self.writeline += "  {:12d}\n".format(int(self.enangdist['head_flags']['NK']))

        self.producterror = False
        for i in range(int(self.enangdist['head_flags']['NK'])):
            if('error' in self.enangdist['product_data'][i]):
                print "error for single product, will not output anything"
                self.producterror = True
            else:
                #product head + tab1
                #  <ZAP>	<AWP> 	<LIP> 	<LAW>
                self.writeline += " {: 8.6e}".format(self.enangdist['product_data'][i]['tab1_flags']['ZAP'])
                self.writeline += " {: 8.6e}".format(self.enangdist['product_data'][i]['tab1_flags']['AWP'])
                self.writeline += "  {:12d}".format(int(self.enangdist['product_data'][i]['tab1_flags']['LIP']))
                self.writeline += "  {:12d}".format(int(self.enangdist['product_data'][i]['tab1_flags']['LAW']))
                self.writeline += "\n"
                #  <QM> <QI> <NP>
                self.writeline += " {: 8.6e}".format(self.enangdist['head_flags']['QM'])
                self.writeline += " {: 8.6e}".format(self.enangdist['head_flags']['QI'])
                self.writeline += "  {:12d}".format(int(self.enangdist['product_data'][i]['tab1_flags']['NP']))
                self.writeline += "\n"
                # <NR>
                self.writeline += "  {:12d}".format(int(self.enangdist['product_data'][i]['tab1_flags']['NR']))
                self.writeline += "\n"
                #<NBT(1)> 	<INT(1)> ... <NBT(NR)> <INT(NR)>
                self.writeline += outputtab(
                    int(self.enangdist['product_data'][i]['tab1_flags']['NR']),
                    map(int, self.enangdist['product_data'][i]['tab1_data']['intpoints']),
                    map(int, self.enangdist['product_data'][i]['tab1_data']['intschemes']),
                    "  {:12d}", "  {:12d}")

                # <x(1)>       <y(1)> ... <x(NP)> <y(NP)>
                self.writeline += outputtab(
                    int(self.enangdist['product_data'][i]['tab1_flags']['NP']),
                    self.enangdist['product_data'][i]['tab1_data']['x'],
                    self.enangdist['product_data'][i]['tab1_data']['y']
                )

                #tab2 + list depending on law
                if(int(self.enangdist['product_data'][i]['tab1_flags']['LAW']) == 0):
                    print("not implemented")
                elif(int(self.enangdist['product_data'][i]['tab1_flags']['LAW']) == 1):
                    self.writecontinuumenergyangledist(i)
                elif(int(self.enangdist['product_data'][i]['tab1_flags']['LAW']) == 2):
                    self.writediscretetwobodyscattering(i)
                elif(int(self.enangdist['product_data'][i]['tab1_flags']['LAW']) == 3):
                    print("not implemented")
                elif(int(self.enangdist['product_data'][i]['tab1_flags']['LAW']) == 4):
                    print("not implemented")
                elif(int(self.enangdist['product_data'][i]['tab1_flags']['LAW']) == 5):
                    print("not implemented")
                elif(int(self.enangdist['product_data'][i]['tab1_flags']['LAW']) == 6):
                    print("not implemented")
                elif(int(self.enangdist['product_data'][i]['tab1_flags']['LAW']) == 7):
                    print("not implemented")
                else:
                    print("Unknown Product Energy-Angle Distribution")

    def writecontinuumenergyangledist(self, idx):
        # <ZA> 	<LANG> 	<LEP> 	<NE>
        self.writeline += "  {:12d}".format(int(self.enangdist['head_flags']['ZA']))
        self.writeline += "  {:12d}".format(int(self.enangdist['product_data'][idx]['tab2_flags']['LANG']))
        self.writeline += "  {:12d}".format(int(self.enangdist['product_data'][idx]['tab2_flags']['LEP']))
        self.writeline += "  {:12d}".format(int(self.enangdist['product_data'][idx]['tab2_flags']['NE']))
        self.writeline += "\n"

        # <NR> (tab2)
        self.writeline += "  {:12d}".format(int(self.enangdist['product_data'][idx]['tab2_flags']['NR']))
        self.writeline += "\n"

        #  	<NBT(1)> 	<INT(1)> ... <NBT(NR)> <INT(NR)>
        # incident energy interpolation ranges
        self.writeline += outputtab(
            int(self.enangdist['product_data'][idx]['tab2_flags']['NR']),
            map(int, self.enangdist['product_data'][idx]['tab2_intdata']['intpoints']),
            map(int, self.enangdist['product_data'][idx]['tab2_intdata']['intschemes']),
            "  {:12d}", "  {:12d}")

        for energy in sorted(self.enangdist['product_data'][idx]['list_data']):
            # <E_1> <NEP_1> <NE_1> <NA_1>
            self.writeline += " {: 8.6e}".format(energy)
            self.writeline += "  {:12d}".format(int(self.enangdist['product_data'][idx]['list_data'][energy]['NEP']))
            self.writeline += "  {:12d}".format(int(self.enangdist['product_data'][idx]['list_data'][energy]['ND']))
            self.writeline += "  {:12d}".format(int(self.enangdist['product_data'][idx]['list_data'][energy]['NA'])+1)
            self.writeline += "\n"
            lastnewline = False
            for j in range(int(self.enangdist['product_data'][idx]['list_data'][energy]['NEP'])):
                self.writeline += " {: 8.6e}".format(self.enangdist['product_data'][idx]['list_data'][energy]['eprimes'][j])

                entrycount = 1
                for k in range(int(self.enangdist['product_data'][idx]['list_data'][energy]['NA'])+1):
                    self.writeline += " {: 8.6e}".format(self.enangdist['product_data'][idx]['list_data'][energy]['coefficientlists'][j][k])
                    entrycount += 1
                    lastnewline = False
                    if(entrycount == 6):
                        self.writeline += "\n"
                        entrycount = 0
                        lastnewline = True
                if(not lastnewline):
                    self.writeline += "\n"
    #END writecontinuumenergyangledist

    def writediscretetwobodyscattering(self, idx):
        # <NE>
        self.writeline += "  {:12d}".format(int(self.enangdist['product_data'][idx]['tab2_flags']['NE']))
        self.writeline += "\n"

        # <NR> (tab2)
        self.writeline += "  {:12d}".format(int(self.enangdist['product_data'][idx]['tab2_flags']['NR']))
        self.writeline += "\n"

        #  	<NBT(1)> 	<INT(1)> ... <NBT(NR)> <INT(NR)>
        # incident energy interpolation ranges
        self.writeline += outputtab(
            int(self.enangdist['product_data'][idx]['tab2_flags']['NR']),
            map(int, self.enangdist['product_data'][idx]['tab2_intdata']['intpoints']),
            map(int, self.enangdist['product_data'][idx]['tab2_intdata']['intschemes']),
            "  {:12d}", "  {:12d}")


        for energy in sorted(self.enangdist['product_data'][idx]['list_data']):
            #   <E_1> <LANG> <NL>
            self.writeline += " {: 8.6e}".format(energy)
            self.writeline += "  {:12d}".format(int(self.enangdist['product_data'][idx]['list_data'][energy]['LANG']))
            self.writeline += "  {:12d}".format(int(self.enangdist['product_data'][idx]['list_data'][energy]['NL']))
            self.writeline += "\n"

            lastnewline = False
            entrycount = 0
            if(int(self.enangdist['product_data'][idx]['list_data'][energy]['LANG']) == 0):
                #legendre
                for j in range(int(self.enangdist['product_data'][idx]['list_data'][energy]['NL'])):
                    self.writeline += " {: 8.6e}".format(self.enangdist['product_data'][idx]['list_data'][energy]['al'][j])
                    entrycount += 1
                    lastnewline = False
                    if(entrycount == 6):
                        self.writeline += "\n"
                        entrycount = 0
                        lastnewline = True
            else:
                for j in range(int(self.enangdist['product_data'][idx]['list_data'][energy]['NL'])):
                    self.writeline += " {: 8.6e}".format(self.enangdist['product_data'][idx]['list_data'][energy]['al'][j][0])
                    self.writeline += " {: 8.6e}".format(self.enangdist['product_data'][idx]['list_data'][energy]['al'][j][1])
                    entrycount += 2
                    lastnewline = False
                    if(entrycount == 6):
                        self.writeline += "\n"
                        entrycount = 0
                        lastnewline = True
                #two coefficients
            if(not lastnewline):
                self.writeline += "\n"
    
    #END writediscretetwobodyscattering
