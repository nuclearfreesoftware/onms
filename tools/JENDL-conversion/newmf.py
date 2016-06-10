#!/usr/bin/python

import numpy as np

#pyne
from pyne.xs.data_source import ENDFDataSource

def _get_tab2(headkeys, lines):
    """Read some lines of the array, treating it as a TAB2 record.

    Parameters
    -----------
    headkeys: iterable, length 6
        An iterable containing the labels for each field in the first
        line. For empty/unassigned fields, use 0.
    lines: two-dimensional array-like
        The lines to be read. Each line should have 6 elements. The first
        line should be the first line of the TAB2 record; since we don't
        know the length of the TAB2 record, the last line should be the last
        line it is plausible for the TAB2 record to end.

    Returns
    --------
    head: dict
        Contains elements of the first card paired with their labels.
    intmeta: dict
        Contains the interpolation meta data.
    total_lines: int
        The number of lines the TAB2 record takes up.
    """
    head = dict(zip(headkeys, lines[0]))
    if 0 in head:
        del head[0]
    nr, nz = int(lines[0][4]), int(lines[0][5])
    
    intmeta_len = (nr*2 - 1)//6 + 1
    intmeta = dict(zip(('intpoints','intschemes'),
                       (lines[1:1+intmeta_len].flat[:nr*2:2],
                        lines[1:1+intmeta_len].flat[1:nr*2:2])))
    total_lines = 1 + intmeta_len
    return head, intmeta, total_lines


class mf6:
    verbose = 0
    totallines = 0
    
    def loadfile(self, filename):
        self.endfds = ENDFDataSource(filename)
        self.endfds.load()

    def setdatasource(self, datasource):
        self.endfds = datasource
        
    def loadmtenangdata(self, nuc, mt):
        self.totallines = 0
        self.enangdist = {}
        self.angrawdata = self.endfds.library.get_rx(nuc, 6, mt).reshape(-1,6)
        # catch exeptions

        head_flags = self.endfds.library._get_head(('ZA', 'AWR', 0, 'LCT', 'NK', 0),
                                                      self.angrawdata[self.totallines])
        self.totallines += 1
        head_flags['QM'] = self.getqm(nuc, mt)
        head_flags['QI'] = self.getqi(nuc, mt)
        self.enangdist['head_flags'] = head_flags

        self.enangdist['product_data'] = []
        self.enangdist['zap_index'] = {}

        # cycle through reaction products
        for subsec in range(0, int(head_flags['NK'])):
            self.readproduct(subsec)
    
    def getdata(self):
        return self.enangdist

    def getqm(self, nuc, mt):
        reactiondata = self.endfds.library.get_xs(nuc, mt, nuc)
        return reactiondata[1]['QM']

    def getqi(self, nuc, mt):
        reactiondata = self.endfds.library.get_xs(nuc, mt, nuc)
        return reactiondata[1]['QI']

    def readproduct(self, subsec):
        int_flags, int_data, int_size = self.endfds.library._get_tab1(
            ('ZAP', 'AWP', 'LIP', 'LAW', 'NR', 'NP'),
            ('x', 'y'),
            self.angrawdata[self.totallines:])
        self.totallines += int_size
        zap = int(int_flags['ZAP'])
        if zap in self.enangdist['zap_index']:
            print "There is already an entry for zap = " + str(zap)
            # handle this!
        else:
            singlezap = {}
            singlezap['tab1_flags'] = int_flags # 
            singlezap['tab1_data'] = int_data
            self.enangdist['product_data'].extend([singlezap])
            insertedidx =  self.enangdist['product_data'].index(singlezap)
            self.enangdist['zap_index'].update({int_flags['ZAP']: insertedidx})

        if int_flags['LAW'] == 0:  # Unknown Distribution (p. 125)
            self.readunknowndis(insertedidx)
        elif int_flags['LAW'] == 1:  # Continuum Energy-Angle Distributions (p. 126)
            self.readcontinuumenergyangledis(insertedidx)
        elif int_flags['LAW'] == 2:  # Discrete Two-Body Scattering (p. 131)
            self.readdiscretetwobodyscattering(insertedidx)
        elif int_flags['LAW'] == 3:  # Isotropic Discrete Emission (p. 132)
            print("LAW={0} is not implemented yet".format(int_flags['LAW']))
        elif int_flags['LAW'] == 4:  # Discrete Two-Body Recoils (p. 132)
            if(int_flags['NP']) > 3:
                print "found one"
                exit(1)
            # law=4 has no substructure.
            # there is nothing to do - option just kept for completeness
            True
        elif int_flags['LAW'] == 5:  # Charged-Particle Elastic Scattering (p. 132)
            print("LAW={0} is not implemented yet".format(int_flags['LAW']))
        elif int_flags['LAW'] == 6:  # N-Body Phase-Space Distributions (p. 136)
            print("LAW={0} is not implemented yet".format(int_flags['LAW']))
        elif int_flags['LAW'] == 7:  # Laboratory Angle-Energy Law (p. 137)
            print("LAW={0} is not implemented yet".format(int_flags['LAW']))
        else:
            print("Unknown Product Energy-Angle Distribution")
    #END readoutproduct

    def readunknowndis(self, idx):
        errormessage = "LAW={0} is not implemented yet".format(self.enangdist['product_data'][idx]['tab1_flags']['LAW'])
        print(errormessage)
        self.enangdist['product_data'][idx] = {}
        self.enangdist['product_data'][idx]['error'] = errormessage

    def readcontinuumenergyangledis(self, idx):
        # read [MAT, 6, MT/ 0.0, 0.0, LANG, LEP, NR, NE/ E int]TAB2
        tab2_head, tab2_meta, tab2_size = _get_tab2(
            (0, 0, 'LANG', 'LEP', 'NR', 'NE'),
            self.angrawdata[self.totallines:])
        self.enangdist['product_data'][idx]['tab2_flags'] = tab2_head
        self.enangdist['product_data'][idx]['tab2_intdata'] = tab2_meta
        self.totallines += tab2_size

        self.enangdist['product_data'][idx]['list_data'] = {}
        for energyidx in range(int(tab2_head['NE'])):
            # read the following structure
            # [MAT, 6, MT/ 0.0,  E1,  ND,  NA,  NW, NEP/
            # E1', b0(E1, E1'), b1(E1, E1'), -------- bNA(E1, E1'),
            # E2', b0(E1, E2'), b1(E1, E2'), -------- bNA(E1, E2'),
            # --------------------------------------------
                            # Enep', b0(E1 , Enep'), b1(E1, E'nep), ---- bNA(E1, Enep')]LIST
            list_head, list_data, list_size = self.endfds.library._get_list(
                (0, 'E1', 'ND', 'NA', 'NW', 'NEP'), ('b'), self.angrawdata[self.totallines:])
            self.totallines += list_size

            # required to flatten data
            list_data = list_data['b']
            
            eprimes = []
            coefficientlists = []
            for energypidx in range(int(list_head['NEP'])):
                eidx = energypidx * (int(list_head['NA']) + 2)
                start = eidx + 1
                stop = start + int(list_head['NA']) + 1
                eprimes.extend([list_data[int(eidx)]])
                coefficientlists.extend([list_data[start:stop]])

            listentry = list_head
            listentry.update({'eprimes': eprimes})
            listentry.update({'coefficientlists': coefficientlists})

            self.enangdist['product_data'][idx]['list_data'].update({list_head['E1']: listentry})
    #END readoutcontinuusenergyangledis
            
    def readdiscretetwobodyscattering(self, idx):
        # read [MAT, 6, MT/ 0.0, 0.0,0, 0, NR, NE/ E int ]TAB2
        tab2_head, tab2_meta, tab2_size = _get_tab2(
            (0, 0, 0, 0, 'NR', 'NE'),
            self.angrawdata[self.totallines:])
        self.enangdist['product_data'][idx]['tab2_flags'] = tab2_head
        self.enangdist['product_data'][idx]['tab2_intdata'] = tab2_meta
        self.totallines += tab2_size
        
        self.enangdist['product_data'][idx]['list_data'] = {}
        for energyidx in range(int(tab2_head['NE'])):
            # read the following structure
            # [MAT, 6, MT/ 0.0, E1, LANG, 0, NW, NL/ A l (E)]LIST
            list_head, list_data, list_size = self.endfds.library._get_list(
                (0, 'E1', 'LANG', 0, 'NW', 'NL'), ('b'), self.angrawdata[self.totallines:])
            self.totallines += list_size

            # required to flatten data
            list_data = list_data['b']

            listentry = list_head

            al = []
            for listidx in range(int(list_head['NL'])):
                if(list_head['LANG'] == 0): # Legendre coefficients
                    al.extend([list_data[listidx]])
                else: # tabulated cosine
                    start = listidx * 2
                    stop = (start + 1) + 1 # upper bound is non-inclusive!
                    al.extend([list_data[start:stop]])
            listentry.update({'al': al})

            self.enangdist['product_data'][idx]['list_data'].update({list_head['E1']: listentry})

        #END readdiscretetwobodyscattering

    def getisotropicdistribution(self, Ei, prodzaid = 0, grid = 'auto'):
        #isotropic in this case: integrated over all possible mu (cosine theta)
        prod = -1
        for zaid in self.enangdist['zap_index']:
            if(prodzaid == int(zaid)):
                prod = self.enangdist['zap_index'][zaid]
        if(prod < 0):
            print("Error: Unknown product index!")
            return {}
        law = self.enangdist['product_data'][prod]['tab1_flags']['LAW']
        # print(law)
        if(law == 0):
            # Unknown Distribution - empty list.
            return {}
        elif(law == 1):
            # 
            lang = int(self.enangdist['product_data'][prod]['tab2_flags']['LANG'])
            if(lang == 1):
                return {}
            elif(lang == 2):
                fulldata = self.interpolateTab2(Ei,
                                           self.enangdist['product_data'][prod]['tab2_intdata'],
                                           self.getElist(prodzaid),
                                           self.enangdist['product_data'][prod]['list_data'],
                                           y = 'eprimes',
                                           z = 'coefficientlists')
                # only first parameter for isotopric / integral
                fdata = [row[0] for row in fulldata['coefficientlists']]
                if(grid != 'auto' and grid != []):
                    fdatai = np.zeros(len(grid))
                    for eidx in range(len(grid)):
                        fdatai[eidx] = self.general1Dinterpol(grid[eidx], fulldata['eprimes'], fdata, 2)
                    return {'e': np.array(grid), 'f': np.array(fdatai)}
                else:
                    return {'e': np.array(fulldata['eprimes']), 'f': np.array(fdata)}
            elif(lang >= 11 and lang <= 15):
                return {}
        elif(law == 2):
            # Discrete Two Body Scattering (has a list)
            return sorted(self.enangdist['product_data'][prod]['list_data'].keys())

        elif(law == 3):
            # Isotropic Discrete Emission - empty list
            return {}
        elif(law == 4):
            # Discrete Two-Body Recoils - empty list
            return {}
        elif(law == 5):
            # Charged-Particle Elastic Scattering
            return sorted(self.enangdist['product_data'][prod]['list_data'].keys())
        elif(law == 6):
            # N-body phase space
            return {}
        elif(law == 7):
            return sorted(self.enangdist['product_data'][prod]['list_data'].keys())

    def getElist(self, prodzaid = 0):
        prod = -1
        for zaid in self.enangdist['zap_index']:
            if(prodzaid == int(zaid)):
                prod = self.enangdist['zap_index'][zaid]
        if(prod < 0):
            print("Error: Unknown product index!")
            return {}
        law = int(self.enangdist['product_data'][prod]['tab1_flags']['LAW'])
        if(law == 0):
            # Unknown Distribution - empty list.
            return []
        elif(law == 1):
            # 
            return sorted(self.enangdist['product_data'][prod]['list_data'].keys())
        elif(law == 2):
            # Discrete Two Body Scattering (has a list)
            return sorted(self.enangdist['product_data'][prod]['list_data'].keys())
        elif(law == 3):
            # Isotropic Discrete Emission - empty list
            return []
        elif(law == 4):
            # Discrete Two-Body Recoils - empty list
            return []
        elif(law == 5):
            # Charged-Particle Elastic Scattering
            return sorted(self.enangdist['product_data'][prod]['list_data'].keys())
        elif(law == 6):
            # N-body phase space
            return []
        elif(law == 7):
            return sorted(self.enangdist['product_data'][prod]['list_data'].keys())

    def interpolateTab1(self, E, tab1data, x = 'x', y = 'y'):
        # get interpolation scheme
        lowint = 0
        for i in range(len(tab1data['intpoints'])):
            highint = tab1data['intpoints'][i] - 1
            if(E <= tab1data[x][highint]):
                break
            lowint = highint
        intscheme = tab1data['intschemes']
        # get low / high E
        index = lowint
        if(E < tab1data[x][index]):
            return tab1data[y][0]
        while(index < highint and E >= tab1data[x][index]):
            if(E == tab1data[x][index]):
                return tab1data[y][index]
            index += 1
        if(index == highint): # E is bigger than last point
            return tab1data[y][-1]
        else:
            return self.simple1Dinterpol(E, tab1data[x][index - 1], tab1data[x][index], tab1data[y][index - 1], tab1data[y][index], intscheme)

    def general1Dinterpol(self, E, xlist, ylist, scheme = 2):
        index = 0
        if(E < xlist[index]):
            return ylist[0]
        while(index < len(ylist) and E >= xlist[index]):
            if(E == xlist[index]):
                return ylist[index]
            index += 1
        if(E > xlist[-1]): # E is bigger than last point
            return ylist[-1]
        else:
            return self.simple1Dinterpol(E, xlist[index - 1], xlist[index],
                                    ylist[index - 1], ylist[index],
                                    scheme)
    
    def simple1Dinterpol(self, E, x1, x2, y1, y2, scheme = 2):
        if(scheme == 1):
            return y1
        elif(scheme == 2):
            return y1 + (E - x1) * ( y2 - y1 ) / (x2 - x1)
        elif(scheme == 3):
            return 0
        elif(scheme == 4):
            return 0
        elif(scheme == 5):
            return 0
        elif(scheme == 6):
            return 0

    def interpolateTab2(self, E, intdata, xdata, listdata, y = 'eprimes', z = 'coefficientlists'):
        #outside area
        if(E < xdata[0]):
            tmp = [0] * len(listdata[xdata[0]][z])
            return {y: [0], z: [tmp]}
        # get interpolation scheme
        lowint = 0
        for i in range(len(intdata['intpoints'])):
            highint = int(intdata['intpoints'][i]) - 1
            if(E <= xdata[highint]):
                break
            lowint = highint
        intscheme = intdata['intschemes']
        index = lowint
        while(index < highint and E >= xdata[index]):
            if(E == xdata[index]):
                return {y: listdata[xdata[index]][y], z: listdata[xdata[index]][z]}
            index += 1
        if(index == highint): # E bigger than last point
            return {y: listdata[xdata[-1]][y], z: listdata[xdata[-1]][z]}
        else:
            tmp = self.simple2Dinterpol(E, xdata[index - 1], xdata[index],
                                        listdata[xdata[index - 1]][y], listdata[xdata[index]][y],
                                        listdata[xdata[index - 1]][z], listdata[xdata[index]][z],
                                        intscheme)
            return {y: tmp['y'], z: tmp['z']}

    def simple2Dinterpol(self, E, x1, x2, y1, y2, z1, z2, scheme = 22):
        if(scheme < 10):
            # cartesian
            return 0
        elif(scheme < 20):
            # mce
            return 0
        else:
            # unit based
            ylowmin = y1[0]
            ylowmax = y1[-1]
            yhighmin = y2[0]
            yhighmax = y2[-1]
            yplow = ylowmin + (E - x1) / (x2 - x1) * (yhighmin - ylowmin)
            yphigh = ylowmax + (E - x1) / (x2 - x1) * (yhighmax - ylowmax)
            # print(ylowmin, ylowmax, yhighmin, yhighmax, yplow, yphigh)

            if(scheme == 21):
                return {'y': y1, 'z': z1}
            elif(scheme == 22):
                ylist = [ (( y - ylowmin)/(ylowmax - ylowmin)) for y in y1 ] + [ (( y - yhighmin)/(yhighmax - yhighmin)) for y in y2 ]
                ylist = sorted(set(ylist))
                zlist = [0] * len(ylist)
                for yidx, yT in enumerate(ylist):
                    yTlow = ylowmin + (ylowmax - ylowmin) * yT
                    yThigh = yhighmin + (yhighmax - yhighmin) * yT
                    parl = [0] * len(z1[0]) # what if z1 / z2 have different parameter length
                    for par in range(len(z1[0])):
                        parlist = [row[par] for row in z1]
                        plow = self.general1Dinterpol(yTlow, y1, parlist, 2)
                        parlist = [row[par] for row in z2]
                        phigh = self.general1Dinterpol(yThigh, y2, parlist, 2)
                        plow = plow * (ylowmax - ylowmin)
                        phigh = phigh * (yhighmax - yhighmin)
                        # ptot = self.simple1Dinterpol((yplow + yT * (yphigh - yplow)),
                        #                              x1, x2, plow, phigh)
                        ptot = self.simple1Dinterpol(E,
                                                     x1, x2, plow, phigh)
                        #print(yT, (yplow + yT * (yphigh - yplow)), plow, phigh, ptot)
                        ptot = ptot / (yphigh - yplow)
                        parl[par] = ptot
                    zlist[yidx] = parl
                ylist = [ (yplow + y * (yphigh - yplow)) for y in ylist ]
                return {'y': ylist, 'z': zlist}
            elif(scheme == 23):
                return {'y': 0, 'z': 0}
            elif(scheme == 24):
                return {'y': 0, 'z': 0}
            elif(scheme == 25):
                return {'y': 0, 'z': 0}
