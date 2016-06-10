/*
 * ONMS - Open Neutron Multiplicity Simulation
 *
 *
 * Copyright (C) 2013-2016 Moritz Kütt
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact: moritz@nuclearfreesoftware.org
 */

#include "NMSMultiplicityResult.hh"
#include <iomanip>

NMSMultiplicityResult::NMSMultiplicityResult(int multiplicitiesInit, int varsplits) {
  multiplicities = multiplicitiesInit;
  splits = varsplits;
  settings.predelay = 4.5;
  settings.adelay = 4096;
  settings.registerLength = 128;
  settings.registerPeriod = 0.5;
  settings.dieaway = 50;
  settings.efficiency = 0.56;

  settings.measurementLength = 100000000;

  settings.doublegatefraction = simplegatefraction();
  settings.triplegatefraction = simplegatefraction() * simplegatefraction();
  
  lastEvent = 0;
  resultsRAptr = new int[multiplicities];
  resultsAptr = new int[multiplicities];
  normresultsRAptr = new double[multiplicities];
  normresultsAptr = new double[multiplicities];
  for(int i = 0; i < multiplicities; i++) {
    resultsRAptr[i] = 0;
    resultsAptr[i] = 0;
    normresultsRAptr[i] = 0;
    normresultsAptr[i] = 0;
  }
  splitresultsRAvec.clear();
  splitresultsAvec.clear();
  splitnormresultsRAvec.clear();
  splitnormresultsAvec.clear();
  if(splits != -1) {
    for(int i = 0; i < splits; i++) {
      std::vector<int> sravec;
      std::vector<int> savec;
      std::vector<double> snravec;
      std::vector<double> snavec;
      for(int j = 0; j < multiplicities; j++) {
	sravec.push_back(0);
	savec.push_back(0);
	snravec.push_back(0);
	snavec.push_back(0);
      }
      splitresultsRAvec.push_back(sravec);
      splitresultsAvec.push_back(savec);
      splitnormresultsRAvec.push_back(snravec);
      splitnormresultsAvec.push_back(snavec);
    }
  }
    

  // Data from SmpSpNuDistData from LLNL Fission Library
  // Matches also data given in P. Santi and M. Miller,
  // NUCLEAR SCIENCE AND ENGINEERING: 160, 190–199 (2008)
  // (except for last value, which is 0.0020406 in latter ref.)
  pu240sf.clear();
  pu240sf.push_back(0.0631852);
  pu240sf.push_back(0.2319644);
  pu240sf.push_back(0.3333230);
  pu240sf.push_back(0.2528207);
  pu240sf.push_back(0.0986461);
  pu240sf.push_back(0.0180199);
  pu240sf.push_back(0.0020407);
  // Zucker Holdren (2 MeV) as tabulated in Verbeke, J. M.; Hagmann, C. & Wright, D. Simulation of Neutron and Gamma Ray Emission from Fission and Photofission Lawrence Livermore National Laboratory, Lawrence Livermore National Laboratory, 2010 (UCRL-AR-228518)
  pu239fif.clear();
  pu239fif.push_back(0.0062555);
  pu239fif.push_back(0.0611921);
  pu239fif.push_back(0.2265608);
  pu239fif.push_back(0.3260637);
  pu239fif.push_back(0.2588354);
  pu239fif.push_back(0.0956070);
  pu239fif.push_back(0.0224705);
  pu239fif.push_back(0.0025946);
  pu239fif.push_back(0.0005205);
}

NMSMultiplicityResult::~NMSMultiplicityResult() {
  delete[] resultsRAptr;
  delete[] resultsAptr;
  delete[] normresultsRAptr;
  delete[] normresultsAptr;
}

NMSMultiplicityResult::NMSMultiplicityResult(const NMSMultiplicityResult& x) {
  settings = x.settings;
  multiplicities = x.multiplicities;
  splits = x.splits;
  lastEvent = x.lastEvent;
  resultsRAptr = new int[multiplicities];
  resultsAptr = new int[multiplicities];
  normresultsRAptr = new double[multiplicities];
  normresultsAptr = new double[multiplicities];
  for(int i = 0; i < multiplicities; i++) {
    resultsRAptr[i] = x.resultsRAptr[i];
    resultsAptr[i] = x.resultsAptr[i];
    normresultsRAptr[i] = x.normresultsRAptr[i];
    normresultsAptr[i] = x.normresultsAptr[i];
  }
  splitresultsRAvec = x.splitresultsRAvec;
  splitresultsAvec = x.splitresultsAvec;
  splitnormresultsRAvec = x.splitnormresultsRAvec;
  splitnormresultsAvec = x.splitnormresultsAvec;

}

NMSMultiplicityResult& NMSMultiplicityResult::operator = (const NMSMultiplicityResult& x) {
  delete[] resultsRAptr;
  delete[] resultsAptr;
  delete[] normresultsRAptr;
  delete[] normresultsAptr;

  settings = x.settings;

  multiplicities = x.multiplicities;
  splits = x.splits;
  lastEvent = x.lastEvent;
  resultsRAptr = new int[multiplicities];
  resultsAptr = new int[multiplicities];
  normresultsRAptr = new double[multiplicities];
  normresultsAptr = new double[multiplicities];
  for(int i = 0; i < multiplicities; i++) {
    resultsRAptr[i] = x.resultsRAptr[i];
    resultsAptr[i] = x.resultsAptr[i];
    normresultsRAptr[i] = x.normresultsRAptr[i];
    normresultsAptr[i] = x.normresultsAptr[i];
  }

  splitresultsRAvec.clear();
  splitresultsAvec.clear();
  splitnormresultsRAvec.clear();
  splitnormresultsAvec.clear();

  splitresultsRAvec = x.splitresultsRAvec;
  splitresultsAvec = x.splitresultsAvec;
  splitnormresultsRAvec = x.splitnormresultsRAvec;
  splitnormresultsAvec = x.splitnormresultsAvec;

  return *this;

}

void NMSMultiplicityResult::SetSettings(NMSMultiplicityMeasurementSetting newset) {
  settings = newset;
}



void NMSMultiplicityResult::setLastEvent(double le) {
  lastEvent = le;
}

double NMSMultiplicityResult::getLastEvent() {
  return lastEvent;
}

int NMSMultiplicityResult::getRA(int i, int split) {
  if( (i >= 0) && (i < multiplicities) ) {
    if(split == -1) {
      return resultsRAptr[i];
    }
    else {
      if(split >= 0 && split < splits) {
	return splitresultsRAvec[split][i];
      }
      else {
	return -1;
      }
    }
  }
  else {
    return -1;
  }
}

int NMSMultiplicityResult::getA(int i, int split) {
  if( (i >= 0) && (i < multiplicities) ) {
    if(split == -1) {
      return resultsAptr[i];
    }
    else {
      if(split >= 0 && split < splits) {
	return splitresultsAvec[split][i];
      }
      else {
	return -1;
      }
    }
  }
  else {
    return -1;
  }
}

void NMSMultiplicityResult::addRA(int i, int increment, int split) {
  if( (i >= 0) && (i < multiplicities) ) {
    if(split == -1) {
      resultsRAptr[i] += increment;
    }
    else {
      if(split >= 0 && split < splits) {
	splitresultsRAvec[split][i] += increment;
      }
    }
  }
}

void NMSMultiplicityResult::addA(int i, int increment, int split) {
  if( (i >= 0) && (i < multiplicities) ) {
    if(split == -1) {
      resultsAptr[i] += increment;
    }
    else {
      if(split >= 0 && split < splits) {
	splitresultsAvec[split][i] += increment;
      }
    }
  }
}

void NMSMultiplicityResult::set(int i, int ra, int a) {
  if( (i >= 0) && (i < multiplicities) )
    {
      if ( ( a >= 0 ) && ( ra >= 0 ) ) {
	resultsRAptr[i] = ra;
	resultsAptr[i] = a;
      }
    }
}

void NMSMultiplicityResult::setRA(int i, int ra) {
  if( (i >= 0) && (i < multiplicities) )
    {
      if ( ra >= 0 ) {
	resultsRAptr[i] = ra;
      }
    }
}

void NMSMultiplicityResult::setA(int i, int a) {
  if( (i >= 0) && (i < multiplicities) )
    {
      if ( a >= 0 ) {
	resultsAptr[i] = a;
      }
    }
}

void NMSMultiplicityResult::normalize() {
  for(int i=0; i < multiplicities; i++) {
    if(getRAsum() != 0) {
      normresultsRAptr[i] = 1.0 * resultsRAptr[i] / getRAsum();
    }
    else {
      normresultsRAptr[i] = 0;
    }
    if(getAsum() != 0) {
      normresultsAptr[i] = 1.0 * resultsAptr[i] / getAsum();
    }
    else {
      normresultsAptr[i] = 0;
    }
  }
  if(splits != -1) {
    for(int split=0; split<splits; split++) {
      for(int i=0; i < multiplicities; i++) {
	if(getRAsum(split) != 0) {
	  splitnormresultsRAvec[split][i] = 1.0 * splitresultsRAvec[split][i] / getRAsum(split);
	}
	else {
	  splitnormresultsRAvec[split][i] = 0;
	}
	if(getAsum(split) != 0) {
	  splitnormresultsAvec[split][i] = 1.0 * splitresultsAvec[split][i] / getAsum(split);
	}
	else {
	  splitnormresultsAvec[split][i] = 0;
	}
      }
    }
  }

}

double NMSMultiplicityResult::getNormRA(int i, int split) {
  if( (i >= 0) && (i < multiplicities) ) {
    if(split == -1) {
      return normresultsRAptr[i];
    }
    else {
      if(split >= 0 and split < splits) {
	return splitnormresultsRAvec[split][i];
      }
      else {
	return -1;
      }
    }
  }
  else {
    return -1;
  }
}

double NMSMultiplicityResult::getNormA(int i, int split) {
  if( (i >= 0) && (i < multiplicities) ) {
    if(split == -1) {
      return normresultsAptr[i];
    }
    else {
      if(split >= 0 and split < splits) {
	return splitnormresultsAvec[split][i];
      }
      else {
	return -1;
      }
    }
  }
  else {
    return -1;
  }
}

int NMSMultiplicityResult::getRAsum(int split) {
  int total = 0;
  for(int i = 0; i < multiplicities; i++) {
    total += getRA(i, split);
  }
  return total;
}

int NMSMultiplicityResult::getAsum(int split) {
  int total = 0;
  for(int i = 0; i < multiplicities; i++) {
    total += getA(i, split);
  }
  return total;
}


double NMSMultiplicityResult::getNuRA(int k, int split) {
  double a = 0;
  double factor = 1;
  if(k >= 0 && k < multiplicities) {
    for(int i = k; i < multiplicities; i++) {
      factor = 1;
      if(i > 0) {
	for(int j = 0; j < k; j++) {
	  factor *= ( i - j);
	}
      }
      a += factor * getNormRA(i, split);
    }
  }
  else {
    a = 0;
  }
  return a;
}

double NMSMultiplicityResult::getNuA(int k, int split) {
  double a = 0;
  double factor = 1;
  if(k >= 0 && k < multiplicities) {
    for(int i = k; i < multiplicities; i++) {
      factor = 1;
      if(i > 0) {
	for(int j = 0; j < k; j++) {
	  factor *= ( i - j);
	}
      }
      a += factor * getNormA(i, split);
    }
  }
  else {
    a = 0;
  }
  return a;
}

double NMSMultiplicityResult::getNu(int k, std::vector<double> &distribution) {
  double a = 0;
  double factor = 1;
  if(k >= 0 && k < distribution.size()) {
    for(int i = k; i < distribution.size(); i++) {
      factor = 1;
      if(i > 0) {
	for(int j = 0; j < k; j++) {
	  factor *= ( i - j);
	}
      }
      a += factor * distribution[i];
    }
  }
  else {
    a = 0;
  }
  return a;
}

double NMSMultiplicityResult::singles(double factor, int split) {
  return getRAsum(split) * factor;
}

double NMSMultiplicityResult::doubles(double factor, int split) {
  // std::streamsize ss = std::cout.precision();
  // std::cout << "singles in doubles" << std::setprecision(12) << singles(factor) << std::endl;
  // std::cout << std::setprecision(ss);
  return singles(factor, split) * (getNuRA(1, split) - getNuA(1, split));
}

double NMSMultiplicityResult::triples(double factor, int split) {
  // std::streamsize ss = std::cout.precision();
  // std::cout << "singles triples" << std::setprecision(12) << singles(factor) << std::endl;
  // std::cout << std::setprecision(ss);
  return singles(1, split) * (getNuRA(2, split) - getNuA(2, split) - 2 * getNuA(1, split) * (getNuRA(1, split) - getNuA(1, split))) / 2 * factor;
}

double NMSMultiplicityResult::avgsingles(double factor) {
  if(splits == -1) {
    return singles(factor);
  }
  else {
    double sum = 0;
    for(int i = 0; i < splits; i++) {
      sum += singles(factor, i) / splits;
    }
    return sum;
  }
}
double NMSMultiplicityResult::avgdoubles(double factor) {
  if(splits == -1) {
    return doubles(factor);
  }
  else {
    double sum = 0;
    for(int i = 0; i < splits; i++) {
      sum += doubles(factor, i) / splits;
    }
    return sum;
  }
}
double NMSMultiplicityResult::avgtriples(double factor) {
  if(splits == -1) {
    return triples(factor);
  }
  else {
    double sum = 0;
    for(int i = 0; i < splits; i++) {
      sum += triples(factor, i) / splits;
    }
    return sum;
  }
}

double NMSMultiplicityResult::uncsingles(double factor) {
  if(splits == -1) {
    return 0;
  }
  else {
    double sum = 0;
    for(int i = 0; i < splits; i++) {
      sum += (singles(factor, i) - avgsingles(factor)) * (singles(factor, i) - avgsingles(factor));
    }
    return sqrt(sum / splits);
  }
}

double NMSMultiplicityResult::uncdoubles(double factor) {
  if(splits == -1) {
    return 0;
  }
  else {
    double sum = 0;
    for(int i = 0; i < splits; i++) {
      sum += (doubles(factor, i) - avgdoubles(factor)) * (doubles(factor, i) - avgdoubles(factor));
    }
    return sqrt(sum / splits);
  }
}
double NMSMultiplicityResult::unctriples(double factor) {
  if(splits == -1) {
    return 0;
  }
  else {
    double sum = 0;
    for(int i = 0; i < splits; i++) {
      sum += (triples(factor, i) - avgtriples(factor)) * (triples(factor, i) - avgtriples(factor));
    }
    return sqrt(sum / splits);
  }
}

void NMSMultiplicityResult::debug() {
  double gate;
  gate = settings.registerLength * settings.registerPeriod;

  std::streamsize ss = std::cout.precision();
  std::cout << std::setprecision(12);
  std::cout << "gate " << gate << std::endl;
  std::cout << "predelay " << settings.predelay << std::endl;
  std::cout << "dieaway " << settings.dieaway << std::endl;
  std::cout << "measurement " << settings.measurementLength << std::endl;

  std::cout << "pu239 nu1 " << getNu(1, pu239fif) << std::endl;
  std::cout << "pu239 nu2 " << getNu(2, pu239fif) << std::endl;
  std::cout << "pu239 nu3 " << getNu(3, pu239fif) << std::endl;
  std::cout << "pu240 nu1 " << getNu(1, pu240sf) << std::endl;
  std::cout << "pu240 nu2 " << getNu(2, pu240sf) << std::endl;
  std::cout << "pu240 nu3 " << getNu(3, pu240sf) << std::endl;
  std::cout << "nuRA(1) " << getNuRA(1) << std::endl;
  std::cout << "nuA(1) " <<  getNuA(1) << std::endl;
  std::cout << "nuRA(2) " << getNuRA(2) << std::endl;
  std::cout << "nuA(2) " <<  getNuA(2) << std::endl;

  std::cout << "mulf fac1 " << exp(- settings.predelay / settings.dieaway) << std::endl;
  std::cout << "mulf fac2 " << (1 - exp( - gate / settings.dieaway)) << std::endl;

  std::cout << "mula " << mula() << std::endl;

  std::cout << std::setprecision(ss);
}

void NMSMultiplicityResult::setGateFractions(double fd, double ft) {
  settings.doublegatefraction = fd;
  settings.triplegatefraction = ft;
}

double NMSMultiplicityResult::simplegatefraction() {
  double gate;
  gate = settings.registerLength * settings.registerPeriod;
  return exp(- settings.predelay / settings.dieaway) * (1 - exp(-gate / settings.dieaway));
  //mulf[predelay_, gate_, dieaway_] :=  Exp[-predelay/dieaway] (1 - Exp[-gate/dieaway])
}

double NMSMultiplicityResult::mula() {
  double secondlength = settings.measurementLength / 1000000;
  double nominator = (-6 * triples(1 / secondlength) * getNu(2, pu240sf) * (getNu(1, pu239fif) - 1));
  double denominator1 = (settings.efficiency * settings.efficiency) * (settings.triplegatefraction) * singles(1 / secondlength);
  //double denominator1 = (settings.efficiency * settings.efficiency) * (mulf() * mulf()) * singles(1 / secondlength);
  double denominatorcommon (getNu(2, pu240sf) * getNu(3, pu239fif) - getNu(3, pu240sf) * getNu(2, pu239fif));
  //std::cout << "mula ml " << secondlength << std::endl;
  // std::cout << "mula nominator " << nominator << std::endl;
  // std::cout << "nu2 sf " << getNu(2, pu240sf) << std::endl;
  // std::cout << "nu1 if " << getNu(1, pu239fif) << std::endl;
  // std::cout << "eff " << settings.efficiency << std::endl;
  // std::cout << "mula denominator " << denominator1 << std::endl;
  return nominator / (denominator1 * denominatorcommon);
}

double NMSMultiplicityResult::mulb() {
  double nominator = 2 * doubles(1 / settings.measurementLength * 1000000) * (getNu(3, pu240sf) * (getNu(1, pu239fif) - 1) - 3 * getNu(2, pu240sf) * getNu(2, pu239fif));
  //  double denominator = settings.efficiency * mulf() * singles(1 / settings.measurementLength * 1000000) * (getNu(2, pu240sf) * getNu(3, pu239fif) - getNu(3, pu240sf) * getNu(2, pu239fif));
  double denominator = settings.efficiency * settings.doublegatefraction * singles(1 / settings.measurementLength * 1000000) * (getNu(2, pu240sf) * getNu(3, pu239fif) - getNu(3, pu240sf) * getNu(2, pu239fif));
  return nominator / denominator;
}

double NMSMultiplicityResult::mulc() {
  double nominator = 6 * doubles(1 / settings.measurementLength * 1000000) * getNu(2, pu240sf) * getNu(2, pu239fif);
  double denominator = settings.efficiency * settings.doublegatefraction * singles(1 / settings.measurementLength * 1000000) * (getNu(2, pu240sf) * getNu(3, pu239fif) - getNu(3, pu240sf) * getNu(2, pu239fif));
  return nominator / denominator - 1;
}

double NMSMultiplicityResult::mul() {
  double* sola = new double;
  double* solb = new double;
  double* solc = new double;

  // std::cout << "Mul a " << mula() << std::endl;
  // std::cout << "Mul b " << mulb() << std::endl;
  // std::cout << "Mul c " << mulc() << std::endl;

  int results = gsl_poly_solve_cubic (mulc(), mulb(), mula(), sola, solb, solc);

  if(results == 1) {
    return *sola;
  }
  else if(results == 3) {
    return *solc;
  }
  else {
    std::cout << "ERROR: Could not solve cubic function to calculate multiplication." << std::endl;
    exit(1);
  }
}

double NMSMultiplicityResult::fissionRate() {
  double part1 = (2 * doubles(1 / settings.measurementLength * 1000000)) / (settings.efficiency * settings.doublegatefraction);
  double part2 = (mul() * (mul() - 1) * singles( 1 / settings.measurementLength * 1000000)) / (getNu(1, pu239fif) - 1);
  double denominator = settings.efficiency * mul() * mul() * getNu(2, pu240sf);
  // std::cout << "F calc: dgf " << settings.doublegatefraction << std::endl;
  // std::cout << "F calc: part1 " << part1 << std::endl;
  return (part1 - part2) / denominator;
}

double NMSMultiplicityResult::alpha() {
  return singles(1 / settings.measurementLength * 1000000) / (fissionRate() * settings.efficiency * getNu(1, pu240sf) * mul()) - 1;
}

double NMSMultiplicityResult::pu240eff(double fissionpergramsecond) {
  return fissionRate() / fissionpergramsecond;
}

double NMSMultiplicityResult::pumass(double pu238, double pu240, double pu242) {
  double pueff = pu240eff();
  return pueff / ( 2.52 * pu238 + pu240 + 1.68 * pu242);
}


void NMSMultiplicityResult::exportToFile(std::string filename) {
  std::ofstream meta;
  std::ofstream table;

  std::string metafilename = filename;
  metafilename.append(".meta");
  std::string tablefilename = filename;
  tablefilename.append(".table");

  table.open(tablefilename.c_str());
  for(int i = 0; i < multiplicities; i++) {
    table << std::setw(13) << i << std::setw(16) << resultsRAptr[i] << std::setw(16) << resultsAptr[i] << std::setw(16) << normresultsRAptr[i] << std::setw(16) << normresultsAptr[i] << std::endl;
  }
  table.close();

  meta.open(metafilename.c_str());
  meta << getLastEvent() << std::endl;
  meta << getRAsum() << std::endl;
  meta << singles(1) << std::endl;
  meta << doubles(1) << std::endl;
  meta << triples(1) << std::endl;
  meta << singles(1000000 / settings.measurementLength) << std::endl;
  meta << doubles(1000000 / settings.measurementLength) << std::endl;
  meta << triples(1000000 / settings.measurementLength) << std::endl;
  meta << pu240eff() << std::endl;

  meta.close();
}

std::ostream& operator<<(std::ostream& os, const NMSMultiplicityResult& mr) {
  os << std::endl;
  os << std::setw(13) << "Multiplicity" << std::setw(16) << "R + A" << std::setw(16) << "A" << std::setw(16) << "R + A (norm)" << std::setw(16) << "A (norm)" << std::endl;;
  for(int i = 0; i < mr.multiplicities && i < 16; i++) {
    os << std::setw(13) << i << std::setw(16) << mr.resultsRAptr[i] << std::setw(16) << mr.resultsAptr[i] << std::setw(16) << mr.normresultsRAptr[i] << std::setw(16) << mr.normresultsAptr[i] << std::endl;
  }
  return os;
}
