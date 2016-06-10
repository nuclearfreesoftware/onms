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

/** \brief Class to store and evaluate multiplicity coincidence results
 *
 *  This class should handle coincidence counting results.
 *  A simple output and different calculational methods are available.
 */

#ifndef NMSMultiplicityResult_h
#define NMSMultiplicityResult_h 1

#include "NMSMultiplicityMeasurementSetting.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h> 

#include <gsl/gsl_poly.h>

class NMSMultiplicityResult {

public:
  NMSMultiplicityResult(int multiplicitiesInit = 128, int splits = -1);

  ~NMSMultiplicityResult();

  NMSMultiplicityResult(const NMSMultiplicityResult& x);
  NMSMultiplicityResult& operator = (const NMSMultiplicityResult& x);

  void SetSettings(NMSMultiplicityMeasurementSetting newset);
  NMSMultiplicityMeasurementSetting GetSettings() { return settings; };

  int getRA(int i, int split = -1);
  int getA(int i, int split = -1);

  void addRA(int i, int increment = 1, int split = -1);
  void addA(int i, int increment = 1, int split = -1);
  void set(int i, int ra, int a);
  void setRA(int i, int ra);
  void setA(int i, int a);

  void reset();

  void normalize();
  double getNormRA(int i, int split = -1);
  double getNormA(int i, int split = -1);

  int getRAsum(int split = -1);
  int getAsum(int split = -1);

  void setLastEvent(double le);
  double getLastEvent();

  double getNuRA(int k, int split = -1);
  double getNuA(int k, int split = -1);

  double getNu(int k, std::vector<double> &distribution);

  double singles(double factor, int split = -1);
  double doubles(double factor, int split = -1);
  double triples(double factor, int split = -1);

  double avgsingles(double factor);
  double avgdoubles(double factor);
  double avgtriples(double factor);
  double uncsingles(double factor);
  double uncdoubles(double factor);
  double unctriples(double factor);

  void setGateFractions(double fd, double ft);
  double simplegatefraction();
  double getdoublegatefraction() { return settings.doublegatefraction; };
  double gettriplegatefraction() { return settings.triplegatefraction; };
  
  double mula();
  double mulb();
  double mulc();
  double mul();
  double fissionRate();
  double alpha();
  double pu240eff(double fissionpergramsecond = 479); // From Thesis! (479(17) fissions/ gram second)
  double pumass(double pu238, double pu240, double pu242);
  
  int getMultiplicities() { return multiplicities; }
  
  void debug();
  void exportToFile(std::string filename);

  friend std::ostream& operator<<(std::ostream& os, const NMSMultiplicityResult& mr);

private:
  int multiplicities;
  int splits;
  int* resultsRAptr;
  int* resultsAptr;
  double* normresultsRAptr;
  double* normresultsAptr;
  double lastEvent;

  std::vector< std::vector < int > > splitresultsRAvec;
  std::vector< std::vector < int > > splitresultsAvec;
  std::vector< std::vector < double > > splitnormresultsRAvec;
  std::vector< std::vector < double > > splitnormresultsAvec;

  NMSMultiplicityMeasurementSetting settings;

  std::vector<double> pu240sf;
  std::vector<double> pu239fif;

};


#endif
