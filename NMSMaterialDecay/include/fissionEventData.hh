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

#ifndef FISSIONEVENTDATA_H
#define FISSIONEVENTDATA_H 1

#define SP_FISSION_ISOTOPES 16
#define SP_FISSION_N 11
#define SP_FISSION_NUBAR_ISOTOPES 18

class fissionEventData 
{
public:
  fissionEventData();
  virtual ~fissionEventData();

  static int sfDataIndex(int isotope, int Cf252option = 0);
  static double fissionEventNu(int isotope, int n, int Cf252option = 0);
  static double getSfNubar(int isotope, int Cf252option = 0);

private:
  static double sfnu[SP_FISSION_ISOTOPES][SP_FISSION_N];

  static int spzaid[SP_FISSION_NUBAR_ISOTOPES];
  static double spnubar[SP_FISSION_NUBAR_ISOTOPES];
};


#endif /* FISSIONEVENTDATA_H */
