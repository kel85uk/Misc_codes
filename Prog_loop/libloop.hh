/* Copyright (C) 2014  kklloh
Lib-loop (Library to loop over system commands) - A Planck-length_scale project for a friend

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#ifndef LIB_LOOP
#define LIB_LOOP
#include <fstream>				/* std::fstream */
#include <iostream>			/* std::cout, std::endl */
#include <stdlib.h>			/* system, NULL, EXIT_FAILURE */

using std::cout;
using std::endl;
using std::fstream;

bool lib_loop(char *program1, char *program2, char *program3, char *program4, char *filename_buffer, double TOL_max, int itermax)
{
  int i;
  cout << "Checking if processor is available... \n";
  fstream fs;
  if (system(NULL)) puts ("Ok");
    else exit (EXIT_FAILURE);
  cout << "From wrapper: Executing command init...\n";
  	i = system(program1); // Change these to the windows system commands, if I remembered correctly it's del instead of rm
	i = system(program2); // Your first initialization program
	if (i != 0) {cout << "Error in initialization... \n"; return false;} // Check for any errors in the command
	double tol = 0.;
  for (int ii = 0; (ii < itermax) && (tol < TOL_max); ++ii){
	  i = system (program3);
		i = system (program4);
		fs.open(filename_buffer);
		fs >> tol;
		fs.close();
		cout << tol << endl;		
		if (i != 0) {cout << "Error in loop execution... \n"; return false;}		
	}
  return true;
}
#endif
