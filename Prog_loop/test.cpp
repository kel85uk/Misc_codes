/* Copyright (C) 2014  kklloh

Test code for libloop library
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

#include "libloop.hh"
#include <iostream>

int main (int argc, char *argv[])
{
	char deletetxt[] = "rm -rvf *.txt";
	char initprog[] = "./init";
	char sumprog[] = "./sum";
	char prodprog[] = "./product";
	char buffer_file[] = "results_prod.txt";
	int	itermax = 100;
	double TOLmax = 100;
	bool return_val = lib_loop(deletetxt,initprog,sumprog,prodprog,buffer_file,TOLmax,itermax);
	if(return_val){
		std::cout << "All ok!\n";
		return 0;
	}
	else{
		std::cout << "Some parts in the wrapper library failed!\n";
		return 1;
	}
}
