#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <algorithm>
using std::tuple;
using std::vector;
using std::get;
using std::sort;
using std::make_tuple;
using std::cout;
using std::endl;
 
typedef tuple<int,int,double> mytuple;
 
bool mycompare (const mytuple &lhs, const mytuple &rhs){
	if(get<0>(lhs) != get<0>(rhs)) return get<0>(lhs) < get<0>(rhs);
  else return get<1>(lhs) < get<1>(rhs);
}

void mat_set(int i, int j, double value, vector<mytuple>& data);

int main(void){
  vector<mytuple> data;
  int i, j; double value;
	mat_set(1,1,3,data);
	mat_set(0,1,2,data);
	mat_set(1,2,3,data);
	mat_set(0,0,6,data);
	mat_set(1,0,2,data);
	vector<int> rowind, colind;

  for(vector<mytuple>::iterator iter = data.begin(); iter != data.end(); iter++){
  		rowind.push_back(get<0>(*iter));
  		colind.push_back(get<1>(*iter));
  }
  for(auto iter1 = rowind.begin(); iter1 != rowind.end(); iter1++){
  		i = iter1 - rowind.begin();
  		cout << (*iter1) << "\t" << colind[i] << "\t" << get<2>(data[i]) << endl;
	}
}

void mat_set(int i, int j, double value, vector<mytuple>& data){
	data.push_back(make_tuple(i,j,value));
	sort(data.begin(),data.end(),mycompare);
}
