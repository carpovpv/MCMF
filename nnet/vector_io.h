
/*
 * Never write like this!!!!
 * Vector binary io for prognosis.
 * 22.12.2007
 */
#ifndef __VECTOR_IO
#define __VECTOR_IO

#include <vector>
#include <iostream>
#include <fstream>
#include <string.h>
#include <typeinfo>

using namespace std;


template <class T>
void 
vector_write(std::ofstream& out, const vector<T> *data)
{
	unsigned long N = data->size();
	out.write(reinterpret_cast<const char*>(&N), sizeof(unsigned long));

	int s = !strcmp(typeid(vector<string>).name(),typeid(*data).name());

	for(unsigned long i=0;i<N;++i)
	{
		if(s){
			string *str = (string *) &data->at(i);
			int n = str->length();
			out.write(reinterpret_cast<const char*>(&n), sizeof(int));
			out.write(str->c_str(),n);
			continue;
		}
		out.write(reinterpret_cast<const char*>(&data->at(i)), sizeof(T));
	}	
}

template <class T>
void
vector_read(std::ifstream& in,vector<T> *data)
{
	unsigned long N;
	in.read(reinterpret_cast<char*>(&N), sizeof(unsigned long));
	int s = !strcmp(typeid(vector<string>).name(),typeid(*data).name());
	
	char *str_buf;
	int str_buf_len = 2048;
	vector<string> *s_data = (vector <string> *) data;
	if (s)	str_buf = new char[str_buf_len];
	for(unsigned long i=0;i<N;i++)
	{
		T buffer;
		if(s)
		{
			int str_length;
			in.read(reinterpret_cast<char*>(&str_length), sizeof(int));
			if(str_length > str_buf_len && str_length <= 2*str_buf_len)
			{
				str_buf_len *= 2;
				delete [] str_buf;
				str_buf = new char[str_buf_len];
			}
			else
			{
				str_buf_len = str_length * 2;
				delete [] str_buf;
				str_buf = new char[str_buf_len];
			}
		
			in.read(str_buf, str_length);
			str_buf[str_length]='\0';
			s_data->push_back(str_buf);
			continue;
		}
			
		in.read(reinterpret_cast<char*>(&buffer), sizeof(T));
		data->push_back(buffer);
	}
	if(s)
		delete [] str_buf;
}
#endif
