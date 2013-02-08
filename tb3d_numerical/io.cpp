// ****************************************************************************
// Input/Output for functions defined on grids
//
// Description:	
//		Additionally to the function the parameters of the grid
//		is saved. Therefore it is possible to read in the function
//		on a later time exactly.
//
// Programm structure: 
//		The first lines of each data file contain the grid parameters.
//		Of course the number of lines depend on the type of grid.
//		These lines are uncommented by '#' so that xmgrace or gnuplot
//		will ignore them.	
// 
// Tobias Stollenwerk, Last Modification: 19.2.2010
// ****************************************************************************
#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<fstream>
#include<vector>
#include<string>
#include<complex>
#include"grid.h"
#include"multigrid.h"
#include"io.h"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

using namespace std;

// ************************************************************************
// *************  create numbered filenames to avoid overwriting  *********
// ************************************************************************
string appFilename(string  folder, string  title)
{
	int count=0;
	string count_str;
	char count_ch[100];
	bool contFlag=true;
	int Nmax=10000;
	string fileName;
	ifstream in;
	while (count < Nmax && contFlag)
	{
		sprintf(count_ch, "%04i", count);
		count_str=count_ch;
		fileName=folder + "/" + title + "_" + count_str + ".dat"; 
		in.open(fileName.c_str());
		if (!(in.fail()))
		{
			count++;
		}
		else
		{
			contFlag=false;
		}
		in.close();
	}
	if (count == Nmax)
	{
		cerr << "Error appFilename. Counter exceeded. Break" << endl;
		exit(1);
	}
	fileName=folder + "/" + title + "_" + count_str;
	return fileName;
}
// ************************************************************************
// *************  input/output for generic grids  *************************
// ************************************************************************
void save (double x, string  folder, string  title)
{
	string datName= folder + "/" + title + ".dat";
	ofstream out;
	out.open(datName.c_str());
	out << scientific << setprecision(17) <<  x << endl;
	out.close();
}
void save ( double value , string  folder ,string  title, int counter)
{
	string datName=folder + "/" + title + ".dat";
	ofstream out;
	out.open(datName.c_str(), ios::app);
	out << scientific << setprecision(17) << counter << "\t" << value << endl;
	out.close();
}
void save (double x, double y,  string  folder, string  title)
{
	string datName= folder + "/" + title + ".dat";
	ofstream out;
	out.open(datName.c_str());
	out << scientific << setprecision(17) <<  x << "\t" << y <<  endl;
	out.close();
}
void save ( vector <double>& A, vector <double>& omega, int N, string  folder,  string  title)
{
	string datName= folder + "/" + title + ".dat";
	ofstream out;
	out.open(datName.c_str());
	for (int i=0; i!=N+1; i++)
	{
		out << scientific << setprecision(17) <<  omega[i] << "\t" << A[i] << endl;
	}
	out.close();
}
void save ( vector <double>& A, grid & dgrid, string  folder, string  title)
{
	string datName= folder + "/" + title + ".dat";
	ofstream out;
	out.open(datName.c_str());
	out << scientific << setprecision(17);
	for (int i=0; i!=dgrid.M+1; i++)
	{
		out << scientific << setprecision(17) <<  dgrid.omega[i] << "\t" << A[i] << endl;
	}
	out.close();
}
void save ( vector <double>& A, grid & dgrid, string  folder, string  title, int counter)
{
	char counter_ch[10];
	string counter_str;
	sprintf(counter_ch, "%04i", counter);
	counter_str=counter_ch;	
	
	string datName=folder + "/" + title + "_" + counter_str + ".dat";
	ofstream out;
	out.open(datName.c_str());
	out << scientific << setprecision(17);
	for (int i=0; i!=dgrid.M+1; i++)
	{
		out << scientific << setprecision(17) <<  dgrid.omega[i] << "\t" << A[i] << endl;
	}
	out.close();
}
void savePseudoPeakPara ( vector <double> & peakPara, string  folder, int counter)
{
	char counter_ch[10];
	string counter_str;
	sprintf(counter_ch, "%04i", counter);
	counter_str=counter_ch;	
	
	string datName=folder + "/pseudoPeakPara_" + counter_str + ".dat";
	ofstream out;
	out.open(datName.c_str());
	out << "\t" << "peakWidth_fu\t\t" << scientific << setprecision(17) <<  peakPara[0] <<  endl;
	out << "\t" << "peakWidth_fd\t\t" << scientific << setprecision(17) <<  peakPara[1] <<  endl;
	out << "\t" << "peakWidth_b \t\t" << scientific << setprecision(17) <<  peakPara[2] <<  endl;
	out << "\t" << "omegaPeak_fu\t\t" << scientific << setprecision(17) <<  peakPara[3] <<  endl;
	out << "\t" << "omegaPeak_fd\t\t" << scientific << setprecision(17) <<  peakPara[4] <<  endl;
	out << "\t" << "omegaPeak_b \t\t" << scientific << setprecision(17) <<  peakPara[5] <<  endl;
	out << "\t" << "omegaPeak_fu_0\t\t" << scientific << setprecision(17) <<  peakPara[6] <<  endl;
	out << "\t" << "omegaPeak_fd_0\t\t" << scientific << setprecision(17) <<  peakPara[7] <<  endl;
	out << "\t" << "omegaPeak_b_0 \t\t" << scientific << setprecision(17) <<  peakPara[8] <<  endl;
	out << "\t" << "omegaPeak_fu_1\t\t" << scientific << setprecision(17) <<  peakPara[9] <<  endl;
	out << "\t" << "omegaPeak_fd_1\t\t" << scientific << setprecision(17) <<  peakPara[10] <<  endl;
	out << "\t" << "omegaPeak_b_1 \t\t" << scientific << setprecision(17) <<  peakPara[11] <<  endl;
	out << "\t" << "NPeak_fu \t\t" << scientific << setprecision(17) <<  peakPara[12] <<  endl;
	out << "\t" << "NPeak_fd \t\t" << scientific << setprecision(17) <<  peakPara[13] <<  endl;
	out << "\t" << "NPeak_b  \t\t" << scientific << setprecision(17) <<  peakPara[14] <<  endl;
	out.close();
}
void savePeakPara ( vector < vector <double> > & peakPara, string  folder, int counter)
{
	char counter_ch[10];
	string counter_str;
	sprintf(counter_ch, "%04i", counter);
	counter_str=counter_ch;	
	
	string datName=folder + "/peakPara_" + counter_str + ".dat";
	ofstream out;
	out.open(datName.c_str());
	out << "#layer\tparameter\tvalue" << endl;
	for (int n=0; n!=peakPara.size(); n++)
	{
		out << n << "\t" << "adu peakWidth:\t" << peakPara[n][0] << endl;
		out << n << "\t" << "adu omegaPeak:\t" << peakPara[n][1] << endl;
		out << n << "\t" << "adu omegaPeak_0:\t" << peakPara[n][2] << endl;
		out << n << "\t" << "adu omegaPeak_1:\t" << peakPara[n][3] << endl;
		out << n << "\t" << "adu NPeak:\t\t" << int(peakPara[n][4]) << endl;
		out << n << "\t" << "add peakWidth:\t" << peakPara[n][5] << endl;
		out << n << "\t" << "add omegaPeak:\t" << peakPara[n][6] << endl;
		out << n << "\t" << "add omegaPeak_0:\t" << peakPara[n][7] << endl;
		out << n << "\t" << "add omegaPeak_1:\t" << peakPara[n][8] << endl;
		out << n << "\t" << "add NPeak:\t\t" << int(peakPara[n][9]) << endl;
	}
	out.close();
}
void read (double & x, string fileName)
{
	char temp[100];
	ifstream in;
	// check if file exist
	in.open(fileName.c_str());
	if (in.fail())
	{
		cerr << endl;
		cerr << "Error while trying to open file: " << fileName << endl;
		cerr << "Break" << endl;
		exit(1);
	}	
	in.close();
	in.open(fileName.c_str());
	in.getline(temp,100);
	x=atof(temp);
	in.close();
}

void read( vector <double>& A, string fileName)
{
	char temp[100];
	ifstream in;
	in.open(fileName.c_str());
	if (in.fail())
	{
		cerr << endl;
		cerr << "Error while trying to open file: " << fileName << endl;
		cerr << "Break" << endl;
		exit(1);
	}	
	in.close();
	in.open(fileName.c_str());
	int n=0;
	while (!in.eof())
	{
		in.getline(temp,100);
		n++;
	}
	if (n-1!=A.size())
	{
		cerr << endl;
		cerr << "Number of data points is incorrect in file: " << fileName << endl;
		cerr << n-1 << " != " << A.size() << endl;
		cerr << "Break" << endl;
		exit(1);
	}
	in.close();

	in.open(fileName.c_str());
	for (int i=0; i!=A.size(); i++)
	{
		in.getline(temp, 100, '\t');
		in.getline(temp, 100);
		A[i]=atof(temp);
	}
	in.close();
}

// read in from a unknow grid and interpolate on grid: mgrid
void readInter( vector <double>& A, string  fileName, grid & mgrid)
{
	char temp[100];
	ifstream in;
	in.open(fileName.c_str());
	if (in.fail())
	{
		cerr << endl;
		cerr << "Error while trying to open file: " << fileName << endl;
		cerr << "Break" << endl;
		exit(1);
	}	
	in.close();
	in.open(fileName.c_str());
	int n_com=0;
	int n_data=0;
	while (!in.eof())
	{
		in.getline(temp,100);
		if (temp[0]=='#')
		{
			n_com++;
		}
		else
		{
			n_data++;
		}
	}
	in.close();
	n_data-=1;

	vector <double> omega1(n_data);
	vector <double> A1(n_data);
	in.open(fileName.c_str());
	for (int i=0; i!=n_com; i++)
	{
		in.ignore(500,'\n');
	}
	for (int i=0; i!=n_data; i++)
	{
		in.getline(temp, 100, '\t');
		omega1[i]=atof(temp);
		in.getline(temp, 100);
		A1[i]=atof(temp);
	}
	in.close();

	// interpolate
	double omegaT, omega_I,omega_Ip, f_I, f_Ip;
	int I;
	for (int j=0; j<=mgrid.M; j++)
	{
		omegaT=mgrid.omega[j];

		// if out of range, set to zero
		if (omegaT<omega1[0])
		{
			A[j]=0;	
		}
		else if (omegaT>omega1[n_data-1])
		{
			A[j]=0;	
		}
		else if (omegaT==omega1[0])
		{
			A[j]=A1[0];
		}
		else if (omegaT==omega1[n_data-1])
		{
			A[j]=A1[n_data-1];
		}
		else
		{
			// calculate inverse of grid (brute force)
			I=0;
			while (I<n_data-1 && omegaT >= omega1[I])
			{
				I++;
			}
			I--;
			omega_I = omega1[I];	
			//cout << omega_I << " <-> " << omegaT << endl;
			f_I=A1[I];
			omega_Ip=omega1[I+1];
			f_Ip=A1[I+1];
			A[j]=f_I+ (f_Ip - f_I) * ((omegaT-omega_I)/(omega_Ip - omega_I));	
		}
	}
}
// read in from a know grid (igrid) and interpolate on other grid (mgrid)
void read( vector <double>& A, string  fileName, grid & mgrid, grid & igrid)
{
	if (A.size()!=mgrid.M+1)
	{
		cerr << endl;
		cerr << "Error while trying to read file: "<< fileName << endl;
		cerr << "Size of container does (" << A.size() << ") not match the given grid length (" << mgrid.M+1 << ")." << endl;
		cerr << "Break" << endl;
		exit(1);
	}
	char temp[100];
	ifstream in;
	in.open(fileName.c_str());
	if (in.fail())
	{
		cerr << endl;
		cerr << "Error while trying to open file: " << fileName << endl;
		cerr << "Break" << endl;
		exit(1);
	}	
	in.close();
	in.open(fileName.c_str());
	int n_com=0;
	int n_data=0;
	while (!in.eof())
	{
		in.getline(temp,100);
		if (temp[0]=='#')
		{
			n_com++;
		}
		else
		{
			n_data++;
		}
	}
	in.close();
	n_data-=1;

	if (n_data!=igrid.M+1)
	{
		cerr << endl;
		cerr << "Read: Number of data points is incorrect in file: " << fileName << endl;
		cerr << n_data << " != " << igrid.M+1 << endl;
		cerr << "Break" << endl;
		exit(1);
	}

	vector <double> A1(n_data);
	in.open(fileName.c_str());
	for (int i=0; i!=n_com; i++)
	{
		in.ignore(500,'\n');
	}
	for (int i=0; i!=n_data; i++)
	{
		in.getline(temp, 100, '\t');
		in.getline(temp, 100);
		A1[i]=atof(temp);
	}
	in.close();

	// interpolate
	double omegaT, omega_I,omega_Ip, f_I, f_Ip;
	int I;
	for (int j=0; j!=mgrid.M+1; j++)
	{
		omegaT=mgrid.omega[j];
		if (omegaT>=igrid.omega_min && omegaT<=igrid.omega_max)
		{
			I=igrid.inverse(omegaT);
			if (igrid.omega[I]<omegaT)
			{
				A[j]=A1[I]+ (A1[I+1] - A1[I])/(igrid.omega[I+1] - igrid.omega[I]) * (omegaT-igrid.omega[I]);	
			}
			else if (igrid.omega[I]>omegaT)
			{
				A[j]=A1[I]+ (A1[I-1] - A1[I])/(igrid.omega[I-1] - igrid.omega[I]) * (omegaT-igrid.omega[I]);	
			}
			else
			{
				A[j]=A1[I];
			}	
		}
		else
		{
			A[j]=0;	
		}
	}
}


// ************************************************************************
// ************* input/output for multigrids ******************************
// ************************************************************************
void save ( vector <double>& A, const multigrid & mgrid, string  folder, string  title)
{
	string datName= folder + "/" + title + ".dat";
	ofstream out;
	out.open(datName.c_str());
	out << scientific << setprecision(17);
	/*
	out << "##MULTGRID-PARAMETERS################" << endl; 
	out << "# N_fgr:    \t" << mgrid.fgridRegions.size() << endl;
	out << "# N_sgr:    \t" << mgrid.sgridRegions.size() << endl;
	out << "#" << endl; 
	out << "# FUNDAMENTAL GRID REGIONS:" << endl; 
	for (int n=0; n<int(mgrid.fgridRegions.size()); n++)
	{
		out << "#" << endl; 
		out << "# type:     \t" << mgrid.fgridRegions[n].type << endl; 
		out << "# id:       \t" << mgrid.fgridRegions[n].id << endl; 
		if (mgrid.fgridRegions[n].type=="log")
		{
			out << "# N_l:      \t" << mgrid.fgridRegions[n].para_i[0] << endl; 
			out << "# N_r:      \t" << mgrid.fgridRegions[n].para_i[1] << endl; 
			out << "# omega_l:  \t" << mgrid.fgridRegions[n].omega_l << endl; 
			out << "# omega_r:  \t" << mgrid.fgridRegions[n].omega_r << endl; 
			out << "# omegak:   \t" << mgrid.fgridRegions[n].para_d[0] << endl; 
			out << "# omegak_0: \t" << mgrid.fgridRegions[n].para_d[1] << endl; 
		}
		else if (mgrid.fgridRegions[n].type=="equi")
		{
			out << "# N_e:      \t" << mgrid.fgridRegions[n].para_i[0] << endl; 
			out << "# omega_l:  \t" << mgrid.fgridRegions[n].omega_l << endl; 
			out << "# omega_r:  \t" << mgrid.fgridRegions[n].omega_r << endl; 
			out << "# omega_c:  \t" << mgrid.fgridRegions[n].omega_c << endl; 
		}
		else if (mgrid.fgridRegions[n].type=="tan")
		{
			out << "# N_t:      \t" << mgrid.fgridRegions[n].para_i[0] << endl; 
			out << "# omega_l:  \t" << mgrid.fgridRegions[n].omega_l << endl; 
			out << "# omega_r:  \t" << mgrid.fgridRegions[n].omega_r << endl; 
			out << "# omega_c:  \t" << mgrid.fgridRegions[n].para_d[0] << endl; 
			out << "# c:        \t" << mgrid.fgridRegions[n].para_d[1] << endl; 
		}
		else 
		{
			cerr << "Error: save multigrid. fgridRegion type is wrong. Break" << endl;
			throw 1;
		}
	}
	out << "#" << endl; 
	out << "# SPECIAL GRID REGIONS:" << endl; 
	for (int n=0; n<int(mgrid.sgridRegions.size()); n++)
	{
		out << "#" << endl; 
		out << "# type:     \t" << mgrid.sgridRegions[n].type << endl; 
		out << "# id:       \t" << mgrid.sgridRegions[n].id << endl; 
		if (mgrid.sgridRegions[n].type=="log")
		{
			out << "# N_l:      \t" << mgrid.sgridRegions[n].para_i[0] << endl; 
			out << "# N_r:      \t" << mgrid.sgridRegions[n].para_i[1] << endl; 
			out << "# omega_l:  \t" << mgrid.sgridRegions[n].omega_l << endl; 
			out << "# omega_r:  \t" << mgrid.sgridRegions[n].omega_r << endl; 
			out << "# omegak:   \t" << mgrid.sgridRegions[n].para_d[0] << endl; 
			out << "# omegak_0: \t" << mgrid.sgridRegions[n].para_d[1] << endl; 
		}
		else if (mgrid.sgridRegions[n].type=="equi")
		{
			out << "# N_e:      \t" << mgrid.sgridRegions[n].para_i[0] << endl; 
			out << "# omega_l:  \t" << mgrid.sgridRegions[n].omega_l << endl; 
			out << "# omega_r:  \t" << mgrid.sgridRegions[n].omega_r << endl; 
			out << "# omega_c:  \t" << mgrid.sgridRegions[n].omega_c << endl; 
		}
		else if (mgrid.sgridRegions[n].type=="tan")
		{
			out << "# N_t:      \t" << mgrid.sgridRegions[n].para_i[0] << endl; 
			out << "# omega_l:  \t" << mgrid.sgridRegions[n].omega_l << endl; 
			out << "# omega_r:  \t" << mgrid.sgridRegions[n].omega_r << endl; 
			out << "# omega_c:  \t" << mgrid.sgridRegions[n].para_d[0] << endl; 
			out << "# c:        \t" << mgrid.sgridRegions[n].para_d[1] << endl; 
		}
		else 
		{
			cerr << "Error: save multigrid. sgridRegion type is wrong. Break" << endl;
			throw 1;
		}
	}
	out << "#####################################" << endl; 
	*/
	for (int i=0; i!=mgrid.M+1; i++)
	{
		out << scientific << setprecision(17) <<  mgrid.omega[i] << "\t" << A[i] << endl;
	}
	out.close();
}
void save ( vector <double>& A, const multigrid & mgrid, string  folder, string  title, int counter)
{
	char counter_ch[10];
	string counter_str;
	sprintf(counter_ch, "%04i", counter);
	counter_str=counter_ch;	
	string datName=folder + "/" + title + "_" + counter_str + ".dat";
	ofstream out;
	out.open(datName.c_str());
	out << scientific << setprecision(17);
	/*
	out << "##MULTGRID-PARAMETERS################" << endl; 
	out << "# N_fgr:    \t" << mgrid.fgridRegions.size() << endl;
	out << "# N_sgr:    \t" << mgrid.sgridRegions.size() << endl;
	out << "#" << endl; 
	out << "# FUNDAMENTAL GRID REGIONS:" << endl; 
	for (int n=0; n<int(mgrid.fgridRegions.size()); n++)
	{
		out << "#" << endl; 
		out << "# type:     \t" << mgrid.fgridRegions[n].type << endl; 
		out << "# id:       \t" << mgrid.fgridRegions[n].id << endl; 
		if (mgrid.fgridRegions[n].type=="log")
		{
			out << "# N_l:      \t" << mgrid.fgridRegions[n].para_i[0] << endl; 
			out << "# N_r:      \t" << mgrid.fgridRegions[n].para_i[1] << endl; 
			out << "# omega_l:  \t" << mgrid.fgridRegions[n].omega_l << endl; 
			out << "# omega_r:  \t" << mgrid.fgridRegions[n].omega_r << endl; 
			out << "# omegak:   \t" << mgrid.fgridRegions[n].para_d[0] << endl; 
			out << "# omegak_0: \t" << mgrid.fgridRegions[n].para_d[1] << endl; 
		}
		else if (mgrid.fgridRegions[n].type=="equi")
		{
			out << "# N_e:      \t" << mgrid.fgridRegions[n].para_i[0] << endl; 
			out << "# omega_l:  \t" << mgrid.fgridRegions[n].omega_l << endl; 
			out << "# omega_r:  \t" << mgrid.fgridRegions[n].omega_r << endl; 
			out << "# omega_c:  \t" << mgrid.fgridRegions[n].omega_c << endl; 
		}
		else if (mgrid.fgridRegions[n].type=="tan")
		{
			out << "# N_t:      \t" << mgrid.fgridRegions[n].para_i[0] << endl; 
			out << "# omega_l:  \t" << mgrid.fgridRegions[n].omega_l << endl; 
			out << "# omega_r:  \t" << mgrid.fgridRegions[n].omega_r << endl; 
			out << "# omega_c:  \t" << mgrid.fgridRegions[n].para_d[0] << endl; 
			out << "# c:        \t" << mgrid.fgridRegions[n].para_d[1] << endl; 
		}
		else 
		{
			cerr << "Error: save multigrid. fgridRegion type is wrong. Break" << endl;
			throw 1;
		}
	}
	out << "#" << endl; 
	out << "# SPECIAL GRID REGIONS:" << endl; 
	for (int n=0; n<int(mgrid.sgridRegions.size()); n++)
	{
		out << "#" << endl; 
		out << "# type:     \t" << mgrid.sgridRegions[n].type << endl; 
		out << "# id:       \t" << mgrid.sgridRegions[n].id << endl; 
		if (mgrid.sgridRegions[n].type=="log")
		{
			out << "# N_l:      \t" << mgrid.sgridRegions[n].para_i[0] << endl; 
			out << "# N_r:      \t" << mgrid.sgridRegions[n].para_i[1] << endl; 
			out << "# omega_l:  \t" << mgrid.sgridRegions[n].omega_l << endl; 
			out << "# omega_r:  \t" << mgrid.sgridRegions[n].omega_r << endl; 
			out << "# omegak:   \t" << mgrid.sgridRegions[n].para_d[0] << endl; 
			out << "# omegak_0: \t" << mgrid.sgridRegions[n].para_d[1] << endl; 
		}
		else if (mgrid.sgridRegions[n].type=="equi")
		{
			out << "# N_e:      \t" << mgrid.sgridRegions[n].para_i[0] << endl; 
			out << "# omega_l:  \t" << mgrid.sgridRegions[n].omega_l << endl; 
			out << "# omega_r:  \t" << mgrid.sgridRegions[n].omega_r << endl; 
			out << "# omega_c:  \t" << mgrid.sgridRegions[n].omega_c << endl; 
		}
		else if (mgrid.sgridRegions[n].type=="tan")
		{
			out << "# N_t:      \t" << mgrid.sgridRegions[n].para_i[0] << endl; 
			out << "# omega_l:  \t" << mgrid.sgridRegions[n].omega_l << endl; 
			out << "# omega_r:  \t" << mgrid.sgridRegions[n].omega_r << endl; 
			out << "# omega_c:  \t" << mgrid.sgridRegions[n].para_d[0] << endl; 
			out << "# c:        \t" << mgrid.sgridRegions[n].para_d[1] << endl; 
		}
		else 
		{
			cerr << "Error: save multigrid. sgridRegion type is wrong. Break" << endl;
			throw 1;
		}
	}
	out << "#####################################" << endl; 
	*/
	for (int i=0; i!=mgrid.M+1; i++)
	{
		out << scientific << setprecision(17) <<  mgrid.omega[i] << "\t" << A[i] << endl;
	}
	out.close();
}
void readGrid(string filename, multigrid & igrid)
{
	// initialize multigrid
	multigrid mgrid;

	// check if file exists
	char temp[100];
	ifstream in;
	in.open(filename.c_str());
	if (in.fail())
	{
		cerr << endl;
		cerr << "Error while trying to open file: " << filename << endl;
		cerr << "Break" << endl;
		exit(1);
	}	
	in.close();

	// read in grid parameter and create grid regions
	in.open(filename.c_str());
	in.ignore(100, '\n');
	in.getline(temp, 100, '\t');
	in.getline(temp, 100);
	int N_fgr=atoi(temp);
	in.getline(temp, 100, '\t');
	in.getline(temp, 100);
	int N_sgr=atoi(temp);

	// fundamental grid regions
	in.ignore(100, '\n');
	in.ignore(100, '\n');
	for (int n=0; n<N_fgr; n++)
	{
		in.ignore(100, '\n');
		in.getline(temp, 100, '\t');
		in.getline(temp, 100);
		string type=temp;
		in.getline(temp, 100, '\t');
		in.getline(temp, 100);
		string id=temp;
		if (type=="log")
		{
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			int N_l=atoi(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			int N_r=atoi(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			double omega_l=atof(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			double omega_r=atof(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			double omegak=atof(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			double omegak_0=atof(temp);

			// create grid region
			mgrid.add_gr_log(N_l, N_r, omega_l, omega_r, omegak, omegak_0, id);
		}
		else if (type=="equi")
		{
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			int N_e=atoi(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			double omega_l=atof(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			double omega_r=atof(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			double omega_c=atof(temp);

			// create grid region
			mgrid.add_gr_equi(N_e, omega_l, omega_r, omega_c, id);
		}
		else if (type=="tan")
		{
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			int N_t=atoi(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			double omega_l=atof(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			double omega_r=atof(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			double omega_c=atof(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			double c=atof(temp);

			// create grid region
			mgrid.add_gr_tan(N_t, omega_l, omega_r, omega_c, c, id);
		}
		else 
		{
			cerr << "Error: read multigrid. fgridregion type '" << type << "' is wrong: . Break" << endl;
			throw 1;
		}
	}

	// special grid regions
	in.ignore(100, '\n');
	in.ignore(100, '\n');
	for (int n=0; n<N_sgr; n++)
	{
		in.ignore(100, '\n');
		in.getline(temp, 100, '\t');
		in.getline(temp, 100);
		string type=temp;
		in.getline(temp, 100, '\t');
		in.getline(temp, 100);
		string id=temp;
		if (type=="log")
		{
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			int N_l=atoi(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			int N_r=atoi(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			double omega_l=atof(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			double omega_r=atof(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			double omegak=atof(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			double omegak_0=atof(temp);

			// create grid region
			mgrid.add_sgr_log(N_l, N_r, omega_l, omega_r, omegak, omegak_0, id);
		}
		else if (type=="equi")
		{
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			int N_e=atoi(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			double omega_l=atof(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			double omega_r=atof(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			double omega_c=atof(temp);

			// create grid region
			mgrid.add_sgr_equi(N_e, omega_l, omega_r, omega_c, id);
		}
		else if (type=="tan")
		{
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			int N_t=atoi(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			double omega_l=atof(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			double omega_r=atof(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			double omega_c=atof(temp);
			in.getline(temp, 100, '\t');
			in.getline(temp, 100);
			double c=atof(temp);

			// create grid region
			mgrid.add_sgr_tan(N_t, omega_l, omega_r, omega_c, c, id);
		}
		else 
		{
			cerr << "Error: read multigrid. sgridregion type '" << type << "' is wrong: . Break" << endl;
			throw 1;
		}
	}
	in.ignore(100, '\n');

	// create multigrid
	mgrid.create();
	// overwrite input grid
	igrid=mgrid;	
}
void readExact ( vector <double>& A, string  fileName, multigrid & mgrid)
{
	char temp[100];
	ifstream in;
	in.open(fileName.c_str());
	if (in.fail())
	{
		cerr << endl;
		cerr << "Error while trying to open file: " << fileName << endl;
		cerr << "Break" << endl;
		exit(1);
	}	
	in.close();
	in.open(fileName.c_str());
	int n_com=0;
	int n_data=0;
	while (!in.eof())
	{
		in.getline(temp,100);
		if (temp[0]=='#')
		{
			n_com++;
		}
		else
		{
			n_data++;
		}
	}
	if (n_data-1!=mgrid.M+1)
	{
		cerr << endl;
		cerr << "Number of data points is incorrect in file: " << fileName << endl;
		cerr << n_data-1 << " != " << mgrid.M+1 << endl;
		cerr << "Break" << endl;
		exit(1);
	}
	in.close();

	in.open(fileName.c_str());
	for (int i=0; i!=n_com; i++)
	{
		in.ignore(500,'\n');
	}
	for (int i=0; i!=mgrid.M+1; i++)
	{
		in.getline(temp, 100, '\t');
		in.getline(temp, 100);
		A[i]=atof(temp);
	}
	in.close();
}

// ************************************************************************
// ************* input/output for multigrids via **************************
// ************* the boost serialization library **************************
// ************************************************************************
void save (multigrid & mg, string folder, string title)
{
	string datName=folder + "/" + title + ".dat";
	ofstream out(datName.c_str());
	boost::archive::text_oarchive oa(out);
	oa << mg;
	out.close();
}
void save (multigrid & mg, string folder, string title, int counter)
{
	char counter_ch[10];
	string counter_str;
	sprintf(counter_ch, "%04i", counter);
	counter_str=counter_ch;	
	string datName=folder + "/" + title + "_" + counter_str + ".dat";
	ofstream out(datName.c_str());
	boost::archive::text_oarchive oa(out);
	oa << mg;
	out.close();
}
void read (multigrid & mg, string filename)
{
	ifstream in;
	in.open(filename.c_str());
	if (in.fail())
	{
		cerr << endl;
		cerr << "Error while trying to open file: " << filename << endl;
		cerr << "Break" << endl;
		exit(1);
	}	
	in.close();
	in.open(filename.c_str());
	boost::archive::text_iarchive ia(in);
	ia >> mg;
	in.close();
}
// ************************************************************************
// ************* input/output for vector of multigrids via ****************
// ************* the boost serialization library **************************
// ************************************************************************
void save (vector<multigrid> & mg, string folder, string title)
{
	string datName=folder + "/" + title + ".dat";
	ofstream out(datName.c_str());
	boost::archive::text_oarchive oa(out);
	oa << mg;
	out.close();
}
void save (vector<multigrid> & mg, string folder, string title, int counter)
{
	char counter_ch[10];
	string counter_str;
	sprintf(counter_ch, "%04i", counter);
	counter_str=counter_ch;	
	string datName=folder + "/" + title + "_" + counter_str + ".dat";
	ofstream out(datName.c_str());
	boost::archive::text_oarchive oa(out);
	oa << mg;
	out.close();
}
void read (vector<multigrid> & mg, string filename)
{
	ifstream in;
	in.open(filename.c_str());
	if (in.fail())
	{
		cerr << endl;
		cerr << "Error while trying to open file: " << filename << endl;
		cerr << "Break" << endl;
		exit(1);
	}	
	in.close();
	in.open(filename.c_str());
	boost::archive::text_iarchive ia(in);
	ia >> mg;
	in.close();
}
// ************************************************************************
// ************* input/output for generic vectors via *********************
// ************* the boost serialization library **************************
// ************************************************************************
void record (vector<double> & pp, string folder, string title)
{
	string datName=folder + "/" + title + ".dat";
	ofstream out(datName.c_str());
	boost::archive::text_oarchive oa(out);
	oa << pp;
	out.close();
}
void record (vector<double> & pp, string folder, string title, int counter)
{
	char counter_ch[10];
	string counter_str;
	sprintf(counter_ch, "%04i", counter);
	counter_str=counter_ch;	
	string datName=folder + "/" + title + "_" + counter_str + ".dat";
	ofstream out(datName.c_str());
	boost::archive::text_oarchive oa(out);
	oa << pp;
	out.close();
}
void load (vector<double> & pp, string filename)
{
	ifstream in;
	in.open(filename.c_str());
	if (in.fail())
	{
		cerr << endl;
		cerr << "Error while trying to open file: " << filename << endl;
		cerr << "Break" << endl;
		exit(1);
	}	
	in.close();
	in.open(filename.c_str());
	boost::archive::text_iarchive ia(in);
	ia >> pp;
	in.close();
}
void record (vector <vector<double> > & pp, string folder, string title)
{
	string datName=folder + "/" + title + ".dat";
	ofstream out(datName.c_str());
	boost::archive::text_oarchive oa(out);
	oa << pp;
	out.close();
}
void record (vector <vector<double> > & pp, string folder, string title, int counter)
{
	char counter_ch[10];
	string counter_str;
	sprintf(counter_ch, "%04i", counter);
	counter_str=counter_ch;	
	string datName=folder + "/" + title + "_" + counter_str + ".dat";
	ofstream out(datName.c_str());
	boost::archive::text_oarchive oa(out);
	oa << pp;
	out.close();
}
void load (vector <vector<double> > & pp, string filename)
{
	ifstream in;
	in.open(filename.c_str());
	if (in.fail())
	{
		cerr << endl;
		cerr << "Error while trying to open file: " << filename << endl;
		cerr << "Break" << endl;
		exit(1);
	}	
	in.close();
	in.open(filename.c_str());
	boost::archive::text_iarchive ia(in);
	ia >> pp;
	in.close();
}
// ************************************************************************
// ************* input/output for layer functions *************************
// ************************************************************************

// save layer and frequency dependent function in a folder for each layer (if it does not exist, create one)
// (each layer has a the same frequency grid)
void save (vector< vector <double> > & A, multigrid & dgrid, string folder, string title)
{
	string modFolder;
	string cmd;
	char n_ch[10];
	string n_str;
	if (A.size()>1)
	{
		for (int n=0; n!=A.size(); n++)
		{
			sprintf(n_ch, "%03i", n);
			n_str=n_ch;
			modFolder=folder + "layer_" + n_str; 
			cmd="mkdir -p " + modFolder;	
			system(cmd.c_str());
			save (A[n], dgrid, modFolder, title);
		}
	}
	else
	{
		save (A[0], dgrid, folder, title);
	}
}
void save (vector< vector <double> > & A, multigrid & dgrid, string folder, string title, int counter)
{
	string modFolder;
	string cmd;
	char n_ch[10];
	string n_str;
	if (A.size()>1)
	{
		for (int n=0; n!=A.size(); n++)
		{
			sprintf(n_ch, "%03i", n);
			n_str=n_ch;
			modFolder=folder + "layer_" + n_str; 
			cmd="mkdir -p " + modFolder;	
			system(cmd.c_str());
			save (A[n], dgrid, modFolder, title, counter);
		}
	}
	else
	{
		save (A[0], dgrid, folder, title, counter);
	}
}
// save layer and frequency dependent function in a folder for each layer (if it does not exist, create one)
// (each layer has a different frequency grid)
void save (vector< vector <double> > & A, vector <multigrid> & dgrid, string folder, string title)
{
	string modFolder;
	string cmd;
	char n_ch[10];
	string n_str;
	if (A.size()>1)
	{
		for (int n=0; n!=A.size(); n++)
		{
			sprintf(n_ch, "%03i", n);
			n_str=n_ch;
			modFolder=folder + "layer_" + n_str; 
			cmd="mkdir -p " + modFolder;	
			system(cmd.c_str());
			save (A[n], dgrid[n], modFolder, title);
		}
	}
	else
	{
		save (A[0], dgrid[0], folder, title);
	}
}
void save (vector< vector <double> > & A, vector <multigrid> & dgrid, string folder, string title, int counter)
{
	string modFolder;
	string cmd;
	char n_ch[10];
	string n_str;
	if (A.size()>1)
	{
		for (int n=0; n!=A.size(); n++)
		{
			sprintf(n_ch, "%03i", n);
			n_str=n_ch;
			modFolder=folder + "layer_" + n_str; 
			cmd="mkdir -p " + modFolder;	
			system(cmd.c_str());
			save (A[n], dgrid[n], modFolder, title, counter);
		}
	}
	else
	{
		save (A[0], dgrid[0], folder, title, counter);
	}
}

// save single values which are layer dependent in a existing file (if it does not exist, create one)
void save ( double value , int layer, string  folder ,string  title)
{
	// save value
	string datName=folder + "/" + title + ".dat";
	ofstream out;
	out.open(datName.c_str(), ios::app);
	out << scientific << setprecision(17) << layer << "\t" << value << endl;
	out.close();
}
void save ( double x, double y , int layer, string  folder ,string  title)
{
	// create layer folder if neccessary
	char layer_ch[10];
	string layer_str;
	sprintf(layer_ch, "%03i", layer);
	layer_str=layer_ch;
	string modFolder=folder + "layer_" + layer_str; 
	string cmd="mkdir -p " + modFolder;	
	system(cmd.c_str());

	// save value
	string datName=modFolder + "/" + title + ".dat";
	ofstream out;
	out.open(datName.c_str(), ios::app);
	out << scientific << setprecision(17) << x << "\t" << y << endl;
	out.close();
}
// save single values which are layer and iteration dependent in a existing file (if it does not exist, create one)
void save ( double value , int layer, string  folder ,string  title, int counter)
{
	// create layer folder if neccessary
	char layer_ch[10];
	string layer_str;
	sprintf(layer_ch, "%03i", layer);
	layer_str=layer_ch;
	string modFolder=folder + "layer_" + layer_str; 
	string cmd="mkdir -p " + modFolder;	
	system(cmd.c_str());

	// save value
	string datName=modFolder + "/" + title + ".dat";
	ofstream out;
	out.open(datName.c_str(), ios::app);
	out << scientific << setprecision(17) << counter << "\t" << value << endl;
	out.close();
}
// save layer dependent functions
void save (vector <double> & A, string folder, string title, int counter)
{
	char counter_ch[10];
	string counter_str;
	sprintf(counter_ch, "%04i", counter);
	counter_str=counter_ch;	
	string datName=folder + "/" + title + "_" + counter_str + ".dat";
	ofstream out;
	out.open(datName.c_str());
	for (int n=0; n!=A.size(); n++)
	{
		out << n << "\t" << scientific << setprecision(17) << A[n] << endl;
	}
	out.close();
}
void save (vector <double> & A, string folder, string title)
{
	string datName=folder + "/" + title + ".dat";
	ofstream out;
	out.open(datName.c_str());
	for (int n=0; n!=A.size(); n++)
	{
		out << n << "\t" << scientific << setprecision(17) << A[n] << endl;
	}
	out.close();
}
void save (vector < vector <double> > & A, string folder, string title, int counter)
{
	char counter_ch[10];
	string counter_str;
	sprintf(counter_ch, "%04i", counter);
	counter_str=counter_ch;	
	string datName=folder + "/" + title + "_" + counter_str + ".dat";
	ofstream out;
	out.open(datName.c_str());
	for (int n=0; n!=A.size(); n++)
	{
		for (int m=0; m!=A[n].size(); m++)
		{
			out << scientific << setprecision(17) << A[n][m] << "\t";
		}
		out << endl;
	}
	out.close();
}
void save (vector < vector <double> > & A, string folder, string title)
{
	string datName=folder + "/" + title + ".dat";
	ofstream out;
	out.open(datName.c_str());
	for (int n=0; n!=A.size(); n++)
	{
		for (int m=0; m!=A[n].size(); m++)
		{
			out << scientific << setprecision(17) << A[n][m] << "\t";
		}
		out << endl;
	}
	out.close();
}
// i/o for x-y-data
void saveXY (vector <vector <double> > & xydata, string folder, string title, int counter)
{
	if (xydata.size()!=2)
	{
		cerr << "Error in save x-y-data routine: Array has wrong size. Break." << endl;
		throw 1;
	}
	else if (xydata[0].size()!=xydata[1].size())
	{
		cerr << "Error in save x-y-data routine: Array has asymmetric size. Break." << endl;
		throw 1;
	}
	char counter_ch[10];
	string counter_str;
	sprintf(counter_ch, "%04i", counter);
	counter_str=counter_ch;	
	string datName=folder + "/" + title + "_" + counter_str + ".dat";
	ofstream out;
	out.open(datName.c_str());
	for (int n=0; n!=xydata[0].size(); n++)
	{
		out << scientific << setprecision(17) << xydata[0][n] << "\t" << xydata[1][n] << endl;
	}
	out.close();
}
void saveXY (vector <vector <double> > & xydata, string folder, string title)
{
	if (xydata.size()!=2)
	{
		cerr << "Error in save x-y-data routine: Array has wrong size. Break." << endl;
		throw 1;
	}
	else if (xydata[0].size()!=xydata[1].size())
	{
		cerr << "Error in save x-y-data routine: Array has asymmetric size. Break." << endl;
		throw 1;
	}
	string datName=folder + "/" + title + ".dat";
	ofstream out;
	out.open(datName.c_str());
	for (int n=0; n!=xydata[0].size(); n++)
	{
		out << scientific << setprecision(17) << xydata[0][n] << "\t" << xydata[1][n] << endl;
	}
	out.close();
}
