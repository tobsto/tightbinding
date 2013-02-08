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
#include<fstream>
#include<vector>
#include<string>
#include<complex>
#include <boost/archive/text_oarchive.hpp> 
#include <boost/archive/text_iarchive.hpp> 

using namespace std;

// ************************************************************************
// *************  create numbered filenames to avoid overwriting  *********
// ************************************************************************
string appFilename(string  folder, string  title);
// ************************************************************************
// *************  input/output for generic grids  *************************
// ************************************************************************
void save ( double x     , string folder, string title);
void save ( double value , string folder ,string title, int counter);
void save (double x, double y,  string  folder, string  title);
void save ( vector <double>& A, vector <double>& omega, int N, string  folder,  string  title);
void save ( vector <double>& A, grid & dgrid, string  folder, string  title);
void save (vector< vector <double> > & A, grid & dgrid, string folder, string title);
void savePeakPara ( vector < vector <double> > & peakPara, string folder, int counter);
void savePseudoPeakPara ( vector <double> & peakPara, string  folder, int counter);
void read (double & x, string fileName);
void read( vector <double>& A, string fileName);
void read( vector <double>& A, string fileName, grid & mgrid, grid & igrid);
void readInter( vector <double>& A, string  fileName, grid & mgrid);

// ************************************************************************
// ************* input/output for multigrids ******************************
// ************************************************************************
void save ( vector <double>& A, const multigrid & mgrid, string  folder, string  title);
void save ( vector <double>& A, const multigrid & mgrid, string  folder, string  title, int counter);
void readGrid(string filename, multigrid & igrid);
void readExact ( vector <double>& A, string  fileName, multigrid & mgrid);

// ************************************************************************
// ************* input/output for layer functions *************************
// ************************************************************************
void save (vector< vector <double> > & A, multigrid & dgrid, string folder, string title);
void save (vector< vector <double> > & A, multigrid & dgrid, string folder, string title, int counter);
void save (vector< vector <double> > & A, vector <multigrid> & dgrid, string folder, string title);
void save (vector< vector <double> > & A, vector <multigrid> & dgrid, string folder, string title, int counter);
void save ( double value , int layer, string  folder ,string  title);
void save ( double value , int layer, string  folder ,string  title, int counter);
void save ( double x, double y , int layer, string  folder ,string  title);
void save (vector <double> & A, string folder, string title, int counter);
void save (vector <double> & A, string folder, string title);
void save (vector < vector <double> > & A, string folder, string title, int counter);
void save (vector < vector <double> > & A, string folder, string title);

// i/o for x-y-data
void saveXY(vector <vector <double> > & xydata, string folder, string title, int counter);
void saveXY (vector <vector <double> > & xydata, string folder, string title);


// ************************************************************************
// ************* input/output for multigrids via **************************
// ************* the boost serialization library **************************
// ************************************************************************
void save (multigrid & mg, string folder, string title);
void save (multigrid & mg, string folder, string title, int counter);
void read (multigrid & mg, string filename);
void save (vector<multigrid> & mg, string folder, string title);
void save (vector<multigrid> & mg, string folder, string title, int counter);
void read (vector<multigrid> & mg, string filename);
namespace boost
{
	namespace serialization
	{
		template<class Archive>
		void serialize(Archive & ar, grid & g, const unsigned int version)
		{
			ar & g.M;
			ar & g.omega_min;
			ar & g.omega_max;
			ar & g.omega;
			ar & g.domega;
		}
		template<class Archive>
		void serialize(Archive & ar, gridRegion & gr, const unsigned int version)
		{
			ar & gr.omega_l;
			ar & gr.omega_c;
			ar & gr.omega_r;
			ar & gr.omega_m;
			ar & gr.omega_p;
			ar & gr.domega_min_l;
			ar & gr.domega_min_r;
			ar & gr.N_d;
			ar & gr.N_i;
			ar & gr.para_d;
			ar & gr.para_i;
			ar & gr.type;
			ar & gr.id;
		}
		template<class Archive>
		void serialize(Archive & ar, subgrid & sg, const unsigned int version)
		{
			ar & sg.omega_l;
			ar & sg.omega_r;
			ar & sg.i_l;
			ar & sg.i_r;
			ar & sg.N_replaced;
			ar & sg.N_d;
			ar & sg.N_i;
			ar & sg.para_d;
			ar & sg.para_i;
			ar & sg.type;
		}	
		template<class Archive>
		void serialize(Archive & ar, multigrid & mg, const unsigned int version)
		{
			ar & boost::serialization::base_object<grid>(mg);
			ar & mg.fsubgrids;
			ar & mg.ssubgrids;
			ar & mg.fgridRegions;
			ar & mg.sgridRegions;
		}
	}
}
// ************************************************************************
// ************* input/output for generic vectors via *********************
// ************* the boost serialization library **************************
// ************************************************************************
void record (vector<double> & pp, string folder, string title);
void record (vector<double> & pp, string folder, string title, int counter);
void load (vector<double> & pp, string filename);
void record (vector <vector<double> > & pp, string folder, string title);
void record (vector <vector<double> > & pp, string folder, string title, int counter);
void load (vector <vector<double> > & pp, string filename);


// ************************************************************************
// ************* columns-class ********************************************
// ************************************************************************
class columns
{
	public:
	int Nr;
	int Nc;
	vector <int> indices;
	vector < vector <double> >* matrix;

	columns(vector < vector <double> > & A, vector<int> & ind)
	{
		this->matrix=&A;
		this->Nr=A.size();
		this->Nc=A[0].size();
		this->indices=ind;
	}
	columns(vector < vector <double> > & A)
	{
		this->matrix=&A;
		this->Nr=A.size();
		this->Nc=A[0].size();
	}
	columns(){};
	void fill (vector < vector <double> > & A, vector<int> & ind)
	{
		this->matrix=&A;
		this->Nr=A.size();
		this->Nc=A[0].size();
		this->indices=ind;
	}
};

// ************************************************************************
// ************* integer columns-class ************************************
// ************************************************************************
class icolumns
{
	public:
	int Nr;
	int Nc;
	vector <int> indices;
	vector < vector <int> >* matrix;

	icolumns(vector < vector <int> > & A, vector<int> & ind)
	{
		this->matrix=&A;
		this->Nr=A.size();
		this->Nc=A[0].size();
		this->indices=ind;
	}
	icolumns(vector < vector <int> > & A)
	{
		this->matrix=&A;
		this->Nr=A.size();
		this->Nc=A[0].size();
	}
	icolumns(){};
	void fill (vector < vector <int> > & A, vector<int> & ind)
	{
		this->matrix=&A;
		this->Nr=A.size();
		this->Nc=A[0].size();
		this->indices=ind;
	}
};

// ************************************************************************
// ************* complex columns-class ************************************
// ************************************************************************
class ccolumns
{
	public:
	int Nr;
	int Nc;
	vector <int> indices;
	vector < vector <complex<double> > >* matrix;

	ccolumns(vector < vector <complex<double> > > & A, vector<int> & ind)
	{
		this->matrix=&A;
		this->Nr=A.size();
		this->Nc=A[0].size();
		this->indices=ind;
	}
	ccolumns(vector < vector <complex<double> > > & A)
	{
		this->matrix=&A;
		this->Nr=A.size();
		this->Nc=A[0].size();
	}
	ccolumns(){};
	void fill (vector < vector <complex<double> > > & A, vector<int> & ind)
	{
		this->matrix=&A;
		this->Nr=A.size();
		this->Nc=A[0].size();
		this->indices=ind;
	}
};

// ************************************************************************
// ************* columns-of-variable-length - class ***********************
// ************************************************************************
class varcolumns
{
	public:
	int Nr;
	int Nc;
	vector < vector <int> > indices;
	vector < vector <double> >* matrix;

	varcolumns(vector < vector <double> > & A, vector < vector<int> > & ind)
	{
		this->matrix=&A;
		this->Nr=A.size();
		this->Nc=A[0].size();
		this->indices=ind;
		/*
		this->indices=vector< vector <int> > (ind.size());
		for (int i=0; i!=ind.size(); i++)
		{
			for (int j=0; j!=ind[i].size(); j++)
			{
				this->indices[i].push_back(ind[i][j]);
			}
		}
		*/
	}
	varcolumns(vector < vector <double> > & A)
	{
		this->matrix=&A;
		this->Nr=A.size();
		this->Nc=A[0].size();
	}
	varcolumns(){};
	void fill (vector < vector <double> > & A, vector < vector<int> >& ind)
	{
		this->matrix=&A;
		this->Nr=A.size();
		this->Nc=A[0].size();
		this->indices=ind;
		/*
		if (Nr!=ind.size())
		{
			this->indices.resize(ind.size());
		}
		for (int i=0; i!=ind.size(); i++)
		{
			this->indices[i].erase(this->indices[i].begin(), this->indices[i].end());
			for (int j=0; j!=ind[i].size(); j++)
			{
				this->indices[i].push_back(ind[i][j]);
			}
		}
		*/
	}
};

// ************************************************************************
// ************* input/output for columns-class via ***********************
// ************* the boost serialization library **************************
// ************************************************************************
namespace boost {
	namespace serialization {

		template<class Archive>
		void serialize(Archive & ar, columns & g, const unsigned int version)
		{
			ar & g.Nr;
			ar & g.Nc;
			ar & g.indices;
		
			int j;
			for (int i=0; i!=g.Nr; i++)
			{		
				for (int jj=0; jj!=g.indices.size(); jj++)
				{		
					j=g.indices[jj];
					ar & (*g.matrix)[i][j];
				}
			}
		}
		template<class Archive>
		void serialize(Archive & ar, ccolumns & g, const unsigned int version)
		{
			ar & g.Nr;
			ar & g.Nc;
			ar & g.indices;
		
			int j;
			for (int i=0; i!=g.Nr; i++)
			{		
				for (int jj=0; jj!=g.indices.size(); jj++)
				{		
					j=g.indices[jj];
					ar & (*g.matrix)[i][j];
				}
			}
		}
		template<class Archive>
		void serialize(Archive & ar, icolumns & g, const unsigned int version)
		{
			ar & g.Nr;
			ar & g.Nc;
			ar & g.indices;
		
			int j;
			for (int i=0; i!=g.Nr; i++)
			{		
				for (int jj=0; jj!=g.indices.size(); jj++)
				{		
					j=g.indices[jj];
					ar & (*g.matrix)[i][j];
				}
			}
		}
		template<class Archive>
		void serialize(Archive & ar, varcolumns & g, const unsigned int version)
		{
			ar & g.Nr;
			ar & g.Nc;
			ar & g.indices;
		
			int j;
			for (int i=0; i!=g.Nr; i++)
			{		
				for (int jj=0; jj!=g.indices[i].size(); jj++)
				{		
					j=g.indices[i][jj];
					ar & (*g.matrix)[i][j];
				}
			}
		}
	}
}

// ******************************************************
// ****** output class for testing purposes *************
// ******************************************************
template<typename T>
class output 
{
	public:
	string filename;
	ofstream out;
	vector<ofstream*> vout;
	// output x-y values with precision 'prec' to 
	// serveral files which filenames differ by 
	// appendix number form 0 to 'N'
	output(string fname, int prec, int N)
	{
		vout.resize(N);
		char n_ch[10];	
		string n_str;
		for (int n=0; n!=N; n++)
		{
			sprintf(n_ch, "%04i", n);
			n_str=n_ch;
			string filename=fname + "_" + n_str + ".dat";
			vout[n]=new ofstream;
			(*this->vout[n]).open(filename.c_str());
			(*this->vout[n]).precision(prec);
		}
	}
	// output x-y values with precision 'prec' to 
	// serveral files which filenames differ by 
	// appendix number form 0 to 'N'. Also 
	// read part of the filename from file named 'codefilename' 
	// and choose if the file is overwritten or written to
	output(string fname, string codefilename, int prec, int N)
	{
		ifstream in(codefilename.c_str());
		if (in.fail())
		{
			cerr << codefilename << " does not exist. Break" << endl;
		}
		string code;
		in >> code;
		in.close();

		vout.resize(N);
		char n_ch[10];	
		string n_str;
		for (int n=0; n!=N; n++)
		{
			sprintf(n_ch, "%04i", n);
			n_str=n_ch;
			string filename=fname + code + "_" + n_str + ".dat";
			vout[n]=new ofstream;
			(*this->vout[n]).open(filename.c_str());
			(*this->vout[n]).precision(prec);
		}
	}
	// output x-y values with precision 'prec' and 
	// append number (up to 1E4) to the filename 
	// if necessary
	output(string fname, int prec)
	{
		char n_ch[10];	
		string n_str;
		for (int n=0; n!=1E4; n++)
		{
			sprintf(n_ch, "%04i", n);
			n_str=n_ch;
			string filename0=fname + "_" + n_str + ".dat";
			ifstream in(filename0.c_str());
			if (in.fail())
			{
				this->filename=filename0;
				break;
			}
		}
		this->out.open(filename.c_str());
		out.precision(prec);
	}
	// output x-y values with precision 'prec' and 
	// choose if the file is overwritten or written to
	output(string fname, int prec, ios_base::openmode mode)
	{
		this->filename=fname + ".dat";
		this->out.open(filename.c_str(), mode);
		out.precision(prec);
	}
	// output x-y values with precision 'prec' and 
	// read part of the filename from file named 'codefilename' 
	// and choose if the file is overwritten or written to
	output(string fname, string codefilename, int prec, ios_base::openmode mode)
	{
		ifstream in(codefilename.c_str());
		if (in.fail())
		{
			cerr << codefilename << " does not exist. Break" << endl;
		}
		string code;
		in >> code;
		in.close();
		this->filename=fname + code + ".dat";
		this->out.open(filename.c_str(), mode);
		out.precision(prec);
	}
	// destructor (close filestream)
	~output()
	{
		this->out.close();
		for (int n=0; n!=this->vout.size(); n++)
		{
			(*this->vout[n]).close();
		}
	}
	// write x-y values
	void write(T & x, T & y)
	{
		this->out << scientific << x << "\t" << y << endl;
	}
	// write x-y values to n-th stream
	void write(T & x, T & y, int n)
	{
		if (n<vout.size())
		{
			(*this->vout[n]) << scientific << x << "\t" << y << endl;
		}
		else
		{
			cerr << "Error: output class: index out of range. break." << endl;
			exit(1);
		}
	}
};
