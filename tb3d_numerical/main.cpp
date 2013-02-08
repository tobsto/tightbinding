// ****************************************************************************
// Tobias Stollenwerk, stollenwerk@th.physik.uni-bonn.de
// ****************************************************************************

#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<fstream>
#include<vector>
#include<string>
#include<complex>
#include<limits>
#include"grid.h"
#include"multigrid.h"
#include"mesh.h"
#include"io.h"
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

void getDos3d_tb(vector<double> & dos3d_tb, multigrid & egrid, int N_k)
{
	// increment in k space
	double delta_k=2*M_PI/double(N_k);
	// hopping matrix element
	double t=1.0/6.0;

	// Calculate band structure
	double kx,ky,kz;
	vector < vector < vector <double> > > epsilon (N_k, vector <vector <double> > (N_k, vector <double> (N_k)));
	for (int i=0; i!=N_k; i++)
	{
		kx=2*M_PI*i/double(N_k)-M_PI;
		for (int j=0; j!=N_k; j++)
		{
			ky=2*M_PI*j/double(N_k)-M_PI;
			for (int l=0; l!=N_k; l++)
			{
				kz=2*M_PI*l/double(N_k)-M_PI;
				//if (kx+ky+kz<=M_PI)
				{
					epsilon[i][j][l]=1.0-2*t*(cos(kx)+cos(ky)+cos(kz));
				}
			}
		}
	}

	for (int n=1; n<=egrid.M-1; n++)
	{
		dos3d_tb[n]=0.0;
		cout << n << " of " << egrid.M << endl;
		for (int i=0; i!=N_k; i++)
		{
			for (int j=0; j!=N_k; j++)
			{
				for (int l=0; l!=N_k; l++)
				{
					if(epsilon[i][j][l]>=egrid.omega[n]-0.5*(egrid.omega[n]-egrid.omega[n-1]) && epsilon[i][j][l] < egrid.omega[n]+0.5*(egrid.omega[n+1]-egrid.omega[n]))
					{
						dos3d_tb[n]+=2.0/(egrid.omega[n+1]-egrid.omega[n-1]);
					}
				}
			}
		}
		// normalization factor
		dos3d_tb[n]*=1.0/double(N_k*N_k*N_k);
	}

	dos3d_tb[0]=0.0;
	dos3d_tb[egrid.M]=0.0;
	for (int i=0; i!=N_k; i++)
	{
		for (int j=0; j!=N_k; j++)
		{
			for (int l=0; l!=N_k; l++)
			{
				if(epsilon[i][j][l]>=egrid.omega[0] && epsilon[i][j][l] < egrid.omega[0]+0.5*(egrid.omega[0]-egrid.omega[1]))
				{
					dos3d_tb[0]+=2.0/(egrid.omega[1]-egrid.omega[0]);
				}

				if(epsilon[i][j][l]>=egrid.omega[egrid.M]-0.5*(egrid.omega[egrid.M]-egrid.omega[egrid.M-1]) && epsilon[i][j][l] <= egrid.omega[egrid.M])
				{
					dos3d_tb[egrid.M]+=2.0/(egrid.omega[egrid.M]-egrid.omega[egrid.M-1]);
				}
			}
		}
	}
	// normalization factor
	dos3d_tb[0]*=1.0/double(N_k*N_k*N_k);
	dos3d_tb[egrid.M]*=1.0/double(N_k*N_k*N_k);
}

double dn(double & epsilon, vector<double> & dos3d_tb, multigrid & egrid, double & n_cr)
{
	double result=0.0;
	for (int n=0; n<=egrid.M && egrid.omega[n]<=epsilon; n++)
	{
		result+=dos3d_tb[n]*egrid.domega[n];
	}
	return 2*result-n_cr;
}

class xTolerance
{
	public:
	double tol;
	xTolerance(double t)
	{
		this->tol=t;
	}
};

void findEpsilon(double & x, double & dntol, double & dmu, vector<double> & dos3d_tb, multigrid & egrid, double & n_cr) 
{
	double ftol=dntol;
	double xa, xb, xc, fa, fb, fc;
	double dx0=dmu;
	double dx=dx0;
	double dx1=0.1;
	// precision
	double xtol=numeric_limits<double>::epsilon()*max(egrid.omega_min, egrid.omega_max);
	int max1=100;
	int max2=100;
	double lowerBound=egrid.omega_min;
	double upperBound=egrid.omega_max;
	int muSearchCount=0;

	xa=x;
	fa=dn(xa, dos3d_tb, egrid, n_cr);
	cout << "\t* mu:\t" << xa << "\t" << "delta_n:\t" << fa << endl;
	muSearchCount++;

	if (fabs(fa)>ftol)
	{
		xb=xa - fa/abs(fa) * dx;
		fb=dn(xb, dos3d_tb, egrid, n_cr);
		cout << "\t* mu:\t" << xb << "\t" << "delta_n:\t" << fb << endl;
		muSearchCount++;
		if (fabs(fb)<=ftol)
		{
			x=xb;
			return;
		}	

		int counter=0;
		bool constFuncFlag=false;
		double factor=1.5;
		while (fa*fb>0.0 && counter<max2 && fabs(fb)>ftol)
		{
			if (fabs(fb-fa)>xtol)
			{
				xc=xb - fb * (xb-xa)/(fb-fa);
				constFuncFlag=false;
			}
			else
			{
				constFuncFlag=true;
			}
			if (constFuncFlag || xc>=upperBound || xc<=lowerBound)
			{
				xc=xb - fb/abs(fb) * dx1;
				while (xc>=upperBound)
				{
					xc-=dx1;
				}
				while (xc<=lowerBound)
				{
					xc+=dx1;
				}
			}
			
			fc=dn(xc, dos3d_tb, egrid, n_cr);
			cout << "\t* mu:\t" << xc << "\t" << "delta_n:\t" << fc << endl;
			muSearchCount++;
	
			xa=xb;
			fa=fb;
			xb=xc;
			fb=fc;
			counter++;
		}
		if (counter==max2)
		{
			// Did not found two points with different signs
			if (fa<0.0)
			{
				x=1.0;
				return;
			}
			else
			{
				x=-1.0;
				return;
			}
		}
		if (fabs(fb)<=ftol)
		{
			x=xb;
			return;
		}
			
		// inital values for bounded secant (bisection) method
		double x1, x2, f1, f2, f;
		if (fa<fb)
		{
			x1=xa;
			f1=fa;
			x2=xb;
			f2=fb;
		}
		else
		{
			x1=xb;
			f1=fb;
			x2=xa;
			f2=fa;
		}

		// container for difference in function values above and below the x-axis 
		vector<double> df;
		// container storing if element of df is below or above x-axis
		vector<bool> below;
		double avdfp;
		double avdfm;
		double fp_before=f2;
		double fm_before=f1;
		bool contFlag=true;
		int N_av=10;
		double df_tol=1E-6;

		// bisection method to find root
		counter=0;
		do
		{
			if (fabs(f2-f1)>xtol)
			{
				// secant method
				x=x1-f1*(x2-x1)/(f2-f1);
			}
			// if there is a division by 0 due to equal function values f1 and f2, use bisection
			else
			{
				x=0.5*(x1+x2);
			}
			// get function value
			f=dn(x, dos3d_tb, egrid, n_cr);
			cout << "\t* mu:\t" << x << "\t" << "delta_n:\t" << f << endl;
			muSearchCount++;


			// Stop iteration for the case of switching between two function values above and below the x-axis
			// this happens, since the numerically evaluated function can not be calculated with arbitrary 
			// accuracy. The algorithm stops if the difference of the function values below and above the x-axis
			// have only tiny variations (<1E-6) over some iterations (10).
			// Therefore one has to keep track of the difference in the function values below (df[0]) and above (df[1])
			// the x-axis and calculate the average over the last iterations (avdfm and avdfp).
			if (f<=0)
			{
				below.push_back(true);
				df.push_back(fabs(f-fm_before));
				fm_before=f;
			}
			else
			{
				df.push_back(fabs(f-fp_before));
				below.push_back(false);
				fp_before=f;
			}
	

			if (counter>=N_av)
			{
				// calculate average of difference in function values over the
				// N_av iterations. Above and below the x-axis. If both values
				// are below df_tol, stop iteration and set root to the last 
				// x value of that branch (above or below) which occured more
				// often in the last N_av iterations. 
				avdfm=0.0;
				avdfp=0.0;
				int N_below=0;
				for (int i=df.size()-N_av; i<df.size(); i++)
				{
					if (below[i])
					{
						avdfm+=df[i];
						N_below++;
					}
					else
					{
						avdfp+=df[i];
					}
				}
				avdfm/=max(N_below,1);
				avdfp/=max(N_av-N_below,1);

				if (avdfm<=df_tol && avdfp<=df_tol)
				{
					contFlag=false;	
				}
			}
			if (contFlag==false)
			{
				//cerr << endl;
				//cerr << "Warning: Error in fixMu routine. Search for root was not successful." << endl;
				//cerr << "Difference in x points fell under machine precision." << endl;
				//throw xTolerance(f);
	double delta_f;
	delta_f=fabs(fp_before-fm_before);
	throw xTolerance(delta_f);
			}
			if (f*f1<0.0)
			{
				x2=x;
				f2=f;
			}
			else
			{
				x1=x;
				f1=f;
			}
			counter++;	
		}		
		while (fabs(f)>ftol && counter<max2 && fabs(x2-x1)>=xtol);

		if (counter==max2)
		{
			cerr << endl;
			cerr << "Warning: Error in fixMu routine. Search for root was not successful."<< endl;
			cerr << "Not successful after " << counter << " iterations. Tolerance " << ftol << " to small?" << endl;
			throw xTolerance(f);
		}
		if (fabs(x2-x1)<xtol)
		{
			//cerr << endl;
			//cerr << "Warning: Error in fixMu routine. Search for root was not successful." << endl;
			//cerr << "Difference in x points fell under machine precision." << endl;
			throw xTolerance(f);
		}
	}
}

int main(int argc, char * argv[])
{
	int N=100;
	string input;
	string output;
	bool inputExists;
	double n_cr;
	try 
	{
		
		// Declare a group of options that will be 
		// allowed only on command line
		po::options_description generic("Generic options");
		generic.add_options()
			("help,h", "produce help message")
		;
		
		// Declare a group of options that will be 
		// allowed both on command line and in
		// config file
		po::options_description options("Standard options");
		options.add_options()
			("N,n"               , po::value<int>(          &N)->default_value(50)          , "Number of sample points in k-space in each direction")
			("filling,c"         , po::value<double>(    &n_cr)->default_value(0.01)        , "Band filling (from 0.0 to 2.0)")
			("output,o"          , po::value<string>( &output)->default_value("output/")   , "Output folder")
			("input,i"           , po::value<string>()                                     , "input folder")
		;

		po::options_description all("Allowed options");
		all.add(generic).add(options);
		
		po::options_description cmdline_options;
		cmdline_options.add(generic).add(options);
		
		po::options_description visible("Allowed options");
		visible.add(generic).add(options);
		
		po::variables_map vm;
		store(parse_command_line(argc, argv, all), vm);
		notify(vm);
		
		inputExists=(!vm["input"].empty()); 

		if (inputExists)
		{
			input=vm["input"].as<string>();
		}
		// help 
		if (vm.count("help")) 
		{
			cout << visible << "\n";
			exit(0);
		}
	}
	catch(exception& e)
	{
		cerr << "Something went wrong during the parameter input. Exception:" << endl;
		cerr << e.what() << endl;
		cerr << "Break." << endl;
		exit(1);
	}

	multigrid egrid;
	vector <double> dos;
	if (inputExists)
	{
		cout << input << endl;
		string pathToGrid=input + "/grid.dat";
		read(egrid, pathToGrid);
		dos=vector <double> (egrid.M+1);
		string pathToDoS=input + "/dos.dat";
		readExact(dos, pathToDoS, egrid);
	}
	else
	{	
		egrid.add_gr_equi(0.0, 2.0, 0.001);
		//egrid.add_sgr_log(-1.0, 0.1, 0.001, 0.01);
		//egrid.add_sgr_log( 1.0, 0.1, 0.001, 0.01);
		//egrid.add_sgr_log(-1.0/3.0, 0.1, 0.001, 0.01);
		//egrid.add_sgr_log( 1.0/3.0, 0.1, 0.001, 0.01);
		//egrid.add_lendpoint(-1.0);
		//egrid.add_rendpoint( 1.0);
		egrid.create();
		dos=vector <double> (egrid.M+1);
	
		getDos3d_tb(dos, egrid, N);
	}
	string cmd="mkdir -p " + output;
	system(cmd.c_str());
	save(dos, egrid, output, "dos");
	save(egrid, output, "grid");

	double epsilon=0.0;
	double dntol=1E-15;
	double depsilon=0.1;
	try
	{
		findEpsilon(epsilon, dntol, depsilon, dos, egrid, n_cr);
	}
	catch(xTolerance & xtol)
	{
		cerr << "Warning: Get Epsilon: Accuracy was only " << fabs(xtol.tol) ;
		cerr << ", relative: " << fabs(xtol.tol)/fabs(epsilon) <<  endl;
		ofstream out("accuracy.dat");
		out << setprecision(15) << scientific << N << "\t" << xtol.tol << endl;
		out.close();
	}
	catch(...)
	{
		cerr << "Error: getDelta_r: Something went wrong. Break." << endl;
		exit(1);
	}
	string occOutput= output + "occNum.dat";
	ofstream out(occOutput.c_str());
	for (double x=0; x<=2.0; x+=1E-3)
	{
		out << setprecision(15) << scientific << x << "\t" << dn(x, dos, egrid, n_cr) << endl;
	}
	out.close();
	cout << "Epsilon= " << scientific << setprecision(17) << epsilon << endl;
	double kf=acos((1.0-epsilon)/(2.0/6.0)-2);
	cout << "Momentum kf (1/a)= " << scientific << setprecision(17) << kf << endl;
	double lambda=2.0*M_PI/(2.0*kf);
	cout << "Friedel wavelength lambda=2*pi/(2*kf) (a)= " << lambda << endl;
}
