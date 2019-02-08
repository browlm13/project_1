/*

					Shallow Water

	L.J. Brown
	Parallel Scientific Computing 001C 1192 

*/

///////////////////////////////////////////////////////
/////                
//      Imports                    

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <chrono>

using namespace std;
using namespace std::chrono;

///////////////////////////////////////////////////////
/////            
//  	Usefull Constants 

const double PI = 3.145926535;

///////////////////////////////////////////////////////
/////               
//      Settings         

/////
// File Names and Settings

// [TODO]: Fix this F I/O stuff
const string LOG_FNAME = "../output_3/shallow_water_log.txt";
const string X_FNAME = "../output_3/xs.csv";
const string Y_FNAME = "../output_3/ys.csv";
const string TIMESTEP_FPREFIX = "../output_3/timesteps/";
const string TIMESTEP_FEXTENSION = ".csv";

const int N_FRAMES = 200;
const bool SILENCE_FRAMES = 1;

////
// Measure Convergence Rate Boolean If True Solution Provided
const bool MEASURE_CONVERGENCE = 0;

////
// CFL Condition
const double CFL = 0.01;

/////
// Mesh Settings

// Mesh Ranges -- x,y,t

const double MIN_X = 0.0;
const double MAX_X = PI;

const double MIN_Y = 0.0;
const double MAX_Y = PI/2;

const double MIN_T = 0.0;
const double MAX_T = 10.0;

// Mesh Size -- x,y

const int NX = 100;
const int NY = 100;

/////
// Forcing Function Settings

double fh(double t, double x, double y){ 	// <cmath> for exp, pow
	//return 0.0;
	return exp(-15.0*pow((t-1.5),2));
}

double fu(double t, double x, double y){
	return 0.0;
}

double fv(double t, double x, double y){
	return 0.0;
}

/////
// Wave Speed

double c(double x, double y){				// <cmath> for cos, pow
	//return 1.0;
	return 0.75 + pow(cos(2*PI*x), 2)*pow(cos(PI*y), 2); //max speed (1.75)
}

/////
// Lambda

double lambda(double delta, double x, double y){
	return CFL*delta/c(x,y); 
}

/////
// Lambda * c^2

double lambda_c2(double delta, double x, double y){
	return CFL*delta*c(x,y);
}

///////////////////////////////////////////////////////
/////
// Helper Functions

void update_lrud(int x_iter,int y_iter, int (& lrud_ref) [4]){

	lrud_ref[0] = x_iter-1;
	lrud_ref[1] = x_iter+1;
	lrud_ref[2] = y_iter+1;
	lrud_ref[3] = y_iter-1;

	if (x_iter == 0){
        lrud_ref[0] = NX-1;
    }
    if (x_iter == NX-1){
        lrud_ref[1] = 0;
    }
    if (y_iter == 0){
        lrud_ref[2] = NY-1;
    }
    if (y_iter == NY-1){
        lrud_ref[3] = 0;
    }
}

// [TODO]: Fix this code alot
bool write_log(string fname){

	ofstream filelog (fname);

	if (filelog.is_open()){

	    filelog << "nx: " << NX << ",\n";
	    filelog << "ny: " << NY << ",\n";

	    filelog << "x_range: [" << MIN_X << ", " << MAX_X << "],\n";
	    filelog << "y_range: [" << MIN_Y << ", " << MAX_Y << "],\n";
	    filelog << "t_range: [" << MIN_T << ", " << MAX_T << "],\n";

	    filelog << "CFL: " << CFL << ",\n";

	    filelog.close();

	    cout << "Completed Writing Program Log to file: " << LOG_FNAME << endl;

	    return 0;
	 }
	 else cout << "Error Opening File: " << fname << endl;
	 return 1;
}

// [TODO]: Fix this code alot
bool write_xy_data(double (& xs_ref) [NX], double (& ys_ref) [NY]){

	// Write xs 
	ofstream xs_file (X_FNAME);

	if (xs_file.is_open()) {
	    for(int i = 0; i < NX; i ++) {
	        xs_file << xs_ref[i] << ", " ;
	    }
	    xs_file.close();
	    cout << "Completed Writing X values to file: " << X_FNAME << endl;
	}
	else {
	 	cout << "Problem Writing X values to file.";
	 	return 1;
	}

	// Write ys 
	ofstream ys_file (Y_FNAME);

	if (ys_file.is_open()) {
	    for(int j = 0; j < NY; j ++) {
	        ys_file << ys_ref[j] << ", " ;
	    }
	    ys_file.close();
	    cout << "Completed Writing Y values to file: " << Y_FNAME << endl;
	}
	else {
	 	cout << "Problem Writing Y values to file.";
	 	return 1;
	}

	return 0;

}

// [TODO]: Fix this code alot
bool write_timestep(int frame_iterator, double t, double (& H_ref) [NX][NY]){

	// Write xs 
	//
	string timestep_fname = TIMESTEP_FPREFIX + to_string(frame_iterator) + TIMESTEP_FEXTENSION;
	
	// console output
	if (!SILENCE_FRAMES){
		cout << timestep_fname << endl;
	}

	ofstream Ht_file (timestep_fname);

	// format -- csv 
	// h00,  ..., h0ny,  t
	// ...
	// hnx0, ..., hnxny, t

	if (Ht_file.is_open()) {

		// iterate over H matrix (transpose)
		for (int j=0; j < NY; j+=1 ){
			for (int i=0; i < NX; i+=1 ){
				Ht_file << H_ref[i][j] << ", " ;
			}

			// add time as final element of row
			Ht_file << t << "\n";
		}
		Ht_file.close();

		if (!SILENCE_FRAMES){
			cout << "Wrote time step: " << frame_iterator << ", t = " << t << ". " << endl;
		}
		return 0;
	}	
	else cout << "Problem Writing Time Step to file.";

	return 1;

}

///////////////////////////////////////////////////////
/////                 					
//      Known Solution for Testing       
//        				
//
////  	Parameters:
//
//////				Wave Speed:         
//						c = 1.0,  "Constant"
//////				Forcing Functions:  
//						fh = fu = fv = 0.0, "None"


/////
// 		Wave Velocity Component Frequencies

const double KX = 1.0;
const double KY = 2.0;
const double KT = sqrt(pow(KX,2) + pow(KY,2)); 		// <cmath> for sqrt, pow

/////
// Wave Height Function

double h(double t, double x, double y){ 			// <cmath> for cos
	// function of (t,x,y)
	//return cos(2*PI*KT*c(x,y) * t) * sin(2*PI*KX * x) * cos(2*PI*KY * y);
	return 0.0;
}

/////
// Wave Velocity Component Functions

double u(double t, double x, double y){				// <cmath> for sin, cos, pow
	// x component of velocity, function of (t,x,y)
	//return pow(c(x,y),2)*(KX/KT) * cos(2*PI*KT * t) * sin(2*PI*KX * x) * cos(2*PI*KY * y);
	return 0.0;
}

// y component of velocity, function of (t,x,y)
double v(double t, double x, double y){				// <cmath> for sin, cos, pow
	//return pow(c(x,y),2)*(KX/KT) * cos(2*PI*KT * t) * cos(2*PI*KX * x) * sin(2*PI*KY * y);
	return 0.0;
}

/////
// Fill True Solution Function

void shallow_water_solution(double t, double (& H_ref) [NX][NY], double (& U_ref) [NX][NY], double (& V_ref) [NX][NY], double (& xs_ref) [NX], double (& ys_ref) [NY]){
	// Fill 2D Arrays With True Solutions For Given Time
	double x = MIN_X;
	double y = MIN_Y;
	for ( int i=0; i<NX; i+=1 ) { // x iterator
		x = xs_ref[i];
		for ( int j=0; j<NY; j+=1 ) { // y iterator
			y = ys_ref[j];
	   		
	   		H_ref[i][j] = h(t,x,y);
	   		U_ref[i][j] = u(t,x,y);
	   		V_ref[i][j] = v(t,x,y);

		}
	}
}

///////////////////////////////////////////////////////
////                 								
//      			Run	Program				         					  
////
///////////////////////////////////////////////////////

int main(){

	// log/cout program start
	cout << "Shallow Water" << endl;

	/* Save basic input information */
	write_log(LOG_FNAME);

	////
	// Initilize x and y arrays, evenly spaced, and determine dx,dy and delta (the minimum)

	double xs [NX] = {0.0};
	double ys [NY] = {0.0};

	double dx = (MAX_X - MIN_X)/NX;
	double dy = (MAX_Y - MIN_Y)/NY;
	double delta = min(dx,dy);

	// fill xs and ys
	double current_x = MIN_X;
	for (int i = 0; i < NX; i += 1){
		xs[i] = current_x;
		current_x += dx;
	}
	double current_y = MIN_Y;
	for (int i = 0; i < NY; i += 1){
		ys[i] = current_y;
		current_y += dy;
	}

	/////
	// Write xs and ys to file
	write_xy_data(xs, ys);

	//find dt and nt
	double dt = (lambda(delta, MIN_X, MIN_Y) + lambda(delta, MAX_X, MAX_Y))/2; //[TODO]: Should be max c(mesh) 
	int nt = (MAX_T - MIN_T)/dt/2;

	cout << "Total number of time steps to run: " << nt << endl;

	////
	// Initilize Matrcies - Initilize Elements To Zero

	double H_old [NX][NY] = {0.0};
	double U_old [NX][NY] = {0.0};
	double V_old [NX][NY] = {0.0};

	double H_new [NX][NY] = {0.0};
	double U_new [NX][NY] = {0.0};
	double V_new [NX][NY] = {0.0};

	double H_swap [NX][NY] = {0.0};
	double U_swap [NX][NY] = {0.0};
	double V_swap [NX][NY] = {0.0};
	
	// [TODO]: Flatten Matrices

	////
	// Initilize t iterator to 2 (Given 2 Time steps
	int t_iter = 2;

	////
	// Initilize Frame Counter to 0
	int frame_iter = 0;

	////
	// For Testing Set Initial 2 Timesteps To True Solutions
	double t = MIN_T;
	shallow_water_solution(t, H_old, U_old, V_old, xs, ys);

	// Write Intital Values to File
	write_timestep(0, t, H_old);

	// Increment Frame Counter
	frame_iter += 1;

	// increment t and take second timestep
	t += dt;
	shallow_water_solution(t, H_old, U_old, V_old, xs, ys);

	// initilize left right up down array
	int lrud [4] = {-1,-1,-1,-1};

	////
	// Begin Timer
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	/////
	// If Measuring Convergence Rate Store Previous error
	if (MEASURE_CONVERGENCE){
		double prev_error = -1;
		double cur_error = -1;
	}


	////
	// Run
	int completion = 0;
	for (t_iter = t_iter; t_iter < nt*nt; t_iter+=1){
		t += dt; 

		// [TODO]: Make forcing funcintions a function of i and j -- precompute x and y values. -- and add to U and V
		
		////
		// leap frog
		double x = MIN_X;
		double y = MIN_Y;
		for (int i = 0; i < NX; i += 1){ // x iterator
			x += dx;
			for (int j = 0; j < NY; j += 1){ // y iterator
				y += dy; 

				// update left right up down array
				update_lrud(i,j,lrud);

	            H_new[i][j] = H_old[i][j] -lambda( delta, i, j ) * ( U_new[lrud[1]][j] - U_new[lrud[0]][j] + V_new[i][lrud[2]] - V_new[i][lrud[3]] ) + 2*dt*fh(t, x, y);
	            U_new[i][j] = U_old[i][j] -lambda_c2( delta, i, j ) * ( H_new[lrud[1]][j] - H_new[lrud[0]][j] ) + 2*dt*fu(t, x, y);
	            V_new[i][j] = V_old[i][j] -lambda_c2( delta, i, j ) * ( H_new[i][lrud[2]] - H_new[i][lrud[3]] ) + 2*dt*fv(t, x, y);

			}
			y = MIN_Y;
		}

		// [TODO]: Move swapping step inside above for loop and remove swapping matrices entirley

		////
		// swap
		for (int i = 0; i < NX; i += 1){ // x iterator
			for (int j = 0; j < NY; j += 1){ // y iterator
				H_swap[i][j] = H_new[i][j];
				U_swap[i][j] = U_new[i][j];
				V_swap[i][j] = V_new[i][j];

				H_old[i][j] = H_new[i][j];
				U_old[i][j] = U_new[i][j];
				V_old[i][j] = V_new[i][j];

				H_new[i][j] = H_swap[i][j];
				U_new[i][j] = U_swap[i][j];
				V_new[i][j] = V_swap[i][j];
			}
		}

		/////
		// Write Timestep to File and Display progress
		int nc = N_FRAMES*(t/(MAX_T - MIN_T)); 
		if ( nc > frame_iter){
			write_timestep(frame_iter, t, H_new);
			frame_iter += 1;
		}

		int pc = 100*(t/(MAX_T - MIN_T)); 
		if ( pc > completion){

			if (completion % 5 == 0){
				cout << "t: " << t << endl;
				cout << "Completion Percentage : " <<  pc << endl;
			
				/*
				// Measure Convergence Rate
				if (MEASURE_CONVERGENCE){

					// time error caluclation and subtract from total

					if (prev_error == -1){
						// compute true solution and current error and store result
						//prev_error = get_error();
					}
					else {
						// compute true solution and display local convergence rate
						// cur_error = get_error();

					}
				}
				*/

			}
			completion = pc;

			// Stop Execution
			if (completion >= 100){

				// Display Program Total Execution Time
				cout << "Maximum time reached" << endl;
				high_resolution_clock::time_point t2 = high_resolution_clock::now();
    			auto duration = duration_cast<microseconds>( t2 - t1 ).count();
    			cout << "Execution Time: " << duration << endl;
    			cout << "Quitting Program. " << endl;
				return 0;
			}
		}


	}

	// Display Program Total Execution Time
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>( t2 - t1 ).count();
	cout << "Execution Time: " << duration << endl;
	cout << "Quitting Program. " << endl;
	return 0;
}