// Temperature_Forecast.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>

#include <fstream>
#include <vector>
#include <list>

#include <string>
#include <sstream>
#include <stdlib.h> 

#include "tnt.h"
#include "jama.h"
#include "jama_lu.h"
#include "jama_cholesky.h"

using namespace std;
using namespace TNT;
using namespace JAMA;

// Function to return temperature.txt
double** read_temperature() {
	ifstream MyFile("temperature.txt");

	// check whether the file opens
	if (MyFile.fail())
	{
		cout << "Error opening file" << endl;
	}

	int year, month, i;
	double T_max, T_min, rain, **input_doubles;
	char sun[100], af[100];

	// Lists to store temperatures, years, months
	list<double> Tmax; list<double> Tmin; list<int> Year; list<int> Month;


	while (MyFile >> year >> month >> T_max >> T_min >> af >> rain >> sun) {
		//cout << year << ", " << month<< ", " << T_max << ", " << T_min << endl;
		Year.push_back(year); Month.push_back(month); Tmax.push_back(T_max); Tmin.push_back(T_min);
	}

	// copy lists into an array
	//-----------------
	double *Years = new double[Year.size()]; //Year array
	copy(Year.begin(), Year.end(), Years);
	double *Months = new double[Month.size()]; //Month array
	copy(Month.begin(), Month.end(), Months);
	double *max_T = new double[Tmax.size()]; //max temperature array
	copy(Tmax.begin(), Tmax.end(), max_T);
	double *min_T = new double[Tmin.size()]; //min temperature array
	copy(Tmin.begin(), Tmin.end(), min_T);
	//------------------

	MyFile.close(); // Done with the file

	static double siz = Year.size();
	double *size = &siz;
	double *returner[5] = { Years, Months, max_T, min_T, size };
	return returner;
}

// Function to return data.txt
double** read_data() {
	ifstream MyFile("data.txt");

	// check whether the file opens
	if (MyFile.fail())
	{
		cout << "Error opening file" << endl;
	}

	double year, population, sunspot, co2;

	// Lists to store temperatures, years, months
	list<double> CO2; list<double> Year; list<double> Population; list<double> Sunspot;


	while (MyFile >> year >> population >> sunspot >> co2) {
		//cout << year << ", " << month<< ", " << T_max << ", " << T_min << endl;
		Year.push_back(year); Population.push_back(population); Sunspot.push_back(sunspot); CO2.push_back(co2);
	}

	// copy lists into an array
	//-----------------
	double *Years2 = new double[Year.size()]; //Year array
	copy(Year.begin(), Year.end(), Years2);
	double *Populations = new double[Population.size()]; //Month array
	copy(Population.begin(), Population.end(), Populations);
	double *Sunspots = new double[Sunspot.size()]; //max temperature array
	copy(Sunspot.begin(), Sunspot.end(), Sunspots);
	double *CO2s = new double[CO2.size()]; //min temperature array
	copy(CO2.begin(), CO2.end(), CO2s);
	//------------------

	MyFile.close(); // Done with the file

	static double sizz = Year.size();
	double *size2 = &sizz;
	double *returner[5] = { Years2, Populations, Sunspots, CO2s, size2 };
	return returner;
}

int main()
{
	// delete at the end
	double **returned = read_temperature();
	double *Years = returned[0];
	double *Months = returned[1];
	double *TMax = returned[2];
	double *TMin = returned[3];
	double *Len_data = returned[4];

	double **returned2 = read_data();
	double *Years2 = returned2[0];
	double *Populations = returned2[1];
	double *Sunspots = returned2[2];
	double *CO2s = returned2[3];
	double *Len_data2 = returned2[4];

	int i, j; //for the loops

	//--- These will be used when data is split into months
	double *temperature_min = new double[int(*Len_data / 12)];
	double *temperature_max = new double[int(*Len_data / 12)];
	double *temperature = new double[int(*Len_data / 12)]; 

	double *year = new double[int(*Len_data / 12)];
	double *month = new double[int(*Len_data / 12)]; 
	//---

	//---User Input Year and Month for prediction---
	int YY;
	cout << "Please enter a year (to predict the temperature anomaly): ";
	cin >> YY;

	int MM;
	cout << "Please enter a Month (1-January, 12-December): " ;
	cin >> MM;
	cout << endl;
	//---

	// Split the data into months
	//--------------
	int count = 0;
	int count2 = 0;
	for (i = 0; i < *Len_data; i++) {

		// get the Nth month
		if (count == (MM - 1))
		{
			temperature_max[count2] = TMax[i];
			temperature_min[count2] = TMin[i];
			year[count2] = Years[i];
			month[count2] = Months[i];

			count2 += 1;
		}
		count += 1;
		if (count == 12) {
			count = 0;
		}
	}

	// define temperature as (t_min - t_max)/2
	for (i = 0; i < *Len_data / 12; i++) {
		temperature[i] = (temperature_max[i] + temperature_min[i]) / 2;
	}
	//--------------

	//-------------------------------Main_prediction algorithms---------------------------------------------
	
	// Model1: Linear Regression (1st order polnomial) using Year and Temperature data only [note temperature.txt from 1948-2018]
	//------
	// Model: t = w0 + w1*x [t: temperature, x: year]
	// where, w1 = (mean(xt) - mean(x)mean(t))/(mean(x^2)-(mean(x))^2) 
	// and, w0 = mean(t) - w1*mean(x)
	//------
	//----------------------------------------------------------------------

	double x_mean, x2_mean, t_mean, xt_mean;
	double *xt = new double[int(*Len_data / 12)];
	double *x_squared = new double[int(*Len_data / 12)]; 

	// get means
	double sum = 0; double sum2 = 0; double sum3 = 0; double sum4 = 0;
	for (i = 0; i < *Len_data / 12; i++) {
		sum += year[i];
		sum2 += temperature[i];

		xt[i] = year[i] * temperature[i];
		sum3 += xt[i];

		x_squared[i] = year[i] * year[i];
		sum4 += x_squared[i];
	}

	x_mean = sum / (*Len_data / 12);
	t_mean = sum2 / (*Len_data / 12);
	xt_mean = sum3 / (*Len_data / 12);
	x2_mean = sum4 / (*Len_data / 12);

	// fitting constants
	double w1 = (xt_mean - x_mean * t_mean) / (x2_mean - x_mean * x_mean);
	double w0 = t_mean - w1 * x_mean;

	// prediction
	double T_predict = w0 + w1 * YY;
	double T_anomaly = T_predict - t_mean;

	//----Calculate standard error of our fittings wrt YY (using the last 20 data points)-----
	double standard_err2 = 0;
	for (i = 0; i < *Len_data/12; i++) {
		standard_err2 += pow(((temperature[i] - w0 - w1 * (i + 1948))), 2);
	}
	double standard_err = pow(standard_err2, 0.5)/pow(70,0.5);

	//-----------------------------------------------------------------------------------------

	// Model2: Multivariable Linear Regression (1st order polnomial) using Year, Temperature, Population, Sunspot #, and global CO2 data [note data.txt from 1961-2018]
	//------
	// Model: t = w0 + w1*x1^0.1  +w2*x2^0.1 +w3*x3^0.1 + w4*x4^0.1) [t: Temperature, x1: Year, x2:Population, x3: Sunspot #, x4: CO2]
	// where, w(optimal) = (X^T*X)^(-1) X^T * t,   [Note: X^T - means transpose(X)] and,
	// w^T=[wo, w1, w2, w3, w4]; X = [1   (x1[0])^0.1   (x2[0])^0.1     (x3[0])^0.1   (x4[0])^0.1    ;   t^T=[t0, t1, ..., tN]; where, ^T:transpose
    //                                .       .             .               .            .
	//                                .       .             .               .            . 
	//                                1   (x1[N])^0.1   (x2[N])^0.1     (x3[N])^0.1   (x4[N])^0.1]
	//------
	//----------------------------------------------------------------------

	double Pp = 0.1; // Power of the model
	// Initialize X
	Array2D<double> XX(int(*Len_data2), 5);
    // Fill the X  Matrix
	for (i = 0; i < *Len_data2; i++) {
		XX[i][0] = 1; XX[i][1] = pow(Years2[i], Pp); XX[i][2] = pow(Populations[i],Pp); XX[i][3] = pow(Sunspots[i],Pp); XX[i][4] = pow(CO2s[i],Pp);
		//XX[i][2] = pow(Years2[i], 1); XX[i][4] = pow(Populations[i],1); XX[i][6] = pow(Sunspots[i],1); XX[i][8] = pow(CO2s[i],1);
		}
	
	// Initialize Transpose(X): X^T
	Array2D<double> XXT(5, int(*Len_data2));
	// Fill the X^T Matrix
	for (i = 0; i < *Len_data2; i++) {
		for (j = 0; j < 5; j++) {
			XXT[j][i] = XX[i][j];
		}
	}

	// Initialize t
	Array2D<double> tt(int(*Len_data2), 1);
	// Fill t Matrix
	for (i = 0; i < *Len_data2; i++) {
		tt[i][0] = temperature[13 + i]; // 13 because temperature.txt start from year 1948 but, data.txt start from 1961
	}

	// set up identity matrix (will be used for inverting matrix (X^T *X))
	Array2D<double> eye(5, 5);
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			eye[i][j] = 0;
			eye[i][i] = 1;
		}
	}

	// Initialize W (optimal parameters for our model)
	Array2D<double> WW(5, 1);
	// Fill in W Matrix
	Array2D<double> WW1(5, 5);
	WW1 = (XXT*XX);
	LU<double> WW1LU(WW1);
	Array2D<double> WW1inv = WW1LU.solve(eye);

	// fitting constants
	WW = WW1inv * XXT *tt;

	// estimates for the future values of CO2, Sunspot #, Heathrow population (needed for future T prediction)
	double P_e = (202347-pow(YY,1.3)) + (1.3*pow(10, -91))*exp(0.109*YY); //CO2 estimate in year YY
	double C_e = -2720 + 1.54*YY; //CO2 estimate in year YY
	double pi = 3.14159265359;
	double S_e = 6.17 + 160.8*pow(sin(pi*(YY - 116.8) / 10.92),2); //Number of Sunspots estimate in year YY

	// Temperature Prediction in year YY (using our multivaruiable 2nd order polynomial)
	double t_predict2 = *WW[0] + *WW[1] * pow(YY, Pp) + *WW[2] * pow(P_e, Pp) + *WW[3] * pow(S_e, Pp) + *WW[4] * pow(C_e, Pp);

	double T_anomaly2 = t_predict2 - t_mean;

	//----------------------------------------------------------------------

	//---------------------Output-----------------------------
		// Output predictions on console
	cout << "Model_1: Temp = w0 +w1*year  where, w0=" << w0 << " and w1=" << w1 << endl;
	cout << "Temperature prediction for year:" << YY << " and Montht:" << MM << " is " << T_predict << " (in degrees Celsius);" << " estimated error on the predicted temperature ~ +-" << standard_err << endl;
	cout << "Temperature anomaly (Temp(year)-mean_Temp(1948-2018)) for year:" << YY << " and month:" << MM << " is " << T_anomaly << "+-" << standard_err << " (in degrees Celsius)" << endl << endl;

	cout << "Model_2: Temp = w0 + w1*year^0.1 + w2*population^0.1 + w3*#_of_sunspots^0.1 + w4*CO2_Conc^0.1  where, w0=" << *WW[0] << " w1=" << *WW[1] << " w2=" << *WW[2] << " w3=" << *WW[3] << " w4=" << *WW[4] << endl;
	cout << "Warning: This is an overfitted model, don't use it for long term predictions (may be used for predictions up to year 2025, for longer term strictly use Model_1)" << endl;
	cout << "Temperature prediction for year:" << YY << " and Montht:" << MM << " is " << t_predict2 << " (in degrees Celsius)" << endl;
	cout << "Temperature anomaly (Temp(year)-mean_Temp(1948-2018)) for year:" << YY << " and month:" << MM << " is " << T_anomaly2 << " (in degrees Celsius)" << endl;

	// Output predictions to the file predictions.txt
	ofstream myfile("prediction.txt");
	if (myfile.is_open())
	{
		myfile << "Model_1: Temp = w0 +w1*year  where, w0=" << w0 << " and w1=" << w1 << endl;
		myfile << "Temperature prediction for year:" << YY << " and Montht:" << MM << " is " << T_predict << " (in degrees Celsius);" << " estimated error on the predicted temperature ~ +-" << standard_err << endl;
		myfile << "Temperature anomaly (Temp(year)-mean_Temp(1948-2018)) for year:" << YY << " and month:" << MM << " is " << T_anomaly << "+-" << standard_err << " (in degrees Celsius)" << endl << endl;
		myfile << "Model_2: Temp = w0 + w1*year^0.1 + w2*population^0.1 + w3*#_of_sunspots^0.1 + w4*CO2_Conc^0.1  where, w0=" << *WW[0] << " w1=" << *WW[1] << " w2=" << *WW[2] << " w3=" << *WW[3] << " w4=" << *WW[4] << endl;
		myfile << "Warning: This is an overfitted model, don't use it for long term predictions (may be used for predictions up to year 2025, for longer term strictly use Model_1)" << endl;
		myfile << "Temperature prediction for year:" << YY << " and Montht:" << MM << " is " << t_predict2 << " (in degrees Celsius)" << endl;
		myfile << "Temperature anomaly (Temp(year)-mean_Temp(1948-2018)) for year:" << YY << " and month:" << MM << " is " << T_anomaly2 << " (in degrees Celsius)" << endl;
		myfile.close();
	}
	else cout << "Error: unable to open file prediction.txt";

	system("pause");
	//--------------------------------------------------------

	// clear the static memory used
	delete temperature_min, temperature_max, temperature, year, month, xt, x_squared;
	delete Years, Months, TMax, TMin, Len_data;
	delete Years2, Populations, Sunspots, CO2s, Len_data2;

	return 0;
}
