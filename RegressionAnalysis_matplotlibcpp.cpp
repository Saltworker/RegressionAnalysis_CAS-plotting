#include <iostream>
#include <cstdlib>      // Get commands // system("pause"); // and // system("cls"); //
#include <fstream>      // Read .txt files
#include <sstream>      // String manipulation

#include <filesystem>   // Utilized to find repository in use
#include <sys/stat.h>   // Utilized for file search function. NOTE: Edit .cpp file name alongside (l. 578) to prevent errors
namespace fs = std::filesystem;

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

class regression { //Default access specifier is // private: //
    double xAxis[1000]{}, yAxis[1000]{};    // Max number of datapoints set to 1000
    double xAxisSorted[1000]{}, yAxisSorted[1000]{}; // Dataset numerically sorted. Utilized under statistical decriptors
    
    const double e = 2.718281828459045;     // Euler's value approximated 

    // Working out values in regression analysis //
    double sumX = 0;
    double sumY = 0;
    double sumXX = 0;
    double sumXY = 0;

    double SSR = 0;
    double SST = 0;

    bool expFejl = false;                   // Can the eksponential regression model be found (without error due to dataset)?
    bool logFejl = false;                   // Can the logarithmic regression model ...?
    bool powFejl = false;                   // Can the power regression model ...?

    std::string RRName[5] = { "linear     ", "exponential", "logarithmic", "power      ", "polynomial " };  // For ranking purposes
    
    std::vector<double> x, y, x_lin, y_lin, x_exp, y_exp, x_log, y_log, x_pow, y_pow, x_poly, y_poly; // Plotting datapoints

public:
    int counter = 0;    // number of datapoints in dataset

    // Linear regression model: y = aLin x + bLin //
    double aLin = 0;    // Slope of linear regression model
    double bLin = 0;    // Intercept of lineær regressionsmodel

    // Exponential regression model: y = aExp e^(cExp * x) //
    double aExp = 0;    // Startværdi af eksponentiel regressionsmodel
    double cExp = 0;    // Fremskrivningsfaktor af eksponentiel regressionsmodel

    // Logarithmic regression model: y = aLog ln(x) + bLog //
    double aLog = 0;    // "Fremskrivningsfaktor"??? af logaritmisk regressionsmodel
    double bLog = 0;    // startværdi af logaritmisk regressionsmodel

    // Power regression model: y = aPot x^cPot //
    double aPow = 0;    // Startværdi af potentiel regressionsmodel
    double cPow = 0;    // Eksponent af potentiel regressionsmodel

    // Second-degree polynomium regression model: a[2] x^2 + a[1] x + a[0] //
    // NOTE: Second-degree polynomium DOES NOT UNDERGO RANKING             //
    double a[3];        // Second-degree polynomium values: a[2] = second-degree term ; a[1] = Slope ; a[0] = Intercept

    double avgX = 0; // Average x-value
    double avgY = 0; // Average y-value

    double varianceY = 0, varianceX = 0;                // Variance VAR-values for x- and y-axis. Used to find (sample) standard deviation(s)

    double RR[5] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };    // Obtained like so => { linRR(), expRR(), logRR(), powRR(), polyRR() };
    double RR_ranked[5];                                // RR[5] but ranked (sorted)
    int RRIndex[5] = { 0, 1, 2, 3, 4 };                 // Exists solely for ranking purposes

    std::string filename;                               // Used to name plots

    /********************************************************************
    *   The constructor does the following:                             *
    *   1. Import datapoints from .txt file                             *
    *       a. Get // xAxis[] // and // xAxis[] //                      *
    *       b. Get number of datapoints as // counter //                *
    *   2. Calculate regression models as // models() //                *
    *   3. Calculate respective R^2-values as // RR_models() //         *
    *   4. Process dataset for use regarding statistical descriptors    *
    ********************************************************************/
    regression(std::string filename) : filename(filename) {
        // PART 1 //
        std::ifstream dataset;      // Initialize ifstream (only reading of the file)
        dataset.open(filename);     // Open dataset file

        if (dataset.is_open()) {    // Does the dataset file exist without the same repository as the .cpp file resides?
            std::string point;      // where "point" is the datapoint itself provided as a singular string

            while (getline(dataset, point)) {           // While there still is a readable datapoint [from up to down]
                std::stringstream s(point);             // std::stringstream separates by INPUT SPACEBAR
                int P = 0;                              // x whilst P=0 and y whilst P=1
                std::string temp, x, y;
                while (s >> temp) {                             // While a readable string value is present within "counter"th line
                    if (P == 0) {                               // Read x-value
                        x = temp;
                        if (x[x.length() - 1] == ',') {         // Commas in between x- and y-values are deleted
                            x.pop_back();
                        }
                        for (int i = 0; i < x.length(); i++) {      // Comma is converted to a period . in x-axis
                            if (x[i] == ',') {
                                x[i] = '.';
                            }
                        }
                        xAxis[counter] = stof(x);
                        if (xAxis[counter] <= 0) {          // Account for potential errors in the following regression models
                            logFejl = true;
                            powFejl = true;
                        }
                    }
                    if (P == 1) {                               //Read y-value
                        y = temp;
                        for (int i = 0; i < y.length(); i++) {      // Comma is converted to a period . in x-axis
                            if (y[i] == ',') {
                                y[i] = '.';
                            }
                        }
                        yAxis[counter] = stof(y);
                        if (yAxis[counter] <= 0) {          // Account for potential errors in the following regression models
                            expFejl = true;
                            powFejl = true;
                        }
                    }
                    P++;
                }
                counter++;                              // count up
                if (counter > 1000) {                   // Account for datapoint limit
                    std::cout << "Attention: Max number of datapoints set to " << counter <<
                                 "... Consider changing that under .cpp (l. 16) & (l. 117) " << "\n";
                    break;                              // Prevent data overflow. Make sure no more than the set datapoint limit is reached
                }
                std::cout << "Datapoint nr." << counter << "\t[" << x << ";" << y << "]\n";
            }
            std::cout << "Dataset length: " << counter << "\n";

            if (counter < 3) {
                std::cout << "ERROR: Dataset must have at least 3 datapoints" << "\n";
                exit(0);
            }

            dataset.close();
        }
        else {
            std::cout << "ERROR: Datafile could not be found" << "\n";
            exit(0);
        }

        // PART 2 //
        models();

        // PART 3 //
        RR_models();

        // PART 4 //
        xAxisSort();
        yAxisSort();
        XYVariance_standardDeviation();
    }

private:        // The following private functions are made so because they are utilized by the constructor alone
    void linReg() { // Make sure linReg goes last to get correct working values under statistical descriptors
        sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;
        for (int i = 0; i < counter; i++) {
            sumX = sumX + xAxis[i];
            sumXX = sumXX + xAxis[i] * xAxis[i];
            sumY = sumY + yAxis[i];
            sumXY = sumXY + xAxis[i] * yAxis[i];
        }
        aLin = (counter * sumXY - sumX * sumY) / (counter * sumXX - sumX * sumX);
        bLin = (sumXX * sumY - sumX * sumXY) / (counter * sumXX - sumX * sumX);
    }

    void expReg() {
        sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;
        for (int i = 0; i < counter; i++) {
            sumX = sumX + xAxis[i];
            sumXX = sumXX + xAxis[i] * xAxis[i];
            sumY = sumY + log(yAxis[i]);
            sumXY = sumXY + xAxis[i] * log(yAxis[i]);
        }
        aExp = (counter * sumXY - sumX * sumY) / (counter * sumXX - sumX * sumX);
        double bExp = (sumXX * sumY - sumX * sumXY) / (counter * sumXX - sumX * sumX);
        cExp = pow(e, bExp); // Back conversion from natural logarithm ln [Yep... log() acutally computes the natural logarithm]
    }

    void logReg() {
        sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;
        for (int i = 0; i < counter; i++) {
            sumX = sumX + log(xAxis[i]);
            sumXX = sumXX + log(xAxis[i]) * log(xAxis[i]);
            sumY = sumY + yAxis[i];
            sumXY = sumXY + log(xAxis[i]) * yAxis[i];
        }
        aLog = (counter * sumXY - sumX * sumY) / (counter * sumXX - sumX * sumX);
        bLog = (sumXX * sumY - sumX * sumXY) / (counter * sumXX - sumX * sumX);
    }
   
    void powReg() {
        sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;
        for (int i = 0; i < counter; i++) {
            sumX = sumX + log10(xAxis[i]);
            sumXX = sumXX + log10(xAxis[i]) * log10(xAxis[i]);
            sumY = sumY + log10(yAxis[i]);
            sumXY = sumXY + log10(xAxis[i]) * log10(yAxis[i]);
        }
        aPow = (counter * sumXY - sumX * sumY) / (counter * sumXX - sumX * sumX);
        double bPot = (sumXX * sumY - sumX * sumXY) / (counter * sumXX - sumX * sumX);
        cPow = pow(10, bPot); // Back conversion from log10
    }

    void secondDegreePolyReg() {    // Made thanks to: https://www.bragitoff.com/2015/09/c-program-for-polynomial-fit-least-squares/
        const int n = 2;            // The degree of the polynomium. Adjust alongside (l. 58), (l. 296), (l. 457) and (l. 649)

        double X[2 * n + 1];
        for (int i = 0; i < 2 * n + 1; i++) {
            X[i] = 0;
            for (int j = 0; j < counter; j++)
                X[i] = X[i] + pow(xAxis[j], i);
        }
        double B[n + 1][n + 2];
        for (int i = 0; i <= n; i++)
            for (int j = 0; j <= n; j++)
                B[i][j] = X[i + j];
        double Y[n + 1];
        for (int i = 0; i < n + 1; i++) {
            Y[i] = 0;
            for (int j = 0; j < counter; j++)
                Y[i] = Y[i] + pow(xAxis[j], i) * yAxis[j];
        }
        for (int i = 0; i <= n; i++)
            B[i][n + 1] = Y[i]; 
        int nn = n + 1;
        
        for (int i = 0; i < nn; i++)
            for (int k = i + 1; k < nn; k++)
                if (B[i][i] < B[k][i])
                    for (int j = 0; j <= nn; j++) {
                        double temp = B[i][j];
                        B[i][j] = B[k][j];
                        B[k][j] = temp;
                    }

        for (int i = 0; i < nn - 1; i++)            
            for (int k = i + 1; k < nn; k++) {
                double t = B[k][i] / B[i][i];
                for (int j = 0; j <= nn; j++)
                    B[k][j] = B[k][j] - t * B[i][j];
            }
        
        for (int i = nn - 1; i >= 0; i--) {
            a[i] = B[i][nn];
            for (int j = 0; j < nn; j++)
                if (j != i)
                    a[i] = a[i] - B[i][j] * a[j];
            a[i] = a[i] / B[i][i];
        }
    }

    void models() { // Calculate all available (possible) regression models
        linReg();
        secondDegreePolyReg();
        if (expFejl == false)
            expReg();
        if (logFejl == false)
            logReg();
        if (powFejl == false)
            powReg();
    }

    void getAverageY() {
        for (int i = 0; i < counter; i++) { // What you think of when you say the "average" or "mean" in a dataset
            avgY = avgY + yAxis[i];
        }
        avgY = avgY / counter;
    }
    void getAverageX() {
        for (int i = 0; i < counter; i++) { // Average/mean x-value. The little brother of the average/mean y-value
            avgX = avgX + xAxis[i];
        }
        avgX = avgX / counter;
    }

    void linRR() {
        SSR = 0, SST = 0; // SSR = Sum Squared Regression; SST = Total Sum of Squares
        for (int i = 0; i < counter; i++) {
            SSR = SSR + pow(yAxis[i] - linModel(xAxis[i]), 2);
            SST = SST + pow(yAxis[i] - avgY, 2);
        }

        RR[0] = 1 - SSR / SST;
    }

    void expRR() {
        SSR = 0, SST = 0;

        for (int i = 0; i < counter; i++) {
            SSR = SSR + pow(yAxis[i] - expModel(xAxis[i]), 2);
            SST = SST + pow(yAxis[i] - avgY, 2);
        }

        RR[1] = 1 - SSR / SST;
    }

    void logRR() {
        SSR = 0, SST = 0;

        for (int i = 0; i < counter; i++) {
            SSR = SSR + pow(yAxis[i] - logModel(xAxis[i]), 2);
            SST = SST + pow(yAxis[i] - avgY, 2);
        }

        RR[2] = 1 - SSR / SST;
    }

    void powRR() {
        SSR = 0, SST = 0;

        for (int i = 0; i < counter; i++) {
            SSR = SSR + pow(yAxis[i] - potModel(xAxis[i]), 2);
            SST = SST + pow(yAxis[i] - avgY, 2);
        }

        RR[3] = 1 - SSR / SST;
    }

    void polyRR() {
        SSR = 0, SST = 0;

        for (int i = 0; i < counter; i++) {
            SSR = SSR + pow(yAxis[i] - polyModel(xAxis[i]), 2);
            SST = SST + pow(yAxis[i] - avgY, 2);
        }

        RR[4] = 1 - SSR / SST;
    }

    void RR_selectionSort(int n) {
        for (int i = 0; i < 5; i++)
            RR_ranked[i] = RR[i];

        int max_index;

        for (int i = 0; i < n - 1; i++) {
            max_index = i;
            for (int j = i + 1; j < n; j++)
                if (RR[j] > RR[max_index])
                    max_index = j;

            // Swap the biggest element with the i-value //
            if (max_index != i) {
                std::swap(RR_ranked[max_index], RR_ranked[i]);
                std::swap(RRName[max_index], RRName[i]);
                std::swap(RRIndex[max_index], RRIndex[i]);
            }
        }
    }

    void RR_models() {          // Calculate all available (possible) R^2 values
        getAverageY();
        getAverageX();          //Can be executed later

        linRR();
        polyRR();
        if (expFejl == false)
            expRR();
        if (logFejl == false)
            logRR();
        if (powFejl == false)
            powRR();
        RR_selectionSort(4);    // Second-degree polynomium DOES NOT UNDERGO RANKING, thus the parameter 4 is used
    }


    std::string RRModel_scientific(int index, bool scientific) { // Regression models using either scientific or fixed (normal) notation
        if (scientific == true) {
            std::cout << std::scientific; // Make sure only the following is printed using scientific notation ...
        }
        else
            std::cout << std::fixed;
        switch (index) {
        case 4:
            std::cout << "y=" << a[2] << "x^2+" << a[1] << "x+" << a[0]; // Polynomial regression model is cramped to make up for the size of the formula
            break;
        case 3:
            std::cout << "y = " << cPow << " x ^ " << aPow;
            break;
        case 2:
            std::cout << "y = " << bLog << " + " << aLog << " ln(x)";
            break;
        case 1:
            std::cout << "y = " << cExp << " e ^ " << aExp << " x";
            break;
        case 0:
            std::cout << "y = " << aLin << " x + " << bLin;
            break;
        }
        std::cout << std::fixed; // ... like this
        return " ";
    }
    void xAxisSort() { //Selection sorting
        for (int i = 0; i < counter; i++)
            xAxisSorted[i] = xAxis[i];

        int min_index;

        for (int i = 0; i < counter - 1; i++) {
            min_index = i;
            for (int j = i + 1; j < counter; j++) {
                if (xAxisSorted[j] < xAxisSorted[min_index])
                    min_index = j;
            }

            // Swap the biggest element with the i-value //
            if (min_index != i)
                std::swap(xAxisSorted[min_index], xAxisSorted[i]);
        }
    }

    void yAxisSort() { //Selection sorting
        for (int i = 0; i < counter; i++)
            yAxisSorted[i] = yAxis[i];

        int min_index;

        for (int i = 0; i < counter - 1; i++) {
            min_index = i;
            for (int j = i + 1; j < counter; j++) {
                if (yAxisSorted[j] < yAxisSorted[min_index])
                    min_index = j;
            }

            // Swap the biggest element with the i-value //
            if (min_index != i)
                std::swap(yAxisSorted[min_index], yAxisSorted[i]);
        }
    }

    void XYVariance_standardDeviation() {
        for (int i = 0; i < counter; i++) {
            varianceX = varianceX + pow(xAxis[i] - avgX, 2);
            varianceY = varianceY + pow(yAxis[i] - avgY, 2);
        }
    }

public:
    void RR_Ranking(int n) { // With highest R^2 value as 1. and lowest as 4.
        std::cout << "REGRESSION MODEL RANKING:" << "\n";
        std::cout << "0. The second-degree polynomium\t\t" << RRModel_scientific(RRIndex[4], true) << "\tR^2 = " << RR_ranked[4] << "\n{ " << RRModel_scientific(RRIndex[4], false) << "}\n\n";

        for (int i = 0; i < n; i++) {
            if (RR[i] != 0.0f)
                std::cout << i + 1 << ". The " << RRName[i] << " regression model\t" << RRModel_scientific(RRIndex[i], true)
                << "\tR^2 = " << RR_ranked[i] << "\n{ " << RRModel_scientific(RRIndex[i], false) << "}\n\n";
        }

        // Check if certain regression models are - can be - used //
        if (expFejl == true)
            std::cout << "WARNING: The exponential regressionsmodel cannot be made (y-axis contains non-postive value)\n";
        if (logFejl == true)
            std::cout << "WARNING: The logarithmic regressionsmodel cannot be made (x-axis contains non-postive value)\n";
        if (powFejl == true)
            std::cout << "WARNING: The power regressionsmodel cannot be made (x- or y-axis contains non-postive value)\n";
    }

    double linModel(double xValue) { return aLin * xValue + bLin; }

    double expModel(double xValue) { return cExp * pow(e, aExp * xValue); }

    double potModel(double xValue) { return cPow * pow(xValue, aPow); }

    double logModel(double xValue) { return aLog * log(xValue) + bLog; }

    double polyModel(double xValue) { return a[2] * pow(xValue, 2) + a[1] * xValue + a[0]; }

    double YQn(float Q) {
        if ((counter * Q == int(counter * Q)))
            return (yAxisSorted[int(Q * counter)] + yAxisSorted[int(Q * counter - 1)]) / 2;
        else
            return yAxisSorted[int(Q * counter)];
    }

    double XQn(float Q) {
        if ((counter * Q == int(counter * Q)))
            return (xAxisSorted[int(Q * counter)] + xAxisSorted[int(Q * counter - 1)]) / 2;
        else
            return xAxisSorted[int(Q * counter)];
    }

    void Xoutliers() { // Find the number of datapoints that lie below 2/3 * Q1 and above 1.5 * Q3
        int outlier_count = 0;
        std::list<double> outliers{};

        double Q1 = XQn(0.25);
        double Q3 = XQn(0.75);
        for (int i = 0; i < counter; i++) {
            if (xAxisSorted[i] < Q1 - 1.5 * (Q3 - Q1)) {
                outlier_count++;
                outliers.push_front(xAxisSorted[i]); // Add the outlier to the list
            }
                
            if (xAxisSorted[i] > Q3 + 1.5 * (Q3 - Q1)) {
                outlier_count++;
                outliers.push_front(xAxisSorted[i]); // Add the outlier to the list
            }
        }
        std::cout << "(X) Outlier count:\t" << std::to_string(outlier_count) << "\t";

        for (auto i : outliers)
            std::cout << "(" << i << ") ";
        std::cout << "\n";
    }
    
    void Youtliers() { // Find the number of datapoints that lie below 2/3 * Q1 and above 1.5 * Q3
        int outlier_count = 0;
        std::list<double> outliers{};

        double Q1 = YQn(0.25);
        double Q3 = YQn(0.75);
        for (int i = 0; i < counter; i++) {
            if (yAxisSorted[i] < Q1 - 1.5 * (Q3 - Q1)) {
                outlier_count++;
                outliers.push_front(yAxisSorted[i]); // Add the outlier to the list
            }

            if (yAxisSorted[i] > Q3 + 1.5 * (Q3 - Q1)) {
                outlier_count++;
                outliers.push_front(yAxisSorted[i]); // Add the outlier to the list
            }
        }
        std::cout << "(Y) Outlier count:\t" << std::to_string(outlier_count) << "\t";

        for (auto i : outliers)
            std::cout << "(" << i << ") ";
        std::cout << "\n";
    }

    double linCorrelationCoefficient() { // r-value for linear regression model
        // Formula: https://www.bmj.com/about-bmj/resources-readers/publications/statistics-square-one/11-correlation-and-regression
        double XY = 0;
        for (int i = 0; i < counter; i++)
            XY = XY + xAxis[i] * yAxis[i];
        return (XY - counter * avgX * avgY) / ( counter * sqrt(varianceX / counter) * sqrt(varianceY / counter) ); 
    }

    void statistical_decriptors() {
        std::cout << "STATISTICAL DESCRIPTORS:\n";
        std::cout << "Number of datapoints: " << counter << "\n\n";

        std::cout << "(X) Domain:\t[" << xAxisSorted[0] << ";" << xAxisSorted[counter - 1] << "]\n";
        std::cout << "(X) Quarters:\t(Q1 = " << XQn(0.25) << ";\tQ2 (median) = " << XQn(0.5) << ";\tQ3 = " << XQn(0.75) << ")\n";
        std::cout << "(X) Spread:\t{Variance omega^2 = " << varianceX << ";\n";
        std::cout << "(X)\t\tStandard deviation omega = " << sqrt(varianceX / counter) << ";\n"; //Lowercase omega, alternatively, goes: \u03A3
        std::cout << "(X)\t\tStandard deviation s = " << sqrt(varianceX / (counter - 1)) << "}\n";
        Xoutliers();
        std::cout << "(X) Average/Mean = " << avgX << "\n\n";

        std::cout << "(Y) Range:\t[" << yAxisSorted[0] << ";" << yAxisSorted[counter - 1] << "]\n";
        std::cout << "(Y) Quarters:\t(Q1 = " << YQn(0.25) << ";\tQ2 (median) = " << YQn(0.5) << ";\tQ3 = " << YQn(0.75) << ")\n";
        std::cout << "(Y) Spread:\t{Variance omega^2 = " << varianceY << ";\n";
        std::cout << "(Y)\t\tStandard deviation omega = " << sqrt(varianceY / counter) << ";\n"; //Lowercase omega, alternatively, goes: \u03A3
        std::cout << "(Y)\t\tSample standard deviation s = " << sqrt(varianceY / (counter - 1)) << "}\n";
        Youtliers();
        std::cout << "(Y) Average/Mean = " << avgY << "\n\n";

        std::cout << "Linear correlation coefficient r:\t" << linCorrelationCoefficient() << "\n\n";
        
        std::cout << "LATEST CALCULATIONS:\n";
        std::cout << "SumX = " << sumX << "\n";
        std::cout << "SumY = " << sumY << "\n";
        std::cout << "SumXX = " << sumXX << "\n";
        std::cout << "SumXY = " << sumXY << "\n";
        std::cout << "SSR = " << SSR << "\n";
        std::cout << "SST = " << SST << "\n\n";
    }

private:
    void plotClear() { x.clear(); y.clear(); x_lin.clear(); y_lin.clear(); x_exp.clear(); y_exp.clear(); 
                       x_log.clear(); y_log.clear(); x_pow.clear(); y_pow.clear(); x_poly.clear(); y_poly.clear(); }

public:
    void plot(std::string lin, std::string exp, std::string log, std::string pow, std::string poly) {
        plotClear();

        for (int i = 0; i < counter; i++) {         // Set dataset up for plotting
            x.push_back(xAxis[i]);
            y.push_back(yAxis[i]);
        }
        for (int i = 0; i < counter; i++) {         // Set linear regression model up for plotting
            x_lin.push_back(xAxis[i]);
            y_lin.push_back(linModel(xAxis[i]));
        }

        if (expFejl == false) {
            for (int i = 0; i < counter; i++) {     // Set exponential regression model up for plotting
                x_exp.push_back(xAxis[i]);
                y_exp.push_back(expModel(xAxis[i]));
            }
        }

        if (logFejl == false) {
            for (int i = 0; i < counter; i++) {     // Set logarithmic regression model up for plotting
                x_log.push_back(xAxis[i]);
                y_log.push_back(logModel(xAxis[i]));
            }
        }

        if (powFejl == false) {
            for (int i = 0; i < counter; i++) {     // Set power regression model up for plotting
                x_pow.push_back(xAxis[i]);
                y_pow.push_back(potModel(xAxis[i]));
            }
        }
        for (int i = 0; i < counter; i++) {         // Set polynomial regression model up for plotting
            x_poly.push_back(xAxis[i]);
            y_poly.push_back(polyModel(xAxis[i]));
        }

        plt::plot(x, y, "bD");
        plt::plot(x_poly, y_poly, { {"color", "#000000"}, {"ls", ":"}, {"label", poly} });
        plt::plot(x_lin, y_lin, { {"color", "#FF0000"}, {"ls", ":"}, {"label", lin} });
        if (expFejl == false)
            plt::plot(x_exp, y_exp, { {"color", "#00FF00"}, {"ls", ":"}, {"label", exp} });
        if (logFejl == false)
            plt::plot(x_log, y_log, { {"color", "#FF00FF"}, {"ls", ":"}, {"label", log} });
        if (powFejl == false)
            plt::plot(x_pow, y_pow, { {"color", "#00FFFF"}, {"ls", ":"}, {"label", pow} });
        plt::xlabel("x");
        plt::ylabel("y");
        plt::legend();
        plt::title(filename);
    }
    void plotResidual(std::string lin, std::string exp, std::string log, std::string pow, std::string poly) {
        plotClear();

        for (int i = 0; i < counter; i++) {         // Set linear regression model up for plotting
            x_lin.push_back(xAxis[i]);
            y_lin.push_back(yAxis[i] - linModel(xAxis[i]));
        }

        if (expFejl == false) {
            for (int i = 0; i < counter; i++) {     // Set exponential regression model up for plotting
                x_exp.push_back(xAxis[i]);
                y_exp.push_back(yAxis[i] - expModel(xAxis[i]));
            }
        }

        if (logFejl == false) {
            for (int i = 0; i < counter; i++) {     // Set logarithmic regression model up for plotting
                x_log.push_back(xAxis[i]);
                y_log.push_back(yAxis[i] - logModel(xAxis[i]));
            }
        }

        if (powFejl == false) {
            for (int i = 0; i < counter; i++) {     // Set power regression model up for plotting
                x_pow.push_back(xAxis[i]);
                y_pow.push_back(yAxis[i] - potModel(xAxis[i]));
            }
        }
        for (int i = 0; i < counter; i++) {         // Set polynomial regression model up for plotting
            x_poly.push_back(xAxis[i]);
            y_poly.push_back(yAxis[i] - polyModel(xAxis[i]));
        }
        std::vector < double> zero = { xAxis[0]};
        plt::plot(zero, zero, { {"color", "#000000"}, {"ls", ":"}, {"label", poly} });
        plt::plot(zero, zero, { {"color", "#FF0000"}, {"ls", ":"}, {"label", lin} });
        if (expFejl == false)
            plt::plot(zero, zero, { {"color", "#00FF00"}, {"ls", ":"}, {"label", exp} });
        if (logFejl == false)
            plt::plot(zero, zero, { {"color", "#FF00FF"}, {"ls", ":"}, {"label", log} });
        if (powFejl == false)
            plt::plot(zero, zero, { {"color", "#00FFFF"}, {"ls", ":"}, {"label", pow} });

        plt::plot(x_poly, y_poly, "k*");
        plt::plot(x_lin, y_lin, "r*");
        if (expFejl == false)
            plt::plot(x_exp, y_exp, "g*");
        if (logFejl == false)
            plt::plot(x_log, y_log, "m*");
        if (powFejl == false)
            plt::plot(x_pow, y_pow, "c*");
        plt::xlabel("x");
        plt::ylabel("Residual");
        plt::legend();
        plt::title(filename);
    }
};

int main()
{
    std::string filename;
    struct stat sb; // (l. 7)

    // Find the absolute path of the following file/program //
    std::filesystem::path filelocation("RegressionAnalysis_matplotlib.cpp");

    int input;
    double xValue;

    while (1) {
        std::cout << "Repository: \n" << fs::absolute(filelocation).remove_filename().string() << "\n";

        // Search for file to use in parameterized constructor //
        while (1) {
            std::cout << "\nType in filename (rememeber .txt):\n";
            std::cin >> filename;
            std::cout << fs::absolute(filelocation).remove_filename().string() << filename << "\n";
            if (stat(filename.c_str(), &sb) == 0 && !(sb.st_mode & S_IFDIR))
                break;
            else
                std::cout << "The file does not exist. Try again\n";
        }

        regression reg(filename); // e.g. reg("dataset.txt");
        std::cout << "Dataset has been imported!\n";
        system("pause");    // Press any key to continue...
        system("cls");      // Clear screen

        do {
            reg.RR_Ranking(4);

            std::cout << "\nTYPE: \n\"1\" to use linear regression model\n";
            std::cout << "\"2\" to use exponential regression model\n";
            std::cout << "\"3\" to use logarithmic regression model\n";
            std::cout << "\"4\" to use power regression model\n";
            std::cout << "\"5\" to use second-degree polynomial regression model\n";
            std::cout << "\"6\" to see statistical descriptors\n";
            std::cout << "\"7\" to see plot of regressionsmodels\t(\"77\" to see residual plot)\t(\"777\" to see both)\n";
            std::cout << "\"8\" to adjust floating point numbers\n";
            std::cout << "\"9\" to select new dataset\n";
            std::cout << "\"0\" to end program\n";
            std::cin >> input;
            switch (input) {
            case 9:
                reg.~regression();      // User-defined destructor to release memory when selecting new dataset
                system("cls");      // Clear screen
                break;
            case 8:
                std::cout << "Enter amount of floating point numbers:\n";
                std::cin >> xValue;
                std::cout.precision(xValue); // Set amount of floating point numbers. Defaults to 6
                system("cls");      // Clear screen
                break;
            case 7: //https://matplotlib-cpp.readthedocs.io/en/latest/docs.html
                reg.plot(std::string("Lin: ") + "y = " + std::to_string(reg.aLin) + " x + " + std::to_string(reg.bLin) + std::string(", R^2 = ") + std::to_string(reg.RR[0]),
                    std::string("Exp: ") + "y = " + std::to_string(reg.cExp) + " e ^ " + std::to_string(reg.aExp) + "x" + std::string(", R^2 = ") + std::to_string(reg.RR[1]),
                    std::string("Log: ") + "y = " + std::to_string(reg.bLog) + " + " + std::to_string(reg.aLog) + " ln(x)" + std::string(", R^2 = ") + std::to_string(reg.RR[2]),
                    std::string("Pow: ") + "y = " + std::to_string(reg.cPow) + " x ^ " + std::to_string(reg.aPow) + std::string(", R^2 = ") + std::to_string(reg.RR[3]),
                    std::string("Poly: ") + "y = " + std::to_string(reg.a[2]) + "x^2 + " + std::to_string(reg.a[1]) + "x + " + std::to_string(reg.a[0]) + std::string(", R^2 = ") + std::to_string(reg.RR[4]));
                plt::show();
                system("pause");    // Press any key to continue...
                system("cls");      // Clear screen
                break;
            case 77:
                reg.plotResidual(std::string("Lin: ") + "y = " + std::to_string(reg.aLin) + " x + " + std::to_string(reg.bLin) + std::string(", R^2 = ") + std::to_string(reg.RR[0]),
                    std::string("Exp: ") + "y = " + std::to_string(reg.cExp) + " e ^ " + std::to_string(reg.aExp) + "x" + std::string(", R^2 = ") + std::to_string(reg.RR[1]),
                    std::string("Log: ") + "y = " + std::to_string(reg.bLog) + " + " + std::to_string(reg.aLog) + " ln(x)" + std::string(", R^2 = ") + std::to_string(reg.RR[2]),
                    std::string("Pow: ") + "y = " + std::to_string(reg.cPow) + " x ^ " + std::to_string(reg.aPow) + std::string(", R^2 = ") + std::to_string(reg.RR[3]),
                    std::string("Poly: ") + "y = " + std::to_string(reg.a[2]) + "x^2 + " + std::to_string(reg.a[1]) + "x + " + std::to_string(reg.a[0]) + std::string(", R^2 = ") + std::to_string(reg.RR[4]));
                plt::show();
                system("pause");    // Press any key to continue...
                system("cls");      // Clear screen
                break;
            case 777:
                reg.plot(std::string("Lin: ") + "y = " + std::to_string(reg.aLin) + " x + " + std::to_string(reg.bLin) + std::string(", R^2 = ") + std::to_string(reg.RR[0]),
                    std::string("Exp: ") + "y = " + std::to_string(reg.cExp) + " e ^ " + std::to_string(reg.aExp) + "x" + std::string(", R^2 = ") + std::to_string(reg.RR[1]),
                    std::string("Log: ") + "y = " + std::to_string(reg.bLog) + " + " + std::to_string(reg.aLog) + " ln(x)" + std::string(", R^2 = ") + std::to_string(reg.RR[2]),
                    std::string("Pow: ") + "y = " + std::to_string(reg.cPow) + " x ^ " + std::to_string(reg.aPow) + std::string(", R^2 = ") + std::to_string(reg.RR[3]),
                    std::string("Poly: ") + "y = " + std::to_string(reg.a[2]) + "x^2 + " + std::to_string(reg.a[1]) + "x + " + std::to_string(reg.a[0]) + std::string(", R^2 = ") + std::to_string(reg.RR[4]));
                plt::figure(2); //Make a seperate plot. Default parameter is 1
                reg.plotResidual(std::string("Lin: ") + "y = " + std::to_string(reg.aLin) + " x + " + std::to_string(reg.bLin) + std::string(", R^2 = ") + std::to_string(reg.RR[0]),
                    std::string("Exp: ") + "y = " + std::to_string(reg.cExp) + " e ^ " + std::to_string(reg.aExp) + "x" + std::string(", R^2 = ") + std::to_string(reg.RR[1]),
                    std::string("Log: ") + "y = " + std::to_string(reg.bLog) + " + " + std::to_string(reg.aLog) + " ln(x)" + std::string(", R^2 = ") + std::to_string(reg.RR[2]),
                    std::string("Pow: ") + "y = " + std::to_string(reg.cPow) + " x ^ " + std::to_string(reg.aPow) + std::string(", R^2 = ") + std::to_string(reg.RR[3]),
                    std::string("Poly: ") + "y = " + std::to_string(reg.a[2]) + "x^2 + " + std::to_string(reg.a[1]) + "x + " + std::to_string(reg.a[0]) + std::string(", R^2 = ") + std::to_string(reg.RR[4]));
                plt::show();
                system("pause");    // Press any key to continue...
                system("cls");      // Clear screen
                break;
            case 6:
                system("cls");      // Clear screen
                reg.statistical_decriptors();
                system("pause");    // Press any key to continue...
                system("cls");      // Clear screen
                break;
            case 5:
                std::cout << "Second-degree polynomial:\nType in x-value:\n";
                std::cin >> xValue;
                std::cout << reg.a[2] << " * " << xValue << "^2 + " << reg.a[1] << " * " << xValue << " + " << reg.a[0] << " = " << reg.polyModel(xValue);
                system("pause");    // Press any key to continue...
                system("cls");      // Clear screen
                break;
            case 4:
                std::cout << "Power \nType in x-value:\n";
                std::cin >> xValue;
                std::cout << reg.cPow << " * " << xValue << " ^ " << reg.aPow << " = " << reg.potModel(xValue) << "\n";
                system("pause");    // Press any key to continue...
                system("cls");      // Clear screen
                break;
            case 3:
                std::cout << "Logarithmic \nType in x-value:\n";
                std::cin >> xValue;
                std::cout << reg.aLog << " * " << "log(" << xValue << ") + " << reg.bLog << " = " << reg.logModel(xValue) << "\n";
                system("pause");    // Press any key to continue...
                system("cls");      // Clear screen
                break;
            case 2:
                std::cout << "Exponential \nType in x-value:\n";
                std::cin >> xValue;
                std::cout << reg.cExp << " * " << reg.aExp << " ^ " << xValue << " = " << reg.expModel(xValue) << "\n";
                system("pause");    // Press any key to continue...
                system("cls");      // Clear screen
                break;
            case 1:
                std::cout << "Linear \nType in x-value:\n";
                std::cin >> xValue;
                std::cout << reg.aLin << " * " << xValue << " + " << reg.bLin << " = " << reg.linModel(xValue) << "\n";
                system("pause");    // Press any key to continue...
                system("cls");      // Clear screen
                break;
            case 0:
                return 0;
            default:
                std::cout << "Case not found. Try again\n";
                system("pause");    // Press any key to continue...
                system("cls");      // Clear screen
            }
        } while (input != 9);
    }
};
