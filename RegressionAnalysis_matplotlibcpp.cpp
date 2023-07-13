/**  UPDATED July 13th 2023
 * 
 * PURPOSE:
 *
 *   Regression analysis tool (class regression) 
 *   utilizing matplotlibcpp.h to graph the results.
 *   Using .txt files as input, the following is found:
 *   - Regression models
 *       ~ Linear
 *       ~ Exponential
 *       ~ Logarithmic
 *       ~ Power
 *       ~ Second-Degree polynomium [Adjustablue in code]
 *   - Correlation coefficient(s) [R^2]
 *   - Graph of regression models
 *   - Other statistical descriptors [MORE PENDING]
 *
 *   TO BE DONE:
 *   1. Get scientific values of "inprecise" fixed values
 *   2. Get additional statistical descriptors
 * 
 * LICENSE:
 * 
 *   MIT License
 * 
 *   Copyright (c) 2023 Mathias Lykholt-Ustrup
 *
 *   Permission is hereby granted, free of charge, to any person obtaining a copy
 *   of this software and associated documentation files (the "Software"), to deal
 *   in the Software without restriction, including without limitation the rights
 *   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *   copies of the Software, and to permit persons to whom the Software is
 *   furnished to do so, subject to the following conditions:

 *   The above copyright notice and this permission notice shall be included in all
 *   copies or substantial portions of the Software.

 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *   SOFTWARE.
 * 
 *   @Author Saltworker
 */

#include <iostream>
#include <cstdlib>      // Get commands // system("pause"); // and // system("cls"); //
#include <fstream>      // Read .txt files
#include <sstream>      // String manipulation

#include <filesystem>   // Utilized to find repository in use. NOTE: Edit .cpp file name alongside (l. 626) to prevent errors
#include <sys/stat.h>   // Utilized for file search function
namespace fs = std::filesystem;

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

class regression {

private:
    double xAxis[1000]{}, yAxis[1000]{};    // Max numer af datapoints set to 1000

    // Working out values in regression analysis //
    double sumX = 0;
    double sumY = 0;
    double sumXX = 0;
    double sumXY = 0;

    double SSR = 0;                         // Sum Squared Regression
    double SST = 0;                         // Total Sum of Squares

    const double e = 2.718281828459045;     // Euler's value approximated 

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

    std::string RRModel[5];                             // Obtained like so => { linReg(), expReg(), logReg(), powReg(), polyReg() };
    double RR[5] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };    // Obtained like so => { linRR(), expRR(), logRR(), powRR(), polyRR() };

    /************************************************************
    *   The constructor does the following:                     *
    *   1. Import datapoints from .txt file                     *
    *       a. Get // xAxis[] // and // xAxis[] //              *
    *       b. Get number of datapoints as // counter //        *
    *   2. Calculate regression models as // models() //        *
    *   3. Calculate respective R^2-values as // RR_models() // *
    *   4. Set datapoints for plotting purposes                 *
    ************************************************************/
    regression(std::string filename) {
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
                        for (int i = 0; i < x.length(); i++) {      // Comma is converted to a period . in x-axis
                            if (x[i] == ',') {
                                x[i] = '.';
                            }
                        }
                        if (x[x.length() - 1] == ',') {         // Commas in between x- and y-values are deleted
                            x.pop_back();
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
                                 "... Consider changing that under .cpp (l. 64) & (l. 165) " << "\n";
                    break;                              // Prevent data overflow. Make sure no more than the set datapoint limit is reached
                }
                std::cout << "Datapoint nr." << counter << "\t[" << x << ";" << y << "]\n";
            }
            std::cout << "Dataset length: " << counter << "\n";

            if (counter < 3) {
                std::cout << "ERROR: Dataset must have at least 3 points" << "\n";
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
        for (int i = 0; i < counter; i++) {         // Set dataset up for plotting
            x.push_back(xAxis[i]);
            y.push_back(yAxis[i]);
        }
        for (int i = 0; i < counter; i++) {         // Set linear regression model up for plotting
            x_lin.push_back(xAxis[i]);
            y_lin.push_back(linModelx(xAxis[i], aLin, bLin));
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
    }


    void reset() { // Working out values for regression analysis purposes
        sumX = 0;
        sumY = 0;
        sumXY = 0;
        sumXX = 0;

        SSR = 0;
        SST = 0;
    }

    void linReg() {
        reset();
        for (int i = 0; i < counter; i++) {
            sumX = sumX + xAxis[i];
            sumXX = sumXX + xAxis[i] * xAxis[i];
            sumY = sumY + yAxis[i];
            sumXY = sumXY + xAxis[i] * yAxis[i];
        }
        aLin = (counter * sumXY - sumX * sumY) / (counter * sumXX - sumX * sumX);
        bLin = (sumXX * sumY - sumX * sumXY) / (counter * sumXX - sumX * sumX);

        RRModel[0] = "y = " + std::to_string(aLin) + "x + " + std::to_string(bLin);
    }

    void expReg() {
        reset();
        for (int i = 0; i < counter; i++) {
            sumX = sumX + xAxis[i];
            sumXX = sumXX + xAxis[i] * xAxis[i];
            sumY = sumY + log(yAxis[i]);
            sumXY = sumXY + xAxis[i] * log(yAxis[i]);
        }
        aExp = (counter * sumXY - sumX * sumY) / (counter * sumXX - sumX * sumX);
        double bExp = (sumXX * sumY - sumX * sumXY) / (counter * sumXX - sumX * sumX);
        cExp = pow(e, bExp); // Back conversion from natural logarithm ln [Yep... log() acutally computes the natural logarithm]

        RRModel[1] = "y = " + std::to_string(cExp) + " e ^ " + std::to_string(aExp) + "x";
    }

    void logReg() {
        reset();
        for (int i = 0; i < counter; i++) {
            sumX = sumX + log(xAxis[i]);
            sumXX = sumXX + log(xAxis[i]) * log(xAxis[i]);
            sumY = sumY + yAxis[i];
            sumXY = sumXY + log(xAxis[i]) * yAxis[i];
        }
        aLog = (counter * sumXY - sumX * sumY) / (counter * sumXX - sumX * sumX);
        bLog = (sumXX * sumY - sumX * sumXY) / (counter * sumXX - sumX * sumX);

        if (bLog >= 0) 
            RRModel[2] = "y = " + std::to_string(aLog) + " ln(x) +" + std::to_string(bLog);
        else
            RRModel[2] = "y = " + std::to_string(aLog) + " ln(x) " + std::to_string(bLog);
    }
   
    void powReg() {
        reset();
        for (int i = 0; i < counter; i++) {
            sumX = sumX + log10(xAxis[i]);
            sumXX = sumXX + log10(xAxis[i]) * log10(xAxis[i]);
            sumY = sumY + log10(yAxis[i]);
            sumXY = sumXY + log10(xAxis[i]) * log10(yAxis[i]);
        }
        aPow = (counter * sumXY - sumX * sumY) / (counter * sumXX - sumX * sumX);
        double bPot = (sumXX * sumY - sumX * sumXY) / (counter * sumXX - sumX * sumX);
        cPow = pow(10, bPot); // Back conversion from log10

        RRModel[3] = "y = " + std::to_string(cPow) + " x ^ " + std::to_string(aPow);
    }

    void secondDegreePolyReg() {    // Made thanks to: https://www.bragitoff.com/2015/09/c-program-for-polynomial-fit-least-squares/
        const int n = 2;            // The degree of the polynomium. Adjust alongside (l. 106), (l. 344), (l. 505) and (l. 697)

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

        RRModel[4] = "y = " + std::to_string(a[2]) + "x^2 + " + std::to_string(a[1]) + "x + " + std::to_string(a[0]);
    }

    void models() { // Calculate all (possible) regression models

        linReg();
        secondDegreePolyReg();
        if (expFejl == false)
            expReg();
        if (logFejl == false)
            logReg();
        if (powFejl == false)
            powReg();
    }

    void linRR() {
        reset();

        for (int i = 0; i < counter; i++) { // Use sumY as average y-value
            sumY = sumY + yAxis[i];
        }
        sumY = sumY / counter;

        for (int i = 0; i < counter; i++) {
            SSR = SSR + pow(yAxis[i] - linModelx(xAxis[i], aLin, bLin), 2);
            SST = SST + pow(yAxis[i] - sumY, 2);
        }

        RR[0] = 1 - SSR / SST;
    }

    void expRR() {
        reset();

        for (int i = 0; i < counter; i++) { // Use sumY as average y-value
            sumY = sumY + yAxis[i];
        }
        sumY = sumY / counter;

        for (int i = 0; i < counter; i++) {
            SSR = SSR + pow(yAxis[i] - expModel(xAxis[i]), 2);
            SST = SST + pow(yAxis[i] - sumY, 2);
        }

        RR[1] = 1 - SSR / SST;
    }

    void logRR() {
        reset();

        for (int i = 0; i < counter; i++) { // Use sumY as average y-value
            sumY = sumY + yAxis[i];
        }
        sumY = sumY / counter;

        for (int i = 0; i < counter; i++) {
            SSR = SSR + pow(yAxis[i] - logModel(xAxis[i]), 2);
            SST = SST + pow(yAxis[i] - sumY, 2);
        }

        RR[2] = 1 - SSR / SST;
    }

    void powRR() {
        reset();

        for (int i = 0; i < counter; i++) { // Use sumY as average y-value
            sumY = sumY + yAxis[i];
        }
        sumY = sumY / counter;

        for (int i = 0; i < counter; i++) {
            SSR = SSR + pow(yAxis[i] - potModel(xAxis[i]), 2);
            SST = SST + pow(yAxis[i] - sumY, 2);
        }

        RR[3] = 1 - SSR / SST;
    }

    void polyRR() {
        reset();

        for (int i = 0; i < counter; i++) { // Use sumY as average y-value
            sumY = sumY + yAxis[i];
        }
        sumY = sumY / counter;

        for (int i = 0; i < counter; i++) {
            SSR = SSR + pow(yAxis[i] - polyModel(xAxis[i]), 2);
            SST = SST + pow(yAxis[i] - sumY, 2);
        }

        RR[4] = 1 - SSR / SST;
    }

    void RR_selectionSort(int n) {
        int max_index;

        for (int i = 0; i < n - 1; i++) {
            max_index = i;
            for (int j = i + 1; j < n; j++) {
                if (RR[j] > RR[max_index])
                    max_index = j;
            }

            // Swap the biggest element with the i-value //
            if (max_index != i) {
                std::swap(RR[max_index], RR[i]);
                std::swap(RRName[max_index], RRName[i]);
                std::swap(RRModel[max_index], RRModel[i]);
            }
        }
    }

    void RR_models() {          // Calculate all (possible) correlation coefficients R^2
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

    void RR_Ranking(int n) {
        std::cout << "REGRESSIONSMODEL RANKING:" << "\n";
        std::cout << "0. The second-degree polynomium\t\t" << RRModel[4] << "\tR^2 = " << RR[4] << "\n\n";

        for (int i = 0; i < n; i++) {
            if (RR[i] != 0.0f)
                std::cout << i + 1 << ". The " << RRName[i] << " regression model\t" << RRModel[i] << "\tR^2 = " << RR[i] << "\n";
        }

        // Check if certain regression models are - can be - used //
        if (expFejl == true)
            std::cout << "WARNING: The exponential regressionsmodel cannot be made (y-axis contains non-postive value)\n";
        if (logFejl == true)
            std::cout << "WARNING: The logarithmic regressionsmodel cannot be made (x-axis contains non-postive value)\n";
        if (powFejl == true)
            std::cout << "WARNING: The power regressionsmodel cannot be made (x- or y-axis contains non-postive value)\n";
    }

    double linModelx(double xValue, double a, double b) { // NOTE: The function accounts for 2-part linear regression
        return a * xValue + b;
    }

    double expModel(double xValue) {
        return cExp * pow(e, aExp * xValue);
    }

    double potModel(double xValue) {
        return cPow * pow(xValue, aPow);
    }

    double logModel(double xValue) {
        return aLog * log(xValue) + bLog;
    }

    double polyModel(double xValue) {
        return a[2] * pow(xValue, 2) + a[1] * xValue + a[0];
    }

    void DMnVM() {   // Print domain and range
        double DM_min = xAxis[0];
        double DM_max = xAxis[counter - 1];
        double VM_min = yAxis[0];
        double VM_max = yAxis[counter - 1];
        std::cout << "Domain" << std::setw(19) << "DM(" << std::to_string(DM_min) << "; " << std::to_string(DM_max) <<
                ")\nRange" << std::setw(20) << "VM(" << std::to_string(VM_min) << "; " << std::to_string(VM_max) << ")\n";
    }

    double xAxisDistance() {                            // The greates distance between datapoints in the x-axis is found
        double xDistance = xAxis[1] - xAxis[0];
        for (int i = 0; i < counter - 1; i++) {
            if (xDistance < xAxis[i + 1] - xAxis[i])    // Account for positive x growth
                xDistance = xAxis[i + 1] - xAxis[i];
            if (xDistance < xAxis[i] - xAxis[i + 1])    // Account for negative x growth
                xDistance = xAxis[i] - xAxis[i + 1];
        }
        return xDistance;
    }

    double yAxisDistance() {                            // The greates distance between datapoints in the y-axis is found
        double yDistance = yAxis[1] - yAxis[0];
        for (int i = 0; i < counter - 1; i++) {
            if (yDistance < yAxis[i + 1] - yAxis[i])    // Account for positive y growth
                yDistance = yAxis[i + 1] - yAxis[i];
            if (yDistance < yAxis[i] - yAxis[i + 1])    // Account for negative y growth
                yDistance = yAxis[i] - yAxis[i + 1];
        }
        return yDistance;
    }

    int outliers() { // Find the number of datapoints that lie below 2/3 * Q1 and above 1.5 * Q3
        int outlier_count = 0;
        double Q1 = yAxis[int(counter / 4)];
        double Q3 = yAxis[int(counter * 3 / 4)];
        for (int i = 0; i < int(counter / 4); i++)
            if (yAxis[i] * 2 / 3 < Q1)
                outlier_count++;
        for (int i = int(counter * 3 / 4); i < counter; i++)
            if (yAxis[i] * 1.5 < Q3)
                outlier_count++;

        return outlier_count;
    }

    std::string linReg2x() { // find 2-part linear regression model and R^2
        if (counter < 6) {
            std::cout << "ERROR: Dataset must have at least 6 points to execute 2-part linear regression" << "\n";
            exit(0);
        }
        
        double aLin1;
        double bLin1;
        double aLin2;
        double bLin2;

        int antal_halv = counter / 2;
        // Linear regression model of the first half of the dataset //
        reset();
        for (int i = 0; i < antal_halv; i++) {
            sumX = sumX + xAxis[i];
            sumXX = sumXX + xAxis[i] * xAxis[i];
            sumY = sumY + yAxis[i];
            sumXY = sumXY + xAxis[i] * yAxis[i];
        }
        aLin1 = (antal_halv * sumXY - sumX * sumY) / (antal_halv * sumXX - sumX * sumX);
        bLin1 = (sumXX * sumY - sumX * sumXY) / (antal_halv * sumXX - sumX * sumX);

        sumY = sumY / antal_halv;

        for (int i = 0; i < antal_halv; i++) {
            SSR = SSR + pow(yAxis[i] - linModelx(xAxis[i], aLin1, bLin1), 2);
            SST = SST + pow(yAxis[i] - sumY, 2);
        }
        double RR_1 = 1 - SSR / SST;

        // Linear regression model of the second half of the dataset //
        reset();
        for (int i = antal_halv; i < counter; i++) {
            sumX = sumX + xAxis[i];
            sumXX = sumXX + xAxis[i] * xAxis[i];
            sumY = sumY + yAxis[i];
            sumXY = sumXY + xAxis[i] * yAxis[i];
        }
        aLin2 = (round((counter + 1) / 2) * sumXY - sumX * sumY) / (round((counter + 1) / 2) * sumXX - sumX * sumX);
        bLin2 = (sumXX * sumY - sumX * sumXY) / (round((counter + 1) / 2) * sumXX - sumX * sumX);

        sumY = sumY / round((counter + 1) / 2);
        for (int i = antal_halv; i < counter; i++) {
            SSR = SSR + pow(yAxis[i] - linModelx(xAxis[i], aLin2, bLin2), 2);
            SST = SST + pow(yAxis[i] - sumY, 2);
        }
        double RR_2 = 1 - SSR / SST;

        return "1) y_1 = " + std::to_string(aLin1) + "x + " + std::to_string(bLin1) + "\tR^2 = " + std::to_string(RR_1) +
            "\n2) y_2 = " + std::to_string(aLin2) + "x + " + std::to_string(bLin2) + "\tR^2 = " + std::to_string(RR_2);
    }

    void plot(std::string lin, std::string exp, std::string log, std::string pot, std::string poly) {

        plt::plot(x, y, "bD");
        plt::plot(x_poly, y_poly, { {"color", "#646464"}, {"ls", ":"}, {"label", poly} });
        plt::plot(x_lin, y_lin, { {"color", "#FF0000"}, {"ls", ":"}, {"label", lin} });
        if (expFejl == false)
            plt::plot(x_exp, y_exp, { {"color", "#00FF00"}, {"ls", ":"}, {"label", exp} });
        if (logFejl == false)
            plt::plot(x_log, y_log, { {"color", "#FF00FF"}, {"ls", ":"}, {"label", log} });
        if (powFejl == false)
            plt::plot(x_pow, y_pow, { {"color", "#00FFFF"}, {"ls", ":"}, {"label", pot} });
    }
};

int main()
{
    std::string filename;
    struct stat sb; // (l. 55)

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

            std::cout << "\nSTATISTICAL DESCRIPTORS:\n";
            std::cout << "Number of datapoints: " << reg.counter << "\n";
            reg.DMnVM(); // Print domain and range
            std::cout << "Greatest x-distance:  " << reg.xAxisDistance() << "\n";
            std::cout << "Greatest y-distance:  " << reg.yAxisDistance() << "\n";
            std::cout << "Outlier count: " << std::setw(8) << reg.outliers() << "\n";

            std::cout << "\nTYPE: \n\"1\" to use linear regression model\n";
            std::cout << "\"2\" to use exponential regression model\n";
            std::cout << "\"3\" to use logarithmic regression model\n";
            std::cout << "\"4\" to use power regression model\n";
            std::cout << "\"5\" to use second-degree polynomial regression model\n";
            std::cout << "\"6\" to use 2-part linear regression model\n";
            std::cout << "\"7\" to see graph of regressionsmodels\n";
            std::cout << "\"9\" to select new dataset\n";
            std::cout << "\"0\" to end program\n";
            std::cin >> input;
            switch (input) {
            case 9:
                reg.~regression();      // User-defined destructor to release memory when selecting new dataset
                system("cls");      // Clear screen
                break;
            case 7: //https://matplotlib-cpp.readthedocs.io/en/latest/docs.html
                reg.plot(std::string("Lin: ") + reg.RRModel[0] + std::string(", R^2 = ") + std::to_string(reg.RR[0]),
                    std::string("Exp: ") + reg.RRModel[1] + std::string(", R^2 = ") + std::to_string(reg.RR[1]),
                    std::string("Log: ") + reg.RRModel[2] + std::string(", R^2 = ") + std::to_string(reg.RR[2]),
                    std::string("Pot: ") + reg.RRModel[3] + std::string(", R^2 = ") + std::to_string(reg.RR[3]),
                    std::string("Poly: ") + reg.RRModel[4] + std::string(", R^2 = ") + std::to_string(reg.RR[4]));
                plt::xlabel("x");
                plt::ylabel("y");
                plt::title(filename);
                plt::legend();
                plt::show();
                system("pause");    // Press any key to continue...
                system("cls");      // Clear screen
                break;
            case 6:
                std::cout << reg.linReg2x() << "\n";
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
                std::cout << reg.aLog << " * " << log(xValue) << " + " << reg.bLog << " = " << reg.logModel(xValue) << "\n";
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
                std::cout << reg.aLin << " * " << xValue << " + " << reg.bLin << " = " << reg.linModelx(xValue, reg.aLin, reg.bLin) << "\n";
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
