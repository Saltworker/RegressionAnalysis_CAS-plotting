/************************************************************
*   Regressionsanalyse v�rkt�j som kan udregne:             *
*   - Line�r, eksponentiel og potentiel regressionsmodel    *
*   - Korrelationskoefficient af regressionsmodeller        *
*   - Outliers inkl. st�rste afstand mellem x- og y-punkter *
*   - x-v�rdi indsat i en af regressionsmodellerne          *
*************************************************************
*   Lavet af: Mathias Lykholt-Ustrup HTXr20                 *
*   Opdateret: 14. juni 2023                                *
************************************************************/

#include <iostream>
#include <cstdlib> //Benyttes kun til ***system("pause");*** og ***system("cls");***
#include <fstream> //Benyttes til at afl�se filer
#include <sstream> //Benyttes til at afl�se strings

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

using namespace std;

class regression {

public:
    double xAkse[999]{}, yAkse[999]{}; //Max 1000 punkter i datas�t

    int antal = 0; // punkter

    //Mellemregninger i regressionsanalyse
    double sumX = 0;
    double sumY = 0;
    double sumXX = 0;
    double sumXY = 0;

    double SSR = 0; //Sum Squared Regression
    double SST = 0; //Total Sum of Squares

    double aLin = 0; //H�ldningskoefficient af line�r regressionsmodel
    double bLin = 0; //Startv�rdi af line�r regressionsmodel

    double aExp = 0; //Startv�rdi af eksponentiel regressionsmodel
    double cExp = 0; //Fremskrivningsfaktor af eksponentiel regressionsmodel

    double aPot = 0; //Startv�rdi af potentiel regressionsmodel
    double cPot = 0; //Eksponent af potentiel regressionsmodel

    const double e = 2.718281828459045; // Tiln�rmet v�rdi for eulers tal

    bool expFejl = false; //Vil der opst� fejl ved at finde eksponentiel regressionsmodel?
    bool potFejl = false; //Vil der opst� fejl ved at finde potentiel regressionsmodel?

    string filnavn;

    regression() : regression("dataset.txt") {}

    regression(string filnavn) : filnavn(filnavn) { //Afl�s punkter fra .txt fil
        ifstream dataset; //Initialis�r ifstream (kun afl�sning af fil)
        dataset.open(filnavn); //�bn datas�t fil

        if (dataset.is_open()) { //�bnes filen? Eksisterer filen i mappe med .cpp fil?
            string point; //string til at behandle given linje

            while (getline(dataset, point)) { //imens der findes en string v�rdi i dataset:
                stringstream s(point); //Start stringstream og seper�r afh�ngigt af mellemrum
                int P = 0; //x n�r P=0 og y n�r P=1
                string temp, x, y; //
                while (s >> temp) { //Mens der afl�ses stringv�rdier i "antal"'te linje
                    if (P == 0) { //x-v�rdi afl�ses
                        x = temp;
                        for (int i = 0; i < x.length(); i++) { // Kommaseperatorer omskrives til punktum i x-akse
                            if (x[i] == ',') {
                                x[i] = '.';
                            }
                        }
                        if (x[x.length() - 1] == ',') { //Fjern eksisterende komma imellem x- og y-v�rdi
                            x.pop_back();
                        }
                        xAkse[antal] = stof(x); //Konvertering fra string til float af x-v�rdi
                        if (xAkse[antal] <= 0) { //Den potentielle regressionsmodel kan ikke benyttes med x=0 eller negative x-v�rdier
                            potFejl = true;
                        }
                    }
                    else { //y-v�rdi afl�ses
                        y = temp;
                        for (int i = 0; i < y.length(); i++) { // Kommaseperator omskrives til punktum i y-akse
                            if (y[i] == ',') {
                                y[i] = '.';
                            }
                        }
                        yAkse[antal] = stof(y); //Konvertering fra string til float af y-v�rdi
                        if (yAkse[antal] <= 0) { //Den potentielle og eksponentielle regressionsmodel kan ikke benyttes med y=0 eller negative y-v�rdier
                            expFejl = true;
                            potFejl = true;
                        }
                    }
                    P++;
                }
                antal++;
                if (antal > 1000) { // Gr�nsen for hvor mange punkter som afl�ses s�ttes
                    cout << "Advarsel: Max antal punkter sat til " << antal << endl; // Hvor "antal" er afh�ngig af linje 14
                    break;
                }
                cout << "Punkt nr." << antal << "\t[" << x << ";" << y << "]" << endl;
            }
            cout << "Laengde af datasaet: " << antal << endl;

            if (antal < 3) {
                cout << "Fejlkode: Mindst 3 punkter i datasaet" << endl;
                exit(0);
            }

            dataset.close();
        }
        else {
            cout << "Fejlkode: Fil kunne ikke aflaeses" << endl;
            exit(0);
        }
    }


    void reset() { //Nulstil regressionsanalyse variabler
        sumX = 0;
        sumY = 0;
        sumXY = 0;
        sumXX = 0;

        SSR = 0;
        SST = 0;
    }

    string linReg() { //Udregner line�r regressionsmodel og returnerer den i string
        reset();
        for (int i = 0; i < antal; i++) {
            sumX = sumX + xAkse[i];
            sumXX = sumXX + xAkse[i] * xAkse[i];
            sumY = sumY + yAkse[i];
            sumXY = sumXY + xAkse[i] * yAkse[i];
        }
        aLin = (antal * sumXY - sumX * sumY) / (antal * sumXX - sumX * sumX);
        bLin = (sumXX * sumY - sumX * sumXY) / (antal * sumXX - sumX * sumX);

        return "y = " + to_string(aLin) + "x + " + to_string(bLin);
    }

    double linModelx(double xValue, double a, double b) { // Funktionen benyttes til flere line�re regressionsmodeller
        return a * xValue + b;
    }

    double LinRR() { // Den bedre metode til at finde R^2 for line�r regressionsmodel
        reset();

        //Find gennemsnitsv�rdi for y-akse og y-v�rdi(er) funet ud fra regressionsmodel
        for (int i = 0; i < antal; i++) {
            sumY = sumY + yAkse[i];
        }
        sumY = sumY / antal;

        for (int i = 0; i < antal; i++) {
            SSR = SSR + pow(yAkse[i] - linModelx(xAkse[i], aLin, bLin), 2);
            SST = SST + pow(yAkse[i] - sumY, 2);
        }

        return 1 - SSR / SST;
    }

    string expReg() { //Udregner eksponentiel regressionsmodel og returnerer den i string
        reset();
        for (int i = 0; i < antal; i++) {
            sumX = sumX + xAkse[i];
            sumXX = sumXX + xAkse[i] * xAkse[i];
            sumY = sumY + log(yAkse[i]);
            sumXY = sumXY + xAkse[i] * log(yAkse[i]);
        }
        aExp = (antal * sumXY - sumX * sumY) / (antal * sumXX - sumX * sumX);
        double bExp = (sumXX * sumY - sumX * sumXY) / (antal * sumXX - sumX * sumX);
        cExp = pow(e, bExp); // Tilbagekonverting fra (naturlig) logaritme

        return "y = " + to_string(cExp) + " e^" + to_string(aExp) + "x";
    }

    double expModel(double xValue) {
        return cExp * pow(e, aExp * xValue);
    }

    double ExpRR() {
        reset();

        //Find gennemsnitsv�rdi for y-akse og y-v�rdi(er) fundet ud fra regressionsmodel
        for (int i = 0; i < antal; i++) {
            sumY = sumY + yAkse[i];
        }
        sumY = sumY / antal;

        for (int i = 0; i < antal; i++) {
            SSR = SSR + pow(yAkse[i] - expModel(xAkse[i]), 2);
            SST = SST + pow(yAkse[i] - sumY, 2);
        }

        return 1 - SSR / SST;
    }

    string potReg() { //Udregner potentiel regressionsmodel og returnerer den i string
        reset();
        for (int i = 0; i < antal; i++) {
            sumX = sumX + log10(xAkse[i]);
            sumXX = sumXX + log10(xAkse[i]) * log10(xAkse[i]);
            sumY = sumY + log10(yAkse[i]);
            sumXY = sumXY + log10(xAkse[i]) * log10(yAkse[i]);
        }
        aPot = (antal * sumXY - sumX * sumY) / (antal * sumXX - sumX * sumX);
        double bPot = (sumXX * sumY - sumX * sumXY) / (antal * sumXX - sumX * sumX);
        cPot = pow(10, bPot); // Tilbagekonertering fra log10

        return "y = " + to_string(cPot) + " x^" + to_string(aPot);
    }

    double potModel(double xValue) {
        return cPot * pow(xValue, aPot);
    }

    double PotRR() {
        reset();

        //Find gennemsnitsv�rdi for y-akse og y-v�rdi(er) fundet ud fra regressionsmodel
        for (int i = 0; i < antal; i++) {
            sumY = sumY + yAkse[i];
        }
        sumY = sumY / antal;

        for (int i = 0; i < antal; i++) {
            SSR = SSR + pow(yAkse[i] - potModel(xAkse[i]), 2);
            SST = SST + pow(yAkse[i] - sumY, 2);
        }

        return 1 - SSR / SST;
    }

    string DMnVM() { //Print definitions- og v�rdim�ngde
        double DM_min = xAkse[0];
        double DM_max = xAkse[antal - 1];
        double VM_min = yAkse[0];
        double VM_max = yAkse[antal - 1];
        return "Definitionsmaengden lyder DM(" + to_string(DM_min) + ";" + to_string(DM_max) +
            ") og vaerdimaengden VM(" + to_string(VM_min) + ";" + to_string(VM_max) + ")";
    }

    double xAkseAfstand() { //Den st�rste afstand mellem punkter ift. x-aksen findes.
        double xAfstand = xAkse[1] - xAkse[0];
        for (int i = 0; i < antal - 1; i++) {
            if (xAfstand < xAkse[i + 1] - xAkse[i]) { //Positiv x tilv�kst tages h�jde for
                xAfstand = xAkse[i + 1] - xAkse[i];
            }
            if (xAfstand < xAkse[i] - xAkse[i + 1]) { //Negativ x tilv�kst tages h�jde for
                xAfstand = xAkse[i] - xAkse[i + 1];
            }
        }

        return xAfstand;
    }

    double yAkseAfstand() { //Den st�rste afstand mellem punkter ift. y-aksen findes.
        double yAfstand = yAkse[1] - yAkse[0];
        for (int i = 0; i < antal - 1; i++) {
            if (yAfstand < yAkse[i + 1] - yAkse[i]) { //Postitiv y tilv�kst tages h�jde for
                yAfstand = yAkse[i + 1] - yAkse[i];
            }
            if (yAfstand < yAkse[i] - yAkse[i + 1]) { //Negativ y tilv�kst tages h�jde for
                yAfstand = yAkse[i] - yAkse[i + 1];
            }
        }

        return yAfstand;
    }

    int outliers() { //Finder antal punkter som ligger under 2/3 * Q1 og over 1.5 * Q3
        int outlier_antal = 0;
        double Q1 = yAkse[int(antal / 4)];
        double Q3 = yAkse[int(antal * 3 / 4)];
        for (int i = 0; i < int(antal / 4); i++) {
            if (yAkse[i] * 2 / 3 < Q1) {
                outlier_antal++;
            }
        }
        for (int i = int(antal * 3 / 4); i < antal; i++) {
            if (yAkse[i] * 1.5 < Q3) {
                outlier_antal++;
            }
        }

        return outlier_antal;
    }

    string linReg2x() { //Udregner line�r regressionsmodel og korrelationskoefficient R^2 for hver den ene og anden halvdel af datas�ttet
        double aLin1;
        double bLin1;
        double aLin2;
        double bLin2;

        int antal_halv = antal / 2;
        //Line�r regressionsmodel for f�rste halvdel af datas�t
        reset();
        for (int i = 0; i < antal_halv; i++) {
            sumX = sumX + xAkse[i];
            sumXX = sumXX + xAkse[i] * xAkse[i];
            sumY = sumY + yAkse[i];
            sumXY = sumXY + xAkse[i] * yAkse[i];
        }
        aLin1 = (antal_halv * sumXY - sumX * sumY) / (antal_halv * sumXX - sumX * sumX);
        bLin1 = (sumXX * sumY - sumX * sumXY) / (antal_halv * sumXX - sumX * sumX);

        sumY = sumY / antal_halv;

        for (int i = 0; i < antal_halv; i++) {
            SSR = SSR + pow(yAkse[i] - linModelx(xAkse[i], aLin1, bLin1), 2);
            SST = SST + pow(yAkse[i] - sumY, 2);
        }
        double RR_1 = 1 - SSR / SST;

        //Line�r regressionsmodel for sidste halvdel af datas�t
        reset();
        for (int i = antal_halv; i < antal; i++) {
            sumX = sumX + xAkse[i];
            sumXX = sumXX + xAkse[i] * xAkse[i];
            sumY = sumY + yAkse[i];
            sumXY = sumXY + xAkse[i] * yAkse[i];
        }
        aLin2 = (round((antal + 1) / 2) * sumXY - sumX * sumY) / (round((antal + 1) / 2) * sumXX - sumX * sumX);
        bLin2 = (sumXX * sumY - sumX * sumXY) / (round((antal + 1) / 2) * sumXX - sumX * sumX);

        sumY = sumY / round((antal + 1) / 2);
        for (int i = antal_halv; i < antal; i++) {
            SSR = SSR + pow(yAkse[i] - linModelx(xAkse[i], aLin2, bLin2), 2);
            SST = SST + pow(yAkse[i] - sumY, 2);
        }
        double RR_2 = 1 - SSR / SST;

        return "1) y_1 = " + to_string(aLin1) + "x + " + to_string(bLin1) + "\tR^2 = " + to_string(RR_1) +
            "\n2) y_2 = " + to_string(aLin2) + "x + " + to_string(bLin2) + "\tR^2 = " + to_string(RR_2);
    }

    void plot(string lin, string exp, string pot) {
        std::vector<double> x, y, x_lin, y_lin, x_exp, y_exp, x_pot, y_pot;

        x.clear();
        y.clear();
        x_lin.clear();
        y_lin.clear();
        x_exp.clear();
        y_exp.clear();
        x_pot.clear();
        y_pot.clear();

        for (int i = 0; i < antal; i++) { //datas�t plot
            x.push_back(xAkse[i]);
            y.push_back(yAkse[i]);
        }
        for (int i = 0; i < antal; i++) { //Line�r regressionsmodel plot
            x_lin.push_back(xAkse[i]);
            y_lin.push_back(linModelx(xAkse[i], aLin, bLin));
        }
        for (int i = 0; i < antal; i++) { //Eksponentiel regressionsmodel plot
            x_exp.push_back(xAkse[i]);
            y_exp.push_back(expModel(xAkse[i]));
        }
        for (int i = 0; i < antal; i++) { //Potentiel regressionsmodel plot
            x_pot.push_back(xAkse[i]);
            y_pot.push_back(potModel(xAkse[i]));
        }

        plt::plot(x, y, "bD");
        plt::plot(x_lin, y_lin, { {"color", "#000064"}, {"ls", ":"}, {"label", lin} });
        if (expFejl == false) {
            plt::plot(x_exp, y_exp, { {"color", "#640064"}, {"ls", ":"}, {"label", exp} });
        }
        if (potFejl == false) {
            plt::plot(x_pot, y_pot, { {"color", "#646464"}, {"ls", ":"}, {"label", pot} });
        }
    }
};

int main()
{    
    string lineaer, eksponentiel, potentiel;
    int rang, input;
    double RR_lin, RR_exp = 0.0f, RR_pot = 0.0f, xV�rdi; //V�rdier som evt. ikke benyttes foruddefineres

    regression reg("dataset_geogebra4.txt"); //reg("dataset_geogebra4.txt");
    cout << "Datasaet er indlaest!" << endl;
    system("pause"); //Press any key to continue...
    system("cls"); //Ryd sk�rm for indl�st datas�t

    lineaer = reg.linReg();
    RR_lin = reg.LinRR();

    if (reg.expFejl == false) {
        eksponentiel = reg.expReg();
        RR_exp = reg.ExpRR();
    }
    if (reg.potFejl == false) {
        potentiel = reg.potReg();
        RR_pot = reg.PotRR();
    }

    double RR[3]{}; // Hvor RR[0] er for line�r, RR[1] for eksponentiel og RR[2] for potentiel regressionsmodel

    while (1) { // Forever loop
        /* Udregn (hvis muligt) regressionsmodeller samt korrelationskoefficienter */

        RR[0] = RR_lin;
        if (reg.expFejl == false) {
            RR[1] = RR_exp;
        }
        else {
            cout << "Advarsel: Eksponentiel regressionsmodel kan ikke laves (y-vaerdi indebaerer ikke-postiv vaerdi)" << endl;
            RR[1] = 0.0f;
        }

        if (reg.potFejl == false) {
            RR[2] = RR_pot;
        }
        else {
            cout << "Advarsel: Potentiel regressionsmodel kan ikke laves (x- eller y-vaerdi indebaerer ikke-postiv vaerdi)" << endl;
            RR[2] = 0.0f;
        }
        cout << "RANGERING AF REGRESSIONSMODELLER: " << endl;

        for (rang = 1; RR[0] + RR[1] + RR[2] > 0; rang++) { // Mens der ikke er regressionsmodeller som er taget h�jde for. SELECTIONSORT GO BRR
            cout << rang << ". ";
            double max = RR[0];
            for (int i = 0; i <= 2; i++) { // Find st�rste korrelationskoefficient
                if (RR[i] > max) {
                    max = RR[i];
                }
            }
            // Den st�rste korrelationskoefficients regressionsmodel findes, printes og s�ttes lig 0
            if (max == RR[2]) {
                cout << "Den potentielle regressionsmodel: " << potentiel << "\t R^2: " << RR[2] << endl;
                RR[2] = 0.0f;
            }
            else {
                if (max == RR[1]) {
                    cout << "Den eksponentielle regressionsmodel: " << eksponentiel << "\t R^2: " << RR[1] << endl;
                    RR[1] = 0.0f;
                }
                else { // if (max = RR[0] {}
                    cout << "Den lineaere regressionsmodel: " << lineaer << "\t R^2: " << RR[0] << endl;
                    RR[0] = 0.0f;
                }
            }
        }

        cout << "\nKOMMENTARER: " << endl;
        cout << "Antal punkter: " << reg.antal << endl;
        cout << reg.DMnVM() << endl; //Definitions- og v�rdim�ngder printes
        cout << "Stoerste x-afstand: " << reg.xAkseAfstand() << endl;
        cout << "Stoerste y-afstand: " << reg.yAkseAfstand() << endl;
        cout << "Antal outliers: " << reg.outliers() << endl;

        cout << "\nINDTAST: \n\"1\" for at benytte lineaer regressionsmodel" << endl;
        cout << "\"2\" for at benytte eksponentiel regressionsmodel" << endl;
        cout << "\"3\" for at benytte potentiel regressionsmodel" << endl;
        cout << "\"4\" for at se 2-delt lineaer regressionsmodel" << endl;
        cout << "\"0\" for at afslutte program" << endl;
        cin >> input;
        switch (input) {
        case 5: //https://matplotlib-cpp.readthedocs.io/en/latest/docs.html
            reg.plot(string("Lin: ") + lineaer + string(", R^2 = ") + to_string(RR_lin), string("Exp: ") + eksponentiel + string(", R^2 = ") + to_string(RR_exp), string("Pot: ") + potentiel + string(", R^2 = ") + to_string(RR_pot));
            plt::xlabel("x");
            plt::ylabel("y");
            plt::title(reg.filnavn);
            plt::legend();
            plt::show();
            system("pause"); //Press any key to continue...
            system("cls"); //Ryd sk�rm for indl�st datas�t
            break;
        case 4:
            cout << reg.linReg2x() << endl;
            system("pause"); //Press any key to continue...
            system("cls"); //Ryd sk�rm for indl�st datas�t
            break;
        case 3:
            cout << "Potentiel \nIndtast x-vaerdi: " << endl;
            cin >> xV�rdi;
            cout << reg.cPot << " * " << xV�rdi << " ^ " << reg.aPot << " = " << reg.potModel(xV�rdi) << endl;
            system("pause"); //Press any key to continue...
            system("cls"); //Ryd sk�rm for indl�st datas�t
            break;
        case 2:
            cout << "Eksponentiel \nIndtast x-vaerdi: " << endl;
            cin >> xV�rdi;
            cout << reg.cExp << " * " << reg.aExp << " ^ " << xV�rdi << " = " << reg.expModel(xV�rdi) << endl;
            system("pause"); //Press any key to continue...
            system("cls"); //Ryd sk�rm for indl�st datas�t
            break;
        case 1:
            cout << "Lineaer \nIndtast x-vaerdi: " << endl;
            cin >> xV�rdi;
            cout << reg.aLin << " * " << xV�rdi << " + " << reg.bLin << " = " << reg.linModelx(xV�rdi, reg.aLin, reg.bLin) << endl;
            system("pause"); //Press any key to continue...
            system("cls"); //Ryd sk�rm for indl�st datas�t
            break;
        case 0:
            return 0;
        default:
            cout << "Kunne ikke finde case. Proev igen." << endl;
            system("pause"); //Press any key to continue...
            system("cls"); //Ryd sk�rm for indl�st datas�t
        }
    }
};