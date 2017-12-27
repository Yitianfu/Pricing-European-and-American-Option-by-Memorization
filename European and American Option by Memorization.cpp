//
//  main.cpp
//   Reference from Class Code IE523
//
//  Created by Yitian Fu on 11/28/17.
//  Copyright © 2017 Yitian Fu. All rights reserved.


#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "normdist.h"
#include "newmat.h"
#include "newmatap.h"

using namespace std;

float up_factor, uptick_prob, downtick_prob, notick_prob, risk_free_rate, strike_price;
float initial_stock_price, expiration_time, volatility, R;
int no_of_divisions;
std::vector<std::vector<double>>eurocall;
std::vector<std::vector<double>>europut;
std::vector<std::vector<double>>americacall;
std::vector<std::vector<double>>americaput;

float max(float a, float b) {
    return (b < a )? a:b;
}
double option_price_put_black_scholes(const double& S,      // spot price
                                      const double& K,      // Strike (exercise) price,
                                      const double& r,      // interest rate
                                      const double& sigma,  // volatility
                                      const double& time){
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    return K*exp(-r*time)*N(-d2) - S*N(-d1);
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
                                       const double& K,       // strike (exercise) price,
                                       const double& r,       // interest rate
                                       const double& sigma,   // volatility
                                       const double& time) {  // time to maturity
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt)+0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    return S*N(d1) - K*exp(-r*time)*N(d2);
};

double N(const double& z) {
    if (z > 6.0) { return 1.0; }; // this guards against overflow
    if (z < -6.0) { return 0.0; };
    double b1 = 0.31938153;
    double b2 = -0.356563782;
    double b3 = 1.781477937;
    double b4 = -1.821255978;
    double b5 = 1.330274429;
    double p = 0.2316419;
    double c2 = 0.3989423;
    double a=fabs(z);
    double t = 1.0/(1.0+a*p);
    double b = c2*exp((-z)*(z/2.0));
    double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
    n = 1.0-b*n;
    if ( z < 0.0 ) n = 1.0 - n;
    return n;
};
//column: k=no of divisions, row: i=2*no of divisions + 1
// Initialize memorization matrix to store value of recursion european_call_option
void Initialize(int no_of_divisions)
{
    for(int j=0;j<no_of_divisions+1;j++)
    {
        std::vector<double>temp;
        eurocall.push_back(temp);
        europut.push_back(temp);
        americacall.push_back(temp);
        americaput.push_back(temp);
    }
    for(int j=0;j<no_of_divisions+1;j++)
    {
        for(int i=0;i<(2*no_of_divisions+1);i++)
        {
            eurocall[j].push_back('a');
            europut[j].push_back('a');
            americacall[j].push_back('a');
            americaput[j].push_back('a');
        }
    }
}




float european_call_option(int k, int i)
{
    
    if (eurocall[k][i+no_of_divisions]!='a')
    {
        return eurocall[k][i+no_of_divisions];
    }
    else
    {
       if (k == no_of_divisions)
       {
           eurocall[k][i+no_of_divisions]=max(0.0, (initial_stock_price*pow(up_factor, ((float) i))) - strike_price);
           return eurocall[k][i+no_of_divisions];
       }
       else
        {
            eurocall[k][i+no_of_divisions]=((uptick_prob*european_call_option(k+1,i+1) +
                            notick_prob*european_call_option(k+1, i)+
                            downtick_prob*european_call_option(k+1, i-1))
                           /R);
            return eurocall[k][i+no_of_divisions];
        }
    }
}


float european_put_option(int k, int i)
{
    
    if (europut[k][i+no_of_divisions]!='a')
    {
        return europut[k][i+no_of_divisions];
    }
    else
    {
        if (k == no_of_divisions)
        {
            europut[k][i+no_of_divisions]=max(0.0,  strike_price-(initial_stock_price*pow(up_factor, ((float) i))));
            return europut[k][i+no_of_divisions];
        }
        else
        {
            europut[k][i+no_of_divisions]=((uptick_prob*european_put_option(k+1,i+1) +
                                             notick_prob*european_put_option(k+1, i)+
                                             downtick_prob*european_put_option(k+1, i-1))
                                            /R);
            return europut[k][i+no_of_divisions];
        }
    }
}



float american_call_option(int k, int i, float current_stock_price)
{
    
    if (americacall[k][i+no_of_divisions]!='a')
    {
        return americacall[k][i+no_of_divisions];
    }
    else
    {
        if (k == no_of_divisions)
        {
            americacall[k][i+no_of_divisions]=max(0.0, (current_stock_price - strike_price));
            return americacall[k][i+no_of_divisions];
        }
    else
       {
           americacall[k][i+no_of_divisions]=max((current_stock_price - strike_price),
                                                 (uptick_prob*american_call_option(k+1, i+1,current_stock_price*up_factor)+
                                                  notick_prob*american_call_option(k+1, i, current_stock_price)+
                                                  downtick_prob*american_call_option(k+1, i-1, current_stock_price/up_factor))/R);
            return americacall[k][i+no_of_divisions];
       }
    }
}

float american_put_option(int k, int i, float current_stock_price)
{
    if (americaput[k][i+no_of_divisions]!='a')
    {
        return americaput[k][i+no_of_divisions];
    }
    else
    {
       if (k == no_of_divisions)
       {
           americaput[k][i+no_of_divisions]=max(0.0, (strike_price - current_stock_price));
           return americaput[k][i+no_of_divisions];
       }
    else
       {
           americaput[k][i+no_of_divisions]=max((strike_price - current_stock_price),
                                                (uptick_prob*american_put_option(k+1, i+1, current_stock_price*up_factor) +
                                                 notick_prob*american_put_option(k+1, i, current_stock_price)+
                                                 downtick_prob*american_put_option(k+1, i-1, current_stock_price/up_factor))/R);
           return americaput[k][i+no_of_divisions];
       }
    }
}






int main (int argc, char* argv[])
{
    
    sscanf (argv[1], "%f", &expiration_time);
    sscanf (argv[2], "%d", &no_of_divisions);
    sscanf (argv[3], "%f", &risk_free_rate);
    sscanf (argv[4], "%f", &volatility);
    sscanf (argv[5], "%f", &initial_stock_price);
    sscanf (argv[6], "%f", &strike_price);
    
    up_factor = exp(volatility*sqrt(2*(expiration_time/((float) no_of_divisions))));
    R = exp(risk_free_rate*expiration_time/((float) no_of_divisions));
    uptick_prob = pow((sqrt(R) - (1/sqrt(up_factor)))/(sqrt(up_factor)-(1/sqrt(up_factor))),2);
    downtick_prob= pow((sqrt(up_factor) - (sqrt(R)))/(sqrt(up_factor)-(1/sqrt(up_factor))),2);
    notick_prob =1-uptick_prob-downtick_prob;
    cout << "Recursive Trinomial European Option Pricing" << endl;
    cout << "Expiration Time (Years) = " << expiration_time << endl;
    cout << "Number of Divisions = " << no_of_divisions << endl;
    cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
    cout << "Initial Stock Price = " << initial_stock_price << endl;
    cout << "Strike Price = " << strike_price << endl;
    cout << "--------------------------------------" << endl;
    cout << "R = "<<R<<endl;
    cout << "Up Factor = " << up_factor << endl;
    cout << "Uptick Probability = " << uptick_prob << endl;
    cout << "Downtick Probability = " << downtick_prob << endl;
    cout << "Notick Probability = " << notick_prob << endl;
    cout << "--------------------------------------" << endl;
    Initialize(no_of_divisions);
   
    double call_price = european_call_option(0, 0);
    cout << "Binomial Price of an European Call Option = " << call_price << endl;
    cout << "Call Price according to Black-Scholes = " <<
    option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate,
                                    volatility, expiration_time) << endl;
    
    double put_price = european_put_option(0, 0);
    cout << "Binomial Price of an European Put Option = " << put_price << endl;
    cout << "Put Price according to Black-Scholes = " <<
    option_price_put_black_scholes(initial_stock_price, strike_price, risk_free_rate,
                                   volatility, expiration_time) << endl;
    
    double ame_call_price = american_call_option(0, 0,initial_stock_price);
    cout << "Trinomial Price of an American Call Option = " << ame_call_price << endl;
    double ame_put_price = american_put_option(0, 0, initial_stock_price);
    cout << "Trinomial Price of an American Put Option = " << ame_put_price << endl;
   
    
    
    
}



