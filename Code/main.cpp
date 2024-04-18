#include <bits/stdc++.h>
using namespace std;
#define rep(i, j, k) for (int i = j; i < k; i++)
#define db double
#define vi vector<db>

const db R = 8.3144618 * 10;
const db Pc = 220.8766; // In bar
const db Tc = 647.096;  // In Kelvin...
const db omega = 0.224;
const vi ar{-7.8595178, 1.8440825, -11.786649, 22.680741, -15.9618719, 1.8012250};
const vi secor{-0.0652869, 1.6790636 * 1e-4,40.838951, 0, 0, -3.9266518 * 1e-2, 0, 2.1157167 * 1e-2, 6.5486487 * 1e-6, 0};
const vi thor{-1.144624 * 1e-2, 2.8274958 * 1e-5, 0, 0, 0, 1.3980876 * 1e-2, 0, -1.4349005 * 1e-2, 0, 0};
const db sigma = 1 + sqrt(2);
const db epsilon = 1 - sqrt(2);


//----------------------------------------------------------------
db spvolh2o(db temp) // Cheked...
{
    return ((1 + (18.1597 * 1e-3 * temp)) / (0.9998 + (18.2249 * 1e-3 * temp) - (7.9222 * 1e-6 * pow(temp, 2)) - (55.4485 * 1e-9 * pow(temp, 3)) + (149.7562 * 1e-12 * pow(temp, 4)) - (393.2952 * 1e-15 * pow(temp, 5))));
}
db rho(db temp, db press) // cheked ...
{
    db b = 19654.32 + 147.037 * temp - 2.2155 * pow(temp, 2) + 1.0478 * 1e-2 * pow(temp, 3) - 2.2789 * 1e-5 * pow(temp, 4);
    db a1 = 3.2891 - 2.391 * 1e-3 * temp + 2.8446 * 1e-4 * pow(temp, 2) - 2.82 * 1e-6 * pow(temp, 3) + 8.477 * 1e-9 * pow(temp, 4);
    db a2 = 6.245 * 1e-5 - 3.913 * 1e-6 * temp - 3.499 * 1e-8 * pow(temp, 2) + 7.942 * 1e-10 * pow(temp, 3) - 3.299 * 1e-12 * pow(temp, 4);
    db v0 = spvolh2o(temp);
    db revrho = v0 - ((v0 * press) / (b + a1 * press + a2 * press));
    return (1 / revrho);
}
//----------------------------------------------------------------

db psh2o(db temp) // checked ....
{
    db t0 = temp + 273.15;
    db gm = (1 - (t0 / Tc));
    db lterm = (Tc / t0) * ((ar[0] * gm) + (ar[1] * (pow(gm, 1.5))) + (ar[2] * pow(gm, 3)) + (ar[3] * pow(gm, 3.5)) + (ar[4] * pow(gm, 4)) + (ar[5] * pow(gm, 7.5)));
    return (exp(lterm) * Pc);
}
//----------------------------------------------------------------
db fugh2o(db temp, db press)
{
    db ps = psh2o(temp);
    return (ps * exp((18.0152 * (press - ps)) / (rho(temp, press) * R * (temp + 273.15))));
}

db kh2o(db temp, db press)
{
    db fh2o = fugh2o(temp, press);
    db eqconst = -2.209+ (3.097*1e-2*temp) - (1.098*1e-4*pow(temp,2)) + (2.048*1e-7*pow(temp,3));
    eqconst = pow(10,eqconst);
    return ((eqconst / (fh2o * press) * exp(((press - 1) * 18.18) / (R * (temp + 273.15)))));
}

db henconst(db temp, db press) // Checked.....
{
    db delb = -5.279063 + 6.187967 * pow((1e3 / (temp + 273.15)), 0.5);
    db neta = -0.114535;
    return (exp((1 - neta) * (log(fugh2o(temp, press))) + neta * (log((R * (temp + 273.15) * rho(temp, press)) / (18.18))) + (2 * rho(temp, press) * delb)));
}

vi parampure(db temp, db press)
{
    db tcritialco2 = 31+273.15;
    db pcriticalco2 = 73.8;
    db t0 = temp + 273.15;
    db asentricfac = 0.228;
    db m = 0.37464 + (1.54226 * asentricfac) - (0.26992 * asentricfac*asentricfac);
    db b = 0.077796 * (R * (tcritialco2 / pcriticalco2));
    db a = 0.457236 * ((pow(R, 2) * pow(tcritialco2, 2)) / pcriticalco2)* pow((1+m*(1-sqrt(t0/tcritialco2))),2);
    db A = (a * press) / (pow((R * t0), 2));
    db B = (b * press) / (R * t0);
    return {a, b, m, A, B};
}

vi pengrobinson(db z, db temp, db press){
    vi param = parampure(temp,press);
    db eqn = pow(z,3)- ((1-param[4])* pow(z,2)) + (param[3]- 2*param[4]- 3*pow(param[4],2))*z - (param[3]*param[4] - pow(param[4],2) -pow(param[4],3));
    db der = 3*pow(z,2) - 2*(1-param[4])*z + (param[3]-2*param[4]-3*pow(param[4],2));
    return {eqn,der};
}

db solvepengrobinson(db zinit,db temp, db press){
    db z = -1;
    db tol = 1e-5;
    db ep = INT_MAX;
    while(ep>tol){
        vi eqn = pengrobinson(zinit,temp,press);
        z= zinit - (eqn[0]/eqn[1]);
        ep = abs(zinit-z);
        zinit = z;
    }
    return z;
}
db alphapr(db temp) // Checked .....
{
    db tcritialco2 = 31+273.15;
    return pow((1 + (0.37464 + 1.54226 * omega - 0.26992 * pow(omega, 2)) * (1 - sqrt(temp+273.15/ tcritialco2))), 2);
}
vi betaqui(db temp, db press)// Checked .......
{ 
    db tcritialco2 = 31+273.15;
    db pcriticalco2 = 73.8;
    db crtemp = (temp + 273.15) / tcritialco2;
    db beta = omega * ((press / pcriticalco2) / (crtemp));
    db qui = (0.45724 * alphapr(crtemp)) / (omega * (crtemp));
    return {beta, qui};
}

db pengrobison(db temp, db press) // Checked........
{
    vi constant = betaqui(temp, press);
    db ep = 1e-5;
    db zinit2 = 1;
    db zg = -1;
    db tol2 = INT_MAX;
    while (tol2 > ep)
    {
        zg = (1 + constant[0]) - (constant[1] * constant[0] * ((zinit2 - constant[0]) / ((zinit2 + epsilon * constant[0]) * (zinit2 + sigma * constant[1]))));
        tol2 = abs(zg - zinit2);
        zinit2 = zg;
    }
    return zg;
}


db param(db temp, db press, vi c) // Checked................................................................
{
    db t0 = temp + 273.15;
    return (c[0] + (c[1] * t0) + (c[2] / t0) + (c[3] * press) + (c[4] / press) + (c[5] * (press / t0)) + (c[6] * (t0 / pow(press, 2))) + ((c[7] * press) / (630 - t0)) + (c[8] * t0 * log(press)) + (c[9] * (press / pow(t0, 2))));
}

db actcoff(db temp, db press, db molality) // Err .....
{
    db term1 = param(temp, press, secor); 
    db term2 = param(temp, press, thor); 
    db finexp = (2 * molality * term1) + (2 * pow(molality, 2) * term2);
    return exp(finexp);
}

db kc02(db temp, db press, db molality) // Checked.....
{
    db henry = henconst(temp, press);
    db coeff = actcoff(temp, press, molality);
    vi pm = parampure(temp, press);
    db tcritialco2 = 31+273.15;
    db pcriticalco2 = 73.8;
    db beta = 0.07779* ((press/pcriticalco2) / ((temp+273.15)/(tcritialco2)));
    db z = solvepengrobinson(1,temp,press);
    // db z = pengrobinson(temp,press);
    // db z = pengrobison(temp,press);
    // cout<<z<<endl<<endl;
    db phi = (z-1) - log(z - pm[4]) - (pm[3]/(pm[4]*(epsilon-sigma))) *log((z+epsilon*pm[4])/(z+sigma*pm[4]));
    phi = exp(phi);
    // db phi =1;
    return ((henry * coeff) / (press * phi));
}

db yh2o(db temp, db press, db molality) // Checked ......
{
    db co2 = kc02(temp, press, molality);
    db h2o = kh2o(temp, press);
    return ((1 - (1 / co2)) / ((1 / h2o) - (1 / co2)));
}

db xco2(db temp, db press, db molality)
{
    db yh20 = yh2o(temp, press, molality);
    db kco2 = kc02(temp, press, molality);
    db yco2n = (1 / (1 + yh20));
    return (yco2n / kco2);
}


db zchecker(db temp, db press)
{
    vi param = parampure(temp, press);
    db zl = solvepengrobinson(1, temp, press);
    db zg = pengrobison(temp, press);
    db check = (zg-zl) + log((zl-param[4])/(zg-param[4])) - (param[3]/(param[4]*(epsilon-sigma)))*log(((zl+sigma*param[4])*(zg+epsilon*param[4]))/((zl+epsilon*param[4])*(zg+sigma*param[4])));  
    return (check>0)?zl:zg;
}
signed main()
{

    // cout<< betaqui(60,83)[0]<<" "<<betaqui(60,83)[1]<<endl;
    // cout<<alphapr(60)<<endl;
    // vi a = parmapure(60, 83);
    // for(auto it:a){
    //     cout<<it<<" ";
    // }
    // cout<<solvepengrobinson(0.1,60,83)<<endl;
    // cout<<endl;
    // cout<<kc02(60,83,1)<<endl;
    // cout<<kc02(60,83,0.0172)<<endl;
    // cout<<yh2o(60,83,1)<<endl;
    // cout<<solvepengrobinson(1,100,150)<<endl;
    // cout<<xco2(100,150,6)<<endl;
    // cout<<"kh2o:- " <<kh2o(60,83)<<endl;
    // cout<<"yh2o:- " <<yh2o(60,83,0.0172)<<endl;
    // cout<<"fh20:- "<<fugh2o(60,83)<<endl;
    // cout<<"gamma:- "<<actcoff(60,83,0.0172)<<endl;
    // cout<<"henconst:- "<<henconst(60,83)<<endl;
    // cout<<"kco2:- "<<kc02(60,83,0.0172)<<endl;
    cout<<"xco2:- "<<xco2(80,83,0.172)*100<<endl;
    // cout<<solvepengrobinson(1,100,150)<<endl;
    // cout<<
    // cout<<pengrobison(100,150)<<endl;
    // cout<<kc02(60,83,0.0172)<<endl;
    // cout<<psh2o(60)<<endl;
    // cout<< pengrobison(100,150)<<endl;
    // cout<<xco2(100,150,6)*100<<endl;
    // cout<<zchecker(60,83)<<endl;
    // cout<<zchecker(100,150)<<endl;
    // db zg = pengrobison(temp, press);
    // cout<<pengrobison(100,150)<<endl;
    // cout<<kh2o(60,83)<<endl;
    // cout<<kh2o(60,83)<<endl;
    // cout<<kc02(60,83)<<endl;
    // cout<<fugh2o(60,83)<<endl;
    // cout<<kc02(80,103,0.0172)<<endl;
}