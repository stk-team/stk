#include "sphere-subdivision.h"

#include <ctime>
#include "./../core/utils.h"

#define yx M_PI_4
int numpts_to_nsubdivsions(const int &N, std::string domain){
    int nsubdivisions = 1;
    int lb_nsubdivisions = sqrt(N/12.0);                       //lower bound on the nSubdivisions

    if(domain == "hemisphere")
        lb_nsubdivisions = sqrt(N/6.0);                       //lower bound on the nSubdivisions

    int ub_nsubdivisions = lb_nsubdivisions+1;                 //upper bound on the nSubdivisions
    int t1_npts = lb_nsubdivisions * lb_nsubdivisions * 12;    //Numpts with lower bound nSubdivs
    int t2_npts = ub_nsubdivisions * ub_nsubdivisions * 12;    //Numpts with upper bound nSubdivs

    //Choose nSubdivision factor that gives Numpts close to N
    if(fabs(N-t1_npts) < fabs(N-t2_npts))
        nsubdivisions = lb_nsubdivisions;
    else
        nsubdivisions = ub_nsubdivisions;
    return nsubdivisions;
}

void xy2phitheta(double &phi, double &theta, const double &x, const double& y){
    double sigma = 10e-10, xc=10e-10, eps=10e-10;
    if(fabs(y) > M_PI_4){ //Polar
        sigma = 2 - fabs(4*y)/M_PI;
        xc = -M_PI + (2*floor((2*(x+M_PI))/M_PI) + 1) * M_PI_4;
        if(fabs(y) > (M_PI_2 - eps)){
            phi = xc + (x-xc);
            theta = sgn(y) * asin(1-((sigma*sigma)/3.0));
        }
        else{
            phi = xc + (x-xc)/sigma;
            theta = sgn(y) * asin(1-((sigma*sigma)/3.0));
        }
    }
    else{ //Equatorial
        phi = x;
        theta = asin((8*y)/(3*M_PI));
    }
}

void phitheta2xyz(double &x, double &y, double &z, const double &phi, const double &theta){
    //here theta is latitude, that is, theta = 90 - colatitude ,
    //where colatitude [0,PI];
    x = cos(theta) * cos(phi);
    y = cos(theta) * sin(phi);
    z = sin(theta);
}

void uv2xy(double &x, double &y, const double &u, const double &v, int face, int nstep){
                double nsteps = nstep;
    if(face == 1){
        x = 1*M_PI_4 + ((u-v)*M_PI_4)/nsteps;
        y = 0*M_PI_4 + ((u+v)*M_PI_4)/nsteps;
        //uv2xy1[{u_,v_}] := {Pi 1/4+(u-v)Pi/4/nsteps,Pi 0/4+(u+v)Pi/4/nsteps};
    }
    else if (face == 2){
        x = 2*M_PI_4 + ((u+v)*M_PI_4)/nsteps;
        y = 1*M_PI_4 + ((-u+v)*M_PI_4)/nsteps;
        //uv2xy2[{u_,v_}] := {Pi 2/4+(u+v)Pi/4/nsteps,Pi 1/4+(-u+v)Pi/4/nsteps};
    }
    else if (face == 3){
        x = 5*M_PI_4 + ((-u+v)*M_PI_4)/nsteps;
        y = 2*M_PI_4 + ((-u-v)*M_PI_4)/nsteps;
        //uv2xy3[{u_,v_}] := {Pi 5/4+(-u+v)Pi/4/nsteps,Pi 2/4+(-u-v)Pi/4/nsteps};
    }
    else if (face == 4){
        x = 8*M_PI_4 + ((-u-v)*M_PI_4)/nsteps;
        y = M_PI_4 + ((u-v)*M_PI_4)/nsteps;
        //uv2xy4[{u_,v_}] := {Pi 8/4+(-u-v)Pi/4/nsteps,Pi 1/4+(u-v)Pi/4/nsteps};
    }
    else if (face == 5){
        x = 5*M_PI_4 + ((-u+v)*M_PI_4)/nsteps;
        y = 0        + ((-u-v)*M_PI_4)/nsteps;
        //uv2xy5[{u_,v_}] := {Pi 5/4+(-u+v)Pi/4/nsteps,-Pi 0/4+(-u-v)Pi/4/nsteps};
    }
    else if (face == 6){
        x = 4*M_PI_4 + ((-u-v)*M_PI_4)/nsteps;
        y = -1*M_PI_4 +((u-v)*M_PI_4)/nsteps;
        //uv2xy6[{u_,v_}] := {Pi 4/4+(-u-v)Pi/4/nsteps,-Pi 1/4+(u-v)Pi/4/nsteps};
    }
    else if(face == 7){
        x =  1*M_PI_4 + ((u-v)*M_PI_4)/nsteps;
        y = -2*M_PI_4 + ((u+v)*M_PI_4)/nsteps;
        //uv2xy7[{u_,v_}] := {Pi 1/4+(u-v)Pi/4/nsteps,-Pi 2/4+(u+v)Pi/4/nsteps};
    }
    else if (face == 8){
        x =  (M_PI *(6.0/4.0)) + ((u+v)*M_PI_4)/nsteps;
        y = -1*M_PI_4 + ((-u+v)*M_PI_4)/nsteps;
        //uv2xy8[{u_,v_}] := {Pi 6/4+(u+v)Pi/4/nsteps,-Pi 1/4+(-u+v)Pi/4/nsteps};
    }
    else if (face == 9){
        x =  2*M_PI_4 + ((u-v)*M_PI_4)/nsteps;
        y = -1*M_PI_4 + ((u+v)*M_PI_4)/nsteps;
        //uv2xy9[{u_,v_}] := {Pi 2/4+(u-v)Pi/4/nsteps,-Pi 1/4+(u+v)Pi/4/nsteps};
    }
    else if (face == 10){
        x =  8*M_PI_4 + ((u-v)*M_PI_4)/nsteps;
        y = -1*M_PI_4 + ((u+v)*M_PI_4)/nsteps;
        //uv2xy10[{u_,v_}] := {Pi 8/4+(u-v)Pi/4/nsteps,-Pi 1/4+(u+v)Pi/4/nsteps};
    }
    else if (face == 11){
        x = 4*M_PI_4 + ((-u+v)*M_PI_4)/nsteps;
        y = 1*M_PI_4 + ((-u-v)*M_PI_4)/nsteps;
        //uv2xy11[{u_,v_}] := {Pi 4/4+(-u+v)Pi/4/nsteps,Pi 1/4+(-u-v)Pi/4/nsteps};
    }
    else if (face == 12){
        x = 6*M_PI_4 + ((-u+v)*M_PI_4)/nsteps;
        y = 1*M_PI_4 + ((-u-v)*M_PI_4)/nsteps;
        //uv2xy12[{u_,v_}] := {Pi 6/4+(-u+v)Pi/4/nsteps,Pi 1/4+(-u-v)Pi/4/nsteps};
    }
}

