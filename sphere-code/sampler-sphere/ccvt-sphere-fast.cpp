/* CCVT Fast Implementation
 * */
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <queue>
#include <exception>
#include <boost/program_options.hpp>
#include <boost/random.hpp>

#include <iomanip>
#include <sys/stat.h>

#include "stk/pointset.hpp"

#include "ccvt-sphere.h"
#include "sampler-sphere.h"
#include "./../domaintransform/domaintransform.h"
#include "./../core/utils-sphere.h"
#include "./../core/constants.h"

/*

//Geodesic distance
template<int DIM, typename POS>
double distance(const stk::Vector<DIM, POS>& a, const stk::Vector<DIM, POS>& b){
    double ax[3],bx[3];

    thetaphi2xyz(ax, a[0], a[1]);
    thetaphi2xyz(bx, b[0], b[1]);
    //normalize_vector(ax);
    //normalize_vector(bx);

    double adotb = 0;
    for(int k=0; k<3;k++)
        adotb += ax[k]*bx[k];

    return acos(adotb);
}

int _ccvtRndGen(int m){
    static boost::mt19937 rng(0);
    static boost::uniform_01<boost::mt19937&> zeroone(rng);
    return (int)(zeroone()*(double)m)%m;
}

template<int DIM, typename POS>
void costList(
        std::map<int, std::map<double, int> >& m,
        double& localDist,
        unsigned int s,
        const stk::Vector<DIM, POS>& sp0,
        const stk::Vector<DIM, POS>& sp1,
        const double& sw0,
        const double& sw1,
        const stk::PointSet<DIM, POS, int>& points,
        std::vector< std::vector<int> >& ptsPerSites)
{
    for(int k=0; k < ptsPerSites[s].size(); k++)
    {
        int pid = ptsPerSites[s][k];
        const stk::Vector<DIM, POS>& p = points[pid].pos();
        double dist0 = distance(p, sp0)/sw0;
        double dist1 = distance(p, sp1)/sw1;
        double energy = -pow(dist1, 2) + pow(dist0, 2);
        m[points[pid].val()][energy] = pid;

        if(localDist < dist0) localDist = dist0;
    }
}


template<int DIM, typename POS>
void ccvt_sphere_fast(stk::PointSet<DIM, POS, double>& sites, //Final point position
                 const stk::PointSet<DIM, POS, int>& points, //Discretization
                 std::vector< std::vector<int> >& ptsPerSites){
    //Init
    std::vector<int> siteRand;
    siteRand.resize(sites.size());
    for(int i=0; i<sites.size(); i++) siteRand[i] = i;

    //Site bound
    std::vector<double> siteBound;
    siteBound.resize(sites.size());
    for(int i=0; i<sites.size(); i++) siteBound[i] = 1.0;

    //std::cout << "init OK" << std::endl;
    //Loop
    int changes = -1;
    int iter = 0;
    while(changes != 0)
    {
        changes = 0;

        std::random_shuffle(siteRand.begin(), siteRand.end(), _ccvtRndGen);

        for(int i=1; i<sites.size(); i++)
        {
            for(int j=0; j<i; j++)
            {
                int s0 = siteRand[i];
                int s1 = siteRand[j];

                const stk::Vector<DIM, POS>& sp0 = sites[s0].pos();
                const stk::Vector<DIM, POS>& sp1 = sites[s1].pos();

                const double& sw0 = sites[s0].val();
                const double& sw1 = sites[s1].val();

                if(siteBound[s0] + siteBound[s1] < distance(sp0, sp1))
                {
                    continue;
                }

                double localDist0 = 0.0;
                double localDist1 = 0.0;

                siteBound[s0] = 0.0;
                std::map<int, std::map<double, int> > pm0;
                costList(pm0, localDist0, s0, sp0, sp1, sw0, sw1, points, ptsPerSites);

                siteBound[s1] = 0.0;
                std::map<int, std::map<double, int> > pm1;
                costList(pm1, localDist1, s1, sp1, sp0, sw1, sw0, points, ptsPerSites);

                int changesInCouple = 0;

                std::map<int, std::map<double, int> >::iterator classIter = pm0.begin();
                while(classIter != pm0.end())
                {
                    int key = classIter->first;

                    std::map<int, std::map<double, int> >::iterator classIter2;
                    classIter2 = pm1.find(key);
                    if(classIter2 != pm1.end())
                    {
                        std::map<double, int>& ppm0 = classIter->second;
                        std::map<double, int>& ppm1 = classIter2->second;

                        std::map<double, int>::reverse_iterator iter0 = ppm0.rbegin();
                        std::map<double, int>::reverse_iterator iter1 = ppm1.rbegin();
                        while(iter0 != ppm0.rend() && iter1 != ppm1.rend() && iter0->first  + iter1->first > 0.0)
                        {
                            std::vector<int>::iterator swap0;
                            swap0 = std::find(
                                        ptsPerSites[s0].begin(),
                                        ptsPerSites[s0].end(),
                                        iter0->second);
                            *swap0 = iter1->second;

                            std::vector<int>::iterator swap1;
                            swap1 = std::find(
                                        ptsPerSites[s1].begin(),
                                        ptsPerSites[s1].end(),
                                        iter1->second);
                            *swap1 = iter0->second;

                            changesInCouple++;

                            iter0++;
                            iter1++;
                        }
                    }
                    classIter++;
                }

                changes += changesInCouple;

                if(changesInCouple)
                {
                    siteBound[s0] = 1.0;
                    siteBound[s1] = 1.0;
                }
                else
                {
                    siteBound[s0] = localDist0*2.0;
                    siteBound[s1] = localDist1*2.0;
                }
            }
        }

        //std::cout << "change : " << changes << std::endl;
        //fprintf(stderr,"\r iter: %d change: %d", iter, changes);
        //iter++;

        //##########################################################################
        //On Sphere: Compute Barycenter of points and update sites to the barycenters
        for(int i=0; i<sites.size(); i++){
            double svec[] = {0,0,0};
            thetaphi2xyz(svec, sites[i].pos()[0], sites[i].pos()[1]);
            //normalize_vector(svec);

            double mean[] = {0,0,0};
            for(int j=0; j < ptsPerSites[i].size(); j++){
                int pid = ptsPerSites[i][j];
                const stk::Vector<DIM, POS>& p = points[pid].pos();

                double pvec[] = {0,0,0};
                thetaphi2xyz(pvec, p[0], p[1]);
                //normalize_vector(pvec);
                for(int k=0; k<3; k++)
                    mean[k] += pvec[k];
            }
            for(int k=0; k<3; k++)
                mean[k] /= (double)ptsPerSites[i].size();

            normalize_vector(mean);

            //Compute theta_phi for the mean position
            double ntheta=0, nphi=0;
            xyz2thetaphi(ntheta, nphi, mean[0], mean[1], mean[2]);
            normalize_theta_phi(ntheta, nphi);
            sites[i].pos()[0] = ntheta;
            sites[i].pos()[1] = nphi;
            //sites[i].pos() += m/(double)ptsPerSites[i].size();
        }
    }
    std::cerr << std::endl;
}

template void ccvt_sphere_fast(stk::PointSet<2, double, double>& sites, //Final point position
                 const stk::PointSet<2, double, int>& points, //Discretization
                 std::vector< std::vector<int> >& ptsPerSites);

*/
void coherent_initialization(stk::PointSet2dd& sites, stk::PointSet2di& points,
                                 std::vector< std::vector<int> >& ptsPerSites, int nsites,
                                 int npts, std::string samplingpattern){
    int sId = 0;
    int pId = 0;
    int totalPts = nsites * npts;
    std::vector<double> P;  //set of allPoints

    //#####################################################################

    if(samplingpattern == "whitenoise"){
        //Complete set of points
        for(int i=0; i < totalPts; i++){
            double p1 = drand48();
            double p2 = drand48();
            double px=0, py=0,pz=0;

            uniformSampleSphere(p1, p2, &px,&py,&pz);

            double ptheta=0.0,pphi=0.0;
            xyz2thetaphi(ptheta, pphi, px,py,pz);

            P.push_back(ptheta);
            P.push_back(pphi);
        }

        //#####################################################################

        //Complete set of sites
        for(int i = 0; i< nsites; i++){
            //sId = sites.size();
            double s1 = drand48();
            double s2 = drand48();
            double sx=0, sy=0,sz=0;

            uniformSampleSphere(s1, s2, &sx, &sy,&sz);

            double stheta=0.0,sphi=0.0;

            xyz2thetaphi(stheta, sphi, sx,sy,sz);
            //std::cout << sx <<" " << sy << "  " << sz << std::endl;
            sites.push_back(stk::Point2dd(stk::Vector2d(stheta,sphi), 1.0));
            //sites.push_back(stk::Point2dd(stk::Vector2d(s1,s2), 1.0));
            ptsPerSites.push_back(std::vector<int>());
        }

        //#####################################################################

        //Assign points to each site;

        for(int i=0; i < P.size();){
            double minDist = 1e10;

            int sId = 0;
            for(int j=0; j<sites.size(); j++){

                stk::Vector<2, double> a(P[i+0],P[i+1]);

                double dist = distance(sites[j].pos(), a);

                if(ptsPerSites[j].size() == npts)
                    continue;
                else{
                    if(minDist > dist){
                        minDist = dist;
                        sId = j;
                    }
                }
            }

            if(ptsPerSites[sId].size() == npts){
                continue;
            }
            else{
                pId = points.size();
                points.push_back(stk::Point2di(stk::Vector2d(P[i+0],P[i+1]), 0));
                ptsPerSites[sId].push_back(pId);
                i+=2;
            }
        }
    }

    //#####################################################################
    /*else if(samplingpattern == "stratified"){

        //Compute all sites:
        std::list<double> S;
        std::vector<double> stratified_sites;
        stratified_sites = spherical_stratified_samples(nsites);
        int temp_nsites = stratified_sites.size() / 3.0;
        for(int i = 0; i < temp_nsites ; i++){
            sId = i;
            double sx = stratified_sites[3*i+0];
            double sy = stratified_sites[3*i+1];
            double sz = stratified_sites[3*i+2];

            double stheta=0.0,sphi=0.0;

            xyz2thetaphi(stheta, sphi, sx,sy,sz);
            //std::cout << sx <<" " << sy << "  " << sz << std::endl;
            sites.push_back(stk::Point2dd(stk::Vector2d(stheta,sphi), 1.0));
            //sites.push_back(stk::Point2dd(stk::Vector2d(s1,s2), 1.0));
            S.push_back(sId);
            S.push_back(stheta);
            S.push_back(sphi);


            ptsPerSites.push_back(std::vector<int>());
        }
        //#####################################################################
        totalPts = sites.size() * npts;

        //Complete set of all points
        for(int i=0; i < totalPts; i++){
            double p1 = drand48();
            double p2 = drand48();
            double px=0, py=0,pz=0;

            uniformSampleSphere(p1, p2, &px,&py,&pz);

            double ptheta=0.0,pphi=0.0;
            xyz2thetaphi(ptheta, pphi, px,py,pz);

            P.push_back(ptheta);
            P.push_back(pphi);
        }

        //#####################################################################
        //Assign points to each site;
        int count=0;
        for(int i=0; i < P.size();i++){
            double minDist = 1e10;
            count++;

            int sId = 0;
            std::list<double>::iterator it1, it2;
            for(auto v = S.begin(); v != S.end();){
                it1 = v;
                int id = *v;
                v++;
                double stheta = *v;
                v++;
                double sphi = *v;
                it2 = ++v;
                stk::Vector<2, double> a(P[i+0],P[i+1]);
                stk::Vector<2,double> s(stheta, sphi);
                double dist = distance(s, a);
                if(minDist > dist){
                    minDist = dist;
                    sId = id;
                }
            }
            pId = points.size();
            points.push_back(stk::Point2di(stk::Vector2d(P[i+0],P[i+1]), 0));
            ptsPerSites[sId].push_back(pId);
            if(ptsPerSites[sId].size() == npts){
                S.erase(it1,it2);
            }
        }
        //#####################################################################

            }
        /*
        //Assign points to each site;
        int count=0;
        for(int i=0; i < P.size();i++){
            double minDist = 1e10;
            count++;

            int sId = 0;
            for(int j=0; j<sites.size(); j++){

                stk::Vector<2, double> a(P[i+0],P[i+1]);

                double dist = distance(sites[j].pos(), a);

                if(ptsPerSites[j].size() == npts)
                    continue;
                else{
                    if(minDist > dist){
                        minDist = dist;
                        sId = j;
                    }
                }
            }

            if(ptsPerSites[sId].size() == npts){
                continue;
            }
            else{
                pId = points.size();
                points.push_back(stk::Point2di(stk::Vector2d(P[i+0],P[i+1]), 0));
                ptsPerSites[sId].push_back(pId);
                i+=2;
            }
        } *
    }
    //#####################################################################
    /*else  if(samplingpattern == "stratified-notused"){
        std::vector<double> stratified_sites;
        stratified_sites = spherical_stratified_samples(nsites);
        int temp_nsites = stratified_sites.size() / 3.0;
        for(int i = 0; i < temp_nsites ; i++){
            sId = i;
            double sx = stratified_sites[3*i+0];
            double sy = stratified_sites[3*i+1];
            double sz = stratified_sites[3*i+2];

            double stheta=0.0,sphi=0.0;

            xyz2thetaphi(stheta, sphi, sx,sy,sz);
            //std::cout << sx <<" " << sy << "  " << sz << std::endl;
            sites.push_back(stk::Point2dd(stk::Vector2d(stheta,sphi), 1.0));
            //sites.push_back(stk::Point2dd(stk::Vector2d(s1,s2), 1.0));

            ptsPerSites.push_back(std::vector<int>());

            for(int j=0; j < npts; j++){
                pId = points.size();

                double p1 = drand48();
                double p2 = drand48();
                double px=0, py=0,pz=0;

                uniformSampleSphere(p1, p2, &px,&py,&pz);

                double ptheta=0.0,pphi=0.0;
                xyz2thetaphi(ptheta, pphi, px,py,pz);

                points.push_back(stk::Point2di(stk::Vector2d(ptheta,pphi), 0));
                //points.push_back(stk::Point2di(stk::Vector2d(p1,p2), 0));

                ptsPerSites[sId].push_back(pId);
            }
        }
    } */
}

/* */

