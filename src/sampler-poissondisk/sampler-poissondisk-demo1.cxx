// $Id: poisson.cpp 790 2011-01-06 19:10:37Z mag $
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include <iomanip>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/timer/timer.hpp>

#include <stk/stk.hpp>

#include "utils.h"
#include <sys/stat.h>

#include "generator.h"

#define NDIM 2

namespace boostPO = boost::program_options;

int main(int argc, char** argv)
{
    /* ARG PARSER *****************************************************/
    std::string fn_output;
    std::string method;
    int nPts;
    int nPatches;

    boostPO::variables_map vm;
    boostPO::options_description desc("Allowed options");
    desc.add_options()
        ("help,h",
            "produce help message")
        ("output,o",
            boostPO::value<std::string>(&fn_output),
            "output filename")
        ("nPts,n",
            boostPO::value<int>(&nPts)->default_value(4096),
            "number of point per point set")
        ("nPatches,m",
            boostPO::value<int>(&nPatches)->default_value(1),
            "number of point sets")
        ("binary,b",
            "write in binary mode")
        ;

    boostPO::positional_options_description p;
    p.add("input", -1);

    try
    {
        boostPO::store(
            boostPO::command_line_parser(argc, argv).
                options(desc).positional(p).run(), vm);
        boostPO::notify(vm);
    }
    catch(boost::program_options::error& e)
    {
        std::cerr << e.what() << std::endl;
        std::cout << desc << std::endl;
        exit(EXIT_FAILURE);
    }

    if(vm.count("help"))
    {
        std::cout << desc << std::endl;
        exit(EXIT_SUCCESS);
    }

    if(vm.count("output") == 0)
    {
        std::cerr << "the option '--output' is required but missing" << std::endl;
        std::cout << desc << std::endl;
        exit(EXIT_FAILURE);
    }

    /* MAIN ***********************************************************/

    try
    {
        srand(time(NULL));

        //Open Output
        stk::io::PointSetOutputStream<2, double, double> stream;
        stream.setValueType(stk::io::PointSetStream::VAL_NONE);
        stream.setPositionType(stk::io::PointSetStream::POS_DOUBLE);
        if(vm.count("binary"))
        {
            stream.setBinary(true);
        }
        stream.open(fn_output);

        double meanTime = 0;
        double meanPts = 0;

        std::cout << "Point set " << 1 << "/" << nPatches << "       \r";
        std::cout.flush();

        float radiusFactor;
        if(nPts <= 1024) radiusFactor = 0.5215;
        else if(nPts >= 4096) radiusFactor = 0.54701;
        else radiusFactor = 0.5215 + (0.54701-0.5215)*((float)nPts-1024.0f)/(4096.0f-1024.0f);

        uint32_t num = nPts*1.05;

        float rad = sqrtf(radiusFactor/static_cast<float>(M_PI*num));

        Sample<NDIM>::setRadius(2.0f*rad);

        //std::cout << "Radius: " << 2.0f*rad << std::endl;

        Counter<NDIM>::init();
        Interval::seed(static_cast<uint32_t>(time(NULL)));
        BigNum<NDIM>::seed(static_cast<uint32_t>(time(NULL)));

        //###############Creating Folders###################################
        boost::filesystem::path source_dir_path( boost::filesystem::current_path() );

        std::cerr << source_dir_path.string() << std::endl;

        std::stringstream ss;
        std::string datafiles, images, graphs;
        ss.str(std::string());
        ss << source_dir_path.string() << "/results/";
        std::string resultFolder = ss.str();
        mkdir(resultFolder.c_str() ,0755);
        ss.str(std::string());
        ss << resultFolder << "pointset-poissodiskT2-n" << nPts << "/";

        mkdir(ss.str().c_str(),0755);

        create_folders(ss.str(), datafiles, images, graphs);
        //##########################################################

        for(int psCur=1; psCur<=nPatches; psCur++)
        {
            stk::PointSet2dd pts;

            boost::timer::cpu_timer timer;

            PDGenerator<NDIM>* gen = 0;

            gen = new PDGenerator<NDIM>(false);

            while (gen->iterate());

            gen->output(pts);

            delete gen;

            Sample<NDIM>::numsamples = 0;
            Sample<NDIM>::sample_pool.clear();

            boost::timer::nanosecond_type timeSystem = timer.elapsed().system;
            boost::timer::nanosecond_type timeUser = timer.elapsed().user;
            boost::timer::nanosecond_type timeWall = timer.elapsed().wall;
            meanTime += (static_cast<double>(timeWall)/1000000000.0 - meanTime)/(psCur+1);

            meanPts += (static_cast<double>(pts.size()) - meanPts)/(psCur+1);

            fprintf(stderr,"\r meanTime %f meanPts %f patch %d",meanTime, meanPts, psCur);
            //std::cout << meanTime << " sec, " << meanPts << " pts" << std::endl;

            ss.str(std::string());
            ss << psCur;
            std::string s1 = ss.str();
            paddedzerosN(s1, nPatches);

            ss.str(std::string());
            ss << datafiles << "pointset-poissodiskT2-n" << nPts << "-" << s1 << ".txt";

            std::ofstream file;
            file.open(ss.str().c_str());

            for(int i = 0; i < pts.size(); i++){
                file << pts.at(i).pos()[0] << " " << pts.at(i).pos()[1] << std::endl;
            }
            file.close();
            //stream.write(pts);
        }

        std::cerr << std::endl;

        std::cout << nPatches << (nPatches > 1 ? " point sets" : " point set") << " generated                             " << std::endl;

        stream.close();
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    exit(EXIT_SUCCESS);
}

