#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/date_time/gregorian/gregorian.hpp"
#include <boost/filesystem.hpp>

#include <sstream>
#include <sys/stat.h>
#include "utils.h"

using namespace boost::gregorian;

void create_folders(std::string homedir, std::string &data, std::string &images, std::string &graphs){

    date d(day_clock::local_day());
    date_facet* facet(new date_facet("%Y-%m-%d"));

//#####################Setting up folders to store data#########################

    std::stringstream ss, ssdate, ssdir, ssdata, ssimg, ssgraph;
    ssdate.imbue(std::locale(ssdate.getloc(), facet));
    ssdate << d;

    //const char* homeDir = getenv ("Home");
    //char final [256];
    //sprintf (final, "%s/Desktop/%s",homeDir, game_name);

    ssdir << homedir << ssdate.str();
    if(!mkdir(ssdir.str().c_str(),0775))
        std::cerr <<"Directory with date!"<< std::endl;

    ssdata << ssdir.str() << "/datafiles/";
    if(!mkdir(ssdata.str().c_str(),0775))
        std::cerr <<"Directory with datafiles!"<< std::endl;

    ssimg << ssdir.str() << "/images/";
    if(!mkdir(ssimg.str().c_str(),0775))
        std::cerr <<"Directory with images!"<< std::endl;

    ssgraph << ssdir.str() << "/graphs/";
    if(!mkdir(ssgraph.str().c_str(),0775))
        std::cerr <<"Directory with graphs!"<< std::endl;

/*
    boost::filesystem::path p1(ssdir.str());

    ssdata << ssdir.str() << "/datafiles/";
    boost::filesystem::path p2(ssdata.str());

    ssimg << ssdir.str() << "/images/";
    boost::filesystem::path p3(ssimg.str());

    ssgraph << ssdir.str() << "/graphs/";
    boost::filesystem::path p4(ssgraph.str());

    try{
        if (boost::filesystem::create_directory(p1)){
            std::cerr <<"Directory with date!"<< std::endl;
        }
/*        if (boost::filesystem::create_directory(p2)){
            std::cerr <<"Directory with datafiles! "<< std::endl;
        }
        if (boost::filesystem::create_directory(p3)){
            std::cerr <<"Directory with images! "<< std::endl;
        }
        if (boost::filesystem::create_directory(p4)){
            std::cerr <<"Directory with graphs! "<< std::endl;
        }

    }
    catch (boost::filesystem::filesystem_error &e){
        std::cerr << e.what() << std::endl;
    }
*/
    data = ssdata.str();
    images = ssimg.str();
    graphs = ssgraph.str();
//############################Folders have been created!###########################

}
