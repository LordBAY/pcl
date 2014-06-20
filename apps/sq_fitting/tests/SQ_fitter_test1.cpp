/**
 * @file SQ_fitter_test1.cpp
 * @brief Test Fit function for Super quadrics
 */
#include <pcl/io/pcd_io.h>
#include <string>
#include <unistd.h>

#include <SQ_fitter.h>

/**
 * @function main
 */
int main( int argc, char* argv[] ) {

  // 0. Variables needed
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud( new pcl::PointCloud<pcl::PointXYZ>() );
  SQ_fitter<pcl::PointXYZ> fitter;
  
  
  // 1. Read pointcloud from input
  int v;
  std::string filename;
  while( (v=getopt( argc, argv, "p:h")) != -1 ) {
	
    switch(v) {
    case 'p': {
      filename = std::string( optarg );
      if( pcl::io::loadPCDFile<pcl::PointXYZ>( filename,
											   *cloud ) == -1 ) {
		std::cout << "\t [ERROR] Could not read pointcloud "<<
filename<< std::endl;
		return 1;
      } else {
		std::cout << "Read pointcloud with "<< cloud->points.size() << " points." << std::endl;
	  }

    } break;
	case 'h': {
	  std::cout <<"Syntax: "<<argv[0]<<" filename.pcd"<< std::endl;
	} break;
    } // end switch
	
  }
  
  if( filename.size() == 0 ) {
	std::cout << "Syntax: "<< argv[0]<< " -p filename.pcd"<< std::endl;
	return 1;
  }


  // 2. Load pointcloud in fitter object
  fitter.setInputCloud( cloud );

  
  double dim[3]; double trans[3]; double rot[3];
  double minDist;
  fitter.getBoundingBox( cloud, dim, trans, rot );
  if( dim[0] <= dim[1] && dim[0] <= dim[2] ) { minDist = dim[0]; }
  if( dim[1] <= dim[0] && dim[1] <= dim[2] ) { minDist = dim[1]; }
  if( dim[2] <= dim[0] && dim[2] <= dim[1] ) { minDist = dim[2]; }
  double smin = minDist*2.0 / 15.0;
  double smax = minDist*2.0 / 5.0;

  // 2. Fit. If successful, visualize and spit out summary
  if( fitter.fit( smax, smin,
				  5, 0.1) ) {

	std::cout << "\t [GOOD] Fit superquadric!"<< std::endl;

	// 3. Spit out result
	fitter.printResults();
	// 4. Visualize
	fitter.visualize();
  
  } else {
	std::cout << "\t [BAD] Did not fit properly" << std::endl;
  }
			
  return 0; 
}


// Local Variables:
// mode: c++
// tab-width: 4
// End:
