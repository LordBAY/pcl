/**
 * @file demo.cpp
 * @brief Demo of super quadric fitting code
 * @author A. Huaman Quispe <ahuaman3@gatech.edu>
 * @date 2014/08/12
 */
#include <iostream>
#include <pcl/io/openni_grabber.h>
#include <pcl/io/pcd_io.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/console/parse.h>

#include <tabletop_segmentation/tabletop_segmentation.h>
#include <tabletop_symmetry/mindGapper.h>
#include <SQ_fitter.h>

#define PointT pcl::PointXYZRGBA

//-- Global variables
boost::shared_ptr<pcl::visualization::CloudViewer> gViewer;

pcl::Grabber* gKinectGrabber;
unsigned int gFilesSaved = 0;
bool gProcessCloud = false;
bool gReadFromFile = false;
std::string gInputFile;
int gType = 0;

TabletopSegmentor<PointT> gTs;

//-- Visualization
std::vector<pcl::PointCloud<PointT>::Ptr> gFittedClouds;


//-- Functions declaration
void printUsage( char* argv0 );
void grabberCallback( const pcl::PointCloud<PointT>::ConstPtr &_cloud );
void keyboardEventOccurred( const pcl::visualization::KeyboardEvent &_event,
			    void *_nothing );
boost::shared_ptr<pcl::visualization::CloudViewer> createViewer( std::string _name );
void processCloud( const pcl::PointCloud<PointT>::Ptr &_cloud );


/**
 * @function main
 * @brief Main function what else?
 */
int main( int argc, char* argv[] ) {

  int c;
  while( (c=getopt(argc, argv,"hr:t:")) != -1 ) {
    switch(c) {
    case 'h': { printUsage(argv[0]); return 1;} break;	
    case 'r': {
      gReadFromFile = true;
      gInputFile = optarg;
    } break;
    case 't': {
      gType = atoi(optarg) == 0? 0 : 1;
    } break;
    case '?': {
      if( optopt == 'r' ) {
	std::cout << "-r Required the input pointcloud file"<< std::endl;
	return 1;
      }
    } break;
    }
  } // end while
  
  /** If reading pointcloud from file */
  if( gReadFromFile ) {
    pcl::PointCloud<PointT>::Ptr cloud( new pcl::PointCloud<PointT>() );
    if( pcl::io::loadPCDFile<PointT>(gInputFile, *cloud ) == -1 ) {
      std::cout << "\t * [ERROR] Reading pointcloud file "<< gInputFile << std::endl;
      return -1;
    } else {
      processCloud( cloud );
    }
  } 
  /** If obtaining pointcloud from OpenNI grabber */
  else {
  
    // Create OpenNI grabber, register key events
    gKinectGrabber = new pcl::OpenNIGrabber();
    if( gKinectGrabber == 0 ) { return 1; }
    boost::function< void (const pcl::PointCloud<PointT>::ConstPtr&) > f = boost::bind( &grabberCallback, _1 );
    gKinectGrabber->registerCallback(f);
    
    // Create viewer and initialize
    gViewer = createViewer("viewOpenNI");
    gKinectGrabber->start();
    
    // Loop expecting capture signal (press ESC)
    while( !gViewer->wasStopped() ) {
      boost::this_thread::sleep( boost::posix_time::seconds(1) );
    }
    
    // If capturing Kinect data, stop the stream before exiting
    gKinectGrabber->stop();  
    
  } // end else
  
  
} // end main



/**
 * @function printUsage
 */
void printUsage( char* argv0 ) {
  std::cout<< "**  Usage:" << std::endl;
  std::cout<< argv0 <<" [ARGS] "<< std::endl;
  std::cout<<" <none>: : Start capturing from a Kinect device, save by pressing SPACE"<<std::endl;
  std::cout<<" -h : Show this help"<<std::endl;
  std::cout<<" -r : Input a pointcloud file"<<std::endl;
  std::cout<<" -t: Type of optimizer (0: CERES, 1: LEVMAR)"<< std::endl;
}

/**
 * @function processCloud
 */
void processCloud( const pcl::PointCloud<PointT>::Ptr &_cloud ) {

  
  /******** STEP 0: INIT PROCESSING FOR CLOUD ********/
    if( !gReadFromFile ) {

      std::stringstream stream;
      stream << "inputCloud" << gFilesSaved << ".pcd";
      std::string filename = stream.str();
      
      if( pcl::io::savePCDFile( filename, *_cloud, true ) == 0 ) {
	gFilesSaved++;
	std::cout << "\t * [INFO] Saved " << filename << "." << std::endl;
      } else {
	std::cout<<"\t * [ERROR] Problem saving "<< filename.c_str() << std::endl;
	return;
      }
    } // end if gReadFromFile

    /************ STEP 1: SEGMENT CLOUD *******************/
    gTs.processCloud( _cloud );


    /*********** STEP 2: MIRROR + FIT ********************/
    // Get camera parameters
    double f = 525;
    int width = 640;
    int height = 480;
    double cx = (double) width / 2;
    double cy = (double) height / 2;

    /** If reading from OpenNIgrabber, get the camera information from the device */
    if( !gReadFromFile ) {
      XnMapOutputMode m = ((pcl::OpenNIGrabber*)gKinectGrabber)->getDevice()->getDepthOutputMode();
      width = (int) m.nXRes;
      height = (int) m.nYRes;
      
      f = (double)((pcl::OpenNIGrabber*)gKinectGrabber)->getDevice()->getDepthFocalLength(0);
      cx = width >> 1;
      cy = height >> 1;
    }
    
    /** Set mindGapper (to generate mirror pointclouds) */
    mindGapper<PointT> mG;
    mG.setTablePlane( gTs.getTableCoeffs() );
    mG.setFittingParams( 6, 5, 0.01,M_PI/18.0 );
    mG.setDeviceParams( width, height, f, (double)cx, (double)cy );
    
    /** Set fitter */
    SQ_fitter<PointT> fitter;
    pcl::PointCloud<PointT>::Ptr fitCloud( new pcl::PointCloud<PointT>() );
    clock_t ts, tf; double dt;
    bool r;
    double smax, smin; int N; double thresh;
    char name[50];

    thresh = 0.1;
    N = 5;
    smax = 0.05; smin = 0.01;
    

    for( int i = 0; i < gTs.getNumClusters(); ++i ) {
      pcl::PointCloud<PointT>::Ptr cluster( new pcl::PointCloud<PointT>() );
      *cluster = gTs.getCluster(i);

      /** Generate mirror pointclouds */ 
      mG.complete( cluster );

      sprintf( name, "mirror_%d.pcd", i );
      pcl::io::savePCDFile( name, *cluster, true );

      /** Fitting SQ to mirrored pointcloud */
      fitter.setInputCloud( cluster );
      ts = clock();
      r = fitter.fit( gType, 
		      smax, smin,
		      N, thresh );
      tf = clock();
      dt = (double)(tf-ts)/CLOCKS_PER_SEC;
      std::cout << "[Cluster "<<i<<"] Fitting time: "<< dt<< std::endl;
      if(r) {
	fitter.printResults();
	fitCloud = fitter.getSampledOutput();
	gFittedClouds.push_back( fitCloud );
       
	sprintf( name, "fit_%d.pcd", i );
	pcl::io::savePCDFile( name, *fitCloud, true );
      } // end if fit was good

    } // end for num clusters


    // View
    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer( new pcl::visualization::PCLVisualizer("input pointcloud") );

    pcl::PointCloud<PointT>::Ptr input( new pcl::PointCloud<PointT>() );
    *input = gTs.getDownFilteredCloud();
    viewer->addPointCloud( _cloud, "inputCloud"  );
    int rc, gc, bc;
    for( int i = 0; i < gFittedClouds.size(); ++i ) {
      sprintf( name, "view_%d", i );
      rc = rand() % 255; gc = rand() % 255; bc = rand() % 255;
      pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZRGBA> col( gFittedClouds[i], rc, gc, bc);
      viewer->addPointCloud( gFittedClouds[i], col, name  );      
    }
    while( !viewer->wasStopped() ) {
	viewer->spinOnce(100);
	boost::this_thread::sleep( boost::posix_time::microseconds(100000) );
    }

    return;
}

/**
 * @function grabberCallback
 */
void grabberCallback( const pcl::PointCloud<PointT>::ConstPtr &_cloud ) {
  
  if( !gViewer->wasStopped() ) {
    gViewer->showCloud( _cloud );
  }
  
  if( gProcessCloud ) {
    pcl::PointCloud<PointT>::Ptr cloud( new pcl::PointCloud<PointT>() );
    *cloud = *_cloud;
    processCloud( cloud );
    gProcessCloud = false;
  }

}

/**
 * @function keyboardEventOccurred
 */
void keyboardEventOccurred( const pcl::visualization::KeyboardEvent &_event,
			   void *_nothing ) {
  
  // Save cloud when pressing space
  if( _event.getKeySym() == "space" && _event.keyDown() ) {
    gProcessCloud = true;
  }
}

/**
 * @function createViewer
 */
boost::shared_ptr<pcl::visualization::CloudViewer> createViewer( std::string _name ) {
  boost::shared_ptr<pcl::visualization::CloudViewer> v( new pcl::visualization::CloudViewer(_name) );
  v->registerKeyboardCallback( keyboardEventOccurred );

  return (v);
}
