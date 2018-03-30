#include <ros/ros.h>
#include <iostream>

// #include <nurbsS.h>
#include <cmath>
// PCL specific includes
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

using namespace std;

void downsampleRow(Eigen::Array<bool,Eigen::Dynamic, 1>& rowFlags, int numRowsDesired){

  // Create linspaced array 
  Eigen::Array<float, Eigen::Dynamic, 1> selectArrayf;
  Eigen::Array<int, Eigen::Dynamic, 1> selectArray;
  selectArrayf.setLinSpaced(numRowsDesired, 0, rowFlags.count()-1);
  selectArray = selectArrayf.round().cast<int>();

  int j = 0;

  for (int i = 0; i < rowFlags.rows(); i++){
    if (rowFlags(i)){
      if (!(selectArray == j).any()){
        rowFlags(i) = false;
      }
      j++;
    }
  }

}

void downsampleCol(Eigen::Array<bool, 1, Eigen::Dynamic>& colFlags, int numcolsDesired){

  // Create linspaced array 
  Eigen::Array<float, Eigen::Dynamic, 1> selectArrayf;
  Eigen::Array<int, Eigen::Dynamic, 1> selectArray;
  selectArrayf.setLinSpaced(numcolsDesired, 0, colFlags.count()-1);
  selectArray = selectArrayf.round().cast<int>();
  
  int j = 0;

  for (int i = 0; i < colFlags.cols(); i++){
    if (colFlags(i)){
      if (!(selectArray == j).any()){
        colFlags(i) = false;
      }
      j++;
    }
  }

}


void regionAverage(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, int i, int j){

  pcl::PointXYZ average(0.0,0.0,0.0);

  int add_count = 0; // to track how many values are averaged. 

  // Back in i
  if (i > 0){
    if (pcl::isFinite(cloud->at(j,i-1))){
      average.getArray3fMap() += cloud->at(j,i-1).getArray3fMap();
      add_count ++;
    }    
  }

  // Forward in i
  if (i < cloud->height-1){
    if (pcl::isFinite(cloud->at(j,i+1))){
      average.getArray3fMap() += cloud->at(j,i+1).getArray3fMap();
      add_count ++;
    }    
  }

  // Back in j
  if (j > 0){
    if (pcl::isFinite(cloud->at(j-1,i))){
      average.getArray3fMap() += cloud->at(j-1,i).getArray3fMap();
      add_count ++;
    }    
  }

  // Forward in j
  if (j < cloud->width-1){
    if (pcl::isFinite(cloud->at(j+1,i))){
      average.getArray3fMap() += cloud->at(j+1,i).getArray3fMap();
      add_count ++;
    }    
  }

  // Divide by count to get the average
  average.getArray3fMap() /= (float)add_count;

  //update cloud
  cloud->at(i,j) = average;

}

void averageOutNans(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, Eigen::Array<int,2,Eigen::Dynamic> nanIndices){
  cout << "Number of cols: " << nanIndices.cols() << endl;
  for (int i = 0; i < nanIndices.cols(); i++){
    cout << "Before: " << cloud->at(nanIndices(0,i),nanIndices(1,i)) << endl;
    regionAverage(cloud,nanIndices(0,i),nanIndices(1,i));
    cout << "After:  " << cloud->at(nanIndices(0,i),nanIndices(1,i)) << endl;
  }

}

void removeRowsNan(Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>& nanArray,Eigen::Array<bool,Eigen::Dynamic, 1>& rowFlags, int maxNanAllowed){

  Eigen::Array<int, Eigen::Dynamic, 1> rowNan = nanArray.rowwise().count().cast<int>();

  // Maximum number of Nans in a row
  int maxNan = rowNan.maxCoeff();

  // If the maximum number Nans is above the set limit
  if (maxNan > maxNanAllowed){
    for (int i = 0; i < rowNan.rows(); i++){
      // If this rows has equal to the maximum number of nans
      if (rowNan(i) == maxNan){
        // Set the row values to zero in the array
        nanArray.row(i) = Eigen::Array<bool, 1, Eigen::Dynamic>::Zero(1,nanArray.cols());

        // Set the rowFlags value to zero
        rowFlags(i) = false;

      }
    }
  }

}

void removeColsNan(Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>& nanArray,Eigen::Array<bool, 1, Eigen::Dynamic>& colFlags, int maxNanAllowed){

  Eigen::Array<int, 1, Eigen::Dynamic> colNan = nanArray.colwise().count().cast<int>();

  // Maximum number of Nans in a col
  int maxNan = colNan.maxCoeff();

  // If the maximum number Nans is above the set limit
  if (maxNan > maxNanAllowed){
    for (int i = 0; i < colNan.cols(); i++){
      // If this cols has equal to the maximum number of nans
      if (colNan(i) == maxNan){
        // Set the col values to zero in the array
        nanArray.col(i) = Eigen::Array<bool, Eigen::Dynamic, 1>::Zero(nanArray.rows(), 1);

        // Set the colFlags value to zero
        colFlags(i) = false;

      }
    }
  }

}

void getNanMatrixFromPointCloud(Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>& nanArray, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud){

  // Expect nanArray to be all false to 

  // Loop through Point cloud 

  for (int i = 0; i < cloud->height; i++){
    for (int j = 0; j < cloud->width; j++){
      // If nan value
      if (!pcl::isFinite(cloud->at(j,i))){
        nanArray(i,j) = true;
      }else{
        nanArray(i,j) = false;
      }
    }
  }
  

}


void pclTest(){
  // pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
  pcl::PCLPointCloud2::Ptr cloud_blob (new pcl::PCLPointCloud2);
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);

  // Fill in the cloud data
  pcl::PCDReader reader;
//   reader.read ("/home/bjm/SpaceCRAFT/ros_ws/src/pcl_testing/data/table_scene_lms400.pcd", *cloud_blob);
  reader.read ("/home/bjm/Dropbox/PhD_Code/Data/3D_Scans/Blensor/Scan01/BlobScan_Data00010.pcd", *cloud_blob);
  // reader.read ("/home/bjm/Dropbox/PhD_Code/NURBS_2018_01_02/test_data.pcd", *cloud_blob);

  // Convert to the templated PointCloud (PointCoud<T> to PointCloud2)
  pcl::fromPCLPointCloud2 (*cloud_blob, *cloud);  

  std::cerr << "PointCloud before filtering: W: " << cloud->width << "H: " << cloud->height << " data points." << std::endl;

  cout << *cloud << endl;
  cout << cloud->points[0].x << endl;
  cout << cloud->points[1].x << endl;

  cout << "Accessing coordinates (0,0): " << cloud->at(0,0) << endl;
  cout << "Checking coordinates (0,0) is finite?: " << pcl::isFinite(cloud->at(0,0)) << endl;
  cout << "Accessing coordinates (0,0).x: " << cloud->at(0,0).x << endl;
  cout << "Accessing coordinates (5,4): " << cloud->at(5,4) << endl;
  cout << "Accessing coordinates (2,9): " << cloud->at(2,9) << endl;


  // cout << "X array" << cloud->points.x << endl;

  // Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic> arrayX = cloud->points.x;

  // Edit the cloud
  cloud->at(0,0).x = NAN;
  cout << "Accessing coordinates (0,0).x: " << cloud->at(0,0).x << endl;

  cloud->at(1,1) = pcl::PointXYZ(1.0,6.0,9.0);
  cout << "Accessing coordinates (1,1): " << cloud->at(1,1) << endl;

  regionAverage(cloud, 0, 0);

  cout << "\n\nNeighbour coords are: " << cloud->at(1,0) << " and " << cloud->at(0,1) << endl;
  cout << "Edited PC is coordinates (0,0): " << cloud->at(0,0) << endl;

  regionAverage(cloud, 5, 5);

  cout << "\n\nNeighbour coords are: " << cloud->at(5,4) << " and " << cloud->at(4,5) << " and " << cloud->at(5,6) << " and " << cloud->at(6,5) << endl;
  cout << "Edited PC is coordinates (0,0): " << cloud->at(5,5) << endl;


  Eigen::Array<int,2,5> nanIndices(2,5); nanIndices << 0, 5, 10, 15, 20,
                                                  10, 15, 20, 10, 15;

  cout << "nanIndices: \n" << nanIndices << endl;

  cout << "\ncloud indices before and after:" << endl;

  averageOutNans(cloud,nanIndices);


  // for (int i = 0; i < cloud->height; i++){
  //   for (int j = 0; j < cloud->width; j++){
  //     // if nan or infinite
  //     if (!pcl::isFinite(cloud->at(j,i))){        
  //       cloud->at(j,i).x = 0.0;// NOTE: at(row,col) so input is at(j, i)
  //       cloud->at(j,i).y = 0.0;
  //       cloud->at(j,i).z = 0.0;
  //       // cout << "replacing Nan at (i,j): (" << i << ", " << j << ")\n" ;
  //     }
  //   }
  // }

  // Eigen::Matrix<float, 45, 45> mask;
  Eigen::MatrixXf mask = Eigen::MatrixXf::Zero(5,5);
  // Eigen::MatrixXf mask = Eigen::MatrixXf::Identity(5,5);
  
  cout << "mask sum: " << mask.sum() << endl;
  cout << "mask row sum: " << mask.rowwise().sum() << endl;
  cout << "mask col sum: " << mask.colwise().sum() << endl;

  Eigen::ArrayXXf table(10, 4);
  table.col(0) = Eigen::ArrayXf::LinSpaced(10, 0, 90);
  table.col(1) = M_PI / 180 * table.col(0);
  table.col(2) = (0*table.col(1).sin())/0;
  table.col(3) = (0*table.col(1))/0;
  table(0,1) = (0*table(0,1))/0;
  table(4,0) = (0*table(4,0))/0;
  std::cout << "  Degrees   Radians      Sine    Cosine\n";
  std::cout << table << std::endl;


  cout << "Is table nan: " << table.isNaN() << endl;
  cout << "Sum of rows isNan:  " << table.isNaN().rowwise().sum() << endl;
  cout << "Sum of columns isNan:  " << table.isNaN().colwise().sum() << endl;
  cout << "Sum of rows  (0) isNan:  " << table.isNaN().row(0).sum() << endl;

  // Eigen::ArrayBase<Derived>::IsNaNReturnType nanMask;

  Eigen::Array<bool,Eigen::Dynamic,Eigen::Dynamic> nanMask(table.rows(),table.cols());
  
  nanMask = table.isNaN();
  

  cout << "IsNaN mask: " << nanMask << endl;

  cout << "IsNaN mask rows " << nanMask.rowwise().count() << endl;

  Eigen::Array<int, Eigen::Dynamic, 1> rowNan0(table.rows(),1);

  rowNan0 = nanMask.rowwise().count().cast<int>();

  cout << "rowNan is: " << rowNan0 << "\nRowNan term 2 is: " << rowNan0(2) << endl;


  nanMask.row(2) = Eigen::Array<bool, 1, Eigen::Dynamic>::Zero(1,nanMask.cols());

  cout << "\n\nNanMask is now: " << nanMask << endl;

  // Eigen::ArrayXf rowNaNSum = nanMask.rowwise().sum();

  // Eigen::isnan(table);
  Eigen::Array<bool,Eigen::Dynamic, 1> rowFlags(10,1); rowFlags << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;
  Eigen::Array<bool,1,Eigen::Dynamic> colFlags(1,4); colFlags << 1, 1, 1, 1;
  cout << "\nRowFlags are now: " << rowFlags << endl;

  removeRowsNan(nanMask, rowFlags, 1);

  cout << "\n\nAfter Remove Nans:\nNanMask is now: " << nanMask << endl;
  cout << "RowFlags are now: " << rowFlags << endl;

  // removeRowsNan(nanMask, rowFlags, 1);
  removeColsNan(nanMask, colFlags, 1);

  cout << "\n\nAfter Remove Nans Col:\nNanMask is now: " << nanMask << endl;
  cout << "ColFlags are now: " << colFlags << endl;  

  // Trial process - given NaNs
  // Initialise
  Eigen::Array<int, 1, Eigen::Dynamic> rowNan;
  
  // Get rows with Nans
  rowNan = nanMask.rowwise().count().cast<int>();

  // cout << rowNan << endl;

  cout << "Rowmax is: " << rowNan.maxCoeff() << endl;

  // Get index where it equals the maximum 



  // get the max of the sums;
  // int maxNan;
  // maxNan = rowNan.cwiseMax(0.0);

  // // Get when the array equals the max
  // Eigen::Array<bool, 1, Eigen::Dynamic> maxNanArr(1,rowNan.cols());
  // maxNanArr = rowNan == maxNan;

  // cout << "MaxNan is " << maxNan << ", and is max: " << maxNanArr << endl;


  Eigen::Array<int, Eigen::Dynamic, 1> selectArray;
  selectArray.setLinSpaced(7, 0, 10);

  cout << "select array is: " << selectArray << " with number of terms: " << selectArray.rows() << endl;

  Eigen::Array<float, Eigen::Dynamic, 1> selectArrayf;
  selectArrayf.setLinSpaced(7, 0, 10);
  cout << "select array is: " << selectArrayf << " with number of terms: " << selectArrayf.rows() << endl;

  selectArray = selectArrayf.round().cast<int>();

  cout << "select array is: " << selectArray << " with number of terms: " << selectArray.rows() << endl;

  cout << "array == 2 " << (selectArray == 2) << endl;
  cout << "(array == 2).any() " << (selectArray == 2).any() << endl;
  cout << "(array == 2).all() " << (selectArray == 2).all() << endl;

  Eigen::Array<bool, Eigen::Dynamic, 1> test(7);
  test << 0, 0, 1, 0, 1, 1, 1;// Initialise with commas
  Eigen::Array<bool, 1, Eigen::Dynamic> test2(7);
  test2 << 0, 0, 1, 0, 1, 1, 1;// Initialise with commas

  cout << "test array is: " << test << endl;

  downsampleRow(test, 3);
  
  cout << "test array after downsampling is: " << test << endl;

  downsampleCol(test2, 3);

  cout << "test2 array after downsampling to 3 is: " << test2 << endl;


  cloud->at(0,0).x = NAN;
  cloud->at(0,0).y = NAN;
  cloud->at(0,0).z = NAN;
  cout << "cloud point 0,0: " << cloud->at(0,0) << "\t is finite: " << pcl::isFinite(cloud->at(0,0)) << endl;

  Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> nanArray(cloud->height,cloud->width);
  nanArray.setZero(cloud->height,cloud->width);

  getNanMatrixFromPointCloud(nanArray, cloud);

  // cout << "Nan array output is: " << nanArray << "count" << nanArray.count() << endl;

  cout << "Number of Nans: " << nanArray.count() << endl;
}

int
main (int argc, char** argv)
{
  // Initialize ROS
  ros::init (argc, argv, "planar_seg_testing");
  ros::NodeHandle nh;

  pclTest();

  // Spin
  // ros::spin ();
}