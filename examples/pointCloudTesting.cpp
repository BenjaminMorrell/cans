#include <ros/ros.h>
#include <iostream>

#include "cans/mapping3D.h"

#include <cmath>
// PCL specific includes
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>



using namespace std;
using namespace PLib;

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

void downsampleCol(Eigen::Array<bool, 1, Eigen::Dynamic>& colFlags, int numColsDesired){

  // Create linspaced array 
  Eigen::Array<float, Eigen::Dynamic, 1> selectArrayf;
  Eigen::Array<int, Eigen::Dynamic, 1> selectArray;
  selectArrayf.setLinSpaced(numColsDesired, 0, colFlags.count()-1);
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
      if (i > 1){
        if (pcl::isFinite(cloud->at(j,i-2))){
          // Step with same delta as neighbours
          average.getArray3fMap() += 2*cloud->at(j,i-1).getArray3fMap() - cloud->at(j,i-2).getArray3fMap();
        }else{
          average.getArray3fMap() += cloud->at(j,i-1).getArray3fMap();
        }

      }else{
        average.getArray3fMap() += cloud->at(j,i-1).getArray3fMap();
      }
      add_count ++;
    }    
  }

  // Forward in i
  if (i < cloud->height-1){
    if (pcl::isFinite(cloud->at(j,i+1))){
      if (i < cloud->height-2){
        if (pcl::isFinite(cloud->at(j,i+2))){
          average.getArray3fMap() += 2*cloud->at(j,i+1).getArray3fMap() - cloud->at(j,i+2).getArray3fMap();    
        }else{
          average.getArray3fMap() += cloud->at(j,i+1).getArray3fMap();
        }
      }else{
        average.getArray3fMap() += cloud->at(j,i+1).getArray3fMap();
      }
      add_count ++;
    }    
  }

  // Back in j
  if (j > 0){
    if (pcl::isFinite(cloud->at(j-1,i))){
      if (j > 1){
        if (pcl::isFinite(cloud->at(j-2,i))){
          average.getArray3fMap() += 2*cloud->at(j-1,i).getArray3fMap() - cloud->at(j-2,i).getArray3fMap();
        }else{
          average.getArray3fMap() += cloud->at(j-1,i).getArray3fMap();
        }
      }else{
        average.getArray3fMap() += cloud->at(j-1,i).getArray3fMap();
      }
      add_count ++;
    }    
  }

  // Forward in j
  if (j < cloud->width-1){
    if (pcl::isFinite(cloud->at(j+1,i))){
      if (j < cloud->width-2){
        if (pcl::isFinite(cloud->at(j+2,i))){
          average.getArray3fMap() += 2*cloud->at(j+1,i).getArray3fMap() - cloud->at(j+2,i).getArray3fMap();
        }else{
          average.getArray3fMap() += cloud->at(j+1,i).getArray3fMap();
        }
      }else{
        average.getArray3fMap() += cloud->at(j+1,i).getArray3fMap();
      }
      add_count ++;
    }    
  }

  // Divide by count to get the average
  average.getArray3fMap() /= (float)add_count;

  //update cloud
  cloud->at(j,i) = average;

}

void averageOutNans(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, Eigen::Array<int,2,Eigen::Dynamic> nanIndices){
  cout << "Number of cols: " << nanIndices.cols() << endl;
  for (int i = 0; i < nanIndices.cols(); i++){
    if (nanIndices(0,i) == -1){
      // Flag meaning that those terms are not valid (not used)
      break;
    }
    cout << "Before: " << cloud->at(nanIndices(1,i),nanIndices(0,i)) << endl;
    regionAverage(cloud,nanIndices(0,i),nanIndices(1,i));
    cout << "After:  " << cloud->at(nanIndices(1,i),nanIndices(0,i)) << endl;
  }

}

bool removeRowsNan(Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>& nanArray,Eigen::Array<bool,Eigen::Dynamic, 1>& rowFlags, int maxNanAllowed){

  Eigen::Array<int, Eigen::Dynamic, 1> rowNan = nanArray.rowwise().count().cast<int>();

  // Maximum number of Nans in a row
  int maxNan = rowNan.maxCoeff();

  int buffer = 5;

  cout << "Max NaN count: " << maxNan << "\tmax allowed is: " << maxNanAllowed << endl;

  // If the maximum number Nans is above the set limit
  if (maxNan > maxNanAllowed){
    for (int i = 0; i < rowNan.rows(); i++){
      // If this rows has equal to the maximum number of nans
      if (rowNan(i) >= maxNan-buffer){
        // Set the row values to zero in the array
        nanArray.row(i) = Eigen::Array<bool, 1, Eigen::Dynamic>::Zero(1,nanArray.cols());

        // Set the rowFlags value to zero
        rowFlags(i) = false;
      }
    }
    return true;
  }else{
    return false;
  }

}

bool removeColsNan(Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>& nanArray,Eigen::Array<bool, 1, Eigen::Dynamic>& colFlags, int maxNanAllowed){

  Eigen::Array<int, 1, Eigen::Dynamic> colNan = nanArray.colwise().count().cast<int>();

  // Maximum number of Nans in a col
  int maxNan = colNan.maxCoeff();

  int buffer = 5;
  // TODO: Make this a user setting - or a member variable of the class in which this function lies 

  cout << "Max NaN count: " << maxNan << "\tmax allowed is: " << maxNanAllowed << endl;

  // If the maximum number Nans is above the set limit
  if (maxNan > maxNanAllowed){
    for (int i = 0; i < colNan.cols(); i++){
      // If this cols has equal to the maximum number of nans
      if (colNan(i) >= maxNan-buffer){
        // Set the col values to zero in the array
        nanArray.col(i) = Eigen::Array<bool, Eigen::Dynamic, 1>::Zero(nanArray.rows(), 1);

        // Set the colFlags value to zero
        colFlags(i) = false;

      }
    }
    return true;
  }else{
    return false;
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
  //       cloud->at(j,i).x = 0.0;// NOTE: at(col,row) so input is at(j, i)
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

void downsampMeshFromScan(){
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

  std::cerr << "PointCloud before filtering: W: " << cloud->width << "\tH: " << cloud->height << "\t data points." << std::endl;

  // Get Nan array
  Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> nanArray(cloud->height,cloud->width);
  nanArray.setZero(cloud->height,cloud->width);

  cout << "nanarray number of rows: " << nanArray.rows() << "\t number of cols: " << nanArray.cols() << endl;

  getNanMatrixFromPointCloud(nanArray, cloud);

  cout << "Number of Nans: " << nanArray.count() << endl;

  // Settings
  int numRowsDesired = 45;
  int numColsDesired = 45;
  int maxNanAllowed = numRowsDesired/5; 

  // Initialise
  Eigen::Array<bool, Eigen::Dynamic, 1> rowFlags(cloud->height, 1); rowFlags.setOnes(cloud->height, 1);
  Eigen::Array<bool, 1, Eigen::Dynamic> colFlags(1, cloud->width); colFlags.setOnes(1, cloud->width);
  // PointCloud (uint32_t width_, uint32_t height_, const PointT& value_ = PointT ())
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloudOut (new pcl::PointCloud<pcl::PointXYZ>(numColsDesired,numRowsDesired,pcl::PointXYZ(0.0,0.0,0.0)));

  cout << "Size of rowFlags: " << rowFlags.rows() << endl;
  cout << "Size of colFlags: " << colFlags.cols() << endl;

  bool exitFlag = false;

  bool nansPresent = true;

  // Loop to remove Nans
  while (!exitFlag){
    nansPresent = removeRowsNan(nanArray, rowFlags, maxNanAllowed);

    if (nansPresent){
      nansPresent = removeColsNan(nanArray, colFlags, maxNanAllowed);
      if (nansPresent){
        // Check dimensions 
        cout << "Nans removed: Row count: " << rowFlags.count() << "\nCol Count: " << colFlags.count() << endl;
        if (rowFlags.count() <= numRowsDesired || colFlags.count() <= numColsDesired){
          exitFlag = true;
          cout << "Exiting because dimensions are too small" << endl;
        }
      }else{
        exitFlag = true;
        cout << "Exiting because there are few enough NaNs" << endl;
      }
    }else{
      exitFlag = true;
      cout << "Exiting because there are few enough NaNs" << endl;
    }
  }

  // Downsample 
  downsampleRow(rowFlags, numRowsDesired);
  downsampleCol(colFlags, numColsDesired);

  // cout << "Row flags are: " << rowFlags << endl;
  // cout << "Col flags are: " << colFlags << endl;
  cout << "Number of nans left: " << nanArray.count() << endl;

  // Extract point cloud
  int ii = 0;
  int jj = 0;
  int ijk = 0;
  Eigen::Array<int, 2, Eigen::Dynamic> nanIndices(2,nanArray.count());
  nanIndices.setConstant(2,nanArray.count(),-1);

  for (int i = 0; i < cloud->height; i++){
    // Reset jj index
    jj = 0; 
    // For each row selected in rowFlags
    if (rowFlags(i)){
      for (int j = 0; j < cloud->width; j++){
        // For each column selected in colFlags
        if (colFlags(j)){
          // Get the data from the cloud
          cloudOut->at(jj,ii) = cloud->at(j,i);

          // Store indices if the value is nan
          if (!pcl::isFinite(cloud->at(j,i))){
            nanIndices(0,ijk) = ii;
            nanIndices(1,ijk) = jj;
            ijk++; 
          }

          // Increment the new column index
          jj++;
        }
      }
      // Increment the new row index
      ii++;
    }
  }

  cout << "ii = " << ii << "\tjj = " << jj << endl;

  // NOTE: this matrix will not be full, as there has been downsampling - so we have less rows and columns  
  cout << "Nan indices are: " << nanIndices << endl;

  // Average Nans
  averageOutNans(cloudOut, nanIndices);

  // Write the downsampled version to disk
  pcl::PCDWriter writer;
  writer.write<pcl::PointXYZ> ("blender_downsampled.pcd", *cloudOut, false);

  Mapping3D mp;

  pcl::PointCloud<pcl::PointXYZ>::Ptr cloudOut2 (new pcl::PointCloud<pcl::PointXYZ>(numColsDesired,numRowsDesired,pcl::PointXYZ(0.0,0.0,0.0)));
  
  mp.meshFromScan(cloudOut2,cloud);
 
  writer.write<pcl::PointXYZ> ("blender_downsampled_2.pcd", *cloudOut2, false);



}

int
main (int argc, char** argv)
{
  // Initialize ROS
  ros::init (argc, argv, "planar_seg_testing");
  ros::NodeHandle nh;

  downsampMeshFromScan();
  // pclTest();

  // Spin
  // ros::spin ();
}