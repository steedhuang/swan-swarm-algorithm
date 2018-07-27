// Authors: Ce MA, Jun Steed HUANG
// Date: 2018-2-2 created, 2018-7-11 updated
// Swan Swarm Selabot avoid a few obstacle with Finsler Geometric Dynamical Window Approach

// Preprocess headers
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <libplayerc++/playerc++.h>
#include <vector>
#include <chrono>
#include <ctime>
#include <math.h>
#include <float.h>
//#include <values.h>

// A couple of constants
#define W_MIN -45 // Minimum of W, W is turning speed
#define W_MAX 45 // Maximum of W
#define X_MIN 0.1 // Minimum of X, X is running speed
#define X_MAX 1   // Maximum of X
#define W_SAMPLE_GRAIN 1 // Sample size of W
#define X_SAMPLE_GRAIN 0.1 // Sample size of X
#define ROBOT_R 0.75 // the logic robot Radius 
#define OBSTACLE_R 0.25 // the logic obstcale Radius
#define ALPHA 0.05 // the weight of heading in evaluation 0.05
#define BETA 0.2 // the weight of dist in evaluation 0.2 
#define GAMMA 0.1  // the weight of velocity in evaluation 0.1
#define ETA 0.1 // the weight of swarm behavior in evaluation 0.1
#define DIST_PARA 100 // the predefined dist value for no obstcale situation
#define TIME_GRAIN 0.1  // time interval size
#define LASER_RANGE 4 // laser probe reach
#define FLOOR_RANGE 8 // floor side reach
#define SWARM_SIZE 60 // number of robots

// Cartesian coordinates and Eular: yaw pitch roll
struct location
{
     double X;
     double Y;
     double YAW;
}
typedef robot_loc_t;
using namespace PlayerCc;
using namespace std;
// flag to quit the simulation
int flag = 0; // flag for shut down
// window 
double time_window;
// converting unit
double degree2radians(double degree){
      return degree*M_PI/180.0;
}
double radians2degree(double radian){
      return radian*180.0/M_PI;
}
// heading
double heading_evaluation(double Yaw){
      return 180-fabs(Yaw);
}
// distance
double dist_evaluation(double traj_end_x, double traj_end_y, RangerProxy &laser, int loc)
{     double dist = 100;
      double angle = loc - 90; // -90 to +90 view 
        //    std::cout<<"angle: "<<angle<<std::endl;  // debug 
        //     usleep(10);
      double ob_x = laser[loc]*cos(angle);
      double ob_y = laser[loc]*sin(angle);
        //    printf("laser[loc]: %f\n", laser[loc]);  // debug 
      double temp_dist = sqrt((ob_x-traj_end_x)*(ob_x-traj_end_x)+(ob_y-traj_end_y)*(ob_y-traj_end_y)) - ROBOT_R - OBSTACLE_R;
      return temp_dist;
}
// speed
double velocity_evaluation(double v){
      return fabs(v);
}
// watch for man
vector<int> obstacle_detection(RangerProxy &laser, double *min_dist)
{
      vector<int> obstacle_location;
      double min_distance = LASER_RANGE;
      for(int i=0;i<laser.GetRangeCount();i++)
      {
            if(laser[i]<LASER_RANGE)
            {     if(min_distance>laser[i])
                  {
                        min_distance = laser[i];
                  }
                  obstacle_location.push_back(i);
            }
      }
      *min_dist = min_distance;
      //      printf("min_distance: %f\n", min_distance);  // debug 
      return obstacle_location;
}
// where I am
vector<robot_loc_t> trajectory_calculation(double X_speed, double W_speed, double *x_end, double *y_end, double *Yaw_status, double time_window)
{     double current_time = 0.0;
      double current_xx = 0.0;
      double current_yy = 0.0;
      double current_yaw = 0.0;
      vector<robot_loc_t> trajs;
      while(current_time<time_window)
      {     current_time = current_time + TIME_GRAIN;
            current_xx = current_xx + X_speed*TIME_GRAIN*cos(degree2radians(fabs(current_yaw)));
            current_yy = current_yy + X_speed*TIME_GRAIN*sin(degree2radians(fabs(current_yaw)));
            current_yaw = current_yaw + TIME_GRAIN*W_speed;
            robot_loc_t lo;
            lo.X = current_xx;
            lo.Y = current_yy;
            lo.YAW = current_yaw;
            trajs.push_back(lo);
      }
      *x_end = current_xx;
      *y_end = current_yy;
      *Yaw_status = *Yaw_status + current_yaw;
      //      printf("current_yy: %f\n", current_yy);  // debug 
      return trajs;
}
/**
 * determine the anarchism-authoritism  distance away from all neighbors
 **/
double follow_neighbour(double x_end, double y_end,\
                        FiducialProxy &neighborFinder, int iD)
{     double total_dist_to_neighbor = 1.0;
      for(int i =0; i< neighborFinder.GetCount(); i++){
            double x_dist = neighborFinder.GetFiducialItem(i).pose.px - x_end;
            double y_dist = neighborFinder.GetFiducialItem(i).pose.py - y_end;
// Assume Hawk is faster than Swan, and Swan is faster than Sparrow, we have:
//            total_dist_to_neighbor = total_dist_to_neighbor + sqrt(x_dist*x_dist*x_dist*x_dist+y_dist*y_dist*y_dist*y_dist)/abs(x_dist*y_dist); // Hawk Finsler distance
            if (fmod(iD, 20) < 18 && fmod(iD, 20) > 1)
              total_dist_to_neighbor = total_dist_to_neighbor + sqrt(x_dist*x_dist*x_dist*x_dist+y_dist*y_dist*y_dist*y_dist-x_dist*y_dist*x_dist*y_dist); // Swan Finsler distance
            else
 total_dist_to_neighbor = total_dist_to_neighbor + sqrt(abs(x_dist*x_dist*x_dist-x_dist*x_dist*y_dist+x_dist*y_dist*y_dist-y_dist*y_dist*y_dist)); // Sparrow Finsler distance
      }
      return total_dist_to_neighbor;
}

// the DWA body
void DynamicWindowApproach(double *X_speed, double *W_speed, double time_window,\
                              RangerProxy &laser,vector<int> obstcale_locs,\
                              Position2dProxy &p2d,\
                              FiducialProxy &neighborFinder, int iD)
{// the sampling window size and settings
      int W_windowSize = (int)((W_MAX-W_MIN)/W_SAMPLE_GRAIN);
      int X_windowSize = (int)((X_MAX-X_MIN)/X_SAMPLE_GRAIN);
      vector<double> Ws(W_windowSize);
      vector<double> Xs(X_windowSize);
      vector<vector<double> > dists(W_windowSize);
      vector<vector<double> > headings(W_windowSize);
      vector<vector<double> > velocitys(W_windowSize);
      vector<vector<double> > swarm_behavior_eva(W_windowSize);
      vector<vector<double> > Evalu(W_windowSize);
      double optimal_eva = -DBL_MAX;
      double total_dist = 0.0;
      double total_heading = 0.0;
      double total_velocity = 0.0;
      double total_swarm_behavior_eva =0.0;
      double optimal_x_speed = *X_speed;
      double optimal_w_speed = *W_speed;
      // printf("\nwindow size: (%d, %d)\n",W_windowSize, X_windowSize);
      // generating all possible v and w
      for(int i=0;i<W_windowSize; i++){
            Ws[i]= W_MIN + i*W_SAMPLE_GRAIN;
      } 
      for(int j=0;j<X_windowSize; j++){
            Xs[j]= X_MIN + j*X_SAMPLE_GRAIN;
      }
      // evaluating all possible v and w
      for(int i=0;i<W_windowSize;i++){
            dists[i]= vector<double>(X_windowSize);
            headings[i] = vector<double>(X_windowSize);
            velocitys[i] = vector<double>(X_windowSize);
            swarm_behavior_eva[i] = vector<double>(X_windowSize);
            for(int j=0;j<X_windowSize;j++){
                  // first generate the trajectory
                  double xx_end = 0.0;
                  double yy_end = 0.0;
                  double yaw_end = radians2degree(p2d.GetYaw());
                  vector<robot_loc_t> trajectory;
                  trajectory = trajectory_calculation(Xs[j], Ws[i], &xx_end, &yy_end, &yaw_end, time_window);
                  double dist_ = DIST_PARA;
                  // distance calculation
                  for(int k=0;k<obstcale_locs.size();k++){
                        double dist_temp = dist_evaluation(xx_end, yy_end, laser, obstcale_locs[k]);
                        if(dist_>dist_temp){
                            dist_ =  dist_temp;
                        }
                  }
                  if(dist_>10*OBSTACLE_R){
                        dist_ = 10*OBSTACLE_R;        // change 2 to 10
                  }
                  if(dist_>1*OBSTACLE_R){
                        dists[i][j] = dist_;
                        total_dist = total_dist + dist_;
                  }else{
                        dists[i][j] = -1;
                  }
                  // heading calculation
                  headings[i][j] = heading_evaluation(yaw_end);
                  // printf("\nhead end: %f\n", headings[i][j]);
                  total_heading = total_heading + headings[i][j];
                  // velocity calculation
                  velocitys[i][j] = velocity_evaluation(Xs[j]);
                  total_velocity = total_velocity + velocitys[i][j];
                  // fixing the middle man
                  swarm_behavior_eva[i][j] = follow_neighbour(xx_end, yy_end, neighborFinder, iD)*fmod(iD+10, 20);
                  total_swarm_behavior_eva = total_swarm_behavior_eva + swarm_behavior_eva[i][j];
                   printf("swarm: %f\n",swarm_behavior_eva[i][j] ); // debug
            }
      }
      // printf("\ntotal head: %f\n", total_heading);
      // find maximum evaluation
      for(int i=0;i<W_windowSize;i++){
            Evalu[i] = vector<double>(X_windowSize);
            for(int j=0;j<X_windowSize;j++){
                  // do normalization to reduce outlier effects
                  Evalu[i][j] = ALPHA*(headings[i][j]/total_heading)+\
                                    BETA*(dists[i][j]/total_dist)+\
                                    GAMMA*(velocitys[i][j]/total_velocity)+\
                                    ETA*(swarm_behavior_eva[i][j]/total_swarm_behavior_eva);
                  if(Evalu[i][j]>optimal_eva&&dists[i][j]!=-1){
                        optimal_eva = Evalu[i][j];
                        optimal_x_speed = Xs[j];
                        optimal_w_speed = Ws[i];
                        //  printf("eva optima: %f", optimal_eva);
                  }
            }
      }
      *X_speed = optimal_x_speed;
      *W_speed = optimal_w_speed;
      // printf("current X: %f ", optimal_x_speed);
      // printf("current W: %f\n", optimal_w_speed);

}
int main(int argc, char *argv[])
{	
      /*need to do this line in c++ only*/
      using namespace PlayerCc;
//===================================INITIALIZATION==================================================//
      PlayerClient    robot60("localhost", 7060);
      Position2dProxy p2dProxy_robot60(&robot60,0);
      RangerProxy      laserProxy_robot60(&robot60,0);
      FiducialProxy neighbor_finderProxy_robot60(&robot60,0);
      PlayerClient    robot61("localhost", 7061);
      Position2dProxy p2dProxy_robot61(&robot61,0);
      RangerProxy      laserProxy_robot61(&robot61,0);
      FiducialProxy neighbor_finderProxy_robot61(&robot61,0);
      PlayerClient    robot62("localhost", 7062);
      Position2dProxy p2dProxy_robot62(&robot62,0);
      RangerProxy      laserProxy_robot62(&robot62,0);
      FiducialProxy neighbor_finderProxy_robot62(&robot62,0);
      PlayerClient    robot63("localhost", 7063);
      Position2dProxy p2dProxy_robot63(&robot63,0);
      RangerProxy      laserProxy_robot63(&robot63,0);
      FiducialProxy neighbor_finderProxy_robot63(&robot63,0);
      PlayerClient    robot64("localhost", 7064);
      Position2dProxy p2dProxy_robot64(&robot64,0);
      RangerProxy      laserProxy_robot64(&robot64,0);
      FiducialProxy neighbor_finderProxy_robot64(&robot64,0);
      PlayerClient    robot65("localhost", 7065);
      Position2dProxy p2dProxy_robot65(&robot65,0);
      RangerProxy      laserProxy_robot65(&robot65,0);
      FiducialProxy neighbor_finderProxy_robot65(&robot65,0);
      PlayerClient    robot66("localhost", 7066);
      Position2dProxy p2dProxy_robot66(&robot66,0);
      RangerProxy      laserProxy_robot66(&robot66,0);
      FiducialProxy neighbor_finderProxy_robot66(&robot66,0);
      PlayerClient    robot67("localhost", 7067);
      Position2dProxy p2dProxy_robot67(&robot67,0);
      RangerProxy      laserProxy_robot67(&robot67,0);
      FiducialProxy neighbor_finderProxy_robot67(&robot67,0);
      PlayerClient    robot68("localhost", 7068);
      Position2dProxy p2dProxy_robot68(&robot68,0);
      RangerProxy      laserProxy_robot68(&robot68,0);
      FiducialProxy neighbor_finderProxy_robot68(&robot68,0);
      PlayerClient    robot69("localhost", 7069);
      Position2dProxy p2dProxy_robot69(&robot69,0);
      RangerProxy      laserProxy_robot69(&robot69,0);
      FiducialProxy neighbor_finderProxy_robot69(&robot69,0);
      PlayerClient    robot70("localhost", 7070);
      Position2dProxy p2dProxy_robot70(&robot70,0);
      RangerProxy      laserProxy_robot70(&robot70,0);
      FiducialProxy neighbor_finderProxy_robot70(&robot70,0);
      PlayerClient    robot71("localhost", 7071);
      Position2dProxy p2dProxy_robot71(&robot71,0);
      RangerProxy      laserProxy_robot71(&robot71,0);
      FiducialProxy neighbor_finderProxy_robot71(&robot71,0);
      PlayerClient    robot72("localhost", 7072);
      Position2dProxy p2dProxy_robot72(&robot72,0);
      RangerProxy      laserProxy_robot72(&robot72,0);
      FiducialProxy neighbor_finderProxy_robot72(&robot72,0);
      PlayerClient    robot73("localhost", 7073);
      Position2dProxy p2dProxy_robot73(&robot73,0);
      RangerProxy      laserProxy_robot73(&robot73,0);
      FiducialProxy neighbor_finderProxy_robot73(&robot73,0);
      PlayerClient    robot74("localhost", 7074);
      Position2dProxy p2dProxy_robot74(&robot74,0);
      RangerProxy      laserProxy_robot74(&robot74,0);
      FiducialProxy neighbor_finderProxy_robot74(&robot74,0);
      PlayerClient    robot75("localhost", 7075);
      Position2dProxy p2dProxy_robot75(&robot75,0);
      RangerProxy      laserProxy_robot75(&robot75,0);
      FiducialProxy neighbor_finderProxy_robot75(&robot75,0);
      PlayerClient    robot76("localhost", 7076);
      Position2dProxy p2dProxy_robot76(&robot76,0);
      RangerProxy      laserProxy_robot76(&robot76,0);
      FiducialProxy neighbor_finderProxy_robot76(&robot76,0);
      PlayerClient    robot77("localhost", 7077);
      Position2dProxy p2dProxy_robot77(&robot77,0);
      RangerProxy      laserProxy_robot77(&robot77,0);
      FiducialProxy neighbor_finderProxy_robot77(&robot77,0);
      PlayerClient    robot78("localhost", 7078);
      Position2dProxy p2dProxy_robot78(&robot78,0);
      RangerProxy      laserProxy_robot78(&robot78,0);
      FiducialProxy neighbor_finderProxy_robot78(&robot78,0);
      PlayerClient    robot79("localhost", 7079);
      Position2dProxy p2dProxy_robot79(&robot79,0);
      RangerProxy      laserProxy_robot79(&robot79,0);
      FiducialProxy neighbor_finderProxy_robot79(&robot79,0);
      PlayerClient    robot80("localhost", 7080);
      Position2dProxy p2dProxy_robot80(&robot80,0);
      RangerProxy      laserProxy_robot80(&robot80,0);
      FiducialProxy neighbor_finderProxy_robot80(&robot80,0);
      PlayerClient    robot81("localhost", 7081);
      Position2dProxy p2dProxy_robot81(&robot81,0);
      RangerProxy      laserProxy_robot81(&robot81,0);
      FiducialProxy neighbor_finderProxy_robot81(&robot81,0);
      PlayerClient    robot82("localhost", 7082);
      Position2dProxy p2dProxy_robot82(&robot82,0);
      RangerProxy      laserProxy_robot82(&robot82,0);
      FiducialProxy neighbor_finderProxy_robot82(&robot82,0);
      PlayerClient    robot83("localhost", 7083);
      Position2dProxy p2dProxy_robot83(&robot83,0);
      RangerProxy      laserProxy_robot83(&robot83,0);
      FiducialProxy neighbor_finderProxy_robot83(&robot83,0);
      PlayerClient    robot84("localhost", 7084);
      Position2dProxy p2dProxy_robot84(&robot84,0);
      RangerProxy      laserProxy_robot84(&robot84,0);
      FiducialProxy neighbor_finderProxy_robot84(&robot84,0);
      PlayerClient    robot85("localhost", 7085);
      Position2dProxy p2dProxy_robot85(&robot85,0);
      RangerProxy      laserProxy_robot85(&robot85,0);
      FiducialProxy neighbor_finderProxy_robot85(&robot85,0);
      PlayerClient    robot86("localhost", 7086);
      Position2dProxy p2dProxy_robot86(&robot86,0);
      RangerProxy      laserProxy_robot86(&robot86,0);
      FiducialProxy neighbor_finderProxy_robot86(&robot86,0);
      PlayerClient    robot87("localhost", 7087);
      Position2dProxy p2dProxy_robot87(&robot87,0);
      RangerProxy      laserProxy_robot87(&robot87,0);
      FiducialProxy neighbor_finderProxy_robot87(&robot87,0);
      PlayerClient    robot88("localhost", 7088);
      Position2dProxy p2dProxy_robot88(&robot88,0);
      RangerProxy      laserProxy_robot88(&robot88,0);
      FiducialProxy neighbor_finderProxy_robot88(&robot88,0);
      PlayerClient    robot89("localhost", 7089);
      Position2dProxy p2dProxy_robot89(&robot89,0);
      RangerProxy      laserProxy_robot89(&robot89,0);
      FiducialProxy neighbor_finderProxy_robot89(&robot89,0);
      PlayerClient    robot90("localhost", 7090);
      Position2dProxy p2dProxy_robot90(&robot90,0);
      RangerProxy      laserProxy_robot90(&robot90,0);
      FiducialProxy neighbor_finderProxy_robot90(&robot90,0);
      PlayerClient    robot91("localhost", 7091);
      Position2dProxy p2dProxy_robot91(&robot91,0);
      RangerProxy      laserProxy_robot91(&robot91,0);
      FiducialProxy neighbor_finderProxy_robot91(&robot91,0);
      PlayerClient    robot92("localhost", 7092);
      Position2dProxy p2dProxy_robot92(&robot92,0);
      RangerProxy      laserProxy_robot92(&robot92,0);
      FiducialProxy neighbor_finderProxy_robot92(&robot92,0);
      PlayerClient    robot93("localhost", 7093);
      Position2dProxy p2dProxy_robot93(&robot93,0);
      RangerProxy      laserProxy_robot93(&robot93,0);
      FiducialProxy neighbor_finderProxy_robot93(&robot93,0);
      PlayerClient    robot94("localhost", 7094);
      Position2dProxy p2dProxy_robot94(&robot94,0);
      RangerProxy      laserProxy_robot94(&robot94,0);
      FiducialProxy neighbor_finderProxy_robot94(&robot94,0);
      PlayerClient    robot95("localhost", 7095);
      Position2dProxy p2dProxy_robot95(&robot95,0);
      RangerProxy      laserProxy_robot95(&robot95,0);
      FiducialProxy neighbor_finderProxy_robot95(&robot95,0);
      PlayerClient    robot96("localhost", 7096);
      Position2dProxy p2dProxy_robot96(&robot96,0);
      RangerProxy      laserProxy_robot96(&robot96,0);
      FiducialProxy neighbor_finderProxy_robot96(&robot96,0);
      PlayerClient    robot97("localhost", 7097);
      Position2dProxy p2dProxy_robot97(&robot97,0);
      RangerProxy      laserProxy_robot97(&robot97,0);
      FiducialProxy neighbor_finderProxy_robot97(&robot97,0);
      PlayerClient    robot98("localhost", 7098);
      Position2dProxy p2dProxy_robot98(&robot98,0);
      RangerProxy      laserProxy_robot98(&robot98,0);
      FiducialProxy neighbor_finderProxy_robot98(&robot98,0);
      PlayerClient    robot99("localhost", 7099);
      Position2dProxy p2dProxy_robot99(&robot99,0);
      RangerProxy      laserProxy_robot99(&robot99,0);
      FiducialProxy neighbor_finderProxy_robot99(&robot99,0);
      PlayerClient    robot100("localhost", 7100);
      Position2dProxy p2dProxy_robot100(&robot100,0);
      RangerProxy      laserProxy_robot100(&robot100,0);
      FiducialProxy neighbor_finderProxy_robot100(&robot100,0);
      PlayerClient    robot101("localhost", 7101);
      Position2dProxy p2dProxy_robot101(&robot101,0);
      RangerProxy      laserProxy_robot101(&robot101,0);
      FiducialProxy neighbor_finderProxy_robot101(&robot101,0);
      PlayerClient    robot102("localhost", 7102);
      Position2dProxy p2dProxy_robot102(&robot102,0);
      RangerProxy      laserProxy_robot102(&robot102,0);
      FiducialProxy neighbor_finderProxy_robot102(&robot102,0);
      PlayerClient    robot103("localhost", 7103);
      Position2dProxy p2dProxy_robot103(&robot103,0);
      RangerProxy      laserProxy_robot103(&robot103,0);
      FiducialProxy neighbor_finderProxy_robot103(&robot103,0);
      PlayerClient    robot104("localhost", 7104);
      Position2dProxy p2dProxy_robot104(&robot104,0);
      RangerProxy      laserProxy_robot104(&robot104,0);
      FiducialProxy neighbor_finderProxy_robot104(&robot104,0);
      PlayerClient    robot105("localhost", 7105);
      Position2dProxy p2dProxy_robot105(&robot105,0);
      RangerProxy      laserProxy_robot105(&robot105,0);
      FiducialProxy neighbor_finderProxy_robot105(&robot105,0);
      PlayerClient    robot106("localhost", 7106);
      Position2dProxy p2dProxy_robot106(&robot106,0);
      RangerProxy      laserProxy_robot106(&robot106,0);
      FiducialProxy neighbor_finderProxy_robot106(&robot106,0);
      PlayerClient    robot107("localhost", 7107);
      Position2dProxy p2dProxy_robot107(&robot107,0);
      RangerProxy      laserProxy_robot107(&robot107,0);
      FiducialProxy neighbor_finderProxy_robot107(&robot107,0);
      PlayerClient    robot108("localhost", 7108);
      Position2dProxy p2dProxy_robot108(&robot108,0);
      RangerProxy      laserProxy_robot108(&robot108,0);
      FiducialProxy neighbor_finderProxy_robot108(&robot108,0);
      PlayerClient    robot109("localhost", 7109);
      Position2dProxy p2dProxy_robot109(&robot109,0);
      RangerProxy      laserProxy_robot109(&robot109,0);
      FiducialProxy neighbor_finderProxy_robot109(&robot109,0);
      PlayerClient    robot110("localhost", 7110);
      Position2dProxy p2dProxy_robot110(&robot110,0);
      RangerProxy      laserProxy_robot110(&robot110,0);
      FiducialProxy neighbor_finderProxy_robot110(&robot110,0);
      PlayerClient    robot111("localhost", 7111);
      Position2dProxy p2dProxy_robot111(&robot111,0);
      RangerProxy      laserProxy_robot111(&robot111,0);
      FiducialProxy neighbor_finderProxy_robot111(&robot111,0);
      PlayerClient    robot112("localhost", 7112);
      Position2dProxy p2dProxy_robot112(&robot112,0);
      RangerProxy      laserProxy_robot112(&robot112,0);
      FiducialProxy neighbor_finderProxy_robot112(&robot112,0);
      PlayerClient    robot113("localhost", 7113);
      Position2dProxy p2dProxy_robot113(&robot113,0);
      RangerProxy      laserProxy_robot113(&robot113,0);
      FiducialProxy neighbor_finderProxy_robot113(&robot113,0);
      PlayerClient    robot114("localhost", 7114);
      Position2dProxy p2dProxy_robot114(&robot114,0);
      RangerProxy      laserProxy_robot114(&robot114,0);
      FiducialProxy neighbor_finderProxy_robot114(&robot114,0);
      PlayerClient    robot115("localhost", 7115);
      Position2dProxy p2dProxy_robot115(&robot115,0);
      RangerProxy      laserProxy_robot115(&robot115,0);
      FiducialProxy neighbor_finderProxy_robot115(&robot115,0);
      PlayerClient    robot116("localhost", 7116);
      Position2dProxy p2dProxy_robot116(&robot116,0);
      RangerProxy      laserProxy_robot116(&robot116,0);
      FiducialProxy neighbor_finderProxy_robot116(&robot116,0);
      PlayerClient    robot117("localhost", 7117);
      Position2dProxy p2dProxy_robot117(&robot117,0);
      RangerProxy      laserProxy_robot117(&robot117,0);
      FiducialProxy neighbor_finderProxy_robot117(&robot117,0);
      PlayerClient    robot118("localhost", 7118);
      Position2dProxy p2dProxy_robot118(&robot118,0);
      RangerProxy      laserProxy_robot118(&robot118,0);
      FiducialProxy neighbor_finderProxy_robot118(&robot118,0);
      PlayerClient    robot119("localhost", 7119);
      Position2dProxy p2dProxy_robot119(&robot119,0);
      RangerProxy      laserProxy_robot119(&robot119,0);
      FiducialProxy neighbor_finderProxy_robot119(&robot119,0);

      // vector<PlayerClient> plyclnts(5);
      // vector<Position2dProxy> p2dProxys(5);
      // vector<RangerProxy> rngProxys(5);
      // vector<FiducialProxy> fidProxys(5);
      // vector<SimulationProxy> SimProxys(5);
// spawn off models
      PlayerClient plyclnts[SWARM_SIZE] ={robot60, robot61, robot62, robot63, robot64, robot65, robot66, robot67, robot68, robot69, robot70, robot71, robot72, robot73, robot74, robot75, robot76, robot77, robot78, robot79, robot80, robot81, robot82, robot83, robot84, robot85, robot86, robot87, robot88, robot89, robot90, robot91, robot92, robot93, robot94, robot95, robot96, robot97, robot98, robot99, robot100, robot101, robot102, robot103, robot104, robot105, robot106, robot107, robot108, robot109, robot110, robot111, robot112, robot113, robot114, robot115, robot116, robot117, robot118, robot119};
      Position2dProxy p2dProxys[SWARM_SIZE] = {p2dProxy_robot60, p2dProxy_robot61, p2dProxy_robot62, p2dProxy_robot63, p2dProxy_robot64, p2dProxy_robot65, p2dProxy_robot66, p2dProxy_robot67, p2dProxy_robot68, p2dProxy_robot69, p2dProxy_robot70, p2dProxy_robot71, p2dProxy_robot72, p2dProxy_robot73, p2dProxy_robot74, p2dProxy_robot75, p2dProxy_robot76, p2dProxy_robot77, p2dProxy_robot78, p2dProxy_robot79, p2dProxy_robot80, p2dProxy_robot81, p2dProxy_robot82, p2dProxy_robot83, p2dProxy_robot84, p2dProxy_robot85, p2dProxy_robot86, p2dProxy_robot87, p2dProxy_robot88, p2dProxy_robot89, p2dProxy_robot90, p2dProxy_robot91, p2dProxy_robot92, p2dProxy_robot93, p2dProxy_robot94, p2dProxy_robot95, p2dProxy_robot96, p2dProxy_robot97, p2dProxy_robot98, p2dProxy_robot99, p2dProxy_robot100, p2dProxy_robot101, p2dProxy_robot102, p2dProxy_robot103, p2dProxy_robot104, p2dProxy_robot105, p2dProxy_robot106, p2dProxy_robot107, p2dProxy_robot108, p2dProxy_robot109, p2dProxy_robot110, p2dProxy_robot111, p2dProxy_robot112, p2dProxy_robot113, p2dProxy_robot114, p2dProxy_robot115, p2dProxy_robot116, p2dProxy_robot117, p2dProxy_robot118, p2dProxy_robot119};
      RangerProxy rngProxys[SWARM_SIZE] = {laserProxy_robot60, laserProxy_robot61, laserProxy_robot62, laserProxy_robot63, laserProxy_robot64, laserProxy_robot65, laserProxy_robot66, laserProxy_robot67, laserProxy_robot68, laserProxy_robot69, laserProxy_robot70, laserProxy_robot71, laserProxy_robot72, laserProxy_robot73, laserProxy_robot74, laserProxy_robot75, laserProxy_robot76, laserProxy_robot77, laserProxy_robot78, laserProxy_robot79, laserProxy_robot80, laserProxy_robot81, laserProxy_robot82, laserProxy_robot83, laserProxy_robot84, laserProxy_robot85, laserProxy_robot86, laserProxy_robot87, laserProxy_robot88, laserProxy_robot89, laserProxy_robot90, laserProxy_robot91, laserProxy_robot92, laserProxy_robot93, laserProxy_robot94, laserProxy_robot95, laserProxy_robot96, laserProxy_robot97, laserProxy_robot98, laserProxy_robot99, laserProxy_robot100, laserProxy_robot101, laserProxy_robot102, laserProxy_robot103, laserProxy_robot104, laserProxy_robot105, laserProxy_robot106, laserProxy_robot107, laserProxy_robot108, laserProxy_robot109, laserProxy_robot110, laserProxy_robot111, laserProxy_robot112, laserProxy_robot113, laserProxy_robot114, laserProxy_robot115, laserProxy_robot116, laserProxy_robot117, laserProxy_robot118, laserProxy_robot119};
      FiducialProxy fidProxys[SWARM_SIZE]={neighbor_finderProxy_robot60, neighbor_finderProxy_robot61, neighbor_finderProxy_robot62, neighbor_finderProxy_robot63, neighbor_finderProxy_robot64, neighbor_finderProxy_robot65, neighbor_finderProxy_robot66, neighbor_finderProxy_robot67, neighbor_finderProxy_robot68, neighbor_finderProxy_robot69, neighbor_finderProxy_robot70, neighbor_finderProxy_robot71, neighbor_finderProxy_robot72, neighbor_finderProxy_robot73, neighbor_finderProxy_robot74, neighbor_finderProxy_robot75, neighbor_finderProxy_robot76, neighbor_finderProxy_robot77, neighbor_finderProxy_robot78, neighbor_finderProxy_robot79, neighbor_finderProxy_robot80, neighbor_finderProxy_robot81, neighbor_finderProxy_robot82, neighbor_finderProxy_robot83, neighbor_finderProxy_robot84, neighbor_finderProxy_robot85, neighbor_finderProxy_robot86, neighbor_finderProxy_robot87, neighbor_finderProxy_robot88, neighbor_finderProxy_robot89, neighbor_finderProxy_robot90, neighbor_finderProxy_robot91, neighbor_finderProxy_robot92, neighbor_finderProxy_robot93, neighbor_finderProxy_robot94, neighbor_finderProxy_robot95, neighbor_finderProxy_robot96, neighbor_finderProxy_robot97, neighbor_finderProxy_robot98, neighbor_finderProxy_robot99, neighbor_finderProxy_robot100, neighbor_finderProxy_robot101, neighbor_finderProxy_robot102, neighbor_finderProxy_robot103, neighbor_finderProxy_robot104, neighbor_finderProxy_robot105, neighbor_finderProxy_robot106, neighbor_finderProxy_robot107, neighbor_finderProxy_robot108, neighbor_finderProxy_robot109, neighbor_finderProxy_robot110, neighbor_finderProxy_robot111, neighbor_finderProxy_robot112, neighbor_finderProxy_robot113, neighbor_finderProxy_robot114, neighbor_finderProxy_robot115, neighbor_finderProxy_robot116, neighbor_finderProxy_robot117, neighbor_finderProxy_robot118, neighbor_finderProxy_robot119};
      // SimulationProxy SimProxys[SWARM_SIZE]={simProxy_robot0, simProxy_robot1, simProxy_robot2, simProxy_robot3, simProxy_robot4, simProxy_robot5, simProxy_robot6, simProxy_robot7, simProxy_robot8, simProxy_robot9, simProxy_robot10, simProxy_robot11};
      vector<double> X_speeds(SWARM_SIZE);
      vector<double> W_speeds(SWARM_SIZE);
      // aiming each robots towards the goal
      for(int i=0;i<SWARM_SIZE;i++){
            X_speeds[i] = 0.25;
            W_speeds[i] = (SWARM_SIZE/2-i)*0.0017; // roughly 3.1415926/180
            p2dProxys[i].SetMotorEnable(1);
            p2dProxys[i].RequestGeom();
            rngProxys[i].RequestGeom();
            fidProxys[i].RequestGeometry();
            rngProxys[i].RequestConfigure();

      }
      // randomize the wheel speed
      srand(time(NULL));
      double Avoid_dist = LASER_RANGE - ROBOT_R*2; //this is twice size of length of robot
      double current_min_dist;
             current_min_dist = Avoid_dist;
             //     printf("current_min_dist: %f\n", current_min_dist);  // debug 
      time_window = 3.0;// every time_window to update once
//======================================LOOP=====================================================================
      while(true)
      {	
	    for(int i = 0;i<SWARM_SIZE;i++){
            // read from the proxies
            // robot.Read();
            plyclnts[i].Read();
            vector<int> locs;
            // Determine obstcales for each robot
            locs = obstacle_detection(rngProxys[i], &current_min_dist);
            DynamicWindowApproach(&X_speeds[i], &W_speeds[i], time_window, rngProxys[i], locs, p2dProxys[i], fidProxys[i], i);
            p2dProxys[i].SetSpeed(X_speeds[i], dtor(W_speeds[i]));
            // Stop when some at next destination
            if ( p2dProxys[i].GetXPos()>FLOOR_RANGE*0.1)   
		flag = flag+1;
            // To make sure no queue overflow
            sleep(time_window/SWARM_SIZE);
            }
            // More than Goldern 61.8% or Pareto 80% next step 10%
            if (flag > SWARM_SIZE*0.1)
            break;
            // sleep(time_window/10);
      }
}
