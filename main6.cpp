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
      PlayerClient    robot300("localhost", 7300);
      Position2dProxy p2dProxy_robot300(&robot300,0);
      RangerProxy      laserProxy_robot300(&robot300,0);
      FiducialProxy neighbor_finderProxy_robot300(&robot300,0);
      PlayerClient    robot301("localhost", 7301);
      Position2dProxy p2dProxy_robot301(&robot301,0);
      RangerProxy      laserProxy_robot301(&robot301,0);
      FiducialProxy neighbor_finderProxy_robot301(&robot301,0);
      PlayerClient    robot302("localhost", 7302);
      Position2dProxy p2dProxy_robot302(&robot302,0);
      RangerProxy      laserProxy_robot302(&robot302,0);
      FiducialProxy neighbor_finderProxy_robot302(&robot302,0);
      PlayerClient    robot303("localhost", 7303);
      Position2dProxy p2dProxy_robot303(&robot303,0);
      RangerProxy      laserProxy_robot303(&robot303,0);
      FiducialProxy neighbor_finderProxy_robot303(&robot303,0);
      PlayerClient    robot304("localhost", 7304);
      Position2dProxy p2dProxy_robot304(&robot304,0);
      RangerProxy      laserProxy_robot304(&robot304,0);
      FiducialProxy neighbor_finderProxy_robot304(&robot304,0);
      PlayerClient    robot305("localhost", 7305);
      Position2dProxy p2dProxy_robot305(&robot305,0);
      RangerProxy      laserProxy_robot305(&robot305,0);
      FiducialProxy neighbor_finderProxy_robot305(&robot305,0);
      PlayerClient    robot306("localhost", 7306);
      Position2dProxy p2dProxy_robot306(&robot306,0);
      RangerProxy      laserProxy_robot306(&robot306,0);
      FiducialProxy neighbor_finderProxy_robot306(&robot306,0);
      PlayerClient    robot307("localhost", 7307);
      Position2dProxy p2dProxy_robot307(&robot307,0);
      RangerProxy      laserProxy_robot307(&robot307,0);
      FiducialProxy neighbor_finderProxy_robot307(&robot307,0);
      PlayerClient    robot308("localhost", 7308);
      Position2dProxy p2dProxy_robot308(&robot308,0);
      RangerProxy      laserProxy_robot308(&robot308,0);
      FiducialProxy neighbor_finderProxy_robot308(&robot308,0);
      PlayerClient    robot309("localhost", 7309);
      Position2dProxy p2dProxy_robot309(&robot309,0);
      RangerProxy      laserProxy_robot309(&robot309,0);
      FiducialProxy neighbor_finderProxy_robot309(&robot309,0);
      PlayerClient    robot310("localhost", 7310);
      Position2dProxy p2dProxy_robot310(&robot310,0);
      RangerProxy      laserProxy_robot310(&robot310,0);
      FiducialProxy neighbor_finderProxy_robot310(&robot310,0);
      PlayerClient    robot311("localhost", 7311);
      Position2dProxy p2dProxy_robot311(&robot311,0);
      RangerProxy      laserProxy_robot311(&robot311,0);
      FiducialProxy neighbor_finderProxy_robot311(&robot311,0);
      PlayerClient    robot312("localhost", 7312);
      Position2dProxy p2dProxy_robot312(&robot312,0);
      RangerProxy      laserProxy_robot312(&robot312,0);
      FiducialProxy neighbor_finderProxy_robot312(&robot312,0);
      PlayerClient    robot313("localhost", 7313);
      Position2dProxy p2dProxy_robot313(&robot313,0);
      RangerProxy      laserProxy_robot313(&robot313,0);
      FiducialProxy neighbor_finderProxy_robot313(&robot313,0);
      PlayerClient    robot314("localhost", 7314);
      Position2dProxy p2dProxy_robot314(&robot314,0);
      RangerProxy      laserProxy_robot314(&robot314,0);
      FiducialProxy neighbor_finderProxy_robot314(&robot314,0);
      PlayerClient    robot315("localhost", 7315);
      Position2dProxy p2dProxy_robot315(&robot315,0);
      RangerProxy      laserProxy_robot315(&robot315,0);
      FiducialProxy neighbor_finderProxy_robot315(&robot315,0);
      PlayerClient    robot316("localhost", 7316);
      Position2dProxy p2dProxy_robot316(&robot316,0);
      RangerProxy      laserProxy_robot316(&robot316,0);
      FiducialProxy neighbor_finderProxy_robot316(&robot316,0);
      PlayerClient    robot317("localhost", 7317);
      Position2dProxy p2dProxy_robot317(&robot317,0);
      RangerProxy      laserProxy_robot317(&robot317,0);
      FiducialProxy neighbor_finderProxy_robot317(&robot317,0);
      PlayerClient    robot318("localhost", 7318);
      Position2dProxy p2dProxy_robot318(&robot318,0);
      RangerProxy      laserProxy_robot318(&robot318,0);
      FiducialProxy neighbor_finderProxy_robot318(&robot318,0);
      PlayerClient    robot319("localhost", 7319);
      Position2dProxy p2dProxy_robot319(&robot319,0);
      RangerProxy      laserProxy_robot319(&robot319,0);
      FiducialProxy neighbor_finderProxy_robot319(&robot319,0);
      PlayerClient    robot320("localhost", 7320);
      Position2dProxy p2dProxy_robot320(&robot320,0);
      RangerProxy      laserProxy_robot320(&robot320,0);
      FiducialProxy neighbor_finderProxy_robot320(&robot320,0);
      PlayerClient    robot321("localhost", 7321);
      Position2dProxy p2dProxy_robot321(&robot321,0);
      RangerProxy      laserProxy_robot321(&robot321,0);
      FiducialProxy neighbor_finderProxy_robot321(&robot321,0);
      PlayerClient    robot322("localhost", 7322);
      Position2dProxy p2dProxy_robot322(&robot322,0);
      RangerProxy      laserProxy_robot322(&robot322,0);
      FiducialProxy neighbor_finderProxy_robot322(&robot322,0);
      PlayerClient    robot323("localhost", 7323);
      Position2dProxy p2dProxy_robot323(&robot323,0);
      RangerProxy      laserProxy_robot323(&robot323,0);
      FiducialProxy neighbor_finderProxy_robot323(&robot323,0);
      PlayerClient    robot324("localhost", 7324);
      Position2dProxy p2dProxy_robot324(&robot324,0);
      RangerProxy      laserProxy_robot324(&robot324,0);
      FiducialProxy neighbor_finderProxy_robot324(&robot324,0);
      PlayerClient    robot325("localhost", 7325);
      Position2dProxy p2dProxy_robot325(&robot325,0);
      RangerProxy      laserProxy_robot325(&robot325,0);
      FiducialProxy neighbor_finderProxy_robot325(&robot325,0);
      PlayerClient    robot326("localhost", 7326);
      Position2dProxy p2dProxy_robot326(&robot326,0);
      RangerProxy      laserProxy_robot326(&robot326,0);
      FiducialProxy neighbor_finderProxy_robot326(&robot326,0);
      PlayerClient    robot327("localhost", 7327);
      Position2dProxy p2dProxy_robot327(&robot327,0);
      RangerProxy      laserProxy_robot327(&robot327,0);
      FiducialProxy neighbor_finderProxy_robot327(&robot327,0);
      PlayerClient    robot328("localhost", 7328);
      Position2dProxy p2dProxy_robot328(&robot328,0);
      RangerProxy      laserProxy_robot328(&robot328,0);
      FiducialProxy neighbor_finderProxy_robot328(&robot328,0);
      PlayerClient    robot329("localhost", 7329);
      Position2dProxy p2dProxy_robot329(&robot329,0);
      RangerProxy      laserProxy_robot329(&robot329,0);
      FiducialProxy neighbor_finderProxy_robot329(&robot329,0);
      PlayerClient    robot330("localhost", 7330);
      Position2dProxy p2dProxy_robot330(&robot330,0);
      RangerProxy      laserProxy_robot330(&robot330,0);
      FiducialProxy neighbor_finderProxy_robot330(&robot330,0);
      PlayerClient    robot331("localhost", 7331);
      Position2dProxy p2dProxy_robot331(&robot331,0);
      RangerProxy      laserProxy_robot331(&robot331,0);
      FiducialProxy neighbor_finderProxy_robot331(&robot331,0);
      PlayerClient    robot332("localhost", 7332);
      Position2dProxy p2dProxy_robot332(&robot332,0);
      RangerProxy      laserProxy_robot332(&robot332,0);
      FiducialProxy neighbor_finderProxy_robot332(&robot332,0);
      PlayerClient    robot333("localhost", 7333);
      Position2dProxy p2dProxy_robot333(&robot333,0);
      RangerProxy      laserProxy_robot333(&robot333,0);
      FiducialProxy neighbor_finderProxy_robot333(&robot333,0);
      PlayerClient    robot334("localhost", 7334);
      Position2dProxy p2dProxy_robot334(&robot334,0);
      RangerProxy      laserProxy_robot334(&robot334,0);
      FiducialProxy neighbor_finderProxy_robot334(&robot334,0);
      PlayerClient    robot335("localhost", 7335);
      Position2dProxy p2dProxy_robot335(&robot335,0);
      RangerProxy      laserProxy_robot335(&robot335,0);
      FiducialProxy neighbor_finderProxy_robot335(&robot335,0);
      PlayerClient    robot336("localhost", 7336);
      Position2dProxy p2dProxy_robot336(&robot336,0);
      RangerProxy      laserProxy_robot336(&robot336,0);
      FiducialProxy neighbor_finderProxy_robot336(&robot336,0);
      PlayerClient    robot337("localhost", 7337);
      Position2dProxy p2dProxy_robot337(&robot337,0);
      RangerProxy      laserProxy_robot337(&robot337,0);
      FiducialProxy neighbor_finderProxy_robot337(&robot337,0);
      PlayerClient    robot338("localhost", 7338);
      Position2dProxy p2dProxy_robot338(&robot338,0);
      RangerProxy      laserProxy_robot338(&robot338,0);
      FiducialProxy neighbor_finderProxy_robot338(&robot338,0);
      PlayerClient    robot339("localhost", 7339);
      Position2dProxy p2dProxy_robot339(&robot339,0);
      RangerProxy      laserProxy_robot339(&robot339,0);
      FiducialProxy neighbor_finderProxy_robot339(&robot339,0);
      PlayerClient    robot340("localhost", 7340);
      Position2dProxy p2dProxy_robot340(&robot340,0);
      RangerProxy      laserProxy_robot340(&robot340,0);
      FiducialProxy neighbor_finderProxy_robot340(&robot340,0);
      PlayerClient    robot341("localhost", 7341);
      Position2dProxy p2dProxy_robot341(&robot341,0);
      RangerProxy      laserProxy_robot341(&robot341,0);
      FiducialProxy neighbor_finderProxy_robot341(&robot341,0);
      PlayerClient    robot342("localhost", 7342);
      Position2dProxy p2dProxy_robot342(&robot342,0);
      RangerProxy      laserProxy_robot342(&robot342,0);
      FiducialProxy neighbor_finderProxy_robot342(&robot342,0);
      PlayerClient    robot343("localhost", 7343);
      Position2dProxy p2dProxy_robot343(&robot343,0);
      RangerProxy      laserProxy_robot343(&robot343,0);
      FiducialProxy neighbor_finderProxy_robot343(&robot343,0);
      PlayerClient    robot344("localhost", 7344);
      Position2dProxy p2dProxy_robot344(&robot344,0);
      RangerProxy      laserProxy_robot344(&robot344,0);
      FiducialProxy neighbor_finderProxy_robot344(&robot344,0);
      PlayerClient    robot345("localhost", 7345);
      Position2dProxy p2dProxy_robot345(&robot345,0);
      RangerProxy      laserProxy_robot345(&robot345,0);
      FiducialProxy neighbor_finderProxy_robot345(&robot345,0);
      PlayerClient    robot346("localhost", 7346);
      Position2dProxy p2dProxy_robot346(&robot346,0);
      RangerProxy      laserProxy_robot346(&robot346,0);
      FiducialProxy neighbor_finderProxy_robot346(&robot346,0);
      PlayerClient    robot347("localhost", 7347);
      Position2dProxy p2dProxy_robot347(&robot347,0);
      RangerProxy      laserProxy_robot347(&robot347,0);
      FiducialProxy neighbor_finderProxy_robot347(&robot347,0);
      PlayerClient    robot348("localhost", 7348);
      Position2dProxy p2dProxy_robot348(&robot348,0);
      RangerProxy      laserProxy_robot348(&robot348,0);
      FiducialProxy neighbor_finderProxy_robot348(&robot348,0);
      PlayerClient    robot349("localhost", 7349);
      Position2dProxy p2dProxy_robot349(&robot349,0);
      RangerProxy      laserProxy_robot349(&robot349,0);
      FiducialProxy neighbor_finderProxy_robot349(&robot349,0);
      PlayerClient    robot350("localhost", 7350);
      Position2dProxy p2dProxy_robot350(&robot350,0);
      RangerProxy      laserProxy_robot350(&robot350,0);
      FiducialProxy neighbor_finderProxy_robot350(&robot350,0);
      PlayerClient    robot351("localhost", 7351);
      Position2dProxy p2dProxy_robot351(&robot351,0);
      RangerProxy      laserProxy_robot351(&robot351,0);
      FiducialProxy neighbor_finderProxy_robot351(&robot351,0);
      PlayerClient    robot352("localhost", 7352);
      Position2dProxy p2dProxy_robot352(&robot352,0);
      RangerProxy      laserProxy_robot352(&robot352,0);
      FiducialProxy neighbor_finderProxy_robot352(&robot352,0);
      PlayerClient    robot353("localhost", 7353);
      Position2dProxy p2dProxy_robot353(&robot353,0);
      RangerProxy      laserProxy_robot353(&robot353,0);
      FiducialProxy neighbor_finderProxy_robot353(&robot353,0);
      PlayerClient    robot354("localhost", 7354);
      Position2dProxy p2dProxy_robot354(&robot354,0);
      RangerProxy      laserProxy_robot354(&robot354,0);
      FiducialProxy neighbor_finderProxy_robot354(&robot354,0);
      PlayerClient    robot355("localhost", 7355);
      Position2dProxy p2dProxy_robot355(&robot355,0);
      RangerProxy      laserProxy_robot355(&robot355,0);
      FiducialProxy neighbor_finderProxy_robot355(&robot355,0);
      PlayerClient    robot356("localhost", 7356);
      Position2dProxy p2dProxy_robot356(&robot356,0);
      RangerProxy      laserProxy_robot356(&robot356,0);
      FiducialProxy neighbor_finderProxy_robot356(&robot356,0);
      PlayerClient    robot357("localhost", 7357);
      Position2dProxy p2dProxy_robot357(&robot357,0);
      RangerProxy      laserProxy_robot357(&robot357,0);
      FiducialProxy neighbor_finderProxy_robot357(&robot357,0);
      PlayerClient    robot358("localhost", 7358);
      Position2dProxy p2dProxy_robot358(&robot358,0);
      RangerProxy      laserProxy_robot358(&robot358,0);
      FiducialProxy neighbor_finderProxy_robot358(&robot358,0);
      PlayerClient    robot359("localhost", 7359);
      Position2dProxy p2dProxy_robot359(&robot359,0);
      RangerProxy      laserProxy_robot359(&robot359,0);
      FiducialProxy neighbor_finderProxy_robot359(&robot359,0);

      // vector<PlayerClient> plyclnts(5);
      // vector<Position2dProxy> p2dProxys(5);
      // vector<RangerProxy> rngProxys(5);
      // vector<FiducialProxy> fidProxys(5);
      // vector<SimulationProxy> SimProxys(5);
// spawn off models
      PlayerClient plyclnts[SWARM_SIZE] ={robot300, robot301, robot302, robot303, robot304, robot305, robot306, robot307, robot308, robot309, robot310, robot311, robot312, robot313, robot314, robot315, robot316, robot317, robot318, robot319, robot320, robot321, robot322, robot323, robot324, robot325, robot326, robot327, robot328, robot329, robot330, robot331, robot332, robot333, robot334, robot335, robot336, robot337, robot338, robot339, robot340, robot341, robot342, robot343, robot344, robot345, robot346, robot347, robot348, robot349, robot350, robot351, robot352, robot353, robot354, robot355, robot356, robot357, robot358, robot359};
      Position2dProxy p2dProxys[SWARM_SIZE] = {p2dProxy_robot300, p2dProxy_robot301, p2dProxy_robot302, p2dProxy_robot303, p2dProxy_robot304, p2dProxy_robot305, p2dProxy_robot306, p2dProxy_robot307, p2dProxy_robot308, p2dProxy_robot309, p2dProxy_robot310, p2dProxy_robot311, p2dProxy_robot312, p2dProxy_robot313, p2dProxy_robot314, p2dProxy_robot315, p2dProxy_robot316, p2dProxy_robot317, p2dProxy_robot318, p2dProxy_robot319, p2dProxy_robot320, p2dProxy_robot321, p2dProxy_robot322, p2dProxy_robot323, p2dProxy_robot324, p2dProxy_robot325, p2dProxy_robot326, p2dProxy_robot327, p2dProxy_robot328, p2dProxy_robot329, p2dProxy_robot330, p2dProxy_robot331, p2dProxy_robot332, p2dProxy_robot333, p2dProxy_robot334, p2dProxy_robot335, p2dProxy_robot336, p2dProxy_robot337, p2dProxy_robot338, p2dProxy_robot339, p2dProxy_robot340, p2dProxy_robot341, p2dProxy_robot342, p2dProxy_robot343, p2dProxy_robot344, p2dProxy_robot345, p2dProxy_robot346, p2dProxy_robot347, p2dProxy_robot348, p2dProxy_robot349, p2dProxy_robot350, p2dProxy_robot351, p2dProxy_robot352, p2dProxy_robot353, p2dProxy_robot354, p2dProxy_robot355, p2dProxy_robot356, p2dProxy_robot357, p2dProxy_robot358, p2dProxy_robot359};
      RangerProxy rngProxys[SWARM_SIZE] = {laserProxy_robot300, laserProxy_robot301, laserProxy_robot302, laserProxy_robot303, laserProxy_robot304, laserProxy_robot305, laserProxy_robot306, laserProxy_robot307, laserProxy_robot308, laserProxy_robot309, laserProxy_robot310, laserProxy_robot311, laserProxy_robot312, laserProxy_robot313, laserProxy_robot314, laserProxy_robot315, laserProxy_robot316, laserProxy_robot317, laserProxy_robot318, laserProxy_robot319, laserProxy_robot320, laserProxy_robot321, laserProxy_robot322, laserProxy_robot323, laserProxy_robot324, laserProxy_robot325, laserProxy_robot326, laserProxy_robot327, laserProxy_robot328, laserProxy_robot329, laserProxy_robot330, laserProxy_robot331, laserProxy_robot332, laserProxy_robot333, laserProxy_robot334, laserProxy_robot335, laserProxy_robot336, laserProxy_robot337, laserProxy_robot338, laserProxy_robot339, laserProxy_robot340, laserProxy_robot341, laserProxy_robot342, laserProxy_robot343, laserProxy_robot344, laserProxy_robot345, laserProxy_robot346, laserProxy_robot347, laserProxy_robot348, laserProxy_robot349, laserProxy_robot350, laserProxy_robot351, laserProxy_robot352, laserProxy_robot353, laserProxy_robot354, laserProxy_robot355, laserProxy_robot356, laserProxy_robot357, laserProxy_robot358, laserProxy_robot359};
      FiducialProxy fidProxys[SWARM_SIZE]={neighbor_finderProxy_robot300, neighbor_finderProxy_robot301, neighbor_finderProxy_robot302, neighbor_finderProxy_robot303, neighbor_finderProxy_robot304, neighbor_finderProxy_robot305, neighbor_finderProxy_robot306, neighbor_finderProxy_robot307, neighbor_finderProxy_robot308, neighbor_finderProxy_robot309, neighbor_finderProxy_robot310, neighbor_finderProxy_robot311, neighbor_finderProxy_robot312, neighbor_finderProxy_robot313, neighbor_finderProxy_robot314, neighbor_finderProxy_robot315, neighbor_finderProxy_robot316, neighbor_finderProxy_robot317, neighbor_finderProxy_robot318, neighbor_finderProxy_robot319, neighbor_finderProxy_robot320, neighbor_finderProxy_robot321, neighbor_finderProxy_robot322, neighbor_finderProxy_robot323, neighbor_finderProxy_robot324, neighbor_finderProxy_robot325, neighbor_finderProxy_robot326, neighbor_finderProxy_robot327, neighbor_finderProxy_robot328, neighbor_finderProxy_robot329, neighbor_finderProxy_robot330, neighbor_finderProxy_robot331, neighbor_finderProxy_robot332, neighbor_finderProxy_robot333, neighbor_finderProxy_robot334, neighbor_finderProxy_robot335, neighbor_finderProxy_robot336, neighbor_finderProxy_robot337, neighbor_finderProxy_robot338, neighbor_finderProxy_robot339, neighbor_finderProxy_robot340, neighbor_finderProxy_robot341, neighbor_finderProxy_robot342, neighbor_finderProxy_robot343, neighbor_finderProxy_robot344, neighbor_finderProxy_robot345, neighbor_finderProxy_robot346, neighbor_finderProxy_robot347, neighbor_finderProxy_robot348, neighbor_finderProxy_robot349, neighbor_finderProxy_robot350, neighbor_finderProxy_robot351, neighbor_finderProxy_robot352, neighbor_finderProxy_robot353, neighbor_finderProxy_robot354, neighbor_finderProxy_robot355, neighbor_finderProxy_robot356, neighbor_finderProxy_robot357, neighbor_finderProxy_robot358, neighbor_finderProxy_robot359};
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
