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
      PlayerClient    robot240("localhost", 7240);
      Position2dProxy p2dProxy_robot240(&robot240,0);
      RangerProxy      laserProxy_robot240(&robot240,0);
      FiducialProxy neighbor_finderProxy_robot240(&robot240,0);
      PlayerClient    robot241("localhost", 7241);
      Position2dProxy p2dProxy_robot241(&robot241,0);
      RangerProxy      laserProxy_robot241(&robot241,0);
      FiducialProxy neighbor_finderProxy_robot241(&robot241,0);
      PlayerClient    robot242("localhost", 7242);
      Position2dProxy p2dProxy_robot242(&robot242,0);
      RangerProxy      laserProxy_robot242(&robot242,0);
      FiducialProxy neighbor_finderProxy_robot242(&robot242,0);
      PlayerClient    robot243("localhost", 7243);
      Position2dProxy p2dProxy_robot243(&robot243,0);
      RangerProxy      laserProxy_robot243(&robot243,0);
      FiducialProxy neighbor_finderProxy_robot243(&robot243,0);
      PlayerClient    robot244("localhost", 7244);
      Position2dProxy p2dProxy_robot244(&robot244,0);
      RangerProxy      laserProxy_robot244(&robot244,0);
      FiducialProxy neighbor_finderProxy_robot244(&robot244,0);
      PlayerClient    robot245("localhost", 7245);
      Position2dProxy p2dProxy_robot245(&robot245,0);
      RangerProxy      laserProxy_robot245(&robot245,0);
      FiducialProxy neighbor_finderProxy_robot245(&robot245,0);
      PlayerClient    robot246("localhost", 7246);
      Position2dProxy p2dProxy_robot246(&robot246,0);
      RangerProxy      laserProxy_robot246(&robot246,0);
      FiducialProxy neighbor_finderProxy_robot246(&robot246,0);
      PlayerClient    robot247("localhost", 7247);
      Position2dProxy p2dProxy_robot247(&robot247,0);
      RangerProxy      laserProxy_robot247(&robot247,0);
      FiducialProxy neighbor_finderProxy_robot247(&robot247,0);
      PlayerClient    robot248("localhost", 7248);
      Position2dProxy p2dProxy_robot248(&robot248,0);
      RangerProxy      laserProxy_robot248(&robot248,0);
      FiducialProxy neighbor_finderProxy_robot248(&robot248,0);
      PlayerClient    robot249("localhost", 7249);
      Position2dProxy p2dProxy_robot249(&robot249,0);
      RangerProxy      laserProxy_robot249(&robot249,0);
      FiducialProxy neighbor_finderProxy_robot249(&robot249,0);
      PlayerClient    robot250("localhost", 7250);
      Position2dProxy p2dProxy_robot250(&robot250,0);
      RangerProxy      laserProxy_robot250(&robot250,0);
      FiducialProxy neighbor_finderProxy_robot250(&robot250,0);
      PlayerClient    robot251("localhost", 7251);
      Position2dProxy p2dProxy_robot251(&robot251,0);
      RangerProxy      laserProxy_robot251(&robot251,0);
      FiducialProxy neighbor_finderProxy_robot251(&robot251,0);
      PlayerClient    robot252("localhost", 7252);
      Position2dProxy p2dProxy_robot252(&robot252,0);
      RangerProxy      laserProxy_robot252(&robot252,0);
      FiducialProxy neighbor_finderProxy_robot252(&robot252,0);
      PlayerClient    robot253("localhost", 7253);
      Position2dProxy p2dProxy_robot253(&robot253,0);
      RangerProxy      laserProxy_robot253(&robot253,0);
      FiducialProxy neighbor_finderProxy_robot253(&robot253,0);
      PlayerClient    robot254("localhost", 7254);
      Position2dProxy p2dProxy_robot254(&robot254,0);
      RangerProxy      laserProxy_robot254(&robot254,0);
      FiducialProxy neighbor_finderProxy_robot254(&robot254,0);
      PlayerClient    robot255("localhost", 7255);
      Position2dProxy p2dProxy_robot255(&robot255,0);
      RangerProxy      laserProxy_robot255(&robot255,0);
      FiducialProxy neighbor_finderProxy_robot255(&robot255,0);
      PlayerClient    robot256("localhost", 7256);
      Position2dProxy p2dProxy_robot256(&robot256,0);
      RangerProxy      laserProxy_robot256(&robot256,0);
      FiducialProxy neighbor_finderProxy_robot256(&robot256,0);
      PlayerClient    robot257("localhost", 7257);
      Position2dProxy p2dProxy_robot257(&robot257,0);
      RangerProxy      laserProxy_robot257(&robot257,0);
      FiducialProxy neighbor_finderProxy_robot257(&robot257,0);
      PlayerClient    robot258("localhost", 7258);
      Position2dProxy p2dProxy_robot258(&robot258,0);
      RangerProxy      laserProxy_robot258(&robot258,0);
      FiducialProxy neighbor_finderProxy_robot258(&robot258,0);
      PlayerClient    robot259("localhost", 7259);
      Position2dProxy p2dProxy_robot259(&robot259,0);
      RangerProxy      laserProxy_robot259(&robot259,0);
      FiducialProxy neighbor_finderProxy_robot259(&robot259,0);
      PlayerClient    robot260("localhost", 7260);
      Position2dProxy p2dProxy_robot260(&robot260,0);
      RangerProxy      laserProxy_robot260(&robot260,0);
      FiducialProxy neighbor_finderProxy_robot260(&robot260,0);
      PlayerClient    robot261("localhost", 7261);
      Position2dProxy p2dProxy_robot261(&robot261,0);
      RangerProxy      laserProxy_robot261(&robot261,0);
      FiducialProxy neighbor_finderProxy_robot261(&robot261,0);
      PlayerClient    robot262("localhost", 7262);
      Position2dProxy p2dProxy_robot262(&robot262,0);
      RangerProxy      laserProxy_robot262(&robot262,0);
      FiducialProxy neighbor_finderProxy_robot262(&robot262,0);
      PlayerClient    robot263("localhost", 7263);
      Position2dProxy p2dProxy_robot263(&robot263,0);
      RangerProxy      laserProxy_robot263(&robot263,0);
      FiducialProxy neighbor_finderProxy_robot263(&robot263,0);
      PlayerClient    robot264("localhost", 7264);
      Position2dProxy p2dProxy_robot264(&robot264,0);
      RangerProxy      laserProxy_robot264(&robot264,0);
      FiducialProxy neighbor_finderProxy_robot264(&robot264,0);
      PlayerClient    robot265("localhost", 7265);
      Position2dProxy p2dProxy_robot265(&robot265,0);
      RangerProxy      laserProxy_robot265(&robot265,0);
      FiducialProxy neighbor_finderProxy_robot265(&robot265,0);
      PlayerClient    robot266("localhost", 7266);
      Position2dProxy p2dProxy_robot266(&robot266,0);
      RangerProxy      laserProxy_robot266(&robot266,0);
      FiducialProxy neighbor_finderProxy_robot266(&robot266,0);
      PlayerClient    robot267("localhost", 7267);
      Position2dProxy p2dProxy_robot267(&robot267,0);
      RangerProxy      laserProxy_robot267(&robot267,0);
      FiducialProxy neighbor_finderProxy_robot267(&robot267,0);
      PlayerClient    robot268("localhost", 7268);
      Position2dProxy p2dProxy_robot268(&robot268,0);
      RangerProxy      laserProxy_robot268(&robot268,0);
      FiducialProxy neighbor_finderProxy_robot268(&robot268,0);
      PlayerClient    robot269("localhost", 7269);
      Position2dProxy p2dProxy_robot269(&robot269,0);
      RangerProxy      laserProxy_robot269(&robot269,0);
      FiducialProxy neighbor_finderProxy_robot269(&robot269,0);
      PlayerClient    robot270("localhost", 7270);
      Position2dProxy p2dProxy_robot270(&robot270,0);
      RangerProxy      laserProxy_robot270(&robot270,0);
      FiducialProxy neighbor_finderProxy_robot270(&robot270,0);
      PlayerClient    robot271("localhost", 7271);
      Position2dProxy p2dProxy_robot271(&robot271,0);
      RangerProxy      laserProxy_robot271(&robot271,0);
      FiducialProxy neighbor_finderProxy_robot271(&robot271,0);
      PlayerClient    robot272("localhost", 7272);
      Position2dProxy p2dProxy_robot272(&robot272,0);
      RangerProxy      laserProxy_robot272(&robot272,0);
      FiducialProxy neighbor_finderProxy_robot272(&robot272,0);
      PlayerClient    robot273("localhost", 7273);
      Position2dProxy p2dProxy_robot273(&robot273,0);
      RangerProxy      laserProxy_robot273(&robot273,0);
      FiducialProxy neighbor_finderProxy_robot273(&robot273,0);
      PlayerClient    robot274("localhost", 7274);
      Position2dProxy p2dProxy_robot274(&robot274,0);
      RangerProxy      laserProxy_robot274(&robot274,0);
      FiducialProxy neighbor_finderProxy_robot274(&robot274,0);
      PlayerClient    robot275("localhost", 7275);
      Position2dProxy p2dProxy_robot275(&robot275,0);
      RangerProxy      laserProxy_robot275(&robot275,0);
      FiducialProxy neighbor_finderProxy_robot275(&robot275,0);
      PlayerClient    robot276("localhost", 7276);
      Position2dProxy p2dProxy_robot276(&robot276,0);
      RangerProxy      laserProxy_robot276(&robot276,0);
      FiducialProxy neighbor_finderProxy_robot276(&robot276,0);
      PlayerClient    robot277("localhost", 7277);
      Position2dProxy p2dProxy_robot277(&robot277,0);
      RangerProxy      laserProxy_robot277(&robot277,0);
      FiducialProxy neighbor_finderProxy_robot277(&robot277,0);
      PlayerClient    robot278("localhost", 7278);
      Position2dProxy p2dProxy_robot278(&robot278,0);
      RangerProxy      laserProxy_robot278(&robot278,0);
      FiducialProxy neighbor_finderProxy_robot278(&robot278,0);
      PlayerClient    robot279("localhost", 7279);
      Position2dProxy p2dProxy_robot279(&robot279,0);
      RangerProxy      laserProxy_robot279(&robot279,0);
      FiducialProxy neighbor_finderProxy_robot279(&robot279,0);
      PlayerClient    robot280("localhost", 7280);
      Position2dProxy p2dProxy_robot280(&robot280,0);
      RangerProxy      laserProxy_robot280(&robot280,0);
      FiducialProxy neighbor_finderProxy_robot280(&robot280,0);
      PlayerClient    robot281("localhost", 7281);
      Position2dProxy p2dProxy_robot281(&robot281,0);
      RangerProxy      laserProxy_robot281(&robot281,0);
      FiducialProxy neighbor_finderProxy_robot281(&robot281,0);
      PlayerClient    robot282("localhost", 7282);
      Position2dProxy p2dProxy_robot282(&robot282,0);
      RangerProxy      laserProxy_robot282(&robot282,0);
      FiducialProxy neighbor_finderProxy_robot282(&robot282,0);
      PlayerClient    robot283("localhost", 7283);
      Position2dProxy p2dProxy_robot283(&robot283,0);
      RangerProxy      laserProxy_robot283(&robot283,0);
      FiducialProxy neighbor_finderProxy_robot283(&robot283,0);
      PlayerClient    robot284("localhost", 7284);
      Position2dProxy p2dProxy_robot284(&robot284,0);
      RangerProxy      laserProxy_robot284(&robot284,0);
      FiducialProxy neighbor_finderProxy_robot284(&robot284,0);
      PlayerClient    robot285("localhost", 7285);
      Position2dProxy p2dProxy_robot285(&robot285,0);
      RangerProxy      laserProxy_robot285(&robot285,0);
      FiducialProxy neighbor_finderProxy_robot285(&robot285,0);
      PlayerClient    robot286("localhost", 7286);
      Position2dProxy p2dProxy_robot286(&robot286,0);
      RangerProxy      laserProxy_robot286(&robot286,0);
      FiducialProxy neighbor_finderProxy_robot286(&robot286,0);
      PlayerClient    robot287("localhost", 7287);
      Position2dProxy p2dProxy_robot287(&robot287,0);
      RangerProxy      laserProxy_robot287(&robot287,0);
      FiducialProxy neighbor_finderProxy_robot287(&robot287,0);
      PlayerClient    robot288("localhost", 7288);
      Position2dProxy p2dProxy_robot288(&robot288,0);
      RangerProxy      laserProxy_robot288(&robot288,0);
      FiducialProxy neighbor_finderProxy_robot288(&robot288,0);
      PlayerClient    robot289("localhost", 7289);
      Position2dProxy p2dProxy_robot289(&robot289,0);
      RangerProxy      laserProxy_robot289(&robot289,0);
      FiducialProxy neighbor_finderProxy_robot289(&robot289,0);
      PlayerClient    robot290("localhost", 7290);
      Position2dProxy p2dProxy_robot290(&robot290,0);
      RangerProxy      laserProxy_robot290(&robot290,0);
      FiducialProxy neighbor_finderProxy_robot290(&robot290,0);
      PlayerClient    robot291("localhost", 7291);
      Position2dProxy p2dProxy_robot291(&robot291,0);
      RangerProxy      laserProxy_robot291(&robot291,0);
      FiducialProxy neighbor_finderProxy_robot291(&robot291,0);
      PlayerClient    robot292("localhost", 7292);
      Position2dProxy p2dProxy_robot292(&robot292,0);
      RangerProxy      laserProxy_robot292(&robot292,0);
      FiducialProxy neighbor_finderProxy_robot292(&robot292,0);
      PlayerClient    robot293("localhost", 7293);
      Position2dProxy p2dProxy_robot293(&robot293,0);
      RangerProxy      laserProxy_robot293(&robot293,0);
      FiducialProxy neighbor_finderProxy_robot293(&robot293,0);
      PlayerClient    robot294("localhost", 7294);
      Position2dProxy p2dProxy_robot294(&robot294,0);
      RangerProxy      laserProxy_robot294(&robot294,0);
      FiducialProxy neighbor_finderProxy_robot294(&robot294,0);
      PlayerClient    robot295("localhost", 7295);
      Position2dProxy p2dProxy_robot295(&robot295,0);
      RangerProxy      laserProxy_robot295(&robot295,0);
      FiducialProxy neighbor_finderProxy_robot295(&robot295,0);
      PlayerClient    robot296("localhost", 7296);
      Position2dProxy p2dProxy_robot296(&robot296,0);
      RangerProxy      laserProxy_robot296(&robot296,0);
      FiducialProxy neighbor_finderProxy_robot296(&robot296,0);
      PlayerClient    robot297("localhost", 7297);
      Position2dProxy p2dProxy_robot297(&robot297,0);
      RangerProxy      laserProxy_robot297(&robot297,0);
      FiducialProxy neighbor_finderProxy_robot297(&robot297,0);
      PlayerClient    robot298("localhost", 7298);
      Position2dProxy p2dProxy_robot298(&robot298,0);
      RangerProxy      laserProxy_robot298(&robot298,0);
      FiducialProxy neighbor_finderProxy_robot298(&robot298,0);
      PlayerClient    robot299("localhost", 7299);
      Position2dProxy p2dProxy_robot299(&robot299,0);
      RangerProxy      laserProxy_robot299(&robot299,0);
      FiducialProxy neighbor_finderProxy_robot299(&robot299,0);

      // vector<PlayerClient> plyclnts(5);
      // vector<Position2dProxy> p2dProxys(5);
      // vector<RangerProxy> rngProxys(5);
      // vector<FiducialProxy> fidProxys(5);
      // vector<SimulationProxy> SimProxys(5);
// spawn off models
      PlayerClient plyclnts[SWARM_SIZE] ={robot240, robot241, robot242, robot243, robot244, robot245, robot246, robot247, robot248, robot249, robot250, robot251, robot252, robot253, robot254, robot255, robot256, robot257, robot258, robot259, robot260, robot261, robot262, robot263, robot264, robot265, robot266, robot267, robot268, robot269, robot270, robot271, robot272, robot273, robot274, robot275, robot276, robot277, robot278, robot279, robot280, robot281, robot282, robot283, robot284, robot285, robot286, robot287, robot288, robot289, robot290, robot291, robot292, robot293, robot294, robot295, robot296, robot297, robot298, robot299};
      Position2dProxy p2dProxys[SWARM_SIZE] = {p2dProxy_robot240, p2dProxy_robot241, p2dProxy_robot242, p2dProxy_robot243, p2dProxy_robot244, p2dProxy_robot245, p2dProxy_robot246, p2dProxy_robot247, p2dProxy_robot248, p2dProxy_robot249, p2dProxy_robot250, p2dProxy_robot251, p2dProxy_robot252, p2dProxy_robot253, p2dProxy_robot254, p2dProxy_robot255, p2dProxy_robot256, p2dProxy_robot257, p2dProxy_robot258, p2dProxy_robot259, p2dProxy_robot260, p2dProxy_robot261, p2dProxy_robot262, p2dProxy_robot263, p2dProxy_robot264, p2dProxy_robot265, p2dProxy_robot266, p2dProxy_robot267, p2dProxy_robot268, p2dProxy_robot269, p2dProxy_robot270, p2dProxy_robot271, p2dProxy_robot272, p2dProxy_robot273, p2dProxy_robot274, p2dProxy_robot275, p2dProxy_robot276, p2dProxy_robot277, p2dProxy_robot278, p2dProxy_robot279, p2dProxy_robot280, p2dProxy_robot281, p2dProxy_robot282, p2dProxy_robot283, p2dProxy_robot284, p2dProxy_robot285, p2dProxy_robot286, p2dProxy_robot287, p2dProxy_robot288, p2dProxy_robot289, p2dProxy_robot290, p2dProxy_robot291, p2dProxy_robot292, p2dProxy_robot293, p2dProxy_robot294, p2dProxy_robot295, p2dProxy_robot296, p2dProxy_robot297, p2dProxy_robot298, p2dProxy_robot299};
      RangerProxy rngProxys[SWARM_SIZE] = {laserProxy_robot240, laserProxy_robot241, laserProxy_robot242, laserProxy_robot243, laserProxy_robot244, laserProxy_robot245, laserProxy_robot246, laserProxy_robot247, laserProxy_robot248, laserProxy_robot249, laserProxy_robot250, laserProxy_robot251, laserProxy_robot252, laserProxy_robot253, laserProxy_robot254, laserProxy_robot255, laserProxy_robot256, laserProxy_robot257, laserProxy_robot258, laserProxy_robot259, laserProxy_robot260, laserProxy_robot261, laserProxy_robot262, laserProxy_robot263, laserProxy_robot264, laserProxy_robot265, laserProxy_robot266, laserProxy_robot267, laserProxy_robot268, laserProxy_robot269, laserProxy_robot270, laserProxy_robot271, laserProxy_robot272, laserProxy_robot273, laserProxy_robot274, laserProxy_robot275, laserProxy_robot276, laserProxy_robot277, laserProxy_robot278, laserProxy_robot279, laserProxy_robot280, laserProxy_robot281, laserProxy_robot282, laserProxy_robot283, laserProxy_robot284, laserProxy_robot285, laserProxy_robot286, laserProxy_robot287, laserProxy_robot288, laserProxy_robot289, laserProxy_robot290, laserProxy_robot291, laserProxy_robot292, laserProxy_robot293, laserProxy_robot294, laserProxy_robot295, laserProxy_robot296, laserProxy_robot297, laserProxy_robot298, laserProxy_robot299};
      FiducialProxy fidProxys[SWARM_SIZE]={neighbor_finderProxy_robot240, neighbor_finderProxy_robot241, neighbor_finderProxy_robot242, neighbor_finderProxy_robot243, neighbor_finderProxy_robot244, neighbor_finderProxy_robot245, neighbor_finderProxy_robot246, neighbor_finderProxy_robot247, neighbor_finderProxy_robot248, neighbor_finderProxy_robot249, neighbor_finderProxy_robot250, neighbor_finderProxy_robot251, neighbor_finderProxy_robot252, neighbor_finderProxy_robot253, neighbor_finderProxy_robot254, neighbor_finderProxy_robot255, neighbor_finderProxy_robot256, neighbor_finderProxy_robot257, neighbor_finderProxy_robot258, neighbor_finderProxy_robot259, neighbor_finderProxy_robot260, neighbor_finderProxy_robot261, neighbor_finderProxy_robot262, neighbor_finderProxy_robot263, neighbor_finderProxy_robot264, neighbor_finderProxy_robot265, neighbor_finderProxy_robot266, neighbor_finderProxy_robot267, neighbor_finderProxy_robot268, neighbor_finderProxy_robot269, neighbor_finderProxy_robot270, neighbor_finderProxy_robot271, neighbor_finderProxy_robot272, neighbor_finderProxy_robot273, neighbor_finderProxy_robot274, neighbor_finderProxy_robot275, neighbor_finderProxy_robot276, neighbor_finderProxy_robot277, neighbor_finderProxy_robot278, neighbor_finderProxy_robot279, neighbor_finderProxy_robot280, neighbor_finderProxy_robot281, neighbor_finderProxy_robot282, neighbor_finderProxy_robot283, neighbor_finderProxy_robot284, neighbor_finderProxy_robot285, neighbor_finderProxy_robot286, neighbor_finderProxy_robot287, neighbor_finderProxy_robot288, neighbor_finderProxy_robot289, neighbor_finderProxy_robot290, neighbor_finderProxy_robot291, neighbor_finderProxy_robot292, neighbor_finderProxy_robot293, neighbor_finderProxy_robot294, neighbor_finderProxy_robot295, neighbor_finderProxy_robot296, neighbor_finderProxy_robot297, neighbor_finderProxy_robot298, neighbor_finderProxy_robot299};
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
