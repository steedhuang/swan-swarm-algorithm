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
      PlayerClient    robot180("localhost", 7180);
      Position2dProxy p2dProxy_robot180(&robot180,0);
      RangerProxy      laserProxy_robot180(&robot180,0);
      FiducialProxy neighbor_finderProxy_robot180(&robot180,0);
      PlayerClient    robot181("localhost", 7181);
      Position2dProxy p2dProxy_robot181(&robot181,0);
      RangerProxy      laserProxy_robot181(&robot181,0);
      FiducialProxy neighbor_finderProxy_robot181(&robot181,0);
      PlayerClient    robot182("localhost", 7182);
      Position2dProxy p2dProxy_robot182(&robot182,0);
      RangerProxy      laserProxy_robot182(&robot182,0);
      FiducialProxy neighbor_finderProxy_robot182(&robot182,0);
      PlayerClient    robot183("localhost", 7183);
      Position2dProxy p2dProxy_robot183(&robot183,0);
      RangerProxy      laserProxy_robot183(&robot183,0);
      FiducialProxy neighbor_finderProxy_robot183(&robot183,0);
      PlayerClient    robot184("localhost", 7184);
      Position2dProxy p2dProxy_robot184(&robot184,0);
      RangerProxy      laserProxy_robot184(&robot184,0);
      FiducialProxy neighbor_finderProxy_robot184(&robot184,0);
      PlayerClient    robot185("localhost", 7185);
      Position2dProxy p2dProxy_robot185(&robot185,0);
      RangerProxy      laserProxy_robot185(&robot185,0);
      FiducialProxy neighbor_finderProxy_robot185(&robot185,0);
      PlayerClient    robot186("localhost", 7186);
      Position2dProxy p2dProxy_robot186(&robot186,0);
      RangerProxy      laserProxy_robot186(&robot186,0);
      FiducialProxy neighbor_finderProxy_robot186(&robot186,0);
      PlayerClient    robot187("localhost", 7187);
      Position2dProxy p2dProxy_robot187(&robot187,0);
      RangerProxy      laserProxy_robot187(&robot187,0);
      FiducialProxy neighbor_finderProxy_robot187(&robot187,0);
      PlayerClient    robot188("localhost", 7188);
      Position2dProxy p2dProxy_robot188(&robot188,0);
      RangerProxy      laserProxy_robot188(&robot188,0);
      FiducialProxy neighbor_finderProxy_robot188(&robot188,0);
      PlayerClient    robot189("localhost", 7189);
      Position2dProxy p2dProxy_robot189(&robot189,0);
      RangerProxy      laserProxy_robot189(&robot189,0);
      FiducialProxy neighbor_finderProxy_robot189(&robot189,0);
      PlayerClient    robot190("localhost", 7190);
      Position2dProxy p2dProxy_robot190(&robot190,0);
      RangerProxy      laserProxy_robot190(&robot190,0);
      FiducialProxy neighbor_finderProxy_robot190(&robot190,0);
      PlayerClient    robot191("localhost", 7191);
      Position2dProxy p2dProxy_robot191(&robot191,0);
      RangerProxy      laserProxy_robot191(&robot191,0);
      FiducialProxy neighbor_finderProxy_robot191(&robot191,0);
      PlayerClient    robot192("localhost", 7192);
      Position2dProxy p2dProxy_robot192(&robot192,0);
      RangerProxy      laserProxy_robot192(&robot192,0);
      FiducialProxy neighbor_finderProxy_robot192(&robot192,0);
      PlayerClient    robot193("localhost", 7193);
      Position2dProxy p2dProxy_robot193(&robot193,0);
      RangerProxy      laserProxy_robot193(&robot193,0);
      FiducialProxy neighbor_finderProxy_robot193(&robot193,0);
      PlayerClient    robot194("localhost", 7194);
      Position2dProxy p2dProxy_robot194(&robot194,0);
      RangerProxy      laserProxy_robot194(&robot194,0);
      FiducialProxy neighbor_finderProxy_robot194(&robot194,0);
      PlayerClient    robot195("localhost", 7195);
      Position2dProxy p2dProxy_robot195(&robot195,0);
      RangerProxy      laserProxy_robot195(&robot195,0);
      FiducialProxy neighbor_finderProxy_robot195(&robot195,0);
      PlayerClient    robot196("localhost", 7196);
      Position2dProxy p2dProxy_robot196(&robot196,0);
      RangerProxy      laserProxy_robot196(&robot196,0);
      FiducialProxy neighbor_finderProxy_robot196(&robot196,0);
      PlayerClient    robot197("localhost", 7197);
      Position2dProxy p2dProxy_robot197(&robot197,0);
      RangerProxy      laserProxy_robot197(&robot197,0);
      FiducialProxy neighbor_finderProxy_robot197(&robot197,0);
      PlayerClient    robot198("localhost", 7198);
      Position2dProxy p2dProxy_robot198(&robot198,0);
      RangerProxy      laserProxy_robot198(&robot198,0);
      FiducialProxy neighbor_finderProxy_robot198(&robot198,0);
      PlayerClient    robot199("localhost", 7199);
      Position2dProxy p2dProxy_robot199(&robot199,0);
      RangerProxy      laserProxy_robot199(&robot199,0);
      FiducialProxy neighbor_finderProxy_robot199(&robot199,0);
      PlayerClient    robot200("localhost", 7200);
      Position2dProxy p2dProxy_robot200(&robot200,0);
      RangerProxy      laserProxy_robot200(&robot200,0);
      FiducialProxy neighbor_finderProxy_robot200(&robot200,0);
      PlayerClient    robot201("localhost", 7201);
      Position2dProxy p2dProxy_robot201(&robot201,0);
      RangerProxy      laserProxy_robot201(&robot201,0);
      FiducialProxy neighbor_finderProxy_robot201(&robot201,0);
      PlayerClient    robot202("localhost", 7202);
      Position2dProxy p2dProxy_robot202(&robot202,0);
      RangerProxy      laserProxy_robot202(&robot202,0);
      FiducialProxy neighbor_finderProxy_robot202(&robot202,0);
      PlayerClient    robot203("localhost", 7203);
      Position2dProxy p2dProxy_robot203(&robot203,0);
      RangerProxy      laserProxy_robot203(&robot203,0);
      FiducialProxy neighbor_finderProxy_robot203(&robot203,0);
      PlayerClient    robot204("localhost", 7204);
      Position2dProxy p2dProxy_robot204(&robot204,0);
      RangerProxy      laserProxy_robot204(&robot204,0);
      FiducialProxy neighbor_finderProxy_robot204(&robot204,0);
      PlayerClient    robot205("localhost", 7205);
      Position2dProxy p2dProxy_robot205(&robot205,0);
      RangerProxy      laserProxy_robot205(&robot205,0);
      FiducialProxy neighbor_finderProxy_robot205(&robot205,0);
      PlayerClient    robot206("localhost", 7206);
      Position2dProxy p2dProxy_robot206(&robot206,0);
      RangerProxy      laserProxy_robot206(&robot206,0);
      FiducialProxy neighbor_finderProxy_robot206(&robot206,0);
      PlayerClient    robot207("localhost", 7207);
      Position2dProxy p2dProxy_robot207(&robot207,0);
      RangerProxy      laserProxy_robot207(&robot207,0);
      FiducialProxy neighbor_finderProxy_robot207(&robot207,0);
      PlayerClient    robot208("localhost", 7208);
      Position2dProxy p2dProxy_robot208(&robot208,0);
      RangerProxy      laserProxy_robot208(&robot208,0);
      FiducialProxy neighbor_finderProxy_robot208(&robot208,0);
      PlayerClient    robot209("localhost", 7209);
      Position2dProxy p2dProxy_robot209(&robot209,0);
      RangerProxy      laserProxy_robot209(&robot209,0);
      FiducialProxy neighbor_finderProxy_robot209(&robot209,0);
      PlayerClient    robot210("localhost", 7210);
      Position2dProxy p2dProxy_robot210(&robot210,0);
      RangerProxy      laserProxy_robot210(&robot210,0);
      FiducialProxy neighbor_finderProxy_robot210(&robot210,0);
      PlayerClient    robot211("localhost", 7211);
      Position2dProxy p2dProxy_robot211(&robot211,0);
      RangerProxy      laserProxy_robot211(&robot211,0);
      FiducialProxy neighbor_finderProxy_robot211(&robot211,0);
      PlayerClient    robot212("localhost", 7212);
      Position2dProxy p2dProxy_robot212(&robot212,0);
      RangerProxy      laserProxy_robot212(&robot212,0);
      FiducialProxy neighbor_finderProxy_robot212(&robot212,0);
      PlayerClient    robot213("localhost", 7213);
      Position2dProxy p2dProxy_robot213(&robot213,0);
      RangerProxy      laserProxy_robot213(&robot213,0);
      FiducialProxy neighbor_finderProxy_robot213(&robot213,0);
      PlayerClient    robot214("localhost", 7214);
      Position2dProxy p2dProxy_robot214(&robot214,0);
      RangerProxy      laserProxy_robot214(&robot214,0);
      FiducialProxy neighbor_finderProxy_robot214(&robot214,0);
      PlayerClient    robot215("localhost", 7215);
      Position2dProxy p2dProxy_robot215(&robot215,0);
      RangerProxy      laserProxy_robot215(&robot215,0);
      FiducialProxy neighbor_finderProxy_robot215(&robot215,0);
      PlayerClient    robot216("localhost", 7216);
      Position2dProxy p2dProxy_robot216(&robot216,0);
      RangerProxy      laserProxy_robot216(&robot216,0);
      FiducialProxy neighbor_finderProxy_robot216(&robot216,0);
      PlayerClient    robot217("localhost", 7217);
      Position2dProxy p2dProxy_robot217(&robot217,0);
      RangerProxy      laserProxy_robot217(&robot217,0);
      FiducialProxy neighbor_finderProxy_robot217(&robot217,0);
      PlayerClient    robot218("localhost", 7218);
      Position2dProxy p2dProxy_robot218(&robot218,0);
      RangerProxy      laserProxy_robot218(&robot218,0);
      FiducialProxy neighbor_finderProxy_robot218(&robot218,0);
      PlayerClient    robot219("localhost", 7219);
      Position2dProxy p2dProxy_robot219(&robot219,0);
      RangerProxy      laserProxy_robot219(&robot219,0);
      FiducialProxy neighbor_finderProxy_robot219(&robot219,0);
      PlayerClient    robot220("localhost", 7220);
      Position2dProxy p2dProxy_robot220(&robot220,0);
      RangerProxy      laserProxy_robot220(&robot220,0);
      FiducialProxy neighbor_finderProxy_robot220(&robot220,0);
      PlayerClient    robot221("localhost", 7221);
      Position2dProxy p2dProxy_robot221(&robot221,0);
      RangerProxy      laserProxy_robot221(&robot221,0);
      FiducialProxy neighbor_finderProxy_robot221(&robot221,0);
      PlayerClient    robot222("localhost", 7222);
      Position2dProxy p2dProxy_robot222(&robot222,0);
      RangerProxy      laserProxy_robot222(&robot222,0);
      FiducialProxy neighbor_finderProxy_robot222(&robot222,0);
      PlayerClient    robot223("localhost", 7223);
      Position2dProxy p2dProxy_robot223(&robot223,0);
      RangerProxy      laserProxy_robot223(&robot223,0);
      FiducialProxy neighbor_finderProxy_robot223(&robot223,0);
      PlayerClient    robot224("localhost", 7224);
      Position2dProxy p2dProxy_robot224(&robot224,0);
      RangerProxy      laserProxy_robot224(&robot224,0);
      FiducialProxy neighbor_finderProxy_robot224(&robot224,0);
      PlayerClient    robot225("localhost", 7225);
      Position2dProxy p2dProxy_robot225(&robot225,0);
      RangerProxy      laserProxy_robot225(&robot225,0);
      FiducialProxy neighbor_finderProxy_robot225(&robot225,0);
      PlayerClient    robot226("localhost", 7226);
      Position2dProxy p2dProxy_robot226(&robot226,0);
      RangerProxy      laserProxy_robot226(&robot226,0);
      FiducialProxy neighbor_finderProxy_robot226(&robot226,0);
      PlayerClient    robot227("localhost", 7227);
      Position2dProxy p2dProxy_robot227(&robot227,0);
      RangerProxy      laserProxy_robot227(&robot227,0);
      FiducialProxy neighbor_finderProxy_robot227(&robot227,0);
      PlayerClient    robot228("localhost", 7228);
      Position2dProxy p2dProxy_robot228(&robot228,0);
      RangerProxy      laserProxy_robot228(&robot228,0);
      FiducialProxy neighbor_finderProxy_robot228(&robot228,0);
      PlayerClient    robot229("localhost", 7229);
      Position2dProxy p2dProxy_robot229(&robot229,0);
      RangerProxy      laserProxy_robot229(&robot229,0);
      FiducialProxy neighbor_finderProxy_robot229(&robot229,0);
      PlayerClient    robot230("localhost", 7230);
      Position2dProxy p2dProxy_robot230(&robot230,0);
      RangerProxy      laserProxy_robot230(&robot230,0);
      FiducialProxy neighbor_finderProxy_robot230(&robot230,0);
      PlayerClient    robot231("localhost", 7231);
      Position2dProxy p2dProxy_robot231(&robot231,0);
      RangerProxy      laserProxy_robot231(&robot231,0);
      FiducialProxy neighbor_finderProxy_robot231(&robot231,0);
      PlayerClient    robot232("localhost", 7232);
      Position2dProxy p2dProxy_robot232(&robot232,0);
      RangerProxy      laserProxy_robot232(&robot232,0);
      FiducialProxy neighbor_finderProxy_robot232(&robot232,0);
      PlayerClient    robot233("localhost", 7233);
      Position2dProxy p2dProxy_robot233(&robot233,0);
      RangerProxy      laserProxy_robot233(&robot233,0);
      FiducialProxy neighbor_finderProxy_robot233(&robot233,0);
      PlayerClient    robot234("localhost", 7234);
      Position2dProxy p2dProxy_robot234(&robot234,0);
      RangerProxy      laserProxy_robot234(&robot234,0);
      FiducialProxy neighbor_finderProxy_robot234(&robot234,0);
      PlayerClient    robot235("localhost", 7235);
      Position2dProxy p2dProxy_robot235(&robot235,0);
      RangerProxy      laserProxy_robot235(&robot235,0);
      FiducialProxy neighbor_finderProxy_robot235(&robot235,0);
      PlayerClient    robot236("localhost", 7236);
      Position2dProxy p2dProxy_robot236(&robot236,0);
      RangerProxy      laserProxy_robot236(&robot236,0);
      FiducialProxy neighbor_finderProxy_robot236(&robot236,0);
      PlayerClient    robot237("localhost", 7237);
      Position2dProxy p2dProxy_robot237(&robot237,0);
      RangerProxy      laserProxy_robot237(&robot237,0);
      FiducialProxy neighbor_finderProxy_robot237(&robot237,0);
      PlayerClient    robot238("localhost", 7238);
      Position2dProxy p2dProxy_robot238(&robot238,0);
      RangerProxy      laserProxy_robot238(&robot238,0);
      FiducialProxy neighbor_finderProxy_robot238(&robot238,0);
      PlayerClient    robot239("localhost", 7239);
      Position2dProxy p2dProxy_robot239(&robot239,0);
      RangerProxy      laserProxy_robot239(&robot239,0);
      FiducialProxy neighbor_finderProxy_robot239(&robot239,0);

      // vector<PlayerClient> plyclnts(5);
      // vector<Position2dProxy> p2dProxys(5);
      // vector<RangerProxy> rngProxys(5);
      // vector<FiducialProxy> fidProxys(5);
      // vector<SimulationProxy> SimProxys(5);
// spawn off models
      PlayerClient plyclnts[SWARM_SIZE] ={robot180, robot181, robot182, robot183, robot184, robot185, robot186, robot187, robot188, robot189, robot190, robot191, robot192, robot193, robot194, robot195, robot196, robot197, robot198, robot199, robot200, robot201, robot202, robot203, robot204, robot205, robot206, robot207, robot208, robot209, robot210, robot211, robot212, robot213, robot214, robot215, robot216, robot217, robot218, robot219, robot220, robot221, robot222, robot223, robot224, robot225, robot226, robot227, robot228, robot229, robot230, robot231, robot232, robot233, robot234, robot235, robot236, robot237, robot238, robot239};
      Position2dProxy p2dProxys[SWARM_SIZE] = {p2dProxy_robot180, p2dProxy_robot181, p2dProxy_robot182, p2dProxy_robot183, p2dProxy_robot184, p2dProxy_robot185, p2dProxy_robot186, p2dProxy_robot187, p2dProxy_robot188, p2dProxy_robot189, p2dProxy_robot190, p2dProxy_robot191, p2dProxy_robot192, p2dProxy_robot193, p2dProxy_robot194, p2dProxy_robot195, p2dProxy_robot196, p2dProxy_robot197, p2dProxy_robot198, p2dProxy_robot199, p2dProxy_robot200, p2dProxy_robot201, p2dProxy_robot202, p2dProxy_robot203, p2dProxy_robot204, p2dProxy_robot205, p2dProxy_robot206, p2dProxy_robot207, p2dProxy_robot208, p2dProxy_robot209, p2dProxy_robot210, p2dProxy_robot211, p2dProxy_robot212, p2dProxy_robot213, p2dProxy_robot214, p2dProxy_robot215, p2dProxy_robot216, p2dProxy_robot217, p2dProxy_robot218, p2dProxy_robot219, p2dProxy_robot220, p2dProxy_robot221, p2dProxy_robot222, p2dProxy_robot223, p2dProxy_robot224, p2dProxy_robot225, p2dProxy_robot226, p2dProxy_robot227, p2dProxy_robot228, p2dProxy_robot229, p2dProxy_robot230, p2dProxy_robot231, p2dProxy_robot232, p2dProxy_robot233, p2dProxy_robot234, p2dProxy_robot235, p2dProxy_robot236, p2dProxy_robot237, p2dProxy_robot238, p2dProxy_robot239};
      RangerProxy rngProxys[SWARM_SIZE] = {laserProxy_robot180, laserProxy_robot181, laserProxy_robot182, laserProxy_robot183, laserProxy_robot184, laserProxy_robot185, laserProxy_robot186, laserProxy_robot187, laserProxy_robot188, laserProxy_robot189, laserProxy_robot190, laserProxy_robot191, laserProxy_robot192, laserProxy_robot193, laserProxy_robot194, laserProxy_robot195, laserProxy_robot196, laserProxy_robot197, laserProxy_robot198, laserProxy_robot199, laserProxy_robot200, laserProxy_robot201, laserProxy_robot202, laserProxy_robot203, laserProxy_robot204, laserProxy_robot205, laserProxy_robot206, laserProxy_robot207, laserProxy_robot208, laserProxy_robot209, laserProxy_robot210, laserProxy_robot211, laserProxy_robot212, laserProxy_robot213, laserProxy_robot214, laserProxy_robot215, laserProxy_robot216, laserProxy_robot217, laserProxy_robot218, laserProxy_robot219, laserProxy_robot220, laserProxy_robot221, laserProxy_robot222, laserProxy_robot223, laserProxy_robot224, laserProxy_robot225, laserProxy_robot226, laserProxy_robot227, laserProxy_robot228, laserProxy_robot229, laserProxy_robot230, laserProxy_robot231, laserProxy_robot232, laserProxy_robot233, laserProxy_robot234, laserProxy_robot235, laserProxy_robot236, laserProxy_robot237, laserProxy_robot238, laserProxy_robot239};
      FiducialProxy fidProxys[SWARM_SIZE]={neighbor_finderProxy_robot180, neighbor_finderProxy_robot181, neighbor_finderProxy_robot182, neighbor_finderProxy_robot183, neighbor_finderProxy_robot184, neighbor_finderProxy_robot185, neighbor_finderProxy_robot186, neighbor_finderProxy_robot187, neighbor_finderProxy_robot188, neighbor_finderProxy_robot189, neighbor_finderProxy_robot190, neighbor_finderProxy_robot191, neighbor_finderProxy_robot192, neighbor_finderProxy_robot193, neighbor_finderProxy_robot194, neighbor_finderProxy_robot195, neighbor_finderProxy_robot196, neighbor_finderProxy_robot197, neighbor_finderProxy_robot198, neighbor_finderProxy_robot199, neighbor_finderProxy_robot200, neighbor_finderProxy_robot201, neighbor_finderProxy_robot202, neighbor_finderProxy_robot203, neighbor_finderProxy_robot204, neighbor_finderProxy_robot205, neighbor_finderProxy_robot206, neighbor_finderProxy_robot207, neighbor_finderProxy_robot208, neighbor_finderProxy_robot209, neighbor_finderProxy_robot210, neighbor_finderProxy_robot211, neighbor_finderProxy_robot212, neighbor_finderProxy_robot213, neighbor_finderProxy_robot214, neighbor_finderProxy_robot215, neighbor_finderProxy_robot216, neighbor_finderProxy_robot217, neighbor_finderProxy_robot218, neighbor_finderProxy_robot219, neighbor_finderProxy_robot220, neighbor_finderProxy_robot221, neighbor_finderProxy_robot222, neighbor_finderProxy_robot223, neighbor_finderProxy_robot224, neighbor_finderProxy_robot225, neighbor_finderProxy_robot226, neighbor_finderProxy_robot227, neighbor_finderProxy_robot228, neighbor_finderProxy_robot229, neighbor_finderProxy_robot230, neighbor_finderProxy_robot231, neighbor_finderProxy_robot232, neighbor_finderProxy_robot233, neighbor_finderProxy_robot234, neighbor_finderProxy_robot235, neighbor_finderProxy_robot236, neighbor_finderProxy_robot237, neighbor_finderProxy_robot238, neighbor_finderProxy_robot239};
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
