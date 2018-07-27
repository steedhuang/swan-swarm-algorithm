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
      PlayerClient    robot360("localhost", 7360);
      Position2dProxy p2dProxy_robot360(&robot360,0);
      RangerProxy      laserProxy_robot360(&robot360,0);
      FiducialProxy neighbor_finderProxy_robot360(&robot360,0);
      PlayerClient    robot361("localhost", 7361);
      Position2dProxy p2dProxy_robot361(&robot361,0);
      RangerProxy      laserProxy_robot361(&robot361,0);
      FiducialProxy neighbor_finderProxy_robot361(&robot361,0);
      PlayerClient    robot362("localhost", 7362);
      Position2dProxy p2dProxy_robot362(&robot362,0);
      RangerProxy      laserProxy_robot362(&robot362,0);
      FiducialProxy neighbor_finderProxy_robot362(&robot362,0);
      PlayerClient    robot363("localhost", 7363);
      Position2dProxy p2dProxy_robot363(&robot363,0);
      RangerProxy      laserProxy_robot363(&robot363,0);
      FiducialProxy neighbor_finderProxy_robot363(&robot363,0);
      PlayerClient    robot364("localhost", 7364);
      Position2dProxy p2dProxy_robot364(&robot364,0);
      RangerProxy      laserProxy_robot364(&robot364,0);
      FiducialProxy neighbor_finderProxy_robot364(&robot364,0);
      PlayerClient    robot365("localhost", 7365);
      Position2dProxy p2dProxy_robot365(&robot365,0);
      RangerProxy      laserProxy_robot365(&robot365,0);
      FiducialProxy neighbor_finderProxy_robot365(&robot365,0);
      PlayerClient    robot366("localhost", 7366);
      Position2dProxy p2dProxy_robot366(&robot366,0);
      RangerProxy      laserProxy_robot366(&robot366,0);
      FiducialProxy neighbor_finderProxy_robot366(&robot366,0);
      PlayerClient    robot367("localhost", 7367);
      Position2dProxy p2dProxy_robot367(&robot367,0);
      RangerProxy      laserProxy_robot367(&robot367,0);
      FiducialProxy neighbor_finderProxy_robot367(&robot367,0);
      PlayerClient    robot368("localhost", 7368);
      Position2dProxy p2dProxy_robot368(&robot368,0);
      RangerProxy      laserProxy_robot368(&robot368,0);
      FiducialProxy neighbor_finderProxy_robot368(&robot368,0);
      PlayerClient    robot369("localhost", 7369);
      Position2dProxy p2dProxy_robot369(&robot369,0);
      RangerProxy      laserProxy_robot369(&robot369,0);
      FiducialProxy neighbor_finderProxy_robot369(&robot369,0);
      PlayerClient    robot370("localhost", 7370);
      Position2dProxy p2dProxy_robot370(&robot370,0);
      RangerProxy      laserProxy_robot370(&robot370,0);
      FiducialProxy neighbor_finderProxy_robot370(&robot370,0);
      PlayerClient    robot371("localhost", 7371);
      Position2dProxy p2dProxy_robot371(&robot371,0);
      RangerProxy      laserProxy_robot371(&robot371,0);
      FiducialProxy neighbor_finderProxy_robot371(&robot371,0);
      PlayerClient    robot372("localhost", 7372);
      Position2dProxy p2dProxy_robot372(&robot372,0);
      RangerProxy      laserProxy_robot372(&robot372,0);
      FiducialProxy neighbor_finderProxy_robot372(&robot372,0);
      PlayerClient    robot373("localhost", 7373);
      Position2dProxy p2dProxy_robot373(&robot373,0);
      RangerProxy      laserProxy_robot373(&robot373,0);
      FiducialProxy neighbor_finderProxy_robot373(&robot373,0);
      PlayerClient    robot374("localhost", 7374);
      Position2dProxy p2dProxy_robot374(&robot374,0);
      RangerProxy      laserProxy_robot374(&robot374,0);
      FiducialProxy neighbor_finderProxy_robot374(&robot374,0);
      PlayerClient    robot375("localhost", 7375);
      Position2dProxy p2dProxy_robot375(&robot375,0);
      RangerProxy      laserProxy_robot375(&robot375,0);
      FiducialProxy neighbor_finderProxy_robot375(&robot375,0);
      PlayerClient    robot376("localhost", 7376);
      Position2dProxy p2dProxy_robot376(&robot376,0);
      RangerProxy      laserProxy_robot376(&robot376,0);
      FiducialProxy neighbor_finderProxy_robot376(&robot376,0);
      PlayerClient    robot377("localhost", 7377);
      Position2dProxy p2dProxy_robot377(&robot377,0);
      RangerProxy      laserProxy_robot377(&robot377,0);
      FiducialProxy neighbor_finderProxy_robot377(&robot377,0);
      PlayerClient    robot378("localhost", 7378);
      Position2dProxy p2dProxy_robot378(&robot378,0);
      RangerProxy      laserProxy_robot378(&robot378,0);
      FiducialProxy neighbor_finderProxy_robot378(&robot378,0);
      PlayerClient    robot379("localhost", 7379);
      Position2dProxy p2dProxy_robot379(&robot379,0);
      RangerProxy      laserProxy_robot379(&robot379,0);
      FiducialProxy neighbor_finderProxy_robot379(&robot379,0);
      PlayerClient    robot380("localhost", 7380);
      Position2dProxy p2dProxy_robot380(&robot380,0);
      RangerProxy      laserProxy_robot380(&robot380,0);
      FiducialProxy neighbor_finderProxy_robot380(&robot380,0);
      PlayerClient    robot381("localhost", 7381);
      Position2dProxy p2dProxy_robot381(&robot381,0);
      RangerProxy      laserProxy_robot381(&robot381,0);
      FiducialProxy neighbor_finderProxy_robot381(&robot381,0);
      PlayerClient    robot382("localhost", 7382);
      Position2dProxy p2dProxy_robot382(&robot382,0);
      RangerProxy      laserProxy_robot382(&robot382,0);
      FiducialProxy neighbor_finderProxy_robot382(&robot382,0);
      PlayerClient    robot383("localhost", 7383);
      Position2dProxy p2dProxy_robot383(&robot383,0);
      RangerProxy      laserProxy_robot383(&robot383,0);
      FiducialProxy neighbor_finderProxy_robot383(&robot383,0);
      PlayerClient    robot384("localhost", 7384);
      Position2dProxy p2dProxy_robot384(&robot384,0);
      RangerProxy      laserProxy_robot384(&robot384,0);
      FiducialProxy neighbor_finderProxy_robot384(&robot384,0);
      PlayerClient    robot385("localhost", 7385);
      Position2dProxy p2dProxy_robot385(&robot385,0);
      RangerProxy      laserProxy_robot385(&robot385,0);
      FiducialProxy neighbor_finderProxy_robot385(&robot385,0);
      PlayerClient    robot386("localhost", 7386);
      Position2dProxy p2dProxy_robot386(&robot386,0);
      RangerProxy      laserProxy_robot386(&robot386,0);
      FiducialProxy neighbor_finderProxy_robot386(&robot386,0);
      PlayerClient    robot387("localhost", 7387);
      Position2dProxy p2dProxy_robot387(&robot387,0);
      RangerProxy      laserProxy_robot387(&robot387,0);
      FiducialProxy neighbor_finderProxy_robot387(&robot387,0);
      PlayerClient    robot388("localhost", 7388);
      Position2dProxy p2dProxy_robot388(&robot388,0);
      RangerProxy      laserProxy_robot388(&robot388,0);
      FiducialProxy neighbor_finderProxy_robot388(&robot388,0);
      PlayerClient    robot389("localhost", 7389);
      Position2dProxy p2dProxy_robot389(&robot389,0);
      RangerProxy      laserProxy_robot389(&robot389,0);
      FiducialProxy neighbor_finderProxy_robot389(&robot389,0);
      PlayerClient    robot390("localhost", 7390);
      Position2dProxy p2dProxy_robot390(&robot390,0);
      RangerProxy      laserProxy_robot390(&robot390,0);
      FiducialProxy neighbor_finderProxy_robot390(&robot390,0);
      PlayerClient    robot391("localhost", 7391);
      Position2dProxy p2dProxy_robot391(&robot391,0);
      RangerProxy      laserProxy_robot391(&robot391,0);
      FiducialProxy neighbor_finderProxy_robot391(&robot391,0);
      PlayerClient    robot392("localhost", 7392);
      Position2dProxy p2dProxy_robot392(&robot392,0);
      RangerProxy      laserProxy_robot392(&robot392,0);
      FiducialProxy neighbor_finderProxy_robot392(&robot392,0);
      PlayerClient    robot393("localhost", 7393);
      Position2dProxy p2dProxy_robot393(&robot393,0);
      RangerProxy      laserProxy_robot393(&robot393,0);
      FiducialProxy neighbor_finderProxy_robot393(&robot393,0);
      PlayerClient    robot394("localhost", 7394);
      Position2dProxy p2dProxy_robot394(&robot394,0);
      RangerProxy      laserProxy_robot394(&robot394,0);
      FiducialProxy neighbor_finderProxy_robot394(&robot394,0);
      PlayerClient    robot395("localhost", 7395);
      Position2dProxy p2dProxy_robot395(&robot395,0);
      RangerProxy      laserProxy_robot395(&robot395,0);
      FiducialProxy neighbor_finderProxy_robot395(&robot395,0);
      PlayerClient    robot396("localhost", 7396);
      Position2dProxy p2dProxy_robot396(&robot396,0);
      RangerProxy      laserProxy_robot396(&robot396,0);
      FiducialProxy neighbor_finderProxy_robot396(&robot396,0);
      PlayerClient    robot397("localhost", 7397);
      Position2dProxy p2dProxy_robot397(&robot397,0);
      RangerProxy      laserProxy_robot397(&robot397,0);
      FiducialProxy neighbor_finderProxy_robot397(&robot397,0);
      PlayerClient    robot398("localhost", 7398);
      Position2dProxy p2dProxy_robot398(&robot398,0);
      RangerProxy      laserProxy_robot398(&robot398,0);
      FiducialProxy neighbor_finderProxy_robot398(&robot398,0);
      PlayerClient    robot399("localhost", 7399);
      Position2dProxy p2dProxy_robot399(&robot399,0);
      RangerProxy      laserProxy_robot399(&robot399,0);
      FiducialProxy neighbor_finderProxy_robot399(&robot399,0);
      PlayerClient    robot400("localhost", 7400);
      Position2dProxy p2dProxy_robot400(&robot400,0);
      RangerProxy      laserProxy_robot400(&robot400,0);
      FiducialProxy neighbor_finderProxy_robot400(&robot400,0);
      PlayerClient    robot401("localhost", 7401);
      Position2dProxy p2dProxy_robot401(&robot401,0);
      RangerProxy      laserProxy_robot401(&robot401,0);
      FiducialProxy neighbor_finderProxy_robot401(&robot401,0);
      PlayerClient    robot402("localhost", 7402);
      Position2dProxy p2dProxy_robot402(&robot402,0);
      RangerProxy      laserProxy_robot402(&robot402,0);
      FiducialProxy neighbor_finderProxy_robot402(&robot402,0);
      PlayerClient    robot403("localhost", 7403);
      Position2dProxy p2dProxy_robot403(&robot403,0);
      RangerProxy      laserProxy_robot403(&robot403,0);
      FiducialProxy neighbor_finderProxy_robot403(&robot403,0);
      PlayerClient    robot404("localhost", 7404);
      Position2dProxy p2dProxy_robot404(&robot404,0);
      RangerProxy      laserProxy_robot404(&robot404,0);
      FiducialProxy neighbor_finderProxy_robot404(&robot404,0);
      PlayerClient    robot405("localhost", 7405);
      Position2dProxy p2dProxy_robot405(&robot405,0);
      RangerProxy      laserProxy_robot405(&robot405,0);
      FiducialProxy neighbor_finderProxy_robot405(&robot405,0);
      PlayerClient    robot406("localhost", 7406);
      Position2dProxy p2dProxy_robot406(&robot406,0);
      RangerProxy      laserProxy_robot406(&robot406,0);
      FiducialProxy neighbor_finderProxy_robot406(&robot406,0);
      PlayerClient    robot407("localhost", 7407);
      Position2dProxy p2dProxy_robot407(&robot407,0);
      RangerProxy      laserProxy_robot407(&robot407,0);
      FiducialProxy neighbor_finderProxy_robot407(&robot407,0);
      PlayerClient    robot408("localhost", 7408);
      Position2dProxy p2dProxy_robot408(&robot408,0);
      RangerProxy      laserProxy_robot408(&robot408,0);
      FiducialProxy neighbor_finderProxy_robot408(&robot408,0);
      PlayerClient    robot409("localhost", 7409);
      Position2dProxy p2dProxy_robot409(&robot409,0);
      RangerProxy      laserProxy_robot409(&robot409,0);
      FiducialProxy neighbor_finderProxy_robot409(&robot409,0);
      PlayerClient    robot410("localhost", 7410);
      Position2dProxy p2dProxy_robot410(&robot410,0);
      RangerProxy      laserProxy_robot410(&robot410,0);
      FiducialProxy neighbor_finderProxy_robot410(&robot410,0);
      PlayerClient    robot411("localhost", 7411);
      Position2dProxy p2dProxy_robot411(&robot411,0);
      RangerProxy      laserProxy_robot411(&robot411,0);
      FiducialProxy neighbor_finderProxy_robot411(&robot411,0);
      PlayerClient    robot412("localhost", 7412);
      Position2dProxy p2dProxy_robot412(&robot412,0);
      RangerProxy      laserProxy_robot412(&robot412,0);
      FiducialProxy neighbor_finderProxy_robot412(&robot412,0);
      PlayerClient    robot413("localhost", 7413);
      Position2dProxy p2dProxy_robot413(&robot413,0);
      RangerProxy      laserProxy_robot413(&robot413,0);
      FiducialProxy neighbor_finderProxy_robot413(&robot413,0);
      PlayerClient    robot414("localhost", 7414);
      Position2dProxy p2dProxy_robot414(&robot414,0);
      RangerProxy      laserProxy_robot414(&robot414,0);
      FiducialProxy neighbor_finderProxy_robot414(&robot414,0);
      PlayerClient    robot415("localhost", 7415);
      Position2dProxy p2dProxy_robot415(&robot415,0);
      RangerProxy      laserProxy_robot415(&robot415,0);
      FiducialProxy neighbor_finderProxy_robot415(&robot415,0);
      PlayerClient    robot416("localhost", 7416);
      Position2dProxy p2dProxy_robot416(&robot416,0);
      RangerProxy      laserProxy_robot416(&robot416,0);
      FiducialProxy neighbor_finderProxy_robot416(&robot416,0);
      PlayerClient    robot417("localhost", 7417);
      Position2dProxy p2dProxy_robot417(&robot417,0);
      RangerProxy      laserProxy_robot417(&robot417,0);
      FiducialProxy neighbor_finderProxy_robot417(&robot417,0);
      PlayerClient    robot418("localhost", 7418);
      Position2dProxy p2dProxy_robot418(&robot418,0);
      RangerProxy      laserProxy_robot418(&robot418,0);
      FiducialProxy neighbor_finderProxy_robot418(&robot418,0);
      PlayerClient    robot419("localhost", 7419);
      Position2dProxy p2dProxy_robot419(&robot419,0);
      RangerProxy      laserProxy_robot419(&robot419,0);
      FiducialProxy neighbor_finderProxy_robot419(&robot419,0);

      // vector<PlayerClient> plyclnts(5);
      // vector<Position2dProxy> p2dProxys(5);
      // vector<RangerProxy> rngProxys(5);
      // vector<FiducialProxy> fidProxys(5);
      // vector<SimulationProxy> SimProxys(5);
// spawn off models
      PlayerClient plyclnts[SWARM_SIZE] ={robot360, robot361, robot362, robot363, robot364, robot365, robot366, robot367, robot368, robot369, robot370, robot371, robot372, robot373, robot374, robot375, robot376, robot377, robot378, robot379, robot380, robot381, robot382, robot383, robot384, robot385, robot386, robot387, robot388, robot389, robot390, robot391, robot392, robot393, robot394, robot395, robot396, robot397, robot398, robot399, robot400, robot401, robot402, robot403, robot404, robot405, robot406, robot407, robot408, robot409, robot410, robot411, robot412, robot413, robot414, robot415, robot416, robot417, robot418, robot419};
      Position2dProxy p2dProxys[SWARM_SIZE] = {p2dProxy_robot360, p2dProxy_robot361, p2dProxy_robot362, p2dProxy_robot363, p2dProxy_robot364, p2dProxy_robot365, p2dProxy_robot366, p2dProxy_robot367, p2dProxy_robot368, p2dProxy_robot369, p2dProxy_robot370, p2dProxy_robot371, p2dProxy_robot372, p2dProxy_robot373, p2dProxy_robot374, p2dProxy_robot375, p2dProxy_robot376, p2dProxy_robot377, p2dProxy_robot378, p2dProxy_robot379, p2dProxy_robot380, p2dProxy_robot381, p2dProxy_robot382, p2dProxy_robot383, p2dProxy_robot384, p2dProxy_robot385, p2dProxy_robot386, p2dProxy_robot387, p2dProxy_robot388, p2dProxy_robot389, p2dProxy_robot390, p2dProxy_robot391, p2dProxy_robot392, p2dProxy_robot393, p2dProxy_robot394, p2dProxy_robot395, p2dProxy_robot396, p2dProxy_robot397, p2dProxy_robot398, p2dProxy_robot399, p2dProxy_robot400, p2dProxy_robot401, p2dProxy_robot402, p2dProxy_robot403, p2dProxy_robot404, p2dProxy_robot405, p2dProxy_robot406, p2dProxy_robot407, p2dProxy_robot408, p2dProxy_robot409, p2dProxy_robot410, p2dProxy_robot411, p2dProxy_robot412, p2dProxy_robot413, p2dProxy_robot414, p2dProxy_robot415, p2dProxy_robot416, p2dProxy_robot417, p2dProxy_robot418, p2dProxy_robot419};
      RangerProxy rngProxys[SWARM_SIZE] = {laserProxy_robot360, laserProxy_robot361, laserProxy_robot362, laserProxy_robot363, laserProxy_robot364, laserProxy_robot365, laserProxy_robot366, laserProxy_robot367, laserProxy_robot368, laserProxy_robot369, laserProxy_robot370, laserProxy_robot371, laserProxy_robot372, laserProxy_robot373, laserProxy_robot374, laserProxy_robot375, laserProxy_robot376, laserProxy_robot377, laserProxy_robot378, laserProxy_robot379, laserProxy_robot380, laserProxy_robot381, laserProxy_robot382, laserProxy_robot383, laserProxy_robot384, laserProxy_robot385, laserProxy_robot386, laserProxy_robot387, laserProxy_robot388, laserProxy_robot389, laserProxy_robot390, laserProxy_robot391, laserProxy_robot392, laserProxy_robot393, laserProxy_robot394, laserProxy_robot395, laserProxy_robot396, laserProxy_robot397, laserProxy_robot398, laserProxy_robot399, laserProxy_robot400, laserProxy_robot401, laserProxy_robot402, laserProxy_robot403, laserProxy_robot404, laserProxy_robot405, laserProxy_robot406, laserProxy_robot407, laserProxy_robot408, laserProxy_robot409, laserProxy_robot410, laserProxy_robot411, laserProxy_robot412, laserProxy_robot413, laserProxy_robot414, laserProxy_robot415, laserProxy_robot416, laserProxy_robot417, laserProxy_robot418, laserProxy_robot419};
      FiducialProxy fidProxys[SWARM_SIZE]={neighbor_finderProxy_robot360, neighbor_finderProxy_robot361, neighbor_finderProxy_robot362, neighbor_finderProxy_robot363, neighbor_finderProxy_robot364, neighbor_finderProxy_robot365, neighbor_finderProxy_robot366, neighbor_finderProxy_robot367, neighbor_finderProxy_robot368, neighbor_finderProxy_robot369, neighbor_finderProxy_robot370, neighbor_finderProxy_robot371, neighbor_finderProxy_robot372, neighbor_finderProxy_robot373, neighbor_finderProxy_robot374, neighbor_finderProxy_robot375, neighbor_finderProxy_robot376, neighbor_finderProxy_robot377, neighbor_finderProxy_robot378, neighbor_finderProxy_robot379, neighbor_finderProxy_robot380, neighbor_finderProxy_robot381, neighbor_finderProxy_robot382, neighbor_finderProxy_robot383, neighbor_finderProxy_robot384, neighbor_finderProxy_robot385, neighbor_finderProxy_robot386, neighbor_finderProxy_robot387, neighbor_finderProxy_robot388, neighbor_finderProxy_robot389, neighbor_finderProxy_robot390, neighbor_finderProxy_robot391, neighbor_finderProxy_robot392, neighbor_finderProxy_robot393, neighbor_finderProxy_robot394, neighbor_finderProxy_robot395, neighbor_finderProxy_robot396, neighbor_finderProxy_robot397, neighbor_finderProxy_robot398, neighbor_finderProxy_robot399, neighbor_finderProxy_robot400, neighbor_finderProxy_robot401, neighbor_finderProxy_robot402, neighbor_finderProxy_robot403, neighbor_finderProxy_robot404, neighbor_finderProxy_robot405, neighbor_finderProxy_robot406, neighbor_finderProxy_robot407, neighbor_finderProxy_robot408, neighbor_finderProxy_robot409, neighbor_finderProxy_robot410, neighbor_finderProxy_robot411, neighbor_finderProxy_robot412, neighbor_finderProxy_robot413, neighbor_finderProxy_robot414, neighbor_finderProxy_robot415, neighbor_finderProxy_robot416, neighbor_finderProxy_robot417, neighbor_finderProxy_robot418, neighbor_finderProxy_robot419};
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
