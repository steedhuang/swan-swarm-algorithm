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
            if (iD >= SWARM_SIZE/3)
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
      PlayerClient    robot0("localhost", 7000);
      Position2dProxy p2dProxy_robot0(&robot0,0);
      RangerProxy      laserProxy_robot0(&robot0,0);
      FiducialProxy neighbor_finderProxy_robot0(&robot0,0);
      // SimulationProxy simProxy_robot0(&robot0,0);

      PlayerClient    robot1("localhost", 7001);
      Position2dProxy p2dProxy_robot1(&robot1,0);
      RangerProxy      laserProxy_robot1(&robot1,0);
      FiducialProxy neighbor_finderProxy_robot1(&robot1,0);
      // SimulationProxy simProxy_robot1(&robot1,0);

      PlayerClient    robot2("localhost", 7002);
      Position2dProxy p2dProxy_robot2(&robot2,0);
      RangerProxy      laserProxy_robot2(&robot2,0);
      FiducialProxy neighbor_finderProxy_robot2(&robot2,0);
      // SimulationProxy simProxy_robot2(&robot2,0);

      PlayerClient    robot3("localhost", 7003);
      Position2dProxy p2dProxy_robot3(&robot3,0);
      RangerProxy      laserProxy_robot3(&robot3,0);
      FiducialProxy neighbor_finderProxy_robot3(&robot3,0);
      // SimulationProxy simProxy_robot3(&robot3,0);

      PlayerClient    robot4("localhost", 7004);
      Position2dProxy p2dProxy_robot4(&robot4,0);
      RangerProxy      laserProxy_robot4(&robot4,0);
      FiducialProxy neighbor_finderProxy_robot4(&robot4,0);
      // SimulationProxy simProxy_robot4(&robot4,0);

      PlayerClient    robot5("localhost", 7005);
      Position2dProxy p2dProxy_robot5(&robot5,0);
      RangerProxy      laserProxy_robot5(&robot5,0);
      FiducialProxy neighbor_finderProxy_robot5(&robot5,0);
      // SimulationProxy simProxy_robot5(&robot5,0);

      PlayerClient    robot6("localhost", 7006);
      Position2dProxy p2dProxy_robot6(&robot6,0);
      RangerProxy      laserProxy_robot6(&robot6,0);
      FiducialProxy neighbor_finderProxy_robot6(&robot6,0);
      // SimulationProxy simProxy_robot6(&robot6,0);

      PlayerClient    robot7("localhost", 7007);
      Position2dProxy p2dProxy_robot7(&robot7,0);
      RangerProxy      laserProxy_robot7(&robot7,0);
      FiducialProxy neighbor_finderProxy_robot7(&robot7,0);
      // SimulationProxy simProxy_robot7(&robot7,0);

      PlayerClient    robot8("localhost", 7008);
      Position2dProxy p2dProxy_robot8(&robot8,0);
      RangerProxy      laserProxy_robot8(&robot8,0);
      FiducialProxy neighbor_finderProxy_robot8(&robot8,0);
      // SimulationProxy simProxy_robot8(&robot8,0);

      PlayerClient    robot9("localhost", 7009);
      Position2dProxy p2dProxy_robot9(&robot9,0);
      RangerProxy      laserProxy_robot9(&robot9,0);
      FiducialProxy neighbor_finderProxy_robot9(&robot9,0);
      // SimulationProxy simProxy_robot9(&robot9,0);

      PlayerClient    robot10("localhost", 7010);
      Position2dProxy p2dProxy_robot10(&robot10,0);
      RangerProxy      laserProxy_robot10(&robot10,0);
      FiducialProxy neighbor_finderProxy_robot10(&robot10,0);
      // SimulationProxy simProxy_robot10(&robot10,0);

      PlayerClient    robot11("localhost", 7011);
      Position2dProxy p2dProxy_robot11(&robot11,0);
      RangerProxy      laserProxy_robot11(&robot11,0);
      FiducialProxy neighbor_finderProxy_robot11(&robot11,0);
      // SimulationProxy simProxy_robot11(&robot11,0);

      PlayerClient    robot12("localhost", 7012);
      Position2dProxy p2dProxy_robot12(&robot12,0);
      RangerProxy      laserProxy_robot12(&robot12,0);
      FiducialProxy neighbor_finderProxy_robot12(&robot12,0);
      // SimulationProxy simProxy_robot2(&robot2,0);

      PlayerClient    robot13("localhost", 7013);
      Position2dProxy p2dProxy_robot13(&robot13,0);
      RangerProxy      laserProxy_robot13(&robot13,0);
      FiducialProxy neighbor_finderProxy_robot13(&robot13,0);
      // SimulationProxy simProxy_robot3(&robot3,0);

      PlayerClient    robot14("localhost", 7014);
      Position2dProxy p2dProxy_robot14(&robot14,0);
      RangerProxy      laserProxy_robot14(&robot14,0);
      FiducialProxy neighbor_finderProxy_robot14(&robot14,0);
      // SimulationProxy simProxy_robot4(&robot4,0);

      PlayerClient    robot15("localhost", 7015);
      Position2dProxy p2dProxy_robot15(&robot15,0);
      RangerProxy      laserProxy_robot15(&robot15,0);
      FiducialProxy neighbor_finderProxy_robot15(&robot15,0);
      // SimulationProxy simProxy_robot5(&robot5,0);

      PlayerClient    robot16("localhost", 7016);
      Position2dProxy p2dProxy_robot16(&robot16,0);
      RangerProxy      laserProxy_robot16(&robot16,0);
      FiducialProxy neighbor_finderProxy_robot16(&robot16,0);
      // SimulationProxy simProxy_robot6(&robot6,0);

      PlayerClient    robot17("localhost", 7017);
      Position2dProxy p2dProxy_robot17(&robot17,0);
      RangerProxy      laserProxy_robot17(&robot17,0);
      FiducialProxy neighbor_finderProxy_robot17(&robot17,0);
      // SimulationProxy simProxy_robot7(&robot7,0);

      PlayerClient    robot18("localhost", 7018);
      Position2dProxy p2dProxy_robot18(&robot18,0);
      RangerProxy      laserProxy_robot18(&robot18,0);
      FiducialProxy neighbor_finderProxy_robot18(&robot18,0);
      // SimulationProxy simProxy_robot8(&robot8,0);

      PlayerClient    robot19("localhost", 7019);
      Position2dProxy p2dProxy_robot19(&robot19,0);
      RangerProxy      laserProxy_robot19(&robot19,0);
      FiducialProxy neighbor_finderProxy_robot19(&robot19,0);
      // SimulationProxy simProxy_robot9(&robot9,0);

      PlayerClient    robot20("localhost", 7020);
      Position2dProxy p2dProxy_robot20(&robot20,0);
      RangerProxy      laserProxy_robot20(&robot20,0);
      FiducialProxy neighbor_finderProxy_robot20(&robot20,0);
      // SimulationProxy simProxy_robot10(&robot10,0);

      PlayerClient    robot21("localhost", 7021);
      Position2dProxy p2dProxy_robot21(&robot21,0);
      RangerProxy      laserProxy_robot21(&robot21,0);
      FiducialProxy neighbor_finderProxy_robot21(&robot21,0);
      // SimulationProxy simProxy_robot11(&robot11,0);

      PlayerClient    robot22("localhost", 7022);
      Position2dProxy p2dProxy_robot22(&robot22,0);
      RangerProxy      laserProxy_robot22(&robot22,0);
      FiducialProxy neighbor_finderProxy_robot22(&robot22,0);
      // SimulationProxy simProxy_robot2(&robot2,0);

      PlayerClient    robot23("localhost", 7023);
      Position2dProxy p2dProxy_robot23(&robot23,0);
      RangerProxy      laserProxy_robot23(&robot23,0);
      FiducialProxy neighbor_finderProxy_robot23(&robot23,0);
      // SimulationProxy simProxy_robot3(&robot3,0);

      PlayerClient    robot24("localhost", 7024);
      Position2dProxy p2dProxy_robot24(&robot24,0);
      RangerProxy      laserProxy_robot24(&robot24,0);
      FiducialProxy neighbor_finderProxy_robot24(&robot24,0);
      // SimulationProxy simProxy_robot4(&robot4,0);

      PlayerClient    robot25("localhost", 7025);
      Position2dProxy p2dProxy_robot25(&robot25,0);
      RangerProxy      laserProxy_robot25(&robot25,0);
      FiducialProxy neighbor_finderProxy_robot25(&robot25,0);
      // SimulationProxy simProxy_robot5(&robot5,0);

      PlayerClient    robot26("localhost", 7026);
      Position2dProxy p2dProxy_robot26(&robot26,0);
      RangerProxy      laserProxy_robot26(&robot26,0);
      FiducialProxy neighbor_finderProxy_robot26(&robot26,0);
      // SimulationProxy simProxy_robot6(&robot6,0);

      PlayerClient    robot27("localhost", 7027);
      Position2dProxy p2dProxy_robot27(&robot27,0);
      RangerProxy      laserProxy_robot27(&robot27,0);
      FiducialProxy neighbor_finderProxy_robot27(&robot27,0);
      // SimulationProxy simProxy_robot7(&robot7,0);

      PlayerClient    robot28("localhost", 7028);
      Position2dProxy p2dProxy_robot28(&robot28,0);
      RangerProxy      laserProxy_robot28(&robot28,0);
      FiducialProxy neighbor_finderProxy_robot28(&robot28,0);
      // SimulationProxy simProxy_robot8(&robot8,0);

      PlayerClient    robot29("localhost", 7029);
      Position2dProxy p2dProxy_robot29(&robot29,0);
      RangerProxy      laserProxy_robot29(&robot29,0);
      FiducialProxy neighbor_finderProxy_robot29(&robot29,0);
      // SimulationProxy simProxy_robot9(&robot9,0);

      PlayerClient    robot30("localhost", 7030);
      Position2dProxy p2dProxy_robot30(&robot30,0);
      RangerProxy      laserProxy_robot30(&robot30,0);
      FiducialProxy neighbor_finderProxy_robot30(&robot30,0);
      PlayerClient    robot31("localhost", 7031);
      Position2dProxy p2dProxy_robot31(&robot31,0);
      RangerProxy      laserProxy_robot31(&robot31,0);
      FiducialProxy neighbor_finderProxy_robot31(&robot31,0);
      PlayerClient    robot32("localhost", 7032);
      Position2dProxy p2dProxy_robot32(&robot32,0);
      RangerProxy      laserProxy_robot32(&robot32,0);
      FiducialProxy neighbor_finderProxy_robot32(&robot32,0);
      PlayerClient    robot33("localhost", 7033);
      Position2dProxy p2dProxy_robot33(&robot33,0);
      RangerProxy      laserProxy_robot33(&robot33,0);
      FiducialProxy neighbor_finderProxy_robot33(&robot33,0);
      PlayerClient    robot34("localhost", 7034);
      Position2dProxy p2dProxy_robot34(&robot34,0);
      RangerProxy      laserProxy_robot34(&robot34,0);
      FiducialProxy neighbor_finderProxy_robot34(&robot34,0);
      PlayerClient    robot35("localhost", 7035);
      Position2dProxy p2dProxy_robot35(&robot35,0);
      RangerProxy      laserProxy_robot35(&robot35,0);
      FiducialProxy neighbor_finderProxy_robot35(&robot35,0);
      PlayerClient    robot36("localhost", 7036);
      Position2dProxy p2dProxy_robot36(&robot36,0);
      RangerProxy      laserProxy_robot36(&robot36,0);
      FiducialProxy neighbor_finderProxy_robot36(&robot36,0);
      PlayerClient    robot37("localhost", 7037);
      Position2dProxy p2dProxy_robot37(&robot37,0);
      RangerProxy      laserProxy_robot37(&robot37,0);
      FiducialProxy neighbor_finderProxy_robot37(&robot37,0);
      PlayerClient    robot38("localhost", 7038);
      Position2dProxy p2dProxy_robot38(&robot38,0);
      RangerProxy      laserProxy_robot38(&robot38,0);
      FiducialProxy neighbor_finderProxy_robot38(&robot38,0);
      PlayerClient    robot39("localhost", 7039);
      Position2dProxy p2dProxy_robot39(&robot39,0);
      RangerProxy      laserProxy_robot39(&robot39,0);
      FiducialProxy neighbor_finderProxy_robot39(&robot39,0);
      PlayerClient    robot40("localhost", 7040);
      Position2dProxy p2dProxy_robot40(&robot40,0);
      RangerProxy      laserProxy_robot40(&robot40,0);
      FiducialProxy neighbor_finderProxy_robot40(&robot40,0);
      PlayerClient    robot41("localhost", 7041);
      Position2dProxy p2dProxy_robot41(&robot41,0);
      RangerProxy      laserProxy_robot41(&robot41,0);
      FiducialProxy neighbor_finderProxy_robot41(&robot41,0);
      PlayerClient    robot42("localhost", 7042);
      Position2dProxy p2dProxy_robot42(&robot42,0);
      RangerProxy      laserProxy_robot42(&robot42,0);
      FiducialProxy neighbor_finderProxy_robot42(&robot42,0);
      PlayerClient    robot43("localhost", 7043);
      Position2dProxy p2dProxy_robot43(&robot43,0);
      RangerProxy      laserProxy_robot43(&robot43,0);
      FiducialProxy neighbor_finderProxy_robot43(&robot43,0);
      PlayerClient    robot44("localhost", 7044);
      Position2dProxy p2dProxy_robot44(&robot44,0);
      RangerProxy      laserProxy_robot44(&robot44,0);
      FiducialProxy neighbor_finderProxy_robot44(&robot44,0);
      PlayerClient    robot45("localhost", 7045);
      Position2dProxy p2dProxy_robot45(&robot45,0);
      RangerProxy      laserProxy_robot45(&robot45,0);
      FiducialProxy neighbor_finderProxy_robot45(&robot45,0);
      PlayerClient    robot46("localhost", 7046);
      Position2dProxy p2dProxy_robot46(&robot46,0);
      RangerProxy      laserProxy_robot46(&robot46,0);
      FiducialProxy neighbor_finderProxy_robot46(&robot46,0);
      PlayerClient    robot47("localhost", 7047);
      Position2dProxy p2dProxy_robot47(&robot47,0);
      RangerProxy      laserProxy_robot47(&robot47,0);
      FiducialProxy neighbor_finderProxy_robot47(&robot47,0);
      PlayerClient    robot48("localhost", 7048);
      Position2dProxy p2dProxy_robot48(&robot48,0);
      RangerProxy      laserProxy_robot48(&robot48,0);
      FiducialProxy neighbor_finderProxy_robot48(&robot48,0);
      PlayerClient    robot49("localhost", 7049);
      Position2dProxy p2dProxy_robot49(&robot49,0);
      RangerProxy      laserProxy_robot49(&robot49,0);
      FiducialProxy neighbor_finderProxy_robot49(&robot49,0);
      PlayerClient    robot50("localhost", 7050);
      Position2dProxy p2dProxy_robot50(&robot50,0);
      RangerProxy      laserProxy_robot50(&robot50,0);
      FiducialProxy neighbor_finderProxy_robot50(&robot50,0);
      PlayerClient    robot51("localhost", 7051);
      Position2dProxy p2dProxy_robot51(&robot51,0);
      RangerProxy      laserProxy_robot51(&robot51,0);
      FiducialProxy neighbor_finderProxy_robot51(&robot51,0);
      PlayerClient    robot52("localhost", 7052);
      Position2dProxy p2dProxy_robot52(&robot52,0);
      RangerProxy      laserProxy_robot52(&robot52,0);
      FiducialProxy neighbor_finderProxy_robot52(&robot52,0);
      PlayerClient    robot53("localhost", 7053);
      Position2dProxy p2dProxy_robot53(&robot53,0);
      RangerProxy      laserProxy_robot53(&robot53,0);
      FiducialProxy neighbor_finderProxy_robot53(&robot53,0);
      PlayerClient    robot54("localhost", 7054);
      Position2dProxy p2dProxy_robot54(&robot54,0);
      RangerProxy      laserProxy_robot54(&robot54,0);
      FiducialProxy neighbor_finderProxy_robot54(&robot54,0);
      PlayerClient    robot55("localhost", 7055);
      Position2dProxy p2dProxy_robot55(&robot55,0);
      RangerProxy      laserProxy_robot55(&robot55,0);
      FiducialProxy neighbor_finderProxy_robot55(&robot55,0);
      PlayerClient    robot56("localhost", 7056);
      Position2dProxy p2dProxy_robot56(&robot56,0);
      RangerProxy      laserProxy_robot56(&robot56,0);
      FiducialProxy neighbor_finderProxy_robot56(&robot56,0);
      PlayerClient    robot57("localhost", 7057);
      Position2dProxy p2dProxy_robot57(&robot57,0);
      RangerProxy      laserProxy_robot57(&robot57,0);
      FiducialProxy neighbor_finderProxy_robot57(&robot57,0);
      PlayerClient    robot58("localhost", 7058);
      Position2dProxy p2dProxy_robot58(&robot58,0);
      RangerProxy      laserProxy_robot58(&robot58,0);
      FiducialProxy neighbor_finderProxy_robot58(&robot58,0);
      PlayerClient    robot59("localhost", 7059);
      Position2dProxy p2dProxy_robot59(&robot59,0);
      RangerProxy      laserProxy_robot59(&robot59,0);
      FiducialProxy neighbor_finderProxy_robot59(&robot59,0);

      // vector<PlayerClient> plyclnts(5);
      // vector<Position2dProxy> p2dProxys(5);
      // vector<RangerProxy> rngProxys(5);
      // vector<FiducialProxy> fidProxys(5);
      // vector<SimulationProxy> SimProxys(5);
// spawn off models
      PlayerClient plyclnts[SWARM_SIZE] ={robot0, robot1, robot2, robot3, robot4, robot5, robot6, robot7, robot8, robot9, robot10, robot11, robot12, robot13, robot14, robot15, robot16, robot17, robot18, robot19, robot20, robot21, robot22, robot23, robot24, robot25, robot26, robot27, robot28, robot29, robot30, robot31, robot32, robot33, robot34, robot35, robot36, robot37, robot38, robot39, robot40, robot41, robot42, robot43, robot44, robot45, robot46, robot47, robot48, robot49, robot50, robot51, robot52, robot53, robot54, robot55, robot56, robot57, robot58, robot59};
      Position2dProxy p2dProxys[SWARM_SIZE] = {p2dProxy_robot0, p2dProxy_robot1, p2dProxy_robot2, p2dProxy_robot3, p2dProxy_robot4, p2dProxy_robot5, p2dProxy_robot6, p2dProxy_robot7, p2dProxy_robot8, p2dProxy_robot9, p2dProxy_robot10, p2dProxy_robot11, p2dProxy_robot12, p2dProxy_robot13, p2dProxy_robot14, p2dProxy_robot15, p2dProxy_robot16, p2dProxy_robot17, p2dProxy_robot18, p2dProxy_robot19, p2dProxy_robot20, p2dProxy_robot21, p2dProxy_robot22, p2dProxy_robot23, p2dProxy_robot24, p2dProxy_robot25, p2dProxy_robot26, p2dProxy_robot27, p2dProxy_robot28, p2dProxy_robot29, p2dProxy_robot30, p2dProxy_robot31, p2dProxy_robot32, p2dProxy_robot33, p2dProxy_robot34, p2dProxy_robot35, p2dProxy_robot36, p2dProxy_robot37, p2dProxy_robot38, p2dProxy_robot39, p2dProxy_robot40, p2dProxy_robot41, p2dProxy_robot42, p2dProxy_robot43, p2dProxy_robot44, p2dProxy_robot45, p2dProxy_robot46, p2dProxy_robot47, p2dProxy_robot48, p2dProxy_robot49, p2dProxy_robot50, p2dProxy_robot51, p2dProxy_robot52, p2dProxy_robot53, p2dProxy_robot54, p2dProxy_robot55, p2dProxy_robot56, p2dProxy_robot57, p2dProxy_robot58, p2dProxy_robot59};
      RangerProxy rngProxys[SWARM_SIZE] = {laserProxy_robot0, laserProxy_robot1, laserProxy_robot2, laserProxy_robot3, laserProxy_robot4, laserProxy_robot5, laserProxy_robot6, laserProxy_robot7, laserProxy_robot8, laserProxy_robot9, laserProxy_robot10, laserProxy_robot11, laserProxy_robot12, laserProxy_robot13, laserProxy_robot14, laserProxy_robot15, laserProxy_robot16, laserProxy_robot17, laserProxy_robot18, laserProxy_robot19, laserProxy_robot20, laserProxy_robot21, laserProxy_robot22, laserProxy_robot23, laserProxy_robot24, laserProxy_robot25, laserProxy_robot26, laserProxy_robot27, laserProxy_robot28, laserProxy_robot29, laserProxy_robot30, laserProxy_robot31, laserProxy_robot32, laserProxy_robot33, laserProxy_robot34, laserProxy_robot35, laserProxy_robot36, laserProxy_robot37, laserProxy_robot38, laserProxy_robot39, laserProxy_robot40, laserProxy_robot41, laserProxy_robot42, laserProxy_robot43, laserProxy_robot44, laserProxy_robot45, laserProxy_robot46, laserProxy_robot47, laserProxy_robot48, laserProxy_robot49, laserProxy_robot50, laserProxy_robot51, laserProxy_robot52, laserProxy_robot53, laserProxy_robot54, laserProxy_robot55, laserProxy_robot56, laserProxy_robot57, laserProxy_robot58, laserProxy_robot59};
      FiducialProxy fidProxys[SWARM_SIZE]={neighbor_finderProxy_robot0, neighbor_finderProxy_robot1, neighbor_finderProxy_robot2, neighbor_finderProxy_robot3, neighbor_finderProxy_robot4, neighbor_finderProxy_robot5, neighbor_finderProxy_robot6, neighbor_finderProxy_robot7, neighbor_finderProxy_robot8, neighbor_finderProxy_robot9, neighbor_finderProxy_robot10, neighbor_finderProxy_robot11, neighbor_finderProxy_robot12, neighbor_finderProxy_robot13, neighbor_finderProxy_robot14, neighbor_finderProxy_robot15, neighbor_finderProxy_robot16, neighbor_finderProxy_robot17, neighbor_finderProxy_robot18, neighbor_finderProxy_robot19, neighbor_finderProxy_robot20, neighbor_finderProxy_robot21, neighbor_finderProxy_robot22, neighbor_finderProxy_robot23, neighbor_finderProxy_robot24, neighbor_finderProxy_robot25, neighbor_finderProxy_robot26, neighbor_finderProxy_robot27, neighbor_finderProxy_robot28, neighbor_finderProxy_robot29, neighbor_finderProxy_robot30, neighbor_finderProxy_robot31, neighbor_finderProxy_robot32, neighbor_finderProxy_robot33, neighbor_finderProxy_robot34, neighbor_finderProxy_robot35, neighbor_finderProxy_robot36, neighbor_finderProxy_robot37, neighbor_finderProxy_robot38, neighbor_finderProxy_robot39, neighbor_finderProxy_robot40, neighbor_finderProxy_robot41, neighbor_finderProxy_robot42, neighbor_finderProxy_robot43, neighbor_finderProxy_robot44, neighbor_finderProxy_robot45, neighbor_finderProxy_robot46, neighbor_finderProxy_robot47, neighbor_finderProxy_robot48, neighbor_finderProxy_robot49, neighbor_finderProxy_robot50, neighbor_finderProxy_robot51, neighbor_finderProxy_robot52, neighbor_finderProxy_robot53, neighbor_finderProxy_robot54, neighbor_finderProxy_robot55, neighbor_finderProxy_robot56, neighbor_finderProxy_robot57, neighbor_finderProxy_robot58, neighbor_finderProxy_robot59};
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
