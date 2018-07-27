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
      PlayerClient    robot120("localhost", 7120);
      Position2dProxy p2dProxy_robot120(&robot120,0);
      RangerProxy      laserProxy_robot120(&robot120,0);
      FiducialProxy neighbor_finderProxy_robot120(&robot120,0);
      PlayerClient    robot121("localhost", 7121);
      Position2dProxy p2dProxy_robot121(&robot121,0);
      RangerProxy      laserProxy_robot121(&robot121,0);
      FiducialProxy neighbor_finderProxy_robot121(&robot121,0);
      PlayerClient    robot122("localhost", 7122);
      Position2dProxy p2dProxy_robot122(&robot122,0);
      RangerProxy      laserProxy_robot122(&robot122,0);
      FiducialProxy neighbor_finderProxy_robot122(&robot122,0);
      PlayerClient    robot123("localhost", 7123);
      Position2dProxy p2dProxy_robot123(&robot123,0);
      RangerProxy      laserProxy_robot123(&robot123,0);
      FiducialProxy neighbor_finderProxy_robot123(&robot123,0);
      PlayerClient    robot124("localhost", 7124);
      Position2dProxy p2dProxy_robot124(&robot124,0);
      RangerProxy      laserProxy_robot124(&robot124,0);
      FiducialProxy neighbor_finderProxy_robot124(&robot124,0);
      PlayerClient    robot125("localhost", 7125);
      Position2dProxy p2dProxy_robot125(&robot125,0);
      RangerProxy      laserProxy_robot125(&robot125,0);
      FiducialProxy neighbor_finderProxy_robot125(&robot125,0);
      PlayerClient    robot126("localhost", 7126);
      Position2dProxy p2dProxy_robot126(&robot126,0);
      RangerProxy      laserProxy_robot126(&robot126,0);
      FiducialProxy neighbor_finderProxy_robot126(&robot126,0);
      PlayerClient    robot127("localhost", 7127);
      Position2dProxy p2dProxy_robot127(&robot127,0);
      RangerProxy      laserProxy_robot127(&robot127,0);
      FiducialProxy neighbor_finderProxy_robot127(&robot127,0);
      PlayerClient    robot128("localhost", 7128);
      Position2dProxy p2dProxy_robot128(&robot128,0);
      RangerProxy      laserProxy_robot128(&robot128,0);
      FiducialProxy neighbor_finderProxy_robot128(&robot128,0);
      PlayerClient    robot129("localhost", 7129);
      Position2dProxy p2dProxy_robot129(&robot129,0);
      RangerProxy      laserProxy_robot129(&robot129,0);
      FiducialProxy neighbor_finderProxy_robot129(&robot129,0);
      PlayerClient    robot130("localhost", 7130);
      Position2dProxy p2dProxy_robot130(&robot130,0);
      RangerProxy      laserProxy_robot130(&robot130,0);
      FiducialProxy neighbor_finderProxy_robot130(&robot130,0);
      PlayerClient    robot131("localhost", 7131);
      Position2dProxy p2dProxy_robot131(&robot131,0);
      RangerProxy      laserProxy_robot131(&robot131,0);
      FiducialProxy neighbor_finderProxy_robot131(&robot131,0);
      PlayerClient    robot132("localhost", 7132);
      Position2dProxy p2dProxy_robot132(&robot132,0);
      RangerProxy      laserProxy_robot132(&robot132,0);
      FiducialProxy neighbor_finderProxy_robot132(&robot132,0);
      PlayerClient    robot133("localhost", 7133);
      Position2dProxy p2dProxy_robot133(&robot133,0);
      RangerProxy      laserProxy_robot133(&robot133,0);
      FiducialProxy neighbor_finderProxy_robot133(&robot133,0);
      PlayerClient    robot134("localhost", 7134);
      Position2dProxy p2dProxy_robot134(&robot134,0);
      RangerProxy      laserProxy_robot134(&robot134,0);
      FiducialProxy neighbor_finderProxy_robot134(&robot134,0);
      PlayerClient    robot135("localhost", 7135);
      Position2dProxy p2dProxy_robot135(&robot135,0);
      RangerProxy      laserProxy_robot135(&robot135,0);
      FiducialProxy neighbor_finderProxy_robot135(&robot135,0);
      PlayerClient    robot136("localhost", 7136);
      Position2dProxy p2dProxy_robot136(&robot136,0);
      RangerProxy      laserProxy_robot136(&robot136,0);
      FiducialProxy neighbor_finderProxy_robot136(&robot136,0);
      PlayerClient    robot137("localhost", 7137);
      Position2dProxy p2dProxy_robot137(&robot137,0);
      RangerProxy      laserProxy_robot137(&robot137,0);
      FiducialProxy neighbor_finderProxy_robot137(&robot137,0);
      PlayerClient    robot138("localhost", 7138);
      Position2dProxy p2dProxy_robot138(&robot138,0);
      RangerProxy      laserProxy_robot138(&robot138,0);
      FiducialProxy neighbor_finderProxy_robot138(&robot138,0);
      PlayerClient    robot139("localhost", 7139);
      Position2dProxy p2dProxy_robot139(&robot139,0);
      RangerProxy      laserProxy_robot139(&robot139,0);
      FiducialProxy neighbor_finderProxy_robot139(&robot139,0);
      PlayerClient    robot140("localhost", 7140);
      Position2dProxy p2dProxy_robot140(&robot140,0);
      RangerProxy      laserProxy_robot140(&robot140,0);
      FiducialProxy neighbor_finderProxy_robot140(&robot140,0);
      PlayerClient    robot141("localhost", 7141);
      Position2dProxy p2dProxy_robot141(&robot141,0);
      RangerProxy      laserProxy_robot141(&robot141,0);
      FiducialProxy neighbor_finderProxy_robot141(&robot141,0);
      PlayerClient    robot142("localhost", 7142);
      Position2dProxy p2dProxy_robot142(&robot142,0);
      RangerProxy      laserProxy_robot142(&robot142,0);
      FiducialProxy neighbor_finderProxy_robot142(&robot142,0);
      PlayerClient    robot143("localhost", 7143);
      Position2dProxy p2dProxy_robot143(&robot143,0);
      RangerProxy      laserProxy_robot143(&robot143,0);
      FiducialProxy neighbor_finderProxy_robot143(&robot143,0);
      PlayerClient    robot144("localhost", 7144);
      Position2dProxy p2dProxy_robot144(&robot144,0);
      RangerProxy      laserProxy_robot144(&robot144,0);
      FiducialProxy neighbor_finderProxy_robot144(&robot144,0);
      PlayerClient    robot145("localhost", 7145);
      Position2dProxy p2dProxy_robot145(&robot145,0);
      RangerProxy      laserProxy_robot145(&robot145,0);
      FiducialProxy neighbor_finderProxy_robot145(&robot145,0);
      PlayerClient    robot146("localhost", 7146);
      Position2dProxy p2dProxy_robot146(&robot146,0);
      RangerProxy      laserProxy_robot146(&robot146,0);
      FiducialProxy neighbor_finderProxy_robot146(&robot146,0);
      PlayerClient    robot147("localhost", 7147);
      Position2dProxy p2dProxy_robot147(&robot147,0);
      RangerProxy      laserProxy_robot147(&robot147,0);
      FiducialProxy neighbor_finderProxy_robot147(&robot147,0);
      PlayerClient    robot148("localhost", 7148);
      Position2dProxy p2dProxy_robot148(&robot148,0);
      RangerProxy      laserProxy_robot148(&robot148,0);
      FiducialProxy neighbor_finderProxy_robot148(&robot148,0);
      PlayerClient    robot149("localhost", 7149);
      Position2dProxy p2dProxy_robot149(&robot149,0);
      RangerProxy      laserProxy_robot149(&robot149,0);
      FiducialProxy neighbor_finderProxy_robot149(&robot149,0);
      PlayerClient    robot150("localhost", 7150);
      Position2dProxy p2dProxy_robot150(&robot150,0);
      RangerProxy      laserProxy_robot150(&robot150,0);
      FiducialProxy neighbor_finderProxy_robot150(&robot150,0);
      PlayerClient    robot151("localhost", 7151);
      Position2dProxy p2dProxy_robot151(&robot151,0);
      RangerProxy      laserProxy_robot151(&robot151,0);
      FiducialProxy neighbor_finderProxy_robot151(&robot151,0);
      PlayerClient    robot152("localhost", 7152);
      Position2dProxy p2dProxy_robot152(&robot152,0);
      RangerProxy      laserProxy_robot152(&robot152,0);
      FiducialProxy neighbor_finderProxy_robot152(&robot152,0);
      PlayerClient    robot153("localhost", 7153);
      Position2dProxy p2dProxy_robot153(&robot153,0);
      RangerProxy      laserProxy_robot153(&robot153,0);
      FiducialProxy neighbor_finderProxy_robot153(&robot153,0);
      PlayerClient    robot154("localhost", 7154);
      Position2dProxy p2dProxy_robot154(&robot154,0);
      RangerProxy      laserProxy_robot154(&robot154,0);
      FiducialProxy neighbor_finderProxy_robot154(&robot154,0);
      PlayerClient    robot155("localhost", 7155);
      Position2dProxy p2dProxy_robot155(&robot155,0);
      RangerProxy      laserProxy_robot155(&robot155,0);
      FiducialProxy neighbor_finderProxy_robot155(&robot155,0);
      PlayerClient    robot156("localhost", 7156);
      Position2dProxy p2dProxy_robot156(&robot156,0);
      RangerProxy      laserProxy_robot156(&robot156,0);
      FiducialProxy neighbor_finderProxy_robot156(&robot156,0);
      PlayerClient    robot157("localhost", 7157);
      Position2dProxy p2dProxy_robot157(&robot157,0);
      RangerProxy      laserProxy_robot157(&robot157,0);
      FiducialProxy neighbor_finderProxy_robot157(&robot157,0);
      PlayerClient    robot158("localhost", 7158);
      Position2dProxy p2dProxy_robot158(&robot158,0);
      RangerProxy      laserProxy_robot158(&robot158,0);
      FiducialProxy neighbor_finderProxy_robot158(&robot158,0);
      PlayerClient    robot159("localhost", 7159);
      Position2dProxy p2dProxy_robot159(&robot159,0);
      RangerProxy      laserProxy_robot159(&robot159,0);
      FiducialProxy neighbor_finderProxy_robot159(&robot159,0);
      PlayerClient    robot160("localhost", 7160);
      Position2dProxy p2dProxy_robot160(&robot160,0);
      RangerProxy      laserProxy_robot160(&robot160,0);
      FiducialProxy neighbor_finderProxy_robot160(&robot160,0);
      PlayerClient    robot161("localhost", 7161);
      Position2dProxy p2dProxy_robot161(&robot161,0);
      RangerProxy      laserProxy_robot161(&robot161,0);
      FiducialProxy neighbor_finderProxy_robot161(&robot161,0);
      PlayerClient    robot162("localhost", 7162);
      Position2dProxy p2dProxy_robot162(&robot162,0);
      RangerProxy      laserProxy_robot162(&robot162,0);
      FiducialProxy neighbor_finderProxy_robot162(&robot162,0);
      PlayerClient    robot163("localhost", 7163);
      Position2dProxy p2dProxy_robot163(&robot163,0);
      RangerProxy      laserProxy_robot163(&robot163,0);
      FiducialProxy neighbor_finderProxy_robot163(&robot163,0);
      PlayerClient    robot164("localhost", 7164);
      Position2dProxy p2dProxy_robot164(&robot164,0);
      RangerProxy      laserProxy_robot164(&robot164,0);
      FiducialProxy neighbor_finderProxy_robot164(&robot164,0);
      PlayerClient    robot165("localhost", 7165);
      Position2dProxy p2dProxy_robot165(&robot165,0);
      RangerProxy      laserProxy_robot165(&robot165,0);
      FiducialProxy neighbor_finderProxy_robot165(&robot165,0);
      PlayerClient    robot166("localhost", 7166);
      Position2dProxy p2dProxy_robot166(&robot166,0);
      RangerProxy      laserProxy_robot166(&robot166,0);
      FiducialProxy neighbor_finderProxy_robot166(&robot166,0);
      PlayerClient    robot167("localhost", 7167);
      Position2dProxy p2dProxy_robot167(&robot167,0);
      RangerProxy      laserProxy_robot167(&robot167,0);
      FiducialProxy neighbor_finderProxy_robot167(&robot167,0);
      PlayerClient    robot168("localhost", 7168);
      Position2dProxy p2dProxy_robot168(&robot168,0);
      RangerProxy      laserProxy_robot168(&robot168,0);
      FiducialProxy neighbor_finderProxy_robot168(&robot168,0);
      PlayerClient    robot169("localhost", 7169);
      Position2dProxy p2dProxy_robot169(&robot169,0);
      RangerProxy      laserProxy_robot169(&robot169,0);
      FiducialProxy neighbor_finderProxy_robot169(&robot169,0);
      PlayerClient    robot170("localhost", 7170);
      Position2dProxy p2dProxy_robot170(&robot170,0);
      RangerProxy      laserProxy_robot170(&robot170,0);
      FiducialProxy neighbor_finderProxy_robot170(&robot170,0);
      PlayerClient    robot171("localhost", 7171);
      Position2dProxy p2dProxy_robot171(&robot171,0);
      RangerProxy      laserProxy_robot171(&robot171,0);
      FiducialProxy neighbor_finderProxy_robot171(&robot171,0);
      PlayerClient    robot172("localhost", 7172);
      Position2dProxy p2dProxy_robot172(&robot172,0);
      RangerProxy      laserProxy_robot172(&robot172,0);
      FiducialProxy neighbor_finderProxy_robot172(&robot172,0);
      PlayerClient    robot173("localhost", 7173);
      Position2dProxy p2dProxy_robot173(&robot173,0);
      RangerProxy      laserProxy_robot173(&robot173,0);
      FiducialProxy neighbor_finderProxy_robot173(&robot173,0);
      PlayerClient    robot174("localhost", 7174);
      Position2dProxy p2dProxy_robot174(&robot174,0);
      RangerProxy      laserProxy_robot174(&robot174,0);
      FiducialProxy neighbor_finderProxy_robot174(&robot174,0);
      PlayerClient    robot175("localhost", 7175);
      Position2dProxy p2dProxy_robot175(&robot175,0);
      RangerProxy      laserProxy_robot175(&robot175,0);
      FiducialProxy neighbor_finderProxy_robot175(&robot175,0);
      PlayerClient    robot176("localhost", 7176);
      Position2dProxy p2dProxy_robot176(&robot176,0);
      RangerProxy      laserProxy_robot176(&robot176,0);
      FiducialProxy neighbor_finderProxy_robot176(&robot176,0);
      PlayerClient    robot177("localhost", 7177);
      Position2dProxy p2dProxy_robot177(&robot177,0);
      RangerProxy      laserProxy_robot177(&robot177,0);
      FiducialProxy neighbor_finderProxy_robot177(&robot177,0);
      PlayerClient    robot178("localhost", 7178);
      Position2dProxy p2dProxy_robot178(&robot178,0);
      RangerProxy      laserProxy_robot178(&robot178,0);
      FiducialProxy neighbor_finderProxy_robot178(&robot178,0);
      PlayerClient    robot179("localhost", 7179);
      Position2dProxy p2dProxy_robot179(&robot179,0);
      RangerProxy      laserProxy_robot179(&robot179,0);
      FiducialProxy neighbor_finderProxy_robot179(&robot179,0);

      // vector<PlayerClient> plyclnts(5);
      // vector<Position2dProxy> p2dProxys(5);
      // vector<RangerProxy> rngProxys(5);
      // vector<FiducialProxy> fidProxys(5);
      // vector<SimulationProxy> SimProxys(5);
// spawn off models
      PlayerClient plyclnts[SWARM_SIZE] ={robot120, robot121, robot122, robot123, robot124, robot125, robot126, robot127, robot128, robot129, robot130, robot131, robot132, robot133, robot134, robot135, robot136, robot137, robot138, robot139, robot140, robot141, robot142, robot143, robot144, robot145, robot146, robot147, robot148, robot149, robot150, robot151, robot152, robot153, robot154, robot155, robot156, robot157, robot158, robot159, robot160, robot161, robot162, robot163, robot164, robot165, robot166, robot167, robot168, robot169, robot170, robot171, robot172, robot173, robot174, robot175, robot176, robot177, robot178, robot179};
      Position2dProxy p2dProxys[SWARM_SIZE] = {p2dProxy_robot120, p2dProxy_robot121, p2dProxy_robot122, p2dProxy_robot123, p2dProxy_robot124, p2dProxy_robot125, p2dProxy_robot126, p2dProxy_robot127, p2dProxy_robot128, p2dProxy_robot129, p2dProxy_robot130, p2dProxy_robot131, p2dProxy_robot132, p2dProxy_robot133, p2dProxy_robot134, p2dProxy_robot135, p2dProxy_robot136, p2dProxy_robot137, p2dProxy_robot138, p2dProxy_robot139, p2dProxy_robot140, p2dProxy_robot141, p2dProxy_robot142, p2dProxy_robot143, p2dProxy_robot144, p2dProxy_robot145, p2dProxy_robot146, p2dProxy_robot147, p2dProxy_robot148, p2dProxy_robot149, p2dProxy_robot150, p2dProxy_robot151, p2dProxy_robot152, p2dProxy_robot153, p2dProxy_robot154, p2dProxy_robot155, p2dProxy_robot156, p2dProxy_robot157, p2dProxy_robot158, p2dProxy_robot159, p2dProxy_robot160, p2dProxy_robot161, p2dProxy_robot162, p2dProxy_robot163, p2dProxy_robot164, p2dProxy_robot165, p2dProxy_robot166, p2dProxy_robot167, p2dProxy_robot168, p2dProxy_robot169, p2dProxy_robot170, p2dProxy_robot171, p2dProxy_robot172, p2dProxy_robot173, p2dProxy_robot174, p2dProxy_robot175, p2dProxy_robot176, p2dProxy_robot177, p2dProxy_robot178, p2dProxy_robot179};
      RangerProxy rngProxys[SWARM_SIZE] = {laserProxy_robot120, laserProxy_robot121, laserProxy_robot122, laserProxy_robot123, laserProxy_robot124, laserProxy_robot125, laserProxy_robot126, laserProxy_robot127, laserProxy_robot128, laserProxy_robot129, laserProxy_robot130, laserProxy_robot131, laserProxy_robot132, laserProxy_robot133, laserProxy_robot134, laserProxy_robot135, laserProxy_robot136, laserProxy_robot137, laserProxy_robot138, laserProxy_robot139, laserProxy_robot140, laserProxy_robot141, laserProxy_robot142, laserProxy_robot143, laserProxy_robot144, laserProxy_robot145, laserProxy_robot146, laserProxy_robot147, laserProxy_robot148, laserProxy_robot149, laserProxy_robot150, laserProxy_robot151, laserProxy_robot152, laserProxy_robot153, laserProxy_robot154, laserProxy_robot155, laserProxy_robot156, laserProxy_robot157, laserProxy_robot158, laserProxy_robot159, laserProxy_robot160, laserProxy_robot161, laserProxy_robot162, laserProxy_robot163, laserProxy_robot164, laserProxy_robot165, laserProxy_robot166, laserProxy_robot167, laserProxy_robot168, laserProxy_robot169, laserProxy_robot170, laserProxy_robot171, laserProxy_robot172, laserProxy_robot173, laserProxy_robot174, laserProxy_robot175, laserProxy_robot176, laserProxy_robot177, laserProxy_robot178, laserProxy_robot179};
      FiducialProxy fidProxys[SWARM_SIZE]={neighbor_finderProxy_robot120, neighbor_finderProxy_robot121, neighbor_finderProxy_robot122, neighbor_finderProxy_robot123, neighbor_finderProxy_robot124, neighbor_finderProxy_robot125, neighbor_finderProxy_robot126, neighbor_finderProxy_robot127, neighbor_finderProxy_robot128, neighbor_finderProxy_robot129, neighbor_finderProxy_robot130, neighbor_finderProxy_robot131, neighbor_finderProxy_robot132, neighbor_finderProxy_robot133, neighbor_finderProxy_robot134, neighbor_finderProxy_robot135, neighbor_finderProxy_robot136, neighbor_finderProxy_robot137, neighbor_finderProxy_robot138, neighbor_finderProxy_robot139, neighbor_finderProxy_robot140, neighbor_finderProxy_robot141, neighbor_finderProxy_robot142, neighbor_finderProxy_robot143, neighbor_finderProxy_robot144, neighbor_finderProxy_robot145, neighbor_finderProxy_robot146, neighbor_finderProxy_robot147, neighbor_finderProxy_robot148, neighbor_finderProxy_robot149, neighbor_finderProxy_robot150, neighbor_finderProxy_robot151, neighbor_finderProxy_robot152, neighbor_finderProxy_robot153, neighbor_finderProxy_robot154, neighbor_finderProxy_robot155, neighbor_finderProxy_robot156, neighbor_finderProxy_robot157, neighbor_finderProxy_robot158, neighbor_finderProxy_robot159, neighbor_finderProxy_robot160, neighbor_finderProxy_robot161, neighbor_finderProxy_robot162, neighbor_finderProxy_robot163, neighbor_finderProxy_robot164, neighbor_finderProxy_robot165, neighbor_finderProxy_robot166, neighbor_finderProxy_robot167, neighbor_finderProxy_robot168, neighbor_finderProxy_robot169, neighbor_finderProxy_robot170, neighbor_finderProxy_robot171, neighbor_finderProxy_robot172, neighbor_finderProxy_robot173, neighbor_finderProxy_robot174, neighbor_finderProxy_robot175, neighbor_finderProxy_robot176, neighbor_finderProxy_robot177, neighbor_finderProxy_robot178, neighbor_finderProxy_robot179};
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
