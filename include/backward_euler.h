#include <Eigen/Dense>

//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//  stiffness(q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters.  
//Output:
//  q - set q to the updated generalized coordinate using Backward Euler time integration
//  qdot - set qdot to the updated generalized velocity using Backward Euler time integration

template<typename FORCE, typename STIFFNESS> 
inline void backward_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force, STIFFNESS &stiffness) {
    Eigen::VectorXd f;                  // 用于存储力f的矩阵
    Eigen::MatrixXd k;                  // 用于存储弹性系数k的矩阵
    Eigen::VectorXd q_new;              // 用于存储经过dt后位移q的值的矩阵
    Eigen::VectorXd qdot_new;           // 用于存储经过dt后速度qdot的值的矩阵   

    force(f, q, qdot);                  // 计算这时的弹簧给球的力f=-kx，存储在f中
    stiffness(k, q, qdot);              // 计算这时的弹簧给球的力k=-k，存储在k中

    qdot_new = (qdot + (f / mass) * dt) / (1 - dt * dt * k(0,0) / mass);  // v1 = [v0 + (a * dt)]/[1-(k*dt^2)/m] 
    q_new = qdot_new * dt + q;                                            // x1 = v1 * dt + x0

    q = q_new;                          // 将新的q通过引用传出
    qdot = qdot_new;                    // 将新的qdot通过引用传出    

}