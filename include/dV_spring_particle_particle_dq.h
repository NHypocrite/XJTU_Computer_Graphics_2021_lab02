#include <Eigen/Dense>
// 使用 k(即stiffness)和 x(即q-l)计算 f ，存储在引用的第一个参数中

//Input:
//  q - the current generalized coordinates for the mass-spring system
//  stiffness - the stiffness parameter (spring constant) for the mass-spring system
//Output:
// dV - the gradient of the potential energy with respect to the generalised coordinates. 
// dV就是势能对x的一阶导数。-> f=kx 胡克定律 (此处还没考虑方向)
void dV_spring_particle_particle_dq(Eigen::VectorXd &dV, const Eigen::VectorXd &q, double stiffness, const Eigen::VectorXd l);