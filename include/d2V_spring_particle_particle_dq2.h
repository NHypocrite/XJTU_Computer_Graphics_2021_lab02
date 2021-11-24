#include <Eigen/Dense>
// 获得 k(即stiffness)，以1x1矩阵的形式存储在引用的第一个参数中

//Input:
//  q - the current generalized coordinates for the mass-spring system
//  stiffness - the stiffness parameter (spring constant) for the mass-spring system
//Output:
// H - the second derivtive of the potential energy 势能二阶导数，k的1x1矩阵存储形式
void d2V_spring_particle_particle_dq2(Eigen::MatrixXd &H, const Eigen::VectorXd &q, double stiffness);