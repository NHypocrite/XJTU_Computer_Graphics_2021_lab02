#include <dV_spring_particle_particle_dq.h>
// 使用 k(即stiffness)和 x(即q-l)计算 f ，存储在引用的第一个参数中
void dV_spring_particle_particle_dq(Eigen::VectorXd &dV, const Eigen::VectorXd &q, double stiffness) {
    //dV.resize(1);
    //compute f
    //dV = q * stiffness;// f = k * x
}