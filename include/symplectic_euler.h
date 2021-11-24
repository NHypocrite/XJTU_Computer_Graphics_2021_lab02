//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//Output:
//  q - set q to the updated generalized coordinate using Symplectic Euler time integration
//  qdot - set qdot to the updated generalized velocity using Symplectic Euler time integration

template<typename FORCE> 
inline void symplectic_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force) {
	Eigen::VectorXd f;                  // 用于存储力f的矩阵
    Eigen::VectorXd q_new;              // 用于存储经过dt后位移q的值的矩阵
    Eigen::VectorXd qdot_new;           // 用于存储经过dt后速度qdot的值的矩阵
	force(f, q, qdot);                  // 计算这时的弹簧给球的力f=-kx，存储在f中

    qdot_new = (f / mass) * dt + qdot;  // v1 = a * dt + v0, 前向欧拉方法
    q_new = qdot_new * dt + q;          // x1 = v1 * dt + x0，后向欧拉方法


    q = q_new;                          // 将新的q通过引用传出
    qdot = qdot_new;                    // 将新的qdot通过引用传出
}