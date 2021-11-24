//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//Output:
//  q - set q to the updated generalized coordinate using Runge-Kutta time integration
//  qdot - set qdot to the updated generalized velocity using Runge-Kutta time integration

template<typename FORCE> 
inline void runge_kutta(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force) {
    // 在（t，t+dt）的区间里取四个点进行加权平均
    Eigen::VectorXd f1, f2, f3, f4;                 // 用于存储力f的矩阵 
    Eigen::VectorXd q_new;                          // 用于存储经过dt后位移q的值的矩阵
    Eigen::VectorXd qdot_new;                       // 用于存储经过dt后速度qdot的值的矩阵
    Eigen::VectorXd q1, q2, q3, q4;                 // 用于存储取的四个点位移q的值的矩阵
    Eigen::VectorXd qdot1, qdot2, qdot3, qdot4;     // 用于存储取的四个点速度qdot的值的矩阵
    
    // 第一个点，位于区间开始处
    q1 = q;                                 
    qdot1 = qdot;
    force(f1, q1, qdot1);

    // 第二个点，位于区间中点
    q2 = qdot1 * dt / 2 + q;                
    qdot2 = f1 / mass * dt / 2 + qdot;
    force(f2, q2, qdot2);

    // 第三个点，位于区间中点
    q3 = qdot2 * dt /2 + q;                 
    qdot3 = f2 / mass * dt / 2 + qdot;
    force(f3, q3, qdot3);

    // 第四个点，位于区间末尾/下一区间开始处
    q4 = qdot3 * dt + q;                    
    qdot4 = f3 / mass * dt + qdot;
    force(f4, q4, qdot4);

    q_new = q + dt * (qdot1 + 2*qdot2 + 2*qdot3 + qdot4) / 6;   // 最终的加权平均新位移
    qdot_new = qdot + dt * ((f1 + 2*f2 + 2*f3 + f4) / mass) / 6;// 最终的加权平均新速度

    q = q_new;                                                  // 将新的q通过引用传出
    qdot = qdot_new;                                            // 将新的qdot通过引用传出
    
}