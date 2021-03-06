
#include <iostream>
#include <visualization.h>

#include <Eigen/Dense>

#include <dV_spring_particle_particle_dq.h>
#include <d2V_spring_particle_particle_dq2.h>
#include <forward_euler.h>
#include <backward_euler.h>
#include <symplectic_euler.h>
#include <runge_kutta.h>

//igl
#include <igl/readOBJ.h>

//simulation parameters
Eigen::VectorXd q;              // q为位移
Eigen::VectorXd q_dot;          // q 的一阶导，可认为是 v

Eigen::VectorXd l;              // 弹簧原长

//VectorXd: 动态长度double型列向量（行数目任意，列数为1的double矩阵）

double mass = 1.0;              // 物体质量
double stiffness = 100.0;       // 硬度，胡克定律的k
double dt = 1e-2;               // 时间微分尺度
int integrator_type = 0;        // 默认积分方法

bool simulate(igl::opengl::glfw::Viewer & viewer) {
    
    //take a time step

    // 函数force：使用 k(即stiffness)和 x(此处即q-l)计算 f=-kx ，以1x1矩阵的形式存储在引用的第一个参数中
    //    ·force以建立一个lambda函数并存储的方式实现
    auto force = [](Eigen::VectorXd &f, const Eigen::VectorXd &q, const Eigen::VectorXd &qdot) { 
        dV_spring_particle_particle_dq(f, q, stiffness, l); 
        f *= -1; 
    };

    // 函数stiff：获得 -k(即stiffness)，以1x1矩阵的形式存储在引用的第一个参数中
    auto stiff = [](Eigen::MatrixXd &k, const Eigen::VectorXd &q, const Eigen::VectorXd &qdot) { 
        d2V_spring_particle_particle_dq2(k, q, stiffness); 
        k *= -1; 
    };

    switch(integrator_type) {
        case 0:
            forward_euler(q,q_dot, dt, mass, force);
            break;
        case 1:
            backward_euler(q,q_dot, dt, mass, force, stiff);
            break;
        case 2:
            symplectic_euler(q,q_dot, dt, mass, force);
            break;
        case 3: 
            runge_kutta(q, q_dot, dt, mass, force);
            break;
    };

    //update mesh positions
    Visualize::rigid_transform_1d(0, q(0));
    Visualize::scale_x(1, q(0));

    return false;
}

int main(int argc, char **argv) {

    std::cout<<"Start A1 \n";

    //check argument for integrator type
    if(argc > 1) {
       if(argv[1] == std::string("be")) { integrator_type = 1; }
       if(argv[1] == std::string("se")) { integrator_type = 2; }
       if(argv[1] == std::string("rk")) { integrator_type = 3; }

     }

    //Load data for animations
    Eigen::MatrixXd V_cow, V_spring;
    Eigen::MatrixXi F_cow, F_spring;

    igl::readOBJ("../data/spot.obj", V_cow, F_cow);
    igl::readOBJ("../data/spring.obj", V_spring, F_spring);

    //setup simulation variables
    q.resize(1);
    q_dot.resize(1);
    l.resize(1);

    q_dot(0) = 0;                   // 初始速度=0m/s
    q(0) = 1;                       // 初始位移=1m
    l(0) = 0;                       // 弹簧原长默认为0m
    if(argc > 2) {
        l(0) = atoi(argv[2]);      // 弹簧原长可由第二个参数传入
    }
    //setup libigl viewer and activate 
    Visualize::setup(q, q_dot);
    Visualize::viewer().callback_post_draw = &simulate;
    Visualize::add_object_to_scene(V_cow, F_cow, Eigen::RowVector3d(244,165,130)/255.);
    Visualize::add_object_to_scene(V_spring, F_spring, Eigen::RowVector3d(200,200,200)/255.);
    Visualize::viewer().launch();
    //std::cout << q(0) << " " << qdot(0);
    //how am I going to organize memory
    //seperate q and q Dot arrays (extra stuff to pass around) <-- this is easier to deal with 
    return 0;
}