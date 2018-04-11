#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H

#include <Eigen/Core>
#include <exception>

//#define NEIRA

/**
 * \mainpage
 *
 * This package implemnts the transformations and their jacobians from \cite Neira1993 and \cite Smith1990. The jacobian of the inverse was not working
 * and has been implemented differently. Please, see <a href="equations.pdf">this document</a> for differentiation details. If Ceres Solver \cite ceres-solver is available a test
 * program comparing the analitical derivatives with the numeric and automatic differentiation from Ceres Solver is compiled.
 *
 * \author Albert Palomer Vila
 *
 */


/**
 * \defgroup TRANSFORMATIONS Transformation implementation of the transformation and its derivatives.
 *
 * \{
 */

/*!
 * \brief Angle between two quaternions. Implemented from <a href="https://math.stackexchange.com/questions/90081/quaternion-distance">Mathematics (Stack Exchange): Quaternion distance</a>
 */
template<typename T>
T quaternionAngle(Eigen::Quaternion<T> q1, Eigen::Quaternion<T> q2)
{
    /// https://math.stackexchange.com/questions/90081/quaternion-distance
    return acos(2*pow(q1.dot(q2),2) - 1);
}

/*!
 * \brief Distance between two quaternions. Implemented from  <a href="https://math.stackexchange.com/questions/90081/quaternion-distance">Mathematics (Stack Exchange): Quaternion distance</a>
 */
template<typename T>
T quaternionDistance(Eigen::Quaternion<T> q1, Eigen::Quaternion<T> q2)
{
    /// https://math.stackexchange.com/questions/90081/quaternion-distance
    return T(1) - pow(q1.dot(q2),T(2));
}


/*!
 * \brief Extract roll, pitch and yaw from a given rotation matrix.
 * \return 3x1 matrix containing  roll, pitch and yaw
 * \callgraph
 * \callergraph
 */
template<typename T>
Eigen::Matrix<T,3,1> getRPY(Eigen::Matrix<T,3,3> rot)
{
    return rot.eulerAngles(2,1,0).colwise().reverse();
}

/*!
 * \brief Extract roll, pitch and yaw from a given quaternion.
 * \return 3x1 matrix containing  roll, pitch and yaw
 * \callgraph
 * \callergraph
 */
template<typename T>
Eigen::Matrix<T,3,1> getRPY(Eigen::Quaternion<T> q)
{
    return getRPY<T>(q.toRotationMatrix());
}

/*!
 * \brief constructs a quaternion from the roll, pitch and yaw rotations. Applied in yaw, pitch and roll order.
 * \return rotation as a quaternion
 * \callgraph
 * \callergraph
 */
template<typename T>
Eigen::Quaternion<T> quaternionFromRPY(Eigen::Matrix<T,3,1> rpy)
{
    // Unit vectors
    Eigen::Matrix<T,3,1> zunit = Eigen::Matrix<T,3,1>::UnitZ();
    Eigen::Matrix<T,3,1> yunit = Eigen::Matrix<T,3,1>::UnitY();
    Eigen::Matrix<T,3,1> xunit = Eigen::Matrix<T,3,1>::UnitX();

    // Rotations
    Eigen::AngleAxis<T> zrot = Eigen::AngleAxis<T>(rpy(2,0),zunit);
    Eigen::AngleAxis<T> yrot = Eigen::AngleAxis<T>(rpy(1,0),yunit);
    Eigen::AngleAxis<T> xrot = Eigen::AngleAxis<T>(rpy(0,0),xunit);

    // Compute quaternion
    Eigen::Quaternion<T> q = zrot*yrot*xrot;

    // Exit
    return q;
}

/*!
 * \brief Rotation matrix from roll pitch yaw.
 * \return 3x3 rotation matrix with the rotations first arround yaw, then pitch and finally roll.
 * \callgraph
 * \callergraph
 */
template<typename T>
Eigen::Matrix<T,3,3> mrot(Eigen::Matrix<T,3,1> rpy)
{
    return quaternionFromRPY<T>(rpy).toRotationMatrix();
}

/*!
 * \brief 3D Compounding operation: x3 = x1 (+) x2.
 * \param t1 3D position (x,y,z)
 * \param q1 3D quaternion orientation
 * \param t2 3D position (x,y,z)
 * \param q2 3D quaternion orientation
 * \param[out] t3 3D position (x,y,z)
 * \param[out] q3 3D quaternion orientation
 * \callgraph
 * \callergraph
 */
template<typename T>
void compound3D(const Eigen::Matrix<T,3,1>& t1, const Eigen::Quaternion<T>& q1, const Eigen::Matrix<T,3,1>& t2, const Eigen::Quaternion<T>& q2, Eigen::Matrix<T,3,1>& t3, Eigen::Quaternion<T>& q3)
{
    t3 = t1 + q1.toRotationMatrix()*t2;
    q3 = q1.toRotationMatrix()*q2.toRotationMatrix();
}

/*!
 * \brief 3D Compounding operation: x3 = x1 (+) x2.
 * \param t1 3D position (x,y,z)
 * \param rpy1 3D orientation (roll,pitch,yaw)
 * \param t2 3D position (x,y,z)
 * \param rpy2 3D orientation (roll,pitch,yaw)
 * \param[out] t3 3D position (x,y,z)
 * \param[out] rpy3 3D orientation (roll,pitch,yaw)
 * \callgraph
 * \callergraph
 */
template<typename T>
void compound3D(const Eigen::Matrix<T,3,1>& t1, const Eigen::Matrix<T,3,1>& rpy1, const Eigen::Matrix<T,3,1>& t2, const Eigen::Matrix<T,3,1>& rpy2, Eigen::Matrix<T,3,1>& t3, Eigen::Matrix<T,3,1>& rpy3)
{
    Eigen::Quaternion<T> q1 = quaternionFromRPY<T>(rpy1);
    Eigen::Quaternion<T> q2 = quaternionFromRPY<T>(rpy2);
    Eigen::Quaternion<T> q3;
    compound3D<T>(t1,q1,t2,q2,t3,q3);
    rpy3 = getRPY<T>(q3);
}

/*!
 * \brief 3D Compounding operation: x3 = x1 (+) x2.
 * \param x1 3D position and orientation (x,y,z,roll,pitch,yaw)
 * \param x2 3D position and orientation (x,y,z,roll,pitch,yaw)
 * \return x3 = x1 (+) x2. Position and orientation (x,y,z,roll,pitch,yaw) of x2 with respect to the frame here x1 is refferenced.
 * \callgraph
 * \callergraph
 */
template<typename T>
Eigen::Matrix<T,6,1> compound3D(const Eigen::Matrix<T,6,1>& x1, const Eigen::Matrix<T,6,1>& x2)
{
    Eigen::Matrix<T,3,1> t1 = x1.block(0,0,3,1);
    Eigen::Matrix<T,3,1> rpy1 = x1.block(3,0,3,1);
    Eigen::Matrix<T,3,1> t2 = x2.block(0,0,3,1);
    Eigen::Matrix<T,3,1> rpy2 = x2.block(3,0,3,1);
    Eigen::Matrix<T,3,1> t3;
    Eigen::Matrix<T,3,1> rpy3;
    compound3D<T>(t1,rpy1,t2,rpy2,t3,rpy3);
    Eigen::Matrix<T,6,1> x3 = Eigen::Matrix<T,6,1>::Zero();
    x3.block(0,0,3,1) = t3;
    x3.block(3,0,3,1) = rpy3;
    return x3;
}

/*!
 * Computes the Jacobians (6x6) of the composition x1 (+) x2 where J1 = d(x1 (+) x2)/dx1 and J2 = d(x1 (+) x2)/dx2
 * \param x1 3D position and orientation (x,y,z,roll,pitch,yaw)
 * \param x2 3D position and orientation (x,y,z,roll,pitch,yaw)
 * \param J1 jacobian with respect to x1
 * \param J2 jacobian with respect to x2
 */
template<typename T>
void compound3DJacobian(const Eigen::Matrix<T,6,1>& x1, const Eigen::Matrix<T,6,1>& x2, Eigen::Matrix<T,6,6>& J1, Eigen::Matrix<T,6,6>& J2)
{
    Eigen::Matrix<T,3,1> t1 = x1.block(0,0,3,1);
    Eigen::Matrix<T,3,1> rpy1 = x1.block(3,0,3,1);
    Eigen::Matrix<T,3,1> t2 = x2.block(0,0,3,1);
    Eigen::Matrix<T,3,1> rpy2 = x2.block(3,0,3,1);
    T t1x = t1(0,0);
    T t1y = t1(1,0);
    T t1z = t1(2,0);
    T t1roll = rpy1(0,0);
    T t1pitch = rpy1(1,0);
    T t1yaw = rpy1(2,0);
    T t2x = t2(0,0);
    T t2y = t2(1,0);
    T t2z = t2(2,0);
    T t2roll = rpy2(0,0);
    T t2pitch = rpy2(1,0);

    // Compute composition
    Eigen::Matrix<T,3,1> t3;
    Eigen::Matrix<T,3,1> rpy3;
    compound3D<T>(t1,rpy1,t2,rpy2,t3,rpy3);
    T t3x = t3(0,0);
    T t3y = t3(1,0);
    T t3z = t3(2,0);
    T t3roll = rpy3(0,0);
    T t3pitch = rpy3(1,0);
    T t3yaw = rpy3(2,0);
    if (cos(t3pitch) == 0)
    {
        throw std::runtime_error("Composition pitch angle is PI/2 or -PI/2 which causes a division by zero in the computation of the jacobian of the composition.");
    }

    // Rotation Matrix
    Eigen::Matrix<T,3,3> R1 = mrot(rpy1);
    T o1x = R1(0,1);
    T o1y = R1(1,1);
    T o1z = R1(2,1);
    T a1x = R1(0,2);
    T a1y = R1(1,2);
    T a1z = R1(2,2);
    Eigen::Matrix<T,3,3> R2 = mrot(rpy2);
    T o2x = R2(0,1);
    T a2x = R2(0,2);

    // Top right block of J1
    Eigen::Matrix<T,3,3> M;
    M(0,0) = t2y *a1x - t2z*o1x;
    M(0,1) = (t3z - t1z)*cos(t1yaw);
    M(0,2) = t1y - t3y;
    M(1,0) = t2y*a1y - t2z*o1y;
    M(1,1) = (t3z - t1z)*sin(t1yaw);
    M(1,2) = t3x - t1x;
    M(2,0) = t2y*a1z-t2z*o1z;
    M(2,1) = -t2x*cos(t1pitch) - t2y*sin(t1pitch)*sin(t1roll) - t2z*sin(t1pitch)*cos(t1roll);
    M(2,2) = 0;

    // Bottom right block of J1
    Eigen::Matrix<T,3,3> K1;
    K1(0,0) = cos(t1pitch)*cos(t3yaw-t1yaw)/cos(t3pitch);
    K1(0,1) = sin(t3yaw - t1yaw)/cos(t3pitch);
    K1(0,2) = 0;
    K1(1,0) = -cos(t1pitch)*sin(t3yaw-t1yaw);
    K1(1,1) = cos(t3yaw-t1yaw);
    K1(1,2) = 0;
    K1(2,0) = (o2x*sin(t3roll)+a2x*cos(t3roll))/cos(t3pitch);
    K1(2,1) = sin(t3pitch)*sin(t3yaw-t1yaw)/cos(t3pitch);
    K1(2,2) = 1;

    // Bottom right block of J2
    Eigen::Matrix<T,3,3> K2;
    K2(0,0) = 1;
    K2(0,1) = sin(t3pitch)*sin(t3roll-t2roll)/cos(t3pitch);
    K2(0,2) = (a1x*cos(t3yaw)+a1y*sin(t3yaw))/cos(t3pitch);
    K2(1,0) = 0;
    K2(1,1) = cos(t3roll - t2roll);
    K2(1,2) = -cos(t2pitch)*sin(t3roll-t2roll);
    K2(2,0) = 0;
    K2(2,1) = sin(t3roll-t2roll)/cos(t3pitch);
    K2(2,2) = cos(t2pitch)*cos(t3roll-t2roll)/cos(t3pitch);

    // Initialize Jacobians
    J1 = Eigen::Matrix<T,6,6>::Identity();
    J2 = Eigen::Matrix<T,6,6>::Identity();
    J1.block(0,3,3,3) = M;
    J1.block(3,3,3,3) = K1;
    J2.block(0,0,3,3) = R1;
    J2.block(3,3,3,3) = K2;
}

/*!
 * \brief 3D Inverse compounding operation: xi = (-) x.
 * \param t 3D position (x,y,z)
 * \param q 3D quaternion orientation
 * \param[out] ti 3D inverted position (x,y,z)
 * \param[out] qi 3D inverted quaternion orientation
 * \callgraph
 * \callergraph
 */
template<typename T>
void inverseCompound3D(const Eigen::Matrix<T,3,1>& t, const Eigen::Quaternion<T>& q, Eigen::Matrix<T,3,1>& ti, Eigen::Quaternion<T>& qi)
{
    qi = q.toRotationMatrix().transpose();
    ti = -qi.toRotationMatrix()*t;
}

/*!
 * \brief 3D Inverse compounding operation: xi = (-) x.
 * \param t 3D position (x,y,z)
 * \param rpy 3D orientation (roll,pitch,yaw)
 * \param[out] ti 3D inverted position (x,y,z)
 * \param[out] rpyi 3D inverted orientation (roll,pitch,yaw)
 * \callgraph
 * \callergraph
 */
template<typename T>
void inverseCompound3D(const Eigen::Matrix<T,3,1>& t, const Eigen::Matrix<T,3,1>& rpy, Eigen::Matrix<T,3,1>& ti, Eigen::Matrix<T,3,1>& rpyi)
{
    Eigen::Quaternion<T> q = quaternionFromRPY<T>(rpy);
    Eigen::Quaternion<T> qi;
    inverseCompound3D<T>(t,q,ti,qi);
    rpyi = getRPY<T>(qi);
}

/*!
 * \brief 3D Inverse compounding operation: xi = (-) x.
 * \return xi = (-) x. Position and orientation (x,y,z,roll,pitch,yaw) of the inverse of x.
 * \callgraph
 * \callergraph
 */
template<typename T>
Eigen::Matrix<T,6,1> inverseCompound3D(const Eigen::Matrix<T,6,1>& x)
{
    Eigen::Matrix<T,3,1> t = x.block(0,0,3,1);
    Eigen::Matrix<T,3,1> rpy = x.block(3,0,3,1);
    Eigen::Matrix<T,3,1> ti;
    Eigen::Matrix<T,3,1> rpyi;
    inverseCompound3D<T>(t,rpy,ti,rpyi);
    Eigen::Matrix<T,6,1> xi;
    xi.block(0,0,3,1) = ti;
    xi.block(3,0,3,1) = rpyi;
    return xi;
}

/*!
 * \brief Computes the Jacobian (6x6) of the pose x.
 */
template<typename T>
Eigen::Matrix<T,6,6> inverseCompound3DJacobian(const Eigen::Matrix<T,6,1>& x)
{
    // Compute inverse
    Eigen::Matrix<T,6,1> xi = inverseCompound3D(x);

    // Split transformation
    Eigen::Matrix<T,3,1> t = x.block(0,0,3,1);
    T t_x = t(0,0);
    T t_y = t(1,0);
    T t_z = t(2,0);
    Eigen::Matrix<T,3,1> rpy = x.block(3,0,3,1);

    // Compute rotation matrix
    Eigen::Matrix<T,3,3> rot = mrot(rpy);

    // Auxiliary computations
    Eigen::Matrix<T,3,1> n = rot.col(0);
    T n_x = n(0,0);
    T n_y = n(1,0);
    T n_z = n(2,0);
    Eigen::Matrix<T,3,1> o = rot.col(1);
    T o_x = o(0,0);
    T o_y = o(1,0);
    T o_z = o(2,0);
    Eigen::Matrix<T,3,1> a = rot.col(2);
    T a_x = a(0,0);
    T a_y = a(1,0);
    T a_z = a(2,0);
    T cr = cos(rpy(0,0));
    T sr = sin(rpy(0,0));
    T cp = cos(rpy(1,0));
    T sp = sin(rpy(1,0));
    T cy = cos(rpy(2,0));
    T sy = sin(rpy(2,0));

    // Compute N
    Eigen::Matrix<T,3,3> N = Eigen::Matrix<T,3,3>::Zero();
    N(0,1) = -n_z*t_x*cy - n_z*t_y*sy + t_z*cp;
    N(0,2) = n_y*t_x - n_x*t_y;
    N(1,0) = xi(2,0);
    N(1,1) = -o_z*t_x*cy - o_z*t_y*sy + t_z*sp*sr;
    N(1,2) = o_y*t_x - o_x*t_y;
    N(2,0) = -xi(1,0);
    N(2,1) = -a_z*t_x*cy - a_z*t_y*sy + t_z*sp*cr;
    N(2,2) = a_y*t_x - a_x*t_y;

#ifdef NEIRA
    // Check a_x
    if (a_x == T(1))
    {
        std::cout<<"Rotation matrix: "<<rot<<std::endl;
        throw std::runtime_error("a_x of the rotation matris is 1 (Rot = (n o a)) and cannot be used to compute the jacobian of the inverse compounding operation.");
    }

    // Compute Q
    Eigen::Matrix<T,3,3> Q = Eigen::Matrix<T,3,3>::Zero();
    T aux = T(1) - pow(a_x,2);
    T aux2 = sqrt(aux);
    Q(0,0) = -n_x/aux;
    Q(0,1) = -o_x*cr/aux;
    Q(0,2) = a_z*a_x/aux;
    Q(1,0) = o_x/aux2;
    Q(1,1) = -a_z*cy/aux2;
    Q(1,2) = a_y/aux2;
    Q(2,0) = n_x*a_x/aux;
    Q(2,1) = -a_y*cy/aux;
    Q(2,2) = -a_z/aux;
#else

    // Compute rotation matrix
    Eigen::Matrix<T,3,3> roti = rot.transpose();

    // Auxiliary computations
    Eigen::Matrix<T,3,1> ni = roti.col(0);
    T ni_x = ni(0,0);
    T ni_y = ni(1,0);
    T ni_z = ni(2,0);
    Eigen::Matrix<T,3,1> oi = roti.col(1);
    T oi_y = oi(1,0);
    T oi_z = oi(2,0);
    Eigen::Matrix<T,3,1> ai = roti.col(2);
    T ai_z = ai(2,0);

    // Derivatives
    Eigen::Matrix<T,3,1> daz;
    daz(0,0) = -cp*sr;
    daz(1,0) = -sp*cr;
    daz(2,0) = 0;
    Eigen::Matrix<T,3,1> doz;
    doz(0,0) = -oi_y;
    doz(1,0) = sy*cp*cr;
    doz(2,0) = ni_z;
    Eigen::Matrix<T,3,1> dnx;
    dnx(0,0) = 0;
    dnx(1,0) = -cy*sp;
    dnx(2,0) = -sy*cp;
    Eigen::Matrix<T,3,1> dny;
    dny(0,0) = ni_z;
    dny(1,0) = cy*cp*sr;
    dny(2,0) = -oi_y;
    Eigen::Matrix<T,3,1> dy;
    dy(0,0) = cy*sp*sr - sy*cr;
    dy(1,0) = -cy*cp*cr;
    dy(2,0) = sy*sp*cr-cy*sr;

    // Compute Q
    Eigen::Matrix<T,3,3> Q = Eigen::Matrix<T,3,3>::Zero();
    T Q000 = pow(ai_z,2)+pow(oi_z,2);
    T Q00 = -oi_z/Q000;
    T Q01 = ai_z/Q000;
    Q(0,0) = Q00*daz(0,0) + Q01*doz(0,0);
    Q(0,1) = Q00*daz(1,0) + Q01*doz(1,0);
    Q(0,2) = Q00*daz(2,0) + Q01*doz(2,0);
    T Q222 = pow(ni_x,2)+pow(ni_y,2);
    T Q20 = -ni_y/Q222;
    T Q21 = ni_x/Q222;
    Q(2,0) = Q20*dnx(0,0) + Q21*dny(0,0);
    Q(2,1) = Q20*dnx(1,0) + Q21*dny(1,0);
    Q(2,2) = Q20*dnx(2,0) + Q21*dny(2,0);

    // Compute derivatives of the new pitch angle
    Eigen::Matrix<T,3,1> dx;
    T phi_ = xi(5,0);
    T sy_ = sin(phi_);
    T cy_ = cos(phi_);
    T dxdphi_ = ni_y*cy_ - cy*cp*sy_;
    dx(0,0) = ni_z*sy_ + dxdphi_*Q(2,0);
    dx(1,0) = cy*cp*sr*sy_ - cy*sp*cy_+ dxdphi_*Q(2,1);
    dx(2,0) = -(oi_y*sy_ + sy*cp*cy_)+ dxdphi_*Q(2,2);

    // Finish jacobian
    T x_ = ni_x*cy_ + ni_y*sy_;
    T y_ = -ni_z;
    T Q111 = pow(x_,2)+pow(y_,2);
    T Q10 = -y_/Q111;
    T Q11 = x_/Q111;
    Q(1,0) = Q10*dx(0,0) + Q11*dy(0,0);
    Q(1,1) = Q10*dx(1,0) + Q11*dy(1,0);
    Q(1,2) = Q10*dx(2,0) + Q11*dy(2,0);
#endif

    // Construct the jacobian
    Eigen::Matrix<T,6,6> J = Eigen::Matrix<T,6,6>::Zero();
    J.block(0,0,3,3) = -rot.transpose();
    J.block(0,3,3,3) = N;
    J.block(3,3,3,3) = Q;

    // Exit
    return J;
}
/** \} */
#endif // TRANSFORMATIONS_H
