#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <Eigen/Dense>

#include <istream>

/**
 * \defgroup TOOLS Other tools used in the library
 *
 * \{
 */

/*!
 * \brief fRand generates a pseudorandom number with a uniform distribution between fMin and fMax.
 * \param fMin minimum value.
 * \param fMax maximum value.
 * \return pseudorandom value.
 * \callgraph
 * \callergraph
 */
template<typename T>
T fRand(T fMin, T fMax)
{
    T f = (T)rand() / T(RAND_MAX);
    return fMin + f * (fMax - fMin);
}

/*!
 * \brief randomPose generates a random pose.
 * \param min_translation
 * \param max_translation
 * \return pose
 * \callgraph
 * \callergraph
 */
template<typename T>
Eigen::Matrix<T,6,1> randomPose(T min_translation, T max_translation)
{
    Eigen::Matrix<T,6,1> x;
    x<<fRand(min_translation,max_translation),fRand(min_translation,max_translation),fRand(min_translation,max_translation),
       fRand(T(-M_PI),T(M_PI)),fRand(T(-M_PI),T(M_PI)),fRand(T(-M_PI),T(M_PI));
    return x;
}

/*!
 * \brief Operator << so quternions can be used with a std::ostream
 */
template <typename T>
std::ostream& operator<<(std::ostream& os, const Eigen::Quaternion<T>& q)
{
    os << "[" << q.w() <<","<< q.x() << "," << q.y() << "," << q.z() << "]";
    return os;
}

/** \} */

#endif // FUNCTIONS_H
