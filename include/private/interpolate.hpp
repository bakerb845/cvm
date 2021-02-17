#ifndef CVM_PRIVATE_INTERPOLATE_HPP
#define CVM_PRIVATE_INTERPOLATE_HPP
namespace
{
/// @brief Simple bilinear interpolation routine.  For a picture see:
///        https://en.wikipedia.org/wiki/Bilinear_interpolation
/// @param[in] x1    The x location of left grid points.
/// @param[in] x2    The x location of right grid points.
/// @param[in] y1    The y location of lower grid points.
/// @param[in] y2    The y location of upper grid points.
/// @param[in] x     The x interpolation location on interval [x1, x2].
/// @param[in] y     The y interpolation location on interval [y1, y2].
/// @param[in] v4    Values at grid points packed [(1,1), (1,2), (2,1), (2,2)]
///                  where (1,1) is the lower left corner, (1,2) is the
///                  the lower right corner, (2,1) is the upper left corner,
///                  and (2,2) is the upper right corner.
/// @result The bilinear interpolatant at (x,y) in the square.
template<class T>
[[maybe_unused]] [[nodiscard]]
T bilinearInterpolation(const T x1, const T x2, 
                        const T y1, const T y2, 
                        const T x,  const T y,
                        const T *v4)
{
    const T one = 1;
    T dx = x2 - x1; 
    T dy = y2 - y1; 
    T det = one/(dx*dy);
    T fq11 = v4[0];
    T fq21 = v4[1];
    T fq12 = v4[2];
    T fq22 = v4[3];
    T fxy = det*( fq11*(x2 - x)*(y2 - y) + fq21*(x - x1)*(y2 - y)
                + fq12*(x2 - x)*(y - y1) + fq22*(x - x1)*(y - y1));
    return fxy;
}
/// @brief Simple trilinear interpolation routine.  For a picture see:
///        https://en.wikipedia.org/wiki/Trilinear_interpolation.
/// @param[in] x0  The x location of the the left grid points.
/// @param[in] y0  The y location of the closest grid points.
/// @param[in] z0  The z location of the bottom grid points.
/// @param[in] x1  The x location of the the right grid points.
/// @param[in] y1  The y location of the furthest grid points.
/// @param[in] z1  The z location of the top grid points.
/// @param[in] x   The x interpolation location on interval [x0, x1].
/// @param[in] y   The y interpolation location on interval [y0, y1].
/// @param[in] z   The z interpolation location on interval [z0, z1].
/// @param[in] v8  The values at
///                [c000, c100, c010, c110, c001, c101, c011, c111]
/// @result The interpolated value at (x,y,z).
template<class T>
[[maybe_unused]] [[nodiscard]]
T trilinearInterpolation(const T x0, const T x1,
                         const T y0, const T y1,
                         const T z0, const T z1,
                         const T x, const T y, const T z,
                         const T *v8)
{
    T dx = x1 - x0; 
    T dy = y1 - y0;
    T dz = z1 - z0;
    T xd = (x - x0)/dx;
    T yd = (y - y0)/dy;
    T zd = (z - z0)/dz;
    T c000 = v8[0];
    T c100 = v8[1];
    T c010 = v8[2];
    T c110 = v8[3];
    T c001 = v8[4];
    T c101 = v8[5];
    T c011 = v8[6];
    T c111 = v8[7];
    // Interpolate in x
    T c00 = c000*(1 - xd) + c100*xd;
    T c01 = c001*(1 - xd) + c101*xd;
    T c10 = c010*(1 - xd) + c110*xd;
    T c11 = c011*(1 - xd) + c111*xd;
    // Interpolate in y
    T c0 = c00*(1 - yd) + c10*yd;
    T c1 = c01*(1 - yd) + c11*yd;
    // Interpolate in z
    T c = c0*(1 - zd) + c1*zd;
    return c; 
}
}
#endif
