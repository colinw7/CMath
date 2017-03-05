#ifndef CMATH_MACROS_H
#define CMATH_MACROS_H

#ifndef M_PI
# define M_PI 3.141592654
#endif

#define PI2       6.283185407
#define PI_DIV_2  1.570796327
#define PI_DIV_4  0.785398163
#define PI3_DIV_2 4.712388980
#define PI_INV    0.318309886

#define PI_F        3.141592654f
#define PI2_F       6.283185407f
#define PI_DIV_2_F  1.570796327f
#define PI_DIV_4_F  0.785398163f
#define PI3_DIV_2_F 4.712388980f
#define PI_INV_F    0.318309886f

#define DEG_TO_RAD(a) (((a)*M_PI)/180.0)
#define RAD_TO_DEG(a) ((a)*(180.0/M_PI))

#define ROUND(r) ((int)((r) + 0.5))

#ifndef REAL_EQ
# define REAL_EQ(r1,r2) (fabs((r1)-(r2)) < 1E-5)

# define FLT_EQ(x,y) (fabs((x)-(y)) < FLT_EPSILON*fabs((x)+(y)) || fabs((x)-(y)) < FLT_MIN)
# define DBL_EQ(x,y) (fabs((x)-(y)) < DBL_EPSILON*fabs((x)+(y)) || fabs((x)-(y)) < DBL_MIN)
#endif

#endif
