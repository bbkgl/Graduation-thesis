#include "udf.h"

#define X_L 400          // 采空区模型的长度
#define L 50             // 基本顶破碎长度
#define Y_L 260          // 工作面宽度
#define DP 0.11          // 平均粒径
#define TIMESTAMP "05121307"   // 时间戳
/****************************************05121307修改*********************************************/
// 修改瓦斯涌出量，将浮煤区和载荷区瓦斯涌出量加大，乘以5
// 修改Z方向粘性阻力，由200---》2000

/*************************************************************************************************/

// 绝对值函数
double my_abs(double x) {
	if (x < 0) return -x;
	else return x;
}

// 单点的孔隙率函数
double porous(double x, double y, double z) {
	double value = ((1 + exp(-0.15 * (Y_L / 2 - my_abs(y - Y_L / 2)))) * (1 - 6 / (9.6 - 3.528 * (1 - exp(-x / (L * 2))))));
	value = sqrt(value);
	if (z >= 60) 
		value = value * (60 - z) / 60;
	return sqrt(value) * 0.6;
}

// 空隙率
DEFINE_PROFILE(porous_profile, thread, position)
{
	real r[ND_ND];
	real x, y, z, value;
	cell_t c;
	printf("Load success!----%s\n", TIMESTAMP);
	begin_c_loop(c, thread)
	{
		C_CENTROID(r, c, thread);
		x = r[0];
		y = r[1];
		z = r[2];
		value = porous(x, y, z);
		C_PROFILE(c, thread, position) = value;
	}
	end_c_loop(c, thread)
}

// 瓦斯涌出量
DEFINE_PROFILE(gas, thread, position)
{
	real r[ND_ND];
	real x, y, z, value;
	cell_t c;
	begin_c_loop(c, thread)
	{
		C_CENTROID(r, c, thread);
		x = r[0];
		y = r[1];
		z = r[2];
		// if (x <= 200)
		// 	value = 5e-8 + 5e-10 * x;
		// else 
		// 	value = 1.5e-7;
		value = 1.5e-7;
		if (z >= 60) value *= 3;
		if (x < 100 && y > 200) 
			value *= pow(1 + (y - 200) / 60, 7);
		if (x <= 150)
			value *= 5;
		C_PROFILE(c, thread, position) = value;
	}
	end_c_loop(c, thread)
}

// 惯性阻力
DEFINE_PROFILE(guanxing, thread, position)
{
	real r[ND_ND];
	real x, y, z, value, n;
	cell_t c;
	begin_c_loop(c, thread)
	{
		C_CENTROID(r, c, thread);
		x = r[0];
		y = r[1];
		z = r[2];
		n = porous(x, y, z);
		value = 3.5 * (1 - n) / (DP * pow(n, 3));
		C_PROFILE(c, thread, position) = value;
	}
	end_c_loop(c, thread)
}

// 粘性阻力系数
DEFINE_PROFILE(nianxing_x, thread, position) 
{
	real r[ND_ND];
	real x, y, z, shentoulv, n;
	cell_t c;
	begin_c_loop(c, thread)
	{
		C_CENTROID(r, c, thread);
		x = r[0];
		y = r[1];
		z = r[2];
		n = porous(x, y, z);
		shentoulv = (DP * DP * pow(n, 3)) / (150 * pow(1 - n, 2));
		C_PROFILE(c, thread, position) = (y + 130) / (2 * shentoulv);
	}
	end_c_loop(c, thread)
}

DEFINE_PROFILE(nianxing_y, thread, position) 
{
	real r[ND_ND];
	real x, y, z, shentoulv, n;
	cell_t c;
	begin_c_loop(c, thread)
	{
		C_CENTROID(r, c, thread);
		x = r[0];
		y = r[1];
		z = r[2];
		n = porous(x, y, z);
		shentoulv = (DP * DP * pow(n, 3)) / (150 * pow(1 - n, 2));
		C_PROFILE(c, thread, position) = 15 / shentoulv;
	}
	end_c_loop(c, thread)
}

DEFINE_PROFILE(nianxing_z, thread, position) 
{
	real r[ND_ND];
	real x, y, z, shentoulv, n;
	cell_t c;
	begin_c_loop(c, thread)
	{
		C_CENTROID(r, c, thread);
		x = r[0];
		y = r[1];
		z = r[2];
		n = porous(x, y, z);
		shentoulv = (DP * DP * pow(n, 3)) / (150 * pow(1 - n, 2));
		C_PROFILE(c, thread, position) = 2000 / shentoulv;
	}
	end_c_loop(c, thread)
}