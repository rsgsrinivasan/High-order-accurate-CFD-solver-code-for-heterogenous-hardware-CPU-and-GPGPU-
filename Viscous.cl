#pragma OPENCL EXTENSION cl_khr_fp64 : enable
typedef struct __attribute__ ((packed)) tag_grid {
	double x, y, z;
	int N0, N1, N2, N3, N4, N5;
	int ID, loc;
} GRID;

typedef struct __attribute__ ((packed)) tag_flow_variables {
	double u, v, w, p, t, rho, mu, a, e;
} FLOW_VARIABLES;

typedef struct __attribute__ ((packed)) tag_Solver_Column_matrix {
	double N0, N1, N2, N3, N4;
} Solver_Column_matrix;

typedef struct __attribute__ ((packed)) tag_TRANSFORMED {
	double zeta_x, zeta_y, zeta_z, eta_x, eta_y, eta_z, xi_x, xi_y, xi_z;
	double zeta_xip, zeta_yip, zeta_zip, eta_xip, eta_yip, eta_zip, xi_xip, xi_yip, xi_zip, zeta_xim, zeta_yim, zeta_zim, eta_xim, eta_yim, eta_zim, xi_xim, xi_yim, xi_zim;
	double zeta_xjp, zeta_yjp, zeta_zjp, eta_xjp, eta_yjp, eta_zjp, xi_xjp, xi_yjp, xi_zjp, zeta_xjm, zeta_yjm, zeta_zjm, eta_xjm, eta_yjm, eta_zjm, xi_xjm, xi_yjm, xi_zjm;
	double zeta_xkp, zeta_ykp, zeta_zkp, eta_xkp, eta_ykp, eta_zkp, xi_xkp, xi_ykp, xi_zkp, zeta_xkm, zeta_ykm, zeta_zkm, eta_xkm, eta_ykm, eta_zkm, xi_xkm, xi_ykm, xi_zkm;
} TRANSFORMED;

__kernel void Viscous(__global Solver_Column_matrix *U_C, __global Solver_Column_matrix *F_C, \
__global Solver_Column_matrix *E_C, __global Solver_Column_matrix *G_C, __global Solver_Column_matrix *Fv_C, \
__global Solver_Column_matrix *Ev_C, __global Solver_Column_matrix *Gv_C, __global const GRID *Grid_P, \
__global const TRANSFORMED *metric, __global const double *det, __global FLOW_VARIABLES *FLOW)
{
    uint Gid = get_global_id(0);

	__private const double Mach = 1.5;
	__private const double Reyl = 500.0;
	__private const double Pr = 0.72;
	__private const double Free_t = 273.0;

    __private int IP1, IP2, IP3, IM1, IM2, IM3;
    __private int JP1, JP2, JP3, JM1, JM2, JM3;
    __private int KP1, KP2, KP3, KM1, KM2, KM3;
    
    __private int LOC_I;
    __private int LOC_IP1, LOC_IP2, LOC_IP3, LOC_IM1, LOC_IM2, LOC_IM3;
    __private int LOC_JP1, LOC_JP2, LOC_JP3, LOC_JM1, LOC_JM2, LOC_JM3;
    __private int LOC_KP1, LOC_KP2, LOC_KP3, LOC_KM1, LOC_KM2, LOC_KM3;

    __private double U0_I;
    __private double U0_IP1, U0_IP2, U0_IP3, U0_IM1, U0_IM2, U0_IM3;
    __private double U0_JP1, U0_JP2, U0_JP3, U0_JM1, U0_JM2, U0_JM3;
    __private double U0_KP1, U0_KP2, U0_KP3, U0_KM1, U0_KM2, U0_KM3;

    __private double U1_I;
    __private double U1_IP1, U1_IP2, U1_IP3, U1_IM1, U1_IM2, U1_IM3;
    __private double U1_JP1, U1_JP2, U1_JP3, U1_JM1, U1_JM2, U1_JM3;
    __private double U1_KP1, U1_KP2, U1_KP3, U1_KM1, U1_KM2, U1_KM3;
    
    __private double U2_I;
    __private double U2_IP1, U2_IP2, U2_IP3, U2_IM1, U2_IM2, U2_IM3;
    __private double U2_JP1, U2_JP2, U2_JP3, U2_JM1, U2_JM2, U2_JM3;
    __private double U2_KP1, U2_KP2, U2_KP3, U2_KM1, U2_KM2, U2_KM3;

    __private double U3_I;
    __private double U3_IP1, U3_IP2, U3_IP3, U3_IM1, U3_IM2, U3_IM3;
    __private double U3_JP1, U3_JP2, U3_JP3, U3_JM1, U3_JM2, U3_JM3;
    __private double U3_KP1, U3_KP2, U3_KP3, U3_KM1, U3_KM2, U3_KM3;

    __private double U4_I;
    __private double U4_IP1, U4_IP2, U4_IP3, U4_IM1, U4_IM2, U4_IM3;
    __private double U4_JP1, U4_JP2, U4_JP3, U4_JM1, U4_JM2, U4_JM3;
    __private double U4_KP1, U4_KP2, U4_KP3, U4_KM1, U4_KM2, U4_KM3;

    __private double FLOW_rho, FLOW_u, FLOW_mu, FLOW_v, FLOW_w, FLOW_p, FLOW_t, FLOW_a, FLOW_e;

    __private double DET_I;
    __private double zeta_x, zeta_y, zeta_z;
    __private double eta_x, eta_y, eta_z;
    __private double xi_x, xi_y, xi_z;

    __private double t_zeta, t_eta, t_xi, u_zeta, u_eta, u_xi;
	__private double v_zeta, v_eta, v_xi, w_zeta, w_eta, w_xi;

	__private double U1_zeta, U2_zeta, U3_zeta, U4_zeta, U5_zeta;
    __private double U1_eta, U2_eta, U3_eta, U4_eta, U5_eta;
    __private double U1_xi, U2_xi, U3_xi, U4_xi, U5_xi;

	__private double d_temp, d_temp1;
	__private double stress_xy, stress_yz, stress_xz, stress_xx, stress_yy, stress_zz;
	__private double KS_x, KS_y, KS_z;
	__private double EDDY_x, EDDY_y, EDDY_z, mod_stress;
	__private double DEL_ZETA, DEL_ETA, DEL_XI;
	__private double Q_x, Q_y, Q_z, DEL_X, DEL_Y, DEL_Z;
	__private double C_MTS, C_T, Rx, Ry, Rz;
	__private double stress_xy_star, stress_yz_star, stress_xz_star;
	__private double divergence_V, cross_V;

    __private double tauzz, tauee, tauxx, tauze, tauzx, tauex, qz, qe, qx;
    __private double TAU_SGS_XY, TAU_SGS_YZ, TAU_SGS_XZ, TAU_SGS_XX, TAU_SGS_YY, TAU_SGS_ZZ;
    __private double H_SGS_X, H_SGS_Y, H_SGS_Z, D_SGS_X, D_SGS_Y, D_SGS_Z;

	__private Solver_Column_matrix tF1_C, tE1_C, tG1_C, tFv1_C, tEv1_C, tGv1_C; 

    IP1 = Grid_P[Gid].N1;
    IP2 = Grid_P[Grid_P[Gid].N1].N1;
    IP3 = Grid_P[Grid_P[Grid_P[Gid].N1].N1].N1;
    IM1 = Grid_P[Gid].N3;
    IM2 = Grid_P[Grid_P[Gid].N3].N3;
    IM3 = Grid_P[Grid_P[Grid_P[Gid].N3].N3].N3;
    
    JP1 = Grid_P[Gid].N0;
    JP2 = Grid_P[Grid_P[Gid].N0].N0;
    JP3 = Grid_P[Grid_P[Grid_P[Gid].N0].N0].N0;
    JM1 = Grid_P[Gid].N2;
    JM2 = Grid_P[Grid_P[Gid].N2].N2;
    JM3 = Grid_P[Grid_P[Grid_P[Gid].N2].N2].N2;
    
    KP1 = Grid_P[Gid].N4;
    KP2 = Grid_P[Grid_P[Gid].N4].N4;
    KP3 = Grid_P[Grid_P[Grid_P[Gid].N4].N4].N4;
    KM1 = Grid_P[Gid].N5;
    KM2 = Grid_P[Grid_P[Gid].N5].N5;
    KM3 = Grid_P[Grid_P[Grid_P[Gid].N5].N5].N5;

    DET_I = det[Gid];
    zeta_x = metric[Gid].zeta_x;
    zeta_y = metric[Gid].zeta_y;
    zeta_z = metric[Gid].zeta_z;
    eta_x = metric[Gid].eta_x;
    eta_y = metric[Gid].eta_y;
    eta_z = metric[Gid].eta_z;
    xi_x = metric[Gid].xi_x;
    xi_y = metric[Gid].xi_y;
    xi_z = metric[Gid].xi_z;

    FLOW_u = FLOW[Gid].u;
    FLOW_v = FLOW[Gid].v;
    FLOW_w = FLOW[Gid].w;
    FLOW_p = FLOW[Gid].p;
    FLOW_t = FLOW[Gid].t;
    FLOW_a = FLOW[Gid].a;
    FLOW_e = FLOW[Gid].e;
    FLOW_mu = FLOW[Gid].mu;
    FLOW_rho = FLOW[Gid].rho;

    //barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

    LOC_I = Grid_P[Gid].loc;
    LOC_IP1 = Grid_P[Grid_P[Gid].N1].loc;
    LOC_IP2 = Grid_P[Grid_P[Grid_P[Gid].N1].N1].loc;
    LOC_IP3 = Grid_P[Grid_P[Grid_P[Grid_P[Gid].N1].N1].N1].loc;
    LOC_IM1 = Grid_P[Grid_P[Gid].N3].loc;
    LOC_IM2 = Grid_P[Grid_P[Grid_P[Gid].N3].N3].loc;
    LOC_IM3 = Grid_P[Grid_P[Grid_P[Grid_P[Gid].N3].N3].N3].loc;
    
    LOC_JP1 = Grid_P[Grid_P[Gid].N0].loc;
    LOC_JP2 = Grid_P[Grid_P[Grid_P[Gid].N0].N0].loc;
    LOC_JP3 = Grid_P[Grid_P[Grid_P[Grid_P[Gid].N0].N0].N0].loc;
    LOC_JM1 = Grid_P[Grid_P[Gid].N2].loc;
    LOC_JM2 = Grid_P[Grid_P[Grid_P[Gid].N2].N2].loc;
    LOC_JM3 = Grid_P[Grid_P[Grid_P[Grid_P[Gid].N2].N2].N2].loc;
    
    LOC_KP1 = Grid_P[Grid_P[Gid].N4].loc;
    LOC_KP2 = Grid_P[Grid_P[Grid_P[Gid].N4].N4].loc;
    LOC_KP3 = Grid_P[Grid_P[Grid_P[Grid_P[Gid].N4].N4].N4].loc;
    LOC_KM1 = Grid_P[Grid_P[Gid].N5].loc;
    LOC_KM2 = Grid_P[Grid_P[Grid_P[Gid].N5].N5].loc;
    LOC_KM3 = Grid_P[Grid_P[Grid_P[Grid_P[Gid].N5].N5].N5].loc;

    //barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

    U0_I = U_C[Gid].N0;
    U1_I = U_C[Gid].N1;
    U2_I = U_C[Gid].N2;
    U3_I = U_C[Gid].N3;
    U4_I = U_C[Gid].N4;

    U0_IP1 = U_C[IP1].N0;
    U1_IP1 = U_C[IP1].N1;
    U2_IP1 = U_C[IP1].N2;
    U3_IP1 = U_C[IP1].N3;
    U4_IP1 = U_C[IP1].N4;

    U0_IP2 = U_C[IP2].N0;
    U1_IP2 = U_C[IP2].N1;
    U2_IP2 = U_C[IP2].N2;
    U3_IP2 = U_C[IP2].N3;
    U4_IP2 = U_C[IP2].N4;

    U0_IP3 = U_C[IP3].N0;
    U1_IP3 = U_C[IP3].N1;
    U2_IP3 = U_C[IP3].N2;
    U3_IP3 = U_C[IP3].N3;
    U4_IP3 = U_C[IP3].N4;

    U0_JP1 = U_C[JP1].N0;
    U1_JP1 = U_C[JP1].N1;
    U2_JP1 = U_C[JP1].N2;
    U3_JP1 = U_C[JP1].N3;
    U4_JP1 = U_C[JP1].N4;

    U0_JP2 = U_C[JP2].N0;
    U1_JP2 = U_C[JP2].N1;
    U2_JP2 = U_C[JP2].N2;
    U3_JP2 = U_C[JP2].N3;
    U4_JP2 = U_C[JP2].N4;

    U0_JP3 = U_C[JP3].N0;
    U1_JP3 = U_C[JP3].N1;
    U2_JP3 = U_C[JP3].N2;
    U3_JP3 = U_C[JP3].N3;
    U4_JP3 = U_C[JP3].N4;

    U0_KP1 = U_C[KP1].N0;
    U1_KP1 = U_C[KP1].N1;
    U2_KP1 = U_C[KP1].N2;
    U3_KP1 = U_C[KP1].N3;
    U4_KP1 = U_C[KP1].N4;

    U0_KP2 = U_C[KP2].N0;
    U1_KP2 = U_C[KP2].N1;
    U2_KP2 = U_C[KP2].N2;
    U3_KP2 = U_C[KP2].N3;
    U4_KP2 = U_C[KP2].N4;

    U0_KP3 = U_C[KP3].N0;
    U1_KP3 = U_C[KP3].N1;
    U2_KP3 = U_C[KP3].N2;
    U3_KP3 = U_C[KP3].N3;
    U4_KP3 = U_C[KP3].N4;

	U0_IM1 = U_C[IM1].N0;
    U1_IM1 = U_C[IM1].N1;
    U2_IM1 = U_C[IM1].N2;
    U3_IM1 = U_C[IM1].N3;
    U4_IM1 = U_C[IM1].N4;

    U0_IM2 = U_C[IM2].N0;
    U1_IM2 = U_C[IM2].N1;
    U2_IM2 = U_C[IM2].N2;
    U3_IM2 = U_C[IM2].N3;
    U4_IM2 = U_C[IM2].N4;

    U0_IM3 = U_C[IM3].N0;
    U1_IM3 = U_C[IM3].N1;
    U2_IM3 = U_C[IM3].N2;
    U3_IM3 = U_C[IM3].N3;
    U4_IM3 = U_C[IM3].N4;

    U0_JM1 = U_C[JM1].N0;
    U1_JM1 = U_C[JM1].N1;
    U2_JM1 = U_C[JM1].N2;
    U3_JM1 = U_C[JM1].N3;
    U4_JM1 = U_C[JM1].N4;

    U0_JM2 = U_C[JM2].N0;
    U1_JM2 = U_C[JM2].N1;
    U2_JM2 = U_C[JM2].N2;
    U3_JM2 = U_C[JM2].N3;
    U4_JM2 = U_C[JM2].N4;

    U0_JM3 = U_C[JM3].N0;
    U1_JM3 = U_C[JM3].N1;
    U2_JM3 = U_C[JM3].N2;
    U3_JM3 = U_C[JM3].N3;
    U4_JM3 = U_C[JM3].N4;

    U0_KM1 = U_C[KM1].N0;
    U1_KM1 = U_C[KM1].N1;
    U2_KM1 = U_C[KM1].N2;
    U3_KM1 = U_C[KM1].N3;
    U4_KM1 = U_C[KM1].N4;

    U0_KM2 = U_C[KM2].N0;
    U1_KM2 = U_C[KM2].N1;
    U2_KM2 = U_C[KM2].N2;
    U3_KM2 = U_C[KM2].N3;
    U4_KM2 = U_C[KM2].N4;

    U0_KM3 = U_C[KM3].N0;
    U1_KM3 = U_C[KM3].N1;
    U2_KM3 = U_C[KM3].N2;
    U3_KM3 = U_C[KM3].N3;
    U4_KM3 = U_C[KM3].N4;



    barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

	if (LOC_I > 0 && LOC_I <= 52 && LOC_IP1 == 100)
	{
		U1_zeta = (1.0 / 6.0)*(11.0*U0_I - 18.0*U0_IM1 + 9.0*U0_IM2 - 2.0*U0_IM3);
		U2_zeta = (1.0 / 6.0)*(11.0*U1_I - 18.0*U1_IM1 + 9.0*U1_IM2 - 2.0*U1_IM3);
		U3_zeta = (1.0 / 6.0)*(11.0*U2_I - 18.0*U2_IM1 + 9.0*U2_IM2 - 2.0*U2_IM3);
		U4_zeta = (1.0 / 6.0)*(11.0*U3_I - 18.0*U3_IM1 + 9.0*U3_IM2 - 2.0*U3_IM3);
		U5_zeta = (1.0 / 6.0)*(11.0*U4_I - 18.0*U4_IM1 + 9.0*U4_IM2 - 2.0*U4_IM3);
	}
	if (LOC_I > 0 && LOC_I <= 52 && LOC_IM1 == 100)
	{
		U1_zeta = (1.0 / 6.0)*(-11.0*U0_I + 18.0*U0_IP1 - 9.0*U0_IP2 + 2.0*U0_IP3);
		U2_zeta = (1.0 / 6.0)*(-11.0*U1_I + 18.0*U1_IP1 - 9.0*U1_IP2 + 2.0*U1_IP3);
		U3_zeta = (1.0 / 6.0)*(-11.0*U2_I + 18.0*U2_IP1 - 9.0*U2_IP2 + 2.0*U2_IP3);
		U4_zeta = (1.0 / 6.0)*(-11.0*U3_I + 18.0*U3_IP1 - 9.0*U3_IP2 + 2.0*U3_IP3);
		U5_zeta = (1.0 / 6.0)*(-11.0*U4_I + 18.0*U4_IP1 - 9.0*U4_IP2 + 2.0*U4_IP3);
	}
	if (LOC_I == 0 && LOC_IP1 > 0 && LOC_IP1 <= 52 && LOC_IP2 == 100)
	{
		U1_zeta = (1.0 / 6.0)*(2.0*U0_IP1 + 3.0*U0_I - 6.0*U0_IM1 + U0_IM2);
		U2_zeta = (1.0 / 6.0)*(2.0*U1_IP1 + 3.0*U1_I - 6.0*U1_IM1 + U1_IM2);
		U3_zeta = (1.0 / 6.0)*(2.0*U2_IP1 + 3.0*U2_I - 6.0*U2_IM1 + U2_IM2);
		U4_zeta = (1.0 / 6.0)*(2.0*U3_IP1 + 3.0*U3_I - 6.0*U3_IM1 + U3_IM2);
		U5_zeta = (1.0 / 6.0)*(2.0*U4_IP1 + 3.0*U4_I - 6.0*U4_IM1 + U4_IM2);
	}
	if (LOC_I == 0 && LOC_IM1 > 0 && LOC_IM1 && LOC_IM2 == 100)
	{
		U1_zeta = (1.0 / 6.0)*((-2.0)*U0_IM1 - 3.0*U0_I + 6.0*U0_IP1 - U0_IP2);
		U2_zeta = (1.0 / 6.0)*((-2.0)*U1_IM1 - 3.0*U1_I + 6.0*U1_IP1 - U1_IP2);
		U3_zeta = (1.0 / 6.0)*((-2.0)*U2_IM1 - 3.0*U2_I + 6.0*U2_IP1 - U2_IP2);
		U4_zeta = (1.0 / 6.0)*((-2.0)*U3_IM1 - 3.0*U3_I + 6.0*U3_IP1 - U3_IP2);
		U5_zeta = (1.0 / 6.0)*((-2.0)*U4_IM1 - 3.0*U4_I + 6.0*U4_IP1 - U4_IP2);
	}
	if (LOC_I == 0 && LOC_IM2 != 100 && LOC_IM1 != 100 && LOC_IP2 != 100 && LOC_IP1 != 100)
	{
		U1_zeta = (1.0 / 12.0)*(8.0*(U0_IP1 - U0_IM1) - (U0_IP2 - U0_IM2));
		U2_zeta = (1.0 / 12.0)*(8.0*(U1_IP1 - U1_IM1) - (U1_IP2 - U1_IM2));
		U3_zeta = (1.0 / 12.0)*(8.0*(U2_IP1 - U2_IM1) - (U2_IP2 - U2_IM2));
		U4_zeta = (1.0 / 12.0)*(8.0*(U3_IP1 - U3_IM1) - (U3_IP2 - U3_IM2));
		U5_zeta = (1.0 / 12.0)*(8.0*(U4_IP1 - U4_IM1) - (U4_IP2 - U4_IM2));
	}
	if (LOC_I > 0 && LOC_IM2 != 100 && LOC_IM1 != 100 && LOC_IP2 != 100 && LOC_IP1 != 100)
	{
		U1_zeta = (1.0 / 12.0)*(8.0*(U0_IP1 - U0_IM1) - (U0_IP2 - U0_IM2));
		U2_zeta = (1.0 / 12.0)*(8.0*(U1_IP1 - U1_IM1) - (U1_IP2 - U1_IM2));
		U3_zeta = (1.0 / 12.0)*(8.0*(U2_IP1 - U2_IM1) - (U2_IP2 - U2_IM2));
		U4_zeta = (1.0 / 12.0)*(8.0*(U3_IP1 - U3_IM1) - (U3_IP2 - U3_IM2));
		U5_zeta = (1.0 / 12.0)*(8.0*(U4_IP1 - U4_IM1) - (U4_IP2 - U4_IM2));
	}
	





	if (LOC_I > 0 && LOC_I <= 52 && LOC_JP1 == 100)
	{
		U1_eta = (1.0 / 6.0)*(11.0*U0_I - 18.0*U0_JM1 + 9.0*U0_JM2 - 2.0*U0_JM3);
		U2_eta = (1.0 / 6.0)*(11.0*U1_I - 18.0*U1_JM1 + 9.0*U1_JM2 - 2.0*U1_JM3);
		U3_eta = (1.0 / 6.0)*(11.0*U2_I - 18.0*U2_JM1 + 9.0*U2_JM2 - 2.0*U2_JM3);
		U4_eta = (1.0 / 6.0)*(11.0*U3_I - 18.0*U3_JM1 + 9.0*U3_JM2 - 2.0*U3_JM3);
		U5_eta = (1.0 / 6.0)*(11.0*U4_I - 18.0*U4_JM1 + 9.0*U4_JM2 - 2.0*U4_JM3);
	}
	if (LOC_I > 0 && LOC_I <= 52 && LOC_JM1 == 100)
	{
		U1_eta = (1.0 / 6.0)*(-11.0*U0_I + 18.0*U0_JP1 - 9.0*U0_JP2 + 2.0*U0_JP3);
		U2_eta = (1.0 / 6.0)*(-11.0*U1_I + 18.0*U1_JP1 - 9.0*U1_JP2 + 2.0*U1_JP3);
		U3_eta = (1.0 / 6.0)*(-11.0*U2_I + 18.0*U2_JP1 - 9.0*U2_JP2 + 2.0*U2_JP3);
		U4_eta = (1.0 / 6.0)*(-11.0*U3_I + 18.0*U3_JP1 - 9.0*U3_JP2 + 2.0*U3_JP3);
		U5_eta = (1.0 / 6.0)*(-11.0*U4_I + 18.0*U4_JP1 - 9.0*U4_JP2 + 2.0*U4_JP3);
	}
	if (LOC_I == 0 && LOC_JP1 > 0 && LOC_JP1 <= 52 && LOC_JP2 == 100)
	{
		U1_eta = (1.0 / 6.0)*(2.0*U0_JP1 + 3.0*U0_I - 6.0*U0_JM1 + U0_JM2);
		U2_eta = (1.0 / 6.0)*(2.0*U1_JP1 + 3.0*U1_I - 6.0*U1_JM1 + U1_JM2);
		U3_eta = (1.0 / 6.0)*(2.0*U2_JP1 + 3.0*U2_I - 6.0*U2_JM1 + U2_JM2);
		U4_eta = (1.0 / 6.0)*(2.0*U3_JP1 + 3.0*U3_I - 6.0*U3_JM1 + U3_JM2);
		U5_eta = (1.0 / 6.0)*(2.0*U4_JP1 + 3.0*U4_I - 6.0*U4_JM1 + U4_JM2);
	}
	if (LOC_I == 0 && LOC_JM1 > 0 && LOC_JM1 && LOC_JM2 == 100)
	{
		U1_eta = (1.0 / 6.0)*((-2.0)*U0_JM1 - 3.0*U0_I + 6.0*U0_JP1 - U0_JP2);
		U2_eta = (1.0 / 6.0)*((-2.0)*U1_JM1 - 3.0*U1_I + 6.0*U1_JP1 - U1_JP2);
		U3_eta = (1.0 / 6.0)*((-2.0)*U2_JM1 - 3.0*U2_I + 6.0*U2_JP1 - U2_JP2);
		U4_eta = (1.0 / 6.0)*((-2.0)*U3_JM1 - 3.0*U3_I + 6.0*U3_JP1 - U3_JP2);
		U5_eta = (1.0 / 6.0)*((-2.0)*U4_JM1 - 3.0*U4_I + 6.0*U4_JP1 - U4_JP2);
	}
	if (LOC_I == 0 && LOC_JM2 != 100 && LOC_JM1 != 100 && LOC_JP2 != 100 && LOC_JP1 != 100)
	{
		U1_eta = (1.0 / 12.0)*(8.0*(U0_JP1 - U0_JM1) - (U0_JP2 - U0_JM2));
		U2_eta = (1.0 / 12.0)*(8.0*(U1_JP1 - U1_JM1) - (U1_JP2 - U1_JM2));
		U3_eta = (1.0 / 12.0)*(8.0*(U2_JP1 - U2_JM1) - (U2_JP2 - U2_JM2));
		U4_eta = (1.0 / 12.0)*(8.0*(U3_JP1 - U3_JM1) - (U3_JP2 - U3_JM2));
		U5_eta = (1.0 / 12.0)*(8.0*(U4_JP1 - U4_JM1) - (U4_JP2 - U4_JM2));
	}
	if (LOC_I > 0 && LOC_JM2 != 100 && LOC_JM1 != 100 && LOC_JP2 != 100 && LOC_JP1 != 100)
	{
		U1_eta = (1.0 / 12.0)*(8.0*(U0_JP1 - U0_JM1) - (U0_JP2 - U0_JM2));
		U2_eta = (1.0 / 12.0)*(8.0*(U1_JP1 - U1_JM1) - (U1_JP2 - U1_JM2));
		U3_eta = (1.0 / 12.0)*(8.0*(U2_JP1 - U2_JM1) - (U2_JP2 - U2_JM2));
		U4_eta = (1.0 / 12.0)*(8.0*(U3_JP1 - U3_JM1) - (U3_JP2 - U3_JM2));
		U5_eta = (1.0 / 12.0)*(8.0*(U4_JP1 - U4_JM1) - (U4_JP2 - U4_JM2));
	}

	if (LOC_I > 0 && LOC_I <= 52 && LOC_KP1 == 100)
	{
		U1_xi = (1.0 / 6.0)*(11.0*U0_I - 18.0*U0_KM1 + 9.0*U0_KM2 - 2.0*U0_KM3);
		U2_xi = (1.0 / 6.0)*(11.0*U1_I - 18.0*U1_KM1 + 9.0*U1_KM2 - 2.0*U1_KM3);
		U3_xi = (1.0 / 6.0)*(11.0*U2_I - 18.0*U2_KM1 + 9.0*U2_KM2 - 2.0*U2_KM3);
		U4_xi = (1.0 / 6.0)*(11.0*U3_I - 18.0*U3_KM1 + 9.0*U3_KM2 - 2.0*U3_KM3);
		U5_xi = (1.0 / 6.0)*(11.0*U4_I - 18.0*U4_KM1 + 9.0*U4_KM2 - 2.0*U4_KM3);
	}
	if (LOC_I > 0 && LOC_I <= 52 && LOC_KM1 == 100)
	{
		U1_xi = (1.0 / 6.0)*(-11.0*U0_I + 18.0*U0_KP1 - 9.0*U0_KP2 + 2.0*U0_KP3);
		U2_xi = (1.0 / 6.0)*(-11.0*U1_I + 18.0*U1_KP1 - 9.0*U1_KP2 + 2.0*U1_KP3);
		U3_xi = (1.0 / 6.0)*(-11.0*U2_I + 18.0*U2_KP1 - 9.0*U2_KP2 + 2.0*U2_KP3);
		U4_xi = (1.0 / 6.0)*(-11.0*U3_I + 18.0*U3_KP1 - 9.0*U3_KP2 + 2.0*U3_KP3);
		U5_xi = (1.0 / 6.0)*(-11.0*U4_I + 18.0*U4_KP1 - 9.0*U4_KP2 + 2.0*U4_KP3);
	}
	if (LOC_I == 0 && LOC_KP1 > 0 && LOC_KP1 <= 52 && LOC_KP2 == 100)
	{
		U1_xi = (1.0 / 6.0)*(2.0*U0_KP1 + 3.0*U0_I - 6.0*U0_KM1 + U0_KM2);
		U2_xi = (1.0 / 6.0)*(2.0*U1_KP1 + 3.0*U1_I - 6.0*U1_KM1 + U1_KM2);
		U3_xi = (1.0 / 6.0)*(2.0*U2_KP1 + 3.0*U2_I - 6.0*U2_KM1 + U2_KM2);
		U4_xi = (1.0 / 6.0)*(2.0*U3_KP1 + 3.0*U3_I - 6.0*U3_KM1 + U3_KM2);
		U5_xi = (1.0 / 6.0)*(2.0*U4_KP1 + 3.0*U4_I - 6.0*U4_KM1 + U4_KM2);
	}
	if (LOC_I == 0 && LOC_KM1 > 0 && LOC_KM1 && LOC_KM2 == 100)
	{
		U1_xi = (1.0 / 6.0)*((-2.0)*U0_KM1 - 3.0*U0_I + 6.0*U0_KP1 - U0_KP2);
		U2_xi = (1.0 / 6.0)*((-2.0)*U1_KM1 - 3.0*U1_I + 6.0*U1_KP1 - U1_KP2);
		U3_xi = (1.0 / 6.0)*((-2.0)*U2_KM1 - 3.0*U2_I + 6.0*U2_KP1 - U2_KP2);
		U4_xi = (1.0 / 6.0)*((-2.0)*U3_KM1 - 3.0*U3_I + 6.0*U3_KP1 - U3_KP2);
		U5_xi = (1.0 / 6.0)*((-2.0)*U4_KM1 - 3.0*U4_I + 6.0*U4_KP1 - U4_KP2);
	}
	if (LOC_I == 0 && LOC_KM2 != 100 && LOC_KM1 != 100 && LOC_KP2 != 100 && LOC_KP1 != 100)
	{
		U1_xi = (1.0 / 12.0)*(8.0*(U0_KP1 - U0_KM1) - (U0_KP2 - U0_KM2));
		U2_xi = (1.0 / 12.0)*(8.0*(U1_KP1 - U1_KM1) - (U1_KP2 - U1_KM2));
		U3_xi = (1.0 / 12.0)*(8.0*(U2_KP1 - U2_KM1) - (U2_KP2 - U2_KM2));
		U4_xi = (1.0 / 12.0)*(8.0*(U3_KP1 - U3_KM1) - (U3_KP2 - U3_KM2));
		U5_xi = (1.0 / 12.0)*(8.0*(U4_KP1 - U4_KM1) - (U4_KP2 - U4_KM2));
	}
	if (LOC_I > 0 && LOC_KM2 != 100 && LOC_KM1 != 100 && LOC_KP2 != 100 && LOC_KP1 != 100)
	{
		U1_xi = (1.0 / 12.0)*(8.0*(U0_KP1 - U0_KM1) - (U0_KP2 - U0_KM2));
		U2_xi = (1.0 / 12.0)*(8.0*(U1_KP1 - U1_KM1) - (U1_KP2 - U1_KM2));
		U3_xi = (1.0 / 12.0)*(8.0*(U2_KP1 - U2_KM1) - (U2_KP2 - U2_KM2));
		U4_xi = (1.0 / 12.0)*(8.0*(U3_KP1 - U3_KM1) - (U3_KP2 - U3_KM2));
		U5_xi = (1.0 / 12.0)*(8.0*(U4_KP1 - U4_KM1) - (U4_KP2 - U4_KM2));
	}
	
   // barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

	u_zeta = (1.0 / U0_I)*(U2_zeta)-(U1_I / ((U0_I * U0_I)))*(U1_zeta);
	u_eta = (1.0 / U0_I)*(U2_eta)-(U1_I / ((U0_I * U0_I)))*(U1_eta);
	u_xi = (1.0 / U0_I)*(U2_xi)-(U1_I / ((U0_I * U0_I)))*(U1_xi);

	v_zeta = (1.0 / U0_I)*(U3_zeta)-(U2_I / ((U0_I * U0_I)))*(U1_zeta);
	v_eta = (1.0 / U0_I)*(U3_eta)-(U2_I / ((U0_I * U0_I)))*(U1_eta);
	v_xi = (1.0 / U0_I)*(U3_xi)-(U2_I / ((U0_I * U0_I)))*(U1_xi);

	w_zeta = (1.0 / U0_I)*(U4_zeta)-(U3_I / ((U0_I * U0_I)))*(U1_zeta);
	w_eta = (1.0 / U0_I)*(U4_eta)-(U3_I / ((U0_I * U0_I)))*(U1_eta);
	w_xi = (1.0 / U0_I)*(U4_xi)-(U3_I / ((U0_I * U0_I)))*(U1_xi);

	t_zeta = (1.0 / U0_I)*(U5_zeta)-(U4_I / ((U0_I * U0_I)))*(U1_zeta)-(U1_I / U0_I)*u_zeta - (U2_I / U0_I)*v_zeta - (U3_I / U0_I)*w_zeta;
	t_eta = (1.0 / U0_I)*(U5_eta)-(U4_I / ((U0_I * U0_I)))*(U1_eta)-(U1_I / U0_I)*u_eta - (U2_I / U0_I)*v_eta - (U3_I / U0_I)*w_eta;
	t_xi = (1.0 / U0_I)*(U5_xi)-(U4_I / ((U0_I * U0_I)))*(U1_xi)-(U1_I / U0_I)*u_xi - (U2_I / U0_I)*v_xi - (U3_I / U0_I)*w_xi;

	d_temp = (1.4*Mach*Mach*0.4*((U4_I / U0_I) - 0.5*(((U1_I * U1_I) / (U0_I * U0_I)) + ((U2_I * U2_I) / (U0_I * U0_I)) + ((U3_I * U3_I) / (U0_I * U0_I)))));
	d_temp1 = sqrt(d_temp)*sqrt(d_temp)*sqrt(d_temp);
	/**************************************SUTHERLAND LAW***************************************************************************/
	FLOW_mu = (d_temp1)*((1.0 + (110.4 / Free_t)) / (d_temp + (110.4 / Free_t)));
	/********************************************************************************************************************************/
	/*****************************************HEAT FLUX******************************************************************************/
	qz = (-1.0)*(FLOW_mu * 1.4 / (Pr*Reyl))*(t_zeta*zeta_x + t_eta * eta_x + t_xi * xi_x);
	qe = (-1.0)*(FLOW_mu * 1.4 / (Pr*Reyl))*(t_zeta*zeta_y + t_eta * eta_y + t_xi * xi_y);
	qx = (-1.0)*(FLOW_mu * 1.4 / (Pr*Reyl))*(t_zeta*zeta_z + t_eta * eta_z + t_xi * xi_z);
	/*********************************************************************************************************************************/
	/***************************************************viscous shear stess terms******************************************************/
	tauzz = (FLOW_mu / Reyl)*((4.0 / 3.0)*(u_zeta*zeta_x + u_eta * eta_x + u_xi * xi_x) - (2.0 / 3.0)*(v_zeta*zeta_y + v_eta * eta_y + v_xi * xi_y) - (2.0 / 3.0)*(zeta_z*w_zeta + eta_z*w_eta + xi_z*w_xi));

	tauee = (FLOW_mu / Reyl)*((4.0 / 3.0)*(v_zeta*zeta_y + v_eta * eta_y + v_xi * xi_y) - (2.0 / 3.0)*(u_zeta*zeta_x + u_eta * eta_x + u_xi * xi_x) - (2.0 / 3.0)*(zeta_z*w_zeta + eta_z*w_eta + xi_z*w_xi));

	tauxx = (FLOW_mu / Reyl)*((4.0 / 3.0)*(w_zeta*zeta_z + w_eta * eta_z + w_xi * xi_z) - (2.0 / 3.0)*(u_zeta*zeta_x + u_eta * eta_x + u_xi * xi_x) - (2.0 / 3.0)*(zeta_y*v_zeta + eta_y*v_eta + xi_y*v_xi));

	tauze = (FLOW_mu / Reyl)*(u_zeta*zeta_y + u_eta * eta_y + u_xi * xi_y + v_zeta * zeta_x + v_eta * eta_x + v_xi * xi_x);

	tauzx = (FLOW_mu / Reyl)*(zeta_z*u_zeta + eta_z*u_eta + xi_z*u_xi + zeta_x*w_zeta + eta_x*w_eta + xi_x*w_xi);

	tauex = (FLOW_mu / Reyl)*(zeta_z*v_zeta + eta_z*v_eta + xi_z*v_xi + zeta_y*w_zeta + eta_y*w_eta + xi_y*w_xi);

	/*******************************************DUCROS SENSOR*********************************************************************/
	divergence_V = u_zeta * zeta_x + u_eta * eta_x + u_xi * xi_x + v_zeta * zeta_y + v_eta * eta_y + v_xi * xi_y + w_zeta * zeta_z + w_eta * eta_z + w_xi * xi_z;
	cross_V = pow(v_zeta*zeta_z + v_eta * eta_z + v_xi * xi_z - (w_zeta*zeta_y + w_eta * eta_y + w_xi * xi_y), 2.0) + \
		pow(u_zeta*zeta_z + u_eta * eta_z + u_xi * xi_z - (w_zeta*zeta_x + w_eta * eta_x + w_xi * xi_x), 2.0) + \
		pow(u_zeta*zeta_y + u_eta * eta_y + u_xi * xi_y - (v_zeta*zeta_x + v_eta * eta_x + v_xi * xi_x), 2.0);

	// DUCROS[Gid] = (divergence_V*divergence_V) / ((divergence_V*divergence_V) + (cross_V)+10e-30);

	/****************************************************************************************************************************/
	/****************************************LARGE EDDY SIMULATION TERMS*********************************************************/
	/***********************************MIXED TIME SCALE EDDY VISCOSITY TERMS****************************************************/
	/*******************************************SGS MOMENTUM TERMS***************************************************************/
	C_MTS = 0.03;
	C_T = 10.0;

	KS_x = pow(((U1_I / U0_I) - ((0.25*U1_IM1 + 0.5*U1_I + 0.25*U1_IP1) / (0.25*U0_IM1 + 0.5*U0_I + 0.25*U0_IP1))), 2.0);
	KS_y = pow(((U2_I / U0_I) - ((0.25*U2_JM1 + 0.5*U2_I + 0.25*U2_JP1) / (0.25*U0_JM1 + 0.5*U0_I + 0.25*U0_JP1))), 2.0);
	KS_z = pow(((U3_I / U0_I) - ((0.25*U3_KM1 + 0.5*U3_I + 0.25*U3_KP1) / (0.25*U0_KM1 + 0.5*U0_I + 0.25*U0_KP1))), 2.0);

	stress_xy_star = (0.5*(u_zeta*zeta_y + u_eta * eta_y + u_xi * xi_y + v_zeta * zeta_x + v_eta * eta_x + v_xi * xi_x));
	stress_yz_star = (0.5*(v_zeta*zeta_z + v_eta * eta_z + v_xi * xi_z + w_zeta * zeta_y + w_eta * eta_y + w_xi * xi_y));
	stress_xz_star = (0.5*(w_zeta*zeta_x + w_eta * eta_x + w_xi * xi_x + u_zeta * zeta_z + u_eta * eta_z + u_xi * xi_z));

	stress_xy = (0.5*(u_zeta*zeta_y + u_eta * eta_y + u_xi * xi_y + v_zeta * zeta_x + v_eta * eta_x + v_xi * xi_x));
	stress_yz = (0.5*(v_zeta*zeta_z + v_eta * eta_z + v_xi * xi_z + w_zeta * zeta_y + w_eta * eta_y + w_xi * xi_y));
	stress_xz = (0.5*(w_zeta*zeta_x + w_eta * eta_x + w_xi * xi_x + u_zeta * zeta_z + u_eta * eta_z + u_xi * xi_z));

	stress_xx = (u_zeta*zeta_x + u_eta * eta_x + u_xi * xi_x);
	stress_yy = (v_zeta*zeta_y + v_eta * eta_y + v_xi * xi_y);
	stress_zz = (w_zeta*zeta_z + w_eta * eta_z + w_xi * xi_z);

	mod_stress = sqrt(2.0*((stress_xy*stress_xy + stress_yz * stress_yz + stress_xz * stress_xz) + (2.0 / 3.0)*(stress_xx*stress_xx + stress_yy * stress_yy + stress_zz * stress_zz)));

	DEL_ZETA = (1.0 / cbrt(DET_I));
	DEL_ETA = (1.0 / cbrt(DET_I));
	DEL_XI = (1.0 / cbrt(DET_I));

	/**********Be carefull with Rx, Ry, Rz******************************/
	/************THE RX Ry and Ry are different from paper to avoid zero divide*************/
	/************STUDY PAPER BEFORE TOUCHING Eddy viscosity ******************************/

	if (mod_stress != 0.0)
	{
		Rx = (DEL_ZETA*(mod_stress)) / sqrt(KS_x);
		Ry = (DEL_ZETA*(mod_stress)) / sqrt(KS_y);
		Rz = (DEL_ZETA*(mod_stress)) / sqrt(KS_z);

		EDDY_x = (C_MTS / (1.0 + (1.0*Rx / (C_T))))*DEL_ZETA*sqrt(KS_x);
		EDDY_y = (C_MTS / (1.0 + (1.0*Ry / (C_T))))*DEL_ZETA*sqrt(KS_y);
		EDDY_z = (C_MTS / (1.0 + (1.0*Rz / (C_T))))*DEL_ZETA*sqrt(KS_z);
	}
	if (mod_stress == 0.0)
	{
		EDDY_x = (C_MTS)*DEL_ZETA*sqrt(KS_x);
		EDDY_y = (C_MTS)*DEL_ZETA*sqrt(KS_y);
		EDDY_z = (C_MTS)*DEL_ZETA*sqrt(KS_z);
	}
	TAU_SGS_XY = -2.0*U0_I * EDDY_x*(stress_xy_star);
	TAU_SGS_YZ = -2.0*U0_I * EDDY_y*(stress_yz_star);
	TAU_SGS_XZ = -2.0*U0_I * EDDY_z*(stress_xz_star);

	TAU_SGS_XX = -2.0*U0_I * EDDY_x*stress_xx;
	TAU_SGS_YY = -2.0*U0_I * EDDY_y*stress_yy;
	TAU_SGS_ZZ = -2.0*U0_I * EDDY_z*stress_zz;
	/**********************************SGS ENERGY FLUX TERMS*******************************************************************/

	Q_x = (1.0 / (1.4*0.4*Mach*Mach))*(-1.0)*(EDDY_x)*(t_zeta*zeta_x + t_eta * eta_x + t_xi * xi_x);
	Q_y = (1.0 / (1.4*0.4*Mach*Mach))*(-1.0)*(EDDY_y)*(t_zeta*zeta_y + t_eta * eta_y + t_xi * xi_y);
	Q_z = (1.0 / (1.4*0.4*Mach*Mach))*(-1.0)*(EDDY_z)*(t_zeta*zeta_z + t_eta * eta_z + t_xi * xi_z);

	DEL_X = TAU_SGS_XX * FLOW_u + TAU_SGS_XY * FLOW_u + TAU_SGS_XZ * FLOW_u;
	DEL_Y = TAU_SGS_XY * FLOW_v + TAU_SGS_YY * FLOW_v + TAU_SGS_YZ * FLOW_v;
	DEL_Z = TAU_SGS_XZ * FLOW_w + TAU_SGS_YZ * FLOW_w + TAU_SGS_ZZ * FLOW_w;

	H_SGS_X = Q_x + DEL_X;
	H_SGS_Y = Q_y + DEL_Y;
	H_SGS_Z = Q_z + DEL_Z;

	D_SGS_X = 0.4*Q_x;
	D_SGS_Y = 0.4*Q_y;
	D_SGS_Z = 0.4*Q_z;
	/********************************************************************************************************/
    tF1_C.N0 = FLOW_rho * FLOW_u;
    tF1_C.N1 = FLOW_rho * FLOW_u * FLOW_u + FLOW_p;
    tF1_C.N2 = FLOW_rho * FLOW_v * FLOW_u;
    tF1_C.N3 = FLOW_rho * FLOW_u * FLOW_w;
    tF1_C.N4 = (((FLOW_p / 0.4) + 0.5*FLOW_rho * (FLOW_u * FLOW_u + FLOW_v * FLOW_v + FLOW_w * FLOW_w)) + FLOW_p)*FLOW_u;

    tE1_C.N0 = FLOW_rho * FLOW_v;
    tE1_C.N1 = FLOW_rho * FLOW_v * FLOW_u;
    tE1_C.N2 = FLOW_rho * FLOW_v * FLOW_v + FLOW_p;
    tE1_C.N3 = FLOW_rho * FLOW_v * FLOW_w;
    tE1_C.N4 = (((FLOW_p / 0.4) + 0.5*FLOW_rho * (FLOW_u * FLOW_u + FLOW_v * FLOW_v + FLOW_w * FLOW_w)) + FLOW_p)*FLOW_v;

    tG1_C.N0 = FLOW_rho * FLOW_w;
    tG1_C.N1 = FLOW_rho * FLOW_u * FLOW_w;
    tG1_C.N2 = FLOW_rho * FLOW_v * FLOW_w;
    tG1_C.N3 = FLOW_rho * FLOW_w * FLOW_w + FLOW_p;
    tG1_C.N4 = (((FLOW_p / 0.4) + 0.5*FLOW_rho * (FLOW_u * FLOW_u + FLOW_v * FLOW_v + FLOW_w * FLOW_w)) + FLOW_p)*FLOW_w;


    tFv1_C.N0 = 0.0;
    tFv1_C.N1 = (-1.0)*tauzz + TAU_SGS_XX;
    tFv1_C.N2 = (-1.0)*tauze + TAU_SGS_XY;
    tFv1_C.N3 = (-1.0)*tauzx + TAU_SGS_XZ;
    tFv1_C.N4 = (((-1.0)*FLOW_u * tauzz - FLOW_v * tauze - FLOW_w * tauzx) + qz + H_SGS_X + D_SGS_X);

    tEv1_C.N0 = 0.0;
    tEv1_C.N1 = (-1.0)*tauze + TAU_SGS_XY;
    tEv1_C.N2 = (-1.0)*tauee + TAU_SGS_YY;
    tEv1_C.N3 = (-1.0)*tauex + TAU_SGS_YZ;
    tEv1_C.N4 = (((-1.0)*FLOW_u * tauze - FLOW_v * tauee - FLOW_w * tauex) + qe + H_SGS_Y + D_SGS_Y);

    tGv1_C.N0 = 0.0;
    tGv1_C.N1 = (-1.0)*tauzx + TAU_SGS_XZ;
    tGv1_C.N2 = (-1.0)*tauex + TAU_SGS_YZ;
    tGv1_C.N3 = (-1.0)*tauxx + TAU_SGS_ZZ;
    tGv1_C.N4 = (((-1.0)*FLOW_u * tauzx - FLOW_v * tauex - FLOW_w * tauxx) + qx + H_SGS_Z + D_SGS_Z);

	F_C[Gid].N0 = DET_I*(zeta_x*tF1_C.N0 + zeta_y*tE1_C.N0 + zeta_z*tG1_C.N0);
	F_C[Gid].N1 = DET_I*(zeta_x*tF1_C.N1 + zeta_y*tE1_C.N1 + zeta_z*tG1_C.N1); 
	F_C[Gid].N2 = DET_I*(zeta_x*tF1_C.N2 + zeta_y*tE1_C.N2 + zeta_z*tG1_C.N2); 
	F_C[Gid].N3 = DET_I*(zeta_x*tF1_C.N3 + zeta_y*tE1_C.N3 + zeta_z*tG1_C.N3); 
	F_C[Gid].N4 = DET_I*(zeta_x*tF1_C.N4 + zeta_y*tE1_C.N4 + zeta_z*tG1_C.N4);  

	E_C[Gid].N0 = DET_I*(eta_x*tF1_C.N0 + eta_y*tE1_C.N0 + eta_z*tG1_C.N0);
	E_C[Gid].N1 = DET_I*(eta_x*tF1_C.N1 + eta_y*tE1_C.N1 + eta_z*tG1_C.N1); 
	E_C[Gid].N2 = DET_I*(eta_x*tF1_C.N2 + eta_y*tE1_C.N2 + eta_z*tG1_C.N2); 
	E_C[Gid].N3 = DET_I*(eta_x*tF1_C.N3 + eta_y*tE1_C.N3 + eta_z*tG1_C.N3); 
	E_C[Gid].N4 = DET_I*(eta_x*tF1_C.N4 + eta_y*tE1_C.N4 + eta_z*tG1_C.N4);  

	G_C[Gid].N0 = DET_I*(xi_x*tF1_C.N0 + xi_y*tE1_C.N0 + xi_z*tG1_C.N0);
	G_C[Gid].N1 = DET_I*(xi_x*tF1_C.N1 + xi_y*tE1_C.N1 + xi_z*tG1_C.N1); 
	G_C[Gid].N2 = DET_I*(xi_x*tF1_C.N2 + xi_y*tE1_C.N2 + xi_z*tG1_C.N2); 
	G_C[Gid].N3 = DET_I*(xi_x*tF1_C.N3 + xi_y*tE1_C.N3 + xi_z*tG1_C.N3); 
	G_C[Gid].N4 = DET_I*(xi_x*tF1_C.N4 + xi_y*tE1_C.N4 + xi_z*tG1_C.N4);  

	Fv_C[Gid].N0 = DET_I*(zeta_x*tFv1_C.N0 + zeta_y*tEv1_C.N0 + zeta_z*tGv1_C.N0);
	Fv_C[Gid].N1 = DET_I*(zeta_x*tFv1_C.N1 + zeta_y*tEv1_C.N1 + zeta_z*tGv1_C.N1); 
	Fv_C[Gid].N2 = DET_I*(zeta_x*tFv1_C.N2 + zeta_y*tEv1_C.N2 + zeta_z*tGv1_C.N2); 
	Fv_C[Gid].N3 = DET_I*(zeta_x*tFv1_C.N3 + zeta_y*tEv1_C.N3 + zeta_z*tGv1_C.N3); 
	Fv_C[Gid].N4 = DET_I*(zeta_x*tFv1_C.N4 + zeta_y*tEv1_C.N4 + zeta_z*tGv1_C.N4);  

	Ev_C[Gid].N0 = DET_I*(eta_x*tFv1_C.N0 + eta_y*tEv1_C.N0 + eta_z*tGv1_C.N0);
	Ev_C[Gid].N1 = DET_I*(eta_x*tFv1_C.N1 + eta_y*tEv1_C.N1 + eta_z*tGv1_C.N1); 
	Ev_C[Gid].N2 = DET_I*(eta_x*tFv1_C.N2 + eta_y*tEv1_C.N2 + eta_z*tGv1_C.N2); 
	Ev_C[Gid].N3 = DET_I*(eta_x*tFv1_C.N3 + eta_y*tEv1_C.N3 + eta_z*tGv1_C.N3); 
	Ev_C[Gid].N4 = DET_I*(eta_x*tFv1_C.N4 + eta_y*tEv1_C.N4 + eta_z*tGv1_C.N4);  

	Gv_C[Gid].N0 = DET_I*(xi_x*tFv1_C.N0 + xi_y*tEv1_C.N0 + xi_z*tGv1_C.N0);
	Gv_C[Gid].N1 = DET_I*(xi_x*tFv1_C.N1 + xi_y*tEv1_C.N1 + xi_z*tGv1_C.N1); 
	Gv_C[Gid].N2 = DET_I*(xi_x*tFv1_C.N2 + xi_y*tEv1_C.N2 + xi_z*tGv1_C.N2); 
	Gv_C[Gid].N3 = DET_I*(xi_x*tFv1_C.N3 + xi_y*tEv1_C.N3 + xi_z*tGv1_C.N3); 
	Gv_C[Gid].N4 = DET_I*(xi_x*tFv1_C.N4 + xi_y*tEv1_C.N4 + xi_z*tGv1_C.N4);  

	FLOW[Gid].u = FLOW_u;
	FLOW[Gid].v = FLOW_v;
	FLOW[Gid].w = FLOW_w;
	FLOW[Gid].p = FLOW_p;
	FLOW[Gid].t = FLOW_t;
	FLOW[Gid].a = FLOW_a;
	FLOW[Gid].e = FLOW_e;
	FLOW[Gid].rho = FLOW_rho;
    FLOW[Gid].mu = FLOW_mu;

}


