#pragma OPENCL EXTENSION cl_khr_fp64 : enable
typedef struct __attribute__ ((packed)) tag_grid {
	double x, y, z;
	int N0, N1, N2, N3, N4, N5;
	int ID, loc;
} GRID;

typedef struct __attribute__ ((packed)) tag_DETERM {
	double ip, im, jp, jm, kp, km;
} DETERM;

typedef struct __attribute__ ((packed)) tag_Solver_Column_matrix {
	double N0, N1, N2, N3, N4;
} Solver_Column_matrix;

typedef struct __attribute__ ((packed)) tag_TRANSFORMED {
	double zeta_x, zeta_y, zeta_z, eta_x, eta_y, eta_z, xi_x, xi_y, xi_z;
	double zeta_xip, zeta_yip, zeta_zip, eta_xip, eta_yip, eta_zip, xi_xip, xi_yip, xi_zip, zeta_xim, zeta_yim, zeta_zim, eta_xim, eta_yim, eta_zim, xi_xim, xi_yim, xi_zim;
	double zeta_xjp, zeta_yjp, zeta_zjp, eta_xjp, eta_yjp, eta_zjp, xi_xjp, xi_yjp, xi_zjp, zeta_xjm, zeta_yjm, zeta_zjm, eta_xjm, eta_yjm, eta_zjm, xi_xjm, xi_yjm, xi_zjm;
	double zeta_xkp, zeta_ykp, zeta_zkp, eta_xkp, eta_ykp, eta_zkp, xi_xkp, xi_ykp, xi_zkp, zeta_xkm, zeta_ykm, zeta_zkm, eta_xkm, eta_ykm, eta_zkm, xi_xkm, xi_ykm, xi_zkm;
} TRANSFORMED;

typedef struct __attribute__ ((packed)) tag_flow_variables {
	double u, v, w, p, t, rho, mu, a, e;
} FLOW_VARIABLES;

typedef struct __attribute__ ((packed)) tag_ROE_AVER {
	double roe_u_ip, roe_v_ip, roe_w_ip, roe_p_ip, roe_h_ip, roe_rho_ip;
	double roe_u_jp, roe_v_jp, roe_w_jp, roe_p_jp, roe_h_jp, roe_rho_jp;
	double roe_u_kp, roe_v_kp, roe_w_kp, roe_p_kp, roe_h_kp, roe_rho_kp;
	double roe_u_im, roe_v_im, roe_w_im, roe_p_im, roe_h_im, roe_rho_im;
	double roe_u_jm, roe_v_jm, roe_w_jm, roe_p_jm, roe_h_jm, roe_rho_jm;
	double roe_u_km, roe_v_km, roe_w_km, roe_p_km, roe_h_km, roe_rho_km;
} ROE_AVERAGING;

__kernel void Weno(const int j, const double Mach, const double del_t, __global const DETERM *deter, __global const double *det,\
__global Solver_Column_matrix *U0, __global Solver_Column_matrix *U1, __global Solver_Column_matrix *U2,\
__global const GRID *Grid_P, __global Solver_Column_matrix *Q_C, __global const TRANSFORMED *metric, __global ROE_AVERAGING *ROE_AVER,\
__global Solver_Column_matrix *Fv, __global Solver_Column_matrix *Ev, __global Solver_Column_matrix *Gv,\
__global Solver_Column_matrix *F, __global Solver_Column_matrix *E, __global Solver_Column_matrix *G, __global FLOW_VARIABLES *FLOW)
{
    uint Gid = get_global_id(0);
    uint Lid = get_local_id(0);

    if (Gid == 1) {
        printf("Hello %d\n", __LINE__);
    }

    __private int k, m;
    __private int IP1, IP2, IP3, IM1, IM2, IM3;
    __private int JP1, JP2, JP3, JM1, JM2, JM3;
    __private int KP1, KP2, KP3, KM1, KM2, KM3;

    __private int LOC_I;
    __private int LOC_IP1, LOC_IP2, LOC_IP3, LOC_IM1, LOC_IM2, LOC_IM3;
    __private int LOC_JP1, LOC_JP2, LOC_JP3, LOC_JM1, LOC_JM2, LOC_JM3;
    __private int LOC_KP1, LOC_KP2, LOC_KP3, LOC_KM1, LOC_KM2, LOC_KM3;

    __private const double Epsilon = 10e-8;
	
	__private double alpha, beta, theta, kk;
	__private double kx_bar, ky_bar, kz_bar;
	__private double phi_sq;
    __private double roe_a;

    
    
    __private double Qi_half_p[3], Qi_half_m[3], Qi_halfn_m[3], Qi_half_np[3];
    __private double Qj_half_p[3], Qj_half_m[3], Qj_halfn_m[3], Qj_half_np[3];
    __private double Qk_half_p[3], Qk_half_m[3], Qk_halfn_m[3], Qk_half_np[3];
        
    __private double IS_Qip[3], IS_Qim[3], IS_Qinp[3], IS_Qinm[3];
    __private double IS_Qjp[3], IS_Qjm[3], IS_Qjnp[3], IS_Qjnm[3];
    __private double IS_Qkp[3], IS_Qkm[3], IS_Qknp[3], IS_Qknm[3];
    
    __private double w_Qip[3], w_Qim[3], w_Qinp[3], w_Qinm[3];
    __private double w_Qjp[3], w_Qjm[3], w_Qjnp[3], w_Qjnm[3];
    __private double w_Qkp[3], w_Qkm[3], w_Qknp[3], w_Qknm[3];
    
    __private double W_Qip[3], W_Qim[3], W_Qinp[3], W_Qinm[3];
    __private double W_Qjp[3], W_Qjm[3], W_Qjnp[3], W_Qjnm[3];
    __private double W_Qkp[3], W_Qkm[3], W_Qknp[3], W_Qknm[3];
    
    __private double Qi_iplus_half_pos_char[5], Qi_iplus_half_neg_char[5], Qi_iminus_half_pos_char[5], Qi_iminus_half_neg_char[5];
    __private double Qj_iplus_half_pos_char[5], Qj_iplus_half_neg_char[5], Qj_iminus_half_pos_char[5], Qj_iminus_half_neg_char[5];
    __private double Qk_iplus_half_pos_char[5], Qk_iplus_half_neg_char[5], Qk_iminus_half_pos_char[5], Qk_iminus_half_neg_char[5];
    
    __private double Qi_iplus_half_pos[5], Qi_iplus_half_neg[5], Qi_iminus_half_pos[5], Qi_iminus_half_neg[5];
    __private double Qj_iplus_half_pos[5], Qj_iplus_half_neg[5], Qj_iminus_half_pos[5], Qj_iminus_half_neg[5];
    __private double Qk_iplus_half_pos[5], Qk_iplus_half_neg[5], Qk_iminus_half_pos[5], Qk_iminus_half_neg[5];

    __private double F_ip_pos[5], F_ip_neg[5], F_im_pos[5], F_im_neg[5];
    __private double E_ip_pos[5], E_ip_neg[5], E_im_pos[5], E_im_neg[5];
    __private double G_ip_pos[5], G_ip_neg[5], G_im_pos[5], G_im_neg[5];
    
    __private double F_jp_pos[5], F_jp_neg[5], F_jm_pos[5], F_jm_neg[5];
    __private double E_jp_pos[5], E_jp_neg[5], E_jm_pos[5], E_jm_neg[5];
    __private double G_jp_pos[5], G_jp_neg[5], G_jm_pos[5], G_jm_neg[5];
    
    __private double F_kp_pos[5], F_kp_neg[5], F_km_pos[5], F_km_neg[5];
    __private double E_kp_pos[5], E_kp_neg[5], E_km_pos[5], E_km_neg[5];
    __private double G_kp_pos[5], G_kp_neg[5], G_km_pos[5], G_km_neg[5];

    __private double F_ip_pos_comp[5], F_ip_neg_comp[5], F_im_pos_comp[5], F_im_neg_comp[5];
    __private double E_ip_pos_comp[5], E_ip_neg_comp[5], E_im_pos_comp[5], E_im_neg_comp[5];
    __private double G_ip_pos_comp[5], G_ip_neg_comp[5], G_im_pos_comp[5], G_im_neg_comp[5];

    __private double eigen_Qip[5], eigen_Qim[5], eigen_Qinp[5], eigen_Qinm[5]; 
    __private double eigen_Qjp[5], eigen_Qjm[5], eigen_Qjnp[5], eigen_Qjnm[5]; 
    __private double eigen_Qkp[5], eigen_Qkm[5], eigen_Qknp[5], eigen_Qknm[5]; 

    __private double DET_I, DET_IP, DET_IM, DET_JP, DET_JM, DET_KP, DET_KM; 

    __private Solver_Column_matrix l_eigen_Q[5];
    __private Solver_Column_matrix r_eigen_Qip[5], r_eigen_Qim[5];
    __private Solver_Column_matrix r_eigen_Qjp[5], r_eigen_Qjm[5];
    __private Solver_Column_matrix r_eigen_Qkp[5], r_eigen_Qkm[5];

    __private double alpha_u_ip, alpha_u_im, alpha_u_jp, alpha_u_jm, alpha_u_kp, alpha_u_km;
    __private double alpha_v_ip, alpha_v_im, alpha_v_jp, alpha_v_jm, alpha_v_kp, alpha_v_km;
    __private double alpha_w_ip, alpha_w_im, alpha_w_jp, alpha_w_jm, alpha_w_kp, alpha_w_km;
    
    __private double final_U[5], dFv[5], dEv[5], dGv[5], dF[5], dE[5], dG[5], L[5]; 
    __private double F_ip[5], F_im[5], E_jp[5], E_jm[5], G_kp[5], G_km[5];

    __private double d2F_d2z_ip[5], d2F_d2z_im[5], d2E_d2e_ip[5], d2E_d2e_im[5], d2G_d2x_ip[5], d2G_d2x_im[5];
    __private double d4F_d4z_ip[5], d4F_d4z_im[5], d4E_d4e_ip[5], d4E_d4e_im[5], d4G_d4x_ip[5], d4G_d4x_im[5];

    __private double F_ip_h[5], F_im_h[5], E_ip_h[5], E_im_h[5], G_ip_h[5], G_im_h[5];

    __private Solver_Column_matrix Q_I;
    __private Solver_Column_matrix Q_IP1, Q_IP2, Q_IP3;
    __private Solver_Column_matrix Q_IM1, Q_IM2, Q_IM3;
    __private Solver_Column_matrix Q_JP1, Q_JP2, Q_JP3;
    __private Solver_Column_matrix Q_JM1, Q_JM2, Q_JM3;
    __private Solver_Column_matrix Q_KP1, Q_KP2, Q_KP3;
    __private Solver_Column_matrix Q_KM1, Q_KM2, Q_KM3;
    __private Solver_Column_matrix Qi_iplus_IP1, Qi_iplus_IP2, Qi_iplus_IP3;
    __private Solver_Column_matrix Qi_iplus_IM1, Qi_iplus_IM2, Qi_iplus_IM3;
    __private Solver_Column_matrix Qi_iminus_IP1, Qi_iminus_IP2, Qi_iminus_IP3;
    __private Solver_Column_matrix Qi_iminus_IM1,  Qi_iminus_IM2,  Qi_iminus_IM3;
    __private Solver_Column_matrix Qj_iplus_JP1,  Qj_iplus_JP2,  Qj_iplus_JP3;
    __private Solver_Column_matrix Qj_iplus_JM1,  Qj_iplus_JM2,  Qj_iplus_JM3;
    __private Solver_Column_matrix Qj_iminus_JP1,  Qj_iminus_JP2,  Qj_iminus_JP3;
    __private Solver_Column_matrix Qj_iminus_JM1,  Qj_iminus_JM2,  Qj_iminus_JM3;
    __private Solver_Column_matrix Qk_iplus_KP1,  Qk_iplus_KP2,  Qk_iplus_KP3;
    __private Solver_Column_matrix Qk_iplus_KM1,  Qk_iplus_KM2,  Qk_iplus_KM3;
    __private Solver_Column_matrix Qk_iminus_KP1,  Qk_iminus_KP2,  Qk_iminus_KP3;
    __private Solver_Column_matrix Qk_iminus_KM1,  Qk_iminus_KM2,  Qk_iminus_KM3;
    __private Solver_Column_matrix Qi_iplus_I,  Qi_iminus_I;
    __private Solver_Column_matrix Qj_iplus_I,  Qj_iminus_I;
    __private Solver_Column_matrix Qk_iplus_I,  Qk_iminus_I;

    if (Gid == 1) {
        printf("Hello %d\n", __LINE__);
    }
    /******************************************************************************************************/

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
    
    DET_IP = deter[Gid].ip;
    DET_IM = deter[Gid].im;
    DET_JP = deter[Gid].jp;
    DET_JM = deter[Gid].jm;
    DET_KP = deter[Gid].kp;
    DET_KM = deter[Gid].km;

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

//    barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

    Q_I.N0 = Q_C[Gid].N0;
    Q_I.N1 = Q_C[Gid].N1;
    Q_I.N2 = Q_C[Gid].N2;
    Q_I.N3 = Q_C[Gid].N3;
    Q_I.N4 = Q_C[Gid].N4;

    Q_IP1.N0 = Q_C[IP1].N0;
    Q_IP1.N1 = Q_C[IP1].N1;
    Q_IP1.N2 = Q_C[IP1].N2;
    Q_IP1.N3 = Q_C[IP1].N3;
    Q_IP1.N4 = Q_C[IP1].N4;

    Q_IP2.N0 = Q_C[IP2].N0;
    Q_IP2.N1 = Q_C[IP2].N1;
    Q_IP2.N2 = Q_C[IP2].N2;
    Q_IP2.N3 = Q_C[IP2].N3;
    Q_IP2.N4 = Q_C[IP2].N4;

    Q_IP3.N0 = Q_C[IP3].N0;
    Q_IP3.N1 = Q_C[IP3].N1;
    Q_IP3.N2 = Q_C[IP3].N2;
    Q_IP3.N3 = Q_C[IP3].N3;
    Q_IP3.N4 = Q_C[IP3].N4;


    Q_JP1.N0 = Q_C[JP1].N0;
    Q_JP1.N1 = Q_C[JP1].N1;
    Q_JP1.N2 = Q_C[JP1].N2;
    Q_JP1.N3 = Q_C[JP1].N3;
    Q_JP1.N4 = Q_C[JP1].N4;

    Q_JP2.N0 = Q_C[JP2].N0;
    Q_JP2.N1 = Q_C[JP2].N1;
    Q_JP2.N2 = Q_C[JP2].N2;
    Q_JP2.N3 = Q_C[JP2].N3;
    Q_JP2.N4 = Q_C[JP2].N4;

    Q_JP3.N0 = Q_C[JP3].N0;
    Q_JP3.N1 = Q_C[JP3].N1;
    Q_JP3.N2 = Q_C[JP3].N2;
    Q_JP3.N3 = Q_C[JP3].N3;
    Q_JP3.N4 = Q_C[JP3].N4;
    
    
    Q_KP1.N0 = Q_C[KP1].N0;
    Q_KP1.N1 = Q_C[KP1].N1;
    Q_KP1.N2 = Q_C[KP1].N2;
    Q_KP1.N3 = Q_C[KP1].N3;
    Q_KP1.N4 = Q_C[KP1].N4;

    Q_KP2.N0 = Q_C[KP2].N0;
    Q_KP2.N1 = Q_C[KP2].N1;
    Q_KP2.N2 = Q_C[KP2].N2;
    Q_KP2.N3 = Q_C[KP2].N3;
    Q_KP2.N4 = Q_C[KP2].N4;

    Q_KP3.N0 = Q_C[KP3].N0;
    Q_KP3.N1 = Q_C[KP3].N1;
    Q_KP3.N2 = Q_C[KP3].N2;
    Q_KP3.N3 = Q_C[KP3].N3;
    Q_KP3.N4 = Q_C[KP3].N4;


    Q_IM1.N0 = Q_C[IM1].N0;
    Q_IM1.N1 = Q_C[IM1].N1;
    Q_IM1.N2 = Q_C[IM1].N2;
    Q_IM1.N3 = Q_C[IM1].N3;
    Q_IM1.N4 = Q_C[IM1].N4;

    Q_IM2.N0 = Q_C[IM2].N0;
    Q_IM2.N1 = Q_C[IM2].N1;
    Q_IM2.N2 = Q_C[IM2].N2;
    Q_IM2.N3 = Q_C[IM2].N3;
    Q_IM2.N4 = Q_C[IM2].N4;

    Q_IM3.N0 = Q_C[IM3].N0;
    Q_IM3.N1 = Q_C[IM3].N1;
    Q_IM3.N2 = Q_C[IM3].N2;
    Q_IM3.N3 = Q_C[IM3].N3;
    Q_IM3.N4 = Q_C[IM3].N4;


    Q_JM1.N0 = Q_C[JM1].N0;
    Q_JM1.N1 = Q_C[JM1].N1;
    Q_JM1.N2 = Q_C[JM1].N2;
    Q_JM1.N3 = Q_C[JM1].N3;
    Q_JM1.N4 = Q_C[JM1].N4;

    Q_JM2.N0 = Q_C[JM2].N0;
    Q_JM2.N1 = Q_C[JM2].N1;
    Q_JM2.N2 = Q_C[JM2].N2;
    Q_JM2.N3 = Q_C[JM2].N3;
    Q_JM2.N4 = Q_C[JM2].N4;

    Q_JM3.N0 = Q_C[JM3].N0;
    Q_JM3.N1 = Q_C[JM3].N1;
    Q_JM3.N2 = Q_C[JM3].N2;
    Q_JM3.N3 = Q_C[JM3].N3;
    Q_JM3.N4 = Q_C[JM3].N4;
    
    
    Q_KM1.N0 = Q_C[KM1].N0;
    Q_KM1.N1 = Q_C[KM1].N1;
    Q_KM1.N2 = Q_C[KM1].N2;
    Q_KM1.N3 = Q_C[KM1].N3;
    Q_KM1.N4 = Q_C[KM1].N4;

    Q_KM2.N0 = Q_C[KM2].N0;
    Q_KM2.N1 = Q_C[KM2].N1;
    Q_KM2.N2 = Q_C[KM2].N2;
    Q_KM2.N3 = Q_C[KM2].N3;
    Q_KM2.N4 = Q_C[KM2].N4;

    Q_KM3.N0 = Q_C[KM3].N0;
    Q_KM3.N1 = Q_C[KM3].N1;
    Q_KM3.N2 = Q_C[KM3].N2;
    Q_KM3.N3 = Q_C[KM3].N3;
    Q_KM3.N4 = Q_C[KM3].N4;
    
    if (Gid == 1) {
        printf("Hello %d\n", __LINE__);
    }

//    barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
    /***********************************************************************************************/
	/****************************************** Qi_iplus********************************************/
    kk = sqrt(metric[Gid].zeta_xip*metric[Gid].zeta_xip + metric[Gid].zeta_yip * metric[Gid].zeta_yip + metric[Gid].zeta_zip * metric[Gid].zeta_zip);
    roe_a = sqrt(0.4*(ROE_AVER[Gid].roe_h_ip - 0.5*(ROE_AVER[Gid].roe_u_ip * ROE_AVER[Gid].roe_u_ip + ROE_AVER[Gid].roe_v_ip * ROE_AVER[Gid].roe_v_ip + ROE_AVER[Gid].roe_w_ip * ROE_AVER[Gid].roe_w_ip)));

	kx_bar = metric[Gid].zeta_xip / kk;
	ky_bar = metric[Gid].zeta_yip / kk;
	kz_bar = metric[Gid].zeta_zip / kk;

	theta = kx_bar * ROE_AVER[Gid].roe_u_ip + ky_bar * ROE_AVER[Gid].roe_v_ip + kz_bar * ROE_AVER[Gid].roe_w_ip;
	phi_sq = 0.5*0.4*(ROE_AVER[Gid].roe_u_ip * ROE_AVER[Gid].roe_u_ip + ROE_AVER[Gid].roe_v_ip * ROE_AVER[Gid].roe_v_ip + ROE_AVER[Gid].roe_w_ip * ROE_AVER[Gid].roe_w_ip);
	alpha = ROE_AVER[Gid].roe_rho_ip / (sqrt(2.0)*roe_a);
	beta = 1.0 / (sqrt(2.0)*ROE_AVER[Gid].roe_rho_ip * roe_a);

	r_eigen_Qip[0].N0 = kx_bar;
	r_eigen_Qip[0].N1 = ky_bar;
	r_eigen_Qip[0].N2 = kz_bar;
	r_eigen_Qip[0].N3 = alpha;
	r_eigen_Qip[0].N4 = alpha;

	r_eigen_Qip[1].N0 = kx_bar * ROE_AVER[Gid].roe_u_ip;
	r_eigen_Qip[1].N1 = ky_bar * ROE_AVER[Gid].roe_u_ip - kz_bar * ROE_AVER[Gid].roe_rho_ip;
	r_eigen_Qip[1].N2 = kz_bar * ROE_AVER[Gid].roe_u_ip + ky_bar * ROE_AVER[Gid].roe_rho_ip;
	r_eigen_Qip[1].N3 = alpha * (ROE_AVER[Gid].roe_u_ip + kx_bar * roe_a);
	r_eigen_Qip[1].N4 = alpha * (ROE_AVER[Gid].roe_u_ip - kx_bar * roe_a);

	r_eigen_Qip[2].N0 = kx_bar * ROE_AVER[Gid].roe_v_ip + kz_bar * ROE_AVER[Gid].roe_rho_ip;
	r_eigen_Qip[2].N1 = ky_bar * ROE_AVER[Gid].roe_v_ip;
	r_eigen_Qip[2].N2 = kz_bar * ROE_AVER[Gid].roe_v_ip - kx_bar * ROE_AVER[Gid].roe_rho_ip;
	r_eigen_Qip[2].N3 = alpha * (ROE_AVER[Gid].roe_v_ip + ky_bar * roe_a);
	r_eigen_Qip[2].N4 = alpha * (ROE_AVER[Gid].roe_v_ip - ky_bar * roe_a);

	r_eigen_Qip[3].N0 = kx_bar * ROE_AVER[Gid].roe_w_ip - ky_bar * ROE_AVER[Gid].roe_rho_ip;
	r_eigen_Qip[3].N1 = ky_bar * ROE_AVER[Gid].roe_w_ip + kx_bar * ROE_AVER[Gid].roe_rho_ip;
	r_eigen_Qip[3].N2 = kz_bar * ROE_AVER[Gid].roe_w_ip;
	r_eigen_Qip[3].N3 = alpha * (ROE_AVER[Gid].roe_w_ip + kz_bar * roe_a);
	r_eigen_Qip[3].N4 = alpha * (ROE_AVER[Gid].roe_w_ip - kz_bar * roe_a);

	r_eigen_Qip[4].N0 = ((kx_bar*phi_sq) / 0.4) + ROE_AVER[Gid].roe_rho_ip * (kz_bar*ROE_AVER[Gid].roe_v_ip - ky_bar * ROE_AVER[Gid].roe_w_ip);
	r_eigen_Qip[4].N1 = ((ky_bar*phi_sq) / 0.4) + ROE_AVER[Gid].roe_rho_ip * (kx_bar*ROE_AVER[Gid].roe_w_ip - kz_bar * ROE_AVER[Gid].roe_u_ip);
	r_eigen_Qip[4].N2 = ((kz_bar*phi_sq) / 0.4) + ROE_AVER[Gid].roe_rho_ip * (ky_bar*ROE_AVER[Gid].roe_u_ip - kx_bar * ROE_AVER[Gid].roe_v_ip);
	r_eigen_Qip[4].N3 = alpha * (((phi_sq + roe_a * roe_a) / (0.4)) + theta * roe_a);
	r_eigen_Qip[4].N4 = alpha * (((phi_sq + roe_a * roe_a) / (0.4)) - theta * roe_a);

    if (Gid == 1) {
        printf("Hello %d\n", __LINE__);
    }
	/*************************************************************************************************/

	l_eigen_Q[0].N0 = kx_bar * (1.0 - (phi_sq / (roe_a * roe_a))) - ((kz_bar*ROE_AVER[Gid].roe_v_ip - ky_bar * ROE_AVER[Gid].roe_w_ip) / ROE_AVER[Gid].roe_rho_ip);
	l_eigen_Q[0].N1 = kx_bar * 0.4*ROE_AVER[Gid].roe_u_ip / (roe_a * roe_a);
	l_eigen_Q[0].N2 = (kx_bar*0.4*ROE_AVER[Gid].roe_v_ip / (roe_a * roe_a)) + (kz_bar / ROE_AVER[Gid].roe_rho_ip);
	l_eigen_Q[0].N3 = (kx_bar*0.4*ROE_AVER[Gid].roe_w_ip / (roe_a * roe_a)) - (ky_bar / ROE_AVER[Gid].roe_rho_ip);
	l_eigen_Q[0].N4 = -kx_bar * 0.4 / (roe_a * roe_a);

	l_eigen_Q[1].N0 = ky_bar * (1.0 - (phi_sq / (roe_a * roe_a))) - ((kx_bar*ROE_AVER[Gid].roe_w_ip - kz_bar * ROE_AVER[Gid].roe_u_ip) / ROE_AVER[Gid].roe_rho_ip);
	l_eigen_Q[1].N1 = (ky_bar*0.4*ROE_AVER[Gid].roe_u_ip / (roe_a * roe_a)) - kz_bar / ROE_AVER[Gid].roe_rho_ip;
	l_eigen_Q[1].N2 = ky_bar * 0.4*ROE_AVER[Gid].roe_v_ip / (roe_a * roe_a);
	l_eigen_Q[1].N3 = (ky_bar*0.4*ROE_AVER[Gid].roe_w_ip / (roe_a * roe_a)) + (kx_bar / ROE_AVER[Gid].roe_rho_ip);
	l_eigen_Q[1].N4 = -ky_bar * 0.4 / (roe_a * roe_a);

	l_eigen_Q[2].N0 = kz_bar * (1.0 - (phi_sq / (roe_a * roe_a))) - ((ky_bar*ROE_AVER[Gid].roe_u_ip - kx_bar * ROE_AVER[Gid].roe_v_ip) / ROE_AVER[Gid].roe_rho_ip);
	l_eigen_Q[2].N1 = kz_bar * 0.4*ROE_AVER[Gid].roe_u_ip / (roe_a * roe_a) + (ky_bar / ROE_AVER[Gid].roe_rho_ip);
	l_eigen_Q[2].N2 = kz_bar * 0.4*ROE_AVER[Gid].roe_v_ip / (roe_a * roe_a) - (kx_bar / ROE_AVER[Gid].roe_rho_ip);
	l_eigen_Q[2].N3 = kz_bar * 0.4*ROE_AVER[Gid].roe_w_ip / (roe_a * roe_a);
	l_eigen_Q[2].N4 = -kz_bar * 0.4 / (roe_a * roe_a);

	l_eigen_Q[3].N0 = beta * (phi_sq - theta * roe_a);
	l_eigen_Q[3].N1 = -beta * (0.4*ROE_AVER[Gid].roe_u_ip - kx_bar * roe_a);
	l_eigen_Q[3].N2 = -beta * (0.4*ROE_AVER[Gid].roe_v_ip - ky_bar * roe_a);
	l_eigen_Q[3].N3 = -beta * (0.4*ROE_AVER[Gid].roe_w_ip - kz_bar * roe_a);
	l_eigen_Q[3].N4 = beta * 0.4;

	l_eigen_Q[4].N0 = beta * (phi_sq + theta * roe_a);
	l_eigen_Q[4].N1 = -beta * (0.4*ROE_AVER[Gid].roe_u_ip + kx_bar * roe_a);
	l_eigen_Q[4].N2 = -beta * (0.4*ROE_AVER[Gid].roe_v_ip + ky_bar * roe_a);
	l_eigen_Q[4].N3 = -beta * (0.4*ROE_AVER[Gid].roe_w_ip + kz_bar * roe_a);
	l_eigen_Q[4].N4 = beta * 0.4;

	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/
    Qi_iplus_IP3.N0 = l_eigen_Q[0].N0 * Q_IP3.N0 + l_eigen_Q[0].N1 * Q_IP3.N1 + l_eigen_Q[0].N2 * Q_IP3.N2 + l_eigen_Q[0].N3 * Q_IP3.N3 + l_eigen_Q[0].N4 * Q_IP3.N4;
    Qi_iplus_IP3.N1 = l_eigen_Q[1].N0 * Q_IP3.N0 + l_eigen_Q[1].N1 * Q_IP3.N1 + l_eigen_Q[1].N2 * Q_IP3.N2 + l_eigen_Q[1].N3 * Q_IP3.N3 + l_eigen_Q[1].N4 * Q_IP3.N4;
    Qi_iplus_IP3.N2 = l_eigen_Q[2].N0 * Q_IP3.N0 + l_eigen_Q[2].N1 * Q_IP3.N1 + l_eigen_Q[2].N2 * Q_IP3.N2 + l_eigen_Q[2].N3 * Q_IP3.N3 + l_eigen_Q[2].N4 * Q_IP3.N4;
    Qi_iplus_IP3.N3 = l_eigen_Q[3].N0 * Q_IP3.N0 + l_eigen_Q[3].N1 * Q_IP3.N1 + l_eigen_Q[3].N2 * Q_IP3.N2 + l_eigen_Q[3].N3 * Q_IP3.N3 + l_eigen_Q[3].N4 * Q_IP3.N4;
    Qi_iplus_IP3.N4 = l_eigen_Q[4].N0 * Q_IP3.N0 + l_eigen_Q[4].N1 * Q_IP3.N1 + l_eigen_Q[4].N2 * Q_IP3.N2 + l_eigen_Q[4].N3 * Q_IP3.N3 + l_eigen_Q[4].N4 * Q_IP3.N4;
    
    Qi_iplus_IP2.N0 = l_eigen_Q[0].N0 * Q_IP2.N0 + l_eigen_Q[0].N1 * Q_IP2.N1 + l_eigen_Q[0].N2 * Q_IP2.N2 + l_eigen_Q[0].N3 * Q_IP2.N3 + l_eigen_Q[0].N4 * Q_IP2.N4;
    Qi_iplus_IP2.N1 = l_eigen_Q[1].N0 * Q_IP2.N0 + l_eigen_Q[1].N1 * Q_IP2.N1 + l_eigen_Q[1].N2 * Q_IP2.N2 + l_eigen_Q[1].N3 * Q_IP2.N3 + l_eigen_Q[1].N4 * Q_IP2.N4;
    Qi_iplus_IP2.N2 = l_eigen_Q[2].N0 * Q_IP2.N0 + l_eigen_Q[2].N1 * Q_IP2.N1 + l_eigen_Q[2].N2 * Q_IP2.N2 + l_eigen_Q[2].N3 * Q_IP2.N3 + l_eigen_Q[2].N4 * Q_IP2.N4;
    Qi_iplus_IP2.N3 = l_eigen_Q[3].N0 * Q_IP2.N0 + l_eigen_Q[3].N1 * Q_IP2.N1 + l_eigen_Q[3].N2 * Q_IP2.N2 + l_eigen_Q[3].N3 * Q_IP2.N3 + l_eigen_Q[3].N4 * Q_IP2.N4;
    Qi_iplus_IP2.N4 = l_eigen_Q[4].N0 * Q_IP2.N0 + l_eigen_Q[4].N1 * Q_IP2.N1 + l_eigen_Q[4].N2 * Q_IP2.N2 + l_eigen_Q[4].N3 * Q_IP2.N3 + l_eigen_Q[4].N4 * Q_IP2.N4;
    
    Qi_iplus_IP1.N0 = l_eigen_Q[0].N0 * Q_IP1.N0 + l_eigen_Q[0].N1 * Q_IP1.N1 + l_eigen_Q[0].N2 * Q_IP1.N2 + l_eigen_Q[0].N3 * Q_IP1.N3 + l_eigen_Q[0].N4 * Q_IP1.N4;
    Qi_iplus_IP1.N1 = l_eigen_Q[1].N0 * Q_IP1.N0 + l_eigen_Q[1].N1 * Q_IP1.N1 + l_eigen_Q[1].N2 * Q_IP1.N2 + l_eigen_Q[1].N3 * Q_IP1.N3 + l_eigen_Q[1].N4 * Q_IP1.N4;
    Qi_iplus_IP1.N2 = l_eigen_Q[2].N0 * Q_IP1.N0 + l_eigen_Q[2].N1 * Q_IP1.N1 + l_eigen_Q[2].N2 * Q_IP1.N2 + l_eigen_Q[2].N3 * Q_IP1.N3 + l_eigen_Q[2].N4 * Q_IP1.N4;
    Qi_iplus_IP1.N3 = l_eigen_Q[3].N0 * Q_IP1.N0 + l_eigen_Q[3].N1 * Q_IP1.N1 + l_eigen_Q[3].N2 * Q_IP1.N2 + l_eigen_Q[3].N3 * Q_IP1.N3 + l_eigen_Q[3].N4 * Q_IP1.N4;
    Qi_iplus_IP1.N4 = l_eigen_Q[4].N0 * Q_IP1.N0 + l_eigen_Q[4].N1 * Q_IP1.N1 + l_eigen_Q[4].N2 * Q_IP1.N2 + l_eigen_Q[4].N3 * Q_IP1.N3 + l_eigen_Q[4].N4 * Q_IP1.N4;
    
    Qi_iplus_I.N0 = l_eigen_Q[0].N0 * Q_I.N0 + l_eigen_Q[0].N1 * Q_I.N1 + l_eigen_Q[0].N2 * Q_I.N2 + l_eigen_Q[0].N3 * Q_I.N3 + l_eigen_Q[0].N4 * Q_I.N4;
    Qi_iplus_I.N1 = l_eigen_Q[1].N0 * Q_I.N0 + l_eigen_Q[1].N1 * Q_I.N1 + l_eigen_Q[1].N2 * Q_I.N2 + l_eigen_Q[1].N3 * Q_I.N3 + l_eigen_Q[1].N4 * Q_I.N4;
    Qi_iplus_I.N2 = l_eigen_Q[2].N0 * Q_I.N0 + l_eigen_Q[2].N1 * Q_I.N1 + l_eigen_Q[2].N2 * Q_I.N2 + l_eigen_Q[2].N3 * Q_I.N3 + l_eigen_Q[2].N4 * Q_I.N4;
    Qi_iplus_I.N3 = l_eigen_Q[3].N0 * Q_I.N0 + l_eigen_Q[3].N1 * Q_I.N1 + l_eigen_Q[3].N2 * Q_I.N2 + l_eigen_Q[3].N3 * Q_I.N3 + l_eigen_Q[3].N4 * Q_I.N4;
    Qi_iplus_I.N4 = l_eigen_Q[4].N0 * Q_I.N0 + l_eigen_Q[4].N1 * Q_I.N1 + l_eigen_Q[4].N2 * Q_I.N2 + l_eigen_Q[4].N3 * Q_I.N3 + l_eigen_Q[4].N4 * Q_I.N4;
    
    Qi_iplus_IM1.N0 = l_eigen_Q[0].N0 * Q_IM1.N0 + l_eigen_Q[0].N1 * Q_IM1.N1 + l_eigen_Q[0].N2 * Q_IM1.N2 + l_eigen_Q[0].N3 * Q_IM1.N3 + l_eigen_Q[0].N4 * Q_IM1.N4;
    Qi_iplus_IM1.N1 = l_eigen_Q[1].N0 * Q_IM1.N0 + l_eigen_Q[1].N1 * Q_IM1.N1 + l_eigen_Q[1].N2 * Q_IM1.N2 + l_eigen_Q[1].N3 * Q_IM1.N3 + l_eigen_Q[1].N4 * Q_IM1.N4;
    Qi_iplus_IM1.N2 = l_eigen_Q[2].N0 * Q_IM1.N0 + l_eigen_Q[2].N1 * Q_IM1.N1 + l_eigen_Q[2].N2 * Q_IM1.N2 + l_eigen_Q[2].N3 * Q_IM1.N3 + l_eigen_Q[2].N4 * Q_IM1.N4;
    Qi_iplus_IM1.N3 = l_eigen_Q[3].N0 * Q_IM1.N0 + l_eigen_Q[3].N1 * Q_IM1.N1 + l_eigen_Q[3].N2 * Q_IM1.N2 + l_eigen_Q[3].N3 * Q_IM1.N3 + l_eigen_Q[3].N4 * Q_IM1.N4;
    Qi_iplus_IM1.N4 = l_eigen_Q[4].N0 * Q_IM1.N0 + l_eigen_Q[4].N1 * Q_IM1.N1 + l_eigen_Q[4].N2 * Q_IM1.N2 + l_eigen_Q[4].N3 * Q_IM1.N3 + l_eigen_Q[4].N4 * Q_IM1.N4;
    
    Qi_iplus_IM2.N0 = l_eigen_Q[0].N0 * Q_IM2.N0 + l_eigen_Q[0].N1 * Q_IM2.N1 + l_eigen_Q[0].N2 * Q_IM2.N2 + l_eigen_Q[0].N3 * Q_IM2.N3 + l_eigen_Q[0].N4 * Q_IM2.N4;
    Qi_iplus_IM2.N1 = l_eigen_Q[1].N0 * Q_IM2.N0 + l_eigen_Q[1].N1 * Q_IM2.N1 + l_eigen_Q[1].N2 * Q_IM2.N2 + l_eigen_Q[1].N3 * Q_IM2.N3 + l_eigen_Q[1].N4 * Q_IM2.N4;
    Qi_iplus_IM2.N2 = l_eigen_Q[2].N0 * Q_IM2.N0 + l_eigen_Q[2].N1 * Q_IM2.N1 + l_eigen_Q[2].N2 * Q_IM2.N2 + l_eigen_Q[2].N3 * Q_IM2.N3 + l_eigen_Q[2].N4 * Q_IM2.N4;
    Qi_iplus_IM2.N3 = l_eigen_Q[3].N0 * Q_IM2.N0 + l_eigen_Q[3].N1 * Q_IM2.N1 + l_eigen_Q[3].N2 * Q_IM2.N2 + l_eigen_Q[3].N3 * Q_IM2.N3 + l_eigen_Q[3].N4 * Q_IM2.N4;
    Qi_iplus_IM2.N4 = l_eigen_Q[4].N0 * Q_IM2.N0 + l_eigen_Q[4].N1 * Q_IM2.N1 + l_eigen_Q[4].N2 * Q_IM2.N2 + l_eigen_Q[4].N3 * Q_IM2.N3 + l_eigen_Q[4].N4 * Q_IM2.N4;
    
    Qi_iplus_IM3.N0 = l_eigen_Q[0].N0 * Q_IM3.N0 + l_eigen_Q[0].N1 * Q_IM3.N1 + l_eigen_Q[0].N2 * Q_IM3.N2 + l_eigen_Q[0].N3 * Q_IM3.N3 + l_eigen_Q[0].N4 * Q_IM3.N4;
    Qi_iplus_IM3.N1 = l_eigen_Q[1].N0 * Q_IM3.N0 + l_eigen_Q[1].N1 * Q_IM3.N1 + l_eigen_Q[1].N2 * Q_IM3.N2 + l_eigen_Q[1].N3 * Q_IM3.N3 + l_eigen_Q[1].N4 * Q_IM3.N4;
    Qi_iplus_IM3.N2 = l_eigen_Q[2].N0 * Q_IM3.N0 + l_eigen_Q[2].N1 * Q_IM3.N1 + l_eigen_Q[2].N2 * Q_IM3.N2 + l_eigen_Q[2].N3 * Q_IM3.N3 + l_eigen_Q[2].N4 * Q_IM3.N4;
    Qi_iplus_IM3.N3 = l_eigen_Q[3].N0 * Q_IM3.N0 + l_eigen_Q[3].N1 * Q_IM3.N1 + l_eigen_Q[3].N2 * Q_IM3.N2 + l_eigen_Q[3].N3 * Q_IM3.N3 + l_eigen_Q[3].N4 * Q_IM3.N4;
    Qi_iplus_IM3.N4 = l_eigen_Q[4].N0 * Q_IM3.N0 + l_eigen_Q[4].N1 * Q_IM3.N1 + l_eigen_Q[4].N2 * Q_IM3.N2 + l_eigen_Q[4].N3 * Q_IM3.N3 + l_eigen_Q[4].N4 * Q_IM3.N4;
    

    /***********************************************************************************************/
	/****************************************** Qi_iminus********************************************/
    kk = sqrt(metric[Gid].zeta_xim*metric[Gid].zeta_xim + metric[Gid].zeta_yim * metric[Gid].zeta_yim + metric[Gid].zeta_zim * metric[Gid].zeta_zim);
    roe_a = sqrt(0.4*(ROE_AVER[Gid].roe_h_im - 0.5*(ROE_AVER[Gid].roe_u_im * ROE_AVER[Gid].roe_u_im + ROE_AVER[Gid].roe_v_im * ROE_AVER[Gid].roe_v_im + ROE_AVER[Gid].roe_w_im * ROE_AVER[Gid].roe_w_im)));

	kx_bar = metric[Gid].zeta_xim / kk;
	ky_bar = metric[Gid].zeta_yim / kk;
	kz_bar = metric[Gid].zeta_zim / kk;

	theta = kx_bar * ROE_AVER[Gid].roe_u_im + ky_bar * ROE_AVER[Gid].roe_v_im + kz_bar * ROE_AVER[Gid].roe_w_im;
	phi_sq = 0.5*0.4*(ROE_AVER[Gid].roe_u_im * ROE_AVER[Gid].roe_u_im + ROE_AVER[Gid].roe_v_im * ROE_AVER[Gid].roe_v_im + ROE_AVER[Gid].roe_w_im * ROE_AVER[Gid].roe_w_im);
	alpha = ROE_AVER[Gid].roe_rho_im / (sqrt(2.0)*roe_a);
	beta = 1.0 / (sqrt(2.0)*ROE_AVER[Gid].roe_rho_im * roe_a);

	r_eigen_Qim[0].N0 = kx_bar;
	r_eigen_Qim[0].N1 = ky_bar;
	r_eigen_Qim[0].N2 = kz_bar;
	r_eigen_Qim[0].N3 = alpha;
	r_eigen_Qim[0].N4 = alpha;

	r_eigen_Qim[1].N0 = kx_bar * ROE_AVER[Gid].roe_u_im;
	r_eigen_Qim[1].N1 = ky_bar * ROE_AVER[Gid].roe_u_im - kz_bar * ROE_AVER[Gid].roe_rho_im;
	r_eigen_Qim[1].N2 = kz_bar * ROE_AVER[Gid].roe_u_im + ky_bar * ROE_AVER[Gid].roe_rho_im;
	r_eigen_Qim[1].N3 = alpha * (ROE_AVER[Gid].roe_u_im + kx_bar * roe_a);
	r_eigen_Qim[1].N4 = alpha * (ROE_AVER[Gid].roe_u_im - kx_bar * roe_a);

	r_eigen_Qim[2].N0 = kx_bar * ROE_AVER[Gid].roe_v_im + kz_bar * ROE_AVER[Gid].roe_rho_im;
	r_eigen_Qim[2].N1 = ky_bar * ROE_AVER[Gid].roe_v_im;
	r_eigen_Qim[2].N2 = kz_bar * ROE_AVER[Gid].roe_v_im - kx_bar * ROE_AVER[Gid].roe_rho_im;
	r_eigen_Qim[2].N3 = alpha * (ROE_AVER[Gid].roe_v_im + ky_bar * roe_a);
	r_eigen_Qim[2].N4 = alpha * (ROE_AVER[Gid].roe_v_im - ky_bar * roe_a);

	r_eigen_Qim[3].N0 = kx_bar * ROE_AVER[Gid].roe_w_im - ky_bar * ROE_AVER[Gid].roe_rho_im;
	r_eigen_Qim[3].N1 = ky_bar * ROE_AVER[Gid].roe_w_im + kx_bar * ROE_AVER[Gid].roe_rho_im;
	r_eigen_Qim[3].N2 = kz_bar * ROE_AVER[Gid].roe_w_im;
	r_eigen_Qim[3].N3 = alpha * (ROE_AVER[Gid].roe_w_im + kz_bar * roe_a);
	r_eigen_Qim[3].N4 = alpha * (ROE_AVER[Gid].roe_w_im - kz_bar * roe_a);

	r_eigen_Qim[4].N0 = ((kx_bar*phi_sq) / 0.4) + ROE_AVER[Gid].roe_rho_im * (kz_bar*ROE_AVER[Gid].roe_v_im - ky_bar * ROE_AVER[Gid].roe_w_im);
	r_eigen_Qim[4].N1 = ((ky_bar*phi_sq) / 0.4) + ROE_AVER[Gid].roe_rho_im * (kx_bar*ROE_AVER[Gid].roe_w_im - kz_bar * ROE_AVER[Gid].roe_u_im);
	r_eigen_Qim[4].N2 = ((kz_bar*phi_sq) / 0.4) + ROE_AVER[Gid].roe_rho_im * (ky_bar*ROE_AVER[Gid].roe_u_im - kx_bar * ROE_AVER[Gid].roe_v_im);
	r_eigen_Qim[4].N3 = alpha * (((phi_sq + roe_a * roe_a) / (0.4)) + theta * roe_a);
	r_eigen_Qim[4].N4 = alpha * (((phi_sq + roe_a * roe_a) / (0.4)) - theta * roe_a);
	/*************************************************************************************************/

	l_eigen_Q[0].N0 = kx_bar * (1.0 - (phi_sq / (roe_a * roe_a))) - ((kz_bar*ROE_AVER[Gid].roe_v_im - ky_bar * ROE_AVER[Gid].roe_w_im) / ROE_AVER[Gid].roe_rho_im);
	l_eigen_Q[0].N1 = kx_bar * 0.4*ROE_AVER[Gid].roe_u_im / (roe_a * roe_a);
	l_eigen_Q[0].N2 = (kx_bar*0.4*ROE_AVER[Gid].roe_v_im / (roe_a * roe_a)) + (kz_bar / ROE_AVER[Gid].roe_rho_im);
	l_eigen_Q[0].N3 = (kx_bar*0.4*ROE_AVER[Gid].roe_w_im / (roe_a * roe_a)) - (ky_bar / ROE_AVER[Gid].roe_rho_im);
	l_eigen_Q[0].N4 = -kx_bar * 0.4 / (roe_a * roe_a);

	l_eigen_Q[1].N0 = ky_bar * (1.0 - (phi_sq / (roe_a * roe_a))) - ((kx_bar*ROE_AVER[Gid].roe_w_im - kz_bar * ROE_AVER[Gid].roe_u_im) / ROE_AVER[Gid].roe_rho_im);
	l_eigen_Q[1].N1 = (ky_bar*0.4*ROE_AVER[Gid].roe_u_im / (roe_a * roe_a)) - kz_bar / ROE_AVER[Gid].roe_rho_im;
	l_eigen_Q[1].N2 = ky_bar * 0.4*ROE_AVER[Gid].roe_v_im / (roe_a * roe_a);
	l_eigen_Q[1].N3 = (ky_bar*0.4*ROE_AVER[Gid].roe_w_im / (roe_a * roe_a)) + (kx_bar / ROE_AVER[Gid].roe_rho_im);
	l_eigen_Q[1].N4 = -ky_bar * 0.4 / (roe_a * roe_a);

	l_eigen_Q[2].N0 = kz_bar * (1.0 - (phi_sq / (roe_a * roe_a))) - ((ky_bar*ROE_AVER[Gid].roe_u_im - kx_bar * ROE_AVER[Gid].roe_v_im) / ROE_AVER[Gid].roe_rho_im);
	l_eigen_Q[2].N1 = kz_bar * 0.4*ROE_AVER[Gid].roe_u_im / (roe_a * roe_a) + (ky_bar / ROE_AVER[Gid].roe_rho_im);
	l_eigen_Q[2].N2 = kz_bar * 0.4*ROE_AVER[Gid].roe_v_im / (roe_a * roe_a) - (kx_bar / ROE_AVER[Gid].roe_rho_im);
	l_eigen_Q[2].N3 = kz_bar * 0.4*ROE_AVER[Gid].roe_w_im / (roe_a * roe_a);
	l_eigen_Q[2].N4 = -kz_bar * 0.4 / (roe_a * roe_a);

	l_eigen_Q[3].N0 = beta * (phi_sq - theta * roe_a);
	l_eigen_Q[3].N1 = -beta * (0.4*ROE_AVER[Gid].roe_u_im - kx_bar * roe_a);
	l_eigen_Q[3].N2 = -beta * (0.4*ROE_AVER[Gid].roe_v_im - ky_bar * roe_a);
	l_eigen_Q[3].N3 = -beta * (0.4*ROE_AVER[Gid].roe_w_im - kz_bar * roe_a);
	l_eigen_Q[3].N4 = beta * 0.4;

	l_eigen_Q[4].N0 = beta * (phi_sq + theta * roe_a);
	l_eigen_Q[4].N1 = -beta * (0.4*ROE_AVER[Gid].roe_u_im + kx_bar * roe_a);
	l_eigen_Q[4].N2 = -beta * (0.4*ROE_AVER[Gid].roe_v_im + ky_bar * roe_a);
	l_eigen_Q[4].N3 = -beta * (0.4*ROE_AVER[Gid].roe_w_im + kz_bar * roe_a);
	l_eigen_Q[4].N4 = beta * 0.4;

	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/
    Qi_iminus_IP3.N0 = l_eigen_Q[0].N0 * Q_IP3.N0 + l_eigen_Q[0].N1 * Q_IP3.N1 + l_eigen_Q[0].N2 * Q_IP3.N2 + l_eigen_Q[0].N3 * Q_IP3.N3 + l_eigen_Q[0].N4 * Q_IP3.N4;
    Qi_iminus_IP3.N1 = l_eigen_Q[1].N0 * Q_IP3.N0 + l_eigen_Q[1].N1 * Q_IP3.N1 + l_eigen_Q[1].N2 * Q_IP3.N2 + l_eigen_Q[1].N3 * Q_IP3.N3 + l_eigen_Q[1].N4 * Q_IP3.N4;
    Qi_iminus_IP3.N2 = l_eigen_Q[2].N0 * Q_IP3.N0 + l_eigen_Q[2].N1 * Q_IP3.N1 + l_eigen_Q[2].N2 * Q_IP3.N2 + l_eigen_Q[2].N3 * Q_IP3.N3 + l_eigen_Q[2].N4 * Q_IP3.N4;
    Qi_iminus_IP3.N3 = l_eigen_Q[3].N0 * Q_IP3.N0 + l_eigen_Q[3].N1 * Q_IP3.N1 + l_eigen_Q[3].N2 * Q_IP3.N2 + l_eigen_Q[3].N3 * Q_IP3.N3 + l_eigen_Q[3].N4 * Q_IP3.N4;
    Qi_iminus_IP3.N4 = l_eigen_Q[4].N0 * Q_IP3.N0 + l_eigen_Q[4].N1 * Q_IP3.N1 + l_eigen_Q[4].N2 * Q_IP3.N2 + l_eigen_Q[4].N3 * Q_IP3.N3 + l_eigen_Q[4].N4 * Q_IP3.N4;
    
    Qi_iminus_IP2.N0 = l_eigen_Q[0].N0 * Q_IP2.N0 + l_eigen_Q[0].N1 * Q_IP2.N1 + l_eigen_Q[0].N2 * Q_IP2.N2 + l_eigen_Q[0].N3 * Q_IP2.N3 + l_eigen_Q[0].N4 * Q_IP2.N4;
    Qi_iminus_IP2.N1 = l_eigen_Q[1].N0 * Q_IP2.N0 + l_eigen_Q[1].N1 * Q_IP2.N1 + l_eigen_Q[1].N2 * Q_IP2.N2 + l_eigen_Q[1].N3 * Q_IP2.N3 + l_eigen_Q[1].N4 * Q_IP2.N4;
    Qi_iminus_IP2.N2 = l_eigen_Q[2].N0 * Q_IP2.N0 + l_eigen_Q[2].N1 * Q_IP2.N1 + l_eigen_Q[2].N2 * Q_IP2.N2 + l_eigen_Q[2].N3 * Q_IP2.N3 + l_eigen_Q[2].N4 * Q_IP2.N4;
    Qi_iminus_IP2.N3 = l_eigen_Q[3].N0 * Q_IP2.N0 + l_eigen_Q[3].N1 * Q_IP2.N1 + l_eigen_Q[3].N2 * Q_IP2.N2 + l_eigen_Q[3].N3 * Q_IP2.N3 + l_eigen_Q[3].N4 * Q_IP2.N4;
    Qi_iminus_IP2.N4 = l_eigen_Q[4].N0 * Q_IP2.N0 + l_eigen_Q[4].N1 * Q_IP2.N1 + l_eigen_Q[4].N2 * Q_IP2.N2 + l_eigen_Q[4].N3 * Q_IP2.N3 + l_eigen_Q[4].N4 * Q_IP2.N4;
    
    Qi_iminus_IP1.N0 = l_eigen_Q[0].N0 * Q_IP1.N0 + l_eigen_Q[0].N1 * Q_IP1.N1 + l_eigen_Q[0].N2 * Q_IP1.N2 + l_eigen_Q[0].N3 * Q_IP1.N3 + l_eigen_Q[0].N4 * Q_IP1.N4;
    Qi_iminus_IP1.N1 = l_eigen_Q[1].N0 * Q_IP1.N0 + l_eigen_Q[1].N1 * Q_IP1.N1 + l_eigen_Q[1].N2 * Q_IP1.N2 + l_eigen_Q[1].N3 * Q_IP1.N3 + l_eigen_Q[1].N4 * Q_IP1.N4;
    Qi_iminus_IP1.N2 = l_eigen_Q[2].N0 * Q_IP1.N0 + l_eigen_Q[2].N1 * Q_IP1.N1 + l_eigen_Q[2].N2 * Q_IP1.N2 + l_eigen_Q[2].N3 * Q_IP1.N3 + l_eigen_Q[2].N4 * Q_IP1.N4;
    Qi_iminus_IP1.N3 = l_eigen_Q[3].N0 * Q_IP1.N0 + l_eigen_Q[3].N1 * Q_IP1.N1 + l_eigen_Q[3].N2 * Q_IP1.N2 + l_eigen_Q[3].N3 * Q_IP1.N3 + l_eigen_Q[3].N4 * Q_IP1.N4;
    Qi_iminus_IP1.N4 = l_eigen_Q[4].N0 * Q_IP1.N0 + l_eigen_Q[4].N1 * Q_IP1.N1 + l_eigen_Q[4].N2 * Q_IP1.N2 + l_eigen_Q[4].N3 * Q_IP1.N3 + l_eigen_Q[4].N4 * Q_IP1.N4;
    
    Qi_iminus_I.N0 = l_eigen_Q[0].N0 * Q_I.N0 + l_eigen_Q[0].N1 * Q_I.N1 + l_eigen_Q[0].N2 * Q_I.N2 + l_eigen_Q[0].N3 * Q_I.N3 + l_eigen_Q[0].N4 * Q_I.N4;
    Qi_iminus_I.N1 = l_eigen_Q[1].N0 * Q_I.N0 + l_eigen_Q[1].N1 * Q_I.N1 + l_eigen_Q[1].N2 * Q_I.N2 + l_eigen_Q[1].N3 * Q_I.N3 + l_eigen_Q[1].N4 * Q_I.N4;
    Qi_iminus_I.N2 = l_eigen_Q[2].N0 * Q_I.N0 + l_eigen_Q[2].N1 * Q_I.N1 + l_eigen_Q[2].N2 * Q_I.N2 + l_eigen_Q[2].N3 * Q_I.N3 + l_eigen_Q[2].N4 * Q_I.N4;
    Qi_iminus_I.N3 = l_eigen_Q[3].N0 * Q_I.N0 + l_eigen_Q[3].N1 * Q_I.N1 + l_eigen_Q[3].N2 * Q_I.N2 + l_eigen_Q[3].N3 * Q_I.N3 + l_eigen_Q[3].N4 * Q_I.N4;
    Qi_iminus_I.N4 = l_eigen_Q[4].N0 * Q_I.N0 + l_eigen_Q[4].N1 * Q_I.N1 + l_eigen_Q[4].N2 * Q_I.N2 + l_eigen_Q[4].N3 * Q_I.N3 + l_eigen_Q[4].N4 * Q_I.N4;
    
    Qi_iminus_IM1.N0 = l_eigen_Q[0].N0 * Q_IM1.N0 + l_eigen_Q[0].N1 * Q_IM1.N1 + l_eigen_Q[0].N2 * Q_IM1.N2 + l_eigen_Q[0].N3 * Q_IM1.N3 + l_eigen_Q[0].N4 * Q_IM1.N4;
    Qi_iminus_IM1.N1 = l_eigen_Q[1].N0 * Q_IM1.N0 + l_eigen_Q[1].N1 * Q_IM1.N1 + l_eigen_Q[1].N2 * Q_IM1.N2 + l_eigen_Q[1].N3 * Q_IM1.N3 + l_eigen_Q[1].N4 * Q_IM1.N4;
    Qi_iminus_IM1.N2 = l_eigen_Q[2].N0 * Q_IM1.N0 + l_eigen_Q[2].N1 * Q_IM1.N1 + l_eigen_Q[2].N2 * Q_IM1.N2 + l_eigen_Q[2].N3 * Q_IM1.N3 + l_eigen_Q[2].N4 * Q_IM1.N4;
    Qi_iminus_IM1.N3 = l_eigen_Q[3].N0 * Q_IM1.N0 + l_eigen_Q[3].N1 * Q_IM1.N1 + l_eigen_Q[3].N2 * Q_IM1.N2 + l_eigen_Q[3].N3 * Q_IM1.N3 + l_eigen_Q[3].N4 * Q_IM1.N4;
    Qi_iminus_IM1.N4 = l_eigen_Q[4].N0 * Q_IM1.N0 + l_eigen_Q[4].N1 * Q_IM1.N1 + l_eigen_Q[4].N2 * Q_IM1.N2 + l_eigen_Q[4].N3 * Q_IM1.N3 + l_eigen_Q[4].N4 * Q_IM1.N4;
    
    Qi_iminus_IM2.N0 = l_eigen_Q[0].N0 * Q_IM2.N0 + l_eigen_Q[0].N1 * Q_IM2.N1 + l_eigen_Q[0].N2 * Q_IM2.N2 + l_eigen_Q[0].N3 * Q_IM2.N3 + l_eigen_Q[0].N4 * Q_IM2.N4;
    Qi_iminus_IM2.N1 = l_eigen_Q[1].N0 * Q_IM2.N0 + l_eigen_Q[1].N1 * Q_IM2.N1 + l_eigen_Q[1].N2 * Q_IM2.N2 + l_eigen_Q[1].N3 * Q_IM2.N3 + l_eigen_Q[1].N4 * Q_IM2.N4;
    Qi_iminus_IM2.N2 = l_eigen_Q[2].N0 * Q_IM2.N0 + l_eigen_Q[2].N1 * Q_IM2.N1 + l_eigen_Q[2].N2 * Q_IM2.N2 + l_eigen_Q[2].N3 * Q_IM2.N3 + l_eigen_Q[2].N4 * Q_IM2.N4;
    Qi_iminus_IM2.N3 = l_eigen_Q[3].N0 * Q_IM2.N0 + l_eigen_Q[3].N1 * Q_IM2.N1 + l_eigen_Q[3].N2 * Q_IM2.N2 + l_eigen_Q[3].N3 * Q_IM2.N3 + l_eigen_Q[3].N4 * Q_IM2.N4;
    Qi_iminus_IM2.N4 = l_eigen_Q[4].N0 * Q_IM2.N0 + l_eigen_Q[4].N1 * Q_IM2.N1 + l_eigen_Q[4].N2 * Q_IM2.N2 + l_eigen_Q[4].N3 * Q_IM2.N3 + l_eigen_Q[4].N4 * Q_IM2.N4;
    
    Qi_iminus_IM3.N0 = l_eigen_Q[0].N0 * Q_IM3.N0 + l_eigen_Q[0].N1 * Q_IM3.N1 + l_eigen_Q[0].N2 * Q_IM3.N2 + l_eigen_Q[0].N3 * Q_IM3.N3 + l_eigen_Q[0].N4 * Q_IM3.N4;
    Qi_iminus_IM3.N1 = l_eigen_Q[1].N0 * Q_IM3.N0 + l_eigen_Q[1].N1 * Q_IM3.N1 + l_eigen_Q[1].N2 * Q_IM3.N2 + l_eigen_Q[1].N3 * Q_IM3.N3 + l_eigen_Q[1].N4 * Q_IM3.N4;
    Qi_iminus_IM3.N2 = l_eigen_Q[2].N0 * Q_IM3.N0 + l_eigen_Q[2].N1 * Q_IM3.N1 + l_eigen_Q[2].N2 * Q_IM3.N2 + l_eigen_Q[2].N3 * Q_IM3.N3 + l_eigen_Q[2].N4 * Q_IM3.N4;
    Qi_iminus_IM3.N3 = l_eigen_Q[3].N0 * Q_IM3.N0 + l_eigen_Q[3].N1 * Q_IM3.N1 + l_eigen_Q[3].N2 * Q_IM3.N2 + l_eigen_Q[3].N3 * Q_IM3.N3 + l_eigen_Q[3].N4 * Q_IM3.N4;
    Qi_iminus_IM3.N4 = l_eigen_Q[4].N0 * Q_IM3.N0 + l_eigen_Q[4].N1 * Q_IM3.N1 + l_eigen_Q[4].N2 * Q_IM3.N2 + l_eigen_Q[4].N3 * Q_IM3.N3 + l_eigen_Q[4].N4 * Q_IM3.N4;
    


    /***********************************************************************************************/
	/****************************************** Qj_iplus********************************************/
    kk = sqrt(metric[Gid].eta_xjp*metric[Gid].eta_xjp + metric[Gid].eta_yjp * metric[Gid].eta_yjp + metric[Gid].eta_zjp * metric[Gid].eta_zjp);
    roe_a = sqrt(0.4*(ROE_AVER[Gid].roe_h_jp - 0.5*(ROE_AVER[Gid].roe_u_jp * ROE_AVER[Gid].roe_u_jp + ROE_AVER[Gid].roe_v_jp * ROE_AVER[Gid].roe_v_jp + ROE_AVER[Gid].roe_w_jp * ROE_AVER[Gid].roe_w_jp)));

	kx_bar = metric[Gid].eta_xjp / kk;
	ky_bar = metric[Gid].eta_yjp / kk;
	kz_bar = metric[Gid].eta_zjp / kk;

	theta = kx_bar * ROE_AVER[Gid].roe_u_jp + ky_bar * ROE_AVER[Gid].roe_v_jp + kz_bar * ROE_AVER[Gid].roe_w_jp;
	phi_sq = 0.5*0.4*(ROE_AVER[Gid].roe_u_jp * ROE_AVER[Gid].roe_u_jp + ROE_AVER[Gid].roe_v_jp * ROE_AVER[Gid].roe_v_jp + ROE_AVER[Gid].roe_w_jp * ROE_AVER[Gid].roe_w_jp);
	alpha = ROE_AVER[Gid].roe_rho_jp / (sqrt(2.0)*roe_a);
	beta = 1.0 / (sqrt(2.0)*ROE_AVER[Gid].roe_rho_jp * roe_a);

	r_eigen_Qjp[0].N0 = kx_bar;
	r_eigen_Qjp[0].N1 = ky_bar;
	r_eigen_Qjp[0].N2 = kz_bar;
	r_eigen_Qjp[0].N3 = alpha;
	r_eigen_Qjp[0].N4 = alpha;

	r_eigen_Qjp[1].N0 = kx_bar * ROE_AVER[Gid].roe_u_jp;
	r_eigen_Qjp[1].N1 = ky_bar * ROE_AVER[Gid].roe_u_jp - kz_bar * ROE_AVER[Gid].roe_rho_jp;
	r_eigen_Qjp[1].N2 = kz_bar * ROE_AVER[Gid].roe_u_jp + ky_bar * ROE_AVER[Gid].roe_rho_jp;
	r_eigen_Qjp[1].N3 = alpha * (ROE_AVER[Gid].roe_u_jp + kx_bar * roe_a);
	r_eigen_Qjp[1].N4 = alpha * (ROE_AVER[Gid].roe_u_jp - kx_bar * roe_a);

	r_eigen_Qjp[2].N0 = kx_bar * ROE_AVER[Gid].roe_v_jp + kz_bar * ROE_AVER[Gid].roe_rho_jp;
	r_eigen_Qjp[2].N1 = ky_bar * ROE_AVER[Gid].roe_v_jp;
	r_eigen_Qjp[2].N2 = kz_bar * ROE_AVER[Gid].roe_v_jp - kx_bar * ROE_AVER[Gid].roe_rho_jp;
	r_eigen_Qjp[2].N3 = alpha * (ROE_AVER[Gid].roe_v_jp + ky_bar * roe_a);
	r_eigen_Qjp[2].N4 = alpha * (ROE_AVER[Gid].roe_v_jp - ky_bar * roe_a);

	r_eigen_Qjp[3].N0 = kx_bar * ROE_AVER[Gid].roe_w_jp - ky_bar * ROE_AVER[Gid].roe_rho_jp;
	r_eigen_Qjp[3].N1 = ky_bar * ROE_AVER[Gid].roe_w_jp + kx_bar * ROE_AVER[Gid].roe_rho_jp;
	r_eigen_Qjp[3].N2 = kz_bar * ROE_AVER[Gid].roe_w_jp;
	r_eigen_Qjp[3].N3 = alpha * (ROE_AVER[Gid].roe_w_jp + kz_bar * roe_a);
	r_eigen_Qjp[3].N4 = alpha * (ROE_AVER[Gid].roe_w_jp - kz_bar * roe_a);

	r_eigen_Qjp[4].N0 = ((kx_bar*phi_sq) / 0.4) + ROE_AVER[Gid].roe_rho_jp * (kz_bar*ROE_AVER[Gid].roe_v_jp - ky_bar * ROE_AVER[Gid].roe_w_jp);
	r_eigen_Qjp[4].N1 = ((ky_bar*phi_sq) / 0.4) + ROE_AVER[Gid].roe_rho_jp * (kx_bar*ROE_AVER[Gid].roe_w_jp - kz_bar * ROE_AVER[Gid].roe_u_jp);
	r_eigen_Qjp[4].N2 = ((kz_bar*phi_sq) / 0.4) + ROE_AVER[Gid].roe_rho_jp * (ky_bar*ROE_AVER[Gid].roe_u_jp - kx_bar * ROE_AVER[Gid].roe_v_jp);
	r_eigen_Qjp[4].N3 = alpha * (((phi_sq + roe_a * roe_a) / (0.4)) + theta * roe_a);
	r_eigen_Qjp[4].N4 = alpha * (((phi_sq + roe_a * roe_a) / (0.4)) - theta * roe_a);
	/*************************************************************************************************/

	l_eigen_Q[0].N0 = kx_bar * (1.0 - (phi_sq / (roe_a * roe_a))) - ((kz_bar*ROE_AVER[Gid].roe_v_jp - ky_bar * ROE_AVER[Gid].roe_w_jp) / ROE_AVER[Gid].roe_rho_jp);
	l_eigen_Q[0].N1 = kx_bar * 0.4*ROE_AVER[Gid].roe_u_jp / (roe_a * roe_a);
	l_eigen_Q[0].N2 = (kx_bar*0.4*ROE_AVER[Gid].roe_v_jp / (roe_a * roe_a)) + (kz_bar / ROE_AVER[Gid].roe_rho_jp);
	l_eigen_Q[0].N3 = (kx_bar*0.4*ROE_AVER[Gid].roe_w_jp / (roe_a * roe_a)) - (ky_bar / ROE_AVER[Gid].roe_rho_jp);
	l_eigen_Q[0].N4 = -kx_bar * 0.4 / (roe_a * roe_a);

	l_eigen_Q[1].N0 = ky_bar * (1.0 - (phi_sq / (roe_a * roe_a))) - ((kx_bar*ROE_AVER[Gid].roe_w_jp - kz_bar * ROE_AVER[Gid].roe_u_jp) / ROE_AVER[Gid].roe_rho_jp);
	l_eigen_Q[1].N1 = (ky_bar*0.4*ROE_AVER[Gid].roe_u_jp / (roe_a * roe_a)) - kz_bar / ROE_AVER[Gid].roe_rho_jp;
	l_eigen_Q[1].N2 = ky_bar * 0.4*ROE_AVER[Gid].roe_v_jp / (roe_a * roe_a);
	l_eigen_Q[1].N3 = (ky_bar*0.4*ROE_AVER[Gid].roe_w_jp / (roe_a * roe_a)) + (kx_bar / ROE_AVER[Gid].roe_rho_jp);
	l_eigen_Q[1].N4 = -ky_bar * 0.4 / (roe_a * roe_a);

	l_eigen_Q[2].N0 = kz_bar * (1.0 - (phi_sq / (roe_a * roe_a))) - ((ky_bar*ROE_AVER[Gid].roe_u_jp - kx_bar * ROE_AVER[Gid].roe_v_jp) / ROE_AVER[Gid].roe_rho_jp);
	l_eigen_Q[2].N1 = kz_bar * 0.4*ROE_AVER[Gid].roe_u_jp / (roe_a * roe_a) + (ky_bar / ROE_AVER[Gid].roe_rho_jp);
	l_eigen_Q[2].N2 = kz_bar * 0.4*ROE_AVER[Gid].roe_v_jp / (roe_a * roe_a) - (kx_bar / ROE_AVER[Gid].roe_rho_jp);
	l_eigen_Q[2].N3 = kz_bar * 0.4*ROE_AVER[Gid].roe_w_jp / (roe_a * roe_a);
	l_eigen_Q[2].N4 = -kz_bar * 0.4 / (roe_a * roe_a);

	l_eigen_Q[3].N0 = beta * (phi_sq - theta * roe_a);
	l_eigen_Q[3].N1 = -beta * (0.4*ROE_AVER[Gid].roe_u_jp - kx_bar * roe_a);
	l_eigen_Q[3].N2 = -beta * (0.4*ROE_AVER[Gid].roe_v_jp - ky_bar * roe_a);
	l_eigen_Q[3].N3 = -beta * (0.4*ROE_AVER[Gid].roe_w_jp - kz_bar * roe_a);
	l_eigen_Q[3].N4 = beta * 0.4;

	l_eigen_Q[4].N0 = beta * (phi_sq + theta * roe_a);
	l_eigen_Q[4].N1 = -beta * (0.4*ROE_AVER[Gid].roe_u_jp + kx_bar * roe_a);
	l_eigen_Q[4].N2 = -beta * (0.4*ROE_AVER[Gid].roe_v_jp + ky_bar * roe_a);
	l_eigen_Q[4].N3 = -beta * (0.4*ROE_AVER[Gid].roe_w_jp + kz_bar * roe_a);
	l_eigen_Q[4].N4 = beta * 0.4;

	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/
    Qj_iplus_JP3.N0 = l_eigen_Q[0].N0 * Q_JP3.N0 + l_eigen_Q[0].N1 * Q_JP3.N1 + l_eigen_Q[0].N2 * Q_JP3.N2 + l_eigen_Q[0].N3 * Q_JP3.N3 + l_eigen_Q[0].N4 * Q_JP3.N4;
    Qj_iplus_JP3.N1 = l_eigen_Q[1].N0 * Q_JP3.N0 + l_eigen_Q[1].N1 * Q_JP3.N1 + l_eigen_Q[1].N2 * Q_JP3.N2 + l_eigen_Q[1].N3 * Q_JP3.N3 + l_eigen_Q[1].N4 * Q_JP3.N4;
    Qj_iplus_JP3.N2 = l_eigen_Q[2].N0 * Q_JP3.N0 + l_eigen_Q[2].N1 * Q_JP3.N1 + l_eigen_Q[2].N2 * Q_JP3.N2 + l_eigen_Q[2].N3 * Q_JP3.N3 + l_eigen_Q[2].N4 * Q_JP3.N4;
    Qj_iplus_JP3.N3 = l_eigen_Q[3].N0 * Q_JP3.N0 + l_eigen_Q[3].N1 * Q_JP3.N1 + l_eigen_Q[3].N2 * Q_JP3.N2 + l_eigen_Q[3].N3 * Q_JP3.N3 + l_eigen_Q[3].N4 * Q_JP3.N4;
    Qj_iplus_JP3.N4 = l_eigen_Q[4].N0 * Q_JP3.N0 + l_eigen_Q[4].N1 * Q_JP3.N1 + l_eigen_Q[4].N2 * Q_JP3.N2 + l_eigen_Q[4].N3 * Q_JP3.N3 + l_eigen_Q[4].N4 * Q_JP3.N4;
    
    Qj_iplus_JP2.N0 = l_eigen_Q[0].N0 * Q_JP2.N0 + l_eigen_Q[0].N1 * Q_JP2.N1 + l_eigen_Q[0].N2 * Q_JP2.N2 + l_eigen_Q[0].N3 * Q_JP2.N3 + l_eigen_Q[0].N4 * Q_JP2.N4;
    Qj_iplus_JP2.N1 = l_eigen_Q[1].N0 * Q_JP2.N0 + l_eigen_Q[1].N1 * Q_JP2.N1 + l_eigen_Q[1].N2 * Q_JP2.N2 + l_eigen_Q[1].N3 * Q_JP2.N3 + l_eigen_Q[1].N4 * Q_JP2.N4;
    Qj_iplus_JP2.N2 = l_eigen_Q[2].N0 * Q_JP2.N0 + l_eigen_Q[2].N1 * Q_JP2.N1 + l_eigen_Q[2].N2 * Q_JP2.N2 + l_eigen_Q[2].N3 * Q_JP2.N3 + l_eigen_Q[2].N4 * Q_JP2.N4;
    Qj_iplus_JP2.N3 = l_eigen_Q[3].N0 * Q_JP2.N0 + l_eigen_Q[3].N1 * Q_JP2.N1 + l_eigen_Q[3].N2 * Q_JP2.N2 + l_eigen_Q[3].N3 * Q_JP2.N3 + l_eigen_Q[3].N4 * Q_JP2.N4;
    Qj_iplus_JP2.N4 = l_eigen_Q[4].N0 * Q_JP2.N0 + l_eigen_Q[4].N1 * Q_JP2.N1 + l_eigen_Q[4].N2 * Q_JP2.N2 + l_eigen_Q[4].N3 * Q_JP2.N3 + l_eigen_Q[4].N4 * Q_JP2.N4;
    
    Qj_iplus_JP1.N0 = l_eigen_Q[0].N0 * Q_JP1.N0 + l_eigen_Q[0].N1 * Q_JP1.N1 + l_eigen_Q[0].N2 * Q_JP1.N2 + l_eigen_Q[0].N3 * Q_JP1.N3 + l_eigen_Q[0].N4 * Q_JP1.N4;
    Qj_iplus_JP1.N1 = l_eigen_Q[1].N0 * Q_JP1.N0 + l_eigen_Q[1].N1 * Q_JP1.N1 + l_eigen_Q[1].N2 * Q_JP1.N2 + l_eigen_Q[1].N3 * Q_JP1.N3 + l_eigen_Q[1].N4 * Q_JP1.N4;
    Qj_iplus_JP1.N2 = l_eigen_Q[2].N0 * Q_JP1.N0 + l_eigen_Q[2].N1 * Q_JP1.N1 + l_eigen_Q[2].N2 * Q_JP1.N2 + l_eigen_Q[2].N3 * Q_JP1.N3 + l_eigen_Q[2].N4 * Q_JP1.N4;
    Qj_iplus_JP1.N3 = l_eigen_Q[3].N0 * Q_JP1.N0 + l_eigen_Q[3].N1 * Q_JP1.N1 + l_eigen_Q[3].N2 * Q_JP1.N2 + l_eigen_Q[3].N3 * Q_JP1.N3 + l_eigen_Q[3].N4 * Q_JP1.N4;
    Qj_iplus_JP1.N4 = l_eigen_Q[4].N0 * Q_JP1.N0 + l_eigen_Q[4].N1 * Q_JP1.N1 + l_eigen_Q[4].N2 * Q_JP1.N2 + l_eigen_Q[4].N3 * Q_JP1.N3 + l_eigen_Q[4].N4 * Q_JP1.N4;
    
    Qj_iplus_I.N0 = l_eigen_Q[0].N0 * Q_I.N0 + l_eigen_Q[0].N1 * Q_I.N1 + l_eigen_Q[0].N2 * Q_I.N2 + l_eigen_Q[0].N3 * Q_I.N3 + l_eigen_Q[0].N4 * Q_I.N4;
    Qj_iplus_I.N1 = l_eigen_Q[1].N0 * Q_I.N0 + l_eigen_Q[1].N1 * Q_I.N1 + l_eigen_Q[1].N2 * Q_I.N2 + l_eigen_Q[1].N3 * Q_I.N3 + l_eigen_Q[1].N4 * Q_I.N4;
    Qj_iplus_I.N2 = l_eigen_Q[2].N0 * Q_I.N0 + l_eigen_Q[2].N1 * Q_I.N1 + l_eigen_Q[2].N2 * Q_I.N2 + l_eigen_Q[2].N3 * Q_I.N3 + l_eigen_Q[2].N4 * Q_I.N4;
    Qj_iplus_I.N3 = l_eigen_Q[3].N0 * Q_I.N0 + l_eigen_Q[3].N1 * Q_I.N1 + l_eigen_Q[3].N2 * Q_I.N2 + l_eigen_Q[3].N3 * Q_I.N3 + l_eigen_Q[3].N4 * Q_I.N4;
    Qj_iplus_I.N4 = l_eigen_Q[4].N0 * Q_I.N0 + l_eigen_Q[4].N1 * Q_I.N1 + l_eigen_Q[4].N2 * Q_I.N2 + l_eigen_Q[4].N3 * Q_I.N3 + l_eigen_Q[4].N4 * Q_I.N4;
    
    Qj_iplus_JM1.N0 = l_eigen_Q[0].N0 * Q_JM1.N0 + l_eigen_Q[0].N1 * Q_JM1.N1 + l_eigen_Q[0].N2 * Q_JM1.N2 + l_eigen_Q[0].N3 * Q_JM1.N3 + l_eigen_Q[0].N4 * Q_JM1.N4;
    Qj_iplus_JM1.N1 = l_eigen_Q[1].N0 * Q_JM1.N0 + l_eigen_Q[1].N1 * Q_JM1.N1 + l_eigen_Q[1].N2 * Q_JM1.N2 + l_eigen_Q[1].N3 * Q_JM1.N3 + l_eigen_Q[1].N4 * Q_JM1.N4;
    Qj_iplus_JM1.N2 = l_eigen_Q[2].N0 * Q_JM1.N0 + l_eigen_Q[2].N1 * Q_JM1.N1 + l_eigen_Q[2].N2 * Q_JM1.N2 + l_eigen_Q[2].N3 * Q_JM1.N3 + l_eigen_Q[2].N4 * Q_JM1.N4;
    Qj_iplus_JM1.N3 = l_eigen_Q[3].N0 * Q_JM1.N0 + l_eigen_Q[3].N1 * Q_JM1.N1 + l_eigen_Q[3].N2 * Q_JM1.N2 + l_eigen_Q[3].N3 * Q_JM1.N3 + l_eigen_Q[3].N4 * Q_JM1.N4;
    Qj_iplus_JM1.N4 = l_eigen_Q[4].N0 * Q_JM1.N0 + l_eigen_Q[4].N1 * Q_JM1.N1 + l_eigen_Q[4].N2 * Q_JM1.N2 + l_eigen_Q[4].N3 * Q_JM1.N3 + l_eigen_Q[4].N4 * Q_JM1.N4;
    
    Qj_iplus_JM2.N0 = l_eigen_Q[0].N0 * Q_JM2.N0 + l_eigen_Q[0].N1 * Q_JM2.N1 + l_eigen_Q[0].N2 * Q_JM2.N2 + l_eigen_Q[0].N3 * Q_JM2.N3 + l_eigen_Q[0].N4 * Q_JM2.N4;
    Qj_iplus_JM2.N1 = l_eigen_Q[1].N0 * Q_JM2.N0 + l_eigen_Q[1].N1 * Q_JM2.N1 + l_eigen_Q[1].N2 * Q_JM2.N2 + l_eigen_Q[1].N3 * Q_JM2.N3 + l_eigen_Q[1].N4 * Q_JM2.N4;
    Qj_iplus_JM2.N2 = l_eigen_Q[2].N0 * Q_JM2.N0 + l_eigen_Q[2].N1 * Q_JM2.N1 + l_eigen_Q[2].N2 * Q_JM2.N2 + l_eigen_Q[2].N3 * Q_JM2.N3 + l_eigen_Q[2].N4 * Q_JM2.N4;
    Qj_iplus_JM2.N3 = l_eigen_Q[3].N0 * Q_JM2.N0 + l_eigen_Q[3].N1 * Q_JM2.N1 + l_eigen_Q[3].N2 * Q_JM2.N2 + l_eigen_Q[3].N3 * Q_JM2.N3 + l_eigen_Q[3].N4 * Q_JM2.N4;
    Qj_iplus_JM2.N4 = l_eigen_Q[4].N0 * Q_JM2.N0 + l_eigen_Q[4].N1 * Q_JM2.N1 + l_eigen_Q[4].N2 * Q_JM2.N2 + l_eigen_Q[4].N3 * Q_JM2.N3 + l_eigen_Q[4].N4 * Q_JM2.N4;
    
    Qj_iplus_JM3.N0 = l_eigen_Q[0].N0 * Q_JM3.N0 + l_eigen_Q[0].N1 * Q_JM3.N1 + l_eigen_Q[0].N2 * Q_JM3.N2 + l_eigen_Q[0].N3 * Q_JM3.N3 + l_eigen_Q[0].N4 * Q_JM3.N4;
    Qj_iplus_JM3.N1 = l_eigen_Q[1].N0 * Q_JM3.N0 + l_eigen_Q[1].N1 * Q_JM3.N1 + l_eigen_Q[1].N2 * Q_JM3.N2 + l_eigen_Q[1].N3 * Q_JM3.N3 + l_eigen_Q[1].N4 * Q_JM3.N4;
    Qj_iplus_JM3.N2 = l_eigen_Q[2].N0 * Q_JM3.N0 + l_eigen_Q[2].N1 * Q_JM3.N1 + l_eigen_Q[2].N2 * Q_JM3.N2 + l_eigen_Q[2].N3 * Q_JM3.N3 + l_eigen_Q[2].N4 * Q_JM3.N4;
    Qj_iplus_JM3.N3 = l_eigen_Q[3].N0 * Q_JM3.N0 + l_eigen_Q[3].N1 * Q_JM3.N1 + l_eigen_Q[3].N2 * Q_JM3.N2 + l_eigen_Q[3].N3 * Q_JM3.N3 + l_eigen_Q[3].N4 * Q_JM3.N4;
    Qj_iplus_JM3.N4 = l_eigen_Q[4].N0 * Q_JM3.N0 + l_eigen_Q[4].N1 * Q_JM3.N1 + l_eigen_Q[4].N2 * Q_JM3.N2 + l_eigen_Q[4].N3 * Q_JM3.N3 + l_eigen_Q[4].N4 * Q_JM3.N4;
    

    /***********************************************************************************************/
	/****************************************** Qj_iminus********************************************/
    kk = sqrt(metric[Gid].eta_xjm*metric[Gid].eta_xjm + metric[Gid].eta_yjm * metric[Gid].eta_yjm + metric[Gid].eta_zjm * metric[Gid].eta_zjm);
    roe_a = sqrt(0.4*(ROE_AVER[Gid].roe_h_jm - 0.5*(ROE_AVER[Gid].roe_u_jm * ROE_AVER[Gid].roe_u_jm + ROE_AVER[Gid].roe_v_jm * ROE_AVER[Gid].roe_v_jm + ROE_AVER[Gid].roe_w_jm * ROE_AVER[Gid].roe_w_jm)));

	kx_bar = metric[Gid].eta_xjm / kk;
	ky_bar = metric[Gid].eta_yjm / kk;
	kz_bar = metric[Gid].eta_zjm / kk;

	theta = kx_bar * ROE_AVER[Gid].roe_u_jm + ky_bar * ROE_AVER[Gid].roe_v_jm + kz_bar * ROE_AVER[Gid].roe_w_jm;
	phi_sq = 0.5*0.4*(ROE_AVER[Gid].roe_u_jm * ROE_AVER[Gid].roe_u_jm + ROE_AVER[Gid].roe_v_jm * ROE_AVER[Gid].roe_v_jm + ROE_AVER[Gid].roe_w_jm * ROE_AVER[Gid].roe_w_jm);
	alpha = ROE_AVER[Gid].roe_rho_jm / (sqrt(2.0)*roe_a);
	beta = 1.0 / (sqrt(2.0)*ROE_AVER[Gid].roe_rho_jm * roe_a);

	r_eigen_Qjm[0].N0 = kx_bar;
	r_eigen_Qjm[0].N1 = ky_bar;
	r_eigen_Qjm[0].N2 = kz_bar;
	r_eigen_Qjm[0].N3 = alpha;
	r_eigen_Qjm[0].N4 = alpha;

	r_eigen_Qjm[1].N0 = kx_bar * ROE_AVER[Gid].roe_u_jm;
	r_eigen_Qjm[1].N1 = ky_bar * ROE_AVER[Gid].roe_u_jm - kz_bar * ROE_AVER[Gid].roe_rho_jm;
	r_eigen_Qjm[1].N2 = kz_bar * ROE_AVER[Gid].roe_u_jm + ky_bar * ROE_AVER[Gid].roe_rho_jm;
	r_eigen_Qjm[1].N3 = alpha * (ROE_AVER[Gid].roe_u_jm + kx_bar * roe_a);
	r_eigen_Qjm[1].N4 = alpha * (ROE_AVER[Gid].roe_u_jm - kx_bar * roe_a);

	r_eigen_Qjm[2].N0 = kx_bar * ROE_AVER[Gid].roe_v_jm + kz_bar * ROE_AVER[Gid].roe_rho_jm;
	r_eigen_Qjm[2].N1 = ky_bar * ROE_AVER[Gid].roe_v_jm;
	r_eigen_Qjm[2].N2 = kz_bar * ROE_AVER[Gid].roe_v_jm - kx_bar * ROE_AVER[Gid].roe_rho_jm;
	r_eigen_Qjm[2].N3 = alpha * (ROE_AVER[Gid].roe_v_jm + ky_bar * roe_a);
	r_eigen_Qjm[2].N4 = alpha * (ROE_AVER[Gid].roe_v_jm - ky_bar * roe_a);

	r_eigen_Qjm[3].N0 = kx_bar * ROE_AVER[Gid].roe_w_jm - ky_bar * ROE_AVER[Gid].roe_rho_jm;
	r_eigen_Qjm[3].N1 = ky_bar * ROE_AVER[Gid].roe_w_jm + kx_bar * ROE_AVER[Gid].roe_rho_jm;
	r_eigen_Qjm[3].N2 = kz_bar * ROE_AVER[Gid].roe_w_jm;
	r_eigen_Qjm[3].N3 = alpha * (ROE_AVER[Gid].roe_w_jm + kz_bar * roe_a);
	r_eigen_Qjm[3].N4 = alpha * (ROE_AVER[Gid].roe_w_jm - kz_bar * roe_a);

	r_eigen_Qjm[4].N0 = ((kx_bar*phi_sq) / 0.4) + ROE_AVER[Gid].roe_rho_jm * (kz_bar*ROE_AVER[Gid].roe_v_jm - ky_bar * ROE_AVER[Gid].roe_w_jm);
	r_eigen_Qjm[4].N1 = ((ky_bar*phi_sq) / 0.4) + ROE_AVER[Gid].roe_rho_jm * (kx_bar*ROE_AVER[Gid].roe_w_jm - kz_bar * ROE_AVER[Gid].roe_u_jm);
	r_eigen_Qjm[4].N2 = ((kz_bar*phi_sq) / 0.4) + ROE_AVER[Gid].roe_rho_jm * (ky_bar*ROE_AVER[Gid].roe_u_jm - kx_bar * ROE_AVER[Gid].roe_v_jm);
	r_eigen_Qjm[4].N3 = alpha * (((phi_sq + roe_a * roe_a) / (0.4)) + theta * roe_a);
	r_eigen_Qjm[4].N4 = alpha * (((phi_sq + roe_a * roe_a) / (0.4)) - theta * roe_a);
	/*************************************************************************************************/

	l_eigen_Q[0].N0 = kx_bar * (1.0 - (phi_sq / (roe_a * roe_a))) - ((kz_bar*ROE_AVER[Gid].roe_v_jm - ky_bar * ROE_AVER[Gid].roe_w_jm) / ROE_AVER[Gid].roe_rho_jm);
	l_eigen_Q[0].N1 = kx_bar * 0.4*ROE_AVER[Gid].roe_u_jm / (roe_a * roe_a);
	l_eigen_Q[0].N2 = (kx_bar*0.4*ROE_AVER[Gid].roe_v_jm / (roe_a * roe_a)) + (kz_bar / ROE_AVER[Gid].roe_rho_jm);
	l_eigen_Q[0].N3 = (kx_bar*0.4*ROE_AVER[Gid].roe_w_jm / (roe_a * roe_a)) - (ky_bar / ROE_AVER[Gid].roe_rho_jm);
	l_eigen_Q[0].N4 = -kx_bar * 0.4 / (roe_a * roe_a);

	l_eigen_Q[1].N0 = ky_bar * (1.0 - (phi_sq / (roe_a * roe_a))) - ((kx_bar*ROE_AVER[Gid].roe_w_jm - kz_bar * ROE_AVER[Gid].roe_u_jm) / ROE_AVER[Gid].roe_rho_jm);
	l_eigen_Q[1].N1 = (ky_bar*0.4*ROE_AVER[Gid].roe_u_jm / (roe_a * roe_a)) - kz_bar / ROE_AVER[Gid].roe_rho_jm;
	l_eigen_Q[1].N2 = ky_bar * 0.4*ROE_AVER[Gid].roe_v_jm / (roe_a * roe_a);
	l_eigen_Q[1].N3 = (ky_bar*0.4*ROE_AVER[Gid].roe_w_jm / (roe_a * roe_a)) + (kx_bar / ROE_AVER[Gid].roe_rho_jm);
	l_eigen_Q[1].N4 = -ky_bar * 0.4 / (roe_a * roe_a);

	l_eigen_Q[2].N0 = kz_bar * (1.0 - (phi_sq / (roe_a * roe_a))) - ((ky_bar*ROE_AVER[Gid].roe_u_jm - kx_bar * ROE_AVER[Gid].roe_v_jm) / ROE_AVER[Gid].roe_rho_jm);
	l_eigen_Q[2].N1 = kz_bar * 0.4*ROE_AVER[Gid].roe_u_jm / (roe_a * roe_a) + (ky_bar / ROE_AVER[Gid].roe_rho_jm);
	l_eigen_Q[2].N2 = kz_bar * 0.4*ROE_AVER[Gid].roe_v_jm / (roe_a * roe_a) - (kx_bar / ROE_AVER[Gid].roe_rho_jm);
	l_eigen_Q[2].N3 = kz_bar * 0.4*ROE_AVER[Gid].roe_w_jm / (roe_a * roe_a);
	l_eigen_Q[2].N4 = -kz_bar * 0.4 / (roe_a * roe_a);

	l_eigen_Q[3].N0 = beta * (phi_sq - theta * roe_a);
	l_eigen_Q[3].N1 = -beta * (0.4*ROE_AVER[Gid].roe_u_jm - kx_bar * roe_a);
	l_eigen_Q[3].N2 = -beta * (0.4*ROE_AVER[Gid].roe_v_jm - ky_bar * roe_a);
	l_eigen_Q[3].N3 = -beta * (0.4*ROE_AVER[Gid].roe_w_jm - kz_bar * roe_a);
	l_eigen_Q[3].N4 = beta * 0.4;

	l_eigen_Q[4].N0 = beta * (phi_sq + theta * roe_a);
	l_eigen_Q[4].N1 = -beta * (0.4*ROE_AVER[Gid].roe_u_jm + kx_bar * roe_a);
	l_eigen_Q[4].N2 = -beta * (0.4*ROE_AVER[Gid].roe_v_jm + ky_bar * roe_a);
	l_eigen_Q[4].N3 = -beta * (0.4*ROE_AVER[Gid].roe_w_jm + kz_bar * roe_a);
	l_eigen_Q[4].N4 = beta * 0.4;

	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/
    Qj_iminus_JP3.N0 = l_eigen_Q[0].N0 * Q_JP3.N0 + l_eigen_Q[0].N1 * Q_JP3.N1 + l_eigen_Q[0].N2 * Q_JP3.N2 + l_eigen_Q[0].N3 * Q_JP3.N3 + l_eigen_Q[0].N4 * Q_JP3.N4;
    Qj_iminus_JP3.N1 = l_eigen_Q[1].N0 * Q_JP3.N0 + l_eigen_Q[1].N1 * Q_JP3.N1 + l_eigen_Q[1].N2 * Q_JP3.N2 + l_eigen_Q[1].N3 * Q_JP3.N3 + l_eigen_Q[1].N4 * Q_JP3.N4;
    Qj_iminus_JP3.N2 = l_eigen_Q[2].N0 * Q_JP3.N0 + l_eigen_Q[2].N1 * Q_JP3.N1 + l_eigen_Q[2].N2 * Q_JP3.N2 + l_eigen_Q[2].N3 * Q_JP3.N3 + l_eigen_Q[2].N4 * Q_JP3.N4;
    Qj_iminus_JP3.N3 = l_eigen_Q[3].N0 * Q_JP3.N0 + l_eigen_Q[3].N1 * Q_JP3.N1 + l_eigen_Q[3].N2 * Q_JP3.N2 + l_eigen_Q[3].N3 * Q_JP3.N3 + l_eigen_Q[3].N4 * Q_JP3.N4;
    Qj_iminus_JP3.N4 = l_eigen_Q[4].N0 * Q_JP3.N0 + l_eigen_Q[4].N1 * Q_JP3.N1 + l_eigen_Q[4].N2 * Q_JP3.N2 + l_eigen_Q[4].N3 * Q_JP3.N3 + l_eigen_Q[4].N4 * Q_JP3.N4;
    
    Qj_iminus_JP2.N0 = l_eigen_Q[0].N0 * Q_JP2.N0 + l_eigen_Q[0].N1 * Q_JP2.N1 + l_eigen_Q[0].N2 * Q_JP2.N2 + l_eigen_Q[0].N3 * Q_JP2.N3 + l_eigen_Q[0].N4 * Q_JP2.N4;
    Qj_iminus_JP2.N1 = l_eigen_Q[1].N0 * Q_JP2.N0 + l_eigen_Q[1].N1 * Q_JP2.N1 + l_eigen_Q[1].N2 * Q_JP2.N2 + l_eigen_Q[1].N3 * Q_JP2.N3 + l_eigen_Q[1].N4 * Q_JP2.N4;
    Qj_iminus_JP2.N2 = l_eigen_Q[2].N0 * Q_JP2.N0 + l_eigen_Q[2].N1 * Q_JP2.N1 + l_eigen_Q[2].N2 * Q_JP2.N2 + l_eigen_Q[2].N3 * Q_JP2.N3 + l_eigen_Q[2].N4 * Q_JP2.N4;
    Qj_iminus_JP2.N3 = l_eigen_Q[3].N0 * Q_JP2.N0 + l_eigen_Q[3].N1 * Q_JP2.N1 + l_eigen_Q[3].N2 * Q_JP2.N2 + l_eigen_Q[3].N3 * Q_JP2.N3 + l_eigen_Q[3].N4 * Q_JP2.N4;
    Qj_iminus_JP2.N4 = l_eigen_Q[4].N0 * Q_JP2.N0 + l_eigen_Q[4].N1 * Q_JP2.N1 + l_eigen_Q[4].N2 * Q_JP2.N2 + l_eigen_Q[4].N3 * Q_JP2.N3 + l_eigen_Q[4].N4 * Q_JP2.N4;
    
    Qj_iminus_JP1.N0 = l_eigen_Q[0].N0 * Q_JP1.N0 + l_eigen_Q[0].N1 * Q_JP1.N1 + l_eigen_Q[0].N2 * Q_JP1.N2 + l_eigen_Q[0].N3 * Q_JP1.N3 + l_eigen_Q[0].N4 * Q_JP1.N4;
    Qj_iminus_JP1.N1 = l_eigen_Q[1].N0 * Q_JP1.N0 + l_eigen_Q[1].N1 * Q_JP1.N1 + l_eigen_Q[1].N2 * Q_JP1.N2 + l_eigen_Q[1].N3 * Q_JP1.N3 + l_eigen_Q[1].N4 * Q_JP1.N4;
    Qj_iminus_JP1.N2 = l_eigen_Q[2].N0 * Q_JP1.N0 + l_eigen_Q[2].N1 * Q_JP1.N1 + l_eigen_Q[2].N2 * Q_JP1.N2 + l_eigen_Q[2].N3 * Q_JP1.N3 + l_eigen_Q[2].N4 * Q_JP1.N4;
    Qj_iminus_JP1.N3 = l_eigen_Q[3].N0 * Q_JP1.N0 + l_eigen_Q[3].N1 * Q_JP1.N1 + l_eigen_Q[3].N2 * Q_JP1.N2 + l_eigen_Q[3].N3 * Q_JP1.N3 + l_eigen_Q[3].N4 * Q_JP1.N4;
    Qj_iminus_JP1.N4 = l_eigen_Q[4].N0 * Q_JP1.N0 + l_eigen_Q[4].N1 * Q_JP1.N1 + l_eigen_Q[4].N2 * Q_JP1.N2 + l_eigen_Q[4].N3 * Q_JP1.N3 + l_eigen_Q[4].N4 * Q_JP1.N4;
    
    Qj_iminus_I.N0 = l_eigen_Q[0].N0 * Q_I.N0 + l_eigen_Q[0].N1 * Q_I.N1 + l_eigen_Q[0].N2 * Q_I.N2 + l_eigen_Q[0].N3 * Q_I.N3 + l_eigen_Q[0].N4 * Q_I.N4;
    Qj_iminus_I.N1 = l_eigen_Q[1].N0 * Q_I.N0 + l_eigen_Q[1].N1 * Q_I.N1 + l_eigen_Q[1].N2 * Q_I.N2 + l_eigen_Q[1].N3 * Q_I.N3 + l_eigen_Q[1].N4 * Q_I.N4;
    Qj_iminus_I.N2 = l_eigen_Q[2].N0 * Q_I.N0 + l_eigen_Q[2].N1 * Q_I.N1 + l_eigen_Q[2].N2 * Q_I.N2 + l_eigen_Q[2].N3 * Q_I.N3 + l_eigen_Q[2].N4 * Q_I.N4;
    Qj_iminus_I.N3 = l_eigen_Q[3].N0 * Q_I.N0 + l_eigen_Q[3].N1 * Q_I.N1 + l_eigen_Q[3].N2 * Q_I.N2 + l_eigen_Q[3].N3 * Q_I.N3 + l_eigen_Q[3].N4 * Q_I.N4;
    Qj_iminus_I.N4 = l_eigen_Q[4].N0 * Q_I.N0 + l_eigen_Q[4].N1 * Q_I.N1 + l_eigen_Q[4].N2 * Q_I.N2 + l_eigen_Q[4].N3 * Q_I.N3 + l_eigen_Q[4].N4 * Q_I.N4;
    
    Qj_iminus_JM1.N0 = l_eigen_Q[0].N0 * Q_JM1.N0 + l_eigen_Q[0].N1 * Q_JM1.N1 + l_eigen_Q[0].N2 * Q_JM1.N2 + l_eigen_Q[0].N3 * Q_JM1.N3 + l_eigen_Q[0].N4 * Q_JM1.N4;
    Qj_iminus_JM1.N1 = l_eigen_Q[1].N0 * Q_JM1.N0 + l_eigen_Q[1].N1 * Q_JM1.N1 + l_eigen_Q[1].N2 * Q_JM1.N2 + l_eigen_Q[1].N3 * Q_JM1.N3 + l_eigen_Q[1].N4 * Q_JM1.N4;
    Qj_iminus_JM1.N2 = l_eigen_Q[2].N0 * Q_JM1.N0 + l_eigen_Q[2].N1 * Q_JM1.N1 + l_eigen_Q[2].N2 * Q_JM1.N2 + l_eigen_Q[2].N3 * Q_JM1.N3 + l_eigen_Q[2].N4 * Q_JM1.N4;
    Qj_iminus_JM1.N3 = l_eigen_Q[3].N0 * Q_JM1.N0 + l_eigen_Q[3].N1 * Q_JM1.N1 + l_eigen_Q[3].N2 * Q_JM1.N2 + l_eigen_Q[3].N3 * Q_JM1.N3 + l_eigen_Q[3].N4 * Q_JM1.N4;
    Qj_iminus_JM1.N4 = l_eigen_Q[4].N0 * Q_JM1.N0 + l_eigen_Q[4].N1 * Q_JM1.N1 + l_eigen_Q[4].N2 * Q_JM1.N2 + l_eigen_Q[4].N3 * Q_JM1.N3 + l_eigen_Q[4].N4 * Q_JM1.N4;
    
    Qj_iminus_JM2.N0 = l_eigen_Q[0].N0 * Q_JM2.N0 + l_eigen_Q[0].N1 * Q_JM2.N1 + l_eigen_Q[0].N2 * Q_JM2.N2 + l_eigen_Q[0].N3 * Q_JM2.N3 + l_eigen_Q[0].N4 * Q_JM2.N4;
    Qj_iminus_JM2.N1 = l_eigen_Q[1].N0 * Q_JM2.N0 + l_eigen_Q[1].N1 * Q_JM2.N1 + l_eigen_Q[1].N2 * Q_JM2.N2 + l_eigen_Q[1].N3 * Q_JM2.N3 + l_eigen_Q[1].N4 * Q_JM2.N4;
    Qj_iminus_JM2.N2 = l_eigen_Q[2].N0 * Q_JM2.N0 + l_eigen_Q[2].N1 * Q_JM2.N1 + l_eigen_Q[2].N2 * Q_JM2.N2 + l_eigen_Q[2].N3 * Q_JM2.N3 + l_eigen_Q[2].N4 * Q_JM2.N4;
    Qj_iminus_JM2.N3 = l_eigen_Q[3].N0 * Q_JM2.N0 + l_eigen_Q[3].N1 * Q_JM2.N1 + l_eigen_Q[3].N2 * Q_JM2.N2 + l_eigen_Q[3].N3 * Q_JM2.N3 + l_eigen_Q[3].N4 * Q_JM2.N4;
    Qj_iminus_JM2.N4 = l_eigen_Q[4].N0 * Q_JM2.N0 + l_eigen_Q[4].N1 * Q_JM2.N1 + l_eigen_Q[4].N2 * Q_JM2.N2 + l_eigen_Q[4].N3 * Q_JM2.N3 + l_eigen_Q[4].N4 * Q_JM2.N4;
    
    Qj_iminus_JM3.N0 = l_eigen_Q[0].N0 * Q_JM3.N0 + l_eigen_Q[0].N1 * Q_JM3.N1 + l_eigen_Q[0].N2 * Q_JM3.N2 + l_eigen_Q[0].N3 * Q_JM3.N3 + l_eigen_Q[0].N4 * Q_JM3.N4;
    Qj_iminus_JM3.N1 = l_eigen_Q[1].N0 * Q_JM3.N0 + l_eigen_Q[1].N1 * Q_JM3.N1 + l_eigen_Q[1].N2 * Q_JM3.N2 + l_eigen_Q[1].N3 * Q_JM3.N3 + l_eigen_Q[1].N4 * Q_JM3.N4;
    Qj_iminus_JM3.N2 = l_eigen_Q[2].N0 * Q_JM3.N0 + l_eigen_Q[2].N1 * Q_JM3.N1 + l_eigen_Q[2].N2 * Q_JM3.N2 + l_eigen_Q[2].N3 * Q_JM3.N3 + l_eigen_Q[2].N4 * Q_JM3.N4;
    Qj_iminus_JM3.N3 = l_eigen_Q[3].N0 * Q_JM3.N0 + l_eigen_Q[3].N1 * Q_JM3.N1 + l_eigen_Q[3].N2 * Q_JM3.N2 + l_eigen_Q[3].N3 * Q_JM3.N3 + l_eigen_Q[3].N4 * Q_JM3.N4;
    Qj_iminus_JM3.N4 = l_eigen_Q[4].N0 * Q_JM3.N0 + l_eigen_Q[4].N1 * Q_JM3.N1 + l_eigen_Q[4].N2 * Q_JM3.N2 + l_eigen_Q[4].N3 * Q_JM3.N3 + l_eigen_Q[4].N4 * Q_JM3.N4;
    
    
    /***********************************************************************************************/
	/****************************************** Qk_iplus********************************************/
    kk = sqrt(metric[Gid].xi_xkp*metric[Gid].xi_xkp + metric[Gid].xi_ykp * metric[Gid].xi_ykp + metric[Gid].xi_zkp * metric[Gid].xi_zkp);
    roe_a = sqrt(0.4*(ROE_AVER[Gid].roe_h_kp - 0.5*(ROE_AVER[Gid].roe_u_kp * ROE_AVER[Gid].roe_u_kp + ROE_AVER[Gid].roe_v_kp * ROE_AVER[Gid].roe_v_kp + ROE_AVER[Gid].roe_w_kp * ROE_AVER[Gid].roe_w_kp)));

	kx_bar = metric[Gid].xi_xkp / kk;
	ky_bar = metric[Gid].xi_ykp / kk;
	kz_bar = metric[Gid].xi_zkp / kk;

	theta = kx_bar * ROE_AVER[Gid].roe_u_kp + ky_bar * ROE_AVER[Gid].roe_v_kp + kz_bar * ROE_AVER[Gid].roe_w_kp;
	phi_sq = 0.5*0.4*(ROE_AVER[Gid].roe_u_kp * ROE_AVER[Gid].roe_u_kp + ROE_AVER[Gid].roe_v_kp * ROE_AVER[Gid].roe_v_kp + ROE_AVER[Gid].roe_w_kp * ROE_AVER[Gid].roe_w_kp);
	alpha = ROE_AVER[Gid].roe_rho_kp / (sqrt(2.0)*roe_a);
	beta = 1.0 / (sqrt(2.0)*ROE_AVER[Gid].roe_rho_kp * roe_a);

	r_eigen_Qkp[0].N0 = kx_bar;
	r_eigen_Qkp[0].N1 = ky_bar;
	r_eigen_Qkp[0].N2 = kz_bar;
	r_eigen_Qkp[0].N3 = alpha;
	r_eigen_Qkp[0].N4 = alpha;

	r_eigen_Qkp[1].N0 = kx_bar * ROE_AVER[Gid].roe_u_kp;
	r_eigen_Qkp[1].N1 = ky_bar * ROE_AVER[Gid].roe_u_kp - kz_bar * ROE_AVER[Gid].roe_rho_kp;
	r_eigen_Qkp[1].N2 = kz_bar * ROE_AVER[Gid].roe_u_kp + ky_bar * ROE_AVER[Gid].roe_rho_kp;
	r_eigen_Qkp[1].N3 = alpha * (ROE_AVER[Gid].roe_u_kp + kx_bar * roe_a);
	r_eigen_Qkp[1].N4 = alpha * (ROE_AVER[Gid].roe_u_kp - kx_bar * roe_a);

	r_eigen_Qkp[2].N0 = kx_bar * ROE_AVER[Gid].roe_v_kp + kz_bar * ROE_AVER[Gid].roe_rho_kp;
	r_eigen_Qkp[2].N1 = ky_bar * ROE_AVER[Gid].roe_v_kp;
	r_eigen_Qkp[2].N2 = kz_bar * ROE_AVER[Gid].roe_v_kp - kx_bar * ROE_AVER[Gid].roe_rho_kp;
	r_eigen_Qkp[2].N3 = alpha * (ROE_AVER[Gid].roe_v_kp + ky_bar * roe_a);
	r_eigen_Qkp[2].N4 = alpha * (ROE_AVER[Gid].roe_v_kp - ky_bar * roe_a);

	r_eigen_Qkp[3].N0 = kx_bar * ROE_AVER[Gid].roe_w_kp - ky_bar * ROE_AVER[Gid].roe_rho_kp;
	r_eigen_Qkp[3].N1 = ky_bar * ROE_AVER[Gid].roe_w_kp + kx_bar * ROE_AVER[Gid].roe_rho_kp;
	r_eigen_Qkp[3].N2 = kz_bar * ROE_AVER[Gid].roe_w_kp;
	r_eigen_Qkp[3].N3 = alpha * (ROE_AVER[Gid].roe_w_kp + kz_bar * roe_a);
	r_eigen_Qkp[3].N4 = alpha * (ROE_AVER[Gid].roe_w_kp - kz_bar * roe_a);

	r_eigen_Qkp[4].N0 = ((kx_bar*phi_sq) / 0.4) + ROE_AVER[Gid].roe_rho_kp * (kz_bar*ROE_AVER[Gid].roe_v_kp - ky_bar * ROE_AVER[Gid].roe_w_kp);
	r_eigen_Qkp[4].N1 = ((ky_bar*phi_sq) / 0.4) + ROE_AVER[Gid].roe_rho_kp * (kx_bar*ROE_AVER[Gid].roe_w_kp - kz_bar * ROE_AVER[Gid].roe_u_kp);
	r_eigen_Qkp[4].N2 = ((kz_bar*phi_sq) / 0.4) + ROE_AVER[Gid].roe_rho_kp * (ky_bar*ROE_AVER[Gid].roe_u_kp - kx_bar * ROE_AVER[Gid].roe_v_kp);
	r_eigen_Qkp[4].N3 = alpha * (((phi_sq + roe_a * roe_a) / (0.4)) + theta * roe_a);
	r_eigen_Qkp[4].N4 = alpha * (((phi_sq + roe_a * roe_a) / (0.4)) - theta * roe_a);
	/*************************************************************************************************/

	l_eigen_Q[0].N0 = kx_bar * (1.0 - (phi_sq / (roe_a * roe_a))) - ((kz_bar*ROE_AVER[Gid].roe_v_kp - ky_bar * ROE_AVER[Gid].roe_w_kp) / ROE_AVER[Gid].roe_rho_kp);
	l_eigen_Q[0].N1 = kx_bar * 0.4*ROE_AVER[Gid].roe_u_kp / (roe_a * roe_a);
	l_eigen_Q[0].N2 = (kx_bar*0.4*ROE_AVER[Gid].roe_v_kp / (roe_a * roe_a)) + (kz_bar / ROE_AVER[Gid].roe_rho_kp);
	l_eigen_Q[0].N3 = (kx_bar*0.4*ROE_AVER[Gid].roe_w_kp / (roe_a * roe_a)) - (ky_bar / ROE_AVER[Gid].roe_rho_kp);
	l_eigen_Q[0].N4 = -kx_bar * 0.4 / (roe_a * roe_a);

	l_eigen_Q[1].N0 = ky_bar * (1.0 - (phi_sq / (roe_a * roe_a))) - ((kx_bar*ROE_AVER[Gid].roe_w_kp - kz_bar * ROE_AVER[Gid].roe_u_kp) / ROE_AVER[Gid].roe_rho_kp);
	l_eigen_Q[1].N1 = (ky_bar*0.4*ROE_AVER[Gid].roe_u_kp / (roe_a * roe_a)) - kz_bar / ROE_AVER[Gid].roe_rho_kp;
	l_eigen_Q[1].N2 = ky_bar * 0.4*ROE_AVER[Gid].roe_v_kp / (roe_a * roe_a);
	l_eigen_Q[1].N3 = (ky_bar*0.4*ROE_AVER[Gid].roe_w_kp / (roe_a * roe_a)) + (kx_bar / ROE_AVER[Gid].roe_rho_kp);
	l_eigen_Q[1].N4 = -ky_bar * 0.4 / (roe_a * roe_a);

	l_eigen_Q[2].N0 = kz_bar * (1.0 - (phi_sq / (roe_a * roe_a))) - ((ky_bar*ROE_AVER[Gid].roe_u_kp - kx_bar * ROE_AVER[Gid].roe_v_kp) / ROE_AVER[Gid].roe_rho_kp);
	l_eigen_Q[2].N1 = kz_bar * 0.4*ROE_AVER[Gid].roe_u_kp / (roe_a * roe_a) + (ky_bar / ROE_AVER[Gid].roe_rho_kp);
	l_eigen_Q[2].N2 = kz_bar * 0.4*ROE_AVER[Gid].roe_v_kp / (roe_a * roe_a) - (kx_bar / ROE_AVER[Gid].roe_rho_kp);
	l_eigen_Q[2].N3 = kz_bar * 0.4*ROE_AVER[Gid].roe_w_kp / (roe_a * roe_a);
	l_eigen_Q[2].N4 = -kz_bar * 0.4 / (roe_a * roe_a);

	l_eigen_Q[3].N0 = beta * (phi_sq - theta * roe_a);
	l_eigen_Q[3].N1 = -beta * (0.4*ROE_AVER[Gid].roe_u_kp - kx_bar * roe_a);
	l_eigen_Q[3].N2 = -beta * (0.4*ROE_AVER[Gid].roe_v_kp - ky_bar * roe_a);
	l_eigen_Q[3].N3 = -beta * (0.4*ROE_AVER[Gid].roe_w_kp - kz_bar * roe_a);
	l_eigen_Q[3].N4 = beta * 0.4;

	l_eigen_Q[4].N0 = beta * (phi_sq + theta * roe_a);
	l_eigen_Q[4].N1 = -beta * (0.4*ROE_AVER[Gid].roe_u_kp + kx_bar * roe_a);
	l_eigen_Q[4].N2 = -beta * (0.4*ROE_AVER[Gid].roe_v_kp + ky_bar * roe_a);
	l_eigen_Q[4].N3 = -beta * (0.4*ROE_AVER[Gid].roe_w_kp + kz_bar * roe_a);
	l_eigen_Q[4].N4 = beta * 0.4;

	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/
    Qk_iplus_KP3.N0 = l_eigen_Q[0].N0 * Q_KP3.N0 + l_eigen_Q[0].N1 * Q_KP3.N1 + l_eigen_Q[0].N2 * Q_KP3.N2 + l_eigen_Q[0].N3 * Q_KP3.N3 + l_eigen_Q[0].N4 * Q_KP3.N4;
    Qk_iplus_KP3.N1 = l_eigen_Q[1].N0 * Q_KP3.N0 + l_eigen_Q[1].N1 * Q_KP3.N1 + l_eigen_Q[1].N2 * Q_KP3.N2 + l_eigen_Q[1].N3 * Q_KP3.N3 + l_eigen_Q[1].N4 * Q_KP3.N4;
    Qk_iplus_KP3.N2 = l_eigen_Q[2].N0 * Q_KP3.N0 + l_eigen_Q[2].N1 * Q_KP3.N1 + l_eigen_Q[2].N2 * Q_KP3.N2 + l_eigen_Q[2].N3 * Q_KP3.N3 + l_eigen_Q[2].N4 * Q_KP3.N4;
    Qk_iplus_KP3.N3 = l_eigen_Q[3].N0 * Q_KP3.N0 + l_eigen_Q[3].N1 * Q_KP3.N1 + l_eigen_Q[3].N2 * Q_KP3.N2 + l_eigen_Q[3].N3 * Q_KP3.N3 + l_eigen_Q[3].N4 * Q_KP3.N4;
    Qk_iplus_KP3.N4 = l_eigen_Q[4].N0 * Q_KP3.N0 + l_eigen_Q[4].N1 * Q_KP3.N1 + l_eigen_Q[4].N2 * Q_KP3.N2 + l_eigen_Q[4].N3 * Q_KP3.N3 + l_eigen_Q[4].N4 * Q_KP3.N4;
    
    Qk_iplus_KP2.N0 = l_eigen_Q[0].N0 * Q_KP2.N0 + l_eigen_Q[0].N1 * Q_KP2.N1 + l_eigen_Q[0].N2 * Q_KP2.N2 + l_eigen_Q[0].N3 * Q_KP2.N3 + l_eigen_Q[0].N4 * Q_KP2.N4;
    Qk_iplus_KP2.N1 = l_eigen_Q[1].N0 * Q_KP2.N0 + l_eigen_Q[1].N1 * Q_KP2.N1 + l_eigen_Q[1].N2 * Q_KP2.N2 + l_eigen_Q[1].N3 * Q_KP2.N3 + l_eigen_Q[1].N4 * Q_KP2.N4;
    Qk_iplus_KP2.N2 = l_eigen_Q[2].N0 * Q_KP2.N0 + l_eigen_Q[2].N1 * Q_KP2.N1 + l_eigen_Q[2].N2 * Q_KP2.N2 + l_eigen_Q[2].N3 * Q_KP2.N3 + l_eigen_Q[2].N4 * Q_KP2.N4;
    Qk_iplus_KP2.N3 = l_eigen_Q[3].N0 * Q_KP2.N0 + l_eigen_Q[3].N1 * Q_KP2.N1 + l_eigen_Q[3].N2 * Q_KP2.N2 + l_eigen_Q[3].N3 * Q_KP2.N3 + l_eigen_Q[3].N4 * Q_KP2.N4;
    Qk_iplus_KP2.N4 = l_eigen_Q[4].N0 * Q_KP2.N0 + l_eigen_Q[4].N1 * Q_KP2.N1 + l_eigen_Q[4].N2 * Q_KP2.N2 + l_eigen_Q[4].N3 * Q_KP2.N3 + l_eigen_Q[4].N4 * Q_KP2.N4;
    
    Qk_iplus_KP1.N0 = l_eigen_Q[0].N0 * Q_KP1.N0 + l_eigen_Q[0].N1 * Q_KP1.N1 + l_eigen_Q[0].N2 * Q_KP1.N2 + l_eigen_Q[0].N3 * Q_KP1.N3 + l_eigen_Q[0].N4 * Q_KP1.N4;
    Qk_iplus_KP1.N1 = l_eigen_Q[1].N0 * Q_KP1.N0 + l_eigen_Q[1].N1 * Q_KP1.N1 + l_eigen_Q[1].N2 * Q_KP1.N2 + l_eigen_Q[1].N3 * Q_KP1.N3 + l_eigen_Q[1].N4 * Q_KP1.N4;
    Qk_iplus_KP1.N2 = l_eigen_Q[2].N0 * Q_KP1.N0 + l_eigen_Q[2].N1 * Q_KP1.N1 + l_eigen_Q[2].N2 * Q_KP1.N2 + l_eigen_Q[2].N3 * Q_KP1.N3 + l_eigen_Q[2].N4 * Q_KP1.N4;
    Qk_iplus_KP1.N3 = l_eigen_Q[3].N0 * Q_KP1.N0 + l_eigen_Q[3].N1 * Q_KP1.N1 + l_eigen_Q[3].N2 * Q_KP1.N2 + l_eigen_Q[3].N3 * Q_KP1.N3 + l_eigen_Q[3].N4 * Q_KP1.N4;
    Qk_iplus_KP1.N4 = l_eigen_Q[4].N0 * Q_KP1.N0 + l_eigen_Q[4].N1 * Q_KP1.N1 + l_eigen_Q[4].N2 * Q_KP1.N2 + l_eigen_Q[4].N3 * Q_KP1.N3 + l_eigen_Q[4].N4 * Q_KP1.N4;
    
    Qk_iplus_I.N0 = l_eigen_Q[0].N0 * Q_I.N0 + l_eigen_Q[0].N1 * Q_I.N1 + l_eigen_Q[0].N2 * Q_I.N2 + l_eigen_Q[0].N3 * Q_I.N3 + l_eigen_Q[0].N4 * Q_I.N4;
    Qk_iplus_I.N1 = l_eigen_Q[1].N0 * Q_I.N0 + l_eigen_Q[1].N1 * Q_I.N1 + l_eigen_Q[1].N2 * Q_I.N2 + l_eigen_Q[1].N3 * Q_I.N3 + l_eigen_Q[1].N4 * Q_I.N4;
    Qk_iplus_I.N2 = l_eigen_Q[2].N0 * Q_I.N0 + l_eigen_Q[2].N1 * Q_I.N1 + l_eigen_Q[2].N2 * Q_I.N2 + l_eigen_Q[2].N3 * Q_I.N3 + l_eigen_Q[2].N4 * Q_I.N4;
    Qk_iplus_I.N3 = l_eigen_Q[3].N0 * Q_I.N0 + l_eigen_Q[3].N1 * Q_I.N1 + l_eigen_Q[3].N2 * Q_I.N2 + l_eigen_Q[3].N3 * Q_I.N3 + l_eigen_Q[3].N4 * Q_I.N4;
    Qk_iplus_I.N4 = l_eigen_Q[4].N0 * Q_I.N0 + l_eigen_Q[4].N1 * Q_I.N1 + l_eigen_Q[4].N2 * Q_I.N2 + l_eigen_Q[4].N3 * Q_I.N3 + l_eigen_Q[4].N4 * Q_I.N4;
    
    Qk_iplus_KM1.N0 = l_eigen_Q[0].N0 * Q_KM1.N0 + l_eigen_Q[0].N1 * Q_KM1.N1 + l_eigen_Q[0].N2 * Q_KM1.N2 + l_eigen_Q[0].N3 * Q_KM1.N3 + l_eigen_Q[0].N4 * Q_KM1.N4;
    Qk_iplus_KM1.N1 = l_eigen_Q[1].N0 * Q_KM1.N0 + l_eigen_Q[1].N1 * Q_KM1.N1 + l_eigen_Q[1].N2 * Q_KM1.N2 + l_eigen_Q[1].N3 * Q_KM1.N3 + l_eigen_Q[1].N4 * Q_KM1.N4;
    Qk_iplus_KM1.N2 = l_eigen_Q[2].N0 * Q_KM1.N0 + l_eigen_Q[2].N1 * Q_KM1.N1 + l_eigen_Q[2].N2 * Q_KM1.N2 + l_eigen_Q[2].N3 * Q_KM1.N3 + l_eigen_Q[2].N4 * Q_KM1.N4;
    Qk_iplus_KM1.N3 = l_eigen_Q[3].N0 * Q_KM1.N0 + l_eigen_Q[3].N1 * Q_KM1.N1 + l_eigen_Q[3].N2 * Q_KM1.N2 + l_eigen_Q[3].N3 * Q_KM1.N3 + l_eigen_Q[3].N4 * Q_KM1.N4;
    Qk_iplus_KM1.N4 = l_eigen_Q[4].N0 * Q_KM1.N0 + l_eigen_Q[4].N1 * Q_KM1.N1 + l_eigen_Q[4].N2 * Q_KM1.N2 + l_eigen_Q[4].N3 * Q_KM1.N3 + l_eigen_Q[4].N4 * Q_KM1.N4;
    
    Qk_iplus_KM2.N0 = l_eigen_Q[0].N0 * Q_KM2.N0 + l_eigen_Q[0].N1 * Q_KM2.N1 + l_eigen_Q[0].N2 * Q_KM2.N2 + l_eigen_Q[0].N3 * Q_KM2.N3 + l_eigen_Q[0].N4 * Q_KM2.N4;
    Qk_iplus_KM2.N1 = l_eigen_Q[1].N0 * Q_KM2.N0 + l_eigen_Q[1].N1 * Q_KM2.N1 + l_eigen_Q[1].N2 * Q_KM2.N2 + l_eigen_Q[1].N3 * Q_KM2.N3 + l_eigen_Q[1].N4 * Q_KM2.N4;
    Qk_iplus_KM2.N2 = l_eigen_Q[2].N0 * Q_KM2.N0 + l_eigen_Q[2].N1 * Q_KM2.N1 + l_eigen_Q[2].N2 * Q_KM2.N2 + l_eigen_Q[2].N3 * Q_KM2.N3 + l_eigen_Q[2].N4 * Q_KM2.N4;
    Qk_iplus_KM2.N3 = l_eigen_Q[3].N0 * Q_KM2.N0 + l_eigen_Q[3].N1 * Q_KM2.N1 + l_eigen_Q[3].N2 * Q_KM2.N2 + l_eigen_Q[3].N3 * Q_KM2.N3 + l_eigen_Q[3].N4 * Q_KM2.N4;
    Qk_iplus_KM2.N4 = l_eigen_Q[4].N0 * Q_KM2.N0 + l_eigen_Q[4].N1 * Q_KM2.N1 + l_eigen_Q[4].N2 * Q_KM2.N2 + l_eigen_Q[4].N3 * Q_KM2.N3 + l_eigen_Q[4].N4 * Q_KM2.N4;
    
    Qk_iplus_KM3.N0 = l_eigen_Q[0].N0 * Q_KM3.N0 + l_eigen_Q[0].N1 * Q_KM3.N1 + l_eigen_Q[0].N2 * Q_KM3.N2 + l_eigen_Q[0].N3 * Q_KM3.N3 + l_eigen_Q[0].N4 * Q_KM3.N4;
    Qk_iplus_KM3.N1 = l_eigen_Q[1].N0 * Q_KM3.N0 + l_eigen_Q[1].N1 * Q_KM3.N1 + l_eigen_Q[1].N2 * Q_KM3.N2 + l_eigen_Q[1].N3 * Q_KM3.N3 + l_eigen_Q[1].N4 * Q_KM3.N4;
    Qk_iplus_KM3.N2 = l_eigen_Q[2].N0 * Q_KM3.N0 + l_eigen_Q[2].N1 * Q_KM3.N1 + l_eigen_Q[2].N2 * Q_KM3.N2 + l_eigen_Q[2].N3 * Q_KM3.N3 + l_eigen_Q[2].N4 * Q_KM3.N4;
    Qk_iplus_KM3.N3 = l_eigen_Q[3].N0 * Q_KM3.N0 + l_eigen_Q[3].N1 * Q_KM3.N1 + l_eigen_Q[3].N2 * Q_KM3.N2 + l_eigen_Q[3].N3 * Q_KM3.N3 + l_eigen_Q[3].N4 * Q_KM3.N4;
    Qk_iplus_KM3.N4 = l_eigen_Q[4].N0 * Q_KM3.N0 + l_eigen_Q[4].N1 * Q_KM3.N1 + l_eigen_Q[4].N2 * Q_KM3.N2 + l_eigen_Q[4].N3 * Q_KM3.N3 + l_eigen_Q[4].N4 * Q_KM3.N4;
    

    /***********************************************************************************************/
	/****************************************** Qk_iminus********************************************/
    kk = sqrt(metric[Gid].xi_xkm*metric[Gid].xi_xkm + metric[Gid].xi_ykm * metric[Gid].xi_ykm + metric[Gid].xi_zkm * metric[Gid].xi_zkm);
    roe_a = sqrt(0.4*(ROE_AVER[Gid].roe_h_km - 0.5*(ROE_AVER[Gid].roe_u_km * ROE_AVER[Gid].roe_u_km + ROE_AVER[Gid].roe_v_km * ROE_AVER[Gid].roe_v_km + ROE_AVER[Gid].roe_w_km * ROE_AVER[Gid].roe_w_km)));

	kx_bar = metric[Gid].xi_xkm / kk;
	ky_bar = metric[Gid].xi_ykm / kk;
	kz_bar = metric[Gid].xi_zkm / kk;

	theta = kx_bar * ROE_AVER[Gid].roe_u_km + ky_bar * ROE_AVER[Gid].roe_v_km + kz_bar * ROE_AVER[Gid].roe_w_km;
	phi_sq = 0.5*0.4*(ROE_AVER[Gid].roe_u_km * ROE_AVER[Gid].roe_u_km + ROE_AVER[Gid].roe_v_km * ROE_AVER[Gid].roe_v_km + ROE_AVER[Gid].roe_w_km * ROE_AVER[Gid].roe_w_km);
	alpha = ROE_AVER[Gid].roe_rho_km / (sqrt(2.0)*roe_a);
	beta = 1.0 / (sqrt(2.0)*ROE_AVER[Gid].roe_rho_km * roe_a);

	r_eigen_Qkm[0].N0 = kx_bar;
	r_eigen_Qkm[0].N1 = ky_bar;
	r_eigen_Qkm[0].N2 = kz_bar;
	r_eigen_Qkm[0].N3 = alpha;
	r_eigen_Qkm[0].N4 = alpha;

	r_eigen_Qkm[1].N0 = kx_bar * ROE_AVER[Gid].roe_u_km;
	r_eigen_Qkm[1].N1 = ky_bar * ROE_AVER[Gid].roe_u_km - kz_bar * ROE_AVER[Gid].roe_rho_km;
	r_eigen_Qkm[1].N2 = kz_bar * ROE_AVER[Gid].roe_u_km + ky_bar * ROE_AVER[Gid].roe_rho_km;
	r_eigen_Qkm[1].N3 = alpha * (ROE_AVER[Gid].roe_u_km + kx_bar * roe_a);
	r_eigen_Qkm[1].N4 = alpha * (ROE_AVER[Gid].roe_u_km - kx_bar * roe_a);

	r_eigen_Qkm[2].N0 = kx_bar * ROE_AVER[Gid].roe_v_km + kz_bar * ROE_AVER[Gid].roe_rho_km;
	r_eigen_Qkm[2].N1 = ky_bar * ROE_AVER[Gid].roe_v_km;
	r_eigen_Qkm[2].N2 = kz_bar * ROE_AVER[Gid].roe_v_km - kx_bar * ROE_AVER[Gid].roe_rho_km;
	r_eigen_Qkm[2].N3 = alpha * (ROE_AVER[Gid].roe_v_km + ky_bar * roe_a);
	r_eigen_Qkm[2].N4 = alpha * (ROE_AVER[Gid].roe_v_km - ky_bar * roe_a);

	r_eigen_Qkm[3].N0 = kx_bar * ROE_AVER[Gid].roe_w_km - ky_bar * ROE_AVER[Gid].roe_rho_km;
	r_eigen_Qkm[3].N1 = ky_bar * ROE_AVER[Gid].roe_w_km + kx_bar * ROE_AVER[Gid].roe_rho_km;
	r_eigen_Qkm[3].N2 = kz_bar * ROE_AVER[Gid].roe_w_km;
	r_eigen_Qkm[3].N3 = alpha * (ROE_AVER[Gid].roe_w_km + kz_bar * roe_a);
	r_eigen_Qkm[3].N4 = alpha * (ROE_AVER[Gid].roe_w_km - kz_bar * roe_a);

	r_eigen_Qkm[4].N0 = ((kx_bar*phi_sq) / 0.4) + ROE_AVER[Gid].roe_rho_km * (kz_bar*ROE_AVER[Gid].roe_v_km - ky_bar * ROE_AVER[Gid].roe_w_km);
	r_eigen_Qkm[4].N1 = ((ky_bar*phi_sq) / 0.4) + ROE_AVER[Gid].roe_rho_km * (kx_bar*ROE_AVER[Gid].roe_w_km - kz_bar * ROE_AVER[Gid].roe_u_km);
	r_eigen_Qkm[4].N2 = ((kz_bar*phi_sq) / 0.4) + ROE_AVER[Gid].roe_rho_km * (ky_bar*ROE_AVER[Gid].roe_u_km - kx_bar * ROE_AVER[Gid].roe_v_km);
	r_eigen_Qkm[4].N3 = alpha * (((phi_sq + roe_a * roe_a) / (0.4)) + theta * roe_a);
	r_eigen_Qkm[4].N4 = alpha * (((phi_sq + roe_a * roe_a) / (0.4)) - theta * roe_a);
	/*************************************************************************************************/

	l_eigen_Q[0].N0 = kx_bar * (1.0 - (phi_sq / (roe_a * roe_a))) - ((kz_bar*ROE_AVER[Gid].roe_v_km - ky_bar * ROE_AVER[Gid].roe_w_km) / ROE_AVER[Gid].roe_rho_km);
	l_eigen_Q[0].N1 = kx_bar * 0.4*ROE_AVER[Gid].roe_u_km / (roe_a * roe_a);
	l_eigen_Q[0].N2 = (kx_bar*0.4*ROE_AVER[Gid].roe_v_km / (roe_a * roe_a)) + (kz_bar / ROE_AVER[Gid].roe_rho_km);
	l_eigen_Q[0].N3 = (kx_bar*0.4*ROE_AVER[Gid].roe_w_km / (roe_a * roe_a)) - (ky_bar / ROE_AVER[Gid].roe_rho_km);
	l_eigen_Q[0].N4 = -kx_bar * 0.4 / (roe_a * roe_a);

	l_eigen_Q[1].N0 = ky_bar * (1.0 - (phi_sq / (roe_a * roe_a))) - ((kx_bar*ROE_AVER[Gid].roe_w_km - kz_bar * ROE_AVER[Gid].roe_u_km) / ROE_AVER[Gid].roe_rho_km);
	l_eigen_Q[1].N1 = (ky_bar*0.4*ROE_AVER[Gid].roe_u_km / (roe_a * roe_a)) - kz_bar / ROE_AVER[Gid].roe_rho_km;
	l_eigen_Q[1].N2 = ky_bar * 0.4*ROE_AVER[Gid].roe_v_km / (roe_a * roe_a);
	l_eigen_Q[1].N3 = (ky_bar*0.4*ROE_AVER[Gid].roe_w_km / (roe_a * roe_a)) + (kx_bar / ROE_AVER[Gid].roe_rho_km);
	l_eigen_Q[1].N4 = -ky_bar * 0.4 / (roe_a * roe_a);

	l_eigen_Q[2].N0 = kz_bar * (1.0 - (phi_sq / (roe_a * roe_a))) - ((ky_bar*ROE_AVER[Gid].roe_u_km - kx_bar * ROE_AVER[Gid].roe_v_km) / ROE_AVER[Gid].roe_rho_km);
	l_eigen_Q[2].N1 = kz_bar * 0.4*ROE_AVER[Gid].roe_u_km / (roe_a * roe_a) + (ky_bar / ROE_AVER[Gid].roe_rho_km);
	l_eigen_Q[2].N2 = kz_bar * 0.4*ROE_AVER[Gid].roe_v_km / (roe_a * roe_a) - (kx_bar / ROE_AVER[Gid].roe_rho_km);
	l_eigen_Q[2].N3 = kz_bar * 0.4*ROE_AVER[Gid].roe_w_km / (roe_a * roe_a);
	l_eigen_Q[2].N4 = -kz_bar * 0.4 / (roe_a * roe_a);

	l_eigen_Q[3].N0 = beta * (phi_sq - theta * roe_a);
	l_eigen_Q[3].N1 = -beta * (0.4*ROE_AVER[Gid].roe_u_km - kx_bar * roe_a);
	l_eigen_Q[3].N2 = -beta * (0.4*ROE_AVER[Gid].roe_v_km - ky_bar * roe_a);
	l_eigen_Q[3].N3 = -beta * (0.4*ROE_AVER[Gid].roe_w_km - kz_bar * roe_a);
	l_eigen_Q[3].N4 = beta * 0.4;

	l_eigen_Q[4].N0 = beta * (phi_sq + theta * roe_a);
	l_eigen_Q[4].N1 = -beta * (0.4*ROE_AVER[Gid].roe_u_km + kx_bar * roe_a);
	l_eigen_Q[4].N2 = -beta * (0.4*ROE_AVER[Gid].roe_v_km + ky_bar * roe_a);
	l_eigen_Q[4].N3 = -beta * (0.4*ROE_AVER[Gid].roe_w_km + kz_bar * roe_a);
	l_eigen_Q[4].N4 = beta * 0.4;

	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/
    Qk_iminus_KP3.N0 = l_eigen_Q[0].N0 * Q_KP3.N0 + l_eigen_Q[0].N1 * Q_KP3.N1 + l_eigen_Q[0].N2 * Q_KP3.N2 + l_eigen_Q[0].N3 * Q_KP3.N3 + l_eigen_Q[0].N4 * Q_KP3.N4;
    Qk_iminus_KP3.N1 = l_eigen_Q[1].N0 * Q_KP3.N0 + l_eigen_Q[1].N1 * Q_KP3.N1 + l_eigen_Q[1].N2 * Q_KP3.N2 + l_eigen_Q[1].N3 * Q_KP3.N3 + l_eigen_Q[1].N4 * Q_KP3.N4;
    Qk_iminus_KP3.N2 = l_eigen_Q[2].N0 * Q_KP3.N0 + l_eigen_Q[2].N1 * Q_KP3.N1 + l_eigen_Q[2].N2 * Q_KP3.N2 + l_eigen_Q[2].N3 * Q_KP3.N3 + l_eigen_Q[2].N4 * Q_KP3.N4;
    Qk_iminus_KP3.N3 = l_eigen_Q[3].N0 * Q_KP3.N0 + l_eigen_Q[3].N1 * Q_KP3.N1 + l_eigen_Q[3].N2 * Q_KP3.N2 + l_eigen_Q[3].N3 * Q_KP3.N3 + l_eigen_Q[3].N4 * Q_KP3.N4;
    Qk_iminus_KP3.N4 = l_eigen_Q[4].N0 * Q_KP3.N0 + l_eigen_Q[4].N1 * Q_KP3.N1 + l_eigen_Q[4].N2 * Q_KP3.N2 + l_eigen_Q[4].N3 * Q_KP3.N3 + l_eigen_Q[4].N4 * Q_KP3.N4;
    
    Qk_iminus_KP2.N0 = l_eigen_Q[0].N0 * Q_KP2.N0 + l_eigen_Q[0].N1 * Q_KP2.N1 + l_eigen_Q[0].N2 * Q_KP2.N2 + l_eigen_Q[0].N3 * Q_KP2.N3 + l_eigen_Q[0].N4 * Q_KP2.N4;
    Qk_iminus_KP2.N1 = l_eigen_Q[1].N0 * Q_KP2.N0 + l_eigen_Q[1].N1 * Q_KP2.N1 + l_eigen_Q[1].N2 * Q_KP2.N2 + l_eigen_Q[1].N3 * Q_KP2.N3 + l_eigen_Q[1].N4 * Q_KP2.N4;
    Qk_iminus_KP2.N2 = l_eigen_Q[2].N0 * Q_KP2.N0 + l_eigen_Q[2].N1 * Q_KP2.N1 + l_eigen_Q[2].N2 * Q_KP2.N2 + l_eigen_Q[2].N3 * Q_KP2.N3 + l_eigen_Q[2].N4 * Q_KP2.N4;
    Qk_iminus_KP2.N3 = l_eigen_Q[3].N0 * Q_KP2.N0 + l_eigen_Q[3].N1 * Q_KP2.N1 + l_eigen_Q[3].N2 * Q_KP2.N2 + l_eigen_Q[3].N3 * Q_KP2.N3 + l_eigen_Q[3].N4 * Q_KP2.N4;
    Qk_iminus_KP2.N4 = l_eigen_Q[4].N0 * Q_KP2.N0 + l_eigen_Q[4].N1 * Q_KP2.N1 + l_eigen_Q[4].N2 * Q_KP2.N2 + l_eigen_Q[4].N3 * Q_KP2.N3 + l_eigen_Q[4].N4 * Q_KP2.N4;
    
    Qk_iminus_KP1.N0 = l_eigen_Q[0].N0 * Q_KP1.N0 + l_eigen_Q[0].N1 * Q_KP1.N1 + l_eigen_Q[0].N2 * Q_KP1.N2 + l_eigen_Q[0].N3 * Q_KP1.N3 + l_eigen_Q[0].N4 * Q_KP1.N4;
    Qk_iminus_KP1.N1 = l_eigen_Q[1].N0 * Q_KP1.N0 + l_eigen_Q[1].N1 * Q_KP1.N1 + l_eigen_Q[1].N2 * Q_KP1.N2 + l_eigen_Q[1].N3 * Q_KP1.N3 + l_eigen_Q[1].N4 * Q_KP1.N4;
    Qk_iminus_KP1.N2 = l_eigen_Q[2].N0 * Q_KP1.N0 + l_eigen_Q[2].N1 * Q_KP1.N1 + l_eigen_Q[2].N2 * Q_KP1.N2 + l_eigen_Q[2].N3 * Q_KP1.N3 + l_eigen_Q[2].N4 * Q_KP1.N4;
    Qk_iminus_KP1.N3 = l_eigen_Q[3].N0 * Q_KP1.N0 + l_eigen_Q[3].N1 * Q_KP1.N1 + l_eigen_Q[3].N2 * Q_KP1.N2 + l_eigen_Q[3].N3 * Q_KP1.N3 + l_eigen_Q[3].N4 * Q_KP1.N4;
    Qk_iminus_KP1.N4 = l_eigen_Q[4].N0 * Q_KP1.N0 + l_eigen_Q[4].N1 * Q_KP1.N1 + l_eigen_Q[4].N2 * Q_KP1.N2 + l_eigen_Q[4].N3 * Q_KP1.N3 + l_eigen_Q[4].N4 * Q_KP1.N4;
    
    Qk_iminus_I.N0 = l_eigen_Q[0].N0 * Q_I.N0 + l_eigen_Q[0].N1 * Q_I.N1 + l_eigen_Q[0].N2 * Q_I.N2 + l_eigen_Q[0].N3 * Q_I.N3 + l_eigen_Q[0].N4 * Q_I.N4;
    Qk_iminus_I.N1 = l_eigen_Q[1].N0 * Q_I.N0 + l_eigen_Q[1].N1 * Q_I.N1 + l_eigen_Q[1].N2 * Q_I.N2 + l_eigen_Q[1].N3 * Q_I.N3 + l_eigen_Q[1].N4 * Q_I.N4;
    Qk_iminus_I.N2 = l_eigen_Q[2].N0 * Q_I.N0 + l_eigen_Q[2].N1 * Q_I.N1 + l_eigen_Q[2].N2 * Q_I.N2 + l_eigen_Q[2].N3 * Q_I.N3 + l_eigen_Q[2].N4 * Q_I.N4;
    Qk_iminus_I.N3 = l_eigen_Q[3].N0 * Q_I.N0 + l_eigen_Q[3].N1 * Q_I.N1 + l_eigen_Q[3].N2 * Q_I.N2 + l_eigen_Q[3].N3 * Q_I.N3 + l_eigen_Q[3].N4 * Q_I.N4;
    Qk_iminus_I.N4 = l_eigen_Q[4].N0 * Q_I.N0 + l_eigen_Q[4].N1 * Q_I.N1 + l_eigen_Q[4].N2 * Q_I.N2 + l_eigen_Q[4].N3 * Q_I.N3 + l_eigen_Q[4].N4 * Q_I.N4;
    
    Qk_iminus_KM1.N0 = l_eigen_Q[0].N0 * Q_KM1.N0 + l_eigen_Q[0].N1 * Q_KM1.N1 + l_eigen_Q[0].N2 * Q_KM1.N2 + l_eigen_Q[0].N3 * Q_KM1.N3 + l_eigen_Q[0].N4 * Q_KM1.N4;
    Qk_iminus_KM1.N1 = l_eigen_Q[1].N0 * Q_KM1.N0 + l_eigen_Q[1].N1 * Q_KM1.N1 + l_eigen_Q[1].N2 * Q_KM1.N2 + l_eigen_Q[1].N3 * Q_KM1.N3 + l_eigen_Q[1].N4 * Q_KM1.N4;
    Qk_iminus_KM1.N2 = l_eigen_Q[2].N0 * Q_KM1.N0 + l_eigen_Q[2].N1 * Q_KM1.N1 + l_eigen_Q[2].N2 * Q_KM1.N2 + l_eigen_Q[2].N3 * Q_KM1.N3 + l_eigen_Q[2].N4 * Q_KM1.N4;
    Qk_iminus_KM1.N3 = l_eigen_Q[3].N0 * Q_KM1.N0 + l_eigen_Q[3].N1 * Q_KM1.N1 + l_eigen_Q[3].N2 * Q_KM1.N2 + l_eigen_Q[3].N3 * Q_KM1.N3 + l_eigen_Q[3].N4 * Q_KM1.N4;
    Qk_iminus_KM1.N4 = l_eigen_Q[4].N0 * Q_KM1.N0 + l_eigen_Q[4].N1 * Q_KM1.N1 + l_eigen_Q[4].N2 * Q_KM1.N2 + l_eigen_Q[4].N3 * Q_KM1.N3 + l_eigen_Q[4].N4 * Q_KM1.N4;
    
    Qk_iminus_KM2.N0 = l_eigen_Q[0].N0 * Q_KM2.N0 + l_eigen_Q[0].N1 * Q_KM2.N1 + l_eigen_Q[0].N2 * Q_KM2.N2 + l_eigen_Q[0].N3 * Q_KM2.N3 + l_eigen_Q[0].N4 * Q_KM2.N4;
    Qk_iminus_KM2.N1 = l_eigen_Q[1].N0 * Q_KM2.N0 + l_eigen_Q[1].N1 * Q_KM2.N1 + l_eigen_Q[1].N2 * Q_KM2.N2 + l_eigen_Q[1].N3 * Q_KM2.N3 + l_eigen_Q[1].N4 * Q_KM2.N4;
    Qk_iminus_KM2.N2 = l_eigen_Q[2].N0 * Q_KM2.N0 + l_eigen_Q[2].N1 * Q_KM2.N1 + l_eigen_Q[2].N2 * Q_KM2.N2 + l_eigen_Q[2].N3 * Q_KM2.N3 + l_eigen_Q[2].N4 * Q_KM2.N4;
    Qk_iminus_KM2.N3 = l_eigen_Q[3].N0 * Q_KM2.N0 + l_eigen_Q[3].N1 * Q_KM2.N1 + l_eigen_Q[3].N2 * Q_KM2.N2 + l_eigen_Q[3].N3 * Q_KM2.N3 + l_eigen_Q[3].N4 * Q_KM2.N4;
    Qk_iminus_KM2.N4 = l_eigen_Q[4].N0 * Q_KM2.N0 + l_eigen_Q[4].N1 * Q_KM2.N1 + l_eigen_Q[4].N2 * Q_KM2.N2 + l_eigen_Q[4].N3 * Q_KM2.N3 + l_eigen_Q[4].N4 * Q_KM2.N4;
    
    Qk_iminus_KM3.N0 = l_eigen_Q[0].N0 * Q_KM3.N0 + l_eigen_Q[0].N1 * Q_KM3.N1 + l_eigen_Q[0].N2 * Q_KM3.N2 + l_eigen_Q[0].N3 * Q_KM3.N3 + l_eigen_Q[0].N4 * Q_KM3.N4;
    Qk_iminus_KM3.N1 = l_eigen_Q[1].N0 * Q_KM3.N0 + l_eigen_Q[1].N1 * Q_KM3.N1 + l_eigen_Q[1].N2 * Q_KM3.N2 + l_eigen_Q[1].N3 * Q_KM3.N3 + l_eigen_Q[1].N4 * Q_KM3.N4;
    Qk_iminus_KM3.N2 = l_eigen_Q[2].N0 * Q_KM3.N0 + l_eigen_Q[2].N1 * Q_KM3.N1 + l_eigen_Q[2].N2 * Q_KM3.N2 + l_eigen_Q[2].N3 * Q_KM3.N3 + l_eigen_Q[2].N4 * Q_KM3.N4;
    Qk_iminus_KM3.N3 = l_eigen_Q[3].N0 * Q_KM3.N0 + l_eigen_Q[3].N1 * Q_KM3.N1 + l_eigen_Q[3].N2 * Q_KM3.N2 + l_eigen_Q[3].N3 * Q_KM3.N3 + l_eigen_Q[3].N4 * Q_KM3.N4;
    Qk_iminus_KM3.N4 = l_eigen_Q[4].N0 * Q_KM3.N0 + l_eigen_Q[4].N1 * Q_KM3.N1 + l_eigen_Q[4].N2 * Q_KM3.N2 + l_eigen_Q[4].N3 * Q_KM3.N3 + l_eigen_Q[4].N4 * Q_KM3.N4;
    
    /******************************************END OF CHARACTERISTIC PLANE TRANSFER******************************************************************************************/
    /*********************************************************************************************************************************************************************/

    /******************************Equation 1 ****************************************/
    /** Qi(i+1/2)_plus**/
    Qi_half_p[0] = (15.0 / 8.0)*Qi_iplus_IP1.N0 - (5.0 / 4.0)*Qi_iplus_IP2.N0 + (3.0 / 8.0)*Qi_iplus_IP3.N0;
    Qi_half_p[1] = (3.0 / 8.0)*Qi_iplus_I.N0 + (3.0 / 4.0)*Qi_iplus_IP1.N0 - (1.0 / 8.0)*Qi_iplus_IP2.N0;
    Qi_half_p[2] = (-1.0 / 8.0)*Qi_iplus_IM1.N0 + (3.0 / 4.0)*Qi_iplus_I.N0 + (3.0 / 8.0)*Qi_iplus_IP1.N0;

    /** Qi(i+1/2)_minus**/
    Qi_half_m[0] = (3.0 / 8.0)*Qi_iplus_I.N0 + (3.0 / 4.0)*Qi_iplus_IP1.N0 - (1.0 / 8.0)*Qi_iplus_IP2.N0;
    Qi_half_m[1] = (-1.0 / 8.0)*Qi_iplus_IM1.N0 + (3.0 / 4.0)*Qi_iplus_I.N0 + (3.0 / 8.0)*Qi_iplus_IP1.N0;
    Qi_half_m[2] = (3.0 / 8.0)*Qi_iplus_IM2.N0 - (5.0 / 4.0)*Qi_iplus_IM1.N0 + (15.0 / 8.0)*Qi_iplus_I.N0;

    /** Qi(i-1/2)_minus**/
    Qi_halfn_m[0] = (3.0 / 8.0)*Qi_iminus_IM1.N0 + (3.0 / 4.0)*Qi_iminus_I.N0 - (1.0 / 8.0)*Qi_iminus_IP1.N0;
    Qi_halfn_m[1] = (-1.0 / 8.0)*Qi_iminus_IM2.N0 + (3.0 / 4.0)*Qi_iminus_IM1.N0 + (3.0 / 8.0)*Qi_iminus_I.N0;
    Qi_halfn_m[2] = (3.0 / 8.0)*Qi_iminus_IM3.N0 - (5.0 / 4.0)*Qi_iminus_IM2.N0 + (15.0 / 8.0)*Qi_iminus_IM1.N0;

    /** Qi(i-1/2)_plus**/
    Qi_half_np[0] = (15.0 / 8.0)*Qi_iminus_I.N0 - (5.0 / 4.0)*Qi_iminus_IP1.N0 + (3.0 / 8.0)*Qi_iminus_IP2.N0;
    Qi_half_np[1] = (3.0 / 8.0)*Qi_iminus_IM1.N0 + (3.0 / 4.0)*Qi_iminus_I.N0 - (1.0 / 8.0)*Qi_iminus_IP1.N0;
    Qi_half_np[2] = (-1.0 / 8.0)*Qi_iminus_IM2.N0 + (3.0 / 4.0)*Qi_iminus_IM1.N0 + (3.0 / 8.0)*Qi_iminus_I.N0;

    /** Qj(i+1/2)_plus**/
    Qj_half_p[0] = (15.0 / 8.0)*Qj_iplus_JP1.N0 - (5.0 / 4.0)*Qj_iplus_JP2.N0 + (3.0 / 8.0)*Qj_iplus_JP3.N0;
    Qj_half_p[1] = (3.0 / 8.0)*Qj_iplus_I.N0 + (3.0 / 4.0)*Qj_iplus_JP1.N0 - (1.0 / 8.0)*Qj_iplus_JP2.N0;
    Qj_half_p[2] = (-1.0 / 8.0)*Qj_iplus_JM1.N0 + (3.0 / 4.0)*Qj_iplus_I.N0 + (3.0 / 8.0)*Qj_iplus_JP1.N0;

    /** Qj(i+1/2)_minus**/
    Qj_half_m[0] = (3.0 / 8.0)*Qj_iplus_I.N0 + (3.0 / 4.0)*Qj_iplus_JP1.N0 - (1.0 / 8.0)*Qj_iplus_JP2.N0;
    Qj_half_m[1] = (-1.0 / 8.0)*Qj_iplus_JM1.N0 + (3.0 / 4.0)*Qj_iplus_I.N0 + (3.0 / 8.0)*Qj_iplus_JP1.N0;
    Qj_half_m[2] = (3.0 / 8.0)*Qj_iplus_JM2.N0 - (5.0 / 4.0)*Qj_iplus_JM1.N0 + (15.0 / 8.0)*Qj_iplus_I.N0;

    /** Qj(i-1/2)_minus**/
    Qj_halfn_m[0] = (3.0 / 8.0)*Qj_iminus_JM1.N0 + (3.0 / 4.0)*Qj_iminus_I.N0 - (1.0 / 8.0)*Qj_iminus_JP1.N0;
    Qj_halfn_m[1] = (-1.0 / 8.0)*Qj_iminus_JM2.N0 + (3.0 / 4.0)*Qj_iminus_JM1.N0 + (3.0 / 8.0)*Qj_iminus_I.N0;
    Qj_halfn_m[2] = (3.0 / 8.0)*Qj_iminus_JM3.N0 - (5.0 / 4.0)*Qj_iminus_JM2.N0 + (15.0 / 8.0)*Qj_iminus_JM1.N0;

    /** Qj(i-1/2)_plus**/
    Qj_half_np[0] = (15.0 / 8.0)*Qj_iminus_I.N0 - (5.0 / 4.0)*Qj_iminus_JP1.N0 + (3.0 / 8.0)*Qj_iminus_JP2.N0;
    Qj_half_np[1] = (3.0 / 8.0)*Qj_iminus_JM1.N0 + (3.0 / 4.0)*Qj_iminus_I.N0 - (1.0 / 8.0)*Qj_iminus_JP1.N0;
    Qj_half_np[2] = (-1.0 / 8.0)*Qj_iminus_JM2.N0 + (3.0 / 4.0)*Qj_iminus_JM1.N0 + (3.0 / 8.0)*Qj_iminus_I.N0;

    /** Qk(i+1/2)_plus**/
    Qk_half_p[0] = (15.0 / 8.0)*Qk_iplus_KP1.N0 - (5.0 / 4.0)*Qk_iplus_KP2.N0 + (3.0 / 8.0)*Qk_iplus_KP3.N0;
    Qk_half_p[1] = (3.0 / 8.0)*Qk_iplus_I.N0 + (3.0 / 4.0)*Qk_iplus_KP1.N0 - (1.0 / 8.0)*Qk_iplus_KP2.N0;
    Qk_half_p[2] = (-1.0 / 8.0)*Qk_iplus_KM1.N0 + (3.0 / 4.0)*Qk_iplus_I.N0 + (3.0 / 8.0)*Qk_iplus_KP1.N0;

    /** Qk(i+1/2)_minus**/
    Qk_half_m[0] = (3.0 / 8.0)*Qk_iplus_I.N0 + (3.0 / 4.0)*Qk_iplus_KP1.N0 - (1.0 / 8.0)*Qk_iplus_KP2.N0;
    Qk_half_m[1] = (-1.0 / 8.0)*Qk_iplus_KM1.N0 + (3.0 / 4.0)*Qk_iplus_I.N0 + (3.0 / 8.0)*Qk_iplus_KP1.N0;
    Qk_half_m[2] = (3.0 / 8.0)*Qk_iplus_KM2.N0 - (5.0 / 4.0)*Qk_iplus_KM1.N0 + (15.0 / 8.0)*Qk_iplus_I.N0;

    /** Qk(i-1/2)_minus**/
    Qk_halfn_m[0] = (3.0 / 8.0)*Qk_iminus_KM1.N0 + (3.0 / 4.0)*Qk_iminus_I.N0 - (1.0 / 8.0)*Qk_iminus_KP1.N0;
    Qk_halfn_m[1] = (-1.0 / 8.0)*Qk_iminus_KM2.N0 + (3.0 / 4.0)*Qk_iminus_KM1.N0 + (3.0 / 8.0)*Qk_iminus_I.N0;
    Qk_halfn_m[2] = (3.0 / 8.0)*Qk_iminus_KM3.N0 - (5.0 / 4.0)*Qk_iminus_KM2.N0 + (15.0 / 8.0)*Qk_iminus_KM1.N0;

    /** Qk(i-1/2)_plus**/
    Qk_half_np[0] = (15.0 / 8.0)*Qk_iminus_I.N0 - (5.0 / 4.0)*Qk_iminus_KP1.N0 + (3.0 / 8.0)*Qk_iminus_KP2.N0;
    Qk_half_np[1] = (3.0 / 8.0)*Qk_iminus_KM1.N0 + (3.0 / 4.0)*Qk_iminus_I.N0 - (1.0 / 8.0)*Qk_iminus_KP1.N0;
    Qk_half_np[2] = (-1.0 / 8.0)*Qk_iminus_KM2.N0 + (3.0 / 4.0)*Qk_iminus_KM1.N0 + (3.0 / 8.0)*Qk_iminus_I.N0;

    /**********************************SMOOTHNESS INDICATOR OR SMOOTHNESS FUNCTION********************************/

    IS_Qim[2] = (13.0 / 12.0)*(pow(Qi_iplus_IM2.N0 - 2.0*Qi_iplus_IM1.N0 + Qi_iplus_I.N0, 2)) + (1.0 / 4.0)*(pow(Qi_iplus_IM2.N0 - 4.0*Qi_iplus_IM1.N0 + 3.0*Qi_iplus_I.N0, 2));
    IS_Qim[1] = (13.0 / 12.0)*(pow(Qi_iplus_IM1.N0 - 2.0*Qi_iplus_I.N0 + Qi_iplus_IP1.N0, 2)) + (1.0 / 4.0)*(pow(Qi_iplus_IM1.N0 - Qi_iplus_IP1.N0, 2));
    IS_Qim[0] = (13.0 / 12.0)*(pow(Qi_iplus_I.N0 - 2.0*Qi_iplus_IP1.N0 + Qi_iplus_IP2.N0, 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iplus_I.N0 - 4.0*Qi_iplus_IP1.N0 + Qi_iplus_IP2.N0, 2));

    IS_Qip[2] = (13.0 / 12.0)*(pow(Qi_iplus_IM1.N0 - 2.0*Qi_iplus_I.N0 + Qi_iplus_IP1.N0, 2)) + (1.0 / 4.0)*(pow(Qi_iplus_IM1.N0 - 4.0*Qi_iplus_I.N0 + 3.0*Qi_iplus_IP1.N0, 2));
    IS_Qip[1] = (13.0 / 12.0)*(pow(Qi_iplus_I.N0 - 2.0*Qi_iplus_IP1.N0 + Qi_iplus_IP2.N0, 2)) + (1.0 / 4.0)*(pow(Qi_iplus_I.N0 - Qi_iplus_IP2.N0, 2));
    IS_Qip[0] = (13.0 / 12.0)*(pow(Qi_iplus_IP1.N0 - 2.0*Qi_iplus_IP2.N0 + Qi_iplus_IP3.N0, 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iplus_IP1.N0 - 4.0*Qi_iplus_IP2.N0 + Qi_iplus_IP3.N0, 2));

    IS_Qinm[2] = (13.0 / 12.0)*(pow(Qi_iminus_IM3.N0 - 2.0*Qi_iminus_IM2.N0 + Qi_iminus_IM1.N0, 2)) + (1.0 / 4.0)*(pow(Qi_iminus_IM3.N0 - 4.0*Qi_iminus_IM2.N0 + 3.0*Qi_iminus_IM1.N0, 2));
    IS_Qinm[1] = (13.0 / 12.0)*(pow(Qi_iminus_IM2.N0 - 2.0*Qi_iminus_IM1.N0 + Qi_iminus_I.N0, 2)) + (1.0 / 4.0)*(pow(Qi_iminus_IM2.N0 - Qi_iminus_I.N0, 2));
    IS_Qinm[0] = (13.0 / 12.0)*(pow(Qi_iminus_IM1.N0 - 2.0*Qi_iminus_I.N0 + Qi_iminus_IP1.N0, 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iminus_IM1.N0 - 4.0*Qi_iminus_I.N0 + Qi_iminus_IP1.N0, 2));

    IS_Qinp[2] = (13.0 / 12.0)*(pow(Qi_iminus_IM2.N0 - 2.0*Qi_iminus_IM1.N0 + Qi_iminus_I.N0, 2)) + (1.0 / 4.0)*(pow(Qi_iminus_IM2.N0 - 4.0*Qi_iminus_IM1.N0 + 3.0*Qi_iminus_I.N0, 2));
    IS_Qinp[1] = (13.0 / 12.0)*(pow(Qi_iminus_IM1.N0 - 2.0*Qi_iminus_I.N0 + Qi_iminus_IP1.N0, 2)) + (1.0 / 4.0)*(pow(Qi_iminus_IM1.N0 - Qi_iminus_IP1.N0, 2));
    IS_Qinp[0] = (13.0 / 12.0)*(pow(Qi_iminus_I.N0 - 2.0*Qi_iminus_IP1.N0 + Qi_iminus_IP2.N0, 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iminus_I.N0 - 4.0*Qi_iminus_IP1.N0 + Qi_iminus_IP2.N0, 2));

    IS_Qjm[2] = (13.0 / 12.0)*(pow(Qj_iplus_JM2.N0 - 2.0*Qj_iplus_JM1.N0 + Qj_iplus_I.N0, 2)) + (1.0 / 4.0)*(pow(Qj_iplus_JM2.N0 - 4.0*Qj_iplus_JM1.N0 + 3.0*Qj_iplus_I.N0, 2));
    IS_Qjm[1] = (13.0 / 12.0)*(pow(Qj_iplus_JM1.N0 - 2.0*Qj_iplus_I.N0 + Qj_iplus_JP1.N0, 2)) + (1.0 / 4.0)*(pow(Qj_iplus_JM1.N0 - Qj_iplus_JP1.N0, 2));
    IS_Qjm[0] = (13.0 / 12.0)*(pow(Qj_iplus_I.N0 - 2.0*Qj_iplus_JP1.N0 + Qj_iplus_JP2.N0, 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iplus_I.N0 - 4.0*Qj_iplus_JP1.N0 + Qj_iplus_JP2.N0, 2));

    IS_Qjp[2] = (13.0 / 12.0)*(pow(Qj_iplus_JM1.N0 - 2.0*Qj_iplus_I.N0 + Qj_iplus_JP1.N0, 2)) + (1.0 / 4.0)*(pow(Qj_iplus_JM1.N0 - 4.0*Qj_iplus_I.N0 + 3.0*Qj_iplus_JP1.N0, 2));
    IS_Qjp[1] = (13.0 / 12.0)*(pow(Qj_iplus_I.N0 - 2.0*Qj_iplus_JP1.N0 + Qj_iplus_JP2.N0, 2)) + (1.0 / 4.0)*(pow(Qj_iplus_I.N0 - Qj_iplus_JP2.N0, 2));
    IS_Qjp[0] = (13.0 / 12.0)*(pow(Qj_iplus_JP1.N0 - 2.0*Qj_iplus_JP2.N0 + Qj_iplus_JP3.N0, 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iplus_JP1.N0 - 4.0*Qj_iplus_JP2.N0 + Qj_iplus_JP3.N0, 2));

    IS_Qjnm[2] = (13.0 / 12.0)*(pow(Qj_iminus_JM3.N0 - 2.0*Qj_iminus_JM2.N0 + Qj_iminus_JM1.N0, 2)) + (1.0 / 4.0)*(pow(Qj_iminus_JM3.N0 - 4.0*Qj_iminus_JM2.N0 + 3.0*Qj_iminus_JM1.N0, 2));
    IS_Qjnm[1] = (13.0 / 12.0)*(pow(Qj_iminus_JM2.N0 - 2.0*Qj_iminus_JM1.N0 + Qj_iminus_I.N0, 2)) + (1.0 / 4.0)*(pow(Qj_iminus_JM2.N0 - Qj_iminus_I.N0, 2));
    IS_Qjnm[0] = (13.0 / 12.0)*(pow(Qj_iminus_JM1.N0 - 2.0*Qj_iminus_I.N0 + Qj_iminus_JP1.N0, 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iminus_JM1.N0 - 4.0*Qj_iminus_I.N0 + Qj_iminus_JP1.N0, 2));

    IS_Qjnp[2] = (13.0 / 12.0)*(pow(Qj_iminus_JM2.N0 - 2.0*Qj_iminus_JM1.N0 + Qj_iminus_I.N0, 2)) + (1.0 / 4.0)*(pow(Qj_iminus_JM2.N0 - 4.0*Qj_iminus_JM1.N0 + 3.0*Qj_iminus_I.N0, 2));
    IS_Qjnp[1] = (13.0 / 12.0)*(pow(Qj_iminus_JM1.N0 - 2.0*Qj_iminus_I.N0 + Qj_iminus_JP1.N0, 2)) + (1.0 / 4.0)*(pow(Qj_iminus_JM1.N0 - Qj_iminus_JP1.N0, 2));
    IS_Qjnp[0] = (13.0 / 12.0)*(pow(Qj_iminus_I.N0 - 2.0*Qj_iminus_JP1.N0 + Qj_iminus_JP2.N0, 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iminus_I.N0 - 4.0*Qj_iminus_JP1.N0 + Qj_iminus_JP2.N0, 2));

    IS_Qkm[2] = (13.0 / 12.0)*(pow(Qk_iplus_KM2.N0 - 2.0*Qk_iplus_KM1.N0 + Qk_iplus_I.N0, 2)) + (1.0 / 4.0)*(pow(Qk_iplus_KM2.N0 - 4.0*Qk_iplus_KM1.N0 + 3.0*Qk_iplus_I.N0, 2));
    IS_Qkm[1] = (13.0 / 12.0)*(pow(Qk_iplus_KM1.N0 - 2.0*Qk_iplus_I.N0 + Qk_iplus_KP1.N0, 2)) + (1.0 / 4.0)*(pow(Qk_iplus_KM1.N0 - Qk_iplus_KP1.N0, 2));
    IS_Qkm[0] = (13.0 / 12.0)*(pow(Qk_iplus_I.N0 - 2.0*Qk_iplus_KP1.N0 + Qk_iplus_KP2.N0, 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iplus_I.N0 - 4.0*Qk_iplus_KP1.N0 + Qk_iplus_KP2.N0, 2));

    IS_Qkp[2] = (13.0 / 12.0)*(pow(Qk_iplus_KM1.N0 - 2.0*Qk_iplus_I.N0 + Qk_iplus_KP1.N0, 2)) + (1.0 / 4.0)*(pow(Qk_iplus_KM1.N0 - 4.0*Qk_iplus_I.N0 + 3.0*Qk_iplus_KP1.N0, 2));
    IS_Qkp[1] = (13.0 / 12.0)*(pow(Qk_iplus_I.N0 - 2.0*Qk_iplus_KP1.N0 + Qk_iplus_KP2.N0, 2)) + (1.0 / 4.0)*(pow(Qk_iplus_I.N0 - Qk_iplus_KP2.N0, 2));
    IS_Qkp[0] = (13.0 / 12.0)*(pow(Qk_iplus_KP1.N0 - 2.0*Qk_iplus_KP2.N0 + Qk_iplus_KP3.N0, 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iplus_KP1.N0 - 4.0*Qk_iplus_KP2.N0 + Qk_iplus_KP3.N0, 2));

    IS_Qknm[2] = (13.0 / 12.0)*(pow(Qk_iminus_KM3.N0 - 2.0*Qk_iminus_KM2.N0 + Qk_iminus_KM1.N0, 2)) + (1.0 / 4.0)*(pow(Qk_iminus_KM3.N0 - 4.0*Qk_iminus_KM2.N0 + 3.0*Qk_iminus_KM1.N0, 2));
    IS_Qknm[1] = (13.0 / 12.0)*(pow(Qk_iminus_KM2.N0 - 2.0*Qk_iminus_KM1.N0 + Qk_iminus_I.N0, 2)) + (1.0 / 4.0)*(pow(Qk_iminus_KM2.N0 - Qk_iminus_I.N0, 2));
    IS_Qknm[0] = (13.0 / 12.0)*(pow(Qk_iminus_KM1.N0 - 2.0*Qk_iminus_I.N0 + Qk_iminus_KP1.N0, 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iminus_KM1.N0 - 4.0*Qk_iminus_I.N0 + Qk_iminus_KP1.N0, 2));

    IS_Qknp[2] = (13.0 / 12.0)*(pow(Qk_iminus_KM2.N0 - 2.0*Qk_iminus_KM1.N0 + Qk_iminus_I.N0, 2)) + (1.0 / 4.0)*(pow(Qk_iminus_KM2.N0 - 4.0*Qk_iminus_KM1.N0 + 3.0*Qk_iminus_I.N0, 2));
    IS_Qknp[1] = (13.0 / 12.0)*(pow(Qk_iminus_KM1.N0 - 2.0*Qk_iminus_I.N0 + Qk_iminus_KP1.N0, 2)) + (1.0 / 4.0)*(pow(Qk_iminus_KM1.N0 - Qk_iminus_KP1.N0, 2));
    IS_Qknp[0] = (13.0 / 12.0)*(pow(Qk_iminus_I.N0 - 2.0*Qk_iminus_KP1.N0 + Qk_iminus_KP2.N0, 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iminus_I.N0 - 4.0*Qk_iminus_KP1.N0 + Qk_iminus_KP2.N0, 2));

    w_Qip[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qip[0]), 2.0));
    w_Qip[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qip[1]), 2.0));
    w_Qip[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qip[2]), 2.0));

    w_Qim[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qim[0]), 2.0));
    w_Qim[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qim[1]), 2.0));
    w_Qim[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qim[2]), 2.0));

    w_Qinp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qinp[0]), 2.0));
    w_Qinp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qinp[1]), 2.0));
    w_Qinp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qinp[2]), 2.0));

    w_Qinm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qinm[0]), 2.0));
    w_Qinm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qinm[1]), 2.0));
    w_Qinm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qinm[2]), 2.0));

    w_Qjp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjp[0]), 2.0));
    w_Qjp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjp[1]), 2.0));
    w_Qjp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjp[2]), 2.0));

    w_Qjm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjm[0]), 2.0));
    w_Qjm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjm[1]), 2.0));
    w_Qjm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjm[2]), 2.0));

    w_Qjnp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjnp[0]), 2.0));
    w_Qjnp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjnp[1]), 2.0));
    w_Qjnp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjnp[2]), 2.0));

    w_Qjnm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjnm[0]), 2.0));
    w_Qjnm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjnm[1]), 2.0));
    w_Qjnm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjnm[2]), 2.0));

    w_Qkp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qkp[0]), 2.0));
    w_Qkp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qkp[1]), 2.0));
    w_Qkp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qkp[2]), 2.0));

    w_Qkm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qkm[0]), 2.0));
    w_Qkm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qkm[1]), 2.0));
    w_Qkm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qkm[2]), 2.0));

    w_Qknp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qknp[0]), 2.0));
    w_Qknp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qknp[1]), 2.0));
    w_Qknp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qknp[2]), 2.0));

    w_Qknm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qknm[0]), 2.0));
    w_Qknm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qknm[1]), 2.0));
    w_Qknm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qknm[2]), 2.0));

    W_Qip[0] = w_Qip[0] / (w_Qip[0] + w_Qip[1] + w_Qip[2]);
    W_Qip[1] = w_Qip[1] / (w_Qip[0] + w_Qip[1] + w_Qip[2]);
    W_Qip[2] = w_Qip[2] / (w_Qip[0] + w_Qip[1] + w_Qip[2]);

    W_Qinp[0] = w_Qinp[0] / (w_Qinp[0] + w_Qinp[1] + w_Qinp[2]);
    W_Qinp[1] = w_Qinp[1] / (w_Qinp[0] + w_Qinp[1] + w_Qinp[2]);
    W_Qinp[2] = w_Qinp[2] / (w_Qinp[0] + w_Qinp[1] + w_Qinp[2]);

    W_Qim[0] = w_Qim[0] / (w_Qim[0] + w_Qim[1] + w_Qim[2]);
    W_Qim[1] = w_Qim[1] / (w_Qim[0] + w_Qim[1] + w_Qim[2]);
    W_Qim[2] = w_Qim[2] / (w_Qim[0] + w_Qim[1] + w_Qim[2]);

    W_Qinm[0] = w_Qinm[0] / (w_Qinm[0] + w_Qinm[1] + w_Qinm[2]);
    W_Qinm[1] = w_Qinm[1] / (w_Qinm[0] + w_Qinm[1] + w_Qinm[2]);
    W_Qinm[2] = w_Qinm[2] / (w_Qinm[0] + w_Qinm[1] + w_Qinm[2]);

    W_Qjp[0] = w_Qjp[0] / (w_Qjp[0] + w_Qjp[1] + w_Qjp[2]);
    W_Qjp[1] = w_Qjp[1] / (w_Qjp[0] + w_Qjp[1] + w_Qjp[2]);
    W_Qjp[2] = w_Qjp[2] / (w_Qjp[0] + w_Qjp[1] + w_Qjp[2]);

    W_Qjnp[0] = w_Qjnp[0] / (w_Qjnp[0] + w_Qjnp[1] + w_Qjnp[2]);
    W_Qjnp[1] = w_Qjnp[1] / (w_Qjnp[0] + w_Qjnp[1] + w_Qjnp[2]);
    W_Qjnp[2] = w_Qjnp[2] / (w_Qjnp[0] + w_Qjnp[1] + w_Qjnp[2]);

    W_Qjm[0] = w_Qjm[0] / (w_Qjm[0] + w_Qjm[1] + w_Qjm[2]);
    W_Qjm[1] = w_Qjm[1] / (w_Qjm[0] + w_Qjm[1] + w_Qjm[2]);
    W_Qjm[2] = w_Qjm[2] / (w_Qjm[0] + w_Qjm[1] + w_Qjm[2]);

    W_Qjnm[0] = w_Qjnm[0] / (w_Qjnm[0] + w_Qjnm[1] + w_Qjnm[2]);
    W_Qjnm[1] = w_Qjnm[1] / (w_Qjnm[0] + w_Qjnm[1] + w_Qjnm[2]);
    W_Qjnm[2] = w_Qjnm[2] / (w_Qjnm[0] + w_Qjnm[1] + w_Qjnm[2]);

    W_Qkp[0] = w_Qkp[0] / (w_Qkp[0] + w_Qkp[1] + w_Qkp[2]);
    W_Qkp[1] = w_Qkp[1] / (w_Qkp[0] + w_Qkp[1] + w_Qkp[2]);
    W_Qkp[2] = w_Qkp[2] / (w_Qkp[0] + w_Qkp[1] + w_Qkp[2]);

    W_Qknp[0] = w_Qknp[0] / (w_Qknp[0] + w_Qknp[1] + w_Qknp[2]);
    W_Qknp[1] = w_Qknp[1] / (w_Qknp[0] + w_Qknp[1] + w_Qknp[2]);
    W_Qknp[2] = w_Qknp[2] / (w_Qknp[0] + w_Qknp[1] + w_Qknp[2]);

    W_Qkm[0] = w_Qkm[0] / (w_Qkm[0] + w_Qkm[1] + w_Qkm[2]);
    W_Qkm[1] = w_Qkm[1] / (w_Qkm[0] + w_Qkm[1] + w_Qkm[2]);
    W_Qkm[2] = w_Qkm[2] / (w_Qkm[0] + w_Qkm[1] + w_Qkm[2]);

    W_Qknm[0] = w_Qknm[0] / (w_Qknm[0] + w_Qknm[1] + w_Qknm[2]);
    W_Qknm[1] = w_Qknm[1] / (w_Qknm[0] + w_Qknm[1] + w_Qknm[2]);
    W_Qknm[2] = w_Qknm[2] / (w_Qknm[0] + w_Qknm[1] + w_Qknm[2]);

    Qi_iplus_half_pos_char[0] = W_Qip[0] * Qi_half_p[0] + W_Qip[1] * Qi_half_p[1] + W_Qip[2] * Qi_half_p[2];
    Qi_iminus_half_pos_char[0] = W_Qinp[0] * Qi_half_np[0] + W_Qinp[1] * Qi_half_np[1] + W_Qinp[2] * Qi_half_np[2];

    Qi_iplus_half_neg_char[0] = W_Qim[0] * Qi_half_m[0] + W_Qim[1] * Qi_half_m[1] + W_Qim[2] * Qi_half_m[2];
    Qi_iminus_half_neg_char[0] = W_Qinm[0] * Qi_halfn_m[0] + W_Qinm[1] * Qi_halfn_m[1] + W_Qinm[2] * Qi_halfn_m[2];

    Qj_iplus_half_pos_char[0] = W_Qjp[0] * Qj_half_p[0] + W_Qjp[1] * Qj_half_p[1] + W_Qjp[2] * Qj_half_p[2];
    Qj_iminus_half_pos_char[0] = W_Qjnp[0] * Qj_half_np[0] + W_Qjnp[1] * Qj_half_np[1] + W_Qjnp[2] * Qj_half_np[2];

    Qj_iplus_half_neg_char[0] = W_Qjm[0] * Qj_half_m[0] + W_Qjm[1] * Qj_half_m[1] + W_Qjm[2] * Qj_half_m[2];
    Qj_iminus_half_neg_char[0] = W_Qjnm[0] * Qj_halfn_m[0] + W_Qjnm[1] * Qj_halfn_m[1] + W_Qjnm[2] * Qj_halfn_m[2];

    Qk_iplus_half_pos_char[0] = W_Qkp[0] * Qk_half_p[0] + W_Qkp[1] * Qk_half_p[1] + W_Qkp[2] * Qk_half_p[2];
    Qk_iminus_half_pos_char[0] = W_Qknp[0] * Qk_half_np[0] + W_Qknp[1] * Qk_half_np[1] + W_Qknp[2] * Qk_half_np[2];

    Qk_iplus_half_neg_char[0] = W_Qkm[0] * Qk_half_m[0] + W_Qkm[1] * Qk_half_m[1] + W_Qkm[2] * Qk_half_m[2];
    Qk_iminus_half_neg_char[0] = W_Qknm[0] * Qk_halfn_m[0] + W_Qknm[1] * Qk_halfn_m[1] + W_Qknm[2] * Qk_halfn_m[2];


    /******************************Equation 2 ****************************************/
    /** Qi(i+1/2)_plus**/
    Qi_half_p[0] = (15.0 / 8.0)*Qi_iplus_IP1.N1 - (5.0 / 4.0)*Qi_iplus_IP2.N1 + (3.0 / 8.0)*Qi_iplus_IP3.N1;
    Qi_half_p[1] = (3.0 / 8.0)*Qi_iplus_I.N1 + (3.0 / 4.0)*Qi_iplus_IP1.N1 - (1.0 / 8.0)*Qi_iplus_IP2.N1;
    Qi_half_p[2] = (-1.0 / 8.0)*Qi_iplus_IM1.N1 + (3.0 / 4.0)*Qi_iplus_I.N1 + (3.0 / 8.0)*Qi_iplus_IP1.N1;

    /** Qi(i+1/2)_minus**/
    Qi_half_m[0] = (3.0 / 8.0)*Qi_iplus_I.N1 + (3.0 / 4.0)*Qi_iplus_IP1.N1 - (1.0 / 8.0)*Qi_iplus_IP2.N1;
    Qi_half_m[1] = (-1.0 / 8.0)*Qi_iplus_IM1.N1 + (3.0 / 4.0)*Qi_iplus_I.N1 + (3.0 / 8.0)*Qi_iplus_IP1.N1;
    Qi_half_m[2] = (3.0 / 8.0)*Qi_iplus_IM2.N1 - (5.0 / 4.0)*Qi_iplus_IM1.N1 + (15.0 / 8.0)*Qi_iplus_I.N1;

    /** Qi(i-1/2)_minus**/
    Qi_halfn_m[0] = (3.0 / 8.0)*Qi_iminus_IM1.N1 + (3.0 / 4.0)*Qi_iminus_I.N1 - (1.0 / 8.0)*Qi_iminus_IP1.N1;
    Qi_halfn_m[1] = (-1.0 / 8.0)*Qi_iminus_IM2.N1 + (3.0 / 4.0)*Qi_iminus_IM1.N1 + (3.0 / 8.0)*Qi_iminus_I.N1;
    Qi_halfn_m[2] = (3.0 / 8.0)*Qi_iminus_IM3.N1 - (5.0 / 4.0)*Qi_iminus_IM2.N1 + (15.0 / 8.0)*Qi_iminus_IM1.N1;

    /** Qi(i-1/2)_plus**/
    Qi_half_np[0] = (15.0 / 8.0)*Qi_iminus_I.N1 - (5.0 / 4.0)*Qi_iminus_IP1.N1 + (3.0 / 8.0)*Qi_iminus_IP2.N1;
    Qi_half_np[1] = (3.0 / 8.0)*Qi_iminus_IM1.N1 + (3.0 / 4.0)*Qi_iminus_I.N1 - (1.0 / 8.0)*Qi_iminus_IP1.N1;
    Qi_half_np[2] = (-1.0 / 8.0)*Qi_iminus_IM2.N1 + (3.0 / 4.0)*Qi_iminus_IM1.N1 + (3.0 / 8.0)*Qi_iminus_I.N1;

    /** Qj(i+1/2)_plus**/
    Qj_half_p[0] = (15.0 / 8.0)*Qj_iplus_JP1.N1 - (5.0 / 4.0)*Qj_iplus_JP2.N1 + (3.0 / 8.0)*Qj_iplus_JP3.N1;
    Qj_half_p[1] = (3.0 / 8.0)*Qj_iplus_I.N1 + (3.0 / 4.0)*Qj_iplus_JP1.N1 - (1.0 / 8.0)*Qj_iplus_JP2.N1;
    Qj_half_p[2] = (-1.0 / 8.0)*Qj_iplus_JM1.N1 + (3.0 / 4.0)*Qj_iplus_I.N1 + (3.0 / 8.0)*Qj_iplus_JP1.N1;

    /** Qj(i+1/2)_minus**/
    Qj_half_m[0] = (3.0 / 8.0)*Qj_iplus_I.N1 + (3.0 / 4.0)*Qj_iplus_JP1.N1 - (1.0 / 8.0)*Qj_iplus_JP2.N1;
    Qj_half_m[1] = (-1.0 / 8.0)*Qj_iplus_JM1.N1 + (3.0 / 4.0)*Qj_iplus_I.N1 + (3.0 / 8.0)*Qj_iplus_JP1.N1;
    Qj_half_m[2] = (3.0 / 8.0)*Qj_iplus_JM2.N1 - (5.0 / 4.0)*Qj_iplus_JM1.N1 + (15.0 / 8.0)*Qj_iplus_I.N1;

    /** Qj(i-1/2)_minus**/
    Qj_halfn_m[0] = (3.0 / 8.0)*Qj_iminus_JM1.N1 + (3.0 / 4.0)*Qj_iminus_I.N1 - (1.0 / 8.0)*Qj_iminus_JP1.N1;
    Qj_halfn_m[1] = (-1.0 / 8.0)*Qj_iminus_JM2.N1 + (3.0 / 4.0)*Qj_iminus_JM1.N1 + (3.0 / 8.0)*Qj_iminus_I.N1;
    Qj_halfn_m[2] = (3.0 / 8.0)*Qj_iminus_JM3.N1 - (5.0 / 4.0)*Qj_iminus_JM2.N1 + (15.0 / 8.0)*Qj_iminus_JM1.N1;

    /** Qj(i-1/2)_plus**/
    Qj_half_np[0] = (15.0 / 8.0)*Qj_iminus_I.N1 - (5.0 / 4.0)*Qj_iminus_JP1.N1 + (3.0 / 8.0)*Qj_iminus_JP2.N1;
    Qj_half_np[1] = (3.0 / 8.0)*Qj_iminus_JM1.N1 + (3.0 / 4.0)*Qj_iminus_I.N1 - (1.0 / 8.0)*Qj_iminus_JP1.N1;
    Qj_half_np[2] = (-1.0 / 8.0)*Qj_iminus_JM2.N1 + (3.0 / 4.0)*Qj_iminus_JM1.N1 + (3.0 / 8.0)*Qj_iminus_I.N1;

    /** Qk(i+1/2)_plus**/
    Qk_half_p[0] = (15.0 / 8.0)*Qk_iplus_KP1.N1 - (5.0 / 4.0)*Qk_iplus_KP2.N1 + (3.0 / 8.0)*Qk_iplus_KP3.N1;
    Qk_half_p[1] = (3.0 / 8.0)*Qk_iplus_I.N1 + (3.0 / 4.0)*Qk_iplus_KP1.N1 - (1.0 / 8.0)*Qk_iplus_KP2.N1;
    Qk_half_p[2] = (-1.0 / 8.0)*Qk_iplus_KM1.N1 + (3.0 / 4.0)*Qk_iplus_I.N1 + (3.0 / 8.0)*Qk_iplus_KP1.N1;

    /** Qk(i+1/2)_minus**/
    Qk_half_m[0] = (3.0 / 8.0)*Qk_iplus_I.N1 + (3.0 / 4.0)*Qk_iplus_KP1.N1 - (1.0 / 8.0)*Qk_iplus_KP2.N1;
    Qk_half_m[1] = (-1.0 / 8.0)*Qk_iplus_KM1.N1 + (3.0 / 4.0)*Qk_iplus_I.N1 + (3.0 / 8.0)*Qk_iplus_KP1.N1;
    Qk_half_m[2] = (3.0 / 8.0)*Qk_iplus_KM2.N1 - (5.0 / 4.0)*Qk_iplus_KM1.N1 + (15.0 / 8.0)*Qk_iplus_I.N1;

    /** Qk(i-1/2)_minus**/
    Qk_halfn_m[0] = (3.0 / 8.0)*Qk_iminus_KM1.N1 + (3.0 / 4.0)*Qk_iminus_I.N1 - (1.0 / 8.0)*Qk_iminus_KP1.N1;
    Qk_halfn_m[1] = (-1.0 / 8.0)*Qk_iminus_KM2.N1 + (3.0 / 4.0)*Qk_iminus_KM1.N1 + (3.0 / 8.0)*Qk_iminus_I.N1;
    Qk_halfn_m[2] = (3.0 / 8.0)*Qk_iminus_KM3.N1 - (5.0 / 4.0)*Qk_iminus_KM2.N1 + (15.0 / 8.0)*Qk_iminus_KM1.N1;

    /** Qk(i-1/2)_plus**/
    Qk_half_np[0] = (15.0 / 8.0)*Qk_iminus_I.N1 - (5.0 / 4.0)*Qk_iminus_KP1.N1 + (3.0 / 8.0)*Qk_iminus_KP2.N1;
    Qk_half_np[1] = (3.0 / 8.0)*Qk_iminus_KM1.N1 + (3.0 / 4.0)*Qk_iminus_I.N1 - (1.0 / 8.0)*Qk_iminus_KP1.N1;
    Qk_half_np[2] = (-1.0 / 8.0)*Qk_iminus_KM2.N1 + (3.0 / 4.0)*Qk_iminus_KM1.N1 + (3.0 / 8.0)*Qk_iminus_I.N1;

    /**********************************SMOOTHNESS INDICATOR OR SMOOTHNESS FUNCTION********************************/

    IS_Qim[2] = (13.0 / 12.0)*(pow(Qi_iplus_IM2.N1 - 2.0*Qi_iplus_IM1.N1 + Qi_iplus_I.N1, 2)) + (1.0 / 4.0)*(pow(Qi_iplus_IM2.N1 - 4.0*Qi_iplus_IM1.N1 + 3.0*Qi_iplus_I.N1, 2));
    IS_Qim[1] = (13.0 / 12.0)*(pow(Qi_iplus_IM1.N1 - 2.0*Qi_iplus_I.N1 + Qi_iplus_IP1.N1, 2)) + (1.0 / 4.0)*(pow(Qi_iplus_IM1.N1 - Qi_iplus_IP1.N1, 2));
    IS_Qim[0] = (13.0 / 12.0)*(pow(Qi_iplus_I.N1 - 2.0*Qi_iplus_IP1.N1 + Qi_iplus_IP2.N1, 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iplus_I.N1 - 4.0*Qi_iplus_IP1.N1 + Qi_iplus_IP2.N1, 2));

    IS_Qip[2] = (13.0 / 12.0)*(pow(Qi_iplus_IM1.N1 - 2.0*Qi_iplus_I.N1 + Qi_iplus_IP1.N1, 2)) + (1.0 / 4.0)*(pow(Qi_iplus_IM1.N1 - 4.0*Qi_iplus_I.N1 + 3.0*Qi_iplus_IP1.N1, 2));
    IS_Qip[1] = (13.0 / 12.0)*(pow(Qi_iplus_I.N1 - 2.0*Qi_iplus_IP1.N1 + Qi_iplus_IP2.N1, 2)) + (1.0 / 4.0)*(pow(Qi_iplus_I.N1 - Qi_iplus_IP2.N1, 2));
    IS_Qip[0] = (13.0 / 12.0)*(pow(Qi_iplus_IP1.N1 - 2.0*Qi_iplus_IP2.N1 + Qi_iplus_IP3.N1, 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iplus_IP1.N1 - 4.0*Qi_iplus_IP2.N1 + Qi_iplus_IP3.N1, 2));

    IS_Qinm[2] = (13.0 / 12.0)*(pow(Qi_iminus_IM3.N1 - 2.0*Qi_iminus_IM2.N1 + Qi_iminus_IM1.N1, 2)) + (1.0 / 4.0)*(pow(Qi_iminus_IM3.N1 - 4.0*Qi_iminus_IM2.N1 + 3.0*Qi_iminus_IM1.N1, 2));
    IS_Qinm[1] = (13.0 / 12.0)*(pow(Qi_iminus_IM2.N1 - 2.0*Qi_iminus_IM1.N1 + Qi_iminus_I.N1, 2)) + (1.0 / 4.0)*(pow(Qi_iminus_IM2.N1 - Qi_iminus_I.N1, 2));
    IS_Qinm[0] = (13.0 / 12.0)*(pow(Qi_iminus_IM1.N1 - 2.0*Qi_iminus_I.N1 + Qi_iminus_IP1.N1, 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iminus_IM1.N1 - 4.0*Qi_iminus_I.N1 + Qi_iminus_IP1.N1, 2));

    IS_Qinp[2] = (13.0 / 12.0)*(pow(Qi_iminus_IM2.N1 - 2.0*Qi_iminus_IM1.N1 + Qi_iminus_I.N1, 2)) + (1.0 / 4.0)*(pow(Qi_iminus_IM2.N1 - 4.0*Qi_iminus_IM1.N1 + 3.0*Qi_iminus_I.N1, 2));
    IS_Qinp[1] = (13.0 / 12.0)*(pow(Qi_iminus_IM1.N1 - 2.0*Qi_iminus_I.N1 + Qi_iminus_IP1.N1, 2)) + (1.0 / 4.0)*(pow(Qi_iminus_IM1.N1 - Qi_iminus_IP1.N1, 2));
    IS_Qinp[0] = (13.0 / 12.0)*(pow(Qi_iminus_I.N1 - 2.0*Qi_iminus_IP1.N1 + Qi_iminus_IP2.N1, 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iminus_I.N1 - 4.0*Qi_iminus_IP1.N1 + Qi_iminus_IP2.N1, 2));

    IS_Qjm[2] = (13.0 / 12.0)*(pow(Qj_iplus_JM2.N1 - 2.0*Qj_iplus_JM1.N1 + Qj_iplus_I.N1, 2)) + (1.0 / 4.0)*(pow(Qj_iplus_JM2.N1 - 4.0*Qj_iplus_JM1.N1 + 3.0*Qj_iplus_I.N1, 2));
    IS_Qjm[1] = (13.0 / 12.0)*(pow(Qj_iplus_JM1.N1 - 2.0*Qj_iplus_I.N1 + Qj_iplus_JP1.N1, 2)) + (1.0 / 4.0)*(pow(Qj_iplus_JM1.N1 - Qj_iplus_JP1.N1, 2));
    IS_Qjm[0] = (13.0 / 12.0)*(pow(Qj_iplus_I.N1 - 2.0*Qj_iplus_JP1.N1 + Qj_iplus_JP2.N1, 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iplus_I.N1 - 4.0*Qj_iplus_JP1.N1 + Qj_iplus_JP2.N1, 2));

    IS_Qjp[2] = (13.0 / 12.0)*(pow(Qj_iplus_JM1.N1 - 2.0*Qj_iplus_I.N1 + Qj_iplus_JP1.N1, 2)) + (1.0 / 4.0)*(pow(Qj_iplus_JM1.N1 - 4.0*Qj_iplus_I.N1 + 3.0*Qj_iplus_JP1.N1, 2));
    IS_Qjp[1] = (13.0 / 12.0)*(pow(Qj_iplus_I.N1 - 2.0*Qj_iplus_JP1.N1 + Qj_iplus_JP2.N1, 2)) + (1.0 / 4.0)*(pow(Qj_iplus_I.N1 - Qj_iplus_JP2.N1, 2));
    IS_Qjp[0] = (13.0 / 12.0)*(pow(Qj_iplus_JP1.N1 - 2.0*Qj_iplus_JP2.N1 + Qj_iplus_JP3.N1, 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iplus_JP1.N1 - 4.0*Qj_iplus_JP2.N1 + Qj_iplus_JP3.N1, 2));

    IS_Qjnm[2] = (13.0 / 12.0)*(pow(Qj_iminus_JM3.N1 - 2.0*Qj_iminus_JM2.N1 + Qj_iminus_JM1.N1, 2)) + (1.0 / 4.0)*(pow(Qj_iminus_JM3.N1 - 4.0*Qj_iminus_JM2.N1 + 3.0*Qj_iminus_JM1.N1, 2));
    IS_Qjnm[1] = (13.0 / 12.0)*(pow(Qj_iminus_JM2.N1 - 2.0*Qj_iminus_JM1.N1 + Qj_iminus_I.N1, 2)) + (1.0 / 4.0)*(pow(Qj_iminus_JM2.N1 - Qj_iminus_I.N1, 2));
    IS_Qjnm[0] = (13.0 / 12.0)*(pow(Qj_iminus_JM1.N1 - 2.0*Qj_iminus_I.N1 + Qj_iminus_JP1.N1, 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iminus_JM1.N1 - 4.0*Qj_iminus_I.N1 + Qj_iminus_JP1.N1, 2));

    IS_Qjnp[2] = (13.0 / 12.0)*(pow(Qj_iminus_JM2.N1 - 2.0*Qj_iminus_JM1.N1 + Qj_iminus_I.N1, 2)) + (1.0 / 4.0)*(pow(Qj_iminus_JM2.N1 - 4.0*Qj_iminus_JM1.N1 + 3.0*Qj_iminus_I.N1, 2));
    IS_Qjnp[1] = (13.0 / 12.0)*(pow(Qj_iminus_JM1.N1 - 2.0*Qj_iminus_I.N1 + Qj_iminus_JP1.N1, 2)) + (1.0 / 4.0)*(pow(Qj_iminus_JM1.N1 - Qj_iminus_JP1.N1, 2));
    IS_Qjnp[0] = (13.0 / 12.0)*(pow(Qj_iminus_I.N1 - 2.0*Qj_iminus_JP1.N1 + Qj_iminus_JP2.N1, 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iminus_I.N1 - 4.0*Qj_iminus_JP1.N1 + Qj_iminus_JP2.N1, 2));

    IS_Qkm[2] = (13.0 / 12.0)*(pow(Qk_iplus_KM2.N1 - 2.0*Qk_iplus_KM1.N1 + Qk_iplus_I.N1, 2)) + (1.0 / 4.0)*(pow(Qk_iplus_KM2.N1 - 4.0*Qk_iplus_KM1.N1 + 3.0*Qk_iplus_I.N1, 2));
    IS_Qkm[1] = (13.0 / 12.0)*(pow(Qk_iplus_KM1.N1 - 2.0*Qk_iplus_I.N1 + Qk_iplus_KP1.N1, 2)) + (1.0 / 4.0)*(pow(Qk_iplus_KM1.N1 - Qk_iplus_KP1.N1, 2));
    IS_Qkm[0] = (13.0 / 12.0)*(pow(Qk_iplus_I.N1 - 2.0*Qk_iplus_KP1.N1 + Qk_iplus_KP2.N1, 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iplus_I.N1 - 4.0*Qk_iplus_KP1.N1 + Qk_iplus_KP2.N1, 2));

    IS_Qkp[2] = (13.0 / 12.0)*(pow(Qk_iplus_KM1.N1 - 2.0*Qk_iplus_I.N1 + Qk_iplus_KP1.N1, 2)) + (1.0 / 4.0)*(pow(Qk_iplus_KM1.N1 - 4.0*Qk_iplus_I.N1 + 3.0*Qk_iplus_KP1.N1, 2));
    IS_Qkp[1] = (13.0 / 12.0)*(pow(Qk_iplus_I.N1 - 2.0*Qk_iplus_KP1.N1 + Qk_iplus_KP2.N1, 2)) + (1.0 / 4.0)*(pow(Qk_iplus_I.N1 - Qk_iplus_KP2.N1, 2));
    IS_Qkp[0] = (13.0 / 12.0)*(pow(Qk_iplus_KP1.N1 - 2.0*Qk_iplus_KP2.N1 + Qk_iplus_KP3.N1, 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iplus_KP1.N1 - 4.0*Qk_iplus_KP2.N1 + Qk_iplus_KP3.N1, 2));

    IS_Qknm[2] = (13.0 / 12.0)*(pow(Qk_iminus_KM3.N1 - 2.0*Qk_iminus_KM2.N1 + Qk_iminus_KM1.N1, 2)) + (1.0 / 4.0)*(pow(Qk_iminus_KM3.N1 - 4.0*Qk_iminus_KM2.N1 + 3.0*Qk_iminus_KM1.N1, 2));
    IS_Qknm[1] = (13.0 / 12.0)*(pow(Qk_iminus_KM2.N1 - 2.0*Qk_iminus_KM1.N1 + Qk_iminus_I.N1, 2)) + (1.0 / 4.0)*(pow(Qk_iminus_KM2.N1 - Qk_iminus_I.N1, 2));
    IS_Qknm[0] = (13.0 / 12.0)*(pow(Qk_iminus_KM1.N1 - 2.0*Qk_iminus_I.N1 + Qk_iminus_KP1.N1, 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iminus_KM1.N1 - 4.0*Qk_iminus_I.N1 + Qk_iminus_KP1.N1, 2));

    IS_Qknp[2] = (13.0 / 12.0)*(pow(Qk_iminus_KM2.N1 - 2.0*Qk_iminus_KM1.N1 + Qk_iminus_I.N1, 2)) + (1.0 / 4.0)*(pow(Qk_iminus_KM2.N1 - 4.0*Qk_iminus_KM1.N1 + 3.0*Qk_iminus_I.N1, 2));
    IS_Qknp[1] = (13.0 / 12.0)*(pow(Qk_iminus_KM1.N1 - 2.0*Qk_iminus_I.N1 + Qk_iminus_KP1.N1, 2)) + (1.0 / 4.0)*(pow(Qk_iminus_KM1.N1 - Qk_iminus_KP1.N1, 2));
    IS_Qknp[0] = (13.0 / 12.0)*(pow(Qk_iminus_I.N1 - 2.0*Qk_iminus_KP1.N1 + Qk_iminus_KP2.N1, 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iminus_I.N1 - 4.0*Qk_iminus_KP1.N1 + Qk_iminus_KP2.N1, 2));

    w_Qip[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qip[0]), 2.0));
    w_Qip[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qip[1]), 2.0));
    w_Qip[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qip[2]), 2.0));

    w_Qim[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qim[0]), 2.0));
    w_Qim[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qim[1]), 2.0));
    w_Qim[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qim[2]), 2.0));

    w_Qinp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qinp[0]), 2.0));
    w_Qinp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qinp[1]), 2.0));
    w_Qinp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qinp[2]), 2.0));

    w_Qinm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qinm[0]), 2.0));
    w_Qinm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qinm[1]), 2.0));
    w_Qinm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qinm[2]), 2.0));

    w_Qjp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjp[0]), 2.0));
    w_Qjp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjp[1]), 2.0));
    w_Qjp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjp[2]), 2.0));

    w_Qjm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjm[0]), 2.0));
    w_Qjm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjm[1]), 2.0));
    w_Qjm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjm[2]), 2.0));

    w_Qjnp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjnp[0]), 2.0));
    w_Qjnp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjnp[1]), 2.0));
    w_Qjnp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjnp[2]), 2.0));

    w_Qjnm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjnm[0]), 2.0));
    w_Qjnm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjnm[1]), 2.0));
    w_Qjnm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjnm[2]), 2.0));

    w_Qkp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qkp[0]), 2.0));
    w_Qkp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qkp[1]), 2.0));
    w_Qkp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qkp[2]), 2.0));

    w_Qkm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qkm[0]), 2.0));
    w_Qkm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qkm[1]), 2.0));
    w_Qkm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qkm[2]), 2.0));

    w_Qknp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qknp[0]), 2.0));
    w_Qknp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qknp[1]), 2.0));
    w_Qknp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qknp[2]), 2.0));

    w_Qknm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qknm[0]), 2.0));
    w_Qknm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qknm[1]), 2.0));
    w_Qknm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qknm[2]), 2.0));

    W_Qip[0] = w_Qip[0] / (w_Qip[0] + w_Qip[1] + w_Qip[2]);
    W_Qip[1] = w_Qip[1] / (w_Qip[0] + w_Qip[1] + w_Qip[2]);
    W_Qip[2] = w_Qip[2] / (w_Qip[0] + w_Qip[1] + w_Qip[2]);

    W_Qinp[0] = w_Qinp[0] / (w_Qinp[0] + w_Qinp[1] + w_Qinp[2]);
    W_Qinp[1] = w_Qinp[1] / (w_Qinp[0] + w_Qinp[1] + w_Qinp[2]);
    W_Qinp[2] = w_Qinp[2] / (w_Qinp[0] + w_Qinp[1] + w_Qinp[2]);

    W_Qim[0] = w_Qim[0] / (w_Qim[0] + w_Qim[1] + w_Qim[2]);
    W_Qim[1] = w_Qim[1] / (w_Qim[0] + w_Qim[1] + w_Qim[2]);
    W_Qim[2] = w_Qim[2] / (w_Qim[0] + w_Qim[1] + w_Qim[2]);

    W_Qinm[0] = w_Qinm[0] / (w_Qinm[0] + w_Qinm[1] + w_Qinm[2]);
    W_Qinm[1] = w_Qinm[1] / (w_Qinm[0] + w_Qinm[1] + w_Qinm[2]);
    W_Qinm[2] = w_Qinm[2] / (w_Qinm[0] + w_Qinm[1] + w_Qinm[2]);

    W_Qjp[0] = w_Qjp[0] / (w_Qjp[0] + w_Qjp[1] + w_Qjp[2]);
    W_Qjp[1] = w_Qjp[1] / (w_Qjp[0] + w_Qjp[1] + w_Qjp[2]);
    W_Qjp[2] = w_Qjp[2] / (w_Qjp[0] + w_Qjp[1] + w_Qjp[2]);

    W_Qjnp[0] = w_Qjnp[0] / (w_Qjnp[0] + w_Qjnp[1] + w_Qjnp[2]);
    W_Qjnp[1] = w_Qjnp[1] / (w_Qjnp[0] + w_Qjnp[1] + w_Qjnp[2]);
    W_Qjnp[2] = w_Qjnp[2] / (w_Qjnp[0] + w_Qjnp[1] + w_Qjnp[2]);

    W_Qjm[0] = w_Qjm[0] / (w_Qjm[0] + w_Qjm[1] + w_Qjm[2]);
    W_Qjm[1] = w_Qjm[1] / (w_Qjm[0] + w_Qjm[1] + w_Qjm[2]);
    W_Qjm[2] = w_Qjm[2] / (w_Qjm[0] + w_Qjm[1] + w_Qjm[2]);

    W_Qjnm[0] = w_Qjnm[0] / (w_Qjnm[0] + w_Qjnm[1] + w_Qjnm[2]);
    W_Qjnm[1] = w_Qjnm[1] / (w_Qjnm[0] + w_Qjnm[1] + w_Qjnm[2]);
    W_Qjnm[2] = w_Qjnm[2] / (w_Qjnm[0] + w_Qjnm[1] + w_Qjnm[2]);

    W_Qkp[0] = w_Qkp[0] / (w_Qkp[0] + w_Qkp[1] + w_Qkp[2]);
    W_Qkp[1] = w_Qkp[1] / (w_Qkp[0] + w_Qkp[1] + w_Qkp[2]);
    W_Qkp[2] = w_Qkp[2] / (w_Qkp[0] + w_Qkp[1] + w_Qkp[2]);

    W_Qknp[0] = w_Qknp[0] / (w_Qknp[0] + w_Qknp[1] + w_Qknp[2]);
    W_Qknp[1] = w_Qknp[1] / (w_Qknp[0] + w_Qknp[1] + w_Qknp[2]);
    W_Qknp[2] = w_Qknp[2] / (w_Qknp[0] + w_Qknp[1] + w_Qknp[2]);

    W_Qkm[0] = w_Qkm[0] / (w_Qkm[0] + w_Qkm[1] + w_Qkm[2]);
    W_Qkm[1] = w_Qkm[1] / (w_Qkm[0] + w_Qkm[1] + w_Qkm[2]);
    W_Qkm[2] = w_Qkm[2] / (w_Qkm[0] + w_Qkm[1] + w_Qkm[2]);

    W_Qknm[0] = w_Qknm[0] / (w_Qknm[0] + w_Qknm[1] + w_Qknm[2]);
    W_Qknm[1] = w_Qknm[1] / (w_Qknm[0] + w_Qknm[1] + w_Qknm[2]);
    W_Qknm[2] = w_Qknm[2] / (w_Qknm[0] + w_Qknm[1] + w_Qknm[2]);

    Qi_iplus_half_pos_char[1] = W_Qip[0] * Qi_half_p[0] + W_Qip[1] * Qi_half_p[1] + W_Qip[2] * Qi_half_p[2];
    Qi_iminus_half_pos_char[1] = W_Qinp[0] * Qi_half_np[0] + W_Qinp[1] * Qi_half_np[1] + W_Qinp[2] * Qi_half_np[2];

    Qi_iplus_half_neg_char[1] = W_Qim[0] * Qi_half_m[0] + W_Qim[1] * Qi_half_m[1] + W_Qim[2] * Qi_half_m[2];
    Qi_iminus_half_neg_char[1] = W_Qinm[0] * Qi_halfn_m[0] + W_Qinm[1] * Qi_halfn_m[1] + W_Qinm[2] * Qi_halfn_m[2];

    Qj_iplus_half_pos_char[1] = W_Qjp[0] * Qj_half_p[0] + W_Qjp[1] * Qj_half_p[1] + W_Qjp[2] * Qj_half_p[2];
    Qj_iminus_half_pos_char[1] = W_Qjnp[0] * Qj_half_np[0] + W_Qjnp[1] * Qj_half_np[1] + W_Qjnp[2] * Qj_half_np[2];

    Qj_iplus_half_neg_char[1] = W_Qjm[0] * Qj_half_m[0] + W_Qjm[1] * Qj_half_m[1] + W_Qjm[2] * Qj_half_m[2];
    Qj_iminus_half_neg_char[1] = W_Qjnm[0] * Qj_halfn_m[0] + W_Qjnm[1] * Qj_halfn_m[1] + W_Qjnm[2] * Qj_halfn_m[2];

    Qk_iplus_half_pos_char[1] = W_Qkp[0] * Qk_half_p[0] + W_Qkp[1] * Qk_half_p[1] + W_Qkp[2] * Qk_half_p[2];
    Qk_iminus_half_pos_char[1] = W_Qknp[0] * Qk_half_np[0] + W_Qknp[1] * Qk_half_np[1] + W_Qknp[2] * Qk_half_np[2];

    Qk_iplus_half_neg_char[1] = W_Qkm[0] * Qk_half_m[0] + W_Qkm[1] * Qk_half_m[1] + W_Qkm[2] * Qk_half_m[2];
    Qk_iminus_half_neg_char[1] = W_Qknm[0] * Qk_halfn_m[0] + W_Qknm[1] * Qk_halfn_m[1] + W_Qknm[2] * Qk_halfn_m[2];


    /******************************Equation 3****************************************/
    /** Qi(i+1/2)_plus**/
    Qi_half_p[0] = (15.0 / 8.0)*Qi_iplus_IP1.N2 - (5.0 / 4.0)*Qi_iplus_IP2.N2 + (3.0 / 8.0)*Qi_iplus_IP3.N2;
    Qi_half_p[1] = (3.0 / 8.0)*Qi_iplus_I.N2 + (3.0 / 4.0)*Qi_iplus_IP1.N2 - (1.0 / 8.0)*Qi_iplus_IP2.N2;
    Qi_half_p[2] = (-1.0 / 8.0)*Qi_iplus_IM1.N2 + (3.0 / 4.0)*Qi_iplus_I.N2 + (3.0 / 8.0)*Qi_iplus_IP1.N2;

    /** Qi(i+1/2)_minus**/
    Qi_half_m[0] = (3.0 / 8.0)*Qi_iplus_I.N2 + (3.0 / 4.0)*Qi_iplus_IP1.N2 - (1.0 / 8.0)*Qi_iplus_IP2.N2;
    Qi_half_m[1] = (-1.0 / 8.0)*Qi_iplus_IM1.N2 + (3.0 / 4.0)*Qi_iplus_I.N2 + (3.0 / 8.0)*Qi_iplus_IP1.N2;
    Qi_half_m[2] = (3.0 / 8.0)*Qi_iplus_IM2.N2 - (5.0 / 4.0)*Qi_iplus_IM1.N2 + (15.0 / 8.0)*Qi_iplus_I.N2;

    /** Qi(i-1/2)_minus**/
    Qi_halfn_m[0] = (3.0 / 8.0)*Qi_iminus_IM1.N2 + (3.0 / 4.0)*Qi_iminus_I.N2 - (1.0 / 8.0)*Qi_iminus_IP1.N2;
    Qi_halfn_m[1] = (-1.0 / 8.0)*Qi_iminus_IM2.N2 + (3.0 / 4.0)*Qi_iminus_IM1.N2 + (3.0 / 8.0)*Qi_iminus_I.N2;
    Qi_halfn_m[2] = (3.0 / 8.0)*Qi_iminus_IM3.N2 - (5.0 / 4.0)*Qi_iminus_IM2.N2 + (15.0 / 8.0)*Qi_iminus_IM1.N2;

    /** Qi(i-1/2)_plus**/
    Qi_half_np[0] = (15.0 / 8.0)*Qi_iminus_I.N2 - (5.0 / 4.0)*Qi_iminus_IP1.N2 + (3.0 / 8.0)*Qi_iminus_IP2.N2;
    Qi_half_np[1] = (3.0 / 8.0)*Qi_iminus_IM1.N2 + (3.0 / 4.0)*Qi_iminus_I.N2 - (1.0 / 8.0)*Qi_iminus_IP1.N2;
    Qi_half_np[2] = (-1.0 / 8.0)*Qi_iminus_IM2.N2 + (3.0 / 4.0)*Qi_iminus_IM1.N2 + (3.0 / 8.0)*Qi_iminus_I.N2;

    /** Qj(i+1/2)_plus**/
    Qj_half_p[0] = (15.0 / 8.0)*Qj_iplus_JP1.N2 - (5.0 / 4.0)*Qj_iplus_JP2.N2 + (3.0 / 8.0)*Qj_iplus_JP3.N2;
    Qj_half_p[1] = (3.0 / 8.0)*Qj_iplus_I.N2 + (3.0 / 4.0)*Qj_iplus_JP1.N2 - (1.0 / 8.0)*Qj_iplus_JP2.N2;
    Qj_half_p[2] = (-1.0 / 8.0)*Qj_iplus_JM1.N2 + (3.0 / 4.0)*Qj_iplus_I.N2 + (3.0 / 8.0)*Qj_iplus_JP1.N2;

    /** Qj(i+1/2)_minus**/
    Qj_half_m[0] = (3.0 / 8.0)*Qj_iplus_I.N2 + (3.0 / 4.0)*Qj_iplus_JP1.N2 - (1.0 / 8.0)*Qj_iplus_JP2.N2;
    Qj_half_m[1] = (-1.0 / 8.0)*Qj_iplus_JM1.N2 + (3.0 / 4.0)*Qj_iplus_I.N2 + (3.0 / 8.0)*Qj_iplus_JP1.N2;
    Qj_half_m[2] = (3.0 / 8.0)*Qj_iplus_JM2.N2 - (5.0 / 4.0)*Qj_iplus_JM1.N2 + (15.0 / 8.0)*Qj_iplus_I.N2;

    /** Qj(i-1/2)_minus**/
    Qj_halfn_m[0] = (3.0 / 8.0)*Qj_iminus_JM1.N2 + (3.0 / 4.0)*Qj_iminus_I.N2 - (1.0 / 8.0)*Qj_iminus_JP1.N2;
    Qj_halfn_m[1] = (-1.0 / 8.0)*Qj_iminus_JM2.N2 + (3.0 / 4.0)*Qj_iminus_JM1.N2 + (3.0 / 8.0)*Qj_iminus_I.N2;
    Qj_halfn_m[2] = (3.0 / 8.0)*Qj_iminus_JM3.N2 - (5.0 / 4.0)*Qj_iminus_JM2.N2 + (15.0 / 8.0)*Qj_iminus_JM1.N2;

    /** Qj(i-1/2)_plus**/
    Qj_half_np[0] = (15.0 / 8.0)*Qj_iminus_I.N2 - (5.0 / 4.0)*Qj_iminus_JP1.N2 + (3.0 / 8.0)*Qj_iminus_JP2.N2;
    Qj_half_np[1] = (3.0 / 8.0)*Qj_iminus_JM1.N2 + (3.0 / 4.0)*Qj_iminus_I.N2 - (1.0 / 8.0)*Qj_iminus_JP1.N2;
    Qj_half_np[2] = (-1.0 / 8.0)*Qj_iminus_JM2.N2 + (3.0 / 4.0)*Qj_iminus_JM1.N2 + (3.0 / 8.0)*Qj_iminus_I.N2;

    /** Qk(i+1/2)_plus**/
    Qk_half_p[0] = (15.0 / 8.0)*Qk_iplus_KP1.N2 - (5.0 / 4.0)*Qk_iplus_KP2.N2 + (3.0 / 8.0)*Qk_iplus_KP3.N2;
    Qk_half_p[1] = (3.0 / 8.0)*Qk_iplus_I.N2 + (3.0 / 4.0)*Qk_iplus_KP1.N2 - (1.0 / 8.0)*Qk_iplus_KP2.N2;
    Qk_half_p[2] = (-1.0 / 8.0)*Qk_iplus_KM1.N2 + (3.0 / 4.0)*Qk_iplus_I.N2 + (3.0 / 8.0)*Qk_iplus_KP1.N2;

    /** Qk(i+1/2)_minus**/
    Qk_half_m[0] = (3.0 / 8.0)*Qk_iplus_I.N2 + (3.0 / 4.0)*Qk_iplus_KP1.N2 - (1.0 / 8.0)*Qk_iplus_KP2.N2;
    Qk_half_m[1] = (-1.0 / 8.0)*Qk_iplus_KM1.N2 + (3.0 / 4.0)*Qk_iplus_I.N2 + (3.0 / 8.0)*Qk_iplus_KP1.N2;
    Qk_half_m[2] = (3.0 / 8.0)*Qk_iplus_KM2.N2 - (5.0 / 4.0)*Qk_iplus_KM1.N2 + (15.0 / 8.0)*Qk_iplus_I.N2;

    /** Qk(i-1/2)_minus**/
    Qk_halfn_m[0] = (3.0 / 8.0)*Qk_iminus_KM1.N2 + (3.0 / 4.0)*Qk_iminus_I.N2 - (1.0 / 8.0)*Qk_iminus_KP1.N2;
    Qk_halfn_m[1] = (-1.0 / 8.0)*Qk_iminus_KM2.N2 + (3.0 / 4.0)*Qk_iminus_KM1.N2 + (3.0 / 8.0)*Qk_iminus_I.N2;
    Qk_halfn_m[2] = (3.0 / 8.0)*Qk_iminus_KM3.N2 - (5.0 / 4.0)*Qk_iminus_KM2.N2 + (15.0 / 8.0)*Qk_iminus_KM1.N2;

    /** Qk(i-1/2)_plus**/
    Qk_half_np[0] = (15.0 / 8.0)*Qk_iminus_I.N2 - (5.0 / 4.0)*Qk_iminus_KP1.N2 + (3.0 / 8.0)*Qk_iminus_KP2.N2;
    Qk_half_np[1] = (3.0 / 8.0)*Qk_iminus_KM1.N2 + (3.0 / 4.0)*Qk_iminus_I.N2 - (1.0 / 8.0)*Qk_iminus_KP1.N2;
    Qk_half_np[2] = (-1.0 / 8.0)*Qk_iminus_KM2.N2 + (3.0 / 4.0)*Qk_iminus_KM1.N2 + (3.0 / 8.0)*Qk_iminus_I.N2;

    /**********************************SMOOTHNESS INDICATOR OR SMOOTHNESS FUNCTION********************************/

    IS_Qim[2] = (13.0 / 12.0)*(pow(Qi_iplus_IM2.N2 - 2.0*Qi_iplus_IM1.N2 + Qi_iplus_I.N2, 2)) + (1.0 / 4.0)*(pow(Qi_iplus_IM2.N2 - 4.0*Qi_iplus_IM1.N2 + 3.0*Qi_iplus_I.N2, 2));
    IS_Qim[1] = (13.0 / 12.0)*(pow(Qi_iplus_IM1.N2 - 2.0*Qi_iplus_I.N2 + Qi_iplus_IP1.N2, 2)) + (1.0 / 4.0)*(pow(Qi_iplus_IM1.N2 - Qi_iplus_IP1.N2, 2));
    IS_Qim[0] = (13.0 / 12.0)*(pow(Qi_iplus_I.N2 - 2.0*Qi_iplus_IP1.N2 + Qi_iplus_IP2.N2, 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iplus_I.N2 - 4.0*Qi_iplus_IP1.N2 + Qi_iplus_IP2.N2, 2));

    IS_Qip[2] = (13.0 / 12.0)*(pow(Qi_iplus_IM1.N2 - 2.0*Qi_iplus_I.N2 + Qi_iplus_IP1.N2, 2)) + (1.0 / 4.0)*(pow(Qi_iplus_IM1.N2 - 4.0*Qi_iplus_I.N2 + 3.0*Qi_iplus_IP1.N2, 2));
    IS_Qip[1] = (13.0 / 12.0)*(pow(Qi_iplus_I.N2 - 2.0*Qi_iplus_IP1.N2 + Qi_iplus_IP2.N2, 2)) + (1.0 / 4.0)*(pow(Qi_iplus_I.N2 - Qi_iplus_IP2.N2, 2));
    IS_Qip[0] = (13.0 / 12.0)*(pow(Qi_iplus_IP1.N2 - 2.0*Qi_iplus_IP2.N2 + Qi_iplus_IP3.N2, 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iplus_IP1.N2 - 4.0*Qi_iplus_IP2.N2 + Qi_iplus_IP3.N2, 2));

    IS_Qinm[2] = (13.0 / 12.0)*(pow(Qi_iminus_IM3.N2 - 2.0*Qi_iminus_IM2.N2 + Qi_iminus_IM1.N2, 2)) + (1.0 / 4.0)*(pow(Qi_iminus_IM3.N2 - 4.0*Qi_iminus_IM2.N2 + 3.0*Qi_iminus_IM1.N2, 2));
    IS_Qinm[1] = (13.0 / 12.0)*(pow(Qi_iminus_IM2.N2 - 2.0*Qi_iminus_IM1.N2 + Qi_iminus_I.N2, 2)) + (1.0 / 4.0)*(pow(Qi_iminus_IM2.N2 - Qi_iminus_I.N2, 2));
    IS_Qinm[0] = (13.0 / 12.0)*(pow(Qi_iminus_IM1.N2 - 2.0*Qi_iminus_I.N2 + Qi_iminus_IP1.N2, 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iminus_IM1.N2 - 4.0*Qi_iminus_I.N2 + Qi_iminus_IP1.N2, 2));

    IS_Qinp[2] = (13.0 / 12.0)*(pow(Qi_iminus_IM2.N2 - 2.0*Qi_iminus_IM1.N2 + Qi_iminus_I.N2, 2)) + (1.0 / 4.0)*(pow(Qi_iminus_IM2.N2 - 4.0*Qi_iminus_IM1.N2 + 3.0*Qi_iminus_I.N2, 2));
    IS_Qinp[1] = (13.0 / 12.0)*(pow(Qi_iminus_IM1.N2 - 2.0*Qi_iminus_I.N2 + Qi_iminus_IP1.N2, 2)) + (1.0 / 4.0)*(pow(Qi_iminus_IM1.N2 - Qi_iminus_IP1.N2, 2));
    IS_Qinp[0] = (13.0 / 12.0)*(pow(Qi_iminus_I.N2 - 2.0*Qi_iminus_IP1.N2 + Qi_iminus_IP2.N2, 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iminus_I.N2 - 4.0*Qi_iminus_IP1.N2 + Qi_iminus_IP2.N2, 2));

    IS_Qjm[2] = (13.0 / 12.0)*(pow(Qj_iplus_JM2.N2 - 2.0*Qj_iplus_JM1.N2 + Qj_iplus_I.N2, 2)) + (1.0 / 4.0)*(pow(Qj_iplus_JM2.N2 - 4.0*Qj_iplus_JM1.N2 + 3.0*Qj_iplus_I.N2, 2));
    IS_Qjm[1] = (13.0 / 12.0)*(pow(Qj_iplus_JM1.N2 - 2.0*Qj_iplus_I.N2 + Qj_iplus_JP1.N2, 2)) + (1.0 / 4.0)*(pow(Qj_iplus_JM1.N2 - Qj_iplus_JP1.N2, 2));
    IS_Qjm[0] = (13.0 / 12.0)*(pow(Qj_iplus_I.N2 - 2.0*Qj_iplus_JP1.N2 + Qj_iplus_JP2.N2, 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iplus_I.N2 - 4.0*Qj_iplus_JP1.N2 + Qj_iplus_JP2.N2, 2));

    IS_Qjp[2] = (13.0 / 12.0)*(pow(Qj_iplus_JM1.N2 - 2.0*Qj_iplus_I.N2 + Qj_iplus_JP1.N2, 2)) + (1.0 / 4.0)*(pow(Qj_iplus_JM1.N2 - 4.0*Qj_iplus_I.N2 + 3.0*Qj_iplus_JP1.N2, 2));
    IS_Qjp[1] = (13.0 / 12.0)*(pow(Qj_iplus_I.N2 - 2.0*Qj_iplus_JP1.N2 + Qj_iplus_JP2.N2, 2)) + (1.0 / 4.0)*(pow(Qj_iplus_I.N2 - Qj_iplus_JP2.N2, 2));
    IS_Qjp[0] = (13.0 / 12.0)*(pow(Qj_iplus_JP1.N2 - 2.0*Qj_iplus_JP2.N2 + Qj_iplus_JP3.N2, 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iplus_JP1.N2 - 4.0*Qj_iplus_JP2.N2 + Qj_iplus_JP3.N2, 2));

    IS_Qjnm[2] = (13.0 / 12.0)*(pow(Qj_iminus_JM3.N2 - 2.0*Qj_iminus_JM2.N2 + Qj_iminus_JM1.N2, 2)) + (1.0 / 4.0)*(pow(Qj_iminus_JM3.N2 - 4.0*Qj_iminus_JM2.N2 + 3.0*Qj_iminus_JM1.N2, 2));
    IS_Qjnm[1] = (13.0 / 12.0)*(pow(Qj_iminus_JM2.N2 - 2.0*Qj_iminus_JM1.N2 + Qj_iminus_I.N2, 2)) + (1.0 / 4.0)*(pow(Qj_iminus_JM2.N2 - Qj_iminus_I.N2, 2));
    IS_Qjnm[0] = (13.0 / 12.0)*(pow(Qj_iminus_JM1.N2 - 2.0*Qj_iminus_I.N2 + Qj_iminus_JP1.N2, 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iminus_JM1.N2 - 4.0*Qj_iminus_I.N2 + Qj_iminus_JP1.N2, 2));

    IS_Qjnp[2] = (13.0 / 12.0)*(pow(Qj_iminus_JM2.N2 - 2.0*Qj_iminus_JM1.N2 + Qj_iminus_I.N2, 2)) + (1.0 / 4.0)*(pow(Qj_iminus_JM2.N2 - 4.0*Qj_iminus_JM1.N2 + 3.0*Qj_iminus_I.N2, 2));
    IS_Qjnp[1] = (13.0 / 12.0)*(pow(Qj_iminus_JM1.N2 - 2.0*Qj_iminus_I.N2 + Qj_iminus_JP1.N2, 2)) + (1.0 / 4.0)*(pow(Qj_iminus_JM1.N2 - Qj_iminus_JP1.N2, 2));
    IS_Qjnp[0] = (13.0 / 12.0)*(pow(Qj_iminus_I.N2 - 2.0*Qj_iminus_JP1.N2 + Qj_iminus_JP2.N2, 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iminus_I.N2 - 4.0*Qj_iminus_JP1.N2 + Qj_iminus_JP2.N2, 2));

    IS_Qkm[2] = (13.0 / 12.0)*(pow(Qk_iplus_KM2.N2 - 2.0*Qk_iplus_KM1.N2 + Qk_iplus_I.N2, 2)) + (1.0 / 4.0)*(pow(Qk_iplus_KM2.N2 - 4.0*Qk_iplus_KM1.N2 + 3.0*Qk_iplus_I.N2, 2));
    IS_Qkm[1] = (13.0 / 12.0)*(pow(Qk_iplus_KM1.N2 - 2.0*Qk_iplus_I.N2 + Qk_iplus_KP1.N2, 2)) + (1.0 / 4.0)*(pow(Qk_iplus_KM1.N2 - Qk_iplus_KP1.N2, 2));
    IS_Qkm[0] = (13.0 / 12.0)*(pow(Qk_iplus_I.N2 - 2.0*Qk_iplus_KP1.N2 + Qk_iplus_KP2.N2, 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iplus_I.N2 - 4.0*Qk_iplus_KP1.N2 + Qk_iplus_KP2.N2, 2));

    IS_Qkp[2] = (13.0 / 12.0)*(pow(Qk_iplus_KM1.N2 - 2.0*Qk_iplus_I.N2 + Qk_iplus_KP1.N2, 2)) + (1.0 / 4.0)*(pow(Qk_iplus_KM1.N2 - 4.0*Qk_iplus_I.N2 + 3.0*Qk_iplus_KP1.N2, 2));
    IS_Qkp[1] = (13.0 / 12.0)*(pow(Qk_iplus_I.N2 - 2.0*Qk_iplus_KP1.N2 + Qk_iplus_KP2.N2, 2)) + (1.0 / 4.0)*(pow(Qk_iplus_I.N2 - Qk_iplus_KP2.N2, 2));
    IS_Qkp[0] = (13.0 / 12.0)*(pow(Qk_iplus_KP1.N2 - 2.0*Qk_iplus_KP2.N2 + Qk_iplus_KP3.N2, 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iplus_KP1.N2 - 4.0*Qk_iplus_KP2.N2 + Qk_iplus_KP3.N2, 2));

    IS_Qknm[2] = (13.0 / 12.0)*(pow(Qk_iminus_KM3.N2 - 2.0*Qk_iminus_KM2.N2 + Qk_iminus_KM1.N2, 2)) + (1.0 / 4.0)*(pow(Qk_iminus_KM3.N2 - 4.0*Qk_iminus_KM2.N2 + 3.0*Qk_iminus_KM1.N2, 2));
    IS_Qknm[1] = (13.0 / 12.0)*(pow(Qk_iminus_KM2.N2 - 2.0*Qk_iminus_KM1.N2 + Qk_iminus_I.N2, 2)) + (1.0 / 4.0)*(pow(Qk_iminus_KM2.N2 - Qk_iminus_I.N2, 2));
    IS_Qknm[0] = (13.0 / 12.0)*(pow(Qk_iminus_KM1.N2 - 2.0*Qk_iminus_I.N2 + Qk_iminus_KP1.N2, 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iminus_KM1.N2 - 4.0*Qk_iminus_I.N2 + Qk_iminus_KP1.N2, 2));

    IS_Qknp[2] = (13.0 / 12.0)*(pow(Qk_iminus_KM2.N2 - 2.0*Qk_iminus_KM1.N2 + Qk_iminus_I.N2, 2)) + (1.0 / 4.0)*(pow(Qk_iminus_KM2.N2 - 4.0*Qk_iminus_KM1.N2 + 3.0*Qk_iminus_I.N2, 2));
    IS_Qknp[1] = (13.0 / 12.0)*(pow(Qk_iminus_KM1.N2 - 2.0*Qk_iminus_I.N2 + Qk_iminus_KP1.N2, 2)) + (1.0 / 4.0)*(pow(Qk_iminus_KM1.N2 - Qk_iminus_KP1.N2, 2));
    IS_Qknp[0] = (13.0 / 12.0)*(pow(Qk_iminus_I.N2 - 2.0*Qk_iminus_KP1.N2 + Qk_iminus_KP2.N2, 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iminus_I.N2 - 4.0*Qk_iminus_KP1.N2 + Qk_iminus_KP2.N2, 2));

    w_Qip[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qip[0]), 2.0));
    w_Qip[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qip[1]), 2.0));
    w_Qip[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qip[2]), 2.0));

    w_Qim[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qim[0]), 2.0));
    w_Qim[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qim[1]), 2.0));
    w_Qim[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qim[2]), 2.0));

    w_Qinp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qinp[0]), 2.0));
    w_Qinp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qinp[1]), 2.0));
    w_Qinp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qinp[2]), 2.0));

    w_Qinm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qinm[0]), 2.0));
    w_Qinm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qinm[1]), 2.0));
    w_Qinm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qinm[2]), 2.0));

    w_Qjp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjp[0]), 2.0));
    w_Qjp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjp[1]), 2.0));
    w_Qjp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjp[2]), 2.0));

    w_Qjm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjm[0]), 2.0));
    w_Qjm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjm[1]), 2.0));
    w_Qjm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjm[2]), 2.0));

    w_Qjnp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjnp[0]), 2.0));
    w_Qjnp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjnp[1]), 2.0));
    w_Qjnp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjnp[2]), 2.0));

    w_Qjnm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjnm[0]), 2.0));
    w_Qjnm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjnm[1]), 2.0));
    w_Qjnm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjnm[2]), 2.0));

    w_Qkp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qkp[0]), 2.0));
    w_Qkp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qkp[1]), 2.0));
    w_Qkp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qkp[2]), 2.0));

    w_Qkm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qkm[0]), 2.0));
    w_Qkm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qkm[1]), 2.0));
    w_Qkm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qkm[2]), 2.0));

    w_Qknp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qknp[0]), 2.0));
    w_Qknp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qknp[1]), 2.0));
    w_Qknp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qknp[2]), 2.0));

    w_Qknm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qknm[0]), 2.0));
    w_Qknm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qknm[1]), 2.0));
    w_Qknm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qknm[2]), 2.0));

    W_Qip[0] = w_Qip[0] / (w_Qip[0] + w_Qip[1] + w_Qip[2]);
    W_Qip[1] = w_Qip[1] / (w_Qip[0] + w_Qip[1] + w_Qip[2]);
    W_Qip[2] = w_Qip[2] / (w_Qip[0] + w_Qip[1] + w_Qip[2]);

    W_Qinp[0] = w_Qinp[0] / (w_Qinp[0] + w_Qinp[1] + w_Qinp[2]);
    W_Qinp[1] = w_Qinp[1] / (w_Qinp[0] + w_Qinp[1] + w_Qinp[2]);
    W_Qinp[2] = w_Qinp[2] / (w_Qinp[0] + w_Qinp[1] + w_Qinp[2]);

    W_Qim[0] = w_Qim[0] / (w_Qim[0] + w_Qim[1] + w_Qim[2]);
    W_Qim[1] = w_Qim[1] / (w_Qim[0] + w_Qim[1] + w_Qim[2]);
    W_Qim[2] = w_Qim[2] / (w_Qim[0] + w_Qim[1] + w_Qim[2]);

    W_Qinm[0] = w_Qinm[0] / (w_Qinm[0] + w_Qinm[1] + w_Qinm[2]);
    W_Qinm[1] = w_Qinm[1] / (w_Qinm[0] + w_Qinm[1] + w_Qinm[2]);
    W_Qinm[2] = w_Qinm[2] / (w_Qinm[0] + w_Qinm[1] + w_Qinm[2]);

    W_Qjp[0] = w_Qjp[0] / (w_Qjp[0] + w_Qjp[1] + w_Qjp[2]);
    W_Qjp[1] = w_Qjp[1] / (w_Qjp[0] + w_Qjp[1] + w_Qjp[2]);
    W_Qjp[2] = w_Qjp[2] / (w_Qjp[0] + w_Qjp[1] + w_Qjp[2]);

    W_Qjnp[0] = w_Qjnp[0] / (w_Qjnp[0] + w_Qjnp[1] + w_Qjnp[2]);
    W_Qjnp[1] = w_Qjnp[1] / (w_Qjnp[0] + w_Qjnp[1] + w_Qjnp[2]);
    W_Qjnp[2] = w_Qjnp[2] / (w_Qjnp[0] + w_Qjnp[1] + w_Qjnp[2]);

    W_Qjm[0] = w_Qjm[0] / (w_Qjm[0] + w_Qjm[1] + w_Qjm[2]);
    W_Qjm[1] = w_Qjm[1] / (w_Qjm[0] + w_Qjm[1] + w_Qjm[2]);
    W_Qjm[2] = w_Qjm[2] / (w_Qjm[0] + w_Qjm[1] + w_Qjm[2]);

    W_Qjnm[0] = w_Qjnm[0] / (w_Qjnm[0] + w_Qjnm[1] + w_Qjnm[2]);
    W_Qjnm[1] = w_Qjnm[1] / (w_Qjnm[0] + w_Qjnm[1] + w_Qjnm[2]);
    W_Qjnm[2] = w_Qjnm[2] / (w_Qjnm[0] + w_Qjnm[1] + w_Qjnm[2]);

    W_Qkp[0] = w_Qkp[0] / (w_Qkp[0] + w_Qkp[1] + w_Qkp[2]);
    W_Qkp[1] = w_Qkp[1] / (w_Qkp[0] + w_Qkp[1] + w_Qkp[2]);
    W_Qkp[2] = w_Qkp[2] / (w_Qkp[0] + w_Qkp[1] + w_Qkp[2]);

    W_Qknp[0] = w_Qknp[0] / (w_Qknp[0] + w_Qknp[1] + w_Qknp[2]);
    W_Qknp[1] = w_Qknp[1] / (w_Qknp[0] + w_Qknp[1] + w_Qknp[2]);
    W_Qknp[2] = w_Qknp[2] / (w_Qknp[0] + w_Qknp[1] + w_Qknp[2]);

    W_Qkm[0] = w_Qkm[0] / (w_Qkm[0] + w_Qkm[1] + w_Qkm[2]);
    W_Qkm[1] = w_Qkm[1] / (w_Qkm[0] + w_Qkm[1] + w_Qkm[2]);
    W_Qkm[2] = w_Qkm[2] / (w_Qkm[0] + w_Qkm[1] + w_Qkm[2]);

    W_Qknm[0] = w_Qknm[0] / (w_Qknm[0] + w_Qknm[1] + w_Qknm[2]);
    W_Qknm[1] = w_Qknm[1] / (w_Qknm[0] + w_Qknm[1] + w_Qknm[2]);
    W_Qknm[2] = w_Qknm[2] / (w_Qknm[0] + w_Qknm[1] + w_Qknm[2]);

    Qi_iplus_half_pos_char[2] = W_Qip[0] * Qi_half_p[0] + W_Qip[1] * Qi_half_p[1] + W_Qip[2] * Qi_half_p[2];
    Qi_iminus_half_pos_char[2] = W_Qinp[0] * Qi_half_np[0] + W_Qinp[1] * Qi_half_np[1] + W_Qinp[2] * Qi_half_np[2];

    Qi_iplus_half_neg_char[2] = W_Qim[0] * Qi_half_m[0] + W_Qim[1] * Qi_half_m[1] + W_Qim[2] * Qi_half_m[2];
    Qi_iminus_half_neg_char[2] = W_Qinm[0] * Qi_halfn_m[0] + W_Qinm[1] * Qi_halfn_m[1] + W_Qinm[2] * Qi_halfn_m[2];

    Qj_iplus_half_pos_char[2] = W_Qjp[0] * Qj_half_p[0] + W_Qjp[1] * Qj_half_p[1] + W_Qjp[2] * Qj_half_p[2];
    Qj_iminus_half_pos_char[2] = W_Qjnp[0] * Qj_half_np[0] + W_Qjnp[1] * Qj_half_np[1] + W_Qjnp[2] * Qj_half_np[2];

    Qj_iplus_half_neg_char[2] = W_Qjm[0] * Qj_half_m[0] + W_Qjm[1] * Qj_half_m[1] + W_Qjm[2] * Qj_half_m[2];
    Qj_iminus_half_neg_char[2] = W_Qjnm[0] * Qj_halfn_m[0] + W_Qjnm[1] * Qj_halfn_m[1] + W_Qjnm[2] * Qj_halfn_m[2];

    Qk_iplus_half_pos_char[2] = W_Qkp[0] * Qk_half_p[0] + W_Qkp[1] * Qk_half_p[1] + W_Qkp[2] * Qk_half_p[2];
    Qk_iminus_half_pos_char[2] = W_Qknp[0] * Qk_half_np[0] + W_Qknp[1] * Qk_half_np[1] + W_Qknp[2] * Qk_half_np[2];

    Qk_iplus_half_neg_char[2] = W_Qkm[0] * Qk_half_m[0] + W_Qkm[1] * Qk_half_m[1] + W_Qkm[2] * Qk_half_m[2];
    Qk_iminus_half_neg_char[2] = W_Qknm[0] * Qk_halfn_m[0] + W_Qknm[1] * Qk_halfn_m[1] + W_Qknm[2] * Qk_halfn_m[2];


    /******************************Equation 4****************************************/
    /** Qi(i+1/2)_plus**/
    Qi_half_p[0] = (15.0 / 8.0)*Qi_iplus_IP1.N3 - (5.0 / 4.0)*Qi_iplus_IP2.N3 + (3.0 / 8.0)*Qi_iplus_IP3.N3;
    Qi_half_p[1] = (3.0 / 8.0)*Qi_iplus_I.N3 + (3.0 / 4.0)*Qi_iplus_IP1.N3 - (1.0 / 8.0)*Qi_iplus_IP2.N3;
    Qi_half_p[2] = (-1.0 / 8.0)*Qi_iplus_IM1.N3 + (3.0 / 4.0)*Qi_iplus_I.N3 + (3.0 / 8.0)*Qi_iplus_IP1.N3;

    /** Qi(i+1/2)_minus**/
    Qi_half_m[0] = (3.0 / 8.0)*Qi_iplus_I.N3 + (3.0 / 4.0)*Qi_iplus_IP1.N3 - (1.0 / 8.0)*Qi_iplus_IP2.N3;
    Qi_half_m[1] = (-1.0 / 8.0)*Qi_iplus_IM1.N3 + (3.0 / 4.0)*Qi_iplus_I.N3 + (3.0 / 8.0)*Qi_iplus_IP1.N3;
    Qi_half_m[2] = (3.0 / 8.0)*Qi_iplus_IM2.N3 - (5.0 / 4.0)*Qi_iplus_IM1.N3 + (15.0 / 8.0)*Qi_iplus_I.N3;

    /** Qi(i-1/2)_minus**/
    Qi_halfn_m[0] = (3.0 / 8.0)*Qi_iminus_IM1.N3 + (3.0 / 4.0)*Qi_iminus_I.N3 - (1.0 / 8.0)*Qi_iminus_IP1.N3;
    Qi_halfn_m[1] = (-1.0 / 8.0)*Qi_iminus_IM2.N3 + (3.0 / 4.0)*Qi_iminus_IM1.N3 + (3.0 / 8.0)*Qi_iminus_I.N3;
    Qi_halfn_m[2] = (3.0 / 8.0)*Qi_iminus_IM3.N3 - (5.0 / 4.0)*Qi_iminus_IM2.N3 + (15.0 / 8.0)*Qi_iminus_IM1.N3;

    /** Qi(i-1/2)_plus**/
    Qi_half_np[0] = (15.0 / 8.0)*Qi_iminus_I.N3 - (5.0 / 4.0)*Qi_iminus_IP1.N3 + (3.0 / 8.0)*Qi_iminus_IP2.N3;
    Qi_half_np[1] = (3.0 / 8.0)*Qi_iminus_IM1.N3 + (3.0 / 4.0)*Qi_iminus_I.N3 - (1.0 / 8.0)*Qi_iminus_IP1.N3;
    Qi_half_np[2] = (-1.0 / 8.0)*Qi_iminus_IM2.N3 + (3.0 / 4.0)*Qi_iminus_IM1.N3 + (3.0 / 8.0)*Qi_iminus_I.N3;

    /** Qj(i+1/2)_plus**/
    Qj_half_p[0] = (15.0 / 8.0)*Qj_iplus_JP1.N3 - (5.0 / 4.0)*Qj_iplus_JP2.N3 + (3.0 / 8.0)*Qj_iplus_JP3.N3;
    Qj_half_p[1] = (3.0 / 8.0)*Qj_iplus_I.N3 + (3.0 / 4.0)*Qj_iplus_JP1.N3 - (1.0 / 8.0)*Qj_iplus_JP2.N3;
    Qj_half_p[2] = (-1.0 / 8.0)*Qj_iplus_JM1.N3 + (3.0 / 4.0)*Qj_iplus_I.N3 + (3.0 / 8.0)*Qj_iplus_JP1.N3;

    /** Qj(i+1/2)_minus**/
    Qj_half_m[0] = (3.0 / 8.0)*Qj_iplus_I.N3 + (3.0 / 4.0)*Qj_iplus_JP1.N3 - (1.0 / 8.0)*Qj_iplus_JP2.N3;
    Qj_half_m[1] = (-1.0 / 8.0)*Qj_iplus_JM1.N3 + (3.0 / 4.0)*Qj_iplus_I.N3 + (3.0 / 8.0)*Qj_iplus_JP1.N3;
    Qj_half_m[2] = (3.0 / 8.0)*Qj_iplus_JM2.N3 - (5.0 / 4.0)*Qj_iplus_JM1.N3 + (15.0 / 8.0)*Qj_iplus_I.N3;

    /** Qj(i-1/2)_minus**/
    Qj_halfn_m[0] = (3.0 / 8.0)*Qj_iminus_JM1.N3 + (3.0 / 4.0)*Qj_iminus_I.N3 - (1.0 / 8.0)*Qj_iminus_JP1.N3;
    Qj_halfn_m[1] = (-1.0 / 8.0)*Qj_iminus_JM2.N3 + (3.0 / 4.0)*Qj_iminus_JM1.N3 + (3.0 / 8.0)*Qj_iminus_I.N3;
    Qj_halfn_m[2] = (3.0 / 8.0)*Qj_iminus_JM3.N3 - (5.0 / 4.0)*Qj_iminus_JM2.N3 + (15.0 / 8.0)*Qj_iminus_JM1.N3;

    /** Qj(i-1/2)_plus**/
    Qj_half_np[0] = (15.0 / 8.0)*Qj_iminus_I.N3 - (5.0 / 4.0)*Qj_iminus_JP1.N3 + (3.0 / 8.0)*Qj_iminus_JP2.N3;
    Qj_half_np[1] = (3.0 / 8.0)*Qj_iminus_JM1.N3 + (3.0 / 4.0)*Qj_iminus_I.N3 - (1.0 / 8.0)*Qj_iminus_JP1.N3;
    Qj_half_np[2] = (-1.0 / 8.0)*Qj_iminus_JM2.N3 + (3.0 / 4.0)*Qj_iminus_JM1.N3 + (3.0 / 8.0)*Qj_iminus_I.N3;

    /** Qk(i+1/2)_plus**/
    Qk_half_p[0] = (15.0 / 8.0)*Qk_iplus_KP1.N3 - (5.0 / 4.0)*Qk_iplus_KP2.N3 + (3.0 / 8.0)*Qk_iplus_KP3.N3;
    Qk_half_p[1] = (3.0 / 8.0)*Qk_iplus_I.N3 + (3.0 / 4.0)*Qk_iplus_KP1.N3 - (1.0 / 8.0)*Qk_iplus_KP2.N3;
    Qk_half_p[2] = (-1.0 / 8.0)*Qk_iplus_KM1.N3 + (3.0 / 4.0)*Qk_iplus_I.N3 + (3.0 / 8.0)*Qk_iplus_KP1.N3;

    /** Qk(i+1/2)_minus**/
    Qk_half_m[0] = (3.0 / 8.0)*Qk_iplus_I.N3 + (3.0 / 4.0)*Qk_iplus_KP1.N3 - (1.0 / 8.0)*Qk_iplus_KP2.N3;
    Qk_half_m[1] = (-1.0 / 8.0)*Qk_iplus_KM1.N3 + (3.0 / 4.0)*Qk_iplus_I.N3 + (3.0 / 8.0)*Qk_iplus_KP1.N3;
    Qk_half_m[2] = (3.0 / 8.0)*Qk_iplus_KM2.N3 - (5.0 / 4.0)*Qk_iplus_KM1.N3 + (15.0 / 8.0)*Qk_iplus_I.N3;

    /** Qk(i-1/2)_minus**/
    Qk_halfn_m[0] = (3.0 / 8.0)*Qk_iminus_KM1.N3 + (3.0 / 4.0)*Qk_iminus_I.N3 - (1.0 / 8.0)*Qk_iminus_KP1.N3;
    Qk_halfn_m[1] = (-1.0 / 8.0)*Qk_iminus_KM2.N3 + (3.0 / 4.0)*Qk_iminus_KM1.N3 + (3.0 / 8.0)*Qk_iminus_I.N3;
    Qk_halfn_m[2] = (3.0 / 8.0)*Qk_iminus_KM3.N3 - (5.0 / 4.0)*Qk_iminus_KM2.N3 + (15.0 / 8.0)*Qk_iminus_KM1.N3;

    /** Qk(i-1/2)_plus**/
    Qk_half_np[0] = (15.0 / 8.0)*Qk_iminus_I.N3 - (5.0 / 4.0)*Qk_iminus_KP1.N3 + (3.0 / 8.0)*Qk_iminus_KP2.N3;
    Qk_half_np[1] = (3.0 / 8.0)*Qk_iminus_KM1.N3 + (3.0 / 4.0)*Qk_iminus_I.N3 - (1.0 / 8.0)*Qk_iminus_KP1.N3;
    Qk_half_np[2] = (-1.0 / 8.0)*Qk_iminus_KM2.N3 + (3.0 / 4.0)*Qk_iminus_KM1.N3 + (3.0 / 8.0)*Qk_iminus_I.N3;

    /**********************************SMOOTHNESS INDICATOR OR SMOOTHNESS FUNCTION********************************/

    IS_Qim[2] = (13.0 / 12.0)*(pow(Qi_iplus_IM2.N3 - 2.0*Qi_iplus_IM1.N3 + Qi_iplus_I.N3, 2)) + (1.0 / 4.0)*(pow(Qi_iplus_IM2.N3 - 4.0*Qi_iplus_IM1.N3 + 3.0*Qi_iplus_I.N3, 2));
    IS_Qim[1] = (13.0 / 12.0)*(pow(Qi_iplus_IM1.N3 - 2.0*Qi_iplus_I.N3 + Qi_iplus_IP1.N3, 2)) + (1.0 / 4.0)*(pow(Qi_iplus_IM1.N3 - Qi_iplus_IP1.N3, 2));
    IS_Qim[0] = (13.0 / 12.0)*(pow(Qi_iplus_I.N3 - 2.0*Qi_iplus_IP1.N3 + Qi_iplus_IP2.N3, 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iplus_I.N3 - 4.0*Qi_iplus_IP1.N3 + Qi_iplus_IP2.N3, 2));

    IS_Qip[2] = (13.0 / 12.0)*(pow(Qi_iplus_IM1.N3 - 2.0*Qi_iplus_I.N3 + Qi_iplus_IP1.N3, 2)) + (1.0 / 4.0)*(pow(Qi_iplus_IM1.N3 - 4.0*Qi_iplus_I.N3 + 3.0*Qi_iplus_IP1.N3, 2));
    IS_Qip[1] = (13.0 / 12.0)*(pow(Qi_iplus_I.N3 - 2.0*Qi_iplus_IP1.N3 + Qi_iplus_IP2.N3, 2)) + (1.0 / 4.0)*(pow(Qi_iplus_I.N3 - Qi_iplus_IP2.N3, 2));
    IS_Qip[0] = (13.0 / 12.0)*(pow(Qi_iplus_IP1.N3 - 2.0*Qi_iplus_IP2.N3 + Qi_iplus_IP3.N3, 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iplus_IP1.N3 - 4.0*Qi_iplus_IP2.N3 + Qi_iplus_IP3.N3, 2));

    IS_Qinm[2] = (13.0 / 12.0)*(pow(Qi_iminus_IM3.N3 - 2.0*Qi_iminus_IM2.N3 + Qi_iminus_IM1.N3, 2)) + (1.0 / 4.0)*(pow(Qi_iminus_IM3.N3 - 4.0*Qi_iminus_IM2.N3 + 3.0*Qi_iminus_IM1.N3, 2));
    IS_Qinm[1] = (13.0 / 12.0)*(pow(Qi_iminus_IM2.N3 - 2.0*Qi_iminus_IM1.N3 + Qi_iminus_I.N3, 2)) + (1.0 / 4.0)*(pow(Qi_iminus_IM2.N3 - Qi_iminus_I.N3, 2));
    IS_Qinm[0] = (13.0 / 12.0)*(pow(Qi_iminus_IM1.N3 - 2.0*Qi_iminus_I.N3 + Qi_iminus_IP1.N3, 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iminus_IM1.N3 - 4.0*Qi_iminus_I.N3 + Qi_iminus_IP1.N3, 2));

    IS_Qinp[2] = (13.0 / 12.0)*(pow(Qi_iminus_IM2.N3 - 2.0*Qi_iminus_IM1.N3 + Qi_iminus_I.N3, 2)) + (1.0 / 4.0)*(pow(Qi_iminus_IM2.N3 - 4.0*Qi_iminus_IM1.N3 + 3.0*Qi_iminus_I.N3, 2));
    IS_Qinp[1] = (13.0 / 12.0)*(pow(Qi_iminus_IM1.N3 - 2.0*Qi_iminus_I.N3 + Qi_iminus_IP1.N3, 2)) + (1.0 / 4.0)*(pow(Qi_iminus_IM1.N3 - Qi_iminus_IP1.N3, 2));
    IS_Qinp[0] = (13.0 / 12.0)*(pow(Qi_iminus_I.N3 - 2.0*Qi_iminus_IP1.N3 + Qi_iminus_IP2.N3, 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iminus_I.N3 - 4.0*Qi_iminus_IP1.N3 + Qi_iminus_IP2.N3, 2));

    IS_Qjm[2] = (13.0 / 12.0)*(pow(Qj_iplus_JM2.N3 - 2.0*Qj_iplus_JM1.N3 + Qj_iplus_I.N3, 2)) + (1.0 / 4.0)*(pow(Qj_iplus_JM2.N3 - 4.0*Qj_iplus_JM1.N3 + 3.0*Qj_iplus_I.N3, 2));
    IS_Qjm[1] = (13.0 / 12.0)*(pow(Qj_iplus_JM1.N3 - 2.0*Qj_iplus_I.N3 + Qj_iplus_JP1.N3, 2)) + (1.0 / 4.0)*(pow(Qj_iplus_JM1.N3 - Qj_iplus_JP1.N3, 2));
    IS_Qjm[0] = (13.0 / 12.0)*(pow(Qj_iplus_I.N3 - 2.0*Qj_iplus_JP1.N3 + Qj_iplus_JP2.N3, 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iplus_I.N3 - 4.0*Qj_iplus_JP1.N3 + Qj_iplus_JP2.N3, 2));

    IS_Qjp[2] = (13.0 / 12.0)*(pow(Qj_iplus_JM1.N3 - 2.0*Qj_iplus_I.N3 + Qj_iplus_JP1.N3, 2)) + (1.0 / 4.0)*(pow(Qj_iplus_JM1.N3 - 4.0*Qj_iplus_I.N3 + 3.0*Qj_iplus_JP1.N3, 2));
    IS_Qjp[1] = (13.0 / 12.0)*(pow(Qj_iplus_I.N3 - 2.0*Qj_iplus_JP1.N3 + Qj_iplus_JP2.N3, 2)) + (1.0 / 4.0)*(pow(Qj_iplus_I.N3 - Qj_iplus_JP2.N3, 2));
    IS_Qjp[0] = (13.0 / 12.0)*(pow(Qj_iplus_JP1.N3 - 2.0*Qj_iplus_JP2.N3 + Qj_iplus_JP3.N3, 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iplus_JP1.N3 - 4.0*Qj_iplus_JP2.N3 + Qj_iplus_JP3.N3, 2));

    IS_Qjnm[2] = (13.0 / 12.0)*(pow(Qj_iminus_JM3.N3 - 2.0*Qj_iminus_JM2.N3 + Qj_iminus_JM1.N3, 2)) + (1.0 / 4.0)*(pow(Qj_iminus_JM3.N3 - 4.0*Qj_iminus_JM2.N3 + 3.0*Qj_iminus_JM1.N3, 2));
    IS_Qjnm[1] = (13.0 / 12.0)*(pow(Qj_iminus_JM2.N3 - 2.0*Qj_iminus_JM1.N3 + Qj_iminus_I.N3, 2)) + (1.0 / 4.0)*(pow(Qj_iminus_JM2.N3 - Qj_iminus_I.N3, 2));
    IS_Qjnm[0] = (13.0 / 12.0)*(pow(Qj_iminus_JM1.N3 - 2.0*Qj_iminus_I.N3 + Qj_iminus_JP1.N3, 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iminus_JM1.N3 - 4.0*Qj_iminus_I.N3 + Qj_iminus_JP1.N3, 2));

    IS_Qjnp[2] = (13.0 / 12.0)*(pow(Qj_iminus_JM2.N3 - 2.0*Qj_iminus_JM1.N3 + Qj_iminus_I.N3, 2)) + (1.0 / 4.0)*(pow(Qj_iminus_JM2.N3 - 4.0*Qj_iminus_JM1.N3 + 3.0*Qj_iminus_I.N3, 2));
    IS_Qjnp[1] = (13.0 / 12.0)*(pow(Qj_iminus_JM1.N3 - 2.0*Qj_iminus_I.N3 + Qj_iminus_JP1.N3, 2)) + (1.0 / 4.0)*(pow(Qj_iminus_JM1.N3 - Qj_iminus_JP1.N3, 2));
    IS_Qjnp[0] = (13.0 / 12.0)*(pow(Qj_iminus_I.N3 - 2.0*Qj_iminus_JP1.N3 + Qj_iminus_JP2.N3, 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iminus_I.N3 - 4.0*Qj_iminus_JP1.N3 + Qj_iminus_JP2.N3, 2));

    IS_Qkm[2] = (13.0 / 12.0)*(pow(Qk_iplus_KM2.N3 - 2.0*Qk_iplus_KM1.N3 + Qk_iplus_I.N3, 2)) + (1.0 / 4.0)*(pow(Qk_iplus_KM2.N3 - 4.0*Qk_iplus_KM1.N3 + 3.0*Qk_iplus_I.N3, 2));
    IS_Qkm[1] = (13.0 / 12.0)*(pow(Qk_iplus_KM1.N3 - 2.0*Qk_iplus_I.N3 + Qk_iplus_KP1.N3, 2)) + (1.0 / 4.0)*(pow(Qk_iplus_KM1.N3 - Qk_iplus_KP1.N3, 2));
    IS_Qkm[0] = (13.0 / 12.0)*(pow(Qk_iplus_I.N3 - 2.0*Qk_iplus_KP1.N3 + Qk_iplus_KP2.N3, 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iplus_I.N3 - 4.0*Qk_iplus_KP1.N3 + Qk_iplus_KP2.N3, 2));

    IS_Qkp[2] = (13.0 / 12.0)*(pow(Qk_iplus_KM1.N3 - 2.0*Qk_iplus_I.N3 + Qk_iplus_KP1.N3, 2)) + (1.0 / 4.0)*(pow(Qk_iplus_KM1.N3 - 4.0*Qk_iplus_I.N3 + 3.0*Qk_iplus_KP1.N3, 2));
    IS_Qkp[1] = (13.0 / 12.0)*(pow(Qk_iplus_I.N3 - 2.0*Qk_iplus_KP1.N3 + Qk_iplus_KP2.N3, 2)) + (1.0 / 4.0)*(pow(Qk_iplus_I.N3 - Qk_iplus_KP2.N3, 2));
    IS_Qkp[0] = (13.0 / 12.0)*(pow(Qk_iplus_KP1.N3 - 2.0*Qk_iplus_KP2.N3 + Qk_iplus_KP3.N3, 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iplus_KP1.N3 - 4.0*Qk_iplus_KP2.N3 + Qk_iplus_KP3.N3, 2));

    IS_Qknm[2] = (13.0 / 12.0)*(pow(Qk_iminus_KM3.N3 - 2.0*Qk_iminus_KM2.N3 + Qk_iminus_KM1.N3, 2)) + (1.0 / 4.0)*(pow(Qk_iminus_KM3.N3 - 4.0*Qk_iminus_KM2.N3 + 3.0*Qk_iminus_KM1.N3, 2));
    IS_Qknm[1] = (13.0 / 12.0)*(pow(Qk_iminus_KM2.N3 - 2.0*Qk_iminus_KM1.N3 + Qk_iminus_I.N3, 2)) + (1.0 / 4.0)*(pow(Qk_iminus_KM2.N3 - Qk_iminus_I.N3, 2));
    IS_Qknm[0] = (13.0 / 12.0)*(pow(Qk_iminus_KM1.N3 - 2.0*Qk_iminus_I.N3 + Qk_iminus_KP1.N3, 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iminus_KM1.N3 - 4.0*Qk_iminus_I.N3 + Qk_iminus_KP1.N3, 2));

    IS_Qknp[2] = (13.0 / 12.0)*(pow(Qk_iminus_KM2.N3 - 2.0*Qk_iminus_KM1.N3 + Qk_iminus_I.N3, 2)) + (1.0 / 4.0)*(pow(Qk_iminus_KM2.N3 - 4.0*Qk_iminus_KM1.N3 + 3.0*Qk_iminus_I.N3, 2));
    IS_Qknp[1] = (13.0 / 12.0)*(pow(Qk_iminus_KM1.N3 - 2.0*Qk_iminus_I.N3 + Qk_iminus_KP1.N3, 2)) + (1.0 / 4.0)*(pow(Qk_iminus_KM1.N3 - Qk_iminus_KP1.N3, 2));
    IS_Qknp[0] = (13.0 / 12.0)*(pow(Qk_iminus_I.N3 - 2.0*Qk_iminus_KP1.N3 + Qk_iminus_KP2.N3, 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iminus_I.N3 - 4.0*Qk_iminus_KP1.N3 + Qk_iminus_KP2.N3, 2));

    w_Qip[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qip[0]), 2.0));
    w_Qip[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qip[1]), 2.0));
    w_Qip[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qip[2]), 2.0));

    w_Qim[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qim[0]), 2.0));
    w_Qim[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qim[1]), 2.0));
    w_Qim[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qim[2]), 2.0));

    w_Qinp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qinp[0]), 2.0));
    w_Qinp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qinp[1]), 2.0));
    w_Qinp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qinp[2]), 2.0));

    w_Qinm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qinm[0]), 2.0));
    w_Qinm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qinm[1]), 2.0));
    w_Qinm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qinm[2]), 2.0));

    w_Qjp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjp[0]), 2.0));
    w_Qjp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjp[1]), 2.0));
    w_Qjp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjp[2]), 2.0));

    w_Qjm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjm[0]), 2.0));
    w_Qjm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjm[1]), 2.0));
    w_Qjm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjm[2]), 2.0));

    w_Qjnp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjnp[0]), 2.0));
    w_Qjnp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjnp[1]), 2.0));
    w_Qjnp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjnp[2]), 2.0));

    w_Qjnm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjnm[0]), 2.0));
    w_Qjnm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjnm[1]), 2.0));
    w_Qjnm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjnm[2]), 2.0));

    w_Qkp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qkp[0]), 2.0));
    w_Qkp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qkp[1]), 2.0));
    w_Qkp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qkp[2]), 2.0));

    w_Qkm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qkm[0]), 2.0));
    w_Qkm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qkm[1]), 2.0));
    w_Qkm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qkm[2]), 2.0));

    w_Qknp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qknp[0]), 2.0));
    w_Qknp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qknp[1]), 2.0));
    w_Qknp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qknp[2]), 2.0));

    w_Qknm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qknm[0]), 2.0));
    w_Qknm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qknm[1]), 2.0));
    w_Qknm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qknm[2]), 2.0));

    W_Qip[0] = w_Qip[0] / (w_Qip[0] + w_Qip[1] + w_Qip[2]);
    W_Qip[1] = w_Qip[1] / (w_Qip[0] + w_Qip[1] + w_Qip[2]);
    W_Qip[2] = w_Qip[2] / (w_Qip[0] + w_Qip[1] + w_Qip[2]);

    W_Qinp[0] = w_Qinp[0] / (w_Qinp[0] + w_Qinp[1] + w_Qinp[2]);
    W_Qinp[1] = w_Qinp[1] / (w_Qinp[0] + w_Qinp[1] + w_Qinp[2]);
    W_Qinp[2] = w_Qinp[2] / (w_Qinp[0] + w_Qinp[1] + w_Qinp[2]);

    W_Qim[0] = w_Qim[0] / (w_Qim[0] + w_Qim[1] + w_Qim[2]);
    W_Qim[1] = w_Qim[1] / (w_Qim[0] + w_Qim[1] + w_Qim[2]);
    W_Qim[2] = w_Qim[2] / (w_Qim[0] + w_Qim[1] + w_Qim[2]);

    W_Qinm[0] = w_Qinm[0] / (w_Qinm[0] + w_Qinm[1] + w_Qinm[2]);
    W_Qinm[1] = w_Qinm[1] / (w_Qinm[0] + w_Qinm[1] + w_Qinm[2]);
    W_Qinm[2] = w_Qinm[2] / (w_Qinm[0] + w_Qinm[1] + w_Qinm[2]);

    W_Qjp[0] = w_Qjp[0] / (w_Qjp[0] + w_Qjp[1] + w_Qjp[2]);
    W_Qjp[1] = w_Qjp[1] / (w_Qjp[0] + w_Qjp[1] + w_Qjp[2]);
    W_Qjp[2] = w_Qjp[2] / (w_Qjp[0] + w_Qjp[1] + w_Qjp[2]);

    W_Qjnp[0] = w_Qjnp[0] / (w_Qjnp[0] + w_Qjnp[1] + w_Qjnp[2]);
    W_Qjnp[1] = w_Qjnp[1] / (w_Qjnp[0] + w_Qjnp[1] + w_Qjnp[2]);
    W_Qjnp[2] = w_Qjnp[2] / (w_Qjnp[0] + w_Qjnp[1] + w_Qjnp[2]);

    W_Qjm[0] = w_Qjm[0] / (w_Qjm[0] + w_Qjm[1] + w_Qjm[2]);
    W_Qjm[1] = w_Qjm[1] / (w_Qjm[0] + w_Qjm[1] + w_Qjm[2]);
    W_Qjm[2] = w_Qjm[2] / (w_Qjm[0] + w_Qjm[1] + w_Qjm[2]);

    W_Qjnm[0] = w_Qjnm[0] / (w_Qjnm[0] + w_Qjnm[1] + w_Qjnm[2]);
    W_Qjnm[1] = w_Qjnm[1] / (w_Qjnm[0] + w_Qjnm[1] + w_Qjnm[2]);
    W_Qjnm[2] = w_Qjnm[2] / (w_Qjnm[0] + w_Qjnm[1] + w_Qjnm[2]);

    W_Qkp[0] = w_Qkp[0] / (w_Qkp[0] + w_Qkp[1] + w_Qkp[2]);
    W_Qkp[1] = w_Qkp[1] / (w_Qkp[0] + w_Qkp[1] + w_Qkp[2]);
    W_Qkp[2] = w_Qkp[2] / (w_Qkp[0] + w_Qkp[1] + w_Qkp[2]);

    W_Qknp[0] = w_Qknp[0] / (w_Qknp[0] + w_Qknp[1] + w_Qknp[2]);
    W_Qknp[1] = w_Qknp[1] / (w_Qknp[0] + w_Qknp[1] + w_Qknp[2]);
    W_Qknp[2] = w_Qknp[2] / (w_Qknp[0] + w_Qknp[1] + w_Qknp[2]);

    W_Qkm[0] = w_Qkm[0] / (w_Qkm[0] + w_Qkm[1] + w_Qkm[2]);
    W_Qkm[1] = w_Qkm[1] / (w_Qkm[0] + w_Qkm[1] + w_Qkm[2]);
    W_Qkm[2] = w_Qkm[2] / (w_Qkm[0] + w_Qkm[1] + w_Qkm[2]);

    W_Qknm[0] = w_Qknm[0] / (w_Qknm[0] + w_Qknm[1] + w_Qknm[2]);
    W_Qknm[1] = w_Qknm[1] / (w_Qknm[0] + w_Qknm[1] + w_Qknm[2]);
    W_Qknm[2] = w_Qknm[2] / (w_Qknm[0] + w_Qknm[1] + w_Qknm[2]);

    Qi_iplus_half_pos_char[3] = W_Qip[0] * Qi_half_p[0] + W_Qip[1] * Qi_half_p[1] + W_Qip[2] * Qi_half_p[2];
    Qi_iminus_half_pos_char[3] = W_Qinp[0] * Qi_half_np[0] + W_Qinp[1] * Qi_half_np[1] + W_Qinp[2] * Qi_half_np[2];

    Qi_iplus_half_neg_char[3] = W_Qim[0] * Qi_half_m[0] + W_Qim[1] * Qi_half_m[1] + W_Qim[2] * Qi_half_m[2];
    Qi_iminus_half_neg_char[3] = W_Qinm[0] * Qi_halfn_m[0] + W_Qinm[1] * Qi_halfn_m[1] + W_Qinm[2] * Qi_halfn_m[2];

    Qj_iplus_half_pos_char[3] = W_Qjp[0] * Qj_half_p[0] + W_Qjp[1] * Qj_half_p[1] + W_Qjp[2] * Qj_half_p[2];
    Qj_iminus_half_pos_char[3] = W_Qjnp[0] * Qj_half_np[0] + W_Qjnp[1] * Qj_half_np[1] + W_Qjnp[2] * Qj_half_np[2];

    Qj_iplus_half_neg_char[3] = W_Qjm[0] * Qj_half_m[0] + W_Qjm[1] * Qj_half_m[1] + W_Qjm[2] * Qj_half_m[2];
    Qj_iminus_half_neg_char[3] = W_Qjnm[0] * Qj_halfn_m[0] + W_Qjnm[1] * Qj_halfn_m[1] + W_Qjnm[2] * Qj_halfn_m[2];

    Qk_iplus_half_pos_char[3] = W_Qkp[0] * Qk_half_p[0] + W_Qkp[1] * Qk_half_p[1] + W_Qkp[2] * Qk_half_p[2];
    Qk_iminus_half_pos_char[3] = W_Qknp[0] * Qk_half_np[0] + W_Qknp[1] * Qk_half_np[1] + W_Qknp[2] * Qk_half_np[2];

    Qk_iplus_half_neg_char[3] = W_Qkm[0] * Qk_half_m[0] + W_Qkm[1] * Qk_half_m[1] + W_Qkm[2] * Qk_half_m[2];
    Qk_iminus_half_neg_char[3] = W_Qknm[0] * Qk_halfn_m[0] + W_Qknm[1] * Qk_halfn_m[1] + W_Qknm[2] * Qk_halfn_m[2];




    /******************************Equation 5****************************************/
    /** Qi(i+1/2)_plus**/
    Qi_half_p[0] = (15.0 / 8.0)*Qi_iplus_IP1.N4 - (5.0 / 4.0)*Qi_iplus_IP2.N4 + (3.0 / 8.0)*Qi_iplus_IP3.N4;
    Qi_half_p[1] = (3.0 / 8.0)*Qi_iplus_I.N4 + (3.0 / 4.0)*Qi_iplus_IP1.N4 - (1.0 / 8.0)*Qi_iplus_IP2.N4;
    Qi_half_p[2] = (-1.0 / 8.0)*Qi_iplus_IM1.N4 + (3.0 / 4.0)*Qi_iplus_I.N4 + (3.0 / 8.0)*Qi_iplus_IP1.N4;

    /** Qi(i+1/2)_minus**/
    Qi_half_m[0] = (3.0 / 8.0)*Qi_iplus_I.N4 + (3.0 / 4.0)*Qi_iplus_IP1.N4 - (1.0 / 8.0)*Qi_iplus_IP2.N4;
    Qi_half_m[1] = (-1.0 / 8.0)*Qi_iplus_IM1.N4 + (3.0 / 4.0)*Qi_iplus_I.N4 + (3.0 / 8.0)*Qi_iplus_IP1.N4;
    Qi_half_m[2] = (3.0 / 8.0)*Qi_iplus_IM2.N4 - (5.0 / 4.0)*Qi_iplus_IM1.N4 + (15.0 / 8.0)*Qi_iplus_I.N4;

    /** Qi(i-1/2)_minus**/
    Qi_halfn_m[0] = (3.0 / 8.0)*Qi_iminus_IM1.N4 + (3.0 / 4.0)*Qi_iminus_I.N4 - (1.0 / 8.0)*Qi_iminus_IP1.N4;
    Qi_halfn_m[1] = (-1.0 / 8.0)*Qi_iminus_IM2.N4 + (3.0 / 4.0)*Qi_iminus_IM1.N4 + (3.0 / 8.0)*Qi_iminus_I.N4;
    Qi_halfn_m[2] = (3.0 / 8.0)*Qi_iminus_IM3.N4 - (5.0 / 4.0)*Qi_iminus_IM2.N4 + (15.0 / 8.0)*Qi_iminus_IM1.N4;

    /** Qi(i-1/2)_plus**/
    Qi_half_np[0] = (15.0 / 8.0)*Qi_iminus_I.N4 - (5.0 / 4.0)*Qi_iminus_IP1.N4 + (3.0 / 8.0)*Qi_iminus_IP2.N4;
    Qi_half_np[1] = (3.0 / 8.0)*Qi_iminus_IM1.N4 + (3.0 / 4.0)*Qi_iminus_I.N4 - (1.0 / 8.0)*Qi_iminus_IP1.N4;
    Qi_half_np[2] = (-1.0 / 8.0)*Qi_iminus_IM2.N4 + (3.0 / 4.0)*Qi_iminus_IM1.N4 + (3.0 / 8.0)*Qi_iminus_I.N4;

    /** Qj(i+1/2)_plus**/
    Qj_half_p[0] = (15.0 / 8.0)*Qj_iplus_JP1.N4 - (5.0 / 4.0)*Qj_iplus_JP2.N4 + (3.0 / 8.0)*Qj_iplus_JP3.N4;
    Qj_half_p[1] = (3.0 / 8.0)*Qj_iplus_I.N4 + (3.0 / 4.0)*Qj_iplus_JP1.N4 - (1.0 / 8.0)*Qj_iplus_JP2.N4;
    Qj_half_p[2] = (-1.0 / 8.0)*Qj_iplus_JM1.N4 + (3.0 / 4.0)*Qj_iplus_I.N4 + (3.0 / 8.0)*Qj_iplus_JP1.N4;

    /** Qj(i+1/2)_minus**/
    Qj_half_m[0] = (3.0 / 8.0)*Qj_iplus_I.N4 + (3.0 / 4.0)*Qj_iplus_JP1.N4 - (1.0 / 8.0)*Qj_iplus_JP2.N4;
    Qj_half_m[1] = (-1.0 / 8.0)*Qj_iplus_JM1.N4 + (3.0 / 4.0)*Qj_iplus_I.N4 + (3.0 / 8.0)*Qj_iplus_JP1.N4;
    Qj_half_m[2] = (3.0 / 8.0)*Qj_iplus_JM2.N4 - (5.0 / 4.0)*Qj_iplus_JM1.N4 + (15.0 / 8.0)*Qj_iplus_I.N4;

    /** Qj(i-1/2)_minus**/
    Qj_halfn_m[0] = (3.0 / 8.0)*Qj_iminus_JM1.N4 + (3.0 / 4.0)*Qj_iminus_I.N4 - (1.0 / 8.0)*Qj_iminus_JP1.N4;
    Qj_halfn_m[1] = (-1.0 / 8.0)*Qj_iminus_JM2.N4 + (3.0 / 4.0)*Qj_iminus_JM1.N4 + (3.0 / 8.0)*Qj_iminus_I.N4;
    Qj_halfn_m[2] = (3.0 / 8.0)*Qj_iminus_JM3.N4 - (5.0 / 4.0)*Qj_iminus_JM2.N4 + (15.0 / 8.0)*Qj_iminus_JM1.N4;

    /** Qj(i-1/2)_plus**/
    Qj_half_np[0] = (15.0 / 8.0)*Qj_iminus_I.N4 - (5.0 / 4.0)*Qj_iminus_JP1.N4 + (3.0 / 8.0)*Qj_iminus_JP2.N4;
    Qj_half_np[1] = (3.0 / 8.0)*Qj_iminus_JM1.N4 + (3.0 / 4.0)*Qj_iminus_I.N4 - (1.0 / 8.0)*Qj_iminus_JP1.N4;
    Qj_half_np[2] = (-1.0 / 8.0)*Qj_iminus_JM2.N4 + (3.0 / 4.0)*Qj_iminus_JM1.N4 + (3.0 / 8.0)*Qj_iminus_I.N4;

    /** Qk(i+1/2)_plus**/
    Qk_half_p[0] = (15.0 / 8.0)*Qk_iplus_KP1.N4 - (5.0 / 4.0)*Qk_iplus_KP2.N4 + (3.0 / 8.0)*Qk_iplus_KP3.N4;
    Qk_half_p[1] = (3.0 / 8.0)*Qk_iplus_I.N4 + (3.0 / 4.0)*Qk_iplus_KP1.N4 - (1.0 / 8.0)*Qk_iplus_KP2.N4;
    Qk_half_p[2] = (-1.0 / 8.0)*Qk_iplus_KM1.N4 + (3.0 / 4.0)*Qk_iplus_I.N4 + (3.0 / 8.0)*Qk_iplus_KP1.N4;

    /** Qk(i+1/2)_minus**/
    Qk_half_m[0] = (3.0 / 8.0)*Qk_iplus_I.N4 + (3.0 / 4.0)*Qk_iplus_KP1.N4 - (1.0 / 8.0)*Qk_iplus_KP2.N4;
    Qk_half_m[1] = (-1.0 / 8.0)*Qk_iplus_KM1.N4 + (3.0 / 4.0)*Qk_iplus_I.N4 + (3.0 / 8.0)*Qk_iplus_KP1.N4;
    Qk_half_m[2] = (3.0 / 8.0)*Qk_iplus_KM2.N4 - (5.0 / 4.0)*Qk_iplus_KM1.N4 + (15.0 / 8.0)*Qk_iplus_I.N4;

    /** Qk(i-1/2)_minus**/
    Qk_halfn_m[0] = (3.0 / 8.0)*Qk_iminus_KM1.N4 + (3.0 / 4.0)*Qk_iminus_I.N4 - (1.0 / 8.0)*Qk_iminus_KP1.N4;
    Qk_halfn_m[1] = (-1.0 / 8.0)*Qk_iminus_KM2.N4 + (3.0 / 4.0)*Qk_iminus_KM1.N4 + (3.0 / 8.0)*Qk_iminus_I.N4;
    Qk_halfn_m[2] = (3.0 / 8.0)*Qk_iminus_KM3.N4 - (5.0 / 4.0)*Qk_iminus_KM2.N4 + (15.0 / 8.0)*Qk_iminus_KM1.N4;

    /** Qk(i-1/2)_plus**/
    Qk_half_np[0] = (15.0 / 8.0)*Qk_iminus_I.N4 - (5.0 / 4.0)*Qk_iminus_KP1.N4 + (3.0 / 8.0)*Qk_iminus_KP2.N4;
    Qk_half_np[1] = (3.0 / 8.0)*Qk_iminus_KM1.N4 + (3.0 / 4.0)*Qk_iminus_I.N4 - (1.0 / 8.0)*Qk_iminus_KP1.N4;
    Qk_half_np[2] = (-1.0 / 8.0)*Qk_iminus_KM2.N4 + (3.0 / 4.0)*Qk_iminus_KM1.N4 + (3.0 / 8.0)*Qk_iminus_I.N4;

    /**********************************SMOOTHNESS INDICATOR OR SMOOTHNESS FUNCTION********************************/

    IS_Qim[2] = (13.0 / 12.0)*(pow(Qi_iplus_IM2.N4 - 2.0*Qi_iplus_IM1.N4 + Qi_iplus_I.N4, 2)) + (1.0 / 4.0)*(pow(Qi_iplus_IM2.N4 - 4.0*Qi_iplus_IM1.N4 + 3.0*Qi_iplus_I.N4, 2));
    IS_Qim[1] = (13.0 / 12.0)*(pow(Qi_iplus_IM1.N4 - 2.0*Qi_iplus_I.N4 + Qi_iplus_IP1.N4, 2)) + (1.0 / 4.0)*(pow(Qi_iplus_IM1.N4 - Qi_iplus_IP1.N4, 2));
    IS_Qim[0] = (13.0 / 12.0)*(pow(Qi_iplus_I.N4 - 2.0*Qi_iplus_IP1.N4 + Qi_iplus_IP2.N4, 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iplus_I.N4 - 4.0*Qi_iplus_IP1.N4 + Qi_iplus_IP2.N4, 2));

    IS_Qip[2] = (13.0 / 12.0)*(pow(Qi_iplus_IM1.N4 - 2.0*Qi_iplus_I.N4 + Qi_iplus_IP1.N4, 2)) + (1.0 / 4.0)*(pow(Qi_iplus_IM1.N4 - 4.0*Qi_iplus_I.N4 + 3.0*Qi_iplus_IP1.N4, 2));
    IS_Qip[1] = (13.0 / 12.0)*(pow(Qi_iplus_I.N4 - 2.0*Qi_iplus_IP1.N4 + Qi_iplus_IP2.N4, 2)) + (1.0 / 4.0)*(pow(Qi_iplus_I.N4 - Qi_iplus_IP2.N4, 2));
    IS_Qip[0] = (13.0 / 12.0)*(pow(Qi_iplus_IP1.N4 - 2.0*Qi_iplus_IP2.N4 + Qi_iplus_IP3.N4, 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iplus_IP1.N4 - 4.0*Qi_iplus_IP2.N4 + Qi_iplus_IP3.N4, 2));

    IS_Qinm[2] = (13.0 / 12.0)*(pow(Qi_iminus_IM3.N4 - 2.0*Qi_iminus_IM2.N4 + Qi_iminus_IM1.N4, 2)) + (1.0 / 4.0)*(pow(Qi_iminus_IM3.N4 - 4.0*Qi_iminus_IM2.N4 + 3.0*Qi_iminus_IM1.N4, 2));
    IS_Qinm[1] = (13.0 / 12.0)*(pow(Qi_iminus_IM2.N4 - 2.0*Qi_iminus_IM1.N4 + Qi_iminus_I.N4, 2)) + (1.0 / 4.0)*(pow(Qi_iminus_IM2.N4 - Qi_iminus_I.N4, 2));
    IS_Qinm[0] = (13.0 / 12.0)*(pow(Qi_iminus_IM1.N4 - 2.0*Qi_iminus_I.N4 + Qi_iminus_IP1.N4, 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iminus_IM1.N4 - 4.0*Qi_iminus_I.N4 + Qi_iminus_IP1.N4, 2));

    IS_Qinp[2] = (13.0 / 12.0)*(pow(Qi_iminus_IM2.N4 - 2.0*Qi_iminus_IM1.N4 + Qi_iminus_I.N4, 2)) + (1.0 / 4.0)*(pow(Qi_iminus_IM2.N4 - 4.0*Qi_iminus_IM1.N4 + 3.0*Qi_iminus_I.N4, 2));
    IS_Qinp[1] = (13.0 / 12.0)*(pow(Qi_iminus_IM1.N4 - 2.0*Qi_iminus_I.N4 + Qi_iminus_IP1.N4, 2)) + (1.0 / 4.0)*(pow(Qi_iminus_IM1.N4 - Qi_iminus_IP1.N4, 2));
    IS_Qinp[0] = (13.0 / 12.0)*(pow(Qi_iminus_I.N4 - 2.0*Qi_iminus_IP1.N4 + Qi_iminus_IP2.N4, 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iminus_I.N4 - 4.0*Qi_iminus_IP1.N4 + Qi_iminus_IP2.N4, 2));

    IS_Qjm[2] = (13.0 / 12.0)*(pow(Qj_iplus_JM2.N4 - 2.0*Qj_iplus_JM1.N4 + Qj_iplus_I.N4, 2)) + (1.0 / 4.0)*(pow(Qj_iplus_JM2.N4 - 4.0*Qj_iplus_JM1.N4 + 3.0*Qj_iplus_I.N4, 2));
    IS_Qjm[1] = (13.0 / 12.0)*(pow(Qj_iplus_JM1.N4 - 2.0*Qj_iplus_I.N4 + Qj_iplus_JP1.N4, 2)) + (1.0 / 4.0)*(pow(Qj_iplus_JM1.N4 - Qj_iplus_JP1.N4, 2));
    IS_Qjm[0] = (13.0 / 12.0)*(pow(Qj_iplus_I.N4 - 2.0*Qj_iplus_JP1.N4 + Qj_iplus_JP2.N4, 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iplus_I.N4 - 4.0*Qj_iplus_JP1.N4 + Qj_iplus_JP2.N4, 2));

    IS_Qjp[2] = (13.0 / 12.0)*(pow(Qj_iplus_JM1.N4 - 2.0*Qj_iplus_I.N4 + Qj_iplus_JP1.N4, 2)) + (1.0 / 4.0)*(pow(Qj_iplus_JM1.N4 - 4.0*Qj_iplus_I.N4 + 3.0*Qj_iplus_JP1.N4, 2));
    IS_Qjp[1] = (13.0 / 12.0)*(pow(Qj_iplus_I.N4 - 2.0*Qj_iplus_JP1.N4 + Qj_iplus_JP2.N4, 2)) + (1.0 / 4.0)*(pow(Qj_iplus_I.N4 - Qj_iplus_JP2.N4, 2));
    IS_Qjp[0] = (13.0 / 12.0)*(pow(Qj_iplus_JP1.N4 - 2.0*Qj_iplus_JP2.N4 + Qj_iplus_JP3.N4, 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iplus_JP1.N4 - 4.0*Qj_iplus_JP2.N4 + Qj_iplus_JP3.N4, 2));

    IS_Qjnm[2] = (13.0 / 12.0)*(pow(Qj_iminus_JM3.N4 - 2.0*Qj_iminus_JM2.N4 + Qj_iminus_JM1.N4, 2)) + (1.0 / 4.0)*(pow(Qj_iminus_JM3.N4 - 4.0*Qj_iminus_JM2.N4 + 3.0*Qj_iminus_JM1.N4, 2));
    IS_Qjnm[1] = (13.0 / 12.0)*(pow(Qj_iminus_JM2.N4 - 2.0*Qj_iminus_JM1.N4 + Qj_iminus_I.N4, 2)) + (1.0 / 4.0)*(pow(Qj_iminus_JM2.N4 - Qj_iminus_I.N4, 2));
    IS_Qjnm[0] = (13.0 / 12.0)*(pow(Qj_iminus_JM1.N4 - 2.0*Qj_iminus_I.N4 + Qj_iminus_JP1.N4, 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iminus_JM1.N4 - 4.0*Qj_iminus_I.N4 + Qj_iminus_JP1.N4, 2));

    IS_Qjnp[2] = (13.0 / 12.0)*(pow(Qj_iminus_JM2.N4 - 2.0*Qj_iminus_JM1.N4 + Qj_iminus_I.N4, 2)) + (1.0 / 4.0)*(pow(Qj_iminus_JM2.N4 - 4.0*Qj_iminus_JM1.N4 + 3.0*Qj_iminus_I.N4, 2));
    IS_Qjnp[1] = (13.0 / 12.0)*(pow(Qj_iminus_JM1.N4 - 2.0*Qj_iminus_I.N4 + Qj_iminus_JP1.N4, 2)) + (1.0 / 4.0)*(pow(Qj_iminus_JM1.N4 - Qj_iminus_JP1.N4, 2));
    IS_Qjnp[0] = (13.0 / 12.0)*(pow(Qj_iminus_I.N4 - 2.0*Qj_iminus_JP1.N4 + Qj_iminus_JP2.N4, 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iminus_I.N4 - 4.0*Qj_iminus_JP1.N4 + Qj_iminus_JP2.N4, 2));

    IS_Qkm[2] = (13.0 / 12.0)*(pow(Qk_iplus_KM2.N4 - 2.0*Qk_iplus_KM1.N4 + Qk_iplus_I.N4, 2)) + (1.0 / 4.0)*(pow(Qk_iplus_KM2.N4 - 4.0*Qk_iplus_KM1.N4 + 3.0*Qk_iplus_I.N4, 2));
    IS_Qkm[1] = (13.0 / 12.0)*(pow(Qk_iplus_KM1.N4 - 2.0*Qk_iplus_I.N4 + Qk_iplus_KP1.N4, 2)) + (1.0 / 4.0)*(pow(Qk_iplus_KM1.N4 - Qk_iplus_KP1.N4, 2));
    IS_Qkm[0] = (13.0 / 12.0)*(pow(Qk_iplus_I.N4 - 2.0*Qk_iplus_KP1.N4 + Qk_iplus_KP2.N4, 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iplus_I.N4 - 4.0*Qk_iplus_KP1.N4 + Qk_iplus_KP2.N4, 2));

    IS_Qkp[2] = (13.0 / 12.0)*(pow(Qk_iplus_KM1.N4 - 2.0*Qk_iplus_I.N4 + Qk_iplus_KP1.N4, 2)) + (1.0 / 4.0)*(pow(Qk_iplus_KM1.N4 - 4.0*Qk_iplus_I.N4 + 3.0*Qk_iplus_KP1.N4, 2));
    IS_Qkp[1] = (13.0 / 12.0)*(pow(Qk_iplus_I.N4 - 2.0*Qk_iplus_KP1.N4 + Qk_iplus_KP2.N4, 2)) + (1.0 / 4.0)*(pow(Qk_iplus_I.N4 - Qk_iplus_KP2.N4, 2));
    IS_Qkp[0] = (13.0 / 12.0)*(pow(Qk_iplus_KP1.N4 - 2.0*Qk_iplus_KP2.N4 + Qk_iplus_KP3.N4, 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iplus_KP1.N4 - 4.0*Qk_iplus_KP2.N4 + Qk_iplus_KP3.N4, 2));

    IS_Qknm[2] = (13.0 / 12.0)*(pow(Qk_iminus_KM3.N4 - 2.0*Qk_iminus_KM2.N4 + Qk_iminus_KM1.N4, 2)) + (1.0 / 4.0)*(pow(Qk_iminus_KM3.N4 - 4.0*Qk_iminus_KM2.N4 + 3.0*Qk_iminus_KM1.N4, 2));
    IS_Qknm[1] = (13.0 / 12.0)*(pow(Qk_iminus_KM2.N4 - 2.0*Qk_iminus_KM1.N4 + Qk_iminus_I.N4, 2)) + (1.0 / 4.0)*(pow(Qk_iminus_KM2.N4 - Qk_iminus_I.N4, 2));
    IS_Qknm[0] = (13.0 / 12.0)*(pow(Qk_iminus_KM1.N4 - 2.0*Qk_iminus_I.N4 + Qk_iminus_KP1.N4, 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iminus_KM1.N4 - 4.0*Qk_iminus_I.N4 + Qk_iminus_KP1.N4, 2));

    IS_Qknp[2] = (13.0 / 12.0)*(pow(Qk_iminus_KM2.N4 - 2.0*Qk_iminus_KM1.N4 + Qk_iminus_I.N4, 2)) + (1.0 / 4.0)*(pow(Qk_iminus_KM2.N4 - 4.0*Qk_iminus_KM1.N4 + 3.0*Qk_iminus_I.N4, 2));
    IS_Qknp[1] = (13.0 / 12.0)*(pow(Qk_iminus_KM1.N4 - 2.0*Qk_iminus_I.N4 + Qk_iminus_KP1.N4, 2)) + (1.0 / 4.0)*(pow(Qk_iminus_KM1.N4 - Qk_iminus_KP1.N4, 2));
    IS_Qknp[0] = (13.0 / 12.0)*(pow(Qk_iminus_I.N4 - 2.0*Qk_iminus_KP1.N4 + Qk_iminus_KP2.N4, 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iminus_I.N4 - 4.0*Qk_iminus_KP1.N4 + Qk_iminus_KP2.N4, 2));

    w_Qip[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qip[0]), 2.0));
    w_Qip[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qip[1]), 2.0));
    w_Qip[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qip[2]), 2.0));

    w_Qim[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qim[0]), 2.0));
    w_Qim[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qim[1]), 2.0));
    w_Qim[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qim[2]), 2.0));

    w_Qinp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qinp[0]), 2.0));
    w_Qinp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qinp[1]), 2.0));
    w_Qinp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qinp[2]), 2.0));

    w_Qinm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qinm[0]), 2.0));
    w_Qinm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qinm[1]), 2.0));
    w_Qinm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qinm[2]), 2.0));

    w_Qjp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjp[0]), 2.0));
    w_Qjp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjp[1]), 2.0));
    w_Qjp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjp[2]), 2.0));

    w_Qjm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjm[0]), 2.0));
    w_Qjm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjm[1]), 2.0));
    w_Qjm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjm[2]), 2.0));

    w_Qjnp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjnp[0]), 2.0));
    w_Qjnp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjnp[1]), 2.0));
    w_Qjnp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjnp[2]), 2.0));

    w_Qjnm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjnm[0]), 2.0));
    w_Qjnm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjnm[1]), 2.0));
    w_Qjnm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjnm[2]), 2.0));

    w_Qkp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qkp[0]), 2.0));
    w_Qkp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qkp[1]), 2.0));
    w_Qkp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qkp[2]), 2.0));

    w_Qkm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qkm[0]), 2.0));
    w_Qkm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qkm[1]), 2.0));
    w_Qkm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qkm[2]), 2.0));

    w_Qknp[0] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qknp[0]), 2.0));
    w_Qknp[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qknp[1]), 2.0));
    w_Qknp[2] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qknp[2]), 2.0));

    w_Qknm[0] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qknm[0]), 2.0));
    w_Qknm[1] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qknm[1]), 2.0));
    w_Qknm[2] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qknm[2]), 2.0));

    W_Qip[0] = w_Qip[0] / (w_Qip[0] + w_Qip[1] + w_Qip[2]);
    W_Qip[1] = w_Qip[1] / (w_Qip[0] + w_Qip[1] + w_Qip[2]);
    W_Qip[2] = w_Qip[2] / (w_Qip[0] + w_Qip[1] + w_Qip[2]);

    W_Qinp[0] = w_Qinp[0] / (w_Qinp[0] + w_Qinp[1] + w_Qinp[2]);
    W_Qinp[1] = w_Qinp[1] / (w_Qinp[0] + w_Qinp[1] + w_Qinp[2]);
    W_Qinp[2] = w_Qinp[2] / (w_Qinp[0] + w_Qinp[1] + w_Qinp[2]);

    W_Qim[0] = w_Qim[0] / (w_Qim[0] + w_Qim[1] + w_Qim[2]);
    W_Qim[1] = w_Qim[1] / (w_Qim[0] + w_Qim[1] + w_Qim[2]);
    W_Qim[2] = w_Qim[2] / (w_Qim[0] + w_Qim[1] + w_Qim[2]);

    W_Qinm[0] = w_Qinm[0] / (w_Qinm[0] + w_Qinm[1] + w_Qinm[2]);
    W_Qinm[1] = w_Qinm[1] / (w_Qinm[0] + w_Qinm[1] + w_Qinm[2]);
    W_Qinm[2] = w_Qinm[2] / (w_Qinm[0] + w_Qinm[1] + w_Qinm[2]);

    W_Qjp[0] = w_Qjp[0] / (w_Qjp[0] + w_Qjp[1] + w_Qjp[2]);
    W_Qjp[1] = w_Qjp[1] / (w_Qjp[0] + w_Qjp[1] + w_Qjp[2]);
    W_Qjp[2] = w_Qjp[2] / (w_Qjp[0] + w_Qjp[1] + w_Qjp[2]);

    W_Qjnp[0] = w_Qjnp[0] / (w_Qjnp[0] + w_Qjnp[1] + w_Qjnp[2]);
    W_Qjnp[1] = w_Qjnp[1] / (w_Qjnp[0] + w_Qjnp[1] + w_Qjnp[2]);
    W_Qjnp[2] = w_Qjnp[2] / (w_Qjnp[0] + w_Qjnp[1] + w_Qjnp[2]);

    W_Qjm[0] = w_Qjm[0] / (w_Qjm[0] + w_Qjm[1] + w_Qjm[2]);
    W_Qjm[1] = w_Qjm[1] / (w_Qjm[0] + w_Qjm[1] + w_Qjm[2]);
    W_Qjm[2] = w_Qjm[2] / (w_Qjm[0] + w_Qjm[1] + w_Qjm[2]);

    W_Qjnm[0] = w_Qjnm[0] / (w_Qjnm[0] + w_Qjnm[1] + w_Qjnm[2]);
    W_Qjnm[1] = w_Qjnm[1] / (w_Qjnm[0] + w_Qjnm[1] + w_Qjnm[2]);
    W_Qjnm[2] = w_Qjnm[2] / (w_Qjnm[0] + w_Qjnm[1] + w_Qjnm[2]);

    W_Qkp[0] = w_Qkp[0] / (w_Qkp[0] + w_Qkp[1] + w_Qkp[2]);
    W_Qkp[1] = w_Qkp[1] / (w_Qkp[0] + w_Qkp[1] + w_Qkp[2]);
    W_Qkp[2] = w_Qkp[2] / (w_Qkp[0] + w_Qkp[1] + w_Qkp[2]);

    W_Qknp[0] = w_Qknp[0] / (w_Qknp[0] + w_Qknp[1] + w_Qknp[2]);
    W_Qknp[1] = w_Qknp[1] / (w_Qknp[0] + w_Qknp[1] + w_Qknp[2]);
    W_Qknp[2] = w_Qknp[2] / (w_Qknp[0] + w_Qknp[1] + w_Qknp[2]);

    W_Qkm[0] = w_Qkm[0] / (w_Qkm[0] + w_Qkm[1] + w_Qkm[2]);
    W_Qkm[1] = w_Qkm[1] / (w_Qkm[0] + w_Qkm[1] + w_Qkm[2]);
    W_Qkm[2] = w_Qkm[2] / (w_Qkm[0] + w_Qkm[1] + w_Qkm[2]);

    W_Qknm[0] = w_Qknm[0] / (w_Qknm[0] + w_Qknm[1] + w_Qknm[2]);
    W_Qknm[1] = w_Qknm[1] / (w_Qknm[0] + w_Qknm[1] + w_Qknm[2]);
    W_Qknm[2] = w_Qknm[2] / (w_Qknm[0] + w_Qknm[1] + w_Qknm[2]);

    Qi_iplus_half_pos_char[4] = W_Qip[0] * Qi_half_p[0] + W_Qip[1] * Qi_half_p[1] + W_Qip[2] * Qi_half_p[2];
    Qi_iminus_half_pos_char[4] = W_Qinp[0] * Qi_half_np[0] + W_Qinp[1] * Qi_half_np[1] + W_Qinp[2] * Qi_half_np[2];

    Qi_iplus_half_neg_char[4] = W_Qim[0] * Qi_half_m[0] + W_Qim[1] * Qi_half_m[1] + W_Qim[2] * Qi_half_m[2];
    Qi_iminus_half_neg_char[4] = W_Qinm[0] * Qi_halfn_m[0] + W_Qinm[1] * Qi_halfn_m[1] + W_Qinm[2] * Qi_halfn_m[2];

    Qj_iplus_half_pos_char[4] = W_Qjp[0] * Qj_half_p[0] + W_Qjp[1] * Qj_half_p[1] + W_Qjp[2] * Qj_half_p[2];
    Qj_iminus_half_pos_char[4] = W_Qjnp[0] * Qj_half_np[0] + W_Qjnp[1] * Qj_half_np[1] + W_Qjnp[2] * Qj_half_np[2];

    Qj_iplus_half_neg_char[4] = W_Qjm[0] * Qj_half_m[0] + W_Qjm[1] * Qj_half_m[1] + W_Qjm[2] * Qj_half_m[2];
    Qj_iminus_half_neg_char[4] = W_Qjnm[0] * Qj_halfn_m[0] + W_Qjnm[1] * Qj_halfn_m[1] + W_Qjnm[2] * Qj_halfn_m[2];

    Qk_iplus_half_pos_char[4] = W_Qkp[0] * Qk_half_p[0] + W_Qkp[1] * Qk_half_p[1] + W_Qkp[2] * Qk_half_p[2];
    Qk_iminus_half_pos_char[4] = W_Qknp[0] * Qk_half_np[0] + W_Qknp[1] * Qk_half_np[1] + W_Qknp[2] * Qk_half_np[2];

    Qk_iplus_half_neg_char[4] = W_Qkm[0] * Qk_half_m[0] + W_Qkm[1] * Qk_half_m[1] + W_Qkm[2] * Qk_half_m[2];
    Qk_iminus_half_neg_char[4] = W_Qknm[0] * Qk_halfn_m[0] + W_Qknm[1] * Qk_halfn_m[1] + W_Qknm[2] * Qk_halfn_m[2];


    /*******************Transfer from characteristic to computational space**********************/

    /**********************************Equation 1 ***********************************************/
    Qi_iplus_half_pos[0] = r_eigen_Qip[0].N0 * (Qi_iplus_half_pos_char[0]) + r_eigen_Qip[0].N1 * (Qi_iplus_half_pos_char[1]) + r_eigen_Qip[0].N2 * (Qi_iplus_half_pos_char[2]) + r_eigen_Qip[0].N3 * (Qi_iplus_half_pos_char[3]) + r_eigen_Qip[0].N4 * (Qi_iplus_half_pos_char[4]);
    Qi_iplus_half_neg[0] = r_eigen_Qip[0].N0 * (Qi_iplus_half_neg_char[0]) + r_eigen_Qip[0].N1 * (Qi_iplus_half_neg_char[1]) + r_eigen_Qip[0].N2 * (Qi_iplus_half_neg_char[2]) + r_eigen_Qip[0].N3 * (Qi_iplus_half_neg_char[3]) + r_eigen_Qip[0].N4 * (Qi_iplus_half_neg_char[4]);

    Qi_iminus_half_pos[0] = r_eigen_Qim[0].N0 * (Qi_iminus_half_pos_char[0]) + r_eigen_Qim[0].N1 * (Qi_iminus_half_pos_char[1]) + r_eigen_Qim[0].N2 * (Qi_iminus_half_pos_char[2]) + r_eigen_Qim[0].N3 * (Qi_iminus_half_pos_char[3]) + r_eigen_Qim[0].N4 * (Qi_iminus_half_pos_char[4]);
    Qi_iminus_half_neg[0] = r_eigen_Qim[0].N0 * (Qi_iminus_half_neg_char[0]) + r_eigen_Qim[0].N1 * (Qi_iminus_half_neg_char[1]) + r_eigen_Qim[0].N2 * (Qi_iminus_half_neg_char[2]) + r_eigen_Qim[0].N3 * (Qi_iminus_half_neg_char[3]) + r_eigen_Qim[0].N4 * (Qi_iminus_half_neg_char[4]);

    Qj_iplus_half_pos[0] = r_eigen_Qjp[0].N0 * (Qj_iplus_half_pos_char[0]) + r_eigen_Qjp[0].N1 * (Qj_iplus_half_pos_char[1]) + r_eigen_Qjp[0].N2 * (Qj_iplus_half_pos_char[2]) + r_eigen_Qjp[0].N3 * (Qj_iplus_half_pos_char[3]) + r_eigen_Qjp[0].N4 * (Qj_iplus_half_pos_char[4]);
    Qj_iplus_half_neg[0] = r_eigen_Qjp[0].N0 * (Qj_iplus_half_neg_char[0]) + r_eigen_Qjp[0].N1 * (Qj_iplus_half_neg_char[1]) + r_eigen_Qjp[0].N2 * (Qj_iplus_half_neg_char[2]) + r_eigen_Qjp[0].N3 * (Qj_iplus_half_neg_char[3]) + r_eigen_Qjp[0].N4 * (Qj_iplus_half_neg_char[4]);

    Qj_iminus_half_pos[0] = r_eigen_Qjm[0].N0 * (Qj_iminus_half_pos_char[0]) + r_eigen_Qjm[0].N1 * (Qj_iminus_half_pos_char[1]) + r_eigen_Qjm[0].N2 * (Qj_iminus_half_pos_char[2]) + r_eigen_Qjm[0].N3 * (Qj_iminus_half_pos_char[3]) + r_eigen_Qjm[0].N4 * (Qj_iminus_half_pos_char[4]);
    Qj_iminus_half_neg[0] = r_eigen_Qjm[0].N0 * (Qj_iminus_half_neg_char[0]) + r_eigen_Qjm[0].N1 * (Qj_iminus_half_neg_char[1]) + r_eigen_Qjm[0].N2 * (Qj_iminus_half_neg_char[2]) + r_eigen_Qjm[0].N3 * (Qj_iminus_half_neg_char[3]) + r_eigen_Qjm[0].N4 * (Qj_iminus_half_neg_char[4]);

    Qk_iplus_half_pos[0] = r_eigen_Qkp[0].N0 * (Qk_iplus_half_pos_char[0]) + r_eigen_Qkp[0].N1 * (Qk_iplus_half_pos_char[1]) + r_eigen_Qkp[0].N2 * (Qk_iplus_half_pos_char[2]) + r_eigen_Qkp[0].N3 * (Qk_iplus_half_pos_char[3]) + r_eigen_Qkp[0].N4 * (Qk_iplus_half_pos_char[4]);
    Qk_iplus_half_neg[0] = r_eigen_Qkp[0].N0 * (Qk_iplus_half_neg_char[0]) + r_eigen_Qkp[0].N1 * (Qk_iplus_half_neg_char[1]) + r_eigen_Qkp[0].N2 * (Qk_iplus_half_neg_char[2]) + r_eigen_Qkp[0].N3 * (Qk_iplus_half_neg_char[3]) + r_eigen_Qkp[0].N4 * (Qk_iplus_half_neg_char[4]);

    Qk_iminus_half_pos[0] = r_eigen_Qkm[0].N0 * (Qk_iminus_half_pos_char[0]) + r_eigen_Qkm[0].N1 * (Qk_iminus_half_pos_char[1]) + r_eigen_Qkm[0].N2 * (Qk_iminus_half_pos_char[2]) + r_eigen_Qkm[0].N3 * (Qk_iminus_half_pos_char[3]) + r_eigen_Qkm[0].N4 * (Qk_iminus_half_pos_char[4]);
    Qk_iminus_half_neg[0] = r_eigen_Qkm[0].N0 * (Qk_iminus_half_neg_char[0]) + r_eigen_Qkm[0].N1 * (Qk_iminus_half_neg_char[1]) + r_eigen_Qkm[0].N2 * (Qk_iminus_half_neg_char[2]) + r_eigen_Qkm[0].N3 * (Qk_iminus_half_neg_char[3]) + r_eigen_Qkm[0].N4 * (Qk_iminus_half_neg_char[4]);


    /**********************************Equation 2 ***********************************************/
    Qi_iplus_half_pos[1] = r_eigen_Qip[1].N0 * (Qi_iplus_half_pos_char[0]) + r_eigen_Qip[1].N1 * (Qi_iplus_half_pos_char[1]) + r_eigen_Qip[1].N2 * (Qi_iplus_half_pos_char[2]) + r_eigen_Qip[1].N3 * (Qi_iplus_half_pos_char[3]) + r_eigen_Qip[1].N4 * (Qi_iplus_half_pos_char[4]);
    Qi_iplus_half_neg[1] = r_eigen_Qip[1].N0 * (Qi_iplus_half_neg_char[0]) + r_eigen_Qip[1].N1 * (Qi_iplus_half_neg_char[1]) + r_eigen_Qip[1].N2 * (Qi_iplus_half_neg_char[2]) + r_eigen_Qip[1].N3 * (Qi_iplus_half_neg_char[3]) + r_eigen_Qip[1].N4 * (Qi_iplus_half_neg_char[4]);

    Qi_iminus_half_pos[1] = r_eigen_Qim[1].N0 * (Qi_iminus_half_pos_char[0]) + r_eigen_Qim[1].N1 * (Qi_iminus_half_pos_char[1]) + r_eigen_Qim[1].N2 * (Qi_iminus_half_pos_char[2]) + r_eigen_Qim[1].N3 * (Qi_iminus_half_pos_char[3]) + r_eigen_Qim[1].N4 * (Qi_iminus_half_pos_char[4]);
    Qi_iminus_half_neg[1] = r_eigen_Qim[1].N0 * (Qi_iminus_half_neg_char[0]) + r_eigen_Qim[1].N1 * (Qi_iminus_half_neg_char[1]) + r_eigen_Qim[1].N2 * (Qi_iminus_half_neg_char[2]) + r_eigen_Qim[1].N3 * (Qi_iminus_half_neg_char[3]) + r_eigen_Qim[1].N4 * (Qi_iminus_half_neg_char[4]);

    Qj_iplus_half_pos[1] = r_eigen_Qjp[1].N0 * (Qj_iplus_half_pos_char[0]) + r_eigen_Qjp[1].N1 * (Qj_iplus_half_pos_char[1]) + r_eigen_Qjp[1].N2 * (Qj_iplus_half_pos_char[2]) + r_eigen_Qjp[1].N3 * (Qj_iplus_half_pos_char[3]) + r_eigen_Qjp[1].N4 * (Qj_iplus_half_pos_char[4]);
    Qj_iplus_half_neg[1] = r_eigen_Qjp[1].N0 * (Qj_iplus_half_neg_char[0]) + r_eigen_Qjp[1].N1 * (Qj_iplus_half_neg_char[1]) + r_eigen_Qjp[1].N2 * (Qj_iplus_half_neg_char[2]) + r_eigen_Qjp[1].N3 * (Qj_iplus_half_neg_char[3]) + r_eigen_Qjp[1].N4 * (Qj_iplus_half_neg_char[4]);

    Qj_iminus_half_pos[1] = r_eigen_Qjm[1].N0 * (Qj_iminus_half_pos_char[0]) + r_eigen_Qjm[1].N1 * (Qj_iminus_half_pos_char[1]) + r_eigen_Qjm[1].N2 * (Qj_iminus_half_pos_char[2]) + r_eigen_Qjm[1].N3 * (Qj_iminus_half_pos_char[3]) + r_eigen_Qjm[1].N4 * (Qj_iminus_half_pos_char[4]);
    Qj_iminus_half_neg[1] = r_eigen_Qjm[1].N0 * (Qj_iminus_half_neg_char[0]) + r_eigen_Qjm[1].N1 * (Qj_iminus_half_neg_char[1]) + r_eigen_Qjm[1].N2 * (Qj_iminus_half_neg_char[2]) + r_eigen_Qjm[1].N3 * (Qj_iminus_half_neg_char[3]) + r_eigen_Qjm[1].N4 * (Qj_iminus_half_neg_char[4]);

    Qk_iplus_half_pos[1] = r_eigen_Qkp[1].N0 * (Qk_iplus_half_pos_char[0]) + r_eigen_Qkp[1].N1 * (Qk_iplus_half_pos_char[1]) + r_eigen_Qkp[1].N2 * (Qk_iplus_half_pos_char[2]) + r_eigen_Qkp[1].N3 * (Qk_iplus_half_pos_char[3]) + r_eigen_Qkp[1].N4 * (Qk_iplus_half_pos_char[4]);
    Qk_iplus_half_neg[1] = r_eigen_Qkp[1].N0 * (Qk_iplus_half_neg_char[0]) + r_eigen_Qkp[1].N1 * (Qk_iplus_half_neg_char[1]) + r_eigen_Qkp[1].N2 * (Qk_iplus_half_neg_char[2]) + r_eigen_Qkp[1].N3 * (Qk_iplus_half_neg_char[3]) + r_eigen_Qkp[1].N4 * (Qk_iplus_half_neg_char[4]);

    Qk_iminus_half_pos[1] = r_eigen_Qkm[1].N0 * (Qk_iminus_half_pos_char[0]) + r_eigen_Qkm[1].N1 * (Qk_iminus_half_pos_char[1]) + r_eigen_Qkm[1].N2 * (Qk_iminus_half_pos_char[2]) + r_eigen_Qkm[1].N3 * (Qk_iminus_half_pos_char[3]) + r_eigen_Qkm[1].N4 * (Qk_iminus_half_pos_char[4]);
    Qk_iminus_half_neg[1] = r_eigen_Qkm[1].N0 * (Qk_iminus_half_neg_char[0]) + r_eigen_Qkm[1].N1 * (Qk_iminus_half_neg_char[1]) + r_eigen_Qkm[1].N2 * (Qk_iminus_half_neg_char[2]) + r_eigen_Qkm[1].N3 * (Qk_iminus_half_neg_char[3]) + r_eigen_Qkm[1].N4 * (Qk_iminus_half_neg_char[4]);


    /**********************************Equation 3 ***********************************************/
    Qi_iplus_half_pos[2] = r_eigen_Qip[2].N0 * (Qi_iplus_half_pos_char[0]) + r_eigen_Qip[2].N1 * (Qi_iplus_half_pos_char[1]) + r_eigen_Qip[2].N2 * (Qi_iplus_half_pos_char[2]) + r_eigen_Qip[2].N3 * (Qi_iplus_half_pos_char[3]) + r_eigen_Qip[2].N4 * (Qi_iplus_half_pos_char[4]);
    Qi_iplus_half_neg[2] = r_eigen_Qip[2].N0 * (Qi_iplus_half_neg_char[0]) + r_eigen_Qip[2].N1 * (Qi_iplus_half_neg_char[1]) + r_eigen_Qip[2].N2 * (Qi_iplus_half_neg_char[2]) + r_eigen_Qip[2].N3 * (Qi_iplus_half_neg_char[3]) + r_eigen_Qip[2].N4 * (Qi_iplus_half_neg_char[4]);

    Qi_iminus_half_pos[2] = r_eigen_Qim[2].N0 * (Qi_iminus_half_pos_char[0]) + r_eigen_Qim[2].N1 * (Qi_iminus_half_pos_char[1]) + r_eigen_Qim[2].N2 * (Qi_iminus_half_pos_char[2]) + r_eigen_Qim[2].N3 * (Qi_iminus_half_pos_char[3]) + r_eigen_Qim[2].N4 * (Qi_iminus_half_pos_char[4]);
    Qi_iminus_half_neg[2] = r_eigen_Qim[2].N0 * (Qi_iminus_half_neg_char[0]) + r_eigen_Qim[2].N1 * (Qi_iminus_half_neg_char[1]) + r_eigen_Qim[2].N2 * (Qi_iminus_half_neg_char[2]) + r_eigen_Qim[2].N3 * (Qi_iminus_half_neg_char[3]) + r_eigen_Qim[2].N4 * (Qi_iminus_half_neg_char[4]);

    Qj_iplus_half_pos[2] = r_eigen_Qjp[2].N0 * (Qj_iplus_half_pos_char[0]) + r_eigen_Qjp[2].N1 * (Qj_iplus_half_pos_char[1]) + r_eigen_Qjp[2].N2 * (Qj_iplus_half_pos_char[2]) + r_eigen_Qjp[2].N3 * (Qj_iplus_half_pos_char[3]) + r_eigen_Qjp[2].N4 * (Qj_iplus_half_pos_char[4]);
    Qj_iplus_half_neg[2] = r_eigen_Qjp[2].N0 * (Qj_iplus_half_neg_char[0]) + r_eigen_Qjp[2].N1 * (Qj_iplus_half_neg_char[1]) + r_eigen_Qjp[2].N2 * (Qj_iplus_half_neg_char[2]) + r_eigen_Qjp[2].N3 * (Qj_iplus_half_neg_char[3]) + r_eigen_Qjp[2].N4 * (Qj_iplus_half_neg_char[4]);

    Qj_iminus_half_pos[2] = r_eigen_Qjm[2].N0 * (Qj_iminus_half_pos_char[0]) + r_eigen_Qjm[2].N1 * (Qj_iminus_half_pos_char[1]) + r_eigen_Qjm[2].N2 * (Qj_iminus_half_pos_char[2]) + r_eigen_Qjm[2].N3 * (Qj_iminus_half_pos_char[3]) + r_eigen_Qjm[2].N4 * (Qj_iminus_half_pos_char[4]);
    Qj_iminus_half_neg[2] = r_eigen_Qjm[2].N0 * (Qj_iminus_half_neg_char[0]) + r_eigen_Qjm[2].N1 * (Qj_iminus_half_neg_char[1]) + r_eigen_Qjm[2].N2 * (Qj_iminus_half_neg_char[2]) + r_eigen_Qjm[2].N3 * (Qj_iminus_half_neg_char[3]) + r_eigen_Qjm[2].N4 * (Qj_iminus_half_neg_char[4]);

    Qk_iplus_half_pos[2] = r_eigen_Qkp[2].N0 * (Qk_iplus_half_pos_char[0]) + r_eigen_Qkp[2].N1 * (Qk_iplus_half_pos_char[1]) + r_eigen_Qkp[2].N2 * (Qk_iplus_half_pos_char[2]) + r_eigen_Qkp[2].N3 * (Qk_iplus_half_pos_char[3]) + r_eigen_Qkp[2].N4 * (Qk_iplus_half_pos_char[4]);
    Qk_iplus_half_neg[2] = r_eigen_Qkp[2].N0 * (Qk_iplus_half_neg_char[0]) + r_eigen_Qkp[2].N1 * (Qk_iplus_half_neg_char[1]) + r_eigen_Qkp[2].N2 * (Qk_iplus_half_neg_char[2]) + r_eigen_Qkp[2].N3 * (Qk_iplus_half_neg_char[3]) + r_eigen_Qkp[2].N4 * (Qk_iplus_half_neg_char[4]);

    Qk_iminus_half_pos[2] = r_eigen_Qkm[2].N0 * (Qk_iminus_half_pos_char[0]) + r_eigen_Qkm[2].N1 * (Qk_iminus_half_pos_char[1]) + r_eigen_Qkm[2].N2 * (Qk_iminus_half_pos_char[2]) + r_eigen_Qkm[2].N3 * (Qk_iminus_half_pos_char[3]) + r_eigen_Qkm[2].N4 * (Qk_iminus_half_pos_char[4]);
    Qk_iminus_half_neg[2] = r_eigen_Qkm[2].N0 * (Qk_iminus_half_neg_char[0]) + r_eigen_Qkm[2].N1 * (Qk_iminus_half_neg_char[1]) + r_eigen_Qkm[2].N2 * (Qk_iminus_half_neg_char[2]) + r_eigen_Qkm[2].N3 * (Qk_iminus_half_neg_char[3]) + r_eigen_Qkm[2].N4 * (Qk_iminus_half_neg_char[4]);


    /**********************************Equation 4 ***********************************************/
    Qi_iplus_half_pos[3] = r_eigen_Qip[3].N0 * (Qi_iplus_half_pos_char[0]) + r_eigen_Qip[3].N1 * (Qi_iplus_half_pos_char[1]) + r_eigen_Qip[3].N2 * (Qi_iplus_half_pos_char[2]) + r_eigen_Qip[3].N3 * (Qi_iplus_half_pos_char[3]) + r_eigen_Qip[3].N4 * (Qi_iplus_half_pos_char[4]);
    Qi_iplus_half_neg[3] = r_eigen_Qip[3].N0 * (Qi_iplus_half_neg_char[0]) + r_eigen_Qip[3].N1 * (Qi_iplus_half_neg_char[1]) + r_eigen_Qip[3].N2 * (Qi_iplus_half_neg_char[2]) + r_eigen_Qip[3].N3 * (Qi_iplus_half_neg_char[3]) + r_eigen_Qip[3].N4 * (Qi_iplus_half_neg_char[4]);

    Qi_iminus_half_pos[3] = r_eigen_Qim[3].N0 * (Qi_iminus_half_pos_char[0]) + r_eigen_Qim[3].N1 * (Qi_iminus_half_pos_char[1]) + r_eigen_Qim[3].N2 * (Qi_iminus_half_pos_char[2]) + r_eigen_Qim[3].N3 * (Qi_iminus_half_pos_char[3]) + r_eigen_Qim[3].N4 * (Qi_iminus_half_pos_char[4]);
    Qi_iminus_half_neg[3] = r_eigen_Qim[3].N0 * (Qi_iminus_half_neg_char[0]) + r_eigen_Qim[3].N1 * (Qi_iminus_half_neg_char[1]) + r_eigen_Qim[3].N2 * (Qi_iminus_half_neg_char[2]) + r_eigen_Qim[3].N3 * (Qi_iminus_half_neg_char[3]) + r_eigen_Qim[3].N4 * (Qi_iminus_half_neg_char[4]);

    Qj_iplus_half_pos[3] = r_eigen_Qjp[3].N0 * (Qj_iplus_half_pos_char[0]) + r_eigen_Qjp[3].N1 * (Qj_iplus_half_pos_char[1]) + r_eigen_Qjp[3].N2 * (Qj_iplus_half_pos_char[2]) + r_eigen_Qjp[3].N3 * (Qj_iplus_half_pos_char[3]) + r_eigen_Qjp[3].N4 * (Qj_iplus_half_pos_char[4]);
    Qj_iplus_half_neg[3] = r_eigen_Qjp[3].N0 * (Qj_iplus_half_neg_char[0]) + r_eigen_Qjp[3].N1 * (Qj_iplus_half_neg_char[1]) + r_eigen_Qjp[3].N2 * (Qj_iplus_half_neg_char[2]) + r_eigen_Qjp[3].N3 * (Qj_iplus_half_neg_char[3]) + r_eigen_Qjp[3].N4 * (Qj_iplus_half_neg_char[4]);

    Qj_iminus_half_pos[3] = r_eigen_Qjm[3].N0 * (Qj_iminus_half_pos_char[0]) + r_eigen_Qjm[3].N1 * (Qj_iminus_half_pos_char[1]) + r_eigen_Qjm[3].N2 * (Qj_iminus_half_pos_char[2]) + r_eigen_Qjm[3].N3 * (Qj_iminus_half_pos_char[3]) + r_eigen_Qjm[3].N4 * (Qj_iminus_half_pos_char[4]);
    Qj_iminus_half_neg[3] = r_eigen_Qjm[3].N0 * (Qj_iminus_half_neg_char[0]) + r_eigen_Qjm[3].N1 * (Qj_iminus_half_neg_char[1]) + r_eigen_Qjm[3].N2 * (Qj_iminus_half_neg_char[2]) + r_eigen_Qjm[3].N3 * (Qj_iminus_half_neg_char[3]) + r_eigen_Qjm[3].N4 * (Qj_iminus_half_neg_char[4]);

    Qk_iplus_half_pos[3] = r_eigen_Qkp[3].N0 * (Qk_iplus_half_pos_char[0]) + r_eigen_Qkp[3].N1 * (Qk_iplus_half_pos_char[1]) + r_eigen_Qkp[3].N2 * (Qk_iplus_half_pos_char[2]) + r_eigen_Qkp[3].N3 * (Qk_iplus_half_pos_char[3]) + r_eigen_Qkp[3].N4 * (Qk_iplus_half_pos_char[4]);
    Qk_iplus_half_neg[3] = r_eigen_Qkp[3].N0 * (Qk_iplus_half_neg_char[0]) + r_eigen_Qkp[3].N1 * (Qk_iplus_half_neg_char[1]) + r_eigen_Qkp[3].N2 * (Qk_iplus_half_neg_char[2]) + r_eigen_Qkp[3].N3 * (Qk_iplus_half_neg_char[3]) + r_eigen_Qkp[3].N4 * (Qk_iplus_half_neg_char[4]);

    Qk_iminus_half_pos[3] = r_eigen_Qkm[3].N0 * (Qk_iminus_half_pos_char[0]) + r_eigen_Qkm[3].N1 * (Qk_iminus_half_pos_char[1]) + r_eigen_Qkm[3].N2 * (Qk_iminus_half_pos_char[2]) + r_eigen_Qkm[3].N3 * (Qk_iminus_half_pos_char[3]) + r_eigen_Qkm[3].N4 * (Qk_iminus_half_pos_char[4]);
    Qk_iminus_half_neg[3] = r_eigen_Qkm[3].N0 * (Qk_iminus_half_neg_char[0]) + r_eigen_Qkm[3].N1 * (Qk_iminus_half_neg_char[1]) + r_eigen_Qkm[3].N2 * (Qk_iminus_half_neg_char[2]) + r_eigen_Qkm[3].N3 * (Qk_iminus_half_neg_char[3]) + r_eigen_Qkm[3].N4 * (Qk_iminus_half_neg_char[4]);


    /**********************************Equation 5 ***********************************************/
    Qi_iplus_half_pos[4] = r_eigen_Qip[4].N0 * (Qi_iplus_half_pos_char[0]) + r_eigen_Qip[4].N1 * (Qi_iplus_half_pos_char[1]) + r_eigen_Qip[4].N2 * (Qi_iplus_half_pos_char[2]) + r_eigen_Qip[4].N3 * (Qi_iplus_half_pos_char[3]) + r_eigen_Qip[4].N4 * (Qi_iplus_half_pos_char[4]);
    Qi_iplus_half_neg[4] = r_eigen_Qip[4].N0 * (Qi_iplus_half_neg_char[0]) + r_eigen_Qip[4].N1 * (Qi_iplus_half_neg_char[1]) + r_eigen_Qip[4].N2 * (Qi_iplus_half_neg_char[2]) + r_eigen_Qip[4].N3 * (Qi_iplus_half_neg_char[3]) + r_eigen_Qip[4].N4 * (Qi_iplus_half_neg_char[4]);

    Qi_iminus_half_pos[4] = r_eigen_Qim[4].N0 * (Qi_iminus_half_pos_char[0]) + r_eigen_Qim[4].N1 * (Qi_iminus_half_pos_char[1]) + r_eigen_Qim[4].N2 * (Qi_iminus_half_pos_char[2]) + r_eigen_Qim[4].N3 * (Qi_iminus_half_pos_char[3]) + r_eigen_Qim[4].N4 * (Qi_iminus_half_pos_char[4]);
    Qi_iminus_half_neg[4] = r_eigen_Qim[4].N0 * (Qi_iminus_half_neg_char[0]) + r_eigen_Qim[4].N1 * (Qi_iminus_half_neg_char[1]) + r_eigen_Qim[4].N2 * (Qi_iminus_half_neg_char[2]) + r_eigen_Qim[4].N3 * (Qi_iminus_half_neg_char[3]) + r_eigen_Qim[4].N4 * (Qi_iminus_half_neg_char[4]);

    Qj_iplus_half_pos[4] = r_eigen_Qjp[4].N0 * (Qj_iplus_half_pos_char[0]) + r_eigen_Qjp[4].N1 * (Qj_iplus_half_pos_char[1]) + r_eigen_Qjp[4].N2 * (Qj_iplus_half_pos_char[2]) + r_eigen_Qjp[4].N3 * (Qj_iplus_half_pos_char[3]) + r_eigen_Qjp[4].N4 * (Qj_iplus_half_pos_char[4]);
    Qj_iplus_half_neg[4] = r_eigen_Qjp[4].N0 * (Qj_iplus_half_neg_char[0]) + r_eigen_Qjp[4].N1 * (Qj_iplus_half_neg_char[1]) + r_eigen_Qjp[4].N2 * (Qj_iplus_half_neg_char[2]) + r_eigen_Qjp[4].N3 * (Qj_iplus_half_neg_char[3]) + r_eigen_Qjp[4].N4 * (Qj_iplus_half_neg_char[4]);

    Qj_iminus_half_pos[4] = r_eigen_Qjm[4].N0 * (Qj_iminus_half_pos_char[0]) + r_eigen_Qjm[4].N1 * (Qj_iminus_half_pos_char[1]) + r_eigen_Qjm[4].N2 * (Qj_iminus_half_pos_char[2]) + r_eigen_Qjm[4].N3 * (Qj_iminus_half_pos_char[3]) + r_eigen_Qjm[4].N4 * (Qj_iminus_half_pos_char[4]);
    Qj_iminus_half_neg[4] = r_eigen_Qjm[4].N0 * (Qj_iminus_half_neg_char[0]) + r_eigen_Qjm[4].N1 * (Qj_iminus_half_neg_char[1]) + r_eigen_Qjm[4].N2 * (Qj_iminus_half_neg_char[2]) + r_eigen_Qjm[4].N3 * (Qj_iminus_half_neg_char[3]) + r_eigen_Qjm[4].N4 * (Qj_iminus_half_neg_char[4]);

    Qk_iplus_half_pos[4] = r_eigen_Qkp[4].N0 * (Qk_iplus_half_pos_char[0]) + r_eigen_Qkp[4].N1 * (Qk_iplus_half_pos_char[1]) + r_eigen_Qkp[4].N2 * (Qk_iplus_half_pos_char[2]) + r_eigen_Qkp[4].N3 * (Qk_iplus_half_pos_char[3]) + r_eigen_Qkp[4].N4 * (Qk_iplus_half_pos_char[4]);
    Qk_iplus_half_neg[4] = r_eigen_Qkp[4].N0 * (Qk_iplus_half_neg_char[0]) + r_eigen_Qkp[4].N1 * (Qk_iplus_half_neg_char[1]) + r_eigen_Qkp[4].N2 * (Qk_iplus_half_neg_char[2]) + r_eigen_Qkp[4].N3 * (Qk_iplus_half_neg_char[3]) + r_eigen_Qkp[4].N4 * (Qk_iplus_half_neg_char[4]);

    Qk_iminus_half_pos[4] = r_eigen_Qkm[4].N0 * (Qk_iminus_half_pos_char[0]) + r_eigen_Qkm[4].N1 * (Qk_iminus_half_pos_char[1]) + r_eigen_Qkm[4].N2 * (Qk_iminus_half_pos_char[2]) + r_eigen_Qkm[4].N3 * (Qk_iminus_half_pos_char[3]) + r_eigen_Qkm[4].N4 * (Qk_iminus_half_pos_char[4]);
    Qk_iminus_half_neg[4] = r_eigen_Qkm[4].N0 * (Qk_iminus_half_neg_char[0]) + r_eigen_Qkm[4].N1 * (Qk_iminus_half_neg_char[1]) + r_eigen_Qkm[4].N2 * (Qk_iminus_half_neg_char[2]) + r_eigen_Qkm[4].N3 * (Qk_iminus_half_neg_char[3]) + r_eigen_Qkm[4].N4 * (Qk_iminus_half_neg_char[4]);

    /**************************END OF TRANSFER FROM CHARCTERISTIC TO COMPUTATIONAL PLANE**************** */

    F_ip_pos[0] = Qi_iplus_half_pos[1];
    F_ip_pos[1] = (pow(Qi_iplus_half_pos[1], 2) / Qi_iplus_half_pos[0]) + 0.4*(Qi_iplus_half_pos[4] - 0.5*((pow(Qi_iplus_half_pos[1], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[2], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[3], 2) / Qi_iplus_half_pos[0])));
    F_ip_pos[2] = (Qi_iplus_half_pos[1] * Qi_iplus_half_pos[2]) / Qi_iplus_half_pos[0];
    F_ip_pos[3] = (Qi_iplus_half_pos[1] * Qi_iplus_half_pos[3]) / Qi_iplus_half_pos[0];
    F_ip_pos[4] = (1.4*Qi_iplus_half_pos[4] - 0.2*((pow(Qi_iplus_half_pos[1], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[2], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[3], 2) / Qi_iplus_half_pos[0])))*(Qi_iplus_half_pos[1] / Qi_iplus_half_pos[0]);

    F_ip_neg[0] = Qi_iplus_half_neg[1];
    F_ip_neg[1] = (pow(Qi_iplus_half_neg[1], 2) / Qi_iplus_half_neg[0]) + 0.4*(Qi_iplus_half_neg[4] - 0.5*((pow(Qi_iplus_half_neg[1], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[2], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[3], 2) / Qi_iplus_half_neg[0])));
    F_ip_neg[2] = (Qi_iplus_half_neg[1] * Qi_iplus_half_neg[2]) / Qi_iplus_half_neg[0];
    F_ip_neg[3] = (Qi_iplus_half_neg[1] * Qi_iplus_half_neg[3]) / Qi_iplus_half_neg[0];
    F_ip_neg[4] = (1.4*Qi_iplus_half_neg[4] - 0.2*((pow(Qi_iplus_half_neg[1], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[2], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[3], 2) / Qi_iplus_half_neg[0])))*(Qi_iplus_half_neg[1] / Qi_iplus_half_neg[0]);

    F_im_pos[0] = Qi_iminus_half_pos[1];
    F_im_pos[1] = (pow(Qi_iminus_half_pos[1], 2) / Qi_iminus_half_pos[0]) + 0.4*(Qi_iminus_half_pos[4] - 0.5*((pow(Qi_iminus_half_pos[1], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[2], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[3], 2) / Qi_iminus_half_pos[0])));
    F_im_pos[2] = (Qi_iminus_half_pos[1] * Qi_iminus_half_pos[2]) / Qi_iminus_half_pos[0];
    F_im_pos[3] = (Qi_iminus_half_pos[1] * Qi_iminus_half_pos[3]) / Qi_iminus_half_pos[0];
    F_im_pos[4] = (1.4*Qi_iminus_half_pos[4] - 0.2*((pow(Qi_iminus_half_pos[1], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[2], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[3], 2) / Qi_iminus_half_pos[0])))*(Qi_iminus_half_pos[1] / Qi_iminus_half_pos[0]);

    F_im_neg[0] = Qi_iminus_half_neg[1];
    F_im_neg[1] = (pow(Qi_iminus_half_neg[1], 2) / Qi_iminus_half_neg[0]) + 0.4*(Qi_iminus_half_neg[4] - 0.5*((pow(Qi_iminus_half_neg[1], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[2], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[3], 2) / Qi_iminus_half_neg[0])));
    F_im_neg[2] = (Qi_iminus_half_neg[1] * Qi_iminus_half_neg[2]) / Qi_iminus_half_neg[0];
    F_im_neg[3] = (Qi_iminus_half_neg[1] * Qi_iminus_half_neg[3]) / Qi_iminus_half_neg[0];
    F_im_neg[4] = (1.4*Qi_iminus_half_neg[4] - 0.2*((pow(Qi_iminus_half_neg[1], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[2], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[3], 2) / Qi_iminus_half_neg[0])))*(Qi_iminus_half_neg[1] / Qi_iminus_half_neg[0]);

    E_ip_pos[0] = Qi_iplus_half_pos[2];
    E_ip_pos[1] = (Qi_iplus_half_pos[1] * Qi_iplus_half_pos[2]) / Qi_iplus_half_pos[0];
    E_ip_pos[2] = (pow(Qi_iplus_half_pos[2], 2) / Qi_iplus_half_pos[0]) + 0.4*(Qi_iplus_half_pos[4] - 0.5*((pow(Qi_iplus_half_pos[1], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[2], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[3], 2) / Qi_iplus_half_pos[0])));
    E_ip_pos[3] = (Qi_iplus_half_pos[2] * Qi_iplus_half_pos[3]) / Qi_iplus_half_pos[0];
    E_ip_pos[4] = (1.4*Qi_iplus_half_pos[4] - 0.2*((pow(Qi_iplus_half_pos[1], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[2], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[3], 2) / Qi_iplus_half_pos[0])))*(Qi_iplus_half_pos[2] / Qi_iplus_half_pos[0]);

    E_ip_neg[0] = Qi_iplus_half_neg[2];
    E_ip_neg[1] = (Qi_iplus_half_neg[1] * Qi_iplus_half_neg[2]) / Qi_iplus_half_neg[0];
    E_ip_neg[2] = (pow(Qi_iplus_half_neg[2], 2) / Qi_iplus_half_neg[0]) + 0.4*(Qi_iplus_half_neg[4] - 0.5*((pow(Qi_iplus_half_neg[1], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[2], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[3], 2) / Qi_iplus_half_neg[0])));
    E_ip_neg[3] = (Qi_iplus_half_neg[2] * Qi_iplus_half_neg[3]) / Qi_iplus_half_neg[0];
    E_ip_neg[4] = (1.4*Qi_iplus_half_neg[4] - 0.2*((pow(Qi_iplus_half_neg[1], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[2], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[3], 2) / Qi_iplus_half_neg[0])))*(Qi_iplus_half_neg[2] / Qi_iplus_half_neg[0]);

    E_im_pos[0] = Qi_iminus_half_pos[2];
    E_im_pos[1] = (Qi_iminus_half_pos[1] * Qi_iminus_half_pos[2]) / Qi_iminus_half_pos[0];
    E_im_pos[2] = (pow(Qi_iminus_half_pos[2], 2) / Qi_iminus_half_pos[0]) + 0.4*(Qi_iminus_half_pos[4] - 0.5*((pow(Qi_iminus_half_pos[1], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[2], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[3], 2) / Qi_iminus_half_pos[0])));
    E_im_pos[3] = (Qi_iminus_half_pos[2] * Qi_iminus_half_pos[3]) / Qi_iminus_half_pos[0];
    E_im_pos[4] = (1.4*Qi_iminus_half_pos[4] - 0.2*((pow(Qi_iminus_half_pos[1], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[2], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[3], 2) / Qi_iminus_half_pos[0])))*(Qi_iminus_half_pos[2] / Qi_iminus_half_pos[0]);

    E_im_neg[0] = Qi_iminus_half_neg[2];
    E_im_neg[1] = (Qi_iminus_half_neg[1] * Qi_iminus_half_neg[2]) / Qi_iminus_half_neg[0];
    E_im_neg[2] = (pow(Qi_iminus_half_neg[2], 2) / Qi_iminus_half_neg[0]) + 0.4*(Qi_iminus_half_neg[4] - 0.5*((pow(Qi_iminus_half_neg[1], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[2], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[3], 2) / Qi_iminus_half_neg[0])));
    E_im_neg[3] = (Qi_iminus_half_neg[2] * Qi_iminus_half_neg[3]) / Qi_iminus_half_neg[0];
    E_im_neg[4] = (1.4*Qi_iminus_half_neg[4] - 0.2*((pow(Qi_iminus_half_neg[1], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[2], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[3], 2) / Qi_iminus_half_neg[0])))*(Qi_iminus_half_neg[2] / Qi_iminus_half_neg[0]);

    G_ip_pos[0] = Qi_iplus_half_pos[3];
    G_ip_pos[1] = (Qi_iplus_half_pos[1] * Qi_iplus_half_pos[3]) / Qi_iplus_half_pos[0];
    G_ip_pos[2] = (Qi_iplus_half_pos[2] * Qi_iplus_half_pos[3]) / Qi_iplus_half_pos[0];
    G_ip_pos[3] = (pow(Qi_iplus_half_pos[3], 2) / Qi_iplus_half_pos[0]) + 0.4*(Qi_iplus_half_pos[4] - 0.5*((pow(Qi_iplus_half_pos[1], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[2], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[3], 2) / Qi_iplus_half_pos[0])));
    G_ip_pos[4] = (1.4*Qi_iplus_half_pos[4] - 0.2*((pow(Qi_iplus_half_pos[1], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[2], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[3], 2) / Qi_iplus_half_pos[0])))*(Qi_iplus_half_pos[3] / Qi_iplus_half_pos[0]);

    G_ip_neg[0] = Qi_iplus_half_neg[3];
    G_ip_neg[1] = (Qi_iplus_half_neg[1] * Qi_iplus_half_neg[3]) / Qi_iplus_half_neg[0];
    G_ip_neg[2] = (Qi_iplus_half_neg[2] * Qi_iplus_half_neg[3]) / Qi_iplus_half_neg[0];
    G_ip_neg[3] = (pow(Qi_iplus_half_neg[3], 2) / Qi_iplus_half_neg[0]) + 0.4*(Qi_iplus_half_neg[4] - 0.5*((pow(Qi_iplus_half_neg[1], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[2], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[3], 2) / Qi_iplus_half_neg[0])));
    G_ip_neg[4] = (1.4*Qi_iplus_half_neg[4] - 0.2*((pow(Qi_iplus_half_neg[1], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[2], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[3], 2) / Qi_iplus_half_neg[0])))*(Qi_iplus_half_neg[3] / Qi_iplus_half_neg[0]);

    G_im_pos[0] = Qi_iminus_half_pos[3];
    G_im_pos[1] = (Qi_iminus_half_pos[1] * Qi_iminus_half_pos[3]) / Qi_iminus_half_pos[0];
    G_im_pos[2] = (Qi_iminus_half_pos[2] * Qi_iminus_half_pos[3]) / Qi_iminus_half_pos[0];
    G_im_pos[3] = (pow(Qi_iminus_half_pos[3], 2) / Qi_iminus_half_pos[0]) + 0.4*(Qi_iminus_half_pos[4] - 0.5*((pow(Qi_iminus_half_pos[1], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[2], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[3], 2) / Qi_iminus_half_pos[0])));
    G_im_pos[4] = (1.4*Qi_iminus_half_pos[4] - 0.2*((pow(Qi_iminus_half_pos[1], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[2], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[3], 2) / Qi_iminus_half_pos[0])))*(Qi_iminus_half_pos[3] / Qi_iminus_half_pos[0]);

    G_im_neg[0] = Qi_iminus_half_neg[3];
    G_im_neg[1] = (Qi_iminus_half_neg[1] * Qi_iminus_half_neg[3]) / Qi_iminus_half_neg[0];
    G_im_neg[2] = (Qi_iminus_half_neg[2] * Qi_iminus_half_neg[3]) / Qi_iminus_half_neg[0];
    G_im_neg[3] = (pow(Qi_iminus_half_neg[3], 2) / Qi_iminus_half_neg[0]) + 0.4*(Qi_iminus_half_neg[4] - 0.5*((pow(Qi_iminus_half_neg[1], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[2], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[3], 2) / Qi_iminus_half_neg[0])));
    G_im_neg[4] = (1.4*Qi_iminus_half_neg[4] - 0.2*((pow(Qi_iminus_half_neg[1], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[2], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[3], 2) / Qi_iminus_half_neg[0])))*(Qi_iminus_half_neg[3] / Qi_iminus_half_neg[0]);


    F_jp_pos[0] = Qj_iplus_half_pos[1];
    F_jp_pos[1] = (pow(Qj_iplus_half_pos[1], 2) / Qj_iplus_half_pos[0]) + 0.4*(Qj_iplus_half_pos[4] - 0.5*((pow(Qj_iplus_half_pos[1], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[2], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[3], 2) / Qj_iplus_half_pos[0])));
    F_jp_pos[2] = (Qj_iplus_half_pos[1] * Qj_iplus_half_pos[2]) / Qj_iplus_half_pos[0];
    F_jp_pos[3] = (Qj_iplus_half_pos[1] * Qj_iplus_half_pos[3]) / Qj_iplus_half_pos[0];
    F_jp_pos[4] = (1.4*Qj_iplus_half_pos[4] - 0.2*((pow(Qj_iplus_half_pos[1], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[2], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[3], 2) / Qj_iplus_half_pos[0])))*(Qj_iplus_half_pos[1] / Qj_iplus_half_pos[0]);

    F_jp_neg[0] = Qj_iplus_half_neg[1];
    F_jp_neg[1] = (pow(Qj_iplus_half_neg[1], 2) / Qj_iplus_half_neg[0]) + 0.4*(Qj_iplus_half_neg[4] - 0.5*((pow(Qj_iplus_half_neg[1], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[2], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[3], 2) / Qj_iplus_half_neg[0])));
    F_jp_neg[2] = (Qj_iplus_half_neg[1] * Qj_iplus_half_neg[2]) / Qj_iplus_half_neg[0];
    F_jp_neg[3] = (Qj_iplus_half_neg[1] * Qj_iplus_half_neg[3]) / Qj_iplus_half_neg[0];
    F_jp_neg[4] = (1.4*Qj_iplus_half_neg[4] - 0.2*((pow(Qj_iplus_half_neg[1], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[2], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[3], 2) / Qj_iplus_half_neg[0])))*(Qj_iplus_half_neg[1] / Qj_iplus_half_neg[0]);

    F_jm_pos[0] = Qj_iminus_half_pos[1];
    F_jm_pos[1] = (pow(Qj_iminus_half_pos[1], 2) / Qj_iminus_half_pos[0]) + 0.4*(Qj_iminus_half_pos[4] - 0.5*((pow(Qj_iminus_half_pos[1], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[2], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[3], 2) / Qj_iminus_half_pos[0])));
    F_jm_pos[2] = (Qj_iminus_half_pos[1] * Qj_iminus_half_pos[2]) / Qj_iminus_half_pos[0];
    F_jm_pos[3] = (Qj_iminus_half_pos[1] * Qj_iminus_half_pos[3]) / Qj_iminus_half_pos[0];
    F_jm_pos[4] = (1.4*Qj_iminus_half_pos[4] - 0.2*((pow(Qj_iminus_half_pos[1], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[2], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[3], 2) / Qj_iminus_half_pos[0])))*(Qj_iminus_half_pos[1] / Qj_iminus_half_pos[0]);

    F_jm_neg[0] = Qj_iminus_half_neg[1];
    F_jm_neg[1] = (pow(Qj_iminus_half_neg[1], 2) / Qj_iminus_half_neg[0]) + 0.4*(Qj_iminus_half_neg[4] - 0.5*((pow(Qj_iminus_half_neg[1], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[2], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[3], 2) / Qj_iminus_half_neg[0])));
    F_jm_neg[2] = (Qj_iminus_half_neg[1] * Qj_iminus_half_neg[2]) / Qj_iminus_half_neg[0];
    F_jm_neg[3] = (Qj_iminus_half_neg[1] * Qj_iminus_half_neg[3]) / Qj_iminus_half_neg[0];
    F_jm_neg[4] = (1.4*Qj_iminus_half_neg[4] - 0.2*((pow(Qj_iminus_half_neg[1], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[2], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[3], 2) / Qj_iminus_half_neg[0])))*(Qj_iminus_half_neg[1] / Qj_iminus_half_neg[0]);

    E_jp_pos[0] = Qj_iplus_half_pos[2];
    E_jp_pos[1] = (Qj_iplus_half_pos[1] * Qj_iplus_half_pos[2]) / Qj_iplus_half_pos[0];
    E_jp_pos[2] = (pow(Qj_iplus_half_pos[2], 2) / Qj_iplus_half_pos[0]) + 0.4*(Qj_iplus_half_pos[4] - 0.5*((pow(Qj_iplus_half_pos[1], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[2], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[3], 2) / Qj_iplus_half_pos[0])));
    E_jp_pos[3] = (Qj_iplus_half_pos[2] * Qj_iplus_half_pos[3]) / Qj_iplus_half_pos[0];
    E_jp_pos[4] = (1.4*Qj_iplus_half_pos[4] - 0.2*((pow(Qj_iplus_half_pos[1], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[2], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[3], 2) / Qj_iplus_half_pos[0])))*(Qj_iplus_half_pos[2] / Qj_iplus_half_pos[0]);

    E_jp_neg[0] = Qj_iplus_half_neg[2];
    E_jp_neg[1] = (Qj_iplus_half_neg[1] * Qj_iplus_half_neg[2]) / Qj_iplus_half_neg[0];
    E_jp_neg[2] = (pow(Qj_iplus_half_neg[2], 2) / Qj_iplus_half_neg[0]) + 0.4*(Qj_iplus_half_neg[4] - 0.5*((pow(Qj_iplus_half_neg[1], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[2], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[3], 2) / Qj_iplus_half_neg[0])));
    E_jp_neg[3] = (Qj_iplus_half_neg[2] * Qj_iplus_half_neg[3]) / Qj_iplus_half_neg[0];
    E_jp_neg[4] = (1.4*Qj_iplus_half_neg[4] - 0.2*((pow(Qj_iplus_half_neg[1], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[2], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[3], 2) / Qj_iplus_half_neg[0])))*(Qj_iplus_half_neg[2] / Qj_iplus_half_neg[0]);

    E_jm_pos[0] = Qj_iminus_half_pos[2];
    E_jm_pos[1] = (Qj_iminus_half_pos[1] * Qj_iminus_half_pos[2]) / Qj_iminus_half_pos[0];
    E_jm_pos[2] = (pow(Qj_iminus_half_pos[2], 2) / Qj_iminus_half_pos[0]) + 0.4*(Qj_iminus_half_pos[4] - 0.5*((pow(Qj_iminus_half_pos[1], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[2], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[3], 2) / Qj_iminus_half_pos[0])));
    E_jm_pos[3] = (Qj_iminus_half_pos[2] * Qj_iminus_half_pos[3]) / Qj_iminus_half_pos[0];
    E_jm_pos[4] = (1.4*Qj_iminus_half_pos[4] - 0.2*((pow(Qj_iminus_half_pos[1], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[2], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[3], 2) / Qj_iminus_half_pos[0])))*(Qj_iminus_half_pos[2] / Qj_iminus_half_pos[0]);

    E_jm_neg[0] = Qj_iminus_half_neg[2];
    E_jm_neg[1] = (Qj_iminus_half_neg[1] * Qj_iminus_half_neg[2]) / Qj_iminus_half_neg[0];
    E_jm_neg[2] = (pow(Qj_iminus_half_neg[2], 2) / Qj_iminus_half_neg[0]) + 0.4*(Qj_iminus_half_neg[4] - 0.5*((pow(Qj_iminus_half_neg[1], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[2], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[3], 2) / Qj_iminus_half_neg[0])));
    E_jm_neg[3] = (Qj_iminus_half_neg[2] * Qj_iminus_half_neg[3]) / Qj_iminus_half_neg[0];
    E_jm_neg[4] = (1.4*Qj_iminus_half_neg[4] - 0.2*((pow(Qj_iminus_half_neg[1], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[2], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[3], 2) / Qj_iminus_half_neg[0])))*(Qj_iminus_half_neg[2] / Qj_iminus_half_neg[0]);

    G_jp_pos[0] = Qj_iplus_half_pos[3];
    G_jp_pos[1] = (Qj_iplus_half_pos[1] * Qj_iplus_half_pos[3]) / Qj_iplus_half_pos[0];
    G_jp_pos[2] = (Qj_iplus_half_pos[2] * Qj_iplus_half_pos[3]) / Qj_iplus_half_pos[0];
    G_jp_pos[3] = (pow(Qj_iplus_half_pos[3], 2) / Qj_iplus_half_pos[0]) + 0.4*(Qj_iplus_half_pos[4] - 0.5*((pow(Qj_iplus_half_pos[1], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[2], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[3], 2) / Qj_iplus_half_pos[0])));
    G_jp_pos[4] = (1.4*Qj_iplus_half_pos[4] - 0.2*((pow(Qj_iplus_half_pos[1], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[2], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[3], 2) / Qj_iplus_half_pos[0])))*(Qj_iplus_half_pos[3] / Qj_iplus_half_pos[0]);

    G_jp_neg[0] = Qj_iplus_half_neg[3];
    G_jp_neg[1] = (Qj_iplus_half_neg[1] * Qj_iplus_half_neg[3]) / Qj_iplus_half_neg[0];
    G_jp_neg[2] = (Qj_iplus_half_neg[2] * Qj_iplus_half_neg[3]) / Qj_iplus_half_neg[0];
    G_jp_neg[3] = (pow(Qj_iplus_half_neg[3], 2) / Qj_iplus_half_neg[0]) + 0.4*(Qj_iplus_half_neg[4] - 0.5*((pow(Qj_iplus_half_neg[1], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[2], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[3], 2) / Qj_iplus_half_neg[0])));
    G_jp_neg[4] = (1.4*Qj_iplus_half_neg[4] - 0.2*((pow(Qj_iplus_half_neg[1], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[2], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[3], 2) / Qj_iplus_half_neg[0])))*(Qj_iplus_half_neg[3] / Qj_iplus_half_neg[0]);

    G_jm_pos[0] = Qj_iminus_half_pos[3];
    G_jm_pos[1] = (Qj_iminus_half_pos[1] * Qj_iminus_half_pos[3]) / Qj_iminus_half_pos[0];
    G_jm_pos[2] = (Qj_iminus_half_pos[2] * Qj_iminus_half_pos[3]) / Qj_iminus_half_pos[0];
    G_jm_pos[3] = (pow(Qj_iminus_half_pos[3], 2) / Qj_iminus_half_pos[0]) + 0.4*(Qj_iminus_half_pos[4] - 0.5*((pow(Qj_iminus_half_pos[1], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[2], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[3], 2) / Qj_iminus_half_pos[0])));
    G_jm_pos[4] = (1.4*Qj_iminus_half_pos[4] - 0.2*((pow(Qj_iminus_half_pos[1], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[2], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[3], 2) / Qj_iminus_half_pos[0])))*(Qj_iminus_half_pos[3] / Qj_iminus_half_pos[0]);

    G_jm_neg[0] = Qj_iminus_half_neg[3];
    G_jm_neg[1] = (Qj_iminus_half_neg[1] * Qj_iminus_half_neg[3]) / Qj_iminus_half_neg[0];
    G_jm_neg[2] = (Qj_iminus_half_neg[2] * Qj_iminus_half_neg[3]) / Qj_iminus_half_neg[0];
    G_jm_neg[3] = (pow(Qj_iminus_half_neg[3], 2) / Qj_iminus_half_neg[0]) + 0.4*(Qj_iminus_half_neg[4] - 0.5*((pow(Qj_iminus_half_neg[1], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[2], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[3], 2) / Qj_iminus_half_neg[0])));
    G_jm_neg[4] = (1.4*Qj_iminus_half_neg[4] - 0.2*((pow(Qj_iminus_half_neg[1], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[2], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[3], 2) / Qj_iminus_half_neg[0])))*(Qj_iminus_half_neg[3] / Qj_iminus_half_neg[0]);


    F_kp_pos[0] = Qk_iplus_half_pos[1];
    F_kp_pos[1] = (pow(Qk_iplus_half_pos[1], 2) / Qk_iplus_half_pos[0]) + 0.4*(Qk_iplus_half_pos[4] - 0.5*((pow(Qk_iplus_half_pos[1], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[2], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[3], 2) / Qk_iplus_half_pos[0])));
    F_kp_pos[2] = (Qk_iplus_half_pos[1] * Qk_iplus_half_pos[2]) / Qk_iplus_half_pos[0];
    F_kp_pos[3] = (Qk_iplus_half_pos[1] * Qk_iplus_half_pos[3]) / Qk_iplus_half_pos[0];
    F_kp_pos[4] = (1.4*Qk_iplus_half_pos[4] - 0.2*((pow(Qk_iplus_half_pos[1], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[2], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[3], 2) / Qk_iplus_half_pos[0])))*(Qk_iplus_half_pos[1] / Qk_iplus_half_pos[0]);

    F_kp_neg[0] = Qk_iplus_half_neg[1];
    F_kp_neg[1] = (pow(Qk_iplus_half_neg[1], 2) / Qk_iplus_half_neg[0]) + 0.4*(Qk_iplus_half_neg[4] - 0.5*((pow(Qk_iplus_half_neg[1], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[2], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[3], 2) / Qk_iplus_half_neg[0])));
    F_kp_neg[2] = (Qk_iplus_half_neg[1] * Qk_iplus_half_neg[2]) / Qk_iplus_half_neg[0];
    F_kp_neg[3] = (Qk_iplus_half_neg[1] * Qk_iplus_half_neg[3]) / Qk_iplus_half_neg[0];
    F_kp_neg[4] = (1.4*Qk_iplus_half_neg[4] - 0.2*((pow(Qk_iplus_half_neg[1], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[2], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[3], 2) / Qk_iplus_half_neg[0])))*(Qk_iplus_half_neg[1] / Qk_iplus_half_neg[0]);

    F_km_pos[0] = Qk_iminus_half_pos[1];
    F_km_pos[1] = (pow(Qk_iminus_half_pos[1], 2) / Qk_iminus_half_pos[0]) + 0.4*(Qk_iminus_half_pos[4] - 0.5*((pow(Qk_iminus_half_pos[1], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[2], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[3], 2) / Qk_iminus_half_pos[0])));
    F_km_pos[2] = (Qk_iminus_half_pos[1] * Qk_iminus_half_pos[2]) / Qk_iminus_half_pos[0];
    F_km_pos[3] = (Qk_iminus_half_pos[1] * Qk_iminus_half_pos[3]) / Qk_iminus_half_pos[0];
    F_km_pos[4] = (1.4*Qk_iminus_half_pos[4] - 0.2*((pow(Qk_iminus_half_pos[1], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[2], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[3], 2) / Qk_iminus_half_pos[0])))*(Qk_iminus_half_pos[1] / Qk_iminus_half_pos[0]);

    F_km_neg[0] = Qk_iminus_half_neg[1];
    F_km_neg[1] = (pow(Qk_iminus_half_neg[1], 2) / Qk_iminus_half_neg[0]) + 0.4*(Qk_iminus_half_neg[4] - 0.5*((pow(Qk_iminus_half_neg[1], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[2], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[3], 2) / Qk_iminus_half_neg[0])));
    F_km_neg[2] = (Qk_iminus_half_neg[1] * Qk_iminus_half_neg[2]) / Qk_iminus_half_neg[0];
    F_km_neg[3] = (Qk_iminus_half_neg[1] * Qk_iminus_half_neg[3]) / Qk_iminus_half_neg[0];
    F_km_neg[4] = (1.4*Qk_iminus_half_neg[4] - 0.2*((pow(Qk_iminus_half_neg[1], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[2], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[3], 2) / Qk_iminus_half_neg[0])))*(Qk_iminus_half_neg[1] / Qk_iminus_half_neg[0]);

    E_kp_pos[0] = Qk_iplus_half_pos[2];
    E_kp_pos[1] = (Qk_iplus_half_pos[1] * Qk_iplus_half_pos[2]) / Qk_iplus_half_pos[0];
    E_kp_pos[2] = (pow(Qk_iplus_half_pos[2], 2) / Qk_iplus_half_pos[0]) + 0.4*(Qk_iplus_half_pos[4] - 0.5*((pow(Qk_iplus_half_pos[1], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[2], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[3], 2) / Qk_iplus_half_pos[0])));
    E_kp_pos[3] = (Qk_iplus_half_pos[2] * Qk_iplus_half_pos[3]) / Qk_iplus_half_pos[0];
    E_kp_pos[4] = (1.4*Qk_iplus_half_pos[4] - 0.2*((pow(Qk_iplus_half_pos[1], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[2], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[3], 2) / Qk_iplus_half_pos[0])))*(Qk_iplus_half_pos[2] / Qk_iplus_half_pos[0]);

    E_kp_neg[0] = Qk_iplus_half_neg[2];
    E_kp_neg[1] = (Qk_iplus_half_neg[1] * Qk_iplus_half_neg[2]) / Qk_iplus_half_neg[0];
    E_kp_neg[2] = (pow(Qk_iplus_half_neg[2], 2) / Qk_iplus_half_neg[0]) + 0.4*(Qk_iplus_half_neg[4] - 0.5*((pow(Qk_iplus_half_neg[1], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[2], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[3], 2) / Qk_iplus_half_neg[0])));
    E_kp_neg[3] = (Qk_iplus_half_neg[2] * Qk_iplus_half_neg[3]) / Qk_iplus_half_neg[0];
    E_kp_neg[4] = (1.4*Qk_iplus_half_neg[4] - 0.2*((pow(Qk_iplus_half_neg[1], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[2], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[3], 2) / Qk_iplus_half_neg[0])))*(Qk_iplus_half_neg[2] / Qk_iplus_half_neg[0]);

    E_km_pos[0] = Qk_iminus_half_pos[2];
    E_km_pos[1] = (Qk_iminus_half_pos[1] * Qk_iminus_half_pos[2]) / Qk_iminus_half_pos[0];
    E_km_pos[2] = (pow(Qk_iminus_half_pos[2], 2) / Qk_iminus_half_pos[0]) + 0.4*(Qk_iminus_half_pos[4] - 0.5*((pow(Qk_iminus_half_pos[1], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[2], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[3], 2) / Qk_iminus_half_pos[0])));
    E_km_pos[3] = (Qk_iminus_half_pos[2] * Qk_iminus_half_pos[3]) / Qk_iminus_half_pos[0];
    E_km_pos[4] = (1.4*Qk_iminus_half_pos[4] - 0.2*((pow(Qk_iminus_half_pos[1], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[2], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[3], 2) / Qk_iminus_half_pos[0])))*(Qk_iminus_half_pos[2] / Qk_iminus_half_pos[0]);

    E_km_neg[0] = Qk_iminus_half_neg[2];
    E_km_neg[1] = (Qk_iminus_half_neg[1] * Qk_iminus_half_neg[2]) / Qk_iminus_half_neg[0];
    E_km_neg[2] = (pow(Qk_iminus_half_neg[2], 2) / Qk_iminus_half_neg[0]) + 0.4*(Qk_iminus_half_neg[4] - 0.5*((pow(Qk_iminus_half_neg[1], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[2], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[3], 2) / Qk_iminus_half_neg[0])));
    E_km_neg[3] = (Qk_iminus_half_neg[2] * Qk_iminus_half_neg[3]) / Qk_iminus_half_neg[0];
    E_km_neg[4] = (1.4*Qk_iminus_half_neg[4] - 0.2*((pow(Qk_iminus_half_neg[1], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[2], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[3], 2) / Qk_iminus_half_neg[0])))*(Qk_iminus_half_neg[2] / Qk_iminus_half_neg[0]);

    G_kp_pos[0] = Qk_iplus_half_pos[3];
    G_kp_pos[1] = (Qk_iplus_half_pos[1] * Qk_iplus_half_pos[3]) / Qk_iplus_half_pos[0];
    G_kp_pos[2] = (Qk_iplus_half_pos[2] * Qk_iplus_half_pos[3]) / Qk_iplus_half_pos[0];
    G_kp_pos[3] = (pow(Qk_iplus_half_pos[3], 2) / Qk_iplus_half_pos[0]) + 0.4*(Qk_iplus_half_pos[4] - 0.5*((pow(Qk_iplus_half_pos[1], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[2], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[3], 2) / Qk_iplus_half_pos[0])));
    G_kp_pos[4] = (1.4*Qk_iplus_half_pos[4] - 0.2*((pow(Qk_iplus_half_pos[1], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[2], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[3], 2) / Qk_iplus_half_pos[0])))*(Qk_iplus_half_pos[3] / Qk_iplus_half_pos[0]);

    G_kp_neg[0] = Qk_iplus_half_neg[3];
    G_kp_neg[1] = (Qk_iplus_half_neg[1] * Qk_iplus_half_neg[3]) / Qk_iplus_half_neg[0];
    G_kp_neg[2] = (Qk_iplus_half_neg[2] * Qk_iplus_half_neg[3]) / Qk_iplus_half_neg[0];
    G_kp_neg[3] = (pow(Qk_iplus_half_neg[3], 2) / Qk_iplus_half_neg[0]) + 0.4*(Qk_iplus_half_neg[4] - 0.5*((pow(Qk_iplus_half_neg[1], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[2], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[3], 2) / Qk_iplus_half_neg[0])));
    G_kp_neg[4] = (1.4*Qk_iplus_half_neg[4] - 0.2*((pow(Qk_iplus_half_neg[1], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[2], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[3], 2) / Qk_iplus_half_neg[0])))*(Qk_iplus_half_neg[3] / Qk_iplus_half_neg[0]);

    G_km_pos[0] = Qk_iminus_half_pos[3];
    G_km_pos[1] = (Qk_iminus_half_pos[1] * Qk_iminus_half_pos[3]) / Qk_iminus_half_pos[0];
    G_km_pos[2] = (Qk_iminus_half_pos[2] * Qk_iminus_half_pos[3]) / Qk_iminus_half_pos[0];
    G_km_pos[3] = (pow(Qk_iminus_half_pos[3], 2) / Qk_iminus_half_pos[0]) + 0.4*(Qk_iminus_half_pos[4] - 0.5*((pow(Qk_iminus_half_pos[1], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[2], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[3], 2) / Qk_iminus_half_pos[0])));
    G_km_pos[4] = (1.4*Qk_iminus_half_pos[4] - 0.2*((pow(Qk_iminus_half_pos[1], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[2], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[3], 2) / Qk_iminus_half_pos[0])))*(Qk_iminus_half_pos[3] / Qk_iminus_half_pos[0]);

    G_km_neg[0] = Qk_iminus_half_neg[3];
    G_km_neg[1] = (Qk_iminus_half_neg[1] * Qk_iminus_half_neg[3]) / Qk_iminus_half_neg[0];
    G_km_neg[2] = (Qk_iminus_half_neg[2] * Qk_iminus_half_neg[3]) / Qk_iminus_half_neg[0];
    G_km_neg[3] = (pow(Qk_iminus_half_neg[3], 2) / Qk_iminus_half_neg[0]) + 0.4*(Qk_iminus_half_neg[4] - 0.5*((pow(Qk_iminus_half_neg[1], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[2], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[3], 2) / Qk_iminus_half_neg[0])));
    G_km_neg[4] = (1.4*Qk_iminus_half_neg[4] - 0.2*((pow(Qk_iminus_half_neg[1], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[2], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[3], 2) / Qk_iminus_half_neg[0])))*(Qk_iminus_half_neg[3] / Qk_iminus_half_neg[0]);

    /******************************************************************************************************************************/

    eigen_Qip[0] = (metric[Gid].zeta_xip*(Qi_iplus_half_pos[1] / Qi_iplus_half_pos[0]) + metric[Gid].zeta_yip * (Qi_iplus_half_pos[2] / Qi_iplus_half_pos[0]) + metric[Gid].zeta_zip * (Qi_iplus_half_pos[3] / Qi_iplus_half_pos[0])) - sqrt(1.4*0.4*((Qi_iplus_half_pos[4] / Qi_iplus_half_pos[0]) - ((pow(Qi_iplus_half_pos[1], 2) + pow(Qi_iplus_half_pos[2], 2) + pow(Qi_iplus_half_pos[3], 2)) / (2.0*pow(Qi_iplus_half_pos[0], 2))))*(metric[Gid].zeta_xip*metric[Gid].zeta_xip + metric[Gid].zeta_yip * metric[Gid].zeta_yip + metric[Gid].zeta_zip * metric[Gid].zeta_zip));
    eigen_Qip[1] = (metric[Gid].zeta_xip*(Qi_iplus_half_pos[1] / Qi_iplus_half_pos[0]) + metric[Gid].zeta_yip * (Qi_iplus_half_pos[2] / Qi_iplus_half_pos[0]) + metric[Gid].zeta_zip * (Qi_iplus_half_pos[3] / Qi_iplus_half_pos[0]));
    eigen_Qip[2] = (metric[Gid].zeta_xip*(Qi_iplus_half_pos[1] / Qi_iplus_half_pos[0]) + metric[Gid].zeta_yip * (Qi_iplus_half_pos[2] / Qi_iplus_half_pos[0]) + metric[Gid].zeta_zip * (Qi_iplus_half_pos[3] / Qi_iplus_half_pos[0]));
    eigen_Qip[3] = (metric[Gid].zeta_xip*(Qi_iplus_half_pos[1] / Qi_iplus_half_pos[0]) + metric[Gid].zeta_yip * (Qi_iplus_half_pos[2] / Qi_iplus_half_pos[0]) + metric[Gid].zeta_zip * (Qi_iplus_half_pos[3] / Qi_iplus_half_pos[0]));
    eigen_Qip[4] = (metric[Gid].zeta_xip*(Qi_iplus_half_pos[1] / Qi_iplus_half_pos[0]) + metric[Gid].zeta_yip * (Qi_iplus_half_pos[2] / Qi_iplus_half_pos[0]) + metric[Gid].zeta_zip * (Qi_iplus_half_pos[3] / Qi_iplus_half_pos[0])) + sqrt(1.4*0.4*((Qi_iplus_half_pos[4] / Qi_iplus_half_pos[0]) - ((pow(Qi_iplus_half_pos[1], 2) + pow(Qi_iplus_half_pos[2], 2) + pow(Qi_iplus_half_pos[3], 2)) / (2.0*pow(Qi_iplus_half_pos[0], 2))))*(metric[Gid].zeta_xip*metric[Gid].zeta_xip + metric[Gid].zeta_yip * metric[Gid].zeta_yip + metric[Gid].zeta_zip * metric[Gid].zeta_zip));

    eigen_Qim[0] = (metric[Gid].zeta_xip*(Qi_iplus_half_neg[1] / Qi_iplus_half_neg[0]) + metric[Gid].zeta_yip * (Qi_iplus_half_neg[2] / Qi_iplus_half_neg[0]) + metric[Gid].zeta_zip * (Qi_iplus_half_neg[3] / Qi_iplus_half_neg[0])) - sqrt(1.4*0.4*((Qi_iplus_half_neg[4] / Qi_iplus_half_neg[0]) - ((pow(Qi_iplus_half_neg[1], 2) + pow(Qi_iplus_half_neg[2], 2) + pow(Qi_iplus_half_neg[3], 2)) / (2.0*pow(Qi_iplus_half_neg[0], 2))))*(metric[Gid].zeta_xip*metric[Gid].zeta_xip + metric[Gid].zeta_yip * metric[Gid].zeta_yip + metric[Gid].zeta_zip * metric[Gid].zeta_zip));
    eigen_Qim[1] = (metric[Gid].zeta_xip*(Qi_iplus_half_neg[1] / Qi_iplus_half_neg[0]) + metric[Gid].zeta_yip * (Qi_iplus_half_neg[2] / Qi_iplus_half_neg[0]) + metric[Gid].zeta_zip * (Qi_iplus_half_neg[3] / Qi_iplus_half_neg[0]));
    eigen_Qim[2] = (metric[Gid].zeta_xip*(Qi_iplus_half_neg[1] / Qi_iplus_half_neg[0]) + metric[Gid].zeta_yip * (Qi_iplus_half_neg[2] / Qi_iplus_half_neg[0]) + metric[Gid].zeta_zip * (Qi_iplus_half_neg[3] / Qi_iplus_half_neg[0]));
    eigen_Qim[3] = (metric[Gid].zeta_xip*(Qi_iplus_half_neg[1] / Qi_iplus_half_neg[0]) + metric[Gid].zeta_yip * (Qi_iplus_half_neg[2] / Qi_iplus_half_neg[0]) + metric[Gid].zeta_zip * (Qi_iplus_half_neg[3] / Qi_iplus_half_neg[0]));
    eigen_Qim[4] = (metric[Gid].zeta_xip*(Qi_iplus_half_neg[1] / Qi_iplus_half_neg[0]) + metric[Gid].zeta_yip * (Qi_iplus_half_neg[2] / Qi_iplus_half_neg[0]) + metric[Gid].zeta_zip * (Qi_iplus_half_neg[3] / Qi_iplus_half_neg[0])) + sqrt(1.4*0.4*((Qi_iplus_half_neg[4] / Qi_iplus_half_neg[0]) - ((pow(Qi_iplus_half_neg[1], 2) + pow(Qi_iplus_half_neg[2], 2) + pow(Qi_iplus_half_neg[3], 2)) / (2.0*pow(Qi_iplus_half_neg[0], 2))))*(metric[Gid].zeta_xip*metric[Gid].zeta_xip + metric[Gid].zeta_yip * metric[Gid].zeta_yip + metric[Gid].zeta_zip * metric[Gid].zeta_zip));

    eigen_Qinp[0] = (metric[Gid].zeta_xim*(Qi_iminus_half_pos[1] / Qi_iminus_half_pos[0]) + metric[Gid].zeta_yim * (Qi_iminus_half_pos[2] / Qi_iminus_half_pos[0]) + metric[Gid].zeta_zim * (Qi_iminus_half_pos[3] / Qi_iminus_half_pos[0])) - sqrt(1.4*0.4*((Qi_iminus_half_pos[4] / Qi_iminus_half_pos[0]) - ((pow(Qi_iminus_half_pos[1], 2) + pow(Qi_iminus_half_pos[2], 2) + pow(Qi_iminus_half_pos[3], 2)) / (2.0*pow(Qi_iminus_half_pos[0], 2))))*(metric[Gid].zeta_xim*metric[Gid].zeta_xim + metric[Gid].zeta_yim * metric[Gid].zeta_yim + metric[Gid].zeta_zim * metric[Gid].zeta_zim));
    eigen_Qinp[1] = (metric[Gid].zeta_xim*(Qi_iminus_half_pos[1] / Qi_iminus_half_pos[0]) + metric[Gid].zeta_yim * (Qi_iminus_half_pos[2] / Qi_iminus_half_pos[0]) + metric[Gid].zeta_zim * (Qi_iminus_half_pos[3] / Qi_iminus_half_pos[0]));
    eigen_Qinp[2] = (metric[Gid].zeta_xim*(Qi_iminus_half_pos[1] / Qi_iminus_half_pos[0]) + metric[Gid].zeta_yim * (Qi_iminus_half_pos[2] / Qi_iminus_half_pos[0]) + metric[Gid].zeta_zim * (Qi_iminus_half_pos[3] / Qi_iminus_half_pos[0]));
    eigen_Qinp[3] = (metric[Gid].zeta_xim*(Qi_iminus_half_pos[1] / Qi_iminus_half_pos[0]) + metric[Gid].zeta_yim * (Qi_iminus_half_pos[2] / Qi_iminus_half_pos[0]) + metric[Gid].zeta_zim * (Qi_iminus_half_pos[3] / Qi_iminus_half_pos[0]));
    eigen_Qinp[4] = (metric[Gid].zeta_xim*(Qi_iminus_half_pos[1] / Qi_iminus_half_pos[0]) + metric[Gid].zeta_yim * (Qi_iminus_half_pos[2] / Qi_iminus_half_pos[0]) + metric[Gid].zeta_zim * (Qi_iminus_half_pos[3] / Qi_iminus_half_pos[0])) + sqrt(1.4*0.4*((Qi_iminus_half_pos[4] / Qi_iminus_half_pos[0]) - ((pow(Qi_iminus_half_pos[1], 2) + pow(Qi_iminus_half_pos[2], 2) + pow(Qi_iminus_half_pos[3], 2)) / (2.0*pow(Qi_iminus_half_pos[0], 2))))*(metric[Gid].zeta_xim*metric[Gid].zeta_xim + metric[Gid].zeta_yim * metric[Gid].zeta_yim + metric[Gid].zeta_zim * metric[Gid].zeta_zim));

    eigen_Qinm[0] = (metric[Gid].zeta_xim*(Qi_iminus_half_neg[1] / Qi_iminus_half_neg[0]) + metric[Gid].zeta_yim * (Qi_iminus_half_neg[2] / Qi_iminus_half_neg[0]) + metric[Gid].zeta_zim * (Qi_iminus_half_neg[3] / Qi_iminus_half_neg[0])) - sqrt(1.4*0.4*((Qi_iminus_half_neg[4] / Qi_iminus_half_neg[0]) - ((pow(Qi_iminus_half_neg[1], 2) + pow(Qi_iminus_half_neg[2], 2) + pow(Qi_iminus_half_neg[3], 2)) / (2.0*pow(Qi_iminus_half_neg[0], 2))))*(metric[Gid].zeta_xim*metric[Gid].zeta_xim + metric[Gid].zeta_yim * metric[Gid].zeta_yim + metric[Gid].zeta_zim * metric[Gid].zeta_zim));
    eigen_Qinm[1] = (metric[Gid].zeta_xim*(Qi_iminus_half_neg[1] / Qi_iminus_half_neg[0]) + metric[Gid].zeta_yim * (Qi_iminus_half_neg[2] / Qi_iminus_half_neg[0]) + metric[Gid].zeta_zim * (Qi_iminus_half_neg[3] / Qi_iminus_half_neg[0]));
    eigen_Qinm[2] = (metric[Gid].zeta_xim*(Qi_iminus_half_neg[1] / Qi_iminus_half_neg[0]) + metric[Gid].zeta_yim * (Qi_iminus_half_neg[2] / Qi_iminus_half_neg[0]) + metric[Gid].zeta_zim * (Qi_iminus_half_neg[3] / Qi_iminus_half_neg[0]));
    eigen_Qinm[3] = (metric[Gid].zeta_xim*(Qi_iminus_half_neg[1] / Qi_iminus_half_neg[0]) + metric[Gid].zeta_yim * (Qi_iminus_half_neg[2] / Qi_iminus_half_neg[0]) + metric[Gid].zeta_zim * (Qi_iminus_half_neg[3] / Qi_iminus_half_neg[0]));
    eigen_Qinm[4] = (metric[Gid].zeta_xim*(Qi_iminus_half_neg[1] / Qi_iminus_half_neg[0]) + metric[Gid].zeta_yim * (Qi_iminus_half_neg[2] / Qi_iminus_half_neg[0]) + metric[Gid].zeta_zim * (Qi_iminus_half_neg[3] / Qi_iminus_half_neg[0])) + sqrt(1.4*0.4*((Qi_iminus_half_neg[4] / Qi_iminus_half_neg[0]) - ((pow(Qi_iminus_half_neg[1], 2) + pow(Qi_iminus_half_neg[2], 2) + pow(Qi_iminus_half_neg[3], 2)) / (2.0*pow(Qi_iminus_half_neg[0], 2))))*(metric[Gid].zeta_xim*metric[Gid].zeta_xim + metric[Gid].zeta_yim * metric[Gid].zeta_yim + metric[Gid].zeta_zim * metric[Gid].zeta_zim));

    eigen_Qjp[0] = (metric[Gid].eta_xjp*(Qj_iplus_half_pos[1] / Qj_iplus_half_pos[0]) + metric[Gid].eta_yjp * (Qj_iplus_half_pos[2] / Qj_iplus_half_pos[0]) + metric[Gid].eta_zjp * (Qj_iplus_half_pos[3] / Qj_iplus_half_pos[0])) - sqrt(1.4*0.4*((Qj_iplus_half_pos[4] / Qj_iplus_half_pos[0]) - ((pow(Qj_iplus_half_pos[1], 2) + pow(Qj_iplus_half_pos[2], 2) + pow(Qj_iplus_half_pos[3], 2)) / (2.0*pow(Qj_iplus_half_pos[0], 2))))*(metric[Gid].eta_xjp*metric[Gid].eta_xjp + metric[Gid].eta_yjp * metric[Gid].eta_yjp + metric[Gid].eta_zjp * metric[Gid].eta_zjp));
    eigen_Qjp[1] = (metric[Gid].eta_xjp*(Qj_iplus_half_pos[1] / Qj_iplus_half_pos[0]) + metric[Gid].eta_yjp * (Qj_iplus_half_pos[2] / Qj_iplus_half_pos[0]) + metric[Gid].eta_zjp * (Qj_iplus_half_pos[3] / Qj_iplus_half_pos[0]));
    eigen_Qjp[2] = (metric[Gid].eta_xjp*(Qj_iplus_half_pos[1] / Qj_iplus_half_pos[0]) + metric[Gid].eta_yjp * (Qj_iplus_half_pos[2] / Qj_iplus_half_pos[0]) + metric[Gid].eta_zjp * (Qj_iplus_half_pos[3] / Qj_iplus_half_pos[0]));
    eigen_Qjp[3] = (metric[Gid].eta_xjp*(Qj_iplus_half_pos[1] / Qj_iplus_half_pos[0]) + metric[Gid].eta_yjp * (Qj_iplus_half_pos[2] / Qj_iplus_half_pos[0]) + metric[Gid].eta_zjp * (Qj_iplus_half_pos[3] / Qj_iplus_half_pos[0]));
    eigen_Qjp[4] = (metric[Gid].eta_xjp*(Qj_iplus_half_pos[1] / Qj_iplus_half_pos[0]) + metric[Gid].eta_yjp * (Qj_iplus_half_pos[2] / Qj_iplus_half_pos[0]) + metric[Gid].eta_zjp * (Qj_iplus_half_pos[3] / Qj_iplus_half_pos[0])) + sqrt(1.4*0.4*((Qj_iplus_half_pos[4] / Qj_iplus_half_pos[0]) - ((pow(Qj_iplus_half_pos[1], 2) + pow(Qj_iplus_half_pos[2], 2) + pow(Qj_iplus_half_pos[3], 2)) / (2.0*pow(Qj_iplus_half_pos[0], 2))))*(metric[Gid].eta_xjp*metric[Gid].eta_xjp + metric[Gid].eta_yjp * metric[Gid].eta_yjp + metric[Gid].eta_zjp * metric[Gid].eta_zjp));

    eigen_Qjm[0] = (metric[Gid].eta_xjp*(Qj_iplus_half_neg[1] / Qj_iplus_half_neg[0]) + metric[Gid].eta_yjp * (Qj_iplus_half_neg[2] / Qj_iplus_half_neg[0]) + metric[Gid].eta_zjp * (Qj_iplus_half_neg[3] / Qj_iplus_half_neg[0])) - sqrt(1.4*0.4*((Qj_iplus_half_neg[4] / Qj_iplus_half_neg[0]) - ((pow(Qj_iplus_half_neg[1], 2) + pow(Qj_iplus_half_neg[2], 2) + pow(Qj_iplus_half_neg[3], 2)) / (2.0*pow(Qj_iplus_half_neg[0], 2))))*(metric[Gid].eta_xjp*metric[Gid].eta_xjp + metric[Gid].eta_yjp * metric[Gid].eta_yjp + metric[Gid].eta_zjp * metric[Gid].eta_zjp));
    eigen_Qjm[1] = (metric[Gid].eta_xjp*(Qj_iplus_half_neg[1] / Qj_iplus_half_neg[0]) + metric[Gid].eta_yjp * (Qj_iplus_half_neg[2] / Qj_iplus_half_neg[0]) + metric[Gid].eta_zjp * (Qj_iplus_half_neg[3] / Qj_iplus_half_neg[0]));
    eigen_Qjm[2] = (metric[Gid].eta_xjp*(Qj_iplus_half_neg[1] / Qj_iplus_half_neg[0]) + metric[Gid].eta_yjp * (Qj_iplus_half_neg[2] / Qj_iplus_half_neg[0]) + metric[Gid].eta_zjp * (Qj_iplus_half_neg[3] / Qj_iplus_half_neg[0]));
    eigen_Qjm[3] = (metric[Gid].eta_xjp*(Qj_iplus_half_neg[1] / Qj_iplus_half_neg[0]) + metric[Gid].eta_yjp * (Qj_iplus_half_neg[2] / Qj_iplus_half_neg[0]) + metric[Gid].eta_zjp * (Qj_iplus_half_neg[3] / Qj_iplus_half_neg[0]));
    eigen_Qjm[4] = (metric[Gid].eta_xjp*(Qj_iplus_half_neg[1] / Qj_iplus_half_neg[0]) + metric[Gid].eta_yjp * (Qj_iplus_half_neg[2] / Qj_iplus_half_neg[0]) + metric[Gid].eta_zjp * (Qj_iplus_half_neg[3] / Qj_iplus_half_neg[0])) + sqrt(1.4*0.4*((Qj_iplus_half_neg[4] / Qj_iplus_half_neg[0]) - ((pow(Qj_iplus_half_neg[1], 2) + pow(Qj_iplus_half_neg[2], 2) + pow(Qj_iplus_half_neg[3], 2)) / (2.0*pow(Qj_iplus_half_neg[0], 2))))*(metric[Gid].eta_xjp*metric[Gid].eta_xjp + metric[Gid].eta_yjp * metric[Gid].eta_yjp + metric[Gid].eta_zjp * metric[Gid].eta_zjp));

    eigen_Qjnp[0] = (metric[Gid].eta_xjm*(Qj_iminus_half_pos[1] / Qj_iminus_half_pos[0]) + metric[Gid].eta_yjm * (Qj_iminus_half_pos[2] / Qj_iminus_half_pos[0]) + metric[Gid].eta_zjm * (Qj_iminus_half_pos[3] / Qj_iminus_half_pos[0])) - sqrt(1.4*0.4*((Qj_iminus_half_pos[4] / Qj_iminus_half_pos[0]) - ((pow(Qj_iminus_half_pos[1], 2) + pow(Qj_iminus_half_pos[2], 2) + pow(Qj_iminus_half_pos[3], 2)) / (2.0*pow(Qj_iminus_half_pos[0], 2))))*(metric[Gid].eta_xjm*metric[Gid].eta_xjm + metric[Gid].eta_yjm * metric[Gid].eta_yjm + metric[Gid].eta_zjm * metric[Gid].eta_zjm));
    eigen_Qjnp[1] = (metric[Gid].eta_xjm*(Qj_iminus_half_pos[1] / Qj_iminus_half_pos[0]) + metric[Gid].eta_yjm * (Qj_iminus_half_pos[2] / Qj_iminus_half_pos[0]) + metric[Gid].eta_zjm * (Qj_iminus_half_pos[3] / Qj_iminus_half_pos[0]));
    eigen_Qjnp[2] = (metric[Gid].eta_xjm*(Qj_iminus_half_pos[1] / Qj_iminus_half_pos[0]) + metric[Gid].eta_yjm * (Qj_iminus_half_pos[2] / Qj_iminus_half_pos[0]) + metric[Gid].eta_zjm * (Qj_iminus_half_pos[3] / Qj_iminus_half_pos[0]));
    eigen_Qjnp[3] = (metric[Gid].eta_xjm*(Qj_iminus_half_pos[1] / Qj_iminus_half_pos[0]) + metric[Gid].eta_yjm * (Qj_iminus_half_pos[2] / Qj_iminus_half_pos[0]) + metric[Gid].eta_zjm * (Qj_iminus_half_pos[3] / Qj_iminus_half_pos[0]));
    eigen_Qjnp[4] = (metric[Gid].eta_xjm*(Qj_iminus_half_pos[1] / Qj_iminus_half_pos[0]) + metric[Gid].eta_yjm * (Qj_iminus_half_pos[2] / Qj_iminus_half_pos[0]) + metric[Gid].eta_zjm * (Qj_iminus_half_pos[3] / Qj_iminus_half_pos[0])) + sqrt(1.4*0.4*((Qj_iminus_half_pos[4] / Qj_iminus_half_pos[0]) - ((pow(Qj_iminus_half_pos[1], 2) + pow(Qj_iminus_half_pos[2], 2) + pow(Qj_iminus_half_pos[3], 2)) / (2.0*pow(Qj_iminus_half_pos[0], 2))))*(metric[Gid].eta_xjm*metric[Gid].eta_xjm + metric[Gid].eta_yjm * metric[Gid].eta_yjm + metric[Gid].eta_zjm * metric[Gid].eta_zjm));

    eigen_Qjnm[0] = (metric[Gid].eta_xjm*(Qj_iminus_half_neg[1] / Qj_iminus_half_neg[0]) + metric[Gid].eta_yjm * (Qj_iminus_half_neg[2] / Qj_iminus_half_neg[0]) + metric[Gid].eta_zjm * (Qj_iminus_half_neg[3] / Qj_iminus_half_neg[0])) - sqrt(1.4*0.4*((Qj_iminus_half_neg[4] / Qj_iminus_half_neg[0]) - ((pow(Qj_iminus_half_neg[1], 2) + pow(Qj_iminus_half_neg[2], 2) + pow(Qj_iminus_half_neg[3], 2)) / (2.0*pow(Qj_iminus_half_neg[0], 2))))*(metric[Gid].eta_xjm*metric[Gid].eta_xjm + metric[Gid].eta_yjm * metric[Gid].eta_yjm + metric[Gid].eta_zjm * metric[Gid].eta_zjm));
    eigen_Qjnm[1] = (metric[Gid].eta_xjm*(Qj_iminus_half_neg[1] / Qj_iminus_half_neg[0]) + metric[Gid].eta_yjm * (Qj_iminus_half_neg[2] / Qj_iminus_half_neg[0]) + metric[Gid].eta_zjm * (Qj_iminus_half_neg[3] / Qj_iminus_half_neg[0]));
    eigen_Qjnm[2] = (metric[Gid].eta_xjm*(Qj_iminus_half_neg[1] / Qj_iminus_half_neg[0]) + metric[Gid].eta_yjm * (Qj_iminus_half_neg[2] / Qj_iminus_half_neg[0]) + metric[Gid].eta_zjm * (Qj_iminus_half_neg[3] / Qj_iminus_half_neg[0]));
    eigen_Qjnm[3] = (metric[Gid].eta_xjm*(Qj_iminus_half_neg[1] / Qj_iminus_half_neg[0]) + metric[Gid].eta_yjm * (Qj_iminus_half_neg[2] / Qj_iminus_half_neg[0]) + metric[Gid].eta_zjm * (Qj_iminus_half_neg[3] / Qj_iminus_half_neg[0]));
    eigen_Qjnm[4] = (metric[Gid].eta_xjm*(Qj_iminus_half_neg[1] / Qj_iminus_half_neg[0]) + metric[Gid].eta_yjm * (Qj_iminus_half_neg[2] / Qj_iminus_half_neg[0]) + metric[Gid].eta_zjm * (Qj_iminus_half_neg[3] / Qj_iminus_half_neg[0])) + sqrt(1.4*0.4*((Qj_iminus_half_neg[4] / Qj_iminus_half_neg[0]) - ((pow(Qj_iminus_half_neg[1], 2) + pow(Qj_iminus_half_neg[2], 2) + pow(Qj_iminus_half_neg[3], 2)) / (2.0*pow(Qj_iminus_half_neg[0], 2))))*(metric[Gid].eta_xjm*metric[Gid].eta_xjm + metric[Gid].eta_yjm * metric[Gid].eta_yjm + metric[Gid].eta_zjm * metric[Gid].eta_zjm));

    eigen_Qkp[0] = (metric[Gid].xi_xkp*(Qk_iplus_half_pos[1] / Qk_iplus_half_pos[0]) + metric[Gid].xi_ykp * (Qk_iplus_half_pos[2] / Qk_iplus_half_pos[0]) + metric[Gid].xi_zkp * (Qk_iplus_half_pos[3] / Qk_iplus_half_pos[0])) - sqrt(1.4*0.4*((Qk_iplus_half_pos[4] / Qk_iplus_half_pos[0]) - ((pow(Qk_iplus_half_pos[1], 2) + pow(Qk_iplus_half_pos[2], 2) + pow(Qk_iplus_half_pos[3], 2)) / (2.0*pow(Qk_iplus_half_pos[0], 2))))*(metric[Gid].xi_xkp*metric[Gid].xi_xkp + metric[Gid].xi_ykp * metric[Gid].xi_ykp + metric[Gid].xi_zkp * metric[Gid].xi_zkp));
    eigen_Qkp[1] = (metric[Gid].xi_xkp*(Qk_iplus_half_pos[1] / Qk_iplus_half_pos[0]) + metric[Gid].xi_ykp * (Qk_iplus_half_pos[2] / Qk_iplus_half_pos[0]) + metric[Gid].xi_zkp * (Qk_iplus_half_pos[3] / Qk_iplus_half_pos[0]));
    eigen_Qkp[2] = (metric[Gid].xi_xkp*(Qk_iplus_half_pos[1] / Qk_iplus_half_pos[0]) + metric[Gid].xi_ykp * (Qk_iplus_half_pos[2] / Qk_iplus_half_pos[0]) + metric[Gid].xi_zkp * (Qk_iplus_half_pos[3] / Qk_iplus_half_pos[0]));
    eigen_Qkp[3] = (metric[Gid].xi_xkp*(Qk_iplus_half_pos[1] / Qk_iplus_half_pos[0]) + metric[Gid].xi_ykp * (Qk_iplus_half_pos[2] / Qk_iplus_half_pos[0]) + metric[Gid].xi_zkp * (Qk_iplus_half_pos[3] / Qk_iplus_half_pos[0]));
    eigen_Qkp[4] = (metric[Gid].xi_xkp*(Qk_iplus_half_pos[1] / Qk_iplus_half_pos[0]) + metric[Gid].xi_ykp * (Qk_iplus_half_pos[2] / Qk_iplus_half_pos[0]) + metric[Gid].xi_zkp * (Qk_iplus_half_pos[3] / Qk_iplus_half_pos[0])) + sqrt(1.4*0.4*((Qk_iplus_half_pos[4] / Qk_iplus_half_pos[0]) - ((pow(Qk_iplus_half_pos[1], 2) + pow(Qk_iplus_half_pos[2], 2) + pow(Qk_iplus_half_pos[3], 2)) / (2.0*pow(Qk_iplus_half_pos[0], 2))))*(metric[Gid].xi_xkp*metric[Gid].xi_xkp + metric[Gid].xi_ykp * metric[Gid].xi_ykp + metric[Gid].xi_zkp * metric[Gid].xi_zkp));

    eigen_Qkm[0] = (metric[Gid].xi_xkp*(Qk_iplus_half_neg[1] / Qk_iplus_half_neg[0]) + metric[Gid].xi_ykp * (Qk_iplus_half_neg[2] / Qk_iplus_half_neg[0]) + metric[Gid].xi_zkp * (Qk_iplus_half_neg[3] / Qk_iplus_half_neg[0])) - sqrt(1.4*0.4*((Qk_iplus_half_neg[4] / Qk_iplus_half_neg[0]) - ((pow(Qk_iplus_half_neg[1], 2) + pow(Qk_iplus_half_neg[2], 2) + pow(Qk_iplus_half_neg[3], 2)) / (2.0*pow(Qk_iplus_half_neg[0], 2))))*(metric[Gid].xi_xkp*metric[Gid].xi_xkp + metric[Gid].xi_ykp * metric[Gid].xi_ykp + metric[Gid].xi_zkp * metric[Gid].xi_zkp));
    eigen_Qkm[1] = (metric[Gid].xi_xkp*(Qk_iplus_half_neg[1] / Qk_iplus_half_neg[0]) + metric[Gid].xi_ykp * (Qk_iplus_half_neg[2] / Qk_iplus_half_neg[0]) + metric[Gid].xi_zkp * (Qk_iplus_half_neg[3] / Qk_iplus_half_neg[0]));
    eigen_Qkm[2] = (metric[Gid].xi_xkp*(Qk_iplus_half_neg[1] / Qk_iplus_half_neg[0]) + metric[Gid].xi_ykp * (Qk_iplus_half_neg[2] / Qk_iplus_half_neg[0]) + metric[Gid].xi_zkp * (Qk_iplus_half_neg[3] / Qk_iplus_half_neg[0]));
    eigen_Qkm[3] = (metric[Gid].xi_xkp*(Qk_iplus_half_neg[1] / Qk_iplus_half_neg[0]) + metric[Gid].xi_ykp * (Qk_iplus_half_neg[2] / Qk_iplus_half_neg[0]) + metric[Gid].xi_zkp * (Qk_iplus_half_neg[3] / Qk_iplus_half_neg[0]));
    eigen_Qkm[4] = (metric[Gid].xi_xkp*(Qk_iplus_half_neg[1] / Qk_iplus_half_neg[0]) + metric[Gid].xi_ykp * (Qk_iplus_half_neg[2] / Qk_iplus_half_neg[0]) + metric[Gid].xi_zkp * (Qk_iplus_half_neg[3] / Qk_iplus_half_neg[0])) + sqrt(1.4*0.4*((Qk_iplus_half_neg[4] / Qk_iplus_half_neg[0]) - ((pow(Qk_iplus_half_neg[1], 2) + pow(Qk_iplus_half_neg[2], 2) + pow(Qk_iplus_half_neg[3], 2)) / (2.0*pow(Qk_iplus_half_neg[0], 2))))*(metric[Gid].xi_xkp*metric[Gid].xi_xkp + metric[Gid].xi_ykp * metric[Gid].xi_ykp + metric[Gid].xi_zkp * metric[Gid].xi_zkp));

    eigen_Qknp[0] = (metric[Gid].xi_xkm*(Qk_iminus_half_pos[1] / Qk_iminus_half_pos[0]) + metric[Gid].xi_ykm * (Qk_iminus_half_pos[2] / Qk_iminus_half_pos[0]) + metric[Gid].xi_zkm * (Qk_iminus_half_pos[3] / Qk_iminus_half_pos[0])) - sqrt(1.4*0.4*((Qk_iminus_half_pos[4] / Qk_iminus_half_pos[0]) - ((pow(Qk_iminus_half_pos[1], 2) + pow(Qk_iminus_half_pos[2], 2) + pow(Qk_iminus_half_pos[3], 2)) / (2.0*pow(Qk_iminus_half_pos[0], 2))))*(metric[Gid].xi_xkm*metric[Gid].xi_xkm + metric[Gid].xi_ykm * metric[Gid].xi_ykm + metric[Gid].xi_zkm * metric[Gid].xi_zkm));
    eigen_Qknp[1] = (metric[Gid].xi_xkm*(Qk_iminus_half_pos[1] / Qk_iminus_half_pos[0]) + metric[Gid].xi_ykm * (Qk_iminus_half_pos[2] / Qk_iminus_half_pos[0]) + metric[Gid].xi_zkm * (Qk_iminus_half_pos[3] / Qk_iminus_half_pos[0]));
    eigen_Qknp[2] = (metric[Gid].xi_xkm*(Qk_iminus_half_pos[1] / Qk_iminus_half_pos[0]) + metric[Gid].xi_ykm * (Qk_iminus_half_pos[2] / Qk_iminus_half_pos[0]) + metric[Gid].xi_zkm * (Qk_iminus_half_pos[3] / Qk_iminus_half_pos[0]));
    eigen_Qknp[3] = (metric[Gid].xi_xkm*(Qk_iminus_half_pos[1] / Qk_iminus_half_pos[0]) + metric[Gid].xi_ykm * (Qk_iminus_half_pos[2] / Qk_iminus_half_pos[0]) + metric[Gid].xi_zkm * (Qk_iminus_half_pos[3] / Qk_iminus_half_pos[0]));
    eigen_Qknp[4] = (metric[Gid].xi_xkm*(Qk_iminus_half_pos[1] / Qk_iminus_half_pos[0]) + metric[Gid].xi_ykm * (Qk_iminus_half_pos[2] / Qk_iminus_half_pos[0]) + metric[Gid].xi_zkm * (Qk_iminus_half_pos[3] / Qk_iminus_half_pos[0])) + sqrt(1.4*0.4*((Qk_iminus_half_pos[4] / Qk_iminus_half_pos[0]) - ((pow(Qk_iminus_half_pos[1], 2) + pow(Qk_iminus_half_pos[2], 2) + pow(Qk_iminus_half_pos[3], 2)) / (2.0*pow(Qk_iminus_half_pos[0], 2))))*(metric[Gid].xi_xkm*metric[Gid].xi_xkm + metric[Gid].xi_ykm * metric[Gid].xi_ykm + metric[Gid].xi_zkm * metric[Gid].xi_zkm));

    eigen_Qknm[0] = (metric[Gid].xi_xkm*(Qk_iminus_half_neg[1] / Qk_iminus_half_neg[0]) + metric[Gid].xi_ykm * (Qk_iminus_half_neg[2] / Qk_iminus_half_neg[0]) + metric[Gid].xi_zkm * (Qk_iminus_half_neg[3] / Qk_iminus_half_neg[0])) - sqrt(1.4*0.4*((Qk_iminus_half_neg[4] / Qk_iminus_half_neg[0]) - ((pow(Qk_iminus_half_neg[1], 2) + pow(Qk_iminus_half_neg[2], 2) + pow(Qk_iminus_half_neg[3], 2)) / (2.0*pow(Qk_iminus_half_neg[0], 2))))*(metric[Gid].xi_xkm*metric[Gid].xi_xkm + metric[Gid].xi_ykm * metric[Gid].xi_ykm + metric[Gid].xi_zkm * metric[Gid].xi_zkm));
    eigen_Qknm[1] = (metric[Gid].xi_xkm*(Qk_iminus_half_neg[1] / Qk_iminus_half_neg[0]) + metric[Gid].xi_ykm * (Qk_iminus_half_neg[2] / Qk_iminus_half_neg[0]) + metric[Gid].xi_zkm * (Qk_iminus_half_neg[3] / Qk_iminus_half_neg[0]));
    eigen_Qknm[2] = (metric[Gid].xi_xkm*(Qk_iminus_half_neg[1] / Qk_iminus_half_neg[0]) + metric[Gid].xi_ykm * (Qk_iminus_half_neg[2] / Qk_iminus_half_neg[0]) + metric[Gid].xi_zkm * (Qk_iminus_half_neg[3] / Qk_iminus_half_neg[0]));
    eigen_Qknm[3] = (metric[Gid].xi_xkm*(Qk_iminus_half_neg[1] / Qk_iminus_half_neg[0]) + metric[Gid].xi_ykm * (Qk_iminus_half_neg[2] / Qk_iminus_half_neg[0]) + metric[Gid].xi_zkm * (Qk_iminus_half_neg[3] / Qk_iminus_half_neg[0]));
    eigen_Qknm[4] = (metric[Gid].xi_xkm*(Qk_iminus_half_neg[1] / Qk_iminus_half_neg[0]) + metric[Gid].xi_ykm * (Qk_iminus_half_neg[2] / Qk_iminus_half_neg[0]) + metric[Gid].xi_zkm * (Qk_iminus_half_neg[3] / Qk_iminus_half_neg[0])) + sqrt(1.4*0.4*((Qk_iminus_half_neg[4] / Qk_iminus_half_neg[0]) - ((pow(Qk_iminus_half_neg[1], 2) + pow(Qk_iminus_half_neg[2], 2) + pow(Qk_iminus_half_neg[3], 2)) / (2.0*pow(Qk_iminus_half_neg[0], 2))))*(metric[Gid].xi_xkm*metric[Gid].xi_xkm + metric[Gid].xi_ykm * metric[Gid].xi_ykm + metric[Gid].xi_zkm * metric[Gid].xi_zkm));

    /********************************************************************************************************************************************************************************************************************************************************************************************************************* */

    F_ip_pos_comp[0] = DET_IP*(metric[Gid].zeta_xip*F_ip_pos[0] + metric[Gid].zeta_yip * E_ip_pos[0] + metric[Gid].zeta_zip * G_ip_pos[0]);
    F_ip_pos_comp[1] = DET_IP*(metric[Gid].zeta_xip*F_ip_pos[1] + metric[Gid].zeta_yip * E_ip_pos[1] + metric[Gid].zeta_zip * G_ip_pos[1]);
    F_ip_pos_comp[2] = DET_IP*(metric[Gid].zeta_xip*F_ip_pos[2] + metric[Gid].zeta_yip * E_ip_pos[2] + metric[Gid].zeta_zip * G_ip_pos[2]);
    F_ip_pos_comp[3] = DET_IP*(metric[Gid].zeta_xip*F_ip_pos[3] + metric[Gid].zeta_yip * E_ip_pos[3] + metric[Gid].zeta_zip * G_ip_pos[3]);
    F_ip_pos_comp[4] = DET_IP*(metric[Gid].zeta_xip*F_ip_pos[4] + metric[Gid].zeta_yip * E_ip_pos[4] + metric[Gid].zeta_zip * G_ip_pos[4]);
    
    F_ip_neg_comp[0] = DET_IP*(metric[Gid].zeta_xip*F_ip_neg[0] + metric[Gid].zeta_yip * E_ip_neg[0] + metric[Gid].zeta_zip * G_ip_neg[0]);
    F_ip_neg_comp[1] = DET_IP*(metric[Gid].zeta_xip*F_ip_neg[1] + metric[Gid].zeta_yip * E_ip_neg[1] + metric[Gid].zeta_zip * G_ip_neg[1]);
    F_ip_neg_comp[2] = DET_IP*(metric[Gid].zeta_xip*F_ip_neg[2] + metric[Gid].zeta_yip * E_ip_neg[2] + metric[Gid].zeta_zip * G_ip_neg[2]);
    F_ip_neg_comp[3] = DET_IP*(metric[Gid].zeta_xip*F_ip_neg[3] + metric[Gid].zeta_yip * E_ip_neg[3] + metric[Gid].zeta_zip * G_ip_neg[3]);
    F_ip_neg_comp[4] = DET_IP*(metric[Gid].zeta_xip*F_ip_neg[4] + metric[Gid].zeta_yip * E_ip_neg[4] + metric[Gid].zeta_zip * G_ip_neg[4]);

    F_im_pos_comp[0] = DET_IM*(metric[Gid].zeta_xim*F_im_pos[0] + metric[Gid].zeta_yim * E_im_pos[0] + metric[Gid].zeta_zim * G_im_pos[0]);
    F_im_pos_comp[1] = DET_IM*(metric[Gid].zeta_xim*F_im_pos[1] + metric[Gid].zeta_yim * E_im_pos[1] + metric[Gid].zeta_zim * G_im_pos[1]);
    F_im_pos_comp[2] = DET_IM*(metric[Gid].zeta_xim*F_im_pos[2] + metric[Gid].zeta_yim * E_im_pos[2] + metric[Gid].zeta_zim * G_im_pos[2]);
    F_im_pos_comp[3] = DET_IM*(metric[Gid].zeta_xim*F_im_pos[3] + metric[Gid].zeta_yim * E_im_pos[3] + metric[Gid].zeta_zim * G_im_pos[3]);
    F_im_pos_comp[4] = DET_IM*(metric[Gid].zeta_xim*F_im_pos[4] + metric[Gid].zeta_yim * E_im_pos[4] + metric[Gid].zeta_zim * G_im_pos[4]);

    F_im_neg_comp[0] = DET_IM*(metric[Gid].zeta_xim*F_im_neg[0] + metric[Gid].zeta_yim * E_im_neg[0] + metric[Gid].zeta_zim * G_im_neg[0]);
    F_im_neg_comp[1] = DET_IM*(metric[Gid].zeta_xim*F_im_neg[1] + metric[Gid].zeta_yim * E_im_neg[1] + metric[Gid].zeta_zim * G_im_neg[1]);
    F_im_neg_comp[2] = DET_IM*(metric[Gid].zeta_xim*F_im_neg[2] + metric[Gid].zeta_yim * E_im_neg[2] + metric[Gid].zeta_zim * G_im_neg[2]);
    F_im_neg_comp[3] = DET_IM*(metric[Gid].zeta_xim*F_im_neg[3] + metric[Gid].zeta_yim * E_im_neg[3] + metric[Gid].zeta_zim * G_im_neg[3]);
    F_im_neg_comp[4] = DET_IM*(metric[Gid].zeta_xim*F_im_neg[4] + metric[Gid].zeta_yim * E_im_neg[4] + metric[Gid].zeta_zim * G_im_neg[4]);

    E_ip_pos_comp[0] = DET_JP*(metric[Gid].eta_xjp*F_jp_pos[0] + metric[Gid].eta_yjp * E_jp_pos[0] + metric[Gid].eta_zjp * G_jp_pos[0]);
    E_ip_pos_comp[1] = DET_JP*(metric[Gid].eta_xjp*F_jp_pos[1] + metric[Gid].eta_yjp * E_jp_pos[1] + metric[Gid].eta_zjp * G_jp_pos[1]);
    E_ip_pos_comp[2] = DET_JP*(metric[Gid].eta_xjp*F_jp_pos[2] + metric[Gid].eta_yjp * E_jp_pos[2] + metric[Gid].eta_zjp * G_jp_pos[2]);
    E_ip_pos_comp[3] = DET_JP*(metric[Gid].eta_xjp*F_jp_pos[3] + metric[Gid].eta_yjp * E_jp_pos[3] + metric[Gid].eta_zjp * G_jp_pos[3]);
    E_ip_pos_comp[4] = DET_JP*(metric[Gid].eta_xjp*F_jp_pos[4] + metric[Gid].eta_yjp * E_jp_pos[4] + metric[Gid].eta_zjp * G_jp_pos[4]);

    E_ip_neg_comp[0] = DET_JP*(metric[Gid].eta_xjp*F_jp_neg[0] + metric[Gid].eta_yjp * E_jp_neg[0] + metric[Gid].eta_zjp * G_jp_neg[0]);
    E_ip_neg_comp[1] = DET_JP*(metric[Gid].eta_xjp*F_jp_neg[1] + metric[Gid].eta_yjp * E_jp_neg[1] + metric[Gid].eta_zjp * G_jp_neg[1]);
    E_ip_neg_comp[2] = DET_JP*(metric[Gid].eta_xjp*F_jp_neg[2] + metric[Gid].eta_yjp * E_jp_neg[2] + metric[Gid].eta_zjp * G_jp_neg[2]);
    E_ip_neg_comp[3] = DET_JP*(metric[Gid].eta_xjp*F_jp_neg[3] + metric[Gid].eta_yjp * E_jp_neg[3] + metric[Gid].eta_zjp * G_jp_neg[3]);
    E_ip_neg_comp[4] = DET_JP*(metric[Gid].eta_xjp*F_jp_neg[4] + metric[Gid].eta_yjp * E_jp_neg[4] + metric[Gid].eta_zjp * G_jp_neg[4]);

    E_im_pos_comp[0] = DET_JM*(metric[Gid].eta_xjm*F_jm_pos[0] + metric[Gid].eta_yjm * E_jm_pos[0] + metric[Gid].eta_zjm * G_jm_pos[0]);
    E_im_pos_comp[1] = DET_JM*(metric[Gid].eta_xjm*F_jm_pos[1] + metric[Gid].eta_yjm * E_jm_pos[1] + metric[Gid].eta_zjm * G_jm_pos[1]);
    E_im_pos_comp[2] = DET_JM*(metric[Gid].eta_xjm*F_jm_pos[2] + metric[Gid].eta_yjm * E_jm_pos[2] + metric[Gid].eta_zjm * G_jm_pos[2]);
    E_im_pos_comp[3] = DET_JM*(metric[Gid].eta_xjm*F_jm_pos[3] + metric[Gid].eta_yjm * E_jm_pos[3] + metric[Gid].eta_zjm * G_jm_pos[3]);
    E_im_pos_comp[4] = DET_JM*(metric[Gid].eta_xjm*F_jm_pos[4] + metric[Gid].eta_yjm * E_jm_pos[4] + metric[Gid].eta_zjm * G_jm_pos[4]);

    E_im_neg_comp[0] = DET_JM*(metric[Gid].eta_xjm*F_jm_neg[0] + metric[Gid].eta_yjm * E_jm_neg[0] + metric[Gid].eta_zjm * G_jm_neg[0]);
    E_im_neg_comp[1] = DET_JM*(metric[Gid].eta_xjm*F_jm_neg[1] + metric[Gid].eta_yjm * E_jm_neg[1] + metric[Gid].eta_zjm * G_jm_neg[1]);
    E_im_neg_comp[2] = DET_JM*(metric[Gid].eta_xjm*F_jm_neg[2] + metric[Gid].eta_yjm * E_jm_neg[2] + metric[Gid].eta_zjm * G_jm_neg[2]);
    E_im_neg_comp[3] = DET_JM*(metric[Gid].eta_xjm*F_jm_neg[3] + metric[Gid].eta_yjm * E_jm_neg[3] + metric[Gid].eta_zjm * G_jm_neg[3]);
    E_im_neg_comp[4] = DET_JM*(metric[Gid].eta_xjm*F_jm_neg[4] + metric[Gid].eta_yjm * E_jm_neg[4] + metric[Gid].eta_zjm * G_jm_neg[4]);

    G_ip_pos_comp[0] = DET_KP*(metric[Gid].xi_xkp*F_kp_pos[0] + metric[Gid].xi_ykp * E_kp_pos[0] + metric[Gid].xi_zkp * G_kp_pos[0]);
    G_ip_pos_comp[1] = DET_KP*(metric[Gid].xi_xkp*F_kp_pos[1] + metric[Gid].xi_ykp * E_kp_pos[1] + metric[Gid].xi_zkp * G_kp_pos[1]);
    G_ip_pos_comp[2] = DET_KP*(metric[Gid].xi_xkp*F_kp_pos[2] + metric[Gid].xi_ykp * E_kp_pos[2] + metric[Gid].xi_zkp * G_kp_pos[2]);
    G_ip_pos_comp[3] = DET_KP*(metric[Gid].xi_xkp*F_kp_pos[3] + metric[Gid].xi_ykp * E_kp_pos[3] + metric[Gid].xi_zkp * G_kp_pos[3]);
    G_ip_pos_comp[4] = DET_KP*(metric[Gid].xi_xkp*F_kp_pos[4] + metric[Gid].xi_ykp * E_kp_pos[4] + metric[Gid].xi_zkp * G_kp_pos[4]);

    G_ip_neg_comp[0] = DET_KP*(metric[Gid].xi_xkp*F_kp_neg[0] + metric[Gid].xi_ykp * E_kp_neg[0] + metric[Gid].xi_zkp * G_kp_neg[0]);
    G_ip_neg_comp[1] = DET_KP*(metric[Gid].xi_xkp*F_kp_neg[1] + metric[Gid].xi_ykp * E_kp_neg[1] + metric[Gid].xi_zkp * G_kp_neg[1]);
    G_ip_neg_comp[2] = DET_KP*(metric[Gid].xi_xkp*F_kp_neg[2] + metric[Gid].xi_ykp * E_kp_neg[2] + metric[Gid].xi_zkp * G_kp_neg[2]);
    G_ip_neg_comp[3] = DET_KP*(metric[Gid].xi_xkp*F_kp_neg[3] + metric[Gid].xi_ykp * E_kp_neg[3] + metric[Gid].xi_zkp * G_kp_neg[3]);
    G_ip_neg_comp[4] = DET_KP*(metric[Gid].xi_xkp*F_kp_neg[4] + metric[Gid].xi_ykp * E_kp_neg[4] + metric[Gid].xi_zkp * G_kp_neg[4]);

    G_im_pos_comp[0] = DET_KM*(metric[Gid].xi_xkm*F_km_pos[0] + metric[Gid].xi_ykm * E_km_pos[0] + metric[Gid].xi_zkm * G_km_pos[0]);
    G_im_pos_comp[1] = DET_KM*(metric[Gid].xi_xkm*F_km_pos[1] + metric[Gid].xi_ykm * E_km_pos[1] + metric[Gid].xi_zkm * G_km_pos[1]);
    G_im_pos_comp[2] = DET_KM*(metric[Gid].xi_xkm*F_km_pos[2] + metric[Gid].xi_ykm * E_km_pos[2] + metric[Gid].xi_zkm * G_km_pos[2]);
    G_im_pos_comp[3] = DET_KM*(metric[Gid].xi_xkm*F_km_pos[3] + metric[Gid].xi_ykm * E_km_pos[3] + metric[Gid].xi_zkm * G_km_pos[3]);
    G_im_pos_comp[4] = DET_KM*(metric[Gid].xi_xkm*F_km_pos[4] + metric[Gid].xi_ykm * E_km_pos[4] + metric[Gid].xi_zkm * G_km_pos[4]);

    G_im_neg_comp[0] = DET_KM*(metric[Gid].xi_xkm*F_km_neg[0] + metric[Gid].xi_ykm * E_km_neg[0] + metric[Gid].xi_zkm * G_km_neg[0]);
    G_im_neg_comp[1] = DET_KM*(metric[Gid].xi_xkm*F_km_neg[1] + metric[Gid].xi_ykm * E_km_neg[1] + metric[Gid].xi_zkm * G_km_neg[1]);
    G_im_neg_comp[2] = DET_KM*(metric[Gid].xi_xkm*F_km_neg[2] + metric[Gid].xi_ykm * E_km_neg[2] + metric[Gid].xi_zkm * G_km_neg[2]);
    G_im_neg_comp[3] = DET_KM*(metric[Gid].xi_xkm*F_km_neg[3] + metric[Gid].xi_ykm * E_km_neg[3] + metric[Gid].xi_zkm * G_km_neg[3]);
    G_im_neg_comp[4] = DET_KM*(metric[Gid].xi_xkm*F_km_neg[4] + metric[Gid].xi_ykm * E_km_neg[4] + metric[Gid].xi_zkm * G_km_neg[4]);

    Qi_iplus_half_pos[0] = DET_IP*Qi_iplus_half_pos[0];
    Qi_iplus_half_pos[1] = DET_IP*Qi_iplus_half_pos[1];
    Qi_iplus_half_pos[2] = DET_IP*Qi_iplus_half_pos[2];
    Qi_iplus_half_pos[3] = DET_IP*Qi_iplus_half_pos[3];
    Qi_iplus_half_pos[4] = DET_IP*Qi_iplus_half_pos[4];

    Qi_iplus_half_neg[0] = DET_IP*Qi_iplus_half_neg[0];
    Qi_iplus_half_neg[1] = DET_IP*Qi_iplus_half_neg[1];
    Qi_iplus_half_neg[2] = DET_IP*Qi_iplus_half_neg[2];
    Qi_iplus_half_neg[3] = DET_IP*Qi_iplus_half_neg[3];
    Qi_iplus_half_neg[4] = DET_IP*Qi_iplus_half_neg[4];

    Qi_iminus_half_pos[0] = DET_IM*Qi_iminus_half_pos[0];
    Qi_iminus_half_pos[1] = DET_IM*Qi_iminus_half_pos[1];
    Qi_iminus_half_pos[2] = DET_IM*Qi_iminus_half_pos[2];
    Qi_iminus_half_pos[3] = DET_IM*Qi_iminus_half_pos[3];
    Qi_iminus_half_pos[4] = DET_IM*Qi_iminus_half_pos[4];

    Qi_iminus_half_neg[0] = DET_IM*Qi_iminus_half_neg[0];
    Qi_iminus_half_neg[1] = DET_IM*Qi_iminus_half_neg[1];
    Qi_iminus_half_neg[2] = DET_IM*Qi_iminus_half_neg[2];
    Qi_iminus_half_neg[3] = DET_IM*Qi_iminus_half_neg[3];
    Qi_iminus_half_neg[4] = DET_IM*Qi_iminus_half_neg[4];

    Qj_iplus_half_pos[0] = DET_JP*Qj_iplus_half_pos[0];
    Qj_iplus_half_pos[1] = DET_JP*Qj_iplus_half_pos[1];
    Qj_iplus_half_pos[2] = DET_JP*Qj_iplus_half_pos[2];
    Qj_iplus_half_pos[3] = DET_JP*Qj_iplus_half_pos[3];
    Qj_iplus_half_pos[4] = DET_JP*Qj_iplus_half_pos[4];

    Qj_iplus_half_neg[0] = DET_JP*Qj_iplus_half_neg[0];
    Qj_iplus_half_neg[1] = DET_JP*Qj_iplus_half_neg[1];
    Qj_iplus_half_neg[2] = DET_JP*Qj_iplus_half_neg[2];
    Qj_iplus_half_neg[3] = DET_JP*Qj_iplus_half_neg[3];
    Qj_iplus_half_neg[4] = DET_JP*Qj_iplus_half_neg[4];

    Qj_iminus_half_pos[0] = DET_JM*Qj_iminus_half_pos[0];
    Qj_iminus_half_pos[1] = DET_JM*Qj_iminus_half_pos[1];
    Qj_iminus_half_pos[2] = DET_JM*Qj_iminus_half_pos[2];
    Qj_iminus_half_pos[3] = DET_JM*Qj_iminus_half_pos[3];
    Qj_iminus_half_pos[4] = DET_JM*Qj_iminus_half_pos[4];

    Qj_iminus_half_neg[0] = DET_JM*Qj_iminus_half_neg[0];
    Qj_iminus_half_neg[1] = DET_JM*Qj_iminus_half_neg[1];
    Qj_iminus_half_neg[2] = DET_JM*Qj_iminus_half_neg[2];
    Qj_iminus_half_neg[3] = DET_JM*Qj_iminus_half_neg[3];
    Qj_iminus_half_neg[4] = DET_JM*Qj_iminus_half_neg[4];

    Qk_iplus_half_pos[0] = DET_KP*Qk_iplus_half_pos[0];
    Qk_iplus_half_pos[1] = DET_KP*Qk_iplus_half_pos[1];
    Qk_iplus_half_pos[2] = DET_KP*Qk_iplus_half_pos[2];
    Qk_iplus_half_pos[3] = DET_KP*Qk_iplus_half_pos[3];
    Qk_iplus_half_pos[4] = DET_KP*Qk_iplus_half_pos[4];

    Qk_iplus_half_neg[0] = DET_KP*Qk_iplus_half_neg[0];
    Qk_iplus_half_neg[1] = DET_KP*Qk_iplus_half_neg[1];
    Qk_iplus_half_neg[2] = DET_KP*Qk_iplus_half_neg[2];
    Qk_iplus_half_neg[3] = DET_KP*Qk_iplus_half_neg[3];
    Qk_iplus_half_neg[4] = DET_KP*Qk_iplus_half_neg[4];

    Qk_iminus_half_pos[0] = DET_KM*Qk_iminus_half_pos[0];
    Qk_iminus_half_pos[1] = DET_KM*Qk_iminus_half_pos[1];
    Qk_iminus_half_pos[2] = DET_KM*Qk_iminus_half_pos[2];
    Qk_iminus_half_pos[3] = DET_KM*Qk_iminus_half_pos[3];
    Qk_iminus_half_pos[4] = DET_KM*Qk_iminus_half_pos[4];

    Qk_iminus_half_neg[0] = DET_KM*Qk_iminus_half_neg[0];
    Qk_iminus_half_neg[1] = DET_KM*Qk_iminus_half_neg[1];
    Qk_iminus_half_neg[2] = DET_KM*Qk_iminus_half_neg[2];
    Qk_iminus_half_neg[3] = DET_KM*Qk_iminus_half_neg[3];
    Qk_iminus_half_neg[4] = DET_KM*Qk_iminus_half_neg[4];

    /******************************************************************************************************************************************************/
    /******************************************************CENTRAL 4th ORDER VISCOUS TERMS*****************************************************************/

    if (LOC_I > 0 && LOC_I <= 52 && (LOC_IP1 == 100))
    {
        dFv[0] = (1.0 / 6.0)*(11.0*Fv[Gid].N0 - 18.0*Fv[IM1].N0 + 9.0*Fv[IM2].N0 - 2.0*Fv[IM3].N0);
        dFv[1] = (1.0 / 6.0)*(11.0*Fv[Gid].N1 - 18.0*Fv[IM1].N1 + 9.0*Fv[IM2].N1 - 2.0*Fv[IM3].N1);
        dFv[2] = (1.0 / 6.0)*(11.0*Fv[Gid].N2 - 18.0*Fv[IM1].N2 + 9.0*Fv[IM2].N2 - 2.0*Fv[IM3].N2);
        dFv[3] = (1.0 / 6.0)*(11.0*Fv[Gid].N3 - 18.0*Fv[IM1].N3 + 9.0*Fv[IM2].N3 - 2.0*Fv[IM3].N3);
        dFv[4] = (1.0 / 6.0)*(11.0*Fv[Gid].N4 - 18.0*Fv[IM1].N4 + 9.0*Fv[IM2].N4 - 2.0*Fv[IM3].N4);
    }
    if (LOC_I > 0 && LOC_I <= 52 && (LOC_IM1 == 100))
    {
        dFv[0] = (1.0 / 6.0)*(-11.0*Fv[Gid].N0 + 18.0*Fv[IP1].N0 - 9.0*Fv[IP2].N0 + 2.0*Fv[IP3].N0);
        dFv[1] = (1.0 / 6.0)*(-11.0*Fv[Gid].N1 + 18.0*Fv[IP1].N1 - 9.0*Fv[IP2].N1 + 2.0*Fv[IP3].N1);
        dFv[2] = (1.0 / 6.0)*(-11.0*Fv[Gid].N2 + 18.0*Fv[IP1].N2 - 9.0*Fv[IP2].N2 + 2.0*Fv[IP3].N2);
        dFv[3] = (1.0 / 6.0)*(-11.0*Fv[Gid].N3 + 18.0*Fv[IP1].N3 - 9.0*Fv[IP2].N3 + 2.0*Fv[IP3].N3);
        dFv[4] = (1.0 / 6.0)*(-11.0*Fv[Gid].N4 + 18.0*Fv[IP1].N4 - 9.0*Fv[IP2].N4 + 2.0*Fv[IP3].N4);
    }
    if (LOC_I == 0 && LOC_IP1 > 0 && LOC_IP1 <= 52 && (LOC_IP2 == 100))
    {
        dFv[0] = (1.0 / 6.0)*((2.0)*Fv[IP1].N0 + 3.0*Fv[Gid].N0 - 6.0*Fv[IM1].N0 + Fv[IM2].N0);
        dFv[1] = (1.0 / 6.0)*((2.0)*Fv[IP1].N1 + 3.0*Fv[Gid].N1 - 6.0*Fv[IM1].N1 + Fv[IM2].N1);
        dFv[2] = (1.0 / 6.0)*((2.0)*Fv[IP1].N2 + 3.0*Fv[Gid].N2 - 6.0*Fv[IM1].N2 + Fv[IM2].N2);
        dFv[3] = (1.0 / 6.0)*((2.0)*Fv[IP1].N3 + 3.0*Fv[Gid].N3 - 6.0*Fv[IM1].N3 + Fv[IM2].N3);
        dFv[4] = (1.0 / 6.0)*((2.0)*Fv[IP1].N4 + 3.0*Fv[Gid].N4 - 6.0*Fv[IM1].N4 + Fv[IM2].N4);
    }
    if (LOC_I == 0 && LOC_IM1 > 0 && LOC_IM1 <= 52 && (LOC_IM2 == 100))
    {
        dFv[0] = (1.0 / 6.0)*((-2.0)*Fv[IM1].N0 - 3.0*Fv[Gid].N0 + 6.0*Fv[IP1].N0 - Fv[IP2].N0);
        dFv[1] = (1.0 / 6.0)*((-2.0)*Fv[IM1].N1 - 3.0*Fv[Gid].N1 + 6.0*Fv[IP1].N1 - Fv[IP2].N1);
        dFv[2] = (1.0 / 6.0)*((-2.0)*Fv[IM1].N2 - 3.0*Fv[Gid].N2 + 6.0*Fv[IP1].N2 - Fv[IP2].N2);
        dFv[3] = (1.0 / 6.0)*((-2.0)*Fv[IM1].N3 - 3.0*Fv[Gid].N3 + 6.0*Fv[IP1].N3 - Fv[IP2].N3);
        dFv[4] = (1.0 / 6.0)*((-2.0)*Fv[IM1].N4 - 3.0*Fv[Gid].N4 + 6.0*Fv[IP1].N4 - Fv[IP2].N4);
    }
    if (LOC_I == 0 && LOC_IM2 != 100 && LOC_IM1 != 100 && LOC_IP2 != 100 && LOC_IP1 != 100)
    {
        dFv[0] = (1.0 / 12.0)*(8.0*(Fv[IP1].N0 - Fv[IM1].N0) - (Fv[IP2].N0 - Fv[IM2].N0));
        dFv[1] = (1.0 / 12.0)*(8.0*(Fv[IP1].N1 - Fv[IM1].N1) - (Fv[IP2].N1 - Fv[IM2].N1));
        dFv[2] = (1.0 / 12.0)*(8.0*(Fv[IP1].N2 - Fv[IM1].N2) - (Fv[IP2].N2 - Fv[IM2].N2));
        dFv[3] = (1.0 / 12.0)*(8.0*(Fv[IP1].N3 - Fv[IM1].N3) - (Fv[IP2].N3 - Fv[IM2].N3));
        dFv[4] = (1.0 / 12.0)*(8.0*(Fv[IP1].N4 - Fv[IM1].N4) - (Fv[IP2].N4 - Fv[IM2].N4));
    }


    if (LOC_I > 0 && LOC_I <= 52 && (LOC_JP1 == 100))
    {
        dEv[0] = (1.0 / 6.0)*(11.0*Ev[Gid].N0 - 18.0*Ev[JM1].N0 + 9.0*Ev[JM2].N0 - 2.0*Ev[JM3].N0);
        dEv[1] = (1.0 / 6.0)*(11.0*Ev[Gid].N1 - 18.0*Ev[JM1].N1 + 9.0*Ev[JM2].N1 - 2.0*Ev[JM3].N1);
        dEv[2] = (1.0 / 6.0)*(11.0*Ev[Gid].N2 - 18.0*Ev[JM1].N2 + 9.0*Ev[JM2].N2 - 2.0*Ev[JM3].N2);
        dEv[3] = (1.0 / 6.0)*(11.0*Ev[Gid].N3 - 18.0*Ev[JM1].N3 + 9.0*Ev[JM2].N3 - 2.0*Ev[JM3].N3);
        dEv[4] = (1.0 / 6.0)*(11.0*Ev[Gid].N4 - 18.0*Ev[JM1].N4 + 9.0*Ev[JM2].N4 - 2.0*Ev[JM3].N4);
    }
    if (LOC_I > 0 && LOC_I <= 52 && (LOC_JM1 == 100))
    {
        dEv[0] = (1.0 / 6.0)*(-11.0*Ev[Gid].N0 + 18.0*Ev[JP1].N0 - 9.0*Ev[JP2].N0 + 2.0*Ev[JP3].N0);
        dEv[1] = (1.0 / 6.0)*(-11.0*Ev[Gid].N1 + 18.0*Ev[JP1].N1 - 9.0*Ev[JP2].N1 + 2.0*Ev[JP3].N1);
        dEv[2] = (1.0 / 6.0)*(-11.0*Ev[Gid].N2 + 18.0*Ev[JP1].N2 - 9.0*Ev[JP2].N2 + 2.0*Ev[JP3].N2);
        dEv[3] = (1.0 / 6.0)*(-11.0*Ev[Gid].N3 + 18.0*Ev[JP1].N3 - 9.0*Ev[JP2].N3 + 2.0*Ev[JP3].N3);
        dEv[4] = (1.0 / 6.0)*(-11.0*Ev[Gid].N4 + 18.0*Ev[JP1].N4 - 9.0*Ev[JP2].N4 + 2.0*Ev[JP3].N4);
    }
    if (LOC_I == 0 && LOC_JP1 > 0 && LOC_JP1 <= 52 && (LOC_JP2 == 100))
    {
        dEv[0] = (1.0 / 6.0)*(2.0*Ev[JP1].N0 + 3.0*Ev[Gid].N0 - 6.0*Ev[JM1].N0 + Ev[JM2].N0);
        dEv[1] = (1.0 / 6.0)*(2.0*Ev[JP1].N1 + 3.0*Ev[Gid].N1 - 6.0*Ev[JM1].N1 + Ev[JM2].N1);
        dEv[2] = (1.0 / 6.0)*(2.0*Ev[JP1].N2 + 3.0*Ev[Gid].N2 - 6.0*Ev[JM1].N2 + Ev[JM2].N2);
        dEv[3] = (1.0 / 6.0)*(2.0*Ev[JP1].N3 + 3.0*Ev[Gid].N3 - 6.0*Ev[JM1].N3 + Ev[JM2].N3);
        dEv[4] = (1.0 / 6.0)*(2.0*Ev[JP1].N4 + 3.0*Ev[Gid].N4 - 6.0*Ev[JM1].N4 + Ev[JM2].N4);
    }
    if (LOC_I == 0 && LOC_JM1 > 0 && LOC_JM1 <= 52 && (LOC_JM2 == 100))
    {
        dEv[0] = (1.0 / 6.0)*(-2.0*Ev[JM1].N0 - 3.0*Ev[Gid].N0 + 6.0*Ev[JP1].N0 - Ev[JP2].N0);
        dEv[1] = (1.0 / 6.0)*(-2.0*Ev[JM1].N1 - 3.0*Ev[Gid].N1 + 6.0*Ev[JP1].N1 - Ev[JP2].N1);
        dEv[2] = (1.0 / 6.0)*(-2.0*Ev[JM1].N2 - 3.0*Ev[Gid].N2 + 6.0*Ev[JP1].N2 - Ev[JP2].N2);
        dEv[3] = (1.0 / 6.0)*(-2.0*Ev[JM1].N3 - 3.0*Ev[Gid].N3 + 6.0*Ev[JP1].N3 - Ev[JP2].N3);
        dEv[4] = (1.0 / 6.0)*(-2.0*Ev[JM1].N4 - 3.0*Ev[Gid].N4 + 6.0*Ev[JP1].N4 - Ev[JP2].N4);
    }
    if (LOC_I == 0 && LOC_JM2 != 100 && LOC_JM1 != 100 && LOC_JP2 != 100 && LOC_JP1 != 100)
    {
        dEv[0] = (1.0 / 12.0)*(8.0*(Ev[JP1].N0 - Ev[JM1].N0) - (Ev[JP2].N0 - Ev[JM2].N0));
        dEv[1] = (1.0 / 12.0)*(8.0*(Ev[JP1].N1 - Ev[JM1].N1) - (Ev[JP2].N1 - Ev[JM2].N1));
        dEv[2] = (1.0 / 12.0)*(8.0*(Ev[JP1].N2 - Ev[JM1].N2) - (Ev[JP2].N2 - Ev[JM2].N2));
        dEv[3] = (1.0 / 12.0)*(8.0*(Ev[JP1].N3 - Ev[JM1].N3) - (Ev[JP2].N3 - Ev[JM2].N3));
        dEv[4] = (1.0 / 12.0)*(8.0*(Ev[JP1].N4 - Ev[JM1].N4) - (Ev[JP2].N4 - Ev[JM2].N4));
    }


    if (LOC_I > 0 && LOC_I <= 52 && (LOC_KP1 == 100))
    {
        dGv[0] = (1.0 / 6.0)*(11.0*Gv[Gid].N0 - 18.0*Gv[KM1].N0 + 9.0*Gv[KM2].N0 - 2.0*Gv[KM3].N0);
        dGv[1] = (1.0 / 6.0)*(11.0*Gv[Gid].N1 - 18.0*Gv[KM1].N1 + 9.0*Gv[KM2].N1 - 2.0*Gv[KM3].N1);
        dGv[2] = (1.0 / 6.0)*(11.0*Gv[Gid].N2 - 18.0*Gv[KM1].N2 + 9.0*Gv[KM2].N2 - 2.0*Gv[KM3].N2);
        dGv[3] = (1.0 / 6.0)*(11.0*Gv[Gid].N3 - 18.0*Gv[KM1].N3 + 9.0*Gv[KM2].N3 - 2.0*Gv[KM3].N3);
        dGv[4] = (1.0 / 6.0)*(11.0*Gv[Gid].N4 - 18.0*Gv[KM1].N4 + 9.0*Gv[KM2].N4 - 2.0*Gv[KM3].N4);
    }
    if (LOC_I > 0 && LOC_I <= 52 && (LOC_KM1 == 100))
    {
        dGv[0] = (1.0 / 6.0)*(-11.0*Gv[Gid].N0 + 18.0*Gv[KP1].N0 - 9.0*Gv[KP2].N0 + 2.0*Gv[KP3].N0);
        dGv[1] = (1.0 / 6.0)*(-11.0*Gv[Gid].N1 + 18.0*Gv[KP1].N1 - 9.0*Gv[KP2].N1 + 2.0*Gv[KP3].N1);
        dGv[2] = (1.0 / 6.0)*(-11.0*Gv[Gid].N2 + 18.0*Gv[KP1].N2 - 9.0*Gv[KP2].N2 + 2.0*Gv[KP3].N2);
        dGv[3] = (1.0 / 6.0)*(-11.0*Gv[Gid].N3 + 18.0*Gv[KP1].N3 - 9.0*Gv[KP2].N3 + 2.0*Gv[KP3].N3);
        dGv[4] = (1.0 / 6.0)*(-11.0*Gv[Gid].N4 + 18.0*Gv[KP1].N4 - 9.0*Gv[KP2].N4 + 2.0*Gv[KP3].N4);
    }
    if (LOC_I == 0 && LOC_KP1 > 0 && LOC_KP1 <= 52 && (LOC_KP2 == 100))
    {
        dGv[0] = (1.0 / 6.0)*(2.0*Gv[KP1].N0 + 3.0*Gv[Gid].N0 - 6.0*Gv[KM1].N0 + Gv[KM2].N0);
        dGv[1] = (1.0 / 6.0)*(2.0*Gv[KP1].N1 + 3.0*Gv[Gid].N1 - 6.0*Gv[KM1].N1 + Gv[KM2].N1);
        dGv[2] = (1.0 / 6.0)*(2.0*Gv[KP1].N2 + 3.0*Gv[Gid].N2 - 6.0*Gv[KM1].N2 + Gv[KM2].N2);
        dGv[3] = (1.0 / 6.0)*(2.0*Gv[KP1].N3 + 3.0*Gv[Gid].N3 - 6.0*Gv[KM1].N3 + Gv[KM2].N3);
        dGv[4] = (1.0 / 6.0)*(2.0*Gv[KP1].N4 + 3.0*Gv[Gid].N4 - 6.0*Gv[KM1].N4 + Gv[KM2].N4);
    }
    if (LOC_I == 0 && LOC_KM1 > 0 && LOC_KM1 <= 52 && (LOC_KM2 == 100))
    {
        dGv[0] = (1.0 / 6.0)*(-2.0*Gv[KM1].N0 - 3.0*Gv[Gid].N0 + 6.0*Gv[KP1].N0 - Gv[KP2].N0);
        dGv[1] = (1.0 / 6.0)*(-2.0*Gv[KM1].N1 - 3.0*Gv[Gid].N1 + 6.0*Gv[KP1].N1 - Gv[KP2].N1);
        dGv[2] = (1.0 / 6.0)*(-2.0*Gv[KM1].N2 - 3.0*Gv[Gid].N2 + 6.0*Gv[KP1].N2 - Gv[KP2].N2);
        dGv[3] = (1.0 / 6.0)*(-2.0*Gv[KM1].N3 - 3.0*Gv[Gid].N3 + 6.0*Gv[KP1].N3 - Gv[KP2].N3);
        dGv[4] = (1.0 / 6.0)*(-2.0*Gv[KM1].N4 - 3.0*Gv[Gid].N4 + 6.0*Gv[KP1].N4 - Gv[KP2].N4);
        
    }
    if (LOC_I == 0 && LOC_KM2 != 100 && LOC_KM1 != 100 && LOC_KP2 != 100 && LOC_KP1 != 100)
    {
        dGv[0] = (1.0 / 12.0)*(8.0*(Gv[KP1].N0 - Gv[KM1].N0) - (Gv[KP2].N0 - Gv[KM2].N0));
        dGv[1] = (1.0 / 12.0)*(8.0*(Gv[KP1].N1 - Gv[KM1].N1) - (Gv[KP2].N1 - Gv[KM2].N1));
        dGv[2] = (1.0 / 12.0)*(8.0*(Gv[KP1].N2 - Gv[KM1].N2) - (Gv[KP2].N2 - Gv[KM2].N2));
        dGv[3] = (1.0 / 12.0)*(8.0*(Gv[KP1].N3 - Gv[KM1].N3) - (Gv[KP2].N3 - Gv[KM2].N3));
        dGv[4] = (1.0 / 12.0)*(8.0*(Gv[KP1].N4 - Gv[KM1].N4) - (Gv[KP2].N4 - Gv[KM2].N4));
    }
    

    if (fabs(eigen_Qip[4]) >= fabs(eigen_Qim[4]))
    {
        alpha_u_ip = fabs(eigen_Qip[4]);
    }
    else if (fabs(eigen_Qim[4]) >= fabs(eigen_Qip[4]))
    {
        alpha_u_ip = fabs(eigen_Qim[4]);
    }

    if (fabs(eigen_Qjp[4]) >= fabs(eigen_Qjm[4]))
    {
        alpha_v_jp = fabs(eigen_Qjp[4]);
    }
    else if (fabs(eigen_Qjm[4]) >= fabs(eigen_Qjp[4]))
    {
        alpha_v_jp = fabs(eigen_Qjm[4]);
    }

    if (fabs(eigen_Qkp[4]) >= fabs(eigen_Qkm[4]))
    {
        alpha_w_kp = fabs(eigen_Qkp[4]);
    }
    else if (fabs(eigen_Qkm[4]) >= fabs(eigen_Qkp[4]))
    {
        alpha_w_kp = fabs(eigen_Qkm[4]);
    }

    if (fabs(eigen_Qinp[4]) >= fabs(eigen_Qinm[4]))
    {
        alpha_u_im = fabs(eigen_Qinp[4]);
    }
    else if (fabs(eigen_Qinm[4]) >= fabs(eigen_Qinp[4]))
    {
        alpha_u_im = fabs(eigen_Qinm[4]);
    }

    if (fabs(eigen_Qjnp[4]) >= fabs(eigen_Qjnm[4]))
    {
        alpha_v_jm = fabs(eigen_Qjnp[4]);
    }
    else if (fabs(eigen_Qjnm[4]) >= fabs(eigen_Qjnp[4]))
    {
        alpha_v_jm = fabs(eigen_Qjnm[4]);
    }

    if (fabs(eigen_Qknp[4]) >= fabs(eigen_Qknm[4]))
    {
        alpha_w_km = fabs(eigen_Qknp[4]);
    }
    else if (fabs(eigen_Qknm[4]) >= fabs(eigen_Qknp[4]))
    {
        alpha_w_km = fabs(eigen_Qknm[4]);
    }

    F_ip[0] = 0.5*(F_ip_pos_comp[0] + F_ip_neg_comp[0] - alpha_u_ip * (Qi_iplus_half_pos[0] - Qi_iplus_half_neg[0]));
    F_ip[1] = 0.5*(F_ip_pos_comp[1] + F_ip_neg_comp[1] - alpha_u_ip * (Qi_iplus_half_pos[1] - Qi_iplus_half_neg[1]));
    F_ip[2] = 0.5*(F_ip_pos_comp[2] + F_ip_neg_comp[2] - alpha_u_ip * (Qi_iplus_half_pos[2] - Qi_iplus_half_neg[2]));
    F_ip[3] = 0.5*(F_ip_pos_comp[3] + F_ip_neg_comp[3] - alpha_u_ip * (Qi_iplus_half_pos[3] - Qi_iplus_half_neg[3]));
    F_ip[4] = 0.5*(F_ip_pos_comp[4] + F_ip_neg_comp[4] - alpha_u_ip * (Qi_iplus_half_pos[4] - Qi_iplus_half_neg[4]));

    F_im[0] = 0.5*(F_im_pos_comp[0] + F_im_neg_comp[0] - alpha_u_im * (Qi_iminus_half_pos[0] - Qi_iminus_half_neg[0]));
    F_im[1] = 0.5*(F_im_pos_comp[1] + F_im_neg_comp[1] - alpha_u_im * (Qi_iminus_half_pos[1] - Qi_iminus_half_neg[1]));
    F_im[2] = 0.5*(F_im_pos_comp[2] + F_im_neg_comp[2] - alpha_u_im * (Qi_iminus_half_pos[2] - Qi_iminus_half_neg[2]));
    F_im[3] = 0.5*(F_im_pos_comp[3] + F_im_neg_comp[3] - alpha_u_im * (Qi_iminus_half_pos[3] - Qi_iminus_half_neg[3]));
    F_im[4] = 0.5*(F_im_pos_comp[4] + F_im_neg_comp[4] - alpha_u_im * (Qi_iminus_half_pos[4] - Qi_iminus_half_neg[4]));

    E_jp[0] = 0.5*(E_ip_pos_comp[0] + E_ip_neg_comp[0] - alpha_v_jp * (Qj_iplus_half_pos[0] - Qj_iplus_half_neg[0]));
    E_jp[1] = 0.5*(E_ip_pos_comp[1] + E_ip_neg_comp[1] - alpha_v_jp * (Qj_iplus_half_pos[1] - Qj_iplus_half_neg[1]));
    E_jp[2] = 0.5*(E_ip_pos_comp[2] + E_ip_neg_comp[2] - alpha_v_jp * (Qj_iplus_half_pos[2] - Qj_iplus_half_neg[2]));
    E_jp[3] = 0.5*(E_ip_pos_comp[3] + E_ip_neg_comp[3] - alpha_v_jp * (Qj_iplus_half_pos[3] - Qj_iplus_half_neg[3]));
    E_jp[4] = 0.5*(E_ip_pos_comp[4] + E_ip_neg_comp[4] - alpha_v_jp * (Qj_iplus_half_pos[4] - Qj_iplus_half_neg[4]));

    E_jm[0] = 0.5*(E_im_pos_comp[0] + E_im_neg_comp[0] - alpha_v_jm * (Qj_iminus_half_pos[0] - Qj_iminus_half_neg[0]));
    E_jm[1] = 0.5*(E_im_pos_comp[1] + E_im_neg_comp[1] - alpha_v_jm * (Qj_iminus_half_pos[1] - Qj_iminus_half_neg[1]));
    E_jm[2] = 0.5*(E_im_pos_comp[2] + E_im_neg_comp[2] - alpha_v_jm * (Qj_iminus_half_pos[2] - Qj_iminus_half_neg[2]));
    E_jm[3] = 0.5*(E_im_pos_comp[3] + E_im_neg_comp[3] - alpha_v_jm * (Qj_iminus_half_pos[3] - Qj_iminus_half_neg[3]));
    E_jm[4] = 0.5*(E_im_pos_comp[4] + E_im_neg_comp[4] - alpha_v_jm * (Qj_iminus_half_pos[4] - Qj_iminus_half_neg[4]));

    G_kp[0] = 0.5*(G_ip_pos_comp[0] + G_ip_neg_comp[0] - alpha_w_kp * (Qk_iplus_half_pos[0] - Qk_iplus_half_neg[0]));
    G_kp[1] = 0.5*(G_ip_pos_comp[1] + G_ip_neg_comp[1] - alpha_w_kp * (Qk_iplus_half_pos[1] - Qk_iplus_half_neg[1]));
    G_kp[2] = 0.5*(G_ip_pos_comp[2] + G_ip_neg_comp[2] - alpha_w_kp * (Qk_iplus_half_pos[2] - Qk_iplus_half_neg[2]));
    G_kp[3] = 0.5*(G_ip_pos_comp[3] + G_ip_neg_comp[3] - alpha_w_kp * (Qk_iplus_half_pos[3] - Qk_iplus_half_neg[3]));
    G_kp[4] = 0.5*(G_ip_pos_comp[4] + G_ip_neg_comp[4] - alpha_w_kp * (Qk_iplus_half_pos[4] - Qk_iplus_half_neg[4]));

    G_km[0] = 0.5*(G_im_pos_comp[0] + G_im_neg_comp[0] - alpha_w_km * (Qk_iminus_half_pos[0] - Qk_iminus_half_neg[0]));
    G_km[1] = 0.5*(G_im_pos_comp[1] + G_im_neg_comp[1] - alpha_w_km * (Qk_iminus_half_pos[1] - Qk_iminus_half_neg[1]));
    G_km[2] = 0.5*(G_im_pos_comp[2] + G_im_neg_comp[2] - alpha_w_km * (Qk_iminus_half_pos[2] - Qk_iminus_half_neg[2]));
    G_km[3] = 0.5*(G_im_pos_comp[3] + G_im_neg_comp[3] - alpha_w_km * (Qk_iminus_half_pos[3] - Qk_iminus_half_neg[3]));
    G_km[4] = 0.5*(G_im_pos_comp[4] + G_im_neg_comp[4] - alpha_w_km * (Qk_iminus_half_pos[4] - Qk_iminus_half_neg[4]));

    
    /****************************************/
    d2F_d2z_ip[0] = (1.0 / 48.0)*(-5.0*F[IM2].N0 + 39.0*F[IM1].N0 - 34.0*F[Gid].N0 - 34.0*F[IP1].N0 + 39.0*F[IP2].N0 - 5.0*F[IP3].N0);
    d2F_d2z_ip[1] = (1.0 / 48.0)*(-5.0*F[IM2].N1 + 39.0*F[IM1].N1 - 34.0*F[Gid].N1 - 34.0*F[IP1].N1 + 39.0*F[IP2].N1 - 5.0*F[IP3].N1);
    d2F_d2z_ip[2] = (1.0 / 48.0)*(-5.0*F[IM2].N2 + 39.0*F[IM1].N2 - 34.0*F[Gid].N2 - 34.0*F[IP1].N2 + 39.0*F[IP2].N2 - 5.0*F[IP3].N2);
    d2F_d2z_ip[3] = (1.0 / 48.0)*(-5.0*F[IM2].N3 + 39.0*F[IM1].N3 - 34.0*F[Gid].N3 - 34.0*F[IP1].N3 + 39.0*F[IP2].N3 - 5.0*F[IP3].N3);
    d2F_d2z_ip[4] = (1.0 / 48.0)*(-5.0*F[IM2].N4 + 39.0*F[IM1].N4 - 34.0*F[Gid].N4 - 34.0*F[IP1].N4 + 39.0*F[IP2].N4 - 5.0*F[IP3].N4);
    
    d2F_d2z_im[0] = (1.0 / 48.0)*(-5.0*F[IM3].N0 + 39.0*F[IM2].N0 - 34.0*F[IM1].N0 - 34.0*F[Gid].N0 + 39.0*F[IP1].N0 - 5.0*F[IP2].N0);
    d2F_d2z_im[1] = (1.0 / 48.0)*(-5.0*F[IM3].N1 + 39.0*F[IM2].N1 - 34.0*F[IM1].N1 - 34.0*F[Gid].N1 + 39.0*F[IP1].N1 - 5.0*F[IP2].N1);
    d2F_d2z_im[2] = (1.0 / 48.0)*(-5.0*F[IM3].N2 + 39.0*F[IM2].N2 - 34.0*F[IM1].N2 - 34.0*F[Gid].N2 + 39.0*F[IP1].N2 - 5.0*F[IP2].N2);
    d2F_d2z_im[3] = (1.0 / 48.0)*(-5.0*F[IM3].N3 + 39.0*F[IM2].N3 - 34.0*F[IM1].N3 - 34.0*F[Gid].N3 + 39.0*F[IP1].N3 - 5.0*F[IP2].N3);
    d2F_d2z_im[4] = (1.0 / 48.0)*(-5.0*F[IM3].N4 + 39.0*F[IM2].N4 - 34.0*F[IM1].N4 - 34.0*F[Gid].N4 + 39.0*F[IP1].N4 - 5.0*F[IP2].N4);

    d2E_d2e_ip[0] = (1.0 / 48.0)*(-5.0*E[JM2].N0 + 39.0*E[JM1].N0 - 34.0*E[Gid].N0 - 34.0*E[JP1].N0 + 39.0*E[JP2].N0 - 5.0*E[JP3].N0);
    d2E_d2e_ip[1] = (1.0 / 48.0)*(-5.0*E[JM2].N1 + 39.0*E[JM1].N1 - 34.0*E[Gid].N1 - 34.0*E[JP1].N1 + 39.0*E[JP2].N1 - 5.0*E[JP3].N1);
    d2E_d2e_ip[2] = (1.0 / 48.0)*(-5.0*E[JM2].N2 + 39.0*E[JM1].N2 - 34.0*E[Gid].N2 - 34.0*E[JP1].N2 + 39.0*E[JP2].N2 - 5.0*E[JP3].N2);
    d2E_d2e_ip[3] = (1.0 / 48.0)*(-5.0*E[JM2].N3 + 39.0*E[JM1].N3 - 34.0*E[Gid].N3 - 34.0*E[JP1].N3 + 39.0*E[JP2].N3 - 5.0*E[JP3].N3);
    d2E_d2e_ip[4] = (1.0 / 48.0)*(-5.0*E[JM2].N4 + 39.0*E[JM1].N4 - 34.0*E[Gid].N4 - 34.0*E[JP1].N4 + 39.0*E[JP2].N4 - 5.0*E[JP3].N4);
    
    d2E_d2e_im[0] = (1.0 / 48.0)*(-5.0*E[JM3].N0 + 39.0*E[JM2].N0 - 34.0*E[JM1].N0 - 34.0*E[Gid].N0 + 39.0*E[JP1].N0 - 5.0*E[JP2].N0);
    d2E_d2e_im[1] = (1.0 / 48.0)*(-5.0*E[JM3].N1 + 39.0*E[JM2].N1 - 34.0*E[JM1].N1 - 34.0*E[Gid].N1 + 39.0*E[JP1].N1 - 5.0*E[JP2].N1);
    d2E_d2e_im[2] = (1.0 / 48.0)*(-5.0*E[JM3].N2 + 39.0*E[JM2].N2 - 34.0*E[JM1].N2 - 34.0*E[Gid].N2 + 39.0*E[JP1].N2 - 5.0*E[JP2].N2);
    d2E_d2e_im[3] = (1.0 / 48.0)*(-5.0*E[JM3].N3 + 39.0*E[JM2].N3 - 34.0*E[JM1].N3 - 34.0*E[Gid].N3 + 39.0*E[JP1].N3 - 5.0*E[JP2].N3);
    d2E_d2e_im[4] = (1.0 / 48.0)*(-5.0*E[JM3].N4 + 39.0*E[JM2].N4 - 34.0*E[JM1].N4 - 34.0*E[Gid].N4 + 39.0*E[JP1].N4 - 5.0*E[JP2].N4);

    d2G_d2x_ip[0] = (1.0 / 48.0)*(-5.0*G[KM2].N0 + 39.0*G[KM1].N0 - 34.0*G[Gid].N0 - 34.0*G[KP1].N0 + 39.0*G[KP2].N0 - 5.0*G[KP3].N0);
    d2G_d2x_ip[1] = (1.0 / 48.0)*(-5.0*G[KM2].N1 + 39.0*G[KM1].N1 - 34.0*G[Gid].N1 - 34.0*G[KP1].N1 + 39.0*G[KP2].N1 - 5.0*G[KP3].N1);
    d2G_d2x_ip[2] = (1.0 / 48.0)*(-5.0*G[KM2].N2 + 39.0*G[KM1].N2 - 34.0*G[Gid].N2 - 34.0*G[KP1].N2 + 39.0*G[KP2].N2 - 5.0*G[KP3].N2);
    d2G_d2x_ip[3] = (1.0 / 48.0)*(-5.0*G[KM2].N3 + 39.0*G[KM1].N3 - 34.0*G[Gid].N3 - 34.0*G[KP1].N3 + 39.0*G[KP2].N3 - 5.0*G[KP3].N3);
    d2G_d2x_ip[4] = (1.0 / 48.0)*(-5.0*G[KM2].N4 + 39.0*G[KM1].N4 - 34.0*G[Gid].N4 - 34.0*G[KP1].N4 + 39.0*G[KP2].N4 - 5.0*G[KP3].N4);
    
    d2G_d2x_im[0] = (1.0 / 48.0)*(-5.0*G[KM3].N0 + 39.0*G[KM2].N0 - 34.0*G[KM1].N0 - 34.0*G[Gid].N0 + 39.0*G[KP1].N0 - 5.0*G[KP2].N0);
    d2G_d2x_im[1] = (1.0 / 48.0)*(-5.0*G[KM3].N1 + 39.0*G[KM2].N1 - 34.0*G[KM1].N1 - 34.0*G[Gid].N1 + 39.0*G[KP1].N1 - 5.0*G[KP2].N1);
    d2G_d2x_im[2] = (1.0 / 48.0)*(-5.0*G[KM3].N2 + 39.0*G[KM2].N2 - 34.0*G[KM1].N2 - 34.0*G[Gid].N2 + 39.0*G[KP1].N2 - 5.0*G[KP2].N2);
    d2G_d2x_im[3] = (1.0 / 48.0)*(-5.0*G[KM3].N3 + 39.0*G[KM2].N3 - 34.0*G[KM1].N3 - 34.0*G[Gid].N3 + 39.0*G[KP1].N3 - 5.0*G[KP2].N3);
    d2G_d2x_im[4] = (1.0 / 48.0)*(-5.0*G[KM3].N4 + 39.0*G[KM2].N4 - 34.0*G[KM1].N4 - 34.0*G[Gid].N4 + 39.0*G[KP1].N4 - 5.0*G[KP2].N4);

    d4F_d4z_ip[0] = (1.0 / 2.0)*(F[IM2].N0 - 3.0*F[IM1].N0 + 2.0*F[Gid].N0 + 2.0*F[IP1].N0 - 3.0*F[IP2].N0 + F[IP3].N0);
    d4F_d4z_ip[1] = (1.0 / 2.0)*(F[IM2].N1 - 3.0*F[IM1].N1 + 2.0*F[Gid].N1 + 2.0*F[IP1].N1 - 3.0*F[IP2].N1 + F[IP3].N1);
    d4F_d4z_ip[2] = (1.0 / 2.0)*(F[IM2].N2 - 3.0*F[IM1].N2 + 2.0*F[Gid].N2 + 2.0*F[IP1].N2 - 3.0*F[IP2].N2 + F[IP3].N2);
    d4F_d4z_ip[3] = (1.0 / 2.0)*(F[IM2].N3 - 3.0*F[IM1].N3 + 2.0*F[Gid].N3 + 2.0*F[IP1].N3 - 3.0*F[IP2].N3 + F[IP3].N3);
    d4F_d4z_ip[4] = (1.0 / 2.0)*(F[IM2].N4 - 3.0*F[IM1].N4 + 2.0*F[Gid].N4 + 2.0*F[IP1].N4 - 3.0*F[IP2].N4 + F[IP3].N4);
    
    d4F_d4z_im[0] = (1.0 / 2.0)*(F[IM3].N0 - 3.0*F[IM2].N0 + 2.0*F[IM1].N0 + 2.0*F[Gid].N0 - 3.0*F[IP1].N0 + F[IP2].N0);
    d4F_d4z_im[1] = (1.0 / 2.0)*(F[IM3].N1 - 3.0*F[IM2].N1 + 2.0*F[IM1].N1 + 2.0*F[Gid].N1 - 3.0*F[IP1].N1 + F[IP2].N1);
    d4F_d4z_im[2] = (1.0 / 2.0)*(F[IM3].N2 - 3.0*F[IM2].N2 + 2.0*F[IM1].N2 + 2.0*F[Gid].N2 - 3.0*F[IP1].N2 + F[IP2].N2);
    d4F_d4z_im[3] = (1.0 / 2.0)*(F[IM3].N3 - 3.0*F[IM2].N3 + 2.0*F[IM1].N3 + 2.0*F[Gid].N3 - 3.0*F[IP1].N3 + F[IP2].N3);
    d4F_d4z_im[4] = (1.0 / 2.0)*(F[IM3].N4 - 3.0*F[IM2].N4 + 2.0*F[IM1].N4 + 2.0*F[Gid].N4 - 3.0*F[IP1].N4 + F[IP2].N4);

    d4E_d4e_ip[0] = (1.0 / 2.0)*(E[JM2].N0 - 3.0*E[JM1].N0 + 2.0*E[Gid].N0 + 2.0*E[JP1].N0 - 3.0*E[JP2].N0 + E[JP3].N0);
    d4E_d4e_ip[1] = (1.0 / 2.0)*(E[JM2].N1 - 3.0*E[JM1].N1 + 2.0*E[Gid].N1 + 2.0*E[JP1].N1 - 3.0*E[JP2].N1 + E[JP3].N1);
    d4E_d4e_ip[2] = (1.0 / 2.0)*(E[JM2].N2 - 3.0*E[JM1].N2 + 2.0*E[Gid].N2 + 2.0*E[JP1].N2 - 3.0*E[JP2].N2 + E[JP3].N2);
    d4E_d4e_ip[3] = (1.0 / 2.0)*(E[JM2].N3 - 3.0*E[JM1].N3 + 2.0*E[Gid].N3 + 2.0*E[JP1].N3 - 3.0*E[JP2].N3 + E[JP3].N3);
    d4E_d4e_ip[4] = (1.0 / 2.0)*(E[JM2].N4 - 3.0*E[JM1].N4 + 2.0*E[Gid].N4 + 2.0*E[JP1].N4 - 3.0*E[JP2].N4 + E[JP3].N4);
    
    d4E_d4e_im[0] = (1.0 / 2.0)*(E[JM3].N0 - 3.0*E[JM2].N0 + 2.0*E[JM1].N0 + 2.0*E[Gid].N0 - 3.0*E[JP1].N0 + E[JP2].N0);
    d4E_d4e_im[1] = (1.0 / 2.0)*(E[JM3].N1 - 3.0*E[JM2].N1 + 2.0*E[JM1].N1 + 2.0*E[Gid].N1 - 3.0*E[JP1].N1 + E[JP2].N1);
    d4E_d4e_im[2] = (1.0 / 2.0)*(E[JM3].N2 - 3.0*E[JM2].N2 + 2.0*E[JM1].N2 + 2.0*E[Gid].N2 - 3.0*E[JP1].N2 + E[JP2].N2);
    d4E_d4e_im[3] = (1.0 / 2.0)*(E[JM3].N3 - 3.0*E[JM2].N3 + 2.0*E[JM1].N3 + 2.0*E[Gid].N3 - 3.0*E[JP1].N3 + E[JP2].N3);
    d4E_d4e_im[4] = (1.0 / 2.0)*(E[JM3].N4 - 3.0*E[JM2].N4 + 2.0*E[JM1].N4 + 2.0*E[Gid].N4 - 3.0*E[JP1].N4 + E[JP2].N4);

    d4G_d4x_ip[0] = (1.0 / 2.0)*(G[KM2].N0 - 3.0*G[KM1].N0 + 2.0*G[Gid].N0 + 2.0*G[KP1].N0 - 3.0*G[KP2].N0 + G[KP3].N0);
    d4G_d4x_ip[1] = (1.0 / 2.0)*(G[KM2].N1 - 3.0*G[KM1].N1 + 2.0*G[Gid].N1 + 2.0*G[KP1].N1 - 3.0*G[KP2].N1 + G[KP3].N1);
    d4G_d4x_ip[2] = (1.0 / 2.0)*(G[KM2].N2 - 3.0*G[KM1].N2 + 2.0*G[Gid].N2 + 2.0*G[KP1].N2 - 3.0*G[KP2].N2 + G[KP3].N2);
    d4G_d4x_ip[3] = (1.0 / 2.0)*(G[KM2].N3 - 3.0*G[KM1].N3 + 2.0*G[Gid].N3 + 2.0*G[KP1].N3 - 3.0*G[KP2].N3 + G[KP3].N3);
    d4G_d4x_ip[4] = (1.0 / 2.0)*(G[KM2].N4 - 3.0*G[KM1].N4 + 2.0*G[Gid].N4 + 2.0*G[KP1].N4 - 3.0*G[KP2].N4 + G[KP3].N4);
    
    d4G_d4x_im[0] = (1.0 / 2.0)*(G[KM3].N0 - 3.0*G[KM2].N0 + 2.0*G[KM1].N0 + 2.0*G[Gid].N0 - 3.0*G[KP1].N0 + G[KP2].N0);
    d4G_d4x_im[1] = (1.0 / 2.0)*(G[KM3].N1 - 3.0*G[KM2].N1 + 2.0*G[KM1].N1 + 2.0*G[Gid].N1 - 3.0*G[KP1].N1 + G[KP2].N1);
    d4G_d4x_im[2] = (1.0 / 2.0)*(G[KM3].N2 - 3.0*G[KM2].N2 + 2.0*G[KM1].N2 + 2.0*G[Gid].N2 - 3.0*G[KP1].N2 + G[KP2].N2);
    d4G_d4x_im[3] = (1.0 / 2.0)*(G[KM3].N3 - 3.0*G[KM2].N3 + 2.0*G[KM1].N3 + 2.0*G[Gid].N3 - 3.0*G[KP1].N3 + G[KP2].N3);
    d4G_d4x_im[4] = (1.0 / 2.0)*(G[KM3].N4 - 3.0*G[KM2].N4 + 2.0*G[KM1].N4 + 2.0*G[Gid].N4 - 3.0*G[KP1].N4 + G[KP2].N4);


    F_ip_h[0] = F_ip[0] - (1.0 / 24.0)*d2F_d2z_ip[0] + (7.0 / 5760.0)*d4F_d4z_ip[0];
    F_ip_h[1] = F_ip[1] - (1.0 / 24.0)*d2F_d2z_ip[1] + (7.0 / 5760.0)*d4F_d4z_ip[1];
    F_ip_h[2] = F_ip[2] - (1.0 / 24.0)*d2F_d2z_ip[2] + (7.0 / 5760.0)*d4F_d4z_ip[2];
    F_ip_h[3] = F_ip[3] - (1.0 / 24.0)*d2F_d2z_ip[3] + (7.0 / 5760.0)*d4F_d4z_ip[3];
    F_ip_h[4] = F_ip[4] - (1.0 / 24.0)*d2F_d2z_ip[4] + (7.0 / 5760.0)*d4F_d4z_ip[4];
    
    F_im_h[0] = F_im[0] - (1.0 / 24.0)*d2F_d2z_im[0] + (7.0 / 5760.0)*d4F_d4z_im[0];
    F_im_h[1] = F_im[1] - (1.0 / 24.0)*d2F_d2z_im[1] + (7.0 / 5760.0)*d4F_d4z_im[1];
    F_im_h[2] = F_im[2] - (1.0 / 24.0)*d2F_d2z_im[2] + (7.0 / 5760.0)*d4F_d4z_im[2];
    F_im_h[3] = F_im[3] - (1.0 / 24.0)*d2F_d2z_im[3] + (7.0 / 5760.0)*d4F_d4z_im[3];
    F_im_h[4] = F_im[4] - (1.0 / 24.0)*d2F_d2z_im[4] + (7.0 / 5760.0)*d4F_d4z_im[4];

    E_ip_h[0] = E_jp[0] - (1.0 / 24.0)*d2E_d2e_ip[0] + (7.0 / 5760.0)*d4E_d4e_ip[0];
    E_ip_h[1] = E_jp[1] - (1.0 / 24.0)*d2E_d2e_ip[1] + (7.0 / 5760.0)*d4E_d4e_ip[1];
    E_ip_h[2] = E_jp[2] - (1.0 / 24.0)*d2E_d2e_ip[2] + (7.0 / 5760.0)*d4E_d4e_ip[2];
    E_ip_h[3] = E_jp[3] - (1.0 / 24.0)*d2E_d2e_ip[3] + (7.0 / 5760.0)*d4E_d4e_ip[3];
    E_ip_h[4] = E_jp[4] - (1.0 / 24.0)*d2E_d2e_ip[4] + (7.0 / 5760.0)*d4E_d4e_ip[4];
    
    E_im_h[0] = E_jm[0] - (1.0 / 24.0)*d2E_d2e_im[0] + (7.0 / 5760.0)*d4E_d4e_im[0];
    E_im_h[1] = E_jm[1] - (1.0 / 24.0)*d2E_d2e_im[1] + (7.0 / 5760.0)*d4E_d4e_im[1];
    E_im_h[2] = E_jm[2] - (1.0 / 24.0)*d2E_d2e_im[2] + (7.0 / 5760.0)*d4E_d4e_im[2];
    E_im_h[3] = E_jm[3] - (1.0 / 24.0)*d2E_d2e_im[3] + (7.0 / 5760.0)*d4E_d4e_im[3];
    E_im_h[4] = E_jm[4] - (1.0 / 24.0)*d2E_d2e_im[4] + (7.0 / 5760.0)*d4E_d4e_im[4];

    G_ip_h[0] = G_kp[0] - (1.0 / 24.0)*d2G_d2x_ip[0] + (7.0 / 5760.0)*d4G_d4x_ip[0];
    G_ip_h[1] = G_kp[1] - (1.0 / 24.0)*d2G_d2x_ip[1] + (7.0 / 5760.0)*d4G_d4x_ip[1];
    G_ip_h[2] = G_kp[2] - (1.0 / 24.0)*d2G_d2x_ip[2] + (7.0 / 5760.0)*d4G_d4x_ip[2];
    G_ip_h[3] = G_kp[3] - (1.0 / 24.0)*d2G_d2x_ip[3] + (7.0 / 5760.0)*d4G_d4x_ip[3];
    G_ip_h[4] = G_kp[4] - (1.0 / 24.0)*d2G_d2x_ip[4] + (7.0 / 5760.0)*d4G_d4x_ip[4];
    
    G_im_h[0] = G_km[0] - (1.0 / 24.0)*d2G_d2x_im[0] + (7.0 / 5760.0)*d4G_d4x_im[0];
    G_im_h[1] = G_km[1] - (1.0 / 24.0)*d2G_d2x_im[1] + (7.0 / 5760.0)*d4G_d4x_im[1];
    G_im_h[2] = G_km[2] - (1.0 / 24.0)*d2G_d2x_im[2] + (7.0 / 5760.0)*d4G_d4x_im[2];
    G_im_h[3] = G_km[3] - (1.0 / 24.0)*d2G_d2x_im[3] + (7.0 / 5760.0)*d4G_d4x_im[3];
    G_im_h[4] = G_km[4] - (1.0 / 24.0)*d2G_d2x_im[4] + (7.0 / 5760.0)*d4G_d4x_im[4];

    dF[0] = F_ip_h[0] - F_im_h[0];
    dF[1] = F_ip_h[1] - F_im_h[1];
    dF[2] = F_ip_h[2] - F_im_h[2];
    dF[3] = F_ip_h[3] - F_im_h[3];
    dF[4] = F_ip_h[4] - F_im_h[4];
    
    dE[0] = E_ip_h[0] - E_im_h[0];
    dE[1] = E_ip_h[1] - E_im_h[1];
    dE[2] = E_ip_h[2] - E_im_h[2];
    dE[3] = E_ip_h[3] - E_im_h[3];
    dE[4] = E_ip_h[4] - E_im_h[4];
    
    dG[0] = G_ip_h[0] - G_im_h[0];
    dG[1] = G_ip_h[1] - G_im_h[1];
    dG[2] = G_ip_h[2] - G_im_h[2];
    dG[3] = G_ip_h[3] - G_im_h[3];
    dG[4] = G_ip_h[4] - G_im_h[4];

    L[0] = (-1.0)*(dF[0] + dE[0] + dG[0] + dFv[0] + dEv[0] + dGv[0]);
    L[1] = (-1.0)*(dF[1] + dE[1] + dG[1] + dFv[1] + dEv[1] + dGv[1]);
    L[2] = (-1.0)*(dF[2] + dE[2] + dG[2] + dFv[2] + dEv[2] + dGv[2]);
    L[3] = (-1.0)*(dF[3] + dE[3] + dG[3] + dFv[3] + dEv[3] + dGv[3]);
    L[4] = (-1.0)*(dF[4] + dE[4] + dG[4] + dFv[4] + dEv[4] + dGv[4]);
    
    /***********************3rd order RK Scheme***************************************************************************/
    if (j == 0)
    {
        U1[Gid].N0 = U0[Gid].N0 + del_t * L[0];
        U1[Gid].N1 = U0[Gid].N1 + del_t * L[1];
        U1[Gid].N2 = U0[Gid].N2 + del_t * L[2];
        U1[Gid].N3 = U0[Gid].N3 + del_t * L[3];
        U1[Gid].N4 = U0[Gid].N4 + del_t * L[4];

        final_U[0] = (1.0 / det[Gid])*U1[Gid].N0;
        final_U[1] = (1.0 / det[Gid])*U1[Gid].N1;
        final_U[2] = (1.0 / det[Gid])*U1[Gid].N2;
        final_U[3] = (1.0 / det[Gid])*U1[Gid].N3;
        final_U[4] = (1.0 / det[Gid])*U1[Gid].N4;
    }
    else if (j == 1)
    {
        U2[Gid].N0 = (3.0 / 4.0)*U0[Gid].N0 + (1.0 / 4.0)*U1[Gid].N0 + 0.25*del_t*L[0];
        U2[Gid].N1 = (3.0 / 4.0)*U0[Gid].N1 + (1.0 / 4.0)*U1[Gid].N1 + 0.25*del_t*L[1];
        U2[Gid].N2 = (3.0 / 4.0)*U0[Gid].N2 + (1.0 / 4.0)*U1[Gid].N2 + 0.25*del_t*L[2];
        U2[Gid].N3 = (3.0 / 4.0)*U0[Gid].N3 + (1.0 / 4.0)*U1[Gid].N3 + 0.25*del_t*L[3];
        U2[Gid].N4 = (3.0 / 4.0)*U0[Gid].N4 + (1.0 / 4.0)*U1[Gid].N4 + 0.25*del_t*L[4];

        final_U[0] = (1.0 / det[Gid])*U2[Gid].N0;
        final_U[1] = (1.0 / det[Gid])*U2[Gid].N1;
        final_U[2] = (1.0 / det[Gid])*U2[Gid].N2;
        final_U[3] = (1.0 / det[Gid])*U2[Gid].N3;
        final_U[4] = (1.0 / det[Gid])*U2[Gid].N4;
    }
    else if (j == 2)
    {
        U0[Gid].N0 = (U0[Gid].N0 / 3.0) + ((2.0 / 3.0)*U2[Gid].N0) + (2.0 / 3.0)*(del_t*L[0]);
        U0[Gid].N1 = (U0[Gid].N1 / 3.0) + ((2.0 / 3.0)*U2[Gid].N1) + (2.0 / 3.0)*(del_t*L[1]);
        U0[Gid].N2 = (U0[Gid].N2 / 3.0) + ((2.0 / 3.0)*U2[Gid].N2) + (2.0 / 3.0)*(del_t*L[2]);
        U0[Gid].N3 = (U0[Gid].N3 / 3.0) + ((2.0 / 3.0)*U2[Gid].N3) + (2.0 / 3.0)*(del_t*L[3]);
        U0[Gid].N4 = (U0[Gid].N4 / 3.0) + ((2.0 / 3.0)*U2[Gid].N4) + (2.0 / 3.0)*(del_t*L[4]);

        final_U[0] = (1.0 / det[Gid])*U0[Gid].N0;
        final_U[1] = (1.0 / det[Gid])*U0[Gid].N1;
        final_U[2] = (1.0 / det[Gid])*U0[Gid].N2;
        final_U[3] = (1.0 / det[Gid])*U0[Gid].N3;
        final_U[4] = (1.0 / det[Gid])*U0[Gid].N4;
    }

    FLOW[Gid].rho = final_U[0];
    FLOW[Gid].u = final_U[1] / final_U[0];
    FLOW[Gid].v = final_U[2] / final_U[0];
    FLOW[Gid].w = final_U[3] / final_U[0];
    FLOW[Gid].e = (final_U[4] / final_U[0]) - 0.5*(pow(FLOW[Gid].u, 2.0) + pow(FLOW[Gid].v, 2.0) + pow(FLOW[Gid].w, 2.0));
    FLOW[Gid].p = FLOW[Gid].e * FLOW[Gid].rho * 0.4;
    FLOW[Gid].t = (1.4*Mach*Mach)*(FLOW[Gid].p / FLOW[Gid].rho);
    FLOW[Gid].a = sqrt(1.4*FLOW[Gid].p / FLOW[Gid].rho);

}


