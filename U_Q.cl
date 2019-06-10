#pragma OPENCL EXTENSION cl_khr_fp64 : enable
typedef struct __attribute__ ((packed)) tag_struct{
    double u_I, v_I, w_I, p_I, t_I, rho_I;
    double u_IP1, v_IP1, w_IP1, p_IP1, t_IP1, rho_IP1;
    double u_JP1, v_JP1, w_JP1, p_JP1, t_JP1, rho_JP1;
    double u_KP1, v_KP1, w_KP1, p_KP1, t_KP1, rho_KP1;
    double u_IM1, v_IM1, w_IM1, p_IM1, t_IM1, rho_IM1;
    double u_JM1, v_JM1, w_JM1, p_JM1, t_JM1, rho_JM1;
    double u_KM1, v_KM1, w_KM1, p_KM1, t_KM1, rho_KM1;
} FLOW_VAR;

typedef struct __attribute__ ((packed)) tag_Solver_Column_matrix {
	double N0, N1, N2, N3, N4;
} Solver_Column_matrix;

typedef struct __attribute__ ((packed)) tag_flow_variables {
	double u, v, w, p, t, rho, mu, a, e;
} FLOW_VARIABLES;

typedef struct __attribute__ ((packed)) tag_grid {
	double x, y, z;
	int N0, N1, N2, N3, N4, N5;
	int ID, loc;
} GRID;

typedef struct __attribute__ ((packed)) tag_ROE_AVER {
	double roe_u_ip, roe_v_ip, roe_w_ip, roe_p_ip, roe_h_ip, roe_rho_ip;
	double roe_u_jp, roe_v_jp, roe_w_jp, roe_p_jp, roe_h_jp, roe_rho_jp;
	double roe_u_kp, roe_v_kp, roe_w_kp, roe_p_kp, roe_h_kp, roe_rho_kp;
	double roe_u_im, roe_v_im, roe_w_im, roe_p_im, roe_h_im, roe_rho_im;
	double roe_u_jm, roe_v_jm, roe_w_jm, roe_p_jm, roe_h_jm, roe_rho_jm;
	double roe_u_km, roe_v_km, roe_w_km, roe_p_km, roe_h_km, roe_rho_km;
} ROE_AVERAGING;

//__local double *lBuffer_flowI, __local FLOW_VAR *flow, __local double *DET_I, __local double *sqr_i, __local double *sqr_b1, __local double *sqr_b3, __local double *sqr_c0, __local double *sqr_c2, __local double *sqr_d4, __local double *sqr_d5
__kernel void U_Q(__global Solver_Column_matrix *U_C, __global Solver_Column_matrix *Q_C, \
__global const GRID *Grid_P, __global const double *det, __global FLOW_VARIABLES *FLOW, \
__global ROE_AVERAGING *ROE_AVER)
{
	uint Gid = get_global_id(0);
  uint Lid = get_local_id(0);

  __private double DET_I;
  __private double sqr_i, sqr_b1, sqr_b3, sqr_c0, sqr_c2, sqr_d4, sqr_d5;
  __private uint IP1, IM1, JP1, JM1, KP1, KM1;

	if (Gid == 1) {
      printf("Hello U_Q %d\n", __LINE__);
  }

  IP1 = Grid_P[Gid].N1;
  IM1 = Grid_P[Gid].N3;
  JP1 = Grid_P[Gid].N0;
  JM1 = Grid_P[Gid].N2;
  KP1 = Grid_P[Gid].N4;
  KM1 = Grid_P[Gid].N5;

  barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

  DET_I = det[Gid];
  sqr_i = sqrt(FLOW[Gid].rho);
	sqr_b1 = sqrt(FLOW[IP1].rho);
	sqr_b3 = sqrt(FLOW[IM1].rho);
	sqr_c0 = sqrt(FLOW[JP1].rho);
	sqr_c2 = sqrt(FLOW[JM1].rho);
	sqr_d4 = sqrt(FLOW[KP1].rho);
	sqr_d5 = sqrt(FLOW[KM1].rho);    

	ROE_AVER[Gid].roe_rho_ip = sqr_b1 * sqr_i;
	ROE_AVER[Gid].roe_u_ip = (sqr_b1*FLOW[IP1].u + sqr_i * FLOW[Gid].u) / (sqr_b1 + sqr_i);
	ROE_AVER[Gid].roe_v_ip = (sqr_b1*FLOW[IP1].v + sqr_i * FLOW[Gid].v) / (sqr_b1 + sqr_i);
	ROE_AVER[Gid].roe_w_ip = (sqr_b1*FLOW[IP1].w + sqr_i * FLOW[Gid].w) / (sqr_b1 + sqr_i);
	ROE_AVER[Gid].roe_h_ip = (sqr_b1*(((FLOW[IP1].p / 0.4) + 0.5*FLOW[IP1].rho * (FLOW[IP1].u * FLOW[IP1].u + FLOW[IP1].v * FLOW[IP1].v + FLOW[IP1].w * FLOW[IP1].w) + FLOW[IP1].p) / FLOW[IP1].rho) + \
		sqr_i*(((FLOW[Gid].p / 0.4) + 0.5*FLOW[Gid].rho * (FLOW[Gid].u * FLOW[Gid].u + FLOW[Gid].v * FLOW[Gid].v + FLOW[Gid].w * FLOW[Gid].w) + FLOW[Gid].p) / FLOW[Gid].rho)) / (sqr_b1 + sqr_i);

  ROE_AVER[Gid].roe_rho_im = sqr_b3 * sqr_i;
	ROE_AVER[Gid].roe_u_im = (sqr_b3*FLOW[IM1].u + sqr_i * FLOW[Gid].u) / (sqr_b3 + sqr_i);
	ROE_AVER[Gid].roe_v_im = (sqr_b3*FLOW[IM1].v + sqr_i * FLOW[Gid].v) / (sqr_b3 + sqr_i);
	ROE_AVER[Gid].roe_w_im = (sqr_b3*FLOW[IM1].w + sqr_i * FLOW[Gid].w) / (sqr_b3 + sqr_i);
	ROE_AVER[Gid].roe_h_im = (sqr_b3*(((FLOW[IM1].p / 0.4) + 0.5*FLOW[IM1].rho * (FLOW[IM1].u * FLOW[IM1].u + FLOW[IM1].v * FLOW[IM1].v + FLOW[IM1].w * FLOW[IM1].w) + FLOW[IM1].p) / FLOW[IM1].rho) + \
		sqr_i*(((FLOW[Gid].p / 0.4) + 0.5*FLOW[Gid].rho * (FLOW[Gid].u * FLOW[Gid].u + FLOW[Gid].v * FLOW[Gid].v + FLOW[Gid].w * FLOW[Gid].w) + FLOW[Gid].p) / FLOW[Gid].rho)) / (sqr_b3 + sqr_i);

  ROE_AVER[Gid].roe_rho_jp = sqr_c0 * sqr_i;
	ROE_AVER[Gid].roe_u_jp = (sqr_c0*FLOW[JP1].u + sqr_i * FLOW[Gid].u) / (sqr_c0 + sqr_i);
	ROE_AVER[Gid].roe_v_jp = (sqr_c0*FLOW[JP1].v + sqr_i * FLOW[Gid].v) / (sqr_c0 + sqr_i);
	ROE_AVER[Gid].roe_w_jp = (sqr_c0*FLOW[JP1].w + sqr_i * FLOW[Gid].w) / (sqr_c0 + sqr_i);
	ROE_AVER[Gid].roe_h_jp = (sqr_c0*(((FLOW[JP1].p / 0.4) + 0.5*FLOW[JP1].rho * (FLOW[JP1].u * FLOW[JP1].u + FLOW[JP1].v * FLOW[JP1].v + FLOW[JP1].w * FLOW[JP1].w) + FLOW[JP1].p) / FLOW[JP1].rho) + \
		sqr_i*(((FLOW[Gid].p / 0.4) + 0.5*FLOW[Gid].rho * (FLOW[Gid].u * FLOW[Gid].u + FLOW[Gid].v * FLOW[Gid].v + FLOW[Gid].w * FLOW[Gid].w) + FLOW[Gid].p) / FLOW[Gid].rho)) / (sqr_c0 + sqr_i);

  ROE_AVER[Gid].roe_rho_jm = sqr_c2 * sqr_i;
	ROE_AVER[Gid].roe_u_jm = (sqr_c2*FLOW[JM1].u + sqr_i * FLOW[Gid].u) / (sqr_c2 + sqr_i);
	ROE_AVER[Gid].roe_v_jm = (sqr_c2*FLOW[JM1].v + sqr_i * FLOW[Gid].v) / (sqr_c2 + sqr_i);
	ROE_AVER[Gid].roe_w_jm = (sqr_c2*FLOW[JM1].w + sqr_i * FLOW[Gid].w) / (sqr_c2 + sqr_i);
	ROE_AVER[Gid].roe_h_jm = (sqr_c2*(((FLOW[JM1].p / 0.4) + 0.5*FLOW[JM1].rho * (FLOW[JM1].u * FLOW[JM1].u + FLOW[JM1].v * FLOW[JM1].v + FLOW[JM1].w * FLOW[JM1].w) + FLOW[JM1].p) / FLOW[JM1].rho) + \
		sqr_i*(((FLOW[Gid].p / 0.4) + 0.5*FLOW[Gid].rho * (FLOW[Gid].u * FLOW[Gid].u + FLOW[Gid].v * FLOW[Gid].v + FLOW[Gid].w * FLOW[Gid].w) + FLOW[Gid].p) / FLOW[Gid].rho)) / (sqr_c2 + sqr_i);

  ROE_AVER[Gid].roe_rho_kp = sqr_d4 * sqr_i;
	ROE_AVER[Gid].roe_u_kp = (sqr_d4*FLOW[KP1].u + sqr_i * FLOW[Gid].u) / (sqr_d4 + sqr_i);
	ROE_AVER[Gid].roe_v_kp = (sqr_d4*FLOW[KP1].v + sqr_i * FLOW[Gid].v) / (sqr_d4 + sqr_i);
	ROE_AVER[Gid].roe_w_kp = (sqr_d4*FLOW[KP1].w + sqr_i * FLOW[Gid].w) / (sqr_d4 + sqr_i);
	ROE_AVER[Gid].roe_h_kp = (sqr_d4*(((FLOW[KP1].p / 0.4) + 0.5*FLOW[KP1].rho * (FLOW[KP1].u * FLOW[KP1].u + FLOW[KP1].v * FLOW[KP1].v + FLOW[KP1].w * FLOW[KP1].w) + FLOW[KP1].p) / FLOW[KP1].rho) + \
		sqr_i*(((FLOW[Gid].p / 0.4) + 0.5*FLOW[Gid].rho * (FLOW[Gid].u * FLOW[Gid].u + FLOW[Gid].v * FLOW[Gid].v + FLOW[Gid].w * FLOW[Gid].w) + FLOW[Gid].p) / FLOW[Gid].rho)) / (sqr_d4 + sqr_i);

  ROE_AVER[Gid].roe_rho_km = sqr_d5 * sqr_i;
	ROE_AVER[Gid].roe_u_km = (sqr_d5*FLOW[KM1].u + sqr_i * FLOW[Gid].u) / (sqr_d5 + sqr_i);
	ROE_AVER[Gid].roe_v_km = (sqr_d5*FLOW[KM1].v + sqr_i * FLOW[Gid].v) / (sqr_d5 + sqr_i);
	ROE_AVER[Gid].roe_w_km = (sqr_d5*FLOW[KM1].w + sqr_i * FLOW[Gid].w) / (sqr_d5 + sqr_i);
	ROE_AVER[Gid].roe_h_km = (sqr_d5*(((FLOW[KM1].p / 0.4) + 0.5*FLOW[KM1].rho * (FLOW[KM1].u * FLOW[KM1].u + FLOW[KM1].v * FLOW[KM1].v + FLOW[KM1].w * FLOW[KM1].w) + FLOW[KM1].p) / FLOW[KM1].rho) + \
		sqr_i*(((FLOW[Gid].p / 0.4) + 0.5*FLOW[Gid].rho * (FLOW[Gid].u * FLOW[Gid].u + FLOW[Gid].v * FLOW[Gid].v + FLOW[Gid].w * FLOW[Gid].w) + FLOW[Gid].p) / FLOW[Gid].rho)) / (sqr_d5 + sqr_i);


  U_C[Gid].N0 = DET_I*FLOW[Gid].rho;
  U_C[Gid].N1 = DET_I*FLOW[Gid].rho * FLOW[Gid].u;
  U_C[Gid].N2 = DET_I*FLOW[Gid].rho * FLOW[Gid].v;
  U_C[Gid].N3 = DET_I*FLOW[Gid].rho * FLOW[Gid].w;
  U_C[Gid].N4 = DET_I*((FLOW[Gid].p / 0.4) + 0.5*FLOW[Gid].rho * (FLOW[Gid].u * FLOW[Gid].u + FLOW[Gid].v * FLOW[Gid].v + FLOW[Gid].w * FLOW[Gid].w));

  Q_C[Gid].N0 = FLOW[Gid].rho;
  Q_C[Gid].N1 = FLOW[Gid].rho * FLOW[Gid].u;
  Q_C[Gid].N2 = FLOW[Gid].rho * FLOW[Gid].v;
  Q_C[Gid].N3 = FLOW[Gid].rho * FLOW[Gid].w;
  Q_C[Gid].N4 = ((FLOW[Gid].p / 0.4) + 0.5*FLOW[Gid].rho * (FLOW[Gid].u * FLOW[Gid].u + FLOW[Gid].v * FLOW[Gid].v + FLOW[Gid].w * FLOW[Gid].w));

    
}
