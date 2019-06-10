//#define __CL_ENABLE_EXCEPTIONS
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#include <iostream>
#include <cmath>
#include <string>
#include "functions.h"
#include <openmpi-x86_64/mpi.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <CL/cl.hpp>

using std::cout;
using std::cin;
using std::string;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::istringstream;

void solver()
{
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status status;

	/*******************OpenCL Initialisation*****************************/
	cl_uint Lsize = 1;
	cl_uint Gsize = (Lsize - (g_node%Lsize)) + g_node;

	cl_uint Weno_Gsize = (Lsize - (sd_node%Lsize)) + sd_node;

	cl::NDRange global(Gsize);
	cl::NDRange Weno_global(Weno_Gsize);
	cl::NDRange local(Lsize);

	cl::Event event;
	cl_int ret;
	cl_uint NumCPUDevices, NumGPUDevices, NumPlatforms;
	size_t A_size;
	std::vector<cl::Device> device_id, GPU_deviceID, CPU_deviceID;
	std::vector<cl::Platform> platforms_id, CPU_platformID, GPU_platformID; 
	cl::CommandQueue queue, queue_viscous, queue_weno;
	
	ret = cl::Platform::get(&platforms_id);
	if (ret != CL_SUCCESS) {
		std::cout << "Error in Obtaining Platform Id" << std::endl;
		exit(0);
	}
	NumPlatforms = platforms_id.size();
/*	
	for (i=0; i< NumPlatforms; i++) {
		ret = platforms_id[i].getDevices(CL_DEVICE_TYPE_CPU, &device_id);
		if (ret == CL_SUCCESS) {
			CPU_platformID.resize(platforms_id.size());
			CPU_deviceID.resize(device_id.size());
			for (j=0; j<platforms_id.size(); j++) {
				CPU_platformID[j] = platforms_id[j];
			}
			for (j=0; j<device_id.size(); j++) {
				CPU_deviceID[j] = device_id[j];
			}
			std::vector<cl::Device>().swap(device_id);
			std::vector<cl::Platform>().swap(platforms_id);
			i = 2;
		}
		else {
			ret = platforms_id[i].getDevices(CL_DEVICE_TYPE_GPU, &device_id);
			if (ret == CL_SUCCESS) {
				GPU_platformID.resize(platforms_id.size());
				GPU_deviceID.resize(device_id.size());
				for (j=0; j<platforms_id.size(); j++) {
					GPU_platformID[j] = platforms_id[j];
				}
				for (j=0; j<device_id.size(); j++) {
					GPU_deviceID[j] = device_id[j];
				}
				std::vector<cl::Device>().swap(device_id);
				std::vector<cl::Platform>().swap(platforms_id);
			}
		}
	}
*/
	ret = platforms_id[0].getDevices(CL_DEVICE_TYPE_CPU, &device_id);
	if (ret == CL_SUCCESS) {
		GPU_platformID.resize(platforms_id.size());
		GPU_deviceID.resize(device_id.size());
		CPU_platformID.resize(platforms_id.size());
		CPU_deviceID.resize(device_id.size());
		GPU_platformID[0] = platforms_id[0];
		CPU_platformID[0] = GPU_platformID[0];
		GPU_deviceID[0] = device_id[0];
		CPU_deviceID[0] = GPU_deviceID[0];
		
		std::vector<cl::Device>().swap(device_id);
		std::vector<cl::Platform>().swap(platforms_id);
	}

	cl_context_properties cps[] = {CL_CONTEXT_PLATFORM, (cl_context_properties)(CPU_platformID[0])(), 0};
   	cl::Context context( CL_DEVICE_TYPE_CPU, cps);	
	//cl::Context context_viscous( CL_DEVICE_TYPE_GPU, cps);	

	std::vector<cl::Device> devices;
	devices = context.getInfo<CL_CONTEXT_DEVICES>();

	queue = cl::CommandQueue(context, CPU_deviceID[0], CL_QUEUE_PROFILING_ENABLE, &ret);
	if (ret != CL_SUCCESS) {
		std::cout << "Error creating Commandqueue for context" << std::endl;
	}
	
	queue_viscous = cl::CommandQueue(context, CPU_deviceID[0], CL_QUEUE_PROFILING_ENABLE, &ret);
	if (ret != CL_SUCCESS) {
		std::cout << "Error creating Commandqueue for context" << std::endl;
	}

	queue_weno = cl::CommandQueue(context, CPU_deviceID[0], CL_QUEUE_PROFILING_ENABLE, &ret);
	if (ret != CL_SUCCESS) {
		std::cout << "Error creating Commandqueue for context" << std::endl;
	}
	/*********************************Read U_Q.cl Kernel file******************************************/
	std::ifstream sourceFile("U_Q.cl");
    if (!sourceFile.good()) {
		std::cerr << "Kernel file not found: " << filename << std::endl;
		exit(1);
	}
	std::string sourceCode(std::istreambuf_iterator<char>(sourceFile),(std::istreambuf_iterator<char>()));
    cl::Program::Sources source(1, std::make_pair(sourceCode.c_str(), sourceCode.length()+1));

	cl::Program program_U_Q;
	program_U_Q = cl::Program(context, source, &ret);
	const char options[] = "-cl-std=CL2.0";
	ret = program_U_Q.build(devices, options);
	cl::STRING_CLASS buildLog;
	program_U_Q.getBuildInfo(devices[0], CL_PROGRAM_BUILD_LOG, &buildLog);
	if (ret != CL_SUCCESS) {
		std::cout << "Failed to build the U_Q program: " << std::endl;
		std::cout << buildLog << std::endl;
		exit(1);
	}

	cl::Kernel kernel_U_Q(program_U_Q, "U_Q", &ret);
	if (ret != CL_SUCCESS) {
		std::cout << "Failed to create kernel for the U_Q program" << std::endl;
		exit(1);
	}
	sourceCode.erase();
	buildLog.erase();
	/**************************************************************************************************/
	
	/*********************************Read Viscous.cl Kernel file******************************************/
	std::ifstream sourceFile2("Viscous.cl");
    if (!sourceFile.good()) {
		std::cerr << "Kernel file not found: " << filename << std::endl;
		exit(1);
	}
	std::string sourceCode2(std::istreambuf_iterator<char>(sourceFile2),(std::istreambuf_iterator<char>()));
    cl::Program::Sources source2(1, std::make_pair(sourceCode2.c_str(), sourceCode2.length()+1));

	cl::Program program_viscous;
	program_viscous = cl::Program(context, source2, &ret);
	ret = program_viscous.build(devices, options);
	program_viscous.getBuildInfo(devices[0], CL_PROGRAM_BUILD_LOG, &buildLog);
	if (ret != CL_SUCCESS) {
		std::cout << "Failed to build the viscous program: " << std::endl;
		std::cout << buildLog << std::endl;
		exit(1);
	}

	cl::Kernel kernel_viscous(program_viscous, "Viscous", &ret);
	if (ret != CL_SUCCESS) {
		std::cout << "Failed to create kernel for the viscous program" << std::endl;
		exit(1);
	}
	/***************************************************************************************************/

	/*********************************Read Weno.cl Kernel file******************************************/
	std::ifstream sourceFile3("Weno.cl");
    if (!sourceFile.good()) {
		std::cerr << "Kernel file not found: " << filename << std::endl;
		exit(1);
	}
	std::string sourceCode3(std::istreambuf_iterator<char>(sourceFile3),(std::istreambuf_iterator<char>()));
    cl::Program::Sources source3(1, std::make_pair(sourceCode3.c_str(), sourceCode3.length()+1));

	cl::Program program_weno;
	program_weno = cl::Program(context, source3, &ret);
	ret = program_weno.build(devices, options);
	program_weno.getBuildInfo(devices[0], CL_PROGRAM_BUILD_LOG, &buildLog);
	if (ret != CL_SUCCESS) {
		std::cout << "Failed to build the weno.cl program: " << std::endl;
		std::cout << buildLog << std::endl;
		exit(1);
	}

	cl::Kernel kernel_weno(program_weno, "Weno", &ret);
	if (ret != CL_SUCCESS) {
		std::cout << "Failed to create kernel for the weno.cl program" << std::endl;
		exit(1);
	}
	/***********************************************************************************************/
	cl::Buffer buffer_U0, buffer_U1, buffer_U2, buffer_Q_C, buffer_Grid, buffer_Det, buffer_FlowVar, buffer_RoeVar;
	cl::Buffer buffer_tzz, buffer_tee, buffer_txx, buffer_tze, buffer_tzx, buffer_tex;
	cl::Buffer buffer_F_C, buffer_E_C, buffer_G_C, buffer_Fv_C, buffer_Ev_C, buffer_Gv_C;        
	cl::Buffer buffer_metric;

	cl::Buffer buffer_j, buffer_Mach, buffer_del_t, buffer_deter ,buffer_Q_I;
	cl::Buffer buffer_Q_IP1, buffer_Q_IP2, buffer_Q_IP3;
	cl::Buffer buffer_Q_IM1, buffer_Q_IM2, buffer_Q_IM3;
	cl::Buffer buffer_Q_JP1, buffer_Q_JP2, buffer_Q_JP3;
	cl::Buffer buffer_Q_JM1, buffer_Q_JM2, buffer_Q_JM3;
	cl::Buffer buffer_Q_KP1, buffer_Q_KP2, buffer_Q_KP3;
	cl::Buffer buffer_Q_KM1, buffer_Q_KM2, buffer_Q_KM3;
	cl::Buffer buffer_Qi_iplus_IP1, buffer_Qi_iplus_IP2, buffer_Qi_iplus_IP3; 
	cl::Buffer buffer_Qi_iplus_IM1, buffer_Qi_iplus_IM2, buffer_Qi_iplus_IM3; 
	cl::Buffer buffer_Qi_iminus_IP1, buffer_Qi_iminus_IP2, buffer_Qi_iminus_IP3; 
	cl::Buffer buffer_Qi_iminus_IM1, buffer_Qi_iminus_IM2, buffer_Qi_iminus_IM3; 
	cl::Buffer buffer_Qj_iplus_JP1, buffer_Qj_iplus_JP2, buffer_Qj_iplus_JP3; 
	cl::Buffer buffer_Qj_iplus_JM1, buffer_Qj_iplus_JM2, buffer_Qj_iplus_JM3; 
	cl::Buffer buffer_Qj_iminus_JP1, buffer_Qj_iminus_JP2, buffer_Qj_iminus_JP3; 
	cl::Buffer buffer_Qj_iminus_JM1, buffer_Qj_iminus_JM2, buffer_Qj_iminus_JM3; 
	cl::Buffer buffer_Qk_iplus_KP1, buffer_Qk_iplus_KP2, buffer_Qk_iplus_KP3; 
	cl::Buffer buffer_Qk_iplus_KM1, buffer_Qk_iplus_KM2, buffer_Qk_iplus_KM3; 
	cl::Buffer buffer_Qk_iminus_KP1, buffer_Qk_iminus_KP2, buffer_Qk_iminus_KP3; 
	cl::Buffer buffer_Qk_iminus_KM1, buffer_Qk_iminus_KM2, buffer_Qk_iminus_KM3; 
	cl::Buffer buffer_Qi_iplus_I, buffer_Qi_iminus_I;
	cl::Buffer buffer_Qj_iplus_I, buffer_Qj_iminus_I;
	cl::Buffer buffer_Qk_iplus_I, buffer_Qk_iminus_I;

	if (USE_GPU == 1) {
		buffer_j = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(j), NULL, NULL);
		buffer_Mach = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(Mach), NULL, NULL);
		buffer_del_t = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(del_t), NULL, NULL);
		buffer_deter = cl::Buffer(context, deter.begin(), deter.end(), true, false, &ret);
	
		buffer_U0 = cl::Buffer(context, U0.begin(), U0.end(), false, false, &ret);
		buffer_U1 = cl::Buffer(context, U1.begin(), U1.end(), false, false, &ret);
		buffer_U2 = cl::Buffer(context, U2.begin(), U2.end(), false, false, &ret);
		buffer_Q_C = cl::Buffer(context, Q_C.begin(), Q_C.end(), false, false, &ret);
		buffer_F_C = cl::Buffer(context, F_C.begin(), F_C.end(), false, false, &ret);
		buffer_E_C = cl::Buffer(context, E_C.begin(), E_C.end(), false, false, &ret);
		buffer_G_C = cl::Buffer(context, G_C.begin(), G_C.end(), false, false, &ret);
		buffer_Fv_C = cl::Buffer(context, Fv_C.begin(), Fv_C.end(), false, false, &ret);
		buffer_Ev_C = cl::Buffer(context, Ev_C.begin(), Ev_C.end(), false, false, &ret);
		buffer_Gv_C = cl::Buffer(context, Gv_C.begin(), Gv_C.end(), false, false, &ret);
		buffer_Grid = cl::Buffer(context, Grid_P.begin(), Grid_P.end(), true, false, &ret);
		buffer_Det = cl::Buffer(context, det.begin(), det.end(), true, false, &ret);
		buffer_FlowVar = cl::Buffer(context, FLOW.begin(), FLOW.end(), false, false, &ret);
		buffer_RoeVar = cl::Buffer(context, ROE_AVER.begin(), ROE_AVER.end(), false, false, &ret);
		buffer_metric = cl::Buffer(context, metric.begin(), metric.end(), true, false, &ret);
	}
	if (USE_GPU == 0) {
		buffer_j = cl::Buffer(context, CL_MEM_USE_HOST_PTR, sizeof(j), NULL, NULL);
		buffer_Mach = cl::Buffer(context, CL_MEM_USE_HOST_PTR, sizeof(Mach), NULL, NULL);
		buffer_del_t = cl::Buffer(context, CL_MEM_USE_HOST_PTR, sizeof(del_t), NULL, NULL);
		buffer_deter = cl::Buffer(context, deter.begin(), deter.end(), true, true, &ret);
	
		buffer_U0 = cl::Buffer(context, U0.begin(), U0.end(), false, true, &ret);
		buffer_U1 = cl::Buffer(context, U1.begin(), U1.end(), false, true, &ret);
		buffer_U2 = cl::Buffer(context, U2.begin(), U2.end(), false, true, &ret);
		buffer_Q_C = cl::Buffer(context, Q_C.begin(), Q_C.end(), false, true, &ret);
		buffer_F_C = cl::Buffer(context, F_C.begin(), F_C.end(), false, true, &ret);
		buffer_E_C = cl::Buffer(context, E_C.begin(), E_C.end(), false, true, &ret);
		buffer_G_C = cl::Buffer(context, G_C.begin(), G_C.end(), false, true, &ret);
		buffer_Fv_C = cl::Buffer(context, Fv_C.begin(), Fv_C.end(), false, true, &ret);
		buffer_Ev_C = cl::Buffer(context, Ev_C.begin(), Ev_C.end(), false, true, &ret);
		buffer_Gv_C = cl::Buffer(context, Gv_C.begin(), Gv_C.end(), false, true, &ret);
		buffer_Grid = cl::Buffer(context, Grid_P.begin(), Grid_P.end(), true, true, &ret);
		buffer_Det = cl::Buffer(context, det.begin(), det.end(), true, true, &ret);
		buffer_FlowVar = cl::Buffer(context, FLOW.begin(), FLOW.end(), false, true, &ret);
		buffer_RoeVar = cl::Buffer(context, ROE_AVER.begin(), ROE_AVER.end(), false, true, &ret);
		buffer_metric = cl::Buffer(context, metric.begin(), metric.end(), true, true, &ret);
	}
	/****************************Kernels arguments********************************************/
	kernel_U_Q.setArg(0, buffer_U0);
	kernel_U_Q.setArg(1, buffer_Q_C);
	kernel_U_Q.setArg(2, buffer_Grid);
	kernel_U_Q.setArg(3, buffer_Det);
	kernel_U_Q.setArg(4, buffer_FlowVar);
	kernel_U_Q.setArg(5, buffer_RoeVar);				

	kernel_viscous.setArg(0, buffer_U0);
	kernel_viscous.setArg(1, buffer_F_C);
	kernel_viscous.setArg(2, buffer_E_C);
	kernel_viscous.setArg(3, buffer_G_C);
	kernel_viscous.setArg(4, buffer_Fv_C);
	kernel_viscous.setArg(5, buffer_Ev_C);
	kernel_viscous.setArg(6, buffer_Gv_C);
	kernel_viscous.setArg(7, buffer_Grid);
	kernel_viscous.setArg(8, buffer_metric);
	kernel_viscous.setArg(9, buffer_Det);
	kernel_viscous.setArg(10, buffer_FlowVar);

	kernel_weno.setArg(0, sizeof(j), &j);
	kernel_weno.setArg(1, sizeof(Mach), &Mach);
	kernel_weno.setArg(2, sizeof(del_t), &del_t);
	kernel_weno.setArg(3, buffer_deter);
	kernel_weno.setArg(4, buffer_Det);
	kernel_weno.setArg(5, buffer_U0);
	kernel_weno.setArg(6, buffer_U1);
	kernel_weno.setArg(7, buffer_U2);
	kernel_weno.setArg(8, buffer_Grid);
	kernel_weno.setArg(9, buffer_Q_C);
	kernel_weno.setArg(10, buffer_metric);
	kernel_weno.setArg(11, buffer_RoeVar);
	kernel_weno.setArg(12, buffer_Fv_C);
	kernel_weno.setArg(13, buffer_Ev_C);
	kernel_weno.setArg(14, buffer_Gv_C);
	kernel_weno.setArg(15, buffer_F_C);
	kernel_weno.setArg(16, buffer_E_C);
	kernel_weno.setArg(17, buffer_G_C);
	kernel_weno.setArg(18, buffer_FlowVar);
	/*****************************************************************************************/

	int k, m, weno_exec, gar, lm, value, u_loc, v_loc, e_loc;
	double epsi, rc;
	int position, memsize, memsize1;
	char *buffer, *buffer2;
	
	//e_inf = e[0][10650];
	epsi = 0.003;

	rc = 0.3;
	epsi = 0.003;

	max_div_um = -1e25;
	max_div_vm = -1e25;
	max_div_wm = -1e25;
	max_div_em = -1e25;

	alpha_u.resize(5);
	alpha_v.resize(5);
	alpha_w.resize(5);
	alpha_u_ip.resize(5);
	alpha_u_im.resize(5);
	alpha_v_ip.resize(5);
	alpha_v_im.resize(5);
	alpha_w_ip.resize(5);
	alpha_w_im.resize(5);
	alpha_u_jp.resize(5);
	alpha_u_jm.resize(5);
	alpha_v_jp.resize(5);
	alpha_v_jm.resize(5);
	alpha_w_jp.resize(5);
	alpha_w_jm.resize(5);
	alpha_u_kp.resize(5);
	alpha_u_km.resize(5);
	alpha_v_kp.resize(5);
	alpha_v_km.resize(5);
	alpha_w_kp.resize(5);
	alpha_w_km.resize(5);

	lm = 0;
	for (iter = itera; iter<iterations; iter++)
	{
		max_div_um = -1e25;
		max_div_vm = -1e25;
		max_div_em = -1e25;
		gar = 0;

		weno_exec = 1;

		for (j = 0; j<3; j++)
		{
			lm++;
			gar++;

			/**************************viscous U_Q kernel********************************************/
			if (USE_GPU == 1) {
				cl::copy(queue, Grid_P.begin(), Grid_P.end(), buffer_Grid);
				cl::copy(queue, det.begin(), det.end(), buffer_Det);
				cl::copy(queue, FLOW.begin(), FLOW.end(), buffer_FlowVar);
			}	
			
			ret = queue.enqueueNDRangeKernel(kernel_U_Q, cl::NullRange, global, local, NULL, &event);
			if (ret != CL_SUCCESS) {
				std::cout << "enqueueNDRangeKernel Failed to execute " << std::endl;
				exit(1);
			}
			event.wait();
			queue.finish();
			
			if (USE_GPU == 1) {
				if (j == 0) {
					cl::copy(queue, buffer_U0, U0.begin(), U0.end());
				}
				if (j == 1) {
					cl::copy(queue, buffer_U1, U1.begin(), U1.end());
				}
				if (j == 2) {
					cl::copy(queue, buffer_U2, U2.begin(), U2.end());
				}
				cl::copy(queue, buffer_Q_C, Q_C.begin(), Q_C.end());
				cl::copy(queue, buffer_RoeVar, ROE_AVER.begin(), ROE_AVER.end());
			}
			/*******************************************************************************************/
/*			cl_ulong time_start;
			cl_ulong time_end;
			event.getProfilingInfo(CL_PROFILING_COMMAND_START, &time_start);
			event.getProfilingInfo(CL_PROFILING_COMMAND_END, &time_end);
			cl_double Exec_time = (time_end-time_start)/float(1e6);

			std::cout << "OpenCl Execution time on Rank " << myrank << " : " << Exec_time << " milliseconds \n";
*/
			/********************************************************************************************/
/*			for (i = 1; i<g_node; i++)
			{

				U[i][j][0] = U0[i].N0;
				U[i][j][1] = U0[i].N1;
				U[i][j][2] = U0[i].N2;
				U[i][j][3] = U0[i].N3;
				U[i][j][4] = U0[i].N4;

				Q[i][0].ip = Q_C[i].N0;
				Q[i][1].ip = Q_C[i].N1;
				Q[i][2].ip = Q_C[i].N2;
				Q[i][3].ip = Q_C[i].N3;
				Q[i][4].ip = Q_C[i].N4;
			}
*/
			/**************************Execute viscous kernel********************************************/
			if (USE_GPU == 1) {			
				cl::copy(queue_viscous, FLOW.begin(), FLOW.end(), buffer_FlowVar);
				if (j == 0) {
					cl::copy(queue_viscous, U0.begin(), U0.end(), buffer_U0);
				}
				if (j == 1) {
					cl::copy(queue_viscous, U1.begin(), U1.end(), buffer_U0);
				}
				if (j == 2) {
					cl::copy(queue_viscous, U2.begin(), U2.end(), buffer_U0);
				}
				cl::copy(queue_viscous, Grid_P.begin(), Grid_P.end(), buffer_Grid);
				cl::copy(queue_viscous, det.begin(), det.end(), buffer_Det);
				cl::copy(queue_viscous, metric.begin(), metric.end(), buffer_metric);
			}
			ret = queue_viscous.enqueueNDRangeKernel(kernel_viscous, cl::NullRange, global, local, NULL, &event);
			if (ret != CL_SUCCESS) {
				std::cout << "enqueueNDRangeKernel Failed to execute for viscous terms " << std::endl;
				exit(1);
			}
			event.wait();
			queue_viscous.finish();

			if (USE_GPU == 1) {	
				cl::copy(queue_viscous, buffer_F_C, F_C.begin(), F_C.end());
				cl::copy(queue_viscous, buffer_E_C, E_C.begin(), E_C.end());
				cl::copy(queue_viscous, buffer_G_C, G_C.begin(), G_C.end());
				cl::copy(queue_viscous, buffer_Fv_C, Fv_C.begin(), Fv_C.end());
				cl::copy(queue_viscous, buffer_Ev_C, Ev_C.begin(), Ev_C.end());
				cl::copy(queue_viscous, buffer_Gv_C, Gv_C.begin(), Gv_C.end());
				cl::copy(queue_viscous, buffer_FlowVar, FLOW.begin(), FLOW.end());
			}
			/*********************************************************************************************/

			for (i = 1; i<g_node; i++)
			{
		
				F[i][0] = F_C[i].N0;
				F[i][1] = F_C[i].N1;
				F[i][2] = F_C[i].N2;
				F[i][3] = F_C[i].N3;
				F[i][4] = F_C[i].N4;

				E[i][0] = E_C[i].N0;
				E[i][1] = E_C[i].N1;
				E[i][2] = E_C[i].N2;
				E[i][3] = E_C[i].N3;
				E[i][4] = E_C[i].N4;

				G[i][0] = G_C[i].N0;
				G[i][1] = G_C[i].N1;
				G[i][2] = G_C[i].N2;
				G[i][3] = G_C[i].N3;
				G[i][4] = G_C[i].N4;

				Fv[i][0] = Fv_C[i].N0;
				Fv[i][1] = Fv_C[i].N1;
				Fv[i][2] = Fv_C[i].N2;
				Fv[i][3] = Fv_C[i].N3;
				Fv[i][4] = Fv_C[i].N4;

				Ev[i][0] = Ev_C[i].N0;
				Ev[i][1] = Ev_C[i].N1;
				Ev[i][2] = Ev_C[i].N2;
				Ev[i][3] = Ev_C[i].N3;
				Ev[i][4] = Ev_C[i].N4;

				Gv[i][0] = Gv_C[i].N0;
				Gv[i][1] = Gv_C[i].N1;
				Gv[i][2] = Gv_C[i].N2;
				Gv[i][3] = Gv_C[i].N3;
				Gv[i][4] = Gv_C[i].N4;

				node[i].val = 0;
			}
			/********************************************************************************************************************************************************************************************/

			/**************************Execute weno kernel********************************************/
			if (USE_GPU == 1) {
				cl::copy(queue_weno, FLOW.begin(), FLOW.end(), buffer_FlowVar);
				cl::copy(queue_weno, U0.begin(), U0.end(), buffer_U0);
				cl::copy(queue_weno, U1.begin(), U1.end(), buffer_U1);
				cl::copy(queue_weno, U2.begin(), U2.end(), buffer_U2);
				cl::copy(queue_weno, Q_C.begin(), Q_C.end(), buffer_Q_C);
				cl::copy(queue_weno, Grid_P.begin(), Grid_P.end(), buffer_Grid);
				cl::copy(queue_weno, det.begin(), det.end(), buffer_Det);
				cl::copy(queue_weno, deter.begin(), deter.end(), buffer_deter);
				cl::copy(queue_weno, metric.begin(), metric.end(), buffer_metric);
				cl::copy(queue_weno, ROE_AVER.begin(), ROE_AVER.end(), buffer_RoeVar);
				cl::copy(queue_weno, F_C.begin(), F_C.end(), buffer_F_C);
				cl::copy(queue_weno, E_C.begin(), E_C.end(), buffer_E_C);
				cl::copy(queue_weno, G_C.begin(), G_C.end(), buffer_G_C);
				cl::copy(queue_weno, Fv_C.begin(), Fv_C.end(), buffer_Fv_C);
				cl::copy(queue_weno, Ev_C.begin(), Ev_C.end(), buffer_Ev_C);
				cl::copy(queue_weno, Gv_C.begin(), Gv_C.end(), buffer_Gv_C);
			}			
			ret = queue_weno.enqueueNDRangeKernel(kernel_weno, cl::NullRange, Weno_global, local, NULL, &event);
			if (ret != CL_SUCCESS) {
				std::cout << "enqueueNDRangeKernel Failed to execute weno kernel " << std::endl;
				exit(1);
			}
			event.wait();
			queue_weno.finish();

			if (USE_GPU == 1) {
			/*	if (j == 0) {
					cl::copy(queue_weno, buffer_U0, U0.begin(), U0.end());
				}*/
				if (j == 1) {
					cl::copy(queue_weno, buffer_U1, U1.begin(), U1.end());
				}
				if (j == 2) {
					cl::copy(queue_weno, buffer_U2, U2.begin(), U2.end());
				}

				cl::copy(queue_weno, buffer_FlowVar, FLOW.begin(), FLOW.end());
			}

			/*************************************************************************************************/

			/***************************Test for OpenCL*****************************************************/
			for (i=1; i<=sd_node; i++) {
				if (j == 2)
				{
					diver[i][1].u = FLOW[i].u;
					diver[i][1].v = FLOW[i].v;
					diver[i][1].e = FLOW[i].e;

					if (max_div_um < fabs(diver[i][1].u - diver[i][0].u))
					{
						max_div_um = fabs(diver[i][1].u - diver[i][0].u);
						u_loc = i;
					}
					if (max_div_vm < fabs(diver[i][1].v - diver[i][0].v))
					{
						max_div_vm = fabs(diver[i][1].v - diver[i][0].v);
						v_loc = i;
					}
					if (max_div_em < fabs(diver[i][1].e - diver[i][0].e))
					{
						max_div_em = fabs(diver[i][1].e - diver[i][0].e);
						e_loc = i;
					}

					diver[i][0].u = diver[i][1].u;
					diver[i][0].v = diver[i][1].v;
					diver[i][0].e = diver[i][1].e;
				}
			}
			/***********************************************************************************************/

			position = 0;
			for (i = 0; i<neigh_pro; i++)
			{
				MPI_Pack_size(recv_c[c[i]] * 9, MPI_DOUBLE, MPI_COMM_WORLD, &memsize1);
				buffer2 = new char[memsize1];                                          /***********carefull with buffer1 and buffer2******************/
				MPI_Pack_size(proc_node[c[i]] * 9, MPI_DOUBLE, MPI_COMM_WORLD, &memsize);
				buffer = new char[memsize];
				position = 0;
				for (k = 0; k<recv_c[c[i]]; k++)
				{
					MPI_Pack(&FLOW[loc_dat[i][k]].u, 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}
				for (k = 0; k<recv_c[c[i]]; k++)
				{
					MPI_Pack(&FLOW[loc_dat[i][k]].v, 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}
				for (k = 0; k<recv_c[c[i]]; k++)
				{
					MPI_Pack(&FLOW[loc_dat[i][k]].w, 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}
				for (k = 0; k<recv_c[c[i]]; k++)
				{
					MPI_Pack(&FLOW[loc_dat[i][k]].p, 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}
				for (k = 0; k<recv_c[c[i]]; k++)
				{
					MPI_Pack(&FLOW[loc_dat[i][k]].t, 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}
				for (k = 0; k<recv_c[c[i]]; k++)
				{
					MPI_Pack(&FLOW[loc_dat[i][k]].rho, 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}
				for (k = 0; k<recv_c[c[i]]; k++)
				{
					MPI_Pack(&FLOW[loc_dat[i][k]].e, 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}
				for (k = 0; k<recv_c[c[i]]; k++)
				{
					MPI_Pack(&FLOW[loc_dat[i][k]].mu, 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}

				MPI_Sendrecv(buffer2, memsize1, MPI_PACKED, c[i], c[i], buffer, memsize, MPI_PACKED, c[i], myrank, MPI_COMM_WORLD, &status);
				position = 0;

				for (k = 0; k<proc_node[c[i]]; k++)
				{
					MPI_Unpack(buffer, memsize, &position, &FLOW[recv_b[i][k]].u, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				}
				for (k = 0; k<proc_node[c[i]]; k++)
				{
					MPI_Unpack(buffer, memsize, &position, &FLOW[recv_b[i][k]].v, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				}
				for (k = 0; k<proc_node[c[i]]; k++)
				{
					MPI_Unpack(buffer, memsize, &position, &FLOW[recv_b[i][k]].w, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				}
				for (k = 0; k<proc_node[c[i]]; k++)
				{
					MPI_Unpack(buffer, memsize, &position, &FLOW[recv_b[i][k]].p, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				}
				for (k = 0; k<proc_node[c[i]]; k++)
				{
					MPI_Unpack(buffer, memsize, &position, &FLOW[recv_b[i][k]].t, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				}
				for (k = 0; k<proc_node[c[i]]; k++)
				{
					MPI_Unpack(buffer, memsize, &position, &FLOW[recv_b[i][k]].rho, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				}
				for (k = 0; k<proc_node[c[i]]; k++)
				{
					MPI_Unpack(buffer, memsize, &position, &FLOW[recv_b[i][k]].e, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				}
				for (k = 0; k<proc_node[c[i]]; k++)
				{
					MPI_Unpack(buffer, memsize, &position, &FLOW[recv_b[i][k]].mu, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				}
				delete[] buffer;
				delete[] buffer2;
				k = 0;
			}

			initial = 1;
			intialise(lm);

			for (k = 1; k <= no_of_tip_send; k++)
			{
				FLOW[tip[k].node].u = 0.0;
				FLOW[tip[k].node].v = 0.0;
				FLOW[tip[k].node].w = 0.0;
				FLOW[tip[k].node].p = FLOW[node[tip[k].node].n_n[3]].p;
				FLOW[tip[k].node].t = FLOW[node[tip[k].node].n_n[3]].t;
				FLOW[tip[k].node].rho = (1.4*Mach*Mach)*(FLOW[tip[k].node].p / FLOW[tip[k].node].t);
				FLOW[tip[k].node].e = FLOW[tip[k].node].p / (0.4*FLOW[tip[k].node].rho);
			}

			for (k = 1; k <= no_of_tip_recv; k++)
			{
				FLOW[tip_recv[k].node].u = 0.0;
				FLOW[tip_recv[k].node].v = 0.0;
				FLOW[tip_recv[k].node].w = 0.0;
				FLOW[tip_recv[k].node].p = FLOW[node[tip_recv[k].node].n_n[3]].p;
				FLOW[tip_recv[k].node].t = FLOW[node[tip_recv[k].node].n_n[3]].t;
				FLOW[tip_recv[k].node].rho = (1.4*Mach*Mach)*(FLOW[tip_recv[k].node].p / FLOW[tip_recv[k].node].t);
				FLOW[tip_recv[k].node].e = FLOW[tip_recv[k].node].p / (0.4*FLOW[tip_recv[k].node].rho);
			}
		}

		position = 0;
		MPI_Pack_size(3, MPI_DOUBLE, MPI_COMM_WORLD, &memsize);
		buffer = new char[memsize];
		MPI_Pack(&max_div_um, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
		MPI_Pack(&max_div_vm, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
		MPI_Pack(&max_div_em, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
		MPI_Send(buffer, position, MPI_PACKED, 0, myrank, MPI_COMM_WORLD);
		delete[] buffer;

		if ((iter % 100) == 0)
		{
			position = 0;
			MPI_Pack_size(sd_node * 10, MPI_DOUBLE, MPI_COMM_WORLD, &memsize);
			buffer = new char[memsize];
			for (i=1; i<= sd_node; i++) {
				MPI_Pack(&FLOW[i].u, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			}
			for (i=1; i<= sd_node; i++) {
				MPI_Pack(&FLOW[i].v, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			}
			for (i=1; i<= sd_node; i++) {
				MPI_Pack(&FLOW[i].w, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			}
			for (i=1; i<= sd_node; i++) {
				MPI_Pack(&FLOW[i].rho, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			}
			for (i=1; i<= sd_node; i++) {
				MPI_Pack(&FLOW[i].p, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			}
			for (i=1; i<= sd_node; i++) {
				MPI_Pack(&FLOW[i].t, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			}
			for (i=1; i<= sd_node; i++) {
				MPI_Pack(&FLOW[i].e, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			}
			for (i=1; i<= sd_node; i++) {
				MPI_Pack(&FLOW[i].mu, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			}
			for (i=1; i<= sd_node; i++) {
				MPI_Pack(&DUCROS[i], 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			}
			MPI_Send(buffer, position, MPI_PACKED, 0, myrank, MPI_COMM_WORLD);
			delete[] buffer;
		}
	}
}
