#ifndef _GQD_API_H_
#define _GQD_API_H_

#include <gqd_type.h>


void GQDInit();
void GQDFree();

void GDDInit();
void GDDFree();

void gpu_dd_add(const GPU_dd* h_a, const GPU_dd* h_b,
                GPU_dd* h_c, const int len);

void gpu_qd_add(const GPU_qd* h_a, const GPU_qd* h_b,
                GPU_qd* h_c, const int len);

void gpu_dd_sub(const GPU_dd* h_a, const GPU_dd* h_b, 
		GPU_dd* h_c, const int len);

void gpu_dd_mul(const GPU_dd* h_a, const GPU_dd* h_b,
                GPU_dd* h_c, const int len);

void gpu_qd_mul(const GPU_qd* h_a, const GPU_qd* h_b,
                GPU_qd* h_c, const int len);


void gpu_dd_div(const GPU_dd* h_a, const GPU_dd* h_b,
                GPU_dd* h_c, const int len);

void gpu_qd_div(const GPU_qd* h_a, const GPU_qd* h_b,
                GPU_qd* h_c, const int len);

void gpu_dd_sqrt(const GPU_dd* h_a, GPU_dd* h_c, const int len);

void gpu_qd_sqrt(const GPU_qd* h_a, GPU_qd* h_c, const int len);


void gpu_dd_exp(const GPU_dd* h_a, GPU_dd* h_c, const int len);

void gpu_qd_exp(const GPU_qd* h_a, GPU_qd* h_c, const int len);


void gpu_dd_sin(const GPU_dd* h_a, GPU_dd* h_c, const int len);

void gpu_qd_sin(const GPU_qd* h_a, GPU_qd* h_c, const int len);


#endif //_GQD_API_H_
