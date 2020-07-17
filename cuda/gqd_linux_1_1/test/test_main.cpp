

#include <stdio.h>
#include <stdlib.h>
#include <qd/qd_real.h>
#include <sys/time.h>
#include <omp.h>

#include "gqd_api.h"
#include "gqd_util.h"


double getSec( struct timeval tvStart, struct timeval tvEnd ) {
        double tStart = (double)tvStart.tv_sec + 1e-6*tvStart.tv_usec;
        double tEnd = (double)tvEnd.tv_sec + 1e-6*tvEnd.tv_usec;
        return (tEnd - tStart);
}



void test_add_dd(const int len) {
	printf("*** addition, double-double, #element = %d ***\n", len);
	dd_real* h_a = (dd_real*)malloc(sizeof(dd_real)*len);
	dd_real* h_b = (dd_real*)malloc(sizeof(dd_real)*len);
	dd_real* gold_c = (dd_real*)malloc(sizeof(dd_real)*len);
	dd_real* ref_c = (dd_real*)malloc(sizeof(dd_real)*len);
	for(int i = 0; i < len; i++) {
		rand(&(h_a[i]));
		rand(&(h_b[i]));
	}

	//GPU
	GPU_dd* gpu_a = (GPU_dd*)malloc(sizeof(GPU_dd)*len);
	GPU_dd* gpu_b = (GPU_dd*)malloc(sizeof(GPU_dd)*len);
	GPU_dd* gpu_c = (GPU_dd*)malloc(sizeof(GPU_dd)*len);
	dd_to_gdd(h_a, gpu_a, len);
	dd_to_gdd(h_b, gpu_b, len);

	gpu_dd_add(gpu_a, gpu_b, gpu_c, len);

	gdd_to_dd(gpu_c, ref_c, len);
	free(gpu_a);
	free(gpu_b);
	free(gpu_c);

	//CPU
	struct timeval start, end;
	gettimeofday(&start, NULL);
#pragma omp parallel for	
	for(int i = 0; i < len; i++) {
		gold_c[i] = h_a[i] + h_b[i];
	}
	gettimeofday(&end, NULL);
	printf("CPU add, dd_real: %f sec\n", getSec(start, end));

	//check result
	checkTwoArray(gold_c, ref_c, len);

	//clean up
	free(h_a);
	free(h_b);
	free(gold_c);
	free(ref_c);
	printf("\n");
}

void test_sub_dd2() {
	const int len = 1;

	printf("*** addition, double-double, #element = %d ***\n", len);
	dd_real* h_a = (dd_real*)malloc(sizeof(dd_real)*len);
	dd_real* h_b = (dd_real*)malloc(sizeof(dd_real)*len);
	dd_real* gold_c = (dd_real*)malloc(sizeof(dd_real)*len);
	dd_real* ref_c = (dd_real*)malloc(sizeof(dd_real)*len);

	h_a[0].x[0] = 5.0013118108914589e-01;
	h_a[0].x[1] = -1.7283934324986953e-17;

	h_b[0].x[0] = 5.001311810891459e-01;
	h_b[0].x[1] = -5.5499622595091138e-17;

	//GPU
	GPU_dd* gpu_a = (GPU_dd*)malloc(sizeof(GPU_dd)*len);
	GPU_dd* gpu_b = (GPU_dd*)malloc(sizeof(GPU_dd)*len);
	GPU_dd* gpu_c = (GPU_dd*)malloc(sizeof(GPU_dd)*len);
	dd_to_gdd(h_a, gpu_a, len);
	dd_to_gdd(h_b, gpu_b, len);

	gpu_dd_sub(gpu_a, gpu_b, gpu_c, len);

	gdd_to_dd(gpu_c, ref_c, len);
	free(gpu_a);
	free(gpu_b);
	free(gpu_c);

	//CPU
	struct timeval start, end;
	gettimeofday(&start, NULL);
#pragma omp parallel for	
	for(int i = 0; i < len; i++) {
		gold_c[i] = h_a[i] - h_b[i];
	}
	gettimeofday(&end, NULL);
	printf("CPU add, dd_real: %f sec\n", getSec(start, end));

	//check result
	checkTwoArray(gold_c, ref_c, len);

	//clean up
	free(h_a);
	free(h_b);
	free(gold_c);
	free(ref_c);
	printf("\n");
}

void test_add_qd(const int len) {
	printf("*** addition, quad-double, #element = %d ***\n", len);
	qd_real* h_a = (qd_real*)malloc(sizeof(qd_real)*len);
	qd_real* h_b = (qd_real*)malloc(sizeof(qd_real)*len);
	qd_real* gold_c = (qd_real*)malloc(sizeof(qd_real)*len);
	qd_real* ref_c = (qd_real*)malloc(sizeof(qd_real)*len);
	for(int i = 0; i < len; i++) {
		rand(&(h_a[i]));
		rand(&(h_b[i]));
	}

	//GPU
	GPU_qd* gpu_a = (GPU_qd*)malloc(sizeof(GPU_qd)*len);
	GPU_qd* gpu_b = (GPU_qd*)malloc(sizeof(GPU_qd)*len);
	GPU_qd* gpu_c = (GPU_qd*)malloc(sizeof(GPU_qd)*len);
	qd_to_gqd(h_a, gpu_a, len);
	qd_to_gqd(h_b, gpu_b, len);

	gpu_qd_add(gpu_a, gpu_b, gpu_c, len);

	gqd_to_qd(gpu_c, ref_c, len);
	free(gpu_a);
	free(gpu_b);
	free(gpu_c);

	//CPU
	struct timeval start, end;
	gettimeofday(&start, NULL);
#pragma omp parallel for
	for(int i = 0; i < len; i++) {
		gold_c[i] = h_a[i] + h_b[i];
	}
	gettimeofday(&end, NULL);
	printf("CPU add, qd_real: %f sec\n", getSec(start, end));

	//check result
	checkTwoArray(gold_c, ref_c, len);

	//clean up
	free(h_a);
	free(h_b);
	free(gold_c);
	free(ref_c);

	printf("\n");
}

void test_mul_dd(const int len) {
	printf("*** multiplication, double-double, #element = %d ***\n", len);
	dd_real* h_a = (dd_real*)malloc(sizeof(dd_real)*len);
	dd_real* h_b = (dd_real*)malloc(sizeof(dd_real)*len);
	dd_real* gold_c = (dd_real*)malloc(sizeof(dd_real)*len);
	dd_real* ref_c = (dd_real*)malloc(sizeof(dd_real)*len);
	for(int i = 0; i < len; i++) {
		rand(&(h_a[i]));
		rand(&(h_b[i]));
	}

	//GPU
	GPU_dd* gpu_a = (GPU_dd*)malloc(sizeof(GPU_dd)*len);
	GPU_dd* gpu_b = (GPU_dd*)malloc(sizeof(GPU_dd)*len);
	GPU_dd* gpu_c = (GPU_dd*)malloc(sizeof(GPU_dd)*len);
	dd_to_gdd(h_a, gpu_a, len);
	dd_to_gdd(h_b, gpu_b, len);

	gpu_dd_mul(gpu_a, gpu_b, gpu_c, len);

	gdd_to_dd(gpu_c, ref_c, len);
	free(gpu_a);
	free(gpu_b);
	free(gpu_c);

	//CPU
	struct timeval start, end;
	gettimeofday(&start, NULL);
#pragma omp parallel for
	for(int i = 0; i < len; i++) {
		gold_c[i] = h_a[i] * h_b[i];
	}
	gettimeofday(&end, NULL);
	printf("CPU dd_real, mul: %f sec.\n", getSec(start, end));

	//check result
	checkTwoArray(gold_c, ref_c, len);

	//clean up
	free(h_a);
	free(h_b);
	free(gold_c);
	free(ref_c);
	printf("\n");
}

void test_mul_qd(const int len) {
	printf("*** multiplication, quad-double, #element = %d ***\n", len);
	qd_real* h_a = (qd_real*)malloc(sizeof(qd_real)*len);
	qd_real* h_b = (qd_real*)malloc(sizeof(qd_real)*len);
	qd_real* gold_c = (qd_real*)malloc(sizeof(qd_real)*len);
	qd_real* ref_c = (qd_real*)malloc(sizeof(qd_real)*len);
	for(int i = 0; i < len; i++) {
		rand(&(h_a[i]));
		rand(&(h_b[i]));
	}

	//GPU
	GPU_qd* gpu_a = (GPU_qd*)malloc(sizeof(GPU_qd)*len);
	GPU_qd* gpu_b = (GPU_qd*)malloc(sizeof(GPU_qd)*len);
	GPU_qd* gpu_c = (GPU_qd*)malloc(sizeof(GPU_qd)*len);
	qd_to_gqd(h_a, gpu_a, len);
	qd_to_gqd(h_b, gpu_b, len);

	gpu_qd_mul(gpu_a, gpu_b, gpu_c, len);

	gqd_to_qd(gpu_c, ref_c, len);
	free(gpu_a);
	free(gpu_b);
	free(gpu_c);

	//CPU
	struct timeval start, end;
	gettimeofday(&start, NULL);
#pragma omp parallel for
	for(int i = 0; i < len; i++) {
		gold_c[i] = h_a[i] * h_b[i];
	}
	gettimeofday(&end, NULL);
	printf("CPU qd_real, mul: %f sec\n", getSec(start, end));

	//check result
	int id = checkTwoArray(gold_c, ref_c, len);
	/*printf("\n***INPUT A***\n");
	printf("%.16e\n", h_a[id].x[0]);
	printf("%.16e\n", h_a[id].x[1]);
	printf("%.16e\n", h_a[id].x[2]);
	printf("%.16e\n", h_a[id].x[3]);

	printf("\n***INPUT B***\n");
	printf("%.16e\n", h_b[id].x[0]);
	printf("%.16e\n", h_b[id].x[1]);
	printf("%.16e\n", h_b[id].x[2]);
	printf("%.16e\n", h_b[id].x[3]);
	*/

	//clean up
	free(h_a);
	free(h_b);
	free(gold_c);
	free(ref_c);

	printf("\n");
}


void test_mul_qd2() {
	const int len = 1;
	printf("*** multiplication, quad-double, #element = %d ***\n", len);
	qd_real* h_a = (qd_real*)malloc(sizeof(qd_real)*len);
	qd_real* h_b = (qd_real*)malloc(sizeof(qd_real)*len);
	qd_real* gold_c = (qd_real*)malloc(sizeof(qd_real)*len);
	qd_real* ref_c = (qd_real*)malloc(sizeof(qd_real)*len);
	for(int i = 0; i < len; i++) {
		rand(&(h_a[i]));
		rand(&(h_b[i]));
	}

	//GPU
	GPU_qd* gpu_a = (GPU_qd*)malloc(sizeof(GPU_qd)*len);
	GPU_qd* gpu_b = (GPU_qd*)malloc(sizeof(GPU_qd)*len);
	GPU_qd* gpu_c = (GPU_qd*)malloc(sizeof(GPU_qd)*len);
	qd_to_gqd(h_a, gpu_a, len);
	qd_to_gqd(h_b, gpu_b, len);

	gpu_qd_mul(gpu_a, gpu_b, gpu_c, len);

	gqd_to_qd(gpu_c, ref_c, len);
	free(gpu_a);
	free(gpu_b);
	free(gpu_c);

	//CPU
	struct timeval start, end;
	gettimeofday(&start, NULL);
#pragma omp parallel for
	for(int i = 0; i < len; i++) {
		gold_c[i] = h_a[i] * h_b[i];
	}
	gettimeofday(&end, NULL);
	printf("CPU qd_real, mul: %f sec\n", getSec(start, end));

	//check result
	int id = 0;
	printf("\n***INPUT A***\n");
	printf("%.16e\n", h_a[id].x[0]);
	printf("%.16e\n", h_a[id].x[1]);
	printf("%.16e\n", h_a[id].x[2]);
	printf("%.16e\n", h_a[id].x[3]);

	printf("\n***INPUT B***\n");
	printf("%.16e\n", h_b[id].x[0]);
	printf("%.16e\n", h_b[id].x[1]);
	printf("%.16e\n", h_b[id].x[2]);
	printf("%.16e\n", h_b[id].x[3]);

	printf("\n*** GOLD RESULT ***\n");
	printf("%.16e\n", gold_c[id].x[0]);
	printf("%.16e\n", gold_c[id].x[1]);
	printf("%.16e\n", gold_c[id].x[2]);
	printf("%.16e\n", gold_c[id].x[3]);

	printf("\n*** GPU RESULT ***\n");
	printf("%.16e\n", ref_c[id].x[0]);
	printf("%.16e\n", ref_c[id].x[1]);
	printf("%.16e\n", ref_c[id].x[2]);
	printf("%.16e\n", ref_c[id].x[3]);


	//clean up
	free(h_a);
	free(h_b);
	free(gold_c);
	free(ref_c);

	printf("\n");
}

void test_div_dd(const int len) {
	printf("*** division, double-double, #element = %d ***\n", len);
	dd_real* h_a = (dd_real*)malloc(sizeof(dd_real)*len);
	dd_real* h_b = (dd_real*)malloc(sizeof(dd_real)*len);
	dd_real* gold_c = (dd_real*)malloc(sizeof(dd_real)*len);
	dd_real* ref_c = (dd_real*)malloc(sizeof(dd_real)*len);
	for(int i = 0; i < len; i++) {
		rand(&(h_a[i]));
		rand(&(h_b[i]));
	}

	//GPU
	GPU_dd* gpu_a = (GPU_dd*)malloc(sizeof(GPU_dd)*len);
	GPU_dd* gpu_b = (GPU_dd*)malloc(sizeof(GPU_dd)*len);
	GPU_dd* gpu_c = (GPU_dd*)malloc(sizeof(GPU_dd)*len);
	dd_to_gdd(h_a, gpu_a, len);
	dd_to_gdd(h_b, gpu_b, len);

	gpu_dd_div(gpu_a, gpu_b, gpu_c, len);

	gdd_to_dd(gpu_c, ref_c, len);
	free(gpu_a);
	free(gpu_b);
	free(gpu_c);

	//CPU
	struct timeval start, end;
	gettimeofday(&start, NULL);
#pragma omp parallel for
	for(int i = 0; i < len; i++) {
		gold_c[i] = h_a[i] / h_b[i];
	}
	gettimeofday(&end, NULL);
	printf("CPU dd_real, div: %f sec.\n", getSec(start, end));

	//check result
	checkTwoArray(gold_c, ref_c, len);

	//clean up
	free(h_a);
	free(h_b);
	free(gold_c);
	free(ref_c);
	printf("\n");
}

void test_div_qd(const int len) {
	printf("*** division, quad-double, #element = %d ***\n", len);
	qd_real* h_a = (qd_real*)malloc(sizeof(qd_real)*len);
	qd_real* h_b = (qd_real*)malloc(sizeof(qd_real)*len);
	qd_real* gold_c = (qd_real*)malloc(sizeof(qd_real)*len);
	qd_real* ref_c = (qd_real*)malloc(sizeof(qd_real)*len);
	for(int i = 0; i < len; i++) {
		rand(&(h_a[i]));
		rand(&(h_b[i]));
	}

	//GPU
	GPU_qd* gpu_a = (GPU_qd*)malloc(sizeof(GPU_qd)*len);
	GPU_qd* gpu_b = (GPU_qd*)malloc(sizeof(GPU_qd)*len);
	GPU_qd* gpu_c = (GPU_qd*)malloc(sizeof(GPU_qd)*len);
	qd_to_gqd(h_a, gpu_a, len);
	qd_to_gqd(h_b, gpu_b, len);

	gpu_qd_div(gpu_a, gpu_b, gpu_c, len);

	gqd_to_qd(gpu_c, ref_c, len);
	free(gpu_a);
	free(gpu_b);
	free(gpu_c);

	//CPU
	struct timeval start, end;
	gettimeofday(&start, NULL);
#pragma omp parallel for
	for(int i = 0; i < len; i++) {
		gold_c[i] = h_a[i] / h_b[i];
	}
	gettimeofday(&end, NULL);
	printf("CPU qd_real, div: %f sec\n", getSec(start, end));

	//check result
	checkTwoArray(gold_c, ref_c, len);

	//clean up
	free(h_a);
	free(h_b);
	free(gold_c);
	free(ref_c);

	printf("\n");
}


void test_sqrt_dd(const int len) {
	printf("*** square root, double-double, #element = %d ***\n", len);
	dd_real* h_a = (dd_real*)malloc(sizeof(dd_real)*len);
	dd_real* gold_c = (dd_real*)malloc(sizeof(dd_real)*len);
	dd_real* ref_c = (dd_real*)malloc(sizeof(dd_real)*len);
	for(int i = 0; i < len; i++) {
		rand(&(h_a[i]));
	}

	//GPU
	GPU_dd* gpu_a = (GPU_dd*)malloc(sizeof(GPU_dd)*len);
	GPU_dd* gpu_c = (GPU_dd*)malloc(sizeof(GPU_dd)*len);
	dd_to_gdd(h_a, gpu_a, len);

	gpu_dd_sqrt(gpu_a, gpu_c, len);

	gdd_to_dd(gpu_c, ref_c, len);
	free(gpu_a);
	free(gpu_c);

	//CPU
	struct timeval start, end;
	gettimeofday(&start, NULL);
#pragma omp parallel for
	for(int i = 0; i < len; i++) {
		gold_c[i] = sqrt(h_a[i]);;
	}
	gettimeofday(&end, NULL);
	printf("CPU dd_real, sqrt: %f sec.\n", getSec(start, end));

	//check result
	int i = checkTwoArray(gold_c, ref_c, len);

	//print the components of element with the max. error
	printf("*** Input components ***\n");
	printf("x1 = %.16e\n", h_a[i].x[0]);
	printf("x2 = %.16e\n", h_a[i].x[1]);	
	printf("***********************\n");

	//clean up
	free(h_a);
	free(gold_c);
	free(ref_c);
	printf("\n");
}


void test_sqrt_dd2() {
	const int len = 1;

	printf("*** square root, double-double, #element = %d ***\n", len);
	dd_real* h_a = (dd_real*)malloc(sizeof(dd_real)*len);
	dd_real* gold_c = (dd_real*)malloc(sizeof(dd_real)*len);
	dd_real* ref_c = (dd_real*)malloc(sizeof(dd_real)*len);
	//for(int i = 0; i < len; i++) {
	//	rand(&(h_a[i]));
	//}
	h_a[0].x[0] = 5.0013118108914589e-01;
	h_a[0].x[1] = -1.7283934324986953e-17;


	//GPU
	GPU_dd* gpu_a = (GPU_dd*)malloc(sizeof(GPU_dd)*len);
	GPU_dd* gpu_c = (GPU_dd*)malloc(sizeof(GPU_dd)*len);
	dd_to_gdd(h_a, gpu_a, len);

	gpu_dd_sqrt(gpu_a, gpu_c, len);

	gdd_to_dd(gpu_c, ref_c, len);
	free(gpu_a);
	free(gpu_c);

	//CPU
	struct timeval start, end;
	gettimeofday(&start, NULL);
#pragma omp parallel for
	for(int i = 0; i < len; i++) {
		gold_c[i] = sqrt(h_a[i]);;
	}
	gettimeofday(&end, NULL);
	printf("CPU dd_real, sqrt: %f sec.\n", getSec(start, end));

	//check result
	printf("*** OUTPUT ***\n");
	printf("GOLD: (%.18e, %.18e)\n", gold_c[0].x[0], gold_c[0].x[1]);
	printf("GPU:  (%.18e, %.18e)\n", ref_c[0].x[0], ref_c[0].x[1]);
	printf("**************\n");

	//clean up
	free(h_a);
	free(gold_c);
	free(ref_c);
	printf("\n");
}


void test_sqrt_qd(const int len) {
	printf("*** square root, quad-double, #element = %d ***\n", len);
	qd_real* h_a = (qd_real*)malloc(sizeof(qd_real)*len);
	qd_real* gold_c = (qd_real*)malloc(sizeof(qd_real)*len);
	qd_real* ref_c = (qd_real*)malloc(sizeof(qd_real)*len);
	for(int i = 0; i < len; i++) {
		rand(&(h_a[i]));
	}

	//GPU
	GPU_qd* gpu_a = (GPU_qd*)malloc(sizeof(GPU_qd)*len);
	GPU_qd* gpu_c = (GPU_qd*)malloc(sizeof(GPU_qd)*len);
	qd_to_gqd(h_a, gpu_a, len);

	gpu_qd_sqrt(gpu_a, gpu_c, len);

	gqd_to_qd(gpu_c, ref_c, len);
	free(gpu_a);
	free(gpu_c);

	//CPU
	struct timeval start, end;
	gettimeofday(&start, NULL);
#pragma omp parallel for
	for(int i = 0; i < len; i++) {
		gold_c[i] = sqrt(h_a[i]);
	}
	gettimeofday(&end, NULL);
	printf("CPU qd_real, sqrt: %f sec\n", getSec(start, end));

	//check result
	checkTwoArray(gold_c, ref_c, len);

	//clean up
	free(h_a);
	free(gold_c);
	free(ref_c);

	printf("\n");
}

void test_exp_dd(const int len) {
	printf("*** exponent, double-double, #element = %d ***\n", len);
	dd_real* h_a = (dd_real*)malloc(sizeof(dd_real)*len);
	dd_real* gold_c = (dd_real*)malloc(sizeof(dd_real)*len);
	dd_real* ref_c = (dd_real*)malloc(sizeof(dd_real)*len);
	for(int i = 0; i < len; i++) {
		rand(&(h_a[i]));
	}

	//GPU
	GPU_dd* gpu_a = (GPU_dd*)malloc(sizeof(GPU_dd)*len);
	GPU_dd* gpu_c = (GPU_dd*)malloc(sizeof(GPU_dd)*len);
	dd_to_gdd(h_a, gpu_a, len);

	gpu_dd_exp(gpu_a, gpu_c, len);

	gdd_to_dd(gpu_c, ref_c, len);
	free(gpu_a);
	free(gpu_c);

	//CPU
	struct timeval start, end;
	gettimeofday(&start, NULL);
#pragma omp parallel for
	for(int i = 0; i < len; i++) {
		gold_c[i] = exp(h_a[i]);;
	}
	gettimeofday(&end, NULL);
	printf("CPU dd_real, exp: %f sec.\n", getSec(start, end));

	//check result
	checkTwoArray(gold_c, ref_c, len);

	//clean up
	free(h_a);
	free(gold_c);
	free(ref_c);
	printf("\n");
}

void test_exp_qd(const int len) {
	printf("*** exponent, quad-double, #element = %d ***\n", len);
	qd_real* h_a = (qd_real*)malloc(sizeof(qd_real)*len);
	qd_real* gold_c = (qd_real*)malloc(sizeof(qd_real)*len);
	qd_real* ref_c = (qd_real*)malloc(sizeof(qd_real)*len);
	for(int i = 0; i < len; i++) {
		rand(&(h_a[i]));
	}

	//GPU
	GPU_qd* gpu_a = (GPU_qd*)malloc(sizeof(GPU_qd)*len);
	GPU_qd* gpu_c = (GPU_qd*)malloc(sizeof(GPU_qd)*len);
	qd_to_gqd(h_a, gpu_a, len);

	gpu_qd_exp(gpu_a, gpu_c, len);

	gqd_to_qd(gpu_c, ref_c, len);
	free(gpu_a);
	free(gpu_c);

	//CPU
	struct timeval start, end;
	gettimeofday(&start, NULL);
#pragma omp parallel for
	for(int i = 0; i < len; i++) {
		gold_c[i] = exp(h_a[i]);
	}
	gettimeofday(&end, NULL);
	printf("CPU qd_real, exp: %f sec\n", getSec(start, end));

	//check result
	checkTwoArray(gold_c, ref_c, len);

	//clean up
	free(h_a);
	free(gold_c);
	free(ref_c);

	printf("\n");
}

void test_sin_dd(const int len) {
	printf("*** sine, double-double, #element = %d ***\n", len);
	dd_real* h_a = (dd_real*)malloc(sizeof(dd_real)*len);
	dd_real* gold_c = (dd_real*)malloc(sizeof(dd_real)*len);
	dd_real* ref_c = (dd_real*)malloc(sizeof(dd_real)*len);
	for(int i = 0; i < len; i++) {
		rand(&(h_a[i]));
	}

	//GPU
	GPU_dd* gpu_a = (GPU_dd*)malloc(sizeof(GPU_dd)*len);
	GPU_dd* gpu_c = (GPU_dd*)malloc(sizeof(GPU_dd)*len);
	dd_to_gdd(h_a, gpu_a, len);

	gpu_dd_sin(gpu_a, gpu_c, len);

	gdd_to_dd(gpu_c, ref_c, len);
	free(gpu_a);
	free(gpu_c);

	//CPU
	struct timeval start, end;
	gettimeofday(&start, NULL);
#pragma omp parallel for
	for(int i = 0; i < len; i++) {
		gold_c[i] = sin(h_a[i]);;
	}
	gettimeofday(&end, NULL);
	printf("CPU dd_real, sin: %f sec.\n", getSec(start, end));

	//check result
	checkTwoArray(gold_c, ref_c, len);

	//clean up
	free(h_a);
	free(gold_c);
	free(ref_c);
	printf("\n");
}

void test_sin_qd(const int len) {
	printf("*** sine, quad-double, #element = %d ***\n", len);
	qd_real* h_a = (qd_real*)malloc(sizeof(qd_real)*len);
	qd_real* gold_c = (qd_real*)malloc(sizeof(qd_real)*len);
	qd_real* ref_c = (qd_real*)malloc(sizeof(qd_real)*len);
	for(int i = 0; i < len; i++) {
		rand(&(h_a[i]));
	}

	//GPU
	GPU_qd* gpu_a = (GPU_qd*)malloc(sizeof(GPU_qd)*len);
	GPU_qd* gpu_c = (GPU_qd*)malloc(sizeof(GPU_qd)*len);
	qd_to_gqd(h_a, gpu_a, len);

	gpu_qd_sin(gpu_a, gpu_c, len);

	gqd_to_qd(gpu_c, ref_c, len);
	free(gpu_a);
	free(gpu_c);

	//CPU
	struct timeval start, end;
	gettimeofday(&start, NULL);
#pragma omp parallel for
	for(int i = 0; i < len; i++) {
		gold_c[i] = sin(h_a[i]);
	}
	gettimeofday(&end, NULL);
	printf("CPU qd_real, sin: %f sec\n", getSec(start, end));

	//check result
	checkTwoArray(gold_c, ref_c, len);

	//clean up
	free(h_a);
	free(gold_c);
	free(ref_c);

	printf("\n");
}


int main() {
	unsigned int old_cw;
	fpu_fix_start(&old_cw);
	const int ompNumThread = 4;
	omp_set_num_threads(4);
	const int len = 8000000;

	printf("GQD TEST SUIT");
	printf("len = %d, ompNumThread = %d\n", len, ompNumThread);
	
	
	printf("/////////////// Start test double-double... ///////////////\n");
	GDDInit();
	test_add_dd(len);
	test_mul_dd(len);
	test_div_dd(len);
	test_sqrt_dd(len);
	test_exp_dd(len);
	test_sin_dd(len);
	GDDFree();
	

	
	printf("/////////////// Start test quad-double... //////////////////\n");
	GQDInit();
	test_add_qd(len);	
	test_mul_qd(len);
	test_div_qd(len);
	test_sqrt_qd(len);
	test_exp_qd(len);
	test_sin_qd(len);	
	GQDFree();
	

	fpu_fix_end(&old_cw);

	return 0;
}

