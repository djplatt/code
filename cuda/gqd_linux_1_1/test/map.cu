#ifndef _GQD_TEST_CU_
#define _GQD_TEST_CU_


/* addition */
template<class T> __global__
void map_add_kernel(const T* d_a, const T* d_b, T* d_c, const int numElement) {
        const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;
        const unsigned int delta = blockDim.x*gridDim.x;

        for( int i = index; i < numElement; i += delta ) {
                d_c[i] = d_a[i] + d_b[i];
        }

}

template<class T> __global__
void map_sub_kernel(const T* d_a, const T* d_b, T* d_c, const int numElement) {
        const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;
        const unsigned int delta = blockDim.x*gridDim.x;

        for( int i = index; i < numElement; i += delta ) {
                d_c[i] = d_a[i] - d_b[i];
        }

}

/* multiplication */
template<class T> __global__
void map_mul_kernel(const T* d_a, const T* d_b, T* d_c, const int numElement) {
        const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;
        const unsigned int delta = blockDim.x*gridDim.x;

        for( int i = index; i < numElement; i += delta ) {
                d_c[i] = d_a[i] * d_b[i];
        }

}

/* division */
template<class T> __global__
void map_div_kernel(const T* d_a, const T* d_b, T* d_c, const int numElement) {
        const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;
        const unsigned int delta = blockDim.x*gridDim.x;

        for( int i = index; i < numElement; i += delta ) {
                d_c[i] = d_a[i] / d_b[i];
        }

}

/* sqrt */
template<class T> __global__
void map_sqrt_kernel(const T* d_a, T* d_c, const int numElement) {
        const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;
        const unsigned int delta = blockDim.x*gridDim.x;

        for( int i = index; i < numElement; i += delta ) {
                d_c[i] = sqrt(d_a[i]);
        }

}



/* exp */
template<class T> __global__
void map_exp_kernel(const T* d_a, T* d_c, const int numElement) {
        const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;
        const unsigned int delta = blockDim.x*gridDim.x;

        for( int i = index; i < numElement; i += delta ) {
                d_c[i] = exp(d_a[i]);
        }

}

/* sin */
template<class T> __global__
void map_sin_kernel(const T* d_a, T* d_c, const int numElement, 
			const T* d_sin_table, const T* d_cos_table) {
        const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;
        const unsigned int delta = blockDim.x*gridDim.x;

        for( int i = index; i < numElement; i += delta ) {
                d_c[i] = sin(d_a[i], d_sin_table, d_cos_table);
        }

}


/* map for GPU_qd */
__global__
void map_kernel( const GPU_qd* d_a, const GPU_qd* d_b, GPU_qd* d_c, const int numElement, 
				const GPU_qd* d_sin_table, const GPU_qd* d_cos_table ) {
	const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int delta = blockDim.x*gridDim.x;

	//GPU_qd a, b, c;

	for( int i = index; i < numElement; i += delta ) {
		//d_c[i] = sin( d_a[i], d_sin_table, d_cos_table );
		//d_c[i] = d_a[i] * d_b[i];
		d_c[i] = exp( d_a[i] );

		//a = d_a[i];
		//b = d_b[i];
		//d_c[i] = exp(a);
		//d_c[i] = d_a[i] + d_b[i];
	}
}


/* map for GPU_dd */
__global__
void map_kernel( const GPU_dd* d_a, const GPU_dd* d_b, GPU_dd* d_c, const int numElement ) {
	const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int delta = blockDim.x*gridDim.x;

	for( int i = index; i < numElement; i += delta ) {
		//d_c[i] = sin( d_a[i], d_sin_table, d_cos_table );
		//d_c[i] = d_a[i] / d_b[i];
		d_c[i] = log( d_a[i] );
	}
}


__global__
void map_kernel( const double* d_a, const double* d_b, double* d_c, const int numElement ) {
	const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int delta = blockDim.x*gridDim.x;

	for( int i = index; i < numElement; i += delta ) {
		//d_c[i] = sin( d_a[i], d_sin_table, d_cos_table );
		//d_c[i] = d_a[i] * d_b[i];
		d_c[i] = exp( d_a[i] );
	}
}



#endif


