#include <stdio.h>
#include <arb.h>

#define prec 300
#define ystep_str "0.125"
#define nsamples  1600
#define npowerseries 28
#define ntaylor   90
#define nasymp    80
#define ymin_str  "0.5"

static arb_t asymp[nasymp],taylor[nsamples][ntaylor],ystep,ymin,ymax,Euler;
#include "wdata.c"

static void Kinit(void) {
	arb_t t1,t2,temp[ntaylor],w[ntaylor];
	int i,j,l;

	arb_init(t1); arb_init(t2);
	for (l=0;l<ntaylor;l++) {
		arb_init(temp[l]);
		arb_init(w[l]);
	}

	arb_init(ystep); arb_set_str(ystep,ystep_str,prec);
	arb_init(ymin); arb_set_str(ymin,ymin_str,prec);
	arb_init(ymax); arb_mul_ui(ymax,ystep,nsamples,prec);
	arb_add(ymax,ymax,ymin,prec);
	arb_init(Euler); arb_const_euler(Euler,prec);
	for (i=0;i<nsamples;i++)
		for (l=0;l<ntaylor;l++)
			arb_init(taylor[i][l]);
	for (i=0;i<nsamples;i++) {
		for (l=0;l<2;l++)
			arb_set_str(w[l],Wdata[i][l],prec);
		arb_mul_ui(t1,ystep,2*i+1,prec);
		arb_mul_2exp_si(t1,t1,-1);
		arb_add(t1,t1,ymin,prec);
		arb_inv(t1,t1,prec);
		arb_neg(t1,t1); // t1 = -1/y
		arb_mul_2exp_si(t2,t1,-2); // t2 = lambda*t1
		for (j=1;j<ntaylor;j++) {
			arb_mul_ui(temp[j],t2,j-1,prec);
			arb_mul(t2,t2,t1,prec);
		}
		for (;l<ntaylor;l++) {
			arb_set(t1,w[l-2]);
			for (j=l-2;j>=0;j--)
				arb_submul(t1,temp[l-j],w[j],prec);
			arb_div_ui(w[l],t1,l*l-l,prec);
		}
		for (l=0;l<ntaylor;l++)
			arb_set(taylor[i][l],w[l]);
	}

	arb_const_pi(t1,prec);
	arb_mul_2exp_si(t1,t1,-1);
	arb_sqrt(t1,t1,prec); // t1 = sqrt(Pi/2)
	for (j=0;j<nasymp;j++) {
		arb_init(asymp[j]);
		arb_set(asymp[j],t1);
		arb_mul_si(t1,t1,-(1+4*j*(j+1)),prec);
		arb_div_si(t1,t1,8*(j+1),prec);
	}

	arb_clear(t1); arb_clear(t2);
	for (l=0;l<ntaylor;l++) {
		arb_clear(temp[l]);
		arb_clear(w[l]);
	}
}

// compute K_0(y)
void K0(arb_t K,const arb_t yin) {
	static int init;
	static arb_t y,t1,t2,y2j;
	int i,j;

	if (!init) {
		Kinit();
		arb_init(y); arb_init(y2j);
		arb_init(t1); arb_init(t2);
		init = 1;
	}
	arb_set(y,yin);
	arb_sub(t1,y,ymin,prec);
	arb_sub(t2,y,ymax,prec);
	if (arb_contains_negative(t1)) {
		arb_mul_2exp_si(t1,y,-1);
		arb_mul(y,t1,t1,prec); // replace y by (y/2)^2
		arb_log(t1,t1,prec);
		arb_add(t1,t1,Euler,prec);
		arb_neg(t1,t1);

		arb_set(K,t1);
		arb_set_ui(y2j,1);
		for (j=1;j<npowerseries;j++) {
			arb_set_ui(t2,j);
			arb_inv(t2,t2,prec);
			arb_add(t1,t1,t2,prec); // t1 += 1.0/j
			arb_mul(y2j,y2j,y,prec);
			arb_div_ui(y2j,y2j,j*j,prec); // y2j *= y/(j*j);
			arb_addmul(K,t1,y2j,prec); // K += t1*y2j;
		}
	} else if (arb_contains_negative(t2)) {
		arb_div(t2,t1,ystep,prec);
		i = arf_get_si(arb_midref(t2),ARF_RND_FLOOR);

		arb_mul_ui(t2,ystep,2*i+1,prec);
		arb_mul_2exp_si(t2,t2,-1);
		arb_add(t2,t2,ymin,prec);
		arb_sub(t1,y,t2,prec); // t1 = y-((i+0.5)*ystep+ymin);
		arb_zero(K);
		for (j=ntaylor-1;j>=0;j--) {
			arb_mul(t2,K,t1,prec);
			arb_add(K,t2,taylor[i][j],prec); // K = K*t1+taylor[i][j];
		}
		arb_sqrt(t1,y,prec);
		arb_div(K,K,t1,prec);
	} else {
		arb_inv(t1,y,prec);
		arb_zero(K);
		for (j=nasymp-1;j>=0;j--) {
			arb_mul(t2,K,t1,prec);
			arb_add(K,t2,asymp[j],prec);
		}
		arb_exp(t1,y,prec);
		arb_sqrt(t2,y,prec);
		arb_mul(t1,t1,t2,prec);
		arb_div(K,K,t1,prec); // K *= exp(-y)/sqrt(y)
	}
}
