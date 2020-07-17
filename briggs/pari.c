/* A simple interface to PARI */
/* Stolen from AR Booker */
#undef ulong
#include <pari/pari.h>
#include <setjmp.h>

static void sigint_handler(void) { exit(0); }
static jmp_buf env;
static const char *inputstring;
static void error_handler(long numerr) {
	fprintf(stderr,"pari error while processing %s\n",inputstring);
	longjmp(env,numerr);
}

static char *pari(const char *s) {
	static int init;
	char *res;
	pari_sp av;

	if (!init) {
		pari_init(10000000,0);
		cb_pari_err_recover = error_handler;
		cb_pari_sigint = sigint_handler;
		init = 1;
	}

	av = avma;
	inputstring = s;
	if (setjmp(env))
		res = (char *)0;
	else
		res = GENtostr(gp_read_str(s));
	avma = av;
	return res;
}

static long parilong(const char *s) {
	char *buf = pari(s);
	long res;
	if (!buf) return 0;
	if (sscanf(buf,"%ld",&res) != 1) res = 0;
	free(buf);
	return res;
}

static char *paristr(char *s) {
	char *buf = pari(s);
	if (!buf) return (char *)0;
	strcpy(s,buf);
	free(buf);
	return s;
}
