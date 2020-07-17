	.file	"int_double13_tests.cpp"
	.text
	.p2align 4,,15
.globl _Z9_fpu_rnddv
	.type	_Z9_fpu_rnddv, @function
_Z9_fpu_rnddv:
.LFB1580:
#APP
# 154 "../includes/int_double13.0.h" 1
	stmxcsr old__SSE_cw(%rip)
# 0 "" 2
#NO_APP
	movl	old__SSE_cw(%rip), %eax
	andl	$8127, %eax
	orb	$32, %ah
	movl	%eax, new__SSE_cw(%rip)
#APP
# 157 "../includes/int_double13.0.h" 1
	ldmxcsr new__SSE_cw(%rip)
# 0 "" 2
#NO_APP
	xorl	%ecx, %ecx
	movabsq	$-9223372036854775808, %rdx
	movq	%rcx, -24(%rsp)
	movq	-24(%rsp), %rax
	movq	%rcx, -16(%rsp)
	movq	%rcx, -40(%rsp)
	movq	%rdx, -32(%rsp)
	movq	%rdx, -56(%rsp)
	movq	%rax, d_zero_zero(%rip)
	movq	-16(%rsp), %rax
	movq	%rdx, -48(%rsp)
	movq	%rdx, -72(%rsp)
	movq	%rcx, -64(%rsp)
	movq	%rax, d_zero_zero+8(%rip)
	movq	-40(%rsp), %rax
	movq	%rax, d_zero(%rip)
	movq	-32(%rsp), %rax
	movq	%rax, d_zero+8(%rip)
	movq	-56(%rsp), %rax
	movq	%rax, d_neg_zero(%rip)
	movq	-48(%rsp), %rax
	movq	%rax, d_neg_zero+8(%rip)
	movq	-72(%rsp), %rax
	movq	%rax, d_neg_neg_zero(%rip)
	movq	-64(%rsp), %rax
	movq	%rax, d_neg_neg_zero+8(%rip)
	ret
.LFE1580:
	.size	_Z9_fpu_rnddv, .-_Z9_fpu_rnddv
	.p2align 4,,15
	.type	_GLOBAL__I_d_zero_zero, @function
_GLOBAL__I_d_zero_zero:
.LFB1586:
	subq	$8, %rsp
.LCFI0:
	movl	$_ZStL8__ioinit, %edi
	call	_ZNSt8ios_base4InitC1Ev
	movl	$__dso_handle, %edx
	movl	$_ZStL8__ioinit, %esi
	movl	$_ZNSt8ios_base4InitD1Ev, %edi
	call	__cxa_atexit
	xorl	%eax, %eax
	movq	%rax, d_zero_zero(%rip)
	movq	%rax, d_zero_zero+8(%rip)
	movq	%rax, d_zero(%rip)
	movq	%rax, d_zero+8(%rip)
	movq	%rax, d_neg_zero(%rip)
	movq	%rax, d_neg_zero+8(%rip)
	movq	%rax, d_neg_neg_zero(%rip)
	movq	%rax, d_neg_neg_zero+8(%rip)
	addq	$8, %rsp
	ret
.LFE1586:
	.size	_GLOBAL__I_d_zero_zero, .-_GLOBAL__I_d_zero_zero
	.section	.ctors,"aw",@progbits
	.align 8
	.quad	_GLOBAL__I_d_zero_zero
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC7:
	.string	"[%30.28e,%30.28e]"
	.text
	.p2align 4,,15
.globl main
	.type	main, @function
main:
.LFB1584:
	pushq	%rbx
.LCFI1:
	subq	$784, %rsp
.LCFI2:
#APP
# 154 "../includes/int_double13.0.h" 1
	stmxcsr old__SSE_cw(%rip)
# 0 "" 2
#NO_APP
	movl	old__SSE_cw(%rip), %eax
	andl	$8127, %eax
	orb	$32, %ah
	movl	%eax, new__SSE_cw(%rip)
#APP
# 157 "../includes/int_double13.0.h" 1
	ldmxcsr new__SSE_cw(%rip)
# 0 "" 2
#NO_APP
	xorl	%ecx, %ecx
	movabsq	$-9223372036854775808, %rdx
	movabsq	$-4616189618054758400, %r8
	movq	%rcx, 48(%rsp)
	movq	48(%rsp), %rax
	movabsq	$4607182418800017408, %rbx
	movq	%rcx, 56(%rsp)
	movq	%rcx, 32(%rsp)
	movabsq	$4158027847206421152, %rdi
	movq	%rdx, 40(%rsp)
	movq	%rdx, 16(%rsp)
	movq	%rax, d_zero_zero(%rip)
	movq	56(%rsp), %rax
	movq	%rdx, 24(%rsp)
	movq	%rdx, (%rsp)
	movabsq	$-5065344189648354656, %rdx
	movq	%r8, 728(%rsp)
	movq	%rcx, 8(%rsp)
	movq	%rax, d_zero_zero+8(%rip)
	movq	32(%rsp), %rax
	movq	728(%rsp), %rcx
	movq	%rbx, 720(%rsp)
	movq	%rdi, 704(%rsp)
	movq	720(%rsp), %rsi
	movl	$.LC7, %edi
	movq	%rax, d_zero(%rip)
	movq	40(%rsp), %rax
	movq	%rdx, 712(%rsp)
	movq	704(%rsp), %rdx
	movq	%rsi, 768(%rsp)
	movq	%rcx, 776(%rsp)
	movq	%rax, d_zero+8(%rip)
	movq	16(%rsp), %rax
	movq	%rax, d_neg_zero(%rip)
	movq	24(%rsp), %rax
	movq	%rax, d_neg_zero+8(%rip)
	movq	(%rsp), %rax
	movq	%rax, d_neg_neg_zero(%rip)
	movq	8(%rsp), %rax
	movq	%rax, d_neg_neg_zero+8(%rip)
	movq	%rdx, 752(%rsp)
	movq	712(%rsp), %rax
	movq	%rsi, 656(%rsp)
	movq	%rcx, 664(%rsp)
	movq	%rdx, 640(%rsp)
	movapd	656(%rsp), %xmm0
	movq	%rax, 648(%rsp)
	movq	%rax, 760(%rsp)
	addpd	640(%rsp), %xmm0
	movapd	%xmm0, 608(%rsp)
	movq	616(%rsp), %rax
	movq	608(%rsp), %rdx
	movq	%rax, 584(%rsp)
	movq	%rdx, 576(%rsp)
	movsd	584(%rsp), %xmm1
	movq	%rax, 600(%rsp)
	movsd	576(%rsp), %xmm0
	movq	%rax, 632(%rsp)
	xorpd	.LC6(%rip), %xmm1
	movq	%rax, 744(%rsp)
	movl	$2, %eax
	movq	%rdx, 592(%rsp)
	movq	%rdx, 624(%rsp)
	movq	%rdx, 736(%rsp)
	call	printf
	movl	$10, %edi
	call	putchar
	movq	768(%rsp), %rax
	movl	$.LC7, %edi
	movq	%rax, 560(%rsp)
	movq	776(%rsp), %rax
	movq	%rax, 568(%rsp)
	movq	752(%rsp), %rax
	movq	%rax, 544(%rsp)
	movq	760(%rsp), %rax
	movq	%rax, 552(%rsp)
	movapd	544(%rsp), %xmm0
	shufpd	$1, %xmm0, %xmm0
	addpd	560(%rsp), %xmm0
	movapd	%xmm0, 512(%rsp)
	movq	520(%rsp), %rax
	movq	512(%rsp), %rdx
	movq	%rax, 488(%rsp)
	movq	%rdx, 480(%rsp)
	movsd	488(%rsp), %xmm1
	movq	%rax, 504(%rsp)
	movsd	480(%rsp), %xmm0
	movq	%rax, 536(%rsp)
	xorpd	.LC6(%rip), %xmm1
	movq	%rax, 744(%rsp)
	movl	$2, %eax
	movq	%rdx, 496(%rsp)
	movq	%rdx, 528(%rsp)
	movq	%rdx, 736(%rsp)
	call	printf
	movl	$10, %edi
	call	putchar
	movq	736(%rsp), %rax
	movl	$.LC7, %edi
	movq	%rax, 464(%rsp)
	movq	744(%rsp), %rax
	movq	%rax, 472(%rsp)
	movapd	464(%rsp), %xmm0
	shufpd	$1, %xmm0, %xmm0
	movapd	%xmm0, 432(%rsp)
	movq	440(%rsp), %rax
	movq	432(%rsp), %rdx
	movq	%rax, 408(%rsp)
	movq	%rdx, 400(%rsp)
	movsd	408(%rsp), %xmm1
	movq	%rax, 424(%rsp)
	movsd	400(%rsp), %xmm0
	movq	%rax, 456(%rsp)
	xorpd	.LC6(%rip), %xmm1
	movq	%rax, 744(%rsp)
	movl	$2, %eax
	movq	%rdx, 416(%rsp)
	movq	%rdx, 448(%rsp)
	movq	%rdx, 736(%rsp)
	call	printf
	movl	$10, %edi
	call	putchar
	movq	768(%rsp), %rax
	movl	$.LC7, %edi
	movapd	d_neg_zero(%rip), %xmm2
	movq	%rax, 384(%rsp)
	movq	776(%rsp), %rax
	movq	%rax, 392(%rsp)
	movq	752(%rsp), %rax
	movapd	384(%rsp), %xmm1
	movq	%rax, 368(%rsp)
	movq	760(%rsp), %rax
	xorpd	d_zero(%rip), %xmm1
	movq	%rax, 376(%rsp)
	movapd	368(%rsp), %xmm3
	movapd	%xmm1, %xmm0
	movapd	%xmm3, %xmm4
	cmpltpd	%xmm2, %xmm0
	shufpd	$1, %xmm3, %xmm4
	movapd	%xmm0, %xmm5
	xorpd	%xmm2, %xmm4
	movapd	%xmm0, %xmm2
	shufpd	$1, %xmm0, %xmm2
	andpd	%xmm4, %xmm5
	andnpd	%xmm3, %xmm0
	movapd	%xmm2, %xmm6
	andpd	%xmm4, %xmm2
	andnpd	%xmm3, %xmm6
	orpd	%xmm0, %xmm5
	movapd	%xmm1, %xmm0
	orpd	%xmm6, %xmm2
	shufpd	$1, %xmm1, %xmm0
	mulpd	%xmm5, %xmm1
	mulpd	%xmm2, %xmm0
	minpd	%xmm0, %xmm1
	movapd	%xmm1, 336(%rsp)
	movq	344(%rsp), %rax
	movq	336(%rsp), %rdx
	movq	%rax, 312(%rsp)
	movq	%rdx, 304(%rsp)
	movsd	312(%rsp), %xmm1
	movq	%rax, 328(%rsp)
	movsd	304(%rsp), %xmm0
	movq	%rax, 360(%rsp)
	xorpd	.LC6(%rip), %xmm1
	movq	%rax, 744(%rsp)
	movl	$2, %eax
	movq	%rdx, 320(%rsp)
	movq	%rdx, 352(%rsp)
	movq	%rdx, 736(%rsp)
	call	printf
	movl	$10, %edi
	call	putchar
	movq	d_zero(%rip), %rax
	movl	$.LC7, %edi
	movq	%rax, 288(%rsp)
	movq	d_zero+8(%rip), %rax
	movsd	288(%rsp), %xmm0
	movq	%rax, 296(%rsp)
	movl	$2, %eax
	movsd	296(%rsp), %xmm1
	xorpd	.LC6(%rip), %xmm1
	call	printf
	movl	$10, %edi
	call	putchar
	movq	d_neg_zero(%rip), %rax
	movl	$.LC7, %edi
	movq	%rax, 272(%rsp)
	movq	d_neg_zero+8(%rip), %rax
	movsd	272(%rsp), %xmm0
	movq	%rax, 280(%rsp)
	movl	$2, %eax
	movsd	280(%rsp), %xmm1
	xorpd	.LC6(%rip), %xmm1
	call	printf
	movl	$10, %edi
	call	putchar
	movq	d_neg_neg_zero(%rip), %rax
	movl	$.LC7, %edi
	movq	%rax, 256(%rsp)
	movq	d_neg_neg_zero+8(%rip), %rax
	movsd	256(%rsp), %xmm0
	movq	%rax, 264(%rsp)
	movl	$2, %eax
	movsd	264(%rsp), %xmm1
	xorpd	.LC6(%rip), %xmm1
	call	printf
	movl	$10, %edi
	call	putchar
	movq	d_zero_zero(%rip), %rax
	movl	$.LC7, %edi
	movq	%rax, 240(%rsp)
	movq	d_zero_zero+8(%rip), %rax
	movsd	240(%rsp), %xmm0
	movq	%rax, 248(%rsp)
	movl	$2, %eax
	movsd	248(%rsp), %xmm1
	xorpd	.LC6(%rip), %xmm1
	call	printf
	movl	$10, %edi
	call	putchar
	movabsq	$-4615739258092021350, %rax
	movq	%rbx, 688(%rsp)
	movq	688(%rsp), %rdx
	movq	%rax, 696(%rsp)
	movq	696(%rsp), %rax
	movl	$.LC7, %edi
	movq	%rdx, 224(%rsp)
	movq	%rdx, 208(%rsp)
	movq	%rax, 232(%rsp)
	movq	%rax, 216(%rsp)
	movapd	224(%rsp), %xmm0
	addpd	208(%rsp), %xmm0
	movapd	%xmm0, 176(%rsp)
	movq	176(%rsp), %rdx
	movq	184(%rsp), %rax
	movq	%rdx, 144(%rsp)
	movq	%rax, 152(%rsp)
	movapd	144(%rsp), %xmm0
	movq	%rdx, 128(%rsp)
	movq	%rax, 136(%rsp)
	addpd	128(%rsp), %xmm0
	movq	%rdx, 160(%rsp)
	movq	%rax, 168(%rsp)
	movq	%rdx, 192(%rsp)
	movq	%rax, 200(%rsp)
	movq	%rdx, 688(%rsp)
	movq	%rax, 696(%rsp)
	movapd	%xmm0, 96(%rsp)
	movq	96(%rsp), %rdx
	movq	104(%rsp), %rax
	movq	%rdx, 80(%rsp)
	movq	%rax, 88(%rsp)
	movq	%rdx, 112(%rsp)
	movq	%rax, 120(%rsp)
	movq	%rdx, 672(%rsp)
	movq	%rax, 680(%rsp)
	movq	%rax, 72(%rsp)
	movq	%rdx, 64(%rsp)
	movl	$2, %eax
	movsd	72(%rsp), %xmm1
	movsd	64(%rsp), %xmm0
	xorpd	.LC6(%rip), %xmm1
	call	printf
	xorl	%eax, %eax
	addq	$784, %rsp
	popq	%rbx
	ret
.LFE1584:
	.size	main, .-main
.globl d_zero_zero
	.bss
	.align 16
	.type	d_zero_zero, @object
	.size	d_zero_zero, 16
d_zero_zero:
	.zero	16
.globl d_zero
	.align 16
	.type	d_zero, @object
	.size	d_zero, 16
d_zero:
	.zero	16
.globl d_neg_zero
	.align 16
	.type	d_neg_zero, @object
	.size	d_neg_zero, 16
d_neg_zero:
	.zero	16
.globl d_neg_neg_zero
	.align 16
	.type	d_neg_neg_zero, @object
	.size	d_neg_neg_zero, 16
d_neg_neg_zero:
	.zero	16
.globl old_cw
	.align 2
	.type	old_cw, @object
	.size	old_cw, 2
old_cw:
	.zero	2
.globl new_cw
	.align 2
	.type	new_cw, @object
	.size	new_cw, 2
new_cw:
	.zero	2
.globl old__SSE_cw
	.align 4
	.type	old__SSE_cw, @object
	.size	old__SSE_cw, 4
old__SSE_cw:
	.zero	4
.globl new__SSE_cw
	.align 4
	.type	new__SSE_cw, @object
	.size	new__SSE_cw, 4
new__SSE_cw:
	.zero	4
	.local	_ZStL8__ioinit
	.comm	_ZStL8__ioinit,1,1
	.weakref	_ZL20__gthrw_pthread_oncePiPFvvE,pthread_once
	.weakref	_ZL27__gthrw_pthread_getspecificj,pthread_getspecific
	.weakref	_ZL27__gthrw_pthread_setspecificjPKv,pthread_setspecific
	.weakref	_ZL22__gthrw_pthread_createPmPK14pthread_attr_tPFPvS3_ES3_,pthread_create
	.weakref	_ZL22__gthrw_pthread_cancelm,pthread_cancel
	.weakref	_ZL26__gthrw_pthread_mutex_lockP15pthread_mutex_t,pthread_mutex_lock
	.weakref	_ZL29__gthrw_pthread_mutex_trylockP15pthread_mutex_t,pthread_mutex_trylock
	.weakref	_ZL28__gthrw_pthread_mutex_unlockP15pthread_mutex_t,pthread_mutex_unlock
	.weakref	_ZL26__gthrw_pthread_mutex_initP15pthread_mutex_tPK19pthread_mutexattr_t,pthread_mutex_init
	.weakref	_ZL30__gthrw_pthread_cond_broadcastP14pthread_cond_t,pthread_cond_broadcast
	.weakref	_ZL25__gthrw_pthread_cond_waitP14pthread_cond_tP15pthread_mutex_t,pthread_cond_wait
	.weakref	_ZL26__gthrw_pthread_key_createPjPFvPvE,pthread_key_create
	.weakref	_ZL26__gthrw_pthread_key_deletej,pthread_key_delete
	.weakref	_ZL30__gthrw_pthread_mutexattr_initP19pthread_mutexattr_t,pthread_mutexattr_init
	.weakref	_ZL33__gthrw_pthread_mutexattr_settypeP19pthread_mutexattr_ti,pthread_mutexattr_settype
	.weakref	_ZL33__gthrw_pthread_mutexattr_destroyP19pthread_mutexattr_t,pthread_mutexattr_destroy
	.section	.rodata.cst16,"aM",@progbits,16
	.align 16
.LC6:
	.long	0
	.long	-2147483648
	.long	0
	.long	0
	.section	.eh_frame,"a",@progbits
.Lframe1:
	.long	.LECIE1-.LSCIE1
.LSCIE1:
	.long	0x0
	.byte	0x1
.globl __gxx_personality_v0
	.string	"zPR"
	.uleb128 0x1
	.sleb128 -8
	.byte	0x10
	.uleb128 0x6
	.byte	0x3
	.long	__gxx_personality_v0
	.byte	0x3
	.byte	0xc
	.uleb128 0x7
	.uleb128 0x8
	.byte	0x90
	.uleb128 0x1
	.align 8
.LECIE1:
.LSFDE1:
	.long	.LEFDE1-.LASFDE1
.LASFDE1:
	.long	.LASFDE1-.Lframe1
	.long	.LFB1580
	.long	.LFE1580-.LFB1580
	.uleb128 0x0
	.align 8
.LEFDE1:
.LSFDE3:
	.long	.LEFDE3-.LASFDE3
.LASFDE3:
	.long	.LASFDE3-.Lframe1
	.long	.LFB1586
	.long	.LFE1586-.LFB1586
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI0-.LFB1586
	.byte	0xe
	.uleb128 0x10
	.align 8
.LEFDE3:
.LSFDE5:
	.long	.LEFDE5-.LASFDE5
.LASFDE5:
	.long	.LASFDE5-.Lframe1
	.long	.LFB1584
	.long	.LFE1584-.LFB1584
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI1-.LFB1584
	.byte	0xe
	.uleb128 0x10
	.byte	0x4
	.long	.LCFI2-.LCFI1
	.byte	0xe
	.uleb128 0x320
	.byte	0x83
	.uleb128 0x2
	.align 8
.LEFDE5:
	.ident	"GCC: (GNU) 4.3.3"
	.section	.note.GNU-stack,"",@progbits
