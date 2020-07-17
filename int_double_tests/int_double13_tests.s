	.file	"int_double13_tests.cpp"
	.text
.globl _Z9_fpu_rnddv
	.type	_Z9_fpu_rnddv, @function
_Z9_fpu_rnddv:
.LFB1564:
#APP
# 80 "../includes/int_double13.0.h" 1
	stmxcsr old__SSE_cw(%rip)
# 0 "" 2
#NO_APP
	movl	old__SSE_cw(%rip), %eax
	andl	$8127, %eax
	orb	$32, %ah
	movl	%eax, new__SSE_cw(%rip)
#APP
# 83 "../includes/int_double13.0.h" 1
	ldmxcsr new__SSE_cw(%rip)
# 0 "" 2
#NO_APP
	movl	$0, %edx
	movq	%rdx, d_zero(%rip)
	movabsq	$-9223372036854775808, %rax
	movq	%rax, d_zero+8(%rip)
	movq	%rdx, d_zero_zero(%rip)
	movq	%rdx, d_zero_zero+8(%rip)
	movq	%rax, d_neg_zero(%rip)
	movq	%rax, d_neg_zero+8(%rip)
	movq	%rdx, d_neg_neg_zero(%rip)
	movq	%rdx, d_neg_neg_zero+8(%rip)
	ret
.LFE1564:
	.size	_Z9_fpu_rnddv, .-_Z9_fpu_rnddv
	.type	_GLOBAL__I_old_cw, @function
_GLOBAL__I_old_cw:
.LFB1589:
	subq	$8, %rsp
.LCFI0:
	movl	$_ZStL8__ioinit, %edi
	call	_ZNSt8ios_base4InitC1Ev
	movl	$__dso_handle, %edx
	movl	$_ZStL8__ioinit, %esi
	movl	$_ZNSt8ios_base4InitD1Ev, %edi
	call	__cxa_atexit
	addq	$8, %rsp
	ret
.LFE1589:
	.size	_GLOBAL__I_old_cw, .-_GLOBAL__I_old_cw
	.section	.ctors,"aw",@progbits
	.align 8
	.quad	_GLOBAL__I_old_cw
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC6:
	.string	"bar    %30.28f,%30.28f\n"
.LC8:
	.string	"foo    %30.28f,%30.28f\n"
.LC10:
	.string	"bletch %30.28f,%30.28f\n"
.LC11:
	.string	"bb     %30.28f,%30.28f\n"
	.text
.globl main
	.type	main, @function
main:
.LFB1587:
	subq	$72, %rsp
.LCFI1:
	call	_Z9_fpu_rnddv
	movsd	.LC2(%rip), %xmm0
	movsd	%xmm0, 48(%rsp)
	movl	$-1717986918, 56(%rsp)
	movl	$-1074685543, 60(%rsp)
	movsd	%xmm0, 32(%rsp)
	movl	$858993459, 40(%rsp)
	movl	$-1074580685, 44(%rsp)
	movapd	d_zero(%rip), %xmm1
	xorpd	.LC5(%rip), %xmm1
	addpd	48(%rsp), %xmm1
	movapd	%xmm1, 16(%rsp)
	movapd	%xmm0, %xmm1
	movl	$.LC6, %edi
	movl	$2, %eax
	call	printf
	movsd	.LC7(%rip), %xmm1
	movsd	32(%rsp), %xmm0
	movl	$.LC8, %edi
	movl	$2, %eax
	call	printf
	movsd	24(%rsp), %xmm1
	xorpd	.LC9(%rip), %xmm1
	movsd	16(%rsp), %xmm0
	movl	$.LC10, %edi
	movl	$2, %eax
	call	printf
	movsd	8(%rsp), %xmm1
	xorpd	.LC9(%rip), %xmm1
	movsd	(%rsp), %xmm0
	movl	$.LC11, %edi
	movl	$2, %eax
	call	printf
	movl	$0, %eax
	addq	$72, %rsp
	ret
.LFE1587:
	.size	main, .-main
.globl old_cw
	.bss
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
.globl d_zero_zero
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
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC2:
	.long	2576980378
	.long	1072798105
	.section	.rodata.cst16,"aM",@progbits,16
	.align 16
.LC5:
	.long	0
	.long	1079574528
	.long	0
	.long	1079574528
	.section	.rodata.cst8
	.align 8
.LC7:
	.long	858993459
	.long	1072902963
	.section	.rodata.cst16
	.align 16
.LC9:
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
	.long	.LFB1564
	.long	.LFE1564-.LFB1564
	.uleb128 0x0
	.align 8
.LEFDE1:
.LSFDE3:
	.long	.LEFDE3-.LASFDE3
.LASFDE3:
	.long	.LASFDE3-.Lframe1
	.long	.LFB1589
	.long	.LFE1589-.LFB1589
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI0-.LFB1589
	.byte	0xe
	.uleb128 0x10
	.align 8
.LEFDE3:
.LSFDE5:
	.long	.LEFDE5-.LASFDE5
.LASFDE5:
	.long	.LASFDE5-.Lframe1
	.long	.LFB1587
	.long	.LFE1587-.LFB1587
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI1-.LFB1587
	.byte	0xe
	.uleb128 0x50
	.align 8
.LEFDE5:
	.ident	"GCC: (GNU) 4.3.3"
	.section	.note.GNU-stack,"",@progbits
