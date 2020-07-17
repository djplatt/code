	.file	"quint1.c"
	.text
	.p2align 4,,15
.globl _Z9issquare2m
	.type	_Z9issquare2m, @function
_Z9issquare2m:
.LFB22:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	movl	%edi, %eax
	andl	$15, %eax
	cmpl	$9, %eax
	jg	.L2
	leal	-2(%rax), %edx
	cmpl	$1, %edx
	jbe	.L2
	cmpl	$5, %eax
	je	.L2
	cmpl	$7, %eax
	je	.L2
	cmpl	$6, %eax
	.p2align 4,,2
	je	.L2
	cmpl	$8, %eax
	.p2align 4,,2
	je	.L2
	movq	%rdi, %rax
	movabsq	$2635249153387078803, %rdx
	mulq	%rdx
	movq	%rdi, %rax
	subq	%rdx, %rax
	shrq	%rax
	addq	%rax, %rdx
	shrq	$2, %rdx
	leaq	0(,%rdx,8), %rax
	subq	%rdx, %rax
	movq	%rdi, %rdx
	subq	%rax, %rdx
	popcntq	%rdx, %rax
	cmpl	$1, %eax
	je	.L10
	.p2align 4,,10
	.p2align 3
.L2:
	xorl	%eax, %eax
	ret
	.p2align 4,,10
	.p2align 3
.L10:
	testq	%rdi, %rdi
	js	.L3
	cvtsi2sdq	%rdi, %xmm0
.L4:
	movsd	.LC0(%rip), %xmm1
	sqrtsd	%xmm0, %xmm0
	ucomisd	%xmm1, %xmm0
	jae	.L5
	cvttsd2siq	%xmm0, %rax
.L6:
	imulq	%rax, %rax
	cmpq	%rdi, %rax
	sete	%al
	movzbl	%al, %eax
	ret
	.p2align 4,,10
	.p2align 3
.L5:
	subsd	%xmm1, %xmm0
	movabsq	$-9223372036854775808, %rdx
	cvttsd2siq	%xmm0, %rax
	xorq	%rdx, %rax
	jmp	.L6
.L3:
	movq	%rdi, %rax
	movq	%rdi, %rdx
	shrq	%rax
	andl	$1, %edx
	orq	%rdx, %rax
	cvtsi2sdq	%rax, %xmm0
	addsd	%xmm0, %xmm0
	jmp	.L4
	.cfi_endproc
.LFE22:
	.size	_Z9issquare2m, .-_Z9issquare2m
	.p2align 4,,15
.globl _Z6isqrt2m
	.type	_Z6isqrt2m, @function
_Z6isqrt2m:
.LFB20:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	testq	%rdi, %rdi
	js	.L12
	cvtsi2sdq	%rdi, %xmm0
.L13:
	sqrtsd	%xmm0, %xmm0
	call	floor
	movsd	.LC0(%rip), %xmm1
	ucomisd	%xmm1, %xmm0
	jae	.L14
	cvttsd2siq	%xmm0, %rax
	addq	$8, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L14:
	.cfi_restore_state
	subsd	%xmm1, %xmm0
	movabsq	$-9223372036854775808, %rdx
	addq	$8, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 8
	cvttsd2siq	%xmm0, %rax
	xorq	%rdx, %rax
	ret
	.p2align 4,,10
	.p2align 3
.L12:
	.cfi_restore_state
	movq	%rdi, %rax
	andl	$1, %edi
	shrq	%rax
	orq	%rdi, %rax
	cvtsi2sdq	%rax, %xmm0
	addsd	%xmm0, %xmm0
	jmp	.L13
	.cfi_endproc
.LFE20:
	.size	_Z6isqrt2m, .-_Z6isqrt2m
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC1:
	.string	"%lu %lu %lu\n"
	.text
	.p2align 4,,15
.globl main
	.type	main, @function
main:
.LFB23:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	pushq	%r14
	.cfi_def_cfa_offset 16
	.cfi_offset 14, -16
	movabsq	$-9223372036854775808, %r14
	pushq	%r13
	.cfi_def_cfa_offset 24
	.cfi_offset 13, -24
	movabsq	$2635249153387078803, %r13
	pushq	%r12
	.cfi_def_cfa_offset 32
	.cfi_offset 12, -32
	movl	$244, %r12d
	pushq	%rbp
	.cfi_def_cfa_offset 40
	.cfi_offset 6, -40
	movl	$730, %ebp
	pushq	%rbx
	.cfi_def_cfa_offset 48
	.cfi_offset 3, -48
	movl	$243, %ebx
	jmp	.L27
	.p2align 4,,10
	.p2align 3
.L18:
	addl	$3, %ebp
	addq	$1, %r12
	addq	$1, %rbx
	cmpl	$-1073741820, %ebp
	je	.L30
.L27:
	movl	%ebp, %eax
	leaq	1(%rbx,%rbx,2), %rcx
	andl	$15, %eax
	cmpl	$9, %eax
	jg	.L18
	leal	-2(%rax), %edx
	cmpl	$1, %edx
	jbe	.L18
	cmpl	$5, %eax
	je	.L18
	cmpl	$7, %eax
	je	.L18
	cmpl	$6, %eax
	.p2align 4,,2
	je	.L18
	cmpl	$8, %eax
	.p2align 4,,2
	je	.L18
	movq	%rcx, %rax
	mulq	%r13
	movq	%rcx, %rax
	subq	%rdx, %rax
	shrq	%rax
	addq	%rax, %rdx
	shrq	$2, %rdx
	leaq	0(,%rdx,8), %rax
	subq	%rdx, %rax
	movq	%rcx, %rdx
	subq	%rax, %rdx
	popcntq	%rdx, %rax
	cmpl	$1, %eax
	jne	.L18
	testq	%rcx, %rcx
	js	.L19
	cvtsi2sdq	%rcx, %xmm0
.L20:
	sqrtsd	%xmm0, %xmm0
	ucomisd	.LC0(%rip), %xmm0
	jae	.L21
	cvttsd2siq	%xmm0, %rax
.L22:
	imulq	%rax, %rax
	cmpq	%rcx, %rax
	jne	.L18
	movl	%r12d, %eax
	andl	$15, %eax
	cmpl	$9, %eax
	jg	.L18
	leal	-2(%rax), %edx
	cmpl	$1, %edx
	jbe	.L18
	cmpl	$5, %eax
	je	.L18
	cmpl	$7, %eax
	je	.L18
	cmpl	$6, %eax
	.p2align 4,,2
	je	.L18
	cmpl	$8, %eax
	.p2align 4,,2
	je	.L18
	movq	%r12, %rax
	mulq	%r13
	movq	%r12, %rax
	subq	%rdx, %rax
	shrq	%rax
	addq	%rax, %rdx
	shrq	$2, %rdx
	leaq	0(,%rdx,8), %rax
	subq	%rdx, %rax
	movq	%r12, %rdx
	subq	%rax, %rdx
	popcntq	%rdx, %rax
	cmpl	$1, %eax
	jne	.L18
	testq	%r12, %r12
	js	.L23
	cvtsi2sdq	%r12, %xmm0
.L24:
	sqrtsd	%xmm0, %xmm0
	ucomisd	.LC0(%rip), %xmm0
	jae	.L25
	cvttsd2siq	%xmm0, %rax
.L26:
	imulq	%rax, %rax
	cmpq	%rax, %r12
	jne	.L18
	movq	%rbx, %rcx
	movl	$3, %edx
	movl	$1, %esi
	movl	$.LC1, %edi
	xorl	%eax, %eax
	call	printf
	jmp	.L18
	.p2align 4,,10
	.p2align 3
.L21:
	subsd	.LC0(%rip), %xmm0
	cvttsd2siq	%xmm0, %rax
	xorq	%r14, %rax
	jmp	.L22
	.p2align 4,,10
	.p2align 3
.L30:
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 40
	popq	%rbp
	.cfi_def_cfa_offset 32
	popq	%r12
	.cfi_def_cfa_offset 24
	popq	%r13
	.cfi_def_cfa_offset 16
	xorl	%eax, %eax
	popq	%r14
	.cfi_def_cfa_offset 8
	ret
.L25:
	.cfi_restore_state
	subsd	.LC0(%rip), %xmm0
	cvttsd2siq	%xmm0, %rax
	xorq	%r14, %rax
	jmp	.L26
.L19:
	movq	%rcx, %rax
	movq	%rcx, %rdx
	shrq	%rax
	andl	$1, %edx
	orq	%rdx, %rax
	cvtsi2sdq	%rax, %xmm0
	addsd	%xmm0, %xmm0
	jmp	.L20
.L23:
	movq	%r12, %rax
	movq	%r12, %rdx
	shrq	%rax
	andl	$1, %edx
	orq	%rdx, %rax
	cvtsi2sdq	%rax, %xmm0
	addsd	%xmm0, %xmm0
	jmp	.L24
	.cfi_endproc
.LFE23:
	.size	main, .-main
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC0:
	.long	0
	.long	1138753536
	.ident	"GCC: (GNU) 4.4.7 20120313 (Red Hat 4.4.7-3)"
	.section	.note.GNU-stack,"",@progbits
