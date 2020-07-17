	.file	"isprime.c"
	.section	.rodata.str1.8,"aMS",@progbits,1
	.align 8
.LC0:
	.string	"There were %lu probable primes in [%ld,%ld]\n"
	.section	.text.startup,"ax",@progbits
	.p2align 4,,15
	.globl	main
	.type	main, @function
main:
.LFB12:
	.cfi_startproc
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	xorl	%esi, %esi
	movabsq	$4611686018327387903, %r9
	xorl	%r10d, %r10d
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	movq	%rsi, %r14
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	movl	$100000001, %r13d
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$40, %rsp
	.cfi_def_cfa_offset 96
	.p2align 4,,7
	.p2align 3
.L2:
	movabsq	$-6148914691236517205, %rax
	movq	%r9, %rcx
	mulq	%r9
	movq	%rdx, 8(%rsp)
	movq	%rax, (%rsp)
	movq	8(%rsp), %rax
	shrq	%rax
	leaq	(%rax,%rax,2), %rax
	cmpq	%rax, %r9
	je	.L17
	movabsq	$-3689348814741910323, %rax
	mulq	%r9
	movq	%rdx, 8(%rsp)
	movq	%rax, (%rsp)
	movq	8(%rsp), %rax
	shrq	$2, %rax
	leaq	(%rax,%rax,4), %rax
	cmpq	%rax, %r9
	je	.L17
	xorl	%edx, %edx
	movl	$7, %esi
	movq	%r9, %rax
	divq	%rsi
	testq	%rdx, %rdx
	je	.L17
	movabsq	$3353953467947191203, %rax
	mulq	%r9
	movq	%rdx, 8(%rsp)
	movq	%rax, (%rsp)
	movq	8(%rsp), %rax
	shrq	%rax
	leaq	(%rax,%rax,4), %rdx
	leaq	(%rax,%rdx,2), %rax
	cmpq	%rax, %r9
	je	.L17
	movabsq	$5675921253449092805, %rax
	mulq	%r9
	movq	%rdx, 8(%rsp)
	movq	%rax, (%rsp)
	movq	8(%rsp), %rax
	shrq	$2, %rax
	leaq	(%rax,%rax,2), %rdx
	leaq	(%rax,%rdx,4), %rax
	cmpq	%rax, %r9
	je	.L17
	movabsq	$-1085102592571150095, %rax
	mulq	%r9
	movq	%rdx, 8(%rsp)
	movq	%rax, (%rsp)
	movq	8(%rsp), %rax
	shrq	$4, %rax
	movq	%rax, %rdx
	salq	$4, %rdx
	addq	%rdx, %rax
	cmpq	%rax, %r9
	je	.L17
	movabsq	$-2912643801112034465, %rax
	mulq	%r9
	movq	%rdx, 8(%rsp)
	movq	%rax, (%rsp)
	movq	8(%rsp), %rax
	shrq	$4, %rax
	leaq	(%rax,%rax,8), %rdx
	leaq	(%rax,%rdx,2), %rax
	cmpq	%rax, %r9
	je	.L17
	xorl	%edx, %edx
	movb	$23, %sil
	movq	%r9, %rax
	divq	%rsi
	testq	%rdx, %rdx
	je	.L17
	xorl	%edx, %edx
	movb	$29, %sil
	movq	%r9, %rax
	divq	%rsi
	testq	%rdx, %rdx
	je	.L17
	xorl	%edx, %edx
	movb	$31, %sil
	movq	%r9, %rax
	divq	%rsi
	testq	%rdx, %rdx
	je	.L17
	movabsq	$-2492803253203993461, %rax
	mulq	%r9
	movq	%rdx, 8(%rsp)
	movq	%rax, (%rsp)
	movq	8(%rsp), %rax
	shrq	$5, %rax
	leaq	(%rax,%rax,8), %rdx
	leaq	(%rax,%rdx,4), %rax
	cmpq	%rax, %r9
	je	.L17
	movabsq	$-4049285284472828403, %rax
	mulq	%r9
	movq	%rdx, 8(%rsp)
	movq	%rax, (%rsp)
	movq	8(%rsp), %rax
	shrq	$5, %rax
	leaq	(%rax,%rax,4), %rdx
	leaq	(%rax,%rdx,8), %rax
	cmpq	%rax, %r9
	je	.L17
	movabsq	$-4718934530483838785, %rax
	mulq	%r9
	movq	%rdx, 8(%rsp)
	movq	%rax, (%rsp)
	movq	8(%rsp), %rax
	shrq	$5, %rax
	leaq	(%rax,%rax,4), %rdx
	leaq	(%rax,%rdx,4), %rdx
	leaq	(%rax,%rdx,2), %rax
	cmpq	%rax, %r9
	je	.L17
	movl	%r9d, %edi
	movl	%r9d, %esi
	imull	%r13d, %edi
	movl	%r9d, %edx
	movl	%r9d, %eax
	movq	%r9, %r8
	movq	%r9, %rbp
	shrq	%r8
	addl	$2, %edi
	xorl	%r11d, %r11d
	imull	%r13d, %edi
	imull	%edi, %esi
	addl	$2, %esi
	imull	%edi, %esi
	imull	%esi, %edx
	addl	$2, %edx
	imull	%esi, %edx
	imull	%edx, %eax
	addl	$2, %eax
	imull	%edx, %eax
	xorl	%edx, %edx
	cltq
	movq	%rax, %rbx
	imulq	%r9, %rbx
	addq	$2, %rbx
	imulq	%rax, %rbx
	orq	$-1, %rax
	divq	%r9
	leaq	1(%rdx), %r12
	subq	%r12, %rbp
	testb	$1, %r8b
	jne	.L16
	.p2align 4,,7
	.p2align 3
.L51:
	shrq	%r8
	incl	%r11d
	testb	$1, %r8b
	je	.L51
.L16:
	leaq	(%r12,%r12), %rdi
	cmpq	%rcx, %rdi
	jb	.L4
	subq	%rcx, %rdi
.L4:
	shrq	%r8
	je	.L5
	movq	%rdi, %rsi
	.p2align 4,,7
	.p2align 3
.L8:
	movq	%rsi, %rax
	mulq	%rsi
	movq	%r10, %rsi
	movq	%rax, (%rsp)
	movq	%rdx, 8(%rsp)
	movq	(%rsp), %rax
	imulq	%rbx, %rax
	imulq	%rax, %rsi
	mulq	%r9
	addq	%rsi, %rdx
	addq	(%rsp), %rax
	adcq	8(%rsp), %rdx
	movq	%rax, 16(%rsp)
	movq	%rdx, 24(%rsp)
	movq	24(%rsp), %rsi
	cmpq	%rcx, %rsi
	jb	.L6
	subq	%rcx, %rsi
.L6:
	testb	$1, %r8b
	je	.L7
	movq	%rsi, %rax
	mulq	%rdi
	movq	%r10, %rdi
	movq	%rax, (%rsp)
	movq	%rdx, 8(%rsp)
	movq	(%rsp), %rax
	imulq	%rbx, %rax
	imulq	%rax, %rdi
	mulq	%r9
	addq	%rdi, %rdx
	addq	(%rsp), %rax
	adcq	8(%rsp), %rdx
	movq	%rax, 16(%rsp)
	movq	%rdx, 24(%rsp)
	movq	24(%rsp), %rdi
	cmpq	%rcx, %rdi
	jb	.L7
	subq	%rcx, %rdi
.L7:
	shrq	%r8
	jne	.L8
.L5:
	cmpq	%rdi, %rbp
	je	.L9
	cmpq	%rdi, %r12
	je	.L9
	.p2align 4,,7
	.p2align 3
.L71:
	decl	%r11d
	js	.L17
	movq	%rdi, %rax
	movq	%r10, %rsi
	mulq	%rdi
	movq	%rax, (%rsp)
	movq	%rdx, 8(%rsp)
	movq	(%rsp), %rax
	imulq	%rbx, %rax
	imulq	%rax, %rsi
	mulq	%r9
	addq	%rsi, %rdx
	addq	(%rsp), %rax
	adcq	8(%rsp), %rdx
	movq	%rax, 16(%rsp)
	movq	%rdx, 24(%rsp)
	movq	24(%rsp), %rdi
	cmpq	%rcx, %rdi
	jb	.L11
	subq	%rcx, %rdi
.L11:
	cmpq	%rdi, %rbp
	jne	.L71
.L9:
	incq	%r14
.L17:
	subl	$2, %r13d
	addq	$2, %r9
	adcq	$0, %r10
	movabsq	$4611686018427387905, %rax
	movq	%r10, %rdx
	xorq	%r9, %rax
	orq	%rax, %rdx
	jne	.L2
	movq	%r14, %rsi
	movabsq	$4611686018427387903, %rcx
	movabsq	$4611686018327387903, %rdx
	movl	$.LC0, %edi
	xorl	%eax, %eax
	call	printf
	addq	$40, %rsp
	.cfi_def_cfa_offset 56
	xorl	%eax, %eax
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE12:
	.size	main, .-main
	.ident	"GCC: (GNU) 4.7.0 20110509 (experimental) 4.7.0"
	.section	.note.GNU-stack,"",@progbits
