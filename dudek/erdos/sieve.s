	.file	"sieve.cpp"
	.text
	.p2align 4,,15
.globl _Z9cross_outm
	.type	_Z9cross_outm, @function
_Z9cross_outm:
.LFB28:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	imulq	%rdi, %rdi
	movq	st(%rip), %rax
	xorl	%edx, %edx
	divq	%rdi
	movq	%rdi, %rax
	subq	%rdx, %rax
	cmpq	%rdi, %rax
	movq	%rax, %rdx
	jne	.L8
	jmp	.L11
	.p2align 4,,10
	.p2align 3
.L6:
	movb	$1, nsfree(%rdx)
	addq	%rdi, %rdx
.L8:
	cmpq	$1073743504, %rdx
	jbe	.L6
	rep
	ret
.L11:
	xorl	%edx, %edx
	jmp	.L6
	.cfi_endproc
.LFE28:
	.size	_Z9cross_outm, .-_Z9cross_outm
	.section	.rodata.str1.8,"aMS",@progbits,1
	.align 8
.LC0:
	.string	"Usage:- %s <sieve size> <start n> <num its>.\n"
	.align 8
.LC1:
	.string	"<sieve size> = L1 cache size per core in KBytes."
	.align 8
.LC2:
	.string	"Expect start to be 0 mod 4. Exiting."
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC3:
	.string	"Need start >=%lu. Exiting.\n"
.LC5:
	.string	"Failed with n=%lu\n"
	.text
	.p2align 4,,15
.globl main
	.type	main, @function
main:
.LFB29:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	pushq	%r14
	.cfi_def_cfa_offset 16
	.cfi_offset 14, -16
	cmpl	$4, %edi
	pushq	%r13
	.cfi_def_cfa_offset 24
	.cfi_offset 13, -24
	pushq	%r12
	.cfi_def_cfa_offset 32
	.cfi_offset 12, -32
	pushq	%rbp
	.cfi_def_cfa_offset 40
	.cfi_offset 6, -40
	pushq	%rbx
	.cfi_def_cfa_offset 48
	.cfi_offset 3, -48
	movq	%rsi, %rbx
	jne	.L31
	movq	8(%rsi), %rdi
	movl	$10, %edx
	xorl	%esi, %esi
	call	strtol
	movl	%eax, %edi
	call	primesieve_set_sieve_size
	movq	16(%rbx), %rdi
	xorl	%esi, %esi
	movl	$10, %edx
	call	strtol
	movq	24(%rbx), %rdi
	movl	$10, %edx
	xorl	%esi, %esi
	movq	%rax, st1(%rip)
	call	strtol
	movq	st1(%rip), %rdx
	movq	%rax, %rbp
	testb	$3, %dl
	jne	.L32
	cmpq	$1680, %rdx
	jbe	.L33
	leaq	-1681(%rdx), %rax
	addq	$1073741824, %rdx
	testq	%rbp, %rbp
	movq	%rdx, en(%rip)
	movq	%rax, st(%rip)
	je	.L16
	xorl	%ebx, %ebx
	movabsq	$-9223372036854775808, %r12
	.p2align 4,,10
	.p2align 3
.L17:
	xorl	%eax, %eax
	.p2align 4,,10
	.p2align 3
.L18:
	movb	$0, nsfree(%rax)
	addq	$1, %rax
	cmpq	$1073743505, %rax
	jne	.L18
	subq	$1, %rdx
	js	.L19
	cvtsi2sdq	%rdx, %xmm0
.L20:
	sqrtsd	%xmm0, %xmm0
	ucomisd	.LC4(%rip), %xmm0
	jae	.L21
	cvttsd2siq	%xmm0, %rsi
.L22:
	movl	$_Z9cross_outm, %edx
	movl	$2, %edi
	call	primesieve_callback_primes
	movq	st1(%rip), %r14
	movq	en(%rip), %rax
	cmpq	%rax, %r14
	jae	.L23
	movl	$1681, %r13d
	.p2align 4,,10
	.p2align 3
.L28:
	cmpb	$0, nsfree-9(%r13)
	je	.L24
	cmpb	$0, nsfree-25(%r13)
	je	.L24
	cmpb	$0, nsfree-49(%r13)
	je	.L24
	cmpb	$0, nsfree-121(%r13)
	je	.L24
	cmpb	$0, nsfree-169(%r13)
	je	.L24
	cmpb	$0, nsfree-289(%r13)
	je	.L24
	cmpb	$0, nsfree-361(%r13)
	je	.L24
	cmpb	$0, nsfree-529(%r13)
	je	.L24
	cmpb	$0, nsfree-841(%r13)
	je	.L24
	cmpb	$0, nsfree-961(%r13)
	je	.L24
	cmpb	$0, nsfree-1369(%r13)
	je	.L24
	cmpb	$0, nsfree-1681(%r13)
	jne	.L34
	.p2align 4,,10
	.p2align 3
.L24:
	cmpb	$0, nsfree-2(%r13)
	leaq	2(%r14), %rsi
	je	.L25
	cmpb	$0, nsfree-7(%r13)
	je	.L25
	cmpb	$0, nsfree-23(%r13)
	je	.L25
	cmpb	$0, nsfree-47(%r13)
	je	.L25
	cmpb	$0, nsfree-119(%r13)
	je	.L25
	cmpb	$0, nsfree-167(%r13)
	je	.L25
	cmpb	$0, nsfree-287(%r13)
	je	.L25
	cmpb	$0, nsfree-359(%r13)
	je	.L25
	cmpb	$0, nsfree-527(%r13)
	je	.L25
	cmpb	$0, nsfree-839(%r13)
	je	.L25
	cmpb	$0, nsfree-959(%r13)
	je	.L25
	cmpb	$0, nsfree-1367(%r13)
	jne	.L35
	.p2align 4,,10
	.p2align 3
.L25:
	cmpb	$0, nsfree-1(%r13)
	leaq	3(%r14), %rsi
	je	.L26
	cmpb	$0, nsfree-6(%r13)
	je	.L26
	cmpb	$0, nsfree-22(%r13)
	je	.L26
	cmpb	$0, nsfree-46(%r13)
	je	.L26
	cmpb	$0, nsfree-118(%r13)
	je	.L26
	cmpb	$0, nsfree-166(%r13)
	je	.L26
	cmpb	$0, nsfree-286(%r13)
	je	.L26
	cmpb	$0, nsfree-358(%r13)
	je	.L26
	cmpb	$0, nsfree-526(%r13)
	je	.L26
	cmpb	$0, nsfree-838(%r13)
	je	.L26
	cmpb	$0, nsfree-958(%r13)
	je	.L26
	cmpb	$0, nsfree-1366(%r13)
	jne	.L36
	.p2align 4,,10
	.p2align 3
.L26:
	movq	en(%rip), %rax
	addq	$4, %r14
	cmpq	%rax, %r14
	jae	.L27
	addq	$4, %r13
	jmp	.L28
.L27:
	movq	st1(%rip), %r14
.L23:
	leaq	1073741824(%rax), %rdx
	addq	$1, %rbx
	addq	$1073741824, %r14
	addq	$1073741824, st(%rip)
	cmpq	%rbx, %rbp
	movq	%r14, st1(%rip)
	movq	%rdx, en(%rip)
	ja	.L17
.L16:
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
.L21:
	.cfi_restore_state
	subsd	.LC4(%rip), %xmm0
	cvttsd2siq	%xmm0, %rsi
	xorq	%r12, %rsi
	jmp	.L22
.L19:
	movq	%rdx, %rax
	andl	$1, %edx
	shrq	%rax
	orq	%rdx, %rax
	cvtsi2sdq	%rax, %xmm0
	addsd	%xmm0, %xmm0
	jmp	.L20
.L36:
	cmpb	$0, nsfree-1678(%r13)
	je	.L26
	movl	$.LC5, %edi
	xorl	%eax, %eax
	call	printf
	jmp	.L26
.L35:
	cmpb	$0, nsfree-1679(%r13)
	je	.L25
	movl	$.LC5, %edi
	xorl	%eax, %eax
	call	printf
	jmp	.L25
.L34:
	movq	%r14, %rsi
	movl	$.LC5, %edi
	xorl	%eax, %eax
	call	printf
	jmp	.L24
.L33:
	movl	$.LC3, %edi
	movl	$1681, %esi
	xorl	%eax, %eax
	call	printf
	xorl	%edi, %edi
	call	exit
.L31:
	movq	(%rsi), %rsi
	movl	$.LC0, %edi
	xorl	%eax, %eax
	call	printf
	movl	$.LC1, %edi
	call	puts
	xorl	%edi, %edi
	call	exit
.L32:
	movl	$.LC2, %edi
	call	puts
	xorl	%edi, %edi
	call	exit
	.cfi_endproc
.LFE29:
	.size	main, .-main
.globl nsfree
	.bss
	.align 32
	.type	nsfree, @object
	.size	nsfree, 1073743505
nsfree:
	.zero	1073743505
.globl st
	.align 8
	.type	st, @object
	.size	st, 8
st:
	.zero	8
.globl en
	.align 8
	.type	en, @object
	.size	en, 8
en:
	.zero	8
.globl st1
	.align 8
	.type	st1, @object
	.size	st1, 8
st1:
	.zero	8
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC4:
	.long	0
	.long	1138753536
	.ident	"GCC: (GNU) 4.4.7 20120313 (Red Hat 4.4.7-3)"
	.section	.note.GNU-stack,"",@progbits
