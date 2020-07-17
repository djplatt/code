	.file	"do_t.cpp"
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC0:
	.string	"Command line:- %s"
.LC1:
	.string	" %s"
.LC2:
	.string	"Usage:- %s <t> <B>\n"
.LC3:
	.string	"1: %ld %ld %ld %ld\n"
.LC4:
	.string	"2: %ld %ld %ld %ld\n"
.LC5:
	.string	"3: %ld %ld %ld %ld\n"
.LC6:
	.string	"4: %ld %ld %ld %ld\n"
	.text
	.p2align 4,,15
.globl main
	.type	main, @function
main:
.LFB25:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	xorl	%eax, %eax
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	movl	%edi, %r13d
	movl	$.LC0, %edi
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	movslq	%r13d, %r12
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	movq	%rsi, %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$168, %rsp
	.cfi_def_cfa_offset 224
	movq	(%rsi), %rsi
	call	printf
	cmpq	$1, %r12
	jbe	.L2
	movl	$1, %ebx
	.p2align 4,,10
	.p2align 3
.L3:
	movq	0(%rbp,%rbx,8), %rsi
	xorl	%eax, %eax
	movl	$.LC1, %edi
	addq	$1, %rbx
	call	printf
	cmpq	%r12, %rbx
	jb	.L3
.L2:
	movl	$10, %edi
	call	putchar
	cmpl	$3, %r13d
	je	.L4
	movq	0(%rbp), %rsi
	movl	$.LC2, %edi
	xorl	%eax, %eax
	call	printf
	movl	$1, %eax
.L5:
	addq	$168, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
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
.L4:
	.cfi_restore_state
	movq	8(%rbp), %rdi
	xorl	%esi, %esi
	movl	$10, %edx
	call	strtol
	movq	16(%rbp), %rdi
	movq	%rax, %rbx
	xorl	%esi, %esi
	movl	$10, %edx
	movq	%rax, 40(%rsp)
	imulq	%rax, %rbx
	call	strtol
	testq	%rax, %rax
	movq	%rax, 32(%rsp)
	je	.L6
	imulq	$3, 40(%rsp), %rax
	leaq	(%rbx,%rbx), %rdx
	movq	$-1, 96(%rsp)
	movl	$1, %r14d
	movq	%rdx, 144(%rsp)
	imulq	$-2, %rbx, %rdx
	movq	%rax, 152(%rsp)
	addq	$2, %rax
	movq	%rax, 88(%rsp)
	movq	152(%rsp), %rax
	movq	%rdx, 136(%rsp)
	movq	152(%rsp), %rdx
	addq	$7, %rax
	imulq	%rbx, %rax
	subq	$3, %rdx
	imulq	%rbx, %rdx
	movq	%rax, 128(%rsp)
	movq	152(%rsp), %rax
	movq	%rdx, 112(%rsp)
	addq	$3, %rax
	imulq	%rbx, %rax
	movq	%rax, 120(%rsp)
	movq	152(%rsp), %rax
	subq	$7, %rax
	imulq	%rbx, %rax
	movq	%rax, 104(%rsp)
	leaq	0(,%rbx,4), %rax
	leaq	(%rax,%rbx), %rdx
	negq	%rax
	subq	%rbx, %rax
	movq	%rdx, 56(%rsp)
	movq	%rax, 48(%rsp)
.L7:
	movq	152(%rsp), %rdx
	leaq	(%r14,%r14), %rax
	movq	88(%rsp), %rbp
	movq	104(%rsp), %r15
	movq	112(%rsp), %r9
	movq	%r14, %r13
	movq	120(%rsp), %r10
	movq	128(%rsp), %r11
	movl	$1, %r12d
	subq	%rax, %rdx
	addq	$5, %rbp
	movq	%rdx, %rax
	addq	88(%rsp), %rax
	movq	96(%rsp), %rdx
	salq	$2, %rdx
	movq	%rax, 64(%rsp)
	movq	88(%rsp), %rax
	movq	%rdx, 80(%rsp)
	addq	%rax, %rax
	movq	%rax, 72(%rsp)
	.p2align 4,,10
	.p2align 3
.L29:
	testq	%r12, %r12
	movq	40(%rsp), %rbx
	je	.L9
	movq	%rbx, %rax
	movq	%r12, %rbx
	jmp	.L10
	.p2align 4,,10
	.p2align 3
.L38:
	movq	%rbx, %rax
	movq	%rdx, %rbx
.L10:
	xorl	%edx, %edx
	divq	%rbx
	testq	%rdx, %rdx
	jne	.L38
.L9:
	testq	%r14, %r14
	je	.L11
	movq	%r14, %rcx
	movq	%rbx, %rax
	jmp	.L13
	.p2align 4,,10
	.p2align 3
.L39:
	movq	%rdx, %rcx
.L13:
	xorl	%edx, %edx
	divq	%rcx
	movq	%rcx, %rax
	testq	%rdx, %rdx
	jne	.L39
	movq	%rcx, %rbx
.L11:
	cmpq	64(%rsp), %rbp
	je	.L14
	movq	%r15, %rdx
	movq	%r15, %rax
	sarq	$63, %rdx
	idivq	%r13
	testq	%rdx, %rdx
	movq	%rax, %rcx
	jne	.L14
	sarq	$63, %rax
	movq	%rax, %rsi
	xorq	%rcx, %rsi
	subq	%rax, %rsi
	cmpq	32(%rsp), %rsi
	jg	.L14
	testq	%rbx, %rbx
	je	.L15
	movq	%rsi, %rax
	movq	%rbx, %rsi
	jmp	.L16
	.p2align 4,,10
	.p2align 3
.L40:
	movq	%rsi, %rax
	movq	%rdx, %rsi
.L16:
	xorl	%edx, %edx
	divq	%rsi
	testq	%rdx, %rdx
	jne	.L40
.L15:
	cmpq	$1, %rsi
	je	.L41
	.p2align 4,,10
	.p2align 3
.L14:
	cmpq	72(%rsp), %rbp
	je	.L17
	movq	%r9, %rdx
	movq	%r9, %rax
	sarq	$63, %rdx
	idivq	%r13
	testq	%rdx, %rdx
	movq	%rax, %rsi
	jne	.L17
	movq	%rax, %rdx
	sarq	$63, %rdx
	movq	%rdx, %rax
	xorq	%rsi, %rax
	subq	%rdx, %rax
	cmpq	32(%rsp), %rax
	jg	.L17
	testq	%rbx, %rbx
	je	.L18
	movq	%rbx, %rcx
	jmp	.L20
	.p2align 4,,10
	.p2align 3
.L42:
	movq	%rcx, %rax
	movq	%rdx, %rcx
.L20:
	xorl	%edx, %edx
	divq	%rcx
	testq	%rdx, %rdx
	jne	.L42
	movq	%rcx, %rax
.L18:
	cmpq	$1, %rax
	je	.L43
	.p2align 4,,10
	.p2align 3
.L17:
	movq	80(%rsp), %rax
	leaq	0(%rbp,%rax), %rax
	testq	%rax, %rax
	je	.L21
	movq	%r10, %rdx
	movq	%r10, %rax
	sarq	$63, %rdx
	idivq	%r13
	testq	%rdx, %rdx
	movq	%rax, %rsi
	jne	.L21
	movq	%rax, %rdx
	sarq	$63, %rdx
	movq	%rdx, %rax
	xorq	%rsi, %rax
	subq	%rdx, %rax
	cmpq	32(%rsp), %rax
	jg	.L21
	testq	%rbx, %rbx
	je	.L22
	movq	%rbx, %rcx
	jmp	.L24
	.p2align 4,,10
	.p2align 3
.L44:
	movq	%rcx, %rax
	movq	%rdx, %rcx
.L24:
	xorl	%edx, %edx
	divq	%rcx
	testq	%rdx, %rdx
	jne	.L44
	movq	%rcx, %rax
.L22:
	cmpq	$1, %rax
	je	.L45
	.p2align 4,,10
	.p2align 3
.L21:
	testq	%rbp, %rbp
	je	.L25
	movq	%r11, %rdx
	movq	%r11, %rax
	sarq	$63, %rdx
	idivq	%r13
	testq	%rdx, %rdx
	movq	%rax, %rcx
	jne	.L25
	movq	%rax, %rdx
	sarq	$63, %rdx
	movq	%rdx, %rax
	xorq	%rcx, %rax
	subq	%rdx, %rax
	cmpq	%rax, 32(%rsp)
	jl	.L25
	testq	%rbx, %rbx
	jne	.L33
	jmp	.L26
	.p2align 4,,10
	.p2align 3
.L46:
	movq	%rbx, %rax
	movq	%rdx, %rbx
.L33:
	xorl	%edx, %edx
	divq	%rbx
	testq	%rdx, %rdx
	jne	.L46
	movq	%rbx, %rax
.L26:
	cmpq	$1, %rax
	je	.L47
	.p2align 4,,10
	.p2align 3
.L25:
	addq	$1, %r12
	addq	$5, %rbp
	addq	56(%rsp), %r11
	addq	%r14, %r13
	addq	56(%rsp), %r10
	addq	48(%rsp), %r9
	addq	48(%rsp), %r15
	cmpq	%r12, 32(%rsp)
	jae	.L29
	movq	144(%rsp), %rdx
	movq	136(%rsp), %rax
	addq	$1, %r14
	subq	$1, 96(%rsp)
	addq	$2, 88(%rsp)
	addq	%rdx, 128(%rsp)
	addq	%rax, 120(%rsp)
	addq	%rdx, 112(%rsp)
	addq	%rax, 104(%rsp)
	cmpq	%r14, 32(%rsp)
	jae	.L7
.L6:
	xorl	%eax, %eax
	jmp	.L5
.L47:
	movq	40(%rsp), %r8
	movq	96(%rsp), %rsi
	movq	%r12, %rdx
	negq	%rdx
	movl	$.LC6, %edi
	xorl	%eax, %eax
	movq	%r9, 16(%rsp)
	movq	%r10, 8(%rsp)
	movq	%r11, 24(%rsp)
	call	printf
	movq	24(%rsp), %r11
	movq	8(%rsp), %r10
	movq	16(%rsp), %r9
	jmp	.L25
.L41:
	movq	40(%rsp), %r8
	movq	%r12, %rdx
	movq	%r14, %rsi
	movl	$.LC3, %edi
	xorl	%eax, %eax
	movq	%r9, 16(%rsp)
	movq	%r10, 8(%rsp)
	movq	%r11, 24(%rsp)
	call	printf
	movq	24(%rsp), %r11
	movq	8(%rsp), %r10
	movq	16(%rsp), %r9
	jmp	.L14
.L43:
	movq	%rsi, %rcx
	movq	40(%rsp), %r8
	movq	96(%rsp), %rsi
	negq	%rcx
	movq	%r12, %rdx
	movl	$.LC4, %edi
	xorl	%eax, %eax
	movq	%r9, 16(%rsp)
	movq	%r10, 8(%rsp)
	movq	%r11, 24(%rsp)
	call	printf
	movq	24(%rsp), %r11
	movq	8(%rsp), %r10
	movq	16(%rsp), %r9
	jmp	.L17
.L45:
	movq	40(%rsp), %r8
	movq	%rsi, %rcx
	movq	%r12, %rdx
	negq	%rcx
	negq	%rdx
	movq	%r14, %rsi
	movl	$.LC5, %edi
	xorl	%eax, %eax
	movq	%r9, 16(%rsp)
	movq	%r10, 8(%rsp)
	movq	%r11, 24(%rsp)
	call	printf
	movq	24(%rsp), %r11
	movq	8(%rsp), %r10
	movq	16(%rsp), %r9
	jmp	.L21
	.cfi_endproc
.LFE25:
	.size	main, .-main
	.ident	"GCC: (GNU) 4.4.7 20120313 (Red Hat 4.4.7-3)"
	.section	.note.GNU-stack,"",@progbits
