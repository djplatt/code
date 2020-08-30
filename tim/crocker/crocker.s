	.file	"crocker.cpp"
	.section	.text.unlikely,"ax",@progbits
.LCOLDB0:
	.text
.LHOTB0:
	.p2align 4,,15
	.globl	_Z8crossoutm
	.type	_Z8crossoutm, @function
_Z8crossoutm:
.LFB28:
	.cfi_startproc
	movq	%rdi, %rax
	movabsq	$-2049638230412172401, %rdx
	mulq	%rdx
	movq	%rdx, %rax
	shrq	$4, %rax
	leaq	(%rax,%rax,8), %rcx
	addq	%rcx, %rcx
	cmpq	%rcx, %rdi
	jne	.L1
	shrq	$7, %rdx
	addq	ns(%rip), %rdx
	andl	$7, %eax
	movq	masks(,%rax,8), %rax
	andb	%al, (%rdx)
.L1:
	rep ret
	.cfi_endproc
.LFE28:
	.size	_Z8crossoutm, .-_Z8crossoutm
	.section	.text.unlikely
.LCOLDE0:
	.text
.LHOTE0:
	.section	.rodata.str1.8,"aMS",@progbits,1
	.align 8
.LC1:
	.string	"Setting up table of potential Crocker Numbers."
	.align 8
.LC2:
	.string	"Crossing out sums of squares plus 0,1,2 powers of 2."
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC3:
	.string	"Found one %u at n=%lu-%lu\n"
	.section	.text.unlikely
.LCOLDB4:
	.text
.LHOTB4:
	.p2align 4,,15
	.globl	_Z10do_crockerv
	.type	_Z10do_crockerv, @function
_Z10do_crockerv:
.LFB29:
	.cfi_startproc
	pushq	%r12
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
	pushq	%rbp
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	movl	$.LC1, %edi
	pushq	%rbx
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	call	puts
	xorl	%eax, %eax
.L5:
	movq	ns(%rip), %rdx
	movb	$-1, (%rdx,%rax)
	addq	$1, %rax
	cmpq	$1048576, %rax
	jne	.L5
	movq	ns(%rip), %rax
	movl	$.LC2, %edi
	movl	$1, %ebx
	movb	$-2, (%rax)
	call	puts
	movl	$1, %r11d
	movabsq	$-2049638230412172401, %rdi
.L8:
	movq	%rbx, %r9
	movq	%rbx, %r10
	imulq	%rbx, %r9
	addq	%r11, %r9
	cmpq	$150994944, %r9
	ja	.L17
.L26:
	movq	%r9, %rax
	mulq	%rdi
	movq	%rdx, %rcx
	shrq	$4, %rcx
	leaq	(%rcx,%rcx,8), %rax
	addq	%rax, %rax
	cmpq	%rax, %r9
	jne	.L9
	shrq	$7, %rdx
	andl	$7, %ecx
	movq	%rdx, %rax
	addq	ns(%rip), %rax
	movq	masks(,%rcx,8), %rdx
	andb	%dl, (%rax)
.L9:
	movl	$1, %r8d
.L13:
	leaq	(%r8,%r9), %rsi
	cmpq	$150994944, %rsi
	ja	.L10
	.p2align 4,,10
	.p2align 3
.L30:
	movq	%rsi, %rax
	mulq	%rdi
	movq	%rdx, %rax
	shrq	$4, %rax
	leaq	(%rax,%rax,8), %rcx
	addq	%rcx, %rcx
	cmpq	%rcx, %rsi
	jne	.L11
	shrq	$7, %rdx
	addq	ns(%rip), %rdx
	andl	$7, %eax
	movq	masks(,%rax,8), %rax
	andb	%al, (%rdx)
.L11:
	addq	%r8, %r8
	leaq	(%rsi,%r8), %rbp
	movq	%r8, %rcx
	cmpq	$150994944, %rbp
	ja	.L13
	.p2align 4,,10
	.p2align 3
.L27:
	movq	%rbp, %rax
	mulq	%rdi
	movq	%rdx, %rax
	shrq	$4, %rax
	leaq	(%rax,%rax,8), %r12
	addq	%r12, %r12
	cmpq	%r12, %rbp
	jne	.L14
	shrq	$7, %rdx
	addq	ns(%rip), %rdx
	andl	$7, %eax
	movq	masks(,%rax,8), %rax
	andb	%al, (%rdx)
.L14:
	addq	%rcx, %rcx
	leaq	(%rsi,%rcx), %rbp
	cmpq	$150994944, %rbp
	jbe	.L27
	leaq	(%r8,%r9), %rsi
	cmpq	$150994944, %rsi
	jbe	.L30
.L10:
	addq	$1, %r10
	movq	%r10, %r9
	imulq	%r10, %r9
	addq	%r11, %r9
	cmpq	$150994944, %r9
	jbe	.L26
.L17:
	addq	$1, %rbx
	movq	%rbx, %r11
	imulq	%rbx, %r11
	cmpq	$150994944, %r11
	jbe	.L8
	xorl	%ebp, %ebp
	xorl	%ebx, %ebx
	jmp	.L20
.L19:
	addq	$144, %rbx
	addq	$1, %rbp
	cmpq	$150994944, %rbx
	je	.L31
.L20:
	movq	ns(%rip), %rax
	movzbl	(%rax,%rbp), %eax
	testb	%al, %al
	je	.L19
	leaq	7(%rbx), %rcx
	movzbl	%al, %esi
	movq	%rbx, %rdx
	movl	$.LC3, %edi
	xorl	%eax, %eax
	call	printf
	jmp	.L19
.L31:
	popq	%rbx
	.cfi_def_cfa_offset 24
	popq	%rbp
	.cfi_def_cfa_offset 16
	popq	%r12
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE29:
	.size	_Z10do_crockerv, .-_Z10do_crockerv
	.section	.text.unlikely
.LCOLDE4:
	.text
.LHOTE4:
	.section	.rodata.str1.8
	.align 8
.LC5:
	.string	"Checking for Crocker numbers <=%lu\n"
	.section	.rodata.str1.1
.LC6:
	.string	"Finished."
	.section	.text.unlikely
.LCOLDB7:
	.section	.text.startup,"ax",@progbits
.LHOTB7:
	.p2align 4,,15
	.globl	main
	.type	main, @function
main:
.LFB30:
	.cfi_startproc
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	movl	$150994944, %esi
	movl	$.LC5, %edi
	xorl	%eax, %eax
	call	printf
	movl	$1048576, %edi
	call	malloc
	movq	%rax, ns(%rip)
	call	_Z10do_crockerv
	movl	$.LC6, %edi
	call	puts
	xorl	%eax, %eax
	addq	$8, %rsp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE30:
	.size	main, .-main
	.section	.text.unlikely
.LCOLDE7:
	.section	.text.startup
.LHOTE7:
	.globl	masks
	.data
	.align 32
	.type	masks, @object
	.size	masks, 64
masks:
	.quad	254
	.quad	253
	.quad	251
	.quad	247
	.quad	239
	.quad	223
	.quad	191
	.quad	127
	.globl	ns
	.bss
	.align 8
	.type	ns, @object
	.size	ns, 8
ns:
	.zero	8
	.ident	"GCC: (GNU) 5.1.0"
	.section	.note.GNU-stack,"",@progbits
