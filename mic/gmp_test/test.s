	.file	"test.c"
	.section	.rodata
.LC0:
	.string	"\nres=%d\n"
	.text
.globl main
	.type	main, @function
main:
.LFB19:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$48, %rsp
	movl	%edi, -36(%rbp)
	movq	%rsi, -48(%rbp)
	leaq	-32(%rbp), %rax
	movq	%rax, %rdi
	call	__gmpz_init
	leaq	-32(%rbp), %rax
	movl	$1, %esi
	movq	%rax, %rdi
	call	__gmpz_set_ui
	leaq	-32(%rbp), %rcx
	leaq	-32(%rbp), %rax
	movl	$100, %edx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	__gmpz_mul_2exp
	movq	stdout(%rip), %rax
	leaq	-32(%rbp), %rdx
	movl	$10, %esi
	movq	%rax, %rdi
	call	__gmpz_out_str
	movl	%eax, -4(%rbp)
	movl	$.LC0, %eax
	movl	-4(%rbp), %edx
	movl	%edx, %esi
	movq	%rax, %rdi
	movl	$0, %eax
	call	printf
	movl	$0, %eax
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE19:
	.size	main, .-main
	.ident	"GCC: (GNU) 4.4.7 20120313 (Red Hat 4.4.7-16)"
	.section	.note.GNU-stack,"",@progbits
