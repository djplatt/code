	.file	"cos_acc.c"
	.comm	old_cw,2,2
	.comm	new_cw,2,2
	.section	.rodata
.LC0:
	.string	"old cw=%x\n"
.LC1:
	.string	"rndd cw=%X\n"
	.text
.globl rndd
	.type	rndd, @function
rndd:
.LFB19:
	.cfi_startproc
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
#APP
# 26 "cos_acc.c" 1
	fnstcw old_cw(%rip)
# 0 "" 2
#NO_APP
	movzwl	old_cw(%rip), %eax
	movswl	%ax,%edx
	movl	$.LC0, %eax
	movl	%edx, %esi
	movq	%rax, %rdi
	movl	$0, %eax
	call	printf
	movzwl	old_cw(%rip), %eax
	andw	$255, %ax
	movw	%ax, new_cw(%rip)
	movzwl	new_cw(%rip), %eax
	orb	$6, %ah
	movw	%ax, new_cw(%rip)
	movzwl	new_cw(%rip), %eax
	movswl	%ax,%edx
	movl	$.LC1, %eax
	movl	%edx, %esi
	movq	%rax, %rdi
	movl	$0, %eax
	call	printf
#APP
# 31 "cos_acc.c" 1
	fldcw new_cw(%rip)
# 0 "" 2
#NO_APP
	addq	$8, %rsp
	ret
	.cfi_endproc
.LFE19:
	.size	rndd, .-rndd
	.section	.rodata
.LC2:
	.string	"rndu cw=%X\n"
	.text
.globl rndu
	.type	rndu, @function
rndu:
.LFB20:
	.cfi_startproc
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
#APP
# 36 "cos_acc.c" 1
	fnstcw old_cw(%rip)
# 0 "" 2
#NO_APP
	movzwl	old_cw(%rip), %eax
	movswl	%ax,%edx
	movl	$.LC0, %eax
	movl	%edx, %esi
	movq	%rax, %rdi
	movl	$0, %eax
	call	printf
	movzwl	old_cw(%rip), %eax
	andw	$255, %ax
	movw	%ax, new_cw(%rip)
	movzwl	new_cw(%rip), %eax
	orb	$10, %ah
	movw	%ax, new_cw(%rip)
	movzwl	new_cw(%rip), %eax
	movswl	%ax,%edx
	movl	$.LC2, %eax
	movl	%edx, %esi
	movq	%rax, %rdi
	movl	$0, %eax
	call	printf
#APP
# 41 "cos_acc.c" 1
	fldcw new_cw(%rip)
# 0 "" 2
#NO_APP
	addq	$8, %rsp
	ret
	.cfi_endproc
.LFE20:
	.size	rndu, .-rndu
	.section	.rodata
.LC3:
	.string	"rndn cw=%X\n"
	.text
.globl rndn
	.type	rndn, @function
rndn:
.LFB21:
	.cfi_startproc
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
#APP
# 46 "cos_acc.c" 1
	fnstcw old_cw(%rip)
# 0 "" 2
#NO_APP
	movzwl	old_cw(%rip), %eax
	movswl	%ax,%edx
	movl	$.LC0, %eax
	movl	%edx, %esi
	movq	%rax, %rdi
	movl	$0, %eax
	call	printf
	movzwl	old_cw(%rip), %eax
	andw	$255, %ax
	movw	%ax, new_cw(%rip)
	movzwl	new_cw(%rip), %eax
	orb	$2, %ah
	movw	%ax, new_cw(%rip)
	movzwl	new_cw(%rip), %eax
	movswl	%ax,%edx
	movl	$.LC3, %eax
	movl	%edx, %esi
	movq	%rax, %rdi
	movl	$0, %eax
	call	printf
#APP
# 51 "cos_acc.c" 1
	fldcw new_cw(%rip)
# 0 "" 2
#NO_APP
	addq	$8, %rsp
	ret
	.cfi_endproc
.LFE21:
	.size	rndn, .-rndn
	.section	.rodata
.LC4:
	.string	"%A\n"
	.text
.globl print_double
	.type	print_double, @function
print_double:
.LFB22:
	.cfi_startproc
	subq	$40, %rsp
	.cfi_def_cfa_offset 48
	movsd	%xmm0, 8(%rsp)
	movsd	8(%rsp), %xmm0
	movl	$.LC4, %eax
	movq	%rax, %rdi
	movl	$1, %eax
	call	printf
	addq	$40, %rsp
	ret
	.cfi_endproc
.LFE22:
	.size	print_double, .-print_double
	.section	.rodata
.LC5:
	.string	"Huge difference. Exiting."
	.text
.globl ulp_diff
	.type	ulp_diff, @function
ulp_diff:
.LFB23:
	.cfi_startproc
	subq	$40, %rsp
	.cfi_def_cfa_offset 48
	movsd	%xmm0, 8(%rsp)
	movsd	%xmm1, (%rsp)
	leaq	8(%rsp), %rax
	movq	%rax, 24(%rsp)
	movq	%rsp, %rax
	movq	%rax, 16(%rsp)
	movq	24(%rsp), %rax
	addq	$4, %rax
	movl	(%rax), %edx
	movq	16(%rsp), %rax
	addq	$4, %rax
	movl	(%rax), %eax
	cmpl	%eax, %edx
	je	.L10
	movl	$.LC5, %edi
	call	puts
	movl	$0, %edi
	call	exit
.L10:
	movq	24(%rsp), %rax
	movl	(%rax), %edx
	movq	16(%rsp), %rax
	movl	(%rax), %eax
	movl	%edx, %ecx
	subl	%eax, %ecx
	movl	%ecx, %eax
	addq	$40, %rsp
	ret
	.cfi_endproc
.LFE23:
	.size	ulp_diff, .-ulp_diff
.globl mycos
	.type	mycos, @function
mycos:
.LFB24:
	.cfi_startproc
	movsd	%xmm0, -24(%rsp)
#APP
# 86 "cos_acc.c" 1
	fldl -24(%rsp)
	fcos
	fstpl -8(%rsp)
# 0 "" 2
#NO_APP
	movq	-8(%rsp), %rax
	movq	%rax, -32(%rsp)
	movsd	-32(%rsp), %xmm0
	ret
	.cfi_endproc
.LFE24:
	.size	mycos, .-mycos
.globl rnd_test
	.type	rnd_test, @function
rnd_test:
.LFB25:
	.cfi_startproc
	movsd	%xmm0, -8(%rsp)
	movsd	%xmm1, -16(%rsp)
	movsd	-8(%rsp), %xmm0
	divsd	-16(%rsp), %xmm0
	ret
	.cfi_endproc
.LFE25:
	.size	rnd_test, .-rnd_test
	.section	.rodata
.LC7:
	.string	"rndn returned %50.48e "
.LC8:
	.string	"rndd returned %50.48e "
.LC9:
	.string	"rndu returned %50.48e "
	.text
.globl main
	.type	main, @function
main:
.LFB26:
	.cfi_startproc
	subq	$328, %rsp
	.cfi_def_cfa_offset 336
	movl	%edi, 44(%rsp)
	movq	%rsi, 32(%rsp)
	movl	$100, %edi
	call	mpfr_set_default_prec
	leaq	208(%rsp), %rax
	movq	%rax, %rdi
	call	mpfr_init
	leaq	176(%rsp), %rax
	movq	%rax, %rdi
	call	mpfr_init
	leaq	112(%rsp), %rax
	movq	%rax, %rdi
	call	mpfi_init
	leaq	48(%rsp), %rax
	movq	%rax, %rdi
	call	mpfi_init
	movl	$0, 308(%rsp)
	movq	32(%rsp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, %rdi
	call	atof
	movsd	%xmm0, 8(%rsp)
	movq	32(%rsp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, %rdi
	call	atof
	movsd	8(%rsp), %xmm1
	call	rnd_test
	movsd	%xmm0, 248(%rsp)
	movl	$0, %eax
	call	rndu
	movq	32(%rsp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, %rdi
	call	atof
	movsd	%xmm0, 16(%rsp)
	movq	32(%rsp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, %rdi
	call	atof
	movsd	16(%rsp), %xmm1
	call	rnd_test
	movsd	%xmm0, 256(%rsp)
	movl	$0, %eax
	call	rndd
	movq	32(%rsp), %rax
	addq	$16, %rax
	movq	(%rax), %rax
	movq	%rax, %rdi
	call	atof
	movsd	%xmm0, 24(%rsp)
	movq	32(%rsp), %rax
	addq	$8, %rax
	movq	(%rax), %rax
	movq	%rax, %rdi
	call	atof
	movsd	24(%rsp), %xmm1
	call	rnd_test
	movsd	%xmm0, 264(%rsp)
	movl	$0, %eax
	call	rndn
	movl	$.LC7, %eax
	movsd	264(%rsp), %xmm0
	movq	%rax, %rdi
	movl	$1, %eax
	call	printf
	movsd	248(%rsp), %xmm0
	call	print_double
	movl	$.LC8, %eax
	movsd	264(%rsp), %xmm0
	movq	%rax, %rdi
	movl	$1, %eax
	call	printf
	movsd	264(%rsp), %xmm0
	call	print_double
	movl	$.LC9, %eax
	movsd	256(%rsp), %xmm0
	movq	%rax, %rdi
	movl	$1, %eax
	call	printf
	movsd	256(%rsp), %xmm0
	call	print_double
	movl	$0, %edi
	call	exit
	.cfi_endproc
.LFE26:
	.size	main, .-main
	.section	.rodata
	.type	__PRETTY_FUNCTION__.6968, @object
	.size	__PRETTY_FUNCTION__.6968, 5
__PRETTY_FUNCTION__.6968:
	.string	"main"
	.ident	"GCC: (Ubuntu 4.4.3-4ubuntu5) 4.4.3"
	.section	.note.GNU-stack,"",@progbits
