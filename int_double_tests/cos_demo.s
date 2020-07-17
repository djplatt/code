	.file	"cos_demo.c"
	.text
	.p2align 4,,15
.globl rndd
	.type	rndd, @function
rndd:
.LFB27:
#APP
# 23 "cos_demo.c" 1
	fnstcw old_cw(%rip)
# 0 "" 2
#NO_APP
	movzbl	old_cw(%rip), %eax
	orb	$7, %ah
	movw	%ax, new_cw(%rip)
#APP
# 28 "cos_demo.c" 1
	fldcw new_cw(%rip)
# 0 "" 2
#NO_APP
	ret
.LFE27:
	.size	rndd, .-rndd
	.p2align 4,,15
.globl rndu
	.type	rndu, @function
rndu:
.LFB28:
#APP
# 33 "cos_demo.c" 1
	fnstcw old_cw(%rip)
# 0 "" 2
#NO_APP
	movzbl	old_cw(%rip), %eax
	orb	$11, %ah
	movw	%ax, new_cw(%rip)
#APP
# 38 "cos_demo.c" 1
	fldcw new_cw(%rip)
# 0 "" 2
#NO_APP
	ret
.LFE28:
	.size	rndu, .-rndu
	.p2align 4,,15
.globl rndn
	.type	rndn, @function
rndn:
.LFB29:
#APP
# 43 "cos_demo.c" 1
	fnstcw old_cw(%rip)
# 0 "" 2
#NO_APP
	movzbl	old_cw(%rip), %eax
	movw	%ax, new_cw(%rip)
#APP
# 48 "cos_demo.c" 1
	fldcw new_cw(%rip)
# 0 "" 2
#NO_APP
	ret
.LFE29:
	.size	rndn, .-rndn
	.p2align 4,,15
.globl mycos
	.type	mycos, @function
mycos:
.LFB32:
	movsd	%xmm0, -24(%rsp)
#APP
# 81 "cos_demo.c" 1
	fldl -24(%rsp)
	fcos
	fstpl -8(%rsp)
# 0 "" 2
#NO_APP
	movsd	-8(%rsp), %xmm0
	ret
.LFE32:
	.size	mycos, .-mycos
	.p2align 4,,15
.globl mysin
	.type	mysin, @function
mysin:
.LFB33:
	movsd	%xmm0, -24(%rsp)
#APP
# 94 "cos_demo.c" 1
	fldl -24(%rsp)
	fsin
	fstpl -8(%rsp)
# 0 "" 2
#NO_APP
	movsd	-8(%rsp), %xmm0
	ret
.LFE33:
	.size	mysin, .-mysin
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC1:
	.string	"%X %X\n"
	.text
	.p2align 4,,15
.globl print_double
	.type	print_double, @function
print_double:
.LFB30:
	subq	$8, %rsp
.LCFI0:
	movl	$.LC1, %edi
	xorl	%eax, %eax
	movl	4(%rsp), %esi
	movl	(%rsp), %edx
	movsd	%xmm0, (%rsp)
	call	printf
	addq	$8, %rsp
	ret
.LFE30:
	.size	print_double, .-print_double
	.section	.rodata.str1.1
.LC2:
	.string	"testing argument x=%30.28e = "
.LC3:
	.string	"cos  "
	.section	.rodata.str1.8,"aMS",@progbits,1
	.align 8
.LC4:
	.string	"Gcc library   %s produced %30.28e = "
	.section	.rodata.str1.1
.LC5:
	.string	"cosl "
	.section	.rodata.str1.8
	.align 8
.LC6:
	.string	"Gcc L library %s produced %30.28e = "
	.section	.rodata.str1.1
.LC7:
	.string	"fcos "
	.section	.rodata.str1.8
	.align 8
.LC8:
	.string	"rounding down %s produced %30.28e = "
	.align 8
.LC9:
	.string	"rounding up   %s produced %30.28e = "
	.align 8
.LC10:
	.string	"rounding near %s produced %30.28e = "
	.text
	.p2align 4,,15
.globl test_val
	.type	test_val, @function
test_val:
.LFB34:
	pushq	%r12
.LCFI1:
	movl	$.LC2, %edi
	movl	$1, %eax
	pushq	%rbx
.LCFI2:
	subq	$104, %rsp
.LCFI3:
	movsd	%xmm0, 24(%rsp)
	call	printf
	movl	88(%rsp), %edx
	movl	$.LC1, %edi
	xorl	%eax, %eax
	movl	92(%rsp), %esi
	fldl	24(%rsp)
	fstpl	88(%rsp)
	call	printf
	movsd	24(%rsp), %xmm0
	call	cos
	movl	$.LC3, %esi
	movl	$.LC4, %edi
	movl	$1, %eax
	movsd	%xmm0, 56(%rsp)
	call	printf
	fldl	56(%rsp)
	movl	$.LC1, %edi
	xorl	%eax, %eax
	movl	88(%rsp), %edx
	movl	92(%rsp), %esi
	fstpl	88(%rsp)
	call	printf
	fldl	24(%rsp)
	fstpt	(%rsp)
	call	cosl
	fstpl	64(%rsp)
	movl	$.LC5, %esi
	movl	$.LC6, %edi
	movl	$1, %eax
	movsd	64(%rsp), %xmm0
	call	printf
	movl	92(%rsp), %esi
	movl	88(%rsp), %edx
	movl	$.LC1, %edi
	movsd	64(%rsp), %xmm0
	xorl	%eax, %eax
	movsd	%xmm0, 88(%rsp)
	call	printf
#APP
# 23 "cos_demo.c" 1
	fnstcw old_cw(%rip)
# 0 "" 2
#NO_APP
	movzwl	old_cw(%rip), %ebx
	movzbl	%bl, %r12d
	movl	%r12d, %eax
	orb	$7, %ah
	movw	%ax, new_cw(%rip)
#APP
# 28 "cos_demo.c" 1
	fldcw new_cw(%rip)
# 0 "" 2
#NO_APP
	fldl	24(%rsp)
	movl	$.LC7, %esi
	movl	$.LC8, %edi
	movl	$1, %eax
	fstpl	88(%rsp)
#APP
# 81 "cos_demo.c" 1
	fldl 88(%rsp)
	fcos
	fstpl 80(%rsp)
# 0 "" 2
#NO_APP
	movsd	80(%rsp), %xmm0
	movsd	%xmm0, 32(%rsp)
	call	printf
	movl	84(%rsp), %esi
	movl	$.LC1, %edi
	xorl	%eax, %eax
	movl	80(%rsp), %edx
	movsd	32(%rsp), %xmm0
	movsd	%xmm0, 80(%rsp)
	call	printf
	movl	%r12d, %eax
	movw	%bx, old_cw(%rip)
	orb	$11, %ah
	movw	%ax, new_cw(%rip)
#APP
# 38 "cos_demo.c" 1
	fldcw new_cw(%rip)
# 0 "" 2
#NO_APP
	fldl	24(%rsp)
	movl	$.LC7, %esi
	movl	$.LC9, %edi
	movl	$1, %eax
	fstpl	80(%rsp)
#APP
# 81 "cos_demo.c" 1
	fldl 80(%rsp)
	fcos
	fstpl 88(%rsp)
# 0 "" 2
#NO_APP
	movsd	88(%rsp), %xmm0
	movsd	%xmm0, 40(%rsp)
	call	printf
	movl	92(%rsp), %esi
	movl	$.LC1, %edi
	xorl	%eax, %eax
	movl	88(%rsp), %edx
	movsd	40(%rsp), %xmm0
	movsd	%xmm0, 88(%rsp)
	call	printf
	movw	%bx, old_cw(%rip)
	movw	%r12w, new_cw(%rip)
#APP
# 48 "cos_demo.c" 1
	fldcw new_cw(%rip)
# 0 "" 2
#NO_APP
	fldl	24(%rsp)
	movl	$.LC7, %esi
	movl	$.LC10, %edi
	movl	$1, %eax
	fstpl	88(%rsp)
#APP
# 81 "cos_demo.c" 1
	fldl 88(%rsp)
	fcos
	fstpl 80(%rsp)
# 0 "" 2
#NO_APP
	movsd	80(%rsp), %xmm0
	movsd	%xmm0, 48(%rsp)
	call	printf
	movl	84(%rsp), %esi
	movl	$.LC1, %edi
	xorl	%eax, %eax
	movl	80(%rsp), %edx
	movsd	48(%rsp), %xmm0
	movsd	%xmm0, 80(%rsp)
	call	printf
	addq	$104, %rsp
	xorl	%eax, %eax
	popq	%rbx
	popq	%r12
	ret
.LFE34:
	.size	test_val, .-test_val
	.p2align 4,,15
.globl main
	.type	main, @function
main:
.LFB35:
	subq	$8, %rsp
.LCFI4:
	xorl	%eax, %eax
	fldl	.LC11(%rip)
	flds	.LC12(%rip)
	.p2align 4,,10
	.p2align 3
.L16:
	addl	$1, %eax
	fmul	%st, %st(1)
	cmpl	$49, %eax
	jne	.L16
	fstp	%st(0)
	fstpl	(%rsp)
	movsd	(%rsp), %xmm0
	call	test_val
	xorl	%eax, %eax
	addq	$8, %rsp
	ret
.LFE35:
	.size	main, .-main
	.section	.rodata.str1.1
.LC13:
	.string	"Huge difference. Exiting."
	.text
	.p2align 4,,15
.globl ulp_diff
	.type	ulp_diff, @function
ulp_diff:
.LFB31:
	subq	$24, %rsp
.LCFI5:
	movl	20(%rsp), %eax
	cmpl	12(%rsp), %eax
	movsd	%xmm0, 16(%rsp)
	movsd	%xmm1, 8(%rsp)
	jne	.L23
	movl	16(%rsp), %eax
	subl	8(%rsp), %eax
	addq	$24, %rsp
	ret
.L23:
	movl	$.LC13, %edi
	call	puts
	xorl	%edi, %edi
	call	exit
.LFE31:
	.size	ulp_diff, .-ulp_diff
	.comm	old_cw,2,2
	.comm	new_cw,2,2
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC11:
	.long	1413326760
	.long	1124671995
	.section	.rodata.cst4,"aM",@progbits,4
	.align 4
.LC12:
	.long	1056964608
	.section	.eh_frame,"a",@progbits
.Lframe1:
	.long	.LECIE1-.LSCIE1
.LSCIE1:
	.long	0x0
	.byte	0x1
	.string	"zR"
	.uleb128 0x1
	.sleb128 -8
	.byte	0x10
	.uleb128 0x1
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
	.long	.LFB27
	.long	.LFE27-.LFB27
	.uleb128 0x0
	.align 8
.LEFDE1:
.LSFDE3:
	.long	.LEFDE3-.LASFDE3
.LASFDE3:
	.long	.LASFDE3-.Lframe1
	.long	.LFB28
	.long	.LFE28-.LFB28
	.uleb128 0x0
	.align 8
.LEFDE3:
.LSFDE5:
	.long	.LEFDE5-.LASFDE5
.LASFDE5:
	.long	.LASFDE5-.Lframe1
	.long	.LFB29
	.long	.LFE29-.LFB29
	.uleb128 0x0
	.align 8
.LEFDE5:
.LSFDE7:
	.long	.LEFDE7-.LASFDE7
.LASFDE7:
	.long	.LASFDE7-.Lframe1
	.long	.LFB32
	.long	.LFE32-.LFB32
	.uleb128 0x0
	.align 8
.LEFDE7:
.LSFDE9:
	.long	.LEFDE9-.LASFDE9
.LASFDE9:
	.long	.LASFDE9-.Lframe1
	.long	.LFB33
	.long	.LFE33-.LFB33
	.uleb128 0x0
	.align 8
.LEFDE9:
.LSFDE11:
	.long	.LEFDE11-.LASFDE11
.LASFDE11:
	.long	.LASFDE11-.Lframe1
	.long	.LFB30
	.long	.LFE30-.LFB30
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI0-.LFB30
	.byte	0xe
	.uleb128 0x10
	.align 8
.LEFDE11:
.LSFDE13:
	.long	.LEFDE13-.LASFDE13
.LASFDE13:
	.long	.LASFDE13-.Lframe1
	.long	.LFB34
	.long	.LFE34-.LFB34
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI1-.LFB34
	.byte	0xe
	.uleb128 0x10
	.byte	0x4
	.long	.LCFI2-.LCFI1
	.byte	0xe
	.uleb128 0x18
	.byte	0x4
	.long	.LCFI3-.LCFI2
	.byte	0xe
	.uleb128 0x80
	.byte	0x83
	.uleb128 0x3
	.byte	0x8c
	.uleb128 0x2
	.align 8
.LEFDE13:
.LSFDE15:
	.long	.LEFDE15-.LASFDE15
.LASFDE15:
	.long	.LASFDE15-.Lframe1
	.long	.LFB35
	.long	.LFE35-.LFB35
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI4-.LFB35
	.byte	0xe
	.uleb128 0x10
	.align 8
.LEFDE15:
.LSFDE17:
	.long	.LEFDE17-.LASFDE17
.LASFDE17:
	.long	.LASFDE17-.Lframe1
	.long	.LFB31
	.long	.LFE31-.LFB31
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI5-.LFB31
	.byte	0xe
	.uleb128 0x20
	.align 8
.LEFDE17:
	.ident	"GCC: (GNU) 4.3.3"
	.section	.note.GNU-stack,"",@progbits
