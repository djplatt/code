/ Speed test routine
/
/ Author: David Harper at obliquity.com
/
/ C prototype:
/
/ void vectorOps(double *a, double *b, double *c, int nsize,
/		 int niters, int mode)
/
/ mode = 0 --> NOP
/	 1 --> ADD
/	 2 --> MULTIPLY
/	 3 --> DIVIDE
/	 4 --> COSINE
/	 5 --> SQRT
/	 6 --> ATAN2
/	 7 --> Y LOG2(X)
/	 8 --> SINCOS
	






	.text

	.align	4
.globl	vectorOps
	.type	vectorOps,@function
vectorOps:
/	Set up the stack frame
	pushl	%ebp
	movl	%esp,%ebp
/	Save registers that we want to restore later
	pushl	%esi
	pushl	%ebx

/	Test niters > 0
	movl	24(%ebp),%ecx
	xorl	%eax,%eax
	cmpl	%ecx,%eax
	jge	.Lexit

/	Test nsize > 0
	movl	20(%ebp),%ecx
	xorl	%eax,%eax
	cmpl	%ecx,%eax
	jge	.Lexit

/	On the value of mode, skip to the relevant section
	movl	28(%ebp),%eax

	
        cmpl    $0,%eax
        je      .Lnop
	
        cmpl    $1,%eax
        je      .Ladd
	
        cmpl    $2,%eax
        je      .Lmultiply
	
        cmpl    $3,%eax
        je      .Ldivide
	
        cmpl    $4,%eax
        je      .Lcosine
	
        cmpl    $5,%eax
        je      .Lsqrt
	
        cmpl    $6,%eax
        je      .Latan
	
        cmpl    $7,%eax
        je      .Lylogx
	
        cmpl    $8,%eax
        je      .Lsincos
	
        cmpl    $9,%eax
        je      .Lexp

/	mode lies outside the range, so exit
	xorl	%eax,%eax
	jmp	.Lexit

	
	.align	4
.Lnop:
	movl	24(%ebp),%ecx

.LnopOuterLoop:
	movl	%ecx,%esi

	movl	20(%ebp),%ecx
	movl	8(%ebp),%eax
	movl	12(%ebp),%ebx
	movl	16(%ebp),%edx

.LnopInnerLoop:
	
	addl	$8,%eax
	addl	$8,%ebx
	addl	$8,%edx

	loop	.LnopInnerLoop

	movl	%esi,%ecx
	loop	.LnopOuterLoop

	movl	24(%ebp),%eax

	jmp	.Lexit

	
	.align	4
.Ladd:
	movl	24(%ebp),%ecx

.LaddOuterLoop:
	movl	%ecx,%esi

	movl	20(%ebp),%ecx
	movl	8(%ebp),%eax
	movl	12(%ebp),%ebx
	movl	16(%ebp),%edx

.LaddInnerLoop:
	fldl	(%eax)
	faddl	(%ebx)
	fstpl	(%edx)
	
	addl	$8,%eax
	addl	$8,%ebx
	addl	$8,%edx

	loop	.LaddInnerLoop

	movl	%esi,%ecx
	loop	.LaddOuterLoop

	movl	24(%ebp),%eax

	jmp	.Lexit

	
	.align	4
.Lmultiply:
	movl	24(%ebp),%ecx

.LmultiplyOuterLoop:
	movl	%ecx,%esi

	movl	20(%ebp),%ecx
	movl	8(%ebp),%eax
	movl	12(%ebp),%ebx
	movl	16(%ebp),%edx

.LmultiplyInnerLoop:
	fldl	(%eax)
	fmull	(%ebx)
	fstpl	(%edx)
	
	addl	$8,%eax
	addl	$8,%ebx
	addl	$8,%edx

	loop	.LmultiplyInnerLoop

	movl	%esi,%ecx
	loop	.LmultiplyOuterLoop

	movl	24(%ebp),%eax

	jmp	.Lexit

	
	.align	4
.Ldivide:
	movl	24(%ebp),%ecx

.LdivideOuterLoop:
	movl	%ecx,%esi

	movl	20(%ebp),%ecx
	movl	8(%ebp),%eax
	movl	12(%ebp),%ebx
	movl	16(%ebp),%edx

.LdivideInnerLoop:
	fldl	(%eax)
	fdivrl	(%ebx)
	fstpl	(%edx)
	
	addl	$8,%eax
	addl	$8,%ebx
	addl	$8,%edx

	loop	.LdivideInnerLoop

	movl	%esi,%ecx
	loop	.LdivideOuterLoop

	movl	24(%ebp),%eax

	jmp	.Lexit

	
	.align	4
.Lcosine:
	movl	24(%ebp),%ecx

.LcosineOuterLoop:
	movl	%ecx,%esi

	movl	20(%ebp),%ecx
	movl	8(%ebp),%eax
	movl	12(%ebp),%ebx
	movl	16(%ebp),%edx

.LcosineInnerLoop:
	fldl	(%eax)
	fcos
	fstpl	(%edx)
	
	addl	$8,%eax
	addl	$8,%ebx
	addl	$8,%edx

	loop	.LcosineInnerLoop

	movl	%esi,%ecx
	loop	.LcosineOuterLoop

	movl	24(%ebp),%eax

	jmp	.Lexit

	
	.align	4
.Lsqrt:
	movl	24(%ebp),%ecx

.LsqrtOuterLoop:
	movl	%ecx,%esi

	movl	20(%ebp),%ecx
	movl	8(%ebp),%eax
	movl	12(%ebp),%ebx
	movl	16(%ebp),%edx

.LsqrtInnerLoop:
	fldl	(%eax)
	fsqrt
	fstpl	(%edx)
	
	addl	$8,%eax
	addl	$8,%ebx
	addl	$8,%edx

	loop	.LsqrtInnerLoop

	movl	%esi,%ecx
	loop	.LsqrtOuterLoop

	movl	24(%ebp),%eax

	jmp	.Lexit

	
	.align	4
.Latan:
	movl	24(%ebp),%ecx

.LatanOuterLoop:
	movl	%ecx,%esi

	movl	20(%ebp),%ecx
	movl	8(%ebp),%eax
	movl	12(%ebp),%ebx
	movl	16(%ebp),%edx

.LatanInnerLoop:
	fldl	(%eax)
	fldl	(%ebx)
	fpatan
	fstpl	(%edx)
	
	addl	$8,%eax
	addl	$8,%ebx
	addl	$8,%edx

	loop	.LatanInnerLoop

	movl	%esi,%ecx
	loop	.LatanOuterLoop

	movl	24(%ebp),%eax

	jmp	.Lexit

	
	.align	4
.Lylogx:
	movl	24(%ebp),%ecx

.LylogxOuterLoop:
	movl	%ecx,%esi

	movl	20(%ebp),%ecx
	movl	8(%ebp),%eax
	movl	12(%ebp),%ebx
	movl	16(%ebp),%edx

.LylogxInnerLoop:
	fldl	(%eax)
	fldl	(%ebx)
	fyl2x
	fstpl	(%edx)
	
	addl	$8,%eax
	addl	$8,%ebx
	addl	$8,%edx

	loop	.LylogxInnerLoop

	movl	%esi,%ecx
	loop	.LylogxOuterLoop

	movl	24(%ebp),%eax

	jmp	.Lexit

	
	.align	4
.Lsincos:
	movl	24(%ebp),%ecx

.LsincosOuterLoop:
	movl	%ecx,%esi

	movl	20(%ebp),%ecx
	movl	8(%ebp),%eax
	movl	12(%ebp),%ebx
	movl	16(%ebp),%edx

.LsincosInnerLoop:
	fldl	(%eax)
	fsincos
	fstpl	(%edx)
	fstpl	(%edx)
	
	addl	$8,%eax
	addl	$8,%ebx
	addl	$8,%edx

	loop	.LsincosInnerLoop

	movl	%esi,%ecx
	loop	.LsincosOuterLoop

	movl	24(%ebp),%eax

	jmp	.Lexit

	
	.align	4
.Lexp:
	movl	24(%ebp),%ecx

.LexpOuterLoop:
	movl	%ecx,%esi

	movl	20(%ebp),%ecx
	movl	8(%ebp),%eax
	movl	12(%ebp),%ebx
	movl	16(%ebp),%edx

.LexpInnerLoop:
	fldl	(%eax)
	f2xm1
	fstpl	(%edx)
	
	addl	$8,%eax
	addl	$8,%ebx
	addl	$8,%edx

	loop	.LexpInnerLoop

	movl	%esi,%ecx
	loop	.LexpOuterLoop

	movl	24(%ebp),%eax

	jmp	.Lexit

.Lexit:
/	Restore the saved registers
	popl	%ebx
	popl	%esi
/	Clear up the stack frame and return
	leave
	ret

.Lend:
	.size	vectorOps,.Lend-vectorOps

	.ident	"David Harper at www.obliquity.com"
