	.file	"check_proth.cpp"
	.text
	.type	_GLOBAL__I__Z7proth_pmmm, @function
_GLOBAL__I__Z7proth_pmmm:
.LFB1535:
	subq	$8, %rsp
.LCFI0:
	movl	$_ZStL8__ioinit, %edi
	call	_ZNSt8ios_base4InitC1Ev
	movl	$__dso_handle, %edx
	movl	$_ZStL8__ioinit, %esi
	movl	$_ZNSt8ios_base4InitD1Ev, %edi
	call	__cxa_atexit
	movl	$_ZN3clnL34cl_random_def_init_helper_instanceE, %edi
	call	_ZN3cln25cl_random_def_init_helperC1Ev
	movl	$__dso_handle, %edx
	movl	$_ZN3clnL34cl_random_def_init_helper_instanceE, %esi
	movl	$_ZN3cln25cl_random_def_init_helperD1Ev, %edi
	call	__cxa_atexit
	movl	$_ZN3clnL27cl_I_classes_dummy_instanceE, %eax
	cmpq	$_ZN3cln15cl_class_fixnumE, %rax
	jne	.L2
	call	abort
.L2:
	movl	$_ZN3clnL36cl_prin_globals_init_helper_instanceE, %edi
	call	_ZN3cln27cl_prin_globals_init_helperC1Ev
	movl	$__dso_handle, %edx
	movl	$_ZN3clnL36cl_prin_globals_init_helper_instanceE, %esi
	movl	$_ZN3cln27cl_prin_globals_init_helperD1Ev, %edi
	addq	$8, %rsp
	jmp	__cxa_atexit
.LFE1535:
	.size	_GLOBAL__I__Z7proth_pmmm, .-_GLOBAL__I__Z7proth_pmmm
	.section	.ctors,"aw",@progbits
	.align 8
	.quad	_GLOBAL__I__Z7proth_pmmm
	.section	.text._ZN3cln4cl_ID1Ev,"axG",@progbits,_ZN3cln4cl_ID1Ev,comdat
	.align 2
	.weak	_ZN3cln4cl_ID1Ev
	.type	_ZN3cln4cl_ID1Ev, @function
_ZN3cln4cl_ID1Ev:
.LFB1339:
	movq	(%rdi), %rdi
	testb	$7, %dil
	jne	.L7
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L8
.L7:
	ret
.L8:
	jmp	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.LFE1339:
	.size	_ZN3cln4cl_ID1Ev, .-_ZN3cln4cl_ID1Ev
.globl _Unwind_Resume
	.section	.text._ZN3cln10cl_I_div_tD1Ev,"axG",@progbits,_ZN3cln10cl_I_div_tD1Ev,comdat
	.align 2
	.weak	_ZN3cln10cl_I_div_tD1Ev
	.type	_ZN3cln10cl_I_div_tD1Ev, @function
_ZN3cln10cl_I_div_tD1Ev:
.LFB1454:
	pushq	%rbp
.LCFI1:
	pushq	%rbx
.LCFI2:
	subq	$8, %rsp
.LCFI3:
	movq	%rdi, %rbp
	movq	8(%rdi), %rdi
	testb	$7, %dil
	jne	.L10
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L15
.L10:
	movq	(%rbp), %rdi
	testb	$7, %dil
	jne	.L13
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L16
.L13:
	addq	$8, %rsp
	popq	%rbx
	popq	%rbp
	ret
.L16:
	addq	$8, %rsp
	popq	%rbx
	popq	%rbp
.LEHB0:
	jmp	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.LEHE0:
.L15:
.LEHB1:
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.LEHE1:
	jmp	.L10
.L14:
	movq	%rax, %rbx
.L11:
	movq	%rbp, %rdi
	call	_ZN3cln4cl_ID1Ev
	movq	%rbx, %rdi
.LEHB2:
	call	_Unwind_Resume
.LEHE2:
.LFE1454:
	.size	_ZN3cln10cl_I_div_tD1Ev, .-_ZN3cln10cl_I_div_tD1Ev
.globl __gxx_personality_v0
	.section	.gcc_except_table,"a",@progbits
.LLSDA1454:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE1454-.LLSDACSB1454
.LLSDACSB1454:
	.uleb128 .LEHB0-.LFB1454
	.uleb128 .LEHE0-.LEHB0
	.uleb128 0x0
	.uleb128 0x0
	.uleb128 .LEHB1-.LFB1454
	.uleb128 .LEHE1-.LEHB1
	.uleb128 .L14-.LFB1454
	.uleb128 0x0
	.uleb128 .LEHB2-.LFB1454
	.uleb128 .LEHE2-.LEHB2
	.uleb128 0x0
	.uleb128 0x0
.LLSDACSE1454:
	.section	.text._ZN3cln10cl_I_div_tD1Ev,"axG",@progbits,_ZN3cln10cl_I_div_tD1Ev,comdat
	.section	.text._Z7pow_modmRN3cln4cl_IES1_,"axG",@progbits,_Z7pow_modmRN3cln4cl_IES1_,comdat
	.weak	_Z7pow_modmRN3cln4cl_IES1_
	.type	_Z7pow_modmRN3cln4cl_IES1_, @function
_Z7pow_modmRN3cln4cl_IES1_:
.LFB1451:
	pushq	%r15
.LCFI4:
	pushq	%r14
.LCFI5:
	pushq	%r13
.LCFI6:
	pushq	%r12
.LCFI7:
	pushq	%rbp
.LCFI8:
	pushq	%rbx
.LCFI9:
	subq	$168, %rsp
.LCFI10:
	movq	%rdi, %r15
	movq	%rsi, %rdi
	movq	%rcx, 8(%rsp)
	movq	(%rdx), %rax
	testb	$7, %al
	jne	.L18
	addl	$1, (%rax)
.L18:
	movq	%rax, 144(%rsp)
	movq	$9, (%r15)
	movq	$1, 128(%rsp)
.LEHB3:
	call	_ZN3cln24cl_I_constructor_from_UQEm
.LEHE3:
	movq	%rax, 128(%rsp)
	movq	$1, 16(%rsp)
	movq	$1, 24(%rsp)
	leaq	144(%rsp), %r13
	leaq	128(%rsp), %rbp
	leaq	64(%rsp), %r14
.L69:
	movq	$1, 112(%rsp)
	leaq	112(%rsp), %rsi
	movq	%r13, %rdi
.LEHB4:
	call	_ZN3cln7compareERKNS_4cl_IES2_
.LEHE4:
	movl	%eax, %ebx
	movq	112(%rsp), %rdi
	testb	$7, %dil
	jne	.L40
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L73
.L40:
	testl	%ebx, %ebx
	jle	.L74
.L42:
	movq	%rbp, %r12
	movq	%r13, %rdi
.LEHB5:
	call	_ZN3cln4oddpERKNS_4cl_IE
	testb	%al, %al
	jne	.L70
	movq	%rbp, %r12
.L22:
	movq	$-1, %rdx
	movq	%r13, %rsi
	leaq	32(%rsp), %rdi
	call	_ZN3cln3ashERKNS_4cl_IEl
	movq	32(%rsp), %rbx
	testb	$7, %bl
	jne	.L29
	addl	$1, (%rbx)
.L29:
	movq	144(%rsp), %rdi
	testb	$7, %dil
	jne	.L30
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L75
.L30:
	movq	%rbx, 144(%rsp)
	movq	32(%rsp), %rdi
	testb	$7, %dil
	jne	.L31
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L76
.L31:
	movq	%rbp, %rdx
	movq	%rbp, %rsi
	movq	%r14, %rdi
	call	_ZN3clnmlERKNS_4cl_IES2_
.LEHE5:
	movq	8(%rsp), %rdx
	movq	%r14, %rsi
	leaq	48(%rsp), %rdi
.LEHB6:
	call	_ZN3cln3modERKNS_4cl_IES2_
.LEHE6:
	movq	48(%rsp), %rbx
	testb	$7, %bl
	jne	.L34
	addl	$1, (%rbx)
.L34:
	movq	128(%rsp), %rdi
	testb	$7, %dil
	jne	.L35
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L77
.L35:
	movq	%rbx, 128(%rsp)
	movq	48(%rsp), %rdi
	testb	$7, %dil
	jne	.L37
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L78
.L37:
	movq	64(%rsp), %rdi
	testb	$7, %dil
	jne	.L69
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	jne	.L69
.LEHB7:
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
	jmp	.L69
.L70:
	movq	%rbp, %r12
	movq	%rbp, %rdx
	movq	%r15, %rsi
	leaq	96(%rsp), %rdi
	call	_ZN3clnmlERKNS_4cl_IES2_
.LEHE7:
	movq	8(%rsp), %rdx
	leaq	96(%rsp), %rsi
	leaq	80(%rsp), %rdi
.LEHB8:
	call	_ZN3cln3modERKNS_4cl_IES2_
.LEHE8:
	movq	80(%rsp), %rbx
	testb	$7, %bl
	jne	.L23
	addl	$1, (%rbx)
.L23:
	movq	(%r15), %rdi
	testb	$7, %dil
	jne	.L24
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L79
.L24:
	movq	%rbx, (%r15)
	movq	80(%rsp), %rdi
	testb	$7, %dil
	jne	.L26
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L80
.L26:
	movq	96(%rsp), %rdi
	testb	$7, %dil
	jne	.L22
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	jne	.L22
.LEHB9:
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
	jmp	.L22
.L73:
	movq	%rbp, %r12
	.p2align 4,,6
	.p2align 3
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.LEHE9:
	testl	%ebx, %ebx
	.p2align 4,,4
	.p2align 3
	jg	.L42
.L74:
	movq	24(%rsp), %rdi
	testb	$7, %dil
	jne	.L45
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L81
.L45:
	movq	16(%rsp), %rdi
	testb	$7, %dil
	jne	.L47
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L82
.L47:
	movq	128(%rsp), %rdi
	testb	$7, %dil
	jne	.L49
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L83
.L49:
	movq	144(%rsp), %rdi
	testb	$7, %dil
	jne	.L17
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L84
.L17:
	movq	%r15, %rax
	addq	$168, %rsp
	popq	%rbx
	popq	%rbp
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	ret
.L78:
.LEHB10:
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.LEHE10:
	jmp	.L37
.L75:
	.p2align 4,,8
	.p2align 3
.LEHB11:
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.LEHE11:
	.p2align 4,,8
	.p2align 3
	jmp	.L30
.L76:
	.p2align 4,,8
	.p2align 3
.LEHB12:
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.LEHE12:
	.p2align 4,,8
	.p2align 3
	jmp	.L31
.L77:
	.p2align 4,,8
	.p2align 3
.LEHB13:
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.LEHE13:
	.p2align 4,,8
	.p2align 3
	jmp	.L35
.L79:
	.p2align 4,,8
	.p2align 3
.LEHB14:
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.LEHE14:
	.p2align 4,,8
	.p2align 3
	jmp	.L24
.L80:
	.p2align 4,,8
	.p2align 3
.LEHB15:
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.LEHE15:
	.p2align 4,,8
	.p2align 3
	jmp	.L26
.L84:
	.p2align 4,,8
	.p2align 3
.LEHB16:
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.LEHE16:
	.p2align 4,,8
	.p2align 3
	jmp	.L17
.L81:
	.p2align 4,,8
	.p2align 3
.LEHB17:
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.LEHE17:
	.p2align 4,,8
	.p2align 3
	jmp	.L45
.L82:
	.p2align 4,,8
	.p2align 3
.LEHB18:
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.LEHE18:
	.p2align 4,,8
	.p2align 3
	jmp	.L47
.L83:
	.p2align 4,,8
	.p2align 3
.LEHB19:
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.LEHE19:
	.p2align 4,,8
	.p2align 3
	jmp	.L49
.L54:
	movq	%rax, %rbx
.L46:
	leaq	16(%rsp), %rdi
	call	_ZN3cln4cl_ID1Ev
	leaq	128(%rsp), %r12
.L48:
	movq	%r12, %rdi
	call	_ZN3cln4cl_ID1Ev
.L50:
	movq	%r15, %rdi
	call	_ZN3cln4cl_ID1Ev
.L65:
.L51:
	movq	%r13, %rdi
	call	_ZN3cln4cl_ID1Ev
	movq	%rbx, %rdi
.LEHB20:
	call	_Unwind_Resume
.LEHE20:
.L58:
	movq	%rax, %rbx
.L38:
	leaq	48(%rsp), %rdi
	call	_ZN3cln4cl_ID1Ev
.L39:
	movq	%r14, %rdi
	call	_ZN3cln4cl_ID1Ev
.L44:
	leaq	16(%rsp), %rdi
	call	_ZN3cln10cl_I_div_tD1Ev
	jmp	.L48
.L56:
	movq	%rax, %rbx
.L33:
	leaq	32(%rsp), %rdi
	call	_ZN3cln4cl_ID1Ev
	.p2align 4,,3
	.p2align 3
	jmp	.L44
.L61:
	movq	%rax, %rbx
.L28:
	leaq	96(%rsp), %rdi
	call	_ZN3cln4cl_ID1Ev
	.p2align 4,,3
	.p2align 3
	jmp	.L44
.L59:
	movq	%rax, %rbx
	jmp	.L39
.L60:
	movq	%rax, %rbx
.L27:
	leaq	80(%rsp), %rdi
	call	_ZN3cln4cl_ID1Ev
	jmp	.L28
.L55:
	movq	%rax, %rbx
.L20:
	movq	128(%rsp), %rdi
	testb	$7, %dil
	jne	.L21
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	jne	.L21
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.L21:
	leaq	144(%rsp), %r13
	jmp	.L50
.L62:
	movq	%rax, %rbx
	jmp	.L44
.L57:
	movq	%rax, %rbx
.L41:
	leaq	112(%rsp), %rdi
	call	_ZN3cln4cl_ID1Ev
	leaq	128(%rsp), %r12
	jmp	.L44
.L64:
	movq	%rax, %rbx
	jmp	.L50
.L63:
	movq	%rax, %rbx
	leaq	128(%rsp), %r12
	jmp	.L48
.LFE1451:
	.size	_Z7pow_modmRN3cln4cl_IES1_, .-_Z7pow_modmRN3cln4cl_IES1_
	.section	.gcc_except_table
.LLSDA1451:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE1451-.LLSDACSB1451
.LLSDACSB1451:
	.uleb128 .LEHB3-.LFB1451
	.uleb128 .LEHE3-.LEHB3
	.uleb128 .L55-.LFB1451
	.uleb128 0x0
	.uleb128 .LEHB4-.LFB1451
	.uleb128 .LEHE4-.LEHB4
	.uleb128 .L57-.LFB1451
	.uleb128 0x0
	.uleb128 .LEHB5-.LFB1451
	.uleb128 .LEHE5-.LEHB5
	.uleb128 .L62-.LFB1451
	.uleb128 0x0
	.uleb128 .LEHB6-.LFB1451
	.uleb128 .LEHE6-.LEHB6
	.uleb128 .L59-.LFB1451
	.uleb128 0x0
	.uleb128 .LEHB7-.LFB1451
	.uleb128 .LEHE7-.LEHB7
	.uleb128 .L62-.LFB1451
	.uleb128 0x0
	.uleb128 .LEHB8-.LFB1451
	.uleb128 .LEHE8-.LEHB8
	.uleb128 .L61-.LFB1451
	.uleb128 0x0
	.uleb128 .LEHB9-.LFB1451
	.uleb128 .LEHE9-.LEHB9
	.uleb128 .L62-.LFB1451
	.uleb128 0x0
	.uleb128 .LEHB10-.LFB1451
	.uleb128 .LEHE10-.LEHB10
	.uleb128 .L59-.LFB1451
	.uleb128 0x0
	.uleb128 .LEHB11-.LFB1451
	.uleb128 .LEHE11-.LEHB11
	.uleb128 .L56-.LFB1451
	.uleb128 0x0
	.uleb128 .LEHB12-.LFB1451
	.uleb128 .LEHE12-.LEHB12
	.uleb128 .L62-.LFB1451
	.uleb128 0x0
	.uleb128 .LEHB13-.LFB1451
	.uleb128 .LEHE13-.LEHB13
	.uleb128 .L58-.LFB1451
	.uleb128 0x0
	.uleb128 .LEHB14-.LFB1451
	.uleb128 .LEHE14-.LEHB14
	.uleb128 .L60-.LFB1451
	.uleb128 0x0
	.uleb128 .LEHB15-.LFB1451
	.uleb128 .LEHE15-.LEHB15
	.uleb128 .L61-.LFB1451
	.uleb128 0x0
	.uleb128 .LEHB16-.LFB1451
	.uleb128 .LEHE16-.LEHB16
	.uleb128 0x0
	.uleb128 0x0
	.uleb128 .LEHB17-.LFB1451
	.uleb128 .LEHE17-.LEHB17
	.uleb128 .L54-.LFB1451
	.uleb128 0x0
	.uleb128 .LEHB18-.LFB1451
	.uleb128 .LEHE18-.LEHB18
	.uleb128 .L63-.LFB1451
	.uleb128 0x0
	.uleb128 .LEHB19-.LFB1451
	.uleb128 .LEHE19-.LEHB19
	.uleb128 .L64-.LFB1451
	.uleb128 0x0
	.uleb128 .LEHB20-.LFB1451
	.uleb128 .LEHE20-.LEHB20
	.uleb128 0x0
	.uleb128 0x0
.LLSDACSE1451:
	.section	.text._Z7pow_modmRN3cln4cl_IES1_,"axG",@progbits,_Z7pow_modmRN3cln4cl_IES1_,comdat
	.text
.globl _Z7proth_pmmm
	.type	_Z7proth_pmmm, @function
_Z7proth_pmmm:
.LFB1455:
	pushq	%r15
.LCFI11:
	pushq	%r14
.LCFI12:
	pushq	%r13
.LCFI13:
	pushq	%r12
.LCFI14:
	pushq	%rbp
.LCFI15:
	pushq	%rbx
.LCFI16:
	subq	$120, %rsp
.LCFI17:
	movq	%rsi, %rbx
	movq	%rdx, %r15
	movq	$1, 96(%rsp)
.LEHB21:
	call	_ZN3cln24cl_I_constructor_from_UQEm
.LEHE21:
	movq	%rax, 96(%rsp)
	leaq	32(%rsp), %rbp
	leaq	96(%rsp), %r13
	movq	%rbx, %rdx
	movq	%r13, %rsi
	movq	%rbp, %rdi
.LEHB22:
	call	_ZN3cln3ashERKNS_4cl_IEl
	movq	32(%rsp), %rbx
	testb	$7, %bl
	jne	.L89
	addl	$1, (%rbx)
.L89:
	movq	96(%rsp), %rdi
	testb	$7, %dil
	jne	.L90
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L122
.L90:
	movq	%rbx, 96(%rsp)
	movq	32(%rsp), %rdi
	testb	$7, %dil
	jne	.L91
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L123
.L91:
	leaq	80(%rsp), %r14
	movq	$-1, %rdx
	movq	%r13, %rsi
	movq	%r14, %rdi
	call	_ZN3cln3ashERKNS_4cl_IEl
.LEHE22:
	leaq	16(%rsp), %rbp
	movq	%r13, %rsi
	movq	%rbp, %rdi
.LEHB23:
	call	_ZN3cln5plus1ERKNS_4cl_IE
	movq	16(%rsp), %rbx
	testb	$7, %bl
	jne	.L94
	addl	$1, (%rbx)
.L94:
	movq	96(%rsp), %rdi
	testb	$7, %dil
	jne	.L95
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L124
.L95:
	movq	%rbx, 96(%rsp)
	movq	16(%rsp), %rdi
	testb	$7, %dil
	jne	.L97
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L125
.L97:
	leaq	64(%rsp), %r12
	movq	%r13, %rcx
	movq	%r14, %rdx
	movq	%r15, %rsi
	movq	%r12, %rdi
	call	_Z7pow_modmRN3cln4cl_IES1_
.LEHE23:
	movq	$9, (%rsp)
	leaq	48(%rsp), %rbp
	movq	%rsp, %rdx
	movq	%r12, %rsi
	movq	%rbp, %rdi
.LEHB24:
	call	_ZN3clnplERKNS_4cl_IES2_
.LEHE24:
	movq	(%rsp), %rdi
	testb	$7, %dil
	jne	.L101
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L126
.L101:
	movq	%r13, %rsi
	movq	%rbp, %rdi
.LEHB25:
	call	_ZN3cln5equalERKNS_4cl_IES2_
.LEHE25:
	movl	%eax, %ebx
	movq	48(%rsp), %rdi
	testb	$7, %dil
	jne	.L102
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L127
.L102:
	movq	64(%rsp), %rdi
	testb	$7, %dil
	jne	.L105
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L128
.L105:
	movq	80(%rsp), %rdi
	testb	$7, %dil
	jne	.L107
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L129
.L107:
	movq	96(%rsp), %rdi
	testb	$7, %dil
	jne	.L109
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	je	.L130
.L109:
	movzbl	%bl, %eax
	addq	$120, %rsp
	popq	%rbx
	popq	%rbp
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	ret
.L130:
.LEHB26:
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.LEHE26:
	jmp	.L109
.L122:
	.p2align 4,,8
	.p2align 3
.LEHB27:
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.LEHE27:
	.p2align 4,,8
	.p2align 3
	jmp	.L90
.L123:
	.p2align 4,,8
	.p2align 3
.LEHB28:
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.LEHE28:
	.p2align 4,,8
	.p2align 3
	jmp	.L91
.L124:
	.p2align 4,,8
	.p2align 3
.LEHB29:
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.LEHE29:
	.p2align 4,,8
	.p2align 3
	jmp	.L95
.L125:
	.p2align 4,,8
	.p2align 3
.LEHB30:
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.LEHE30:
	.p2align 4,,8
	.p2align 3
	jmp	.L97
.L126:
	.p2align 4,,8
	.p2align 3
.LEHB31:
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
	.p2align 4,,8
	.p2align 3
	jmp	.L101
.L127:
	.p2align 4,,8
	.p2align 3
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.LEHE31:
	.p2align 4,,8
	.p2align 3
	jmp	.L102
.L128:
	.p2align 4,,8
	.p2align 3
.LEHB32:
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.LEHE32:
	.p2align 4,,8
	.p2align 3
	jmp	.L105
.L129:
	.p2align 4,,8
	.p2align 3
.LEHB33:
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.LEHE33:
	.p2align 4,,8
	.p2align 3
	jmp	.L107
.L114:
	movq	%rax, %r15
.L93:
	movq	%rbp, %rdi
	call	_ZN3cln4cl_ID1Ev
.L108:
	movq	%r13, %rdi
	call	_ZN3cln4cl_ID1Ev
	movq	%r15, %rdi
.LEHB34:
	call	_Unwind_Resume
.LEHE34:
.L117:
	movq	%rax, %r15
.L106:
	movq	%r14, %rdi
	call	_ZN3cln4cl_ID1Ev
	jmp	.L108
.L115:
	movq	%rax, %r15
.L103:
	movq	%rbp, %rdi
	call	_ZN3cln4cl_ID1Ev
.L104:
	movq	%r12, %rdi
	call	_ZN3cln4cl_ID1Ev
	jmp	.L106
.L112:
	movq	%rax, %r15
.L100:
	movq	%rsp, %rdi
	call	_ZN3cln4cl_ID1Ev
	.p2align 4,,4
	.p2align 3
	jmp	.L104
.L118:
	movq	%rax, %r15
	.p2align 4,,2
	.p2align 3
	jmp	.L108
.L111:
	movq	%rax, %rbx
.L87:
	movq	96(%rsp), %rdi
	testb	$7, %dil
	.p2align 4,,2
	.p2align 3
	jne	.L88
	movl	(%rdi), %eax
	subl	$1, %eax
	movl	%eax, (%rdi)
	testl	%eax, %eax
	jne	.L88
	call	_ZN3cln19cl_free_heap_objectEPNS_7cl_heapE
.L88:
	movq	%rbx, %rdi
.LEHB35:
	call	_Unwind_Resume
.LEHE35:
.L116:
	movq	%rax, %r15
	jmp	.L104
.L113:
	movq	%rax, %r15
.L98:
	movq	%rbp, %rdi
	call	_ZN3cln4cl_ID1Ev
	.p2align 4,,2
	.p2align 3
	jmp	.L106
.LFE1455:
	.size	_Z7proth_pmmm, .-_Z7proth_pmmm
	.section	.gcc_except_table
.LLSDA1455:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE1455-.LLSDACSB1455
.LLSDACSB1455:
	.uleb128 .LEHB21-.LFB1455
	.uleb128 .LEHE21-.LEHB21
	.uleb128 .L111-.LFB1455
	.uleb128 0x0
	.uleb128 .LEHB22-.LFB1455
	.uleb128 .LEHE22-.LEHB22
	.uleb128 .L118-.LFB1455
	.uleb128 0x0
	.uleb128 .LEHB23-.LFB1455
	.uleb128 .LEHE23-.LEHB23
	.uleb128 .L117-.LFB1455
	.uleb128 0x0
	.uleb128 .LEHB24-.LFB1455
	.uleb128 .LEHE24-.LEHB24
	.uleb128 .L112-.LFB1455
	.uleb128 0x0
	.uleb128 .LEHB25-.LFB1455
	.uleb128 .LEHE25-.LEHB25
	.uleb128 .L115-.LFB1455
	.uleb128 0x0
	.uleb128 .LEHB26-.LFB1455
	.uleb128 .LEHE26-.LEHB26
	.uleb128 0x0
	.uleb128 0x0
	.uleb128 .LEHB27-.LFB1455
	.uleb128 .LEHE27-.LEHB27
	.uleb128 .L114-.LFB1455
	.uleb128 0x0
	.uleb128 .LEHB28-.LFB1455
	.uleb128 .LEHE28-.LEHB28
	.uleb128 .L118-.LFB1455
	.uleb128 0x0
	.uleb128 .LEHB29-.LFB1455
	.uleb128 .LEHE29-.LEHB29
	.uleb128 .L113-.LFB1455
	.uleb128 0x0
	.uleb128 .LEHB30-.LFB1455
	.uleb128 .LEHE30-.LEHB30
	.uleb128 .L117-.LFB1455
	.uleb128 0x0
	.uleb128 .LEHB31-.LFB1455
	.uleb128 .LEHE31-.LEHB31
	.uleb128 .L116-.LFB1455
	.uleb128 0x0
	.uleb128 .LEHB32-.LFB1455
	.uleb128 .LEHE32-.LEHB32
	.uleb128 .L117-.LFB1455
	.uleb128 0x0
	.uleb128 .LEHB33-.LFB1455
	.uleb128 .LEHE33-.LEHB33
	.uleb128 .L118-.LFB1455
	.uleb128 0x0
	.uleb128 .LEHB34-.LFB1455
	.uleb128 .LEHE34-.LEHB34
	.uleb128 0x0
	.uleb128 0x0
	.uleb128 .LEHB35-.LFB1455
	.uleb128 .LEHE35-.LEHB35
	.uleb128 0x0
	.uleb128 0x0
.LLSDACSE1455:
	.text
	.section	.rodata.str1.8,"aMS",@progbits,1
	.align 8
.LC0:
	.string	"Usgage:- check_proth <infile>."
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC1:
	.string	"rb"
	.section	.rodata.str1.8
	.align 8
.LC2:
	.string	"Failed to open %s for binary input. Exiting.\n"
	.align 8
.LC3:
	.string	"First h chosen=%lu was too far from starting h=%lu,\n"
	.align 8
.LC4:
	.string	"Proth number at %lu*2^%lu+1 with a=%lu failed.\n"
	.section	.rodata.str1.1
.LC5:
	.string	"Error reading infile."
	.section	.rodata.str1.8
	.align 8
.LC6:
	.string	"Gap from h=%lu to h=%lu too large for n=%lu.\n"
	.align 8
.LC7:
	.string	"Final prime at or after h=%lu was too far from end at h=%lu for n=%lu.\n"
	.text
.globl main
	.type	main, @function
main:
.LFB1456:
	pushq	%r15
.LCFI18:
	pushq	%r14
.LCFI19:
	pushq	%r13
.LCFI20:
	pushq	%r12
.LCFI21:
	pushq	%rbp
.LCFI22:
	pushq	%rbx
.LCFI23:
	subq	$56, %rsp
.LCFI24:
	cmpl	$2, %edi
	jne	.L144
	movq	8(%rsi), %rdi
	movl	$.LC1, %esi
	call	fopen
	movq	%rax, %r12
	testq	%rax, %rax
	je	.L145
	leaq	48(%rsp), %rdi
	movq	%rax, %rcx
	movl	$1, %edx
	movl	$8, %esi
	call	fread
	leaq	40(%rsp), %rdi
	movq	%r12, %rcx
	movl	$1, %edx
	movl	$8, %esi
	call	fread
	leaq	32(%rsp), %rdi
	movq	%r12, %rcx
	movl	$1, %edx
	movl	$8, %esi
	call	fread
	leaq	24(%rsp), %r15
	movq	%r12, %rcx
	movl	$1, %edx
	movl	$8, %esi
	movq	%r15, %rdi
	call	fread
	leaq	16(%rsp), %r14
	movq	%r12, %rcx
	movl	$1, %edx
	movl	$8, %esi
	movq	%r14, %rdi
	call	fread
	movq	48(%rsp), %rcx
	movabsq	$4000000000000000000, %r13
	sarq	%cl, %r13
	movq	%r13, %rax
	shrq	%rax
	movq	%rax, (%rsp)
	movq	16(%rsp), %rbp
	movq	40(%rsp), %rdx
	movq	%rbp, %rax
	subq	%rdx, %rax
	cmpq	%rax, (%rsp)
	jb	.L146
	movq	24(%rsp), %rdx
	testq	%rdx, %rdx
	je	.L142
	jmp	.L152
.L138:
	movq	%rbx, %rbp
.L142:
	movq	%r12, %rcx
	movl	$1, %edx
	movl	$8, %esi
	movq	%r15, %rdi
	call	fread
	subq	$1, %rax
	jne	.L148
	movq	%r12, %rcx
	movl	$1, %edx
	movl	$8, %esi
	movq	%r14, %rdi
	call	fread
	subq	$1, %rax
	jne	.L149
	movq	16(%rsp), %rbx
	movq	%rbx, %rax
	subq	%rbp, %rax
	cmpq	%rax, %r13
	jb	.L150
	movq	24(%rsp), %rdx
	testq	%rdx, %rdx
	je	.L138
	movq	48(%rsp), %rsi
	movq	%rbx, %rdi
	call	_Z7proth_pmmm
	testl	%eax, %eax
	jne	.L138
.L143:
	movq	24(%rsp), %rcx
	movq	48(%rsp), %rdx
	movq	16(%rsp), %rsi
	movl	$.LC4, %edi
	call	printf
	xorl	%edi, %edi
	call	exit
.L148:
	movq	32(%rsp), %rdx
	movq	%rdx, %rax
	subq	(%rsp), %rax
	cmpq	%rbp, %rax
	ja	.L151
	movq	%r12, %rdi
	call	fclose
	xorl	%eax, %eax
	addq	$56, %rsp
	popq	%rbx
	popq	%rbp
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	ret
.L152:
	movq	%rcx, %rsi
	movq	%rbp, %rdi
	call	_Z7proth_pmmm
	testl	%eax, %eax
	jne	.L142
	.p2align 4,,2
	.p2align 3
	jmp	.L143
.L149:
	movl	$.LC5, %edi
	call	puts
	xorl	%edi, %edi
	call	exit
.L150:
	movq	48(%rsp), %rcx
	movq	%rbx, %rdx
	movq	%rbp, %rsi
	movl	$.LC6, %edi
	xorl	%eax, %eax
	call	printf
	xorl	%edi, %edi
	call	exit
.L144:
	movl	$.LC0, %edi
	call	puts
	xorl	%edi, %edi
	call	exit
.L146:
	movq	%rbp, %rsi
	movl	$.LC3, %edi
	xorl	%eax, %eax
	call	printf
	xorl	%edi, %edi
	call	exit
.L145:
	movl	$.LC2, %edi
	xorl	%eax, %eax
	call	printf
	xorl	%edi, %edi
	call	exit
.L151:
	movq	48(%rsp), %rcx
	movq	%rbp, %rsi
	movl	$.LC7, %edi
	xorl	%eax, %eax
	call	printf
	xorl	%edi, %edi
	call	exit
.LFE1456:
	.size	main, .-main
	.local	_ZStL8__ioinit
	.comm	_ZStL8__ioinit,1,1
	.local	_ZN3clnL34cl_random_def_init_helper_instanceE
	.comm	_ZN3clnL34cl_random_def_init_helper_instanceE,1,1
	.local	_ZN3clnL27cl_I_classes_dummy_instanceE
	.comm	_ZN3clnL27cl_I_classes_dummy_instanceE,1,1
	.local	_ZN3clnL36cl_prin_globals_init_helper_instanceE
	.comm	_ZN3clnL36cl_prin_globals_init_helper_instanceE,1,1
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
	.section	.eh_frame,"a",@progbits
.Lframe1:
	.long	.LECIE1-.LSCIE1
.LSCIE1:
	.long	0x0
	.byte	0x1
	.string	"zPLR"
	.uleb128 0x1
	.sleb128 -8
	.byte	0x10
	.uleb128 0x7
	.byte	0x3
	.long	__gxx_personality_v0
	.byte	0x3
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
	.long	.LFB1535
	.long	.LFE1535-.LFB1535
	.uleb128 0x4
	.long	0x0
	.byte	0x4
	.long	.LCFI0-.LFB1535
	.byte	0xe
	.uleb128 0x10
	.align 8
.LEFDE1:
.LSFDE3:
	.long	.LEFDE3-.LASFDE3
.LASFDE3:
	.long	.LASFDE3-.Lframe1
	.long	.LFB1339
	.long	.LFE1339-.LFB1339
	.uleb128 0x4
	.long	0x0
	.align 8
.LEFDE3:
.LSFDE5:
	.long	.LEFDE5-.LASFDE5
.LASFDE5:
	.long	.LASFDE5-.Lframe1
	.long	.LFB1454
	.long	.LFE1454-.LFB1454
	.uleb128 0x4
	.long	.LLSDA1454
	.byte	0x4
	.long	.LCFI1-.LFB1454
	.byte	0xe
	.uleb128 0x10
	.byte	0x4
	.long	.LCFI2-.LCFI1
	.byte	0xe
	.uleb128 0x18
	.byte	0x4
	.long	.LCFI3-.LCFI2
	.byte	0xe
	.uleb128 0x20
	.byte	0x83
	.uleb128 0x3
	.byte	0x86
	.uleb128 0x2
	.align 8
.LEFDE5:
.LSFDE7:
	.long	.LEFDE7-.LASFDE7
.LASFDE7:
	.long	.LASFDE7-.Lframe1
	.long	.LFB1451
	.long	.LFE1451-.LFB1451
	.uleb128 0x4
	.long	.LLSDA1451
	.byte	0x4
	.long	.LCFI4-.LFB1451
	.byte	0xe
	.uleb128 0x10
	.byte	0x4
	.long	.LCFI5-.LCFI4
	.byte	0xe
	.uleb128 0x18
	.byte	0x4
	.long	.LCFI6-.LCFI5
	.byte	0xe
	.uleb128 0x20
	.byte	0x4
	.long	.LCFI7-.LCFI6
	.byte	0xe
	.uleb128 0x28
	.byte	0x4
	.long	.LCFI8-.LCFI7
	.byte	0xe
	.uleb128 0x30
	.byte	0x4
	.long	.LCFI9-.LCFI8
	.byte	0xe
	.uleb128 0x38
	.byte	0x4
	.long	.LCFI10-.LCFI9
	.byte	0xe
	.uleb128 0xe0
	.byte	0x83
	.uleb128 0x7
	.byte	0x86
	.uleb128 0x6
	.byte	0x8c
	.uleb128 0x5
	.byte	0x8d
	.uleb128 0x4
	.byte	0x8e
	.uleb128 0x3
	.byte	0x8f
	.uleb128 0x2
	.align 8
.LEFDE7:
.LSFDE9:
	.long	.LEFDE9-.LASFDE9
.LASFDE9:
	.long	.LASFDE9-.Lframe1
	.long	.LFB1455
	.long	.LFE1455-.LFB1455
	.uleb128 0x4
	.long	.LLSDA1455
	.byte	0x4
	.long	.LCFI11-.LFB1455
	.byte	0xe
	.uleb128 0x10
	.byte	0x4
	.long	.LCFI12-.LCFI11
	.byte	0xe
	.uleb128 0x18
	.byte	0x4
	.long	.LCFI13-.LCFI12
	.byte	0xe
	.uleb128 0x20
	.byte	0x4
	.long	.LCFI14-.LCFI13
	.byte	0xe
	.uleb128 0x28
	.byte	0x4
	.long	.LCFI15-.LCFI14
	.byte	0xe
	.uleb128 0x30
	.byte	0x4
	.long	.LCFI16-.LCFI15
	.byte	0xe
	.uleb128 0x38
	.byte	0x4
	.long	.LCFI17-.LCFI16
	.byte	0xe
	.uleb128 0xb0
	.byte	0x83
	.uleb128 0x7
	.byte	0x86
	.uleb128 0x6
	.byte	0x8c
	.uleb128 0x5
	.byte	0x8d
	.uleb128 0x4
	.byte	0x8e
	.uleb128 0x3
	.byte	0x8f
	.uleb128 0x2
	.align 8
.LEFDE9:
.LSFDE11:
	.long	.LEFDE11-.LASFDE11
.LASFDE11:
	.long	.LASFDE11-.Lframe1
	.long	.LFB1456
	.long	.LFE1456-.LFB1456
	.uleb128 0x4
	.long	0x0
	.byte	0x4
	.long	.LCFI18-.LFB1456
	.byte	0xe
	.uleb128 0x10
	.byte	0x4
	.long	.LCFI19-.LCFI18
	.byte	0xe
	.uleb128 0x18
	.byte	0x4
	.long	.LCFI20-.LCFI19
	.byte	0xe
	.uleb128 0x20
	.byte	0x4
	.long	.LCFI21-.LCFI20
	.byte	0xe
	.uleb128 0x28
	.byte	0x4
	.long	.LCFI22-.LCFI21
	.byte	0xe
	.uleb128 0x30
	.byte	0x4
	.long	.LCFI23-.LCFI22
	.byte	0xe
	.uleb128 0x38
	.byte	0x4
	.long	.LCFI24-.LCFI23
	.byte	0xe
	.uleb128 0x70
	.byte	0x83
	.uleb128 0x7
	.byte	0x86
	.uleb128 0x6
	.byte	0x8c
	.uleb128 0x5
	.byte	0x8d
	.uleb128 0x4
	.byte	0x8e
	.uleb128 0x3
	.byte	0x8f
	.uleb128 0x2
	.align 8
.LEFDE11:
	.ident	"GCC: (GNU) 4.3.3"
	.section	.note.GNU-stack,"",@progbits
