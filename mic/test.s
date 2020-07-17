# mark_description "Intel(R) C Intel(R) 64 Compiler XE for applications running on Intel(R) MIC Architecture, Version 14.0.1.106";
# mark_description " Build 20131008";
# mark_description "-S -mmic";
	.file "test.cpp"
	.text
..TXTST0:
# -- Begin  main
# mark_begin;
# Threads 4
        .align    16,0x90
	.globl main
main:
# parameter 1: %edi
# parameter 2: %rsi
..B1.1:                         # Preds ..B1.0 Latency 21
..___tag_value_main.1:                                          #5.1
        pushq     %rbp                                          #5.1
..___tag_value_main.3:                                          #
        movq      %rsp, %rbp                                    #5.1
..___tag_value_main.4:                                          #
        andq      $-128, %rsp                                   #5.1
        pushq     %r14                                          #5.1 c1
        pushq     %r15                                          #5.1 c5
        subq      $112, %rsp                                    #5.1 c5
..___tag_value_main.6:                                          #
        movq      %rsi, %r15                                    #5.1 c9
        movl      %edi, %r14d                                   #5.1 c9
        movq      $0x004000002, %rsi                            #5.1 c13
        movl      $3, %edi                                      #5.1 c17
        call      __intel_new_feature_proc_init                 #5.1 c21
                                # LOE rbx r12 r13 r15 r14d
..B1.12:                        # Preds ..B1.1 Latency 7
        stmxcsr   (%rsp)                                        #5.1 c1
        orl       $32832, (%rsp)                                #5.1 c2
        ldmxcsr   (%rsp)                                        #5.1 c6
        cmpl      $2, %r14d                                     #6.12 c7
        jne       ..B1.9        # Prob 1%                       #6.12 c7
                                # LOE rbx r12 r13 r15
..B1.2:                         # Preds ..B1.12 Latency 5
        movq      8(%r15), %rdi                                 #9.16 c1
        call      atol                                          #9.16 c5
                                # LOE rax rbx r12 r13
..B1.13:                        # Preds ..B1.2 Latency 1
        movq      %rax, %r8                                     #9.16 c1
                                # LOE rbx r8 r12 r13
..B1.3:                         # Preds ..B1.13 Latency 5
        xorl      %esi, %esi                                    #10.24 c1
        xorl      %edi, %edi                                    #11.17 c1
        testq     %r8, %r8                                      #11.22 c5
        jle       ..B1.7        # Prob 10%                      #11.22 c5
                                # LOE rbx rsi rdi r8 r12 r13
..B1.4:                         # Preds ..B1.3 Latency 1
        movq      $0xd96ed1192687859b, %rcx                     #12.5 c1
                                # LOE rcx rbx rsi rdi r8 r12 r13
..B1.5:                         # Preds ..B1.5 ..B1.4 Latency 43
        addq      %rdi, %rsi                                    #12.14 c1
        incq      %rdi                                          #11.26 c1
        movq      %rsi, %rax                                    #12.17 c5
        mulq      %rcx                                          #12.17 c9
        shrq      $20, %rdx                                     #12.17 c21
        imulq     $-1234567, %rdx, %r9                          #12.17 c25
        addq      %r9, %rsi                                     #12.17 c39
        cmpq      %r8, %rdi                                     #11.22 c43
        jl        ..B1.5        # Prob 82%                      #11.22 c43
                                # LOE rcx rbx rsi rdi r8 r12 r13
..B1.7:                         # Preds ..B1.5 ..B1.3 Latency 5
        movl      $.L_2__STRING.0, %edi                         #14.3 c1
        xorl      %eax, %eax                                    #14.3 c1
..___tag_value_main.8:                                          #14.3
        call      printf                                        #14.3
..___tag_value_main.9:                                          #
                                # LOE rbx r12 r13
..B1.8:                         # Preds ..B1.7 Latency 13
        xorl      %eax, %eax                                    #16.1 c1
        addq      $112, %rsp                                    #16.1 c1
..___tag_value_main.10:                                         #16.1
        popq      %r15                                          #16.1
..___tag_value_main.11:                                         #16.1
        popq      %r14                                          #16.1
        movq      %rbp, %rsp                                    #16.1 c13
        popq      %rbp                                          #16.1
..___tag_value_main.12:                                         #
        ret                                                     #16.1
..___tag_value_main.14:                                         #
                                # LOE
..B1.9:                         # Preds ..B1.12                 # Infreq Latency 5
        xorl      %edi, %edi                                    #7.5 c1
        call      exit                                          #7.5 c5
        .align    16,0x90
..___tag_value_main.18:                                         #
                                # LOE
# mark_end;
	.type	main,@function
	.size	main,.-main
	.data
# -- End  main
	.section .rodata.str1.4, "aMS",@progbits,1
	.align 4
	.align 4
.L_2__STRING.0:
	.long	1970496850
	.long	1998615660
	.long	622883681
	.long	681068
	.type	.L_2__STRING.0,@object
	.size	.L_2__STRING.0,16
	.data
	.section .note.GNU-stack, ""
// -- Begin DWARF2 SEGMENT .eh_frame
	.section .eh_frame,"a",@progbits
.eh_frame_seg:
	.align 8
	.4byte 0x0000001c
	.8byte 0x00507a0100000000
	.4byte 0x09107801
	.byte 0x00
	.8byte __gxx_personality_v0
	.4byte 0x9008070c
	.2byte 0x0001
	.byte 0x00
	.4byte 0x00000094
	.4byte 0x00000024
	.8byte ..___tag_value_main.1
	.8byte ..___tag_value_main.18-..___tag_value_main.1
	.2byte 0x0400
	.4byte ..___tag_value_main.3-..___tag_value_main.1
	.2byte 0x100e
	.byte 0x04
	.4byte ..___tag_value_main.4-..___tag_value_main.3
	.4byte 0x8610060c
	.2byte 0x0402
	.4byte ..___tag_value_main.6-..___tag_value_main.4
	.8byte 0xff800d1c380e0e10
	.8byte 0xfffffff80d1affff
	.8byte 0x800d1c380e0f1022
	.8byte 0xfffff00d1affffff
	.2byte 0x22ff
	.byte 0x04
	.4byte ..___tag_value_main.10-..___tag_value_main.6
	.2byte 0x04cf
	.4byte ..___tag_value_main.11-..___tag_value_main.10
	.2byte 0x04ce
	.4byte ..___tag_value_main.12-..___tag_value_main.11
	.4byte 0xc608070c
	.byte 0x04
	.4byte ..___tag_value_main.14-..___tag_value_main.12
	.8byte 0x0e0e10028610060c
	.8byte 0x1affffff800d1c38
	.8byte 0x0f1022fffffff80d
	.8byte 0xffffff800d1c380e
	.8byte 0x0022fffffff00d1a
	.4byte 0x00000000
	.byte 0x00
# End
