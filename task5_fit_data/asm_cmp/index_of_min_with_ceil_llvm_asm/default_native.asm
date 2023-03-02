.LCPI0_0:
        .quad   0x7fffffffffffffff
.LCPI0_1:
        .quad   0x7ff0000000000000
<alloc::vec::Vec<f64> as example::IndexOfMinWithFloor<f64>>::index_of_min_with_ceil:
        push    rax
        mov     rsi, qword ptr [rdi + 16]
        test    rsi, rsi
        je      .LBB0_1
        mov     r8, qword ptr [rdi + 8]
        mov     r10b, 1
        xor     eax, eax
        vmovddup        xmm1, qword ptr [rip + .LCPI0_0]
        vmovsd  xmm2, qword ptr [rip + .LCPI0_1]
        xor     ecx, ecx
.LBB0_3:
        cmp     rcx, rsi
        mov     r9, rsi
        cmova   r9, rcx
        jmp     .LBB0_4
.LBB0_10:
        inc     rcx
        cmp     rcx, rsi
        jae     .LBB0_11
.LBB0_4:
        cmp     r9, rcx
        je      .LBB0_16
        vmovsd  xmm3, qword ptr [r8 + 8*rcx]
        vucomisd        xmm3, xmm0
        jae     .LBB0_10
        vandpd  xmm4, xmm3, xmm1
        vucomisd        xmm2, xmm4
        jbe     .LBB0_10
        mov     rdx, rcx
        test    r10b, 1
        jne     .LBB0_14
        cmp     rdi, rsi
        jae     .LBB0_9
        vmovsd  xmm4, qword ptr [r8 + 8*rdi]
        vucomisd        xmm4, xmm3
        mov     rdx, rdi
        jbe     .LBB0_14
        mov     rdx, rcx
.LBB0_14:
        inc     rcx
        mov     eax, 1
        xor     r10d, r10d
        mov     rdi, rdx
        cmp     rcx, rsi
        jb      .LBB0_3
        pop     rcx
        ret
.LBB0_11:
        mov     rdx, rdi
        pop     rcx
        ret
.LBB0_1:
        xor     eax, eax
        pop     rcx
        ret
.LBB0_16:
        lea     rdx, [rip + .L__unnamed_1]
        mov     rdi, r9
        call    qword ptr [rip + core::panicking::panic_bounds_check@GOTPCREL]
        ud2
.LBB0_9:
        lea     rdx, [rip + .L__unnamed_2]
        call    qword ptr [rip + core::panicking::panic_bounds_check@GOTPCREL]
        ud2

.L__unnamed_3:
        .ascii  "/app/example.rs"

.L__unnamed_1:
        .quad   .L__unnamed_3
        .asciz  "\017\000\000\000\000\000\000\000\b\000\000\000\020\000\000"

.L__unnamed_2:
        .quad   .L__unnamed_3
        .asciz  "\017\000\000\000\000\000\000\000\r\000\000\0001\000\000"
