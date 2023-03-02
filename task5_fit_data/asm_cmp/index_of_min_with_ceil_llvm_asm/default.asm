.LCPI0_0:
        .quad   0x7fffffffffffffff
        .quad   0x7fffffffffffffff
.LCPI0_1:
        .quad   0x7ff0000000000000
<alloc::vec::Vec<f64> as example::IndexOfMinWithFloor<f64>>::index_of_min_with_ceil:
        push    rax
        mov     rsi, qword ptr [rdi + 16]
        test    rsi, rsi
        je      .LBB0_1
        mov     r9, qword ptr [rdi + 8]
        mov     r10b, 1
        xor     eax, eax
        movapd  xmm1, xmmword ptr [rip + .LCPI0_0]
        movsd   xmm2, qword ptr [rip + .LCPI0_1]
        xor     r8d, r8d
        mov     rcx, r8
        jmp     .LBB0_3
.LBB0_9:
        inc     rcx
        cmp     rcx, rsi
        jae     .LBB0_10
.LBB0_3:
        cmp     rcx, rsi
        jae     .LBB0_15
        movsd   xmm3, qword ptr [r9 + 8*rcx]
        ucomisd xmm3, xmm0
        jae     .LBB0_9
        movapd  xmm4, xmm3
        andpd   xmm4, xmm1
        ucomisd xmm2, xmm4
        jbe     .LBB0_9
        mov     rdx, rcx
        test    r10b, 1
        jne     .LBB0_13
        cmp     rdi, rsi
        jae     .LBB0_8
        movsd   xmm4, qword ptr [r9 + 8*rdi]
        ucomisd xmm4, xmm3
        mov     rdx, rdi
        jbe     .LBB0_13
        mov     rdx, rcx
.LBB0_13:
        inc     rcx
        mov     eax, 1
        xor     r10d, r10d
        mov     rdi, rdx
        mov     r8, rcx
        cmp     rcx, rsi
        jb      .LBB0_3
        pop     rcx
        ret
.LBB0_10:
        mov     rdx, rdi
        pop     rcx
        ret
.LBB0_1:
        xor     eax, eax
        pop     rcx
        ret
.LBB0_15:
        cmp     r8, rsi
        cmovbe  r8, rsi
        lea     rdx, [rip + .L__unnamed_1]
        mov     rdi, r8
        call    qword ptr [rip + core::panicking::panic_bounds_check@GOTPCREL]
        ud2
.LBB0_8:
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
