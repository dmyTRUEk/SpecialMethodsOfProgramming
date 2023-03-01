<alloc::vec::Vec<T,A> as core::ops::index::Index<I>>::index:
        cmp     rdx, rsi
        jae     .LBB0_1
        lea     rax, [rdi + 8*rdx]
        ret
.LBB0_1:
        push    rax
        mov     rdi, rdx
        mov     rdx, rcx
        call    qword ptr [rip + core::panicking::panic_bounds_check@GOTPCREL]
        ud2

.LCPI1_0:
        .quad   0x7fffffffffffffff
        .quad   0x7fffffffffffffff
.LCPI1_1:
        .quad   0x7ff0000000000000
<alloc::vec::Vec<f64> as example::IndexOfMinWithFloor<f64>>::index_of_min_with_ceil:
        push    rax
        mov     rsi, qword ptr [rdi + 16]
        test    rsi, rsi
        je      .LBB1_6
        mov     rdi, qword ptr [rdi + 8]
        xor     edx, edx
        movapd  xmm1, xmmword ptr [rip + .LCPI1_0]
        movsd   xmm2, qword ptr [rip + .LCPI1_1]
        jmp     .LBB1_2
.LBB1_5:
        inc     rdx
        cmp     rsi, rdx
        je      .LBB1_6
.LBB1_2:
        cmp     rsi, rdx
        je      .LBB1_8
        movsd   xmm3, qword ptr [rdi + 8*rdx]
        ucomisd xmm3, xmm0
        jae     .LBB1_5
        andpd   xmm3, xmm1
        ucomisd xmm2, xmm3
        jbe     .LBB1_5
        lea     rcx, [rip + .L__unnamed_1]
        call    <alloc::vec::Vec<T,A> as core::ops::index::Index<I>>::index
        lea     rdi, [rip + .L__unnamed_2]
        lea     rdx, [rip + .L__unnamed_3]
        mov     esi, 43
        call    qword ptr [rip + core::panicking::panic@GOTPCREL]
        ud2
.LBB1_6:
        xor     eax, eax
        pop     rcx
        ret
.LBB1_8:
        lea     rdx, [rip + .L__unnamed_4]
        mov     rdi, rsi
        call    qword ptr [rip + core::panicking::panic_bounds_check@GOTPCREL]
        ud2

.L__unnamed_5:
        .ascii  "/app/example.rs"

.L__unnamed_4:
        .quad   .L__unnamed_5
        .asciz  "\017\000\000\000\000\000\000\000\b\000\000\000\020\000\000"

.L__unnamed_1:
        .quad   .L__unnamed_5
        .asciz  "\017\000\000\000\000\000\000\000\016\000\000\000\037\000\000"

.L__unnamed_2:
        .ascii  "called `Option::unwrap()` on a `None` value"

.L__unnamed_3:
        .quad   .L__unnamed_5
        .asciz  "\017\000\000\000\000\000\000\000\016\000\000\000B\000\000"
