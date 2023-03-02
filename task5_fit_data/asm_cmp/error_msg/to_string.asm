example::e:
        push    rbx
        mov     rbx, rdi
        mov     edi, 20
        mov     esi, 1
        call    qword ptr [rip + __rust_alloc@GOTPCREL]
        test    rax, rax
        je      .LBB0_1
        vmovups xmm0, xmmword ptr [rip + .L__unnamed_1]
        vmovups xmmword ptr [rax], xmm0
        mov     dword ptr [rax + 16], 1735617824
        mov     qword ptr [rbx], 20
        mov     qword ptr [rbx + 8], rax
        mov     qword ptr [rbx + 16], 20
        mov     rax, rbx
        pop     rbx
        ret
.LBB0_1:
        mov     edi, 20
        mov     esi, 1
        call    qword ptr [rip + alloc::alloc::handle_alloc_error@GOTPCREL]
        ud2

.L__unnamed_1:
        .ascii  "this is my error msg"
