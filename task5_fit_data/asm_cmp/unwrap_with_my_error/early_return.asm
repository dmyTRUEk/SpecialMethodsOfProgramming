example::unwrap_with_my_error:
        push    rbx
        mov     rbx, rdi
        cmp     rsi, 1      ; check if `rsi` is 1
        jne     .LBB0_1
        movabs  rax, -3689348814741910323
        mulx    rax, rax, rax
        shr     rax, 3
        mov     qword ptr [rbx], rax
        mov     qword ptr [rbx + 8], 0
        mov     rax, rbx
        pop     rbx
        ret
.LBB0_1:
        mov     edi, 15
        mov     esi, 1
        call    qword ptr [rip + __rust_alloc@GOTPCREL]
        test    rax, rax
        je      .LBB0_5
        movabs  rcx, 8245935278387129711
        mov     qword ptr [rax + 7], rcx
        movabs  rcx, 8031170983519877485
        mov     qword ptr [rax], rcx
        mov     qword ptr [rbx], 15
        mov     qword ptr [rbx + 8], rax
        mov     qword ptr [rbx + 16], 15
        mov     rax, rbx
        pop     rbx
        ret
.LBB0_5:
        mov     edi, 15
        mov     esi, 1
        call    qword ptr [rip + alloc::alloc::handle_alloc_error@GOTPCREL]
        ud2
