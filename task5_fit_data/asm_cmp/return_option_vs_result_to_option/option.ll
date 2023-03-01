.LCPI0_0:
        .quad   0x3ff60221426fe719
example::o:
        xor     eax, eax
        vucomisd        xmm0, qword ptr [rip + .LCPI0_0]
        seta    al
        ret
