.LCPI0_0:
        .quad   0x3ff60221426fe719
.LCPI0_1:
        .quad   0x0000000000000014
example::o:
        xor     eax, eax
        vucomisd        xmm0, qword ptr [rip + .LCPI0_0]
        seta    al
        kmovd   k1, eax
        vmovsd  xmm1, qword ptr [rip + .LCPI0_1]
        vmovsd  xmm1 {k1}, xmm1, xmm0
        vmovapd xmm0, xmm1
        ret
