*! dwlsmv -- Module 11: the fit indices of one fit, as locals for posting
*!
*! Written 2026-10-21 (option A; F5 names frozen by the user 2026-10-21).
*!
*!     void dwlsmv_m11_fitidx(struct dwlsmv_stage1 vector st,
*!                            real scalar T3, real scalar D)
*!
*! Sets the caller's locals rmsea_s, cfi_s, tli_s, T3_b, df_T3_b, which
*! dwlsmv.ado posts as e() scalars of the same names. n = N - G.
*! Nothing is printed here; _dwlsmv_header prints the indices.
*!
*! SPDX-License-Identifier: GPL-3.0-or-later
*! SPDX-FileCopyrightText: 2026 Erik Ruzek

version 16

mata:

void dwlsmv_m11_fitidx(struct dwlsmv_stage1 vector st, real scalar T3, real scalar D)
{
    real scalar G, g, N
    real rowvector b, fi

    G = length(st)
    N = 0
    for (g = 1; g <= G; g++) {
        N = N + st[g].N
    }
    b  = dwlsmv_fitidx_base(st)
    fi = dwlsmv_fitidx(T3, D, N - G, G, b[1], b[2])
    st_local("rmsea_s", strofreal(fi[1], "%21.0g"))
    st_local("cfi_s",   strofreal(fi[2], "%21.0g"))
    st_local("tli_s",   strofreal(fi[3], "%21.0g"))
    st_local("T3_b",    strofreal(b[1],  "%21.0g"))
    st_local("df_T3_b", strofreal(b[2],  "%21.0g"))
}

end
