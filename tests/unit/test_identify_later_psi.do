*! Later-group factor variance and the redundant-parameter count (2026-11-02,
*! user). Spec 4.8 G5: a later group's psi F F is fixed at 1 or free, and any
*! other fixed value is refused, rc 3498. The rank refusal's count states each
*! point's deficiency when the second point's rank is higher. Both came from
*! the Module 3 review.
clear all
set more off
set linesize 255

local srcs : dir "src/mata" files "*.mata"
local srcs : list clean srcs
local types dwlsmv_types.mata
local srcs : list srcs - types
local load src/mata/dwlsmv_types.mata
foreach f of local srcs {
    local load `load' src/mata/`f'
}
foreach f of local load {
    quietly do "`f'"
}

mata:
string scalar ttt_text(real scalar G, string scalar tie, | string scalar x3)
{
    string scalar s
    real scalar g

    s = "groups " + strofreal(G) + char(10)
    for (g = 1; g <= G; g++) {
        s = s + "group " + strofreal(g) + char(10) + "factor F" + char(10) +
            "ind y1 y2 y3 y4" + char(10) + "lambda F y1 free" + char(10) +
            "lambda F y2 free" + char(10) + "lambda F y3 free" + char(10) +
            "lambda F y4 free" + char(10)
        if (g == 1) {
            s = s + "metric unitvariance F" + char(10)
        }
        else {
            s = s + "psi F F fixed 1" + char(10)
            if (g == 3 & x3 != "") {
                s = s + x3
            }
        }
    }
    return(s + tie + char(10))
}

struct dwlsmv_stage1 vector ttt_st(real rowvector K1)
{
    struct dwlsmv_stage1 vector st
    struct dwlsmv_statmap scalar m
    real scalar g, i

    st = dwlsmv_stage1(1, cols(K1))
    for (g = 1; g <= cols(K1); g++) {
        m         = dwlsmv_statmap_build(K1[g] \ J(3, 1, 3))
        st[g].map = m
        st[g].s   = J(m.D, 1, 0.3)
        for (i = 1; i <= m.T; i++) {
            st[g].s[i] = invnormal(m.thr_k[i] / m.K[m.var_i[i]]) + 0.02 * m.var_i[i] + 0.05 * (g - 1)
        }
    }
    return(st)
}

struct dwlsmv_parmap scalar ttt_id(real rowvector K1, string scalar tie)
{
    struct dwlsmv_stage1 vector st
    struct dwlsmv_statmap vector m
    real scalar g

    st = ttt_st(K1)
    m  = dwlsmv_statmap(1, cols(K1))
    for (g = 1; g <= cols(K1); g++) {
        m[g] = st[g].map
    }
    return(dwls_identify(ttt_text(cols(K1), tie), m, st))
}

void tlp_txt(string scalar s)
{
    struct dwlsmv_parmap scalar p
    struct dwlsmv_statmap vector m
    struct dwlsmv_stage1 vector st
    real scalar g

    st = ttt_st((3,3))
    m  = dwlsmv_statmap(1, 2)
    for (g = 1; g <= 2; g++) {
        m[g] = st[g].map
    }
    p = dwls_identify(s, m, st)
    st_numscalar("tlp_npar", p.npar)
}
void tlp_val(string scalar v)
{
    tlp_txt(subinstr(ttt_text(2, ""), "psi F F fixed 1", "psi F F " + v))
}
void tlp_two()
{
    string scalar s, n
    real scalar g

    n = char(10)
    s = "groups 2" + n
    for (g = 1; g <= 2; g++) {
        s = s + "group " + strofreal(g) + n + "factor F1 F2" + n + "ind y1 y2 y3 y4" + n +
            "lambda F1 y1 free" + n + "lambda F1 y2 free" + n +
            "lambda F2 y3 free" + n + "lambda F2 y4 free" + n
        if (g == 1) {
            s = s + "metric unitvariance F1" + n + "metric unitvariance F2" + n
        }
        else {
            s = s + "psi F1 F1 fixed 1" + n + "psi F2 F2 fixed 1" + n
        }
        s = s + "psi F2 F1 free" + n + "beta F2 F1 free" + n
    }
    tlp_txt(s)
}
end

local nfail 0
program tlp_run
    args v
    tempfile lg
    quietly log using "`lg'", text replace name(tlp)
    capture noisily mata: tlp_val("`v'")
    local rc = _rc
    quietly log close tlp
    mata: st_local("hit", strofreal(strpos(invtokens(cat(st_local("lg"))'), "a later group's factor variance") > 0))
    c_local tlp_rc `rc'
    c_local tlp_hit `hit'
end
foreach v in "fixed -1" "fixed .5" "fixed 0" {
    tlp_run "`v'"
    if `tlp_rc' != 3498 | !`tlp_hit' {
        display as error "  psi F F `v' FAILED (rc `tlp_rc', message `tlp_hit')"
        local ++nfail
    }
    else display as text "  psi F F `v': refused by the rule"
}
* free here is rank-refused (loadings untied), so only the rule's silence is asserted
foreach v in "fixed 1" "free" {
    tlp_run "`v'"
    if `tlp_hit' {
        display as error "  psi F F `v' FAILED: the rule fired"
        local ++nfail
    }
    else display as text "  psi F F `v': the rule is silent"
}

tempfile lg
quietly log using "`lg'", text replace name(tlp)
capture noisily mata: tlp_two()
local rc = _rc
quietly log close tlp
mata: st_local("hit", strofreal(strpos(invtokens(cat(st_local("lg"))'), "6 parameter(s) are redundant here and 2 at a second point") > 0))
if `rc' != 3498 | !`hit' {
    display as error "  redundant count FAILED (rc `rc', message `hit')"
    local ++nfail
}
else display as text "  redundant count: each point's deficiency stated"

if `nfail' {
    display as error "test_identify_later_psi: `nfail' failure(s)"
    exit 9
}
display as result "test_identify_later_psi: all checks passed"
exit 0
