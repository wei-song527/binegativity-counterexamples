"""Supplementary verification for
"Two-qubit binegativity can increase under local amplitude damping".

Run: python3 verify_binegativity.py
Dependencies: numpy, matplotlib. Rational checks use fractions.Fraction;
eigensolver, grid, and plotting checks use floating-point arithmetic.
The figure and verification_results.json are written beside this script.
See README.txt for the verification scope and reference results.
"""

import sys

if sys.flags.optimize:
    raise RuntimeError("Run without -O, -OO, or PYTHONOPTIMIZE; assertions are required.")

from fractions import Fraction as F
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent


def matrix(rows, denominator=1):
    return np.array([[F(x, denominator) for x in row] for row in rows], dtype=object)


def pt(rho):
    return rho.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)


def negative_part(x):
    vals, vecs = np.linalg.eigh(np.asarray(x, dtype=complex))
    return (vecs * np.maximum(-vals, 0)) @ vecs.conj().T


def n2(rho):
    x = negative_part(pt(rho))
    return float(np.trace(x).real + 2*np.trace(negative_part(pt(x))).real)


def n2_vector(rho):
    vals, vecs = np.linalg.eigh(np.asarray(pt(rho), dtype=complex))
    if vals[0] >= 0:
        return 0.0
    return float(-vals[0] * (1+2*abs(np.linalg.det(vecs[:, 0].reshape(2, 2)))))


def negativity(rho):
    return float(2*np.trace(negative_part(pt(rho))).real)


def concurrence_x(rho):
    return float(2*max(0, abs(rho[0, 3])-np.sqrt(rho[1, 1]*rho[2, 2])))


def damping_x(rho, z):
    """Exact for Fraction-valued rho and z; eta=1-z*z."""
    a, b, c, d = [rho[i, i] for i in range(4)]
    w = rho[0, 3]
    out = np.zeros((4, 4), dtype=object)
    out[0, 0] = a+(1-z*z)*c
    out[1, 1] = b+(1-z*z)*d
    out[2, 2] = z*z*c
    out[3, 3] = z*z*d
    out[0, 3] = out[3, 0] = z*w
    return out


def reset_x(rho, z):
    a, b, c, d = [rho[i, i] for i in range(4)]
    w = rho[0, 3]
    out = np.zeros((4, 4), dtype=object)
    out[0, 0] = a+(1-z*z)*(c+d)
    out[1, 1] = b
    out[2, 2] = z*z*c
    out[3, 3] = z*z*d
    out[0, 3] = out[3, 0] = z*w
    return out


def eigenpair_exact(rho, t, ratio):
    # The normalized negative eigenvector is (0,ratio,-1,0)/sqrt(ratio^2+1).
    v = np.array([F(0), F(ratio), F(-1), F(0)], dtype=object)
    assert np.array_equal(pt(rho) @ v, -t*v)
    cv = F(2*ratio, ratio*ratio+1)
    return t*(1+cv)


def state_checks(rho, full_rank):
    assert sum(rho[i, i] for i in range(4)) == 1
    assert np.array_equal(rho, rho.T)
    # Exact positive-definiteness/semidefiniteness of this block structure.
    assert all(rho[i, i] >= 0 for i in range(4))
    determinant = rho[0, 0]*rho[3, 3]-rho[0, 3]**2
    assert determinant >= 0
    if full_rank:
        assert all(rho[i, i] > 0 for i in range(4)) and determinant > 0
    return float(np.linalg.eigvalsh(np.asarray(rho, dtype=float)).min())


def apply(rho, ks):
    return sum(k @ rho @ k.conj().T for k in ks)


def check_kraus(ks):
    assert np.allclose(sum(k.conj().T@k for k in ks), np.eye(4), atol=2e-15, rtol=0)


def check_case(rho, sigma, z, reset, input_t, output_t, input_ratio, output_ratio,
               expected_input, expected_output, full_rank):
    actual = reset_x(rho, z) if reset else damping_x(rho, z)
    assert np.array_equal(actual, sigma)
    minima = [state_checks(x, full_rank) for x in (rho, sigma)]
    exact_in = eigenpair_exact(rho, input_t, input_ratio)
    exact_out = eigenpair_exact(sigma, output_t, output_ratio)
    assert exact_in == expected_input and exact_out == expected_output
    assert exact_out > exact_in
    f, s = np.asarray(rho, dtype=float), np.asarray(sigma, dtype=float)
    eye = np.eye(2)
    jump = np.array([[0., 1.], [0., 0.]])
    zf = float(z)
    ks = [np.kron(np.diag([1., zf]), eye)]
    if reset:
        ks += [np.sqrt(1-zf*zf)*np.kron(jump, b)
               for b in (np.diag([1., 0.]), jump)]
    else:
        ks += [np.sqrt(1-zf*zf)*np.kron(jump, eye)]
    check_kraus(ks)
    assert np.allclose(apply(f, ks), s, atol=2e-15, rtol=0)
    for x, ex in ((f, exact_in), (s, exact_out)):
        assert abs(n2(x)-float(ex)) < 2e-14
        assert abs(n2(x)-n2_vector(x)) < 2e-14
    assert negativity(s) <= negativity(f)
    assert concurrence_x(s) <= concurrence_x(f)
    return dict(input_N2=str(exact_in), output_N2=str(exact_out),
                gain=str(exact_out-exact_in), ratio=str(exact_out/exact_in),
                minimum_eigenvalues=minima)


def family(p):
    rho = np.zeros((4, 4))
    rho[0, 0] = rho[3, 3] = rho[0, 3] = rho[3, 0] = p/2
    rho[2, 2] = 1-p
    return rho


def h(u):
    s = np.hypot(1, u)
    return u/(s+u)*(1+1/s)  # stable form of sqrt(1+u^2)-u


def main():
    rho = matrix([[8000,0,0,88],[0,3,0,0],[0,0,696,0],[88,0,0,1]],8700)
    sigma = matrix([[622022,0,0,2244],[0,280,0,0],
                    [0,0,6264,0],[2244,0,0,9]],628575)
    local = check_case(rho,sigma,F(6,17),False,F(2,2175),F(156,209525),8,3,
                       F(54,47125),F(1248,1047625),True)
    assert F(local['gain']) == F(618,13619125)
    assert F(local['ratio']) == F(2704,2601)
    assert F(156,209525)/F(2,2175) == F(6786,8381)

    rho0 = matrix([[5,0,0,5],[0,0,0,0],[0,0,24,0],[5,0,0,5]],34)
    sigma0 = matrix([[2029,0,0,225],[0,0,0,0],
                     [0,0,600,0],[225,0,0,125]],2754)
    reset = check_case(rho0,sigma0,F(5,9),True,F(1,34),F(25,918),5,3,
                       F(9,221),F(20,459),False)
    assert F(reset['gain']) == F(1,351)
    d = np.diag([F(1),F(1),F(5,9),F(5,9)])
    retained = d@rho0@d
    q0 = sum(retained[i,i] for i in range(4))
    assert q0 == F(565,1377) and 1-q0 == F(812,1377)
    assert eigenpair_exact(retained,F(25,918),3) == F(20,459)
    reset['probabilities'] = [str(q0),str(1-q0)]

    full_rho = matrix([[990,0,0,495],[0,32,0,0],
                       [0,0,2408,0],[495,0,0,990]],4420)
    full_sigma = reset_x(full_rho,F(4,7))
    # Explicit output for the appendix example, retained here for verification.
    assert np.array_equal(full_sigma,matrix([
        [160644,0,0,13860],[0,1568,0,0],
        [0,0,38528,0],[13860,0,0,15840]],216580))
    full_reset = check_case(full_rho,full_sigma,F(4,7),True,
                            F(67,4420),F(109,7735),5,3,
                            F(603,28730),F(872,38675),True)
    assert F(full_reset['gain']) == F(1567,1005550)
    m0=np.kron(np.diag([1.,4/7]),np.eye(2))
    m1=np.kron(np.diag([0.,np.sqrt(33)/7]),np.eye(2))
    rr=np.asarray(full_rho,dtype=float)
    assert abs(n2(m1@rr@m1.T))<2e-14
    assert abs(n2(m0@rr@m0.T)-float(F(872,38675)))<2e-14

    errors = []
    for p in np.linspace(.001,3/7-.001,301):
        r = (1-p)/p
        z = 4/(3*r)
        a = family(p)
        b = np.asarray(reset_x(a,z),dtype=float)
        u, v = n2(a), n2(b)
        pred = [p/(2*r)*h(r), p/(2*r)*h(r*z)]
        errors.extend([abs(u-pred[0]),abs(v-pred[1])])
        assert v > u
        ms = [np.kron(np.diag([1.,z]),np.eye(2)),
              np.kron(np.diag([0.,np.sqrt(1-z*z)]),np.eye(2))]
        branches = [m@a@m.T for m in ms]
        assert abs(n2(branches[1])) < 2e-14
        assert abs(sum(n2(x) for x in branches)-v) < 2e-14
    assert max(errors) < 2e-14
    # Exact derivative sign is justified in the manuscript. Check its formula.
    for v in np.linspace(.001,1/3,101):
        derivative=(1+v)**2*(1-3*v+v*v-v**3)/(1+v*v)**2
        assert derivative > 0

    # Restriction used in the earlier noise study: equal middle diagonals.
    ew_checks=0
    for p in np.linspace(0,1,11):
        for theta in np.linspace(0,np.pi/2,13):
            psi=np.array([np.cos(theta),0,0,np.sin(theta)])
            a=(1-p)*np.eye(4)/4+p*np.outer(psi,psi)
            assert abs(n2(a)-negativity(a)) < 2e-14
            for eta in (.2,.6,.9):
                b=np.asarray(damping_x(a,np.sqrt(1-eta)),dtype=float)
                assert n2(b) <= n2(a)+2e-14
                ew_checks+=1

    # Figure: normalized spectral quantities under ordinary local damping.
    plt.rcParams.update({'font.family':'serif','font.size':9,
                         'mathtext.fontset':'stix','axes.linewidth':.7,
                         'pdf.fonttype':42,'ps.fonttype':42})
    fig, ax = plt.subplots(figsize=(3.38,2.55))
    eta = np.linspace(0,1,1001)
    initial=np.asarray(rho,dtype=float)
    vals=np.array([[n2(a),negativity(a),concurrence_x(a)]
                   for a in (np.asarray(damping_x(initial,np.sqrt(1-e)),dtype=float)
                             for e in eta)])
    ratios=vals/vals[0]
    ax.axhline(1,color='.75',lw=.7,zorder=0)
    for j,(label,color,style) in enumerate([
          (r'$N_2$','#1565A7','-'),(r'$N$','#C47724','--'),(r'$C$','.35',':')]):
        ax.plot(eta,ratios[:,j],label=label,color=color,ls=style,lw=1.6)
    ex=float(F(253,289)); ey=float(F(2704,2601))
    ax.plot(ex,ey,'o',ms=4,mfc='white',mec='#1565A7',mew=1.2)
    ax.annotate(r'$\eta_*=253/289$',xy=(ex,ey),xytext=(.38,1.065),
                arrowprops=dict(arrowstyle='-',lw=.65,color='.35'),fontsize=8)
    ax.set(xlabel=r'Damping probability $\eta$',ylabel='Value / initial value',
           xlim=(0,1),ylim=(0,1.16))
    ax.set_xticks(np.linspace(0,1,6)); ax.set_yticks(np.linspace(0,1,6))
    ax.legend(frameon=False,loc='lower left',ncol=3,handlelength=1.8,
              columnspacing=1.1,borderaxespad=.5)
    ax.tick_params(direction='in',top=True,right=True,length=3)
    fig.tight_layout(pad=.55)
    fig.savefig(ROOT/'binegativity_damping.pdf')
    fig.savefig(ROOT/'binegativity_damping.png',dpi=220)
    plt.close(fig)
    result={'local_amplitude_damping':local,'one_way_LOCC':reset,
            'full_rank_one_way_LOCC':full_reset,
            'continuous_family_grid_points':301,'max_family_formula_error':max(errors),
            'earlier_input_family_checks':ew_checks,
            'convention':'N=||rho^Gamma||_1-1; N2=Tr[(rho^Gamma)_-]+2Tr[(((rho^Gamma)_-)^Gamma)_-]'}
    (ROOT/'verification_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
