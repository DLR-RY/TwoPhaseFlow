

// Update alpha1
#include "alphaSuSp.H"
advector.advect(Sp,(Su + divU*min(alpha1(), scalar(1)))());

// Update rhoPhi
rhoPhi = advector.getRhoPhi(rho1, rho2);
alphaPhi10 = advector.alphaPhi();
rhoCpLPhi = advector.getRhoPhi(rho1*tCpl, 0*rho1*tCpl);
rhoCpVPhi = advector.getRhoPhi(0*rho2*tCpv, rho2*tCpv);

alpha2 = 1.0 - alpha1;

Info<< "Phase-1 volume fraction = "
    << alpha1.weightedAverage(mesh.Vsc()).value()
    << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
    << "  Max(" << alpha1.name() << ") - 1 = " << max(alpha1).value() - 1
    << endl;
