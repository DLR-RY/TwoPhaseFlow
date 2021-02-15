{
    // autoPtr<phaseSystem::massTransferTable>
    //     massTransferPtr(fluid.massTransfer());

    // phaseSystem::massTransferTable&
    //     massTransfer(massTransferPtr());

    PtrList<volScalarField>& Y1 = mixture.phase1().Y();
    PtrList<volScalarField>& Y2 = mixture.phase2().Y();
    // why alphaEff ?
    volScalarField muEff = turbulence.muEff(); //mixture.phase1().rho();
    // volScalarField muEff2 = turbulence.alphaEff(); //mixture.phase2().rho();
    if (!Y1.empty())
    {
        label inertIndex1 = mixture.phase1().inertIndex();
        volScalarField Yt(0.0*Y1[0]);

        forAll(Y1, i)
        {
            tmp<fvScalarMatrix> Y1iEqn(mixture.phase1().YiEqn(Y1[i],muEff));

            if (Y1iEqn.valid())
            {
                Y1iEqn =
                (
                    Y1iEqn
                //  ==
                //    *massTransfer[Y1[i].name()]
                //   + fvOptions(alpha1, rho1, Y1[i])
                );

                Y1iEqn->relax();
                Y1iEqn->solve(mesh.solver("Yi"));

                Y1[i].max(0.0);
                Yt += Y1[i];
            }
        }

        Y1[inertIndex1] = scalar(1) - Yt;
        Y1[inertIndex1].max(0.0);
         
    }
 

    if (!Y2.empty())
    {
        label inertIndex2 = mixture.phase2().inertIndex();
        volScalarField Yt(0.0*Y2[0]);

        forAll(Y2, i)
        {
            tmp<fvScalarMatrix> Y2iEqn(mixture.phase2().YiEqn(Y2[i],muEff));

            if (Y2iEqn.valid())
            {
                // Y2iEqn =
                // (
                //     Y2iEqn
                // //  ==
                // //    *massTransfer[Y2[i].name()]
                // //   + fvOptions(alpha2, rho2, Y2[i])
                // );

                Y2iEqn->relax();
                Y2iEqn->solve(mesh.solver("Yi"));
                Y2[i].max(0.0);
                Yt += Y2[i];
            }
        }

        Y2[inertIndex2] = scalar(1) - Yt;
        Y2[inertIndex2].max(0.0);
    }

    // fluid.massTransfer(); // updates interfacial mass flow rates
}
