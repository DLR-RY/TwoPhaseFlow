
scalar Lp = 0;
scalar L1 = 0;
scalar L2 = 0;
scalar U1 = 0;
scalar U2 = 0;

forAll(alpha1, cellI)
{
    if (Foam::abs(mesh.C()[cellI].x() < Foam::cbrt(mesh.V()[cellI])) && (Foam::abs(mesh.C()[cellI].y()) < Foam::cbrt(mesh.V()[cellI]))
    {
        const unallocLabelList& neigh1 = mesh.cellCells()[cellI];
        scalar alpha_0 = alpha1[cellI];
        scalar alpha_1 = alpha1[neigh1[0]];
        scalar dalpha_z = 0;

        forAll(neigh1, neigh)
        {
            dz = (mesh.C[cellI].z() - mesh.C[neigh1[neigh]].z());
            if(Foam::abs(dz) >= Foam::cbrt(mesh.V()[cellI])/2)
            {
                dalpha_z = (alpha_0 - alpha1[neigh1[neigh]])/dz;

                if(dalpha_z > 0)
                {
                    L1 = mesh.C()[cellI].z();
                    U1 = U[cellI].z();
                }

                else if(dalpha_z < 0)
                {
                    L2 = mesh.C()[cellI].z();
                    U2 = U[cellI].z();
                }

            }
        }
    }
}

Lp = L2 - L1;

