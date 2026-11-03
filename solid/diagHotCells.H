// DIAGNOSI TEMPORANEA (P10): bilancio di energia e di acqua delle celle del solido piu' calde di
// diagHotCellsT [K] (controlDict, 0 = spento). Per ogni sotto-passo accettato stampa i termini di
// TsEqn e pcEqn moltiplicati per V/Amax, cioe' in W/m2 e kg/(m2 s) della faccia piu' grande della cella.
{
    const scalar diagT =
        runTime.controlDict().lookupOrDefault<scalar>("diagHotCellsT", 0);

    if (diagT > 0 && gMax(Ts.primitiveField()) > diagT)
    {
        const scalar dt = runTime.deltaT().value();
        const dimensionedScalar cap_l("cap_l", dimensionSet(0,2,-2,-1,0,0,0), 4182);
        const dimensionedScalar cap_v("cap_v", dimensionSet(0,2,-2,-1,0,0,0), 1880);
        const dimensionedScalar Tref("Tref", dimensionSet(0,0,0,1,0,0,0), 273.15);
        const dimensionedScalar L_v("L_v", dimensionSet(0,2,-2,0,0,0,0), 2.5e6);
        const dimensionedScalar rhol("rhol", dimMass/dimVolume, 1000);

        const surfaceScalarField phiGd
        (
            (fvc::interpolate(Krel, "Krel")*rhol*g) & mesh.Sf()
        );
        const surfaceScalarField phiGTd
        (
            (cap_l*fvc::interpolate(Ts - Tref)*fvc::interpolate(Krel, "Krel")*rhol*g) & mesh.Sf()
        );

        const volScalarField hCond(fvc::laplacian(lambda_m, Ts, "laplacian(Krel,pc)"));
        const volScalarField hKpt(fvc::laplacian(((Ts - Tref)*cap_v + L_v)*K_pt, Ts, "laplacian(Krel,pc)"));
        const volScalarField hLiq(fvc::laplacian((Ts - Tref)*cap_l*Krel, pc, "laplacian(Krel,pc)"));
        const volScalarField hVap(fvc::laplacian(((Ts - Tref)*cap_v + L_v)*K_v, pc, "laplacian(Krel,pc)"));
        const volScalarField hGrav(-fvc::div(phiGTd));
        const volScalarField mPc(fvc::laplacian(Krel + K_v, pc, "laplacian(Krel,pc)"));
        const volScalarField mKpt(fvc::laplacian(K_pt, Ts, "laplacian(Krel,pc)"));
        const volScalarField mGrav(-fvc::div(phiGd));

        // le 4 celle piu' calde del processore sopra la soglia
        DynamicList<label> hot;
        forAll(Ts, celli)
        {
            if (Ts[celli] > diagT) hot.append(celli);
        }
        labelList order;
        sortedOrder(scalarField(UIndirectList<scalar>(Ts.primitiveField(), hot)), order);

        const scalarField& V = mesh.V();
        const surfaceScalarField& magSf = mesh.magSf();
        for (label k = order.size() - 1; k >= max(order.size() - 4, 0); k--)
        {
            const label c = hot[order[k]];

            scalar Amax = 0;
            forAll(mesh.cells()[c], fi)
            {
                const label f = mesh.cells()[c][fi];
                if (mesh.isInternalFace(f))
                {
                    Amax = max(Amax, magSf[f]);
                }
                else
                {
                    const label patchi = mesh.boundaryMesh().whichPatch(f);
                    Amax = max(Amax, magSf.boundaryField()[patchi][f - mesh.boundaryMesh()[patchi].start()]);
                }
            }
            const scalar s = V[c]/Amax;

            const scalar Ct = rho_m[c]*cap_m[c] + ws[c]*4182.0;
            const scalar stor = Ct*(Ts[c] - Ts_old[c])/dt;
            const scalar sumH = hCond[c] + hKpt[c] + hLiq[c] + hVap[c] + hGrav[c];

            Pout<< "diagHot t " << runTime.value() << " dt " << dt << " cella " << c
                << " C " << mesh.C()[c] << " Ts " << Ts[c] << " Ts_old " << Ts_old[c]
                << " ws " << ws[c] << " ws_old " << ws_old[c] << " pc " << pc[c] << nl
                << "    calore W/m2: accumulo " << stor*s
                << " cond " << hCond[c]*s << " Kpt " << hKpt[c]*s
                << " liquido " << hLiq[c]*s << " vapore " << hVap[c]*s
                << " gravita " << hGrav[c]*s << " accumulo-somma " << (stor - sumH)*s << nl
                << "    acqua kg/m2s: dws/dt " << (ws[c] - ws_old[c])/dt*s
                << " Crel*dpc/dt " << Crel[c]*(pc[c] - pc_old[c])/dt*s
                << " flussi pc " << mPc[c]*s << " Kpt " << mKpt[c]*s << " gravita " << mGrav[c]*s << nl
                << "    coeff: lambda " << lambda_m[c] << " K_pt " << K_pt[c] << " K_v " << K_v[c]
                << " Krel " << Krel[c] << " Crel " << Crel[c] << " V/Amax " << s << endl;
        }
    }
}
