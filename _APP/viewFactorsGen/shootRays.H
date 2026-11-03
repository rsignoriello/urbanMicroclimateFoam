// All rays expressed as start face (local) index and end face (global).
// Rays are accumulated in batches of at most maxDynListLength pairs and
// flushed to surfacesMesh.findLine() each time the batch is full.
// This prevents MPI_ERR_COUNT (int32 overflow) on large meshes.
//
// Safe upper bound: point = 3 doubles = 24 bytes.
// INT_MAX / 24 ~ 89M.  Default 10M gives ~240 MB per findLine call.
const label maxDynListLength
(
    viewFactorDict.lookupOrDefault<label>("maxDynListLength", 10000000)
);

const Switch useProcessorAABBCulling
(
    viewFactorDict.lookupOrDefault<Switch>("useProcessorAABBCulling", true)
);

const Switch useProcessorFaceGrid
(
    viewFactorDict.lookupOrDefault<Switch>("useProcessorFaceGrid", true)
);

const label processorFaceGridN
(
    max(viewFactorDict.lookupOrDefault<label>("processorFaceGridN", 12), 1)
);

List<point> remoteMin(Pstream::nProcs(), point(GREAT, GREAT, GREAT));
List<point> remoteMax(Pstream::nProcs(), point(-GREAT, -GREAT, -GREAT));

List<vector> remoteCellDelta(Pstream::nProcs(), Zero);
List<List<DynamicList<label>>> remoteFaceBins(Pstream::nProcs());

for (label procj = 0; procj < Pstream::nProcs(); procj++)
{
    const pointField& fc = remoteCoarseCf[procj];
    if (fc.size())
    {
        point pMin(fc[0]);
        point pMax(fc[0]);
        forAll(fc, k)
        {
            const point& p = fc[k];
            pMin.x() = min(pMin.x(), p.x());
            pMin.y() = min(pMin.y(), p.y());
            pMin.z() = min(pMin.z(), p.z());
            pMax.x() = max(pMax.x(), p.x());
            pMax.y() = max(pMax.y(), p.y());
            pMax.z() = max(pMax.z(), p.z());
        }
        remoteMin[procj] = pMin;
        remoteMax[procj] = pMax;

        if (useProcessorFaceGrid && processorFaceGridN > 1)
        {
            const vector span = pMax - pMin;
            const vector cellDelta
            (
                max(span.x()/processorFaceGridN, SMALL),
                max(span.y()/processorFaceGridN, SMALL),
                max(span.z()/processorFaceGridN, SMALL)
            );
            const vector invCellDelta
            (
                1.0/cellDelta.x(),
                1.0/cellDelta.y(),
                1.0/cellDelta.z()
            );

            remoteCellDelta[procj] = cellDelta;

            const label nCells =
                processorFaceGridN*processorFaceGridN*processorFaceGridN;
            remoteFaceBins[procj].setSize(nCells);

            forAll(fc, facei)
            {
                const point& p = fc[facei];

                const label ix = min
                (
                    max(label((p.x() - pMin.x())*invCellDelta.x()), 0),
                    processorFaceGridN - 1
                );
                const label iy = min
                (
                    max(label((p.y() - pMin.y())*invCellDelta.y()), 0),
                    processorFaceGridN - 1
                );
                const label iz = min
                (
                    max(label((p.z() - pMin.z())*invCellDelta.z()), 0),
                    processorFaceGridN - 1
                );

                const label cellId =
                    ix + processorFaceGridN*(iy + processorFaceGridN*iz);
                remoteFaceBins[procj][cellId].append(facei);
            }
        }
    }
}

for (label proci = 0; proci < Pstream::nProcs(); proci++)
{
    // Shoot rays from me to proci. Note that even if processor has
    // 0 faces we still need to call findLine to keep calls synced.

    const pointField& myFc   = remoteCoarseCf[Pstream::myProcNo()];
    const vectorField& myArea = remoteCoarseSf[Pstream::myProcNo()];
    const labelField& myAgg   = remoteCoarseAgg[Pstream::myProcNo()];

    const pointField& remoteFc   = remoteCoarseCf[proci];
    const pointField& remoteArea = remoteCoarseSf[proci];
    const labelField& remoteAgg  = remoteCoarseAgg[proci];

    // i,j persist across batches so we resume where we left off.
    label i = 0;
    label j = 0;
    label candidateFacei = -1;
    label candidateFacePos = 0;
    DynamicList<label> candidateFaces;

    bool procNotDone = (myFc.size() > 0);

    do
    {
        // --- Accumulate up to maxDynListLength candidate ray pairs ---
        DynamicField<point> start(min(maxDynListLength, coarseMesh.nFaces()));
        DynamicField<point> end(start.capacity());
        DynamicList<label>  startIndex(start.capacity());
        DynamicList<label>  endIndex(start.capacity());
        DynamicList<label>  startAgg(start.capacity());
        DynamicList<label>  endAgg(start.capacity());

        bool batchFull = false;
        while (i < myFc.size() && !batchFull)
        {
            const point&  fc   = myFc[i];
            const vector& fA   = myArea[i];
            const label&  fAgg = myAgg[i];

            if (useProcessorAABBCulling && remoteFc.size())
            {
                const point& pMin = remoteMin[proci];
                const point& pMax = remoteMax[proci];

                scalar minDot = GREAT;
                for (label cx = 0; cx < 2; cx++)
                {
                    const scalar x = (cx ? pMax.x() : pMin.x());
                    for (label cy = 0; cy < 2; cy++)
                    {
                        const scalar y = (cy ? pMax.y() : pMin.y());
                        for (label cz = 0; cz < 2; cz++)
                        {
                            const scalar z = (cz ? pMax.z() : pMin.z());
                            const scalar dDotN = ((point(x, y, z) - fc) & fA);
                            minDot = min(minDot, dDotN);
                        }
                    }
                }

                // Need (remFc-fc)&fA < 0 for any candidate. If even the
                // minimum box-corner projection is >= 0, skip this face.
                if (minDot >= 0)
                {
                    j = 0;
                    i++;
                    candidateFacei = -1;
                    continue;
                }
            }

            if
            (
                useProcessorFaceGrid
             && processorFaceGridN > 1
             && remoteFaceBins[proci].size()
            )
            {
                if (candidateFacei != i)
                {
                    candidateFaces.clear();
                    candidateFacei = i;
                    candidateFacePos = 0;

                    const point& pMin = remoteMin[proci];
                    const vector& cellDelta = remoteCellDelta[proci];

                    for (label cx = 0; cx < processorFaceGridN; cx++)
                    {
                        const scalar x0 = pMin.x() + cx*cellDelta.x();
                        const scalar x1 = x0 + cellDelta.x();
                        for (label cy = 0; cy < processorFaceGridN; cy++)
                        {
                            const scalar y0 = pMin.y() + cy*cellDelta.y();
                            const scalar y1 = y0 + cellDelta.y();
                            for (label cz = 0; cz < processorFaceGridN; cz++)
                            {
                                const scalar z0 = pMin.z() + cz*cellDelta.z();
                                const scalar z1 = z0 + cellDelta.z();

                                scalar minDotCell = GREAT;
                                for (label px = 0; px < 2; px++)
                                {
                                    const scalar x = (px ? x1 : x0);
                                    for (label py = 0; py < 2; py++)
                                    {
                                        const scalar y = (py ? y1 : y0);
                                        for (label pz = 0; pz < 2; pz++)
                                        {
                                            const scalar z = (pz ? z1 : z0);
                                            const scalar dDotN =
                                                ((point(x, y, z) - fc) & fA);
                                            minDotCell = min(minDotCell, dDotN);
                                        }
                                    }
                                }

                                if (minDotCell < 0)
                                {
                                    const label cellId =
                                        cx + processorFaceGridN
                                       *(cy + processorFaceGridN*cz);
                                    const DynamicList<label>& binFaces =
                                        remoteFaceBins[proci][cellId];
                                    forAll(binFaces, bI)
                                    {
                                        candidateFaces.append(binFaces[bI]);
                                    }
                                }
                            }
                        }
                    }
                }

                while
                (
                    candidateFacePos < candidateFaces.size()
                 && !batchFull
                )
                {
                    j = candidateFaces[candidateFacePos++];

                    if (proci != Pstream::myProcNo() || i != j)
                    {
                        const point&  remFc  = remoteFc[j];
                        const vector& remA   = remoteArea[j];
                        const label&  remAgg = remoteAgg[j];

                        const vector d = remFc - fc;

                        if (((d & fA) < 0.) && ((d & remA) > 0))
                        {
                            start.append(fc + 0.001*d);
                            startIndex.append(i);
                            startAgg.append
                            (
                                globalNumbering.toGlobal(proci, fAgg)
                            );
                            end.append(fc + 0.999*d);
                            endIndex.append
                            (
                                globalNumbering.toGlobal(proci, j)
                            );
                            endAgg.append
                            (
                                globalNumbering.toGlobal(proci, remAgg)
                            );

                            if (startIndex.size() >= maxDynListLength)
                            {
                                batchFull = true;
                            }
                        }
                    }
                }

                if (!batchFull && candidateFacePos >= candidateFaces.size())
                {
                    j = 0;
                    i++;
                    candidateFacei = -1;
                }
            }
            else
            {
                while (j < remoteFc.size() && !batchFull)
                {
                    if (proci != Pstream::myProcNo() || i != j)
                    {
                        const point&  remFc  = remoteFc[j];
                        const vector& remA   = remoteArea[j];
                        const label&  remAgg = remoteAgg[j];

                        const vector d = remFc - fc;

                        if (((d & fA) < 0.) && ((d & remA) > 0))
                        {
                            start.append(fc + 0.001*d);
                            startIndex.append(i);
                            startAgg.append
                            (
                                globalNumbering.toGlobal(proci, fAgg)
                            );
                            end.append(fc + 0.999*d);
                            endIndex.append
                            (
                                globalNumbering.toGlobal(proci, j)
                            );
                            endAgg.append
                            (
                                globalNumbering.toGlobal(proci, remAgg)
                            );

                            if (startIndex.size() >= maxDynListLength)
                            {
                                // Batch full: advance j so next batch
                                // resumes at j+1, same i.
                                batchFull = true;
                                j++;
                            }
                        }
                    }
                    if (!batchFull) j++;
                }

                if (!batchFull)
                {
                    j = 0;
                    i++;
                }
            }
        }

        procNotDone = (i < myFc.size());

        // --- Process this batch (collective across all procs) ---
        List<pointIndexHit> hitInfo(startIndex.size());
        surfacesMesh.findLine(start, end, hitInfo);

        // Return hit global agglo index
        labelList aggHitIndex;
        surfacesMesh.getField(hitInfo, aggHitIndex);

        DynamicList<label> dRayIs;

        // Collect the rays which don't have obstacle in between
        // rayStartFace and rayEndFace. Rays hitting themselves are
        // stored in dRayIs for further tracing.
        forAll(hitInfo, rayI)
        {
            if (!hitInfo[rayI].hit())
            {
                rayStartFace.append(startIndex[rayI]);
                rayEndFace.append(endIndex[rayI]);
            }
            else if (aggHitIndex[rayI] == startAgg[rayI])
            {
                dRayIs.append(rayI);
            }
        }

        start.clear();

        // Continue rays which hit themselves.
        bool firstLoop = true;
        DynamicField<point> startHitItself;
        DynamicField<point> endHitItself;
        label iter = 0;

        do
        {
            labelField rayIs;
            rayIs.transfer(dRayIs);
            dRayIs.clear();

            forAll(rayIs, rayI)
            {
                const label rayID = rayIs[rayI];
                label hitIndex = firstLoop ? rayIs[rayI] : rayI;

                if (hitInfo[hitIndex].hit())
                {
                    if (aggHitIndex[hitIndex] == startAgg[rayID])
                    {
                        const vector& endP   = end[rayID];
                        const vector& startP = hitInfo[hitIndex].hitPoint();
                        const vector  d      = endP - startP;

                        startHitItself.append(startP + 0.01*d);
                        endHitItself.append(startP + 1.01*d);

                        dRayIs.append(rayID);
                    }
                    else if (aggHitIndex[hitIndex] == endAgg[rayID])
                    {
                        rayStartFace.append(startIndex[rayID]);
                        rayEndFace.append(endIndex[rayID]);
                    }
                }
            }

            hitInfo.clear();
            hitInfo.resize(dRayIs.size());

            surfacesMesh.findLine(startHitItself, endHitItself, hitInfo);
            surfacesMesh.getField(hitInfo, aggHitIndex);

            endHitItself.clear();
            startHitItself.clear();
            firstLoop = false;
            iter++;

        } while (returnReduce(hitInfo.size(), orOp<bool>()) > 0 && iter < 10);

        // --- End of batch ---

    } while (returnReduce(procNotDone, orOp<bool>()));
}
