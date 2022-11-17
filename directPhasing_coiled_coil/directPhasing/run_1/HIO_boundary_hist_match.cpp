/*
icpc *.cpp -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -O
*/

#include "../crysPara.hpp"

int main() {

    //compute resolution of each unique structure factor, output uniqStruFactReso
    computeUniqStruFactReso();
    
    /***************************************************************************

                compute true phases for all allowed origin translation

    ***************************************************************************/
    //read in Fmodel got from Phenix (and CCP4) with bulk solvent correction, output struFactAmplCal and struFactPhase
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << "../" << PDB_CODE << "_fmodel.hkl,.HKL,.sca";
    inputFmodel(fileNameStream.str());

    //generate true phases for all allowed origin trans, input uniqStruFactAmplCal and uniqStruFactPhase, output truePhaseOrigTran
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << PDB_CODE << "_true_phase_orig_tran_"; //.txt
    generateTruePhaseForAllOrigTran("not output true phases", fileNameStream.str());


    /***************************************************************************

                compute true prot mask for all allowed origin translation

    ***************************************************************************/

    //read in Fmodel got from Phenix (and CCP4) with bulk solvent correction, output struFactAmplCal and struFactPhase
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << "../" << PDB_CODE << "_fatom.hkl,.HKL,.sca";
    inputFmodel(fileNameStream.str());

    //generate true prot mask for all origin trans, input uniqStruFactAmplCal and uniqStruFactPhase, output trueProtMaskOrigTran
    generateTrueProtMaskForAllOrigTran();
    

    /***************************************************************************

                input observed structure factor amplitude

    ***************************************************************************/
    //read in observerd diffraction data in experiment, reserved in uniqStruFactStat, uniqStruFactAmplObs, uniqStruFactAmplObsSigm
    //data format (3%5d, %2d, %9.1f, %7.1f, %8.3f) for h, k, l, obsStatus, obssf, obssfSigma, resolution
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << "../" << PDB_CODE << "_uniq_sf.txt";
    inputUniqStruFactAmplObs(fileNameStream.str());


    /***************************************************************************

                input observed reference histogram

    ***************************************************************************/
    
    //read in the boundary values of reference protein hist, reserved in stdBoundLargeBin
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << "../" << PDB_CODE_REFE_HIST << "_hist.txt";
    inputRefeHist(fileNameStream.str());


    /***************************************************************************

                initialize random density in asu

    ***************************************************************************/
    /*
    //read in density in asu from deposited pdb file, each atom is represented by density 1.0
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << "../" << PDB_CODE << ".pdb";
    unsigned long seed = mixClockTimeGetid();
    ofstream outputSeed("rand_seed_init_dens_asu.txt");
    outputSeed << seed << endl;    outputSeed.close();
    srand (seed); //initialize random seed
    inputDensAsuFromPdb(fileNameStream.str(), 0.05); //input random 2% atom positions
    */
    
    //generate random density in asu
    unsigned long seed = mixClockTimeGetid();
    ofstream outputSeed("rand_seed_init_dens_asu.txt");
    outputSeed << seed << endl;    outputSeed.close();
    generateRandDensAsu(seed);
    

    /***************************************************************************

                input SigmWeigObs

    ***************************************************************************/
    
    //read in the sigmweigobs
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << "../" << PDB_CODE << "_iter_sigm.txt";
    inputSigmWeigObsAndSigmWeigAvg(fileNameStream.str());
    
    //apply symmetric operation to get density in unit cell from densAsu, input densAsu, output densUnitCell
    applySymmOperDensAsu();
    
    //prepare for fft
    setHalfGridPhase(); //set half grid phase shift for fft, initialize halfGridPosi and halfGridNega

    //apply backward fft, input densUnitCell, output allStruFactAmpl and allStruFactPhase
    applyFFTDensToStruFact();

    //catch uniqStruFactAmplCal and uniqStruFactPhase from allStruFactAmpl and allStruFactPhase
    catchUniqStruFactAmplCalAndPhase();


    /***************************************************************************

                start iteration on density and phase operation

    ***************************************************************************/
    int iorig0 = 0; //initial origin choice

    //record the convergence status
    ofstream outputConvStat("log_conv_stat.txt");

    for(int iter = 0; iter <= numIter - 1; iter++) {
    
        //vary solvent content, input solvContIni, solvContFina, output new solvCont
        solvCont = varySolvContLinearly(iter);

        //vary sigmWeigAvg
        sigmWeigAvg = varySigmWeigAvgLinearly(iter);
        
        //vary sigmaWeigObs
        sigmWeigObs = varySigmWeigObsNonLinearly(iter);

        //vary Fobs, input uniqStruFactAmplOrigObs, output uniqStruFactAmplObs
        varyWeightUniqStruFactAmplObs(sigmWeigObs);
        
        //fill up missing reflections, input uniqStruFactAmplCal and uniqStruFactAmplObs, output missStruFactAmpl
        fillMissRefl(iter);

        //compute allStruFactAmpl from uniqStruFactAmplObs(and missStruFactAmpl)
        computeAllStruFactAmplFromUniqStruFactAmplObs();

        //apply forward fft, input allStruFactAmpl and allStruFactPhase, output densUnitCell
        applyFFTStruFactToDens();

        //get densAsuTemp from densUnitCell, input densUnitCell, output densAsuTemp
        catchDensAsuTempFromDensUnitCell();

        //compute weighted average density in unit cell, input densUnitCell, output weigAvgDensUnitCell
        computeWeigAvgDensUnitCell(sigmWeigAvg);

        //get weigAvgDensAsu from weigAvgDensUnitCell, input weigAvgDensUnitCell, output weigAvgDensAsu
        catchWeigAvgDensAsuFromWeigAvgDensUnitCell();

        //compute a densCutoff value on weigAvgDensAsu by bins, input weigAvgDensAsu, return densCutoff
        densCutoff = computeWeigAvgDensCutoffAsu(1.0-solvCont);

        //make a protein mask in asu, input weigAvgDensAsu, output protMaskAsu
        makeProtMaskAsu(densCutoff);
        
        //apply HIO in asu to update densAsu, input previous densAsu and current densAsuTemp, output new densAsu
        applyHIO(iter); //HIODensLimit

        //apply histogram matching on densAsu, input densAsu, output new densAsu
        applyHistMatch();

        //apply solvent flattening on densAsu during the last iterations, input densAsu, output new densAsu
        if(numIter-iter < numIterSolvFlat) applySolvFlat();

         //apply symmetric operation to get densUniCell from updated densAsu, input densAsu, output densUnitCell
        applySymmOperDensAsu();

        //apply backward fft, input densUnitCell, output allStruFactAmpl and allStruFactPhase
        applyFFTDensToStruFact();

        //catch uniqStruFactAmplCal and uniqStruFactPhase from allStruFactAmpl and allStruFactPhase
        catchUniqStruFactAmplCalAndPhase();

        //compute scaleFact between Fobs and Fcal for working data set
        scaleFact[iter] = computeScaleFact();

        //compute R value
        Rwork[iter] = computeRvalue(1); //obsStatus = 1
        Rfree[iter] = computeRvalue(2); //obsStatus = 2
        
        //compute mean phase error for all allowed origin translations in resolution sphere
        float meanPhaseErrorResoSphere[numResoShell][numOrigTran] = {};
        computeMeanPhaseErrorResoSphere(meanPhaseErrorResoSphere);
        for(int ishell = 0; ishell <= numResoShell-1; ++ishell)
            for(int iorig = 0; iorig <= numOrigTran-1; ++iorig)
                meanPhaseErrorOrigTranResoSphere[iter][ishell][iorig]= meanPhaseErrorResoSphere[ishell][iorig];

        //compute mean phase error for all allowed origin translations in resolution shell
        float meanPhaseErrorResoShell[numResoShell][numOrigTran] = {};
        computeMeanPhaseErrorResoShell(meanPhaseErrorResoShell);
        for(int ishell = 0; ishell <= numResoShell-1; ++ishell)
            for(int iorig = 0; iorig <= numOrigTran-1; ++iorig)
                meanPhaseErrorOrigTranResoShell[iter][ishell][iorig]= meanPhaseErrorResoShell[ishell][iorig];

        //compute correlation coefficient for all allowed origin translations
        float corrCoef[numResoShell][numOrigTran] = {};
        computeCorrCoefInSphere(corrCoef);
        for(int ishell = 0; ishell <= numResoShell-1; ++ishell)
            for(int iorig = 0; iorig <= numOrigTran-1; ++iorig)
                corrCoefOrigTran[iter][ishell][iorig]= corrCoef[ishell][iorig];
        
        //compute prot mask error for all allowed origin translations
        float protMaskMatch[numOrigTran] = {};
        computeProtMaskMatch(protMaskMatch);
        for(int iorig = 0; iorig <= numOrigTran-1; ++iorig)
            protMaskMatchOrigTran[iter][iorig]= protMaskMatch[iorig];

        //check origin choice
        iorig0 = checkOrigChoiceByProtMaskMatch(protMaskMatch);
        outputConvStat << setfill(' ') << setw(8) << fixed << iter << ',' << setfill(' ') << setw(5) << fixed << iorig0;
        //outputConvStat << ',' << setfill(' ') << setw(6) << setprecision(3) << fixed << protMaskMatch[iorig0] << endl;
        for(int iorig = 0; iorig <= numOrigTran-1; ++iorig)
            outputConvStat << ',' << setfill(' ') << setw(6) << setprecision(3) << fixed << protMaskMatch[iorig];
        outputConvStat << endl;
        
        //check convergence
        checkConvByMeanPhaseError(meanPhaseErrorResoSphere); //actually by meanPhaseErrorResoSphere[0][iorig]
       
        //output the progress
        cout << iter << " / " << numIter << endl;
    }
    
    outputConvStat.close();

    /***************************************************************************

                output calculated results

    ***************************************************************************/
    //output R vlaue
    fileNameString = "R_work.txt";
    outputRvalue(fileNameString, Rwork);

    fileNameString = "R_free.txt";
    outputRvalue(fileNameString, Rfree);
    
    //--------------------------------------------------------------------------
    //output mean phase error in resolution sphere
    fileNameString = "mean_phase_error_reso_sphere.txt";
    outputMeanPhaseErrorResoSphere(fileNameString, iorig0);

    //--------------------------------------------------------------------------
    //output mean phase error in resolution shell
    fileNameString = "mean_phase_error_reso_shell.txt";
    outputMeanPhaseErrorResoShell(fileNameString, iorig0);

    //--------------------------------------------------------------------------
    //output correlation coefficient
    fileNameString = "corr_coef.txt";
    outputCorrCoef(fileNameString, iorig0);
    
    //--------------------------------------------------------------------------
    //output protein mask match
    fileNameString = "prot_mask_match.txt";
    outputProtMaskMatch(fileNameString, iorig0);
    
    //--------------------------------------------------------------------------
    /*
    //output the final protMask in unit cell to a pdb file
    fileNameString = "prot_mask_fina.pdb";
    densCutoff = computeWeigAvgDensCutoffAsu(1.0-solvCont);
    outputWeigAvgDensUnitCell(fileNameString, densCutoff);

    //--------------------------------------------------------------------------
    //output the final density in unit cell to a pdb file
    fileNameString = "densUnitCell_fina.pdb";
    densCutoff = computeDensCutoffUnitCell((1.0-solvCont)*0.05);
    outputDensUnitCell(fileNameString, densCutoff);
    */
    
    //--------------------------------------------------------------------------
    //output uniqStruFactPhase and missStruFactAmpl
    fileNameStream.str( string() );
    fileNameStream.clear();
    fileNameStream << PDB_CODE << "_phase.cif";
    outputCifFileForCoot(fileNameStream.str());
    
    //--------------------------------------------------------------------------

    return 0;
}
