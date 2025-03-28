clear;

sigLoc = "D:\Errol's Quests\UBC_20250319_SNR\2025-03-18_15-53-54_CHAN01_sig.mat";
refLoc = "D:\Errol's Quests\UBC_20250319_SNR\2025-03-18_15-53-54_CHAN01_ref.mat";
metaLoc = "D:\Errol's Quests\UBC_20250319_SNR\2025-03-18_15-53-54_CHAN01_metadata.mat";

saveLoc = "D:\Errol's Quests\UBC_20250319_SNR\DataAnalysisCode\tempSave1";

b = postProcessingObject_SNR(sigLoc,refLoc,metaLoc);

b.setAnaRange([2,5000]);
b.setApoLen(95);
% b.plotTimeDomSig();
% b.showApoWindow();
b.applyApoWindow();
b.removeSplineBackground();
b.blackmanHarris();
b.zeroPad(200);
b.fftAndMertz();
b.freqAvg(10);
%b.plotAbsorbance()
%b.SNRFinding([7.86e12 7.875e12],[7.88e12 8e12],0.05,"newTechTest");
%b.SNRFinding([6.096e12 6.106e12],[6.28e12 6.4e12],0.05,newTechTest2);
%b.SNRFinding([1.868e12 1.878e12],[1.74e12 1.85e12],0.05);


