//PAPER 1 ATLAS
    
#include <cmath>
#include <algorithm>
#include <iterator>
#include <iostream>
    
    cout << "##################################################\n";
    cout << "##                   New Event                  ##\n";
    cout << "##################################################\n\n";
    
    cout << "--------------------------------------------------\n";
    cout << "Electrons Tight:\n\n";
    for(int i=0; i < electronsTight.size(); i++){
        cout << "ElectronID: " << electronsTight[i] << "\n";
        cout << "Momentum: " << electronsTight[i]->PT << "\n";
        cout << "ETA: " << (electronsTight[i]->P4()).Eta() << "\n\n" ;
    }
    
    cout << "--------------------------------------------------\n";
    cout << "Muons:\n\n";
    for(int i=0; i < muonsCombined.size(); i++){
        cout << "MuonID: " << muonsCombined[i] << "\n";
        cout << "Momentum: " << muonsCombined[i]->PT << "\n";
        cout << "ETA: " << (muonsCombined[i]->P4()).Eta() << "\n\n";
    }
    
    cout << "--------------------------------------------------\n";
    cout << "Jets:\n\n";
    for(int i=0; i < jets.size(); i++){
        cout << "JetID: " << jets[i] << "\n";
        cout << "Momentum: " << jets[i]->PT << "\n";
        cout << "ETA: " << (jets[i]->P4()).Eta() << "\n\n";
    }
    
    cout << "--------------------------------------------------\n";
    cout << "Photons\n\n";
    for(int i=0; i < photons.size(); i++){
        cout << "PhotonID: " << photons[i] << "\n";
        cout << "Momentum: " << photons[i]->PT << "\n";
        cout << "ETA: " << (photons[i]->P4()).Eta() << "\n\n";
    }
    
    //////////////////////////
    // EVENT RECONSTRUCTION //
    //////////////////////////
    
    electronsTight = filterPhaseSpace(electronsTight, 25, -2.47, 2.47, true);
    std::vector<Electron*> isoElectrons = filterIsolation(electronsTight); /* as defined in AnalysisManager */
    
    cout << "AFTER ISOLATION\n";
    
    cout << "--------------------------------------------------\n";
    cout << "Electrons Tight\n\n";
    for(int i=0; i < isoElectrons.size(); i++){
        cout << "ElectronID: " << isoElectrons[i] << "\n";
        cout << "Momentum: " << isoElectrons[i]->PT << "\n";
        cout << "ETA: " << (isoElectrons[i]->P4()).Eta() << "\n\n" ;
    }
    
    muonsCombined = filterPhaseSpace(muonsCombined, 25, -2.5, 2.5);
    // Defining the isolation in AnalysisManager.
    std::vector<Muon*> isoMuons = filterIsolation(muonsCombined); /* as defined in AnalysisManager */
    
    cout << "Muons\n\n";
    for(int i=0; i < isoMuons.size(); i++){
        cout << "MuonID: " << isoMuons[i] << "\n";
        cout << "Momentum: " << isoMuons[i]->PT << "\n";
        cout << "ETA: " << (isoMuons[i]->P4()).Eta() << "\n\n";
    }
    
    photons = filterPhaseSpace(photons, 15, -2.37, 2.37, true);
    
    cout << "--------------------------------------------------\n";
    cout << "Photons\n\n";
    for(int i=0; i < photons.size(); i++){
        cout << "PhotonID: " << photons[i] << "\n";
        cout << "Momentum: " << photons[i]->PT << "\n";
        cout << "ETA: " << (photons[i]->P4()).Eta() << "\n\n";
    }
    
    jets = filterPhaseSpace(jets, 25, -2.5, 2.5);
    
    cout << "--------------------------------------------------\n";
    cout << "Jets\n\n";
    for(int i=0; i < jets.size(); i++){
        cout << "JetID: " << jets[i] << "\n";
        cout << "Momentum: " << jets[i]->PT << "\n";
        cout << "ETA: " << (jets[i]->P4()).Eta() << "\n\n";
    }
    
    

    cout << "--------------------------------------------------\n";
    cout << "BEFORE MUON-JET OVERLAP\n\n";
    if(isoMuons.size() != 0){
        for(int j=0; j < isoMuons.size(); j++){
            cout << "MuonID: " << isoMuons[j] << "\n\n";
            for(int i=0; i < jets.size(); i++){
                cout << "JetID: " << jets[i] << "\n";
                cout << "DR: " << ((jets[i]->P4()).DeltaR(isoMuons[j]->P4())) << "\n\n" ;
            }
        }
    }
    
    isoMuons = overlapRemoval(isoMuons, jets, 0.4);
    
    cout << "AFTER MUON-JET OVERLAP\n";
    for(int i=0; i < isoMuons.size(); i++){
        cout << "MuonID: " << isoMuons[i] << "\n";
        cout << "Momentum: " << isoMuons[i]->PT << "\n";
        cout << "ETA: " << (isoMuons[i]->P4()).Eta() << "\n\n";
    }
    
    
    
    cout << "--------------------------------------------------\n";
    cout << "BEFORE JET-ELECTRON OVERLAP\n\n";
    if(isoElectrons.size() != 0){
        for(int j=0; j < isoElectrons.size(); j++){
            cout << "ElectronID: " << isoElectrons[j] << "\n\n";
            for(int i=0; i < jets.size(); i++){
                cout << "JetID: "<< jets[i] << "\n";
                cout << "DR: " << ((jets[i]->P4()).DeltaR(isoElectrons[j]->P4())) << "\n\n" ;
            }
        }
    }
    
    jets = overlapRemoval(jets, isoElectrons, 0.2);
    
    cout << "AFTER JET-ELECTRON OVERLAP\n";
    for(int i=0; i < jets.size(); i++){
        cout << "JetID: " << jets[i] << "\n";
        cout << "Momentum: " << jets[i]->PT << "\n";
        cout << "ETA: " << (jets[i]->P4()).Eta() << "\n\n";
    }
    
    
    
    cout << "--------------------------------------------------\n";
    cout << "BEFORE JET-PHOTON OVERLAP\n\n";
    if(photons.size() != 0){
        for(int j=0; j < photons.size(); j++){
            cout << "PhotonID: " << photons[j] << "\n\n";
            for(int i=0; i < jets.size(); i++){
                cout << "JetID: " << jets[i] << "\n";
                cout << "DR: " << ((jets[i]->P4()).DeltaR(photons[j]->P4())) << "\n\n" ;
            }
        }
    }
    
    jets = overlapRemoval(jets, photons, 0.1);
    
    cout << "AFTER JET-PHOTON OVERLAP\n";
    cout << "Jets\n\n";
    for(int i=0; i < jets.size(); i++){
        cout << "JetID: " << jets[i] << "\n";
        cout << "Momentum: " << jets[i]->PT << "\n";
        cout << "ETA: " << (jets[i]->P4()).Eta() << "\n\n";
    }
    
    
    
    
    cout << "--------------------------------------------------\n";
    cout << "BEFORE ELECTRON-JET OVERLAP\n\n";
    if(isoElectrons.size() != 0){
        for(int j=0; j < isoElectrons.size(); j++){
            cout << "ElectronID: " << isoElectrons[j] << "\n\n";
            for(int i=0; i < jets.size(); i++){
                cout << "JetID: " << jets[i] << "\n";
                cout << "DR: " << ((jets[i]->P4()).DeltaR(isoElectrons[j]->P4())) << "\n\n" ;
            }
        }
    }
    
    isoElectrons = overlapRemoval(isoElectrons, jets, 0.4);
    
    cout << "AFTER ELECTRON-JET OVERLAP\n";
    cout << "Electrons Tight\n\n";
    for(int i=0; i < isoElectrons.size(); i++){
        cout << "ElectronID: " << isoElectrons[i] << "\n";
        cout << "Momentum: " << electronsTight[i]->PT << "\n";
        cout << "ETA: " << (electronsTight[i]->P4()).Eta() << "\n\n" ;
    }
    
    missingET->addMuons(isoMuons);
    
    /////////////////////
    // EVENT SELECTION //
    /////////////////////
    
    /* SINGLE-LEPTON TRIGGER REQUIREMENTS */
    
    cout << "--------------------------------------------------\n";
    cout << "SINGLE-LEPTON TRIGGER REQUIREMENTS\n\n";
    
    /* For electrons */
    
    cout << "Electrons before\n\n";
    for(int i=0; i < isoElectrons.size(); i++){
        cout << "ElectronID: " << isoElectrons[i] << "\n";
        cout << "Momentum: " << electronsTight[i]->PT << "\n";
        cout << "ETA: " << (electronsTight[i]->P4()).Eta() << "\n\n" ;
    }
    
    std::vector<Electron*> AllElectrons;
    for (int i = 0; i < electronsTight.size(); i++) {
        if (electronsTight[i]->PT < 60) {
            if ((std::find(isoElectrons.begin(), isoElectrons.end(), electronsTight[i]) != isoElectrons.end()) && (electronsTight[i]->PT > 24)) {
                AllElectrons.push_back(electronsTight[i]);
            }
        }
        else {
            AllElectrons.push_back(electronsTight[i]);
        }
    }
    
    cout << "Electrons after\n\n";
    for(int i=0; i < AllElectrons.size(); i++){
        cout << "ElectronID: " << AllElectrons[i] << "\n";
        cout << "Momentum: " << AllElectrons[i]->PT << "\n";
        cout << "ETA: " << (AllElectrons[i]->P4()).Eta() << "\n\n" ;
    }
    
    /* For muons */
    
    cout << "Muons before\n\n";
    for(int i=0; i < isoMuons.size(); i++){
        cout << "MuonID: " << isoMuons[i] << "\n";
        cout << "Momentum: " << isoMuons[i]->PT << "\n";
        cout << "ETA: " << (isoMuons[i]->P4()).Eta() << "\n\n";
    }
    
    std::vector<Muon*> AllMuons;
    for (int j = 0; j < muonsCombined.size(); j++) {
        if (muonsCombined[j]->PT < 34) {
            if ((std::find(isoMuons.begin(), isoMuons.end(), muonsCombined[j]) != isoMuons.end()) && (muonsCombined[j]->PT > 24)) {
                AllMuons.push_back(muonsCombined[j]);
            }
        }
        else {
            AllMuons.push_back(muonsCombined[j]);
        }
    }
    
    cout << "Muons after\n\n";
    for(int i=0; i < AllMuons.size(); i++){
        cout << "MuonID: " << AllMuons[i] << "\n";
        cout << "Momentum: " << AllMuons[i]->PT << "\n";
        cout << "ETA: " << (AllMuons[i]->P4()).Eta() << "\n\n";
    }
    
    /* EVENTS MUST HAVE EXACTLY 1 LEPTON */
    
    cout << "--------------------------------------------------\n";
    cout << "EXACTLY 1 LEPTON\n\n";
    
    cout << "Electrons\n\n";
    for(int i=0; i < AllElectrons.size(); i++){
        cout << "ElectronID: " << AllElectrons[i] << "\n";
        cout << "Momentum: " << AllElectrons[i]->PT << "\n";
        cout << "ETA: " << (AllElectrons[i]->P4()).Eta() << "\n\n" ;
    }
    
    cout << "Muons\n\n";
    for(int i=0; i < AllMuons.size(); i++){
        cout << "MuonID: " << AllMuons[i] << "\n";
        cout << "Momentum: " << AllMuons[i]->PT << "\n";
        cout << "ETA: " << (AllMuons[i]->P4()).Eta() << "\n\n";
    }
    
    if (AllElectrons.size() + AllMuons.size() != 1) {
        return;
    }
    
    /* EVENTS MUST HAVE AT LEAST 4 JETS */
    
    cout << "--------------------------------------------------\n";
    cout << "AT LEAST 4 LEPTONS\n\n";
    
    cout << "Jets\n\n";
    for(int i=0; i < jets.size(); i++){
        cout << "JetID: " << jets[i] << "\n";
        cout << "Momentum: " << jets[i]->PT << "\n";
        cout << "ETA: " << (jets[i]->P4()).Eta() << "\n\n";
    }
    
    if (jets.size() < 4) {
        return;
    }
    
    /* EVENTS MUST HAVE AT LEAST 1 B-JET */
    
    cout << "--------------------------------------------------\n";
    cout << "AT LEAST 1 B-JET\n\n";
    
    std::vector<Jet*> BJets;
    
    for (int i = 0; i < jets.size(); i++) {
        if (checkBTag(jets[i])){
            BJets.push_back(jets[i]);
        }
    }
    
    if (BJets.empty()) {
        return;
    }
    
    cout << "B-Jets\n\n";
    for(int i=0; i < BJets.size(); i++){
        cout << "B-JetID: " << BJets[i] << "\n";
        cout << "Momentum: " << BJets[i]->PT << "\n";
        cout << "ETA: " << (BJets[i]->P4()).Eta() << "\n\n";
    }
    
    
    /* ADDITIONAL REQUIREMENTS ON THE MISSING TRANSVERSE MOMENTUM
     AND TRANSVERSE MASS OF THE W BOSON CANDIDATE. */
    
    cout << "--------------------------------------------------\n";
    cout << "MISSING ET AND TRANSVERSE MASS OF W BOSON CANDIDATE\n\n";
    cout << "Missing ET: " << missingET->PT << "\n\n";
    
    /* For electrons */
    
    cout << "Electrons\n\n";
    for(int i=0; i < AllElectrons.size(); i++){
        cout << "ElectronID: " << AllElectrons[i] << "\n";
        cout << "Momentum: " << AllElectrons[i]->PT << "\n";
        cout << "ETA: " << (AllElectrons[i]->P4()).Eta() << "\n\n" ;
        cout << "MissingET and electron transverse mass: " << (missingET->P4()+AllElectrons[i]->P4()).Mt() << "\n\n";
    }
    
    if (AllElectrons.size() == 1) {
        if ((missingET->PT < 30) || ((missingET->P4()+AllElectrons[0]->P4()).Mt() < 30)) {
            return;
        }
    }
    
    /* For muons */
    
    cout << "Muons\n\n";
    for(int i=0; i < AllMuons.size(); i++){
        cout << "MuonID: " << AllMuons[i] << "\n";
        cout << "Momentum: " << AllMuons[i]->PT << "\n";
        cout << "ETA: " << (AllMuons[i]->P4()).Eta() << "\n\n";
        cout << "MissingET and muon transverse mass: " << (missingET->P4()+AllMuons[i]->P4()).Mt() << "\n\n";
    }
    
    if (AllMuons.size() == 1) {
        if ((missingET->PT < 20) || (missingET->PT + (missingET->P4()+AllMuons[0]->P4()).Mt() < 60)) {
            return;
        }
    }
    
    /* MORE EVENT REQUIREMENTS */
    
    /* We need exactly 1 photon */
    
    cout << "--------------------------------------------------\n";
    cout << "EXACTLY 1 PHOTON \n\n";
    for(int i=0; i < photons.size(); i++){
        cout << "PhotonID: " << photons[i] << "\n";
        cout << "Momentum: " << photons[i]->PT << "\n";
        cout << "ETA: " << (photons[i]->P4()).Eta() << "\n\n";
    }
    
    if (photons.size() != 1) {
        return;
    }
    
    /* Events with a jet within a cone of dR = 0.5 around the selected photon are rejected to remove photon radiation from quarks. */
    
    cout << "--------------------------------------------------\n";
    cout << "JET AROUND PHOTON\n\n";
    
    for(int j=0; j < photons.size(); j++){
        cout << "PhotonID: " << photons[j] << "\n\n";
        for(int i=0; i < jets.size(); i++){
            cout << "JetID: " << jets[i] << "\n";
            cout << "DR: " << ((jets[i]->P4()).DeltaR(photons[j]->P4())) << "\n\n" ;
        }
    }
    
    
    for(int j = 0; j < jets.size(); j++){
        if (abs(jets[j]->P4().DeltaR(photons[0]->P4())) <= 0.5){
            return;
        }
    }
    
    /* Photon near lepton */
    
    cout << "--------------------------------------------------\n";
    cout << "PHOTON AROUND LEPTON\n\n";
    
    for(int j=0; j < AllElectrons.size(); j++){
        cout << "ElectronID: " << AllElectrons[j] << "\n";
        for(int i=0; i < photons.size(); i++){
            cout << "PhotonID: " << photons[i] << "\n";
            cout << "DR: " << ((photons[i]->P4()).DeltaR(AllElectrons[j]->P4())) << "\n\n" ;
        }
    }
    
    if((AllElectrons.size() == 1) && (photons[0]->P4().DeltaR(AllElectrons[0]->P4())) <= 0.7){
        return;
    }
    
    for(int j=0; j < AllMuons.size(); j++){
        cout << "MuonID: " << AllMuons[j] << "\n";
        for(int i=0; i < photons.size(); i++){
            cout << "PhotonID: " << photons[i] << "\n";
            cout << "DR: "<< ((photons[i]->P4()).DeltaR(AllMuons[j]->P4())) << "\n\n" ;
        }
    }
    
    if((AllMuons.size() == 1) && (photons[0]->P4().DeltaR(AllMuons[0]->P4())) <= 0.7){
        return;
    }
    
    /* Invariant mass of the electron and the photon has to be outside a 5GeV mass window around the Z boson mass. */
    
    cout << "--------------------------------------------------\n";
    cout << "ELECTRON + PHOTON IM \n\n";
    
    for(int j=0; j < AllElectrons.size(); j++){
        cout << "ElectronID: " << AllElectrons[j] << "\n\n";
        for(int i=0; i < photons.size(); i++){
            cout << "PhotonID: " << photons[i] << "\n";
            cout << "Invariant Mass: " << (photons[i]->P4() + AllElectrons[j]->P4()).M() << "\n\n" ;
        }
    }
    
    if (AllElectrons.size() == 1) {
        int mass = (photons[0]->P4() + AllElectrons[0]->P4()).M();
        if ((mass > 86) && (mass < 96)) {
            return;
        }
    }
    
    ///////////////////////////////////
    // SIGNAL REGIONS CATEGORIZATION //
    ///////////////////////////////////
    
    cout << "--------------------------------------------------\n";
    cout << "Event passed all requirements\n";
    cout << "--------------------------------------------------\n\n";
    
    countSignalEvent("tt+photon");
    