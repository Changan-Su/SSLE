//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file analysis/AnaEx01/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.



HistoManager::HistoManager()
{
  // Create or get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("csv");
  // the default file type can be overriden in run macro
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  if (!fFactoryOn) {
    analysisManager->SetVerboseLevel(1);

#ifdef G4MULTITHREADED
    analysisManager->SetNtupleMerging(true);
#endif

    analysisManager->SetHistoDirectoryName("histo");
    analysisManager->SetNtupleDirectoryName("ntuple");
  }

  G4bool fileOpen = analysisManager->OpenFile("AnaEx01");
  if (!fileOpen) {
    G4cerr << "\n---> HistoManager::Book(): cannot open " << analysisManager->GetFileName()
           << G4endl;
    return;
  }

  if (!fFactoryOn) {
    // Create histograms
    analysisManager->CreateH1("EAbs", "Edep in absorber (MeV)", 100, 0., 800 * MeV);
    analysisManager->CreateH1("EGap", "Edep in gap (MeV)", 100, 0., 100 * MeV);
    analysisManager->CreateH1("LAbs", "trackL in absorber (mm)", 100, 0., 1 * m);
    analysisManager->CreateH1("LGap", "trackL in gap (mm)", 100, 0., 50 * cm);
    analysisManager->CreateH1("PhotonCounts", "Number of optical photons in PMT", 200, 0, 2000);
    analysisManager->CreateH1("PhotonDepth","Depth index of captured photons", 20, -0.5, 19.5);
    // Create Ntuples
    analysisManager->CreateNtuple("Ntuple1", "Edep");
    analysisManager->CreateNtupleDColumn("Eabs");  // column Id = 0
    analysisManager->CreateNtupleDColumn("Egap");  // column Id = 1
    analysisManager->FinishNtuple();

    analysisManager->CreateNtuple("Ntuple2", "TrackL");
    analysisManager->CreateNtupleDColumn("Labs");  // column Id = 0
    analysisManager->CreateNtupleDColumn("Lgap");  // column Id = 1
    analysisManager->FinishNtuple();

    analysisManager->CreateNtuple("PhotonNtuple", "PhotonCounts");
    analysisManager->CreateNtupleIColumn("PhotonCounts");
    analysisManager->FinishNtuple();

    analysisManager->CreateNtuple("PhotonGenerated", "GeneratedCounts");
    analysisManager->CreateNtupleIColumn("GeneratedCounts");
    analysisManager->FinishNtuple();

    // 假设最深 20 层够用；如果想用实际层数，可在 DetectorConstruction 构造好后把层数传进来

    analysisManager->CreateNtuple("PhotonDepthNtuple", "Depth index of each captured photon");
    analysisManager->CreateNtupleIColumn("DepthIndex");
    analysisManager->FinishNtuple();

    fFactoryOn = true;
  }

  G4cout << "\n----> Output file is open in " << analysisManager->GetFileName() << "."
         << analysisManager->GetFileType() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......




void HistoManager::FillPhotonNtuple(G4int photonCounts)
{
  
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(2, 0, photonCounts); 
    analysisManager->AddNtupleRow(2);
}

void HistoManager::FillPhotonGeneratedNtuple(G4int photonGenerated)
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(3, 0, photonGenerated); // 新增ntuple编号为3
    analysisManager->AddNtupleRow(3);
}

void HistoManager::FillPhotonHisto(G4int photonCounts)
  {G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillH1(4, photonCounts); }

void HistoManager::Save()
{
  if (!fFactoryOn) {
    return;
  }

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();

  G4cout << "\n----> Histograms and ntuples are saved\n" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillHisto(G4int ih, G4double xbin, G4double weight)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillH1(ih, xbin, weight);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Normalize(G4int ih, G4double fac)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  auto h1 = analysisManager->GetH1(ih);
  if (h1 != nullptr) {
    h1->scale(fac);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void HistoManager::FillNtuple(G4double energyAbs, G4double energyGap, G4double trackLAbs,G4double trackLGap)
// {
//   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
//   // Fill 1st ntuple ( id = 0)
//   analysisManager->FillNtupleDColumn(0, 0, energyAbs);
//   analysisManager->FillNtupleDColumn(0, 1, energyGap);
//   analysisManager->AddNtupleRow(0);
//   // Fill 2nd ntuple ( id = 1)
//   analysisManager->FillNtupleDColumn(1, 0, trackLAbs);
//   analysisManager->FillNtupleDColumn(1, 1, trackLGap);
//   analysisManager->AddNtupleRow(1);
// }



void HistoManager::FillNtuple(G4double energyAbs, G4double trackLAbs)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Fill energy (Ntuple1)
  analysisManager->FillNtupleDColumn(0, 0, energyAbs);
  analysisManager->FillNtupleDColumn(0, 1, 0.);  // 默认 gap 为 0
  analysisManager->AddNtupleRow(0);

  // Fill track length (Ntuple2)
  analysisManager->FillNtupleDColumn(1, 0, trackLAbs);
  analysisManager->FillNtupleDColumn(1, 1, 0.);  // 默认 gap 为 0
  analysisManager->AddNtupleRow(1);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void HistoManager::PrintStatistic()
{

  if (!fFactoryOn) {
    return;
  }

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  G4cout << "\n ----> print histograms statistic \n" << G4endl;
  for (G4int i = 0; i < analysisManager->GetNofH1s(); ++i) {
    G4String name = analysisManager->GetH1Name(i);
    auto h1 = analysisManager->GetH1(i);

    G4String unitCategory;
    if (name[0U] == 'E') {
      unitCategory = "Energy";
    }
    if (name[0U] == 'L') {
      unitCategory = "Length";
    }
    // we use an explicit unsigned int type for operator [] argument
    // to avoid problems with windows compiler

    G4cout << name << ": mean = " << G4BestUnit(h1->mean(), unitCategory)
           << " rms = " << G4BestUnit(h1->rms(), unitCategory) << G4endl;
  }
}

void HistoManager::FillPhotonDepth(G4int depthIdx)
{
    G4AnalysisManager::Instance()->FillH1(5, depthIdx);   // 5 是上面那条 H1 的序号
    G4AnalysisManager::Instance()->FillNtupleIColumn(4,0,depthIdx); // 4 是新 ntuple
    G4AnalysisManager::Instance()->AddNtupleRow(4);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
