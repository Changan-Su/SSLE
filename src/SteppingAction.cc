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
//
/// \file B1/src/SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SteppingAction.hh"
#include "G4OpticalPhoton.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Event.hh"
#pragma message ("!!! Using EventAction.hh from: " __FILE__)
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction) : fEventAction(eventAction) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void SteppingAction::UserSteppingAction(const G4Step* step)
// {
//   if (!fScoringVolume||!fPMTVolume) {
//     const auto detConstruction = static_cast<const DetectorConstruction*>(
//       G4RunManager::GetRunManager()->GetUserDetectorConstruction());
//       fScoringVolume = detConstruction->GetScoringVolume();
//       fPMTVolume = detConstruction->GetPMTVolume();  // 获取 PMT 逻辑体积
//   }

//   // get volume of the current step
//   G4LogicalVolume* volume =
//     step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
//      if (volume != fScoringVolume && volume != fPMTVolume) return;
//      const auto* secondaries = step->GetSecondaryInCurrentStep();// 获取当前步的次级粒子列表
//      int opticalPhotonCount = 0;// 统计光子数量
//      for (const auto& track : *secondaries) {
//        if (track->GetDefinition()->GetParticleName() == "opticalphoton") {
//          opticalPhotonCount++;
//        }
//      }
//      if (opticalPhotonCount > 0) {
//         if (volume == fScoringVolume){
//           fEventAction->AddPhotonInDetector(opticalPhotonCount);// 统计光子数量
//         }
//         else if(volume == fPMTVolume) {
//           fEventAction->AddPhotonInPMT(opticalPhotonCount);// 统计 PMT 中的光子数量
//         }
//      }



//     //   // collect energy deposited in this step
//     G4double edepStep = step->GetTotalEnergyDeposit();//本次步能量沉积
//     G4double stepLength = step->GetStepLength();//统计步长

//     //统计所有粒子能量沉积
//      if (edepStep > 0.) {
//         fEventAction->AddEdep(edepStep);
//         fEventAction->AddAbsorption(edepStep, step->GetStepLength());// 统计能量吸收和步长
//     }

//   }

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringVolume || !fPMTVolume) {
    const auto detConstruction = static_cast<const DetectorConstruction*>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detConstruction->GetScoringVolume();
    fPMTVolume    = detConstruction->GetPMTVolume();
  }

  G4Track* track = step->GetTrack();


  if (track->GetDefinition() == G4OpticalPhoton::Definition()) {
    // 只在刚生成的那一刻+1
    if (track->GetCurrentStepNumber() == 1) {
        fEventAction->AddPhotonGenerated(1);
    }
    G4LogicalVolume* volume =
        step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
    G4String volName = volume->GetName();
    if (volName.contains("PMT")) {
        fEventAction->AddPhotonInPMT(1);
        track->SetTrackStatus(fStopAndKill);  // 防止重复计数
    }
    return;
}

  G4double edepStep   = step->GetTotalEnergyDeposit();
  G4double stepLength = step->GetStepLength();

  if (edepStep > 0.) {
    fEventAction->AddEdep(edepStep);
    fEventAction->AddAbsorption(edepStep, stepLength);
    //  G4cout << "Edep: " << edepStep / keV << " keV in volume: " 
    //        << step->GetPreStepPoint()->GetPhysicalVolume()->GetName() << G4endl;
      const auto* secondaries = step->GetSecondaryInCurrentStep();
    int nOptical = 0;
    for (const auto& s : *secondaries) {
        if (s->GetDefinition() == G4OpticalPhoton::Definition()) {
            nOptical++;
        }
    }
    // G4cout << "[ScintDebug] Edep: " << edepStep / keV 
    //        << " keV  => OpticalPhotons: " << nOptical << G4endl;
  
  
          }
}  

}  // namespace B1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

