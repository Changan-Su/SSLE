#include "G4RunManagerFactory.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4UImanager.hh"

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "MyPhysicsList.hh"  

#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalParameters.hh"
#include "G4OpticalPhysics.hh"
using namespace B1;  // 确保使用你的命名空间


int main(int argc, char** argv)
{
    G4UIExecutive* ui = nullptr;
    if (argc == 1) {
        ui = new G4UIExecutive(argc, argv);
    }

    auto runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);

    runManager->SetUserInitialization(new DetectorConstruction());
    // runManager->SetUserInitialization(new MyPhysicsList()); 

        
    G4VModularPhysicsList* physicsList = new FTFP_BERT;
    physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());

    auto opticalPhysics = new G4OpticalPhysics();
    auto opticalParams = G4OpticalParameters::Instance();

    opticalParams->SetWLSTimeProfile("delta");

    opticalParams->SetScintTrackSecondariesFirst(true);

    opticalParams->SetCerenkovMaxPhotonsPerStep(100);
    opticalParams->SetCerenkovMaxBetaChange(10.0);
    opticalParams->SetCerenkovTrackSecondariesFirst(true);

    physicsList->RegisterPhysics(opticalPhysics);
    runManager->SetUserInitialization(physicsList);

    runManager->SetUserInitialization(new ActionInitialization());



    auto visManager = new G4VisExecutive();
    visManager->Initialize();

    auto UImanager = G4UImanager::GetUIpointer();

    if (!ui) {
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    } else {
        UImanager->ApplyCommand("/control/execute init_vis.mac");
        ui->SessionStart();
        delete ui;
    }

    delete visManager;
    delete runManager;
}
