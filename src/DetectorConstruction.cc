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
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Trd.hh"

#include "G4Region.hh"
#include "G4ProductionCuts.hh"


namespace B1
{

  
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  G4VPhysicalVolume* DetectorConstruction::Construct()
  {
    // Get nist material manager
    G4NistManager* nist = G4NistManager::Instance();
  
    // Envelope parameters
    //
    G4double env_sizeXY = 20 * cm, env_sizeZ = 50 * cm;
    G4Material* env_mat = nist->FindOrBuildMaterial("G4_Galactic");
  
    // Option to switch on/off checking of volumes overlaps
    //
    G4bool checkOverlaps = true;
  
    //
    // World
    //
    G4double world_sizeXY = env_sizeXY;
    G4double world_sizeZ =  env_sizeZ;
    G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
  
    auto solidWorld =
      new G4Box("World",  // its name
                0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);  // its size
  
    auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
                                          world_mat,  // its material
                                          "World");  // its name
  
    auto physWorld = new G4PVPlacement(nullptr,  // no rotation
                                       G4ThreeVector(),  // at (0,0,0)
                                       logicWorld,  // its logical volume
                                       "World",  // its name
                                       nullptr,  // its mother  volume
                                       false,  // no boolean operation
                                       0,  // copy number
                                       checkOverlaps);  // overlaps checking
  
    //
    // Envelope
    //
    auto solidEnv = new G4Box("Envelope",  // its name
                              0.5 * env_sizeXY, 0.5 * env_sizeXY, 0.5 * env_sizeZ);  // its size
  
    auto logicEnv = new G4LogicalVolume(solidEnv,  // its solid
                                        env_mat,  // its material
                                        "Envelope");  // its name
  
    new G4PVPlacement(nullptr,  // no rotation
                      G4ThreeVector(),  // at (0,0,0)
                      logicEnv,  // its logical volume
                      "Envelope",  // its name
                      logicWorld,  // its mother  volume
                      false,  // no boolean operation
                      0,  // copy number
                      checkOverlaps);  // overlaps checking
  

    G4Element* elNa = G4NistManager::Instance()->FindOrBuildElement("Na");
    G4Element* elI  = G4NistManager::Instance()->FindOrBuildElement("I");
    G4double density = 3.67 * g/cm3;
    G4Material* NaI_Tl = new G4Material("NaI_Tl", density, 2);
    NaI_Tl->AddElement(elNa, 1);
    NaI_Tl->AddElement(elI, 1);

    // 材料属性表用LXe风格
    std::vector<G4double> nai_Energy = {2.07 * eV, 2.34 * eV, 2.62 * eV, 2.89 * eV, 3.10 * eV};
    std::vector<G4double> nai_SCINT = {1.0, 1.0, 1.0, 1.0, 1.0}; // 保持能量范围覆盖NaI发射区间
    std::vector<G4double> nai_RIND = {1.85, 1.85, 1.85, 1.85, 1.85};
    // std::vector<G4double> nai_ABSL = {38. * cm, 38. * cm, 38. * cm, 38. * cm, 38. * cm};
    std::vector<G4double> nai_ABSL = {100. * cm, 100. * cm, 100. * cm, 100. * cm, 100. * cm}; // 吸收长度设置为100cm

    auto nai_mt = new G4MaterialPropertiesTable();
    nai_mt->AddProperty("SCINTILLATIONCOMPONENT1", nai_Energy, nai_SCINT);
    nai_mt->AddProperty("SCINTILLATIONCOMPONENT2", nai_Energy, nai_SCINT); // 没slow就全0也行
    nai_mt->AddProperty("RINDEX", nai_Energy, nai_RIND);
    nai_mt->AddProperty("ABSLENGTH", nai_Energy, nai_ABSL);
    nai_mt->AddConstProperty("SCINTILLATIONYIELD", 38000. / MeV); // 或者先试12000看看能否产光
    nai_mt->AddConstProperty("RESOLUTIONSCALE", 1.0);
    nai_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 250. * ns);
    nai_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 0. * ns);
    nai_mt->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
    nai_mt->AddConstProperty("SCINTILLATIONYIELD2", 0.0);
    NaI_Tl->SetMaterialPropertiesTable(nai_mt);
    // G4Element* elXe = G4NistManager::Instance()->FindOrBuildElement("Xe");
    // G4double density = 2.96 * g/cm3;
    // G4Material* LXe = new G4Material("LXe", density, 1);
    // LXe->AddElement(elXe, 1);
    // G4MaterialPropertiesTable* LXe_mt = new G4MaterialPropertiesTable();
    // std::vector<G4double> lxe_Energy = {7.0 * eV, 7.07 * eV, 7.14 * eV};
    // std::vector<G4double> lxe_SCINT = {0.1, 1.0, 0.1};
    // std::vector<G4double> lxe_RIND = {1.59, 1.57, 1.54};
    // std::vector<G4double> lxe_ABSL = {35. * cm, 35. * cm, 35. * cm};
    // LXe_mt = new G4MaterialPropertiesTable();
    // LXe_mt->AddProperty("SCINTILLATIONCOMPONENT1", lxe_Energy, lxe_SCINT);
    // LXe_mt->AddProperty("SCINTILLATIONCOMPONENT2", lxe_Energy, lxe_SCINT);
    // LXe_mt->AddProperty("RINDEX", lxe_Energy, lxe_RIND);
    // LXe_mt->AddProperty("ABSLENGTH", lxe_Energy, lxe_ABSL);
    // LXe_mt->AddConstProperty("SCINTILLATIONYIELD", 12000. / MeV);
    // LXe_mt->AddConstProperty("RESOLUTIONSCALE", 1.0);
    // LXe_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 20. * ns);
    // LXe_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 45. * ns);
    // LXe_mt->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
    // LXe_mt->AddConstProperty("SCINTILLATIONYIELD2", 0.0);
    // LXe->SetMaterialPropertiesTable(LXe_mt);

    // G4Material* NaI_Tl = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");
    G4double crystal_rmin = 0. * cm, crystal_rmax = 5. * cm;
    G4double crystal_hz = 5 * cm;
    G4double crystal_phimin = 0. * deg, crystal_phimax = 360. * deg;
    auto solidCrystal = new G4Tubs("Crystal", crystal_rmin, crystal_rmax, crystal_hz/2, crystal_phimin, crystal_phimax);
    auto logicCrystal = new G4LogicalVolume(solidCrystal, NaI_Tl, "Crystal");
    // auto logicCrystal = new G4LogicalVolume(solidCrystal, LXe, "Crystal");
    G4ThreeVector posCrystal = G4ThreeVector(0, 0, 0);
    auto physCrystal = new G4PVPlacement(nullptr, posCrystal, logicCrystal, "Crystal", logicEnv, false, 0, checkOverlaps);
  //set the crystal as scoring volume

          // ================= 设置 NaI 区域的 Production Cut =================
    G4Region* regionNaI = new G4Region("NaIRegion");
    logicCrystal->SetRegion(regionNaI);
    regionNaI->AddRootLogicalVolume(logicCrystal);

    // G4ProductionCuts* cutsNaI = new G4ProductionCuts();
    // cutsNaI->SetProductionCut(10.0 *mm, G4ProductionCuts::GetIndex("gamma"));
    // cutsNaI->SetProductionCut(10.0 *mm, G4ProductionCuts::GetIndex("e-"));
    // cutsNaI->SetProductionCut(10.0 *mm, G4ProductionCuts::GetIndex("e+"));
    // regionNaI->SetProductionCuts(cutsNaI);
    // =============================================================



    fScoringVolume = logicCrystal;


      
    //PMT 
    // G4Material* pmt_mat = nist->FindOrBuildMaterial("G4_Galactic");
    // auto pmt_mpt = new G4MaterialPropertiesTable();
    // G4double energy[] = {2.0*eV, 3.5*eV};  // 350nm到600nm
    // G4double rindex_pmt[] = {1.0, 1.0};    // 真空 折射率 1.0
    // pmt_mpt->AddProperty("RINDEX", energy, rindex_pmt, 2);
    // pmt_mat->SetMaterialPropertiesTable(pmt_mpt);

    // --- 定义 PMT 材料 ---
    G4Material* pmt_mat = nist->FindOrBuildMaterial("G4_GLASS_PLATE");

    // --- 先定义属性数组 ---
    G4double energy[2] = {2.0 * eV, 3.5 * eV};          // Photon energy range
    G4double rindex_pmt[2] = {1.52, 1.52};              // PMT glass refractive index
    G4double efficiency[2] = {0.9, 0.9};              // Quantum efficiency 25%

    // --- 定义材质属性表 ---
    G4MaterialPropertiesTable* pmt_mpt = new G4MaterialPropertiesTable();
    pmt_mpt->AddProperty("RINDEX", energy, rindex_pmt, 2);
    pmt_mpt->AddProperty("EFFICIENCY", energy, efficiency, 2);

    // --- 材料挂属性表
    pmt_mat->SetMaterialPropertiesTable(pmt_mpt);

    G4double mptTheckness = 3 * mm ;
    G4Tubs* solidPMT = new G4Tubs("PMT", 0., crystal_rmax, mptTheckness/2, crystal_phimin, crystal_phimax);
    auto logicPMT = new G4LogicalVolume(solidPMT, pmt_mat, "PMT");
    fPMTVolume = logicPMT;  // Store PMT logical volume for optical photon tracking
    auto physPMT = new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, crystal_hz/2 + mptTheckness / 2),  // Position PMT above the crystal
                      logicPMT,  // its logical volume
                      "PMT",  // its name
                      logicEnv,  // its mother volume
                      false,  // no boolean operation
                      0,  // copy number
                      checkOverlaps);  // overlaps checking
    auto opticalSurface = new G4OpticalSurface("CrystalToPMTSurface");
    opticalSurface->SetType(dielectric_dielectric); // 介质-介质
    opticalSurface->SetModel(unified);
    // opticalSurface->SetFinish(polished); // 镜面反射，抛光面
    // opticalSurface->SetFinish(polishedfrontpainted);
    opticalSurface->SetFinish(groundfrontpainted);
    opticalSurface->SetSigmaAlpha(0.1); // 表面粗糙度
    
    fPhysCrystal = physCrystal;
    fPhysPMT = physPMT;

    new G4LogicalBorderSurface("CrystalToPMT",
    fPhysCrystal,    // 晶体的物理体
    fPhysPMT,        // PMT的物理体
    opticalSurface  // 光学表面
    );
    
 
// Reflector Solid
G4double reflectorThickness = 0.1 * mm;
auto solidOuter = new G4Tubs("Outer",
    0.,
    crystal_rmax + reflectorThickness,
    (crystal_hz/2 + reflectorThickness),
    0.*deg, 360.*deg);

auto solidInner = new G4Tubs("Inner",
    0.,
    crystal_rmax,
    (crystal_hz/2),
    0.*deg, 360.*deg);

auto solidReflector = new G4SubtractionSolid("Reflector", solidOuter, solidInner);

// Reflector Material
G4Material* PTFE = nist->FindOrBuildMaterial("G4_TEFLON");// 聚四氟乙烯 (PTFE) 材料
G4LogicalVolume* logicReflector = new G4LogicalVolume(solidReflector, PTFE, "Reflector");

// 放置 Reflector
auto physReflector = new G4PVPlacement(
    nullptr,
    G4ThreeVector(0, 0, 0), // 中心和晶体对齐
    logicReflector,
    "Reflector",
    logicEnv,
    false,
    0,
    checkOverlaps);

// Reflector Optical Surface
G4OpticalSurface* reflectorSurface = new G4OpticalSurface("ReflectorSurface");
reflectorSurface->SetType(dielectric_metal);
reflectorSurface->SetFinish(groundfrontpainted);
reflectorSurface->SetModel(unified);

// Reflector Surface Properties
G4MaterialPropertiesTable* mptReflector = new G4MaterialPropertiesTable();
G4double ephoton[2] = {2.0 * eV, 3.5 * eV};
G4double reflectivity[2] = {0.98, 0.98};
G4double r_efficiency[2] = {0.0, 0.0};

mptReflector->AddProperty("REFLECTIVITY", ephoton, reflectivity, 2);
mptReflector->AddProperty("EFFICIENCY", ephoton, r_efficiency, 2);
reflectorSurface->SetMaterialPropertiesTable(mptReflector);

new G4LogicalSkinSurface("ReflectorSkin", logicReflector, reflectorSurface);


    //
    // Shape 2
    //
    // G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
    // G4ThreeVector pos2 = G4ThreeVector(0, -1 * cm, 7 * cm);
  
    // // Trapezoid shape
    // G4double shape2_dxa = 12 * cm, shape2_dxb = 12 * cm;
    // G4double shape2_dya = 10 * cm, shape2_dyb = 16 * cm;
    // G4double shape2_dz = 6 * cm;
    // auto solidShape2 =
    //   new G4Trd("Shape2",  // its name
    //             0.5 * shape2_dxa, 0.5 * shape2_dxb, 0.5 * shape2_dya, 0.5 * shape2_dyb,
    //             0.5 * shape2_dz);  // its size
  
    // auto logicShape2 = new G4LogicalVolume(solidShape2,  // its solid
    //                                        shape2_mat,  // its material
    //                                        "Shape2");  // its name
  
    // new G4PVPlacement(nullptr,  // no rotation
    //                   pos2,  // at position
    //                   logicShape2,  // its logical volume
    //                   "Shape2",  // its name
    //                   logicEnv,  // its mother  volume
    //                   false,  // no boolean operation
    //                   0,  // copy number
    //                   checkOverlaps);  // overlaps checking
  
    // // Set Shape2 as scoring volume
    // //
    // fScoringVolume = logicShape2;
  
    //
    // always return the physical World
    //
    return physWorld;
  }
}