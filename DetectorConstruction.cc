#include "DetectorConstruction.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Exception.hh"
#include "G4UserLimits.hh"
#include "G4GenericMessenger.hh"
#include "G4RunManager.hh"

DetectorConstruction::DetectorConstruction() {
    DefineCommands();
}

DetectorConstruction::~DetectorConstruction() {
    delete fMessenger;
}

void DetectorConstruction::DefineCommands() {
    fMessenger = new G4GenericMessenger(this, "/detector/", "Detector control");

    auto& matCmd = fMessenger->DeclareMethod("targetMaterial", &DetectorConstruction::SetTargetMaterial, "Set target material.");
    matCmd.SetParameterName("material", false);
    matCmd.SetDefaultValue("Target_30Si");

    fMessenger->DeclareProperty("customZ", fCustomZ, "Z of custom material");
    fMessenger->DeclareProperty("customA", fCustomA, "A of custom material");
    fMessenger->DeclarePropertyWithUnit("customDensity", "g/cm3", fCustomDensity, "Density of custom material");
    fMessenger->DeclareProperty("customSymbol", fCustomSymbol, "Chemical symbol of custom material");
}

void DetectorConstruction::SetTargetMaterial(const G4String& materialName) {
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* mat = nist->FindOrBuildMaterial(materialName);
    
    if (!mat) {
        G4cerr << "Material " << materialName << " not found!" << G4endl;
        return;
    }
    
    if (fTargetLV) {
        fTargetLV->SetMaterial(mat);
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
}



G4VPhysicalVolume* DetectorConstruction::Construct() {
    G4NistManager* nist = G4NistManager::Instance();
    if (!nist) {
        G4Exception("DetectorConstruction::Construct", "ERR01", FatalException, "NistManager not found.");
    }

    // Build target material dynamically using the custom properties set by macros
    G4String name = "Custom_Z" + std::to_string(fCustomZ) + "_A" + std::to_string(fCustomA);
    G4Isotope* iso = new G4Isotope(name + "_iso", fCustomZ, fCustomA, fCustomA * g/mole);
    G4Element* el = new G4Element(name + "_el", fCustomSymbol, 1);
    el->AddIsotope(iso, 100.0 * perCent);

    // If fCustomDensity is exactly 1.5 (internal units), we should multiply it by g/cm3 manually in case 
    // it was NOT set via messenger (i.e. default value was used)
    G4double finalDensity = fCustomDensity;
    if (finalDensity < 100.0) { 
        // fallback if it wasn't set by macro
        finalDensity = finalDensity * g/cm3;
    }
    
    G4Material* matTarget = new G4Material(name + "_mat", finalDensity, 1);
    matTarget->AddElement(el, 1);

    G4Material* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
    G4Box* solidWorld = new G4Box("World", 10*cm, 10*cm, 10*cm);
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, vacuum, "World");
    G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, true);

    G4Box* solidTarget = new G4Box("Target", 1*cm, 1*cm, 0.5*mm);
    fTargetLV = new G4LogicalVolume(solidTarget, matTarget, "Target");
    
    // Set step limits (e.g., 0.1 micrometer maximum step) to ensure accurate dE/dx integration
    fTargetLV->SetUserLimits(new G4UserLimits(0.1 * um));
    
    new G4PVPlacement(0, G4ThreeVector(0,0,0), fTargetLV, "Target", logicWorld, false, 0, true);

    if (!physWorld) {
        G4Exception("DetectorConstruction::Construct", "ERR02", FatalException, "Physical volume creation failed.");
    }

    return physWorld;
}
