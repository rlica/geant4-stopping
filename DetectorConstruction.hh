#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class G4GenericMessenger;
class G4LogicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    DetectorConstruction();
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;
    void SetTargetMaterial(const G4String& materialName);

private:
    void DefineCommands();

    G4GenericMessenger* fMessenger = nullptr;
    G4LogicalVolume* fTargetLV = nullptr;
    
    G4int fCustomZ = 14;
    G4int fCustomA = 30;
    G4double fCustomDensity = 1.5;
    G4String fCustomSymbol = "Si";

    // Full target thickness along z, in G4 internal length units (1.0 = 1 mm).
    // Settable via /detector/targetThickness or /detector/targetArealDensity.
    G4double fTargetThickness = 1.0;
    // Areal density in internal units (set in mg/cm2 via macro). 0 = not set.
    // When > 0 it overrides fTargetThickness (converted using the target density).
    G4double fTargetArealDensity = 0.0;
};

#endif
