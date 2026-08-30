#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class G4GenericMessenger;
#include <vector>

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* anEvent) override;

    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

private:
    void DefineCommands();

    G4ParticleGun* fParticleGun;
    G4GenericMessenger* fMessenger = nullptr;

    G4int fTargetZ = 14;
    G4int fTargetA = 30;
    G4int fRecoilZ = 14;
    G4int fRecoilA = 32;
    G4int fEjectileZ = 8;
    G4int fEjectileA = 16;
    
    G4double fThetaMin = 0.0;
    G4double fThetaMax = 180.0; // degrees
    
    G4bool fIsThinTarget = false;

    // Cross-section data
    void LoadCrossSectionData();
    G4double GetCrossSection(G4double energy);
    std::vector<G4double> fEnergyData;
    std::vector<G4double> fCrossSectionData;

    // Angular distribution data
    void LoadAngularDistribution();
    G4double SampleThetaCM();
    std::vector<G4double> fThetaData;
    std::vector<G4double> fAngularCDF;

    // Fast lookup for slowing down
    G4bool fIsInitialized = false;
    void InitializeSlowingDown();
    std::vector<G4double> fDepthTable;
    std::vector<G4double> fEnergyAtDepth;
    std::vector<G4double> fSigmaAtDepth;
    G4double fMaxSigma = 0.0;
    G4double fMaxDepth = 0.0;
    G4double fTargetThickness; // initialized in .cc
};

#endif
