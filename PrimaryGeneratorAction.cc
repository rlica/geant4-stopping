#include "PrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4Exception.hh"
#include "Randomize.hh"
#include "G4GenericMessenger.hh"
#include <fstream>
#include <sstream>
#include "G4EmCalculator.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include <algorithm>

PrimaryGeneratorAction::PrimaryGeneratorAction() {
    fParticleGun = new G4ParticleGun(1);
    fTargetThickness = 1.0 * mm;
    DefineCommands();
    LoadCrossSectionData();
    LoadAngularDistribution();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fParticleGun;
    delete fMessenger;
}

void PrimaryGeneratorAction::DefineCommands() {
    fMessenger = new G4GenericMessenger(this, "/reaction/", "Transfer reaction kinematics control");
    fMessenger->DeclareProperty("targetZ", fTargetZ, "Target Atomic Number");
    fMessenger->DeclareProperty("targetA", fTargetA, "Target Mass Number");
    fMessenger->DeclareProperty("recoilZ", fRecoilZ, "Recoil Atomic Number");
    fMessenger->DeclareProperty("recoilA", fRecoilA, "Recoil Mass Number");
    fMessenger->DeclareProperty("ejectileZ", fEjectileZ, "Ejectile Atomic Number");
    fMessenger->DeclareProperty("ejectileA", fEjectileA, "Ejectile Mass Number");
    fMessenger->DeclareProperty("thetaMin", fThetaMin, "Min CM angle of ejectile (deg)");
    fMessenger->DeclareProperty("thetaMax", fThetaMax, "Max CM angle of ejectile (deg)");
    fMessenger->DeclareProperty("thinTarget", fIsThinTarget, "Set to true to neglect beam energy loss and interact at the surface.");
}

void PrimaryGeneratorAction::LoadCrossSectionData() {
    std::ifstream file("../cross_sections_18O_30Si_2n_Estar_5p5_6p5.txt");
    if (!file.is_open()) {
        G4cerr << "WARNING: Failed to open ../cross_sections_18O_30Si_2n_Estar_5p5_6p5.txt" << G4endl;
        return;
    }
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        G4double E, tot_mb, est_mb;
        std::string tot_str, est_str;
        if (iss >> E >> tot_str >> est_str) {
            size_t pos = tot_str.find('D');
            if (pos != std::string::npos) tot_str[pos] = 'e';
            tot_mb = std::stod(tot_str);
            fEnergyData.push_back(E * MeV);
            fCrossSectionData.push_back(tot_mb);
        }
    }
}

G4double PrimaryGeneratorAction::GetCrossSection(G4double energy) {
    if (fEnergyData.empty()) return 1.0; // fallback uniform probability if file missing
    if (energy <= fEnergyData.front()) return fCrossSectionData.front();
    if (energy >= fEnergyData.back()) return fCrossSectionData.back();
    
    for (size_t i = 0; i < fEnergyData.size() - 1; ++i) {
        if (energy >= fEnergyData[i] && energy <= fEnergyData[i+1]) {
            G4double t = (energy - fEnergyData[i]) / (fEnergyData[i+1] - fEnergyData[i]);
            return fCrossSectionData[i] + t * (fCrossSectionData[i+1] - fCrossSectionData[i]);
        }
    }
    return 0.0;
}

void PrimaryGeneratorAction::LoadAngularDistribution() {
    std::ifstream file("../angular_distribution.txt");
    if (!file.is_open()) {
        G4cerr << "WARNING: Failed to open ../angular_distribution.txt. Falling back to isotropic (uniform in cos(theta))." << G4endl;
        return;
    }

    std::vector<G4double> thetas;
    std::vector<G4double> dSigma_dOmega;

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        G4double theta_deg, ds_do;
        if (iss >> theta_deg >> ds_do) {
            thetas.push_back(theta_deg * deg);
            dSigma_dOmega.push_back(ds_do);
        }
    }

    if (thetas.size() < 2) {
        G4cerr << "WARNING: Not enough points in angular_distribution.txt" << G4endl;
        return;
    }

    // Build the CDF
    // P(theta) ~ int_0^theta (dSigma/dOmega) * sin(theta') d(theta')
    fThetaData = thetas;
    fAngularCDF.push_back(0.0);
    G4double integral = 0.0;
    
    for (size_t i = 0; i < thetas.size() - 1; ++i) {
        G4double d_theta = thetas[i+1] - thetas[i];
        G4double theta_mid = 0.5 * (thetas[i+1] + thetas[i]);
        G4double ds_do_mid = 0.5 * (dSigma_dOmega[i+1] + dSigma_dOmega[i]);
        
        // Element of solid angle is ~ sin(theta) d(theta)
        integral += ds_do_mid * std::sin(theta_mid) * d_theta;
        fAngularCDF.push_back(integral);
    }

    // Normalize CDF to 1.0
    G4double total = fAngularCDF.back();
    if (total > 0.0) {
        for (auto& val : fAngularCDF) {
            val /= total;
        }
    }
}

G4double PrimaryGeneratorAction::SampleThetaCM() {
    // If no data loaded or invalid, fall back to true 3D isotropic sampling (uniform in cos(theta))
    if (fAngularCDF.empty() || fAngularCDF.back() == 0.0) {
        G4double cosTheta = -1.0 + 2.0 * G4UniformRand();
        return std::acos(cosTheta);
    }

    // Sample from CDF
    G4double u = G4UniformRand();
    auto it = std::lower_bound(fAngularCDF.begin(), fAngularCDF.end(), u);
    
    if (it == fAngularCDF.begin()) return fThetaData.front();
    if (it == fAngularCDF.end()) return fThetaData.back();
    
    size_t idx = std::distance(fAngularCDF.begin(), it);
    
    // Linear interpolation
    G4double cdf1 = fAngularCDF[idx - 1];
    G4double cdf2 = fAngularCDF[idx];
    G4double t1 = fThetaData[idx - 1];
    G4double t2 = fThetaData[idx];
    
    G4double t = (u - cdf1) / (cdf2 - cdf1);
    return t1 + t * (t2 - t1);
}

void PrimaryGeneratorAction::InitializeSlowingDown() {
    fIsInitialized = true;
    G4EmCalculator emCalc;
    G4ParticleDefinition* beamDef = fParticleGun->GetParticleDefinition();
    G4double E0 = fParticleGun->GetParticleEnergy();
    
    G4LogicalVolume* targetLV = G4LogicalVolumeStore::GetInstance()->GetVolume("Target", false);
    if (!targetLV || !targetLV->GetMaterial()) return;
    G4Material* mat = targetLV->GetMaterial();
    
    G4double R0 = emCalc.GetRange(E0, beamDef, mat);
    G4double R22 = emCalc.GetRange(22.0*MeV, beamDef, mat);
    
    fMaxDepth = R0 - R22;
    if (fMaxDepth < 0) fMaxDepth = 0.0;
    if (fMaxDepth > fTargetThickness) fMaxDepth = fTargetThickness;
    
    G4int nBins = 1000;
    fMaxSigma = 0.0;
    
    // Instead of using GetKinEnergy which is inconsistent for generic ions,
    // we build the table by varying Energy from E0 down to 22 MeV,
    // compute Range, and then invert it to get Energy vs z.
    
    std::vector<G4double> e_table;
    std::vector<G4double> z_table;
    
    G4double dE = (E0 - 22.0*MeV) / nBins;
    for (int i = 0; i <= nBins; ++i) {
        G4double E = E0 - i * dE;
        G4double R = emCalc.GetRange(E, beamDef, mat);
        G4double z = R0 - R;
        e_table.push_back(E);
        z_table.push_back(z);
    }
    
    fMaxDepth = z_table.back();
    if (fMaxDepth < 0) fMaxDepth = 0.0;
    if (fMaxDepth > fTargetThickness) fMaxDepth = fTargetThickness;
    
    G4double dz = fMaxDepth / nBins;
    for (int i = 0; i <= nBins; ++i) {
        G4double z = i * dz;
        
        // Find E by interpolating z_table
        G4double E_z = 22.0*MeV; // default
        if (z <= 0) E_z = E0;
        else if (z >= fMaxDepth) E_z = 22.0*MeV;
        else {
            for (size_t j = 0; j < z_table.size() - 1; ++j) {
                if (z >= z_table[j] && z <= z_table[j+1]) {
                    G4double t = (z - z_table[j]) / (z_table[j+1] - z_table[j]);
                    E_z = e_table[j] + t * (e_table[j+1] - e_table[j]);
                    break;
                }
            }
        }
        
        G4double sig = GetCrossSection(E_z);
        fDepthTable.push_back(z);
        fEnergyAtDepth.push_back(E_z);
        fSigmaAtDepth.push_back(sig);
        if (sig > fMaxSigma) fMaxSigma = sig;
    }
    
    std::ofstream outFile("18O_energy_vs_depth.txt");
    if (outFile.is_open()) {
        outFile << "# Depth(um) Energy(MeV)\n";
        for (size_t i = 0; i < fDepthTable.size(); ++i) {
            outFile << fDepthTable[i]/um << " " << fEnergyAtDepth[i]/MeV << "\n";
        }
        outFile.close();
    }
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
    // If not set by macro, set default particle and energy on first event
    if (fParticleGun->GetParticleDefinition() == nullptr || 
        fParticleGun->GetParticleDefinition()->GetParticleName() == "geantino") {
        G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(8, 18, 0.0);
        fParticleGun->SetParticleDefinition(ion);
        fParticleGun->SetParticleEnergy(26.0 * MeV);
    }
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    
    if (!fIsThinTarget && !fIsInitialized) {
        InitializeSlowingDown();
    }

    G4double zInteraction = 0.0;
    G4double E_interaction = fParticleGun->GetParticleEnergy();

    if (!fIsThinTarget) {
        // Rejection sampling for interaction depth z
        G4bool accepted = false;
        G4int maxTries = 10000;
        G4int tries = 0;
        
        while (!accepted && tries < maxTries) {
            G4double zCand = G4UniformRand() * fMaxDepth;
            G4double sig_cand = 0.0;
            
            G4int nBins = fDepthTable.size() - 1;
            G4int bin = -1;
            for (int i = 0; i < nBins; ++i) {
                if (zCand >= fDepthTable[i] && zCand < fDepthTable[i+1]) {
                    bin = i;
                    break;
                }
            }
            if (bin >= 0) {
                sig_cand = fSigmaAtDepth[bin];
                if (G4UniformRand() * fMaxSigma < sig_cand) {
                    zInteraction = zCand;
                    E_interaction = fEnergyAtDepth[bin];
                    accepted = true;
                }
            }
            tries++;
        }
        
        if (!accepted) {
            zInteraction = 0.0;
            E_interaction = fParticleGun->GetParticleEnergy();
        }
    }

    G4double targetFrontFace = -fTargetThickness / 2.0; 
    G4double zVertex = targetFrontFace + zInteraction;
    
    fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., zVertex));

    // Save beam state
    G4ParticleDefinition* beamDef = fParticleGun->GetParticleDefinition();
    G4double beamEnergy = fParticleGun->GetParticleEnergy();
    G4ThreeVector beamDir = fParticleGun->GetParticleMomentumDirection();

    // Masses
    G4double ma = beamDef->GetPDGMass();
    G4double mA = G4IonTable::GetIonTable()->GetIonMass(fTargetZ, fTargetA);
    G4double mB = G4IonTable::GetIonTable()->GetIonMass(fRecoilZ, fRecoilA);
    G4double mb = G4IonTable::GetIonTable()->GetIonMass(fEjectileZ, fEjectileA);

    // Lab Kinematics of Beam at vertex
    G4double Ta = E_interaction;
    G4double Ea = Ta + ma;
    G4double pA_mom = std::sqrt(Ta*Ta + 2.0*Ta*ma);

    // Center of Mass Energy
    G4double s = ma*ma + mA*mA + 2.0*Ea*mA;
    G4double Ecm = std::sqrt(s);

    // Adjust Q-value for excitation energy (E* = 6.0 MeV)
    G4double E_star = 6.0 * MeV;
    G4double mB_star = mB + E_star;

    // Recoil (B) CM Energy and Momentum
    G4double EB_cm = (s + mB_star*mB_star - mb*mb) / (2.0*Ecm);
    G4double pb_cm2 = EB_cm*EB_cm - mB_star*mB_star;
    
    if (pb_cm2 <= 0.0) {
        G4Exception("PrimaryGeneratorAction", "ERR_KIN", JustWarning, "Reaction is kinematically forbidden.");
        fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., zVertex));
        fParticleGun->GeneratePrimaryVertex(anEvent);
        return;
    }

    G4double pb_cm = std::sqrt(pb_cm2);

    // Sample Ejectile CM Angle using the angular distribution
    G4double theta_cm = SampleThetaCM();
    
    // Ensure we respect min/max limits if provided (though FRESCO usually covers 0-180 anyway)
    if (theta_cm < fThetaMin*deg) theta_cm = fThetaMin*deg;
    if (theta_cm > fThetaMax*deg) theta_cm = fThetaMax*deg;
    
    G4double phi = G4UniformRand() * 2.0 * CLHEP::pi;

    // Ejectile momentum in CM
    G4double pbx = pb_cm * std::sin(theta_cm) * std::cos(phi);
    G4double pby = pb_cm * std::sin(theta_cm) * std::sin(phi);
    G4double pbz = pb_cm * std::cos(theta_cm);

    // Recoil momentum in CM (opposite to ejectile)
    G4double pBx = -pbx;
    G4double pBy = -pby;
    G4double pBz = -pbz;

    // Lorentz Boost to Lab frame
    G4double Etot = Ea + mA;
    G4double beta_z = pA_mom / Etot;
    G4double gamma = Etot / Ecm;

    G4double pBx_lab = pBx;
    G4double pBy_lab = pBy;
    G4double pBz_lab = gamma * (pBz + beta_z * EB_cm); // Use actual CM energy of recoil
    G4double EB_lab = gamma * (EB_cm + beta_z * pBz);
    
    // The generated kinetic energy MUST be calculated relative to the excited mass
    // because the 6.0 MeV excitation is internal energy, not translational kinetic energy!
    G4double TB_lab = EB_lab - mB_star;

    G4ThreeVector dirB(pBx_lab, pBy_lab, pBz_lab);

    // Set particle gun to shoot recoil (32Si)
    fParticleGun->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(fRecoilZ, fRecoilA, 0.0));
    fParticleGun->SetParticleEnergy(TB_lab);
    fParticleGun->SetParticleMomentumDirection(dirB.unit());
    fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., zVertex));
    fParticleGun->GeneratePrimaryVertex(anEvent);

    // Restore beam state for the next event / UI query
    fParticleGun->SetParticleDefinition(beamDef);
    fParticleGun->SetParticleEnergy(beamEnergy);
    fParticleGun->SetParticleMomentumDirection(beamDir);
}
