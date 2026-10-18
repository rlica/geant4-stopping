#include "RunAction.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
#include <cmath>

#include "PrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4IonTable.hh"
#include <ctime>
#include <fstream>
#include <sstream>

RunAction::RunAction() : G4UserRunAction() {
    G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
    accumulableManager->RegisterAccumulable(fSumTime);
    accumulableManager->RegisterAccumulable(fSumTime2);
    accumulableManager->RegisterAccumulable(fSumRecoilEnergy);
    accumulableManager->RegisterAccumulable(fSumInteractionDepth);
    accumulableManager->RegisterAccumulable(fSumStoppingDepth);
    accumulableManager->RegisterAccumulable(fStopCount);
    accumulableManager->RegisterAccumulable(fTransmitCount);
    accumulableManager->RegisterAccumulable(fBackCount);
    accumulableManager->RegisterAccumulable(fSumTInitE);
    accumulableManager->RegisterAccumulable(fSumTInitE2);
    accumulableManager->RegisterAccumulable(fSumTExitE);
    accumulableManager->RegisterAccumulable(fSumTExitE2);
    accumulableManager->RegisterAccumulable(fSumTInitBeta);
    accumulableManager->RegisterAccumulable(fSumTExitBeta);

    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("csv");
    analysisManager->SetVerboseLevel(0);
    analysisManager->SetFileName("stopping_data");
    analysisManager->CreateNtuple("RecoilData", "Recoil Stopping Data");
    analysisManager->CreateNtupleDColumn("StoppingTime_ps");
    analysisManager->CreateNtupleDColumn("RecoilEnergy_MeV");
    analysisManager->CreateNtupleDColumn("RecoilInitBeta");
    analysisManager->CreateNtupleDColumn("ExitEnergy_MeV");
    analysisManager->CreateNtupleDColumn("ExitBeta");
    analysisManager->FinishNtuple();
}

void RunAction::BeginOfRunAction(const G4Run*) {
    G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
    accumulableManager->Reset();

    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->OpenFile();
}

void RunAction::AddStoppingData(G4double time, G4double initialEnergy, G4double interactionDepth, G4double stoppingDepth) {
    fSumTime += time;
    fSumTime2 += time * time;
    fSumRecoilEnergy += initialEnergy;
    fSumInteractionDepth += interactionDepth;
    fSumStoppingDepth += stoppingDepth;
    fStopCount += 1;
}

void RunAction::AddExitData(G4double initEnergy, G4double initBeta, G4double exitEnergy, G4double exitBeta, G4bool forward) {
    if (forward) {
        fTransmitCount += 1;
        fSumTInitE += initEnergy;
        fSumTInitE2 += initEnergy * initEnergy;
        fSumTExitE += exitEnergy;
        fSumTExitE2 += exitEnergy * exitEnergy;
        fSumTInitBeta += initBeta;
        fSumTExitBeta += exitBeta;
    } else {
        // Backward exits are recorded in the ntuple (negative energies) but do
        // not enter the forward-transmitted averages.
        fBackCount += 1;
    }
}

void RunAction::EndOfRunAction(const G4Run* run) {
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();

    G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
    accumulableManager->Merge();

    const G4int count = fStopCount.GetValue();
    const G4int nForward = fTransmitCount.GetValue();
    const G4int nBack = fBackCount.GetValue();

    if ((count + nForward + nBack) > 0 && IsMaster()) {
        G4String particleName = "Unknown";
        G4double particleEnergy = 0.0;
        const PrimaryGeneratorAction* pga = dynamic_cast<const PrimaryGeneratorAction*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
        if (pga) {
            const G4ParticleGun* gun = pga->GetParticleGun();
            if (gun && gun->GetParticleDefinition()) {
                particleName = gun->GetParticleDefinition()->GetParticleName();
                particleEnergy = gun->GetParticleEnergy();
            }
        }

        G4String materialName = "Unknown";
        G4double density = 0.0;
        G4double thickness = 0.0; // full layer thickness along z
        G4LogicalVolume* targetLV = G4LogicalVolumeStore::GetInstance()->GetVolume("Target", false);
        if (targetLV && targetLV->GetMaterial()) {
            materialName = targetLV->GetMaterial()->GetName();
            density = targetLV->GetMaterial()->GetDensity();
            if (auto* solidTarget = dynamic_cast<G4Box*>(targetLV->GetSolid())) {
                thickness = 2.0 * solidTarget->GetZHalfLength();
            }
        }
        const G4bool thinTarget = pga && pga->IsThinTarget();

        // Reaction target as defined in the kinematics (targetZ/targetA).
        // In thin-target mode it is not simulated as a material; in
        // thick-target mode the simulated layer is the target itself.
        G4String reactionTargetName = "Unknown";
        if (pga) {
            auto* targetIon = G4IonTable::GetIonTable()->GetIon(pga->GetTargetZ(), pga->GetTargetA(), 0.0);
            if (targetIon) reactionTargetName = targetIon->GetParticleName();
        }

        std::ostringstream oss;
        const std::time_t now = std::time(nullptr);
        oss << "\n=============================================\n";
        oss << "               Run Summary                   \n";
        oss << "=============================================\n";
        oss << " Run finished:           " << std::ctime(&now); // ctime() appends '\n'
        oss << " Primary Particle:       " << particleName << "\n";
        oss << " Beam Energy:            " << particleEnergy / MeV << " MeV\n";
        oss << " Reaction Target:        " << reactionTargetName
            << " (Z=" << (pga ? pga->GetTargetZ() : 0)
            << ", A=" << (pga ? pga->GetTargetA() : 0) << ")";
        if (pga) {
            oss << (thinTarget
                ? " - not simulated (reaction at the target surface)"
                : " - simulated (also the stopping medium)");
        }
        oss << "\n";
        oss << " Stopping Mode:          " << (thinTarget
            ? "thin target (simulated layer acts as stopper/backing)"
            : "thick target (target is also the stopping medium)") << "\n";
        oss << " Stopper Material:       " << materialName << "\n";
        if (density > 0) {
            oss << " Stopper Density:        " << density / (g/cm3) << " g/cm3\n";
        }
        if (thickness > 0) {
            const G4double thicknessUm = thickness / um;
            oss << " Stopper Thickness:      " << thicknessUm << " um\n";
            if (density > 0) {
                // areal density [mg/cm2] = rho [g/cm3] * t [um] / 10
                oss << " Stopper Areal Density:  "
                       << (density / (g/cm3)) * thicknessUm / 10.0 << " mg/cm2\n";
            }
        }
        oss << " Total number of events: " << run->GetNumberOfEvent() << "\n";

        if (count > 0) {
            G4double sum = fSumTime.GetValue();
            G4double sum2 = fSumTime2.GetValue();

            G4double mean = sum / count;
            G4double variance = (sum2 / count) - (mean * mean);
            if (variance < 0) variance = 0;
            G4double stdDev = std::sqrt(variance);
            G4double stdErr = stdDev / std::sqrt(count);

            G4double meanEnergy = fSumRecoilEnergy.GetValue() / count;
            G4double meanInteractionDepth = fSumInteractionDepth.GetValue() / count;
            G4double meanStoppingDepth = fSumStoppingDepth.GetValue() / count;

            oss << " Particles stopped:      " << count << "\n";
            oss << " Avg Recoil Init Energy: " << meanEnergy / MeV << " MeV\n";
            oss << " Avg Interaction Depth:  " << meanInteractionDepth / um << " um\n";
            oss << " Avg Stopping Depth:     " << meanStoppingDepth / um << " um\n";
            oss << " Average stopping time:  " << mean / ps << " ps\n";
            oss << " Standard Deviation:     " << stdDev / ps << " ps\n";
            oss << " Error Bar (Std Error):  " << stdErr / ps << " ps\n";
        }

        if (nForward > 0) {
            const G4int n = nForward;
            const G4int total = count + nForward + nBack;

            const G4double meanInitE = fSumTInitE.GetValue() / n;
            const G4double meanInitE2 = fSumTInitE2.GetValue() / n;
            G4double varInitE = meanInitE2 - meanInitE * meanInitE;
            if (varInitE < 0) varInitE = 0;
            const G4double meanExitE = fSumTExitE.GetValue() / n;
            const G4double meanExitE2 = fSumTExitE2.GetValue() / n;
            G4double varExitE = meanExitE2 - meanExitE * meanExitE;
            if (varExitE < 0) varExitE = 0;
            const G4double meanInitBeta = fSumTInitBeta.GetValue() / n;
            const G4double meanExitBeta = fSumTExitBeta.GetValue() / n;

            oss << " Particles transmitted:  " << n << " of " << total << " (forward exits)";
            if (nBack > 0) {
                oss << ", " << nBack << " (backward exits)";
            }
            oss << "\n";
            oss << " Avg Recoil Energy before medium: " << meanInitE / MeV << " MeV  (v/c = " << meanInitBeta << ")\n";
            oss << " Avg Recoil Energy after medium:  " << meanExitE / MeV << " MeV  (v/c = " << meanExitBeta << ")\n";
            oss << " Avg Recoil Energy Loss:          " << (meanInitE - meanExitE) / MeV << " MeV\n";
            oss << " StdDev Recoil Energy after:      " << std::sqrt(varExitE) / MeV << " MeV\n";
            oss << " StdErr Recoil Energy after:      " << std::sqrt(varExitE / n) / MeV << " MeV\n";
        }
        oss << "=============================================\n";

        // Print to stdout and save next to the CSV (the current working
        // directory, i.e. the per-macro outputs folder for runner-script runs).
        G4cout << oss.str() << G4endl;
        std::ofstream summaryFile("run_summary.txt");
        if (summaryFile) {
            summaryFile << oss.str();
        } else {
            G4cout << "WARNING: could not write run_summary.txt in the current directory." << G4endl;
        }
    }
}
