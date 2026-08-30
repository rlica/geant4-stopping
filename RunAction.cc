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

RunAction::RunAction() : G4UserRunAction() {
    G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
    accumulableManager->RegisterAccumulable(fSumTime);
    accumulableManager->RegisterAccumulable(fSumTime2);
    accumulableManager->RegisterAccumulable(fSumRecoilEnergy);
    accumulableManager->RegisterAccumulable(fSumInteractionDepth);
    accumulableManager->RegisterAccumulable(fSumStoppingDepth);
    accumulableManager->RegisterAccumulable(fStopCount);

    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("csv");
    analysisManager->SetVerboseLevel(0);
    analysisManager->SetFileName("stopping_data");
    analysisManager->CreateNtuple("RecoilData", "Recoil Stopping Data");
    analysisManager->CreateNtupleDColumn("StoppingTime_ps");
    analysisManager->CreateNtupleDColumn("RecoilEnergy_MeV");
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

void RunAction::EndOfRunAction(const G4Run* run) {
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();

    G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
    accumulableManager->Merge();

    G4int count = fStopCount.GetValue();
    if (count > 0 && IsMaster()) {
        G4double sum = fSumTime.GetValue();
        G4double sum2 = fSumTime2.GetValue();

        G4double mean = sum / count;
        G4double variance = (sum2 / count) - (mean * mean);
        if (variance < 0) variance = 0;
        G4double stdDev = std::sqrt(variance);
        G4double stdErr = stdDev / std::sqrt(count);

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
        G4LogicalVolume* targetLV = G4LogicalVolumeStore::GetInstance()->GetVolume("Target", false);
        if (targetLV && targetLV->GetMaterial()) {
            materialName = targetLV->GetMaterial()->GetName();
            density = targetLV->GetMaterial()->GetDensity();
        }

        G4double meanEnergy = fSumRecoilEnergy.GetValue() / count;
        G4double meanInteractionDepth = fSumInteractionDepth.GetValue() / count;
        G4double meanStoppingDepth = fSumStoppingDepth.GetValue() / count;

        G4cout << "\n=============================================" << G4endl;
        G4cout << "               Run Summary                   " << G4endl;
        G4cout << "=============================================" << G4endl;
        G4cout << " Primary Particle:       " << particleName << G4endl;
        G4cout << " Beam Energy:            " << particleEnergy / MeV << " MeV" << G4endl;
        G4cout << " Target Material:        " << materialName << G4endl;
        if (density > 0) {
            G4cout << " Target Density:         " << density / (g/cm3) << " g/cm3" << G4endl;
        }
        G4cout << " Total number of events: " << run->GetNumberOfEvent() << G4endl;
        G4cout << " Particles stopped:      " << count << G4endl;
        G4cout << " Avg Recoil Init Energy: " << meanEnergy / MeV << " MeV" << G4endl;
        G4cout << " Avg Interaction Depth:  " << meanInteractionDepth / um << " um" << G4endl;
        G4cout << " Avg Stopping Depth:     " << meanStoppingDepth / um << " um" << G4endl;
        G4cout << " Average stopping time:  " << mean / ps << " ps" << G4endl;
        G4cout << " Standard Deviation:     " << stdDev / ps << " ps" << G4endl;
        G4cout << " Error Bar (Std Error):  " << stdErr / ps << " ps" << G4endl;
        G4cout << "=============================================\n" << G4endl;
    }
}
