#include "SteppingAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4Exception.hh"
#include "RunAction.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"

void SteppingAction::UserSteppingAction(const G4Step* step) {
    if (!step) {
        G4Exception("SteppingAction::UserSteppingAction", "ERR_STEP_01", JustWarning, "Null step encountered.");
        return;
    }

    G4Track* track = step->GetTrack();
    if (!track) {
        G4Exception("SteppingAction::UserSteppingAction", "ERR_STEP_02", JustWarning, "Null track encountered.");
        return;
    }
    
    if (track->GetDefinition()->GetAtomicMass() == 32 && track->GetDefinition()->GetAtomicNumber() == 14) {
        if (track->GetKineticEnergy() <= 0.0) {
            G4double stoppingTime = track->GetGlobalTime();
            G4double initialEnergy = track->GetVertexKineticEnergy();
            
            // Target half-length is 0.5 mm, so front face is at -0.5 mm
            G4double targetFrontFace = -0.5 * mm;
            G4double interactionDepth = track->GetVertexPosition().z() - targetFrontFace;
            G4double stoppingDepth = track->GetPosition().z() - targetFrontFace;
            
            RunAction* runAction = (RunAction*)G4RunManager::GetRunManager()->GetUserRunAction();
            if (runAction) {
                runAction->AddStoppingData(stoppingTime, initialEnergy, interactionDepth, stoppingDepth);
            }

            auto analysisManager = G4AnalysisManager::Instance();
            analysisManager->FillNtupleDColumn(0, stoppingTime / ps);
            analysisManager->FillNtupleDColumn(1, initialEnergy / MeV);
            analysisManager->AddNtupleRow();
            
            track->SetTrackStatus(fStopAndKill);
        }
    }
}
