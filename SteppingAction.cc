#include "SteppingAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4Exception.hh"
#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include <cmath>
#include <algorithm>

// Half thickness of the Target box along z, from the actual geometry (cached).
static G4double TargetHalfThickness() {
    static G4double halfZ = 0.5 * mm;
    static G4bool cached = false;
    if (!cached) {
        G4LogicalVolume* lv = G4LogicalVolumeStore::GetInstance()->GetVolume("Target", false);
        if (lv) {
            G4Box* box = dynamic_cast<G4Box*>(lv->GetSolid());
            if (box) halfZ = box->GetZHalfLength();
        }
        cached = true;
    }
    return halfZ;
}

// Relativistic v/c from kinetic energy and rest mass.
static G4double BetaFromEnergy(G4double kineticEnergy, G4double mass) {
    G4double gamma = 1.0 + kineticEnergy / mass;
    return std::sqrt(std::max(0.0, 1.0 - 1.0 / (gamma * gamma)));
}

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

    // Only the recoil of interest, as configured via /reaction/recoilZ|A
    // (falls back to the 32Si default when not set).
    static G4int recoilA = -1;
    static G4int recoilZ = -1;
    if (recoilA < 0) {
        auto* pga = static_cast<const PrimaryGeneratorAction*>(
            G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
        recoilA = pga ? pga->GetRecoilA() : 32;
        recoilZ = pga ? pga->GetRecoilZ() : 14;
    }
    if (track->GetDefinition()->GetAtomicMass() != recoilA || track->GetDefinition()->GetAtomicNumber() != recoilZ) {
        return;
    }

    const G4VPhysicalVolume* preVol = step->GetPreStepPoint()->GetPhysicalVolume();
    const G4bool inTarget = (preVol != nullptr && preVol->GetName() == "Target");

    // Record the state right after the reaction ("before the medium"): energy and v/c.
    if (inTarget && track->GetCurrentStepNumber() == 1) {
        fRecoilInitEnergy = step->GetPreStepPoint()->GetKineticEnergy();
        fRecoilInitBeta = step->GetPreStepPoint()->GetVelocity() / CLHEP::c_light;
        RunAction* runAction = (RunAction*)G4RunManager::GetRunManager()->GetUserRunAction();
        if (runAction) {
            runAction->AddRecoilCreated(fRecoilInitEnergy);
        }
    }

    const G4StepPoint* postPoint = step->GetPostStepPoint();

    // --- Case 1: recoil comes to rest inside the medium (thick target) ---
    if (track->GetKineticEnergy() <= 0.0) {
        G4double stoppingTime = track->GetGlobalTime();
        G4double initialEnergy = track->GetVertexKineticEnergy();

        G4double targetFrontFace = -TargetHalfThickness();
        G4double interactionDepth = track->GetVertexPosition().z() - targetFrontFace;
        G4double stoppingDepth = track->GetPosition().z() - targetFrontFace;

        RunAction* runAction = (RunAction*)G4RunManager::GetRunManager()->GetUserRunAction();
        if (runAction) {
            runAction->AddStoppingData(stoppingTime, initialEnergy, interactionDepth, stoppingDepth);
        }

        G4double initBeta = (fRecoilInitBeta > 0.0)
            ? fRecoilInitBeta
            : BetaFromEnergy(initialEnergy, track->GetDefinition()->GetPDGMass());

        auto analysisManager = G4AnalysisManager::Instance();
        analysisManager->FillNtupleDColumn(0, stoppingTime / ps);
        analysisManager->FillNtupleDColumn(1, initialEnergy / MeV);
        analysisManager->FillNtupleDColumn(2, initBeta);
        analysisManager->FillNtupleDColumn(3, 0.0);
        analysisManager->FillNtupleDColumn(4, 0.0);
        analysisManager->AddNtupleRow();

        track->SetTrackStatus(fStopAndKill);
        return;
    }

    // --- Case 2: recoil leaves the Target volume through a boundary (thin target) ---
    if (inTarget && postPoint->GetStepStatus() == fGeomBoundary) {
        const G4VPhysicalVolume* postVol = postPoint->GetPhysicalVolume();
        if (!postVol || postVol->GetName() != "Target") {
            G4double exitEnergy = postPoint->GetKineticEnergy();
            G4double exitBeta = postPoint->GetVelocity() / CLHEP::c_light;
            // Sign encodes the exit face: positive = forward (downstream) exit,
            // negative = backward (upstream) exit.
            const G4bool forward = (postPoint->GetPosition().z() > 0.0);
            if (!forward) {
                exitEnergy = -exitEnergy;
                exitBeta = -exitBeta;
            }

            RunAction* runAction = (RunAction*)G4RunManager::GetRunManager()->GetUserRunAction();
            if (runAction) {
                runAction->AddExitData(fRecoilInitEnergy, fRecoilInitBeta, exitEnergy, exitBeta, forward);
            }

            auto analysisManager = G4AnalysisManager::Instance();
            analysisManager->FillNtupleDColumn(0, track->GetGlobalTime() / ps);
            analysisManager->FillNtupleDColumn(1, fRecoilInitEnergy / MeV);
            analysisManager->FillNtupleDColumn(2, fRecoilInitBeta);
            analysisManager->FillNtupleDColumn(3, exitEnergy / MeV);
            analysisManager->FillNtupleDColumn(4, exitBeta);
            analysisManager->AddNtupleRow();

            // Nothing happens in vacuum: stop tracking the transmitted recoil.
            track->SetTrackStatus(fStopAndKill);
        }
        return;
    }

    // --- Case 3: recoil started on the target surface but outside the Target volume
    // (edge case of /reaction/thinTarget true): classify by flight direction.
    if (!inTarget && track->GetCurrentStepNumber() == 1) {
        const G4bool forward = (track->GetMomentumDirection().z() > 0.0);
        G4double exitEnergy = step->GetPreStepPoint()->GetKineticEnergy();
        G4double exitBeta = step->GetPreStepPoint()->GetVelocity() / CLHEP::c_light;
        if (!forward) {
            exitEnergy = -exitEnergy;
            exitBeta = -exitBeta;
        }

        RunAction* runAction = (RunAction*)G4RunManager::GetRunManager()->GetUserRunAction();
        if (runAction) {
            runAction->AddRecoilCreated(std::fabs(exitEnergy));
            runAction->AddExitData(std::fabs(exitEnergy), std::fabs(exitBeta), exitEnergy, exitBeta, forward);
        }

        auto analysisManager = G4AnalysisManager::Instance();
        analysisManager->FillNtupleDColumn(0, track->GetGlobalTime() / ps);
        analysisManager->FillNtupleDColumn(1, std::fabs(exitEnergy) / MeV);
        analysisManager->FillNtupleDColumn(2, std::fabs(exitBeta));
        analysisManager->FillNtupleDColumn(3, exitEnergy / MeV);
        analysisManager->FillNtupleDColumn(4, exitBeta);
        analysisManager->AddNtupleRow();

        track->SetTrackStatus(fStopAndKill);
    }
}
